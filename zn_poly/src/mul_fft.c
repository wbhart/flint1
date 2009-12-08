/*
   mul_fft.c:  polynomial multiplication and and middle product via
               Schonhage/Nussbaumer FFT
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


/*
   The multiplication algorithm is essentially that of [Sch77].
   
   We map the problem to S[Z]/(Z^K - 1), where S = R[Y]/(Y^M + 1), M and K are
   powers of two, and K <= 2M (this ensures we have enough roots of unity in
   S for the FFTs). The inputs are split into pieces of size M/2. We need K to
   be large enough that the product can be resolved unambiguously in
   S[Z]/(Z^K - 1), and we want M minimal subject to these conditions (we end
   up with M and K around sqrt(n1 + n2)).


   Our middle product algorithm conceptually has two pieces:
      (1) reducing the problem from a middle product in R[X] to a middle
          product in S[Z], and
      (2) computing the middle product in S[Z] via transposed truncated
          Fourier transforms, according to the transposition principle
          (see [BLS03]).
   

   (1) Reduction from R[X] to S[Z].

   Suppose we want to compute the middle product of op1[0, n1) and op2[0, n2).
   We first pad op1 by p zeroes on the left, where p satisfies:
      * n2 + p - 1 is divisible by M/2,
      * 1 <= p <= M/2.
   Call this op1', which has length n1' = n1 + p.
   
   Now split op1' into m1 chunks of length M/2 and op2 into m2 chunks of
   length M/2 (both zero-padded on the right to get up to a full chunk of
   length M/2), i.e.
      m1 = ceil(n1' / (M/2))
      m2 = ceil(n2 / (M/2)).
   
   Compute the middle product of the resulting polynomials of length m1 and m2
   over S. (For this we require that m1 <= K, so that the FFTs work.) The
   result has length m1 - m2 + 1; we must show that the result has enough
   information to reconstruct the middle product of the original op1 and op2
   in R[X].
   
   The first coefficient of the middle product of op1 and op2 would usually be
   at index n2 - 1 into the full product op1 * op2. Therefore it appears at
   index n2 + p - 1 into the full product op1' * op2. This appears at index
   (n2 + p - 1) / (M/2) = m2 of the middle product we performed over S, or in
   other words, the *second* coefficient (since the first coefficient would be
   the one at index m2 - 1). This is good since we also need the overlapping
   data from the previous coefficient.
   
   The last coefficient of the middle product of op1 and op2 would usually be
   at index n1 - 1 into the full product op1 * op2. That's at index n1 + p - 1
   into the full product of op1' * op2. This appears at index
   floor((n1 + p - 1) / (M/2)) <= m1 - 1 of the middle product over S, which
   is good since that's the last one we computed.
   
   
   (2) Middle product in S[Z].
   
   Fix a polynomial F in S[Z] of length n, and let m >= n.
   
   Consider the linear map that sends the length m - n + 1 polynomial G to 
   G * F (which has length m).
   
   The transpose of this map sends a length m polynomial H to the reversal
   of the middle product of F and H' (of length m - n + 1), where H' is the
   reversal of H.
   
   Therefore our algorithm is:
      * compute usual FFT of F
      * compute transposed IFFT of the reversal of H
      * multiply coefficients pointwise in S
      * compute transposed FFT of product, and reverse the output

*/


#include <stdio.h>
#include "zn_poly_internal.h"


/* ============================================================================

     splitting and combining routines

============================================================================ */


void
fft_split (pmfvec_t res, const ulong* op, size_t n, size_t k, ulong x,
           ulong b)
{
   const zn_mod_struct* mod = res->mod;
   ulong M = res->M;
   pmf_t dest = res->data;

   // handle completely zero blocks from leading zeroes
   for (; k >= M/2; k -= M/2, dest += res->skip)
   {
      dest[0] = b;
      zn_array_zero (dest + 1, M);
   }
   
   // handle block with partially leading zeroes
   if (k)
   {
      dest[0] = b;
      zn_array_zero (dest + 1, k);
      size_t left = M/2 - k;
      
      if (n < left)
      {
         zn_array_scalar_mul_or_copy (dest + 1 + k, op, n, x, mod);
         zn_array_zero (dest + 1 + k + n, M - n - k);
         return;
      }

      zn_array_scalar_mul_or_copy (dest + 1 + k, op, left, x, mod);
      zn_array_zero (dest + 1 + M/2, M/2);
      n -= left;
      op += left;
      dest += res->skip;
   }

   // handle complete blocks of length M/2
   for (; n >= M/2; n -= M/2, op += M/2, dest += res->skip)
   {
      dest[0] = b;
      zn_array_scalar_mul_or_copy (dest + 1, op, M/2, x, mod);
      zn_array_zero (dest + 1 + M/2, M/2);
   }
   
   // last block of fractional length
   if (n)
   {
      dest[0] = b;
      zn_array_scalar_mul_or_copy (dest + 1, op, n, x, mod);
      zn_array_zero (dest + 1 + n, M - n);
   }
}


/*
   If neg == 0, copies op[0, n) into res[0, n).
   If neg == 1, copies the negative of op[0, n) into res[0, n).
*/
#define zn_array_signed_copy \
    ZNP_zn_array_signed_copy
ZNP_INLINE void
zn_array_signed_copy (ulong* res, const ulong* op, ulong n, int neg,
                      const zn_mod_t mod)
{
   if (neg)
      zn_array_neg (res, op, n, mod);
   else
      zn_array_copy (res, op, n);
}


/*
   This routine adds the last M/2 coefficients of op1 to the first M/2
   coefficients of op2, and writes them to res[0, M/2). If n < M/2, it
   only writes the first n coefficients, and ignores the rest.
   
   If op1 is NULL, it is treated as being zero. Ditto for op2.
   
   The main complication in this routine is dealing with the bias fields
   of op1 and op2, so some segments need to be added and some subtracted.
   We still do everything in a single pass.
*/
#define fft_combine_chunk \
    ZNP_fft_combine_chunk
void
fft_combine_chunk (ulong* res, size_t n, pmf_const_t op1,
                   pmf_const_t op2, ulong M, const zn_mod_t mod)
{
   n = ZNP_MIN (n, M/2);

   if (op1 == NULL && op2 == NULL)
   {
      // both inputs are zero; just writes zeroes to the output
      zn_array_zero (res, n);
      return;
   }

   // We want to start reading from the Y^(M/2) coefficient of op1, which
   // is located at index M/2 - bias(op1) mod 2M. But really there are only
   // M coefficients, so we reduce the index mod M, and put neg1 = 1 if we
   // have to negate the coefficients (i.e. they wrapped around
   // negacyclically). If op1 is zero, just put s1 = ULONG_MAX.
   ulong s1 = ULONG_MAX;
   int neg1;
   if (op1)
   {
      s1 = (M/2 - op1[0]) & (2*M - 1);
      neg1 = (s1 >= M);
      if (neg1)
         s1 -= M;
   }

   // Similarly for op2, but we want to start reading from the Y^(M/2)
   // coefficient.
   ulong s2 = ULONG_MAX;
   int neg2;
   if (op2)
   {
      s2 = (-op2[0]) & (2*M - 1);
      neg2 = (s2 >= M);
      if (neg2)
         s2 -= M;
   }
   
   // Swap the inputs so that s1 <= s2.
   if (s1 > s2)
   {
      pmf_const_t op_temp = op1; op1 = op2; op2 = op_temp;
      ulong s_temp = s1; s1 = s2; s2 = s_temp;
      int neg_temp = neg1; neg1 = neg2; neg2 = neg_temp;
   }
   
   // advance beyond bias fields
   op1++;
   op2++;
   
   if (s2 == ULONG_MAX)
   {
      // One of the inputs is zero; may assume it's op2. We only need to
      // work with op1.

      // op1 looks like this:
      //
      //      0      s1              M
      // op1: BBBBBBBBAAAAAAAAAAAAAAAA
      //
      // The A parts need to be copied with the same sign; the B parts need
      // to have the sign flipped.

      if (n <= M - s1)
         // Only need part of AAAA up to n.
         zn_array_signed_copy (res, op1 + s1, n, neg1, mod);
      else
      {
         // Copy AAAAA
         zn_array_signed_copy (res, op1 + s1, M - s1, neg1, mod);
         // Negate BBBBB
         zn_array_signed_copy (res + M - s1, op1, n - M + s1, !neg1, mod);
      }

      return;
   }

   // Neither op1 nor op2 are zero.
   
   // The picture looks like this:
   //
   //      0      s1                             M
   // op1: CCCCCCCCAAAAAAAAAAAAABBBBBBBBBBBBBBBBBB
   //
   //                              s2           M
   // op2: BBBBBBBBBBBBBBBBBBCCCCCCCCAAAAAAAAAAAAA

   // Combine the portions marked AAAA
   // (bail out if we reach n)
   if (n <= M - s2)
   {
      zn_skip_array_signed_add (res, 1, n, op2 + s2, neg2,
                                op1 + s1, neg1, mod);
      return;
   }

   res = zn_skip_array_signed_add (res, 1, M - s2, op2 + s2, neg2,
                                   op1 + s1, neg1, mod);
   n -= (M - s2);

   // Combine the portions marked BBBB
   // (bail out if we reach n)
   if (n <= s2 - s1)
   {
      zn_skip_array_signed_add (res, 1, n, op2, !neg2,
                                op1 + s1 + M - s2, neg1, mod);
      return;
   }
   
   res = zn_skip_array_signed_add (res, 1, s2 - s1, op2, !neg2,
                                   op1 + s1 + M - s2, neg1, mod);
   n -= (s2 - s1);

   // Combine the portions marked CCCC
   zn_skip_array_signed_add (res, 1, (n >= s1) ? s1 : n,
                             op2 + s2 - s1, !neg2, op1, !neg1, mod);
}



void
fft_combine (ulong* res, size_t n, const pmfvec_t op, ulong z, int skip_first)
{
   if (z == 0)
   {
      // zero it out
      zn_array_zero (res, n);
      return;
   }
   
   if (!skip_first)
   {
      // Copy the relevant part of the first coefficient
      size_t k = ZNP_MIN (n, op->M/2);
      fft_combine_chunk (res, k, NULL, op->data, op->M, op->mod);
      res += k;
      n -= k;
   }
   
   // In the loop below, ptr1 = (i-1)-th coefficient, ptr2 = i-th coefficient
   pmf_const_t ptr1 = op->data;
   pmf_const_t ptr2 = op->data + op->skip;

   ulong i;
   for (i = 1; i < z && n >= op->M/2;
        i++, n -= op->M/2, res += op->M/2, ptr1 += op->skip, ptr2 += op->skip)
   {
      // Add first half of i-th coefficient to second half of (i-1)-th
      // coefficient
      fft_combine_chunk (res, n, ptr1, ptr2, op->M, op->mod);
   }

   if (i < z)
   {
      // Ran out of output space before getting to last pmf_t.
      // Do the same add operation as above, but stop when the buffer is full.
      fft_combine_chunk (res, n, ptr1, ptr2, op->M, op->mod);
      return;
   }

   // Arrived at last coefficient, still haven't exhausted output buffer.
   // Copy second half of last coefficient, and zero-pad to the end.
   fft_combine_chunk (res, n, ptr1, NULL, op->M, op->mod);
   if (n > op->M/2)
      zn_array_zero (res + op->M/2, n - op->M/2);
}



/* ============================================================================

     multiplication routine

============================================================================ */


void
mul_fft_params (unsigned* lgK, unsigned* lgM, ulong* m1, ulong* m2,
                size_t n1, size_t n2)
{
   unsigned _lgM;
   size_t _m1, _m2, _m3;
   ulong M;

   // increase lgM until all the conditions are satisfied
   for (_lgM = 1; ; _lgM++)
   {
      _m1 = CEIL_DIV_2EXP (n1, _lgM - 1);      // = ceil(n1 / (M/2))
      _m2 = CEIL_DIV_2EXP (n2, _lgM - 1);      // = ceil(n2 / (M/2))
      _m3 = _m1 + _m2 - 1;

      M = 1UL << _lgM;
      if (_m3 <= 2 * M)
         break;
   }

   *lgM = _lgM;
   *lgK = (_m3 > M) ? (_lgM + 1) : _lgM;
   *m1 = _m1;
   *m2 = _m2;
}



ulong
zn_array_mul_fft_fudge (size_t n1, size_t n2, int sqr, const zn_mod_t mod)
{
   unsigned lgK, lgM;
   ulong m1, m2;
   mul_fft_params (&lgK, &lgM, &m1, &m2, n1, n2);
   
   // need to divide by 2^lgK coming from FFT
   ulong fudge1 = zn_mod_pow2 (-lgK, mod);
   // and take into account fudge from pointwise multiplies
   ulong fudge2 = pmfvec_mul_fudge (lgM, sqr, mod);
   
   return zn_mod_mul (fudge1, fudge2, mod);
}



void zn_array_mul_fft (ulong* res,
                       const ulong* op1, size_t n1,
                       const ulong* op2, size_t n2,
                       ulong x, const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   unsigned lgK, lgM;
   
   // number of pmf_t coefficients for each input poly
   ulong m1, m2;

   // figure out how big the transform needs to be
   mul_fft_params (&lgK, &lgM, &m1, &m2, n1, n2);
   
   // number of pmf_t coefficients for output poly
   ulong m3 = m1 + m2 - 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ptrdiff_t skip = M + 1;
   
   pmfvec_t vec1, vec2;
   
   int sqr = (op1 == op2  &&  n1 == n2);

   if (!sqr)
   {
      // multiplying two distinct inputs

      // split inputs into pmf_t's and perform FFTs
      pmfvec_init (vec1, lgK, skip, lgM, mod);
      fft_split (vec1, op1, n1, 0, 1, 0);
      pmfvec_fft (vec1, m3, m1, 0);

      // note: we apply the fudge factor here, because the second input is
      // shorter than both the first input and the output :-)
      pmfvec_init (vec2, lgK, skip, lgM, mod);
      fft_split (vec2, op2, n2, 0, x, 0);
      pmfvec_fft (vec2, m3, m2, 0);

      // pointwise multiplication
      pmfvec_mul (vec1, vec1, vec2, m3, 1);

      pmfvec_clear (vec2);
   }
   else
   {
      // squaring a single input
   
      // split input into pmf_t's and perform FFTs
      pmfvec_init (vec1, lgK, skip, lgM, mod);
      fft_split (vec1, op1, n1, 0, 1, 0);
      pmfvec_fft (vec1, m3, m1, 0);

      // pointwise multiplication
      pmfvec_mul (vec1, vec1, vec1, m3, 1);
   }

   // inverse FFT, and write output
   pmfvec_ifft (vec1, m3, 0, m3, 0);
   size_t n3 = n1 + n2 - 1;
   fft_combine (res, n3, vec1, m3, 0);

   pmfvec_clear (vec1);
   
   // if we're squaring, then we haven't applied the fudge factor yet,
   // so do it now
   if (sqr)
      zn_array_scalar_mul_or_copy (res, res, n3, x, mod);
}



/* ============================================================================

     middle product routines

============================================================================ */


void
mulmid_fft_params (unsigned* lgK, unsigned* lgM, ulong* m1, ulong* m2,
                   ulong* p, size_t n1, size_t n2)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   unsigned _lgM;
   size_t _m1;
   ulong M, _p;

   // increase lgM until all the conditions are satisfied
   for (_lgM = 1; ; _lgM++)
   {
      M = 1UL << _lgM;
      _p = ((-n2) & (M/2 - 1)) + 1;
      _m1 = CEIL_DIV_2EXP (n1 + _p, _lgM - 1);
      if (_m1 <= 2 * M)
         break;
   }

   *lgM = _lgM;
   *lgK = (_m1 > M) ? (_lgM + 1) : _lgM;
   *p = _p;
   *m1 = _m1;
   *m2 = CEIL_DIV_2EXP (n2, _lgM - 1);
}



ulong
zn_array_mulmid_fft_precomp1_fudge (size_t n1, size_t n2, const zn_mod_t mod)
{
   unsigned lgK, lgM;
   ulong m1, m2, p;
   mulmid_fft_params (&lgK, &lgM, &m1, &m2, &p, n1, n2);
   
   // need to divide by 2^lgK coming from FFT
   ulong fudge1 = zn_mod_pow2 (-lgK, mod);
   // and take into account fudge from pointwise multiplies
   ulong fudge2 = pmfvec_mul_fudge (lgM, 0, mod);
   
   return zn_mod_mul (fudge1, fudge2, mod);
}


ulong
zn_array_mulmid_fft_fudge (size_t n1, size_t n2, const zn_mod_t mod)
{
   return zn_array_mulmid_fft_precomp1_fudge (n1, n2, mod);
}


void
zn_array_mulmid_fft_precomp1_init (zn_array_mulmid_fft_precomp1_t res,
                                   const ulong* op1, size_t n1, size_t n2,
                                   ulong x, const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   res->n1 = n1;
   res->n2 = n2;
   
   unsigned lgK, lgM;
   
   mulmid_fft_params (&lgK, &lgM, &res->m1, &res->m2, &res->p, n1, n2);

   ulong M = 1UL << lgM;
   ptrdiff_t skip = M + 1;

   // allocate space for transposed IFFT
   pmfvec_init (res->vec1, lgK, skip, lgM, mod);
   
   // split input, with padding, in reversed order, and apply requested
   // scaling factor
   pmfvec_reverse (res->vec1, res->m1);
   fft_split (res->vec1, op1, n1, res->p, x, 0);
   pmfvec_reverse (res->vec1, res->m1);
   
   // transposed IFFT first input
   pmfvec_tpifft (res->vec1, res->m1, 0, res->m1, 0);
}


void
zn_array_mulmid_fft_precomp1_execute
                    (ulong* res, const ulong* op2, ulong x,
                     const zn_array_mulmid_fft_precomp1_t precomp)
{
   const pmfvec_struct* vec1 = precomp->vec1;
   size_t n1 = precomp->n1;
   size_t n2 = precomp->n2;
   ulong m1 = precomp->m1;
   ulong m2 = precomp->m2;

   pmfvec_t vec2;
   pmfvec_init (vec2, vec1->lgK, vec1->skip, vec1->lgM, vec1->mod);

   // split and compute FFT of second input (with requested scaling factor)
   fft_split (vec2, op2, n2, 0, x, 0);
   pmfvec_fft (vec2, m1, m2, 0);
   
   // pointwise multiply against precomputed transposed IFFT of first input
   pmfvec_mul (vec2, vec1, vec2, m1, 0);

   // transposed FFT
   ulong m3 = m1 - m2 + 1;
   pmfvec_tpfft (vec2, m1, m3, 0);
   
   // reverse output and combine
   pmfvec_reverse (vec2, m3);
   fft_combine (res, n1 - n2 + 1, vec2, m3, 1);
   pmfvec_reverse (vec2, m3);

   pmfvec_clear (vec2);
}


void
zn_array_mulmid_fft_precomp1_clear (zn_array_mulmid_fft_precomp1_t op)
{
   pmfvec_clear (op->vec1);
}


void
zn_array_mulmid_fft (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     ulong x, const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   
   // re-use the precomp1 code
   zn_array_mulmid_fft_precomp1_t precomp;
   zn_array_mulmid_fft_precomp1_init (precomp, op1, n1, n2, x, mod);
   zn_array_mulmid_fft_precomp1_execute (res, op2, 1, precomp);
   zn_array_mulmid_fft_precomp1_clear (precomp);
}


// end of file ****************************************************************
