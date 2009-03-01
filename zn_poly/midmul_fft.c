/*
   midmul_fft.c:  middle product by Schonhage FFT
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
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
   The notation below is the same as in mul_fft.c.
   
   Our middle product algorithm conceptually has two pieces:
      (1) reducing the problem from a middle product in R[X] to a middle
          product in S[Z], and
      (2) computing the middle product in S[Z] via transposed truncated
          Fourier transforms, according to the transposition principle
          (see [BLS03]).
   

   (1) Reduction from R[X] to S[Z].

   Suppose we want to compute the middle product of op1[0, len1) and
   op2[0, len2). We first pad op1 by _pad_ zeroes on the left, where _pad_
   satisfies:
      * len2 + pad - 1 is divisible by M/2,
      * 1 <= pad <= M/2.
   Call this op1', which has length len1' = len1 + pad.
   
   Now split op1' into _coeffs1_ chunks of length M/2 and op2 into _coeffs2_
   chunks of length M/2 (both zero-padded on the right to get up to a full
   chunk of length M/2), i.e.
      coeffs1 = ceil(len1' / (M/2))
      coeffs2 = ceil(len2 / (M/2)).
   
   Compute the middle product of the resulting polynomials of length _coeffs1_
   and _coeffs2_ over S. (For this we require that coeffs1 <= K, so that the
   FFTs work.) The result has length coeffs1 - coeffs2 + 1; we must show that
   the result has enough information to reconstruct the middle product of the
   original op1 and op2 in R[X].
   
   The first coefficient of the middle product of op1 and op2 would usually be
   at index len2 - 1 into the full product op1 * op2. Therefore it appears at
   index len2 + pad - 1 into the full product op1' * op2. This appears at
   index (len2 + pad - 1) / (M/2) = coeffs2 of the middle product we performed
   over S, or in other words, the *second* coefficient (since the first
   coefficient would be the one at index coeffs2 - 1). This is good since we
   also need the overlapping data from the previous coefficient.
   
   The last coefficient of the middle product of op1 and op2 would usually be
   at index len1 - 1 into the full product op1 * op2. That's at index
   len1 + pad - 1 into the full product of op1' * op2. This appears at index
   floor((len1 + pad - 1) / (M/2)) <= coeffs1 - 1 of the middle product over
   S, which is good since that's the last one we computed.
   
   
   (2) Middle product in S[Z].
   
   Fix a polynomial F in S[Z] of length n, and let m >= n.
   
   Consider the linear map that sends the length m-n+1 polynomial G to G*F
   (which has length m).
   
   The transpose of this map sends a length m polynomial H to the reversal
   of the middle product of F and H' (of length m-n+1), where H' is the
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

     transposed FFT and IFFT routines

============================================================================ */


/*
   The following routines, i.e.

      zn_pmf_vec_fft_transposed_notrunc_iterative(),
      zn_pmf_vec_fft_transposed_small(),
      zn_pmf_vec_fft_transposed_factor(),
      zn_pmf_vec_fft_transposed(),
      zn_pmf_vec_ifft_transposed_notrunc_iterative(),
      zn_pmf_vec_ifft_transposed_small(),
      zn_pmf_vec_ifft_transposed_factor(),
      zn_pmf_vec_ifft_transposed(),
   
   are *transposed* versions of the corresponding routines from mul_fft.c.
   They are used in the FFT middle product code.
   
   The FFT routines in mul_fft.c are S-linear maps from S^nonzero to S^length.
   The transposed versions are S-linear maps from S^length to S^nonzero whose
   matrices are the transpose of the mul_fft.c versions.
   
   The IFFT routines in mul_fft.c are S-linear maps from S^nonzero to
   S^(length + forward). The transposed versions are S-linear maps from
   S^(length + forward) to S^nonzero whose matrices are the transpose of the
   mul_fft.c versions.
   
   The algorithms are transposed essentially by reversing them, and
   transposing every step of the algorithm; see for example [BLS03] for how
   this is done. We don't have comments on these routines; see the comments
   on the corresponding routines in mul_fft.c.
*/


void zn_pmf_vec_fft_transposed_notrunc_iterative(zn_pmf_vec_t op, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);

   if (op->lgK == 0)
      return;

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M;
   ulong r_last = M >> (op->lgK - 1);
   twist <<= (op->lgK - 1);
   ptrdiff_t half_skip = op->skip;
   ulong* op_end = op->data + (op->skip << op->lgK);
   ulong* op_ptr;
   ulong* start;

   for (; r >= r_last; r >>= 1, half_skip <<= 1, twist >>= 1)
   {
      for (start = op->data, s = twist; s < M; s += r, start += op->skip)
      {
         for (op_ptr = start; op_ptr < op_end; op_ptr += 2*half_skip)
         {
            zn_pmf_rotate(op_ptr + half_skip, M + s);
            zn_pmf_bfly(op_ptr + half_skip, op_ptr, M, mod);
         }
      }
   }
}



void zn_pmf_vec_fft_transposed_small(zn_pmf_vec_t op, ulong length,
                                     ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(length <= op->K);
   ZNP_ASSERT(nonzero <= op->K);

   if (length == op->K && nonzero == op->K)
   {
      zn_pmf_vec_fft_transposed_notrunc_iterative(op, twist);
      return;
   }
   
   if (op->K == 1)
   {
      if (length == 0 && nonzero == 1)
         zn_pmf_zero(op->data, op->M);
      return;
   }

   const zn_mod_struct* mod = op->mod;

   op->lgK--;
   op->K >>= 1;
   
   long i;
   ulong M = op->M;
   ulong U = op->K;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half_skip = skip << op->lgK;

   if (length <= U)
   {
      zn_pmf_vec_fft_transposed_small(op, length, ZNP_MIN(nonzero, U),
                                      twist << 1);

      ulong* op_ptr = op->data;
      for (i = 0; i < (long)(nonzero - U); i++, op_ptr += skip)
         zn_pmf_set(op_ptr + half_skip, op_ptr, M);
   }
   else
   {
      ulong nonzero2 = ZNP_MIN(nonzero, U);

      op->data += half_skip;
      zn_pmf_vec_fft_transposed_small(op, length - U, nonzero2, twist << 1);
      op->data -= half_skip;
      zn_pmf_vec_fft_transposed_small(op, U, nonzero2, twist << 1);

      i = nonzero2 - 1;
      ulong* op_ptr = op->data + i * skip;
      ulong r = M >> op->lgK;
      ulong s = twist + i * r;
      
      for (; i >= (long)(nonzero - U) && i >= 0; i--, op_ptr -= skip, s -= r)
      {
         zn_pmf_rotate(op_ptr + half_skip, s);
         zn_pmf_add(op_ptr, op_ptr + half_skip, M, mod);
      }
      
      for (; i >= 0; i--, op_ptr -= skip, s -= r)
      {
         zn_pmf_rotate(op_ptr + half_skip, M + s);
         zn_pmf_bfly(op_ptr + half_skip, op_ptr, M, mod);
      }
   }

   op->K <<= 1;
   op->lgK++;
}



void zn_pmf_vec_fft_transposed_factor(zn_pmf_vec_t op, unsigned lgT,
                                      ulong length, ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(lgT > 0  &&  lgT < op->lgK);
   ZNP_ASSERT(length <= op->K);
   ZNP_ASSERT(nonzero <= op->K);

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;
   
   ulong length_U = length & (U - 1);
   ulong length_T = length >> lgU;
   
   ulong length_T_ceil = length_T + (length_U > 0);
   
   ulong nonzero_T = nonzero >> lgU;
   ulong nonzero_U = nonzero & (U - 1);
   ulong nonzero_U2 = nonzero_T ? U : nonzero_U;

   ulong r = op->M >> (lgK - 1);
   ulong s, i;
   
   op->K = U;
   op->lgK = lgU;
   twist <<= lgT;

   for (i = 0; i < length_T; i++, op->data += skip_U)
      zn_pmf_vec_fft_transposed(op, U, nonzero_U2, twist);
   
   if (i < T)
      zn_pmf_vec_fft_transposed(op, length_U, nonzero_U2, twist);

   op->data = data;
   op->K = T;
   op->lgK = lgT;
   op->skip = skip_U;
   twist >>= lgT;

   for (i = 0, s = twist; i < nonzero_U; i++, op->data += skip, s += r)
      zn_pmf_vec_fft_transposed(op, length_T_ceil, nonzero_T + 1, s);
      
   if (nonzero_T)
   {
      for (; i < U; i++, op->data += skip, s += r)
         zn_pmf_vec_fft_transposed(op, length_T_ceil, nonzero_T, s);
   }
   
   op->data = data;
   op->skip = skip;
   op->K = K;
   op->lgK = lgK;
}



void zn_pmf_vec_fft_transposed(zn_pmf_vec_t op, ulong length, ulong nonzero,
                               ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(length <= op->K);
   ZNP_ASSERT(nonzero <= op->K);

   if (op->K <= 2  ||  2 * op->K * op->M * sizeof(ulong) <= ZNP_CACHE_SIZE)
   {
      zn_pmf_vec_fft_transposed_small(op, length, nonzero, twist);
   }
   else
   {
      zn_pmf_vec_fft_transposed_factor(op, op->lgK / 2, length,
                                       nonzero, twist);
   }
}



void zn_pmf_vec_ifft_transposed_notrunc_iterative(zn_pmf_vec_t op, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);

   if (op->lgK == 0)
      return;

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M >> (op->lgK - 1);
   ptrdiff_t half_skip = op->skip << (op->lgK - 1);
   ulong* op_end = op->data + (op->skip << op->lgK);
   ulong* op_ptr;
   ulong* start;
   
   for (; r <= M; r <<= 1, half_skip >>= 1, twist <<= 1)
   {
      for (start = op->data, s = twist; s < M; s += r, start += op->skip)
      {
         for (op_ptr = start; op_ptr < op_end; op_ptr += 2*half_skip)
         {
            zn_pmf_bfly(op_ptr, op_ptr + half_skip, M, mod);
            zn_pmf_rotate(op_ptr + half_skip, M - s);
         }
      }
   }
}



void zn_pmf_vec_ifft_transposed_small(zn_pmf_vec_t op, ulong length,
                                      int forward, ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(nonzero <= op->K);
   ZNP_ASSERT(length <= nonzero);
   ZNP_ASSERT(length + forward <= op->K);

   if (length == op->K)
   {
      zn_pmf_vec_ifft_transposed_notrunc_iterative(op, twist);
      return;
   }

   if (op->K == 1)
   {
      if (nonzero && !forward)
         zn_pmf_zero(op->data, op->M);
      return;
   }

   const zn_mod_struct* mod = op->mod;

   op->lgK--;
   op->K >>= 1;
   
   long i;
   ulong M = op->M;
   ulong U = op->K;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half_skip = skip << op->lgK;

   if (length + forward <= U)
   {
      i = ZNP_MIN(nonzero, U);
      long last_zero_fwd_bfly = ZNP_MAX(nonzero - i, length);
      long last_zero_cross_bfly = ZNP_MIN(nonzero - i, length);

      zn_pmf_t op_ptr = op->data;

      for (i = 0; i < last_zero_cross_bfly; i++, op_ptr += skip)
      {
         zn_pmf_set(op_ptr + half_skip, op_ptr, M);
         zn_pmf_rotate(op_ptr + half_skip, M);
         zn_pmf_add(op_ptr, op_ptr, M, mod);
      }
      
      for (; i < length; i++, op_ptr += skip)
         zn_pmf_add(op_ptr, op_ptr, M, mod);
      
      zn_pmf_vec_ifft_transposed_small(op, length, forward,
                                       ZNP_MIN(nonzero, U), twist << 1);
      
      for (; i < last_zero_fwd_bfly; i++, op_ptr += skip)
      {
         zn_pmf_divby2(op_ptr, M, mod);
         zn_pmf_set(op_ptr + half_skip, op_ptr, M);
      }
      
      for (; i < ZNP_MIN(nonzero, U); i++, op_ptr += skip)
         zn_pmf_divby2(op_ptr, M, mod);
   }
   else
   {
      long last_zero_cross_bfly = nonzero - U;
      long last_cross_bfly = length - U;
      ulong r = M >> op->lgK;
      ulong s = twist;
      zn_pmf_t op_ptr = op->data;

      for (i = 0; i < last_cross_bfly; i++, s += r, op_ptr += skip)
      {
         zn_pmf_bfly(op_ptr, op_ptr + half_skip, M, mod);
         zn_pmf_rotate(op_ptr + half_skip, M - s);
      }

      op->data += half_skip;
      zn_pmf_vec_ifft_transposed_small(op, length - U, forward, U, twist << 1);
      op->data -= half_skip;
      
      for (; i < last_zero_cross_bfly; i++, s += r, op_ptr += skip)
      {
         zn_pmf_rotate(op_ptr + half_skip, M + s);
         zn_pmf_sub(op_ptr + half_skip, op_ptr, M, mod);
         zn_pmf_sub(op_ptr, op_ptr + half_skip, M, mod);
      }
      
      for (; i < U; i++, s += r, op_ptr += skip)
      {
         zn_pmf_add(op_ptr, op_ptr, M, mod);
         zn_pmf_rotate(op_ptr + half_skip, s);
         zn_pmf_add(op_ptr, op_ptr + half_skip, M, mod);
      }
      
      zn_pmf_vec_ifft_transposed_notrunc_iterative(op, twist << 1);
   }

   op->K <<= 1;
   op->lgK++;
}



void zn_pmf_vec_ifft_transposed_factor(zn_pmf_vec_t op, unsigned lgT,
                                       ulong length, int forward,
                                       ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(nonzero <= op->K);
   ZNP_ASSERT(length <= nonzero);
   ZNP_ASSERT(length + forward <= op->K);
   ZNP_ASSERT(lgT > 0  &&  lgT < op->lgK);

   if (nonzero == 0)
      return;

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;

   ulong length_U = length & (U - 1);
   ulong length_T = length >> lgU;

   ulong nonzero_U = nonzero & (U - 1);
   ulong nonzero_T = nonzero >> lgU;

   ulong r = op->M >> (lgK - 1);
   ulong s, i;
   ulong twist_T = twist << lgT;

   if (length_U || forward)
   {
      op->lgK = lgT;
      op->K = T;
      op->skip = skip_U;
      
      for (i = 0, op->data = data, s = twist; i < length_U && i < nonzero_U;
           i++, op->data += skip, s += r)
      {
         zn_pmf_vec_ifft_transposed(op, length_T + 1, 0, nonzero_T + 1, s);
      }
      if (nonzero_T)
      {
         for (; i < length_U; i++, op->data += skip, s += r)
            zn_pmf_vec_ifft_transposed(op, length_T + 1, 0, nonzero_T, s);
      }

      op->data = data + length_T * skip_U;
      op->lgK = lgU;
      op->K = U;
      op->skip = skip;
      zn_pmf_vec_ifft_transposed(op, length_U, forward,
                                 nonzero_T ? U : nonzero_U, twist_T);
   }

   op->lgK = lgT;
   op->K = T;
   op->skip = skip_U;

   for (i = length_U, op->data = data + (skip * length_U),
        s = twist + (r * length_U); i < nonzero_U;
        i++, op->data += skip, s += r)
   {
      zn_pmf_vec_ifft_transposed(op, length_T, length_U || forward,
                                 nonzero_T + 1, s);
   }
   if (nonzero_T)
   {
      for (; i < U; i++, op->data += skip, s += r)
         zn_pmf_vec_ifft_transposed(op, length_T, length_U || forward,
                                    nonzero_T, s);
   }

   op->data = data;
   op->skip = skip;
   op->lgK = lgU;
   op->K = U;
   for (i = 0; i < length_T; i++, op->data += skip_U)
      zn_pmf_vec_ifft_transposed(op, U, 0, U, twist_T);

   op->data = data;
   op->lgK = lgK;
   op->K = K;
}



void zn_pmf_vec_ifft_transposed(zn_pmf_vec_t op, ulong length, int forward,
                                ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(nonzero <= op->K);
   ZNP_ASSERT(length <= nonzero);
   ZNP_ASSERT(length + forward <= op->K);

   if (op->K <= 2  ||  2 * op->K * op->M * sizeof(ulong) <= ZNP_CACHE_SIZE)
   {
      zn_pmf_vec_ifft_transposed_small(op, length, forward, nonzero, twist);
   }
   else
   {
      zn_pmf_vec_ifft_transposed_factor(op, op->lgK / 2, length, forward,
                                        nonzero, twist);
   }
}



/* ============================================================================

     main array middle product routines

============================================================================ */


void midmul_fft_params(unsigned* lgK, unsigned* lgM,
                       ulong* coeffs1, ulong* coeffs2, ulong* pad,
                       size_t len1, size_t len2)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   unsigned _lgM;
   size_t _coeffs1;
   ulong M, _pad;

   // increase lgM until all the conditions are satisfied
   for (_lgM = 1; ; _lgM++)
   {
      M = 1UL << _lgM;
      _pad = ((-len2) & (M/2 - 1)) + 1;
      _coeffs1 = CEIL_DIV_2EXP(len1 + _pad, _lgM - 1);
      if (_coeffs1 <= 2*M)
         break;
   }

   *lgM = _lgM;
   *lgK = (_coeffs1 > M) ? (_lgM + 1) : _lgM;
   *pad = _pad;
   *coeffs1 = _coeffs1;
   *coeffs2 = CEIL_DIV_2EXP(len2, _lgM - 1);
}



ulong zn_array_midmul_fft_precomp1_get_fudge(size_t len1, size_t len2,
                                             const zn_mod_t mod)
{
   unsigned lgK, lgM;
   ulong coeffs1, coeffs2, pad;
   midmul_fft_params(&lgK, &lgM, &coeffs1, &coeffs2, &pad, len1, len2);
   
   // need to divide by 2^lgK coming from FFT
   ulong fudge1 = zn_mod_pow2(-lgK, mod);
   // and take into account fudge from pointwise multiplies
   ulong fudge2 = zn_pmf_vec_mul_get_fudge(lgM, 0, mod);
   
   return zn_mod_mul(fudge1, fudge2, mod);
}


ulong zn_array_midmul_fft_get_fudge(size_t len1, size_t len2,
                                    const zn_mod_t mod)
{
   return zn_array_midmul_fft_precomp1_get_fudge(len1, len2, mod);
}


void zn_array_midmul_fft_precomp1_init(
           zn_array_midmul_fft_precomp1_t res, const ulong* op1,
           size_t len1, size_t len2, ulong scale, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   res->len1 = len1;
   res->len2 = len2;
   
   unsigned lgK, lgM;
   
   midmul_fft_params(&lgK, &lgM, &res->coeffs1, &res->coeffs2, &res->pad,
                     len1, len2);

   ulong M = 1UL << lgM;
   ptrdiff_t skip = M + 1;

   // allocate space for transposed IFFT
   zn_pmf_vec_init(res->vec1, lgK, skip, lgM, mod);
   
   // split input, with padding, in reversed order, and apply requested
   // scaling factor
   zn_pmf_vec_reverse(res->vec1, res->coeffs1);
   fft_split(res->vec1, op1, len1, res->pad, scale, 0);
   zn_pmf_vec_reverse(res->vec1, res->coeffs1);
   
   // transposed IFFT first input
   zn_pmf_vec_ifft_transposed(res->vec1, res->coeffs1, 0, res->coeffs1, 0);
}


void zn_array_midmul_fft_precomp1_execute(
            ulong* res, const ulong* op2, ulong scale,
            const zn_array_midmul_fft_precomp1_t precomp)
{
   const zn_pmf_vec_struct* vec1 = precomp->vec1;
   size_t len1 = precomp->len1;
   size_t len2 = precomp->len2;
   ulong coeffs1 = precomp->coeffs1;
   ulong coeffs2 = precomp->coeffs2;

   zn_pmf_vec_t vec2;
   zn_pmf_vec_init(vec2, vec1->lgK, vec1->skip, vec1->lgM, vec1->mod);

   // split and compute FFT of second input (with requested scaling factor)
   fft_split(vec2, op2, len2, 0, scale, 0);
   zn_pmf_vec_fft(vec2, coeffs1, coeffs2, 0);
   
   // pointwise multiply against precomputed transposed IFFT of first input
   zn_pmf_vec_mul(vec2, vec1, vec2, coeffs1, 0);

   // transposed FFT
   ulong coeffs3 = coeffs1 - coeffs2 + 1;
   zn_pmf_vec_fft_transposed(vec2, coeffs1, coeffs3, 0);
   
   // reverse output and combine
   zn_pmf_vec_reverse(vec2, coeffs3);
   fft_combine(res, len1 - len2 + 1, vec2, coeffs3, 1);
   zn_pmf_vec_reverse(vec2, coeffs3);

   zn_pmf_vec_clear(vec2);
}


void zn_array_midmul_fft_precomp1_clear(zn_array_midmul_fft_precomp1_t op)
{
   zn_pmf_vec_clear(op->vec1);
}


void zn_array_midmul_fft(ulong* res, const ulong* op1, size_t len1,
                         const ulong* op2, size_t len2,
                         ulong scale, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   
   // re-use the precomp1 code
   zn_array_midmul_fft_precomp1_t precomp;
   zn_array_midmul_fft_precomp1_init(precomp, op1, len1, len2, scale, mod);
   zn_array_midmul_fft_precomp1_execute(res, op2, 1, precomp);
   zn_array_midmul_fft_precomp1_clear(precomp);
}


// end of file ****************************************************************
