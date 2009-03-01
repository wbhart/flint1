/*
   mul_fft.c:  multiplication by Schonhage FFT
   
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
   The algorithm is essentially that of [Sch77].
   
   Basically we map the problem to S[Z]/(Z^K - 1), where S = R[Y]/(Y^M + 1),
   M and K are powers of two, and K <= 2M (this ensures we have enough roots
   of unity in S for the FFTs). The inputs are split into pieces of size M/2.
   We need K to be large enough that the product can be resolved unambiguously
   in S[Z]/(Z^K - 1), and we want M minimal subject to these conditions
   (we end up with M and K around sqrt(len1 + len2)).
*/


#include <stdio.h>
#include "zn_poly_internal.h"


/* ============================================================================

     splitting and combining routines

============================================================================ */


/*
   A convenience function used in fft_split below. It behaves just like
   zn_array_scalar_mul, except it uses the obvious optimisation if x == 1.
*/
#define zn_array_scalar_mul_or_copy \
    ZNP_zn_array_scalar_mul_or_copy
void zn_array_scalar_mul_or_copy(ulong* res, const ulong* op, size_t len,
                                 ulong x, const zn_mod_t mod)
{
   if (x == 1)
   {
      if (res == op)
         return;
      zn_array_copy(res, op, len);
   }
   else
      zn_array_scalar_mul(res, op, len, x, mod);
}


void fft_split(zn_pmf_vec_t res, const ulong* op, size_t len, size_t lead,
               ulong scale, ulong bias)
{
   const zn_mod_struct* mod = res->mod;
   ulong M = res->M;
   zn_pmf_t dest = res->data;

   // handle completely zero blocks from leading zeroes
   for (; lead >= M/2; lead -= M/2, dest += res->skip)
   {
      dest[0] = bias;
      zn_array_zero(dest + 1, M);
   }
   
   // handle block with partially leading zeroes
   if (lead)
   {
      dest[0] = bias;
      zn_array_zero(dest + 1, lead);
      size_t left = M/2 - lead;
      
      if (len < left)
      {
         zn_array_scalar_mul_or_copy(dest + 1 + lead, op, len, scale, mod);
         zn_array_zero(dest + 1 + lead + len, M - len - lead);
         return;
      }

      zn_array_scalar_mul_or_copy(dest + 1 + lead, op, left, scale, mod);
      zn_array_zero(dest + 1 + M/2, M/2);
      len -= left;
      op += left;
      dest += res->skip;
   }

   // handle complete blocks of length M/2
   for (; len >= M/2; len -= M/2, op += M/2, dest += res->skip)
   {
      dest[0] = bias;
      zn_array_scalar_mul_or_copy(dest + 1, op, M/2, scale, mod);
      zn_array_zero(dest + 1 + M/2, M/2);
   }
   
   // last block of fractional length
   if (len)
   {
      dest[0] = bias;
      zn_array_scalar_mul_or_copy(dest + 1, op, len, scale, mod);
      zn_array_zero(dest + 1 + len, M - len);
   }
}


/*
   If neg == 0, copies op[0, len) into res[0, len).
   If neg == 1, copies the negative of op[0, len) into res[0, len).
*/
#define zn_array_signed_copy \
    ZNP_zn_array_signed_copy
ZNP_INLINE
void zn_array_signed_copy(ulong* res, ulong len,
                          const ulong* op, int neg, const zn_mod_t mod)
{
   if (neg)
      zn_array_neg(res, op, len, mod);
   else
      zn_array_copy(res, op, len);
}


/*
   This routine adds the last M/2 coefficients of op1 to the first M/2
   coefficients of op2, and writes them to res[0, M/2). If len < M/2, it
   only writes the first _len_ coefficients, and ignores the rest.
   
   If op1 is NULL, it is treated as being zero. Ditto for op2.
   
   The main complication in this routine is dealing with the bias fields
   of op1 and op2, so some segments need to be added and some subtracted.
   We still do everything in a single pass.
*/
#define fft_combine_chunk \
    ZNP_fft_combine_chunk
void fft_combine_chunk(ulong* res, size_t len, zn_pmf_const_t op1,
                       zn_pmf_const_t op2, ulong M, const zn_mod_t mod)
{
   len = ZNP_MIN(len, M/2);

   if (op1 == NULL && op2 == NULL)
   {
      // both inputs are zero; just writes zeroes to the output
      zn_array_zero(res, len);
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
      const ulong* op_temp = op1; op1 = op2; op2 = op_temp;
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

      if (len <= M - s1)
         // Only need part of AAAA up to len.
         zn_array_signed_copy(res, len, op1 + s1, neg1, mod);
      else
      {
         // Copy AAAAA
         zn_array_signed_copy(res, M - s1, op1 + s1, neg1, mod);
         // Negate BBBBB
         zn_array_signed_copy(res + M - s1, len - M + s1, op1, !neg1, mod);
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
   // (bail out if we reach len)
   if (len <= M - s2)
   {
      zn_skip_array_signed_add(res, 1, len, op2 + s2, neg2, 
                               op1 + s1, neg1, mod);
      return;
   }

   res = zn_skip_array_signed_add(res, 1, M - s2, op2 + s2, neg2,
                                  op1 + s1, neg1, mod);
   len -= (M - s2);

   // Combine the portions marked BBBB
   // (bail out if we reach len)
   if (len <= s2 - s1)
   {
      zn_skip_array_signed_add(res, 1, len, op2, !neg2,
                               op1 + s1 + M - s2, neg1, mod);
      return;
   }
   
   res = zn_skip_array_signed_add(res, 1, s2 - s1, op2, !neg2,
                                  op1 + s1 + M - s2, neg1, mod);
   len -= (s2 - s1);

   // Combine the portions marked CCCC
   zn_skip_array_signed_add(res, 1, (len >= s1) ? s1 : len,
                            op2 + s2 - s1, !neg2, op1, !neg1, mod);
}



void fft_combine(ulong* res, size_t len, const zn_pmf_vec_t op, ulong nonzero,
                 int skip_first)
{
   if (nonzero == 0)
   {
      // zero it out
      zn_array_zero(res, len);
      return;
   }
   
   if (!skip_first)
   {
      // Copy the relevant part of the first coefficient
      size_t amount = ZNP_MIN(len, op->M/2);
      fft_combine_chunk(res, amount, NULL, op->data, op->M, op->mod);
      res += amount;
      len -= amount;
   }
   
   // In the loop below, ptr1 = (i-1)-th coefficient, ptr2 = i-th coefficient
   const ulong* ptr1 = op->data;
   const ulong* ptr2 = op->data + op->skip;

   ulong i;
   for (i = 1; i < nonzero && len >= op->M/2; i++, len -= op->M/2,
        res += op->M/2, ptr1 += op->skip, ptr2 += op->skip)
   {
      // Add first half of i-th coefficient to second half of (i-1)-th
      // coefficient
      fft_combine_chunk(res, len, ptr1, ptr2, op->M, op->mod);
   }

   if (i < nonzero)
   {
      // Ran out of output space before getting to last zn_pmf_t.
      // Do the same add operation as above, but stop when the buffer is full.
      fft_combine_chunk(res, len, ptr1, ptr2, op->M, op->mod);
      return;
   }

   // Arrived at last coefficient, still haven't exhausted output buffer.
   // Copy second half of last coefficient, and zero-pad to the end.
   fft_combine_chunk(res, len, ptr1, NULL, op->M, op->mod);
   if (len > op->M/2)
      zn_array_zero(res + op->M/2, len - op->M/2);
}


/* ============================================================================

     FFT routines

============================================================================ */


/*
   The various FFT functions below operate on a zn_pmf_vec_t, and compute
   inplace:
   
      b_k = u^{k'} \sum_{i=0}^{K-1} w^{i k'} a_i,

   where w = Y^(2M/k) (the standard K-th root of unity), and u = Y^twist.
   Must have 0 <= twist < 2M/K. The notation k' indicates the bit-reversal
   of k of length lgK.
   
   Several versions have a _length_ parameter; they only compute the first
   _length_ coefficients. (The remaining coefficients are used in intermediate
   computations, and contain junk at the end.) If there is no _length_
   parameter, all K coefficients are computed.
   
   Several versions have a _nonzero_ parameter; they assume that the input
   coefficients are zero from index _nonzero_ and beyond. They never read from
   those coefficients. If there is no _nonzero_ parameter, all the inputs
   are used.
   
   There are four versions of the FFT routine:
   
      * zn_pmf_vec_fft(): main entry point, delegates to one of the other
        routines based on the size of the transform.
        
      * zn_pmf_vec_fft_notrunc_iterative(): low-overhead iterative FFT,
        without any truncation.
        
      * zn_pmf_vec_fft_small(): handles the top layer of butterflies, and then
        recurses into the two halves. This is intended for fairly small
        transforms, where locality is not a big issue. The algorithm
        implemented here is essentially van der Hoeven's "truncated Fourier
        transform" [vdH04], [vdH05].
        
      * zn_pmf_vec_fft_factor(): intended for large transforms, where locality
        is an issue. It factors the FFT into U transforms of length T followed
        by T transforms of length U, where K = T * U. This is done recursively
        until we hope we're in L1 cache.
        
        The algorithm is straightforward, but I believe it to be new. It is
        simultaneously a generalisation of van der Hoeven's truncated
        transform and Bailey's FFT algorithm [Bai89]. (I used something
        similar in the ZmodF_poly module in FLINT.)

        (Note: zn_pmf_vec_fft_small() is essentially equivalent to
        zn_pmf_vec_fft_factor() with lgT = 1.)
        
        The factoring version is called recursively until we hope we're in L1
        cache, at which point we switch to zn_pmf_vec_fft_small().


   NOTE: our approach to improving cache performance is certainly not ideal.
   The main problem is that we can get address conflicts, especially since
   everything gets spread out by powers of two. Mitigating factors:
   associativity in the cache; the extra bias word scrambles the addresses
   somewhat; when the transforms gets large, so do the coefficients, so we
   don't expect to fit that many in cache anyway.

   NOTE: these functions are not thread-safe. Apart from modifying the input
   inplace, they also temporarily modify the zn_pmf_vec_t structs themselves.

*/


void zn_pmf_vec_fft_notrunc_iterative(zn_pmf_vec_t op, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   
   if (op->lgK == 0)
      return;

   // just plain butterfly loop

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M >> (op->lgK - 1);    // 2M/K = index for K-th root of unity
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
            zn_pmf_rotate(op_ptr + half_skip, M + s);
         }
      }
   }
}



void zn_pmf_vec_fft_small(zn_pmf_vec_t op, ulong length, ulong nonzero,
                          ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(length <= op->K);
   ZNP_ASSERT(nonzero <= op->K);

   if (length == 0)
      return;

   if (length == op->K && nonzero == op->K)
   {
      // No truncation requested; use iterative version
      zn_pmf_vec_fft_notrunc_iterative(op, twist);
      return;
   }
   
   if (op->K == 1)
   {
      // Base case length 1 transform
      if (length == 1 && nonzero == 0)
         zn_pmf_zero(op->data, op->M);
      return;
   }

   const zn_mod_struct* mod = op->mod;

   // We treat the input as two rows and U columns, in row-major order.
   
   // descend to first row (first half of op)
   op->lgK--;
   op->K >>= 1;
   
   long i;
   ulong M = op->M;
   ulong U = op->K;
   ulong* op_ptr = op->data;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half_skip = skip << op->lgK;
   ulong nonzero2 = ZNP_MIN(nonzero, U);

   if (length <= U)
   {
      // Only need the first output of the first layer of butterflies.
      for (i = 0; i < (long)(nonzero - U); i++, op_ptr += skip)
         zn_pmf_add(op_ptr, op_ptr + half_skip, M, mod);
      
      // Recurse into top row
      zn_pmf_vec_fft_small(op, length, nonzero2, twist << 1);
   }
   else
   {
      // Need both outputs from the first layer of butterflies.
      ulong s = twist;
      ulong r = M >> op->lgK;

      for (i = 0; i < (long)(nonzero - U); i++, op_ptr += skip, s += r)
      {
         zn_pmf_bfly(op_ptr, op_ptr + half_skip, M, mod);
         zn_pmf_rotate(op_ptr + half_skip, M + s);
      }

      // Butterflies where second input is zero
      for (; i < nonzero2; i++, op_ptr += skip, s += r)
      {
         zn_pmf_set(op_ptr + half_skip, op_ptr, M);
         zn_pmf_rotate(op_ptr + half_skip, s);
      }
      
      // Recurse into top row...
      zn_pmf_vec_fft_small(op, U, nonzero2, twist << 1);

      // ... and recurse into bottom row
      op->data += half_skip;
      zn_pmf_vec_fft_small(op, length - U, nonzero2, twist << 1);
      op->data -= half_skip;
   }

   // pop back to whole transform
   op->K <<= 1;
   op->lgK++;
}



/*
   As described above, this splits the length K transform into T rows by
   U columns, where K = U * T, T = 2^lgT, U = 2^lgU.
   
   Must have 0 < lgT < lgK.
*/
void zn_pmf_vec_fft_factor(zn_pmf_vec_t op, unsigned lgT,
                           ulong length, ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(lgT > 0  &&  lgT < op->lgK);
   ZNP_ASSERT(length <= op->K);
   ZNP_ASSERT(nonzero <= op->K);
   
   if (length == 0)
      return;

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;
   
   // We need _length_ output coefficients, starting from the top-left,
   // in row-major order.

   // Write length = U * length_T + length_U, where 0 <= length_U < U
   ulong length_U = length & (U - 1);
   ulong length_T = length >> lgU;
   
   // length_T_ceil = number of rows of output, including the last partial row
   ulong length_T_ceil = length_T + (length_U > 0);
   
   // Write nonzero = U * nonzero_T + nonzero_U, where 0 <= nonzero_U < U
   ulong nonzero_T = nonzero >> lgU;
   ulong nonzero_U = nonzero & (U - 1);

   ulong r = op->M >> (lgK - 1);     // 2M/K = index for K-th root of unity
   ulong s, i;
   
   // --------------- FFTs along columns

   op->K = T;
   op->lgK = lgT;
   op->skip = skip_U;

   // First handle the columns with nonzero_T + 1 input coefficients.
   for (i = 0, s = twist; i < nonzero_U; i++, op->data += skip, s += r)
      zn_pmf_vec_fft(op, length_T_ceil, nonzero_T + 1, s);
      
   if (nonzero_T)
   {
      // Handle the remaining columns, which only have nonzero_T input
      // coefficients.
      for (; i < U; i++, op->data += skip, s += r)
         zn_pmf_vec_fft(op, length_T_ceil, nonzero_T, s);

      // Update nonzero_U to reflect the fact that the last row doesn't have
      // any zeroes any more.
      nonzero_U = U;
   }
   
   // --------------- FFTs along rows

   op->data = data;
   op->K = U;
   op->lgK = lgU;
   op->skip = skip;
   twist <<= lgT;

   // Handle the first length_T rows.
   for (i = 0; i < length_T; i++, op->data += skip_U)
      zn_pmf_vec_fft(op, U, nonzero_U, twist);
   
   // For the last row, we only need the first length_U outputs:
   if (i < T)
      zn_pmf_vec_fft(op, length_U, nonzero_U, twist);
   
   // --------------- restore parameters

   op->data = data;
   op->K = K;
   op->lgK = lgK;
}



void zn_pmf_vec_fft(zn_pmf_vec_t op, ulong length, ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(length <= op->K);
   ZNP_ASSERT(nonzero <= op->K);
   
   if (op->K <= 2  ||  2 * op->K * op->M * sizeof(ulong) <= ZNP_CACHE_SIZE)
   {
      // FFT is pretty small; use small version
      zn_pmf_vec_fft_small(op, length, nonzero, twist);
   }
   else
   {
      // FFT is relatively big; use factoring algorithm instead
      zn_pmf_vec_fft_factor(op, op->lgK / 2, length, nonzero, twist);
   }
}



/* ============================================================================

     inverse FFT routines

============================================================================ */

/*
   The IFFTs are a little more complicated than the FFTs, mostly because of
   the truncation.
   
   Let a_i and b_k be as described above for the FFTs. The IFFT functions
   take as input the array

      b_{0'}, b_{1'}, ..., b_{(m-1)'}, K*a_m, K*a_{m+1}, ..., K*a_{K-1},
      
   where m = _length_. If the _forward_ flag is zero, the output of the
   IFFT is:
   
      K*a_0, K*a_1, ..., K*a_{m-1}, ... (then K-m junk coefficients)
      
   If the _forward_ flag is set, the output is

      K*a_0, K*a_1, ..., K*a_{m-1}, b_m, ... (then K-m-1 junk coefficients)
   
   i.e. it also computes one coefficient of the *forward* FFT.
   
   If there is no _length_ parameter, then we take m = k. In this case we
   require that forward = 0, and the routine becomes the (non-truncated)
   IFFT as usually understood, with inputs in bit-reversed order and outputs
   in normal order.

   Several versions have a _nonzero_ parameter; they assume that the input
   coefficients are zero from index _nonzero_ and beyond. They never read from
   those coefficients. Must have nonzero >= length. If there is no _nonzero_
   parameter, all the inputs are used.
   
   There are four versions of the IFFT routine:
   
      * zn_pmf_vec_ifft(): main entry point, delegates to one of the other
        routines based on the size of the transform.
        
      * zn_pmf_vec_ifft_notrunc_iterative(): low-overhead iterative IFFT,
        without any truncation.
        
      * zn_pmf_vec_ifft_small(): recurses into the two halves, and handles the
        top layer of butterflies. This is intended for fairly small transforms,
        where locality is not a big issue. The algorithm implemented here is
        essentially van der Hoeven's "truncated inverse Fourier transform".
        
      * zn_pmf_vec_ifft_factor(): intended for large transforms, where locality
        is an issue. It factors the FFT into U transforms of length T and
        T transforms of length U, where K = T * U. This is done recursively
        until we hope we're in L1 cache.
        
        The algorithm is not as simple as the FFT version; it is necessary
        to alternate between "row" and "column" transforms in a slightly
        complicated way. I believe the algorithm to be new.

        (Note: zn_pmf_vec_ifft_small() is essentially equivalent to
        zn_pmf_vec_ifft_factor() with lgT = 1.)
        

   NOTE: the remarks made for the FFTs concerning cache-friendliness and
   thread-safety apply equally to these functions.

*/

void zn_pmf_vec_ifft_notrunc_iterative(zn_pmf_vec_t op, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);

   if (op->lgK == 0)
      return;

   // just plain butterfly loop

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M;
   ulong r_last = M >> (op->lgK - 1);    // 2M/K = index for K-th root of unity
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
            zn_pmf_rotate(op_ptr + half_skip, M - s);
            zn_pmf_bfly(op_ptr + half_skip, op_ptr, M, mod);
         }
      }
   }
}



void zn_pmf_vec_ifft_small(zn_pmf_vec_t op, ulong length, int forward,
                           ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(nonzero <= op->K);
   ZNP_ASSERT(length <= nonzero);
   ZNP_ASSERT(length + forward <= op->K);

   if (length == op->K)
   {
      // No truncation requested; use iterative version
      zn_pmf_vec_ifft_notrunc_iterative(op, twist);
      return;
   }

   if (op->K == 1)
   {
      // Base case length 1 transform; may assume that length == 0
      if (forward && nonzero == 0)
         zn_pmf_zero(op->data, op->M);
      return;
   }
   
   const zn_mod_struct* mod = op->mod;

   // We treat the input as two rows and U columns, in row-major order.
   
   // descend to first row (first half of op)
   op->K >>= 1;
   op->lgK--;

   ulong M = op->M;
   ulong U = op->K;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half_skip = skip << op->lgK;

   // symbols in the following diagrams:
   // A = fully untransformed coefficient (one of the a_i)
   // B = intermediate coefficient
   // C = fully transformed coefficient (one of the b_k)
   // a, b, c = same as three above, but implied zero
   // ? = garbage that we don't care about
   // * = the "forward" C coefficient, or "?" if no forward coefficient
   //     requested
   
   // The horizontal transforms convert between B and C.
   // The vertical butterflies convert between A and B.
   
   if (length + forward <= U)
   {
      // The input could look like one of the following:
      // CCCCAAAA      CCCCAAAA      CCCCAAaa      CCCCaaaa
      // AAAAAAaa  or  AAaaaaaa  or  aaaaaaaa  or  aaaaaaaa

      long i = ZNP_MIN(nonzero, U);
      long last_zero_fwd_bfly = ZNP_MAX(nonzero - i, length);
      long last_zero_cross_bfly = ZNP_MIN(nonzero - i, length);

      zn_pmf_t op_ptr = op->data + skip * (--i);

      // First some forward butterflies ("Aa" => "B?") to make them look like:
      // CCCCAABB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAAAA??  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= last_zero_fwd_bfly; i--, op_ptr -= skip)
      {
         // (2*a0, ?) -> (a0, ?)   = (b0, ?)
         zn_pmf_divby2(op_ptr, M, mod);
      }

      // Then some forward butterflies ("AA" => "B?") to make them look like:
      // CCCCBBBB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= (long) length; i--, op_ptr -= skip)
      {
         // (2*a0, 2*a1) -> (a0 + a1, ?)   = (b0, ?)
         zn_pmf_add(op_ptr, op_ptr + half_skip, M, mod);
         zn_pmf_divby2(op_ptr, M, mod);
      }

      // Transform the first row to make them look like:
      // BBBB*???      BBBB*???      BBBB*???      BBBB*???
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      zn_pmf_vec_ifft_small(op, length, forward, ZNP_MIN(nonzero, U),
                            twist << 1);

      // Cross butterflies ("Ba" => "A?") to make them look like:
      // BBBB*???      BBAA*???      AAAA*???      AAAA*???
      // AAAA????  or  AA??????  or  ??????aa  or  ????aaaa
      for (; i >= last_zero_cross_bfly; i--, op_ptr -= skip)
      {
         // (b0, ?) -> (2*b0, ?)    = (2*a0, ?)
         zn_pmf_add(op_ptr, op_ptr, M, mod);
      }
      
      // Cross butterflies ("BA" => "A?") to make them look like:
      // AAAA*???      AAAA*???      AAAA*???      AAAA*???
      // ????????  or  ????????  or  ??????aa  or  ????aaaa
      for (; i >= 0; i--, op_ptr -= skip)
      {
         // (b0, 2*a1) -> (2*b0 - 2*a1, ?)     = (2*a0, ?)
         zn_pmf_add(op_ptr, op_ptr, M, mod);
         zn_pmf_sub(op_ptr, op_ptr + half_skip, M, mod);
      }
   }
   else
   {
      // The input looks like one of these:
      // CCCCCCCC                     CCCCCCCC
      // AAAAaaaa (forward == 1)  or  CCCAAAaa
   
      // Transform first row (no truncation necessary) to make them look like:
      // BBBBBBBB                     BBBBBBBB
      // AAAAaaaa (forward == 1)  or  CCCAAAaa
      zn_pmf_vec_ifft_notrunc_iterative(op, twist << 1);

      long i = U - 1;
      ulong r = M >> op->lgK;     // 2M/K = index for K-th root of unity
      ulong s = twist + r * i;
      zn_pmf_t op_ptr = op->data + skip * i;
      
      long last_zero_cross_bfly = nonzero - U;
      long last_cross_bfly = length - U;
      
      // Cross butterflies ("Ba" => "AB") to make them look like:
      // BBBBAAAA                     BBBBBBAA
      // AAAABBBB (forward == 1)  or  CCCAAABB
      for (; i >= last_zero_cross_bfly; i--, s -= r, op_ptr -= skip)
      {
         // (b0, ?) -> (2*b0, w*b0)     = (2*a0, b1)
         zn_pmf_set(op_ptr + half_skip, op_ptr, M);
         zn_pmf_rotate(op_ptr + half_skip, s);
         zn_pmf_add(op_ptr, op_ptr, M, mod);
      }

      // Cross butterflies ("BA" => "AB") to make them look like:
      // AAAAAAAA                     BBBAAAAA
      // BBBBBBBB (forward == 1)  or  CCCBBBBB
      for (; i >= last_cross_bfly; i--, s -= r, op_ptr -= skip)
      {
         // (b0, 2*a1) -> (2*(b0-a1), w*(b0-2*a1))    = (2*a0, b1)
         zn_pmf_sub(op_ptr + half_skip, op_ptr, M, mod);
         zn_pmf_sub(op_ptr, op_ptr + half_skip, M, mod);
         zn_pmf_rotate(op_ptr + half_skip, M + s);
      }

      // Transform second row to make them look like:
      // AAAAAAAA                     BBBAAAAA
      // *??????? (forward == 1)  or  BBB*????
      op->data += half_skip;
      zn_pmf_vec_ifft_small(op, length - U, forward, U, twist << 1);
      op->data -= half_skip;

      // Inverse butterflies ("BB" => "AA") to make them look like:
      // AAAAAAAA                     AAAAAAAA
      // *??????? (forward == 1)  or  AAA*????
      for (; i >= 0; i--, s -= r, op_ptr -= skip)
      {
         // (b0, b1) -> (b0 + w*b1, b0 - w*b1)    = (2*a0, 2*a1)
         zn_pmf_rotate(op_ptr + half_skip, M - s);
         zn_pmf_bfly(op_ptr + half_skip, op_ptr, M, mod);
      }
   }

   // pop back to full size
   op->K <<= 1;
   op->lgK++;
}



void zn_pmf_vec_ifft_factor(zn_pmf_vec_t op, unsigned lgT, ulong length,
                            int forward, ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(nonzero <= op->K);
   ZNP_ASSERT(length <= nonzero);
   ZNP_ASSERT(length + forward <= op->K);
   ZNP_ASSERT(lgT > 0  &&  lgT < op->lgK);

   if (nonzero == 0)
   {
      if (forward)
         zn_pmf_zero(op->data, op->M);
      return;
   }

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;

   // Write length = U * length_T + length_U, where 0 <= length_U < U
   ulong length_U = length & (U - 1);
   ulong length_T = length >> lgU;

   // Write nonzero = U * nonzero_T + nonzero_U, where 0 <= nonzero_U < U
   ulong nonzero_U = nonzero & (U - 1);
   ulong nonzero_T = nonzero >> lgU;

   ulong r = op->M >> (lgK - 1);     // 2M/K = index for K-th root of unity
   ulong s, i;
   ulong twist_T = twist << lgT;

   // where:
   // symbols in the following diagrams:
   // A = fully untransformed coefficient (one of the a_i)
   // B = intermediate coefficient
   // C = fully transformed coefficient (one of the b_k)
   // ? = garbage that we don't care about
   // * = the "forward" C coefficient, or "?" if no forward coefficient
   //     requested

   // The input looks something like this:
   //
   // CCCCCCCC
   // CCCCCCCC
   // CCCAAAAA
   // AAAAAAAA
   //
   // (we won't bother marking in the locations of zeroes on the diagrams)
   //
   // The horizontal transforms convert between B and C.
   // The vertical transforms convert between A and B.

   // First do row transforms for complete rows, to make it look like:
   // BBBBBBBB
   // BBBBBBBB
   // CCCAAAAA
   // AAAAAAAA
   op->lgK = lgU;
   op->K = U;
   for (i = 0; i < length_T; i++, op->data += skip_U)
      zn_pmf_vec_ifft(op, U, 0, U, twist_T);

   // Column transforms for the rightmost columns, to obtain
   // BBBAAAAA
   // BBBAAAAA
   // CCCBBBBB
   // AAA?????
   op->lgK = lgT;
   op->K = T;
   op->skip = skip_U;

   for (i = length_U, op->data = data + (skip * length_U),
        s = twist + (r * length_U); i < nonzero_U;
        i++, op->data += skip, s += r)
   {
      zn_pmf_vec_ifft(op, length_T, length_U || forward, nonzero_T + 1, s);
   }
   if (nonzero_T)
   {
      for (; i < U; i++, op->data += skip, s += r)
         zn_pmf_vec_ifft(op, length_T, length_U || forward, nonzero_T, s);
   }

   // If there is still a partial row to deal with....
   if (length_U || forward)
   {
      // Transform the partial row to obtain
      // BBBAAAAA
      // BBBAAAAA
      // BBB*????
      // AAA?????
      op->data = data + length_T * skip_U;
      op->lgK = lgU;
      op->K = U;
      op->skip = skip;
      zn_pmf_vec_ifft(op, length_U, forward,
                      nonzero_T ? U : nonzero_U, twist_T);

      // Column transforms for the leftmost columns, to obtain
      // AAAAAAAA
      // AAAAAAAA
      // AAA*????
      // ????????
      op->lgK = lgT;
      op->K = T;
      op->skip = skip_U;
      
      for (i = 0, op->data = data, s = twist; i < length_U && i < nonzero_U;
           i++, op->data += skip, s += r)
      {
         zn_pmf_vec_ifft(op, length_T + 1, 0, nonzero_T + 1, s);
      }
      if (nonzero_T)
      {
         for (; i < length_U; i++, op->data += skip, s += r)
            zn_pmf_vec_ifft(op, length_T + 1, 0, nonzero_T, s);
      }
   }
   
   // restore parameters
   op->lgK = lgK;
   op->K = K;
   op->skip = skip;
   op->data = data;
}



void zn_pmf_vec_ifft(zn_pmf_vec_t op, ulong length, int forward,
                     ulong nonzero, ulong twist)
{
   ZNP_ASSERT(op->lgK <= op->lgM + 1);
   ZNP_ASSERT(twist * op->K < 2*op->M);
   ZNP_ASSERT(nonzero <= op->K);
   ZNP_ASSERT(length <= nonzero);
   ZNP_ASSERT(length + forward <= op->K);

   if (op->K <= 2  ||  2 * op->K * op->M * sizeof(ulong) <= ZNP_CACHE_SIZE)
   {
      // IFFT is pretty small; use small version
      zn_pmf_vec_ifft_small(op, length, forward, nonzero, twist);
   }
   else
   {
      // IFFT is relatively big; use factoring algorithm instead
      zn_pmf_vec_ifft_factor(op, op->lgK / 2, length, forward, nonzero, twist);
   }
}



/* ============================================================================

     main array multiplication routine

============================================================================ */


void mul_fft_params(unsigned* lgK, unsigned* lgM,
                    ulong* coeffs1, ulong* coeffs2,
                    size_t len1, size_t len2)
{
   unsigned _lgM;
   size_t _coeffs1, _coeffs2, _coeffs3;
   ulong M;

   // increase lgM until all the conditions are satisfied
   for (_lgM = 1; ; _lgM++)
   {
      _coeffs1 = CEIL_DIV_2EXP(len1, _lgM - 1);      // = ceil(len1 / (M/2))
      _coeffs2 = CEIL_DIV_2EXP(len2, _lgM - 1);      // = ceil(len2 / (M/2))
      _coeffs3 = _coeffs1 + _coeffs2 - 1;

      M = 1UL << _lgM;
      if (_coeffs3 <= 2*M)
         break;
   }

   *lgM = _lgM;
   *lgK = (_coeffs3 > M) ? (_lgM + 1) : _lgM;
   *coeffs1 = _coeffs1;
   *coeffs2 = _coeffs2;
}



ulong zn_array_mul_fft_get_fudge(size_t len1, size_t len2, int squaring,
                                 const zn_mod_t mod)
{
   unsigned lgK, lgM;
   ulong coeffs1, coeffs2;
   mul_fft_params(&lgK, &lgM, &coeffs1, &coeffs2, len1, len2);
   
   // need to divide by 2^lgK coming from FFT
   ulong fudge1 = zn_mod_pow2(-lgK, mod);
   // and take into account fudge from pointwise multiplies
   ulong fudge2 = zn_pmf_vec_mul_get_fudge(lgM, squaring, mod);
   
   return zn_mod_mul(fudge1, fudge2, mod);
}



void zn_array_mul_fft(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, ulong scale,
                      const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   unsigned lgK, lgM;
   
   // number of zn_pmf_t coefficients for each input poly
   ulong coeffs1, coeffs2;

   // figure out how big the transform needs to be
   mul_fft_params(&lgK, &lgM, &coeffs1, &coeffs2, len1, len2);
   
   // number of zn_pmf_t coefficients for output poly
   ulong length = coeffs1 + coeffs2 - 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ptrdiff_t skip = M + 1;
   
   zn_pmf_vec_t vec1;
   
   int squaring = (op1 == op2  &&  len1 == len2);

   if (!squaring)
   {
      // multiplying two distinct inputs

      // split inputs into zn_pmf_t's and perform FFTs
      zn_pmf_vec_init(vec1, lgK, skip, lgM, mod);
      fft_split(vec1, op1, len1, 0, 1, 0);
      zn_pmf_vec_fft(vec1, length, coeffs1, 0);

      zn_pmf_vec_t vec2;
      
      // note: we apply the fudge factor here, because the second input is
      // shorter than both the first input and the output :-)
      zn_pmf_vec_init(vec2, lgK, skip, lgM, mod);
      fft_split(vec2, op2, len2, 0, scale, 0);
      zn_pmf_vec_fft(vec2, length, coeffs2, 0);

      // pointwise multiplication
      zn_pmf_vec_mul(vec1, vec1, vec2, length, 1);

      zn_pmf_vec_clear(vec2);
   }
   else
   {
      // squaring a single input
   
      // split input into zn_pmf_t's and perform FFTs
      zn_pmf_vec_init(vec1, lgK, skip, lgM, mod);
      fft_split(vec1, op1, len1, 0, 1, 0);
      zn_pmf_vec_fft(vec1, length, coeffs1, 0);

      // pointwise multiplication
      zn_pmf_vec_mul(vec1, vec1, vec1, length, 1);
   }

   // inverse FFT, and write output
   zn_pmf_vec_ifft(vec1, length, 0, length, 0);
   size_t len_res = len1 + len2 - 1;
   fft_combine(res, len_res, vec1, length, 0);

   zn_pmf_vec_clear(vec1);
   
   // if we're squaring, then we haven't applied the fudge factor yet,
   // so do it now
   if (scale != 1 && squaring)
      zn_array_scalar_mul(res, res, len_res, scale, mod);
}


// end of file ****************************************************************
