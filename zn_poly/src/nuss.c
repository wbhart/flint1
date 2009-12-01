/*
   nuss.c:  negacyclic multiplications via Nussbaumer's algorithm
   
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
   The main routine exported from this module is nuss_mul(). This
   takes two arrays of length L = 2^lgL and computes their negacyclic
   convolution, using Nussbaumer's convolution algorithm [Nus80]. These are
   used primarily for the pointwise multiplications in the Schonhage FFT
   (see mul_fft.c).
   
   It is optimised for *small* problems; we don't worry much about locality.

   For example, on my development machine, for 63-bit coefficients it beats
   KS multiplication at length 256; for 8-bit coefficients it wins at length
   2048.
   
   The algorithm is as follows. Let R = Z/mZ. The input polynomials live in
   R[X]/(X^L + 1). We map the convolution problem to S[Z]/(Z^K - 1),
   where S = R[Y]/(Y^M + 1), M and K are powers of two, MK = 2L, K <= 2M,
   and M is minimal subject to these conditions (in other words M and K are
   around sqrt(L), and S has K-th roots of unity). An input polynomial 

      \sum_{j=0}^{M-1} \sum_{i=0}^{K/2 - 1} a_{i + j*K/2} X^{i + j*K/2}

   is mapped to

      \sum_i \sum_j (a_{i + j*K/2} Y^j) Z^i.
   
   Note that this map involves *transposing* the input data, unlike the
   splitting step in the main Schonhage convolution routine. The nice thing
   about this map is that the original negacyclic property of R[X]/(X^L + 1)
   gets reflected neatly in S. The inverse map sends Z -> X, Y -> X^(K/2).

*/

#include "zn_poly_internal.h"


/*
   The functions nuss_split() and nuss_fft() together perform the splitting
   and FFT stages of Nussbaumer's algorithm. 
   
   The input array op has length M*K/2 (where M = res->M, K = res->K).

   In effect the two routines accomplish the following. The input is split
   into M chunks of length K/2, which are then transposed into the first K/2
   coefficients of res. The last K/2 coefficients are set to zero. Then we
   compute the DFT:

      b_k = sum_{i=0}^{K-1} w^{ik} a_i,
   
   where w is the standard K-th root of unity (i.e. Y^(2M/K)). The FFT is
   inplace, and outputs are in bit-reversed order.
   
   The nuss_split() function actually incorporates the first two layers of the
   FFT (to avoid unnecessary operations on zero coefficients). The nuss_fft()
   function handles the remaining lgK - 2 layers.

   We require that 2M >= K >= 4.
*/


#define nuss_split \
    ZNP_nuss_split
void
nuss_split (pmfvec_t res, const ulong* op)
{
   ZNP_ASSERT (res->lgK >= 2);
   ZNP_ASSERT (res->lgM + 1 >= res->lgK);

   // Let b[0], ..., b[K/2 - 1] be the fourier coefficients obtained by
   // performing the split. (The coefficients b[K/2], ..., b[K - 1] are zero.)

   // After the first FFT pass, the coefficients would be
   //     b[0], b[1], ..., b[K/2 - 1],
   // and then
   //     b[0], w * b[1], ..., w^(K/2 - 1) * b[K/2 - 1].
   
   // After the second FFT pass, the coefficients would be
   //                (b[j] +     b[j + K/4])   for 0 <= j < K/4,
   //         w^(2j) (b[j] -     b[j + K/4])   for 0 <= j < K/4,
   //            w^j (b[j] + I * b[j + K/4])   for 0 <= j < K/4,
   //         w^(3j) (b[j] - I * b[j + K/4])   for 0 <= j < K/4,
   // where I = w^(K/4) is the fourth root of unity.
   
   // We do all this in one pass, computing the four combinations
   //   b[j] +     b[j + K/4]
   //   b[j] -     b[j + K/4]
   //   b[j] + I * b[j + K/4]
   //   b[j] - I * b[j + K/4]
   // directly from the input, simultaneously with the transposition, and
   // then throw the w^j adjustments into the bias fields.
   
   ulong M = res->M, K = res->K;
   const zn_mod_struct* mod = res->mod;

   // dest points to j-th output coefficient
   ulong* dest = res->data + 1;
   // dest + half points to (j + K/4)-th output coefficient
   ptrdiff_t half = res->skip << (res->lgK - 2);

   // w = Y^r is the primitive K-th root of unity
   ulong r = M >> (res->lgK - 1);
   ulong i, j, s = 0;
   
   for (j = 0; j < K/4; j++, dest += res->skip, s += r)
   {
      const ulong* src = op + j;
      
      // apply twists
      dest[-1] = 0;
      dest[-1 + half] = 2 * s;
      dest[-1 + 2 * half] = s;
      dest[-1 + 3 * half] = 3 * s;
      
      // do quadruple butterfly and transposition
      if (zn_mod_is_slim (mod))
      {
         // slim version
         for (i = 0; i < M/2; i++, src += K/2)
         {
            ulong x0 = src[0];
            ulong x1 = src[K/4];
            ulong x2 = src[M*K/4];
            ulong x3 = src[M*K/4 + K/4];
            
            dest[i]                   = zn_mod_add_slim (x0, x1, mod);
            dest[i + half]            = zn_mod_sub_slim (x0, x1, mod);
            dest[i + 2 * half]        = zn_mod_sub_slim (x0, x3, mod);
            dest[i + 3 * half]        = zn_mod_add_slim (x0, x3, mod);
            dest[i + M/2]             = zn_mod_add_slim (x2, x3, mod);
            dest[i + half + M/2]      = zn_mod_sub_slim (x2, x3, mod);
            dest[i + 2 * half + M/2]  = zn_mod_add_slim (x2, x1, mod);
            dest[i + 3 * half + M/2]  = zn_mod_sub_slim (x2, x1, mod);
         }
      }
      else
      {
         // non-slim version
         for (i = 0; i < M/2; i++, src += K/2)
         {
            ulong x0 = src[0];
            ulong x1 = src[K/4];
            ulong x2 = src[M*K/4];
            ulong x3 = src[M*K/4 + K/4];
            
            dest[i]                   = zn_mod_add (x0, x1, mod);
            dest[i + half]            = zn_mod_sub (x0, x1, mod);
            dest[i + 2 * half]        = zn_mod_sub (x0, x3, mod);
            dest[i + 3 * half]        = zn_mod_add (x0, x3, mod);
            dest[i + M/2]             = zn_mod_add (x2, x3, mod);
            dest[i + half + M/2]      = zn_mod_sub (x2, x3, mod);
            dest[i + 2 * half + M/2]  = zn_mod_add (x2, x1, mod);
            dest[i + 3 * half + M/2]  = zn_mod_sub (x2, x1, mod);
         }
      }
   }
}



#define nuss_fft \
    ZNP_nuss_fft
void
nuss_fft (pmfvec_t op)
{
   ZNP_ASSERT (op->lgK >= 2);
   ZNP_ASSERT (op->lgM + 1 >= op->lgK);
   
   if (op->lgK == 2)
      return;

   const zn_mod_struct* mod = op->mod;
   ulong M = op->M;
   ulong s, r = op->M >> (op->lgK - 3);
   ptrdiff_t half = op->skip << (op->lgK - 3);
   ulong* end = op->data + (op->skip << op->lgK);
   ulong* start;
   ulong* p;

   for (; r <= M; r <<= 1, half >>= 1)
   for (s = 0, start = op->data; s < M; s += r, start += op->skip)
   for (p = start; p < end; p += 2 * half)
   {
      pmf_bfly (p, p + half, M, mod);
      pmf_rotate (p + half, M + s);
   }
}


/*
   Inverse FFT, i.e. computes

      a_i = sum_{k=0}^{K-1} w^{-ik} b_k.

   The inputs are in bit-reversed order, and outputs in usual order.
*/
#define nuss_ifft \
    ZNP_nuss_ifft
void
nuss_ifft (pmfvec_t op)
{
   const zn_mod_struct* mod = op->mod;
   ulong M = op->M;
   ulong s, r = M;
   ulong r_last = op->M >> (op->lgK - 1);
   ptrdiff_t half = op->skip;
   ulong* end = op->data + (op->skip << op->lgK);
   ulong* p;
   ulong* start;

   for (; r >= r_last; r >>= 1, half <<= 1)
   for (s = 0, start = op->data; s < M; s += r, start += op->skip)
   for (p = start; p < end; p += 2 * half)
   {
      pmf_rotate (p + half, M - s);
      pmf_bfly (p + half, p, M, mod);
   }
}


/*
   This routine performs the reverse Nussbaumer substitution, i.e.
   maps Z -> X, Y -> X^(K/2), performing appropriate additions/subtractions
   to combine the overlapping coefficients, stores results in the res array,
   of length M*K/2.
   
   The main complication in this routine is the bias field in the fourier
   coefficients, i.e. the coefficients might be rotated by random angles
   (and the data needs to be transposed too). We still do everything in a
   single pass.
*/
#define nuss_combine \
    ZNP_nuss_combine
void
nuss_combine (ulong* res, const pmfvec_t op)
{
   ulong i, j;
   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   // src1 points to i-th fourier coefficient
   // src2 points to (i + K/2)-th fourier coefficient
   // These get added/subtracted into res[i + j*K/2], 0 <= j < M.
   ulong* src1 = op->data + 1;
   ulong* src2 = op->data + op->skip * op->K/2 + 1;
   
   for (i = 0; i < op->K/2; i++, src1 += op->skip, src2 += op->skip)
   {
      // We'll be writing to dest[j*K/2], 0 <= j < M
      ulong* dest = res + i;

      // We want to start reading from the Y^i coefficient of src1, which
      // is located at index -bias(src1) mod 2M. But really there are only
      // M coefficients, so we reduce the index mod M, and put neg1 = 1 if
      // we have to negate the coefficients (i.e. they wrapped around
      // negacyclically).
      ulong s1 = (-src1[-1]) & (2*M - 1);
      int neg1 = (s1 >= M);
      if (neg1)
         s1 -= M;

      // Ditto for s2 and neg2 with respect to src2, except that we want
      // to start reading from the coefficient of Y^(i-1).
      ulong s2 = (-1 - src2[-1]) & (2*M - 1);
      int neg2 = (s2 >= M);
      if (neg2)
         s2 -= M;

      // Swap the inputs so that s1 <= s2. (Actually we don't want to disturb
      // src1 and src2, so put the pointers into x1 and x2 instead.)
      ulong* x1;
      ulong* x2;
      if (s1 < s2)
      {
         x1 = src1;
         x2 = src2;
      }
      else
      {
         x1 = src2;
         x2 = src1;
         ulong s_temp = s1; s1 = s2; s2 = s_temp;
         int neg_temp = neg1; neg1 = neg2; neg2 = neg_temp;
      }

      // Okay, now the picture looks like this:
      //
      //     0      s1                             M
      // x1: CCCCCCCCAAAAAAAAAAAAABBBBBBBBBBBBBBBBBB
      //
      //                              s2           M
      // x2: BBBBBBBBBBBBBBBBBBCCCCCCCCAAAAAAAAAAAAA

      // Combine the portions marked AAAA
      dest = zn_skip_array_signed_add (dest, op->K/2, M - s2, x2 + s2, neg2,
                                       x1 + s1, neg1, mod);

      // Combine the portions marked BBBB (x2 stuff gets negated)
      dest = zn_skip_array_signed_add (dest, op->K/2, s2 - s1, x2, !neg2,
                                       x1 + s1 + M - s2, neg1, mod);

      // Combine the portions marked CCCC (both inputs get negated)
      zn_skip_array_signed_add (dest, op->K/2, s1, x2 + s2 - s1, !neg2,
                                x1, !neg1, mod);
   }
}



#define nuss_pointwise_mul_fudge \
    ZNP_nuss_pointwise_mul_fudge
ulong
nuss_pointwise_mul_fudge (unsigned lgM, int sqr, const zn_mod_t mod)
{
   ulong M = 1UL << lgM;
   return _zn_array_mul_fudge (M, M, sqr, mod);
}



/*
   Multiplies fourier coefficients in op1 by those in op2, stores results
   in res. Inplace operation is okay. Automatically uses faster squaring
   version if inputs are the same pmfvec_t object.
   
   NOTE: for now this routine never recurses into Nussbaumer multiplication.
   Doing so would start to become relevant when the original multiplication
   problem has length around 10^9.
   
   The result comes out divided by a fudge factor, which can be recovered via
   nuss_pointwise_mul_fudge().
*/
#define nuss_pointwise_mul \
    ZNP_nuss_pointwise_mul
void
nuss_pointwise_mul (pmfvec_t res, const pmfvec_t op1, const pmfvec_t op2)
{
   ZNP_ASSERT (pmfvec_compatible (res, op1));
   ZNP_ASSERT (pmfvec_compatible (res, op2));

   ulong i, M = res->M;
   pmf_t dest = res->data;
   pmf_const_t src1 = op1->data;
   pmf_const_t src2 = op2->data;

   ZNP_FASTALLOC (temp, ulong, 6624, 2 * M);

   temp[2*M - 1] = 0;

   for (i = 0; i < res->K; i++, dest += res->skip,
        src1 += op1->skip, src2 += op2->skip)
   {
      // add biases
      dest[0] = src1[0] + src2[0];

      // plain multiplication...
      _zn_array_mul (temp, src1 + 1, M, src2 + 1, M, 1, res->mod);
      // ... negacyclic reduction.
      zn_array_sub (dest + 1, temp, temp + M, M, res->mod);
   }

   ZNP_FASTFREE (temp);
}


void
nuss_params (unsigned* lgK, unsigned* lgM, unsigned lgL)
{
   *lgK = (lgL / 2) + 1;
   *lgM = lgL - *lgK + 1;
}


ulong
nuss_mul_fudge (unsigned lgL, int sqr, const zn_mod_t mod)
{
   unsigned lgK, lgM;
   nuss_params (&lgK, &lgM, lgL);
   
   // need to divide by 2^lgK coming from FFT
   ulong fudge1 = zn_mod_pow2 (-lgK, mod);
   // and take into account fudge from pointwise multiplies
   ulong fudge2 = nuss_pointwise_mul_fudge (lgM, sqr, mod);
   
   return zn_mod_mul (fudge1, fudge2, mod);
}


void
nuss_mul (ulong* res, const ulong* op1, const ulong* op2,
                pmfvec_t vec1, pmfvec_t vec2)
{
   ZNP_ASSERT (vec1->lgM + 1 >= vec1->lgK);
   
   if (op1 != op2)
   {
      ZNP_ASSERT (pmfvec_compatible (vec1, vec2));
      
      // split inputs into fourier coefficients and perform FFTs
      nuss_split (vec1, op1);
      nuss_fft (vec1);
      nuss_split (vec2, op2);
      nuss_fft (vec2);
      
      // multiply fourier coefficients into vec1
      nuss_pointwise_mul (vec1, vec1, vec2);
   }
   else
   {
      // split input into fourier coefficients and perform FFT
      nuss_split (vec1, op1);
      nuss_fft (vec1);
      
      // square fourier coefficients into vec1
      nuss_pointwise_mul (vec1, vec1, vec1);
   }

   // inverse FFT
   nuss_ifft (vec1);

   // recombine into result
   nuss_combine (res, vec1);
}


// end of file ****************************************************************
