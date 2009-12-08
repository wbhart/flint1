/*
   pmfvec_fft.c:  FFT/IFFT and transposed FFT/IFFT routines for pmfvec_t
   
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

#include <stdio.h>
#include "zn_poly_internal.h"


/* ============================================================================

     FFT routines

============================================================================ */


void
pmfvec_fft_basecase (pmfvec_t op, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   
   if (op->lgK == 0)
      return;

   // just plain butterfly loop

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M >> (op->lgK - 1);
   ptrdiff_t half = op->skip << (op->lgK - 1);
   ulong* end = op->data + (op->skip << op->lgK);
   ulong* p;
   ulong* start;
   
   for (; r <= M; r <<= 1, half >>= 1, t <<= 1)
   for (start = op->data, s = t; s < M; s += r, start += op->skip)
   for (p = start; p < end; p += 2 * half)
   {
      pmf_bfly (p, p + half, M, mod);
      pmf_rotate (p + half, M + s);
   }
}



void
pmfvec_fft_dc (pmfvec_t op, ulong n, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (n >= 1 && n <= op->K);
   ZNP_ASSERT (z >= 1 && z <= op->K);

   if (op->K == 1)
      return;

   if (n == op->K  &&  z == op->K)
   {
      // No truncation requested; use iterative version
      pmfvec_fft_basecase (op, t);
      return;
   }

   const zn_mod_struct* mod = op->mod;

   // We treat the input as two rows and U = K/2 columns, in row-major order.
   
   // descend to first row (first half of op)
   op->lgK--;
   op->K >>= 1;
   
   long i;
   ulong M = op->M;
   ulong U = op->K;
   ulong* p = op->data;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half = skip << op->lgK;
   ulong z2 = ZNP_MIN (z, U);

   if (n <= U)
   {
      // Only need the first output of the first layer of butterflies.
      for (i = 0; i < (long)(z - U); i++, p += skip)
         pmf_add (p, p + half, M, mod);
      
      // Recurse into top row
      pmfvec_fft_dc (op, n, z2, t << 1);
   }
   else
   {
      // Need both outputs from the first layer of butterflies.
      ulong s = t;
      ulong r = M >> op->lgK;
      
      for (i = 0; i < (long)(z - U); i++, p += skip, s += r)
      {
         pmf_bfly (p, p + half, M, mod);
         pmf_rotate (p + half, M + s);
      }

      // Butterflies where second input is zero
      for (; i < z2; i++, p += skip, s += r)
      {
         pmf_set (p + half, p, M);
         pmf_rotate (p + half, s);
      }
      
      // Recurse into top row...
      pmfvec_fft_dc (op, U, z2, t << 1);

      // ... and recurse into bottom row
      op->data += half;
      pmfvec_fft_dc (op, n - U, z2, t << 1);
      op->data -= half;
   }

   // pop back to whole transform
   op->K <<= 1;
   op->lgK++;
}



/*
   As described above, this splits the length K transform into T = 2^lgT rows
   by U = 2^lgU columns, where K = U * T.
   
   Must have 0 < lgT < lgK.
*/
void
pmfvec_fft_huge (pmfvec_t op, unsigned lgT, ulong n, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (lgT > 0  &&  lgT < op->lgK);
   ZNP_ASSERT (n >= 1 && n <= op->K);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   
   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;
   
   // We need n output coefficients, starting from the top-left, in row-major
   // order.

   // Write n = U * nT + nU, where 0 <= nU < U
   ulong nU = n & (U - 1);
   ulong nT = n >> lgU;
   
   // nT_ceil = number of rows of output, including the last partial row
   ulong nT_ceil = nT + (nU > 0);
   
   // Write z = U * zT + zU, where 0 <= zU < U
   ulong zT = z >> lgU;
   ulong zU = z & (U - 1);
   ulong zU2 = zT ? U : zU;

   ulong r = op->M >> (lgK - 1);
   ulong s, i;
   
   // --------------- FFTs along columns

   op->K = T;
   op->lgK = lgT;
   op->skip = skip_U;

   // First handle the columns with zT + 1 input coefficients.
   for (i = 0, s = t; i < zU; i++, op->data += skip, s += r)
      pmfvec_fft (op, nT_ceil, zT + 1, s);
      
   // Handle the remaining columns, which only have zT input coefficients.
   for (; i < zU2; i++, op->data += skip, s += r)
      pmfvec_fft (op, nT_ceil, zT, s);
   
   // --------------- FFTs along rows

   op->data = data;
   op->K = U;
   op->lgK = lgU;
   op->skip = skip;
   t <<= lgT;

   // Handle the first nT rows.
   for (i = 0; i < nT; i++, op->data += skip_U)
      pmfvec_fft (op, U, zU2, t);
   
   // For the last row, we only need the first nU outputs:
   if (nU)
      pmfvec_fft (op, nU, zU2, t);
   
   // --------------- restore parameters

   op->data = data;
   op->K = K;
   op->lgK = lgK;
}



void
pmfvec_fft (pmfvec_t op, ulong n, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (n >= 1 && n <= op->K);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   
   if (op->K <= 2  ||  2 * op->K * op->M * sizeof (ulong) <= ZNP_CACHE_SIZE)
   {
      // FFT is pretty small; use divide-and-conquer
      pmfvec_fft_dc (op, n, z, t);
   }
   else
   {
      // FFT is relatively big; use factoring algorithm instead
      pmfvec_fft_huge (op, op->lgK / 2, n, z, t);
   }
}



/* ============================================================================

     inverse FFT routines

============================================================================ */

void
pmfvec_ifft_basecase (pmfvec_t op, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);

   if (op->lgK == 0)
      return;

   // just plain butterfly loop

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M;
   ulong r_last = M >> (op->lgK - 1);
   t <<= (op->lgK - 1);
   ptrdiff_t half = op->skip;
   ulong* end = op->data + (op->skip << op->lgK);
   ulong* p;
   ulong* start;

   for (; r >= r_last; r >>= 1, half <<= 1, t >>= 1)
   for (start = op->data, s = t; s < M; s += r, start += op->skip)
   for (p = start; p < end; p += 2 * half)
   {
      pmf_rotate (p + half, M - s);
      pmf_bfly (p + half, p, M, mod);
   }
}



void
pmfvec_ifft_dc (pmfvec_t op, ulong n, int fwd, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   ZNP_ASSERT (n + fwd >= 1 && n + fwd <= op->K);
   ZNP_ASSERT (n <= z);

   if (op->K == 1)
      return;
   
   if (n == op->K)
   {
      // No truncation requested; use iterative version
      pmfvec_ifft_basecase (op, t);
      return;
   }

   const zn_mod_struct* mod = op->mod;

   // We treat the input as two rows and U = K / 2 columns, in row-major order.
   
   // descend to first row (first half of op)
   op->K >>= 1;
   op->lgK--;

   ulong M = op->M;
   ulong U = op->K;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half = skip << op->lgK;

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
   
   if (n + fwd <= U)
   {
      // The input could look like one of the following:
      // CCCCAAAA      CCCCAAAA      CCCCAAaa      CCCCaaaa
      // AAAAAAaa  or  AAaaaaaa  or  aaaaaaaa  or  aaaaaaaa

      long zU2 = ZNP_MIN (z, U);
      long last_zero_fwd_bfly = ZNP_MAX (z - zU2, n);
      long last_zero_cross_bfly = ZNP_MIN (z - zU2, n);

      long i = zU2 - 1;
      pmf_t p = op->data + skip * i;

      // First some forward butterflies ("Aa" => "B?") to make them look like:
      // CCCCAABB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAAAA??  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= last_zero_fwd_bfly; i--, p -= skip)
      {
         // (2*a0, ?) -> (a0, ?)   = (b0, ?)
         pmf_divby2 (p, M, mod);
      }

      // Then some forward butterflies ("AA" => "B?") to make them look like:
      // CCCCBBBB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= (long) n; i--, p -= skip)
      {
         // (2*a0, 2*a1) -> (a0 + a1, ?)   = (b0, ?)
         pmf_add (p, p + half, M, mod);
         pmf_divby2 (p, M, mod);
      }

      // Transform the first row to make them look like:
      // BBBB*???      BBBB*???      BBBB*???      BBBB*???
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      pmfvec_ifft_dc (op, n, fwd, zU2, t << 1);

      // Cross butterflies ("Ba" => "A?") to make them look like:
      // BBBB*???      BBAA*???      AAAA*???      AAAA*???
      // AAAA????  or  AA??????  or  ??????aa  or  ????aaaa
      for (; i >= last_zero_cross_bfly; i--, p -= skip)
      {
         // (b0, ?) -> (2*b0, ?)    = (2*a0, ?)
         pmf_add (p, p, M, mod);
      }
      
      // Cross butterflies ("BA" => "A?") to make them look like:
      // AAAA*???      AAAA*???      AAAA*???      AAAA*???
      // ????????  or  ????????  or  ??????aa  or  ????aaaa
      for (; i >= 0; i--, p -= skip)
      {
         // (b0, 2*a1) -> (2*b0 - 2*a1, ?)     = (2*a0, ?)
         pmf_add (p, p, M, mod);
         pmf_sub (p, p + half, M, mod);
      }
   }
   else
   {
      // The input looks like one of these:
      // CCCCCCCC                 CCCCCCCC
      // AAAAaaaa (fwd == 1)  or  CCCAAAaa
   
      // Transform first row (no truncation necessary) to make them look like:
      // BBBBBBBB                 BBBBBBBB
      // AAAAaaaa (fwd == 1)  or  CCCAAAaa
      pmfvec_ifft_basecase (op, t << 1);

      long i = U - 1;
      ulong r = M >> op->lgK;
      ulong s = t + r * i;
      pmf_t p = op->data + skip * i;
      
      long last_zero_cross_bfly = z - U;
      long last_cross_bfly = n - U;
      
      // Cross butterflies ("Ba" => "AB") to make them look like:
      // BBBBAAAA                 BBBBBBAA
      // AAAABBBB (fwd == 1)  or  CCCAAABB
      for (; i >= last_zero_cross_bfly; i--, s -= r, p -= skip)
      {
         // (b0, ?) -> (2*b0, w*b0)     = (2*a0, b1)
         pmf_set (p + half, p, M);
         pmf_rotate (p + half, s);
         pmf_add (p, p, M, mod);
      }

      // Cross butterflies ("BA" => "AB") to make them look like:
      // AAAAAAAA                 BBBAAAAA
      // BBBBBBBB (fwd == 1)  or  CCCBBBBB
      for (; i >= last_cross_bfly; i--, s -= r, p -= skip)
      {
         // (b0, 2*a1) -> (2*(b0-a1), w*(b0-2*a1))    = (2*a0, b1)
         pmf_sub (p + half, p, M, mod);
         pmf_sub (p, p + half, M, mod);
         pmf_rotate (p + half, M + s);
      }

      // Transform second row to make them look like:
      // AAAAAAAA                 BBBAAAAA
      // *??????? (fwd == 1)  or  BBB*????
      op->data += half;
      pmfvec_ifft_dc (op, n - U, fwd, U, t << 1);
      op->data -= half;

      // Inverse butterflies ("BB" => "AA") to make them look like:
      // AAAAAAAA                 AAAAAAAA
      // *??????? (fwd == 1)  or  AAA*????
      for (; i >= 0; i--, s -= r, p -= skip)
      {
         // (b0, b1) -> (b0 + w*b1, b0 - w*b1)    = (2*a0, 2*a1)
         pmf_rotate (p + half, M - s);
         pmf_bfly (p + half, p, M, mod);
      }
   }

   // pop back to full size
   op->K <<= 1;
   op->lgK++;
}



void
pmfvec_ifft_huge (pmfvec_t op, unsigned lgT, ulong n, int fwd, ulong z,
                  ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   ZNP_ASSERT (n + fwd >= 1 && n + fwd <= op->K);
   ZNP_ASSERT (n <= z);
   ZNP_ASSERT (lgT > 0  &&  lgT < op->lgK);

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;

   // Write n = U * nT + nU, where 0 <= nU < U
   ulong nU = n & (U - 1);
   ulong nT = n >> lgU;

   // Write z = U * zT + zU, where 0 <= zU < U
   ulong zU = z & (U - 1);
   ulong zT = z >> lgU;
   ulong zU2 = zT ? U : zU;
   
   ulong mU1 = ZNP_MIN (zU, nU);
   ulong mU2 = ZNP_MAX (zU, nU);

   int fwd2 = nU || fwd;

   ulong r = op->M >> (lgK - 1);
   ulong s, i;
   ulong tT = t << lgT;

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
   for (i = 0; i < nT; i++, op->data += skip_U)
      pmfvec_ifft (op, U, 0, U, tT);

   // Column transforms for the rightmost columns, to obtain
   // BBBAAAAA
   // BBBAAAAA
   // CCCBBBBB
   // AAA?????
   op->lgK = lgT;
   op->K = T;
   op->skip = skip_U;

   for (i = nU, op->data = data + (skip * nU), s = t + (r * nU); i < mU2;
        i++, op->data += skip, s += r)
   {
      pmfvec_ifft (op, nT, fwd2, zT + 1, s);
   }
   for (; i < zU2; i++, op->data += skip, s += r)
      pmfvec_ifft (op, nT, fwd2, zT, s);

   // If there is still a partial row to deal with....
   if (fwd2)
   {
      // Transform the partial row to obtain
      // BBBAAAAA
      // BBBAAAAA
      // BBB*????
      // AAA?????
      op->data = data + nT * skip_U;
      op->lgK = lgU;
      op->K = U;
      op->skip = skip;
      pmfvec_ifft (op, nU, fwd, zU2, tT);

      // Column transforms for the leftmost columns, to obtain
      // AAAAAAAA
      // AAAAAAAA
      // AAA*????
      // ????????
      op->lgK = lgT;
      op->K = T;
      op->skip = skip_U;
      
      for (i = 0, op->data = data, s = t; i < mU1;
           i++, op->data += skip, s += r)
      {
         pmfvec_ifft (op, nT + 1, 0, zT + 1, s);
      }
      for (; i < nU; i++, op->data += skip, s += r)
         pmfvec_ifft (op, nT + 1, 0, zT, s);
   }
   
   // restore parameters
   op->lgK = lgK;
   op->K = K;
   op->skip = skip;
   op->data = data;
}



void
pmfvec_ifft (pmfvec_t op, ulong n, int fwd, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (z <= op->K);
   ZNP_ASSERT (n <= z);
   ZNP_ASSERT (n + fwd <= op->K);

   if (op->K <= 2  ||  2 * op->K * op->M * sizeof(ulong) <= ZNP_CACHE_SIZE)
   {
      // IFFT is pretty small; use use divide-and-conquer
      pmfvec_ifft_dc (op, n, fwd, z, t);
   }
   else
   {
      // IFFT is relatively big; use factoring algorithm instead
      pmfvec_ifft_huge (op, op->lgK / 2, n, fwd, z, t);
   }
}



/* ============================================================================

     transposed FFT routines

============================================================================ */


void
pmfvec_tpfft_basecase (pmfvec_t op, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);

   if (op->lgK == 0)
      return;

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M;
   ulong r_last = M >> (op->lgK - 1);
   t <<= (op->lgK - 1);
   ptrdiff_t half = op->skip;
   ulong* end = op->data + (op->skip << op->lgK);
   ulong* p;
   ulong* start;

   for (; r >= r_last; r >>= 1, half <<= 1, t >>= 1)
   for (start = op->data, s = t; s < M; s += r, start += op->skip)
   for (p = start; p < end; p += 2 * half)
   {
      pmf_rotate (p + half, M + s);
      pmf_bfly (p + half, p, M, mod);
   }
}



void
pmfvec_tpfft_dc (pmfvec_t op, ulong n, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (n >= 1 && n <= op->K);
   ZNP_ASSERT (z >= 1 && z <= op->K);

   if (op->K == 1)
      return;

   if (n == op->K  &&  z == op->K)
   {
      pmfvec_tpfft_basecase (op, t);
      return;
   }
   
   const zn_mod_struct* mod = op->mod;

   op->lgK--;
   op->K >>= 1;
   
   long i;
   ulong M = op->M;
   ulong U = op->K;
   ulong* p = op->data;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half = skip << op->lgK;
   ulong z2 = ZNP_MIN (z, U);

   if (n <= U)
   {
      pmfvec_tpfft_dc (op, n, z2, t << 1);

      for (i = 0; i < (long)(z - U); i++, p += skip)
         pmf_set (p + half, p, M);
   }
   else
   {
      op->data += half;
      pmfvec_tpfft_dc (op, n - U, z2, t << 1);
      op->data -= half;
      pmfvec_tpfft_dc (op, U, z2, t << 1);

      ulong s = t;
      ulong r = M >> op->lgK;
      
      for (i = 0; i < (long)(z - U); i++, p += skip, s += r)
      {
         pmf_rotate (p + half, M + s);
         pmf_bfly (p + half, p, M, mod);
      }

      for (; i < z2; i++, p += skip, s += r)
      {
         pmf_rotate (p + half, s);
         pmf_add (p, p + half, M, mod);
      }
   }

   op->K <<= 1;
   op->lgK++;
}



void
pmfvec_tpfft_huge (pmfvec_t op, unsigned lgT, ulong n, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (lgT > 0  &&  lgT < op->lgK);
   ZNP_ASSERT (n >= 1 && n <= op->K);
   ZNP_ASSERT (z >= 1 && z <= op->K);

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;
   
   ulong nU = n & (U - 1);
   ulong nT = n >> lgU;
   
   ulong nT_ceil = nT + (nU > 0);
   
   ulong zT = z >> lgU;
   ulong zU = z & (U - 1);
   ulong zU2 = zT ? U : zU;

   ulong r = op->M >> (lgK - 1);
   ulong s, i;
   
   op->K = U;
   op->lgK = lgU;
   t <<= lgT;

   for (i = 0; i < nT; i++, op->data += skip_U)
      pmfvec_tpfft (op, U, zU2, t);
   
   if (nU)
      pmfvec_tpfft (op, nU, zU2, t);

   op->data = data;
   op->K = T;
   op->lgK = lgT;
   op->skip = skip_U;
   t >>= lgT;

   for (i = 0, s = t; i < zU; i++, op->data += skip, s += r)
      pmfvec_tpfft (op, nT_ceil, zT + 1, s);
      
   for (; i < zU2; i++, op->data += skip, s += r)
      pmfvec_tpfft (op, nT_ceil, zT, s);
   
   op->data = data;
   op->skip = skip;
   op->K = K;
   op->lgK = lgK;
}



void
pmfvec_tpfft (pmfvec_t op, ulong n, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (n >= 1 && n <= op->K);
   ZNP_ASSERT (z >= 1 && z <= op->K);

   if (op->K <= 2  ||  2 * op->K * op->M * sizeof (ulong) <= ZNP_CACHE_SIZE)
   {
      pmfvec_tpfft_dc (op, n, z, t);
   }
   else
   {
      pmfvec_tpfft_huge (op, op->lgK / 2, n, z, t);
   }
}



/* ============================================================================

     transposed inverse IFFT routines

============================================================================ */

void
pmfvec_tpifft_basecase (pmfvec_t op, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2*op->M);

   if (op->lgK == 0)
      return;

   ulong M = op->M;
   const zn_mod_struct* mod = op->mod;

   ulong s, r = M >> (op->lgK - 1);
   ptrdiff_t half = op->skip << (op->lgK - 1);
   ulong* end = op->data + (op->skip << op->lgK);
   ulong* p;
   ulong* start;
   
   for (; r <= M; r <<= 1, half >>= 1, t <<= 1)
   for (start = op->data, s = t; s < M; s += r, start += op->skip)
   for (p = start; p < end; p += 2 * half)
   {
      pmf_bfly (p, p + half, M, mod);
      pmf_rotate (p + half, M - s);
   }
}



void
pmfvec_tpifft_dc (pmfvec_t op, ulong n, int fwd, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   ZNP_ASSERT (n + fwd >= 1 && n + fwd <= op->K);
   ZNP_ASSERT (n <= z);

   if (op->K == 1)
      return;

   if (n == op->K)
   {
      pmfvec_tpifft_basecase (op, t);
      return;
   }

   const zn_mod_struct* mod = op->mod;

   op->lgK--;
   op->K >>= 1;

   long i;
   ulong M = op->M;
   ulong U = op->K;
   pmf_t p = op->data;
   ptrdiff_t skip = op->skip;
   ptrdiff_t half = skip << op->lgK;

   if (n + fwd <= U)
   {
      long zU2 = ZNP_MIN (z, U);
      long last_zero_fwd_bfly = ZNP_MAX (z - zU2, n);
      long last_zero_cross_bfly = ZNP_MIN (z - zU2, n);

      for (i = 0; i < last_zero_cross_bfly; i++, p += skip)
      {
         pmf_set (p + half, p, M);
         pmf_rotate (p + half, M);
         pmf_add (p, p, M, mod);
      }
      
      for (; i < n; i++, p += skip)
         pmf_add (p, p, M, mod);
      
      pmfvec_tpifft_dc (op, n, fwd, zU2, t << 1);
      
      for (; i < last_zero_fwd_bfly; i++, p += skip)
      {
         pmf_divby2 (p, M, mod);
         pmf_set (p + half, p, M);
      }
      
      for (; i < zU2; i++, p += skip)
         pmf_divby2 (p, M, mod);
   }
   else
   {
      long last_zero_cross_bfly = z - U;
      long last_cross_bfly = n - U;
      ulong r = M >> op->lgK;
      ulong s = t;

      for (i = 0; i < last_cross_bfly; i++, s += r, p += skip)
      {
         pmf_bfly (p, p + half, M, mod);
         pmf_rotate (p + half, M - s);
      }

      op->data += half;
      pmfvec_tpifft_dc (op, n - U, fwd, U, t << 1);
      op->data -= half;
      
      for (; i < last_zero_cross_bfly; i++, s += r, p += skip)
      {
         pmf_rotate (p + half, M + s);
         pmf_sub (p + half, p, M, mod);
         pmf_sub (p, p + half, M, mod);
      }
      
      for (; i < U; i++, s += r, p += skip)
      {
         pmf_add (p, p, M, mod);
         pmf_rotate (p + half, s);
         pmf_add (p, p + half, M, mod);
      }
      
      pmfvec_tpifft_basecase (op, t << 1);
   }

   op->K <<= 1;
   op->lgK++;
}



void
pmfvec_tpifft_huge (pmfvec_t op, unsigned lgT, ulong n, int fwd, ulong z,
                    ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   ZNP_ASSERT (n + fwd >= 1 && n + fwd <= op->K);
   ZNP_ASSERT (n <= z);
   ZNP_ASSERT (lgT > 0  &&  lgT < op->lgK);

   unsigned lgK = op->lgK;
   unsigned lgU = lgK - lgT;
   
   ulong K = op->K;
   ulong T = 1UL << lgT;
   ulong U = 1UL << lgU;
   
   ptrdiff_t skip = op->skip;
   ptrdiff_t skip_U = skip << lgU;
   
   ulong* data = op->data;

   ulong nU = n & (U - 1);
   ulong nT = n >> lgU;

   ulong zU = z & (U - 1);
   ulong zT = z >> lgU;
   ulong zU2 = zT ? U : zU;

   ulong mU1 = ZNP_MIN (zU, nU);
   ulong mU2 = ZNP_MAX (zU, nU);

   int fwd2 = nU || fwd;

   ulong r = op->M >> (lgK - 1);
   ulong s, i;
   ulong tT = t << lgT;

   if (fwd2)
   {
      op->lgK = lgT;
      op->K = T;
      op->skip = skip_U;
      
      for (i = 0, op->data = data, s = t; i < mU1;
           i++, op->data += skip, s += r)
      {
         pmfvec_tpifft (op, nT + 1, 0, zT + 1, s);
      }
      for (; i < nU; i++, op->data += skip, s += r)
         pmfvec_tpifft (op, nT + 1, 0, zT, s);

      op->data = data + nT * skip_U;
      op->lgK = lgU;
      op->K = U;
      op->skip = skip;
      pmfvec_tpifft (op, nU, fwd, zU2, tT);
   }

   op->lgK = lgT;
   op->K = T;
   op->skip = skip_U;

   for (i = nU, op->data = data + (skip * nU), s = t + (r * nU); i < mU2;
        i++, op->data += skip, s += r)
   {
      pmfvec_tpifft (op, nT, fwd2, zT + 1, s);
   }
   for (; i < zU2; i++, op->data += skip, s += r)
      pmfvec_tpifft (op, nT, fwd2, zT, s);

   op->data = data;
   op->skip = skip;
   op->lgK = lgU;
   op->K = U;
   for (i = 0; i < nT; i++, op->data += skip_U)
      pmfvec_tpifft (op, U, 0, U, tT);

   op->data = data;
   op->lgK = lgK;
   op->K = K;
}



void
pmfvec_tpifft (pmfvec_t op, ulong n, int fwd, ulong z, ulong t)
{
   ZNP_ASSERT (op->lgK <= op->lgM + 1);
   ZNP_ASSERT (t * op->K < 2 * op->M);
   ZNP_ASSERT (z >= 1 && z <= op->K);
   ZNP_ASSERT (n + fwd >= 1 && n + fwd <= op->K);
   ZNP_ASSERT (n <= z);

   if (op->K <= 2  ||  2 * op->K * op->M * sizeof (ulong) <= ZNP_CACHE_SIZE)
   {
      pmfvec_tpifft_dc (op, n, fwd, z, t);
   }
   else
   {
      pmfvec_tpifft_huge (op, op->lgK / 2, n, fwd, z, t);
   }
}



// end of file ****************************************************************
