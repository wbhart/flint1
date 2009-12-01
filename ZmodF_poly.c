/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

ZmodF_poly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "ZmodF_mul.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"
#include "fmpz.h"

/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodF_poly_init(ZmodF_poly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;

   poly->storage = (mp_limb_t*) flint_heap_alloc(bufs * (n+1));
     
   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) flint_heap_alloc_bytes(bufs*sizeof(ZmodF_t));
   
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   unsigned long i;
   for (i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}


void ZmodF_poly_clear(ZmodF_poly_t poly)
{
   flint_heap_free(poly->coeffs);
   flint_heap_free(poly->storage);
}

void ZmodF_poly_stack_init(ZmodF_poly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;

   poly->storage = (mp_limb_t*) flint_stack_alloc(bufs * (n+1));
     
   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) flint_stack_alloc_bytes(bufs*sizeof(ZmodF_t));
   
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   unsigned long i;
   for (i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}

void ZmodF_poly_stack_clear(ZmodF_poly_t poly)
{
   flint_stack_release();
   flint_stack_release();
}

/****************************************************************************

   Basic Arithmetic Routines
   
****************************************************************************/


void ZmodF_poly_set(ZmodF_poly_t x, ZmodF_poly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->n == y->n);

   unsigned long i;
   for (i = 0; i < y->length; i++)
      ZmodF_set(x->coeffs[i], y->coeffs[i], x->n);

   x->length = y->length;
}


void ZmodF_poly_pointwise_mul(ZmodF_poly_t res, ZmodF_poly_t x, ZmodF_poly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   unsigned long j;

   ZmodF_mul_info_t info;
   ZmodF_mul_info_init(info, x->n, x == y);
   
   ulong i;
   if (x != y)
      for (i = 0; i < x->length; i++)
      {
         if (i+8 < x->length)
         {
            for (j = 0; j < x->n; j += 8) FLINT_PREFETCH(x->coeffs[i+8], j);
            for (j = 0; j < y->n; j += 8) FLINT_PREFETCH(y->coeffs[i+8], j);
         }
         ZmodF_mul_info_mul(info, res->coeffs[i], x->coeffs[i], y->coeffs[i]);
      }
   else
   {
      unsigned long i;
      for (i = 0; i < x->length; i++)
      {
         if (i+8 < x->length)
         {
            for (j = 0; j < x->n; j += 8) FLINT_PREFETCH(x->coeffs[i+8], j);
         }
         ZmodF_mul_info_mul(info, res->coeffs[i], x->coeffs[i], x->coeffs[i]);
      }
   }

   ZmodF_mul_info_clear(info);

   res->length = x->length;
}


void ZmodF_poly_add(ZmodF_poly_t res, ZmodF_poly_t x, ZmodF_poly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   unsigned long i;
   for (i = 0; i < x->length; i++)
      ZmodF_add(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodF_poly_sub(ZmodF_poly_t res, ZmodF_poly_t x, ZmodF_poly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   unsigned long i;
   for (i = 0; i < x->length; i++)
      ZmodF_sub(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodF_poly_normalise(ZmodF_poly_t poly)
{
   unsigned long i;
   for (i = 0; i < poly->length; i++)
      ZmodF_normalise(poly->coeffs[i], poly->n);
}


void ZmodF_poly_rescale(ZmodF_poly_t poly)
{
   if (poly->depth == 0)
      return;

   unsigned long i;
   for (i = 0; i < poly->length; i++)
      ZmodF_short_div_2exp(poly->coeffs[i], poly->coeffs[i],
                           poly->depth, poly->n);
}

void ZmodF_poly_rescale_range(ZmodF_poly_t poly, unsigned long start, unsigned long n)
{
   if (poly->depth == 0)
      return;
      
   unsigned long length = FLINT_MIN(n, poly->length);
   
   unsigned long i;
   for (i = start; i < length; i++)
      ZmodF_short_div_2exp(poly->coeffs[i], poly->coeffs[i],
                           poly->depth, poly->n);
}


/****************************************************************************

   Forward fourier transforms (internal code)

****************************************************************************/


void _ZmodF_poly_FFT_iterative(
            ZmodF_t* x, unsigned long depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   FLINT_ASSERT(depth >= 1);

   unsigned long i, s, start;
   ZmodF_t* y, * z;

   // root is the (2^depth)-th root of unity for the current layer,
   // measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS) >> depth;
   FLINT_ASSERT(twist < root);

   // half = half the current block length
   unsigned long half = 1UL << (depth - 1);
   unsigned long half_skip = half * skip;

   unsigned long layer;

   // =========================================================================
   // Special case for first layer, if root and/or twist involve sqrt2
   // rotations.

   if ((root | twist) & 1)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length <= half)
      {
         // Only need first output for each butterfly,
         // i.e. (a, b) -> (a + b, ?)
         if (nonzero > half)
            for (i = 0, y = x; i < nonzero - half; i++, y += skip)
               ZmodF_add(y[0], y[0], y[half_skip], n);
      }
      else
      {
         // Need both outputs for each butterfly.
         if (nonzero <= half)
         {
            // The second half of each butterfly input are zeroes, so we just
            // computing (a, 0) -> (a, ra), where r is the appropriate root
            // of unity.
            for (i = 0, s = twist, y = x; i < nonzero;
                 i++, s += root, y += skip)
            {
               ZmodF_mul_sqrt2exp(y[half_skip], y[0], s, n);
            }
         }
         else
         {
            // If nonzero > half, then we need some full butterflies...
            for (i = 0, s = twist, y = x; i < nonzero - half;
                 i++, s += root, y += skip)
            {
               ZmodF_forward_butterfly_sqrt2exp(y, y + half_skip, 
                                                scratch, s, n);
            }
            // and also some partial butterflies (a, 0) -> (a, ra).
            for (; i < half; i++, s += root, y += skip)
               ZmodF_mul_sqrt2exp(y[half_skip], y[0], s, n);
         }
      }

      // Here we switch to measuring roots as powers of 2, but we also need
      // to update to the next layer's roots, and these two actions cancel
      // each other out :-)
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
      if (nonzero > 2*half)
         nonzero = 2*half;
   
      layer = 1;
   }
   else
   {
      // no special case for first layer
      layer = 0;

      // switch to measuring roots as powers of 2
      root >>= 1;
      twist >>= 1;
   }

   // =========================================================================
   // This section handles the layers where there are still zero coefficients
   // to take advantage of. In most cases, this will only happen for the
   // first layer or two, so we don't bother with specialised limbshift-only
   // code for these layers.

   // Note: from here on there are no sqrt2 rotations, and we measure all
   // roots as powers of 2.
   
   for (; (layer < depth) && (nonzero < 2*half); layer++)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length_remainder > half)
      {
         // If length overhangs by more than half the block, then we need to
         // perform full butterflies on the last block (i.e. the last block
         // doesn't get any special treatment).
         length_quantised += 2*half;
      }
      else if (length_remainder)
      {
         // If length overhangs the block by at most half the block size,
         // then we only need to compute the first output of each butterfly
         // for this block, i.e. (a, b) -> (a + b)
         if (nonzero > half)
         {
            y = x + skip * length_quantised;
            for (i = 0; i < nonzero - half; i++, y += skip)
               ZmodF_add(y[0], y[0], y[half_skip], n);
         }
      }

      if (nonzero <= half)
      {
         // If nonzero <= half, then the second half of each butterfly input
         // are zeroes, so we just computing (a, 0) -> (a, ra), where r is the
         // appropriate root of unity.
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, s = twist, z = y; i < nonzero;
                 i++, s += root, z += skip)
            {
               ZmodF_mul_2exp(z[half_skip], z[0], s, n);
            }
         }
      }
      else
      {
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            // If nonzero > half, then we need some full butterflies...
            for (i = 0, s = twist, z = y; i < nonzero - half;
                 i++, s += root, z += skip)
            {
               ZmodF_forward_butterfly_2exp(z, z + half_skip, scratch, s, n);
            }
            // and also some partial butterflies (a, 0) -> (a, ra).
            for (; i < half; i++, s += root, z += skip)
               ZmodF_mul_2exp(z[half_skip], z[0], s, n);
         }
      }

      // Update roots of unity
      twist <<= 1;
      root <<= 1;
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
      
      if (nonzero > 2*half)
         // no more zero coefficients to take advantage of:
         nonzero = 2*half;
   }

   // =========================================================================
   // Now we may assume there are no more zero coefficients.

   for (; layer < depth; layer++)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length_remainder > half)
      {
         // If length overhangs by more than half the block, then we need to
         // perform full butterflies on the last block (i.e. the last block
         // doesn't get any special treatment).
         length_quantised += 2*half;
      }
      else if (length_remainder)
      {
         // If length overhangs the block by at most half the block size,
         // then we only need to compute the first output of each butterfly
         // for this block, i.e. (a, b) -> (a + b)
         y = x + skip * length_quantised;
         for (i = 0; i < half; i++, y += skip)
            ZmodF_add(y[0], y[0], y[half_skip], n);
      }
      
      // To keep the inner loops long, we have two versions of the next loop.
      if (layer < depth/2)
      {
         // Version 1: only a few relatively long blocks.
         
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, s = twist, z = y; i < half; i++, s += root, z += skip)
               ZmodF_forward_butterfly_2exp(z, z + half_skip, scratch, s, n);
         }
      }
      else
      {
         // Version 2: lots of short blocks.
         
         // Two sub-versions, depending on whether the rotations are all by
         // a whole number of limbs.
         if ((root | twist) & (FLINT_BITS - 1))
         {
            // Version 2a: rotations still involve bitshifts.
            for (i = 0, s = twist, y = x; i < half; i++, s += root, y += skip)
               for (start = 0, z = y; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_forward_butterfly_2exp(z, z + half_skip,
                                               scratch, s, n);
               }
         }
         else
         {
            // Version 2b: rotations involve only limbshifts.
            unsigned long root_limbs = root >> FLINT_LG_BITS_PER_LIMB;

            if (twist == 0)
            {
               // special case, since ZmodF_forward_butterfly_Bexp doesn't
               // allow zero rotation count
               for (start = 0, z = x; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_simple_butterfly(z, z + half_skip, scratch, n);
               }
               i = 1;
               y = x + skip;
               s = root_limbs;
            }
            else
            {
               i = 0;
               y = x;
               s = twist >> FLINT_LG_BITS_PER_LIMB;
            }
            
            for (; i < half; i++, s += root_limbs, y += skip)
               for (start = 0, z = y; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_forward_butterfly_Bexp(z, z + half_skip,
                                               scratch, s, n);
               }
         }
      }

      // Update roots of unity
      twist <<= 1;
      root <<= 1;
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
   }
}


/*
Factors FFT of length 2^depth into length 2^rows_depth and length 2^cols_depth
transforms
*/
void _ZmodF_poly_FFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(rows_depth >= 1);
   FLINT_ASSERT(cols_depth >= 1);
   
   unsigned long depth = rows_depth + cols_depth;
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS) >> depth;
   FLINT_ASSERT(twist < root);

   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long length_whole_rows = length_cols ?
                                     (length_rows + 1) : length_rows;
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // column transforms
   for (i = 0, y = x, j = twist; i < nonzero_cols; i++, y += skip, j += root)
      _ZmodF_poly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                     length_whole_rows, j, n, scratch);

   if (nonzero_rows)
   {
      for (; i < cols; i++, y += skip, j += root)
         _ZmodF_poly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                        length_whole_rows, j, n, scratch);
      nonzero_cols = cols;
   }
   
   // row transforms
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodF_poly_FFT(y, cols_depth, skip, nonzero_cols, cols,
                     twist << rows_depth, n, scratch);

   if (length_cols)
      // The relevant portion of the last row:
      _ZmodF_poly_FFT(y, cols_depth, skip, nonzero_cols, length_cols,
                     twist << rows_depth, n, scratch);
}



/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of buffers assumed to be nonzero
length = number of fourier coefficients requested
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodF_poly_FFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                    unsigned long nonzero, unsigned long length,
                    unsigned long twist, unsigned long n,
                    ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   FLINT_ASSERT(depth >= 1);

   // If the data fits in L1 (2^depth coefficients of length n+1, plus a
   // scratch buffer), then use the iterative transform. Otherwise factor the
   // FFT into two chunks.
   if (depth == 1 ||
       ((1 << depth) + 1) * (n+1) <= ZMODFPOLY_FFT_FACTOR_THRESHOLD)
   {
      _ZmodF_poly_FFT_iterative(x, depth, skip, nonzero, length,
                               twist, n, scratch);
   }
   else
   {
      unsigned long rows_depth = depth >> 1;
      unsigned long cols_depth = depth - rows_depth;
      _ZmodF_poly_FFT_factor(x, rows_depth, cols_depth, skip, nonzero, length,
                            twist, n, scratch);
   }
}



/****************************************************************************

   Inverse fourier transforms (internal code)

****************************************************************************/


/*
This one is for when there is no truncation.
*/
void _ZmodF_poly_IFFT_iterative(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(depth >= 1);

   // root is the (2^(layer+1))-th root unity for each layer,
   // measured as a power of sqrt2
   long root = 2*n*FLINT_BITS;
   twist <<= (depth - 1);
   FLINT_ASSERT(twist < root);

   unsigned long half = 1;
   unsigned long half_skip = skip;
   unsigned long size = 1UL << depth;
   unsigned long layer, start, i, s;
   ZmodF_t* y, * z;
   
   // First group of layers; lots of small blocks.

   for (layer = 0; layer < depth/2; layer++)
   {
      // no sqrt2 should be involved here
      FLINT_ASSERT(!((twist | root) & 1));

      // change roots to be measured as powers of 2
      // (also this updates for the next layer in advance)
      root >>= 1;
      twist >>= 1;
      
      if ((root | twist) & (FLINT_BITS-1))
      {
         // This version allows bitshifts
         for (i = 0, y = x, s = twist; i < half; i++, s += root, y += skip)
            for (start = 0, z = y; start < size;
                 start += 2*half, z += 2*half_skip)
            {
               ZmodF_inverse_butterfly_2exp(z, z + half_skip, scratch, s, n);
            }
      }
      else
      {
         // This version is limbshifts only
         unsigned long root_limbs = root >> FLINT_LG_BITS_PER_LIMB;

         if (twist == 0)
         {
            // special case since ZmodF_inverse_butterfly_Bexp doesn't allow
            // zero rotation count
            for (start = 0, z = x; start < size;
                 start += 2*half, z += 2*half_skip)
            {
               ZmodF_simple_butterfly(z, z + half_skip, scratch, n);
            }
         
            i = 1;
            s = root_limbs;
            y = x + skip;
         }
         else
         {
            i = 0;
            s = twist >> FLINT_LG_BITS_PER_LIMB;
            y = x;
         }
         
         for (; i < half; i++, s += root_limbs, y += skip)
            for (start = 0, z = y; start < size;
                 start += 2*half, z += 2*half_skip)
            {
               ZmodF_inverse_butterfly_Bexp(z, z + half_skip, scratch, s, n);
            }
      }
      
      half <<= 1;
      half_skip <<= 1;
   }


   // Second group of layers; just a few large blocks.
   
   for (; layer < depth; layer++)
   {
      if ((root | twist) & 1)
      {
         // sqrt2 is involved. This had better be the last layer.
         FLINT_ASSERT(layer == depth - 1);
         
         for (i = 0, z = x, s = twist; i < half; i++, s += root, z += skip)
            ZmodF_inverse_butterfly_sqrt2exp(z, z + half_skip, scratch, s, n);
         
         return;
      }
      else
      {
         // Only bitshifts.

         // change roots to be measured as powers of 2
         // (also this updates for the next layer in advance)
         twist >>= 1;
         root >>= 1;
         
         for (start = 0, y = x; start < size;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, z = y, s = twist; i < half; i++, s += root, z += skip)
               ZmodF_inverse_butterfly_2exp(z, z + half_skip, scratch, s, n);
         }
      }
   
      half <<= 1;
      half_skip <<= 1;
   }
}



/*
This one's for working in L1 when truncation is involved. It splits into
two halves.
*/
void _ZmodF_poly_IFFT_recursive(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long nonzero, unsigned long length, int extra,
               unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   FLINT_ASSERT(depth >= 1);

   long size = 1UL << depth;

   if (length == size)
   {
      // no truncation necessary
      _ZmodF_poly_IFFT_iterative(x, depth, skip, twist, n, scratch);
      return;
   }

   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   long root = (4*n*FLINT_BITS) >> depth;
   FLINT_ASSERT(twist < root);

   long cols = size >> 1;
   long half = skip << (depth - 1);

   // symbols in the following diagrams:
   // A = fully untransformed coefficient
   // a = fully untransformed coefficient (implied zero)
   // B = intermediate coefficient
   // b = intermediate coefficient (implied zero)
   // C = fully transformed coefficient
   // c = fully transformed coefficient (implied zero)
   // ? = garbage that we don't care about
   // * = the extra C coefficient, or "?" if no extra coefficient requested
   
   // the horizontal transforms convert between B and C
   // the vertical butterflies convert between A and B

   if ((length < cols) || (length == cols && !extra))
   {
      // The input could look like one of the following:
      // CCCCAAAA      CCCCAAAA      CCCCAAaa      CCCCaaaa
      // AAAAAAaa  or  AAaaaaaa  or  aaaaaaaa  or  aaaaaaaa

      long i, last_zero_forward_butterfly, last_zero_cross_butterfly;

      if (nonzero <= cols)
      {
         i = nonzero - 1;
         last_zero_forward_butterfly = length;
         last_zero_cross_butterfly = 0;
      }
      else
      {
         i = cols - 1;
         if (nonzero > length + cols)
         {
            last_zero_forward_butterfly = nonzero - cols;
            last_zero_cross_butterfly = length;
         }
         else
         {
            last_zero_forward_butterfly = length;
            last_zero_cross_butterfly = nonzero - cols;
         }
      }
      
      ZmodF_t* y = x + skip*i;

      // First some forward butterflies ("Aa" => "B?") to make them look like:
      // CCCCAABB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAAAA??  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= last_zero_forward_butterfly; i--, y -= skip)
      {
         // (2*a0, ?) -> (a0, ?)   = (b0, ?)
         ZmodF_short_div_2exp(y[0], y[0], 1, n);
      }

      // Then some forward butterflies ("AA" => "B?") to make them look like:
      // CCCCBBBB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= (long)length; i--, y -= skip)
      {
         // (2*a0, 2*a1) -> (a0 + a1, ?)   = (b0, ?)
         ZmodF_add(y[0], y[0], y[half], n);
         ZmodF_short_div_2exp(y[0], y[0], 1, n);
      }

      // Transform the first row to make them look like:
      // BBBB*???      BBBB*???      BBBB*???      BBBB*???
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      if (depth > 1)
         _ZmodF_poly_IFFT_recursive(x, depth - 1, skip,
                                   (nonzero < cols) ? nonzero : cols,
                                   length, extra, twist << 1, n, scratch);
      
      // Cross butterflies ("Ba" => "A?") to make them look like:
      // BBBB*???      BBAA*???      AAAA*???      AAAA*???
      // AAAA????  or  AA??????  or  ??????aa  or  ????aaaa
      for (; i >= last_zero_cross_butterfly; i--, y -= skip)
      {
         // (b0, ?) -> (2*b0, ?)    = (2*a0, ?)
         ZmodF_add(y[0], y[0], y[0], n);
      }
         
      // Cross butterflies ("BA" => "A?") to make them look like:
      // AAAA*???      AAAA*???      AAAA*???      AAAA*???
      // ????????  or  ????????  or  ??????aa  or  ????aaaa
      for (; i >= 0; i--, y -= skip)
      {
         // (b0, 2*a1) -> (2*b0 - 2*a1, ?)     = (2*a0, ?)
         ZmodF_add(y[0], y[0], y[0], n);
         ZmodF_sub(y[0], y[0], y[half], n);
      }
   }
   else
   {
      // The input looks like one of these:
      // CCCCCCCC                   CCCCCCCC
      // AAAAaaaa (extra == 1)  or  CCCAAAaa
   
      // Transform first row (no truncation necessary) to make them look like:
      // BBBBBBBB                   BBBBBBBB
      // AAAAaaaa (extra == 1)  or  CCCAAAaa
      if (depth > 1)
         _ZmodF_poly_IFFT_iterative(x, depth - 1, skip, twist << 1, n, scratch);

      long i = cols - 1;
      unsigned long s = twist + root*i;
      ZmodF_t* y = x + skip*i;
      
      long last_zero_cross_butterfly = nonzero - cols;
      long last_cross_butterfly = length - cols;
   
      // Cross butterflies ("Ba" => "AB") to make them look like:
      // BBBBAAAA                   BBBBBBAA
      // AAAABBBB (extra == 1)  or  CCCAAABB
      for (; i >= last_zero_cross_butterfly; i--, s -= root, y -= skip)
      {
         // (b0, ?) -> (2*b0, w*b0)     = (2*a0, b1)
         ZmodF_mul_sqrt2exp(y[half], y[0], s, n);
         ZmodF_add(y[0], y[0], y[0], n);
      }
         
      // Cross butterflies ("BA" => "AB") to make them look like:
      // AAAAAAAA                   BBBAAAAA
      // BBBBBBBB (extra == 1)  or  CCCBBBBB
      for (; i >= last_cross_butterfly; i--, s -= root, y -= skip)
      {
         // (b0, 2*a1) -> (2*(b0-a1), w*(b0-2*a1))    = (2*a0, b1)
         ZmodF_sub(scratch[0], y[0], y[half], n);
         ZmodF_add(y[0], y[0], scratch[0], n);
         ZmodF_mul_sqrt2exp(y[half], scratch[0], s, n);
      }
      
      // Transform second row to make them look like:
      // AAAAAAAA                   BBBAAAAA
      // *??????? (extra == 1)  or  BBB*????
      if (depth > 1)
         _ZmodF_poly_IFFT_recursive(x + skip*cols, depth - 1, skip, cols,
                                   length - cols, extra, twist << 1, n,
                                   scratch);

      // Inverse butterflies ("BB" => "AA") to make them look like:
      // AAAAAAAA                   AAAAAAAA
      // *??????? (extra == 1)  or  AAA*????
      for (; i >= 0; i--, s -= root, y -= skip)
      {
         // (b0, b1) -> (b0 + w*b1, b0 - w*b1)    = (2*a0, 2*a1)
         ZmodF_inverse_butterfly_sqrt2exp(y, y + half, scratch, s, n);
      }
   }
}



/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodF_poly_IFFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            int extra, unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(rows_depth >= 1);
   FLINT_ASSERT(cols_depth >= 1);

   unsigned long depth = rows_depth + cols_depth;
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS) >> depth;
   FLINT_ASSERT(twist < root);
   
   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // row transforms for the rows where we have all fourier coefficients
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodF_poly_IFFT(y, cols_depth, skip, cols, cols, 0,
                      twist << rows_depth, n, scratch);

   // column transforms where we have enough information
   for (i = length_cols, y = x + (skip * length_cols),
        j = twist + (root*length_cols);
        i < nonzero_cols; i++, y += skip, j += root)
   {
      _ZmodF_poly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                      length_rows, length_cols ? 1 : extra, j, n, scratch);
   }
   if (nonzero_rows)
      for (; i < cols; i++, y += skip, j += root)
         _ZmodF_poly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                         length_rows, length_cols ? 1 : extra, j, n, scratch);

   if (length_cols)
   {
      // a single switcheroo row transform
      _ZmodF_poly_IFFT(x + length_rows * (skip << cols_depth), cols_depth,
                      skip, (nonzero_rows ? cols : nonzero_cols),
                      length_cols, extra, twist << rows_depth, n, scratch);

      // remaining column transforms
      for (i = 0, y = x, j = twist; i < length_cols && i < nonzero_cols;
           i++, y += skip, j += root)
      {
         _ZmodF_poly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                         length_rows + 1, 0, j, n, scratch);
      }
      if (nonzero_rows)
      {
         for (; i < length_cols; i++, y += skip, j += root)
            _ZmodF_poly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                            length_rows + 1, 0, j, n, scratch);
      }
   }
   else if (extra)
   {
      // need one extra trivial fourier coefficient
      x += length_rows * (skip << cols_depth);
      for (i = 1, y = x + skip; i < (nonzero_rows ? cols : nonzero_cols);
           i++, y += skip)
      {
         ZmodF_add(x[0], x[0], y[0], n);
      }
      ZmodF_short_div_2exp(x[0], x[0], cols_depth, n);
   }
}


/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodF_poly_IFFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                     unsigned long nonzero, unsigned long length, int extra,
                     unsigned long twist, unsigned long n,
                     ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   FLINT_ASSERT(depth >= 1);

   // If the data fits in L1 (2^depth coefficients of length n+1, plus a
   // scratch buffer), then use the iterative transform. Otherwise factor the
   // FFT into two chunks.
   if (depth == 1 ||
       ((1 << depth) + 1) * (n+1) <= ZMODFPOLY_FFT_FACTOR_THRESHOLD)
   {
      _ZmodF_poly_IFFT_recursive(x, depth, skip, nonzero, length, extra,
                                twist, n, scratch);
   }
   else
   {
      unsigned long rows_depth = depth >> 1;
      unsigned long cols_depth = depth - rows_depth;
      _ZmodF_poly_IFFT_factor(x, rows_depth, cols_depth, skip, nonzero, length,
                             extra, twist, n, scratch);
   }
}


/****************************************************************************

   Forward "dual" fourier transforms (internal code)

(twists are applied *before* the transform instead of afterwards, so these
are used for e.g. negacyclic transforms)

****************************************************************************/

/*
Let M = 2^depth
2^root = Mth root of unity
input is assumed to be mod x^M - a^M, where a = 2^twist

assumes twist nonzero
*/
void _ZmodF_poly_FFT_dual_recursive(
            ZmodF_t* x, unsigned long depth,
            unsigned long twist, unsigned long root,
            unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(twist);
   FLINT_ASSERT(twist < root);

   // =========================================================================
   // special cases for length <= 4

   if (depth == 2)
   {
      // length == 4
      
      // ----------------------------------------------------------------------
      // Do the outer layer of two butterflies first. This is basically an
      // unrolled version of the length >= 8 case below.

      unsigned long bits = (2*twist) & (FLINT_BITS-1);
      unsigned long limbs = n - (twist >> (FLINT_LG_BITS_PER_LIMB-1));
      
      if (bits)
      {
         // each butterfly needs a bitshift
         bits = FLINT_BITS - bits;
         if (--limbs)
         {
            ZmodF_short_div_2exp(*scratch, x[2], bits, n);
            ZmodF_div_Bexp_add(x[2], x[0], *scratch, limbs, n);
            ZmodF_div_Bexp_sub(x[0], x[0], *scratch, limbs, n);

            ZmodF_short_div_2exp(*scratch, x[3], bits, n);
            ZmodF_div_Bexp_add(x[3], x[1], *scratch, limbs, n);
            ZmodF_div_Bexp_sub(x[1], x[1], *scratch, limbs, n);
         }
         else
         {
            ZmodF_short_div_2exp(*scratch, x[2], bits, n);
            ZmodF_add(x[2], x[0], *scratch, n);
            ZmodF_sub(x[0], x[0], *scratch, n);

            ZmodF_short_div_2exp(*scratch, x[3], bits, n);
            ZmodF_add(x[3], x[1], *scratch, n);
            ZmodF_sub(x[1], x[1], *scratch, n);
         }
      }
      else
      {
         // no bitshifts needed
         ZmodF_div_Bexp_add(*scratch, x[0], x[2], limbs, n);
         ZmodF_swap(scratch, x+2);
         ZmodF_div_Bexp_sub(x[0], x[0], *scratch, limbs, n);

         ZmodF_div_Bexp_add(*scratch, x[1], x[3], limbs, n);
         ZmodF_swap(scratch, x+3);
         ZmodF_div_Bexp_sub(x[1], x[1], *scratch, limbs, n);
      }

      // ----------------------------------------------------------------------
      // Now do the bottom layer, two "blocks" of one butterfly each.

      twist = n*FLINT_BITS - twist;
      ZmodF_inverse_butterfly_2exp(x, x+1, scratch, twist, n);
      ZmodF_swap(x, x+1);
      ZmodF_inverse_butterfly_2exp(x+2, x+3, scratch, twist - root, n);
      ZmodF_swap(x+2, x+3);

      return;
   }

   if (depth <= 1)
   {
      // length == 1 or 2
      if (depth == 1)
      {
         ZmodF_inverse_butterfly_2exp(x, x+1, scratch,
                                      n*FLINT_BITS - twist, n);
         ZmodF_swap(x, x+1);
      }
      return;
   }
   
   // =========================================================================
   // general case for length >= 8

   // butterflies (a, b) -> (a + w*b, a - w*b), where w = 2^(amount).
   unsigned long half = 1 << (depth - 1);
   ZmodF_t* y = x + half;
   unsigned long amount = twist << (depth - 1);
   unsigned long bits = amount & (FLINT_BITS-1);
   unsigned long limbs = n - (amount >> FLINT_LG_BITS_PER_LIMB);
   
   if (bits)
   {
      // each butterfly needs a bitshift
      bits = FLINT_BITS - bits;
      if (--limbs)
      {
         unsigned long i;
         for (i = 0; i < half; i++)
         {
            ZmodF_short_div_2exp(*scratch, y[i], bits, n);
            ZmodF_div_Bexp_add(y[i], x[i], *scratch, limbs, n);
            ZmodF_div_Bexp_sub(x[i], x[i], *scratch, limbs, n);
         }
      }
      else
      {
         unsigned long i;
         for (i = 0; i < half; i++)
         {
            ZmodF_short_div_2exp(*scratch, y[i], bits, n);
            ZmodF_add(y[i], x[i], *scratch, n);
            ZmodF_sub(x[i], x[i], *scratch, n);
         }
      }
   }
   else
   {
      // all butterflies are limbshifts only
      unsigned long i;
      for (i = 0; i < half; i++)
      {
         ZmodF_div_Bexp_add(*scratch, x[i], y[i], limbs, n);
         ZmodF_swap(scratch, y+i);
         ZmodF_div_Bexp_sub(x[i], x[i], *scratch, limbs, n);
      }
   }
   
   // =========================================================================
   // recurse into two halves

   _ZmodF_poly_FFT_dual_recursive(x, depth-1, twist, root << 1, n, scratch);
   _ZmodF_poly_FFT_dual_recursive(x + half, depth-1, twist + root, root << 1,
                                 n, scratch);
}



void _ZmodF_poly_IFFT_dual_recursive(
            ZmodF_t* x, unsigned long depth,
            unsigned long twist, unsigned long root,
            unsigned long n, ZmodF_t* scratch)
{
   // =========================================================================
   // special cases for length <= 4
   
   if (depth == 2)
   {
      // ----------------------------------------------------------------------
      // Do the inner layer of two "blocks" of one butterfly each.

      unsigned long temp = n*FLINT_BITS - twist;
      ZmodF_forward_butterfly_2exp(x+3, x+2, scratch, temp - root, n);
      ZmodF_swap(x+2, x+3);
      ZmodF_forward_butterfly_2exp(x+1, x, scratch, temp, n);
      ZmodF_swap(x, x+1);

      // ----------------------------------------------------------------------
      // Now do the outer layer of two butterflies. This is basically an
      // unrolled version of the length >= 8 case below.

      unsigned long amount = 2*twist;
      unsigned long bits = amount & (FLINT_BITS-1);
      unsigned long limbs = n - (amount >> FLINT_LG_BITS_PER_LIMB);

      if (bits)
      {
         // each butterfly needs a bitshift
         if (limbs != n)
         {
            ZmodF_sub_mul_Bexp(*scratch, x[2], x[0], limbs, n);
            ZmodF_add(x[0], x[0], x[2], n);
            ZmodF_short_div_2exp(x[2], *scratch, bits, n);

            ZmodF_sub_mul_Bexp(*scratch, x[3], x[1], limbs, n);
            ZmodF_add(x[1], x[1], x[3], n);
            ZmodF_short_div_2exp(x[3], *scratch, bits, n);
         }
         else
         {
            ZmodF_sub(*scratch, x[0], x[2], n);
            ZmodF_add(x[0], x[0], x[2], n);
            ZmodF_short_div_2exp(x[2], *scratch, bits, n);

            ZmodF_sub(*scratch, x[1], x[3], n);
            ZmodF_add(x[1], x[1], x[3], n);
            ZmodF_short_div_2exp(x[3], *scratch, bits, n);
         }
      }
      else
      {
         // no bitshifts required
         ZmodF_sub_mul_Bexp(*scratch, x[2], x[0], limbs, n);
         ZmodF_add(x[0], x[0], x[2], n);
         ZmodF_swap(x+2, scratch);

         ZmodF_sub_mul_Bexp(*scratch, x[3], x[1], limbs, n);
         ZmodF_add(x[1], x[1], x[3], n);
         ZmodF_swap(x+3, scratch);
      }

      return;
   }

   if (depth <= 1)
   {
      if (depth == 1)
      {
         ZmodF_forward_butterfly_2exp(x+1, x, scratch, twist, n);
         ZmodF_swap(x, x+1);
      }
      return;
   }

   unsigned long half = 1 << (depth - 1);
   
   // =========================================================================
   // recurse into two halves

   _ZmodF_poly_IFFT_dual_recursive(x, depth-1, twist, root << 1, n, scratch);
   _ZmodF_poly_IFFT_dual_recursive(x + half, depth-1, twist + root, root << 1,
                                  n, scratch);

   // =========================================================================
   // general case for length >= 8

   // butterflies (a, b) -> (a + b, w*(a - b)), where w = 2^(-amount).
   ZmodF_t* y = x + half;
   unsigned long amount = twist << (depth - 1);
   unsigned long bits = amount & (FLINT_BITS-1);
   unsigned long limbs = n - (amount >> FLINT_LG_BITS_PER_LIMB);

   if (bits)
   {
      // each butterfly needs a bitshift
      if (limbs != n)
      {
         unsigned long i;
         for (i = 0; i < half; i++)
         {
            ZmodF_sub_mul_Bexp(*scratch, y[i], x[i], limbs, n);
            ZmodF_add(x[i], x[i], y[i], n);
            ZmodF_short_div_2exp(y[i], *scratch, bits, n);
         }
      }
      else
      {
         unsigned long i;
         for (i = 0; i < half; i++)
         {
            ZmodF_sub(*scratch, x[i], y[i], n);
            ZmodF_add(x[i], x[i], y[i], n);
            ZmodF_short_div_2exp(y[i], *scratch, bits, n);
         }
      }
   }
   else
   {
      // all butterflies are limbshifts only
      unsigned long i;
      for (i = 0; i < half; i++)
      {
         ZmodF_sub_mul_Bexp(*scratch, y[i], x[i], limbs, n);
         ZmodF_add(x[i], x[i], y[i], n);
         ZmodF_swap(y+i, scratch);
      }
   }
}


/****************************************************************************

   Fourier Transform Routines

****************************************************************************/


void ZmodF_poly_FFT(ZmodF_poly_t poly, unsigned long length)
{
   FLINT_ASSERT(length <= (1UL << poly->depth));
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         unsigned long i;
         for (i = 0; i < length; i++)
            ZmodF_zero(poly->coeffs[i], poly->n);
      }
      else
      {
         if (poly->depth >= 1)
            _ZmodF_poly_FFT(poly->coeffs, poly->depth, 1, poly->length,
                           length, 0, poly->n, poly->scratch);
      }
   }

   poly->length = length;
}


void ZmodF_poly_IFFT(ZmodF_poly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (poly->length && poly->depth)
      _ZmodF_poly_IFFT(poly->coeffs, poly->depth, 1, poly->length,
                      poly->length, 0, 0, poly->n, poly->scratch);
}


// res may alias x or y
// x and y may alias each other
void ZmodF_poly_convolution(ZmodF_poly_t res, ZmodF_poly_t x, ZmodF_poly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodF_poly_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodF_poly_FFT(y, length);
      
   ZmodF_poly_pointwise_mul(res, x, y);
   ZmodF_poly_IFFT(res);
   ZmodF_poly_rescale(res);
}

// res may alias x or y
// x and y may alias each other
// only computes the coefficients in the range 
// [start, n) the rest are rubbish
void ZmodF_poly_convolution_range(ZmodF_poly_t res, ZmodF_poly_t x, 
                                        ZmodF_poly_t y, unsigned long start, unsigned long n)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodF_poly_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodF_poly_FFT(y, length);
      
   ZmodF_poly_pointwise_mul(res, x, y);
   ZmodF_poly_IFFT(res);
   ZmodF_poly_rescale_range(res, start, n);
}


/****************************************************************************

   Negacyclic Fourier Transform Routines
   
****************************************************************************/


/*
ignores length of poly
*/
void ZmodF_poly_negacyclic_FFT(ZmodF_poly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   unsigned long twist = (poly->n * FLINT_BITS) >> poly->depth;

   _ZmodF_poly_FFT_dual_recursive(poly->coeffs, poly->depth, twist, 2*twist, poly->n, poly->scratch);
   poly->length = 1 << poly->depth;
}


void ZmodF_poly_negacyclic_IFFT(ZmodF_poly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   unsigned long twist = (poly->n * FLINT_BITS) >> poly->depth;
   _ZmodF_poly_IFFT_dual_recursive(poly->coeffs, poly->depth, twist, 2*twist, poly->n, poly->scratch);
   poly->length = 1 << poly->depth;
}


void ZmodF_poly_negacyclic_convolution(ZmodF_poly_t res,
                                      ZmodF_poly_t x, ZmodF_poly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long size = 1UL << res->depth;
   
   ZmodF_poly_negacyclic_FFT(x);
   if (x != y)    // take care of aliasing
      ZmodF_poly_negacyclic_FFT(y);
      
   ZmodF_poly_pointwise_mul(res, x, y);
   ZmodF_poly_negacyclic_IFFT(res);
   ZmodF_poly_rescale(res);
   res->length = size;
}


// end of file ****************************************************************
