/******************************************************************************

Development code for Schonhage-Strassen polynomial multiplication over Z.

 (C) 2006 William Hart and David Harvey


This implementation of SS multiplication uses the following optimisations.
These optimisations only make a difference for fairly small coefficients.
For large coefficients, the running time is completely dominated by the
pointwise multiplications anyway.

* Coefficients are represented as raw arrays of limbs. Arithmetic is performed
mainly with GMP's mpn layer. We don't use mpz_t; they are too slow.

* We reduce mod p as rarely as possible. Instead, we let bits build up in the
last limb.

* Butterflies are accomplished using as few passes over the data as possible.
In particular, for those operations that can't be performed in-place, we use
indirection techniques to make sure that we don't have to waste time moving
the data back into the "correct" place.

* We use a cache-friendly algorithm for large transform lengths. It
decomposes the transform into many much smaller transforms which utilise the
cache much better. We also supply some cache hints (for compilers/processors
that understand them).

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gmp.h>
#include "flint.h"
#include "Z-ssmul.h"
#include "mpn_extras.h"
#include "ssfft.h"
#include "flint-manager.h"


/*
This flag is just here temporarily while we evaluate the new cache-friendly
FFT code. If set, it uses the new code.
 */
#define USE_NEW_FFT 1

/*============================================================================

Various low-level conversion functions for operating on the coefficients of 
our FFT.

B = number of bits per limb = FLINT_BITS_PER_LIMB.

We avoid using mpz_t's to represent coefficients whilst doing the FFT, IFFT,
pointwise multiplications and scaling because of the added overhead of 
managing the sign and the buffer length.

=============================================================================*/

void Z_convert_raw_to_mpz(mpz_t output, mp_limb_t* input, unsigned long n)
{
   mpz_import(output, n, -1, sizeof(mp_limb_t), 0, 0, input);
}

/*==============================================================================

Arithmetic Support Functions for the FFT and IFFT. 

Everywhere we are working modulo a fermat number p = 2^(nB) + 1, where n >= 1.

All coefficients are stored in arrays of limbs of length n+1 (least significant
limb first). The limbs represent a single signed integer in 2's complement;
the high bit of x[n] is the sign bit.

In the FFT, we try to do reductions mod p as rarely as possible; we just let
bits build up in the extra limb as long as we can get away with it. To ensure
correctness, each function's documentation describes exactly how much of the
last limb can get chewed up.

===============================================================================*/

/*
Divides a coefficient by 2^s mod p. The shift s MUST be in the range 0 < s < B;
i.e. this is a rotation by a fraction of a limb.

src and dest can be the same buffer. If they're not, they must be disjoint.

Total bitlength of coefficient is not increased. That is, if
    0 <= abs(src) < 2^(nB + k),
then also
    0 <= abs(dest) < 2^(nB + k).

*/
void Z_rotate_right_bits(mp_limb_t* dest, mp_limb_t* src,
                       unsigned long s, unsigned long n)
{
   mp_limb_t high = src[n];
   mp_limb_t overflow = mpn_rshift(dest, src, n + 1, s);

   if ((mp_limb_signed_t)(high) < 0)
      // If src was negative, we have to add extra 1's at the top to
      // make the shifted version negative (because mpn_rshift doesn't
      // respect sign)
      dest[n] += (mp_limb_signed_t)(-1L) << (FLINT_BITS_PER_LIMB - s);
   
   // rotate the overflow back around to the top limb
   mpn_sub_1(dest + n - 1, dest + n - 1, 2, overflow);
}


/*
Reduces the input into the canonical range 0 <= x < p. Result is inplace.

In most cases, this involves simply distributing the top limb down to the
bottom.

In rare cases, it will need to propagate a carry to the next limb. (How rare
depends on how full the top limb is. The probability will be lower for
smaller FFT transform lengths.)

In extremely rare cases it will need to propagate the carry further, possibly
doing a whole pass over the data.

In pathologically rare cases it will perform TWO passes over the data. For
example if the input is exactly p-1, it will first subtract p to obtain -1,
and then add p back to get p-1.

*/
void Z_reduce_mod_p_exact(mp_limb_t* x, unsigned long n)
{
   mp_limb_t hi = x[n];

   if ((mp_limb_signed_t) hi < 0)
   {
      // If top limb (hi) is negative, we add -hi multiples of p
      x[n] = 0;
      mpn_add_1(x, x, n + 1, -hi);

      // If the result is >= p (very unlikely)...
      if (x[n] && x[0])
      {
         // ... need to subtract off p.
         x[n] = 0;
         x[0]--;
      }
   }
   else
   {
      // If top limb (hi) is non-negative, we subtract hi multiples of p
      x[n] = 0;
      mpn_sub_1(x, x, n + 1, hi);

      // If the result is negative (very unlikely)...
      if (x[n])
      {
         // ... need to add back p.
         x[n] = 0;
         mpn_add_1(x, x, n + 1, 1);
      }
   }
}


/*
Computes dest = src * 2^(Bs) mod p.

i.e. shifts left by a whole number of limbs.

MUST have 0 <= s < n.

dest and src may NOT overlap.

The top limb is distributed down, so the result is guaranteed to satisfy
    -2^B + 1 <= dest <= 2^(nB) + 2^B - 2
In other words it doesn't use more than one bit of the top limb.

*/
void Z_rotate_mod_p_limbs(mp_limb_t* dest, mp_limb_t* src,
                        unsigned long s, unsigned long n)
{
   // put negative of high limbs of src into low limbs of dest
   negate_limbs(dest, src + n - s, s + 1);
   mp_limb_t carry = dest[s];

   // put low limbs of src into high limbs of dest
   copy_limbs(dest + s, src, n - s);
   dest[n] = 0;

   // propagate carry
   signed_add_1(dest + s, n - s + 1, carry);
}


/*
Computes dest = src * 2^s mod p.

i.e. shifts left by s bits.

MUST have 0 <= s < Bn.

dest and src may NOT overlap.

The top limb is distributed down, so the result is guaranteed to satisfy
    -2^B + 1 <= dest <= 2^(nB) + 2^B - 2
In other words it doesn't use more than one bit of the top limb.

*/
void Z_rotate_mod_p_bits(mp_limb_t* dest, mp_limb_t* src,
                       unsigned long s, unsigned long n)
{
   // split rotation into a whole number of limbs and leftover number of bits

   // dirty: word size must be a power of two
   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
    
   if (bits == 0)
   {
      // need special case here for bits == 0 because Z_rotate_right_bits
      // does not allow shift count to be zero
      Z_rotate_mod_p_limbs(dest, src, s, n);
      return;
   }

   // Instead of shifting left by s limbs and bits bits, we will
   // shift left by s+1 limbs and RIGHT by (B - bits) bits.
   s++;

   // put the negative of the high limbs of src into the low limbs of dest
   negate_limbs(dest, src + n - s, s + 1);
   mp_limb_t carry = dest[s];
   
   // put the low limbs of src into the high limbs of dest
   if (n != s)
      copy_limbs(dest + s, src, n - s);
   dest[n] = 0;
   
   // propagate carry
   signed_add_1(dest + s, n - s + 1, carry);
   
   // do the fractional limb rotation
   Z_rotate_right_bits(dest, dest, FLINT_BITS_PER_LIMB - bits, n);
}


/*============================================================================

FFT butterflies.

We have separate functions for doing butterflies with rotations by an
arbitrary number of bits and by a whole number of limbs; this saves on some
arithmetic and branching in the deeper recursion levels.

All of these functions operate on raw coefficient buffers (n+1 limbs as
described above). However, they accept POINTERS to mp_limb_t*. The reason is
that it might not be efficient to store the output in the same buffers as
the input, so the functions are permitted to store one of the results in the
scratch buffer, and then permute the pointers so that *a is the first output
and *b is the second output (and *scratch points to whatever the third
remaining buffer is).

============================================================================*/

/* 
Swaps arrays of limbs efficiently

*/
inline void Z_swap_limb_ptrs(mp_limb_t** x, mp_limb_t** y)
{
    mp_limb_t* temp;
    temp = *x;
    *x = *y;
    *y = temp;
}


/*
Computes  (a, b) -> (a + b, 2^(Bs) (a - b)).

i.e. performs a forward FFT butterfly with a rotation by a whole number
of limbs.

MUST have 0 < s < n.

Suppose the inputs satisfy
   0 <= abs(a) < 2^(nB + k),
   0 <= abs(b) < 2^(nB + k),
for some integer 0 <= k < B-1.

Then the outputs have the following guarantees:
   0 <= abs(a) < 2^(nB + k + 1),
   0 <= abs(b) < 2^(nB + 1).

In other words, a is at most one bit bigger than either of the inputs, and b
does not use more than one bit of the last limb.

*/
void Z_fft_butterfly_limbs(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                         unsigned long s, unsigned long n)
{
   // get low limbs of a - b into high limbs of output
   (*scratch)[n] = -mpn_sub_n(*scratch + s, *a, *b, n - s);
    
   // get high limbs of b - a into low limbs of output
   mp_limb_t overflow =
         (*b)[n] - (*a)[n] - mpn_sub_n(*scratch, *b + n - s, *a + n - s, s);
   signed_add_1(*scratch + s, n + 1 - s, overflow);
    
   // a = a + b
   mpn_add_n(*a, *a, *b, n+1);
    
   // swap rotated version of a - b back into b
   Z_swap_limb_ptrs(b, scratch);
}


/*
Computes  (a, b) -> (a + b, 2^s (a - b)).

i.e. performs a forward FFT butterfly with a rotation by an arbitrary number
of bits.

MUST have 0 < s < Bn.

Output guarantees are the same as for Z_fft_butterfly_limbs.

*/
void Z_fft_butterfly_bits(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                        unsigned long s, unsigned long n)
{
   // split rotation into a whole number of limbs and leftover number of bits

   // dirty: word size must be a power of two
   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
    
   if (bits == 0)
   {
      // need special case here for bits == 0 because Z_rotate_right_bits
      // does not allow shift count to be zero
      Z_fft_butterfly_limbs(a, b, scratch, s, n);
      return;
   }
   
   // Instead of shifting left by s limbs and bits bits, we will
   // shift left by s+1 limbs and RIGHT by (B - bits) bits.
   s++;
   
   // First do the subtraction and whole-limb rotation in a single pass.
   // We use pretty much the same code sequence as in Z_fft_butterfly_limbs,
   // except we also need to cover the case s == n.
   if (s == n)
      (*scratch)[n] = 0;
   else
      (*scratch)[n] = -mpn_sub_n(*scratch + s, *a, *b, n - s);
    
   mp_limb_t overflow =
           (*b)[n] - (*a)[n] - mpn_sub_n(*scratch, *b + n - s, *a + n - s, s);
   signed_add_1(*scratch + s, n + 1 - s, overflow);
   
   // a = a + b
   mpn_add_n(*a, *a, *b, n + 1);
   
   // Now do the fractional-limb rotation
   Z_rotate_right_bits(*b, *scratch, FLINT_BITS_PER_LIMB - bits, n);
}


/*
Computes  (a, b) -> (a - 2^(Bs) b, a + 2^(Bs) b).

i.e. performs an inverse FFT butterfly with a rotation by a whole number
of limbs.

MUST have 0 < s < n.

The output guarantees for this function are very tight: the high limbs of
both outputs differ by at most 2 from the high limb of the input a.

*/
void Z_ifft_butterfly_limbs(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                          unsigned long s, unsigned long n)
{
   // add low limbs of b to high limbs of a
   (*scratch)[n] = (*a)[n] + mpn_add_n(*scratch + s, *a + s, *b, n - s);
    
   // subtract high limbs of b from low limbs of a
   mp_limb_t overflow = -(*b)[n] - mpn_sub_n(*scratch, *a, *b + n - s, s);
   signed_add_1(*scratch + s, n + 1 - s, overflow);
   
   // subtract low limbs of b from high limbs of a
   (*a)[n] -= mpn_sub_n(*a + s, *a + s, *b, n - s);
   
   // add low limbs of a to high limbs of b
   overflow = (*b)[n] + mpn_add_n(*a, *a, *b + n - s, s);
   signed_add_1(*a + s, n + 1 - s, overflow);
   
   // swap answer into b
   Z_swap_limb_ptrs(b, scratch);
}


/*
Computes  (a, b) -> (a - 2^s b, a + 2^s b).

i.e. performs an inverse FFT butterfly with a rotation by an arbitrary number
of bits.

MUST have 0 < s < nB.

Output guarantees are the same as for Z_ifft_butterfly_limbs.

*/
void Z_ifft_butterfly_bits(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                         unsigned long s, unsigned long n)
{
   // split shift into a whole number of limbs and leftover number of bits

   // dirty: word size must be a power of two
   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
   
   if (bits == 0)
   {
      // need special case here for bits == 0 because Z_rotate_right_bits
      // does not allow shift count to be zero
      Z_ifft_butterfly_limbs(a, b, scratch, s, n);
      return;
   }
   
   // Instead of shifting left by s limbs and bits bits, we will
   // shift left by s+1 limbs and RIGHT by (B - bits) bits.
   s++;
   
   // First do the fractional-limb rotation.
   Z_rotate_right_bits(*scratch, *b, FLINT_BITS_PER_LIMB - bits, n);
   
   // Now basically use the same code sequence as in Z_ifft_butterfly_limbs,
   // except we also need to cover the case s == n.
   
   if (n == s)
      (*b)[n] = (*a)[n];
   else
      (*b)[n] = (*a)[n] + mpn_add_n(*b + s, *a + s, *scratch, n - s);
   
   mp_limb_t overflow =
               -(*scratch)[n] - mpn_sub_n(*b, *a, *scratch + n - s, s);
   signed_add_1(*b + s, n + 1 - s, overflow);
   
   if (n != s)
      (*a)[n] -= mpn_sub_n(*a + s, *a + s, *scratch, n - s);
   
   overflow = (*scratch)[n] + mpn_add_n(*a, *a, *scratch + n - s, s);
   signed_add_1(*a + s, n + 1 - s, overflow);
}


/*============================================================================

Main FFT functions. These functions call the butterflies in the right order to
perform FFTs and IFFTs.

=============================================================================*/

#if USE_NEW_FFT

/*
Basecase of the FFT. Once the FFT has been broken up into small enough pieces
this function actually performs the FFT recursively for each of those pieces.

*/
void Z_fft_recursive(mp_limb_t** start, unsigned long skip,
                   unsigned long start_r, unsigned long skip_r,
                   unsigned long depth, mp_limb_t** scratch,
                   unsigned long n, int first)
{
   int fractional_limb_shift = (skip_r | start_r) & (FLINT_BITS_PER_LIMB - 1);

   unsigned long half = skip << (depth - 1);
   mp_limb_t** middle = start + half;
   unsigned long offset, r;

   if (first)
   {
      if (fractional_limb_shift)
      {
         for (offset = 0, r = start_r; offset < half;
              offset += skip, r += skip_r)
            Z_rotate_mod_p_bits(middle[offset], start[offset], r, n);
      }
      else
      {
         unsigned long skip_r_limbs = skip_r / FLINT_BITS_PER_LIMB;

         for (offset = 0, r = start_r / FLINT_BITS_PER_LIMB; offset < half;
              offset += skip, r += skip_r_limbs)
         {
            Z_rotate_mod_p_limbs(middle[offset], start[offset], r, n);
         }
      }

      if (depth == 1)
         return;
   }
   else
   {
      if (start_r == 0)
      {
         // scratch = a - b
         mpn_sub_n(*scratch, *start, *middle, n + 1);
         
         // a += b
         mpn_add_n(*start, *start, *middle, n + 1);
         
         Z_swap_limb_ptrs(middle, scratch);
      }
      else
         Z_fft_butterfly_bits(start, middle, scratch, start_r, n);

      if (depth == 1)
         return;
      
      if (fractional_limb_shift)
      {
         for (offset = skip, r = start_r + skip_r; offset < half;
              offset += skip, r += skip_r)
         {
             Z_fft_butterfly_bits(start + offset, middle + offset, scratch, r, n);
         } 

      }
      else
      {
         unsigned long skip_r_limbs = skip_r / FLINT_BITS_PER_LIMB;
         unsigned long start_r_limbs = start_r / FLINT_BITS_PER_LIMB;
         
         for (offset = skip, r = start_r_limbs + skip_r_limbs; offset < half;
              offset += skip, r += skip_r_limbs)
         {
             Z_fft_butterfly_limbs(start + offset, middle + offset,
                                 scratch, r, n);
         }
      }
   }

   Z_fft_recursive(start, skip, start_r << 1, skip_r << 1, depth - 1,
                 scratch, n, 0);
   Z_fft_recursive(middle, skip, start_r << 1, skip_r << 1, depth - 1,
                 scratch, n, 0);
}


/*
transform of length 2^depth

coefficients are at start, start + skip, ... start + (length-1)*skip

corresponding values of r are:
 start_r, start_r + skip_r, ... start_r + (length/2 - 1)*skip_r

the idea of this routine is to get the transforms down to fit in cache
cache if possible, then switch to a simple recursive algorithm
 
*/
void Z_fft_main(mp_limb_t** start, unsigned long skip,
              unsigned long start_r, unsigned long skip_r,
              unsigned long depth, mp_limb_t** scratch,
              unsigned long n, int first, int crossover)
{
   if (crossover == -1)
   {
      // work out crossover = optimal depth to switch over to plain
      // recursive algorithm
      unsigned long cache_length =
                        2*FLINT_CACHE_SIZE / ((n+1) * sizeof(mp_limb_t)) - 1;
      for (crossover = 0; cache_length > 0; cache_length >>= 1, crossover++);
      if (crossover == 0)
         crossover = 1;
   }

   if (depth <= crossover)
   {
      Z_fft_recursive(start, skip, start_r, skip_r, depth, scratch, n, first);
      return;
   }

   // Factor FFT into two pieces.
   // Write depth = depth1 + depth2.
   // First piece does the first depth1 layers.
   // Second piece does the last depth2 layers.
   
   // In other words we do
   // length2 transforms of length length1, followed by
   // length1 transforms of length length2,
   // where length1 * length2 = 2^depth.

   unsigned long depth2 = crossover;
   unsigned long depth1 = depth - depth2;

   unsigned long length1 = 1 << depth1;
   unsigned long length2 = 1 << depth2;

   mp_limb_t** next_start = start;
   unsigned long next_skip = skip << depth2;
   unsigned long next_start_r = start_r;
   unsigned long next_skip_r = skip_r << depth2;

   for (unsigned long i = 0; i < length2; i++, next_start += skip,
        next_start_r += skip_r)
      Z_fft_main(next_start, next_skip, next_start_r, next_skip_r,
               depth1, scratch, n, first, crossover);

   next_start = start;
   next_start_r = start_r << depth1;
   next_skip_r = skip_r << depth1;

   for (unsigned long i = 0; i < length1; i++, next_start += next_skip)
      Z_fft_main(next_start, skip, next_start_r, next_skip_r,
               depth2, scratch, n, 0, crossover);
}

/*
Basecase for the IFFT. Just the inverse of Z_fft_recursive.

*/
void Z_ifft_recursive(mp_limb_t** start, unsigned long skip,
                    unsigned long start_r, unsigned long skip_r,
                    unsigned long depth, mp_limb_t** scratch,
                    unsigned long n)
{
   int fractional_limb_shift = (skip_r | start_r) & (FLINT_BITS_PER_LIMB - 1);

   unsigned long half = skip << (depth - 1);
   mp_limb_t** middle = start + half;
   unsigned long offset, r;

   if (depth > 1)
   {
      Z_ifft_recursive(start, skip, start_r << 1, skip_r << 1, depth - 1,
                     scratch, n);
      Z_ifft_recursive(middle, skip, start_r << 1, skip_r << 1, depth - 1,
                     scratch, n);
   }

   if (start_r == 0)
   {
      // scratch = a - b
      mpn_sub_n(*scratch, *start, *middle, n + 1);
      
      // a += b
      mpn_add_n(*start, *start, *middle, n + 1);
      
      Z_swap_limb_ptrs(middle, scratch);
   }
   else
      Z_ifft_butterfly_bits(start, middle, scratch, 
           n*FLINT_BITS_PER_LIMB - start_r, n);
   
   if (depth == 1)
      return;
   
   if (fractional_limb_shift)
   {
      for (offset = skip, r = n*FLINT_BITS_PER_LIMB - start_r - 
           skip_r; offset < half; offset += skip, r -= skip_r)
      {
          Z_ifft_butterfly_bits(start + offset, middle + offset, scratch, r, n);
      }
   }
   else
   {
      unsigned long skip_r_limbs = skip_r / FLINT_BITS_PER_LIMB;
      unsigned long start_r_limbs = start_r / FLINT_BITS_PER_LIMB;
      
      for (offset = skip, r = n - start_r_limbs - skip_r_limbs;
             offset < half; offset += skip, r -= skip_r_limbs)
      {
         Z_ifft_butterfly_limbs(start + offset, middle + offset,
                              scratch, r, n);
      }
   }
}

/*
Main IFFT routine. It breaks the IFFT up into pieces and calls 
Z_ifft_recursive on the pieces. It is the inverse of Z_fft_main.

*/
void Z_ifft_main(mp_limb_t** start, unsigned long skip,
               unsigned long start_r, unsigned long skip_r,
               unsigned long depth, mp_limb_t** scratch,
               unsigned long n, int crossover)
{
   if (crossover == -1)
   {
      // work out crossover = optimal depth to switch over to plain
      // recursive algorithm
      unsigned long cache_length =
                        2*FLINT_CACHE_SIZE / ((n+1) * sizeof(mp_limb_t)) - 1;
      for (crossover = 0; cache_length > 0; cache_length >>= 1, crossover++);
      if (crossover == 0)
         crossover = 1;
   }

   // do the forward FFT backwards

   if (depth <= crossover)
   {
      Z_ifft_recursive(start, skip, start_r, skip_r, depth, scratch, n);
      return;
   }

   unsigned long depth2 = crossover;
   unsigned long depth1 = depth - depth2;

   unsigned long length1 = 1 << depth1;
   unsigned long length2 = 1 << depth2;

   mp_limb_t** next_start = start;
   unsigned long next_skip = skip << depth2;
   unsigned long next_start_r = start_r << depth1;
   unsigned long next_skip_r = skip_r << depth1;

   for (unsigned long i = 0; i < length1; i++, next_start += next_skip)
      Z_ifft_main(next_start, skip, next_start_r, next_skip_r,
                depth2, scratch, n, crossover);

   next_start = start;
   next_start_r = start_r;
   next_skip_r = skip_r << depth2;

   for (unsigned long i = 0; i < length2; i++, next_start += skip,
        next_start_r += skip_r)
   {
      Z_ifft_main(next_start, next_skip, next_start_r, next_skip_r,
                depth1, scratch, n, crossover);
   }
}


#else

// todo: remove these older FFT functions at some point,
// if we're convinced the new ones are unconditionally better

/* we are doing a recursive, inplace transform of length 2^depth,
twiddling by 2^r. Output is in bit-reversed order as usual.
data[i*(n+1)] is the beginning of the i-th coefficient.

r should be 0 <= r < B*n

*/

void Z_fft(mp_limb_t** start, unsigned long length, mp_limb_t** scratch,
         unsigned long r, unsigned long n, int first)
{
   unsigned long half = length >> 1;
   mp_limb_t** middle = start + half;
   unsigned long offset, r_mult;
    
   if (first)
   {
      copy_limbs(middle[0], start[0], n + 1);
      
      if (r & (FLINT_BITS_PER_LIMB - 1))
      {
         // not shifting by a multiple of the limb length
         for (offset = 1, r_mult = r; offset < half; offset++, r_mult += r)
            Z_rotate_mod_p_bits(middle[offset], start[offset], r_mult, n);
      }
      else
      {
         // shifting by a multiple of the limb length
         unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
         for (offset = 1, r_mult = rdivB; offset < half;
              offset++, r_mult += rdivB)
            Z_rotate_mod_p_limbs(middle[offset], start[offset], r_mult, n);
      }
   }
   else
   {
      // scratch = a - b
      mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
      
      // a += b
      mpn_add_n(start[0], start[0], middle[0], n + 1);
      
      Z_swap_limb_ptrs(middle, scratch);
      
      if (half == 1)
         return;
      
      if (r & (FLINT_BITS_PER_LIMB - 1))
      {
         // not shifting by a multiple of the limb length
         for (offset = 1, r_mult = r; offset < half; offset++, r_mult += r)
            Z_fft_butterfly_bits(start + offset, middle + offset,
                               scratch, r_mult, n);
      }
      else
      {
         // shifting by a multiple of the limb length
         unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
         for (offset = 1, r_mult = rdivB; offset < half;
              offset++, r_mult += rdivB)
            Z_fft_butterfly_limbs(start + offset, middle + offset,
                                scratch, r_mult, n);
      }
   }
   
   Z_fft(middle, length >> 1, scratch, r << 1, n, 0);
   Z_fft(start, length >> 1, scratch, r << 1, n, 0);
}

/*
Inverse of Z_fft.

*/
void Z_ifft(mp_limb_t** start, unsigned long length, mp_limb_t** scratch,
          unsigned long r, unsigned long n)
{
   unsigned long half = length >> 1;
   mp_limb_t** middle = start + half;
   unsigned long offset, r_mult = 0;
   unsigned long nB = FLINT_BITS_PER_LIMB * n;
   
   if (half == 1)
   {
      // scratch = a - b
      mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
      
      // a += b
      mpn_add_n(start[0], start[0], middle[0], n + 1);
      
      Z_swap_limb_ptrs(middle, scratch);
      
      return;
   }
   
   Z_ifft(middle, half, scratch, r << 1, n);
   Z_ifft(start, half, scratch, r << 1, n);
   
   // scratch = a - b
   mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
   
   // a += b
   mpn_add_n(start[0], start[0], middle[0], n + 1);
   
   Z_swap_limb_ptrs(middle, scratch);
   
   if (r & (FLINT_BITS_PER_LIMB - 1))
   {
      // not shifting by a multiple of the limb length
      for (offset = 1, r_mult = nB - r; offset < half; offset++, r_mult -= r)
      {
         Z_ifft_butterfly_bits(start + offset, middle + offset,
                             scratch, r_mult, n);
      }
   }
   else
   {
      // shifting by a multiple of the limb length
      unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
      for (offset = 1, r_mult = n - rdivB; offset < half;
           offset++, r_mult -= rdivB)
         Z_ifft_butterfly_limbs(start + offset, middle + offset,
                              scratch, r_mult, n);
   }
}

#endif

void Z_ssmul_main(mp_limb_t** array1, unsigned long length1,
                 mp_limb_t** array2, unsigned long length2,
                 mp_limb_t** scratch, unsigned long log_length,
                 unsigned long n)
{
    unsigned long i, j;

    // root of unity is sqrt2^r
    // todo: add assertion here to check divisibility by 2^log_length
    // todo: change this to shift by FLINT_LG_BITS_PER_LIMB
    unsigned long r = (4*n*FLINT_BITS_PER_LIMB) >> log_length;

    // number of fourier coefficients required (i.e. length of output)
    // todo: add assertion here to check output_length <= 2^log_length
    unsigned long output_length = length1 + length2 - 1;

    // do FFTs

    ssfft_fft(array1, 1, log_length, length1, output_length, 0, r, n, scratch);
    ssfft_fft(array2, 1, log_length, length2, output_length, 0, r, n, scratch);

    // pointwise multiplies

    mp_limb_t* mul_scratch = (mp_limb_t*) flint_malloc_limbs(2*n);

    for (i = 0; i < output_length; i++)
    {
       for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
       for (j = 0; j < n; j += 8) FLINT_PREFETCH(array2[i+1], j);

       Z_reduce_mod_p_exact(array1[i], n);
       Z_reduce_mod_p_exact(array2[i], n);

       // We've reduced the coefficients to the range 0 <= x < p.
       // If they are < p-1, then they fit into n limbs. Otherwise we have
       // special cases as dealt with in the following code.

       if (array1[i][n])
       {
          // special case: array1[i] == -1 (mod p),
          // so answer is (-1) * array2[i]
          negate_limbs(array1[i], array2[i], n+1);
          Z_reduce_mod_p_exact(array1[i], n);
       }
       else if (array2[i][n])
       {
          // special case: array2[i] == -1 (mod p),
          // so answer is (-1) * array1[i]
          negate_limbs(array1[i], array1[i], n+1);
          Z_reduce_mod_p_exact(array1[i], n);
       }
       else
       {
          // usual case, just do the multiplication
          if (n >= 1796) Z_SSMul(mul_scratch, array1[i], array2[i], n, n);
          else mpn_mul_n(mul_scratch, array1[i], array2[i], n);
          // reduce the result mod p
          array1[i][n] = -mpn_sub_n(array1[i], mul_scratch, mul_scratch + n, n);
       }
    }

    flint_free(mul_scratch);

    // do inverse FFT
    ssfft_ifft(array1, 1, log_length, output_length, output_length, 0,
               0, r, n, scratch);
}


/*============================================================================

Main Schonhage-Strassen routine exported from this module.

=============================================================================*/

// todo: write an SSSqr for squaring instead of multiplying

/* 
Main polynomial multiplication routine. It performs a Schoenhage-Strassen
algorithm, calling the above Z_fft and Z_ifft functions. res is set to the 
polynomial data1*data2.

At approximately the
point when the degree of the polynomials exceeds the input coefficient bit
size, it starts bundling coefficients together using the Kronecker-Schoenhage
technique as refined by Paul Zimmerman and David Harvey.

log_length is the ceiling of the log_2 of the length of the *input* 
polynomial (length = deg+1).  coeff_bits is the number of bits in the 
maximimum *input* coefficient. We assume that the length of the output
polynomial will always be less than 2^B where B = FLINT_BITS_PER_LIMB

todo: add code to accept Zpoly's and compute the other two parameters
from the given poly inputs

*/
void Z_SSMul(mp_limb_t* res, mp_limb_t* data1, mp_limb_t* data2,
              unsigned long limbs, unsigned long limbs2)
{
   unsigned long twk;
   
   twk = 2;   
      
   unsigned long i, j, skip;
   unsigned long input_limbs = limbs;
   
   unsigned long coeff_bits = limbs*FLINT_BITS_PER_LIMB;
   // number of bits that output coefficients would normally require:
   unsigned long output_bits = 2*coeff_bits;
   
   unsigned long log_length = 0; 
   unsigned long length = 1;
   
   while ((twk*length <= 2*output_bits)&&(input_limbs>1)) 
   {
         length<<=1;
         log_length++;
         input_limbs = (input_limbs+1)>>1;
         output_bits = 2*FLINT_BITS_PER_LIMB*input_limbs+log_length;
   }
   // round up to a multiple of the bits that will support the convolution length
   output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   coeff_bits = (output_bits-log_length)/2;
   input_limbs = (coeff_bits)/FLINT_BITS_PER_LIMB;
   coeff_bits = FLINT_BITS_PER_LIMB*input_limbs;
   //printf("%ld %ld %ld %ld \n",coeff_bits,input_limbs,output_bits,log_length);
   
   // make nB larger than the output coefficient size (B = FLINT_BITS_PER_LIMB)
   // this effectively rounds the output coefficients up to a multiple of the 
   // limb length
   unsigned long n = (output_bits - 1) / FLINT_BITS_PER_LIMB + 1;
   
   // now we set length to the *output* polynomial length (deg+1)
   // since the input polynomial length is not needed further
   length = 1 << (log_length + 1);
   unsigned long trunc_length = (limbs-1)/input_limbs + 1;
   unsigned long trunc_length2 = (limbs2-1)/input_limbs + 1;
   
   // We want p = 2^rm+1 to be just larger than the largest output
   // coefficient and where m=2^k is bigger than B = FLINT_BITS_PER_LIMB
   // and in such a way that p will have 2l roots of unity where l is the
   // length of the input polynomials (i.e. l = deg+1 of input polys), so that
   // we can do a convolution of length 2l. In other words: k >= l and 2^k > B. 
   // But the output coefficients are smaller than 2^nB+1 by definition, 
   // and we already ensured nB was a multiple of l, so these two 
   // conditions are satisfied if we set r = nB / (length/2).
   // (Recall l = length/2).
   
   //unsigned long r = n*FLINT_BITS_PER_LIMB * 2 / length;

   // We need the following memory allocated:
   //
   // * Space for 2*length + 1 coefficients of size n+1 limbs
   //   (the extra 1 is for a scratch buffer)
   //
   // * Space for 2*length + 1 pointers to coefficients
   //
   // * Working space for a pointwise multiplication of length 2n
   //   We use n+1 limbs instead of n so that we can allow
   //   carries to build up without reducing mod p all the time

   // dirty: this code assumes that a pointer is the same length as a limb
   // on account of limb_alloc just taking the number of limbs that the caller
   // wants as a parameter
   
   mp_limb_t* array = (mp_limb_t*) flint_malloc_limbs((2*length + 1) * (n+1) + 2*n);
   mp_limb_t** pointarr = (mp_limb_t**) flint_malloc_limbs(2*length + 1);

   // partition up array into coefficient buffers of length n+1...
   
   pointarr[0] = array;
   for (i = 0; i < 2*length; i++)
      pointarr[i+1] = pointarr[i] + (n+1);
      
   // ...and pointwise multiplication working space
   
   mp_limb_t** array1 = pointarr;
   mp_limb_t** array2 = pointarr + length;
   mp_limb_t** scratch = pointarr + 2*length;
   // convert data for first FFT
   for (skip = 0, i = 0; (i < length/2) && (skip+input_limbs <= limbs); i++, skip+=input_limbs)
   {
      clear_limbs(array1[i],n+1);
      // prefetch entire bundled coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
      // convert a coefficient
      copy_limbs(array1[i], data1+skip, input_limbs);
   }
   if (i < length/2) clear_limbs(array1[i],n+1);
   if (limbs>skip) copy_limbs(array1[i], data1+skip, limbs-skip);
   i++;
   for (; i < length/2; i++) clear_limbs(array1[i],n+1);
      
   // convert data for second FFT
   for (skip = 0, i = 0; (i < length/2) && (skip+input_limbs <= limbs2); i++, skip+=input_limbs)
   {
      clear_limbs(array2[i],n+1);
      // prefetch entire bundled coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array2[i+1], j);
      // convert a coefficients
      copy_limbs(array2[i], data2+skip, input_limbs);
   }
   if (i < length/2) clear_limbs(array2[i],n+1);
   if (limbs2>skip) copy_limbs(array2[i], data2+skip, limbs2-skip);
   i++;
   for (; i < length/2; i++) clear_limbs(array2[i],n+1);
     
   Z_ssmul_main(array1, trunc_length, array2, trunc_length2, scratch, log_length+1, n);
   
   // We have to clear the output polynomial, since we have to add to its
   // coefficients rather than just copy into them
   
   clear_limbs(res,2*limbs);
    
   for (skip = 0, i = 0; (i < trunc_length + trunc_length2 - 1) && (skip+n <= 2*limbs); i++, skip+=input_limbs)
   { 
      // divide by appropriate normalising power of 2
      // (I'm assuming here that the transform length will always be 
      // less than 2^B)
      
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
      Z_rotate_right_bits(array1[i], array1[i], log_length+1, n);
      Z_reduce_mod_p_exact(array1[i], n);
      
      // convert all the coefficients out 
      mpn_add(res+skip, res+skip, n+1, array1[i], n);
      
   } 
   
   while (skip < 2*limbs)
   {
      if ((2*limbs > skip) && (i < trunc_length+trunc_length2 - 1))
      {
         for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
         Z_rotate_right_bits(array1[i], array1[i], log_length+1, n);
         Z_reduce_mod_p_exact(array1[i], n);
      
         // convert all the coefficients out 
         mpn_add(res+skip, res+skip, 2*limbs - skip, array1[i], 2*limbs - skip);
      }  
      i++;
      skip+=input_limbs;
   }  
   
   // free allocated space
   flint_free(pointarr);
   flint_free(array);
}

void Z_fast_mul(mpz_t res, mpz_t a, mpz_t b)
{
   unsigned long int limbs;
   if (a->_mp_size > (1796*64)/FLINT_BITS_PER_LIMB)
   {
      if (a->_mp_size >= b->_mp_size) limbs = a->_mp_size; 
      else limbs = b->_mp_size; 
      mp_limb_t* output = (mp_limb_t*) limb_alloc(2*limbs,0);
      if (a->_mp_size > b->_mp_size) Z_SSMul(output, a->_mp_d, b->_mp_d, a->_mp_size, b->_mp_size);
      else Z_SSMul(output, b->_mp_d, a->_mp_d, b->_mp_size, a->_mp_size);
      Z_convert_raw_to_mpz(res,output,2*limbs);
      if (mpz_sgn(res) != mpz_sgn(a)*mpz_sgn(b)) mpz_neg(res,res);
      limb_release();
   } else mpz_mul(res, a, b);
}

// end of file ****************************************************************
