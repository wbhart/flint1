/******************************************************************************

 ZmodF_mul.c

 Copyright (C) 2007, David Harvey and William Hart
 
 Routines for multiplication of elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS.
 
******************************************************************************/

#include <math.h>
#include "ZmodF.h"
#include "ZmodF_poly.h"
#include "ZmodF_mul.h"
#include "mpn_extras.h"
#include "ZmodF_mul-tuning.h"
#include "longlong_wrapper.h"
#include "longlong.h"


/******************************************************************************

   Plain multiplication

******************************************************************************/

/*
   Normalises a and b, and then attempts to multiply them mod p, putting result
   in res. It only succeeds if one of the inputs is exactly -1 mod p, in which
   case it returns 1. Otherwise it just returns 0.
*/
long _ZmodF_mul_handle_minus1(ZmodF_t res, ZmodF_t a, ZmodF_t b,
                              unsigned long n)
{
   ZmodF_normalise(a, n);
   ZmodF_normalise(b, n);
      
   if (a[n])
   {
      // a = -1 mod p
      ZmodF_neg(res, b, n);
      return 1;
   }

   if (b[n])
   {
      // b = -1 mod p
      ZmodF_neg(res, a, n);
      return 1;
   }
   
   return 0;
}


/*
   Same as _ZmodF_mul_handle_minus1, but for squaring.
*/
long _ZmodF_sqr_handle_minus1(ZmodF_t res, ZmodF_t a, unsigned long n)
{
   ZmodF_normalise(a, n);

   if (a[n])
   {
      // a = -1 mod p
      if (a == res)
         // handle aliasing
         res[n] = 0;
      else
         ZmodF_zero(res, n);
      res[0] = 1;

      return 1;
   }

   return 0;
}


/*
   Computes res := a * b mod p = B^n + 1, where a and b are of length n.
   scratch must be of length 2n, and not overlap any of a, b, res.
   a and b must be normalised, and have zero overflow limbs.
   Any combinations of a, b, res may be aliased.
   The output is not necessarily normalised.
*/
inline
void _ZmodF_mul(ZmodF_t res, ZmodF_t a, ZmodF_t b, mp_limb_t* scratch,
                unsigned long n)
{
   FLINT_ASSERT(a[n] == 0);
   FLINT_ASSERT(b[n] == 0);

   // Detect zero limbs at the top of a and b, and reduce multiplication size
   // appropriately.

   // (Note: the reason we do this is that in the FFTs we often *do* get
   // Fourier coefficients with plenty of zeroes; for example the first and
   // second coefficients are x0 + x1 + ... + xM and x0 - x1 + x2 - ... - xM,
   // which are about half the size of the other coefficients. The overhead
   // from performing the checks is negligible if n is large enough, and we
   // expect say n >= 16 if the FFTs are tuned correctly.)

   unsigned long limbs_out = 2*n;
   unsigned long limbs1 = n;
   while (limbs1 && !a[limbs1-1])
   {
      scratch[--limbs_out] = 0;
      limbs1--;
   }
   unsigned long limbs2 = n;
   while (limbs2 && !b[limbs2-1])
   {
      scratch[--limbs_out] = 0;
      limbs2--;
   }
   if ((limbs1 == 0) || (limbs2 == 0)) 
   {
      F_mpn_clear(res, n+1);
      return;
   }
   
   // do the product into scratch
   if (limbs1 >= limbs2)
      F_mpn_mul(scratch, a, limbs1, b, limbs2);
   else
      F_mpn_mul(scratch, b, limbs2, a, limbs1);
      
   // reduce mod p
   res[n] = -mpn_sub_n(res, scratch, scratch + n, n);
}


void ZmodF_mul(ZmodF_t res, ZmodF_t a, ZmodF_t b, mp_limb_t* scratch,
               unsigned long n)
{
   // try special cases a = -1 or b = -1 mod p
   if (_ZmodF_mul_handle_minus1(res, a, b, n))
      return;
      
   // that didn't work, run ordinary multiplication
   _ZmodF_mul(res, a, b, scratch, n);
}


void ZmodF_sqr(ZmodF_t res, ZmodF_t a, mp_limb_t* scratch, unsigned long n)
{
   // try special case a = -1 mod p
   if (_ZmodF_sqr_handle_minus1(res, a, n))
      return;

   // that didn't work, run ordinary multiplication
   _ZmodF_mul(res, a, a, scratch, n);
}


/******************************************************************************

   ZmodF_mul_info initialisation routines

******************************************************************************/

/*
   Initialises info to use plain mpn_mul_n for multiplication.
*/
void ZmodF_mul_info_init_plain(ZmodF_mul_info_t info, unsigned long n,
                               int squaring)
{
   info->n = n;
   info->squaring = squaring;
   info->algo = ZMODF_MUL_ALGO_PLAIN;
   info->scratch = (mp_limb_t*) flint_stack_alloc(2*n);
}


/*
   Initialises info to use the threeway algorithm for multiplication.

   PRECONDITIONS:
      n must be divisible by 3
*/
void ZmodF_mul_info_init_threeway(ZmodF_mul_info_t info, unsigned long n,
                                  int squaring)
{
   FLINT_ASSERT(n % 3 == 0);

   info->n = n;
   info->squaring = squaring;
   info->m = n/3;
   info->algo = ZMODF_MUL_ALGO_THREEWAY;
   // todo: maybe can use less memory here when squaring:
   info->scratch = (mp_limb_t*) flint_stack_alloc(3*n + 1);
}


/*
   Initialises info to use FFT algorithm for multiplication.

   If m == 0, it will automatically choose a suitable value for m, and will
   choose k = 0. (This won't always be optimal, but it should be pretty close.)

   PRECONDITIONS:
      n*FLINT_BITS must be divisible by 2^depth (so that the input can be
      broken up into 2^depth pieces)
      
      0 <= k <= 2
      k <= m
      
      m*FLINT_BITS must be divisible by 2^depth (so that the coefficients
      support a negacyclic FFT of the given depth).
      
      (m+k)*FLINT_BITS >= 2*n*FLINT_BITS/2^depth + 1 + depth (so that getting
      the convolution mod (B^m + 1)*B^k determines the result uniquely)
*/
void ZmodF_mul_info_init_fft(ZmodF_mul_info_t info, unsigned long n,
                             unsigned long depth, unsigned long m,
                             unsigned long k, int squaring)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(k <= m);

   info->algo = ZMODF_MUL_ALGO_FFT;
   info->n = n;
   info->squaring = squaring;

   if (!m)
   {
      // automatically determine reasonable values for m and k
      
      unsigned long input_bits = (n*FLINT_BITS) >> depth;
      unsigned long output_bits = 2*input_bits + 1 + depth;
      m = ((output_bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1;
      k = 0;
      
      // m needs to be divisible by 2^shift
      unsigned long shift = 0;
      if (depth > FLINT_LG_BITS_PER_LIMB)
         shift = depth - FLINT_LG_BITS_PER_LIMB;

      // if m is in the feasible range for the threeway algorithm, try to
      // round it up/down to a multiple of 3
      if (m < ZmodF_mul_threeway_fft_threshold)
      {
         // first try rounding down to a multiple of 3, using k to compensate
         // (this only works if the smaller m satisfies the right divisibility
         // conditions)
         unsigned long smaller_m = (m / 3) * 3;
         if (((smaller_m >> shift) << shift) == smaller_m)
         {
            k = m - smaller_m;
            m = smaller_m;
         }
         else
         {
            // that didn't work; just round up to a multiple of 3*2^shift
            unsigned long round = 3 << shift;
            m = (((m-1) / round) + 1) * round;
         }
      }
      else
      {
         // threeway not feasible.

         // first try rounding down to a multiple of 2^shift, using
         // k to compensate
         unsigned long smaller_m = (m >> shift) << shift;
         if (m - smaller_m <= 2)
         {
            k = m - smaller_m;
            m = smaller_m;
         }
         else
         {
            // that didn't work; just round up to a multiple of 2^shift
            m = (((m - 1) >> shift) + 1) << shift;
         }
      }
   }

   FLINT_ASSERT((m*FLINT_BITS) % (1 << depth) == 0);
   FLINT_ASSERT(k <= 2);
   FLINT_ASSERT((m+k)*FLINT_BITS >= ((2*n*FLINT_BITS) >> depth) + 1 + depth);
   info->m = m;
   info->k = k;

   // For the ZmodF_poly routines, we'll only need m+1 limbs (i.e. to work
   // mod B^m + 1). But later on we'll need m+k+1 limbs to store the result
   // mod (B^m + 1)*B^k, so we allocate these spare limbs in advance.
   ZmodF_poly_init(info->polys[0], depth, m + k, 1);
   ZmodF_poly_decrease_n(info->polys[0], m);
   
   if (!squaring)
   {
      ZmodF_poly_init(info->polys[1], depth, m + k, 1);
      ZmodF_poly_decrease_n(info->polys[1], m);
   }

   // todo: maybe can use less memory here when squaring:
   if (k)
      info->scratch = (mp_limb_t*) flint_stack_alloc((3*k) << depth);
   else
      info->scratch = NULL;
}


/*
  Looks up n in the tuning table and returns optimal negacyclic FFT depth.
*/
unsigned long _ZmodF_mul_best_fft_depth(unsigned long n, int squaring)
{
   unsigned long* table = squaring ? ZmodF_sqr_fft_table : ZmodF_mul_fft_table;

   unsigned long i;
   for (i = 0; table[i]; i++)
      if (n < table[i])
         return i + 3;

   // We've gone beyond the end of the table; need to choose a value
   // somewhat heuristically. We extrapolate from the last table entry,
   // assuming that the convolution length should be proportional to
   // sqrt(total bitsize).
   unsigned long depth = i + 3 +
                 (unsigned long) floor(log(1.0 * n / table[i-1]) / log(4.0));

   // need n*FLINT_BITS divisible by 2^depth
   while ((n*FLINT_BITS) & ((1 << depth) - 1))
      depth--;
   
   return depth;
}



void ZmodF_mul_info_init(ZmodF_mul_info_t info, unsigned long n, int squaring)
{
   if (!squaring)
   {
      if (n < ZmodF_mul_plain_threeway_threshold)
      {
         ZmodF_mul_info_init_plain(info, n, 0);
         return;
      }

      if (n % 3 == 0)
      {
         if (n < ZmodF_mul_threeway_fft_threshold)
         {
            ZmodF_mul_info_init_threeway(info, n, 0);
            return;
         }
      }
      else
      {
         if (n < ZmodF_mul_plain_fft_threshold)
         {
            ZmodF_mul_info_init_plain(info, n, 0);
            return;
         }
      }
   }
   else
   {
      if (n < ZmodF_sqr_plain_threeway_threshold)
      {
         ZmodF_mul_info_init_plain(info, n, 1);
         return;
      }

      if (n % 3 == 0)
      {
         if (n < ZmodF_sqr_threeway_fft_threshold)
         {
            ZmodF_mul_info_init_threeway(info, n, 1);
            return;
         }
      }
      else
      {
         if (n < ZmodF_sqr_plain_fft_threshold)
         {
            ZmodF_mul_info_init_plain(info, n, 1);
            return;
         }
      }
   }

   unsigned long depth = _ZmodF_mul_best_fft_depth(n, squaring);
   ZmodF_mul_info_init_fft(info, n, depth, 0, 0, squaring);
}


void ZmodF_mul_info_clear(ZmodF_mul_info_t info)
{
   if (info->scratch)
      flint_stack_release();

   if (info->algo == ZMODF_MUL_ALGO_FFT)
   {
      if (!info->squaring)
         ZmodF_poly_clear(info->polys[1]);
      ZmodF_poly_clear(info->polys[0]);
   }
}


/******************************************************************************

   FFT multiplication routines

******************************************************************************/

/*
Splits x into equally sized pieces.

That is, let R = (B^n)^(1/M), where M = transform length of "poly".
Let x = \sum_{i=0}^{M-1} c_i R^i, where each c_i is in [0, R).
This function computes the c_i and stores them as coefficients of "poly".

PRECONDITIONS:
   x must be normalised and of length n, and have zero overflow limb
   (i.e. is != -1 mod p)
   
   n*FLINT_BITS must be divisible by M

*/
void _ZmodF_mul_fft_split(ZmodF_poly_t poly, ZmodF_t x, unsigned long n)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << poly->depth) == 0);
   FLINT_ASSERT(x[n] == 0);

   unsigned long size = 1UL << poly->depth;
   
   // we'll split x into chunks each having "bits" bits
   unsigned long bits = (n * FLINT_BITS) >> poly->depth;
   // round it up to a whole number of limbs
   unsigned long limbs = ((bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1;
   
   // last_mask is applied to the last limb of each target coefficient to
   // zero out the bits that don't belong there
   unsigned long last_mask = (1UL << (bits & (FLINT_BITS-1))) - 1;
   if (!last_mask)
      last_mask = -1UL;

   // [start, end) are the bit-indices into x of the current chunk
   unsigned long start, end, i;
   for (i = 0, start = 0, end = bits; i < size; i++, start = end, end += bits)
   {
      // figure out which limbs contain the data for this chunk
      unsigned long start_limb = start >> FLINT_LG_BITS_PER_LIMB;
      unsigned long end_limb = ((end-1) >> FLINT_LG_BITS_PER_LIMB) + 1;
      
      // shift/copy the limbs containing the chunk into the target coefficient
      unsigned long start_bits = start & (FLINT_BITS-1);
      if (start_bits)
         mpn_rshift(poly->coeffs[i], x + start_limb, end_limb - start_limb,
                    start_bits);
      else
         F_mpn_copy(poly->coeffs[i], x + start_limb, end_limb - start_limb);

      // zero out the high bits that shouldn't contribute to this coefficient
      poly->coeffs[i][limbs-1] &= last_mask;

      // zero out remaining limbs
      F_mpn_clear(poly->coeffs[i] + limbs, poly->n + 1 - limbs);
   }
}



/*
Combines coefficients of poly into a ZmodF_t.

More precisely, let the coefficients of "poly" be c_0, ..., c_{M-1}, where
M = transform length of poly. The coefficients are assumed to be normalised,
and lie in the range [0, q), where q = (B^m + 1)*B^k. (Here "normalised" means
that the top m+1 limbs are normalised in the usual ZmodF_t sense; the bottom
k limbs are arbitrary.)

However the coefficients are interpreted as signed integers in the range
(-E, E) where E = B^(m+k)/2. (This is necessary since the output of a
negacyclic convolution can have negative coefficients.)

The function computes res := \sum_{i=0}^{M-1} c_i R^i   (mod p)
(not necessarily normalised), where R = (B^n)^(1/M).

NOTE: the input poly is destroyed.

*/
void _ZmodF_mul_fft_combine(ZmodF_t res, ZmodF_poly_t poly,
                            unsigned long m, unsigned long k, unsigned long n)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->n > 0);

   ZmodF_zero(res, n);

   unsigned long size = 1 << poly->depth;

   // "bits" is the number of bits apart that each coefficient must be stored
   unsigned long bits = (n * FLINT_BITS) >> poly->depth;
   // "start" is the bit-index into x where the current coeff will be stored
   unsigned long start;
   long i;

   // first loop: twiddle bits to deal with signs of coefficients
   // (do this in reverse order to avoid quadratic running time)
   for (i = size-1, start = (size-1)*bits; i >= 0; i--, start -= bits)
   {
      // start_limb, start_bit indicate where the current coefficient will
      // be stored; end_limb points beyond the last limb
      unsigned long start_limb = start >> FLINT_LG_BITS_PER_LIMB;
      unsigned long start_bit = start & (FLINT_BITS-1);
      unsigned long end_limb = start_limb + m + k + 1;

      if (poly->coeffs[i][m + k] ||
          ((mp_limb_signed_t) poly->coeffs[i][m + k - 1] < 0))
      {
         // coefficient is negative, so it was stored as (B^m + 1)*B^k plus
         // the real coefficient, so we need to correct for that
         unsigned long fiddle = 1UL << start_bit;
         
         // first correct for B^k
         unsigned long index = start_limb + k;
         if (index < n)
            mpn_sub_1(res + index, res + index, n + 1 - index, fiddle);
         else
         {
            // handle negacyclic wraparound
            index -= n;
            mpn_add_1(res + index, res + index, n + 1 - index, fiddle);
         }
         
         // now correct for B^(m+k)
         index = start_limb + m + k;
         if (index < n)
            mpn_sub_1(res + index, res + index, n + 1 - index, fiddle);
         else
         {
            // handle negacyclic wraparound
            index -= n;
            mpn_add_1(res + index, res + index, n + 1 - index, fiddle);
         }
      }
   }

   // main loop: add in bulk of coefficients
   for (i = 0, start = 0; i < size; i++, start += bits)
   {
      // start_limb, start_bit indicate where the current coefficient will
      // be stored; end_limb points beyond the last limb
      unsigned long start_limb = start >> FLINT_LG_BITS_PER_LIMB;
      unsigned long start_bit = start & (FLINT_BITS-1);
      unsigned long end_limb = start_limb + m + k + 1;

      // shift coefficient to the left to line it up with the spot where
      // it's going to get added in
      if (start_bit)
         mpn_lshift(poly->coeffs[i], poly->coeffs[i], m + k + 1, start_bit);
                    
      if (end_limb <= n)
      {
         // the coefficient fits nicely into the output, so add it in
         mpn_add(res + start_limb, res + start_limb, n + 1 - start_limb,
                 poly->coeffs[i], m + k + 1);
      }
      else
      {
         // the coefficient needs to wrap around negacyclically
         res[n] += mpn_add_n(res + start_limb, res + start_limb,
                             poly->coeffs[i], n - start_limb);
         mpn_sub(res, res, n + 1, poly->coeffs[i] + n - start_limb,
                 end_limb - n);
      }
   }
}



/*
in1 and in2 are two polynomials of length len, with coefficients mod B.
This function computes their negacyclic convolution (mod B), stores result
at out, also length len.

out must not alias in1 or in2.

Only naive convolution is implemented.
*/
void _ZmodF_mul_fft_convolve_modB(unsigned long* out, unsigned long* in1,
                                  unsigned long* in2, unsigned long len)
{
   unsigned long i, j;
   
   for (i = 0; i < len; i++)
      out[i] = in1[0] * in2[i];

   for (i = 1; i < len; i++)
   {
      for (j = 0; j < len - i; j++)
         out[i+j] += in1[i] * in2[j];
      for (; j < len; j++)
         out[i+j-len] -= in1[i] * in2[j];
   }
}



/*
   Multiplies (in1_hi*B + in1_lo) by (in2_hi*B + in2_lo) modulo B^2;
   puts result into (out_hi*B + out_lo).
*/
inline
void mul_modB2(unsigned long* out_hi, unsigned long* out_lo,
               unsigned long in1_hi, unsigned long in1_lo,
               unsigned long in2_hi, unsigned long in2_lo)
{
   umul_ppmm(*out_hi, *out_lo, in1_lo, in2_lo);
   *out_hi += (in1_hi * in2_lo) + (in2_hi * in1_lo);
}



/*
in1 and in2 are two polynomials of length len, with coefficients mod B^2.
This function computes their negacyclic convolution (mod B^2), stores result
at out, also length len. Each coefficient uses exactly 2 limbs.

out must not alias in1 or in2.

Only naive convolution is implemented.
*/
void _ZmodF_mul_fft_convolve_modB2(unsigned long* out, unsigned long* in1,
                                   unsigned long* in2, unsigned long len)
{
   unsigned long i, j;
   unsigned long hi, lo;

   for (i = 0; i < len; i++)
      mul_modB2(out + 2*i+1, out + 2*i, in1[1], in1[0], in2[2*i+1], in2[2*i]);

   for (i = 1; i < len; i++)
   {
      for (j = 0; j < len - i; j++)
      {
         mul_modB2(&hi, &lo, in1[2*i+1], in1[2*i], in2[2*j+1], in2[2*j]);
         add_ssaaaa(out[2*(i+j)+1], out[2*(i+j)],
                    out[2*(i+j)+1], out[2*(i+j)], hi, lo);
      }
      for (; j < len; j++)
      {
         mul_modB2(&hi, &lo, in1[2*i+1], in1[2*i], in2[2*j+1], in2[2*j]);
         sub_ddmmss(out[2*(i+j-len)+1], out[2*(i+j-len)],
                    out[2*(i+j-len)+1], out[2*(i+j-len)], hi, lo);
      }
   }
}


/*
   Reduces an array of ZmodF_t's mod B, stores result at out.
*/
void _ZmodF_mul_fft_reduce_modB(unsigned long* out, ZmodF_t* in,
                                unsigned long len)
{
   for (unsigned long i = 0; i < len; i++)
      out[i] = in[i][0];
}


/*
   Reduces an array of ZmodF_t's mod B^2, stores result at out (two limbs
   per coefficient).
*/
void _ZmodF_mul_fft_reduce_modB2(unsigned long* out, ZmodF_t* in,
                                 unsigned long len)
{
   for (unsigned long i = 0; i < len; i++)
   {
      out[2*i] = in[i][0];
      out[2*i+1] = in[i][1];
   }
}


/*
Computes res = a*b using FFT algorithm.
Output is not necessarily normalised.
If a == b, this function automatically uses a faster squaring algorithm.

Expects both inputs to be normalised and != -1 mod p.
*/
inline
void _ZmodF_mul_info_mul_fft(ZmodF_mul_info_t info, ZmodF_t res,
                             ZmodF_t a, ZmodF_t b)
{
   FLINT_ASSERT(info->algo == ZMODF_MUL_ALGO_FFT);
   FLINT_ASSERT(!a[info->n]);
   FLINT_ASSERT(!b[info->n]);

   unsigned long len = 1UL << info->polys[0]->depth;
   unsigned long n = info->n;
   unsigned long m = info->m;
   unsigned long k = info->k;

   if (a != b)
   {
      // distinct operands case
      FLINT_ASSERT(!info->squaring);

      // split inputs into 2^depth pieces
      _ZmodF_mul_fft_split(info->polys[0], a, n);
      _ZmodF_mul_fft_split(info->polys[1], b, n);

      // negacyclic convolution mod B^k
      if (k == 1)
      {
         _ZmodF_mul_fft_reduce_modB(info->scratch + len,
                                    info->polys[0]->coeffs, len);
         _ZmodF_mul_fft_reduce_modB(info->scratch + 2*len,
                                    info->polys[1]->coeffs, len);
         _ZmodF_mul_fft_convolve_modB(info->scratch, info->scratch + len,
                                      info->scratch + 2*len, len);
      }
      else if (k == 2)
      {
         _ZmodF_mul_fft_reduce_modB2(info->scratch + 2*len,
                                     info->polys[0]->coeffs, len);
         _ZmodF_mul_fft_reduce_modB2(info->scratch + 4*len,
                                     info->polys[1]->coeffs, len);
         _ZmodF_mul_fft_convolve_modB2(info->scratch, info->scratch + 2*len,
                                       info->scratch + 4*len, len);
      }

      // negacyclic convolution mod B^m + 1 using FFT
      ZmodF_poly_negacyclic_convolution(info->polys[0], info->polys[0],
                                        info->polys[1]);
   }
   else
   {
      // squaring case
      
      // split input into 2^depth pieces
      _ZmodF_mul_fft_split(info->polys[0], a, n);
      
      // negacyclic convolution mod B^k
      if (k == 1)
      {
         _ZmodF_mul_fft_reduce_modB(info->scratch + len,
                                    info->polys[0]->coeffs, len);
         _ZmodF_mul_fft_convolve_modB(info->scratch, info->scratch + len,
                                      info->scratch + len, len);
      }
      else if (k == 2)
      {
         _ZmodF_mul_fft_reduce_modB2(info->scratch + 2*len,
                                     info->polys[0]->coeffs, len);
         _ZmodF_mul_fft_convolve_modB2(info->scratch, info->scratch + 2*len,
                                       info->scratch + 2*len, len);
      }

      // negacyclic convolution mod B^m + 1 using FFT
      ZmodF_poly_negacyclic_convolution(info->polys[0], info->polys[0],
                                        info->polys[0]);
   }

   // use CRT to determine coefficients of convolution mod (B^m + 1)*B^k.
   // Basically the idea is: adjust bottom k limbs to make them agree with the
   // known mod B^k result, and compensate the top k limbs accordingly (i.e.
   // add a multiple of B^m + 1).
   if (k)
   {
      for (unsigned long i = 0; i < len; i++)
      {
         ZmodF_t coeff = info->polys[0]->coeffs[i];
         ZmodF_fast_reduce(coeff, m);

         // zero out the top limbs
         for (unsigned long j = 0; j < k; j++)
            coeff[m+1+j] = 0;

         // compute the amount we need to adjust by
         mp_limb_t adjust[2];
         mpn_sub_n(adjust, info->scratch + k*i, coeff, k);
         
         // perform adjustment
         mpn_add(coeff, coeff, m+k+1, adjust, k);
         mpn_add(coeff + m, coeff + m, k+1, adjust, k);

         // normalise result
         ZmodF_normalise(coeff + k, m);
      }
   }
   else
      ZmodF_poly_normalise(info->polys[0]);

   // substitute back to get integer mod B^n + 1
   _ZmodF_mul_fft_combine(res, info->polys[0], m, k, n);
}



/******************************************************************************

   Threeway multiplication routines

******************************************************************************/

/*
   Assume a is length 3m, and normalised, and != -1 mod p.
   Reduces a mod B^m + 1, stores result at res, in usual ZmodF_t format (m+1
   limbs), not necessarily normalised.
*/
inline
void _ZmodF_mul_threeway_reduce1(ZmodF_t res, ZmodF_t a, unsigned long m)
{
   FLINT_ASSERT(a[3*m] == 0);
   res[m] = mpn_add_n(res, a, a+2*m, m);
   res[m] -= mpn_sub_n(res, res, a+m, m);
}


/*
   Assume a is length 3m, and normalised, and != -1 mod p.
   Reduces a mod B^2m - B^m + 1, stores result at res, exactly 2m limbs.
   res must have room for 2m+1 limbs, even though the last limb is not used
   for the answer.

   Note: in some cases there are two possible answers, since B^2m - B^m + 1
   is less than B^2m. In these cases either answer may be produced.
*/
inline
void _ZmodF_mul_threeway_reduce2(mp_limb_t* res, ZmodF_t a, unsigned long m)
{
   res[2*m] = mpn_add_n(res+m, a+m, a+2*m, m);
   long borrow = mpn_sub_n(res, a, a+2*m, m);
   mpn_sub_1(res+m, res+m, m+1, borrow);
   if (res[2*m])
   {
      FLINT_ASSERT(res[2*m] == 1);
      // subtract B^2m - B^m + 1, then it's guaranteed to be normalised
      mpn_sub_1(res, res, 2*m, 1);
      mpn_add_1(res+m, res+m, m, 1);
   }
}


/*
   Computes res = a*b using threeway algorithm.
   Output is not necessarily normalised.
   If a == b, this function automatically uses a faster squaring algorithm.

   Expects both inputs to be normalised and != -1 mod p.
*/
inline
void _ZmodF_mul_info_mul_threeway(ZmodF_mul_info_t info,
                                  ZmodF_t res, ZmodF_t a, ZmodF_t b)
{
   FLINT_ASSERT(info->algo == ZMODF_MUL_ALGO_THREEWAY);
   FLINT_ASSERT(!a[info->n]);
   FLINT_ASSERT(!b[info->n]);

   unsigned long m = info->m;
   mp_limb_t* buf1 = info->scratch;
   mp_limb_t* buf2 = buf1 + m+1;
   mp_limb_t* buf3 = buf2 + 4*m;

   if (a != b)
   {
      // distinct operands case
      FLINT_ASSERT(!info->squaring);
      
      // reduce a and b mod B^m + 1
      _ZmodF_mul_threeway_reduce1(buf1, a, m);
      _ZmodF_mul_threeway_reduce1(buf2, b, m);

      // buf1 := a*b  mod B^m + 1
      ZmodF_mul(buf1, buf1, buf2, buf3, m);

      // reduce inputs mod B^2m - B^m + 1
      _ZmodF_mul_threeway_reduce2(buf2, a, m);
      _ZmodF_mul_threeway_reduce2(buf2 + 2*m, b, m);

      // multiply mod B^2m - B^m + 1
      mpn_mul_n(buf3, buf2, buf2 + 2*m, 2*m);
   }
   else
   {
      // squaring case

      // reduce a mod B^m + 1
      _ZmodF_mul_threeway_reduce1(buf1, a, m);

      // buf1 := a*a  mod B^m + 1
      ZmodF_mul(buf1, buf1, buf1, buf3, m);

      // reduce a mod B^2m - B^m + 1
      _ZmodF_mul_threeway_reduce2(buf2, a, m);

      // square mod B^2m - B^m + 1
      mpn_mul_n(buf3, buf2, buf2, 2*m);
   }

   // reduce mod B^3m + 1 (inplace)
   buf3[3*m] = -mpn_sub(buf3, buf3, 3*m, buf3 + 3*m, m);

   // Now buf3 (length 3m+1) is congruent to a*b mod B^2m - B^m + 1,
   // and buf1 (length m+1) is congruent to a*b mod B^m + 1.

   // Need to adjust buf3 to make it congruent to buf1 mod B^m + 1,
   // without modifying it mod B^2m - B^m + 1.
   // Strategy is:
   // Let X = (buf1 - buf3)/3 mod B^m + 1.
   // Add (B^2m - B^m + 1)*X to buf3.

   // buf2 := 3*X
   ZmodF_normalise(buf3, 3*m);
   if (buf3[3*m])
      // special case: buf3 = -1 mod B^3m + 1
      mpn_add_1(buf2, buf1, m+1, 1);
   else
   {
      _ZmodF_mul_threeway_reduce1(buf2, buf3, m);
      ZmodF_sub(buf2, buf1, buf2, m);
   }

   // buf2 := X
   ZmodF_divby3(buf2, buf2, m);
   ZmodF_normalise(buf2, m);
   
   // res := buf3 + (B^2m - B^m + 1)*X.
   if (buf2[m])
   {
      // special case: X = -1 mod B^m + 1
      ZmodF_set(res, buf3, 3*m);
      mpn_sub_1(res, res, 3*m+1, 1);
      mpn_add_1(res + m, res + m, 2*m+1, 1);
      mpn_sub_1(res + 2*m, res + 2*m, m+1, 1);
   }
   else
   {
      // usual case
      mp_limb_t carry1 = mpn_add_n(res, buf3, buf2, m);
      mp_limb_t carry2 = mpn_sub_n(res + m, buf3 + m, buf2, m);
      res[3*m] = buf3[3*m] + mpn_add_n(res + 2*m, buf3 + 2*m, buf2, m);
      mpn_add_1(res + m, res + m, 2*m+1, carry1);
      mpn_sub_1(res + 2*m, res + 2*m, m+1, carry2);
   }
}



/******************************************************************************

   Main ZmodF_mul_info multiplication routine

******************************************************************************/


void ZmodF_mul_info_mul(ZmodF_mul_info_t info,
                        ZmodF_t res, ZmodF_t a, ZmodF_t b)
{
   // try special cases a = -1 or b = -1 mod p
   if (a != b)
   {
      if (_ZmodF_mul_handle_minus1(res, a, b, info->n))
         return;
   }
   else
   {
      if (_ZmodF_sqr_handle_minus1(res, a, info->n))
         return;
   }

   if (info->algo == ZMODF_MUL_ALGO_PLAIN)
   {
      _ZmodF_mul(res, a, b, info->scratch, info->n);
      return;
   }

   else if (info->algo == ZMODF_MUL_ALGO_THREEWAY)
      _ZmodF_mul_info_mul_threeway(info, res, a, b);

   else
   {
      FLINT_ASSERT(info->algo == ZMODF_MUL_ALGO_FFT);
      _ZmodF_mul_info_mul_fft(info, res, a, b);
   }
}


// end of file ****************************************************************
