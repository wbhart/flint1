/******************************************************************************

 ZmodF_mul.c

 Copyright (C) 2007, David Harvey and William Hart
 
 Routines for multiplication of elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS.
 
******************************************************************************/

#include "ZmodF.h"
#include "ZmodF_poly.h"
#include "ZmodF_mul.h"
#include "Z_mpn.h"
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
      clear_limbs(res, n+1);
      return;
   }
   
   // do the product into scratch
   if (limbs1 >= limbs2)
      Z_mpn_mul(scratch, a, limbs1, b, limbs2);
   else
      Z_mpn_mul(scratch, b, limbs2, a, limbs1);
      
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
initialises info to use plain mpn_mul_n for multiplication
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
initialises info to use the "three-way split" algorithm for multiplication

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
initialises info to use a negacyclic FFT algorithm for multiplication,
with a transform depth of "depth".

n*FLINT_BITS must be divisible by 2^depth (so that the coefficients can be
broken up into 2^depth pieces)
*/
void ZmodF_mul_info_init_negacyclic(ZmodF_mul_info_t info, unsigned long n,
                                    unsigned long depth, int squaring)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << depth) == 0);

   info->algo = ZMODF_MUL_ALGO_NEGACYCLIC;
   info->n = n;
   info->squaring = squaring;

   // work out how many limbs the small coefficients need to have
   unsigned long input_bits = (n*FLINT_BITS) >> depth;
   unsigned long output_bits = 2*input_bits + 1 + depth;
   info->m = ((output_bits-1) >> FLINT_LG_BITS_PER_LIMB) + 1;

   // round up m so that 2^(depth+1) roots of unity are available,
   // i.e. need m*FLINT_BITS divisible by 2^depth
   if (depth > FLINT_LG_BITS_PER_LIMB)
   {
      unsigned long shift = depth - FLINT_LG_BITS_PER_LIMB;
      info->m = (((info->m - 1) >> shift) + 1) << shift;
   }
   
   ZmodF_poly_init(info->polys[0], depth, info->m, 1);
   if (!squaring)
      ZmodF_poly_init(info->polys[1], depth, info->m, 1);
}



/*
initialises info to use 2nd negacyclic FFT algorithm for multiplication,
with a transform depth of "depth".

n*FLINT_BITS must be divisible by 2^depth (so that the coefficients can be
broken up into 2^depth pieces)
*/
void ZmodF_mul_info_init_negacyclic2(ZmodF_mul_info_t info, unsigned long n, 
                                     unsigned long depth, int squaring)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << depth) == 0);

   info->algo = ZMODF_MUL_ALGO_NEGACYCLIC2;
   info->n = n;
   info->squaring = squaring;

   // work out how many limbs the small coefficients need to have
   // (remember one extra limb will be supplied by working mod B)
   unsigned long input_bits = (n*FLINT_BITS) >> depth;
   unsigned long output_bits = 2*input_bits + 1 + depth;
   info->m = ((output_bits-1) >> FLINT_LG_BITS_PER_LIMB);

   // round up m so that 2^(depth+1) roots of unity are available,
   // i.e. need m*FLINT_BITS divisible by 2^depth
   if (depth > FLINT_LG_BITS_PER_LIMB)
   {
      unsigned long shift = depth - FLINT_LG_BITS_PER_LIMB;
      info->m = (((info->m - 1) >> shift) + 1) << shift;
   }
   
   // we allocate m+2 limbs for each coefficient, but for the ZmodF_poly
   // routines we only use m+1 limbs (i.e. work mod B^m + 1). Then, later on
   // we'll need the full m+2 limbs to store the result mod B^(m+1) + B.
   ZmodF_poly_init(info->polys[0], depth, info->m+1, 1);
   ZmodF_poly_decrease_n(info->polys[0], info->m);
   if (!squaring)
   {
      ZmodF_poly_init(info->polys[1], depth, info->m+1, 1);
      ZmodF_poly_decrease_n(info->polys[1], info->m);
   }
   // todo: maybe can use less memory here when squaring:
   info->scratch = (mp_limb_t*) flint_stack_alloc(3UL << depth);
}


/*
Looks up n in the tuning table and returns optimal negacyclic FFT depth.
*/
unsigned long _ZmodF_mul_get_depth(unsigned long n)
{
   for (unsigned long depth = 0; ; depth++)
   {
      unsigned long value = ZMODF_MUL_NEGACYCLIC_THRESHOLD[depth];
      // bail out if either
      // (a) value == 0: means that n was too big for the table, or
      // (b) value > n: found the right value in the table

      // todo: for really really humungous input, (a) is the wrong strategy,
      // won't even be close to n log(n) anymore.... need a heuristic for this
      if (value == 0 || value > n)
         return depth + ZMODF_MUL_MIN_NEGACYCLIC_DEPTH;
   }
}


// same as above, but for squaring
unsigned long _ZmodF_sqr_get_depth(unsigned long n)
{
   for (unsigned long depth = 0; ; depth++)
   {
      unsigned long value = ZMODF_SQR_NEGACYCLIC_THRESHOLD[depth];
      if (value == 0 || value > n)
         return depth + ZMODF_SQR_MIN_NEGACYCLIC_DEPTH;
   }
}


void ZmodF_mul_info_init(ZmodF_mul_info_t info, unsigned long n, int squaring)
{
   if (!squaring)
   {
      if (n < ZMODF_MUL_PLAIN_THREEWAY_THRESHOLD)
      {
         // n is tiny, just use plain algorithm
         ZmodF_mul_info_init_plain(info, n, 0);
         return;
      }

#if ZMODF_MUL_ENABLE_THREEWAY
      if (n % 3 == 0)
      {
         if (n < ZMODF_MUL_THREEWAY_NEGACYCLIC_THRESHOLD)
         {
            // n is larger, use threeway algorithm
            ZmodF_mul_info_init_threeway(info, n, 0);
            return;
         }
      }
      else
#endif
      {
         if (n < ZMODF_MUL_PLAIN_NEGACYCLIC_THRESHOLD)
         {
            // n is small, just use plain algorithm
            ZmodF_mul_info_init_plain(info, n, 0);
            return;
         }
      }

#if ZMODF_MUL_ENABLE_NEGACYCLIC
      // neither plain nor threeway is appropriate;
      // need to select FFT of appropriate depth
      unsigned long depth = _ZmodF_mul_get_depth(n);
      
      // need to check that n supports an FFT of that depth, and if not,
      // take the largest possible depth that works
      while ((n * FLINT_BITS) & ((1 << depth) - 1))
         depth--;
         
//    switched off during development:
//      ZmodF_mul_info_init_negacyclic2(info, n, depth, 0);
      ZmodF_mul_info_init_negacyclic(info, n, depth, 0);

#else
      ZmodF_mul_info_init_plain(info, n, 0);
#endif
   }
   else
   {
      // squaring

      if (n < ZMODF_SQR_PLAIN_THREEWAY_THRESHOLD)
      {
         // n is tiny, just use plain algorithm
         ZmodF_mul_info_init_plain(info, n, 1);
         return;
      }

#if ZMODF_MUL_ENABLE_THREEWAY
      if (n % 3 == 0)
      {
         if (n < ZMODF_SQR_THREEWAY_NEGACYCLIC_THRESHOLD)
         {
            // n is larger, use threeway algorithm
            ZmodF_mul_info_init_threeway(info, n, 1);
            return;
         }
      }
      else
#endif
      {
         if (n < ZMODF_SQR_PLAIN_NEGACYCLIC_THRESHOLD)
         {
            // n is small, just use plain algorithm
            ZmodF_mul_info_init_plain(info, n, 1);
            return;
         }
      }
      
#if ZMODF_MUL_ENABLE_NEGACYCLIC
      // neither plain nor threeway is appropriate;
      // need to select FFT of appropriate depth
      unsigned long depth = _ZmodF_sqr_get_depth(n);
      
      // need to check that n supports an FFT of that depth, and if not,
      // take the largest possible depth that works
      while ((n * FLINT_BITS) & ((1 << depth) - 1))
         depth--;

//    switched off during development:
//      ZmodF_mul_info_init_negacyclic2(info, n, depth, 1);
      ZmodF_mul_info_init_negacyclic(info, n, depth, 1);

#else
      ZmodF_mul_info_init_plain(info, n, 1);
#endif
   }
}


void ZmodF_mul_info_clear(ZmodF_mul_info_t info)
{
   if (info->algo != ZMODF_MUL_ALGO_NEGACYCLIC)
   {
      flint_stack_release();
   }

   if (info->algo == ZMODF_MUL_ALGO_NEGACYCLIC ||
       info->algo == ZMODF_MUL_ALGO_NEGACYCLIC2)
   {
      if (!info->squaring)
         ZmodF_poly_clear(info->polys[1]);
      ZmodF_poly_clear(info->polys[0]);
   }
}


/******************************************************************************

   Negacyclic multiplication splitting/combining routines

******************************************************************************/

/*
Splits x into equally sized pieces.
Number of pieces = transform length of "poly".
Number of bits of x (i.e. n*FLINT_BITS) must be divisible by
the transform length.

Assumes x is normalised and of length n, and has zero overflow limb.
*/
void _ZmodF_mul_negacyclic_split(ZmodF_poly_t poly, ZmodF_t x, unsigned long n)
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
         copy_limbs(poly->coeffs[i], x + start_limb, end_limb - start_limb);

      // zero out the high bits that shouldn't contribute to this coefficient
      poly->coeffs[i][limbs-1] &= last_mask;

      // zero out remaining limbs
      clear_limbs(poly->coeffs[i] + limbs, poly->n + 1 - limbs);
   }
}


/*
Combines coefficients of poly into a ZmodF_t.

More precisely, let the coefficients of "poly" be c_0, ..., c_{M-1}, where
M = transform length of poly. The coefficients are assumed to be normalised,
and lie in the range [0, q), where q = B^m + 1, where m is the coefficient
length associated to poly. However this function interprets the coefficients
as representing integers in the range [-(q+1)/2, (q-3)/2)] since the output
of a negacyclic convolution can have *negative* coefficients.

Then it computes x := \sum_{i=0}^{M-1} c_i R^i   (mod p) (not necessarily
normalised), where R = (B^n)^(1/M).

*/
void _ZmodF_mul_negacyclic_combine(ZmodF_t x, ZmodF_poly_t poly,
                                   unsigned long n)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << poly->depth) == 0);

   ZmodF_zero(x, n);

   unsigned long i, carry;
   unsigned long size = 1 << poly->depth;
   unsigned long m = poly->n;

   // "bits" is the number of bits apart that each coefficient must be stored
   unsigned long bits = (n * FLINT_BITS) >> poly->depth;
   // "start" is the bit-index into x where the current coeff will be stored
   unsigned long start;
   
   for (i = 0, start = 0; i < size; i++, start += bits)
   {
      // start_limb, start_bit indicate where the current coefficient will
      // be stored; end_limb points beyond the last limb
      unsigned long start_limb = start >> FLINT_LG_BITS_PER_LIMB;
      unsigned long start_bit = start & (FLINT_BITS-1);
      unsigned long end_limb = start_limb + m + 1;

      if (poly->coeffs[i][m])
      {
         // special case: coefficient is -1 mod p
         mpn_sub_1(x + start_limb, x + start_limb, n + 1 - start_limb,
                   1UL << start_bit);
      }
      else
      {
         long negative = ((mp_limb_signed_t) poly->coeffs[i][m - 1] < 0);
      
         // shift coefficient to the left to line it up with the spot where
         // it's going to get added in
         if (start_bit)
            mpn_lshift(poly->coeffs[i], poly->coeffs[i], m + 1, start_bit);
                       
         if (end_limb <= n)
         {
            // the coefficient fits nicely into the output, so add it in
            mpn_add(x + start_limb, x + start_limb, n + 1 - start_limb,
                    poly->coeffs[i], m + 1);

            if (negative)
            {
               // If the coefficient is negative, then it's stored as
               // B^m + 1 plus the real coefficient, so we need to correct
               // for that
               
               // todo: are we possibly inducing quadratic runtime here!!!!??
               // i.e. with long borrow propagations?
               mpn_sub_1(x + start_limb, x + start_limb, n + 1 - start_limb,
                         1UL << start_bit);
               mpn_sub_1(x + end_limb - 1, x + end_limb - 1,
                         n + 2 - end_limb, 1UL << start_bit);
            }
         }
         else
         {
            // the coefficient needs to wrap around negacyclically
            x[n] += mpn_add_n(x + start_limb, x + start_limb, poly->coeffs[i],
                              n - start_limb);
            mpn_sub(x, x, n + 1, poly->coeffs[i] + n - start_limb,
                    end_limb - n);

            if (negative)
            {
               // handle negative coefficient as above
               mpn_sub_1(x + start_limb, x + start_limb, n + 1 - start_limb,
                         1UL << start_bit);
               end_limb -= n;
               mpn_add_1(x + end_limb - 1, x + end_limb - 1,
                         n + 2 - end_limb, 1UL << start_bit);
            }
         }
      }
   }
}


/******************************************************************************

   Negacyclic multiplication (2nd version) splitting/combining routines

******************************************************************************/


/*
in1 and in2 are two polynomials of length len, with coefficients mod B.
This function computes their negacyclic convolution (mod B), stores result
at out, also length len.

out must not alias in1 or in2.

Currently only naive convolution is implemented.
*/
void _ZmodF_mul_negacyclic2_convolve_modB(
            mp_limb_t* out, mp_limb_t* in1, mp_limb_t* in2, unsigned long len)
{
   unsigned long i, j;

   for (i = 0; i < len; i++)
      out[i] = in1[i] * in2[0];

   for (i = 1; i < len; i++)
   {
      for (j = 0; j < len - i; j++)
         out[i+j] += in1[i] * in2[j];
      for (; j < len; j++)
         out[i+j-len] -= in1[i] * in2[j];
   }
}



/*
Combines coefficients of poly into a ZmodF_t.

More precisely, let the coefficients of "poly" be c_0, ..., c_{M-1}, where
M = transform length of poly. The coefficients are assumed to be normalised,
and lie in the range [0, q), where q = B^(m+1) + B, where m is the coefficient
length associated to poly. (Note that we've abused the ZmodF_poly slightly,
since we are storing m+2 limbs per coefficient instead of m+1 limbs.)

However this function interprets the coefficients as representing integers in
the range (-B^(m+1)/2, B^(m+1)/2) since the output of a negacyclic convolution
can have *negative* coefficients.

Then it computes x := \sum_{i=0}^{M-1} c_i R^i   (mod p) (not necessarily
normalised), where R = (B^n)^(1/M).

*/
void _ZmodF_mul_negacyclic2_combine(ZmodF_t x, ZmodF_poly_t poly,
                                    unsigned long n)
{
   FLINT_ASSERT((n * FLINT_BITS) % (1 << poly->depth) == 0);

   ZmodF_zero(x, n);

   unsigned long i, carry;
   unsigned long size = 1 << poly->depth;
   unsigned long m = poly->n;

   // "bits" is the number of bits apart that each coefficient must be stored
   unsigned long bits = (n * FLINT_BITS) >> poly->depth;
   // "start" is the bit-index into x where the current coeff will be stored
   unsigned long start;
   
   for (i = 0, start = 0; i < size; i++, start += bits)
   {
      // start_limb, start_bit indicate where the current coefficient will
      // be stored; end_limb points beyond the last limb
      unsigned long start_limb = start >> FLINT_LG_BITS_PER_LIMB;
      unsigned long start_bit = start & (FLINT_BITS-1);
      unsigned long end_limb = start_limb + m + 2;

      if (poly->coeffs[i][m + 1])
      {
         abort();
         /*
         // special case: coefficient has overflow limb set....
         mpn_sub_1(x + start_limb, x + start_limb, n + 1 - start_limb,
                   1UL << start_bit);
         */
      }
      else
      {
         long negative = ((mp_limb_signed_t) poly->coeffs[i][m] < 0);
      
         // shift coefficient to the left to line it up with the spot where
         // it's going to get added in
         if (start_bit)
            mpn_lshift(poly->coeffs[i], poly->coeffs[i], m + 2, start_bit);
                       
         if (end_limb <= n)
         {
            // the coefficient fits nicely into the output, so add it in
            mpn_add(x + start_limb, x + start_limb, n + 1 - start_limb,
                    poly->coeffs[i], m + 2);

            if (negative)
            {
               // If the coefficient is negative, then it's stored as
               // B^(m+1) + B plus the real coefficient, so we need to correct
               // for that
               mpn_sub_1(x + start_limb + 1, x + start_limb + 1,
                         n - start_limb, 1UL << start_bit);
               mpn_sub_1(x + end_limb - 1, x + end_limb - 1,
                         n + 2 - end_limb, 1UL << start_bit);
            }
         }
         else
         {
            // the coefficient needs to wrap around negacyclically
            x[n] += mpn_add_n(x + start_limb, x + start_limb, poly->coeffs[i],
                              n - start_limb);
            mpn_sub(x, x, n + 1, poly->coeffs[i] + n - start_limb,
                    end_limb - n);

            if (negative)
            {
               // handle negative coefficient as above
               mpn_sub_1(x + start_limb + 1, x + start_limb + 1,
                         n + 1 - start_limb, 1UL << start_bit);
               end_limb -= n;
               mpn_add_1(x + end_limb - 1, x + end_limb - 1,
                         n + 2 - end_limb, 1UL << start_bit);
            }
         }
      }
   }
}


/******************************************************************************

   Threeway multiplication splitting/combining routines

******************************************************************************/

/*
Assume a is length 3m, and normalised, and != -1 mod p.
Reduces a mod B^m + 1, stores result at res, in usual ZmodF_t format (m+1
limbs), in normalised form.
 */
inline
void _ZmodF_mul_threeway_reduce1(ZmodF_t res, ZmodF_t a, unsigned long m)
{
   res[m] = mpn_add_n(res, a, a+2*m, m);
   res[m] -= mpn_sub_n(res, res, a+m, m);
   ZmodF_normalise(res, m);
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
given a of length m+1 (normalised, ZmodF format) and b of length 2*m,
computes res of length 3*m (ZmodF format, not normalised)
such that res = a mod B^m + 1, and res = b mod B^2m - B^m + 1
 */
inline
void _ZmodF_mul_threeway_crt(mp_limb_t* res, ZmodF_t a, mp_limb_t* b,
                             unsigned long m)
{
   unsigned long carry0, carry1;

   // --------- idempotent for b is (-B^2m + B^m + 2)/3

   // let b = b0 + B^m*b1,
   // so (-B^2m + B^m + 2)*(b0 + B^m*b1)
   // = (2*b0 + b1) + B^m*(2*b1 + b0) + B^2m*(-b0 + b1)

   // put b0 + b1 into first m limbs
   carry0 = mpn_add_n(res, b, b+m, m);
   // put 2*b1 + b0 into second m limbs
   carry1 = mpn_add_n(res+m, res, b+m, m) + carry0;
   // put 2*b0 + b1 into first m limbs
   carry0 += mpn_add_n(res, res, b, m);
   // put b1 - b0 into third m limbs
   res[3*m] = -mpn_sub_n(res+2*m, b+m, b, m);

   // propagate carries
   mpn_add_1(res+m, res+m, 2*m+1, carry0);
   mpn_add_1(res+2*m, res+2*m, m+1, carry1);

   // --------- idempotent for a is (B^2m - B^m + 1)/3

   if (a[m])
   {
      // special case if a = -1 mod B^m + 1
      mpn_sub_1(res, res, 3*m+1, 1);
      mpn_add_1(res+m, res+m, 2*m+1, 1);
      mpn_sub_1(res+2*m, res+2*m, m+1, 1);
   }
   else
   {
      // usual case
      mpn_add(res, res, 3*m+1, a, m);
      mpn_sub(res+m, res+m, 2*m+1, a, m);
      mpn_add(res+2*m, res+2*m, m+1, a, m);
   }

   // --------- correct for idempotents being off by factor of 3

   // make overflow limb nonnegative
   ZmodF_fast_reduce(res, 3*m);

   // compute a "total" which is congruent to res mod 3
   unsigned long total = 0;
   for (unsigned long i = 0; i < 3*m+1; i++)
   {
      total += (res[i] & ((1UL << (FLINT_BITS/2)) - 1));
      total += (res[i] >> (FLINT_BITS/2));
   }

   // add "total" times B^3n + 1 (the latter is 2 mod 3)
   mpn_add_1(res, res, 3*m+1, total);
   res[3*m] += total;
   
   unsigned long rem = mpn_divexact_by3(res, res, 3*m+1);
   FLINT_ASSERT(!rem);

   ZmodF_fast_reduce(res, 3*m);
}



/******************************************************************************

   Main ZmodF_mul_info multiplication routines

******************************************************************************/


void ZmodF_mul_info_mul(ZmodF_mul_info_t info,
                        ZmodF_t res, ZmodF_t a, ZmodF_t b)
{
   FLINT_ASSERT(!info->squaring);

   // try special cases a = -1 or b = -1 mod p
   if (_ZmodF_mul_handle_minus1(res, a, b, info->n))
      return;

   if (info->algo == ZMODF_MUL_ALGO_PLAIN)
   {
      // plain mpn_mul_n multiplication
      _ZmodF_mul(res, a, b, info->scratch, info->n);
      return;
   }
   else if (info->algo == ZMODF_MUL_ALGO_THREEWAY)
   {
      // threeway split multiplication
      unsigned long m = info->m;
      mp_limb_t* buf1 = info->scratch;
      mp_limb_t* buf2 = buf1 + m+1;
      mp_limb_t* buf3 = buf2 + 4*m;

      // reduce inputs mod B^m + 1
      _ZmodF_mul_threeway_reduce1(buf1, a, m);
      _ZmodF_mul_threeway_reduce1(buf2, b, m);
      // multiply mod B^m + 1
      ZmodF_mul(buf1, buf1, buf2, buf3, m);
      ZmodF_normalise(buf1, m);

      // reduce inputs mod B^2m - B^m + 1
      _ZmodF_mul_threeway_reduce2(buf2, a, m);
      _ZmodF_mul_threeway_reduce2(buf2 + 2*m, b, m);

      // multiply
      mpn_mul_n(buf3, buf2, buf2 + 2*m, 2*m);
      // reduce mod B^3m + 1 (inplace)
      buf3[3*m] = -mpn_sub(buf3, buf3, 3*m, buf3 + 3*m, m);

      ZmodF_normalise(buf3, 3*m);
      // reduce mod B^2m - B^m + 1
      _ZmodF_mul_threeway_reduce2(buf2, buf3, m);
      
      // combine results to get answer mod B^3m + 1
      _ZmodF_mul_threeway_crt(res, buf1, buf2, m);
   }
   else if (info->algo == ZMODF_MUL_ALGO_NEGACYCLIC)
   {
      // negacyclic FFT multiplication
      _ZmodF_mul_negacyclic_split(info->polys[0], a, info->n);
      _ZmodF_mul_negacyclic_split(info->polys[1], b, info->n);

      ZmodF_poly_negacyclic_convolution(info->polys[0], info->polys[0],
                                        info->polys[1]);
      ZmodF_poly_normalise(info->polys[0]);

      _ZmodF_mul_negacyclic_combine(res, info->polys[0], info->n);
   }
   else
   {
      // negacyclic FFT multiplication, version 2
      FLINT_ASSERT(info->algo == ZMODF_MUL_ALGO_NEGACYCLIC2);

      // split inputs into 2^depth pieces
      _ZmodF_mul_negacyclic_split(info->polys[0], a, info->n);
      _ZmodF_mul_negacyclic_split(info->polys[1], b, info->n);
      
      // reduce all coefficient mod B, store results at scratch + len and
      // scratch + 2*len
      unsigned long len = 1UL << info->polys[0]->depth;
      for (unsigned long i = 0; i < len; i++)
         info->scratch[len + i] = info->polys[0]->coeffs[i][0];
      for (unsigned long i = 0; i < len; i++)
         info->scratch[2*len + i] = info->polys[1]->coeffs[i][0];

      // do convolution mod B
      _ZmodF_mul_negacyclic2_convolve_modB(
             info->scratch, info->scratch + len, info->scratch + 2*len, len);
      
      // do the convolution mod B^m + 1
      ZmodF_poly_negacyclic_convolution(info->polys[0], info->polys[0],
                                        info->polys[1]);

      // use CRT to determine coefficients of convolution mod B^(m+1) + B
      // (store them as m+2 limbs each, in 2's complement, similar to ZmodF_t,
      // results are normalised)
      // todo: factor out this code into a function
      for (unsigned long i = 0; i < len; i++)
      {
         mp_limb_t* coeff = info->polys[0]->coeffs[i];
         ZmodF_normalise(coeff, info->m);
         
         if (coeff[info->m])
         {
            // special case for -1 mod B^m + 1
            // todo: write this
            abort();
         }
         else
         {
            // usual case, when not -1 mod B^m + 1
            sub_ddmmss(coeff[info->m+1], coeff[info->m], 0, info->scratch[i],
                       0, coeff[0]);
            coeff[0] = info->scratch[i];
         }

         // normalise
         ZmodF_normalise(coeff+1, info->m);
      }

      // combine results
      _ZmodF_mul_negacyclic2_combine(res, info->polys[0], info->n);
   }
}


void ZmodF_mul_info_sqr(ZmodF_mul_info_t info, ZmodF_t res, ZmodF_t a)
{
   // try special case a = -1 mod p
   if (_ZmodF_sqr_handle_minus1(res, a, info->n))
      return;

   if (info->algo == ZMODF_MUL_ALGO_PLAIN)
   {
      // plain mpn_mul_n multiplication
      _ZmodF_mul(res, a, a, info->scratch, info->n);
      return;
   }
   else if (info->algo == ZMODF_MUL_ALGO_THREEWAY)
   {
      // threeway split multiplication
      unsigned long m = info->m;
      mp_limb_t* buf1 = info->scratch;
      mp_limb_t* buf2 = buf1 + m+1;
      mp_limb_t* buf3 = buf2 + 4*m;

      // reduce input mod B^m + 1
      _ZmodF_mul_threeway_reduce1(buf1, a, m);
      // multiply mod B^m + 1
      ZmodF_sqr(buf1, buf1, buf3, m);
      ZmodF_normalise(buf1, m);

      // reduce input mod B^2m - B^m + 1
      _ZmodF_mul_threeway_reduce2(buf2, a, m);

      // multiply
      mpn_mul_n(buf3, buf2, buf2, 2*m);
      // reduce mod B^3m + 1 (inplace)
      buf3[3*m] = -mpn_sub(buf3, buf3, 3*m, buf3 + 3*m, m);

      ZmodF_normalise(buf3, 3*m);
      // reduce mod B^2m - B^m + 1
      _ZmodF_mul_threeway_reduce2(buf2, buf3, m);
      
      // combine results to get answer mod B^3m + 1
      _ZmodF_mul_threeway_crt(res, buf1, buf2, m);
   }
   else if (info->algo == ZMODF_MUL_ALGO_NEGACYCLIC)
   {
      // negacyclic FFT multiplication
      _ZmodF_mul_negacyclic_split(info->polys[0], a, info->n);

      ZmodF_poly_negacyclic_convolution(info->polys[0], info->polys[0],
                                       info->polys[0]);
      ZmodF_poly_normalise(info->polys[0]);

      _ZmodF_mul_negacyclic_combine(res, info->polys[0], info->n);
   }
   else
   {
      // negacyclic FFT multiplication, version 2
      FLINT_ASSERT(info->algo == ZMODF_MUL_ALGO_NEGACYCLIC2);

      // split inputs into 2^depth pieces
      _ZmodF_mul_negacyclic_split(info->polys[0], a, info->n);
      
      // reduce all coefficient mod B, store results at scratch + len
      unsigned long len = 1UL << info->polys[0]->depth;
      for (unsigned long i = 0; i < len; i++)
         info->scratch[len + i] = info->polys[0]->coeffs[i][0];

      // do convolution mod B
      _ZmodF_mul_negacyclic2_convolve_modB(
             info->scratch, info->scratch + len, info->scratch + len, len);
      
      // now do the convolution mod B^m + 1
      ZmodF_poly_negacyclic_convolution(info->polys[0], info->polys[0],
                                        info->polys[0]);

      // use CRT to determine coefficients of convolution mod B^(m+1) + B
      // (store them as m+2 limbs each, in 2's complement, similar to ZmodF_t,
      // results are normalised)
      for (unsigned long i = 0; i < len; i++)
      {
         mp_limb_t* coeff = info->polys[0]->coeffs[i];
         ZmodF_normalise(coeff, info->m);
         
         if (coeff[info->m])
         {
            // special case for -1 mod B^m + 1
            // todo: write this
            abort();
         }
         else
         {
            // usual case, when not -1 mod B^m + 1
            sub_ddmmss(coeff[info->m+1], coeff[info->m], 0, info->scratch[i],
                       0, coeff[0]);
            coeff[0] = info->scratch[i];
         }

         // normalise
         ZmodF_normalise(coeff+1, info->m);
      }

      // combine results
      _ZmodF_mul_negacyclic2_combine(res, info->polys[0], info->n);
   }
}


// end of file ****************************************************************
