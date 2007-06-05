/******************************************************************************

 ZmodF_mul.c

 Copyright (C) 2007, David Harvey
 
 Routines for multiplication of elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
******************************************************************************/

#include "ZmodF.h"
#include "ZmodFpoly.h"
#include "ZmodF_mul.h"


#define ENABLE_NEGACYCLIC_MULTS 1


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
   assumes a and b are normalised, and that their overflow limbs are 0
*/
static inline
void _ZmodF_mul(ZmodF_t res, ZmodF_t a, ZmodF_t b, mp_limb_t* scratch,
                unsigned long n)
{
   FLINT_ASSERT(a[n] == 0);
   
   // do the product
   mpn_mul_n(scratch, a, b, n);
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


// when the coefficient size is between negacyclic_threshold[i] and
// negacyclic_threshold[i+1] limbs, we use a negacyclic FFT of depth
// i + min_negacyclic_depth. If it's less than the first entry, we just use
// plain old mpn_mul_n

// todo: the first four of these values are tuned for sage.math, the rest
// I pulled from my rear end... we need automatic tuning code for this
// (it's very difficult because the performance is not very smooth yet)
unsigned long negacyclic_threshold_table[] = 
   {320, 500, 1000, 2200, 5000, 10000, 20000, 40000};

unsigned long negacyclic_threshold_table_size =
   sizeof(negacyclic_threshold_table) / sizeof(negacyclic_threshold_table[0]);

unsigned long min_negacyclic_depth = 4;


unsigned long ZmodF_mul_precomp_get_feasible_n(unsigned long *depth,
                                               unsigned long n)
{
   // n is small enough to just use mpn_mul_n directly
   if (n < negacyclic_threshold_table[0])
      return n;
   
   // n is big enough, search the table for d = best convolution depth
   unsigned long d = 0;
   for (unsigned long i = 1; i < negacyclic_threshold_table_size; i++)
   {
      if (n < negacyclic_threshold_table[i])
      {
         d = i + min_negacyclic_depth - 1;
         break;
      }
   }
   
   if (d == 0)
   {
      // n was too big for the table, we need to estimate a value for d

      // todo: for now just use the maximum table size, since we aren't
      // planning on multiplying anything that humungous anyway...?
      d = negacyclic_threshold_table_size + min_negacyclic_depth;
   }

   // Round up n to a multiple of 2^d/FLINT_BITS_PER_LIMB, so that n can
   // be broken up into 2^d chunks.
   if (d > FLINT_LG_BITS_PER_LIMB)
   {
      unsigned long shift = d - FLINT_LG_BITS_PER_LIMB;
      n = (((n-1) >> shift) + 1) << shift;
   }
      
   if (depth != NULL)
      *depth = d;
   
   return n;
}


void ZmodF_mul_precomp_init(ZmodF_mul_precomp_t info, unsigned long n,
                            int squaring)
{
   unsigned long depth;
   info->n = n;
   
#if ENABLE_NEGACYCLIC_MULTS
   if ((n < negacyclic_threshold_table[0]) ||
       (n != ZmodF_mul_precomp_get_feasible_n(&depth, n)))
#else
   if (1)
#endif
   {
      // use mpn_mul_n
      info->use_fft = 0;

      // todo: this should use stack-based memory manager; and the
      // documentation for this function NEEDS TO MENTION THIS
      info->scratch = (mp_limb_t*) malloc(2*n * sizeof(mp_limb_t));
   }
   else
   {
      // use negacyclic FFTs
      info->use_fft = 1;

      // work out how many limbs the small coefficients need to have
      unsigned long input_bits = (n*FLINT_BITS_PER_LIMB) >> depth;
      unsigned long output_bits = 2*input_bits + 1 + depth;
      unsigned long next_n = ((output_bits-1) >> FLINT_LG_BITS_PER_LIMB) + 1;

      // round up next_n so that enough roots of unity are available,
      // i.e. need FLINT_BITS_PER_LIMB*next_n divisible by 2^(depth-1)
      if (depth-1 > FLINT_LG_BITS_PER_LIMB)
      {
         unsigned long shift = depth - 1 - FLINT_LG_BITS_PER_LIMB;
         next_n = (((next_n - 1) >> shift) + 1) << shift;
      }

      ZmodFpoly_init(info->polys[0], depth, next_n, 1);
      ZmodFpoly_init(info->polys[1], depth, next_n, 1);
   }
}


void ZmodF_mul_precomp_clear(ZmodF_mul_precomp_t info)
{
   if (info->use_fft)
   {
      ZmodFpoly_clear(info->polys[1]);
      ZmodFpoly_clear(info->polys[0]);
   }
   else
   {
      free(info->scratch);
   }
}


/*
Splits coefficient x into a ZmodFpoly_t.
Assumes x is normalised and of length n.
*/
void _ZmodF_mul_split(ZmodFpoly_t poly, ZmodF_t x, unsigned long n)
{
   FLINT_ASSERT((n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);

   unsigned long size = 1UL << poly->depth;
   
   // we'll split x into chunks each having "bits" bits
   unsigned long bits = (n * FLINT_BITS_PER_LIMB) >> poly->depth;
   // round it up to a whole number of limbs
   unsigned long limbs = ((bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1;
   
   // last_mask is applied to the last limb of each target coefficient to
   // zero out the bits that don't belong there
   unsigned long last_mask = (1UL << (bits & (FLINT_BITS_PER_LIMB-1))) - 1;
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
      unsigned long start_bits = start & (FLINT_BITS_PER_LIMB-1);
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
Combines coefficients of poly into a ZmodF_t, doing things appropriately mod p.
Coefficients of poly must be normalised.
*/
void _ZmodF_mul_combine(ZmodF_t x, ZmodFpoly_t poly, unsigned long n)
{
   ZmodF_zero(x, n);

   unsigned long i, carry;
   unsigned long size = 1 << poly->depth;

   // "bits" is the number of bits apart that each coefficient must be stored
   unsigned long bits = (n * FLINT_BITS_PER_LIMB) >> poly->depth;
   // "start" is the bit-index into x where the current coeff will be stored
   unsigned long start;
   
   for (i = 0, start = 0; i < size; i++, start += bits)
   {
      // start_limb, start_bit indicate where the current coefficient will
      // be stored; end_limb points beyond the last limb
      unsigned long start_limb = start >> FLINT_LG_BITS_PER_LIMB;
      unsigned long start_bit = start & (FLINT_BITS_PER_LIMB-1);
      unsigned long end_limb = start_limb + poly->n + 1;

      if (poly->coeffs[i][poly->n])
      {
         // special case: coefficient is -1 mod p
         mpn_sub_1(x + start_limb, x + start_limb, n + 1 - start_limb,
                   1UL << start_bit);
      }
      else
      {
         long negative = ((mp_limb_signed_t) poly->coeffs[i][poly->n - 1] < 0);
      
         // shift coefficient to the left to line it up with the spot where
         // it's going to get added in
         if (start_bit)
            mpn_lshift(poly->coeffs[i], poly->coeffs[i], poly->n + 1,
                       start_bit);
                       
         if (end_limb <= n)
         {
            // the coefficient fits nicely into the output, so add it in
            mpn_add(x + start_limb, x + start_limb, n + 1 - start_limb,
                    poly->coeffs[i], poly->n + 1);

            if (negative)
            {
               // If the coefficient is negative, then it's stored as
               // 2^(FLINT_BITS_PER_LIMB*poly->n) + 1 plus the real
               // coefficient, so we need to correct for that
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


void ZmodF_mul_precomp(ZmodF_mul_precomp_t info,
                       ZmodF_t res, ZmodF_t a, ZmodF_t b)
{
   // try special cases a = -1 or b = -1 mod p
   if (_ZmodF_mul_handle_minus1(res, a, b, info->n))
      return;
   
   if (!info->use_fft)
   {
      // plain mpn_mul_n multiplication
      _ZmodF_mul(res, a, b, info->scratch, info->n);
      return;
   }
   
   // negacyclic FFT multiplication
   _ZmodF_mul_split(info->polys[0], a, info->n);
   _ZmodF_mul_split(info->polys[1], b, info->n);

   ZmodFpoly_negacyclic_convolution(info->polys[0], info->polys[0],
                                    info->polys[1]);
   ZmodFpoly_normalise(info->polys[0]);

   _ZmodF_mul_combine(res, info->polys[0], info->n);
}


void ZmodF_sqr_precomp(ZmodF_mul_precomp_t info, ZmodF_t res, ZmodF_t a)
{
   // try special case a = -1 mod p
   if (_ZmodF_sqr_handle_minus1(res, a, info->n))
      return;

   if (!info->use_fft)
   {
      // plain mpn_mul_n multiplication
      _ZmodF_mul(res, a, a, info->scratch, info->n);
      return;
   }
   
   // negacyclic FFT multiplication
   _ZmodF_mul_split(info->polys[0], a, info->n);

   ZmodFpoly_negacyclic_convolution(info->polys[0], info->polys[0],
                                    info->polys[0]);
   ZmodFpoly_normalise(info->polys[0]);

   _ZmodF_mul_combine(res, info->polys[0], info->n);
}


// end of file ****************************************************************
