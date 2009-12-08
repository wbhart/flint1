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
/******************************************************************************

 mpn_extras.c
 
 Extra functions for manipulating mpn's and limbs.

 Copyright (C) 2006, William Hart

 mp_limb_t mpn_divmod_1_preinv was adapted from GMP, (C) Free Software Foundation

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "long_extras.h"
#include "memory-manager.h"
#include "mpn_extras.h"
#include "ZmodF_poly.h"
#include "ZmodF_mul.h"
#include "F_mpn_mul-tuning.h"

#define DEBUG2 1

/*=======================================================================================

    Performs division by a limb d and places the quotient in qp and returns the 
    remainder. Requires a single limb approximation to 1/d as input. If the most
    significant bit of d is not 1 it expects d to be shifted left (by norm bits)
    until the most significant bit is 1 before the inverse is computed. However 
    the original d should be supplied to the function, not the shifted d. 
    
    This code has been adapted from code found in the GMP package version 4.2.1
    (divrem_1.c) (C) Free Software Foundation
*/

mp_limb_t F_mpn_divrem_ui_precomp(mp_limb_t * qp, mp_limb_t * up, 
                                  unsigned long un, mp_limb_t d, mp_limb_t dinv)
{
  mp_size_t  n;
  mp_size_t  i;
  mp_limb_t  n1, n0;
  mp_limb_t  r = 0;
 
  unsigned long norm;

  n = un;
  if (n == 0)
    return 0;
  
  count_lead_zeros(norm, d);

  qp += (n - 1);   /* Make qp point at most significant quotient limb */

  if ((d & (1L<<(FLINT_BITS-1))) != 0)
  {
     if (un != 0)
     {
        /* High quotient limb is 0 or 1, skip a divide step. */
	    mp_limb_t q;
	    r = up[un - 1];
	    q = (r >= d);
	    *qp-- = q;
	    r -= (d & -q);
	    n--;
	    un--;
	 }

     /* Multiply-by-inverse, divisor already normalized. */
     for (i = un - 1; i >= 0; i--)
     {
        n0 = up[i];
        udiv_qrnnd_preinv (*qp, r, r, n0, d, dinv);
        qp--;
     }
     return r;
  } else
  {
     /* Most significant bit of divisor == 0.  */
     
     /* Skip a division if high < divisor (high quotient 0).  Testing here
	 before normalizing will still skip as often as possible.  */
     if (un != 0)
	 {
	    n1 = up[un - 1];
	    if (n1 < d)
        {
           r = n1;
	       *qp-- = 0;
	       n--;
	       if (n == 0) return r;
	       un--;
        }
	 }  

     d <<= norm;
     r <<= norm;

     if (un != 0)
     {
        n1 = up[un - 1];
        r |= (n1 >> (FLINT_BITS - norm));
        for (i = un - 2; i >= 0; i--)
		{
		  n0 = up[i];
		  udiv_qrnnd_preinv (*qp, r, r, 
				     ((n1 << norm) | (n0 >> (FLINT_BITS - norm))), d, dinv);
		  qp--;
		  n1 = n0;
		}
        udiv_qrnnd_preinv (*qp, r, r, n1 << norm, d, dinv);
        qp--;
     }
     
     return r >> norm;
  }
}

mp_limb_t F_mpn_addmul(mp_limb_t * rp, mp_limb_t * s1p, unsigned long s1n, 
                                      mp_limb_t * s2p, unsigned long s2n)
{
   if (s2n == 0) return 0;
   
   mp_limb_t carry;
   
   carry = mpn_addmul_1(rp, s1p, s1n, s2p[0]);
   unsigned long i;
   for (i = 1; i < s2n; i++)
   {
      carry = mpn_add_1(rp+i+s1n-1, rp+i+s1n-1, 1, carry); 
      if (s2p[i]) carry += mpn_addmul_1(rp+i, s1p, s1n, s2p[i]);
   }
   carry = mpn_add_1(rp+s2n+s1n-1, rp+s2n+s1n-1, 1, carry); 
   return carry;
}

/*=====================================================================================

   Fast Integer Multiplication Code
   
=====================================================================================*/

unsigned long MUL_TWK_VALS[MUL_TWK_COUNT][3] = 
{
   {2000, 2140, 1024},
   {2140, 2430, 64},
   {2430, 2580, 1024},
   {2580, 2700, 64},
   {2700, 2880, 4096},
   {2880, 3850, 16},
   {3850, 4220, 4},
   {4220, 4400, 1024},
   {4400, 4850, 16},
   {4850, 5700, 1024},
   {5700, 7900, 4},
   {7900, 8900, 1024},
   {8900, 97000, 4},
   {97000, 127000, 1},
   {127000, 262000, 4},
   {262000, 517000, 1},
   {517000, 1050000, 4},
   {1050000, 2060000, 1},
   {2060000, 4230000, 4},
   {4230000, 8350000, 1}
};

unsigned long SQR_TWK_VALS[SQR_TWK_COUNT][3] = 
{
   {1564, 1994, 16},
   {1994, 2952, 64},
   {2952, 5921, 16},
   {5921, 32575, 4},
   {32575, 40006, 16},
   {40006, 66526, 4},
   {66526, 127370, 1},
   {127370, 257473, 4},
   {257473, 520507, 1},
   {520507, 1050000, 4},
   {1050000, 2060000, 1},
   {2060000, 4230000, 4},
   {4230000, 8350000, 1}
};

unsigned long FFT_MUL_TWK[FFT_MUL_COUNT][2] = 
{
   {1734, 8},
   {3450, 9},
   {3710, 8},
   {3840, 9},
   {4080, 8},
   {4220, 9},
   {4480, 8},
   {4600, 9},
   {4870, 8},
   {4980, 9},
   {5230, 8},
   {5380, 9},
   {6140, 10},
   {7150, 9},
   {7700, 10},
   {8700, 9},
   {9200, 10},
   {12300, 11},
   {34800, 10},
   {36900, 12},
   {41000, 11},
   {49200, 12},
   {99000, 13},
   {130000, 12},
   {198000, 13},
   {396000, 14},
   {526000, 13},
   {1040000, 14},
   {3150000, 15},
   {6300000, 16},
   {8400000, 15},
   {12700000, 17},
   {15900000, 16}
};

unsigned long FFT_SQR_TWK[FFT_SQR_COUNT][2] = 
{
   {1300, 8},
   {2700, 9},
   {2950, 8},
   {3080, 9},
   {4100, 8},
   {4240, 9},
   {5100, 10},
   {5630, 9},
   {6150, 10},
   {8700, 9},
   {9200, 10},
   {12300, 11},
   {24700, 12},
   {29000, 11},
   {35000, 10},
   {36900, 11},
   {49000, 12},
   {99000, 13},
   {130000, 12},
   {198000, 13},
   {396000, 14},
   {660000, 13},
   {780000, 14},
   {170000, 15},
   {210000, 14},
   {420000, 15},
   {880000, 16},
   {1060000, 15},
   {1260000, 17},
   {1680000, 16}
};

/*
   Splits an mpn into segments of length coeff_limbs and stores in a ZmodF_poly
   in zero padded coefficients of length output_limbs, for use in FFT 
   convolution code. Assumes that the input is total_limbs in length. 
   Used by the large integer multiplication code 
   (F_mpn_mul and F_mpn_mul_precache and F_mpn_mul_trunc)
*/

void F_mpn_FFT_split(ZmodF_poly_t poly, const mp_limb_t * limbs, const unsigned long total_limbs,
                               unsigned long coeff_limbs, unsigned long output_limbs)
{
   unsigned long length = (total_limbs-1)/coeff_limbs + 1;
   unsigned long i, j, skip;
   
   for (skip = 0, i = 0; skip+coeff_limbs <= total_limbs; skip+=coeff_limbs, i++)
   {
      if (i + 1 < length)
		 for (j = 0; j + 8 < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      
      F_mpn_clear(poly->coeffs[i], output_limbs+1);
      // convert a coefficient
      F_mpn_copy(poly->coeffs[i], limbs+skip, coeff_limbs);
   }
   if (i < length) F_mpn_clear(poly->coeffs[i], output_limbs+1);
   if (total_limbs > skip) F_mpn_copy(poly->coeffs[i], limbs+skip, total_limbs-skip);
   
   poly->length = length;
}

/*
   Splits an mpn into segments of length _bits_ and stores in a ZmodF_poly
   in zero padded coefficients of length output_limbs, for use in FFT 
   convolution code. Assumes that the input is total_limbs in length. 
   Used by the large integer multiplication code 
   (F_mpn_mul)
   
   It is assumed that bits is not divisible by FLINT_BITS
*/

void F_mpn_FFT_split_bits(ZmodF_poly_t poly, const mp_limb_t * limbs, const unsigned long total_limbs,
                               unsigned long bits, unsigned long output_limbs)
{
   unsigned long length = (FLINT_BITS*total_limbs-1)/bits + 1;
   unsigned long i, j;
   
   unsigned long top_bits = ((FLINT_BITS-1)&bits);
   if (top_bits == 0)
   {
      F_mpn_FFT_split(poly, limbs, total_limbs, bits >> FLINT_LG_BITS_PER_LIMB, output_limbs);
      return;
   }
   unsigned long coeff_limbs = (bits>>FLINT_LG_BITS_PER_LIMB) + 1;
   unsigned long mask = (1L<<top_bits)-1L;
   unsigned long shift_bits = 0L;
   unsigned long const * limb_ptr = limbs;                      
    
   for (i = 0; i < length - 1; i++)
   {
      for (j = 0; j + 8 < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      
      F_mpn_clear(poly->coeffs[i], output_limbs+1);
      // convert a coefficient
      if (!shift_bits)
      {
         F_mpn_copy(poly->coeffs[i], limb_ptr, coeff_limbs);
         poly->coeffs[i][coeff_limbs-1] &= mask;
         limb_ptr += (coeff_limbs-1);
         shift_bits += top_bits;
      } else
      {
         mpn_rshift(poly->coeffs[i], limb_ptr, coeff_limbs, shift_bits);
         limb_ptr += (coeff_limbs-1);
         shift_bits += top_bits;
         if (shift_bits >= FLINT_BITS)
         {
            limb_ptr++;
            poly->coeffs[i][coeff_limbs-1] += (limb_ptr[0] << (FLINT_BITS - (shift_bits - top_bits)));
            shift_bits -= FLINT_BITS; 
         }
         poly->coeffs[i][coeff_limbs-1] &= mask;
      }                      
   }
   
   F_mpn_clear(poly->coeffs[i], output_limbs+1);
   unsigned long limbs_left = total_limbs - (limb_ptr - limbs);
   if (!shift_bits)
   {
      F_mpn_copy(poly->coeffs[i], limb_ptr, limbs_left);
   } else
   {
      mpn_rshift(poly->coeffs[i], limb_ptr, limbs_left, shift_bits);
   }                      
      
   poly->length = length;
}

/*
   Recombines coefficients of a ZmodF_poly after doing a convolution. Assumes 
   each of the coefficients of the ZmodF_poly is output_limbs long, that each 
   of the coefficients is being shifted by a multiple of coeff_limbs and added
   to an mpn which is total_limbs long. It is assumed that the mpn has been 
   zeroed in advance.
   Used by the large integer multiplication code (F_mpn_mul and F_mpn_mul_precache
   and F_mpn_mul_trunc)
*/

void F_mpn_FFT_combine(mp_limb_t * res, ZmodF_poly_t poly, unsigned long coeff_limbs, 
                             unsigned long output_limbs, unsigned long total_limbs)
{
   unsigned long skip, i, j;
   unsigned long length = poly->length;
   
   for (skip = 0, i = 0; (i < length) && (skip+output_limbs <= total_limbs); i++, skip+=coeff_limbs)
   { 
      for (j = 0; j < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      mpn_add(res+skip, res+skip, output_limbs+1, poly->coeffs[i], output_limbs);      
   } 
   while ((skip < total_limbs) && (i < length))
   {
      mpn_add(res+skip, res+skip, total_limbs - skip, poly->coeffs[i], FLINT_MIN(total_limbs - skip, output_limbs));
      i++;
      skip+=coeff_limbs;
   }  

}

/*
   Recombines coefficients of a ZmodF_poly after doing a convolution. Assumes 
   each of the coefficients of the ZmodF_poly is output_limbs long, that each 
   of the coefficients is being shifted by a multiple of _bits_ and added
   to an mpn which is total_limbs long. It is assumed that the mpn has been 
   zeroed in advance.
   Used by the large integer multiplication code (F_mpn_mul)
   
   It is assumed that bits is not divisible by FLINT_BITS
*/

void F_mpn_FFT_combine_bits(mp_limb_t * res, ZmodF_poly_t poly, unsigned long bits, 
                             unsigned long output_limbs, unsigned long total_limbs)
{
   unsigned long top_bits = ((FLINT_BITS-1)&bits);
   if (top_bits == 0)
   {
      F_mpn_FFT_combine(res, poly, bits >> FLINT_LG_BITS_PER_LIMB, output_limbs, total_limbs);
      return;
   }
   
   unsigned long coeff_limbs = (bits>>FLINT_LG_BITS_PER_LIMB) + 1;
   unsigned long i, j;
   unsigned long length = poly->length;
   unsigned long * temp = (unsigned long *) flint_heap_alloc(output_limbs+1);
   unsigned long shift_bits = 0;
   unsigned long * limb_ptr = res;
   unsigned long * end = res + total_limbs;
   
   for (i = 0; (i < length) && (limb_ptr + output_limbs < end); i++)
   { 
      for (j = 0; j < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      if (shift_bits)
      {
         mpn_lshift(temp, poly->coeffs[i], output_limbs+1, shift_bits);
         mpn_add_n(limb_ptr, limb_ptr, temp, output_limbs+1);
      } else
      {
         mpn_add(limb_ptr, limb_ptr, output_limbs+1, poly->coeffs[i], output_limbs);
      }
      shift_bits += top_bits;
      limb_ptr += (coeff_limbs - 1);
      if (shift_bits >= FLINT_BITS)
      {
         limb_ptr++;
         shift_bits -= FLINT_BITS;
      }      
   } 
   while ((limb_ptr < end) && (i < length))
   {
      if (shift_bits)
      {
         mpn_lshift(temp, poly->coeffs[i], output_limbs+1, shift_bits);
         mpn_add_n(limb_ptr, limb_ptr, temp, end - limb_ptr);
      } else
      {
         mpn_add_n(limb_ptr, limb_ptr, poly->coeffs[i], end - limb_ptr);
      }
      shift_bits += top_bits;
      limb_ptr += (coeff_limbs - 1);
      if (shift_bits >= FLINT_BITS)
      {
         limb_ptr++;
         shift_bits -= FLINT_BITS;
      }  
      i++;    
   }
   
   flint_heap_free(temp);      
}

/*
   Compute optimal lengths for the polynomials that coeff1 and coeff2 are broken into
   in the convolution based long integer code.
   We want the sum of the two lengths to satisfy the SS condition (with sqrt2):
      
      if 2^l1 < length1 + length2 <= 2^l2 then 2^(l2-1) divides output_bits
      
   Requires limbs1 and limbs2 are at least 1, ensures length1 and length2 are at least 1
*/

#define F_mpn_mul_ADJUST \
do { \
   /* Compute the coefficient size for breaking the two long integers up */ \
   coeff_limbs = (limbs1+limbs2-1)/(length)+1; \
   if (coeff_limbs == 1L) /* This is as far as we can go */ \
   \
   { \
      length1 = limbs1; \
      length2 = limbs2; \
      done = 1; \
   } \
   while ((limbs1-1)/(coeff_limbs)+(limbs2-1)/(coeff_limbs)+2 > length) coeff_limbs++; \
   /* Compute the number of bits for the output coefficients */ \
   output_bits = (2*coeff_limbs+1)*FLINT_BITS; \
   output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1); \
   /* Try and compute a more optimal coefficient size to break up inputs */ \
   coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2; \
   if ((long) coeff_limbs <= 1L) coeff_limbs = 1; \
   /* Compute the lengths of the polys the coefficients will be broken into with this coeff size */ \
   length1 = (limbs1-1)/coeff_limbs+1; \
   length2 = (limbs2-1)/coeff_limbs+1; \
} while (0)

#define F_mpn_mul_TUNING \
do { \
   if (twk > 64) \
   { \
      length = 2; \
      log_length = 1; \
      \
      int done = 0; \
      \
      while ((length < 2*output_bits) && !done) \
      { \
         /* We are outside the optimal SS region, so double the length */ \
         length<<=1; \
         log_length++; \
         F_mpn_mul_ADJUST; \
      } \
      \
      while ((twk > 64) && (length >= 4)) \
      { \
         log_length--; \
         length>>=1; \
         twk>>=2; \
      } \
      \
      F_mpn_mul_ADJUST; \
      \
   } else \
   { \
      int done = 0; \
      \
      while ((twk*length < 2*output_bits) && !done) \
      { \
         /* We are outside the optimal SS region, so double the length */ \
         length<<=1; \
         log_length++; \
         F_mpn_mul_ADJUST; \
      } \
   } \
} while (0)

mp_limb_t __F_mpn_mul(mp_limb_t * res, const mp_limb_t * data1, const unsigned long limbs1, 
                                      const mp_limb_t * data2, const unsigned long limbs2, unsigned long log_length)
{
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long s1 = (FLINT_BIT_COUNT(data1[limbs1-1]) + FLINT_BIT_COUNT(data2[limbs2-1]) <= FLINT_BITS);
   unsigned long total_limbs = coeff_limbs - s1;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 1;

   unsigned long bits;
      
   do
   {
      bits = (((limbs1 << FLINT_LG_BITS_PER_LIMB)-1) >> (log_length-1)) + 1;
      output_bits = 2*bits + log_length2;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   
      bits = (output_bits - log_length2)/2;
      length1 = ((limbs1 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      length2 = ((limbs2 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      log_length2++;
   } while ((length2 > (1L<<(log_length2-1))) || (length1 > (1L<<(log_length-1))));
   
   n = (output_bits-1)/FLINT_BITS+1;
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld, %ld, %ld\n", bits, length1, length2, output_bits, coeff_limbs, n, log_length);
#endif   
   ZmodF_poly_t poly1;
   ZmodF_poly_init(poly1, log_length, n, 1);
   F_mpn_FFT_split_bits(poly1, data1, limbs1, bits, n);
   
   ulong length = length1 + length2 - 1;
   ulong size = 1UL << log_length;
   if (length > size)
   length = size;

   ZmodF_poly_FFT(poly1, length);
   
	if ((data1 == data2) && (limbs1 == limbs2))
   {
      ZmodF_poly_pointwise_mul(poly1, poly1, poly1);
	} else
	{
		ZmodF_poly_t poly2;
      ZmodF_poly_init(poly2, log_length, n, 1);
      F_mpn_FFT_split_bits(poly2, data2, limbs2, bits, n);

      ZmodF_poly_FFT(poly2, length);
      
      ZmodF_poly_pointwise_mul(poly1, poly1, poly2);
	   ZmodF_poly_clear(poly2);
	}	
		
   ZmodF_poly_IFFT(poly1);
   ZmodF_poly_rescale(poly1);   
   
   ZmodF_poly_normalise(poly1);
   
   F_mpn_clear(res, limbs1+limbs2);
   
   F_mpn_FFT_combine_bits(res, poly1, bits, n, total_limbs);
   ZmodF_poly_clear(poly1);
   
   return res[limbs1+limbs2-1];
}

mp_limb_t __F_mpn_mul_trunc(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2, 
                                      unsigned long log_length, unsigned long trunc)
{
   unsigned long length = 1;
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long s1 = (FLINT_BIT_COUNT(data1[limbs1-1]) + FLINT_BIT_COUNT(data2[limbs2-1]) <= FLINT_BITS);
   unsigned long total_limbs = coeff_limbs - s1;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 1;

   unsigned long bits;
      
   do
   {
      bits = (((limbs1 << FLINT_LG_BITS_PER_LIMB)-1) >> (log_length-1)) + 1;
      output_bits = 2*bits + log_length2;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   
      bits = (output_bits - log_length2)/2;
      length1 = ((limbs1 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      length2 = ((limbs2 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      log_length2++;
   } while ((length2 > (1L<<(log_length2-1))) || (length1 > (1L<<(log_length-1))));
   
   n = (output_bits-1)/FLINT_BITS+1;
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld, %ld\n", bits, length1, length2, output_bits, coeff_limbs, n);
#endif   
   ZmodF_poly_t poly1;
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   F_mpn_FFT_split_bits(poly1, data1, limbs1, bits, n);
   
   if ((data1 == data2) && (limbs1 == limbs2))
   {
      // identical operands case
      ZmodF_poly_convolution_range(poly1, poly1, poly1, 0, (trunc*FLINT_BITS-1)/bits+1);
   }
   else
   {
      // distinct operands case
      ZmodF_poly_t poly2;
      ZmodF_poly_stack_init(poly2, log_length, n, 1);
      F_mpn_FFT_split_bits(poly2, data2, limbs2, bits, n);

      ZmodF_poly_convolution_range(poly1, poly1, poly2, 0, (trunc*FLINT_BITS-1)/bits+1);

      ZmodF_poly_stack_clear(poly2);
   }
   
   poly1->length = FLINT_MIN(poly1->length, (trunc*FLINT_BITS-1)/bits+1);
   ZmodF_poly_normalise(poly1);
   
   F_mpn_clear(res, trunc);
   
   F_mpn_FFT_combine_bits(res, poly1, bits, n, trunc);
   ZmodF_poly_stack_clear(poly1);
   
   return res[trunc-1];
}

/*
   Multiply two integers in mpn format
   
   WARNING: This function requires limbs1+limbs2 output limbs when limbs1+limbs2
   < FLINT_FFT_LIMBS_CROSSOVER but may require one less limb otherwise. The function
   will return 0 if it did not require (and indeed did not zero) the extra limb,
   otherwise it returns the (non zero) value of this high limb after multiplication.
   
   Assumes neither of limbs1, limbs2 is zero. 
*/

mp_limb_t F_mpn_mul(mp_limb_t * res, const mp_limb_t * data1, const unsigned long limbs1, 
                                      const mp_limb_t * data2, const unsigned long limbs2)
{
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long twk;
   
   if (coeff_limbs/2 > FFT_TUNE_CUTOFF)
   {
      twk = 0;
      while ((1L<<(2*twk)) < FLINT_BITS*coeff_limbs)
      {
         twk++;
      }   
   } else if ((data1 != data2) || (limbs1 != limbs2))
   {
      if (coeff_limbs/2 < FFT_MUL_TWK[0][0]) 
         return mpn_mul(res, data1, limbs1, data2, limbs2);
      else
      {
         unsigned long i = 0;
         while ((i < FFT_MUL_COUNT-1) && (coeff_limbs/2 > FFT_MUL_TWK[i+1][0])) i++;
         twk = FFT_MUL_TWK[i][1];
      }
   } else
   {
      if (coeff_limbs/2 < FFT_SQR_TWK[0][0]) 
         return mpn_mul(res, data1, limbs1, data1, limbs1);
      else
      {
         unsigned long i = 0;
         while ((i < FFT_SQR_COUNT-1) && (coeff_limbs/2 > FFT_SQR_TWK[i+1][0])) i++;
         twk = FFT_SQR_TWK[i][1];
      }
   }


   return __F_mpn_mul(res, data1, limbs1, data2, limbs2, twk);
}

/*
   Multiply two integers in mpn format truncating to _trunc_ output limbs
   Assumes none of limbs1, limbs2 and trunc is zero. 
   
   WARNING: the output "res" needs to have "limbs1+limbs2"
   limbs allocated and we must have limbs1 >= limbs2 >= 1.
*/

mp_limb_t F_mpn_mul_trunc(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                        mp_limb_t * data2, unsigned long limbs2, unsigned long trunc)
{
   unsigned long coeff_limbs = limbs1 + limbs2;
   if (trunc > coeff_limbs) trunc = coeff_limbs;
   unsigned long twk;
   
   if (coeff_limbs/2 > FFT_TUNE_CUTOFF)
   {
      twk = 0;
      while ((1L<<(2*twk)) < FLINT_BITS*coeff_limbs)
      {
         twk++;
      }   
   } else if ((data1 != data2) || (limbs1 != limbs2))
   {
      if (coeff_limbs/2 < FFT_MUL_TWK[0][0]) 
      {
         mpn_mul(res, data1, limbs1, data2, limbs2);
         return res[trunc-1];
      } else
      {
         unsigned long i = 0;
         while ((i < FFT_MUL_COUNT-1) && (coeff_limbs/2 > FFT_MUL_TWK[i+1][0])) i++;
         twk = FFT_MUL_TWK[i][1];
      }
   } else
   {
      if (coeff_limbs/2 < FFT_SQR_TWK[0][0]) 
      {
         mpn_mul(res, data1, limbs1, data1, limbs1);
         return res[trunc-1];
      } else
      {
         unsigned long i = 0;
         while ((i < FFT_SQR_COUNT-1) && (coeff_limbs/2 > FFT_SQR_TWK[i+1][0])) i++;
         twk = FFT_SQR_TWK[i][1];
      }
   }


   return __F_mpn_mul_trunc(res, data1, limbs1, data2, limbs2, twk, trunc);
}

/*   
   Precompute an FFT for integer multiplication.
   Assumes neither of limbs1, limbs2 is zero. 
*/

void F_mpn_mul_precache_init(F_mpn_precache_t precache, mp_limb_t * data1, unsigned long limbs1, unsigned long limbs2)
{
   if (limbs1 == 0)
   {
      precache->poly = NULL;
      return;
   }

   int swapped = 0;
   if (limbs2 > limbs1)
   {
      unsigned long temp = limbs1;
      limbs1 = limbs2;
      limbs2 = temp;
      swapped = 1;
   }

   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long log_length;
   
   if (coeff_limbs/2 > FFT_TUNE_CUTOFF)
   {
      log_length = 0;
      while ((1L<<(2*log_length)) < FLINT_BITS*coeff_limbs)
      {
         log_length++;
      }   
   } else
   {
      unsigned long i = 0;
      while ((i < FFT_SQR_COUNT-1) && (coeff_limbs/2 > FFT_SQR_TWK[i+1][0])) i++;
      log_length = FFT_SQR_TWK[i][1];
   }
   
   unsigned long length = 1;
   
   unsigned long total_limbs = coeff_limbs;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 1;
   
   unsigned long bits;
      
   do
   {
      bits = (((limbs1 << FLINT_LG_BITS_PER_LIMB)-1) >> (log_length-1)) + 1;
      output_bits = 2*bits + log_length2;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   
      bits = (output_bits - log_length2)/2;
      length1 = ((limbs1 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      length2 = ((limbs2 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      log_length2++;
   } while ((length2 > (1L<<(log_length2-1))) || (length1 > (1L<<(log_length-1))));
   
   n = (output_bits-1)/FLINT_BITS+1;

   if (swapped)
   {
      unsigned long temp = limbs1;
      limbs1 = limbs2;
      limbs2 = temp;
      temp = length1;
      length1 = length2;
      length2 = temp; 
   }
      
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodF_poly_p poly1;
   poly1 = (ZmodF_poly_p) malloc(sizeof(ZmodF_poly_struct));
   ZmodF_poly_init(poly1, log_length, n, 1);
   F_mpn_FFT_split_bits(poly1, data1, limbs1, bits, n);
   
   unsigned long size = (1L<<poly1->depth);
   ZmodF_poly_FFT(poly1, size);
   precache->type = FFT_PRE;
   precache->bits = bits;
   precache->length = length1;
   precache->length2 = length2;
   precache->coeff_limbs = coeff_limbs;
   precache->limbs1 = limbs1;
   precache->limbs2 = limbs2;
   precache->poly = poly1; 
   precache->msl_bits = FLINT_BIT_COUNT(data1[limbs1-1]);
}

void F_mpn_mul_precache_clear(F_mpn_precache_t precache)
{
   if (precache->type == FFT_PRE) 
   {
      if (precache->poly) 
      {
         ZmodF_poly_clear(precache->poly);
         free(precache->poly);
      }
   }   
}

/*   
   Compute an integer multiplication given a precomputed FFT for 
   one of the integers.
   Assumes neither of limbs1, limbs2 is zero. 
*/

mp_limb_t F_mpn_mul_precache(mp_limb_t * res, mp_limb_t * data2, unsigned long limbs2, F_mpn_precache_t precache)
{
   ZmodF_poly_t poly2;
   ZmodF_poly_stack_init(poly2, precache->poly->depth, precache->poly->n, 1);
   int s1 = (FLINT_BIT_COUNT(data2[limbs2-1]) + precache->msl_bits <= FLINT_BITS);
   
   F_mpn_FFT_split_bits(poly2, data2, limbs2, precache->bits, precache->poly->n);
   
   ZmodF_poly_FFT(poly2, precache->length+poly2->length-1);
   ZmodF_poly_pointwise_mul(poly2, poly2, precache->poly);
   ZmodF_poly_IFFT(poly2);
   ZmodF_poly_rescale(poly2);
   
   ZmodF_poly_normalise(poly2);
   F_mpn_clear(res, precache->limbs1 + limbs2 - s1);
   
   F_mpn_FFT_combine_bits(res, poly2, precache->bits, precache->poly->n, precache->limbs1 + limbs2 - s1);
   
   ZmodF_poly_stack_clear(poly2);
   
   if (s1) return 0;
   else return res[precache->limbs1+limbs2-1];
}

mp_limb_t F_mpn_mul_precache_trunc(mp_limb_t * res, mp_limb_t * data2, unsigned long limbs2, F_mpn_precache_t precache, unsigned long trunc)
{
   if (trunc == 0) return 0;
   ZmodF_poly_t poly2;
   ZmodF_poly_stack_init(poly2, precache->poly->depth, precache->poly->n, 1);
   int s1 = (FLINT_BIT_COUNT(data2[limbs2-1]) + precache->msl_bits <= FLINT_BITS);
   if (trunc > precache->limbs1+limbs2 - s1) trunc = precache->limbs1+limbs2 - s1;
   F_mpn_FFT_split_bits(poly2, data2, limbs2, precache->bits, precache->poly->n);
   
   ZmodF_poly_FFT(poly2, precache->length+poly2->length-1);
   ZmodF_poly_pointwise_mul(poly2, poly2, precache->poly);
   ZmodF_poly_IFFT(poly2);
   ZmodF_poly_rescale_range(poly2, 0, (trunc*FLINT_BITS-1)/precache->bits+1);
   poly2->length = FLINT_MIN(poly2->length,(trunc*FLINT_BITS-1)/precache->bits+1);
   ZmodF_poly_normalise(poly2);
   F_mpn_clear(res, precache->limbs1 + limbs2);
   
   F_mpn_FFT_combine_bits(res, poly2, precache->bits, precache->poly->n, trunc);
   
   ZmodF_poly_stack_clear(poly2);
   
   return res[trunc-1];
}

/*
   Multiply integers data1 of length limbs1 by data2 of length limbs2 and return limbs
	[start... trunc) of the output assuming that limbs2 <= limbs1/2
*/

mp_limb_t __F_mpn_mul_middle(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2, 
                                      unsigned long start, unsigned long trunc)
{
   unsigned long coeff_limbs = trunc;
   unsigned long twk;
   
   if (coeff_limbs/2 > FFT_TUNE_CUTOFF)
   {
      twk = 0;
      while ((1L<<(2*twk)) < FLINT_BITS*coeff_limbs)
      {
         twk++;
      }   
   } else if ((data1 != data2) || (limbs1 != limbs2))
   {
      if (coeff_limbs/2 < FFT_MUL_TWK[0][0]) 
      {
         mpn_mul(res, data1, limbs1, data2, limbs2);
         return res[trunc-1];
      } else
      {
         unsigned long i = 0;
         while ((i < FFT_MUL_COUNT-1) && (coeff_limbs/2 > FFT_MUL_TWK[i+1][0])) i++;
         twk = FFT_MUL_TWK[i][1];
      }
   } else
   {
      if (coeff_limbs/2 < FFT_SQR_TWK[0][0]) 
      {
         mpn_mul(res, data1, limbs1, data1, limbs1);
         return res[trunc-1];
      } else
      {
         unsigned long i = 0;
         while ((i < FFT_SQR_COUNT-1) && (coeff_limbs/2 > FFT_SQR_TWK[i+1][0])) i++;
         twk = FFT_SQR_TWK[i][1];
      }
   }

   unsigned long log_length = twk;
   unsigned long length = 1;
   
   unsigned long total_limbs = coeff_limbs;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 1;

   unsigned long bits;
      
   do
   {
      bits = (((limbs1 << FLINT_LG_BITS_PER_LIMB)-1) >> (log_length)) + 1;
      output_bits = 2*bits + log_length2;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   
      bits = (output_bits - log_length2)/2;
      length1 = ((limbs1 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      length2 = ((limbs2 << FLINT_LG_BITS_PER_LIMB)-1)/bits + 1;
      log_length2++;
   } while ((length2 > (1L<<(log_length2))) || (length1 > (1L<<(log_length))) || (length1 + length2 > (1L<<log_length) + (1L<<(log_length-1))));
   
   n = (output_bits-1)/FLINT_BITS+1;

	// We break limbs1 and limbs2 into length1 and length2 chunks of _bits_ bits 
	// such that length1 <= 2^log_length and we set output_bits so that it is
   // divisible by 2^(log_length-1) and >= 2*bits + log_length2
	// we make nB >= output_bits
	// we also ensure that no more than 1/3 of the product wraps around

#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld, %ld, %ld\n", bits, length1, length2, output_bits, coeff_limbs, n, log_length);
#endif   
   ZmodF_poly_t poly1;
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   F_mpn_FFT_split_bits(poly1, data1, limbs1, bits, n);
   
	// we perform a convolution of length 2^log_length
	// if length2 <= 2^(log_length-1) then the total length of a full product would be 
	// < 2^log_length + 2^(log_length-1)
	// as the FFT will wrap around, the bottom 2^(log_length-1)-1 terms may be messed up

	// we are only interested in limbs [start, trunc) of the output and as each of the FFT
	// coefficients will be staggered by _bits_ bits when added to the output, we are only
	// interested in the terms (start*FLINT_BITS - output_bits)/bits-1 ... 
	// (trunc*FLINT_BITS-1)/bits+1 of the output of the FFT
	// note it is not (start*FLINT_BITS)/bits as the previous terms are output_bits wide 
	// and so may overlap the limb we want 

   long first_term = (start*FLINT_BITS - output_bits)/bits;
	if (first_term < 0L) first_term = 0L;
	
	if ((data1 == data2) && (limbs1 == limbs2))
   {
      // identical operands case
      ZmodF_poly_convolution_range(poly1, poly1, poly1, first_term, (trunc*FLINT_BITS-1)/bits+1);
   }
   else
   {
      // distinct operands case
      ZmodF_poly_t poly2;
      ZmodF_poly_stack_init(poly2, log_length, n, 1);
      F_mpn_FFT_split_bits(poly2, data2, limbs2, bits, n);

      ZmodF_poly_convolution_range(poly1, poly1, poly2, first_term, (trunc*FLINT_BITS-1)/bits+1);

      ZmodF_poly_stack_clear(poly2);
   }
   
   poly1->length = (trunc*FLINT_BITS-1)/bits+1;
   ZmodF_poly_normalise(poly1);
   
   F_mpn_clear(res, trunc);
   
   // When we combine the FFT terms it is possible that there is a carry out of the 
	// top term which wrapped around, which may add 1 to the first limb we are interested in
	// However suppose we are in the special case where the integers being multiplied are 
	// Kronecker segmented versions of polynomials of length 2n and n respectively, where
	// each coefficient is packed into _bits_ bits
	// Then the top _bits_ bits of the product of the integers is necessarily 0, as the length 
	// of the product is only 3n-1
	// Thus to prove that a carry does not occur into the terms we are interested in, it is
	// sufficient to note that the n-th term of the product cannot be 2^bits - 1 (the only
	// situation where a carry from earlier overlaps could cause a carry out of the n-th term)
	// Suppose the number of bits of the original coefficients of the polynomials is b and 
	// that n <= 2^k so that 2*b + k <= bits
	// Then the n-th coefficient of the product polynomial can be at most 2^k * (2^b-1)^2. QED.

	F_mpn_FFT_combine_bits(res, poly1, bits, n, trunc);
   ZmodF_poly_stack_clear(poly1);
   
   return res[trunc-1];
}

/*
   This function is intended to be used with F_mpn_mul_precache_init
	It computes the product as per __F_mpn_mul_middle, however there is no wrap around
	of the FFT as F_mpn_mul_precache_init precomputed the FFT for use with a full 
	precomputed product
	The saving here is that only the required coefficients are computed, and of course
	the precomputed FFT can be reused
*/

mp_limb_t __F_mpn_mul_middle_precache(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      F_mpn_precache_t precache, 
                                      unsigned long start, unsigned long trunc)
{
   ZmodF_poly_t poly1;
   ZmodF_poly_stack_init(poly1, precache->poly->depth, precache->poly->n, 1);
   F_mpn_FFT_split_bits(poly1, data1, limbs1, precache->bits, precache->poly->n);
   
   unsigned long length = precache->poly->length + poly1->length - 1;
	unsigned long log_length2 = ceil_log2(poly1->length);
	unsigned long output_bits = 2*precache->bits + log_length2;
   unsigned long size = (1L<<precache->poly->depth);
   if (length > size) length = size;
   ZmodF_poly_FFT(poly1, length);
   ZmodF_poly_pointwise_mul(poly1, poly1, precache->poly);
   ZmodF_poly_IFFT(poly1);
   ZmodF_poly_rescale_range(poly1, (start*FLINT_BITS - output_bits)/precache->bits, (trunc*FLINT_BITS-1)/precache->bits+1);
   
   poly1->length = FLINT_MIN(poly1->length, (trunc*FLINT_BITS-1)/precache->bits+1);
   ZmodF_poly_normalise(poly1);
   
   F_mpn_clear(res, trunc);
   
   F_mpn_FFT_combine_bits(res, poly1, precache->bits, precache->poly->n, trunc);
   ZmodF_poly_stack_clear(poly1);
   
   return res[trunc-1];
}



