/****************************************************************************

Zpoly_mpn.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint-manager.h"
#include "mpn_extras.h"
#include "ZmodFpoly.h"
#include "Z_mpn.h"


void Z_split_limbs(ZmodFpoly_t poly, mp_limb_t * limbs, unsigned long total_limbs,
                               unsigned long coeff_limbs, unsigned long output_limbs)
{
   unsigned long length = (total_limbs-1)/coeff_limbs + 1;
   unsigned long i, j, skip;
   
   for (skip = 0, i = 0; skip+coeff_limbs <= total_limbs; skip+=coeff_limbs, i++)
   {
      for (j = 0; j < output_limbs; j += 8) FLINT_PREFETCH(poly->coeffs[i+1], j);
      
      clear_limbs(poly->coeffs[i], output_limbs+1);
      // convert a coefficient
      copy_limbs(poly->coeffs[i], limbs+skip, coeff_limbs);
   }
   if (i < length) clear_limbs(poly->coeffs[i], output_limbs+1);
   if (total_limbs > skip) copy_limbs(poly->coeffs[i], limbs+skip, total_limbs-skip);
   i++;
   
   poly->length = length;
}

void Z_combine_limbs(mp_limb_t * res, ZmodFpoly_t poly, unsigned long coeff_limbs, 
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

void Z_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2)
{
   unsigned long length = 1;
   unsigned long log_length = 0;
   
   unsigned long output_limbs = limbs1+limbs2;
   
   unsigned long twk = 8;
   
   while (twk*length < output_limbs*FLINT_BITS_PER_LIMB)
   {
      length<<=1;
      log_length++;
      output_limbs = 2*(limbs1+limbs2)/(length+1)+1;
   }
   
   unsigned long coeff_limbs = (limbs1 + limbs2-1)/length+1;
   if (limbs1%coeff_limbs) coeff_limbs++;
   
   unsigned long length1 = (limbs1-1)/coeff_limbs + 1;
   unsigned long length2 = (limbs2-1)/coeff_limbs + 1;
   
   // Recompute output_limbs
   unsigned log_length2 = 0;
   while ((1<<log_length2) < length2) log_length2++;
    
   unsigned long output_bits = 2*coeff_limbs*FLINT_BITS_PER_LIMB + log_length2;
   output_limbs = (output_bits-1)/FLINT_BITS_PER_LIMB+1;
   
   // Round to a number of bits supported by the convolution length
   output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   unsigned long n = (output_bits - 1) / FLINT_BITS_PER_LIMB + 1;
   
   ZmodFpoly_t poly1;
   ZmodFpoly_init(poly1, log_length, n, 1);
   Z_split_limbs(poly1, data1, limbs1, coeff_limbs, n);
   
   if (data1 == data2 && limbs1 == limbs2)
   {
      // identical operands case
      ZmodFpoly_convolution(poly1, poly1, poly1);
   }
   else
   {
      // distinct operands case
      ZmodFpoly_t poly2;
      ZmodFpoly_init(poly2, log_length, n, 1);
      Z_split_limbs(poly2, data2, limbs2, coeff_limbs, n);

      ZmodFpoly_convolution(poly1, poly1, poly2);

      ZmodFpoly_clear(poly2);
   }
   
   ZmodFpoly_normalise(poly1);
   
   clear_limbs(res, limbs1 + limbs2);
   
   Z_combine_limbs(res, poly1, coeff_limbs, output_limbs, limbs1 + limbs2);
   ZmodFpoly_clear(poly1);
}
