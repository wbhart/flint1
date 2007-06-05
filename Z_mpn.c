/****************************************************************************

Zpoly_mpn.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory-manager.h"
#include "mpn_extras.h"
#include "ZmodFpoly.h"
#include "Z_mpn.h"
#include "ZmodF_mul.h"
#include "Z_mpn_mul-tuning.h"

#define DEBUG 0

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
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long output_bits = coeff_limbs*FLINT_BITS_PER_LIMB;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 0;
   
   unsigned long twk;
   
   if (data1 != data2)
   {
      if (coeff_limbs < MUL_TWK_SMALL_CUTOFF) twk = MUL_TWK_SMALL_DEFAULT;
      else if (coeff_limbs > MUL_TWK_LARGE_CUTOFF) twk = MUL_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < MUL_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs < MUL_TWK_VALS[twk_count][0]) || (coeff_limbs < MUL_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = MUL_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   } else
   {
      if (coeff_limbs < SQR_TWK_SMALL_CUTOFF) twk = SQR_TWK_SMALL_DEFAULT;
      else if (coeff_limbs > SQR_TWK_LARGE_CUTOFF) twk = SQR_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < SQR_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs < SQR_TWK_VALS[twk_count][0]) || (coeff_limbs < SQR_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = SQR_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   }
   
   while (twk*length < 2*output_bits)
   {
      length<<=1;
      log_length++;
      coeff_limbs = (limbs1+limbs2-1)/length+1;
      while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
      output_bits = (2*coeff_limbs+1)*FLINT_BITS_PER_LIMB;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      coeff_limbs = ((output_bits - FLINT_BITS_PER_LIMB)/FLINT_BITS_PER_LIMB)/2;
      if ((long) coeff_limbs < 1) coeff_limbs = 1;
      length1 = (limbs1-1)/coeff_limbs+1;
      length2 = (limbs2-1)/coeff_limbs+1;
   }
      
   n = output_bits/FLINT_BITS_PER_LIMB;
   n = ZmodF_mul_precomp_get_feasible_n(NULL, n);
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodFpoly_t poly1;
   ZmodFpoly_stack_init(poly1, log_length, n, 1);
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
      ZmodFpoly_stack_init(poly2, log_length, n, 1);
      Z_split_limbs(poly2, data2, limbs2, coeff_limbs, n);

      ZmodFpoly_convolution(poly1, poly1, poly2);

      ZmodFpoly_stack_clear(poly2);
   }
   
   ZmodFpoly_normalise(poly1);
   
   clear_limbs(res, limbs1 + limbs2);
   
   Z_combine_limbs(res, poly1, coeff_limbs, 2*coeff_limbs+1, limbs1 + limbs2);
   ZmodFpoly_stack_clear(poly1);
}

void Z_mul(mpz_t res, mpz_t a, mpz_t b)
{
   unsigned long int limbs;
   if (a->_mp_size + b->_mp_size > 128000/FLINT_BITS_PER_LIMB)
   {
      if (a->_mp_size >= b->_mp_size) limbs = a->_mp_size;
      else limbs = b->_mp_size;
      mp_limb_t* output = (mp_limb_t*) flint_stack_alloc(a->_mp_size+b->_mp_size);
      Z_mpn_mul(output, a->_mp_d, a->_mp_size, b->_mp_d, b->_mp_size);
      mpz_import(res, a->_mp_size+b->_mp_size, -1, sizeof(mp_limb_t), 0, 0, output);
      if (mpz_sgn(res) != mpz_sgn(a)*mpz_sgn(b)) mpz_neg(res,res);
      flint_stack_release();
   } else mpz_mul(res, a, b);
}
