/****************************************************************************

Z_mpn.c: Z arithmetic

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory-manager.h"
#include "mpn_extras.h"
#include "ZmodF_poly.h"
#include "Z_mpn.h"
#include "ZmodF_mul.h"
#include "Z_mpn_mul-tuning.h"

#define DEBUG 0

void Z_split_limbs(ZmodF_poly_t poly, mp_limb_t * limbs, unsigned long total_limbs,
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

void Z_combine_limbs(mp_limb_t * res, ZmodF_poly_t poly, unsigned long coeff_limbs, 
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

mp_limb_t Z_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2)
{
   unsigned long length = 1;
   unsigned long log_length = 0;
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 0;
   
   unsigned long twk;
   
   if (coeff_limbs/2 < 2300) 
   {
      return mpn_mul(res, data1, limbs1, data2, limbs2);
   }
   
   if (data1 != data2)
   {
      if (coeff_limbs/2 < MUL_TWK_SMALL_CUTOFF) twk = MUL_TWK_SMALL_DEFAULT;
      else if (coeff_limbs/2 > MUL_TWK_LARGE_CUTOFF) twk = MUL_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < MUL_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs/2 < MUL_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < MUL_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = MUL_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   } else
   {
      if (coeff_limbs/2 < SQR_TWK_SMALL_CUTOFF) twk = SQR_TWK_SMALL_DEFAULT;
      else if (coeff_limbs/2 > SQR_TWK_LARGE_CUTOFF) twk = SQR_TWK_LARGE_DEFAULT;
      else 
      {  
         for (unsigned long twk_count = 0; twk_count < SQR_TWK_COUNT; twk_count++)
         {
            if ((coeff_limbs/2 < SQR_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < SQR_TWK_VALS[twk_count][1])) continue;
            else 
            {
               twk = SQR_TWK_VALS[twk_count][2];
               break;
            }
         }
      }
   }
   
   if (twk > 64)
   {
      length = 2;
      log_length = 1;
      while ((1<<(log_length-1)) < output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
      while (twk > 64)
      {
         log_length--;
         length>>=1;
         twk>>=2;
      }
      coeff_limbs = (limbs1+limbs2-1)/length+1;
      while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
      output_bits = (2*coeff_limbs+1)*FLINT_BITS;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      while ((output_bits%3) != 0) output_bits+=(1<<(log_length-1));
      coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
      if ((long) coeff_limbs < 1) coeff_limbs = 1;
      length1 = (limbs1-1)/coeff_limbs+1;
      length2 = (limbs2-1)/coeff_limbs+1;
      log_length = 1;
      while ((1<<log_length) < length1 + length2) log_length++;
      length = (1<<log_length);        
   }
   else
   {
      while (twk*length < 2*output_bits)
      {
         length<<=1;
         log_length++;
         coeff_limbs = (limbs1+limbs2-1)/length+1;
         while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
         output_bits = (2*coeff_limbs+1)*FLINT_BITS;
         output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
         coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
         if ((long) coeff_limbs < 1) coeff_limbs = 1;
         length1 = (limbs1-1)/coeff_limbs+1;
         length2 = (limbs2-1)/coeff_limbs+1;
      }
   }
         
   n = output_bits/FLINT_BITS;
   //printf("n= %ld\n",n);
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodF_poly_t poly1;
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   Z_split_limbs(poly1, data1, limbs1, coeff_limbs, n);
   
   if (data1 == data2 && limbs1 == limbs2)
   {
      // identical operands case
      ZmodF_poly_convolution(poly1, poly1, poly1);
   }
   else
   {
      // distinct operands case
      ZmodF_poly_t poly2;
      ZmodF_poly_stack_init(poly2, log_length, n, 1);
      Z_split_limbs(poly2, data2, limbs2, coeff_limbs, n);

      ZmodF_poly_convolution(poly1, poly1, poly2);

      ZmodF_poly_stack_clear(poly2);
   }
   
   ZmodF_poly_normalise(poly1);
   
   clear_limbs(res, limbs1 + limbs2);
   
   Z_combine_limbs(res, poly1, coeff_limbs, 2*coeff_limbs+1, limbs1 + limbs2);
   ZmodF_poly_stack_clear(poly1);
   
   return res[limbs1+limbs2-1];
}

void Z_mpn_mul_precomp_init(Z_mpn_precomp_t precomp, mp_limb_t * data1, unsigned long limbs1, unsigned long limbs2)
{
   unsigned long length = 1;
   unsigned long log_length = 0;
   
   unsigned long coeff_limbs = limbs1 + limbs2;
   unsigned long output_bits = coeff_limbs*FLINT_BITS;
   unsigned long n = coeff_limbs;
 
   unsigned long length1 = 1;
   unsigned long length2 = 1;
   
   unsigned log_length2 = 0;
   
   unsigned long twk;
   
   if (coeff_limbs/2 < MUL_TWK_SMALL_CUTOFF) twk = MUL_TWK_SMALL_DEFAULT;
   else if (coeff_limbs/2 > MUL_TWK_LARGE_CUTOFF) twk = MUL_TWK_LARGE_DEFAULT;
   else 
   {  
      for (unsigned long twk_count = 0; twk_count < MUL_TWK_COUNT; twk_count++)
      {
         if ((coeff_limbs/2 < MUL_TWK_VALS[twk_count][0]) || (coeff_limbs/2 < MUL_TWK_VALS[twk_count][1])) continue;
         else 
         {
            twk = MUL_TWK_VALS[twk_count][2];
            break;
         }
      }
   }
   
   while (twk*length < 2*output_bits)
   {
      length<<=1;
      log_length++;
      coeff_limbs = (limbs1+limbs2-1)/length+1;
      while ((limbs1-1)/coeff_limbs+(limbs2-1)/coeff_limbs+2 > length) coeff_limbs++;
      output_bits = (2*coeff_limbs+1)*FLINT_BITS;
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      coeff_limbs = ((output_bits - FLINT_BITS)/FLINT_BITS)/2;
      if ((long) coeff_limbs < 1) coeff_limbs = 1;
      length1 = (limbs1-1)/coeff_limbs+1;
      length2 = (limbs2-1)/coeff_limbs+1;
   }
      
   n = output_bits/FLINT_BITS;
#if DEBUG
   printf("%ld, %ld, %ld, %ld, %ld\n", length1, length2, output_bits, coeff_limbs, log_length);
#endif   
   ZmodF_poly_p poly1;
   poly1 = (ZmodF_poly_p) malloc(sizeof(ZmodF_poly_struct));
   ZmodF_poly_stack_init(poly1, log_length, n, 1);
   Z_split_limbs(poly1, data1, limbs1, coeff_limbs, n);
   
   ZmodF_poly_FFT(poly1, length1 + length2 - 1);
   
   precomp->type = FFT_PRE;
   precomp->length = length1 + length2 - 1;
   precomp->length2 = length2;
   precomp->coeff_limbs = coeff_limbs;
   precomp->limbs1 = limbs1;
   precomp->limbs2 = limbs2;
   precomp->poly = poly1; 
}

void Z_mpn_mul_precomp_clear(Z_mpn_precomp_t precomp)
{
   if (precomp->type == FFT_PRE) ZmodF_poly_stack_clear(precomp->poly);
}

mp_limb_t Z_mpn_mul_precomp(mp_limb_t * res, mp_limb_t * data2, unsigned long limbs2, Z_mpn_precomp_t precomp)
{
   ZmodF_poly_t poly2;
   ZmodF_poly_stack_init(poly2, precomp->poly->depth, precomp->poly->n, 1);
   
   Z_split_limbs(poly2, data2, limbs2, precomp->coeff_limbs, precomp->poly->n);
   
   for (unsigned long i = poly2->length; i < precomp->length2; i++)
   {
      clear_limbs(poly2->coeffs[i], poly2->n+1);
   }
   poly2->length = precomp->length2;
   
   ZmodF_poly_FFT(poly2, precomp->length);
   ZmodF_poly_pointwise_mul(poly2, poly2, precomp->poly);
   ZmodF_poly_IFFT(poly2);
   ZmodF_poly_rescale(poly2);
   
   ZmodF_poly_normalise(poly2);
   
   clear_limbs(res, precomp->limbs1 + precomp->limbs2);
   
   Z_combine_limbs(res, poly2, precomp->coeff_limbs, 2*precomp->coeff_limbs+1, precomp->limbs1 + precomp->limbs2);
   
   ZmodF_poly_stack_clear(poly2);
   
   return res[precomp->limbs1+precomp->limbs2-1];
}

void Z_mul(mpz_t res, mpz_t a, mpz_t b)
{
   unsigned long int limbs;
   if (a->_mp_size + b->_mp_size > 128000/FLINT_BITS) 
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

