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

Z_mpn-test.c: test module for Z_mpn module

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "mpn_extras.h"
#include "mpz_poly.h"
#include "test-support.h"
#include "F_mpn_mul-tuning.h"

#define VARY_BITS 1
#define SIGNS 1

#define DEBUG 0    // prints debug information
#define DEBUG2 1 


/****************************************************************************

   Test code for Conversion Routines
   
****************************************************************************/


unsigned long randint(unsigned long randsup) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(mpz_poly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_zero(pol);

   for (unsigned long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, randstate, bits);
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
       
   }
   
   mpz_clear(temp);
} 

void randpoly_unsigned(mpz_poly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);

   mpz_poly_zero(pol);
   
   for (unsigned long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, randstate, bits);
       }
       mpz_poly_set_coeff(pol, i, temp);
       
   }
   
   mpz_clear(temp);
} 

int test_F_mpn_splitcombine_bits()
{
    mp_limb_t * int1, * int2;
    ZmodF_poly_t poly;
    int result = 1;
    
    for (unsigned long count = 0; (count < 30000) && (result == 1); count++)
    {
        unsigned long limbs = randint(300)+1;
        unsigned long bits = randint(500)+1;
        unsigned long coeff_limbs = randint(100) + (bits-1)/FLINT_BITS + 1;
        unsigned long length = (FLINT_BITS*limbs - 1)/bits + 1;
        unsigned long log_length = 0;
        while ((1L << log_length) < length) log_length++;
        
#if DEBUG
        printf("limbs = %ld, bits = %ld, coeff_limbs = %ld\n", limbs, bits, coeff_limbs);
#endif
        
        int1 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs);
        int2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs);
        ZmodF_poly_init(poly, log_length, coeff_limbs, 0);
        
        mpn_random2(int1, limbs);
        F_mpn_FFT_split_bits(poly, int1, limbs, bits, coeff_limbs);
        F_mpn_clear(int2, limbs);
        F_mpn_FFT_combine_bits(int2, poly, bits, coeff_limbs, limbs);

#if DEBUG
        F_mpn_printx(int1, limbs); printf("\n\n");
        for (unsigned long i = 0; i < length; i++) { F_mpn_printx(poly->coeffs[i], coeff_limbs); printf("\n");}
        printf("\n");
        F_mpn_printx(int2, limbs); printf("\n\n");
#endif

        for (unsigned long j = 0; j < limbs; j++)
        {
           if (int1[j] != int2[j]) result = 0;
        }
        
        ZmodF_poly_clear(poly);
        free(int2);
        free(int1);
    }
    
    return result;
}

int test_F_mpn_mul_precache()
{
   mp_limb_t * int1, * int2, * product, * product2;
   F_mpn_precache_t precache;
   mp_limb_t msl;
   int result = 1;
   
   for (unsigned long count = 0; (count < 30) && (result == 1); count++)
   {
      unsigned long limbs2 = randint(2*FLINT_FFT_LIMBS_CROSSOVER)+1;
      unsigned long limbs1 = randint(2*FLINT_FFT_LIMBS_CROSSOVER)+1;
   
      int1 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs1);

      mpn_random2(int1, limbs1);
   
      F_mpn_mul_precache_init(precache, int1, limbs1, limbs2);   
           
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {    
#if DEBUG
         printf("%ld, %ld\n",limbs1, limbs2);
#endif

         unsigned long limbs3 = randint(limbs2)+1;
         int2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs3);
         product = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         product2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         
         F_mpn_clear(int2, limbs3);
         mpn_random2(int2, limbs3);
      
         F_mpn_mul_precache(product, int2, limbs3, precache);
         
         if (limbs1 > limbs3) msl = mpn_mul(product2, int1, limbs1, int2, limbs3);
         else msl = mpn_mul(product2, int2, limbs3, int1, limbs1);
      
         for (unsigned long j = 0; j < limbs1+limbs3 - (msl == 0); j++)
         {
            if (product[j] != product2[j]) result = 0;
         }
      
         free(product2);
         free(product);
         free(int2);
      }
   
      F_mpn_mul_precache_clear(precache);
      
      free(int1);
   }   
   
   return result;
}

int test_F_mpn_mul_precache_trunc()
{
   mp_limb_t * int1, * int2, * product, * product2;
   F_mpn_precache_t precache;
   mp_limb_t msl;
   int result = 1;
   
   for (unsigned long count = 0; (count < 30) && (result == 1); count++)
   {
      unsigned long limbs2 = randint(2*FLINT_FFT_LIMBS_CROSSOVER)+1;
      unsigned long limbs1 = randint(2*FLINT_FFT_LIMBS_CROSSOVER)+1;
   
      int1 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs1);

      mpn_random2(int1, limbs1);
   
      F_mpn_mul_precache_init(precache, int1, limbs1, limbs2);   
           
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {    
         unsigned long limbs3 = randint(limbs2)+1;
         unsigned long trunc = randint(2*(limbs1+limbs3));
#if DEBUG
         printf("limbs1 = %ld, limbs3 = %ld, trunc = %ld\n", limbs1, limbs3, trunc);
#endif

         int2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs3);
         product = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         product2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         
         F_mpn_clear(int2, limbs3);
         mpn_random2(int2, limbs3);
      
         if (limbs1 > limbs3) F_mpn_mul_trunc(product2, int1, limbs1, int2, limbs3, trunc);
         else F_mpn_mul_trunc(product2, int2, limbs3, int1, limbs1, trunc);
         F_mpn_mul_precache_trunc(product, int2, limbs3, precache, trunc);
      
         for (unsigned long j = 0; j < FLINT_MIN(trunc, limbs1+limbs3); j++)
         {
            if (product[j] != product2[j]) 
            {
               printf("Failure at %ld\n", j);
               result = 0;
            }
         }
      
         free(product2);
         free(product);
         free(int2);
      }
   
      F_mpn_mul_precache_clear(precache);
      
      free(int1);
   }   
   
   return result;
}

int test_F_mpn_mul()
{
   mp_limb_t * int1, * int2, * product, * product2;
   mp_limb_t msl, msl2;
   int result = 1;
   
   for (unsigned long count = 0; (count < 30) && (result == 1); count++)
   {
      unsigned long limbs2 = randint(2*FLINT_FFT_LIMBS_CROSSOVER)+1;
      unsigned long limbs1 = limbs2 + randint(1000);
   
      int1 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs1);

      mpn_random2(int1, limbs1);
        
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {    
#if DEBUG
         printf("%ld, %ld\n",limbs1, limbs2);
#endif

         int2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs2);
         product = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         product2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         
         F_mpn_clear(int2, limbs2);
         mpn_random2(int2, randint(limbs2-1)+1);
      
         msl = F_mpn_mul(product, int1, limbs1, int2, limbs2);
         
         msl2 = mpn_mul(product2, int1, limbs1, int2, limbs2);
      
         for (unsigned long j = 0; j < limbs1+limbs2 - (msl == 0); j++)
         {
            if (product[j] != product2[j]) result = 0;
         }
         
         result &= (msl == msl2);
         
         free(product2);
         free(product);
         free(int2);
      }
   
      free(int1);
   }   
   
   return result;
}

int test_F_mpn_mul_trunc()
{
   mp_limb_t * int1, * int2, * product, * product2;
   mp_limb_t msl;
   int result = 1;
   
   for (unsigned long count = 0; (count < 30) && (result == 1); count++)
   {
      unsigned long limbs2 = randint(2*FLINT_FFT_LIMBS_CROSSOVER)+1;
      unsigned long limbs1 = limbs2 + randint(1000);
      
      int1 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs1);

      mpn_random2(int1, limbs1);
        
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {    
#if DEBUG
         printf("%ld, %ld\n",limbs1, limbs2);
#endif

         unsigned long trunc = randint(limbs1 + limbs2 - 1)+1;
         int2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs2);
         product = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         product2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         
         F_mpn_clear(int2, limbs2);
         mpn_random2(int2, limbs2);
      
         F_mpn_mul_trunc(product, int1, limbs1, int2, limbs2, trunc);
         
         mpn_mul(product2, int1, limbs1, int2, limbs2);
      
         for (unsigned long j = 0; j < trunc; j++)
         {
            if (product[j] != product2[j]) result = 0;
         }
         
         free(product2);
         free(product);
         free(int2);
      }
   
      free(int1);
   }   
   
   return result;
}

/****************************************************************************

   Main test functions

****************************************************************************/

void F_mpn_test_all()
{
   int success, all_success = 1;

   RUN_TEST(F_mpn_splitcombine_bits);
   RUN_TEST(F_mpn_mul);
   RUN_TEST(F_mpn_mul_trunc);
   RUN_TEST(F_mpn_mul_precache);
   RUN_TEST(F_mpn_mul_precache_trunc);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpn_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}



// end of file ****************************************************************
