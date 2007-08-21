/****************************************************************************

Z_mpn-test.c: test module for Z_mpn module

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "Z_mpn.h"
#include "mpz_poly.h"
#include "test-support.h"

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

int test_Z_mpn_mul_precomp()
{
   mp_limb_t * int1, * int2, * product, * product2;
   Z_mpn_precomp_t precomp;
   mp_limb_t msl;
   int result = 1;
   
   for (unsigned long count = 0; (count < 100) && (result == 1); count++)
   {
      unsigned long limbs2 = randint(100) + 100;
      unsigned long limbs1 = limbs2 + randint(100);
   
      int1 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs1);

      mpn_random2(int1, limbs1);
   
      Z_mpn_mul_precomp_init(precomp, int1, limbs1, limbs2);   
           
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {    
#if DEBUG
         printf("%ld, %ld\n",limbs1, limbs2);
#endif

         int2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*limbs2);
         product = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         product2 = (mp_limb_t *) malloc(sizeof(mp_limb_t)*(limbs1+limbs2));
         
         clear_limbs(int2, limbs2);
         mpn_random2(int2, randint(limbs2-1)+1);
      
         Z_mpn_mul_precomp(product, int2, limbs2, precomp);
         
         msl = mpn_mul(product2, int1, limbs1, int2, limbs2);
      
         for (unsigned long j = 0; j < limbs1+limbs2 - (msl == 0); j++)
         {
            if (product[j] != product2[j]) result = 0;
         }
      
         free(product2);
         free(product);
         free(int2);
      }
   
      Z_mpn_mul_precomp_clear(precomp);
      
      free(int1);
   }   
   
   return result;
}


/****************************************************************************

   Main test functions

****************************************************************************/

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void Z_mpn_test_all()
{
   int success, all_success = 1;

   RUN_TEST(Z_mpn_mul_precomp);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   test_support_init();
   Z_mpn_test_all();
   test_support_cleanup();

   return 0;
}



// end of file ****************************************************************
