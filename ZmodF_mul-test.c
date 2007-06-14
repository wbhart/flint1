/****************************************************************************

ZmodF_mul-test.c: test module for ZmodF_mul module

Copyright (C) 2007, David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "ZmodFpoly.h"
#include "ZmodF_mul.h"

gmp_randstate_t randstate;


int test_ZmodF_mul()
{
   mp_limb_t scratch[2*MAX_N];
   
   for (unsigned long n = 1; n <= 5; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
         for (int outbuf = 0; outbuf <= 2; outbuf++)
            for (int inbuf1 = 0; inbuf1 <= 2; inbuf1++)
               for (int inbuf2 = 0; inbuf2 <= 2; inbuf2++)
               {
                  setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));
            
                  ZmodF_mul(coeffs[outbuf], coeffs[inbuf1], coeffs[inbuf2],
                            scratch, n);

                  if (!check_coeffs(3, n))
                     return 0;
                     
                  mpz_mul(global_mpz, coeffs_mpz_in[inbuf1],
                          coeffs_mpz_in[inbuf2]);
                  mpz_mod(global_mpz, global_mpz, global_p);
                  
                  if (mpz_cmp(coeffs_mpz_out[outbuf], global_mpz))
                     return 0;
               }

   return 1;
}


int test_ZmodF_sqr()
{
   mp_limb_t scratch[2*MAX_N];
   
   for (unsigned long n = 1; n <= 5; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
         for (int outbuf = 0; outbuf <= 1; outbuf++)
            for (int inbuf = 0; inbuf <= 1; inbuf++)
            {
               setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 2));
         
               ZmodF_sqr(coeffs[outbuf], coeffs[inbuf], scratch, n);

               if (!check_coeffs(2, n))
                  return 0;
                  
               mpz_mul(global_mpz, coeffs_mpz_in[inbuf],
                       coeffs_mpz_in[inbuf]);
               mpz_mod(global_mpz, global_mpz, global_p);
               
               if (mpz_cmp(coeffs_mpz_out[outbuf], global_mpz))
                  return 0;
            }

   return 1;
}


int test_ZmodF_mul()
{
}


/****************************************************************************

   Main test functions

****************************************************************************/


void ZmodF_mul_test_all()
{
   int success, all_success = 1;

   RUN_TEST(ZmodF_mul);
   RUN_TEST(ZmodF_sqr);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(randstate);
   
   ZmodF_mul_test_all();

   return 0;
}



// end of file ****************************************************************
