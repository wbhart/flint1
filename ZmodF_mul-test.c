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


#define DEBUG 0    // prints debug information


/*   // disabled temporarily
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
*/


int test__ZmodF_mul_negacyclic_split()
{
   int success = 1;

   mpz_t x, y, z;
   mpz_init(x);
   mpz_init(y);
   mpz_init(z);
   mp_limb_t buf[300];

   for (unsigned long n = 1; n < 200 && success; n++)
   {
      for (unsigned long depth = 0;
           ((n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0) && success; depth++)
      {
         ZmodFpoly_t poly;
         ZmodFpoly_init(poly, depth, n, 1);

#if DEBUG
         printf("n = %d, depth = %d\n", n, depth);
#endif
         
         for (unsigned long trial = 0; trial < 30; trial++)
         {
            mpz_rrandomb(x, randstate, n*FLINT_BITS_PER_LIMB);
            memset(buf, 0, n+1);
            mpz_export(buf, NULL, -1, sizeof(mp_limb_t), 0, 0, x);
            
            _ZmodF_mul_negacyclic_split(poly, buf, n);

            for (unsigned long i = 0; i < (1 << depth); i++)
            {
               mpz_tdiv_r_2exp(y, x, (n*FLINT_BITS_PER_LIMB) >> depth);
               mpz_tdiv_q_2exp(x, x, (n*FLINT_BITS_PER_LIMB) >> depth);
               mpz_import(z, poly->n+1, -1, sizeof(mp_limb_t), 0, 0, poly->coeffs[i]);
               if (mpz_cmp(z, y))
                  success = 0;
            }
         }
         
         ZmodFpoly_clear(poly);
      }
   }
   
   mpz_clear(x);
   mpz_clear(y);
   mpz_clear(z);
   return success;
}


int test__ZmodF_mul_negacyclic_combine()
{
   return 0;
}


int test__ZmodF_mul_threeway_reduce1()
{
   return 0;
}


int test__ZmodF_mul_threeway_reduce2()
{
   return 0;
}


int test__ZmodF_mul_threeway_crt()
{
   return 0;
}


int test_ZmodF_mul_info_mul_plain()
{
   return 0;
}


int test_ZmodF_mul_info_mul_threeway()
{
   return 0;
}


int test_ZmodF_mul_info_mul_negacyclic()
{
   return 0;
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


void ZmodF_mul_test_all()
{
   int success, all_success = 1;

/*
   RUN_TEST(ZmodF_mul);
   RUN_TEST(ZmodF_sqr);
*/
   RUN_TEST(_ZmodF_mul_negacyclic_split);
   RUN_TEST(_ZmodF_mul_negacyclic_combine);
   RUN_TEST(_ZmodF_mul_threeway_reduce1);
   RUN_TEST(_ZmodF_mul_threeway_reduce2);
   RUN_TEST(_ZmodF_mul_threeway_crt);
   RUN_TEST(ZmodF_mul_info_mul_plain);
   RUN_TEST(ZmodF_mul_info_mul_threeway);
   RUN_TEST(ZmodF_mul_info_mul_negacyclic);

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
