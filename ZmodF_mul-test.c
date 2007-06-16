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


/*
Prints the ZmodF_t to stdout in hex, each limb in a separate block,
most significant limb (i.e. the overflow limb) first.
*/
void ZmodF_print(ZmodF_t x, unsigned long n)
{
   for (long i = n; i >= 0; i--)
#if FLINT_BITS == 64
      printf("%016lx ", x[i]);
#else
      printf("%08lx ", x[i]);
#endif
}



/*
Prints each coefficient of the polynomial on a separate line.
*/
void ZmodFpoly_print(ZmodFpoly_t x)
{
   for (unsigned long k = 0; k < (1UL << x->depth); k++)
   {
      ZmodF_print(x->coeffs[k], x->n);
      printf("\n");
   }
}



unsigned long random_ulong(unsigned long max)
{
   return gmp_urandomm_ui(randstate, max);
}


// uses mpz_rrandomb to write random bits of length k limbs to x
void random_limbs(mp_limb_t* x, unsigned long k)
{
   mpz_t y;
   mpz_init(y);
   
   mpz_rrandomb(y, randstate, k*FLINT_BITS);
   if (random_ulong(2))
   {
      mpz_t power;
      mpz_init(power);
      
      mpz_set_ui(power, 1);
      mpz_mul_2exp(power, power, k*FLINT_BITS);
      mpz_sub_ui(power, power, 1);
      mpz_sub(y, power, y);
      
      mpz_clear(power);
   }

   memset(x, 0, k * sizeof(mp_limb_t));
   mpz_export(x, NULL, -1, sizeof(mp_limb_t), 0, 0, y);
   
   mpz_clear(y);
}


/*
int test_ZmodF_mul()
{
   mp_limb_t scratch[2*MAX_N];
   
   for (unsigned long n = 1; n <= 5; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
         for (int outbuf = 0; outbuf <= 2; outbuf++)
            for (int inbuf1 = 0; inbuf1 <= 2; inbuf1++)
               for (int inbuf2 = 0; inbuf2 <= 2; inbuf2++)
               {
                  setup_coeffs(3, n, random_ulong(FLINT_BITS - 2));
            
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
               setup_coeffs(2, n, random_ulong(FLINT_BITS - 2));
         
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
           ((n*FLINT_BITS) % (1 << depth) == 0) && success; depth++)
      {
         unsigned long bits = (n*FLINT_BITS) >> depth;
         unsigned long m = (bits-1)/FLINT_BITS + 1;
      
         ZmodFpoly_t poly;
         ZmodFpoly_init(poly, depth, m, 1);

#if DEBUG
         printf("n = %d, depth = %d, m = %d\n", n, depth, m);
#endif
         
         for (unsigned long trial = 0; trial < 120; trial++)
         {
            mpz_rrandomb(x, randstate, n*FLINT_BITS);
            memset(buf, 0, (n+1) * sizeof(mp_limb_t));
            mpz_export(buf, NULL, -1, sizeof(mp_limb_t), 0, 0, x);
            
            _ZmodF_mul_negacyclic_split(poly, buf, n);

            for (unsigned long i = 0; i < (1 << depth); i++)
            {
               mpz_tdiv_r_2exp(y, x, bits);
               mpz_tdiv_q_2exp(x, x, bits);
               mpz_import(z, m+1, -1, sizeof(mp_limb_t), 0, 0, poly->coeffs[i]);
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
   int success = 1;
   
   mpz_t x, y, p, q, half_q, total;
   mpz_init(x);
   mpz_init(y);
   mpz_init(q);
   mpz_init(p);
   mpz_init(half_q);
   mpz_init(total);

   mp_limb_t buf[300];

   for (unsigned long n = 1; n < 200 && success; n++)
   {
      for (unsigned long depth = 0;
           ((n*FLINT_BITS) % (1 << depth) == 0) && success; depth++)
      {
         unsigned long bits = (n*FLINT_BITS) >> depth;
         unsigned long m = (bits-1)/FLINT_BITS + 1;

         ZmodFpoly_t poly;
         ZmodFpoly_init(poly, depth, m, 1);

#if DEBUG
         printf("n = %d, depth = %d, m = %d\n", n, depth, m);
#endif

         // p := B^n + 1
         mpz_set_ui(p, 1);
         mpz_mul_2exp(p, p, n*FLINT_BITS);
         mpz_add_ui(p, p, 1);
         
         // q := B^m + 1
         // half_q := B^m / 2  = (q-1)/2
         mpz_set_ui(half_q, 1);
         mpz_mul_2exp(half_q, half_q, m*FLINT_BITS - 1);
         mpz_add(q, half_q, half_q);
         mpz_add_ui(q, q, 1);
         
         for (unsigned long trial = 0; trial < 100; trial++)
         {
            mpz_set_ui(total, 0);
         
            for (long i = (1 << depth) - 1; i >= 0; i--)
            {
               memset(poly->coeffs[i], 0, (m+1) * sizeof(mp_limb_t));
                  
               // select a random x in [0, q)
               if (random_ulong(5) == 0)
               {
                  // try the special case -1 mod q every now and then
                  mpz_sub_ui(x, q, 1);
               }
               else
               {
                  mpz_rrandomb(x, randstate, m*FLINT_BITS);
                  if (random_ulong(2))
                  {
                     // flip leading bit
                     mpz_sub(x, q, x);
                     mpz_sub_ui(x, x, 2);
                  }
               }
               
               // store in poly coefficient
               mpz_export(poly->coeffs[i], NULL, -1, sizeof(mp_limb_t), 0, 0, x);
               
               // convert to range [-(q+1)/2, (q-3)/2]
               mpz_add(x, x, half_q);
               mpz_add_ui(x, x, 1);
               mpz_mod(x, x, q);
               mpz_sub(x, x, half_q);
               mpz_sub_ui(x, x, 1);
               
               // add to running total
               mpz_mul_2exp(total, total, bits);
               mpz_add(total, total, x);
            }

            // reduce total mod p and compare to result of target function
            mpz_mod(total, total, p);
            
            _ZmodF_mul_negacyclic_combine(buf, poly, n);
            ZmodF_normalise(buf, n);
            mpz_import(y, n+1, -1, sizeof(mp_limb_t), 0, 0, buf);

            if (mpz_cmp(y, total))
            {
               success = 0;
               abort();
            }
         }
         
         ZmodFpoly_clear(poly);
      }
   }

   mpz_clear(x);
   mpz_clear(y);
   mpz_clear(q);
   mpz_clear(p);
   mpz_clear(half_q);
   mpz_clear(total);

   return success;
}


int test__ZmodF_mul_threeway_reduce()
{
   int success = 1;

   mp_limb_t in[2000];
   mp_limb_t out1[2000];
   mp_limb_t out2[2000];
   mp_limb_t test[2000];

   mpz_t x, y, power, power2, mod1, mod2;
   mpz_init(x);
   mpz_init(y);
   mpz_init(power);
   mpz_init(power2);
   mpz_init(mod1);
   mpz_init(mod2);

   for (unsigned long n = 3; n < 1000 && success; n += 3)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif
      
      // power = B^n
      mpz_set_ui(power, 1);
      mpz_mul_2exp(power, power, n*FLINT_BITS);

      // power2 = B^(2n/3)
      mpz_set_ui(power2, 1);
      mpz_mul_2exp(power2, power2, 2*n/3*FLINT_BITS);

      // mod1 = B^(n/3) + 1
      mpz_set_ui(mod1, 1);
      mpz_mul_2exp(mod1, mod1, n/3*FLINT_BITS);
      mpz_add_ui(mod1, mod1, 1);

      // mod2 = B^(2n/3) - B^(n/3) + 1
      mpz_set(mod2, mod1);
      mpz_mul_2exp(mod2, mod2, n/3*FLINT_BITS);
      mpz_sub(mod2, mod2, mod1);
      mpz_sub(mod2, mod2, mod1);
      mpz_add_ui(mod2, mod2, 3);

      for (unsigned long trial = 0; trial < 250 && success; trial++)
      {
         ZmodF_zero(in, n);
         mpz_rrandomb(x, randstate, n*FLINT_BITS);
         if (random_ulong(2))
         {
            // flip bits with probability 1/2
            mpz_sub(x, power, x);
            mpz_sub_ui(x, x, 1);
         }
         mpz_export(in, NULL, -1, sizeof(mp_limb_t), 0, 0, x);

         _ZmodF_mul_threeway_reduce1(out1, in, n/3);
         mpz_mod(y, x, mod1);
         memset(test, 0, (n/3 + 1) * sizeof(mp_limb_t));
         mpz_export(test, NULL, -1, sizeof(mp_limb_t), 0, 0, y);
         if (memcmp(test, out1, (n/3 + 1) * sizeof(mp_limb_t)))
             success = 0;

         _ZmodF_mul_threeway_reduce2(out2, in, n/3);
         mpz_mod(y, x, mod2);
         memset(test, 0, (2*n/3) * sizeof(mp_limb_t));
         mpz_export(test, NULL, -1, sizeof(mp_limb_t), 0, 0, y);
         if (memcmp(test, out2, 2*n/3 * sizeof(mp_limb_t)))
         {
            // didn't work... check if the "other answer" is correct
            mpz_add(y, y, mod2);
            if (mpz_cmp(y, power2) >= 0)
               success = 0;
            else
            {
               memset(test, 0, (2*n/3) * sizeof(mp_limb_t));
               mpz_export(test, NULL, -1, sizeof(mp_limb_t), 0, 0, y);
               if (memcmp(test, out2, 2*n/3 * sizeof(mp_limb_t)))
                  return 0;
            }
         }
      }
   }

   mpz_clear(mod2);
   mpz_clear(mod1);
   mpz_clear(power2);
   mpz_clear(power);
   mpz_clear(y);
   mpz_clear(x);

   return success;
}


int test__ZmodF_mul_threeway_crt()
{
   int success = 1;
   
   mp_limb_t in1[2000];
   mp_limb_t in2[2000];
   mp_limb_t out[2000];
   mpz_t x1, x2, y1, y2, z, mod1, mod2;
   
   mpz_init(x1);
   mpz_init(x2);
   mpz_init(y1);
   mpz_init(y2);
   mpz_init(z);
   mpz_init(mod1);
   mpz_init(mod2);

   for (unsigned long n = 3; n < 1000 && success; n += 3)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif

      // mod1 = B^(n/3) + 1
      mpz_set_ui(mod1, 1);
      mpz_mul_2exp(mod1, mod1, n/3*FLINT_BITS);
      mpz_add_ui(mod1, mod1, 1);

      // mod2 = B^(2n/3) - B^(n/3) + 1
      mpz_set(mod2, mod1);
      mpz_mul_2exp(mod2, mod2, n/3*FLINT_BITS);
      mpz_sub(mod2, mod2, mod1);
      mpz_sub(mod2, mod2, mod1);
      mpz_add_ui(mod2, mod2, 3);

      for (unsigned long trial = 0; trial < 100 && success; trial++)
      {
         if (random_ulong(5) == 0)
         {
            // every now and then generate -1 mod B^(n/3) + 1
            ZmodF_zero(in1, n/3);
            in1[n/3] = 1;
         }
         else
         {
            random_limbs(in1, n/3);
            in1[n/3] = 0;
         }
         
         random_limbs(in2, 2*n/3);
         
         mpz_import(x1, n/3 + 1, -1, sizeof(mp_limb_t), 0, 0, in1);
         mpz_import(x2, 2*n/3, -1, sizeof(mp_limb_t), 0, 0, in2);
         mpz_mod(x2, x2, mod2);
         
         _ZmodF_mul_threeway_crt(out, in1, in2, n/3);
         ZmodF_normalise(out, n);
         mpz_import(z, n+1, -1, sizeof(mp_limb_t), 0, 0, out);

         mpz_mod(y1, z, mod1);
         mpz_mod(y2, z, mod2);
         if (mpz_cmp(y1, x1) || mpz_cmp(y2, x2))
            success = 0;
      }
   }

   mpz_clear(mod1);
   mpz_clear(mod2);
   mpz_clear(x1);
   mpz_clear(x2);
   mpz_clear(y1);
   mpz_clear(y2);
   mpz_clear(z);

   return success;
}


int test_ZmodF_mul_info_mul_plain()
{
   return 0;
}


int test_ZmodF_mul_info_mul_threeway()
{
   int success = 1;

   mp_limb_t in1[2000];
   mp_limb_t in2[2000];
   mp_limb_t out_plain[2000];
   mp_limb_t out_threeway[2000];

   mpz_t x;
   mpz_init(x);

   for (unsigned long n = 3; n < 1000 && success; n += 3)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif

      ZmodF_mul_info_t info_plain, info_threeway;
      ZmodF_mul_info_init_threeway(info_threeway, n);
      ZmodF_mul_info_init_plain(info_plain, n);

      for (unsigned long trial = 0; trial < 250 && success; trial++)
      {
         ZmodF_zero(in1, n);
         mpz_rrandomb(x, randstate, n*FLINT_BITS + FLINT_BITS/2);
         mpz_export(in1, NULL, -1, sizeof(mp_limb_t), 0, 0, x);
         
         ZmodF_zero(in2, n);
         mpz_rrandomb(x, randstate, n*FLINT_BITS + FLINT_BITS/2);
         mpz_export(in2, NULL, -1, sizeof(mp_limb_t), 0, 0, x);
         
         ZmodF_mul_info_mul(info_plain, out_plain, in1, in2);
         ZmodF_mul_info_mul(info_threeway, out_threeway, in1, in2);
         
         ZmodF_normalise(out_plain, n);
         ZmodF_normalise(out_threeway, n);
         
         if (memcmp(out_plain, out_threeway, (n+1) * sizeof(mp_limb_t)))
            success = 0;
      }
      
      ZmodF_mul_info_clear(info_plain);
      ZmodF_mul_info_clear(info_threeway);
   }

   mpz_clear(x);

   return success;
}


int test_ZmodF_mul_info_mul_negacyclic()
{
   int success = 1;

   mp_limb_t in1[2000];
   mp_limb_t in2[2000];
   mp_limb_t out_plain[2000];
   mp_limb_t out_negacyclic[2000];

   mpz_t x;
   mpz_init(x);

   for (unsigned long n = 1; n < 1000 && success; n++)
   {
      for (unsigned long depth = 1;
           (n*FLINT_BITS) % (1 << depth) == 0
           && (depth <= FLINT_LG_BITS_PER_LIMB + 4)
           && success; depth++)
      {

#if DEBUG
         printf("n = %d, depth = %d\n", n, depth);
#endif

         ZmodF_mul_info_t info_plain, info_negacyclic;
         ZmodF_mul_info_init_negacyclic(info_negacyclic, n, depth);
         ZmodF_mul_info_init_plain(info_plain, n);

         for (unsigned long trial = 0; trial < 10 && success; trial++)
         {
            ZmodF_zero(in1, n);
            mpz_rrandomb(x, randstate, n*FLINT_BITS + FLINT_BITS/2);
            mpz_export(in1, NULL, -1, sizeof(mp_limb_t), 0, 0, x);

            ZmodF_zero(in2, n);
            mpz_rrandomb(x, randstate, n*FLINT_BITS + FLINT_BITS/2);
            mpz_export(in2, NULL, -1, sizeof(mp_limb_t), 0, 0, x);

            ZmodF_mul_info_mul(info_plain, out_plain, in1, in2);
            ZmodF_mul_info_mul(info_negacyclic, out_negacyclic, in1, in2);

            ZmodF_normalise(out_plain, n);
            ZmodF_normalise(out_negacyclic, n);

            if (memcmp(out_plain, out_negacyclic, (n+1) * sizeof(mp_limb_t)))
               success = 0;
         }

         ZmodF_mul_info_clear(info_plain);
         ZmodF_mul_info_clear(info_negacyclic);
      }
   }

   mpz_clear(x);

   return success;
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

//   RUN_TEST(ZmodF_mul);
//   RUN_TEST(ZmodF_sqr);

//   RUN_TEST(_ZmodF_mul_negacyclic_split);
//   RUN_TEST(_ZmodF_mul_negacyclic_combine);
//   RUN_TEST(_ZmodF_mul_threeway_reduce);
   RUN_TEST(_ZmodF_mul_threeway_crt);
   RUN_TEST(ZmodF_mul_info_mul_plain);
//   RUN_TEST(ZmodF_mul_info_mul_threeway);
//   RUN_TEST(ZmodF_mul_info_mul_negacyclic);

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
