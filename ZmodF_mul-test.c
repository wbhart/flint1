/****************************************************************************

ZmodF_mul-test.c: test module for ZmodF_mul module

Copyright (C) 2007, David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "ZmodF_mul.h"
#include "test-support.h"



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
            random_limbs(buf, n);
            buf[n] = 0;
            mpn_to_mpz(x, buf, n);
            
            _ZmodF_mul_negacyclic_split(poly, buf, n);

            for (unsigned long i = 0; i < (1 << depth); i++)
            {
               mpz_tdiv_r_2exp(y, x, bits);
               mpz_tdiv_q_2exp(x, x, bits);
               mpn_to_mpz(z, poly->coeffs[i], m+1);
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
               // select a random x in [0, q)
               if (random_ulong(5) == 0)
               {
                  // try the special case -1 mod q every now and then
                  memset(poly->coeffs[i], 0, (m+1) * sizeof(mp_limb_t));
                  poly->coeffs[i][m] = 1;
               }
               else
               {
                  random_limbs(poly->coeffs[i], m);
                  poly->coeffs[i][m] = 0;
               }
               
               // store in poly coefficient
               mpn_to_mpz(x, poly->coeffs[i], m+1);
               
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
            mpn_to_mpz(y, buf, n+1);

            if (mpz_cmp(y, total))
               success = 0;
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

   for (unsigned long n = 3; n < 300 && success; n += 3)
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
         random_limbs(in, n);
         in[n] = 0;
         mpn_to_mpz(x, in, n+1);

         _ZmodF_mul_threeway_reduce1(out1, in, n/3);
         mpz_mod(y, x, mod1);
         mpz_to_mpn(test, n/3 + 1, y);
         if (mpn_cmp(test, out1, n/3 + 1))
             success = 0;

         _ZmodF_mul_threeway_reduce2(out2, in, n/3);
         mpz_mod(y, x, mod2);
         mpz_to_mpn(test, 2*n/3, y);
         if (mpn_cmp(test, out2, 2*n/3))
         {
            // didn't work... check if the "other answer" is correct
            mpz_add(y, y, mod2);
            if (mpz_cmp(y, power2) >= 0)
               success = 0;
            else
            {
               mpz_to_mpn(test, 2*n/3, y);
               if (mpn_cmp(test, out2, 2*n/3))
                  success = 0;
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

   for (unsigned long n = 3; n < 300 && success; n += 3)
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

         mpn_to_mpz(x1, in1, n/3 + 1);
         mpn_to_mpz(x2, in2, 2*n/3);
         mpz_mod(x2, x2, mod2);
         
         _ZmodF_mul_threeway_crt(out, in1, in2, n/3);
         ZmodF_normalise(out, n);
         mpn_to_mpz(z, out, n+1);

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
   int success = 1;

   mp_limb_t in1[2000];
   mp_limb_t in2[2000];
   mp_limb_t out[2000];

   mpz_t x1, x2, y, z, p;
   mpz_init(x1);
   mpz_init(x2);
   mpz_init(y);
   mpz_init(z);
   mpz_init(p);

   for (unsigned long n = 1; n < 100 && success; n++)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif

      // p = B^n + 1
      mpz_set_ui(p, 1);
      mpz_mul_2exp(p, p, n*FLINT_BITS);
      mpz_add_ui(p, p, 1);

      ZmodF_mul_info_t info;
      ZmodF_mul_info_init_plain(info, n);

      for (unsigned long trial = 0; trial < 1000 && success; trial++)
      {
         if (random_ulong(4) == 0)
         {
            // put in -1 mod p every now and then
            ZmodF_zero(in1, n);
            in1[n] = 1;
         }
         else
         {
            random_limbs(in1, n);
            in1[n] = 0;
         }

         if (random_ulong(4) == 0)
         {
            // put in -1 mod p every now and then
            ZmodF_zero(in2, n);
            in2[n] = 1;
         }
         else
         {
            random_limbs(in2, n);
            in2[n] = 0;
         }

         // test multiplication
         
         mpn_to_mpz(x1, in1, n+1);
         mpn_to_mpz(x2, in2, n+1);
         mpz_mul(z, x1, x2);
         mpz_mod(z, z, p);
         
         ZmodF_mul_info_mul(info, out, in1, in2);
         ZmodF_normalise(out, n);
         mpn_to_mpz(y, out, n+1);
         
         if (mpz_cmp(y, z))
            success = 0;

         // test squaring

         mpz_mul(z, x1, x1);
         mpz_mod(z, z, p);
         ZmodF_mul_info_sqr(info, out, in1);
         ZmodF_normalise(out, n);
         mpn_to_mpz(y, out, n+1);
         
         if (mpz_cmp(y, z))
            success = 0;

      }
      
      ZmodF_mul_info_clear(info);
   }

   mpz_clear(x1);
   mpz_clear(x2);
   mpz_clear(y);
   mpz_clear(z);
   mpz_clear(p);

   return success;
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

   for (unsigned long n = 3; n < 300 && success; n += 3)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif

      ZmodF_mul_info_t info_plain, info_threeway;
      ZmodF_mul_info_init_threeway(info_threeway, n);
      ZmodF_mul_info_init_plain(info_plain, n);

      for (unsigned long trial = 0; trial < 250 && success; trial++)
      {
         if (random_ulong(4) == 0)
         {
            // put in -1 mod p every now and then
            ZmodF_zero(in1, n);
            in1[n] = 1;
         }
         else
         {
            random_limbs(in1, n);
            in1[n] = 0;
         }

         if (random_ulong(4) == 0)
         {
            // put in -1 mod p every now and then
            ZmodF_zero(in2, n);
            in2[n] = 1;
         }
         else
         {
            random_limbs(in2, n);
            in2[n] = 0;
         }

         // test multiplication
         
         ZmodF_mul_info_mul(info_plain, out_plain, in1, in2);
         ZmodF_mul_info_mul(info_threeway, out_threeway, in1, in2);
         
         ZmodF_normalise(out_plain, n);
         ZmodF_normalise(out_threeway, n);
         
         if (mpn_cmp(out_plain, out_threeway, n+1))
            success = 0;

         // test squaring

         ZmodF_mul_info_sqr(info_plain, out_plain, in1);
         ZmodF_mul_info_sqr(info_threeway, out_threeway, in1);

         ZmodF_normalise(out_plain, n);
         ZmodF_normalise(out_threeway, n);
         
         if (mpn_cmp(out_plain, out_threeway, n+1))
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
            if (random_ulong(4) == 0)
            {
               // put in -1 mod p every now and then
               ZmodF_zero(in1, n);
               in1[n] = 1;
            }
            else
            {
               random_limbs(in1, n);
               in1[n] = 0;
            }

            if (random_ulong(4) == 0)
            {
               // put in -1 mod p every now and then
               ZmodF_zero(in2, n);
               in2[n] = 1;
            }
            else
            {
               random_limbs(in2, n);
               in2[n] = 0;
            }

            // test multiplication

            ZmodF_mul_info_mul(info_plain, out_plain, in1, in2);
            ZmodF_mul_info_mul(info_negacyclic, out_negacyclic, in1, in2);

            ZmodF_normalise(out_plain, n);
            ZmodF_normalise(out_negacyclic, n);

            if (mpn_cmp(out_plain, out_negacyclic, n+1))
               success = 0;

            // test squaring
            
            ZmodF_mul_info_sqr(info_plain, out_plain, in1);
            ZmodF_mul_info_sqr(info_negacyclic, out_negacyclic, in1);

            ZmodF_normalise(out_plain, n);
            ZmodF_normalise(out_negacyclic, n);

            if (mpn_cmp(out_plain, out_negacyclic, n+1))
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

   RUN_TEST(_ZmodF_mul_negacyclic_split);
   RUN_TEST(_ZmodF_mul_negacyclic_combine);
   RUN_TEST(_ZmodF_mul_threeway_reduce);
   RUN_TEST(_ZmodF_mul_threeway_crt);
   RUN_TEST(ZmodF_mul_info_mul_plain);
   RUN_TEST(ZmodF_mul_info_mul_threeway);
   RUN_TEST(ZmodF_mul_info_mul_negacyclic);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   test_support_init();
   
   ZmodF_mul_test_all();

   test_support_cleanup();

   return 0;
}



// end of file ****************************************************************
