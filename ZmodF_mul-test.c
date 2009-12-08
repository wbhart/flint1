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

ZmodF_mul-test.c: test module for ZmodF_mul module

Copyright (C) 2007, David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "ZmodF_mul.h"
#include "test-support.h"
#include "memory-manager.h"

#define DEBUG 0    // prints debug information

/*
Prints the ZmodF_t to stdout in hex, each limb in a separate block,
most significant limb (i.e. the overflow limb) first.
*/
void ZmodF_print(ZmodF_t x, unsigned long n)
{
   long i;
   for (i = n; i >= 0; i--)
#if FLINT_BITS == 64
      printf("%016lx ", x[i]);
#else
      printf("%08lx ", x[i]);
#endif
}



/*
Prints each coefficient of the polynomial on a separate line.
*/
void ZmodF_poly_print(ZmodF_poly_t x)
{
   unsigned long k;
   for (k = 0; k < (1UL << x->depth); k++)
   {
      ZmodF_print(x->coeffs[k], x->n);
      printf("\n");
   }
}



int test__ZmodF_mul_fft_split()
{
   int success = 1;

   mpz_t x, y, z;
   mpz_init(x);
   mpz_init(y);
   mpz_init(z);
   mp_limb_t buf[300];

   unsigned long n;
   for (n = 1; n < 200 && success; n++)
   {
      unsigned long depth;
      for (depth = 0;
           ((n*FLINT_BITS) % (1 << depth) == 0) && success; depth++)
      {
         unsigned long bits = (n*FLINT_BITS) >> depth;
         unsigned long m = (bits-1)/FLINT_BITS + 1;
      
         ZmodF_poly_t poly;
         ZmodF_poly_init(poly, depth, m, 1);

#if DEBUG
         printf("n = %d, depth = %d, m = %d\n", n, depth, m);
#endif
         
         unsigned long trial;
         for (trial = 0; trial < 120; trial++)
         {
            random_limbs(buf, n);
            buf[n] = 0;
            mpn_to_mpz(x, buf, n);
            
            _ZmodF_mul_fft_split(poly, buf, n);

            unsigned long i;
            for (i = 0; i < (1 << depth); i++)
            {
               mpz_tdiv_r_2exp(y, x, bits);
               mpz_tdiv_q_2exp(x, x, bits);
               mpn_to_mpz(z, poly->coeffs[i], m+1);
               if (mpz_cmp(z, y))
                  success = 0;
            }
         }
         
         ZmodF_poly_clear(poly);
      }
   }
   
   mpz_clear(x);
   mpz_clear(y);
   mpz_clear(z);

   return success;
}


int test__ZmodF_mul_fft_combine()
{
   int success = 1;
   
   mpz_t x, y, p, q, r, s, total;
   mpz_init(x);
   mpz_init(y);
   mpz_init(s);
   mpz_init(r);
   mpz_init(q);
   mpz_init(p);
   mpz_init(total);

   mp_limb_t buf[300];

   unsigned long n;
   for (n = 1; n < 80 && success; n++)
   {
      unsigned long depth;
      for (depth = 0;
           ((n*FLINT_BITS) % (1 << depth) == 0) && success; depth++)
      {
         unsigned long m;
         for (m = 1; m < n/4 && success; m++)
         {
            unsigned long k;
            for (k = 0; k < 5 && success; k++)
            {
#if DEBUG
               printf("n = %ld, depth = %ld, m = %ld, k = %ld\n", n, depth, m, k);
#endif

               ZmodF_poly_t poly;
               ZmodF_poly_init(poly, depth, m+k, 1);

               // p := B^n + 1
               mpz_set_ui(p, 1);
               mpz_mul_2exp(p, p, n*FLINT_BITS);
               mpz_add_ui(p, p, 1);
               
               // q := (B^m + 1)*B^k
               mpz_set_ui(q, 1);
               mpz_mul_2exp(q, q, m*FLINT_BITS);
               mpz_add_ui(q, q, 1);
               mpz_mul_2exp(q, q, k*FLINT_BITS);

               // r := B^(m+k) - 1
               mpz_set_ui(r, 1);
               mpz_mul_2exp(r, r, (m+k)*FLINT_BITS);
               mpz_sub_ui(r, r, 1);
               
               // s := B^(m+k)/2
               mpz_set_ui(s, 1);
               mpz_mul_2exp(s, s, (m+k)*FLINT_BITS - 1);
                  
               unsigned long trial;
               for (trial = 0; trial < 20 && success; trial++)
               {
                  mpz_set_ui(total, 0);
               
                  long i;
                  for (i = (1 << depth) - 1; i >= 0; i--)
                  {
                     // select random x in (0, B^(m+k))
                     mpz_set_ui(x, 0);
                     while (!mpz_sgn(x))
                     {
                        mpz_rrandomb(x, randstate, (m+k)*FLINT_BITS);
                        if (random_ulong(2))    // to get high bit 0 sometimes
                           mpz_sub(x, r, x);
                     }
                        
                     // push it down to (-B^(m+k)/2, B^(m+k)/2)
                     mpz_sub(x, x, s);
                     
                     // add it to running total
                     mpz_mul_2exp(total, total, (n*FLINT_BITS) >> depth);
                     mpz_add(total, total, x);
                     
                     // normalise it into [0, q), and store in polynomial
                     mpz_mod(x, x, q);
                     mpz_to_mpn(poly->coeffs[i], m+k+1, x);
                  }
                  
                  // compare result to target function
                  _ZmodF_mul_fft_combine(buf, poly, m, k, n);
                  ZmodF_normalise(buf, n);
                  mpn_to_mpz(y, buf, n+1);
                  
                  mpz_mod(total, total, p);
                  if (mpz_cmp(total, y))
                     success = 0;
               }

               ZmodF_poly_clear(poly);
            }
         }
      }
   }

   mpz_clear(x);
   mpz_clear(y);
   mpz_clear(s);
   mpz_clear(r);
   mpz_clear(q);
   mpz_clear(p);
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

   unsigned long n;
   for (n = 3; n < 300 && success; n += 3)
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

      unsigned long trial;
      for (trial = 0; trial < 250 && success; trial++)
      {
         random_limbs(in, n);
         in[n] = 0;
         mpn_to_mpz(x, in, n+1);

         _ZmodF_mul_threeway_reduce1(out1, in, n/3);
         ZmodF_normalise(out1, n/3);
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

   unsigned long n;
   for (n = 1; n < 100 && success; n++)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif

      // p = B^n + 1
      mpz_set_ui(p, 1);
      mpz_mul_2exp(p, p, n*FLINT_BITS);
      mpz_add_ui(p, p, 1);

      ZmodF_mul_info_t info;
      ZmodF_mul_info_init_plain(info, n, 0);

      unsigned long trial;
      for (trial = 0; trial < 1000 && success; trial++)
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
         
         ZmodF_mul_info_mul(info, out, in1, in1);
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

   unsigned long n;
   for (n = 3; n < 100 && success; n += 3)
   {
#if DEBUG
      printf("n = %d\n", n);
#endif

      ZmodF_mul_info_t info_plain, info_threeway;
      ZmodF_mul_info_init_threeway(info_threeway, n, 0);
      ZmodF_mul_info_init_plain(info_plain, n, 0);

      unsigned long trial;
      for (trial = 0; trial < 50000 && success; trial++)
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
         
         ZmodF_mul_info_mul(info_plain, out_plain, in1, in1);
         ZmodF_mul_info_mul(info_threeway, out_threeway, in1, in1);
         
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



int test_ZmodF_mul_info_mul_fft()
{
   int success = 1;

   mp_limb_t in1[1000];
   mp_limb_t in2[1000];
   mp_limb_t out_plain[1000];
   mp_limb_t out_fft[1000];

   mpz_t x;
   mpz_init(x);

   unsigned long n;
   for (n = 1; n < 300 && success; n++)
   {
      unsigned long depth;
      for (depth = 1;
           (n*FLINT_BITS) % (1 << depth) == 0
           && (depth <= FLINT_LG_BITS_PER_LIMB + 4)
           && success; depth++)
      {
         unsigned long input_bits = (n*FLINT_BITS) >> depth;
         unsigned long output_bits = 2*input_bits + 1 + depth;
         unsigned long target_m =
                            ((output_bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1;
      
         unsigned long m;
         for (m = target_m - 2; m <= target_m + 3 && success; m++)
         {
            if ((m*FLINT_BITS) % (1 << depth) != 0)
               continue;

            unsigned long k;
            for (k = 0; k <= 2 && success; k++)
            {
               if (m + k < target_m)
                  continue;
			   if (m + k > n)
			      continue;
               if (k > m)
                  continue;
            
#if DEBUG
               printf("n = %ld, depth = %ld, m = %ld, k = %ld\n", n, depth, m, k);
#endif

               ZmodF_mul_info_t info_plain, info_fft;
               ZmodF_mul_info_init_plain(info_plain, n, 0);
               ZmodF_mul_info_init_fft(info_fft, n, depth, m, k, 0);

               unsigned long trial;
               for (trial = 0; trial < 10 && success; trial++)
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
                  ZmodF_mul_info_mul(info_fft, out_fft, in1, in2);

                  ZmodF_normalise(out_plain, n);
                  ZmodF_normalise(out_fft, n);

                  if (mpn_cmp(out_plain, out_fft, n+1))
                     success = 0;

                  // test squaring
         
                  ZmodF_mul_info_mul(info_plain, out_plain, in1, in1);
                  ZmodF_mul_info_mul(info_fft, out_fft, in1, in1);

                  ZmodF_normalise(out_plain, n);
                  ZmodF_normalise(out_fft, n);

                  if (mpn_cmp(out_plain, out_fft, n+1))
                     success = 0;
               }

               ZmodF_mul_info_clear(info_fft);
               ZmodF_mul_info_clear(info_plain);
            }
         }
      }
   }

   mpz_clear(x);

   return success;
}



/****************************************************************************

   Main test functions

****************************************************************************/

void ZmodF_mul_test_all()
{
   int success, all_success = 1;

   RUN_TEST(_ZmodF_mul_fft_split);
   RUN_TEST(_ZmodF_mul_fft_combine);
   RUN_TEST(_ZmodF_mul_threeway_reduce);
   RUN_TEST(ZmodF_mul_info_mul_plain);
   RUN_TEST(ZmodF_mul_info_mul_threeway);
   RUN_TEST(ZmodF_mul_info_mul_fft);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   test_support_init();
   ZmodF_mul_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();
   
   return 0;
}


// end of file ****************************************************************
