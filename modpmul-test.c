/******************************************************************************

 modpmul-test.c
 Test routines for modpmul.c.

 Copyright (C) 2006, David Harvey

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <gmp.h>
#include "flint.h"
#include "profiler.h"
#include "modpmul.h"


void print_array(unsigned long* data, unsigned long length)
{
   for (int i = 0; i < length; i++)
      printf("%016lx  %ld\n", data[i], data[i]);
   printf("\n");
}


/*
 stores a prime p with FLINT_REDC_BITS bits at *p,
 and a root of unity of order 2^n at *w.
 
 These are just for testing purposes in this module.
*/
void get_test_prime(unsigned long* p, unsigned long *w, unsigned long n)
{
#if FLINT_REDC_BITS == 63
   *p = 0x5700000000000001UL;   // a 63 bit prime
   *w = 0x3e003802f87a45ebUL;   // has order 2^56 mod p
   unsigned long order = 56;
#elif FLINT_REDC_BITS == 31
   *p = 0x78000001UL;           // a 31 bit prime
   *w = 0x1a427a41UL;           // has order 2^27 mod p
   unsigned long order = 27;
#else
#error get_test_prime only works for 32-bit or 64-bit machines
#endif
   assert(n <= order);
   
   // power it up to the requested root of unity
   for (unsigned long i = 0; i < order - n; i++)
      *w = mod_mul(*w, *w, *p);
}


void test_modp_fft_basecase(unsigned long max_n, gmp_randstate_t* state)
{
   for (unsigned long n = 3; n <= max_n; n++)
   {
      printf("testing modp_fft_basecase (length 2^%ld)... ", n);
      fflush(stdout);
   
      unsigned long p;
      unsigned long w;
      get_test_prime(&p, &w, n);

      struct redc_precomp_t redc_info;
      redc_precomp_init(&redc_info, p);
      
      unsigned long Rw = convert_to_redc2(w, &redc_info);

      unsigned long length = 1 << n;
      unsigned long i, j, k;
      unsigned long* data1;
      unsigned long* data2;
      unsigned long* data3;
      
      data1 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data2 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data3 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      
      // generate random data to transform
      for (i = 0; i < length; i++)
         data1[i] = gmp_urandomm_ui(*state, p);
      
      // do naive transform (data1 => data2)
      mpz_t accum, w_pow, w_pow_pow;
      mpz_init(accum);
      mpz_init(w_pow);
      mpz_init(w_pow_pow);
      mpz_set_ui(w_pow, 1);
      for (i = 0; i < length; i++)
      {
         mpz_set_ui(accum, 0);
         mpz_set_ui(w_pow_pow, 1);
         for (j = 0; j < length; j++)
         {
            mpz_addmul_ui(accum, w_pow_pow, data1[j]);
            mpz_mul(w_pow_pow, w_pow_pow, w_pow);
            mpz_mod_ui(w_pow_pow, w_pow_pow, p);
         }
         data2[i] = mpz_fdiv_ui(accum, p);
         mpz_mul_ui(w_pow, w_pow, w);
         mpz_mod_ui(w_pow, w_pow, p);
      }
      mpz_clear(accum);
      mpz_clear(w_pow);
      mpz_clear(w_pow_pow);
      
      // bit-reverse output (data2 => data3)
      for (i = 0; i < length; i++)
      {
         j = 0;
         for (k = 0; k < n; k++)
            j = (j << 1) + ((i >> k) & 1);
         data3[j] = data2[i];
      }
      
      // do fast transform (data1 => data1)
      struct modp_fft_precomp_t info;
      modp_fft_precomp_init(&info, &redc_info, n, Rw, 12345); // force basecase
      modp_fft(&info, data1, 0);
      modp_fft_precomp_clear(&info);
      
      // check FFT results == naive DFT results
      for (i = 0; i < length; i++)
         assert(data1[i] == data3[i]);
      
      free(data1);
      free(data2);
      free(data3);
      
      printf("ok\n");
   }
}


void test_modp_fft(unsigned long max_n, gmp_randstate_t* state)
{
   for (unsigned long n = 3; n <= max_n; n++)
   {
      printf("testing modp_fft (length 2^%ld)... ", n);
      fflush(stdout);
   
      unsigned long p;
      unsigned long w;
      get_test_prime(&p, &w, n);

      struct redc_precomp_t redc_info;
      redc_precomp_init(&redc_info, p);
      
      unsigned long Rw = convert_to_redc2(w, &redc_info);

      unsigned long length = 1 << n;
      unsigned long i, j, k;
      unsigned long* data1;
      unsigned long* data2;
      unsigned long* data3;
      struct modp_fft_precomp_t info;
      
      data1 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data2 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data3 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      
      // generate random data to transform
      for (i = 0; i < length; i++)
         data1[i] = gmp_urandomm_ui(*state, p);
      
      // try transform using basecase (data2 => data2)
      for (i = 0; i < length; i++)
         data2[i] = data1[i];
      modp_fft_precomp_init(&info, &redc_info, n, Rw, 12345); // force basecase
      modp_fft(&info, data2, 0);
      modp_fft_precomp_clear(&info);
   
      // now using cache-friendly fft routine (data1 => data1)
      modp_fft_precomp_init(&info, &redc_info, n, Rw, 5); // force factoring
      modp_fft(&info, data1, 0);
      modp_fft_precomp_clear(&info);
      
      // compare results of two algorithms
      for (i = 0; i < length; i++)
         assert(data1[i] == data2[i]);
      
      free(data1);
      free(data2);
      free(data3);
      
      printf("ok\n");
   }
}


// tries FFT and IFFT on random data and checks if we got back to where
// we started
void test_modp_ifft(unsigned long max_n, gmp_randstate_t* state, int basecase)
{
   for (unsigned long n = 3; n <= max_n; n++)
   {
      if (basecase)
         printf("testing modp_ifft_basecase (length 2^%ld)... ", n);
      else
         printf("testing modp_ifft (length 2^%ld)... ", n);
      fflush(stdout);
   
      unsigned long p;
      unsigned long w;
      get_test_prime(&p, &w, n);

      struct redc_precomp_t redc_info;
      redc_precomp_init(&redc_info, p);
      
      unsigned long Rw = convert_to_redc2(w, &redc_info);

      unsigned long length = 1 << n;
      unsigned long i, j, k;
      unsigned long* data1;
      unsigned long* data2;
      struct modp_fft_precomp_t info;
      
      data1 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data2 = (unsigned long*) malloc(sizeof(unsigned long) * length);

      // generate random data to transform
      for (i = 0; i < length; i++)
         data1[i] = gmp_urandomm_ui(*state, p);
      
      for (i = 0; i < length; i++)
         data2[i] = data1[i];
      
      // forward transform
      modp_fft_precomp_init(&info, &redc_info, n, Rw,
                            FLINT_MODP_FFT_BASECASE_THRESHOLD);
      modp_fft(&info, data2, 0);
      modp_fft_precomp_clear(&info);

      // inverse transform
      modp_fft_precomp_init(&info, &redc_info, n, Rw, basecase ? 12345 : 5);
      modp_ifft(&info, data2);
      modp_fft_precomp_clear(&info);

      // remove normalising factor, and normalise ifft output
      for (i = 0; i < length; i++)
      {
         data1[i] = mod_mul(data1[i], 1 << n, p);
         if (data2[i] >= p)
            data2[i] -= p;
      }
      
      // check it's equal to the input
      for (i = 0; i < length; i++)
         assert(data1[i] == data2[i]);
      
      free(data1);
      free(data2);

      printf("ok\n");
   }
}


void test_modp_convolution(unsigned long max_n, gmp_randstate_t* state)
{
   for (unsigned long n = 3; n <= max_n; n++)
   {
      printf("testing modp_convolution (length 2^%ld)... ", n);
      fflush(stdout);
   
      unsigned long p;
      unsigned long w;
      get_test_prime(&p, &w, n);

      struct redc_precomp_t redc_info;
      redc_precomp_init(&redc_info, p);
      
      unsigned long Rw = convert_to_redc2(w, &redc_info);

      unsigned long length = 1 << n;
      unsigned long i, j, k;
      unsigned long* data1;
      unsigned long* data2;
      unsigned long* data3;
      
      data1 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data2 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      data3 = (unsigned long*) malloc(sizeof(unsigned long) * length);
      
      // generate random data to convolve
      for (i = 0; i < length; i++)
      {
         data1[i] = gmp_urandomm_ui(*state, p);
         data2[i] = gmp_urandomm_ui(*state, p);
         data3[i] = 0;
      }

      // do naive convolution (data1 * data2 => data3)
      for (i = 0; i < length; i++)
      {
         for (j = 0; j < length; j++)
         {
            unsigned long k = i + j;
            if (k >= length)
               k -= length;
            data3[k] += mod_mul(data1[i], data2[j], p);
            if (data3[k] >= p)
               data3[k] -= p;
         }
      }

      // do convolution via FFTs (data1 * data2 => data1)
      modp_convolution(&redc_info, n, Rw, data1, data2);

      // check FFT results == naive result
      for (i = 0; i < length; i++)
         assert(data1[i] == data3[i]);
      
      free(data1);
      free(data2);
      free(data3);
      
      printf("ok\n");
   }
}


void profile_modp_fft(gmp_randstate_t* state, int inverse)
{
   unsigned long p;
   unsigned long w;
   get_test_prime(&p, &w, 1);

   struct redc_precomp_t redc_info;
   redc_precomp_init(&redc_info, p);
   
   unsigned long max_n = 22;

   unsigned long* data =
           (unsigned long*) malloc(sizeof(unsigned long) * (1 << max_n));

   for (unsigned long i = 0; i < (1 << max_n); i++)
      data[i] = gmp_urandomm_ui(*state, p);

   for (unsigned long n = 3; n <= max_n; n++)
   {
      get_test_prime(&p, &w, n);
      unsigned long Rw = convert_to_redc2(w, &redc_info);

      unsigned long trials = 20000000 / (1 << n) / n;
      if (trials == 0)
         trials = 1;
         
      init_all_clocks();

      start_clock(0);

      for (unsigned long i = 0; i < trials; i++)
      {
         struct modp_fft_precomp_t info;
         modp_fft_precomp_init(&info, &redc_info, n, Rw,
                               FLINT_MODP_FFT_BASECASE_THRESHOLD);
         if (inverse)
            modp_ifft(&info, data);
         else
            modp_fft(&info, data, 0);
         modp_fft_precomp_clear(&info);
      }
      
      stop_clock(0);

      //      printf("%ld %lf\n", n, get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0 / trials);

      printf("%6.3lf\n", get_clock(0) / trials / n / (1 << n) * 50.0);
   }

   free(data);
}


void profile_modp_convolution(gmp_randstate_t* state)
{
   unsigned long p;
   unsigned long w;
   get_test_prime(&p, &w, 1);

   struct redc_precomp_t redc_info;
   redc_precomp_init(&redc_info, p);
   
   unsigned long max_n = 22;

   unsigned long* data1 =
           (unsigned long*) malloc(sizeof(unsigned long) * (1 << max_n));
   unsigned long* data2 =
           (unsigned long*) malloc(sizeof(unsigned long) * (1 << max_n));

   for (unsigned long i = 0; i < (1 << max_n); i++)
   {
      data1[i] = gmp_urandomm_ui(*state, p);
      data2[i] = gmp_urandomm_ui(*state, p);
   }

   for (unsigned long n = 3; n <= max_n; n++)
   {
      get_test_prime(&p, &w, n);
      unsigned long Rw = convert_to_redc2(w, &redc_info);

      unsigned long trials = 20000000 / (1 << n) / n;
      if (trials == 0)
         trials = 1;
         
      init_all_clocks();

      start_clock(0);

      for (unsigned long i = 0; i < trials; i++)
         modp_convolution(&redc_info, n, Rw, data1, data2);
      
      stop_clock(0);

      //      printf("%ld %lf\n", n, get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0 / trials);

      //      printf("%6.3lf\n", get_clock(0) / trials / n / (1 << n) * 50.0);
      printf("%ld %lf\n", n, get_clock(0) / trials / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0);
   }

   free(data1);
   free(data2);
}


int main()
{
   gmp_randstate_t state;
   gmp_randinit_default(state);

#if 1
   profile_modp_convolution(&state);
#elif 0
   test_modp_fft_basecase(10, &state);
   test_modp_ifft(10, &state, 1);
   test_modp_fft(20, &state);
   test_modp_ifft(20, &state, 0);
   test_modp_convolution(12, &state);
#else
   profile_modp_fft(&state, 0);
#endif

   return 0;
}
