/******************************************************************************

 ZmodF-test.c: test module for ZmodF module

 Copyright (C) 2007, David Harvey


TODO: establish and test overflow bit guarantees

 
******************************************************************************/

#include <stdio.h>
#include "ZmodF.h"

gmp_randstate_t ZmodF_test_randstate;
mpz_t global_mpz;   // to avoid frequent mpz_init calls


// ============================================================================
// Test case support code


/*
Prints the ZmodF_t to the given stream in hex, each limb in a separate block,
most significant limb (i.e. the overflow limb) first.
*/
void ZmodF_print(ZmodF_t x, unsigned long n)
{
   for (long i = n; i >= 0; i--)
#if FLINT_BITS_PER_LIMB == 64
      printf("%016lx ", x[i]);
#else
      printf("%08lx ", x[i]);
#endif
}


unsigned long random_ulong(unsigned long max)
{
   return gmp_urandomm_ui(ZmodF_test_randstate, max);
}


/*
Generates a random ZmodF_t with at most overflow_bits used in the
overflow limb. More precisely, the top (FLINT_BITS_PER_LIMB - overflow_bits)
bits will all be equal to the sign bit. It uses mpz_rrandomb to get long
strings of 0's and 1's.
*/
void ZmodF_random(ZmodF_t x, unsigned long n, unsigned long overflow_bits)
{
   ZmodF_zero(x, n);

   mpz_rrandomb(global_mpz, ZmodF_test_randstate, (n+1)*FLINT_BITS_PER_LIMB);
   mpz_export(x, NULL, -1, sizeof(mp_limb_t), 0, 0, global_mpz);

   // GMP has a "bug" where the top bit of the output of mpz_rrandomb
   // is always set. So we flip everything with probability 1/2.
   if (random_ulong(2))
      for (unsigned long i = 0; i <= n; i++)
         x[i] = ~x[i];

   // Now copy the sign bit downwards so that only overflow_bits bits are used.
   if ((mp_limb_signed_t) x[n] >= 0)
      x[n] &= (1UL << overflow_bits) - 1;
   else
      x[n] |= ~((1UL << overflow_bits) - 1);
}


#if FLINT_BITS_PER_LIMB == 64
#define SENTRY_LIMB 0x0123456789abcdefUL
#else
#define SETNRY_LIMB 0x01234567UL
#endif


#define MAX_COEFFS 5
#define MAX_N 30

mp_limb_t global_buf[MAX_COEFFS * (MAX_N + 3)];

ZmodF_t coeffs[MAX_COEFFS];
mpz_t coeffs_mpz_in[MAX_COEFFS];
mpz_t coeffs_mpz_out[MAX_COEFFS];
unsigned long global_n = 0;
mpz_t global_p;


/*
Converts given ZmodF_t into mpz_t format, reduced into [0, p) range.
Assumes global_n and global_p are set correctly.
*/
void ZmodF_convert_out(mpz_t output, ZmodF_t input)
{
   int negative = ((mp_limb_signed_t) input[global_n] < 0);
   
   if (negative)
      for (int i = 0; i <= global_n; i++)
         input[i] = ~input[i];
         
   mpz_import(output, global_n+1, -1, sizeof(mp_limb_t), 0, 0, input);
   
   if (negative)
   {
      mpz_add_ui(output, output, 1);
      mpz_neg(output, output);
      for (int i = 0; i <= global_n; i++)
         input[i] = ~input[i];
   }

   mpz_mod(output, output, global_p);
}


/*
y := x * 2^(s/2)  mod p    (using a very naive algorithm)
y may alias x
Assumes global_n and global_p are set correctly.
*/
void naive_mul_sqrt2exp(mpz_t y, mpz_t x, unsigned long s)
{
   static mpz_t temp;
   static int init = 0;
   
   if (!init)
   {
      mpz_init(temp);
      init = 1;
   }

   if (s & 1)
   {
      mpz_mul_2exp(y, x, s/2 + global_n*FLINT_BITS_PER_LIMB/4);
      mpz_mul_2exp(temp, y, global_n*FLINT_BITS_PER_LIMB/2);
      mpz_sub(y, temp, y);
      mpz_mod(y, y, global_p);
   }
   else
   {
      mpz_mul_2exp(y, x, s/2);
      mpz_mod(y, y, global_p);
   }
}


/*
y := x * 2^(-s/2)  mod p    (using a very naive algorithm)
y may alias x
Assumes global_n and global_p are set correctly.
*/
void naive_div_sqrt2exp(mpz_t y, mpz_t x, unsigned long s)
{
   naive_mul_sqrt2exp(y, x, 4*global_n*FLINT_BITS_PER_LIMB - s);
}



/*
Sets up "count" ZmodF_t's in the global array with random data.
Makes coeffs[0], ..., coeffs[count-1] point to those buffers. Adds sentry
limbs at the ends of each buffer. Converts each coefficient to mpz_t form
in coeffs_mpz_in. Sets global_n := n and global_p := B^n + 1.
*/
void setup_coeffs(unsigned long count, unsigned long n, unsigned long overflow_bits)
{
   assert(n <= MAX_N);
   assert(count <= MAX_COEFFS);

   // update global_p only if n has changed since last time
   if (n != global_n)
   {
      global_n = n;
      mpz_set_ui(global_p, 1);
      mpz_mul_2exp(global_p, global_p, n*FLINT_BITS_PER_LIMB);
      mpz_add_ui(global_p, global_p, 1);
   }

   // make pointers point to the right place
   coeffs[0] = global_buf + 1;
   for (int i = 1; i < count; i++)
      coeffs[i] = coeffs[i-1] + (n+3);
      
   // add sentry limbs
   for (int i = 0; i < count; i++)
      coeffs[i][-1] = coeffs[i][n+1] = SENTRY_LIMB;
      
   // generate random coefficients
   for (int i = 0; i < count; i++)
      ZmodF_random(coeffs[i], n, overflow_bits);
      
   // convert coefficients to coeffs_mpz_in
   for (int i = 0; i < count; i++)
      ZmodF_convert_out(coeffs_mpz_in[i], coeffs[i]);
}


/*
Checks that the sentries have not been overwritten, and that the first "count"
pointers in "coeffs" points to correct distinct buffers. Converts each
coefficient to mpz_t form in coeffs_mpz_out.

Returns 1 on success.
*/
int check_coeffs(unsigned long count, unsigned long n)
{
   // check sentry limbs
   for (int i = 0; i < count; i++)
   {
      if (coeffs[i][-1] != SENTRY_LIMB)
         return 0;
      if (coeffs[i][n+1] != SENTRY_LIMB)
         return 0;
   }

   // check pointers point to valid buffers
   for (int i = 0; i < count; i++)
   {
      unsigned long offset = coeffs[i] - global_buf;
      if (offset % (n+3) != 1)
         return 0;
      if ((offset - 1) / (n+3) >= count)
         return 0;
   }
   
   // check pointers point to distinct buffers
   for (int i = 0; i < count; i++)
      for (int j = i+1; j < count; j++)
         if (coeffs[i] == coeffs[j])
            return 0;

   // convert coefficients to coeffs_mpz_out
   for (int i = 0; i < count; i++)
      ZmodF_convert_out(coeffs_mpz_out[i], coeffs[i]);
   
   return 1;
}



// ============================================================================
// Actual test cases for ZmodF functions


int test_ZmodF_normalise()
{
   for (unsigned long n = 1; n <= 5; n++)
      for (unsigned long trial = 0; trial < 100000; trial++)
      {
         setup_coeffs(1, n, random_ulong(FLINT_BITS_PER_LIMB - 2));
         
         ZmodF_normalise(coeffs[0], n);

         if (!check_coeffs(1, n))
            return 0;

         // check normalised value is in [0, p)
         if (coeffs[0][n])
         {
            if (coeffs[0][n] != 1)
               return 0;
            for (unsigned long i = 0; i < n; i++)
               if (coeffs[0][i])
                  return 0;
         }
            
         // check output actually equals input mod p
         if (mpz_cmp(coeffs_mpz_in[0], coeffs_mpz_out[0]))
            return 0;
      }
   
   return 1;
}


int test_ZmodF_fast_reduce()
{
   for (unsigned long n = 1; n <= 5; n++)
      for (unsigned long trial = 0; trial < 100000; trial++)
      {
         setup_coeffs(1, n, random_ulong(FLINT_BITS_PER_LIMB - 2));
         
         ZmodF_fast_reduce(coeffs[0], n);

         if (!check_coeffs(1, n))
            return 0;

         // check high limb of normalised value is in [0, 2]
         if (coeffs[0][n] > 2)
            return 0;
            
         // check output actually equals input mod p
         if (mpz_cmp(coeffs_mpz_in[0], coeffs_mpz_out[0]))
            return 0;
      }
   
   return 1;
}


int test_ZmodF_neg()
{
   for (unsigned long n = 1; n <= 5; n++)
      for (unsigned long trial = 0; trial < 50000; trial++)
         for (int inplace = 0; inplace <= 1; inplace++)
         {
            setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 2));
            
            ZmodF_neg(coeffs[inplace], coeffs[0], n);

            if (!check_coeffs(2, n))
               return 0;

            mpz_neg(global_mpz, coeffs_mpz_out[inplace]);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_in[0], global_mpz))
               return 0;
         }

   return 1;
}


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


int test_ZmodF_short_div_2exp()
{
   for (unsigned long n = 1; n <= 3; n++)
      for (unsigned long trial = 0; trial < 2000; trial++)
         for (unsigned long s = 1; s < FLINT_BITS_PER_LIMB; s++)
            for (int inplace = 0; inplace <= 1; inplace++)
            {
               setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

               ZmodF_short_div_2exp(coeffs[inplace], coeffs[0], s, n);

               if (!check_coeffs(2, n))
                  return 0;

               naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[0], 2*s);
               if (mpz_cmp(coeffs_mpz_out[inplace], global_mpz))
                  return 0;
            }

   return 1;
}


int test_ZmodF_mul_Bexp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 20000; trial++)
         for (unsigned long s = 1; s < n; s++)
         {
            setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_mul_Bexp(coeffs[1], coeffs[0], s, n);

            if (!check_coeffs(2, n))
               return 0;

            naive_mul_sqrt2exp(global_mpz, coeffs_mpz_in[0],
                               2*FLINT_BITS_PER_LIMB*s);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_div_Bexp_sub()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
         for (unsigned long s = 1; s < n; s++)
            for (int inbuf1 = 0; inbuf1 <= 2; inbuf1++)
               for (int inbuf2 = 1; inbuf2 <= 2; inbuf2++)
               {
                  setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

                  ZmodF_div_Bexp_sub(coeffs[0], coeffs[inbuf1], coeffs[inbuf2],
                                     s, n);

                  if (!check_coeffs(3, n))
                     return 0;

                  naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[inbuf2],
                                     2*FLINT_BITS_PER_LIMB*s);
                  mpz_sub(global_mpz, coeffs_mpz_in[inbuf1], global_mpz);
                  mpz_mod(global_mpz, global_mpz, global_p);
                  if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
                     return 0;
               }

   return 1;
}


int test_ZmodF_div_Bexp_add()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
         for (unsigned long s = 1; s < n; s++)
            for (int inbuf1 = 0; inbuf1 <= 2; inbuf1++)
               for (int inbuf2 = 1; inbuf2 <= 2; inbuf2++)
               {
                  setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

                  ZmodF_div_Bexp_add(coeffs[0], coeffs[inbuf1], coeffs[inbuf2],
                                     s, n);

                  if (!check_coeffs(3, n))
                     return 0;

                  naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[inbuf2],
                                     2*FLINT_BITS_PER_LIMB*s);
                  mpz_add(global_mpz, coeffs_mpz_in[inbuf1], global_mpz);
                  mpz_mod(global_mpz, global_mpz, global_p);
                  if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
                     return 0;
               }

   return 1;
}


int test_ZmodF_sub_mul_Bexp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
         for (unsigned long s = 1; s < n; s++)
            for (int inbuf1 = 1; inbuf1 <= 2; inbuf1++)
               for (int inbuf2 = 1; inbuf2 <= 2; inbuf2++)
               {
                  setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

                  ZmodF_sub_mul_Bexp(coeffs[0], coeffs[inbuf1], coeffs[inbuf2],
                                     s, n);

                  if (!check_coeffs(3, n))
                     return 0;

                  mpz_sub(global_mpz, coeffs_mpz_in[inbuf1],
                          coeffs_mpz_in[inbuf2]);
                  naive_mul_sqrt2exp(global_mpz, global_mpz,
                                     2*FLINT_BITS_PER_LIMB*s);
                  if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
                     return 0;
               }

   return 1;
}


int test_ZmodF_mul_pseudosqrt2_n_odd()
{
   for (unsigned long n = 1; n <= 9; n += 2)
      for (unsigned long trial = 0; trial < 8000; trial++)
         for (unsigned long s = 0; s < 2*n; s++)
         {
            setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 4));

            ZmodF_mul_pseudosqrt2_n_odd(coeffs[1], coeffs[0], s, n);

            if (!check_coeffs(2, n))
               return 0;

            mpz_mul_2exp(global_mpz, coeffs_mpz_in[0], n*FLINT_BITS_PER_LIMB/2);
            mpz_sub(global_mpz, coeffs_mpz_in[0], global_mpz);
            mpz_mul_2exp(global_mpz, global_mpz, s*FLINT_BITS_PER_LIMB);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_mul_pseudosqrt2_n_even()
{
   for (unsigned long n = 2; n <= 10; n += 2)
      for (unsigned long trial = 0; trial < 8000; trial++)
         for (unsigned long s = 0; s < 2*n; s++)
         {
            setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_mul_pseudosqrt2_n_even(coeffs[1], coeffs[0], s, n);

            if (!check_coeffs(2, n))
               return 0;

            mpz_mul_2exp(global_mpz, coeffs_mpz_in[0], n*FLINT_BITS_PER_LIMB/2);
            mpz_sub(global_mpz, coeffs_mpz_in[0], global_mpz);
            mpz_mul_2exp(global_mpz, global_mpz, s*FLINT_BITS_PER_LIMB);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_mul_2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 500; trial++)
         for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
         {
            setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 3));

            ZmodF_mul_2exp(coeffs[1], coeffs[0], s, n);

            if (!check_coeffs(2, n))
               return 0;

            naive_mul_sqrt2exp(global_mpz, coeffs_mpz_in[0], 2*s);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_mul_sqrt2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 500; trial++)
         for (unsigned long s = 0; s < 2*n*FLINT_BITS_PER_LIMB; s++)
         {
            setup_coeffs(2, n, random_ulong(FLINT_BITS_PER_LIMB - 6));

            ZmodF_mul_sqrt2exp(coeffs[1], coeffs[0], s, n);

            if (!check_coeffs(2, n))
               return 0;

            naive_mul_sqrt2exp(global_mpz, coeffs_mpz_in[0], s);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_sub_mul_2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 100; trial++)
         for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
            for (int inbuf1 = 1; inbuf1 <= 2; inbuf1++)
               for (int inbuf2 = 1; inbuf2 <= 2; inbuf2++)
               {
                  setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

                  ZmodF_sub_mul_2exp(coeffs[0], coeffs[inbuf1], coeffs[inbuf2],
                                     s, n);

                  if (!check_coeffs(3, n))
                     return 0;

                  mpz_sub(global_mpz, coeffs_mpz_in[inbuf1],
                          coeffs_mpz_in[inbuf2]);
                  naive_mul_sqrt2exp(global_mpz, global_mpz, 2*s);
                  if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
                     return 0;
               }

   return 1;
}


int test_ZmodF_forward_butterfly_Bexp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 25000; trial++)
         for (unsigned long s = 1; s < n; s++)
         {
            setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_forward_butterfly_Bexp(&coeffs[0], &coeffs[1], &coeffs[2],
                                         s, n);

            if (!check_coeffs(3, n))
               return 0;

            mpz_sub(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
            naive_mul_sqrt2exp(global_mpz, global_mpz,
                               2*FLINT_BITS_PER_LIMB*s);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
            
            mpz_add(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_forward_butterfly_2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 400; trial++)
         for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
         {
            setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_forward_butterfly_2exp(&coeffs[0], &coeffs[1], &coeffs[2],
                                         s, n);

            if (!check_coeffs(3, n))
               return 0;

            mpz_sub(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
            naive_mul_sqrt2exp(global_mpz, global_mpz, 2*s);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
            
            mpz_add(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_forward_butterfly_sqrt2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 400; trial++)
         for (unsigned long s = 0; s < 2*n*FLINT_BITS_PER_LIMB; s++)
         {
            setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_forward_butterfly_sqrt2exp(&coeffs[0], &coeffs[1],
                                             &coeffs[2], s, n);

            if (!check_coeffs(3, n))
               return 0;

            mpz_sub(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
            naive_mul_sqrt2exp(global_mpz, global_mpz, s);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
            
            mpz_add(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_inverse_butterfly_2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 400; trial++)
         for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
         {
            setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_inverse_butterfly_2exp(&coeffs[0], &coeffs[1], &coeffs[2],
                                         s, n);

            if (!check_coeffs(3, n))
               return 0;

            naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[1], 2*s);
            mpz_add(global_mpz, coeffs_mpz_in[0], global_mpz);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
               return 0;
            
            naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[1], 2*s);
            mpz_sub(global_mpz, coeffs_mpz_in[0], global_mpz);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_inverse_butterfly_sqrt2exp()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 400; trial++)
         for (unsigned long s = 0; s < 2*n*FLINT_BITS_PER_LIMB; s++)
         {
            setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

            ZmodF_inverse_butterfly_sqrt2exp(&coeffs[0], &coeffs[1],
                                             &coeffs[2], s, n);

            if (!check_coeffs(3, n))
               return 0;

            naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[1], s);
            mpz_add(global_mpz, coeffs_mpz_in[0], global_mpz);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
               return 0;
            
            naive_div_sqrt2exp(global_mpz, coeffs_mpz_in[1], s);
            mpz_sub(global_mpz, coeffs_mpz_in[0], global_mpz);
            mpz_mod(global_mpz, global_mpz, global_p);
            if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
               return 0;
         }

   return 1;
}


int test_ZmodF_simple_butterfly()
{
   for (unsigned long n = 1; n <= 6; n++)
      for (unsigned long trial = 0; trial < 4000; trial++)
      {
         setup_coeffs(3, n, random_ulong(FLINT_BITS_PER_LIMB - 2));

         ZmodF_simple_butterfly(&coeffs[0], &coeffs[1], &coeffs[2], n);

         if (!check_coeffs(3, n))
            return 0;

         mpz_sub(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
         mpz_mod(global_mpz, global_mpz, global_p);
         if (mpz_cmp(coeffs_mpz_out[1], global_mpz))
            return 0;
         
         mpz_add(global_mpz, coeffs_mpz_in[0], coeffs_mpz_in[1]);
         mpz_mod(global_mpz, global_mpz, global_p);
         if (mpz_cmp(coeffs_mpz_out[0], global_mpz))
            return 0;
      }

   return 1;
}



#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void ZmodF_test_all()
{
   int success, all_success = 1;

   RUN_TEST(ZmodF_normalise);
   RUN_TEST(ZmodF_fast_reduce);
   RUN_TEST(ZmodF_neg);
   RUN_TEST(ZmodF_mul);
   RUN_TEST(ZmodF_sqr);
   RUN_TEST(ZmodF_short_div_2exp);
   RUN_TEST(ZmodF_mul_Bexp);
   RUN_TEST(ZmodF_div_Bexp_sub);
   RUN_TEST(ZmodF_div_Bexp_add);
   RUN_TEST(ZmodF_sub_mul_Bexp);
   RUN_TEST(ZmodF_mul_pseudosqrt2_n_odd);
   RUN_TEST(ZmodF_mul_pseudosqrt2_n_even);
   RUN_TEST(ZmodF_mul_2exp);
   RUN_TEST(ZmodF_mul_sqrt2exp);
   RUN_TEST(ZmodF_sub_mul_2exp);
   RUN_TEST(ZmodF_forward_butterfly_Bexp);
   RUN_TEST(ZmodF_forward_butterfly_2exp);
   RUN_TEST(ZmodF_forward_butterfly_sqrt2exp);
   RUN_TEST(ZmodF_inverse_butterfly_2exp);
   RUN_TEST(ZmodF_inverse_butterfly_sqrt2exp);
   RUN_TEST(ZmodF_simple_butterfly);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(ZmodF_test_randstate);
   mpz_init(global_mpz);
   mpz_init(global_p);
   for (int i = 0; i < MAX_COEFFS; i++)
   {
      mpz_init(coeffs_mpz_in[i]);
      mpz_init(coeffs_mpz_out[i]);
   }
   
   ZmodF_test_all();
   
   return 0;
}



// end of file ****************************************************************
