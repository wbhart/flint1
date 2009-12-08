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

ZmodF_poly-test.c: test module for ZmodF_poly module

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "fmpz_poly.h"
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
   
   unsigned long i;
   for (i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits+1);
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
   
   unsigned long i;
   for (i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits+1);
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

/****************************************************************************

   Test code for Fourier Transform Routines

****************************************************************************/


/*
Prints the ZmodF_t, each limb in a separate block, most significant limb
(i.e. the overflow limb) first.
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


/*
Generates a random ZmodF_poly_t with at most overflow_bits used in the
overflow limb for each coefficient.

The ZmodF_poly_t should already be initialised. This function ignores the
"length" attribute.
*/
void ZmodF_poly_random(ZmodF_poly_t x, unsigned long overflow_bits)
{
   unsigned long n = x->n;
   
   mpz_t temp;
   mpz_init(temp);

   unsigned long k;
   for (k = 0; k < (1UL << x->depth); k++)
   {
      ZmodF_t y = x->coeffs[k];
   
      ZmodF_zero(y, n);
      mpz_rrandomb(temp, randstate, (n+1)*FLINT_BITS);
      mpz_export(y, NULL, -1, sizeof(mp_limb_t), 0, 0, temp);

      // GMP has a "bug" where the top bit of the output of mpz_rrandomb
      // is always set. So we flip everything with probability 1/2.
      unsigned long i;
      if (random_ulong(2))
         for (i = 0; i <= n; i++)
            y[i] = ~y[i];

      // Copy the sign bit downwards so that only overflow_bits bits are used.
      if ((mp_limb_signed_t) y[n] >= 0)
         y[n] &= (1UL << overflow_bits) - 1;
      else
         y[n] |= ~((1UL << overflow_bits) - 1);
   }

   mpz_clear(temp);
}


mpz_t global_p;
unsigned long global_n = 0;


// Sets:
// global_n := n,
// global_p = 2^(FLINT_BITS*n) + 1
void set_global_n(unsigned long n)
{
   if (n != global_n)
   {
      global_n = n;
      mpz_set_ui(global_p, 1);
      mpz_mul_2exp(global_p, global_p, n*FLINT_BITS);
      mpz_add_ui(global_p, global_p, 1);
   }
}


/*
Converts given ZmodF_t into mpz_t format, reduced into [0, p) range.
Assumes global_n and global_p are set correctly.
*/
void ZmodF_convert_out(mpz_t output, ZmodF_t input)
{
   int negative = ((mp_limb_signed_t) input[global_n] < 0);
   
   int i;
   if (negative)
      for (i = 0; i <= global_n; i++)
         input[i] = ~input[i];
         
   mpz_import(output, global_n+1, -1, sizeof(mp_limb_t), 0, 0, input);
   
   if (negative)
   {
      mpz_add_ui(output, output, 1);
      mpz_neg(output, output);
      int i;
      for (i = 0; i <= global_n; i++)
         input[i] = ~input[i];
   }

   mpz_mod(output, output, global_p);
}


/*
Converts input polynomial to mpz_poly format. Each output coefficient is
normalised into [0, p). All 2^depth coefficients are converted.

Assumes that output is already initialised.
*/
void ZmodF_poly_convert_out(mpz_poly_t output, ZmodF_poly_t input)
{
   unsigned long size = 1UL << input->depth;
   unsigned long n = input->n;
   
   mpz_poly_ensure_alloc(output, size);
   set_global_n(n);
   
   unsigned long k;
   for (k = 0; k < size; k++)
      ZmodF_convert_out(output->coeffs[k], input->coeffs[k]);
      
   output->length = size;
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
      mpz_mul_2exp(y, x, s/2 + global_n*FLINT_BITS/4);
      mpz_mul_2exp(temp, y, global_n*FLINT_BITS/2);
      mpz_sub(y, temp, y);
      mpz_mod(y, y, global_p);
   }
   else
   {
      mpz_mul_2exp(y, x, s/2);
      mpz_mod(y, y, global_p);
   }
}


// root and twist are powers of sqrt2
void naive_FFT(mpz_poly_t x, unsigned long depth, unsigned long root,
               unsigned long twist, unsigned long n)
{
   static mpz_t temp;
   static int init = 0;
   
   if (!init)
   {
      mpz_init(temp);
      init = 1;
   }

   unsigned long size = 1UL << depth;
   
   unsigned long d;
   for (d = 0; d < depth; d++)
   {
      unsigned long half = 1UL << (depth - d - 1);
      unsigned long start;
      for (start = 0; start < size; start += 2*half)
      {
         unsigned long i;
         for (i = 0; i < half; i++)
         {
            mpz_t* a = &x->coeffs[start + i];
            mpz_t* b = &x->coeffs[start + half + i];
            mpz_add(temp, *a, *b);
            mpz_sub(*b, *a, *b);
            naive_mul_sqrt2exp(*b, *b, twist + i*root);
            mpz_mod(*a, temp, global_p);
         }
      }
      root <<= 1;
      twist <<= 1;
   }
}



int test__ZmodF_poly_FFT_iterative_case(
         unsigned long depth, unsigned long nonzero, unsigned long length,
         unsigned long twist, unsigned long n)
{
   mpz_poly_t poly1, poly2;
   ZmodF_poly_t f;

   unsigned long size = 1UL << depth;
   unsigned long root = 4*n*FLINT_BITS / size;
                  
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   ZmodF_poly_init(f, depth, n, 1);

   int success = 1;
   set_global_n(n);
         
   ZmodF_poly_random(f, 4);
   ZmodF_poly_convert_out(poly1, f);
   unsigned long i;
   for (i = nonzero; i < size; i++)
      mpz_set_ui(poly1->coeffs[i], 0);

   naive_FFT(poly1, depth, root, twist, n);

   _ZmodF_poly_FFT_iterative(f->coeffs, depth, 1, nonzero, length,
                            twist, n, f->scratch);
   ZmodF_poly_convert_out(poly2, f);

   for (i = 0; i < length; i++)
      if (mpz_cmp(poly1->coeffs[i], poly2->coeffs[i]))
         success = 0;
   
   ZmodF_poly_clear(f);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);

   return success;
}


int test__ZmodF_poly_FFT_iterative()
{
   int success = 1;

   unsigned long depth;
   for (depth = 1; depth <= 11 && success; depth++)
   {
      unsigned long size = 1UL << depth;
   
      // need 4*n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         unsigned long num_trials = 40000 / (1 << depth);
         unsigned long trial;
         for (trial = 0; trial < num_trials && success; trial++)
         {
            unsigned long nonzero, length, twist, root;
            
            if (depth == 0)
               nonzero = length = 1;
            else
            {
               nonzero = random_ulong(size-1) + 1;
               length = random_ulong(size-1) + 1;
            }

            twist = random_ulong(4*n*FLINT_BITS / size);
            success = success && test__ZmodF_poly_FFT_iterative_case(
                                           depth, nonzero, length, twist, n);
         }
      }
   }

   return success;
}


int test__ZmodF_poly_FFT_factor_case(
         unsigned long rows_depth, unsigned long cols_depth,
         unsigned long nonzero, unsigned long length,
         unsigned long twist, unsigned long n)
{
   mpz_poly_t poly1, poly2;
   ZmodF_poly_t f;

   unsigned long depth = rows_depth + cols_depth;
   unsigned long size = 1UL << depth;
   unsigned long root = 4*n*FLINT_BITS / size;
                  
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   ZmodF_poly_init(f, depth, n, 1);

   int success = 1;
   set_global_n(n);
   
   ZmodF_poly_random(f, 4);
   ZmodF_poly_convert_out(poly1, f);
   unsigned long i;
   for (i = nonzero; i < size; i++)
      mpz_set_ui(poly1->coeffs[i], 0);

   naive_FFT(poly1, depth, root, twist, n);

   _ZmodF_poly_FFT_factor(f->coeffs, rows_depth, cols_depth, 1, nonzero,
                         length, twist, n, f->scratch);
   ZmodF_poly_convert_out(poly2, f);

   for (i = 0; i < length; i++)
      if (mpz_cmp(poly1->coeffs[i], poly2->coeffs[i]))
         success = 0;
   
   ZmodF_poly_clear(f);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);

   return success;
}


int test__ZmodF_poly_FFT_factor()
{
   int success = 1;
   
   unsigned long depth, depth1;
   for (depth = 2; depth <= 6 && success; depth++)
      for (depth1 = 1; depth1 < depth && success; depth1++)
      {
         unsigned long depth2 = depth - depth1;
         unsigned long size = 1UL << depth;
      
         // need 4*n*FLINT_BITS divisible by 2^depth
         unsigned long n = size / (4*FLINT_BITS);
         if (n == 0)
            n = 1;
         
#if DEBUG
         printf("depth1 = %d, depth2 = %d, n = %d\n", depth1, depth2, n);
#endif

         unsigned long length, nonzero;
         for (length = 1; length <= size; length++)
            for (nonzero = 1; nonzero <= size; nonzero++)
            {
               unsigned long num_trials = 1000000 / (1 << (3*depth));
               if (num_trials == 0)
                  num_trials = 1;
               unsigned long trial;
               for (trial = 0; trial < num_trials; trial++)
               {
                  unsigned long twist = random_ulong(
                                             4*n*FLINT_BITS / size);
                  success = success && test__ZmodF_poly_FFT_factor_case(
                                 depth1, depth2, nonzero, length, twist, n);
               }
            }
      }

   return success;
}



// root and twist are powers of sqrt2
void naive_IFFT(mpz_poly_t x, unsigned long depth, unsigned long root,
                unsigned long twist, unsigned long n)
{
   static mpz_t temp;
   static int init = 0;
   
   if (!init)
   {
      mpz_init(temp);
      init = 1;
   }

   unsigned long size = 1UL << depth;
   root <<= (depth - 1);
   twist <<= (depth - 1);
   
   unsigned long d;
   for (d = 0; d < depth; d++)
   {
      unsigned long half = 1UL << d;
      unsigned long start;
      for (start = 0; start < size; start += 2*half)
      {
         unsigned long i;
         for (i = 0; i < half; i++)
         {
            mpz_t* a = &x->coeffs[start + i];
            mpz_t* b = &x->coeffs[start + half + i];
            naive_mul_sqrt2exp(*b, *b, 4*n*FLINT_BITS - (twist + i*root));
            mpz_add(temp, *a, *b);
            mpz_sub(*b, *a, *b);
            mpz_mod(*a, temp, global_p);
            mpz_mod(*b, *b, global_p);
         }
      }
      root >>= 1;
      twist >>= 1;
   }
}


int test__ZmodF_poly_IFFT_recursive_case(
         unsigned long depth, unsigned long nonzero, unsigned long length,
         int extra, unsigned long twist, unsigned long n)
{
   mpz_poly_t poly1, poly2;
   ZmodF_poly_t f;
   mpz_t extra_coeff;
   mpz_init(extra_coeff);

   unsigned long size = 1UL << depth;
   unsigned long root = 4*n*FLINT_BITS / size;
                  
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   ZmodF_poly_init(f, depth, n, 1);

   int success = 1;
   set_global_n(n);

   // run truncated inverse transform on random data
   ZmodF_poly_random(f, 4);
   ZmodF_poly_convert_out(poly1, f);
   _ZmodF_poly_IFFT_recursive(f->coeffs, depth, 1, nonzero, length, extra,
                             twist, n, f->scratch);
   
   // reassemble the untransformed coefficients
   ZmodF_poly_convert_out(poly2, f);
   if (extra)
      // save extra coefficient if necessary
      mpz_set(extra_coeff, poly2->coeffs[length]);
   unsigned long i;
   for (i = length; i < nonzero; i++)
      mpz_set(poly2->coeffs[i], poly1->coeffs[i]);
   for (i = nonzero; i < size; i++)
      mpz_set_ui(poly2->coeffs[i], 0);
   
   // run forward transform on proposed untransformed coefficients
   naive_FFT(poly2, depth, root, twist, n);
   // rescale
   for (i = 0; i < size; i++)
      naive_mul_sqrt2exp(poly2->coeffs[i], poly2->coeffs[i],
                         2*(2*n*FLINT_BITS - depth));
   // check the first few agree with input
   for (i = 0; i < length; i++)
      if (mpz_cmp(poly2->coeffs[i], poly1->coeffs[i]))
         success = 0;
   // check the extra coefficient is correct too
   if (extra)
      if (mpz_cmp(poly2->coeffs[length], extra_coeff))
         success = 0;

   mpz_clear(extra_coeff);
   ZmodF_poly_clear(f);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);

   return success;
}



int test__ZmodF_poly_IFFT_recursive()
{
   int success = 1;

   unsigned long depth;
   for (depth = 1; depth <= 6 && success; depth++)
   {
      unsigned long size = 1UL << depth;
      
      // need 4*n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         unsigned long nonzero, length;
         int extra;
         for (nonzero = 1; nonzero <= size; nonzero++)
            for (extra = 0; extra < 2; extra++)
               for (length = 1-extra; length <= nonzero; length++)
               {
                  if (extra && length == size)
                     continue;

                  unsigned long num_trials = 100000 / (1 << (3*depth));
                  if (num_trials == 0)
                     num_trials = 1;
                  unsigned long trial;
                  for (trial = 0; trial < num_trials; trial++)
                  {
                     unsigned long twist = random_ulong(
                                               4*n*FLINT_BITS / size);
                     success = success && test__ZmodF_poly_IFFT_recursive_case(
                                    depth, nonzero, length, extra, twist, n);
                  }
               }
      }
   }

   return success;
}



int test__ZmodF_poly_IFFT_iterative_case(unsigned long depth,
                                        unsigned long twist, unsigned long n)
{
   mpz_poly_t poly1, poly2;
   ZmodF_poly_t f;
   mpz_t extra_coeff;
   mpz_init(extra_coeff);

   unsigned long size = 1UL << depth;
   unsigned long root = 4*n*FLINT_BITS / size;
                  
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   ZmodF_poly_init(f, depth, n, 1);

   int success = 1;
   set_global_n(n);

   ZmodF_poly_random(f, 4);
   ZmodF_poly_convert_out(poly1, f);
   naive_IFFT(poly1, depth, root, twist, n);
   _ZmodF_poly_IFFT_iterative(f->coeffs, depth, 1, twist, n, f->scratch);
   ZmodF_poly_convert_out(poly2, f);

   unsigned long i;
   for (i = 0; i < size; i++)
      if (mpz_cmp(poly1->coeffs[i], poly2->coeffs[i]))
         success = 0;

   mpz_clear(extra_coeff);
   ZmodF_poly_clear(f);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);

   return success;
}



int test__ZmodF_poly_IFFT_iterative()
{
   int success = 1;

   unsigned long depth;
   for (depth = 1; depth <= 9 && success; depth++)
   {
      unsigned long size = 1UL << depth;
      
      // need 4*n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         unsigned long num_trials = 100000 / (1 << (depth));
         if (num_trials == 0)
            num_trials = 1;
         unsigned long trial;
         for (trial = 0; trial < num_trials; trial++)
         {
            unsigned long twist = random_ulong(
                                      4*n*FLINT_BITS / size);
            success = success && test__ZmodF_poly_IFFT_iterative_case(
                                                            depth, twist, n);
         }
      }
   }

   return success;
}



int test__ZmodF_poly_IFFT()
{
   mpz_poly_t poly1, poly2;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_t extra_coeff;
   mpz_init(extra_coeff);

   int success = 1;

   unsigned long depth;
   for (depth = 1; depth <= 11 && success; depth++)
   {
      unsigned long size = 1UL << depth;
   
      // need 4*n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
         ZmodF_poly_t f;
         ZmodF_poly_init(f, depth, n, 1);

#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         set_global_n(n);
         
         unsigned long num_trials = 40000 / (1 << depth);
         unsigned long trial;
         for (trial = 0; trial < num_trials; trial++)
         {
            unsigned long nonzero, length, twist, root;
            int extra = random_ulong(2);
            
            if (depth == 0)
            {
               nonzero = 1;
               length = 1 - extra;
            }
            else
            {
               nonzero = random_ulong(size-1) + 1;
               length = random_ulong(nonzero) + 1 - extra;
            }

            root = 4*n*FLINT_BITS / size;
            twist = random_ulong(root);
            
            // run truncated inverse transform on random data
            ZmodF_poly_random(f, 4);
            ZmodF_poly_convert_out(poly1, f);
            _ZmodF_poly_IFFT(f->coeffs, depth, 1, nonzero, length, extra,
                            twist, n, f->scratch);
            
            // reassemble the untransformed coefficients
            ZmodF_poly_convert_out(poly2, f);
            if (extra)
               // save extra coefficient if necessary
               mpz_set(extra_coeff, poly2->coeffs[length]);
            unsigned long i;
            for (i = length; i < nonzero; i++)
               mpz_set(poly2->coeffs[i], poly1->coeffs[i]);
            for (i = nonzero; i < size; i++)
               mpz_set_ui(poly2->coeffs[i], 0);
            
            // run forward transform on proposed untransformed coefficients
            naive_FFT(poly2, depth, root, twist, n);
            // rescale
            for (i = 0; i < size; i++)
               naive_mul_sqrt2exp(poly2->coeffs[i], poly2->coeffs[i],
                                  2*(2*n*FLINT_BITS - depth));
            // check the first few agree with input
            for (i = 0; i < length; i++)
               if (mpz_cmp(poly2->coeffs[i], poly1->coeffs[i]))
                  success = 0;
            // check the extra coefficient is correct too
            if (extra)
               if (mpz_cmp(poly2->coeffs[length], extra_coeff))
                  success = 0;
         }
         
         ZmodF_poly_clear(f);
      }
   }

   mpz_clear(extra_coeff);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
   return success;
}


// x and y should both have length 2^depth
// this version just multiplies out the convolution
void really_naive_convolution(mpz_poly_t res, mpz_poly_t x, mpz_poly_t y,
                              unsigned long depth)
{
   unsigned long size = 1UL << depth;
   mpz_poly_ensure_alloc(res, size);
   res->length = size;
   
   unsigned long i;
   for (i = 0; i < size; i++)
      mpz_set_ui(res->coeffs[i], 0);
   
   unsigned long j;
   for (i = 0; i < size; i++)
      for (j = 0; j < size; j++)
         mpz_addmul(res->coeffs[(i+j) % size], x->coeffs[i], y->coeffs[j]);
   
   for (i = 0; i < size; i++)
      mpz_mod(res->coeffs[i], res->coeffs[i], global_p);
}


// x and y should both have length 2^depth
// this version uses naive_FFT and naive_IFFT mod p
void naive_convolution(mpz_poly_t res, mpz_poly_t x, mpz_poly_t y,
                       unsigned long depth, unsigned long n)
{
   unsigned long size = 1UL << depth;
   mpz_poly_t xt, yt;
   mpz_poly_init(xt);
   mpz_poly_init(yt);
   
   mpz_poly_set(xt, x);
   mpz_poly_set(yt, y);
   
   naive_FFT(xt, depth, (4*n*FLINT_BITS) >> depth, 0, n);
   naive_FFT(yt, depth, (4*n*FLINT_BITS) >> depth, 0, n);

   mpz_poly_ensure_alloc(res, size);
   unsigned long i;
   for (i = 0; i < (1 << depth); i++)
   {
      mpz_mul(res->coeffs[i], xt->coeffs[i], yt->coeffs[i]);
      mpz_mul_2exp(res->coeffs[i], res->coeffs[i],
                   2*n*FLINT_BITS - depth);
      mpz_mod(res->coeffs[i], res->coeffs[i], global_p);
   }
   res->length = size;
   
   naive_IFFT(res, depth, (4*n*FLINT_BITS) >> depth, 0, n);
   
   mpz_poly_clear(xt);
   mpz_poly_clear(yt);
}


int test_ZmodF_poly_convolution()
{
   mpz_poly_t poly1, poly2, poly3, poly4;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init(poly3);
   mpz_poly_init(poly4);
   int success = 1;

   unsigned long depth;
   for (depth = 0; depth <= 11 && success; depth++)
   {
      unsigned long size = 1UL << depth;
   
      // need 4*n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
         ZmodF_poly_t f1, f2, f3;
         ZmodF_poly_init(f1, depth, n, 1);
         ZmodF_poly_init(f2, depth, n, 1);
         ZmodF_poly_init(f3, depth, n, 1);

#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         set_global_n(n);
         
         // switch to FFT-based convolution even for the test code, otherwise
         // tests get too slow
         int use_really_naive = (depth <= 5);
         
         unsigned long num_trials = (use_really_naive ? 50000 : 20000) /
                                    ((1 << depth) * n);
         if (num_trials == 0)
            num_trials = 1;
         
         unsigned long trial;
         for (trial = 0; trial < num_trials && success; trial++)
         {
            unsigned long len1 = random_ulong(size+1);
            unsigned long len2 = random_ulong(size+1);

            ZmodF_poly_random(f1, 4);
            ZmodF_poly_random(f2, 4);
            f1->length = len1;
            f2->length = len2;

            ZmodF_poly_convert_out(poly1, f1);
            unsigned long i;
            for (i = len1; i < size; i++)
               mpz_set_ui(poly1->coeffs[i], 0);
            ZmodF_poly_convert_out(poly2, f2);
            for (i = len2; i < size; i++)
               mpz_set_ui(poly2->coeffs[i], 0);

            ZmodF_poly_convolution(f3, f1, f2);

            unsigned long out_len = len1 + len2 - 1;
            if (out_len > size)
               out_len = size;
            
			for (i = out_len; i < size; i++)
				ZmodF_zero(f3->coeffs[i], n);

            ZmodF_poly_convert_out(poly3, f3);
            if (use_really_naive)
               really_naive_convolution(poly4, poly1, poly2, depth);
            else
               naive_convolution(poly4, poly1, poly2, depth, n);
            
            for (i = 0; i < out_len; i++)
               if (mpz_cmp(poly3->coeffs[i], poly4->coeffs[i]))
                  success = 0;
         }
         
         ZmodF_poly_clear(f3);
         ZmodF_poly_clear(f2);
         ZmodF_poly_clear(f1);
      }
   }

   mpz_poly_clear(poly4);
   mpz_poly_clear(poly3);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
   return success;
}

int test_ZmodF_poly_convolution_range()
{
   mpz_poly_t poly1, poly2, poly3, poly4;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init(poly3);
   mpz_poly_init(poly4);
   int success = 1;

   unsigned long depth;
   for (depth = 0; depth <= 11 && success; depth++)
   {
      unsigned long size = 1UL << depth;
   
      // need 4*n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
         ZmodF_poly_t f1, f2, f3;
         ZmodF_poly_init(f1, depth, n, 1);
         ZmodF_poly_init(f2, depth, n, 1);
         ZmodF_poly_init(f3, depth, n, 1);

#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         set_global_n(n);
         
         // switch to FFT-based convolution even for the test code, otherwise
         // tests get too slow
         int use_really_naive = (depth <= 5);
         
         unsigned long num_trials = (use_really_naive ? 50000 : 20000) /
                                    ((1 << depth) * n);
         if (num_trials == 0)
            num_trials = 1;
         
         unsigned long trial;
         for (trial = 0; trial < num_trials && success; trial++)
         {
            unsigned long len1 = random_ulong(size+1);
            unsigned long len2 = random_ulong(size+1);
            
            unsigned long out_len = len1 + len2 - 1;
            if (out_len > size)
               out_len = size;

            unsigned long trunc;
            if (out_len) trunc = random_ulong(out_len)+1;
            else trunc = 0;
#if DEBUG
            printf("len1 = %ld, len2 = %ld, trunc = %ld\n", len1, len2, trunc); 
#endif
            ZmodF_poly_random(f1, 4);
            ZmodF_poly_random(f2, 4);
            f1->length = len1;
            f2->length = len2;

            ZmodF_poly_convert_out(poly1, f1);
            unsigned long i;
            for (i = len1; i < size; i++)
               mpz_set_ui(poly1->coeffs[i], 0);
            ZmodF_poly_convert_out(poly2, f2);
            for (i = len2; i < size; i++)
               mpz_set_ui(poly2->coeffs[i], 0);

            ZmodF_poly_convolution_range(f3, f1, f2, 0, trunc);

            for (i = trunc; i < size; i++)
				ZmodF_zero(f3->coeffs[i], n);
            ZmodF_poly_convert_out(poly3, f3);
            
			if (use_really_naive)
               really_naive_convolution(poly4, poly1, poly2, depth);
            else
               naive_convolution(poly4, poly1, poly2, depth, n);
            
            for (i = 0; i < trunc; i++)
               if (mpz_cmp(poly3->coeffs[i], poly4->coeffs[i]))
                  success = 0;
         }
         
         ZmodF_poly_clear(f3);
         ZmodF_poly_clear(f2);
         ZmodF_poly_clear(f1);
      }
   }

   mpz_poly_clear(poly4);
   mpz_poly_clear(poly3);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
   return success;
}


// x and y should both have length 2^depth
// this version just multiplies out the convolution
void really_naive_negacyclic_convolution(mpz_poly_t res, mpz_poly_t x, mpz_poly_t y,
                                         unsigned long depth)
{
   unsigned long size = 1UL << depth;
   mpz_poly_ensure_alloc(res, size);
   res->length = size;
   
   unsigned long i;
   for (i = 0; i < size; i++)
      mpz_set_ui(res->coeffs[i], 0);
   
   unsigned long j;
   for (i = 0; i < size; i++)
      for (j = 0; j < size; j++)
      {
         unsigned long k = i + j;
         if (k < size)
            mpz_addmul(res->coeffs[k], x->coeffs[i], y->coeffs[j]);
         else
            mpz_submul(res->coeffs[k-size], x->coeffs[i], y->coeffs[j]);
      }
   
   for (i = 0; i < size; i++)
      mpz_mod(res->coeffs[i], res->coeffs[i], global_p);
}


// x and y should both have length 2^depth
// this version uses naive_FFT and naive_IFFT mod p
void naive_negacyclic_convolution(mpz_poly_t res, mpz_poly_t x, mpz_poly_t y,
                                  unsigned long depth, unsigned long n)
{
   unsigned long size = 1UL << depth;
   mpz_poly_t xt, yt;
   mpz_poly_init(xt);
   mpz_poly_init(yt);
   
   mpz_poly_set(xt, x);
   mpz_poly_set(yt, y);

   unsigned long i;
   for (i = 0; i < (1 << depth); i++)
   {
      naive_mul_sqrt2exp(xt->coeffs[i], xt->coeffs[i], (2*i*n*FLINT_BITS) >> depth);
      naive_mul_sqrt2exp(yt->coeffs[i], yt->coeffs[i], (2*i*n*FLINT_BITS) >> depth);
   }
   naive_FFT(xt, depth, (4*n*FLINT_BITS) >> depth, 0, n);
   naive_FFT(yt, depth, (4*n*FLINT_BITS) >> depth, 0, n);

   mpz_poly_ensure_alloc(res, size);
   for (i = 0; i < (1 << depth); i++)
   {
      mpz_mul(res->coeffs[i], xt->coeffs[i], yt->coeffs[i]);
      mpz_mul_2exp(res->coeffs[i], res->coeffs[i],
                   2*n*FLINT_BITS - depth);
      mpz_mod(res->coeffs[i], res->coeffs[i], global_p);
   }
   res->length = size;
   
   naive_IFFT(res, depth, (4*n*FLINT_BITS) >> depth, 0, n);
   for (i = 0; i < (1 << depth); i++)
   {
      naive_mul_sqrt2exp(res->coeffs[i], res->coeffs[i], (4*n*FLINT_BITS) - ((2*i*n*FLINT_BITS) >> depth));
   }

   mpz_poly_clear(xt);
   mpz_poly_clear(yt);
}


int test_ZmodF_poly_negacyclic_convolution()
{
   mpz_poly_t poly1, poly2, poly3, poly4;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init(poly3);
   mpz_poly_init(poly4);
   int success = 1;

   unsigned long depth;
   for (depth = 0; depth <= 10 && success; depth++)
   {
      unsigned long size = 1UL << depth;
   
      // need n*FLINT_BITS divisible by 2^depth
      unsigned long n_skip = size / (FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;
         
      unsigned long n;
      for (n = n_skip; n < 6*n_skip && success; n += n_skip)
      {
         ZmodF_poly_t f1, f2, f3;
         ZmodF_poly_init(f1, depth, n, 1);
         ZmodF_poly_init(f2, depth, n, 1);
         ZmodF_poly_init(f3, depth, n, 1);

#if DEBUG
         printf("depth = %d, n = %d\n", depth, n);
#endif

         set_global_n(n);
         
         // switch to FFT-based convolution even for the test code, otherwise
         // tests get too slow
         int use_really_naive = (depth <= 5);
         
         unsigned long num_trials = (use_really_naive ? 50000 : 20000) /
                                    ((1 << depth) * n);
         if (num_trials == 0)
            num_trials = 1;
         
         unsigned long trial;
         for (trial = 0; trial < num_trials && success; trial++)
         {
            ZmodF_poly_random(f1, 4);
            ZmodF_poly_random(f2, 4);
            f1->length = size;
            f2->length = size;

            ZmodF_poly_convert_out(poly1, f1);
            ZmodF_poly_convert_out(poly2, f2);

            ZmodF_poly_negacyclic_convolution(f3, f1, f2);

            ZmodF_poly_convert_out(poly3, f3);
            if (use_really_naive)
               really_naive_negacyclic_convolution(poly4, poly1, poly2, depth);
            else
               naive_negacyclic_convolution(poly4, poly1, poly2, depth, n);
            
            unsigned long i;
            for (i = 0; i < size; i++)
               if (mpz_cmp(poly3->coeffs[i], poly4->coeffs[i]))
                  success = 0;
         }
         
         ZmodF_poly_clear(f3);
         ZmodF_poly_clear(f2);
         ZmodF_poly_clear(f1);
      }
   }

   mpz_poly_clear(poly4);
   mpz_poly_clear(poly3);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
   return success;
}


/****************************************************************************

   Main test functions

****************************************************************************/

void ZmodF_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(_ZmodF_poly_FFT_iterative);
   RUN_TEST(_ZmodF_poly_FFT_factor);
   RUN_TEST(_ZmodF_poly_IFFT_recursive);
   RUN_TEST(_ZmodF_poly_IFFT_iterative);
   RUN_TEST(_ZmodF_poly_IFFT);
   RUN_TEST(ZmodF_poly_convolution);
   RUN_TEST(ZmodF_poly_convolution_range);
   RUN_TEST(ZmodF_poly_negacyclic_convolution);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   test_support_init();
   mpz_init(global_p);
   
   ZmodF_poly_test_all();

   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}



// end of file ****************************************************************
