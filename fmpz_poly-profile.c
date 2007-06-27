/****************************************************************************

fmpz_poly-profile.c

Profiling for fmpz_poly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include <math.h>
#include <gmp.h>
#include "profiler-main.h"
#include "flint.h"
#include "memory-manager.h"
#include "fmpz_poly.h"
#include "mpz_poly.h"
#include "test-support.h"

//=============================================================================

// whether to generate signed or unsigned random polys
#define SIGNS 0


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
       bits = maxbits;
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

// ============================================================================


/*
Calls prof2d_sample(length, bits) for all length, bits combinations
such that length*bits < max_bits, with length and bits spaced out by
the given ratio
 */
void run_triangle(unsigned long max_bits, double ratio)
{
   int max_iter = (int) ceil(log((double) max_bits) / log(ratio));

   unsigned long last_length = 0;
   for (unsigned long i = 0; i <= max_iter; i++)
   {
      unsigned long length = (unsigned long) floor(pow(ratio, i));
      if (length != last_length)
      {
         last_length = length;

         unsigned long last_bits = 0;
         for (unsigned long j = 0; j <= max_iter; j++)
         {
            unsigned long bits = (unsigned long) floor(pow(ratio, j));
            if (bits != last_bits)
            {
               last_bits = bits;

               if (bits * length < max_bits)
                  prof2d_sample(length, bits);
            }
         }
      }
   }
}


// ============================================================================


void sample_fmpz_poly_mul_KS(unsigned long length, unsigned long bits,
                             unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   mpz_poly_t r_poly, r_poly2;  
   
   mpz_poly_init(r_poly); 
   mpz_poly_init(r_poly2); 
   mpz_poly_realloc(r_poly, length);
   mpz_poly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, (output_bits-1)/FLINT_BITS+1);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   for (unsigned long i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
         randpoly(r_poly, length, bits);
         mpz_poly_to_fmpz_poly(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         mpz_poly_to_fmpz_poly(poly2, r_poly2);
      }
      prof2d_start();
      _fmpz_poly_mul_KS(poly3, poly1, poly2);
      prof2d_stop();
   }
   
   mpz_poly_clear(r_poly);
   mpz_poly_clear(r_poly2);
   
   _fmpz_poly_stack_clear(poly3);
   _fmpz_poly_stack_clear(poly2);
   _fmpz_poly_stack_clear(poly1);
   
}


char* prof2dDriverString_fmpz_poly_mul_KS(char* params)
{
   return "fmpz_poly_mul_KS over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}

char* prof2dDriverDefaultParams_fmpz_poly_mul_KS()
{
   return "16000000 1.2";
}


void prof2dDriver_fmpz_poly_mul_KS(char* params)
{
   int max_bits;
   double ratio;

   sscanf(params, "%ld %lf", &max_bits, &ratio);

   test_support_init();
   prof2d_set_sampler(sample_fmpz_poly_mul_KS);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// ****************************************************************************

void sample_fmpz_poly_mul_SS(unsigned long length, unsigned long bits,
                             unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   mpz_poly_t r_poly, r_poly2;  
   
   mpz_poly_init(r_poly); 
   mpz_poly_init(r_poly2); 
   mpz_poly_realloc(r_poly, length);
   mpz_poly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, (output_bits-1)/FLINT_BITS+1);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
 
   
   for (unsigned long i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
         randpoly(r_poly, length, bits);
         mpz_poly_to_fmpz_poly(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         mpz_poly_to_fmpz_poly(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul_SS(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   mpz_poly_clear(r_poly);
   mpz_poly_clear(r_poly2);
   
   _fmpz_poly_stack_clear(poly3);
   _fmpz_poly_stack_clear(poly2);
   _fmpz_poly_stack_clear(poly1);
   
}


char* prof2dDriverString_fmpz_poly_mul_SS(char* params)
{
   return "fmpz_poly_mul_SS over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}

char* prof2dDriverDefaultParams_fmpz_poly_mul_SS()
{
   return "16000000 1.2";
}

void prof2dDriver_fmpz_poly_mul_SS(char* params)
{
   int max_bits;
   double ratio;

   sscanf(params, "%ld %lf", &max_bits, &ratio);
   
   test_support_init();
   prof2d_set_sampler(sample_fmpz_poly_mul_SS);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// ============================================================================


void sample_fmpz_poly_mul_karatsuba(unsigned long length, unsigned long bits,
                                    unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   mpz_poly_t r_poly, r_poly2;  
   
   mpz_poly_init(r_poly); 
   mpz_poly_init(r_poly2); 
   mpz_poly_realloc(r_poly, length);
   mpz_poly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, poly1->limbs+poly2->limbs+1);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
 
   
   for (unsigned long i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
         randpoly(r_poly, length, bits);
         mpz_poly_to_fmpz_poly(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         mpz_poly_to_fmpz_poly(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul_karatsuba(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   mpz_poly_clear(r_poly);
   mpz_poly_clear(r_poly2);
   
   _fmpz_poly_stack_clear(poly3);
   _fmpz_poly_stack_clear(poly2);
   _fmpz_poly_stack_clear(poly1);
   
}


char* prof2dDriverString_fmpz_poly_mul_karatsuba(char* params)
{
   return "fmpz_poly_mul_karatsuba over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}

char* prof2dDriverDefaultParams_fmpz_poly_mul_karatsuba()
{
   return "300000 1.2";
}

void prof2dDriver_fmpz_poly_mul_karatsuba(char* params)
{
   int max_bits;
   double ratio;

   sscanf(params, "%ld %lf", &max_bits, &ratio);

   test_support_init();
   prof2d_set_sampler(sample_fmpz_poly_mul_karatsuba);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// ============================================================================


void sample_fmpz_poly_mul(unsigned long length, unsigned long bits,
                          unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   mpz_poly_t r_poly, r_poly2;  
   
   mpz_poly_init(r_poly); 
   mpz_poly_init(r_poly2); 
   mpz_poly_realloc(r_poly, length);
   mpz_poly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, (output_bits-1)/FLINT_BITS+1);
   
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
 
   
   for (unsigned long i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
         randpoly(r_poly, length, bits);
         mpz_poly_to_fmpz_poly(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         mpz_poly_to_fmpz_poly(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   mpz_poly_clear(r_poly);
   mpz_poly_clear(r_poly2);
   
   _fmpz_poly_stack_clear(poly3);
   _fmpz_poly_stack_clear(poly2);
   _fmpz_poly_stack_clear(poly1);
   
}


char* prof2dDriverString_fmpz_poly_mul(char* params)
{
   return "fmpz_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}

char* prof2dDriverDefaultParams_fmpz_poly_mul()
{
   return "16000000 1.2";
}

void prof2dDriver_fmpz_poly_mul(char* params)
{
   int max_bits;
   double ratio;

   sscanf(params, "%ld %lf", &max_bits, &ratio);

   test_support_init();
   prof2d_set_sampler(sample_fmpz_poly_mul);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// end of file ****************************************************************
