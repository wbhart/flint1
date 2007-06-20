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

//================================================================================

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
    
   unsigned long r_count;
   
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
   return "fmpz_poly_mul_KS over various lengths and various bit sizes";
}


/*
Parameters for this target are:
   length_min: minimum length
   length_max: maximum length
   ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_fmpz_poly_mul_KS(char* params)
{
   int length_min, length_max, bits_min;
   double ratio;

   if (strlen(params) == 0)
   {
      // default parameters:
      length_min = 1;
      length_max = 1000000;
      bits_min = 1;
      
      ratio = 1.2;
   }
   else
   {
      sscanf(params, "%ld %ld %lf %ld", &length_min, &length_max,
                                        &ratio, &bits_min);
   }

   test_support_init();

   prof2d_set_sampler(sample_fmpz_poly_mul_KS);

   for (unsigned long length = length_min; length < length_max;
        length = (int)(ceil(ratio * (float)length)))
   {
      for (unsigned long bits = bits_min; (bits <= length); 
                                  bits = (int)(ceil(ratio * bits)))
      {
         if (bits == 444) bits = 445;
         if (length == 444) length = 445;
         if (bits == 924) bits = 925;
         if (length == 924) length = 925;
         if (bits == 51144) bits = 51145;
         if (length == 51144) length = 51145;
         if (bits == 61374) bits = 61375;
         if (length == 61374) length = 61375;
         if (bits == 106056) bits = 106057;
         if (length == 106056) length = 106057;
         if (bits*length > 1000000) continue;
         prof2d_sample(length, bits);
      }
   }

   test_support_cleanup();
}

// **************************************************************************************

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
    
   unsigned long r_count;
   
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
   return "fmpz_poly_mul_SS over various lengths and various bit sizes";
}


/*
Parameters for this target are:
   length_min: minimum length
   length_max: maximum length
   ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_fmpz_poly_mul_SS(char* params)
{
   int length_min, length_max, bits_min, bits_max;
   double ratio;

   if (strlen(params) == 0)
   {
      // default parameters:
      length_min = 1;
      length_max = 1000000;
      bits_max = 1000000;
      bits_min = 1;
      
      ratio = 1.2;
   }
   else
   {
      sscanf(params, "%ld %ld %lf %ld", &length_min, &length_max,
                                        &ratio, &bits_min);
   }

   test_support_init();

   prof2d_set_sampler(sample_fmpz_poly_mul_SS);

   for (unsigned long length = length_min; length <= length_max;
        length = (int)(ceil(ratio * (float)length)))
   {
      for (unsigned long bits = bits_min; bits <= bits_max; 
                                  bits = (int)(ceil(ratio * bits)))
      {
         if (bits * length > 1000000) continue;
         
         if (bits == 444) bits = 445;
         if (length == 444) length = 445;
         if (bits == 924) bits = 925;
         if (length == 924) length = 925;
         if (bits == 51144) bits = 51145;
         if (length == 51144) length = 51145;
         if (bits == 61374) bits = 61375;
         if (length == 61374) length = 61375;
         if (bits == 106056) bits = 106057;
         if (length == 106056) length = 106057;
         if (bits*length > 1000000) continue;
         prof2d_sample(length, bits);
      }
   }

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
    
   unsigned long r_count;
   
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
   return "fmpz_poly_mul_karatsuba over various lengths and various bit sizes";
}


/*
Parameters for this target are:
   length_min: minimum length
   length_max: maximum length
   ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_fmpz_poly_mul_karatsuba(char* params)
{
   int length_min, length_max, bits_min, bits_max;
   double ratio;

   if (strlen(params) == 0)
   {
      // default parameters:
      length_min = 1;
      length_max = 1000000;
      bits_max = 1000000;
      bits_min = 1;
      
      ratio = 1.2;
   }
   else
   {
      sscanf(params, "%ld %ld %lf %ld", &length_min, &length_max,
                                        &ratio, &bits_min);
   }

   test_support_init();

   prof2d_set_sampler(sample_fmpz_poly_mul_karatsuba);

   for (unsigned long length = length_min; length <= length_max;
        length = (int)(ceil(ratio * (float)length)))
   {
      for (unsigned long bits = bits_min; bits <= bits_max; 
                                  bits = (int)(ceil(ratio * bits)))
      {
         if (bits * length > 1000000) continue;
         
         if (bits == 444) bits = 445;
         if (length == 444) length = 445;
         if (bits == 924) bits = 925;
         if (length == 924) length = 925;
         if (bits == 51144) bits = 51145;
         if (length == 51144) length = 51145;
         if (bits == 61374) bits = 61375;
         if (length == 61374) length = 61375;
         if (bits == 106056) bits = 106057;
         if (length == 106056) length = 106057;
         if (bits*length > 1000000) continue;
         prof2d_sample(length, bits);
      }
   }

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
   
   unsigned long r_count;
   
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
   return "fmpz_poly_mul over various lengths and various bit sizes";
}


/*
Parameters for this target are:
   length_min: minimum length
   length_max: maximum length
   ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_fmpz_poly_mul(char* params)
{
   int length_min, length_max, bits_min, bits_max;
   double ratio;

   if (strlen(params) == 0)
   {
      // default parameters:
      length_min = 1;
      length_max = 16000000;
      bits_max = 16000000;
      bits_min = 1;
      
      ratio = 1.2;
   }
   else
   {
      sscanf(params, "%ld %ld %lf %ld", &length_min, &length_max,
                                        &ratio, &bits_min);
   }

   test_support_init();

   prof2d_set_sampler(sample_fmpz_poly_mul);

   for (unsigned long length = length_min; length <= length_max;
        length = (int)(ceil(ratio * (float)length)))
   {
      for (unsigned long bits = bits_min; bits <= bits_max; 
                                  bits = (int)(ceil(ratio * bits)))
      {
         if (bits * length > 16000000) continue;
         
         if (bits == 444) bits = 445;
         if (length == 444) length = 445;
         if (bits == 924) bits = 925;
         if (length == 924) length = 925;
         if (bits == 51144) bits = 51145;
         if (length == 51144) length = 51145;
         if (bits == 61374) bits = 61375;
         if (length == 61374) length = 61375;
         if (bits == 106056) bits = 106057;
         if (length == 106056) length = 106057;
         if (bits*length > 1000000) continue;
         prof2d_sample(length, bits);
      }
   }

   test_support_cleanup();
}



// end of file ****************************************************************
