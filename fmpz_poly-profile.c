/****************************************************************************

fmpz_poly-profile.c

Profiling for ZmodFpoly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include <math.h>
#include <gmp.h>
#include "profiler-main.h"
#include "flint.h"
#include "memory-manager.h"
#include "fmpz_poly.h"
#include "Zpoly.h"

//================================================================================

#define SIGNS 0

gmp_randstate_t Zpoly_test_randstate;

unsigned long randint(unsigned long randsup) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(Zpoly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   //if (pol->coeffs) Zpoly_clear(pol);
   //Zpoly_init3(pol, length, maxbits);
   for (unsigned long i = 0; i < length; i++)
   {
       bits = maxbits;
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, Zpoly_test_randstate, bits);
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
       }
       Zpoly_set_coeff(pol, i, temp);
       
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
   Zpoly_t r_poly, r_poly2;  
   
   Zpoly_init(r_poly); 
   Zpoly_init(r_poly2); 
   Zpoly_realloc(r_poly, length);
   Zpoly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, (output_bits-1)/FLINT_BITS_PER_LIMB+1);
    
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
         _fmpz_poly_convert_in(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         _fmpz_poly_convert_in(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul_KS(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   Zpoly_clear(r_poly);
   Zpoly_clear(r_poly2);
   
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
      sscanf(params, "%d %d %lf %d", &length_min, &length_max,
                                     &ratio, &bits_min);
   }

   gmp_randinit_default(Zpoly_test_randstate);

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
}

// **************************************************************************************

void sample_fmpz_poly_mul_SS(unsigned long length, unsigned long bits,
                          unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   Zpoly_t r_poly, r_poly2;  
   
   Zpoly_init(r_poly); 
   Zpoly_init(r_poly2); 
   Zpoly_realloc(r_poly, length);
   Zpoly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, (output_bits-1)/FLINT_BITS_PER_LIMB+1);
    
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
         _fmpz_poly_convert_in(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         _fmpz_poly_convert_in(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul_SS(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   Zpoly_clear(r_poly);
   Zpoly_clear(r_poly2);
   
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
      sscanf(params, "%d %d %lf %d", &length_min, &length_max,
                                     &ratio, &bits_min);
   }

   gmp_randinit_default(Zpoly_test_randstate);

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
}

// ============================================================================


void sample_fmpz_poly_mul_karatsuba(unsigned long length, unsigned long bits,
                          unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   Zpoly_t r_poly, r_poly2;  
   
   Zpoly_init(r_poly); 
   Zpoly_init(r_poly2); 
   Zpoly_realloc(r_poly, length);
   Zpoly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly3, 2*length-1, poly1->limbs+poly2->limbs+2);
    
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
         _fmpz_poly_convert_in(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         _fmpz_poly_convert_in(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul_karatsuba(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   Zpoly_clear(r_poly);
   Zpoly_clear(r_poly2);
   
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
      sscanf(params, "%d %d %lf %d", &length_min, &length_max,
                                     &ratio, &bits_min);
   }

   gmp_randinit_default(Zpoly_test_randstate);

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
}


// ============================================================================


void sample_fmpz_poly_mul(unsigned long length, unsigned long bits,
                          unsigned long count)
{
   unsigned long m = ceil_log2(length);
   unsigned long output_bits = 2*bits+m;
   
   fmpz_poly_t poly1, poly2, poly3;
   Zpoly_t r_poly, r_poly2;  
   
   Zpoly_init(r_poly); 
   Zpoly_init(r_poly2); 
   Zpoly_realloc(r_poly, length);
   Zpoly_realloc(r_poly2, length);
  
   _fmpz_poly_stack_init(poly1, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   _fmpz_poly_stack_init(poly2, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
   if ((bits > 256) && (length <=16) && (length > 2)) _fmpz_poly_stack_init(poly3, 2*length-1, poly1->limbs+poly2->limbs+2);
   else  _fmpz_poly_stack_init(poly3, 2*length-1, (output_bits-1)/FLINT_BITS_PER_LIMB+1);
   
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
         _fmpz_poly_convert_in(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         _fmpz_poly_convert_in(poly2, r_poly2);
      }
       prof2d_start();
       _fmpz_poly_mul(poly3, poly1, poly2);
       prof2d_stop();
   }
   
   Zpoly_clear(r_poly);
   Zpoly_clear(r_poly2);
   
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
      sscanf(params, "%d %d %lf %d", &length_min, &length_max,
                                     &ratio, &bits_min);
   }

   gmp_randinit_default(Zpoly_test_randstate);

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
}



// end of file ****************************************************************
