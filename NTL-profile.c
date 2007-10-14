/****************************************************************************

NTL-profile.c

Profiling for NTL

Copyright (C) 2007, Tomasz Lechowski, William Hart and David Harvey

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

#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <stdio.h>


//=============================================================================

NTL_CLIENT

// whether to generate signed or unsigned random polys
#define SIGNS 0


unsigned long randint(unsigned long randsup) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

// ============================================================================


/*
Calls prof2d_sample(length, bits, NULL) for all length, bits combinations
such that length*bits < max_bits, with length and bits spaced out by
the given ratio
*/
 
void run_triangle(unsigned long max_bits, double ratio)
{
   int max_iter = (int) ceil(log((double) max_bits) / log(ratio));

   unsigned long last_length = 0;
   for (unsigned long i = 0; i <= max_iter; i++)
   {
      unsigned long length = (unsigned long) floor(powl(ratio, i));
      if (length != last_length)
      {
         last_length = length;

         unsigned long last_bits = 0;
         for (unsigned long j = 0; j <= max_iter; j++)
         {
            unsigned long bits = (unsigned long) floor(powl(ratio, j));
            if (bits != last_bits)
            {
               last_bits = bits;

               if (bits * length < max_bits)
                  prof2d_sample(length, bits, NULL);
            }
         }
      }
   }
}


// ============================================================================

void sample_NTL_Kar_poly_mul(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{

       
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length+1);
    poly2.SetMaxLength(length+1);
    poly3.SetMaxLength((length)+(length)+1);


   
   
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
	    for (unsigned long j = 0; j<length+1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		RandomBits(a,bits);
		SetCoeff(poly2,j,a);
		}
      }
       prof_start();
       KarMul(poly3, poly1, poly2);
       prof_stop();
   }
   
   
}




char* profDriverString_NTL_Kar_poly_mul(char* params)
{
   return "NTL_Kar_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_Kar_poly_mul()
{
   return "16000000 1.2";
}
 
void profDriver_NTL_Kar_poly_mul(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_Kar_poly_mul);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}



// ============================================================================


void sample_NTL_Hom_poly_mul(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length+1);
    poly2.SetMaxLength(length+1);
    poly3.SetMaxLength((length)+(length)+1);


   
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
	    for (unsigned long j = 0; j<length+1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		RandomBits(a,bits);
		SetCoeff(poly2,j,a);
		}
      }
       prof_start();
       HomMul(poly3, poly1, poly2);
       prof_stop();
   }
   

}




char* profDriverString_NTL_Hom_poly_mul(char* params)
{
   return "NTL_Hom_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_Hom_poly_mul()
{
   return "16000000 1.2";
}
 
void profDriver_NTL_Hom_poly_mul(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_Hom_poly_mul);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}
// ============================================================================

void sample_NTL_SS_poly_mul(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{
       
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length+1);
    poly2.SetMaxLength(length+1);
    poly3.SetMaxLength((length)+(length)+1);

   
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
	    for (unsigned long j = 0; j<length+1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		RandomBits(a,bits);
		SetCoeff(poly2,j,a);
		}
      }
       prof_start();
       SSMul(poly3, poly1, poly2);
       prof_stop();
   }
   

   
}




char* profDriverString_NTL_SS_poly_mul(char* params)
{
   return "NTL_SS_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_SS_poly_mul()
{
   return "16000000 1.2";
}
 
void profDriver_NTL_SS_poly_mul(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_SS_poly_mul);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// ============================================================================

void sample_NTL_Plain_poly_mul(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{
       
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length+1);
    poly2.SetMaxLength(length+1);
    poly3.SetMaxLength((length)+(length)+1);

   
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
	    for (unsigned long j = 0; j<length+1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		RandomBits(a,bits);
		SetCoeff(poly2,j,a);
		}
      }
       prof_start();
       PlainMul(poly3, poly1, poly2);
       prof_stop();
   }
   

   
}




char* profDriverString_NTL_Plain_poly_mul(char* params)
{
   return "NTL_Plain_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_Plain_poly_mul()
{
   return "16000000 1.2";
}
 
void profDriver_NTL_Plain_poly_mul(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_Plain_poly_mul);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// ============================================================================



void sample_NTL_poly_mul(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{
       
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length+1);
    poly2.SetMaxLength(length+1);
    poly3.SetMaxLength((length)+(length)+1);

   
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
	    for (unsigned long j = 0; j<length+1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		RandomBits(a,bits);
		SetCoeff(poly2,j,a);
		}
      }
       prof_start();
       mul(poly3, poly1, poly2);
       prof_stop();
   }
   

   
}




char* profDriverString_NTL_poly_mul(char* params)
{
   return "NTL_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_poly_mul()
{
   return "16000000 1.2";
}
 
void profDriver_NTL_poly_mul(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_poly_mul);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}


// ============================================================================

void sample_NTL_poly_div1(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{

       
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length);
    poly2.SetMaxLength(length);
    poly3.SetMaxLength((length)+(length)-1);


   
   
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 10000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 4;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   for (unsigned long i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
	    do
        {
           for (unsigned long j = 0; j < length; j++)
		   {
		      RandomBits(a,bits);
		      SetCoeff(poly1,j,a);
		   }
        } while (IsZero(poly1));
        for (unsigned long j = 0; j < length; j++)
		{
           RandomBits(a,bits);
		   SetCoeff(poly2,j,a);
		}
      }
      
      mul(poly3, poly1, poly2);
      prof_start();
      for (unsigned long count2 = 0; count2 < r_count; count2++)
      {
         divide(poly2, poly3, poly1);
      }
      prof_stop();
      
      i += (r_count-1);
   }  
}




char* profDriverString_NTL_poly_div1(char* params)
{
   return "NTL_poly_div1 over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_poly_div1()
{
   return "1000000 1.2";
}
 
void profDriver_NTL_poly_div1(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_poly_div1);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}



// ============================================================================



void sample_NTL_poly_div2(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{

       
    ZZX poly1;
    ZZX poly2;
    ZZX poly3;
    ZZ a;
    poly1.SetMaxLength(length);
    poly2.SetMaxLength(length);
    poly3.SetMaxLength(2*length-1);


   
   
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
	    for (unsigned long j = 0; j<length-1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		}
		SetCoeff(poly1,length-1,1);
	    for (unsigned long j = 0; j<2*length-1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly3,j,a);
		}
      }
       prof_start();
       div(poly2, poly3, poly1);
       prof_stop();
   }
   
   
}




char* profDriverString_NTL_poly_div2(char* params)
{
   return "NTL_poly_div2 over various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_poly_div2()
{
   return "1100000 1.77827941";
}
 
void profDriver_NTL_poly_div2(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_poly_div2);
   run_triangle(max_bits, ratio);
   test_support_cleanup();
}

// end of file ****************************************************************


