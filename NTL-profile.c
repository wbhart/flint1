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
#include <NTL/ZZXFactoring.h>
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
   unsigned long i;
   for (i = 0; i <= max_iter; i++)
   {
      unsigned long length = (unsigned long) floor(powl(ratio, i));
      if (length != last_length)
      {
         last_length = length;

         unsigned long last_bits = 0;
         unsigned long j;
         for (j = 0; j <= max_iter; j++)
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

void sample_NTL_factor(unsigned long length, unsigned long bits,
                          void* arg, unsigned long count)
{
    ZZX poly1, poly2, poly3;
    ZZ a, c;
    vec_pair_ZZX_long factors;
    
    poly1.SetMaxLength(length);
    //poly2.SetMaxLength(length);
    //poly3.SetMaxLength(2*length-1);


   
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   unsigned long i;
   for (i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
	    unsigned long j;
	    for (j = 0; j<length; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		//RandomBits(a,bits);
		//SetCoeff(poly2,j,a);
		}
		//mul(poly3, poly1, poly2);
      }
       prof_start();
       factor(c, factors, poly1);
       prof_stop();
   }
   

}




char* profDriverString_NTL_factor(char* params)
{
   return "NTL_factor factors polynomials of various lengths and various bit sizes.\n"
   "Parameters are: max bitsize; ratio between consecutive lengths/bitsizes.";
}
 
char* profDriverDefaultParams_NTL_factor()
{
   return "10000 1.2";
}
 
void profDriver_NTL_factor(char* params)
{
   unsigned long max_bits;
   double ratio;
    
   sscanf(params, "%ld %lf", &max_bits, &ratio);
 
   test_support_init();
   prof2d_set_sampler(sample_NTL_factor);
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
    poly1.SetMaxLength(length);
    poly2.SetMaxLength(length);
    poly3.SetMaxLength(2*length-1);

   
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   unsigned long i;
   for (i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
	    unsigned long j;
	    for (j = 0; j<length; j++)
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
    poly3.SetMaxLength(2*length-1);


   
   
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 10000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 4;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   unsigned long i;
   for (i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
	    do
        {
           unsigned long j;
           for (j = 0; j < length; j++)
		   {
		      RandomBits(a,bits);
		      SetCoeff(poly1,j,a);
		   }
        } while (IsZero(poly1));
        unsigned long j;
        for (j = 0; j < length; j++)
		{
           RandomBits(a,bits);
		   SetCoeff(poly2,j,a);
		}
      }
      
      mul(poly3, poly1, poly2);
      prof_start();
      unsigned long count2;
      for (count2 = 0; count2 < r_count; count2++)
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
   
   unsigned long i;
   for (i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
	    unsigned long j;
	    for (j = 0; j<length-1; j++)
		{
		RandomBits(a,bits);
		SetCoeff(poly1,j,a);
		}
		SetCoeff(poly1,length-1,1);
	    unsigned long j;
	    for (j = 0; j<2*length-1; j++)
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


