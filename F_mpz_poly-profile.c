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

F_mpz_poly-profile.c

Profiling for F_mpz_poly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include <math.h>
#include <gmp.h>
#include "profiler-main.h"
#include "flint.h"
#include "memory-manager.h"
#include "F_mpz_poly.h"
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
   
   unsigned long i;
   for (i = 0; i < length; i++)
   {
       bits = maxbits;
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_urandomb(temp, randstate, bits);
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
      unsigned long length = (unsigned long) floor(pow(ratio, i));
      if (length != last_length)
      {
         last_length = length;

         unsigned long last_bits = 0;
         unsigned long j;
         for (j = 0; j <= max_iter; j++)
         {
            unsigned long bits = (unsigned long) floor(pow(ratio, j));
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

/*
   Calls prof2d_sample(length, bits, NULL) for all length, bits combinations
   such that the output coefficients of the polynomials being multiplied fit
	into FLINT_BITS-2 bits.
 */

void run_block(ulong max_length, double ratio)
{
	int max_iter = (int) ceil(log((double) max_length) / log(ratio));
	ulong last_length = 0;
	
	ulong i;
	for (i = 0; i <= max_iter; i++)
	{
		ulong length = (ulong) floor(pow(ratio, i));
      if (length != last_length)
      {
         last_length = length;
         ulong m = ceil_log2(length);
	         
         ulong bits;
         for (bits = 1; bits <= 31; bits++)
	      {
		      if (m + 2*bits <= FLINT_BITS - 2) prof2d_sample(length, bits, NULL);
	      }
		}
	}
}

// ============================================================================


void sample_F_mpz_poly_mul_KS(ulong length, ulong bits,
                             void* arg, ulong count)
{
   ulong m = ceil_log2(length);
   ulong output_bits = 2*bits+m;
   
   F_mpz_poly_t poly1, poly2, poly3;
   mpz_poly_t r_poly, r_poly2;  
   
   mpz_poly_init(r_poly); 
   mpz_poly_init(r_poly2); 
   mpz_poly_realloc(r_poly, length);
   mpz_poly_realloc(r_poly2, length);
  
   F_mpz_poly_init2(poly1, length);
   F_mpz_poly_init2(poly2, length);
   F_mpz_poly_init2(poly3, 2*length-1);
   
   ulong r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   ulong i;
   for (i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
         randpoly(r_poly, length, bits);
         mpz_poly_to_F_mpz_poly(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         mpz_poly_to_F_mpz_poly(poly2, r_poly2);
      }
      prof_start();
      F_mpz_poly_mul_KS(poly3, poly1, poly2);
      prof_stop();
   }
   
   mpz_poly_clear(r_poly);
   mpz_poly_clear(r_poly2);
   
   F_mpz_poly_clear(poly3);
   F_mpz_poly_clear(poly2);
   F_mpz_poly_clear(poly1);
   
}


char* profDriverString_F_mpz_poly_mul_KS(char* params)
{
   return "F_mpz_poly_mul_KS over various lengths and various bit sizes.\n"
   "Parameters are: max length; ratio between consecutive lengths.";
}

char* profDriverDefaultParams_F_mpz_poly_mul_KS()
{
   return "4000000 1.2";
}


void profDriver_F_mpz_poly_mul_KS(char* params)
{
   ulong max_length;
   double ratio;

   sscanf(params, "%ld %lf", &max_length, &ratio);

   test_support_init();
   prof2d_set_sampler(sample_F_mpz_poly_mul_KS);
   run_triangle(max_length, ratio);
   test_support_cleanup();
}

// ============================================================================


void sample_F_mpz_poly_mul(ulong length, ulong bits,
                             void* arg, ulong count)
{
   ulong m = ceil_log2(length);
   ulong output_bits = 2*bits+m;
   
   F_mpz_poly_t poly1, poly2, poly3;
   mpz_poly_t r_poly, r_poly2;  
   
   mpz_poly_init(r_poly); 
   mpz_poly_init(r_poly2); 
   mpz_poly_realloc(r_poly, length);
   mpz_poly_realloc(r_poly2, length);
  
   F_mpz_poly_init2(poly1, length);
   F_mpz_poly_init2(poly2, length);
   F_mpz_poly_init2(poly3, 2*length-1);
   
   ulong r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
   
   ulong i;
   for (i = 0; i < count; i++)
   {
      if (i%r_count == 0)
      {
         randpoly(r_poly, length, bits);
         mpz_poly_to_F_mpz_poly(poly1, r_poly);
         randpoly(r_poly2, length, bits);
         mpz_poly_to_F_mpz_poly(poly2, r_poly2);
      }
      prof_start();
      F_mpz_poly_mul(poly3, poly1, poly2);
      prof_stop();
   }
   
   mpz_poly_clear(r_poly);
   mpz_poly_clear(r_poly2);
   
   F_mpz_poly_clear(poly3);
   F_mpz_poly_clear(poly2);
   F_mpz_poly_clear(poly1);
   
}


char* profDriverString_F_mpz_poly_mul(char* params)
{
   return "F_mpz_poly_mul over various lengths and various bit sizes.\n"
   "Parameters are: max length; ratio between consecutive lengths.";
}

char* profDriverDefaultParams_F_mpz_poly_mul()
{
   return "4000000 1.2";
}


void profDriver_F_mpz_poly_mul(char* params)
{
   ulong max_length;
   double ratio;

   sscanf(params, "%ld %lf", &max_length, &ratio);

   test_support_init();
   prof2d_set_sampler(sample_F_mpz_poly_mul);
   run_triangle(max_length, ratio);
   test_support_cleanup();
}

// ============================================================================

// end of file ****************************************************************
