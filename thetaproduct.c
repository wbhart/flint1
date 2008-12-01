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

/*
   Demo FLINT program for computing products of theta functions.
   
   (C) 2008 William Hart
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "theta.h"

#ifdef HAVE_ZNPOLY
#include "zn_poly.h"
#endif

#define MIKE_LIMIT 100000000L
#define MIKE_LIMIT2 MIKE_LIMIT/4
#define BLOCK 100000
#define COUNT MIKE_LIMIT/BLOCK
#define K 1

/*int main(void)
{
	long * array1 = (long *) flint_heap_alloc(1000);
	long quad[4] = {1, -1, 1, -1};
	
	theta_1d_quadchar(quad, 3, 1, 0, array1, 1000, 1000);

	theta_print(array1, 1000, 1000);

	return 0;
}*/

/*int main(void)
{
   unsigned long p = 65521;

   zn_mod_t mod;
   zn_mod_init(mod, p);
   
   //------------------------------------------------------------
	
   unsigned long * theta_A1 = (unsigned long *) flint_heap_alloc(MIKE_LIMIT2);
   
   long * array1 = (long *) flint_heap_alloc(BLOCK);
   
   for (unsigned long start = 0; start < MIKE_LIMIT; start += BLOCK)
   {   
	  theta_2d_A1(array1, start, BLOCK);
   
      unsigned long start2 = start/8;
	  
	  for (unsigned long i = 0; i < BLOCK/8; i++)
      {
         theta_A1[start2 + i] = array1[8*i+K];
      }
   }
   
   flint_heap_free(array1);

   //-------------------------------------------------------------
   
   unsigned long * theta_B = (unsigned long *) flint_heap_alloc(MIKE_LIMIT2);
   
   long * array2 = (long *) flint_heap_alloc(BLOCK);
   
   for (unsigned long start = 0; start < MIKE_LIMIT; start += BLOCK)
   {   
	  theta_1d_B(array2, start, BLOCK);
   
      unsigned long start2 = start/8;
	  
	  long coeff;
	  
	  for (unsigned long i = 0; i < BLOCK/8; i++)
      {
         coeff = array2[8*i];
		 
		 if (coeff < 0L) theta_B[start2 + i] = p + coeff;
		 else theta_B[start2 + i] = coeff;
      }
   }
   
   flint_heap_free(array2);
   
   //----------------------------------------------------------------------
   unsigned long * theta_prod = flint_heap_alloc(2*MIKE_LIMIT2-1);

   printf("here1\n");
   zn_array_mul_fft_dft(theta_prod, theta_A1, MIKE_LIMIT2, theta_B, MIKE_LIMIT2, 4, mod);

   flint_heap_free(theta_A1);
   flint_heap_free(theta_B);
  
   printf("here2\n");
   
   unsigned long s = 0;
   
   long coeff;
   for(unsigned long j = 0; j < MIKE_LIMIT2; j++)
   {
       coeff = theta_prod[j];
	   if (!coeff) 
	   {
		  //printf("%ld ", 8*j+K);
		  s++;
	   }
   }
   printf("\n%ld zeroes\n", s);
   
   flint_heap_free(theta_prod);
   
   
   return 0;
}*/

int main(void)
{
   fmpz_poly_t theta_1, theta_C, theta_prod;
   
   //--------------------------------------------------------------
   
   long * array1 = (long *) flint_heap_alloc(BLOCK);
   
   long character[4] = {1, -1, 1, -1};
	
   fmpz_poly_init2(theta_1, MIKE_LIMIT2, 1); 
   
   for (unsigned long start = 0; start < MIKE_LIMIT; start += BLOCK)
   {   
	  theta_1d_quadchar_0(character, 4, 4, 1, 4, array1, start, BLOCK);
   
      unsigned long start2 = start/4;
	  
	  for (unsigned long i = 0; i < BLOCK/4; i++)
      {
         _fmpz_poly_set_coeff_si(theta_1, start2 + i, array1[4*i+K]);
      }
   }
   
   theta_1->length = MIKE_LIMIT2;
   
   flint_heap_free(array1);

   printf("Computed first theta function.\n");

   //-------------------------------------------------------------
   
   long * array2 = (long *) flint_heap_alloc(BLOCK);
  
   fmpz_poly_init2(theta_C, MIKE_LIMIT2, 1); 
   
   long character2[4] = {1, -1, 1, -1};
	
   for (unsigned long start = 0; start < MIKE_LIMIT; start += BLOCK)
   {   
	  unsigned long start2 = start/4;
	  
	  theta_2d_C(array2, start2, BLOCK/4);
   
	  for (unsigned long i = 0; i < BLOCK/4; i++)
      {
         _fmpz_poly_set_coeff_si(theta_C, start2 + i, array2[i]);
      }
   }
   
   theta_C->length = MIKE_LIMIT2;
   
   flint_heap_free(array2);
   
   printf("Computed second theta function.\n");

   //----------------------------------------------------------------------

   fmpz_poly_init2(theta_prod, MIKE_LIMIT2, 1);
   
   _fmpz_poly_mul_KS_trunc(theta_prod, theta_1, theta_C, MIKE_LIMIT2, -24);
   
   fmpz_poly_clear(theta_1);
   fmpz_poly_clear(theta_C);
  
   printf("Completed multiplication\n");

   unsigned long s = 0;
   long max = 0;
   long coeff;
   for(unsigned long j = 0; j < MIKE_LIMIT2; j++)
   {
       coeff = fmpz_poly_get_coeff_si(theta_prod, j);
	   if (coeff > max) max = coeff;
	   if (-coeff > max) max = -coeff;
	   if (!coeff) 
	   {
		  s++;
	   }
   }
   printf("max = %ld\n", max);
   printf("\n%ld zeroes\n", s);
   
   fmpz_poly_clear(theta_prod);
   
   
   return 0;
}
