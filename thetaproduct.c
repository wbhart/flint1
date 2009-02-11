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

#define LIMIT 416700000L
#define BLOCK    100000
#define COUNT (LIMIT/BLOCK)

int main(void)
{
   fmpz_poly_t theta_1, theta_2, theta_3, theta_prod;
   
   //--------------------------------------------------------------
   
   long * array1 = (long *) flint_heap_alloc(BLOCK);
   
   fmpz_poly_init2(theta_1, LIMIT, 1); 
   
   for (unsigned long start = 0; start < LIMIT*2; start += BLOCK)
   {   
	  theta_1d_0(1, 1, 0, array1, start, BLOCK);
   
     unsigned long start2 = start/2;
	  
	  for (unsigned long i = 0; i < BLOCK/2; i++)
      {
         _fmpz_poly_set_coeff_si(theta_1, start2 + i, array1[2*i]);
      }
   }
   
   theta_1->length = LIMIT;
   
   flint_heap_free(array1);

   printf("Computed first theta function.\n");

	//fmpz_poly_print_pretty(theta_1, "x"); printf("\n");

   //-------------------------------------------------------------
   
   array1 = (long *) flint_heap_alloc(BLOCK);
	
	long character[4] = {1, -1, 1, -1};
	
   fmpz_poly_init2(theta_2, LIMIT, 1); 
   
   for (unsigned long start = 0; start < LIMIT; start += BLOCK)
   {   
	  theta_1d_quadchar(character, 3, 1, 0, array1, start, BLOCK);
   
     unsigned long start2 = start;
	  
	  for (unsigned long i = 0; i < BLOCK/2; i++)
      {
         _fmpz_poly_set_coeff_si(theta_2, start2 + 2*i, array1[2*i]);
      }
   }
   
   theta_2->length = LIMIT;
   
   flint_heap_free(array1);

   printf("Computed second theta function.\n");

	//fmpz_poly_print_pretty(theta_2, "x"); printf("\n");

   //----------------------------------------------------------------------

	array1 = (long *) flint_heap_alloc(BLOCK);
	
	long character2[4] = {1, -1, 1, -1};
	
   fmpz_poly_init2(theta_3, LIMIT+BLOCK, 1); 
   
   for (unsigned long start=0; start < 1+(LIMIT-1)/3; start += BLOCK)
   {   
	  theta_1d_quadchar(character2, 1, 0, 0, array1, start, BLOCK);
   
     unsigned long start2 = start*3;
	  
	  for (unsigned long i = 0; i < BLOCK; i++)
      {
         _fmpz_poly_set_coeff_si(theta_3, start2 + 3*i, array1[i]);
      }
   }
   
   theta_3->length = LIMIT;
   
   flint_heap_free(array1);

   printf("Computed third theta function.\n");

	//fmpz_poly_print_pretty(theta_3, "x"); printf("\n");

	//----------------------------------------------------------------------

   fmpz_poly_init2(theta_prod, LIMIT, 1);
   
   _fmpz_poly_mul_KS_trunc(theta_prod, theta_1, theta_2, LIMIT, -24);
   
	fmpz_poly_clear(theta_1);
   fmpz_poly_clear(theta_2);
   
	_fmpz_poly_mul_KS_trunc(theta_prod, theta_prod, theta_3, LIMIT, -24);
   
   fmpz_poly_clear(theta_3);
  
   printf("Completed multiplication\n");

   #define OFFSET (1L<<17)
   #define LEN (1L<<18)

   // sieve out non-squarefree coefficients
   timeit_start(t0);

   for(long a = 1; a < limit ; a++) {
      long ab = a*2; // a*b
      long ab2 = ab*2; // a*b*b
      while (ab2 < limit) {
          // check the coefficient is in our series
          if ( (ab2-k) % mod == 0 )
                  // flag non-squarefree coefficients as -OFFSET
                  fmpz_poly_set_coeff_si(theta_prod, (ab2-k)/mod, -OFFSET);

          // iterate b++, ab = a*b, ab2 = a*b^2
          ab2 += ab;
          ab += a;
          ab2 += ab;
      }
   }

   timeit_stop(t0);
   fprintf(stderr, "Sieve out non-squarefree coefficients: cpu = %ld ms  wall = %ld ms\n", t0->cpu, t0->wall);

   timeit_start(t0);
   unsigned long arr[LEN];
   for(long i = 0; i < LEN; i++)
      arr[i] = 0;

   unsigned long s = 0;
   long maxneg = 0;
   long maxpos = 0;
   long coeff;
   for(unsigned long j = k; j < limit; j += mod)
   {
      coeff = fmpz_poly_get_coeff_si(theta_prod, j/mod);

      // skip non-squarefree coefficients
      if(coeff == -OFFSET) {
          //printf("arr[%ld] : non-squarefree\n", j);
          continue;
      }

      //printf("arr[%ld] = %ld\n", j, coeff);

      arr[OFFSET+coeff]++;
      if (coeff > maxpos) maxpos = coeff;
      if (coeff < maxneg) maxneg = coeff;
      if (!coeff)
      {
         printf("%ld\n", j);
         s++;
      }
   }
   printf("\n\nmaxneg = %ld, maxpos = %ld\nnumzeros = %ld\n\n", maxneg, maxpos, s);
   for(long i = maxneg; i <= maxpos; i++) {
      printf("VALUE = %ld, count = %ld\n", i, arr[i+OFFSET]);
   }
   timeit_stop(t0);
   fprintf(stderr, "Counting zeroes: cpu = %ld ms  wall = %ld ms\n", t0->cpu, t0->wall);
   fprintf(stderr, "\n%ld zeroes\n", s);
   
	theta_prod->length = 30000;
	//	fmpz_poly_print_pretty(theta_prod, "x");

   fmpz_poly_clear(theta_prod);
   
   
   return 0;
}
