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
#include "F_mpz.h"
#include "F_mpz_poly.h"
#include "theta.h"
#include "profiler.h"

#define LIMIT 2000000000L
#define BLOCK  100000
#define BUNDLE 100

#define COUNT (LIMIT/BLOCK)

int main(void)
{
   F_mpz_poly_t theta_1, theta_2, theta_3, theta_prod, p1, p2, out;
   
   //--------------------------------------------------------------
   
   long * array1 = (long *) flint_heap_alloc(BLOCK);
   
   F_mpz_poly_init(theta_1);
   F_mpz_poly_fit_length(theta_1, LIMIT);
   theta_1->length = LIMIT;
   
   for (ulong start = 0; start < LIMIT*2; start += BLOCK)
   {   
	  theta_1d_0(1, 1, 0, array1, start, BLOCK);
   
      ulong start2 = start/2;
	  
	  for (ulong i = 0; i < BLOCK/2; i++)
      {
         F_mpz_set_si(theta_1->coeffs + start2 + i, array1[2*i]);
      }
   }
   
   _F_mpz_poly_normalise(theta_1);

   flint_heap_free(array1);

   F_mpz_poly_init(p1);

   F_mpz_poly_pack_bytes(p1, theta_1, BUNDLE, 3);

   F_mpz_poly_clear(theta_1);
    

   printf("Computed and packed first theta function.\n");

   //-------------------------------------------------------------
   
   array1 = (long *) flint_heap_alloc(BLOCK);
	
   long character[4] = {1, -1, 1, -1};
	
   F_mpz_poly_init(theta_2);
   F_mpz_poly_fit_length(theta_2, LIMIT); 
   
   theta_2->length = LIMIT;
   
   for (ulong start = 0; start < LIMIT; start += BLOCK)
   {   
	  theta_1d_quadchar(character, 3, 1, 0, array1, start, BLOCK);
   
      ulong start2 = start;
	  
	  for (ulong i = 0; i < BLOCK/2; i++)
      {
         F_mpz_set_si(theta_2->coeffs + start2 + 2*i, array1[2*i]);
      }
   }
   
   _F_mpz_poly_normalise(theta_2);

   flint_heap_free(array1);

   F_mpz_poly_init(p2);

   F_mpz_poly_pack_bytes(p2, theta_2, BUNDLE, 3);

   F_mpz_poly_clear(theta_2);

   printf("Computed and packed second theta function.\n");

   //----------------------------------------------------------------------

   F_mpz_poly_init(out);
   
   F_mpz_poly_mul_modular_trunc(out, p1, p2, 0, (LIMIT - 1)/BUNDLE + 1);
	
   printf("First product computed\n");

   F_mpz_poly_init(theta_prod);
   F_mpz_poly_fit_length(theta_prod, LIMIT + BUNDLE);
   
   F_mpz_poly_unpack_bytes(theta_prod, out, BUNDLE, 3);

   F_mpz_poly_clear(out);
   _F_mpz_cleanup2();
    
   F_mpz_poly_truncate(theta_prod, LIMIT);

   printf("First unpacking computed\n");

   F_mpz_poly_init(p1);

   F_mpz_poly_pack_bytes(p1, theta_prod, BUNDLE, 3);

   F_mpz_poly_clear(theta_prod);

   printf("Repacking completed\n");

   //----------------------------------------------------------------------

   array1 = (long *) flint_heap_alloc(BLOCK);
	
   long character2[4] = {1, -1, 1, -1};
	
   F_mpz_poly_init(theta_3);
   F_mpz_poly_fit_length(theta_3, LIMIT + 3*BLOCK);
   
   theta_3->length = LIMIT;
   
   for (ulong start = 0; start < LIMIT/3; start += BLOCK)
   {   
	  theta_1d_quadchar(character2, 1, 0, 0, array1, start, BLOCK);
    
      ulong start2 = start*3;
	  
	  for (ulong i = 0; i < BLOCK; i++)
      {
         F_mpz_set_si(theta_3->coeffs + start2 + 3*i, array1[i]);
      }
   }
   
   _F_mpz_poly_normalise(theta_3);

   flint_heap_free(array1);

   F_mpz_poly_init(p2);

   F_mpz_poly_pack_bytes(p2, theta_3, BUNDLE, 3);

   F_mpz_poly_clear(theta_3);

   printf("Computed and packed third theta function.\n");

   //----------------------------------------------------------------------

   
   F_mpz_poly_init(out);
   
   F_mpz_poly_mul_modular_trunc(out, p1, p2, 0, (LIMIT - 1)/BUNDLE + 1);
	
   printf("Completed multiplication\n");

   F_mpz_poly_init(theta_prod);
   F_mpz_poly_fit_length(theta_prod, LIMIT + BUNDLE);
   
   F_mpz_poly_unpack_bytes(theta_prod, out, BUNDLE, 3);

   F_mpz_poly_clear(out);
   _F_mpz_cleanup2();
   
   F_mpz_poly_truncate(theta_prod, LIMIT);

   printf("Unpacking complete\n");

   /*timeit_start(t0);
   #define OFFSET (1L<<17)
   #define LEN (1L<<18)

   // sieve out non-squarefree coefficients
   
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

   timeit_t t0;
   
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
   //	fmpz_poly_print_pretty(theta_prod, "x"); */

   F_mpz_poly_clear(theta_prod);
   
   
   return 0;
}
