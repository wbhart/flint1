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
#include <omp.h>
#include "flint.h"
#include "F_mpz.h"
#include "F_mpz_poly.h"
#include "theta.h"
#include "profiler.h"

#define LIMIT 3000000000L
#define BLOCK  100000L
#define BUNDLE 200L
#define BYTES 2L

#define MOD 8L
#define K 1L

#define COUNT (LIMIT/BLOCK)

int main(void)
{
   FILE * myfile, * myfile2;
   myfile = fopen("zeros1mod8", "w"); 
   myfile2 = fopen("stats1mod8", "w"); 
   
   F_mpz_poly_t theta_1, theta_2, theta_3, theta_prod, p1, p2, out;
   
   //--------------------------------------------------------------
   
   long * array1 = (long *) flint_heap_alloc(16*BLOCK);
   
   F_mpz_poly_init(theta_1);
   F_mpz_poly_fit_length(theta_1, LIMIT);
   theta_1->length = LIMIT;
   
#pragma omp parallel
	{
#pragma omp for
	for (long start = 0; start < LIMIT*MOD; start += BLOCK)
   {   
	  ulong j = omp_get_thread_num();
	  theta_2d_A1(array1 + j*BLOCK, start, BLOCK);
   
      ulong start2 = start/MOD;
	  
	  for (ulong i = 0; i < BLOCK/MOD; i++)
      {
         F_mpz_set_si(theta_1->coeffs + start2 + i, array1[j*BLOCK + MOD*i+K]);
      }
   }
	}
   
   _F_mpz_poly_normalise(theta_1);

	flint_heap_free(array1);

   F_mpz_poly_init(p1);

   F_mpz_poly_pack_bytes(p1, theta_1, BUNDLE, BYTES);

   F_mpz_poly_clear(theta_1);
    
   printf("Computed and packed first theta function.\n");
  
   //-------------------------------------------------------------
   
   array1 = (long *) flint_heap_alloc(16*BLOCK);
	
   F_mpz_poly_init(theta_2);
   F_mpz_poly_fit_length(theta_2, LIMIT); 
   
   theta_2->length = LIMIT;
   
#pragma omp parallel
	{
#pragma omp for
	for (long start = 0; start < LIMIT*MOD; start += BLOCK)
   {   
	  ulong j = omp_get_thread_num();
	  theta_2d_B(array1 + j*BLOCK, start, BLOCK);
   
      ulong start2 = start/MOD;
	  
	  for (ulong i = 0; i < BLOCK/MOD; i++)
      {
         F_mpz_set_si(theta_2->coeffs + start2 + i, array1[j*BLOCK + MOD*i]);
      }
   }
	}
   
   _F_mpz_poly_normalise(theta_2);

   flint_heap_free(array1);

   F_mpz_poly_init(p2);

   F_mpz_poly_pack_bytes(p2, theta_2, BUNDLE, BYTES);

   F_mpz_poly_clear(theta_2);

   printf("Computed and packed second theta function.\n");

   //----------------------------------------------------------------------

   F_mpz_poly_init(out);
   
   F_mpz_poly_mul_modular_trunc(out, p1, p2, 0, (LIMIT - 1)/BUNDLE + 1);
	
   printf("First product computed\n");

   F_mpz_poly_init(theta_prod);
   
   //F_mpz_poly_unpack_bytes(theta_prod, out, BUNDLE, 2);

	ulong limbs = ((2*BUNDLE)*BYTES*8 - 1)/FLINT_BITS + 2; // max number of limbs of each large coeff
	ulong extras = (FLINT_BYTES_PER_LIMB - 1)/BYTES + 1;

	ulong length_max = BUNDLE*out->length + BUNDLE + extras;

	F_mpz_poly_fit_length(theta_prod, length_max);
	
	// zero coeffs, as we will be adding to them
	theta_prod->length = length_max;
   for (ulong i = 0; i < length_max; i++)
		F_mpz_zero(theta_prod->coeffs + i);

	mp_limb_t * arr1 = flint_heap_alloc(limbs);
	mp_limb_t * arr2 = flint_heap_alloc(limbs);
	mp_limb_t * arr3 = flint_heap_alloc(limbs);

	int neg1 = 0;
	int neg2 = 0;
	int neg3 = 0;

   mp_limb_t * temp;
	int temp_n;
	long limb;
	long carry = 0L;
	int borrow = 0;
   int borrow2 = 0;

	F_mpn_clear(arr1, limbs);
	F_mpn_clear(arr2, limbs);
	F_mpn_clear(arr3, limbs);

	short int * s1 = (short int *) arr1;
	short int * s2 = (short int *) arr2;
	short int * s3 = (short int *) arr3;
		
   F_mpz_get_limbs(arr1, out->coeffs); // initialise with first big coefficient
	if (F_mpz_sgn(out->coeffs) < 0) neg1 = 1;

	for (ulong i = 0; i < length_max; i+=BUNDLE)
	{
		ulong j;

	   for (j = 0; (j < extras) && (i + j < length_max); j++)
	   {
         if (neg1) 
			{
				borrow2 = -(s1[j] < 0);
			   carry -= (long) s1[j];
			} else 
			{
				borrow2 = (s1[j] < 0);
			   carry += (long) s1[j];
			}

			if (neg2) 
			{
				borrow2 -= (s2[BUNDLE + j] < 0);
			   carry -= (long) s2[BUNDLE + j];
			} else 
			{
				borrow2 += (s2[BUNDLE + j] < 0);
			   carry += (long) s2[BUNDLE + j];
			}

			if (neg3) 
			{
				borrow2 -= (s3[2*BUNDLE + j] < 0);
			   carry -= (long) s3[2*BUNDLE + j];
			} else 
			{
				borrow2 += (s3[2*BUNDLE + j] < 0);
			   carry += (long) s3[2*BUNDLE + j];
			}

			carry += borrow;

			theta_prod->coeffs[i + j] = (long) ((short int) carry);
			carry -= theta_prod->coeffs[i + j];
			carry >>= (BYTES*8);

			borrow = borrow2;
		}

		for ( ; (j < BUNDLE) && (i + j < length_max); j++)
		{
         if (neg1) 
			{
				borrow2 = -(s1[j] < 0);
			   carry -= (long) s1[j];
			} else 
			{
				borrow2 = (s1[j] < 0);
			   carry += (long) s1[j];
			}

			if (neg2) 
			{
				borrow2 -= (s2[BUNDLE + j] < 0);
			   carry -= (long) s2[BUNDLE + j];
			} else 
			{
				borrow2 += (s2[BUNDLE + j] < 0);
			   carry += (long) s2[BUNDLE + j];
			}

			carry += borrow;

			theta_prod->coeffs[i + j] = (long) ((short int) carry);
			carry -= theta_prod->coeffs[i + j];
			carry >>= (BYTES*8);

			borrow = borrow2;
		}

		temp = arr3;
		arr3 = arr2;
		arr2 = arr1;
		arr1 = temp;

		s1 = (short int *) arr1;
	   s2 = (short int *) arr2;
	   s3 = (short int *) arr3;
		
		neg3 = neg2;
		neg2 = neg1;

		F_mpn_clear(arr1, limbs);
	   ulong k = i/BUNDLE + 1;
		if (k < out->length)
		{
			F_mpz_get_limbs(arr1, out->coeffs + k); // initialise with first big coefficient
		   
			if (F_mpz_sgn(out->coeffs + k) < 0) neg1 = 1;
		   else neg1 = 0;
	   }
	}

   flint_heap_free(arr1);
	flint_heap_free(arr2);
	flint_heap_free(arr3);

   F_mpz_poly_clear(out);
   _F_mpz_cleanup2();

   F_mpz_poly_truncate(theta_prod, LIMIT);

	printf("First unpacking computed, theta_prod has length %ld\n", theta_prod->length);

   #define OFFSET (1L<<17)
   #define LEN (1L<<18)

   // sieve out non-squarefree coefficients

   for(long a = 1; a < LIMIT*MOD + K ; a++) {
      long ab = a*2L; // a*b
      long ab2 = ab*2L; // a*b*b
      while (ab2 < LIMIT*MOD + K) {
          // check the coefficient is in our series
          if ( (ab2-K) % MOD == 0L )
                  // flag non-squarefree coefficients as -OFFSET
                  F_mpz_set_si(theta_prod->coeffs + (ab2-K)/MOD, -OFFSET);

          // iterate b++, ab = a*b, ab2 = a*b^2
          ab2 += ab;
          ab += a;
          ab2 += ab;
      }
   }

   fprintf(stderr, "Sieve out non-squarefree coefficients done\n");

   unsigned long arr[LEN];
   for(long i = 0; i < LEN; i++)
      arr[i] = 0L;

   unsigned long s = 0L;
   long maxneg = 0L;
   long maxpos = 0L;
   long coeff;
   for(unsigned long j = 0; j < LIMIT; j ++)
   {
      coeff = F_mpz_get_si(theta_prod->coeffs + j);

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
         fprintf(myfile, "%ld ", MOD*j+K);
         s++;
      }
   }

   fprintf(myfile, "\n");
   fclose(myfile);

   for(long i = maxneg; i <= maxpos; i++) {
      fprintf(myfile2, "VALUE = %ld, count = %ld\n", i, arr[i+OFFSET]);
   }
   
   fclose(myfile2);

   printf("\n\nmaxneg = %ld, maxpos = %ld\nnumzeros = %ld\n\n", maxneg, maxpos, s);

   F_mpz_poly_clear(theta_prod);
   
   return 0;
}
