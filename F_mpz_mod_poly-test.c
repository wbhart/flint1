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

F_mpz_mod_poly-test.c: Test code for F_mpz_mod_poly.c and F_mpz_mod_poly.h

Copyright (C) 2009, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "flint.h"
#include "mpz_poly.h"
#include "F_mpz_mod_poly.h"
#include "memory-manager.h"
#include "test-support.h"

#define VARY_BITS 1 // random coefficients have random number of bits up to the limit given
#define SIGNS 1 // random coefficients will be randomly signed
#define SPARSE 1 // polynomials are spares (triggers more corner cases)
#define ITER 1 // if you want all tests to run longer, increase this

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // allows easy switching of debugging code on and off when debugging (if inserted)
#define DEBUG2 1 

// generate a random mpz_poly_t with up to the given length and number of bits per coefficient
void mpz_randpoly(mpz_poly_t pol, long length, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_ensure_alloc(pol, length);
	mpz_poly_zero(pol);
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = z_randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
#if SPARSE
          if (z_randint(10) == 1) mpz_rrandomb(temp, randstate, bits);
          else mpz_set_ui(temp, 0);
#else
          mpz_rrandomb(temp, randstate, bits);
#endif
#if SIGNS
          if (z_randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

// generate a random mpz_poly_t with up to the given length and number of bits per coefficient
void mpz_randpoly_unsigned(mpz_poly_t pol, long length, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_ensure_alloc(pol, length);
	mpz_poly_zero(pol);
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = z_randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
#if SPARSE
          if (z_randint(10) == 1) mpz_rrandomb(temp, randstate, bits);
          else mpz_set_ui(temp, 0);
#else
          mpz_rrandomb(temp, randstate, bits);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

// generate a dense random mpz_poly_t with up to the given length and number of bits per coefficient
void mpz_randpoly_dense(mpz_poly_t pol, long length, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_ensure_alloc(pol, length);
	mpz_poly_zero(pol);
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = z_randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, randstate, bits);
#if SIGNS
          if (z_randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

// same as for mpz_randpoly above, except it creates an F_mpz_mod_poly
// WARNING: do not use for testing of conversion between the two formats
void F_mpz_mod_randpoly(F_mpz_mod_poly_t poly, ulong length, ulong bits)
{
	mpz_poly_t m_poly;
	mpz_poly_init(m_poly);
	mpz_randpoly(m_poly, length, bits);
	mpz_poly_to_F_mpz_mod_poly(poly, m_poly);
	mpz_poly_clear(m_poly);
}

void F_mpz_random_modulus(F_mpz_t P, ulong bits)
{
   F_mpz_random(P, bits);
   if (F_mpz_is_zero(P)) F_mpz_add_ui(P, P, 1);
}

void mpz_poly_reduce(mpz_poly_t m_poly, F_mpz_t P)
{
   mpz_t mP;
   mpz_init(mP);
   
   F_mpz_get_mpz(mP, P);
   for (ulong i = 0; i < m_poly->length; i++)
      mpz_mod(m_poly->coeffs[i], m_poly->coeffs[i], mP);
   mpz_poly_normalise(m_poly);

   mpz_clear(mP);
}

int test_F_mpz_mod_poly_convert()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_mod_poly_t F_poly;
   F_mpz_t P;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_init(P);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
      
      F_mpz_random(P, bits);
      F_mpz_add_ui(P, P, 1);
      F_mpz_mod_poly_init(F_poly, P);
      
      mpz_poly_to_F_mpz_mod_poly(F_poly, m_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly);
      mpz_poly_reduce(m_poly1, P);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly);
      F_mpz_clear(P);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_mod_poly_mul()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_mod_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits, length1, length2;
   F_mpz_t P;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   F_mpz_init(P);
     
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_mul(res, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_mul(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
		F_mpz_mod_poly_clear(res);
   }
   
   // alias poly1 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_mul(F_poly1, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly1);

      mpz_poly_mul(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly2 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_mul(F_poly2, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly2);

      mpz_poly_mul(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly1 and poly2
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      
      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      
		F_mpz_mod_poly_mul(res, F_poly1, F_poly1);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_mul(res1, m_poly1, m_poly1);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(res);
   }
   
   F_mpz_clear(P);
   mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_mod_poly_mul_trunc_left()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_mod_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits, length1, length2, trunc;
   F_mpz_t P;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   F_mpz_init(P);
     
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);

      trunc = z_randint(length1 + length2);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_mul_trunc_left(res, F_poly1, F_poly2, trunc);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_mul(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);
      for (ulong i = 0; i < FLINT_MIN(res1->length, trunc); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
      mpz_poly_normalise(res2);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld, trunc = %ld\n", m_poly1->length, bits, m_poly2->length, trunc);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
		F_mpz_mod_poly_clear(res);
   }
   
   // alias poly1 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      trunc = z_randint(length1 + length2);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_mul_trunc_left(F_poly1, F_poly1, F_poly2, trunc);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly1);

      mpz_poly_mul(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);
      for (ulong i = 0; i < FLINT_MIN(res1->length, trunc); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
      mpz_poly_normalise(res2);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld, trunc = %ld\n", m_poly1->length, bits, m_poly2->length, trunc);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly2 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      trunc = z_randint(length1 + length2);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_mul_trunc_left(F_poly2, F_poly1, F_poly2, trunc);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly2);

      mpz_poly_mul(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);
      for (ulong i = 0; i < FLINT_MIN(res1->length, trunc); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
      mpz_poly_normalise(res2);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld, trunc = %ld\n", m_poly1->length, bits, m_poly2->length, trunc);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly1 and poly2
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		
      trunc = z_randint(2*length1);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      
      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      
		F_mpz_mod_poly_mul_trunc_left(res, F_poly1, F_poly1, trunc);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_mul(res1, m_poly1, m_poly1);
		mpz_poly_reduce(res1, P);
      for (ulong i = 0; i < FLINT_MIN(res1->length, trunc); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
      mpz_poly_normalise(res2);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld, trunc = %ld\n", m_poly1->length, bits, m_poly2->length, trunc);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(res);
   }
   
   F_mpz_clear(P);
   mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_mod_poly_add()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_mod_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits, length1, length2;
   F_mpz_t P;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   F_mpz_init(P);
     
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_add(res, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_add(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
		F_mpz_mod_poly_clear(res);
   }
   
   // alias poly1 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_add(F_poly1, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly1);

      mpz_poly_add(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly2 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_add(F_poly2, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly2);

      mpz_poly_add(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly1 and poly2
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      
      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      
		F_mpz_mod_poly_add(res, F_poly1, F_poly1);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_add(res1, m_poly1, m_poly1);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(res);
   }
   
   F_mpz_clear(P);
   mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_mod_poly_sub()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_mod_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits, length1, length2;
   F_mpz_t P;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   F_mpz_init(P);
     
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_sub(res, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_sub(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
		F_mpz_mod_poly_clear(res);
   }
   
   // alias poly1 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_sub(F_poly1, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly1);

      mpz_poly_sub(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly2 and res
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      F_mpz_mod_randpoly(F_poly1, length1, bits);
      F_mpz_mod_randpoly(F_poly2, length2, bits);

      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      
		F_mpz_mod_poly_sub(F_poly2, F_poly1, F_poly2);
		F_mpz_mod_poly_to_mpz_poly(res2, F_poly2);

      mpz_poly_sub(res1, m_poly1, m_poly2);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
   }
   
   // alias poly1 and poly2
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length1 = z_randint(500);
		
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(res, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      
      F_mpz_mod_poly_to_mpz_poly(m_poly1, F_poly1);
      
		F_mpz_mod_poly_sub(res, F_poly1, F_poly1);
		F_mpz_mod_poly_to_mpz_poly(res2, res);

      mpz_poly_sub(res1, m_poly1, m_poly1);
		mpz_poly_reduce(res1, P);

      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits = %ld, m_poly2->length = %ld\n", m_poly1->length, bits, m_poly2->length);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(res);
   }
   
   F_mpz_clear(P);
   mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

void F_mpz_mod_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
	
   RUN_TEST(F_mpz_mod_poly_convert); 
   RUN_TEST(F_mpz_mod_poly_add); 
   RUN_TEST(F_mpz_mod_poly_sub); 
   RUN_TEST(F_mpz_mod_poly_mul); 
   RUN_TEST(F_mpz_mod_poly_mul_trunc_left); 

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_mod_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();
	_F_mpz_cleanup();

   return 0;
}
