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
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

// same as for mpz_randpoly above, except it creates an F_mpz_mod_poly
// WARNING: do not use for testing of conversion between the two formats
void F_mpz_mod_randpoly(F_mpz_mod_poly_t poly, ulong length, ulong bits)
{
	if (length == 0) 
   {
      F_mpz_mod_poly_zero(poly);
      return;
   }
   
   mpz_poly_t m_poly;
	mpz_poly_init(m_poly);
	mpz_randpoly(m_poly, length, bits);
	mpz_poly_to_F_mpz_mod_poly(poly, m_poly);
	mpz_poly_clear(m_poly);
}

void F_mpz_randpoly(F_mpz_poly_t poly, ulong length, ulong bits)
{
	mpz_poly_t m_poly;
	mpz_poly_init(m_poly);
	mpz_randpoly(m_poly, length, bits);
	mpz_poly_to_F_mpz_poly(poly, m_poly);
	mpz_poly_clear(m_poly);
}

void F_mpz_random_modulus(F_mpz_t P, ulong bits)
{
   F_mpz_random(P, bits);

   if (F_mpz_is_zero(P)) F_mpz_add_ui(P, P, 1);
}

void F_mpz_random_prime_modulus(F_mpz_t P, ulong bits)
{
   do
   {
      F_mpz_random(P, bits);
      if (COEFF_IS_MPZ(*P))
      {
         __mpz_struct * mpz_ptr = F_mpz_ptr_mpz(*P);
         mpz_nextprime(mpz_ptr, mpz_ptr);
         if (mpz_sizeinbase(mpz_ptr, 2) > bits) *P = 0L;
      } else
      {
         *P = z_nextprime(*P, 0);
         if ((FLINT_BIT_COUNT(*P) > bits) || (FLINT_BIT_COUNT(*P) > FLINT_BITS - 2)) 
            *P = 0L;
      }
   } while (*P == 0L);
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

void F_mpz_poly_reduce(F_mpz_poly_t m_poly, F_mpz_t P)
{
   for (ulong i = 0; i < m_poly->length; i++)
      F_mpz_mod(m_poly->coeffs + i, m_poly->coeffs + i, P);
   _F_mpz_poly_normalise(m_poly);
}

int test_F_mpz_mod_poly_init_realloc_clear()
{
   F_mpz_mod_poly_t F_poly;
   F_mpz_t P;
   int result = 1;
   ulong bits, length;
      
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_init(P);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      
      F_mpz_random(P, bits);
      F_mpz_add_ui(P, P, 1);
      F_mpz_mod_poly_init(F_poly, P);
      
      F_mpz_mod_randpoly(F_poly, length, bits);
      length = z_randint(100);
      F_mpz_mod_poly_realloc(F_poly, length);
      bits = z_randint(200) + 1;
      F_mpz_mod_randpoly(F_poly, length, bits);
           
      F_mpz_mod_poly_clear(F_poly);
      F_mpz_clear(P);
   }
   
   return result;
}

int test_F_mpz_mod_poly_to_mpz_poly()
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

int test_F_mpz_mod_poly_to_F_mpz_poly()
{
   F_mpz_poly_t m_poly1, m_poly2;
   F_mpz_mod_poly_t F_poly;
   F_mpz_t P;
   int result = 1;
   ulong bits, length;
   
   F_mpz_poly_init(m_poly1); 
   F_mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_init(P);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      F_mpz_randpoly(m_poly1, length, bits);
      
      F_mpz_random(P, bits);
      F_mpz_add_ui(P, P, 1);
      F_mpz_mod_poly_init(F_poly, P);
      
      F_mpz_poly_to_F_mpz_mod_poly(F_poly, m_poly1);
      F_mpz_mod_poly_to_F_mpz_poly(m_poly2, F_poly);
      F_mpz_poly_reduce(m_poly1, P);
          
      result = F_mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         F_mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         F_mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly);
      F_mpz_clear(P);
   }
   
   F_mpz_poly_clear(m_poly1);
   F_mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_mod_poly_to_zmod_poly()
{
   zmod_poly_t m_poly;
   F_mpz_mod_poly_t F_poly1, F_poly2;
   F_mpz_t P;
   int result = 1;
   ulong bits, length;
     
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_init(P);

      bits = z_randint(FLINT_BITS) + 1;
      length = z_randint(100);
      
      F_mpz_random(P, bits);
      F_mpz_add_ui(P, P, 1);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      zmod_poly_init(m_poly, F_mpz_get_ui(P)); 

      F_mpz_mod_randpoly(F_poly1, length, bits);
      
      F_mpz_mod_poly_to_zmod_poly(m_poly, F_poly1);
      zmod_poly_to_F_mpz_mod_poly(F_poly2, m_poly);
          
      result = F_mpz_mod_poly_equal(F_poly1, F_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, F_poly1->length, F_poly2->length);
         zmod_poly_print(m_poly); printf("\n");
 		}
          
      zmod_poly_clear(m_poly);
      F_mpz_mod_poly_clear(F_poly1);
      F_mpz_mod_poly_clear(F_poly2);
      F_mpz_clear(P);
   }
    
   return result;
}

int test_F_mpz_mod_poly_scalar_mul()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_mod_poly_t F_poly, F_poly2;
   F_mpz_t P, x;
   mpz_t mx;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_init(mx);

   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_init(P);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
      
      F_mpz_random(P, bits);
      F_mpz_add_ui(P, P, 1);
      F_mpz_random(x, bits);
      F_mpz_mod(x, x, P);
      F_mpz_get_mpz(mx, x);
      F_mpz_mod_poly_init(F_poly, P);
      F_mpz_mod_poly_init(F_poly2, P);
      
      mpz_poly_to_F_mpz_mod_poly(F_poly, m_poly1);
      F_mpz_mod_poly_scalar_mul(F_poly2, F_poly, x);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly2);
      for (ulong i = 0; i < m_poly1->length; i++)
         mpz_mul(m_poly1->coeffs[i], m_poly1->coeffs[i], mx);
      mpz_poly_reduce(m_poly1, P);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
         F_mpz_print(x); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly);
      F_mpz_mod_poly_clear(F_poly2);
      F_mpz_clear(P);
   }

   // alias operands
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_init(P);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
      
      F_mpz_random(P, bits);
      F_mpz_add_ui(P, P, 1);
      F_mpz_random(x, bits);
      F_mpz_mod(x, x, P);
      F_mpz_get_mpz(mx, x);
      F_mpz_mod_poly_init(F_poly, P);
      
      mpz_poly_to_F_mpz_mod_poly(F_poly, m_poly1);
      F_mpz_mod_poly_scalar_mul(F_poly, F_poly, x);
      F_mpz_mod_poly_to_mpz_poly(m_poly2, F_poly);
      for (ulong i = 0; i < m_poly1->length; i++)
         mpz_mul(m_poly1->coeffs[i], m_poly1->coeffs[i], mx);
      mpz_poly_reduce(m_poly1, P);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
         F_mpz_print(x); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly);
      F_mpz_clear(P);
   }

   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   mpz_clear(mx);

   return result;
}

int test_F_mpz_mod_poly_setequal()
{
   F_mpz_mod_poly_t F_poly1, F_poly2;
   F_mpz_t P;
   int result = 1;
   ulong bits, length;

   // check equal polys
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length = z_randint(100);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
 
      F_mpz_mod_randpoly(F_poly1, length, bits);
      F_mpz_mod_poly_set(F_poly2, F_poly1);
          
      result = F_mpz_mod_poly_equal(F_poly1, F_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, len1 = %ld, len2 = %ld\n", length, bits, F_poly1->length, F_poly2->length);
		}
          
      F_mpz_mod_poly_clear(F_poly1);
      F_mpz_mod_poly_clear(F_poly2);
      F_mpz_clear(P);
   }

   // check unequal polys
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length = z_randint(100);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
 
      F_mpz_mod_randpoly(F_poly1, length, bits);
      F_mpz_mod_poly_set(F_poly2, F_poly1);
      if (F_poly2->length == 0)
      {
         F_mpz_mod_poly_fit_length(F_poly2, 1);
         F_mpz_set_ui(F_poly2->coeffs, 1);
         F_poly2->length = 1;
      } else
         F_mpz_add_ui(F_poly2->coeffs, F_poly2->coeffs, 1);
          
      result = !F_mpz_mod_poly_equal(F_poly1, F_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, len1 = %ld, len2 = %ld\n", length, bits, F_poly1->length, F_poly2->length);
		}
          
      F_mpz_mod_poly_clear(F_poly1);
      F_mpz_mod_poly_clear(F_poly2);
      F_mpz_clear(P);
   }

   return result;
}

int test_F_mpz_mod_poly_swap()
{
   F_mpz_mod_poly_t F_poly1, F_poly2, F_poly3;
   F_mpz_t P;
   int result = 1;
   ulong bits, length;

   // check equal polys
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 1;
      length = z_randint(100);
     
      F_mpz_random_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(F_poly3, P);
 
      F_mpz_mod_randpoly(F_poly1, length, bits);
      F_mpz_mod_randpoly(F_poly2, length, bits);
      F_mpz_mod_poly_set(F_poly3, F_poly1);
      F_mpz_mod_poly_swap(F_poly2, F_poly1);
          
      result = F_mpz_mod_poly_equal(F_poly3, F_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, len3 = %ld, len2 = %ld\n", length, bits, F_poly3->length, F_poly2->length);
		}
          
      F_mpz_mod_poly_clear(F_poly1);
      F_mpz_mod_poly_clear(F_poly2);
      F_mpz_mod_poly_clear(F_poly3);
      F_mpz_clear(P);
   }

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

int test_F_mpz_mod_poly_divrem_basecase()
{
   F_mpz_mod_poly_t F_poly1, F_poly2, F_poly3, Q, R;
   int result = 1;
   ulong bits, length1, length2;
   F_mpz_t P;
   
   F_mpz_init(P);
     
   // test exact division
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 2;
      length1 = z_randint(500);
		length2 = z_randint(500) + 1;
     
      F_mpz_random_prime_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(F_poly3, P);
      F_mpz_mod_poly_init(Q, P);
      F_mpz_mod_poly_init(R, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      do {F_mpz_mod_randpoly(F_poly2, length2, bits);} while (F_poly2->length == 0);
      
		F_mpz_mod_poly_mul(F_poly3, F_poly1, F_poly2);
		F_mpz_mod_poly_divrem_basecase(Q, R, F_poly3, F_poly2);

      result = F_mpz_mod_poly_equal(Q, F_poly1); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits = %ld, length2 = %ld, bits\n", F_poly1->length, bits, F_poly2->length, bits);
         F_mpz_mod_poly_print(Q); printf("\n");
         F_mpz_mod_poly_print(F_poly1); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
		F_mpz_mod_poly_clear(F_poly3);
		F_mpz_mod_poly_clear(Q);
		F_mpz_mod_poly_clear(R);
   }
      
   F_mpz_clear(P);
   
   return result;
}

int test_F_mpz_mod_poly_divrem_divconquer()
{
   F_mpz_mod_poly_t F_poly1, F_poly2, F_poly3, Q, R;
   int result = 1;
   ulong bits, length1, length2;
   F_mpz_t P;
   
   F_mpz_init(P);
     
   // test exact division
   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
		bits = z_randint(200) + 2;
      length1 = z_randint(500);
		length2 = z_randint(500) + 1;
     
      F_mpz_random_prime_modulus(P, bits);
      F_mpz_mod_poly_init(F_poly1, P);
      F_mpz_mod_poly_init(F_poly2, P);
      F_mpz_mod_poly_init(F_poly3, P);
      F_mpz_mod_poly_init(Q, P);
      F_mpz_mod_poly_init(R, P);

      F_mpz_mod_randpoly(F_poly1, length1, bits);
      do {F_mpz_mod_randpoly(F_poly2, length2, bits);} while (F_poly2->length == 0);
      
		F_mpz_mod_poly_mul(F_poly3, F_poly1, F_poly2);
		F_mpz_mod_poly_divrem_divconquer(Q, R, F_poly3, F_poly2);
      
      result = F_mpz_mod_poly_equal(Q, F_poly1); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits = %ld, length2 = %ld, bits\n", F_poly1->length, bits, F_poly2->length, bits);
         F_mpz_mod_poly_print(Q); printf("\n");
         F_mpz_mod_poly_print(F_poly1); printf("\n");
		}
          
      F_mpz_mod_poly_clear(F_poly1);
		F_mpz_mod_poly_clear(F_poly2);
		F_mpz_mod_poly_clear(F_poly3);
		F_mpz_mod_poly_clear(Q);
		F_mpz_mod_poly_clear(R);
   }
      
   F_mpz_clear(P);
   
   return result;
}

void F_mpz_mod_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
	
   RUN_TEST(F_mpz_mod_poly_init_realloc_clear); 
   RUN_TEST(F_mpz_mod_poly_to_mpz_poly); 
   RUN_TEST(F_mpz_mod_poly_to_F_mpz_poly); 
   RUN_TEST(F_mpz_mod_poly_to_zmod_poly); 
   RUN_TEST(F_mpz_mod_poly_setequal); 
   RUN_TEST(F_mpz_mod_poly_swap); 
   RUN_TEST(F_mpz_mod_poly_add); 
   RUN_TEST(F_mpz_mod_poly_sub); 
   RUN_TEST(F_mpz_mod_poly_scalar_mul); 
   RUN_TEST(F_mpz_mod_poly_mul); 
   RUN_TEST(F_mpz_mod_poly_mul_trunc_left); 
   RUN_TEST(F_mpz_mod_poly_divrem_basecase); 
   RUN_TEST(F_mpz_mod_poly_divrem_divconquer); 
   
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
