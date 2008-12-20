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

F_mpz_poly-test.c: Test code for F_mpz_poly.c and F_mpz_poly.h

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "flint.h"
#include "mpz_poly.h"
#include "F_mpz_poly.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "test-support.h"
#include "zmod_poly.h"

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

// same as for mpz_randpoly above, except it creates an F_mpz_poly
// WARNING: do not use for testing of conversion between the two formats
void F_mpz_randpoly(F_mpz_poly_t poly, ulong length, ulong bits)
{
	mpz_poly_t m_poly;
	mpz_poly_init(m_poly);
	mpz_randpoly(m_poly, length, bits);
	mpz_poly_to_F_mpz_poly(poly, m_poly);
	mpz_poly_clear(m_poly);
}

int test_F_mpz_poly_convert()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_add()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_add(res, F_poly1, F_poly2);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_add(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_add(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_add(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_add(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_add(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_sub()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_sub(res, F_poly1, F_poly2);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_sub(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_sub(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_sub(res1, m_poly2, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_sub(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_sub(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_getset_coeff_si()
{
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   long coeff, coeff2;
   ulong coeff_bits, coeff_num;
   
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      // set random coeffs in the poly
		for (ulong count2 = 0; (count2 < 500) && result == 1; count2++)
      {
          coeff_bits = z_randint(FLINT_BITS);
          coeff = z_randbits(coeff_bits);
          coeff_num = z_randint(F_poly->length);
              
			 if (z_randint(2)) coeff = -coeff;
              
			 if (F_poly->length)
          {
              F_mpz_poly_set_coeff_si(F_poly, coeff_num, coeff);
              coeff2 = F_mpz_poly_get_coeff_si(F_poly, coeff_num);
				  
				  result = (coeff2 == coeff);
				  if (!result)
				  {
					  printf("Error: length = %ld, coeff_num = %ld, coeff = %ld, coeff2 = %ld\n", length, coeff_num, coeff, coeff2);
				  }
          }
      }

      F_mpz_poly_clear(F_poly);
   }
   
   return result; 
}

int test_F_mpz_poly_getset_coeff_ui()
{
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   ulong coeff, coeff2;
   ulong coeff_bits, coeff_num;
   
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      // set random coeffs in the poly
		for (ulong count2 = 0; (count2 < 500) && result == 1; count2++)
      {
          coeff_bits = z_randint(FLINT_BITS+1);
          coeff = z_randbits(coeff_bits);
          coeff_num = z_randint(F_poly->length);
                  
			 if (F_poly->length)
          {
              F_mpz_poly_set_coeff_ui(F_poly, coeff_num, coeff);
              coeff2 = F_mpz_poly_get_coeff_ui(F_poly, coeff_num);
				  
				  result = (coeff2 == coeff);
				  if (!result)
				  {
					  printf("Error: length = %ld, coeff_num = %ld, coeff = %u, coeff2 = %u\n", length, coeff_num, coeff, coeff2);
				  }
          }
      }

      F_mpz_poly_clear(F_poly);
   }
   
   return result; 
}

int test_F_mpz_poly_getset_coeff_mpz()
{
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   mpz_t coeff, coeff2;
   ulong coeff_bits, coeff_num;
	mpz_init(coeff);
	mpz_init(coeff2);
   
   for (ulong count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      // set random coeffs in the poly
		for (ulong count2 = 0; (count2 < 300) && result == 1; count2++)
      {
          coeff_bits = z_randint(200);
          mpz_rrandomb(coeff, randstate, coeff_bits);
          coeff_num = z_randint(F_poly->length);
              
			 if (z_randint(2)) mpz_neg(coeff, coeff);
              
			 if (F_poly->length)
          {
              F_mpz_poly_set_coeff_mpz(F_poly, coeff_num, coeff);
              F_mpz_poly_get_coeff_mpz(coeff2, F_poly, coeff_num);
				  
				  result = (mpz_cmp(coeff2, coeff) == 0);
				  if (!result)
				  {
					  gmp_printf("Error: length = %ld, coeff_num = %ld, coeff = %Zd, coeff2 = %Zd\n", length, coeff_num, coeff, coeff2);
				  }
          }
      }

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_clear(coeff);
	mpz_clear(coeff2);
	return result; 
}

int test_F_mpz_poly_set()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      F_mpz_poly_set(F_poly2, F_poly1);
		F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
   }
   
	// aliasing is trivial for set and doesn't need testing

   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_equal()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   
   // random polys unlikely to be equal, test against mpz_poly_equal
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
          
      result = (mpz_poly_equal(m_poly1, m_poly2) == F_mpz_poly_equal(F_poly1, F_poly2)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

	// polys are equal
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      F_mpz_poly_set(F_poly2, F_poly1);
          
      result = (F_mpz_poly_equal(F_poly1, F_poly2)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

	// polys are aliased
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
          
      result = (F_mpz_poly_equal(F_poly1, F_poly1)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
   }

	mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_swap()
{
   F_mpz_poly_t F_poly, F_poly2, F_poly3, F_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   // reverse and reverse back, second time aliased
	for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      bits2 = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);   
      F_mpz_poly_init(F_poly3);   
      F_mpz_poly_init(F_poly4);   
      
      length = z_randint(100);
      length2 = z_randint(100);
      
      F_mpz_randpoly(F_poly, length, bits); 
      F_mpz_randpoly(F_poly2, length2, bits2); 
                
      F_mpz_poly_set(F_poly3, F_poly);
      F_mpz_poly_set(F_poly4, F_poly2);

      F_mpz_poly_swap(F_poly2, F_poly);
		
		result = (F_mpz_poly_equal(F_poly4, F_poly) && F_mpz_poly_equal(F_poly3, F_poly2));
      if (!result)
		{
			printf("Error: length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
			printf("F_poly->length = %ld, F_poly2->length = %ld, F_poly3->length = %ld, F_poly4->length = %ld\n", F_poly->length, F_poly2->length, F_poly3->length, F_poly4->length);
		}

      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly3);
      F_mpz_poly_clear(F_poly4);
   }

	// Aliasing is trivial for swap so doesn't need testing

   return result; 
}

int test_F_mpz_poly_max_bits()
{
   mpz_poly_t m_poly1;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length, next_bits;
	long sign, mpz_bits, test_bits;
   
   mpz_poly_init(m_poly1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_bits = 0;
      sign = 1L;
      for (ulong i = 0; i < m_poly1->length; i++)
      {
         next_bits = mpz_sizeinbase(m_poly1->coeffs[i], 2);
         if (next_bits > mpz_bits) mpz_bits = next_bits;
         if (mpz_sgn(m_poly1->coeffs[i]) < 0L) sign = -1L;
      }
      mpz_bits = sign*mpz_bits;

		test_bits = F_mpz_poly_max_bits(F_poly);
          
      result = (mpz_bits == test_bits); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, m_poly1->length = %ld\n", length, bits, m_poly1->length);
         printf("mpz_bits = %ld, test_bits = %ld\n", mpz_bits, test_bits);
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   
   return result;
}

int test_F_mpz_poly_max_bits1()
{
   mpz_poly_t m_poly1;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length, next_bits;
	long sign, mpz_bits, test_bits;
   
   mpz_poly_init(m_poly1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(FLINT_BITS-2) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_bits = 0;
      sign = 1L;
      for (ulong i = 0; i < m_poly1->length; i++)
      {
         next_bits = mpz_sizeinbase(m_poly1->coeffs[i], 2);
         if (next_bits > mpz_bits) mpz_bits = next_bits;
         if (mpz_sgn(m_poly1->coeffs[i]) < 0L) sign = -1L;
      }
      mpz_bits = sign*mpz_bits;

		test_bits = F_mpz_poly_max_bits1(F_poly);
          
      result = (mpz_bits == test_bits); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, m_poly1->length = %ld\n", length, bits, m_poly1->length);
         printf("mpz_bits = %ld, test_bits = %ld\n", mpz_bits, test_bits);
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   
   return result;
}

int test_F_mpz_poly_max_limbs()
{
   mpz_poly_t m_poly1;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length, next_limbs;
	long sign, mpz_limbs, test_limbs;
   
   mpz_poly_init(m_poly1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(500) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_limbs = 0;
      for (ulong i = 0; i < m_poly1->length; i++)
      {
         next_limbs = mpz_size(m_poly1->coeffs[i]);
         if (next_limbs > mpz_limbs) mpz_limbs = next_limbs;
      }
      
		test_limbs = F_mpz_poly_max_limbs(F_poly);
          
      result = (mpz_limbs == test_limbs); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, m_poly1->length = %ld\n", length, bits, m_poly1->length);
         printf("mpz_limbs = %ld, test_limbs = %ld\n", mpz_limbs, test_limbs);
		}
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly1);
   
   return result;
}

int test_F_mpz_poly_neg()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, F_poly4;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // negate and negate back
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      F_mpz_randpoly(F_poly1, length1, bits1);
      
		F_mpz_poly_neg(F_poly2, F_poly1);
      F_mpz_poly_neg(F_poly2, F_poly2);
          
      result = (F_mpz_poly_equal(F_poly1, F_poly2)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

	// sub equals negate and add, included aliased negation
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(F_poly4);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      bits2 = z_randint(200) + 1;
      length2 = z_randint(100);
      F_mpz_randpoly(F_poly1, length1, bits1);
      F_mpz_randpoly(F_poly2, length2, bits2);
      
		F_mpz_poly_sub(F_poly3, F_poly1, F_poly2);
      F_mpz_poly_neg(F_poly2, F_poly2);
      F_mpz_poly_add(F_poly4, F_poly1, F_poly2);
          
      result = (F_mpz_poly_equal(F_poly3, F_poly4)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(F_poly4);
   }

   return result;
}

int test_F_mpz_poly_reverse()
{
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   unsigned long bits, length, length2;
   
   // reverse and reverse back, second time aliased
	for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);   
      
      length = z_randint(100);
      length2 = length + z_randint(200);
      
      F_mpz_randpoly(F_poly, length, bits); 
                
      F_mpz_poly_reverse(F_poly2, F_poly, length2);
		F_mpz_poly_reverse(F_poly2, F_poly2, length2);
           
      result = F_mpz_poly_equal(F_poly2, F_poly);
      if (!result)
		{
			printf("Error: length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
			printf("F_poly->length = %ld, F_poly2->length = %ld\n", F_poly->length, F_poly2->length);
		}

      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   // completely aliased reverse and reverse back
	for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);   
      
      length = z_randint(100);
      length2 = length + z_randint(200);
       
      F_mpz_randpoly(F_poly, length, bits); 
         
      F_mpz_poly_set(F_poly2, F_poly);
      F_mpz_poly_reverse(F_poly, F_poly, length2);
      F_mpz_poly_reverse(F_poly, F_poly, length2);
           
      result = F_mpz_poly_equal(F_poly2, F_poly);
		if (!result)
		{
			printf("Error: length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
			printf("F_poly->length = %ld, F_poly2->length = %ld\n", F_poly->length, F_poly2->length);
		}
      
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }

   return result; 
}

int test_F_mpz_poly_shift()
{
   F_mpz_poly_t F_poly, F_poly2, F_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long shift;
   
   F_mpz_poly_init(F_poly);
   F_mpz_poly_init(F_poly2);
   F_mpz_poly_init(F_poly3);
      
	// left shift followed by right shift
	for (unsigned long count1 = 0; (count1 < 700) && (result == 1) ; count1++)
   {
      bits = z_randint(500)+ 1;
      length = z_randint(20);        
      shift = z_randint(20);
		    
		F_mpz_randpoly(F_poly, length, bits); 
      F_mpz_poly_set(F_poly3, F_poly);
      
		F_mpz_poly_left_shift(F_poly2, F_poly, shift); 
      F_mpz_poly_right_shift(F_poly, F_poly2, shift);
      
      result = F_mpz_poly_equal(F_poly3, F_poly);
		if (!result) printf("Error: bits = %ld, length = %ld, shift = %ld, F_poly->length = %ld, F_poly3->length = %ld\n", bits, length, shift, F_poly->length, F_poly3->length);
   }

	// F_poly3 is used uninitialised after this point so clear it
	F_mpz_poly_clear(F_poly3);

   // left shift followed by right shift, completely aliased
	for (unsigned long count1 = 0; (count1 < 700) && (result == 1) ; count1++)
   {
      bits = z_randint(500)+ 1;
      length = z_randint(500);        
      shift = z_randint(100);
		    
		F_mpz_randpoly(F_poly, length, bits); 
      F_mpz_poly_set(F_poly2, F_poly);
      
		F_mpz_poly_left_shift(F_poly, F_poly, shift); 
      F_mpz_poly_right_shift(F_poly, F_poly, shift);
      
      result = F_mpz_poly_equal(F_poly2, F_poly);
		if (!result) printf("Error: bits = %ld, length = %ld, shift = %ld\n", bits, length, shift);
   }
   
   // explicit check of right shift
	for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
   { 
       bits = z_randint(500)+ 1;
       length = z_randint(500);        
       if (length) shift = z_randint(length);
		 else shift = 0;
		
       do F_mpz_randpoly(F_poly, length, bits); 
       while (F_poly->length < length);

       F_poly3->length = F_poly->length - shift;
       F_poly3->coeffs = F_poly->coeffs + shift;
		 
		 F_mpz_poly_right_shift(F_poly2, F_poly, shift);      
          
		 result = F_mpz_poly_equal(F_poly3, F_poly2);
		 if (!result) printf("Error: bits = %ld, length = %ld, shift = %ld\n", bits, length, shift);
   }

   F_mpz_poly_clear(F_poly);
   F_mpz_poly_clear(F_poly2);
  
   return result; 
}

int test_F_mpz_poly_scalar_mul_ui()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, length;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_poly_scalar_mul_ui(F_poly2, F_poly, mult);
          
      mpz_poly_init(m_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mult == 0L) result = (F_poly2->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased multiply
	for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
          
		F_mpz_poly_scalar_mul_ui(F_poly, F_poly, mult);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      if (mult == 0L) result = (F_poly->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_clear(temp);
   
   return result; 
}

int test_F_mpz_poly_scalar_mul_si()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, bits2, length;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_poly_scalar_mul_si(F_poly2, F_poly, mult);
          
      mpz_poly_init(m_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mult == 0L) result = (F_poly2->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            if (mult < 0L) 
			   {
				   mpz_mul_ui(temp, m_poly->coeffs[i], -mult);
				   mpz_neg(temp, temp);
			   } else mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased multiply
	for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
          
		F_mpz_poly_scalar_mul_si(F_poly, F_poly, mult);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      if (mult == 0L) result = (F_poly->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            if (mult < 0L) 
			   {
				   mpz_mul_ui(temp, m_poly->coeffs[i], -mult);
				   mpz_neg(temp, temp);
			   } else mpz_mul_ui(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }          
	   }
		
	   if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_clear(temp);
   
   return result; 
}

int test_F_mpz_poly_scalar_mul()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
	F_mpz_t x;
   int result = 1;
   ulong bits, bits2, length;
   mpz_t temp, mult;
   mpz_init(temp);
   mpz_init(mult);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
		F_mpz_init(x);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(200);
		mpz_rrandomb(mult, randstate, bits2);
		if (z_randint(2)) mpz_neg(mult, mult);

		F_mpz_set_mpz(x, mult);
      
		F_mpz_poly_scalar_mul(F_poly2, F_poly, x);
          
      mpz_poly_init(m_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mpz_sgn(mult) == 0) result = (F_poly2->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

		if (!result) 
		{
			gmp_printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %Zd\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_clear(x);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased multiply
	for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      F_mpz_init(x);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      bits2 = z_randint(200);
		mpz_rrandomb(mult, randstate, bits2);
		if (z_randint(2)) mpz_neg(mult, mult);
          
		F_mpz_set_mpz(x, mult);
      
		F_mpz_poly_scalar_mul(F_poly, F_poly, x);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      if (mpz_sgn(mult) == 0) result = (F_poly->length == 0);
		else
		{
			for (ulong i = 0; i < m_poly->length; i++)
         {
            mpz_mul(temp, m_poly->coeffs[i], mult);
            result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
         }
		}

	   if (!result) 
		{
			gmp_printf("Error: length = %ld, bits = %ld, bits2 = %ld, mult = %Zd\n", length, bits, bits2, mult);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_clear(x);
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_clear(temp);
   mpz_clear(mult);
   
   return result; 
}

int test_F_mpz_poly_mul_classical()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
		length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_classical(res, F_poly1, F_poly2);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: m_poly1->length = %ld, bits1 = %ld, m_poly2->length = %ld, bits2 = %ld\n", m_poly1->length, bits1, m_poly2->length, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_classical(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_classical(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of poly1 and poly2
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_classical(res, F_poly1, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive(res1, m_poly1, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_mul_karatsuba()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
		mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_karatsuba(res, F_poly1, F_poly2);			
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_karatsuba(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_karatsuba(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_karatsuba(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_karatsuba(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_karatsuba(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of poly1 and poly2
	for (ulong count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_karatsuba(res, F_poly1, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
		mpz_poly_set(m_poly2, m_poly1);
      mpz_poly_mul_karatsuba(res1, m_poly2, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_bit_pack()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   mp_limb_t * array;
	int result = 1;
   ulong bits, length, depth;
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 0; (count1 < 500) && (result == 1) ; count1++)
   {
      bits = z_randint(FLINT_BITS - 4) + 2;
      
      length = z_randint(1000)+1;
      
		F_mpz_poly_init2(F_poly, length);
      F_mpz_poly_init2(F_poly2, length);
   
      do mpz_randpoly(m_poly, length, bits - 1);
      while (m_poly->length < length);
          
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
      ulong n = (bits*length - 1)/FLINT_BITS + 1;
		array = flint_heap_alloc(n);
          
      F_mpz_poly_bit_pack(array, n, F_poly, bits, length, 0L);
      F_poly2->length = length;
          
      F_mpz_poly_bit_unpack(F_poly2, array, length, bits); 
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
      flint_heap_free(array);          
          
		result = mpz_poly_equal(m_poly, m_poly2);
      if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
      
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   // try negating the coefficients
	for (ulong count1 = 0; (count1 < 500) && (result == 1) ; count1++)
   {
      bits = z_randint(FLINT_BITS - 4) + 2;
      length = z_randint(1000)+1;
      
      F_mpz_poly_init2(F_poly, length);
      F_mpz_poly_init2(F_poly2, length);
    
      do mpz_randpoly(m_poly, length, bits - 1);
      while (m_poly->length < length);
          
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
      ulong n = (bits*length - 1)/FLINT_BITS + 1;
		array = flint_heap_alloc(n);
          
      F_mpz_poly_bit_pack(array, n, F_poly, bits, length, -1L);
      F_poly2->length = length;
          
      F_mpz_poly_bit_unpack(F_poly2, array, length, bits); 
      F_mpz_poly_neg(F_poly2, F_poly2);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
      flint_heap_free(array);          
          
		result = mpz_poly_equal(m_poly, m_poly2);
      if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
      
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_bit_pack_unsigned()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   mp_limb_t * array;
	int result = 1;
   ulong bits, length, depth;
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 0; (count1 < 500) && (result == 1) ; count1++)
   {
      bits = z_randint(FLINT_BITS - 4) + 2;
      
      length = z_randint(1000)+1;
      
		F_mpz_poly_init2(F_poly, length);
      F_mpz_poly_init2(F_poly2, length);
   
      do mpz_randpoly_unsigned(m_poly, length, bits - 1);
      while (m_poly->length < length);
          
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
      ulong n = (bits*length - 1)/FLINT_BITS + 1;
		array = flint_heap_alloc(n);
          
      F_mpz_poly_bit_pack_unsigned(array, n, F_poly, bits, length);
      F_poly2->length = length;
                
      F_mpz_poly_bit_unpack_unsigned(F_poly2, array, length, bits); 
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
      flint_heap_free(array);          
          
		result = mpz_poly_equal(m_poly, m_poly2);
      if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}
      
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_bit_pack2()
{
   mpz_poly_t m_poly, m_poly2, m_poly3;
   F_mpz_poly_t F_poly, F_poly2, F_poly3;
   mp_limb_t * array1, * array2;
	int result = 1;
   ulong bits, length, depth, bundle;
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(m_poly3); 
   
   for (ulong count1 = 0; (count1 < 500) && (result == 1) ; count1++)
   {
      bits = z_randint(FLINT_BITS-4)+ 2;
      
      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      
      length = z_randint(1000)+1;
      
		F_mpz_poly_fit_length(F_poly, length);
      F_mpz_poly_fit_length(F_poly2, length);
      F_mpz_poly_fit_length(F_poly3, length);
         
      do mpz_randpoly(m_poly, length, bits - 1);
      while (m_poly->length < length);
      
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
      
      ulong n = (bits*length - 1)/FLINT_BITS + 1;
		array1 = flint_heap_alloc(n);
      array2 = flint_heap_alloc(n);
          
      F_mpz_poly_bit_pack2(array1, array2, n, F_poly, bits, length, 0L, 0L);
             
      F_mpz_poly_bit_unpack(F_poly2, array1, length, bits); 
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
      F_mpz_poly_bit_unpack(F_poly3, array2, length, bits); 
      F_mpz_poly_to_mpz_poly(m_poly3, F_poly3);

		for (ulong i = 0; i < m_poly3->length; i++)
			if (i%2 == 1) mpz_neg(m_poly3->coeffs[i], m_poly3->coeffs[i]);
          
      flint_heap_free(array1);          
      flint_heap_free(array2);          
          
		result &= mpz_poly_equal(m_poly, m_poly2);
      result &= mpz_poly_equal(m_poly, m_poly3);

      if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n\n");
         mpz_poly_print_pretty(m_poly3, "x"); printf("\n\n");
		}
      
      F_mpz_poly_clear(F_poly);
      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly3);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   mpz_poly_clear(m_poly3);
   
   return result;
}

int test_F_mpz_poly_mul_KS()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
		mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_KS(res, F_poly1, F_poly2);			
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive_KS(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// try unsigned coefficients
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly_unsigned(m_poly1, length1, bits1);
		mpz_randpoly_unsigned(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_KS(res, F_poly1, F_poly2);			
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive_KS(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_KS(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_karatsuba(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_KS(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive_KS(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of poly1 and poly2
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_KS(res, F_poly1, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
		mpz_poly_set(m_poly2, m_poly1);
      mpz_poly_mul_naive_KS(res1, m_poly2, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_mul_KS2()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
		mpz_randpoly(m_poly1, length1, bits1);
		mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_KS2(res, F_poly1, F_poly2);			
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive_KS(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// try unsigned coefficients
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly_unsigned(m_poly1, length1, bits1);
		mpz_randpoly_unsigned(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_KS2(res, F_poly1, F_poly2);			
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive_KS(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(res);
   }
   
	// test aliasing of res and poly1
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_KS2(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_karatsuba(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
   // test aliasing of res and poly2
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      bits2 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_KS2(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_naive_KS(res1, m_poly1, m_poly2);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
	
	// test aliasing of poly1 and poly2
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint((FLINT_BITS - 9)/2) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_KS2(res, F_poly1, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
		mpz_poly_set(m_poly2, m_poly1);
      mpz_poly_mul_naive_KS(res1, m_poly2, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(res);
   }
   
	mpz_poly_clear(res1);
   mpz_poly_clear(res2);
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

void F_mpz_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(F_mpz_poly_convert); 
   RUN_TEST(F_mpz_poly_getset_coeff_si); 
   RUN_TEST(F_mpz_poly_getset_coeff_ui); 
   RUN_TEST(F_mpz_poly_getset_coeff_mpz); 
   RUN_TEST(F_mpz_poly_set); 
   RUN_TEST(F_mpz_poly_equal); 
   RUN_TEST(F_mpz_poly_swap); 
	RUN_TEST(F_mpz_poly_max_bits1);
   RUN_TEST(F_mpz_poly_max_bits);
   RUN_TEST(F_mpz_poly_max_limbs); 
   RUN_TEST(F_mpz_poly_neg);
	RUN_TEST(F_mpz_poly_reverse); 
   RUN_TEST(F_mpz_poly_add); 
   RUN_TEST(F_mpz_poly_sub); 
   RUN_TEST(F_mpz_poly_shift); 
   RUN_TEST(F_mpz_poly_scalar_mul_ui); 
   RUN_TEST(F_mpz_poly_scalar_mul_si); 
   RUN_TEST(F_mpz_poly_scalar_mul); 
   RUN_TEST(F_mpz_poly_mul_classical); 
   RUN_TEST(F_mpz_poly_mul_karatsuba); 
	RUN_TEST(F_mpz_poly_bit_pack);
   RUN_TEST(F_mpz_poly_bit_pack_unsigned);
   RUN_TEST(F_mpz_poly_bit_pack2);
   RUN_TEST(F_mpz_poly_mul_KS); 
	RUN_TEST(F_mpz_poly_mul_KS2);
	
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();
	_F_mpz_cleanup();

   return 0;
}
