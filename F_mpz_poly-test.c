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
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		    
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
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		    
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
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		    
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
      mpz_poly_mul_classical(res1, m_poly1, m_poly1);		
		    
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
			printf("Error2: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
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
			printf("Error3: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
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
			printf("Error4: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
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
			printf("Error5: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
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

		bits1 = z_randint(200) + 4;
      bits2 = z_randint(200) + 4;
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

int test_F_mpz_poly_byte_pack_unsigned()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2;
   mp_limb_t * array;
   int result = 1;
   ulong bits, length, bytes;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 1; (count1 < 200*ITER) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 2;
      bytes = ((bits - 1)>>3) + 1;
      
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      for (ulong count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = z_randint(1000) + 1;
#if DEBUG
          printf("%ld, %ld\n", length, bits);
#endif
          F_mpz_poly_fit_length(F_poly1, length);
          F_mpz_poly_fit_length(F_poly2, length);
          
          do mpz_randpoly_unsigned(m_poly1, length, bits/2);
          while (mpz_poly_length(m_poly1) < length);

#if DEBUG
          for (ulong j = 0; j < m_poly1->length; j++)
             gmp_printf("%Zx, ", m_poly1->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
          array = flint_heap_alloc(((bytes*length-1)>>FLINT_LG_BYTES_PER_LIMB) + 2);
          
          F_mpz_poly_byte_pack(array, F_poly1, length, bytes, 1L);
             
		    F_mpz_poly_zero(F_poly2);
		    F_mpz_poly_fit_length(F_poly2, length);
			 F_poly2->length = length;
			 F_mpn_clear(F_poly2->coeffs, length);
          
          F_mpz_poly_byte_unpack_unsigned(F_poly2, array, length, bytes);  
          
			 F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);
          
          flint_heap_free(array);
          
#if DEBUG
          for (ulong j = 0; j < m_poly2->length; j++)
             gmp_printf("%Zx, ", m_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(m_poly1, m_poly2);
			 if (!result)
			 {
				 mpz_poly_print(m_poly1); printf("\n\n");
				 mpz_poly_print(m_poly2); printf("\n\n");
			 }
      }   
      
		F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_byte_pack()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2;
   mp_limb_t * array;
   int result = 1;
   ulong bits, length, bytes;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 1; (count1 < 200*ITER) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 2;
      bytes = ((bits - 1)>>3) + 1;
      
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      for (ulong count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = z_randint(1000) + 1;
#if DEBUG
          printf("%ld, %ld\n", length, bits);
#endif
          F_mpz_poly_fit_length(F_poly1, length);
          F_mpz_poly_fit_length(F_poly2, length);
          
          do mpz_randpoly(m_poly1, length, bits/2);
          while (mpz_poly_length(m_poly1) < length);

#if DEBUG
          for (ulong j = 0; j < m_poly1->length; j++)
             gmp_printf("%Zx, ", m_poly1->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
          array = flint_heap_alloc(((bytes*length-1)>>FLINT_LG_BYTES_PER_LIMB) + 2);
          
          F_mpz_poly_byte_pack(array, F_poly1, length, bytes, 1L);
             
		    F_mpz_poly_zero(F_poly2);
		    F_mpz_poly_fit_length(F_poly2, length);
			 F_poly2->length = length;
			 F_mpn_clear(F_poly2->coeffs, length);
          
          F_mpz_poly_byte_unpack(F_poly2, array, length, bytes);  
          
			 F_mpz_poly_to_mpz_poly(m_poly2, F_poly2);

          flint_heap_free(array);
          
#if DEBUG
          for (ulong j = 0; j < m_poly2->length; j++)
             gmp_printf("%Zx, ", m_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(m_poly1, m_poly2);
			 if (!result)
			 {
				 mpz_poly_print(m_poly1); printf("\n\n");
				 mpz_poly_print(m_poly2); printf("\n\n");
			 }
      }   
      
		F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
   }
   
   mpz_poly_clear(m_poly1);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_mul_SS()
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

		bits1 = z_randint(200) + 4;
      bits2 = z_randint(200) + 4;
      length1 = z_randint(100);
      length2 = z_randint(100);
      
		mpz_randpoly(m_poly1, length1, bits1);
		mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_SS(res, F_poly1, F_poly2);			
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
      
		F_mpz_poly_mul_SS(res, F_poly1, F_poly2);			
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
      
		F_mpz_poly_mul_SS(res, res, F_poly1);
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
      
		F_mpz_poly_mul_SS(res, F_poly1, res);
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
      
		F_mpz_poly_mul_SS(res, F_poly1, F_poly1);
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

int test_F_mpz_poly_pack_bytes()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly1, F_poly2, F_poly3;
   int result = 1;
   ulong bits, length, bytes;
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 

   for (ulong count1 = 1; (count1 < 10*ITER) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ FLINT_BITS;
      bytes = ((bits-1)>>3)+1;
      
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      
	  length = z_randint(1000) + 80;
	  ulong n = z_randint(length/5 + 1) + 1;
	  ulong limbs = ((2*n - 1)*bytes*8 - 1)/FLINT_BITS + 1;
			 
#if DEBUG
      printf("%ld, %ld, %ld, %ld, %ld\n", length, bits, n, limbs, bytes);
#endif

      do mpz_randpoly(m_poly, length, bits/2);
      while (m_poly->length < length);

#if DEBUG
      for (unsigned j = 0; j < m_poly->length; j++)
         gmp_printf("%Zx, ", m_poly->coeffs[j]);
      printf("\n\n");
#endif

      mpz_poly_to_F_mpz_poly(F_poly1, m_poly);
          
      F_mpz_poly_pack_bytes(F_poly2, F_poly1, n, bytes);
	  F_mpz_poly_unpack_bytes(F_poly3, F_poly2, n, bytes);
          
	  F_mpz_poly_to_mpz_poly(m_poly2, F_poly3);
                   
#if DEBUG
	  for (unsigned j = 0; j < m_poly2->length; j++)
         gmp_printf("%Zx, ", m_poly2->coeffs[j]);
      printf("\n\n");
#endif
          
      result = F_mpz_poly_equal(F_poly1, F_poly3);
      if (!result) 
      {
	     mpz_poly_print_pretty(m_poly, "x"); printf("\n\n");
		 mpz_poly_print_pretty(m_poly2, "x"); printf("\n\n");
      }   

      F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly3);
   }
   
   mpz_poly_clear(m_poly);
   mpz_poly_clear(m_poly2);
   
   return result;
}

int test_F_mpz_poly_tofromstring()
{
   mpz_poly_t test_poly;
   F_mpz_poly_t test_F_mpz_poly, test_F_mpz_poly2;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = z_randint(100)+ 1;
      
      F_mpz_poly_init2(test_F_mpz_poly, (bits-1)/FLINT_BITS+1);
      F_mpz_poly_init2(test_F_mpz_poly2, (bits-1)/FLINT_BITS+1+z_randint(30));
      FILE * testfile; 
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = z_randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          F_mpz_poly_fit_length(test_F_mpz_poly, length);
          F_mpz_poly_fit_length(test_F_mpz_poly2, length);
          mpz_randpoly(test_poly, length, bits); 

          mpz_poly_to_F_mpz_poly(test_F_mpz_poly, test_poly);
          
          char * strbuf = F_mpz_poly_to_string(test_F_mpz_poly);
          int OK = F_mpz_poly_from_string(test_F_mpz_poly2, strbuf);
          free(strbuf);
//          F_mpz_poly_check_normalisation(test_F_mpz_poly2);
          result = F_mpz_poly_equal(test_F_mpz_poly2, test_F_mpz_poly) && OK;
           
      }
            
      F_mpz_poly_clear(test_F_mpz_poly);
      F_mpz_poly_clear(test_F_mpz_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_F_mpz_poly_to_zmod_poly()
{
   mpz_poly_t test_poly;
   F_mpz_t temp;
   F_mpz_poly_t test_F_mpz_poly, test_F_mpz_poly2;
	zmod_poly_t test_zmod_poly;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly);
   F_mpz_init(temp); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = z_randint(FLINT_BITS-3) + 1;
      
      F_mpz_poly_init(test_F_mpz_poly);
      F_mpz_poly_init(test_F_mpz_poly2);

      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
          length = z_randint(20);
			 ulong p = z_nextprime(1L<<(bits+1), 0);
			 zmod_poly_init(test_zmod_poly, p);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          F_mpz_poly_fit_length(test_F_mpz_poly, length);
          mpz_randpoly(test_poly, length, bits);
           
#if DEBUG
          mpz_poly_print_pretty(test_poly, "x");
          printf("\n\n");
#endif
          mpz_poly_to_F_mpz_poly(test_F_mpz_poly, test_poly);
//          F_mpz_poly_check_normalisation(test_F_mpz_poly);
          
			 F_mpz_poly_to_zmod_poly(test_zmod_poly, test_F_mpz_poly);
			 zmod_poly_to_F_mpz_poly(test_F_mpz_poly2, test_zmod_poly);

// Now reduce poly mod p

          for (unsigned long c = 0; c < test_F_mpz_poly->length; c++){
             F_mpz_sub(temp, test_F_mpz_poly->coeffs+c, test_F_mpz_poly2->coeffs+c);
             result = (F_mpz_mod_ui(temp, temp, p) == 0L);
             if (!result)
               break;
          }

			 if (!result)
			 {
				 F_mpz_poly_print(test_F_mpz_poly); printf("\n\n");
				 F_mpz_poly_print(test_F_mpz_poly2); printf("\n\n");
				 printf("p = %ld\n", p);
			 }

			 zmod_poly_clear(test_zmod_poly);
      }   
          
      F_mpz_poly_clear(test_F_mpz_poly);
      F_mpz_poly_clear(test_F_mpz_poly2);
   }

   F_mpz_clear(temp);  
   mpz_poly_clear(test_poly);
   
   return result;
}


/******************************************************************************

Square-Free Factorization

******************************************************************************/

void F_mpz_poly_squarefree(F_mpz_poly_factor_t fac, F_mpz_t content, F_mpz_poly_t F)
{

   F_mpz_poly_content(content, F);

   F_mpz_poly_t f;

   F_mpz_poly_init(f);

   F_mpz_poly_scalar_div_exact(f, F, content);

   F_mpz_poly_factor_clear(fac);

   F_mpz_poly_factor_init(fac);

   if (f->length == 1)
      return;

   F_mpz_poly_t d, v, w, s, t1;
   F_mpz_poly_init(d);
   F_mpz_poly_init(v);
   F_mpz_poly_init(w);
   F_mpz_poly_init(s);
   F_mpz_poly_init(t1);


   F_mpz_poly_derivative(t1, f);
   F_mpz_poly_gcd(d, f, t1);

   if (d->length == 1){
      F_mpz_poly_factor_insert(fac, f, 1);

      F_mpz_poly_clear(d);
      F_mpz_poly_clear(v);
      F_mpz_poly_clear(w);
      F_mpz_poly_clear(s);
      F_mpz_poly_clear(t1);
      F_mpz_poly_clear(f);

      return;
   }

   F_mpz_poly_div(v, f, d);
   F_mpz_poly_div(w, t1, d);

   long i = 0;

   for ( ; ;){

      i = i + 1;


      F_mpz_poly_derivative(t1, v);
      F_mpz_poly_sub(s, w, t1);

      if (s->length == 0){
         if (v->length > 1)
            F_mpz_poly_factor_insert(fac, v, i);

         F_mpz_poly_clear(d);
         F_mpz_poly_clear(v);
         F_mpz_poly_clear(w);
         F_mpz_poly_clear(s);
         F_mpz_poly_clear(t1);
         F_mpz_poly_clear(f);
         return;
      }

      F_mpz_poly_gcd(d, v, s);
      F_mpz_poly_div(v, v, d);
      F_mpz_poly_div(w, s, d);

      if (d->length > 1)
         F_mpz_poly_factor_insert(fac, d, i);
   } 

   F_mpz_poly_clear(f);
   F_mpz_poly_clear(d);
   F_mpz_poly_clear(v);
   F_mpz_poly_clear(w);
   F_mpz_poly_clear(s);
   F_mpz_poly_clear(t1);

   return;
}

/******************************************************************************

Naive modp F_mpz functions

******************************************************************************/


void F_mpz_poly_rem_modp_naive(F_mpz_poly_t R, F_mpz_poly_t A, F_mpz_poly_t B, F_mpz_t p){

//Here we want to work as if F_mpz_poly's are modded out by an F_mpz p.  Perhaps we assume that A,B are already reduced by p...

   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      F_mpz_poly_set(R, A);      
      return;
   }

   F_mpz_poly_set(R,A);
   long coeff = A->length - 1;

   F_mpz_poly_t pre_inv_B, Bm1, qB;

   F_mpz_poly_init2(Bm1, B->length);

   F_mpz_poly_set(Bm1, B);
//   F_mpz_poly_init2(qB, Bm1->length);

   //need a copy of B for the many multiplication/subtractions

   //Also want to precompute 1/leadcoeff(B) mod P

   F_mpz_t B_lead_inv;
   F_mpz_init(B_lead_inv);
   F_mpz_t coeff_Q;
   F_mpz_init( coeff_Q);

   F_mpz_invert(B_lead_inv, B->coeffs + (B->length - 1), p);

   long R_length;

   F_mpz_t temp;
   F_mpz_init( temp );

   while (coeff >= (long) Bm1->length - 1){

      while( (coeff >= (long) Bm1->length - 1) && F_mpz_is_zero( R->coeffs + coeff ) ){
         coeff--;
      }

      if (coeff >= (long) B->length - 1){

         //coeff_Q here
         F_mpz_mul2(temp, R->coeffs + coeff, B_lead_inv);
         F_mpz_mod(coeff_Q, temp, p);

         F_mpz_poly_init(qB);

         F_mpz_poly_fit_length(qB, Bm1->length);

         qB->length = Bm1->length;

         for (long i = 0; i < Bm1->length; i++){
            F_mpz_mul2(temp, Bm1->coeffs + i , coeff_Q );
            F_mpz_mod(qB->coeffs + i, temp, p);
         }
 
         //for each coeff do mulmod by coeff_Q and write to temp_B

         // shift qB but left_shift doesn't exist it's now just shift I think

         F_mpz_poly_left_shift(qB, qB, coeff - B->length + 1);

         R_length = coeff;

         F_mpz_poly_sub(R, R, qB);

         R->length = R_length;

         F_mpz_poly_clear(qB);

         for (long i = 0; i < R->length; i++){
            F_mpz_mod(R->coeffs + i, R->coeffs + i, p);
         }


         //subtract a shifted temp_B from R

         coeff--;

      }

   }

   R->length = B->length - 1;
   _F_mpz_poly_normalise(R);

   F_mpz_clear(B_lead_inv);

   F_mpz_poly_clear(Bm1);
   F_mpz_clear(coeff_Q);
   F_mpz_clear(temp);


}


void F_mpz_poly_mulmod_modp_naive(F_mpz_poly_t R, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t B, F_mpz_t p){

//Want to compute f*g mod B where all polynomials are considered modulo a F_mpz modulus p.  Need that the leading coeff of H is invertible mod p.

   if (B->length == 0){
      printf("FLINT Exception: Divide by zero\n");
      abort();
   }

   if ( (B->length ==1) || (f->length == 0) || (g->length == 0) ){
      F_mpz_poly_zero(R);
      return;
   }

   F_mpz_poly_t prod;
   F_mpz_poly_init(prod);
   F_mpz_poly_mul(prod, f, g);
   F_mpz_poly_rem_modp_naive(R, prod, B, p);
   F_mpz_poly_clear(prod);

}




/******************************************************************************************************************************************************

Hensel Lifting procedures

****************************************************************************************************************************************************/



void _Build_Hensel_Tree(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, zmod_poly_factor_t fac){

   long r;
   unsigned long p;

   r = fac->num_factors;
   p = (fac->factors[0])->p;

   long i, j, s;
   long minp, mind;
   long tmp;

   zmod_poly_t V[2*r-2];
   zmod_poly_t W[2*r-2];

   for (i = 0; i < 2*r-2; i++)
      zmod_poly_init(V[i],p);

   for (i = 0; i < 2*r-2; i++)
      zmod_poly_init(W[i],p);



/* We will have five arrays: a v of fmpz_polys and a V of zmod_polys also a w and a W and link.  Here's the idea, we will sort each leaf and node of a factor tree by degree, in fact choosing to multiply the two smallest factors, then the next two smallest (factors or products) until a tree is made.  The tree will be stored in the v's.  The first two elements of v will be the smallest modular factors, the last two elements of v will multiply to form F itself.  Since v will be rearranging the original factors we will need to be able to recover the original order.  For this we use link which has nonnegative even numbers and negative numbers.  Link is an array of longs which aligns with V/v  if link has a negative number in spot j that means V[j] is an original modular factor which has been lifted, if link[j] is a nonnegative even number then V[j] stores a product of the two entries at V[ link[j] ] and V[ link[j]+1 ].  W/w plays the role of the extended GCD, at V[0], V[2], V[4], etc we have a new product, W[0], W[2], W[4], etc are the XGCD compliments of the V's.  So V[0]*W[0]+V[1]*W[1] = 1 mod p^(something)  these will be lifted along with the entries in V.  It's not enough to just lift each factor we have to lift the entire tree and the tree of XGCD inverses.*/

   for (i = 0; i < r; i++){
      zmod_poly_set(V[i], fac->factors[i]);
      link[i] = -(i+1);
   }

   for (j = 0; j < 2*r - 4; j += 2){

      minp = j;
      mind = zmod_poly_degree(V[j]);

      for (s = j+1; s < i; s++){
         if (zmod_poly_degree(V[s]) < mind){
            minp = s;
            mind = zmod_poly_degree(V[s]);
         }
      }

      zmod_poly_swap(V[j],V[minp]);

      tmp = link[j];
      link[j] = link[minp];
      link[minp] = tmp; 
      //swap link[j] and V[minp]

      minp = j+1;
      mind = zmod_poly_degree(V[j+1]);

      for ( s = j+2; s < i; s++){
         if (zmod_poly_degree(V[s]) < mind ){
            minp = s;
            mind = zmod_poly_degree(V[s]);
         }
      }

      zmod_poly_swap(V[j+1],V[minp]);

      tmp = link[j+1];
      link[j+1] = link[minp];
      link[minp] = tmp; 
      //swap link[j+1] and V[minp]

      zmod_poly_mul(V[i], V[j], V[j+1]);
      link[i] = j;
      i++;

   }

   zmod_poly_t d;
   zmod_poly_init(d, p);

   for (j = 0; j < 2*r - 2; j += 2){
      zmod_poly_xgcd_euclidean(d, W[j], W[j+1], V[j], V[j+1]);
      //Make a check for d!=1
   }

   for (j = 0; j < 2*r-2; j++){
      zmod_poly_to_F_mpz_poly(v[j], V[j]);
      zmod_poly_to_F_mpz_poly(w[j], W[j]);
   }

   for (i = 0; i < 2*r-2; i++)
      zmod_poly_clear(V[i]);

   for (i = 0; i < 2*r-2; i++)
      zmod_poly_clear(W[i]);

   zmod_poly_clear(d);

}


//const long p should be fmpz_t p, but right now we dont support giant moduli

void _Hensel_Lift(F_mpz_poly_t Gout, F_mpz_poly_t Hout, F_mpz_poly_t Aout, F_mpz_poly_t Bout, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t h, F_mpz_poly_t a, F_mpz_poly_t b, F_mpz_t p, F_mpz_t p1){

   F_mpz_poly_t c, g1, h1, G, H, A, B;
   F_mpz_poly_init(c);
   F_mpz_poly_init(g1);
   F_mpz_poly_init(h1);
   F_mpz_poly_init(G);
   F_mpz_poly_init(H);
   F_mpz_poly_init(A);
   F_mpz_poly_init(B);

   F_mpz_poly_mul(c, g, h);

   F_mpz_poly_sub(c, f, c);

   F_mpz_poly_scalar_div_exact(c, c, p);
   //Make a check that c is divisible by p

//When I make a precomputing function use GG, HH instead

   F_mpz_poly_mulmod_modp_naive(h1, c, a, h, p1);
   F_mpz_poly_mulmod_modp_naive(g1, c, b, g, p1);
//   F_mpz_poly_smod(g1, g1, p1);
//   F_mpz_poly_smod(h1, h1, p1);

   F_mpz_poly_scalar_mul(g1, g1, p);
   F_mpz_poly_scalar_mul(h1, h1, p);

   F_mpz_poly_add(G, g, g1);
   F_mpz_poly_add(H, h, h1);
   
//Lifting the inverses now

   F_mpz_poly_t a1, b1, t1, t2, unity;
   F_mpz_poly_init(a1);
   F_mpz_poly_init(b1);
   F_mpz_poly_init(t1);
   F_mpz_poly_init(t2);
   F_mpz_poly_init(unity);

   F_mpz_poly_set_coeff_si(unity, 0, -1);

   F_mpz_poly_mul(t1, a, G);
   F_mpz_poly_mul(t2, b, H);

   F_mpz_poly_add(t1, t1, t2);
   F_mpz_poly_add(t1, t1, unity);
   F_mpz_poly_neg(t1, t1);

   F_mpz_poly_scalar_div_exact(t1, t1, p);
//Make a check that t1 is divisible by p

   F_mpz_poly_mulmod_modp_naive(a1, t1, a, h, p1);
   F_mpz_poly_mulmod_modp_naive(b1, t1, b, g, p1);
//   F_mpz_poly_smod(b1, b1, p1);
//   F_mpz_poly_smod(a1, a1, p1);

   F_mpz_poly_scalar_mul(a1, a1, p);
   F_mpz_poly_add(A, a, a1);

   F_mpz_poly_scalar_mul(b1, b1, p);
   F_mpz_poly_add(B, b, b1);

   F_mpz_poly_set(Gout, G);
   F_mpz_poly_set(Hout, H);
   F_mpz_poly_set(Aout, A);
   F_mpz_poly_set(Bout, B);

   F_mpz_poly_clear(a1);
   F_mpz_poly_clear(b1);
   F_mpz_poly_clear(t1);
   F_mpz_poly_clear(t2);
   F_mpz_poly_clear(unity);

   F_mpz_poly_clear(c);
   F_mpz_poly_clear(g1);
   F_mpz_poly_clear(h1);
   F_mpz_poly_clear(G);
   F_mpz_poly_clear(H);
   F_mpz_poly_clear(A);
   F_mpz_poly_clear(B);

}




void _Rec_Tree_Hensel_Lift(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, F_mpz_t p, F_mpz_poly_t f, long j, long inv, F_mpz_t p1){

   if (j < 0) return;

   if (inv)
      _Hensel_Lift(v[j], v[j+1], w[j], w[j+1], f, v[j], v[j+1], w[j], w[j+1], p, p1);
//   else
//      _Hensel_Lift1(v[j], v[j+1], f, v[j], v[j+1], w[j], w[j+1], p, p1);
//altered to check a bug, should be Hensel_Lift1

   _Rec_Tree_Hensel_Lift(link, v, w, p, v[j],   link[j],   inv, p1);
   _Rec_Tree_Hensel_Lift(link, v, w, p, v[j+1], link[j+1], inv, p1);

}

void _Tree_Hensel_Lift(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, long e0, long e1, F_mpz_poly_t f, long inv, long p, long r, F_mpz_t P){

   //Want to make these into fmpz's before too long
   F_mpz_t temp, p0, p1;
   F_mpz_init(p0);
   F_mpz_init(p1);
   F_mpz_init(temp);

   F_mpz_set_ui(temp, p);

   F_mpz_pow_ui(p0, temp, e0);
   F_mpz_pow_ui(p1, temp, e1 - e0);

   _Rec_Tree_Hensel_Lift(link, v, w, p0, f, 2*r-4, inv, p1);

   F_mpz_mul2(P, p0, p1);

   F_mpz_clear(temp);
   F_mpz_clear(p0);
   F_mpz_clear(p1);

}

/***************************************************

Naive Zassenhaus

***************************/

void F_mpz_poly_zassenhaus_naive(F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc){

   ulong r = lifted_fac->num_factors;

   F_mpz_poly_t f;
   F_mpz_poly_init(f);
   F_mpz_poly_set(f, F);

   //Here is where we can check the factors
   F_mpz_poly_t Q,R;
   F_mpz_poly_init(Q);
   F_mpz_poly_init(R);
   F_mpz_poly_t tryme;
   F_mpz_poly_init(tryme);

   F_mpz_t temp_lc;

//F_mpz_poly_divrem(Q, R, f, tryme)
//if R==0 then insert tryme with exponents[j] and do some updating or something
//also replace f by Q then



   int k;
   int l;
   int indx;
   int used_arr[r];
   for(l = 0; l < r; l++){
      used_arr[l] = 0;
   }


   for (k = 1; k < r; k++){


      ulong count = 0;

      ulong sub_arr[k];

      for(l = 0; l < k; l++){
         sub_arr[l] = l;
      }

      indx = k-1;

      sub_arr[indx]--;

      while ((indx >= 0)){
         sub_arr[indx] = sub_arr[indx] + 1;
         for (l = indx + 1; l < k; l++){
            sub_arr[l] = sub_arr[l-1] + 1UL;
         }

         if (sub_arr[k-1] > r-1UL ){
            indx--;
         }
         else{


            for(l = 0; l < k; l++){
               if (used_arr[sub_arr[l]] == 1)
                  break;
            }
//Need to involve lc, perhaps set coeff 0 to lc and do lc * rest and check if under M_bits... here I'm using a trial division... hmm
            F_mpz_poly_fit_length(tryme, 1UL);
            tryme->length = 1UL;
            F_mpz_set(tryme->coeffs + 0, lc);
            for(l = 0; l < k; l++){
               F_mpz_poly_mul(tryme, tryme, lifted_fac->factors[sub_arr[l]]);
            }

            F_mpz_poly_smod(tryme, tryme, P);
            F_mpz_init(temp_lc);
            F_mpz_poly_content(temp_lc, tryme);
//FINDME
            F_mpz_poly_scalar_div_exact(tryme, tryme, temp_lc);

            F_mpz_poly_divrem(Q, R, f, tryme);

            if (R->length == 0){

//start here
               for(l = 0; l < k; l++){
                  printf("%d, ", sub_arr[l]);
               }
               printf("\n");
//      to here        FOUND ONE!!!!!
               F_mpz_poly_factor_insert(final_fac, tryme, exp);

//and here
               F_mpz_poly_print(tryme); printf(" tryme with lc\n");

               tryme->length = 1UL;
               F_mpz_set_ui(tryme->coeffs + 0, 1);
               for(l = 0; l < k; l++){
                  F_mpz_poly_mul(tryme, tryme, lifted_fac->factors[sub_arr[l]]);
               }

               F_mpz_poly_smod(tryme, tryme, P);
               F_mpz_poly_print(tryme); printf(" is factor without lc\n");

// to here
               
               for(l = 0; l < k; l++){
                  used_arr[sub_arr[l]] = 1;
                  count++;
               }
               F_mpz_poly_set(f, Q);
               F_mpz_set(lc, Q->coeffs + Q->length - 1 );

               //If r-count = k then the rest are irreducible.

            }

            F_mpz_clear(temp_lc);
            indx = k-1;
         }
      }

//This is where we switch to the next loop for k.
//So we will have found all factors using <= k local factors
//We should/could update f to be the rest divided away (or multiply the remaining)
// could also adjust r.  It is the number of remaining factors  
//so if you update then test if r = k or k+1 in which case the remaining f is irreducible.

   }

   ulong test = 0;

   for (l = 0; l < r; l++){
      test = test + used_arr[l];
   }

   if (test == 0)
      F_mpz_poly_factor_insert(final_fac, f, exp);

   F_mpz_poly_clear(f);
   F_mpz_poly_clear(tryme);
   F_mpz_poly_clear(Q);
   F_mpz_poly_clear(R);

   return;

}


/*********

   Factoring wrapper after square free part

*********/

void F_mpz_poly_factor_sq_fr_prim( F_mpz_poly_factor_t final_fac, ulong exp, F_mpz_poly_t f ){


   if (f->length <= 1){
      return;
   }

   if (f->length == 2){
      F_mpz_poly_factor_insert( final_fac, f, exp);
      return;
   }

   ulong len = f->length;

   F_mpz_t lc;

   F_mpz_init(lc);   

   F_mpz_set(lc, f->coeffs + len - 1);

   ulong M_bits = 0;

   M_bits = M_bits + F_mpz_bits(lc) + FLINT_ABS(F_mpz_poly_max_bits(f)) + len + (long)ceil(log2((double) len));

   zmod_poly_t F, F_d, F_sbo;

   int tryme = 1;

   ulong p = 2UL;

   long i;

   for (i = 0; (i < 200) && (tryme == 1); i++){
      zmod_poly_init(F, p);

      zmod_poly_init(F_d, p);

      zmod_poly_init(F_sbo, p);

      F_mpz_poly_to_zmod_poly(F, f);

      if (F->length < f->length){
         p = z_nextprime( p, 0);
         zmod_poly_clear(F_d);
         zmod_poly_clear(F_sbo);
         zmod_poly_clear(F);
         continue;
      }

      zmod_poly_derivative(F_d, F);      

      zmod_poly_gcd(F_sbo, F, F_d);

      if (zmod_poly_is_one(F_sbo)){
         tryme = 0;
      }
      else{
         p = z_nextprime( p, 0);

         zmod_poly_clear(F);
      }


      zmod_poly_clear(F_d);

      zmod_poly_clear(F_sbo);

   }

   if (i == 200){
      printf("wasn't square_free after 100 primes, maybe an error\n");

      zmod_poly_clear(F);

      F_mpz_clear(lc);

      return;
   }

   ulong a;

   a = (long) ceil( (double) M_bits / log2( (double)p ) );

   a = (long) pow( (double) 2, ceil( log2( (double) a ) ) );

   zmod_poly_factor_t fac;

   zmod_poly_factor_init(fac);

   unsigned long lead_coeff;

   lead_coeff = zmod_poly_factor(fac, F);

   long r;

   r = fac->num_factors;

   if (r > 20){

//In the future this is where we will call hoeij/Novocin method

      printf("r larger than 20, might take too long, I'm stoping\n");

      zmod_poly_clear(F);
      zmod_poly_factor_clear(fac);
      F_mpz_clear(lc);

      return;
   }

   if (r == 0){
      printf("FLINT-exception: something broke\n");
      zmod_poly_clear(F);
      zmod_poly_factor_clear(fac);
      F_mpz_clear(lc);
      abort();
   }

   if (r == 1){

      F_mpz_poly_factor_insert(final_fac, f, exp);
      F_mpz_poly_print(f); printf(" irreducible\n");

      zmod_poly_clear(F);
      zmod_poly_factor_clear(fac);
      F_mpz_clear(lc);

      return;
   }

//Begin Hensel Lifting phase, make the tree in v, w, and link

   F_mpz_poly_t v[2*r-2];
   F_mpz_poly_t w[2*r-2];   

   long link[2*r-2];

//P will be the F_mpz modulus

   F_mpz_t P, big_P;
   F_mpz_init(P);
   F_mpz_init(big_P);
   F_mpz_set_ui(P, p);
   F_mpz_pow_ui(big_P, P, a);

//Make a copy of f, f1, such that f1 monic and equiv to f mod p^a

   F_mpz_poly_t f1;

   F_mpz_poly_init(f1);

   F_mpz_t lc_inv;

   F_mpz_init( lc_inv );

   if (F_mpz_is_one(lc)){
      F_mpz_poly_set(f1, f);
   }
   else if (F_mpz_is_m1(lc)){
      F_mpz_poly_neg(f1, f);
   }
   else{
      F_mpz_mod(lc_inv, lc, big_P);

      int OK = F_mpz_invert(lc_inv, lc_inv, big_P);

      if (OK == 0){
         printf(" some problem\n");
         abort();
      }

      F_mpz_poly_scalar_mul(f1, f, lc_inv);

      F_mpz_print(lc_inv); printf(" is it big\n");

      F_mpz_poly_smod(f1, f1, big_P);
   }

   for (i = 0; i < 2*r-2; i++)
      F_mpz_poly_init(v[i]);

   for (i = 0; i < 2*r-2; i++)
      F_mpz_poly_init(w[i]);

   _Build_Hensel_Tree(link, v, w, fac);

//clearing fac early, don't need them anymore

   zmod_poly_factor_clear(fac);

   ulong e = 1;

   while( (2 * e) <= a ){
      _Tree_Hensel_Lift(link, v, w, e, 2*e, f1, 1, p, r, P);
      e = 2*e;
   }

   F_mpz_poly_print(f1); printf("\n was the special monic f equiv mod p^e\n");
   F_mpz_print(big_P); printf("\n was big_P\n");

   F_mpz_poly_clear( f1 );
   F_mpz_clear( lc_inv );
//Have now Hensel lifted to p^a for the precalculated a, in the optimized version we will lift even less

//Here let's make a list of Hensel lifted factors for grabbing information and trial testing.
   F_mpz_poly_factor_t lifted_fac;
   F_mpz_poly_factor_init(lifted_fac);

//If r>5 then we'll need to allocate enough room for r factors

   if (r > lifted_fac->alloc)
   {
      lifted_fac->factors = (F_mpz_poly_t *) flint_heap_realloc_bytes(lifted_fac->factors, sizeof(F_mpz_poly_t)*r);
      lifted_fac->exponents = (unsigned long *) flint_heap_realloc(lifted_fac->exponents, r);
      for (unsigned long i = lifted_fac->alloc; i < r; i++)
         F_mpz_poly_init(lifted_fac->factors[i]);
      lifted_fac->alloc = r;
   }


//Now we should undo the mystical link part to find the original local factors in original order, see the long explanation in the Hensel code
   for(i = 0; i < 2*r -2; i++){
      if (link[i] < 0){
         F_mpz_poly_smod(v[i], v[i], P); 
         F_mpz_poly_set(lifted_fac->factors[-link[i]-1],v[i]);
         lifted_fac->exponents[-link[i]-1] = 1L; 
      }
   }
   lifted_fac->num_factors = r;

//Now we are ready to to the Zassenhaus testing... later the r > 20 (or even 10) test could go here

   F_mpz_poly_factor_print(lifted_fac);

   F_mpz_poly_zassenhaus_naive(final_fac, lifted_fac, f, P, exp, lc);

//Done factoring, just gotta clean house

   for (i = 0; i < 2*r-2; i++){
      F_mpz_poly_clear(v[i]);
      F_mpz_poly_clear(w[i]);
   }

   F_mpz_poly_factor_clear(lifted_fac);
   zmod_poly_clear(F);
   F_mpz_clear(lc);
   F_mpz_clear(big_P);
   F_mpz_clear(P);
   return;
}


void F_mpz_poly_factor(F_mpz_poly_factor_t final_fac, F_mpz_t cong, F_mpz_poly_t G){

//Need to store the square free factors

   if (G->length == 0){
      F_mpz_set_ui(cong, 0UL);
      return;
   }

   if (G->length == 1){
      F_mpz_set(cong, G->coeffs);
      return;
   }

   F_mpz_poly_t g;
   F_mpz_poly_init(g);

   if (G->length == 2){
      F_mpz_poly_content(cong, G);
      F_mpz_poly_scalar_div_exact(g, G, cong);
      F_mpz_poly_factor_insert( final_fac, g, 1UL);
      return;
   }

   ulong x_pow = 0;
   while (F_mpz_is_zero(G->coeffs + x_pow)){
      x_pow++; 
   }


   if (x_pow != 0){
      F_mpz_poly_t temp_x;
      F_mpz_poly_init(temp_x);
      F_mpz_poly_set_coeff_ui(temp_x, 1, 1);
      F_mpz_poly_factor_insert(final_fac, temp_x, x_pow);
      F_mpz_poly_clear(temp_x);
   }

   F_mpz_poly_right_shift(g, G, x_pow);

   F_mpz_poly_factor_t sq_fr_fac;

   F_mpz_poly_factor_init( sq_fr_fac );

   F_mpz_poly_squarefree(sq_fr_fac, cong, g);

//Now we can go through and factor each square free guy and add it to final factors.

   for(ulong j = 0; j < sq_fr_fac->num_factors; j++){
      F_mpz_poly_factor_sq_fr_prim(final_fac, sq_fr_fac->exponents[j], sq_fr_fac->factors[j]);
   }


   F_mpz_poly_factor_clear( sq_fr_fac );
   F_mpz_poly_clear(g);

}



int test_F_mpz_poly_factor()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   F_mpz_poly_factor_t F_factors;
   F_mpz_t content;
   int result = 1;
   long num_facs;
   ulong bits1, bits2, length1, length2;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   for (ulong count1 = 0; (count1 < 200*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);
      F_mpz_poly_factor_init(F_factors);
      F_mpz_init(content);

		bits1 = z_randint(100) + 1;
      bits2 = z_randint(100) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
		mpz_randpoly(m_poly2, length2, bits2);

		mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
//      printf(" floating not here\n");

//      F_mpz_poly_print(F_poly1); printf(" poly1\n");
//      F_mpz_poly_print(F_poly2); printf(" poly2\n");

		F_mpz_poly_mul(res, F_poly1, F_poly2);			

//      F_mpz_poly_print(res); printf(" res\n");

      F_mpz_poly_factor(F_factors, content, res);

		F_mpz_poly_to_mpz_poly(res2, res);

//      printf(" floating not here\n");

      F_mpz_poly_clear(res);
      F_mpz_poly_init(res);

      F_mpz_poly_set_coeff_ui(res, 0, 1);

      for (num_facs = 0; num_facs < F_factors->num_factors; num_facs++){
         for (ulong pow = 0; pow < F_factors->exponents[num_facs]; pow++){
            F_mpz_poly_mul(res, res, F_factors->factors[num_facs]);
         }
      }

      F_mpz_poly_scalar_mul( res, res, content);
		F_mpz_poly_to_mpz_poly(res1, res);
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error2: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         mpz_poly_print_pretty(m_poly1, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
         F_mpz_poly_factor_print(F_factors);
         F_mpz_print(content); printf(" content \n");
         mpz_poly_print_pretty(res1, "x"); printf("\n");
         mpz_poly_print_pretty(res2, "x"); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(res);
      F_mpz_poly_factor_clear(F_factors);
      F_mpz_clear(content);
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
	
   RUN_TEST(F_mpz_poly_factor);
	
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
