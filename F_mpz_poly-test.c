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

double random_d()
{
   if (z_randint(2)) return rand()/((double) RAND_MAX + 1);
   else return -rand()/((double) RAND_MAX + 1);
}

void F_mpz_test_random(F_mpz_t f, ulong bits)
{
	if (bits == 0)
	{
		F_mpz_zero(f);
      return;
	}
	
	mpz_t temp;
	mpz_init(temp);
	
	mpz_rrandomb(temp, randstate, bits);
#if SIGNS
	if (z_randint(2)) mpz_neg(temp, temp);
#endif
   
	F_mpz_set_mpz(f, temp);

   mpz_clear(temp);
}

// generate a random mpz_poly_t with up to the given length and number of bits per coefficient
void mpz_randpoly(mpz_poly_t pol, long length, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_ensure_alloc(pol, length);
	mpz_poly_zero(pol);
   
   long i;
   for (i = 0; i < length; i++)
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
   
   long i;
   for (i = 0; i < length; i++)
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
   
   long i;
   for (i = 0; i < length; i++)
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

int test_F_mpz_poly_to_mpz_poly()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   ulong count1;
   for (count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
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

int test_F_mpz_poly_to_fmpz_poly()
{
   mpz_poly_t m_poly;
   fmpz_poly_t f_poly;
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits, length;
   
   mpz_poly_init(m_poly); 

   ulong count1;
   for (count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      fmpz_poly_init(f_poly);
      
      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly);
      F_mpz_poly_to_fmpz_poly(f_poly, F_poly1);
      fmpz_poly_to_F_mpz_poly(F_poly2, f_poly);
          
      result = F_mpz_poly_equal(F_poly1, F_poly2); 
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, F_poly1->length, F_poly2->length);
         F_mpz_poly_print_pretty(F_poly1, "x"); printf("\n");
         F_mpz_poly_print_pretty(F_poly2, "x"); printf("\n");
		}
          
      fmpz_poly_clear(f_poly);
      
      F_mpz_poly_clear(F_poly1);
      F_mpz_poly_clear(F_poly2);
   }
   
   mpz_poly_clear(m_poly);
  
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

   ulong count1;
   for (count1 = 0; (count1 < 20000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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

   ulong count1;
   for (count1 = 0; (count1 < 20000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      // set random coeffs in the poly
		ulong count2;
		for (count2 = 0; (count2 < 500) && result == 1; count2++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      // set random coeffs in the poly
		ulong count2;
		for (count2 = 0; (count2 < 500) && result == 1; count2++)
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
					  printf("Error: length = %ld, coeff_num = %ld, coeff = %ld, coeff2 = %ld\n", length, coeff_num, coeff, coeff2);
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
   
   ulong count1;
   for (count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      
      F_mpz_poly_init(F_poly);

      length = z_randint(100)+1;        
      
		F_mpz_randpoly(F_poly, length, bits); 
		    
      // set random coeffs in the poly
		ulong count2;
		for (count2 = 0; (count2 < 300) && result == 1; count2++)
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

   ulong count1;
   for (count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
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
	ulong count1;
	for (count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
	for (count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
	for (count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
	unsigned long count1;
	for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(200) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_bits = 0;
      sign = 1L;
      ulong i;
      for (i = 0; i < m_poly1->length; i++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(FLINT_BITS-2) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_bits = 0;
      sign = 1L;
      ulong i;
      for (i = 0; i < m_poly1->length; i++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 100000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly);

      bits = z_randint(500) + 1;
      length = z_randint(100);
      mpz_randpoly(m_poly1, length, bits);
           
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      
		mpz_limbs = 0;
      ulong i;
      for (i = 0; i < m_poly1->length; i++)
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

int test_F_mpz_poly_tofromstring()
{
   mpz_poly_t test_poly;
   F_mpz_poly_t test_F_mpz_poly, test_F_mpz_poly2;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   unsigned long count1;
   for (count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = z_randint(100)+ 1;
      
      F_mpz_poly_init2(test_F_mpz_poly, (bits-1)/FLINT_BITS+1);
      F_mpz_poly_init2(test_F_mpz_poly2, (bits-1)/FLINT_BITS+1+z_randint(30));
      FILE * testfile; 
      
      unsigned long count2;
      for (count2 = 0; (count2 < 10) && (result == 1); count2++)
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

int test_F_mpz_poly_neg()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, F_poly4;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // negate and negate back
	ulong count1;
	for (count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
	for (count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
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
	unsigned long count1;
	for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
	unsigned long count1, count2;
	for (count1 = 0; (count1 < 700) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 700) && (result == 1) ; count1++)
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
	for (count2 = 0; (count2 < 1000) && (result == 1); count2++)
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
   unsigned long count1;
   for (count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = z_randint(FLINT_BITS-3) + 1;
      
      F_mpz_poly_init(test_F_mpz_poly);
      F_mpz_poly_init(test_F_mpz_poly2);

      unsigned long count2;
      for (count2 = 0; (count2 < 100) && (result == 1); count2++)
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
          unsigned long c;
          for (c = 0; c < test_F_mpz_poly->length; c++){
             F_mpz_sub(temp, test_F_mpz_poly->coeffs + c, test_F_mpz_poly2->coeffs + c);
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
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
          
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mult == 0L) result = (F_poly2->length == 0);
		else
		{
			ulong i;
			for (i = 0; i < m_poly->length; i++)
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
	for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
			ulong i;
			for (i = 0; i < m_poly->length; i++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
          
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mult == 0L) result = (F_poly2->length == 0);
		else
		{
			ulong i;
			for (i = 0; i < m_poly->length; i++)
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
	for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
			ulong i;
			for (i = 0; i < m_poly->length; i++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
          
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      if (mpz_sgn(mult) == 0) result = (F_poly2->length == 0);
		else
		{
			ulong i;
			for (i = 0; i < m_poly->length; i++)
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
	for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
			ulong i;
			for (i = 0; i < m_poly->length; i++)
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

int test_F_mpz_poly_scalar_divexact()
{
   F_mpz_poly_t F_poly, F_poly2, F_poly3;
	F_mpz_t x;
   int result = 1;
   ulong bits, bits2, length;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
		F_mpz_init(x);
      
      F_mpz_randpoly(F_poly, length, bits); 
          
      F_mpz_test_random(x, 200);
     
		F_mpz_poly_scalar_mul(F_poly2, F_poly, x);
		F_mpz_poly_scalar_divexact(F_poly3, F_poly2, x);
          
      result = (F_mpz_poly_equal(F_poly, F_poly3));

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         F_mpz_print(x); printf("\n");
         F_mpz_poly_print_pretty(F_poly, "x"); printf("\n");
         F_mpz_poly_print_pretty(F_poly3, "x"); printf("\n");
		}

      F_mpz_clear(x);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased 
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
		F_mpz_init(x);
      
      F_mpz_randpoly(F_poly, length, bits); 
          
      F_mpz_test_random(x, 200);
     
		F_mpz_poly_scalar_mul(F_poly2, F_poly, x);
		F_mpz_poly_scalar_divexact(F_poly2, F_poly2, x);
          
      result = (F_mpz_poly_equal(F_poly, F_poly2));

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         F_mpz_print(x); printf("\n");
         F_mpz_poly_print_pretty(F_poly, "x"); printf("\n");
         F_mpz_poly_print_pretty(F_poly2, "x"); printf("\n");
		}

      F_mpz_clear(x);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
      
   return result; 
}

int test_F_mpz_poly_scalar_smod()
{
   F_mpz_poly_t F_poly, F_poly2, F_poly3;
	F_mpz_t x, y;
   int result = 1;
   ulong bits, bits2, length, i;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ 1;
      bits2 = z_randint(300)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      F_mpz_init(x);
      F_mpz_init(y);
      
      F_mpz_randpoly(F_poly, length, bits); 
          
      do { F_mpz_test_random(x, bits2); }
      while (F_mpz_is_zero(x));

		F_mpz_poly_scalar_smod(F_poly2, F_poly, x);
		
      for (i = 0; i < F_poly2->length; i++)
      {
         F_mpz_smod(y, F_poly->coeffs + i, x);
         result &= (F_mpz_equal(y, F_poly2->coeffs + i));
      }

      for ( ; i < F_poly->length; i++)
      {
         F_mpz_smod(y, F_poly->coeffs + i, x);
         result &= F_mpz_is_zero(y);
      }

      if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         F_mpz_print(x); printf("\n");
         F_mpz_poly_print_pretty(F_poly, "x"); printf("\n");
         F_mpz_poly_print_pretty(F_poly2, "x"); printf("\n");
		}

      F_mpz_clear(x);
		F_mpz_clear(y);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased 
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(300)+ 1;
      bits2 = z_randint(300)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      F_mpz_init(x);
      F_mpz_init(y);
      
      F_mpz_randpoly(F_poly, length, bits); 
          
      do { F_mpz_test_random(x, bits2); }
      while (F_mpz_is_zero(x));

		F_mpz_poly_scalar_smod(F_poly2, F_poly, x);
		F_mpz_poly_scalar_smod(F_poly, F_poly, x);
		
      result = (F_mpz_poly_equal(F_poly, F_poly2));

      if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         F_mpz_print(x); printf("\n");
         F_mpz_poly_print_pretty(F_poly, "x"); printf("\n");
         F_mpz_poly_print_pretty(F_poly2, "x"); printf("\n");
		}

      F_mpz_clear(x);
		F_mpz_clear(y);
		F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
      
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

   ulong count1;
   for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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

int test_F_mpz_poly_mul_classical_trunc_left()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2, trunc;
   long len_out;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   ulong count1;
   for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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

      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_classical_trunc_left(res, F_poly1, F_poly2, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);
      long i;
      for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_classical_trunc_left(res, res, F_poly1, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_classical_trunc_left(res, F_poly1, res, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      len_out = 2*m_poly1->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_classical_trunc_left(res, F_poly1, F_poly1, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly1);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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

   ulong count1;
   for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 2500*ITER) && (result == 1) ; count1++)
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

int test_F_mpz_poly_mul_karatsuba_trunc_left()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2, trunc;
   long len_out;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   ulong count1;
   for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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

      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_karatsuba_trunc_left(res, F_poly1, F_poly2, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);
      long i;
      for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_karatsuba_trunc_left(res, res, F_poly1, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(100);
      length2 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_karatsuba_trunc_left(res, F_poly1, res, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly2);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      mpz_randpoly(m_poly1, length1, bits1);
           
      len_out = 2*m_poly1->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_karatsuba_trunc_left(res, F_poly1, F_poly1, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul_classical(res1, m_poly1, m_poly1);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
         mpz_set_ui(res1->coeffs[i], 0);
      mpz_poly_normalise(res1);
		    
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

int test_F_mpz_poly_mul()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2, trunc;
   long len_out;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   ulong count1;
   for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
     
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul(res, F_poly1, F_poly2);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly2);
		    
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
	for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(500);
      length2 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul(res, res, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly2);		
		    
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
	for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(500);
      length2 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul(res, F_poly1, res);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly2);		
		    
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
	for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      length1 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul(res, F_poly1, F_poly1);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly1);		
		    
      result = mpz_poly_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
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

int test_F_mpz_poly_mul_trunc_left()
{
   mpz_poly_t m_poly1, m_poly2, res1, res2;
   F_mpz_poly_t F_poly1, F_poly2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2, trunc;
   long len_out;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 
   mpz_poly_init(res1); 
   mpz_poly_init(res2); 

   ulong count1;
   for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(500);
		length2 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);

      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
		F_mpz_poly_mul_trunc_left(res, F_poly1, F_poly2, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly2);
      long i;
      for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(500);
      length2 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_trunc_left(res, res, F_poly1, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly2);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      length1 = z_randint(500);
      length2 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
      mpz_randpoly(m_poly2, length2, bits2);
           
      len_out = m_poly1->length + m_poly2->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      mpz_poly_to_F_mpz_poly(res, m_poly2);
      
		F_mpz_poly_mul_trunc_left(res, F_poly1, res, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly2);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
		    
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
	for (count1 = 0; (count1 < 500*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(res);

		bits1 = z_randint(200) + 1;
      length1 = z_randint(500);
      mpz_randpoly(m_poly1, length1, bits1);
           
      len_out = 2*m_poly1->length - 1;
      trunc = z_randint(len_out + 1);
           
      mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
      
		F_mpz_poly_mul_trunc_left(res, F_poly1, F_poly1, trunc);
		F_mpz_poly_to_mpz_poly(res2, res);
      mpz_poly_mul(res1, m_poly1, m_poly1);		
		long i;
		for (i = 0; (i < len_out) && (i < trunc) && (i < res1->length); i++)
      {
         mpz_set_ui(res1->coeffs[i], 0);
         mpz_set_ui(res2->coeffs[i], 0);
      }
      mpz_poly_normalise(res1);
		    
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

   ulong count1;
   for (count1 = 0; (count1 < 500) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 500) && (result == 1) ; count1++)
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

   ulong count1;
   for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(res);
      F_mpz_poly_factor_init(F_factors);
      F_mpz_init(content);

	  do {
	    bits1 = z_randint(50) + 1;
        bits2 = z_randint(50) + 1;
        length1 = z_randint(50);
        length2 = z_randint(50);
        mpz_randpoly(m_poly1, length1, bits1);
		mpz_randpoly(m_poly2, length2, bits2);

		mpz_poly_to_F_mpz_poly(F_poly1, m_poly1);
        mpz_poly_to_F_mpz_poly(F_poly2, m_poly2);
      
//      printf(" floating not here\n");

//      F_mpz_poly_print(F_poly1); printf(" poly1\n");
//      F_mpz_poly_print(F_poly2); printf(" poly2\n");

		F_mpz_poly_mul(res, F_poly1, F_poly2);			
	  } while (res->length == 0);

      F_mpz_poly_print(res); printf(" res\n");

      F_mpz_poly_factor(F_factors, content, res);

		F_mpz_poly_to_mpz_poly(res2, res);

//      printf(" floating not here\n");

      F_mpz_poly_clear(res);
      F_mpz_poly_init(res);

      F_mpz_poly_set_coeff_ui(res, 0, 1);

      for (num_facs = 0; num_facs < F_factors->num_factors; num_facs++){
         ulong pow;
         for (pow = 0; pow < F_factors->exponents[num_facs]; pow++){
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

int test_F_mpz_poly_bit_pack_unsigned()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   mp_limb_t * array;
	int result = 1;
   ulong bits, length, depth;
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 

   ulong count1;
   for (count1 = 0; (count1 < 500) && (result == 1) ; count1++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 500) && (result == 1) ; count1++)
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

		ulong i;
		for (i = 0; i < m_poly3->length; i++)
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

   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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

   ulong count1;
   for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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

   ulong count1;
   for (count1 = 1; (count1 < 200*ITER) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 2;
      bytes = ((bits - 1)>>3) + 1;
      
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      ulong count2;
      for (count2 = 0; (count2 < 10) && (result == 1); count2++)
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
          ulong j;
          for (j = 0; j < m_poly1->length; j++)
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
          ulong j;
          for (j = 0; j < m_poly2->length; j++)
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

   ulong count1;
   for (count1 = 1; (count1 < 200*ITER) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 2;
      bytes = ((bits - 1)>>3) + 1;
      
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);

      ulong count2;
      for (count2 = 0; (count2 < 10) && (result == 1); count2++)
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
          ulong j;
          for (j = 0; j < m_poly1->length; j++)
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
          ulong j;
          for (j = 0; j < m_poly2->length; j++)
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

   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1) ; count1++)
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

   ulong count1;
   for (count1 = 1; (count1 < 1000*ITER) && (result == 1) ; count1++)
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

int test_F_mpz_poly_divrem_basecase()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_basecase(Q, R, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, F_poly2) && (R->length == 0)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_divrem_basecase(Q, R, F_poly2, F_poly1);
      F_mpz_poly_mul(F_poly3, Q, F_poly1);
      F_mpz_poly_add(F_poly3, F_poly3, R);

      result = (F_mpz_poly_equal(F_poly3, F_poly2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(F_poly3); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_div_basecase()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_basecase(Q, R, F_poly3, F_poly1);
      F_mpz_poly_div_basecase(Q2, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_divrem_basecase(Q, R, F_poly2, F_poly1);
      F_mpz_poly_div_basecase(Q2, F_poly2, F_poly1);
      
      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_div_divconquer_recursive()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, Qb, Qb2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(Qb);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_div_divconquer_recursive(Q, Qb, F_poly3, F_poly1);
      F_mpz_poly_div_basecase(Q2, F_poly3, F_poly1);
      
      result = (F_mpz_poly_equal(Q, Q2) && F_mpz_poly_equal(Qb, F_poly3)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(F_poly3); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
         F_mpz_poly_print(Qb); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(Qb);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Qb);
      F_mpz_poly_init(Qb2);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_div_divconquer_recursive(Q, Qb, F_poly2, F_poly1);
      F_mpz_poly_mul(Qb2, Q, F_poly1);

      result = (F_mpz_poly_equal(Qb, Qb2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Qb); printf("\n");
         F_mpz_poly_print(Qb2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Qb);
		F_mpz_poly_clear(Qb2);
   }

	return result;
}

int test_F_mpz_poly_divrem_divconquer()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_divconquer(Q, R, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, F_poly2) && (R->length == 0)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_divrem_divconquer(Q, R, F_poly2, F_poly1);
      F_mpz_poly_mul(F_poly3, Q, F_poly1);
      F_mpz_poly_add(F_poly3, F_poly3, R);

      result = (F_mpz_poly_equal(F_poly3, F_poly2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(F_poly3); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_divrem_basecase_low()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_basecase_low(Q, R, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, F_poly2) && (R->length == 0)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_divrem_basecase_low(Q, R, F_poly2, F_poly1);
      F_mpz_poly_mul(F_poly3, Q, F_poly1);
      F_mpz_poly_add(F_poly3, F_poly3, R);
      F_mpz_poly_truncate(F_poly3, F_poly1->length - 1);
      F_mpz_poly_truncate(F_poly2, F_poly1->length - 1);

      result = ((F_mpz_poly_equal(F_poly3, F_poly2)) && (R->length <= F_poly1->length - 1)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(F_poly3); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_div_divconquer_recursive_low()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, Qb, Qb2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(Qb);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_div_divconquer_recursive_low(Q, Qb, F_poly3, F_poly1);
      F_mpz_poly_div_basecase(Q2, F_poly3, F_poly1);
      F_mpz_poly_truncate(F_poly3, F_poly1->length - 1);
      
      result = (F_mpz_poly_equal(Q, Q2) && F_mpz_poly_equal(Qb, F_poly3)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(F_poly3); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
         F_mpz_poly_print(Qb); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(Qb);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Qb);
      F_mpz_poly_init(Qb2);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_div_divconquer_recursive_low(Q, Qb, F_poly2, F_poly1);
      F_mpz_poly_mul(Qb2, Q, F_poly1);
      F_mpz_poly_truncate(Qb2, F_poly1->length - 1);

      result = (F_mpz_poly_equal(Qb, Qb2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Qb); printf("\n");
         F_mpz_poly_print(Qb2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Qb);
		F_mpz_poly_clear(Qb2);
   }

	return result;
}

int test_F_mpz_poly_div_divconquer()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_basecase(Q, R, F_poly3, F_poly1);
      F_mpz_poly_div_divconquer(Q2, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }
   
   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_divrem_basecase(Q, R, F_poly2, F_poly1);
      F_mpz_poly_div_divconquer(Q2, F_poly2, F_poly1);
      
      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_div_hensel()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_basecase(Q, R, F_poly3, F_poly1);
      F_mpz_poly_div_hensel(Q2, F_poly3, F_poly3->length, F_poly1, F_poly1->length);

      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }

	return result;
}

int test_F_mpz_poly_divexact()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_divrem_basecase(Q, R, F_poly3, F_poly1);
      F_mpz_poly_divexact(Q2, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }

	return result;
}

int test_F_mpz_poly_pseudo_divrem_basecase()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, R;
   int result = 1;
   ulong bits1, bits2, length1, length2, d;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_pseudo_divrem_basecase(Q, R, &d, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, F_poly2) && (R->length == 0)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }

   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_pseudo_divrem_basecase(Q, R, &d, F_poly2, F_poly1);
      F_mpz_poly_mul(F_poly3, Q, F_poly1);
      F_mpz_poly_add(F_poly3, F_poly3, R);
      
      F_mpz * B_lead = F_poly1->coeffs + F_poly1->length - 1;
      F_mpz_t pow;
      F_mpz_init(pow);
      F_mpz_pow_ui(pow, B_lead, d);
      F_mpz_poly_scalar_mul(F_poly2, F_poly2, pow);
      F_mpz_clear(pow);

      result = (F_mpz_poly_equal(F_poly3, F_poly2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(F_poly3); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(R); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_pseudo_div_basecase()
{
   F_mpz_poly_t F_poly1, F_poly2, F_poly3, Q, Q2, R;
   int result = 1;
   ulong bits1, bits2, length1, length2, d;
   
   // test exact division
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(F_poly3);
      F_mpz_poly_init(Q);
      
		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_mul(F_poly3, F_poly1, F_poly2);			
		F_mpz_poly_pseudo_div_basecase(Q, &d, F_poly3, F_poly1);

      result = (F_mpz_poly_equal(Q, F_poly2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(F_poly3);
		F_mpz_poly_clear(Q);
   }
   
   // test inexact division
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      F_mpz_poly_init(Q);
      F_mpz_poly_init(Q2);
      F_mpz_poly_init(R);

		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200);
      length1 = z_randint(100) + 1;
      length2 = z_randint(100);
      
      do F_mpz_randpoly(F_poly1, length1, bits1);
      while (F_poly1->length == 0);
      F_mpz_randpoly(F_poly2, length2, bits2);
           
		F_mpz_poly_pseudo_div_basecase(Q, &d, F_poly2, F_poly1);
      F_mpz_poly_pseudo_divrem_basecase(Q2, R, &d, F_poly2, F_poly1);

      result = (F_mpz_poly_equal(Q, Q2)); 
		
      if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
         F_mpz_poly_print(F_poly1); printf("\n");
         F_mpz_poly_print(F_poly2); printf("\n");
         F_mpz_poly_print(Q); printf("\n");
         F_mpz_poly_print(Q2); printf("\n");
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
		F_mpz_poly_clear(Q);
		F_mpz_poly_clear(Q2);
		F_mpz_poly_clear(R);
   }
      
	return result;
}

int test_F_mpz_poly_derivative()
{
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   F_mpz_t t;
   F_mpz_init(t);
            
   // check coeffs of derivative
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      F_mpz_poly_init(F_poly2);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100);
      F_mpz_randpoly(F_poly1, length1, bits1);
      
		F_mpz_poly_derivative(F_poly2, F_poly1);
      
      if (F_poly1->length <= 1)
         result = (F_poly2->length == 0);
      else
      {
         long j, coeff;

         for (j = 0; j < 100; j++)
         {
            coeff = z_randint(F_poly1->length - 1) + 1;
            F_mpz_mul_ui(t, F_poly1->coeffs + coeff, coeff);
            result &= F_mpz_equal(t, F_poly2->coeffs + coeff - 1);
         }
      } 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
		F_mpz_poly_clear(F_poly2);
   }

   F_mpz_clear(t);

   return result;
}

int test_F_mpz_poly_content()
{
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   F_mpz_t c1, c2;
   F_mpz_init(c1);
   F_mpz_init(c2);
            
   // check giving a content to a content free polynomial yields the correct content
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100) + 1;
      
      do {
         F_mpz_randpoly(F_poly1, length1, bits1);
         F_mpz_poly_content(c1, F_poly1);
      } while (!F_mpz_is_one(c1));

      F_mpz_test_random(c1, 100);
      F_mpz_poly_scalar_mul(F_poly1, F_poly1, c1);

      F_mpz_poly_content(c2, F_poly1);

		result = F_mpz_equal(c1, c2);

		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_poly_clear(F_poly1);
   }

   F_mpz_clear(c1);
   F_mpz_clear(c2);

   return result;
}

int test_F_mpz_poly_eval_horner_d()
{
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   double d1, d2, val;
   
   // check giving a content to a content free polynomial yields the correct content
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100) + 1;
      
      do {
         F_mpz_randpoly(F_poly1, length1, bits1);
      } while (F_poly1->length < 1);

      val = random_d();

      F_poly1->coeffs++;
      F_poly1->length--;

      d1 = F_mpz_poly_eval_horner_d(F_poly1, val);
      
      F_poly1->coeffs--;
      F_poly1->length++;

      d1 *= val;
      d1 += F_mpz_get_d(F_poly1->coeffs);

      d2 = F_mpz_poly_eval_horner_d(F_poly1, val);
      
      result = (d1 == d2);

		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, d1 = %lf, d2 = %lf\n", 
                                                    length1, bits1, d1, d2);
		}
          
      F_mpz_poly_clear(F_poly1);
   }

   return result;
}

int test_F_mpz_poly_eval_horner_d_2exp()
{
   F_mpz_poly_t F_poly1, F_poly2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   double d1, d2, val;
   long exp1, exp2, exp;
   
   // check that negating a poly returns minus the evaluation
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100) + 1;
      
      F_mpz_randpoly(F_poly1, length1, bits1);
      
      val = random_d();

      d1 = F_mpz_poly_eval_horner_d_2exp(&exp1, F_poly1, val);
      
      F_mpz_poly_neg(F_poly1, F_poly1); 
      d2 = F_mpz_poly_eval_horner_d_2exp(&exp2, F_poly1, val);
      
      result = ((d1 == -d2) && (exp1 == exp2));

		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, d1 = %lf, d2 = %lf, "
            "exp1 = %ld, exp2 = %ld\n", length1, bits1, d1, d2, exp1, exp2);
		}
          
      F_mpz_poly_clear(F_poly1);
   }

   // check that multiplying by a power of 2 works
	for (ulong count1 = 0; (count1 < 20000*ITER) && (result == 1); count1++)
   {
      F_mpz_poly_init(F_poly1);
      
		bits1 = z_randint(200) + 1;
      length1 = z_randint(100) + 1;
      
      F_mpz_randpoly(F_poly1, length1, bits1);
      
      val = random_d();

      d1 = F_mpz_poly_eval_horner_d_2exp(&exp1, F_poly1, val);
      
      exp = z_randint(FLINT_BITS);
      F_mpz_poly_scalar_mul_ui(F_poly1, F_poly1, 1UL<<exp); 
      d2 = F_mpz_poly_eval_horner_d_2exp(&exp2, F_poly1, val);
      
      result = ((fabs(d1 - d2) < 0.000000000000001f) && ((d1 == 0) || (exp1 + exp == exp2)));

		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, d1 = %lf, d2 = %lf, "
            "exp1 = %ld, exp2 = %ld, exp = %ld\n", length1, bits1, d1, 
            d2, exp1, exp2, exp);
		}
          
      F_mpz_poly_clear(F_poly1);
   }

   return result;
}

int test_F_mpz_poly_scalar_abs()
{
   mpz_poly_t m_poly, m_poly2;
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, length, i;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(m_poly); 
   mpz_poly_init(m_poly2); 
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);  

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      F_mpz_poly_scalar_abs(F_poly2, F_poly);
          
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly2); 
          
      for (i = 0; i < m_poly->length; i++)
      {
         mpz_abs(temp, m_poly->coeffs[i]);
         result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
      }

		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
         mpz_poly_print_pretty(m_poly, "x"); printf("\n");
         mpz_poly_print_pretty(m_poly2, "x"); printf("\n");
		}

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
   
   // aliased abs
	for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      length = z_randint(200);        

      F_mpz_poly_init(F_poly);
      
		mpz_randpoly(m_poly, length, bits); 
      mpz_poly_to_F_mpz_poly(F_poly, m_poly);
          
      F_mpz_poly_scalar_abs(F_poly, F_poly);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly); 
          
      for (i = 0; i < m_poly->length; i++)
      {
         mpz_abs(temp, m_poly->coeffs[i]);
         result &= (mpz_cmp(temp, m_poly2->coeffs[i]) == 0);
      }
		
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld\n", length, bits);
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

int test_F_mpz_poly_CLD_bound()
{
   F_mpz_poly_t F_poly, F_poly2;
   int result = 1;
   ulong bits, length, i;
   F_mpz_t sum, bound;
   
   /* We check that CLD_bound is between the absolutely value of the n-th 
   coeff of f' and the sum of the absolute values of the coeffs of f' */
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(20)+ 1;
      length = z_randint(20);  
      
      printf("bits = %ld, length = %ld\n", bits, length);

      F_mpz_poly_init(F_poly);
      F_mpz_poly_init(F_poly2);
      
      F_mpz_init(sum);
      F_mpz_init(bound);
      
      F_mpz_randpoly(F_poly, length, bits); 
          
      F_mpz_poly_derivative(F_poly2, F_poly);
      F_mpz_poly_scalar_abs(F_poly2, F_poly2);
               
      for (i = 0; i < F_poly2->length; i++)
         F_mpz_add(sum, sum, F_poly2->coeffs + i);
      printf("bits = %ld, length = %ld\n", bits, length);

      for (i = 0; i < F_poly2->length && result == 1; i++)
      {
         printf("i = %ld, len = %ld\n", i, F_poly->length);
         F_mpz_poly_print(F_poly); printf("\n");
         F_mpz_poly_CLD_bound(bound, F_poly, i);
         printf("done = %ld\n", i);
         result &= (F_mpz_cmp(F_poly2->coeffs + i, bound) <= 0);
         result &= (F_mpz_cmp(sum, bound) >= 0);
      }
       printf("bits = %ld, length = %ld\n", bits, length);
  
		if (!result) 
		{
			printf("Error: length = %ld, bits = %ld, i = %ld\n", length, bits, i);
         F_mpz_print(sum);
         F_mpz_print(F_poly2->coeffs + i);
         F_mpz_print(bound);   
		}

      F_mpz_clear(sum);
      F_mpz_clear(bound);

      F_mpz_poly_clear(F_poly2);
      F_mpz_poly_clear(F_poly);
   }
       
   return result; 
}

void F_mpz_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %d\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(F_mpz_poly_derivative); 
   RUN_TEST(F_mpz_poly_content); 
   RUN_TEST(F_mpz_poly_eval_horner_d); 
   RUN_TEST(F_mpz_poly_eval_horner_d_2exp); 
   RUN_TEST(F_mpz_poly_scalar_abs); 
   RUN_TEST(F_mpz_poly_to_mpz_poly); 
   RUN_TEST(F_mpz_poly_to_fmpz_poly); 
   RUN_TEST(F_mpz_poly_CLD_bound); 
   RUN_TEST(F_mpz_poly_getset_coeff_si); 
   RUN_TEST(F_mpz_poly_getset_coeff_ui); 
   RUN_TEST(F_mpz_poly_getset_coeff_mpz); 
   RUN_TEST(F_mpz_poly_set); 
   RUN_TEST(F_mpz_poly_equal); 
   RUN_TEST(F_mpz_poly_swap); 
   RUN_TEST(F_mpz_poly_max_bits1);
   RUN_TEST(F_mpz_poly_max_bits);
   RUN_TEST(F_mpz_poly_max_limbs);
   RUN_TEST(F_mpz_poly_tofromstring);
   RUN_TEST(F_mpz_poly_neg);
   RUN_TEST(F_mpz_poly_reverse); 
   RUN_TEST(F_mpz_poly_add); 
   RUN_TEST(F_mpz_poly_sub); 
   RUN_TEST(F_mpz_poly_shift);
   RUN_TEST(F_mpz_poly_to_zmod_poly);
   RUN_TEST(F_mpz_poly_scalar_mul_ui); 
   RUN_TEST(F_mpz_poly_scalar_mul_si); 
   RUN_TEST(F_mpz_poly_scalar_mul);
   RUN_TEST(F_mpz_poly_scalar_divexact);
   RUN_TEST(F_mpz_poly_scalar_smod);
   RUN_TEST(F_mpz_poly_mul_classical); 
   RUN_TEST(F_mpz_poly_mul_classical_trunc_left); 
   RUN_TEST(F_mpz_poly_mul_karatsuba); 
   RUN_TEST(F_mpz_poly_mul_karatsuba_trunc_left); 
   RUN_TEST(F_mpz_poly_bit_pack);
   RUN_TEST(F_mpz_poly_bit_pack_unsigned);
   RUN_TEST(F_mpz_poly_bit_pack2);
   RUN_TEST(F_mpz_poly_byte_pack_unsigned); 
   RUN_TEST(F_mpz_poly_byte_pack); 
   RUN_TEST(F_mpz_poly_mul_KS); 
   RUN_TEST(F_mpz_poly_mul_KS2);
   RUN_TEST(F_mpz_poly_mul_SS); 
   RUN_TEST(F_mpz_poly_mul); 
   RUN_TEST(F_mpz_poly_mul_trunc_left); 
   RUN_TEST(F_mpz_poly_pack_bytes); 
	RUN_TEST(F_mpz_poly_divrem_basecase); 
	RUN_TEST(F_mpz_poly_div_basecase); 
	RUN_TEST(F_mpz_poly_div_divconquer_recursive); 
	RUN_TEST(F_mpz_poly_divrem_divconquer); 
	RUN_TEST(F_mpz_poly_divrem_basecase_low); 
	RUN_TEST(F_mpz_poly_div_divconquer_recursive_low); 
	RUN_TEST(F_mpz_poly_div_divconquer); 
	RUN_TEST(F_mpz_poly_div_hensel); 
	RUN_TEST(F_mpz_poly_divexact); 
	RUN_TEST(F_mpz_poly_pseudo_divrem_basecase); 
	RUN_TEST(F_mpz_poly_pseudo_div_basecase); 
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
