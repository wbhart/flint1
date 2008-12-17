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

fmpz_poly-test.c: Test code for fmpz_poly.c and fmpz_poly.h

Copyright (C) 2007, David Harvey
Copyright (C) 2007, 2008, William Hart
Copyright (C) 2008, Richard Howell-Peak

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "test-support.h"
#include "zmod_poly.h"
#include "long_extras.h"
#include "mpz_poly.h"

#define VARY_BITS 1
#define SPARSE 1

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

/* 
   Generate a random integer in the range [0, limit) 
   If limit == 0, return a random limb
*/

unsigned long randint(unsigned long limit) 
{
#if FLINT_BITS == 32
    static uint64_t randval = 4035456057U;
    randval = ((uint64_t)randval*(uint64_t)1025416097U+(uint64_t)286824430U)%(uint64_t)4294967311U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)randval%limit;
#else
    static unsigned long randval = 4035456057U;
    static unsigned long randval2 = 6748392731U;
    randval = ((unsigned long)randval*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    randval2 = ((unsigned long)randval2*(unsigned long)1647637699U+(unsigned long)286824428U)%(unsigned long)4294967357U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)(randval+(randval2<<32))%limit;
#endif
}

/*
   Generate a random integer with up to the given number of bits [0, FLINT_BITS]
*/

unsigned long randbits(unsigned long bits)
{
   return randint(l_shift(1L, bits));
}

/* Return a random prime of (upto) the given number of bits [2, FLINT_BITS] */

unsigned long randprime(unsigned long bits)
{
   unsigned long limit, rand;
   
   if (bits < 2)
   {
      printf("FLINT Exception: attempt to generate prime < 2!\n");
      abort();
   }
   
   if (bits == FLINT_BITS)
   {
      do
      {
         rand = randbits(bits);
      
#if FLINT_BITS == 32
      }  while (rand > 4294967290UL);
#else
      }  while (rand > 18446744073709551556UL);
#endif
      rand = z_nextprime(rand);

   } else
   {
      do
      {
         rand = randbits(bits);
         rand = z_nextprime(rand);
      } while ((rand >> bits) > 0L);
   }
   
   return rand;
}

/* Generate a random zmod polynomial with the modulus n of the given length with 
   normalised coefficients */

void randpoly(zmod_poly_t poly, long length, unsigned long n)
{
   if (length == 0) 
   {
      zmod_poly_fit_length(poly, 1);
      poly->length = 0;
      return;
   }
              
   zmod_poly_fit_length(poly, length);
   
   for (unsigned long i = 0; i < length; i++)
      poly->coeffs[i] = randint(n);
   poly->length = length;
      
   __zmod_poly_normalise(poly);
} 

int test___zmod_poly_normalise()
{
   zmod_poly_t poly;
   int result = 1;
   unsigned long bits, length, nz_coeff;
   
   for (unsigned long count1 = 0; (count1 < 100000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS - 1) + 2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      
      length = randint(100) + 1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      randpoly(poly, length, modulus); 
      
		nz_coeff = z_randint(length);
		poly->coeffs[nz_coeff] = 1L;

		for (ulong i = nz_coeff + 1; i < length; i++) poly->coeffs[i] = 0L;
      poly->length = length;

		if (length != nz_coeff + 1) result &= (!__zmod_poly_normalised(poly));
      __zmod_poly_normalise(poly);
      result &= (__zmod_poly_normalised(poly));

      result = (poly->length == nz_coeff + 1);
      if (!result)
		{
			printf("nz_coeff = %ld, poly->length = %ld\n", nz_coeff, poly->length);
		}

      zmod_poly_clear(poly);
   }
         
   return result; 
}
  
int test_zmod_poly_truncate()
{
   zmod_poly_t poly;
   int result = 1;
   unsigned long bits, length, trunc;
   
   for (unsigned long count1 = 0; (count1 < 100000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS - 1) + 2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      
      length = randint(100) + 1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      randpoly(poly, length, modulus); 
      
		trunc = z_randint(length);
		if (trunc) poly->coeffs[trunc - 1] = 1L;

		if (poly->length < trunc) poly->length = trunc;

		zmod_poly_truncate(poly, trunc);

      result = (poly->length == trunc);
      if (!result)
		{
			printf("trunc = %ld, poly->length = %ld\n", trunc, poly->length);
		}

      zmod_poly_clear(poly);
   }
         
   return result; 
}
  
int test_zmod_poly_reverse()
{
   zmod_poly_t poly, poly2;
   int result = 1;
   unsigned long bits, length, length2;
   
   for (unsigned long count1 = 0; (count1 < 100000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);   
      
      length = randint(100);
      length2 = length + randint(200);
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
#endif

      randpoly(poly, length, modulus); 
                
      zmod_poly_reverse(poly2, poly, length2);
      zmod_poly_reverse(poly2, poly2, length2);
           
      result = zmod_poly_equal(poly2, poly);
      
      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);   
      
      length = randint(100);
      length2 = length + randint(200);
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
#endif

      randpoly(poly, length, modulus); 
          
      zmod_poly_set(poly2, poly);
      zmod_poly_reverse(poly, poly, length2);
      zmod_poly_reverse(poly, poly, length2);
           
      result = zmod_poly_equal(poly2, poly);
      
      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
   }
      
   return result; 
}
  
int test_zmod_poly_to_from_string()
{
   zmod_poly_t poly, poly2;
   int result = 1;
   unsigned long bits, length;
	char * str;
   
   for (unsigned long count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);   
      
      length = randint(100);

#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      randpoly(poly, length, modulus); 
                
      str = zmod_poly_to_string(poly);
		zmod_poly_from_string(poly2, str);
           
      result = zmod_poly_equal(poly2, poly);
		if (!result)
		{
			zmod_poly_print(poly); printf("\n\n");
         zmod_poly_print(poly2); printf("\n\n");
		}
      
      free(str);
		zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
   }
         
   return result; 
}
  
int test_zmod_poly_fprint_fread()
{
   zmod_poly_t poly, poly2;
   int result = 1;
   unsigned long bits, length;
	
   for (unsigned long count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);   
      
      length = randint(100);

#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      randpoly(poly, length, modulus); 
                
      FILE * file = fopen("tmp", "w");
		zmod_poly_fprint(poly, file);
		fflush(file);
		fclose(file);
      file = fopen("tmp", "r");
		zmod_poly_fread(poly2, file);
		fclose(file);
           
      result = zmod_poly_equal(poly2, poly);
		if (!result)
		{
			zmod_poly_print(poly); printf("\n\n");
         zmod_poly_print(poly2); printf("\n\n");
		}
      
      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
   }
         
   return result; 
}
  
int test_zmod_poly_addsub()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_add(res, pol1, pol2);
         zmod_poly_sub(res, res, pol2);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res);  
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);

         zmod_poly_set(res,pol1);
         zmod_poly_add(pol1, pol1, pol2);
         zmod_poly_sub(pol1, pol1, pol2);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res);  
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);

         zmod_poly_set(res,pol2);
         zmod_poly_add(pol2, pol1, pol2);
         zmod_poly_sub(pol1, pol2, pol1);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res);  
   }   
   return result;
}

int test_zmod_poly_add_no_red()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_add_no_red(res, pol1, pol2);
			for (ulong i = 0; i < res->length; i++)
				res->coeffs[i] %= res->p;
         zmod_poly_sub(res, res, pol2);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res);  
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_set(res,pol1);
         zmod_poly_add_no_red(pol1, pol1, pol2);
			for (ulong i = 0; i < pol1->length; i++)
				pol1->coeffs[i] %= pol1->p;
         zmod_poly_sub(pol1, pol1, pol2);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res);  
   }   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_set(res,pol2);
         zmod_poly_add_no_red(pol2, pol1, pol2);
			for (ulong i = 0; i < pol2->length; i++)
				pol2->coeffs[i] %= pol2->p;
         zmod_poly_sub(pol2, pol2, pol1);
         
         result &= zmod_poly_equal(res, pol2);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res);  
   } 

   return result;
}

int test_zmod_poly_neg()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_sub(res1, pol1, pol2);
         zmod_poly_neg(res2, pol2);
         zmod_poly_add(res2, res2, pol1);
         
         result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2);  
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_sub(res1, pol1, pol2);
         zmod_poly_neg(pol2, pol2);
         zmod_poly_add(res2, pol2, pol1);
         
         result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2);  
   }   
   return result;
}

int test_zmod_poly_bits()
{
   int result = 1;
   zmod_poly_t pol1;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
			randpoly(pol1, length1, modulus);
         
         ulong res1 = zmod_poly_bits(pol1);
			ulong res2 = 0;
			
			mpz_t temp;
			mpz_init(temp);
			for (ulong i = 0; i < pol1->length; i++)
			{
				mpz_set_ui(temp, pol1->coeffs[i]);
				ulong bits_t = mpz_sizeinbase(temp, 2);
				if (bits_t > res2) res2 = bits_t;
			}
         mpz_clear(temp);

         result &= (res1 == res2);
         
         if (!result)
         {
            printf("res1 = %ld, res2 = %ld\n", res1, res2); printf("\n");
         }
      }
      
      zmod_poly_clear(pol1);
   }

	return result;
}

int test_zmod_poly_shift()
{
   int result = 1;
   zmod_poly_t pol1, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long shift = randint(100);
         
         randpoly(pol1, length1, modulus);
         
         zmod_poly_left_shift(res, pol1, shift);
         zmod_poly_right_shift(res, res, shift);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res); 
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long shift = randint(100);
         
         randpoly(pol1, length1, modulus);

         zmod_poly_set(res,pol1);         
         zmod_poly_left_shift(pol1, pol1, shift);
         zmod_poly_right_shift(pol1, pol1, shift);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res); 
   }   
   return result;
}

int test_zmod_poly_swap()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         unsigned long shift = randint(100);
         
         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_sub(res1, pol1, pol2);
         zmod_poly_swap(pol1, pol2);
         zmod_poly_sub(res2, pol2, pol1);
         
         result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2);
   }
   
   return result;
}

int test_zmod_poly_setequal()
{
   int result = 1;
   zmod_poly_t pol1, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
         randpoly(pol1, length1, modulus);
         
         zmod_poly_set(res, pol1);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res); 
   }

   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
         randpoly(pol1, length1, modulus);
         
         zmod_poly_set(pol1, pol1);

         zmod_poly_set(res, pol1);

         result &= zmod_poly_equal(pol1, pol1);
         
         result &= zmod_poly_equal(res, res);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res); 
   }
   
   return result;
}

int test__zmod_poly_attach()
{
   int result = 1;
   zmod_poly_t pol1, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
     
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
         randpoly(pol1, length1, modulus);
         
         _zmod_poly_attach(res, pol1);
         
         result &= zmod_poly_equal(res, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
   }
   
   return result;
}

int test__zmod_poly_attach_truncate()
{
   int result = 1;
   zmod_poly_t pol1, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
     
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100) + 1;
         
         do randpoly(pol1, length1, modulus);
			while (!pol1->length);
         
         ulong trunc = z_randint(pol1->length + 1);
			_zmod_poly_attach_truncate(res, pol1, trunc);
         
			zmod_poly_truncate(pol1, trunc);

         result &= zmod_poly_equal(res, pol1);
         
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
      }
      
      zmod_poly_clear(pol1);
   }
   
   return result;
}

int test__zmod_poly_attach_shift()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
     
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100) + 1;
         
         do randpoly(pol1, length1, modulus);
			while (!pol1->length);
         
         ulong shift = z_randint(pol1->length + 1);
			_zmod_poly_attach_shift(res, pol1, shift);
         
			zmod_poly_right_shift(pol2, pol1, shift);

         result &= zmod_poly_equal(res, pol2);
         
         if (!result)
         {
            zmod_poly_print(res); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
         }
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
   }
   
   return result;
}

int test_zmod_poly_getset_coeff()
{
   int result = 1;
   zmod_poly_t pol1, pol2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long num = randint(200);
         unsigned long coeff = randint(modulus);
         
         randpoly(pol1, length1, modulus);
         zmod_poly_set(pol2, pol1);
         zmod_poly_set_coeff_ui(pol1, num, coeff);
         
         result &= (coeff == zmod_poly_get_coeff_ui(pol1, num));
         
         if (num + 1 > length1) 
         {
            zmod_poly_set_coeff_ui(pol1, num, 0);
            result &= zmod_poly_equal(pol1, pol2);
         }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1); 
      zmod_poly_clear(pol2);
   }
   
   return result;
}

int test_zmod_poly_mul_classicalKS()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            for (unsigned long i = 0; i < 10; i++)
               zmod_poly_mul_KS(res2, pol1, pol2, 0);
            
            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            for (unsigned long i = 0; i < 10; i++)
               zmod_poly_mul(res2, pol1, pol2);
            
            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }

	   zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }

   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      
		for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            zmod_poly_mul(pol1, pol1, pol2);
            
            result &= zmod_poly_equal(res1, pol1);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
            }
#endif
      
      }
	      
		zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
   }


   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            zmod_poly_mul(pol2, pol1, pol2);
            
            result &= zmod_poly_equal(res1, pol2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
            }
#endif
      
      }

      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
  }
   
   return result;
}

int test_zmod_poly_sqr_classicalKS()
{
   int result = 1;
   zmod_poly_t pol1, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            
            zmod_poly_sqr_classical(res1, pol1);
            zmod_poly_mul_KS(res2, pol1, pol1, 0);
         
            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_sqr()
{
   int result = 1;
   zmod_poly_t pol1, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            
            zmod_poly_sqr(res1, pol1);
            zmod_poly_sqr_classical(res2, pol1);
         
            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            
            zmod_poly_sqr_classical(res2, pol1);
            zmod_poly_sqr(pol1, pol1);
            
            result &= zmod_poly_equal(pol1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_classical_trunc()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
         unsigned long trunc;
         if (length1 + length2 > 1) trunc = randint(2*(length1 + length2 - 1));
         else trunc = 0;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            zmod_poly_truncate(res1, trunc);
            zmod_poly_mul_classical_trunc(res2, pol1, pol2, trunc);
            
            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_KS_trunc()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
         unsigned long trunc;
         if (length1 + length2 > 1) trunc = randint(2*(length1 + length2 - 1));
         else trunc = 0;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_KS(res1, pol1, pol2, 0);
            zmod_poly_truncate(res1, trunc);
            zmod_poly_mul_KS_trunc(res2, pol1, pol2, 0, trunc);
            
            result &= zmod_poly_equal(res1, res2);
            
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_trunc_n()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
         unsigned long trunc;
         if (length1 + length2 > 1) trunc = randint(2*(length1 + length2 - 1));
         else trunc = 0;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_KS(res1, pol1, pol2, 0);
            zmod_poly_truncate(res1, trunc);
            zmod_poly_mul_trunc_n(res2, pol1, pol2, trunc);
            
            result &= zmod_poly_equal(res1, res2);
            
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_trunc_n_precache()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 30) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(2000)+2000;
         unsigned long length2 = randint(2000)+2000;
         
         unsigned long trunc;
         trunc = randint(length1 + length2 - 2000)+2000;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_precache_t pre;
            zmod_poly_mul_trunc_n_precache_init(pre, pol2, 0, trunc);
            zmod_poly_mul_trunc_n_precache(res1, pol1, pre, trunc);
            if (pol1->length > pol2->length) zmod_poly_mul_KS_trunc(res2, pol1, pol2, 0, trunc);
            else zmod_poly_mul_KS_trunc(res2, pol2, pol1, 0, trunc);
            zmod_poly_mul_precache_clear(pre);

            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test__zmod_poly_mul_KS_precache()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 30) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(2000)+2000;
         unsigned long length2 = randint(2000)+2000;
                  
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_precache_t pre;
            zmod_poly_mul_precache_init(pre, pol2, 0, length1);
				zmod_poly_fit_length(res1, length1 + length2 - 1);
            _zmod_poly_mul_KS_precache(res1, pol1, pre, 0);
            if (pol1->length > pol2->length) zmod_poly_mul_KS(res2, pol1, pol2, 0);
            else zmod_poly_mul_KS(res2, pol2, pol1, 0);
            zmod_poly_mul_precache_clear(pre);

            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_precache()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 30) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(2000)+2000;
         unsigned long length2 = randint(2000)+2000;
                  
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_precache_t pre;
            zmod_poly_mul_precache_init(pre, pol2, 0, length1);
            zmod_poly_mul_precache(res1, pol1, pre);
            if (pol1->length > pol2->length) zmod_poly_mul_KS(res2, pol1, pol2, 0);
            else zmod_poly_mul_KS(res2, pol2, pol1, 0);
            zmod_poly_mul_precache_clear(pre);

            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

#if USE_MIDDLE_PRODUCT
int test_zmod_poly_mul_KS_middle()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 30) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(1000)+1000;
         unsigned long length2 = (length1+1)/2;
         
         unsigned long trunc = length1;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_KS_trunc(res1, pol1, pol2, 0, trunc);
            for (unsigned long i = 0; i < trunc/2; i++)
               res1->coeffs[i] = 0L;
            zmod_poly_mul_KS_middle(res2, pol1, pol2, 0, trunc);
            
            result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_middle_precache()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 30) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(1000)+1000;
         unsigned long length2 = (length1+1)/2;
         
         unsigned long trunc = length1;
			zmod_poly_precache_t pre;
         
#if DEBUG
         printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

         randpoly(pol1, length1, modulus);
         randpoly(pol2, length2, modulus);
         
         zmod_poly_mul_KS_trunc(res1, pol1, pol2, 0, trunc);
         for (unsigned long i = 0; i < trunc/2; i++)
            res1->coeffs[i] = 0L;
		   zmod_poly_mul_trunc_n_precache_init(pre, pol2, 0, length2);
         zmod_poly_mul_middle_precache(res2, pol1, pre, trunc);
         zmod_poly_mul_precache_clear(pre);
            
         result &= zmod_poly_equal(res1, res2);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
         }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}
#endif

int test_zmod_poly_mul_classical_trunc_left()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
         unsigned long trunc;
         if (length1 + length2 > 1) trunc = randint(2*(length1 + length2 - 1));
         else trunc = 0;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            zmod_poly_mul_classical_trunc_left(res2, pol1, pol2, trunc);
            
            for (unsigned long i = trunc; i < res1->length; i++)
               if (res1->coeffs[i] != res2->coeffs[i]) result = 0;
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_mul_trunc_left_n()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         
         unsigned long trunc;
         if (length1 + length2 > 1) trunc = randint(2*(length1 + length2 - 1));
         else trunc = 0;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld, trunc = %ld\n", bits, length1, length2, modulus, trunc);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_classical(res1, pol1, pol2);
            zmod_poly_mul_trunc_left_n(res2, pol1, pol2, trunc);
            
            for (unsigned long i = trunc; i < res1->length; i++)
               if (res1->coeffs[i] != res2->coeffs[i]) result = 0;
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
            }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_scalar_mul()
{
   int result = 1;
   zmod_poly_t pol1, res1;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
         unsigned long scalar = randint(modulus-1) + 1;
         unsigned long scalar_inv = z_invert(scalar, modulus);
         
#if DEBUG
         printf("length1 = %ld, bits = %ld, modulus = %ld, scalar = %ld, scalar_inv = %ld\n", length1, bits, modulus, scalar, scalar_inv);
#endif
         
         randpoly(pol1, length1, modulus);
         
         zmod_poly_scalar_mul(res1, pol1, scalar);
         zmod_poly_scalar_mul(res1, res1, scalar_inv);
         
         result &= zmod_poly_equal(res1, pol1);
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1); 
   }
   
   return result;
}

int test_zmod_poly_make_monic()
{
   int result = 1;
   zmod_poly_t pol1, res1;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
         unsigned long scalar;
         
#if DEBUG
         printf("length1 = %ld, bits = %ld, modulus = %ld\n", length1, bits, modulus);
#endif
         
         randpoly(pol1, length1, modulus);
         
         if (pol1->length)
			{
				scalar = pol1->coeffs[pol1->length - 1];
				zmod_poly_make_monic(res1, pol1);

				result &= (res1->coeffs[res1->length - 1] == 1);

				zmod_poly_scalar_mul(res1, res1, scalar);
            
            result &= zmod_poly_equal(res1, pol1);
			}

#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1); 
   }

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         
         unsigned long scalar;
         
#if DEBUG
         printf("length1 = %ld, bits = %ld, modulus = %ld\n", length1, bits, modulus);
#endif
         
         randpoly(pol1, length1, modulus);
         
         if (pol1->length)
			{
				scalar = pol1->coeffs[pol1->length - 1];
            zmod_poly_set(res1, pol1);
				zmod_poly_make_monic(pol1, pol1);

				result &= (pol1->coeffs[pol1->length - 1] == 1);

				zmod_poly_scalar_mul(pol1, pol1, scalar);
            
            result &= zmod_poly_equal(res1, pol1);
			}

#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1); 
   }
   return result;
}

int test_zmod_poly_divrem_classical()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q, R;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      zmod_poly_init(R, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_divrem_classical(Q, R, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
            zmod_poly_print(R); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
      zmod_poly_clear(R); 
   }
   
   return result;
}

int test_zmod_poly_divrem_newton()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q, R;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      zmod_poly_init(R, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_divrem_newton(Q, R, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
            zmod_poly_print(R); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
      zmod_poly_clear(R); 
   }
   
   return result;
}

int test_zmod_poly_divrem()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q, R;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      zmod_poly_init(R, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_divrem(Q, R, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
            zmod_poly_print(R); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
      zmod_poly_clear(R); 
   }
   
   return result;
}

int test_zmod_poly_div_classical()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_div_classical(Q, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
   }
   
   return result;
}

int test_zmod_poly_div()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_div(Q, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
   }
   
   return result;
}

int test_zmod_poly_divrem_divconquer()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q, R;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      zmod_poly_init(R, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(500);
         unsigned long length2 = randint(500);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_divrem_divconquer(Q, R, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
            zmod_poly_print(R); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
      zmod_poly_clear(R); 
   }
   
   return result;
}

int test_zmod_poly_div_divconquer()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_div_divconquer(Q, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
   }
   
   return result;
}

int test_zmod_poly_div_divconquer_recursive()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q, BQ;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 400) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      zmod_poly_init(BQ, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100);
         unsigned long length2 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_div_divconquer(Q, res1, pol2);
					zmod_poly_mul(pol2, BQ, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
               result &= zmod_poly_equal(BQ, pol2);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
            zmod_poly_print(BQ); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
      zmod_poly_clear(BQ); 
   }
   
   return result;
}

int test_zmod_poly_newton_invert_basecase()
{
   zmod_poly_t poly, poly2, poly3;
   int result = 1;
   unsigned long bits, length, n;
   
   for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);
      zmod_poly_init(poly3, modulus);
           
      length = random_ulong(64)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

         do randpoly(poly, length, modulus); 
         while (poly->length == 0);
      
         zmod_poly_set_coeff_ui(poly, poly->length - 1, 1L);
      
         n = randint(poly->length) + 1;
      
         zmod_poly_newton_invert_basecase(poly2, poly, n);
      
         zmod_poly_mul(poly3, poly, poly2);
           
         for (unsigned long i = 0; i < n - 1; i++)
         {
            result &= (poly3->coeffs[i+poly3->length-n] == 0L);
         }
         result &= (poly3->coeffs[poly3->length-1] == 1L);
      
#if DEBUG
      if (!result)
      {
         zmod_poly_print(poly); printf("\n");
         zmod_poly_print(poly2); printf("\n");
         zmod_poly_print(poly3); printf("\n");
      }
#endif
      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
      zmod_poly_clear(poly3);
   }
      
   return result; 
}

int test_zmod_poly_newton_invert()
{
   zmod_poly_t poly, poly2, poly3;
   int result = 1;
   unsigned long bits, length;
   
   for (unsigned long count1 = 0; (count1 < 30) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);
      zmod_poly_init(poly3, modulus);
            
      length = random_ulong(5000)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         do randpoly(poly, length, modulus); 
         while ((poly->length == 0) || (poly->coeffs[0] == 0L));
            
         zmod_poly_newton_invert(poly2, poly, length);
      
         zmod_poly_mul_trunc_n(poly3, poly, poly2, length);
            
         result &= (poly3->length == 1);
         result &= (poly3->coeffs[0] == 1L);
      }

      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
      zmod_poly_clear(poly3);
   }
   
   return result; 
}

int test_zmod_poly_div_series()
{
   zmod_poly_t poly, poly2, poly3, poly4;
   int result = 1;
   unsigned long bits, length;
   
   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);
      zmod_poly_init(poly3, modulus);
      zmod_poly_init(poly4, modulus);
            
      length = randint(200)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

         do randpoly(poly, length, modulus); 
         while ((poly->length == 0) || (poly->coeffs[0] == 0L));
      
         randpoly(poly2, length, modulus); 
                
         zmod_poly_div_series(poly3, poly2, poly, length);
      
         zmod_poly_mul_trunc_n(poly4, poly3, poly, length);
            
         result = zmod_poly_equal(poly4, poly2);
      
#if DEBUG
      if (!result)
      {
         zmod_poly_print(poly); printf("\n");
         zmod_poly_print(poly2); printf("\n");
         zmod_poly_print(poly3); printf("\n");
         zmod_poly_print(poly4); printf("\n");
      }
#endif
      
      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
      zmod_poly_clear(poly3);
      zmod_poly_clear(poly4);
   }
   
   return result; 
}

int test_zmod_poly_div_newton()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, Q;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(Q, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(200);
         unsigned long length2 = randint(200);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
         unsigned log_length = 0L;
         while ((1L<<log_length) < FLINT_MAX(length1, length2)) log_length++;
         
            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul(res1, pol1, pol2);
            if (pol2->length)
            {
               zmod_poly_div_newton(Q, res1, pol2);
         
               result &= zmod_poly_equal(Q, pol1);
            }
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(Q); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1); 
      zmod_poly_clear(Q); 
   }
   
   return result;
}

int test_zmod_poly_gcd_euclidean()
{
   int result = 1;
   zmod_poly_t pol1, pol2, pol3, res1, res2, res3;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd_euclidean(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_gcd_euclidean(res1, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (res1->length == pol3->length));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd_euclidean(res1, pol1, pol2);
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_gcd_euclidean(pol1, pol1, pol2);//alias testing
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol1, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (pol1->length == pol3->length));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd_euclidean(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_gcd_euclidean(pol2, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol2, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (pol2->length == pol3->length));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
   return result;
}

int test_zmod_poly_gcd_invert()
{
   int result = 1;
   zmod_poly_t pol1, pol2, pol3, res1, res2, res3, res4;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      zmod_poly_init(res4, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+2;
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               if (pol2->length != 0) zmod_poly_divrem_newton(res2, pol1, pol1, pol2);
               else zmod_poly_zero(pol1);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while ((res1->length != 1) || (pol1->length == 0));

            zmod_poly_gcd_invert(res1, pol1, pol2);
            zmod_poly_mul(res2, res1, pol1);
            zmod_poly_divrem_newton(res4, res3, res2, pol2);
            
            result &= (res3->length == 1);
         
#if DEBUG2
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
            zmod_poly_print(res3); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
      zmod_poly_clear(res4); 
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      zmod_poly_init(res4, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+2;
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               if (pol2->length != 0) zmod_poly_divrem_newton(res2, pol1, pol1, pol2);
               else zmod_poly_zero(pol1);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while ((res1->length != 1) || (pol1->length == 0));

            zmod_poly_set(res1,pol1);//alias testing
            zmod_poly_gcd_invert(pol1, pol1, pol2);
            zmod_poly_mul(res2, res1, pol1);
            zmod_poly_divrem_newton(res4, res3, res2, pol2);
            
            result &= (res3->length == 1);
         
#if DEBUG2
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
            zmod_poly_print(res3); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
      zmod_poly_clear(res4); 
   }   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      zmod_poly_init(res4, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+2;
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               if (pol2->length != 0) zmod_poly_divrem_newton(res2, pol1, pol1, pol2);
               else zmod_poly_zero(pol1);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while ((res1->length != 1) || (pol1->length == 0));

            zmod_poly_set(res1,pol2);
            zmod_poly_gcd_invert(pol2, pol1, pol2);
            zmod_poly_mul(res2, pol2, pol1);
            zmod_poly_divrem_newton(res4, res3, res2, res1);
            
            result &= (res3->length == 1);
         
#if DEBUG2
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(res2); printf("\n\n");
            zmod_poly_print(res3); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
      zmod_poly_clear(res4); 
   }
   return result;
}

int test_zmod_poly_xgcd()
{
   int result = 1;
   zmod_poly_t s, t, pol1, pol2, pol3, res1, res2, res3;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_xgcd(res1, s, t, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(s, s, pol1);
            zmod_poly_mul(t, t, pol2);
            zmod_poly_add(s, s, t);

            result &= ((res3->length == 0) && (res1->length == pol3->length) && zmod_poly_equal(res1, s));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
            
            zmod_poly_set(res1,pol1);
            zmod_poly_xgcd(pol1, s, t, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol1, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(s, s, res1);
            zmod_poly_mul(t, t, pol2);
            zmod_poly_add(s, s, t);

            result &= ((res3->length == 0) && (pol1->length == pol3->length) && zmod_poly_equal(pol1, s));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
            
            zmod_poly_set(res1,pol2);
            zmod_poly_xgcd(pol2, s, t, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol2, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(s, s, pol1);
            zmod_poly_mul(t, t, res1);
            zmod_poly_add(s, s, t);

            result &= ((res3->length == 0) && (pol2->length == pol3->length) && zmod_poly_equal(pol2, s));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_set(s,pol1);
            zmod_poly_xgcd(res1, pol1, t, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(pol1, pol1, s);
            zmod_poly_mul(t, t, pol2);
            zmod_poly_add(pol1, pol1, t);

            result &= ((res3->length == 0) && (res1->length == pol3->length) && zmod_poly_equal(res1, pol1));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_set(s,pol2);
            zmod_poly_xgcd(res1, pol2, t, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(pol2, pol2, pol1);
            zmod_poly_mul(t, t, s);
            zmod_poly_add(pol2, pol2, t);

            result &= ((res3->length == 0) && (res1->length == pol3->length) && zmod_poly_equal(res1, pol2));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_set(t,pol1);
            zmod_poly_xgcd(res1, s, pol1, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(s, s, t);
            zmod_poly_mul(pol1, pol1, pol2);
            zmod_poly_add(s, s, pol1);

            result &= ((res3->length == 0) && (res1->length == pol3->length) && zmod_poly_equal(res1, s));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(s, modulus);
      zmod_poly_init(t, modulus);
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_set(t,pol2);
            zmod_poly_xgcd(res1, s, pol2, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            zmod_poly_mul(s, s, pol1);
            zmod_poly_mul(pol2, pol2, t);
            zmod_poly_add(s, s, pol2);

            result &= ((res3->length == 0) && (res1->length == pol3->length) && zmod_poly_equal(res1, s));         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(res1); printf("\n\n");
            zmod_poly_print(s); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(s);
      zmod_poly_clear(t);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
   return result;
}

int test_zmod_poly_resultant_euclidean()
{
   int result = 1;
   zmod_poly_t pol1, pol2, lin;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 500) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(lin, modulus);

      unsigned long r1 = randint(FLINT_MIN(10, modulus)); 
      unsigned long r2 = randint(FLINT_MIN(10, modulus)); 
      unsigned long * roots1 = flint_stack_alloc(r1+1);
      unsigned long * roots2 = flint_stack_alloc(r2+1);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
#if DEBUG
            printf("r1 = %ld, r2 = %ld, modulus = %ld\n", r1, r2, modulus);
#endif

            int exists;

            for (unsigned long i = 0; i < r1; )
            {
               exists = 0;
               unsigned long n = randint(modulus);
               for (unsigned long j = 0; j < i; j++)
                  if (roots1[j] == n) exists = 1;
               if (!exists) 
               {
                  roots1[i] = n;
                  i++;
               }
            }
            
            for (unsigned long i = 0; i < r2; )
            {
               exists = 0;
               unsigned long n = randint(modulus);
               for (unsigned long j = 0; j < i; j++)
                  if (roots2[j] == n) exists = 1;
               if (!exists) 
               {
                  roots2[i] = n;
                  i++;
               }
            }
            
            zmod_poly_set_coeff_ui(pol1, 0, 1);
            pol1->length = 1;
            zmod_poly_set_coeff_ui(pol2, 0, 1);
            pol2->length = 1;

            zmod_poly_set_coeff_ui(lin, 1, 1L);
            lin->length = 2;
            
            for (unsigned long i = 0; i < r1; i++)
            {
               zmod_poly_set_coeff_ui(lin, 0, z_submod(0, roots1[i], modulus));
               zmod_poly_mul(pol1, pol1, lin);
            }

            for (unsigned long i = 0; i < r2; i++)
            {
               zmod_poly_set_coeff_ui(lin, 0, z_submod(0, roots2[i], modulus));
               zmod_poly_mul(pol2, pol2, lin);
            }

            unsigned long res1, res2;

            res1 = 1;
            for (unsigned long i = 0; i < r1; i++)
            {
               for (unsigned long j = 0; j < r2; j++)
               {
                  res1 = z_mulmod2_precomp(res1, z_submod(roots1[i], roots2[j], modulus), modulus, pol1->p_inv);
               }
            }
 
            res2 = zmod_poly_resultant_euclidean(pol1, pol2);

            result = (res1 == res2);
         
#if DEBUG
            if (!result)
            {
               printf("res1 = %ld, res2 = %ld\n", res1, res2);
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               for (unsigned long i = 0; i < r1; i++) printf("%ld, ", roots1[i]); 
               printf("\n");
               for (unsigned long i = 0; i < r2; i++) printf("%ld, ", roots2[i]); 
               printf("\n");
            }
#endif
      
      }
      
      flint_stack_release();
      flint_stack_release();
      zmod_poly_clear(lin);
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
   }
   
   return result;
}

//simple theoretically failsafe derivative function
void simple_derivative(zmod_poly_t x_primed, zmod_poly_t x)
{
	unsigned long length = zmod_poly_length(x);
	unsigned long p = zmod_poly_modulus(x);
	mpz_t num, index;
	mpz_init(index);
	mpz_init(num);
	pre_inv_t p_inv = z_precompute_inverse(p);

	//First of all, we make a fmpz_poly of the polynomial
	mpz_poly_t f;
	mpz_poly_init(f);
	
	for(unsigned long i = 0; i < length; ++i)
	   mpz_poly_set_coeff_ui(f, i, zmod_poly_get_coeff_ui(x, i));
	
	//Now we take the derivative of that:
	mpz_init(index);
	mpz_init(num);
	for(unsigned long i = 0; i < length; ++i)
	{
	   mpz_poly_get_coeff(num, f, i+1); // mpn_poly returns 0 if i+1 > length - 1 
	   mpz_set_ui(index, i+1);		
	   mpz_mul(num, num, index);
	   mpz_poly_set_coeff(f, i, num);
	}
	mpz_clear(index);
	mpz_clear(num);	
	   
	//zero the polynomial
	zmod_poly_zero(x_primed);

	//now we reduce all the coefficients mod p and reconstruct the polynomial.
	mpz_t coeff, new_coeff;
	//final poly will have degree one less than original so we must ignore the last coefficient
	if(length !=0) length--;
	mpz_init(coeff);
	mpz_init(new_coeff);
    for (unsigned long i = 0; i < length; i++)
	{
	   mpz_poly_get_coeff (coeff , f , i);
	   zmod_poly_set_coeff_ui(x_primed, i, mpz_mod_ui(new_coeff, coeff, p));		
	}
    mpz_clear(coeff);
	mpz_clear(new_coeff);
	mpz_poly_clear(f);
}

int test_zmod_poly_derivative()
{
	int result = 1;
	zmod_poly_t poly1, res1, res2;
	unsigned long bits;
    unsigned long modulus;

	//test random polynomials
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
	{
		bits = randint(FLINT_BITS-2)+2;
		
		do {modulus = randprime(bits);} while (modulus < 2);
		
		for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
		{

			unsigned long length = randint(100)+1;
#if DEBUG
            printf("length = %ld, bits = %ld, modulus = %ld\n", length, bits, modulus);
#endif			
	 		//init
			zmod_poly_init(poly1, modulus);
			zmod_poly_init(res1, modulus);
			zmod_poly_init(res2, modulus);
			
			//random poly test
			randpoly(poly1, length, modulus);
			zmod_poly_derivative(res1, poly1);
			simple_derivative(res2, poly1);
			result &= zmod_poly_equal(res1, res2);

			if (!result)
			{
				printf("\npoly1 = ");zmod_poly_print(poly1); printf("\n\n");
				printf("res1 = "); zmod_poly_print(res1); printf("\n\n");
				printf("res2 = "); zmod_poly_print(res2); printf("\n\n");
			}
			
			zmod_poly_clear(poly1);
			zmod_poly_clear(res1);
			zmod_poly_clear(res2);
		}
	}
	
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
	{
	   bits = randint(FLINT_BITS-2)+2;
		
	   do {modulus = randprime(bits);} while (modulus < 2);
	   //create zero polynomials
#if DEBUG
            printf("bits = %ld, modulus = %ld\n", bits, modulus);
#endif			
	   zmod_poly_init(poly1, modulus);
	   zmod_poly_zero(poly1);
	   zmod_poly_init(res1, modulus);
	   zmod_poly_init(res2, modulus);
	   zmod_poly_zero(res2);
			
	   zmod_poly_derivative(res1, poly1);
	   result &= zmod_poly_equal(res1, res2);
	   if(!result) printf("Failed on zero test\n");
			
	   zmod_poly_clear(poly1);
       zmod_poly_clear(res1);
	   zmod_poly_clear(res2);
	}
			
	//test special case poly1 == poly2
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
	{
		bits = randint(FLINT_BITS-2)+2;
		
		do {modulus = randprime(bits);} while (modulus < 2);
		
		for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
		{

			unsigned long length = randint(100)+1;
#if DEBUG
            printf("length = %ld, bits = %ld, modulus = %ld\n", length, bits, modulus);
#endif			
			zmod_poly_init(poly1, modulus);
			zmod_poly_init(res1, modulus);
			randpoly(poly1, length, modulus);
			
			zmod_poly_derivative(res1, poly1);
			zmod_poly_derivative(poly1, poly1);
			result &= zmod_poly_equal(res1, poly1);
			if(!result)
				printf("failed on equal \n");
			
			zmod_poly_clear(poly1);
			zmod_poly_clear(res1);
		}	
	}
	
	return result;
}

int test_zmod_poly_mulmod()
{
   int result = 1;
   zmod_poly_t pol1, pol2, res1, res2, quot, rem, f;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
      zmod_poly_init(f, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         unsigned long length3 = randint(400)+1;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
            do randpoly(f, length3, modulus);
			while (zmod_poly_is_zero(f));
         
            zmod_poly_mul(res1, pol1, pol2);
            zmod_poly_mulmod(res2, pol1, pol2, f);
			zmod_poly_sub(res1, res1, res2);
			zmod_poly_divrem(quot, rem, res1, f);
            result &= zmod_poly_is_zero(rem);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
               zmod_poly_print(rem); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(quot);  
      zmod_poly_clear(rem);  
      zmod_poly_clear(f);  
   }

   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
      zmod_poly_init(f, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         unsigned long length3 = randint(400)+1;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
            do randpoly(f, length3, modulus);
			while (zmod_poly_is_zero(f));
         
            zmod_poly_mul(res1, pol1, pol2);
            zmod_poly_mulmod(pol1, pol1, pol2, f);
			zmod_poly_sub(res1, res1, pol1);
			zmod_poly_divrem(quot, rem, res1, f);
            result &= zmod_poly_is_zero(rem);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
               zmod_poly_print(rem); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(quot);  
      zmod_poly_clear(rem);  
      zmod_poly_clear(f);  
   }

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
      zmod_poly_init(f, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         unsigned long length3 = randint(400)+1;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
            do randpoly(f, length3, modulus);
			while (zmod_poly_is_zero(f));
         
            zmod_poly_mul(res1, pol1, pol2);
            zmod_poly_mulmod(pol2, pol1, pol2, f);
			zmod_poly_sub(res1, res1, pol2);
			zmod_poly_divrem(quot, rem, res1, f);
            result &= zmod_poly_is_zero(rem);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
               zmod_poly_print(rem); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(quot);  
      zmod_poly_clear(rem);  
      zmod_poly_clear(f);  
   }


   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
      zmod_poly_init(f, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(400);
         unsigned long length2 = randint(400);
         unsigned long length3 = randint(400)+1;
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
            do randpoly(f, length3, modulus);
			while (zmod_poly_is_zero(f));
         
            zmod_poly_set(res2,f);
            zmod_poly_mul(res1, pol1, pol2);
            zmod_poly_mulmod(f, pol1, pol2, f);
			zmod_poly_sub(res1, res1, f);
			zmod_poly_divrem(quot, rem, res1, res2);
            result &= zmod_poly_is_zero(rem);
         
#if DEBUG
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(pol2); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
               zmod_poly_print(rem); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(quot);  
      zmod_poly_clear(rem);  
      zmod_poly_clear(f);  
   }

   
   return result;
}

int test_zmod_poly_powmod()
{
   int result = 1;
   zmod_poly_t pol1, res1, res2, f, temp;
   unsigned long bits;
   long exp;
   
   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(f, modulus);
      zmod_poly_init(temp, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length3 = randint(100)+1;
         
		 exp = randint(30);
         if (randint(2) && (length3 != 1)) exp = -exp;
		 if ((exp == 0) && ((length3 == 1) || (length1 == 0))) exp++;
               
#if DEBUG
            printf("exp = %ld, bits = %ld, length1 = %ld, length3 = %ld, modulus = %ld\n", exp, bits, length1, length3, modulus);
#endif

            do
			{
			   randpoly(pol1, length1, modulus);
               do randpoly(f, length3, modulus);
			   while (zmod_poly_is_zero(f));
			   zmod_poly_divrem(temp, pol1, pol1, f);
			   zmod_poly_gcd(temp, pol1, f);
			} while (((pol1->length == 0) && (exp <= 0)) || ((pol1->length != 0) && (exp < 0L) && (temp->length != 1)));
            
            zmod_poly_powmod(res1, pol1, exp, f);
			zmod_poly_set_coeff_ui(res2, 0, 1L);
			res2->length = 1;
			for (unsigned long i = 0; i < FLINT_ABS(exp); i++)
			{
			   zmod_poly_mulmod(res2, res2, pol1, f);
			}
            
            if (exp >= 0L) result &= zmod_poly_equal(res1, res2);
			else
			{
			   zmod_poly_mulmod(temp, res1, res2, f);
			   result &= (temp->length == 1);
			}
            
#if DEBUG2
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
			   zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
               if (exp < 0L) zmod_poly_print(temp); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(f);  
      zmod_poly_clear(temp);  
   }

   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(f, modulus);
      zmod_poly_init(temp, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length3 = randint(100)+1;
         
		 exp = randint(30);
         if (randint(2) && (length3 != 1)) exp = -exp;
		 if ((exp == 0) && ((length3 == 1) || (length1 == 0))) exp++;
               
#if DEBUG
            printf("exp = %ld, bits = %ld, length1 = %ld, length3 = %ld, modulus = %ld\n", exp, bits, length1, length3, modulus);
#endif

            do
			{
			   randpoly(pol1, length1, modulus);
               do randpoly(f, length3, modulus);
			   while (zmod_poly_is_zero(f));
			   zmod_poly_divrem(temp, pol1, pol1, f);
			   zmod_poly_gcd(temp, pol1, f);
			} while (((pol1->length == 0) && (exp <= 0)) || ((pol1->length != 0) && (exp < 0L) && (temp->length != 1)));
            
            zmod_poly_set(res1,pol1);
            zmod_poly_powmod(pol1, pol1, exp, f);

			zmod_poly_set_coeff_ui(res2, 0, 1L);
			res2->length = 1;
			for (unsigned long i = 0; i < FLINT_ABS(exp); i++)
			{
			   zmod_poly_mulmod(res2, res2, res1, f);
			}
            
            if (exp >= 0L) result &= zmod_poly_equal(pol1, res2);
			else
			{
			   zmod_poly_mulmod(temp, pol1, res2, f);
			   result &= (temp->length == 1);
			}
            
#if DEBUG2
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
			   zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
               if (exp < 0L) zmod_poly_print(temp); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(f);  
      zmod_poly_clear(temp);  
   }

   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(f, modulus);
      zmod_poly_init(temp, modulus);
      
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length3 = randint(100)+1;
         
		 exp = randint(30);
         if (randint(2) && (length3 != 1)) exp = -exp;
		 if ((exp == 0) && ((length3 == 1) || (length1 == 0))) exp++;
               
#if DEBUG
            printf("exp = %ld, bits = %ld, length1 = %ld, length3 = %ld, modulus = %ld\n", exp, bits, length1, length3, modulus);
#endif

            do
			{
			   randpoly(pol1, length1, modulus);
               do randpoly(f, length3, modulus);
			   while (zmod_poly_is_zero(f));
			   zmod_poly_divrem(temp, pol1, pol1, f);
			   zmod_poly_gcd(temp, pol1, f);
			} while (((pol1->length == 0) && (exp <= 0)) || ((pol1->length != 0) && (exp < 0L) && (temp->length != 1)));
            
            zmod_poly_set(res1,f);
            zmod_poly_powmod(f, pol1, exp, f);

			zmod_poly_set_coeff_ui(res2, 0, 1L);
			res2->length = 1;
			for (unsigned long i = 0; i < FLINT_ABS(exp); i++)
			{
			   zmod_poly_mulmod(res2, res2, pol1, res1);
			}
            
            if (exp >= 0L) result &= zmod_poly_equal(f, res2);
			else
			{
			   zmod_poly_mulmod(temp, f, res2, res1);
			   result &= (temp->length == 1);
			}
            
#if DEBUG2
            if (!result)
            {
               zmod_poly_print(pol1); printf("\n\n");
               zmod_poly_print(f); printf("\n\n");
			   zmod_poly_print(res1); printf("\n\n");
               zmod_poly_print(res2); printf("\n\n");
               if (exp < 0L) zmod_poly_print(temp); printf("\n\n");
            }
#endif
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
      zmod_poly_clear(f);  
      zmod_poly_clear(temp);  
   }
   
   return result;
}

int test_zmod_poly_isirreducible()
{
   zmod_poly_t poly, poly2, poly3;
   int result = 1;
   unsigned long bits, length, length2;
   
   for (unsigned long count1 = 0; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(poly, modulus);
      zmod_poly_init(poly2, modulus);
      zmod_poly_init(poly3, modulus);
      
		ulong length = randint(10) + 2;
      do 
		{
			randpoly(poly, length, modulus); 
			zmod_poly_make_monic(poly, poly);
		}
      while ((!zmod_poly_isirreducible(poly)) || (poly->length < 2));
      
      zmod_poly_factor_t factors;
		zmod_poly_factor_init(factors);
		zmod_poly_factor_berlekamp(factors, poly);
		result &= (factors->num_factors == 1);
		if (!result)
		{
			printf("Error: irreducible polynomial should not have non-trivial factors!\n");
			zmod_poly_print(poly); printf("\n");
		}
      zmod_poly_factor_clear(factors);

		ulong length2 = randint(10) + 2;
      
		do 
		{
			randpoly(poly2, length2, modulus); 
         zmod_poly_make_monic(poly2, poly2);
		} while ((!zmod_poly_isirreducible(poly2)) || (poly2->length < 2));

      zmod_poly_mul(poly3, poly, poly2);

      result &= !zmod_poly_isirreducible(poly3);
      if (!result)
		{
			printf("Error: reducible polynomial declared irreducible!\n");
			zmod_poly_print(poly3); printf("\n");
		}

      zmod_poly_clear(poly);
      zmod_poly_clear(poly2);
      zmod_poly_clear(poly3);
   }

   return result; 
}

int test_zmod_poly_factor_square_free()
{
   int result = 1;
   zmod_poly_t pol1, poly, quot, rem;
   zmod_poly_factor_t res;
   unsigned long bits;
   unsigned long modulus;
   ulong exponents[5];

   for (unsigned long count1 = 0; (count1 < 300) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
 	   zmod_poly_init(pol1, modulus);
      zmod_poly_init(poly, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
     
	   zmod_poly_zero(pol1);
		zmod_poly_set_coeff_ui(pol1, 0, 1);
		
		ulong length = randint(7) + 2;
      do 
	   {
	      randpoly(poly, length, modulus); 
		   zmod_poly_make_monic(poly, poly);
	   }
      while ((!zmod_poly_isirreducible(poly)) || (poly->length < 2));
      exponents[0] = z_randprime(5);
		ulong prod1 = exponents[0];
		for (ulong i = 0; i < exponents[0]; i++) zmod_poly_mul(pol1, pol1, poly);
		
	   ulong num_factors = z_randint(5)+1;
	   for (ulong i = 1; i < num_factors; i++)
	   {
		   do 
	      {
	         length = randint(7) + 2;
            randpoly(poly, length, modulus); 
		      zmod_poly_make_monic(poly, poly);
			   if (poly->length) zmod_poly_divrem(quot, rem, pol1, poly);
	      }
         while ((!zmod_poly_isirreducible(poly)) || (poly->length < 2) || (rem->length == 0));
		   do exponents[i] = z_randprime(5);
			while (prod1 % exponents[i] == 0);
         prod1 *= exponents[i];
		   for (ulong j = 0; j < exponents[i]; j++) zmod_poly_mul(pol1, pol1, poly);
	   }
     
		zmod_poly_factor_init(res);
      zmod_poly_factor_square_free(res, pol1);

	   result &= (res->num_factors == num_factors);
		if (result)
		{
			ulong prod2 = 1;
			for (ulong i = 0; i < num_factors; i++)
				prod2 *= res->exponents[i];
         result &= (prod1 == prod2);
		}
		if (!result)
		{
			printf("Error: exponents don't match. Modulus = %ld\n", modulus);
			for (ulong i = 0; i < res->num_factors; i++) printf("%ld ", res->exponents[i]);
			printf("\n");
         for (ulong i = 0; i < num_factors; i++) printf("%ld ", exponents[i]);
			printf("\n");
		}

	   zmod_poly_clear(quot);
	   zmod_poly_clear(rem);
	   zmod_poly_clear(pol1);
	   zmod_poly_factor_clear(res);
   }
}

int test_zmod_poly_factor_berlekamp()
{
   int result = 1;
   zmod_poly_t pol1, poly, quot, rem;
   zmod_poly_factor_t res;
   unsigned long bits;
   unsigned long modulus;

   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
 	   zmod_poly_init(pol1, modulus);
      zmod_poly_init(poly, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
     
	   ulong length = randint(10) + 2;
      do 
	   {
	      randpoly(pol1, length, modulus); 
		   zmod_poly_make_monic(pol1, pol1);
	   }
      while ((!zmod_poly_isirreducible(pol1)) || (pol1->length < 2));
		  
	   ulong num_factors = randint(5) + 1;
	   for (ulong i = 1; i < num_factors; i++)
	   {
		   do 
	      {
	         length = randint(10) + 2;
            randpoly(poly, length, modulus); 
		      zmod_poly_make_monic(poly, poly);
			   if (poly->length) zmod_poly_divrem(quot, rem, pol1, poly);
	      }
         while ((!zmod_poly_isirreducible(poly)) || (poly->length < 2) || (rem->length == 0));
		   zmod_poly_mul(pol1, pol1, poly);
	   }
     
	   zmod_poly_factor_init(res);
      zmod_poly_factor_berlekamp(res, pol1);

	   result = (res->num_factors == num_factors);
      
	   zmod_poly_clear(quot);
	   zmod_poly_clear(rem);
	   zmod_poly_clear(pol1);
	   zmod_poly_factor_clear(res);
   }
}

int test_zmod_poly_factor()
{
   int result = 1;
   zmod_poly_t pol1, poly, quot, rem;
   zmod_poly_factor_t res;
   unsigned long bits;
   unsigned long modulus;
   ulong exponents[5];

   for (unsigned long count1 = 0; (count1 < 300) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
 	   zmod_poly_init(pol1, modulus);
      zmod_poly_init(poly, modulus);
      zmod_poly_init(quot, modulus);
      zmod_poly_init(rem, modulus);
     
	   zmod_poly_zero(pol1);
		zmod_poly_set_coeff_ui(pol1, 0, 1);
		
		ulong length = randint(7) + 2;
      do 
	   {
	      randpoly(poly, length, modulus); 
		   zmod_poly_make_monic(poly, poly);
	   }
      while ((!zmod_poly_isirreducible(poly)) || (poly->length < 2));
      exponents[0] = z_randint(30)+1;
		ulong prod1 = exponents[0];
		for (ulong i = 0; i < exponents[0]; i++) zmod_poly_mul(pol1, pol1, poly);
		
	   ulong num_factors = z_randint(5)+1;
	   for (ulong i = 1; i < num_factors; i++)
	   {
		   do 
	      {
	         length = randint(7) + 2;
            randpoly(poly, length, modulus); 
		      zmod_poly_make_monic(poly, poly);
			   if (poly->length) zmod_poly_divrem(quot, rem, pol1, poly);
	      }
         while ((!zmod_poly_isirreducible(poly)) || (poly->length < 2) || (rem->length == 0));
		   exponents[i] = z_randint(30)+1;
			prod1 *= exponents[i];
		   for (ulong j = 0; j < exponents[i]; j++) zmod_poly_mul(pol1, pol1, poly);
	   }
     
		zmod_poly_factor_init(res);
      ulong lead = zmod_poly_factor(res, pol1);
      
	   result &= (res->num_factors == num_factors);
		if (!result)
		   printf("Error: number of factors incorrect, %ld, %ld\n", res->num_factors, num_factors);
		

		zmod_poly_t product;
		zmod_poly_init(product, pol1->p);
		zmod_poly_set_coeff_ui(product, 0, 1);
		for (ulong i = 0; i < res->num_factors; i++)
			for (ulong j = 0; j < res->exponents[i]; j++)
				zmod_poly_mul(product, product, res->factors[i]);
		zmod_poly_scalar_mul(product, product, lead);
		result &= zmod_poly_equal(pol1, product);
		if (!result)
		{
			printf("Error: product of factors does not equal original polynomial\n");
			zmod_poly_print(pol1); printf("\n");
			zmod_poly_print(product); printf("\n");
		}
		zmod_poly_clear(product);
      
	   zmod_poly_clear(quot);
	   zmod_poly_clear(rem);
	   zmod_poly_clear(pol1);
	   zmod_poly_factor_clear(res);
   }
}

int test_zmod_poly_2x2_mat_mul_classical_strassen()
{
   int result = 1;
   zmod_poly_2x2_mat_t A, B, R1, R2;
   unsigned long bits, length;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_2x2_mat_init(A, modulus);
      zmod_poly_2x2_mat_init(B, modulus);
      zmod_poly_2x2_mat_init(R1, modulus);
      zmod_poly_2x2_mat_init(R2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length = randint(100)+1;
         
         randpoly(A->a, randint(length), modulus);
         randpoly(A->b, randint(length), modulus);
         randpoly(A->c, randint(length), modulus);
         randpoly(A->d, randint(length), modulus);
         
         randpoly(B->a, randint(length), modulus);
         randpoly(B->b, randint(length), modulus);
         randpoly(B->c, randint(length), modulus);
         randpoly(B->d, randint(length), modulus);
         
			zmod_poly_2x2_mat_mul_classical(R1, A, B);
         zmod_poly_2x2_mat_mul_strassen(R2, A, B);

         result &= zmod_poly_equal(R1->a, R2->a);
         result &= zmod_poly_equal(R1->b, R2->b);
         result &= zmod_poly_equal(R1->c, R2->c);
         result &= zmod_poly_equal(R1->d, R2->d);
         
         if (!result)
         {
				printf("%ld, %ld, %ld, %ld\n", A->a->length, A->b->length, A->c->length, A->d->length);
				printf("%ld, %ld, %ld, %ld\n", B->a->length, B->b->length, B->c->length, B->d->length);
         }
      }
      
      zmod_poly_2x2_mat_clear(A);
      zmod_poly_2x2_mat_clear(B);
      zmod_poly_2x2_mat_clear(R1);
      zmod_poly_2x2_mat_clear(R2);
   }
   
   return result;
}

int test_zmod_poly_2x2_mat_mul()
{
   int result = 1;
   zmod_poly_2x2_mat_t A, B, R1, R2;
   unsigned long bits, length;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-1)+2;
      unsigned long modulus;
      
      do {modulus = randbits(bits);} while (modulus < 2);
      
      zmod_poly_2x2_mat_init(A, modulus);
      zmod_poly_2x2_mat_init(B, modulus);
      zmod_poly_2x2_mat_init(R1, modulus);
      zmod_poly_2x2_mat_init(R2, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length = randint(100)+1;
         
         randpoly(A->a, randint(length), modulus);
         randpoly(A->b, randint(length), modulus);
         randpoly(A->c, randint(length), modulus);
         randpoly(A->d, randint(length), modulus);
         
         randpoly(B->a, randint(length), modulus);
         randpoly(B->b, randint(length), modulus);
         randpoly(B->c, randint(length), modulus);
         randpoly(B->d, randint(length), modulus);
         
			zmod_poly_2x2_mat_mul_strassen(R2, A, B);
         zmod_poly_2x2_mat_mul(R1, A, B);
         
         result &= zmod_poly_equal(R1->a, R2->a);
         result &= zmod_poly_equal(R1->b, R2->b);
         result &= zmod_poly_equal(R1->c, R2->c);
         result &= zmod_poly_equal(R1->d, R2->d);
         
         if (!result)
         {
				printf("%ld, %ld, %ld, %ld\n", A->a->length, A->b->length, A->c->length, A->d->length);
				printf("%ld, %ld, %ld, %ld\n", B->a->length, B->b->length, B->c->length, B->d->length);
         }
      }
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length = randint(100)+1;
         
         randpoly(A->a, randint(length), modulus);
         randpoly(A->b, randint(length), modulus);
         randpoly(A->c, randint(length), modulus);
         randpoly(A->d, randint(length), modulus);
         
         randpoly(B->a, randint(length), modulus);
         randpoly(B->b, randint(length), modulus);
         randpoly(B->c, randint(length), modulus);
         randpoly(B->d, randint(length), modulus);
         
			zmod_poly_2x2_mat_mul_strassen(R2, A, B);
         zmod_poly_2x2_mat_mul(A, A, B);
         
         result &= zmod_poly_equal(A->a, R2->a);
         result &= zmod_poly_equal(A->b, R2->b);
         result &= zmod_poly_equal(A->c, R2->c);
         result &= zmod_poly_equal(A->d, R2->d);
         
         if (!result)
         {
				printf("%ld, %ld, %ld, %ld\n", A->a->length, A->b->length, A->c->length, A->d->length);
				printf("%ld, %ld, %ld, %ld\n", B->a->length, B->b->length, B->c->length, B->d->length);
         }
      }
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length = randint(100)+1;
         
         randpoly(A->a, randint(length), modulus);
         randpoly(A->b, randint(length), modulus);
         randpoly(A->c, randint(length), modulus);
         randpoly(A->d, randint(length), modulus);
         
         randpoly(B->a, randint(length), modulus);
         randpoly(B->b, randint(length), modulus);
         randpoly(B->c, randint(length), modulus);
         randpoly(B->d, randint(length), modulus);
         
			zmod_poly_2x2_mat_mul_strassen(R2, A, B);
         zmod_poly_2x2_mat_mul(B, A, B);
         
         result &= zmod_poly_equal(B->a, R2->a);
         result &= zmod_poly_equal(B->b, R2->b);
         result &= zmod_poly_equal(B->c, R2->c);
         result &= zmod_poly_equal(B->d, R2->d);
         
         if (!result)
         {
				printf("%ld, %ld, %ld, %ld\n", A->a->length, A->b->length, A->c->length, A->d->length);
				printf("%ld, %ld, %ld, %ld\n", B->a->length, B->b->length, B->c->length, B->d->length);
         }
      }
      
      zmod_poly_2x2_mat_clear(A);
      zmod_poly_2x2_mat_clear(B);
      zmod_poly_2x2_mat_clear(R1);
      zmod_poly_2x2_mat_clear(R2);
   }
   
   return result;
}

int test_zmod_poly_half_gcd()
{
   int result = 1;
   zmod_poly_t pol1, pol2, pol3, res1;
   zmod_poly_2x2_mat_t Q;
	unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_2x2_mat_init(Q, modulus);
				
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length2 = randint(50)+1;
         unsigned long length1 = length2 + randint(50)+1;
         unsigned long length3 = randint(50)+1;
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while ((res1->length != 1) || (pol2->length >= pol1->length) || (pol2->length == 0));

            do
				{
					randpoly(pol3, length3, modulus);
				} while (pol3->length == 0);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
           
            ulong d_red = pol1->length/2 + 1;
				
				zmod_poly_half_gcd(Q, pol1, pol2);
            zmod_poly_swap(Q->a, Q->d);
				zmod_poly_neg(Q->b, Q->b);
				zmod_poly_neg(Q->c, Q->c);

				zmod_poly_mul(Q->a, Q->a, pol1);
            zmod_poly_mul(Q->b, Q->b, pol2);
            zmod_poly_add(Q->a, Q->a, Q->b);

            zmod_poly_mul(Q->c, Q->c, pol1);
            zmod_poly_mul(Q->d, Q->d, pol2);
            zmod_poly_add(Q->c, Q->c, Q->d);

            result &= (Q->a->length >= pol1->length/2 + 1);
				result &= (pol1->length/2 + 1 > Q->c->length);
#if DEBUG2
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
				printf("deg(c') = %ld, ceil(deg(a)/2) = %ld, deg(d') = %ld\n", Q->a->length, pol1->length/2 + 1, Q->c->length);
         }
#endif
      }
      
      zmod_poly_2x2_mat_clear(Q);
		zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
   }
   
   return result;
}

int test_zmod_poly_half_gcd_iter()
{
   int result = 1;
   zmod_poly_t pol1, pol2, pol3, res1;
   zmod_poly_2x2_mat_t Q;
	unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_2x2_mat_init(Q, modulus);
				
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length2 = randint(50)+1;
         unsigned long length1 = length2 + randint(50)+1;
         unsigned long length3 = randint(50)+1;
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while ((res1->length != 1) || (pol2->length >= pol1->length) || (pol2->length == 0));

            do
				{
					randpoly(pol3, length3, modulus);
				} while (pol3->length == 0);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
           
            ulong d_red = pol1->length/2 + 1;
				
				zmod_poly_half_gcd_iter(Q, pol1, pol2);
            zmod_poly_swap(Q->a, Q->d);
				zmod_poly_neg(Q->b, Q->b);
				zmod_poly_neg(Q->c, Q->c);

				zmod_poly_mul(Q->a, Q->a, pol1);
            zmod_poly_mul(Q->b, Q->b, pol2);
            zmod_poly_add(Q->a, Q->a, Q->b);

            zmod_poly_mul(Q->c, Q->c, pol1);
            zmod_poly_mul(Q->d, Q->d, pol2);
            zmod_poly_add(Q->c, Q->c, Q->d);

            result &= (Q->a->length >= pol1->length/2 + 1);
				result &= (pol1->length/2 + 1 > Q->c->length);
#if DEBUG2
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
				printf("deg(c') = %ld, ceil(deg(a)/2) = %ld, deg(d') = %ld\n", Q->a->length, pol1->length/2 + 1, Q->c->length);
         }
#endif
      }
      
      zmod_poly_2x2_mat_clear(Q);
		zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
   }
   
   return result;
}

int test_zmod_poly_gcd_hgcd()
{
   int result = 1;
   zmod_poly_t pol1, pol2, pol3, res1, res2, res3;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length2 = randint(100)+1;
         unsigned long length1 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
				zmod_poly_gcd_hgcd(res1, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (res1->length == pol3->length));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length2 = randint(100)+1;
         unsigned long length1 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);

            zmod_poly_set(res1,pol1);         
				zmod_poly_gcd_hgcd(pol1, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol1, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (pol1->length == pol3->length));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }

   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         unsigned long length2 = randint(100)+1;
         unsigned long length1 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);

            zmod_poly_set(res1,pol2);         
				zmod_poly_gcd_hgcd(pol2, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol2, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (pol2->length == pol3->length));
         
#if DEBUG
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
#endif
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
   return result;
}

int test_zmod_poly_gcd()
{
   int result = 1;
   zmod_poly_t pol1, pol2, pol3, res1, res2, res3;
   unsigned long bits;
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd_euclidean(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_gcd(res1, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (res1->length == pol3->length));
         
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
   // alias both inputs
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            randpoly(pol1, length1, modulus);
            
            zmod_poly_gcd(res1, pol1, pol1);
            if (pol1->length != 0) zmod_poly_divrem_newton(res2, res3, res1, pol1);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (res1->length == pol1->length));
         
         if (!result)
         {
            zmod_poly_print(pol1); printf("\n\n");
            zmod_poly_print(res1); printf("\n\n");
         }
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }

	// alias first input and output
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd_euclidean(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_gcd(pol1, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol1, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (pol1->length == pol3->length));
         
         if (!result)
         {
            zmod_poly_print(pol3); printf("\n\n");
            zmod_poly_print(pol1); printf("\n\n");
         }
      }

		zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }

      
	// alias second input and output
	for (unsigned long count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = randint(FLINT_BITS-2)+2;
      unsigned long modulus;
      
      do {modulus = randprime(bits);} while (modulus < 2);
      
      zmod_poly_init(pol1, modulus);
      zmod_poly_init(pol2, modulus);
      zmod_poly_init(pol3, modulus);
      zmod_poly_init(res1, modulus);
      zmod_poly_init(res2, modulus);
      zmod_poly_init(res3, modulus);
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      {
         unsigned long length1 = randint(100)+1;
         unsigned long length2 = randint(100)+1;
         unsigned long length3 = randint(100);
         
#if DEBUG
         printf("length1 = %ld, length2 = %ld, bits = %ld, modulus = %ld\n", length1, length2, bits, modulus);
#endif
         
            do 
            {
               randpoly(pol1, length1, modulus);
               randpoly(pol2, length2, modulus);
               zmod_poly_gcd_euclidean(res1, pol1, pol2);          
            } while (res1->length != 1);

            randpoly(pol3, length3, modulus);
            zmod_poly_mul(pol1, pol1, pol3);
            zmod_poly_mul(pol2, pol2, pol3);
         
            zmod_poly_gcd(pol2, pol1, pol2);
            if (pol3->length != 0) zmod_poly_divrem_newton(res2, res3, pol2, pol3);
            else zmod_poly_zero(res3);

            result &= ((res3->length == 0) && (pol2->length == pol3->length));
         
         if (!result)
         {
            zmod_poly_print(pol3); printf("\n\n");
            zmod_poly_print(pol2); printf("\n\n");
         }
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(pol3);
      zmod_poly_clear(res1); 
      zmod_poly_clear(res2); 
      zmod_poly_clear(res3); 
   }
   
   return result;
}

void zmod_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
	RUN_TEST(zmod_poly_fprint_fread); 
#endif
   RUN_TEST(zmod_poly_to_from_string); 
   RUN_TEST(__zmod_poly_normalise); 
   RUN_TEST(zmod_poly_truncate); 
   RUN_TEST(zmod_poly_reverse); 
   RUN_TEST(zmod_poly_add_no_red); 
   RUN_TEST(zmod_poly_addsub); 
   RUN_TEST(zmod_poly_neg); 
   RUN_TEST(zmod_poly_bits); 
   RUN_TEST(zmod_poly_shift); 
   RUN_TEST(zmod_poly_swap); 
   RUN_TEST(zmod_poly_setequal); 
   RUN_TEST(_zmod_poly_attach); 
   RUN_TEST(_zmod_poly_attach_truncate); 
   RUN_TEST(_zmod_poly_attach_shift); 
   RUN_TEST(zmod_poly_derivative); 
   RUN_TEST(zmod_poly_getset_coeff); 
   RUN_TEST(zmod_poly_mul); 
	RUN_TEST(zmod_poly_mul_trunc_n); 
   RUN_TEST(zmod_poly_mul_classicalKS); 
   RUN_TEST(zmod_poly_sqr); 
   RUN_TEST(zmod_poly_sqr_classicalKS); 
   RUN_TEST(zmod_poly_mul_classical_trunc); 
   RUN_TEST(zmod_poly_mul_classical_trunc_left); 
   RUN_TEST(zmod_poly_mul_trunc_left_n); 
   RUN_TEST(zmod_poly_mul_KS_trunc); 
#if USE_MIDDLE_PRODUCT
   RUN_TEST(zmod_poly_mul_KS_middle); 
   RUN_TEST(zmod_poly_mul_middle_precache); 
#endif
   RUN_TEST(_zmod_poly_mul_KS_precache); 
   RUN_TEST(zmod_poly_mul_precache); 
   RUN_TEST(zmod_poly_mul_trunc_n_precache); 
   RUN_TEST(zmod_poly_scalar_mul); 
   RUN_TEST(zmod_poly_make_monic); 
   RUN_TEST(zmod_poly_divrem_classical); 
   RUN_TEST(zmod_poly_divrem_newton); 
   RUN_TEST(zmod_poly_divrem); 
   RUN_TEST(zmod_poly_div_classical); 
   RUN_TEST(zmod_poly_div); 
   RUN_TEST(zmod_poly_divrem_divconquer); 
   RUN_TEST(zmod_poly_div_divconquer); 
   RUN_TEST(zmod_poly_div_divconquer_recursive); 
   RUN_TEST(zmod_poly_newton_invert_basecase);
   RUN_TEST(zmod_poly_newton_invert); 
   RUN_TEST(zmod_poly_div_series); 
   RUN_TEST(zmod_poly_div_newton); 
   RUN_TEST(zmod_poly_gcd_euclidean); 
	RUN_TEST(zmod_poly_half_gcd_iter); 
   RUN_TEST(zmod_poly_gcd_hgcd); 
   RUN_TEST(zmod_poly_half_gcd); 
   RUN_TEST(zmod_poly_gcd); 
   RUN_TEST(zmod_poly_gcd_invert); 
   RUN_TEST(zmod_poly_xgcd); 
	RUN_TEST(zmod_poly_resultant_euclidean);
   RUN_TEST(zmod_poly_mulmod); 
   RUN_TEST(zmod_poly_powmod); 
   RUN_TEST(zmod_poly_isirreducible); 
   RUN_TEST(zmod_poly_factor_berlekamp); 
   RUN_TEST(zmod_poly_factor_square_free); 
   RUN_TEST(zmod_poly_factor); 
   RUN_TEST(zmod_poly_2x2_mat_mul_classical_strassen); 
   RUN_TEST(zmod_poly_2x2_mat_mul); 
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   zmod_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}
