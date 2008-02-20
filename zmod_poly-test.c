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

Copyright (C) 2007, William Hart and David Harvey

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

#define VARY_BITS 0
#define SPARSE 0

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
      zmod_poly_ensure_alloc(poly, 1);
      poly->length = 0;
      return;
   }
              
   zmod_poly_ensure_alloc(poly, length);
   
   for (unsigned long i = 0; i < length; i++)
      poly->coeffs[i] = randint(n);
   poly->length = length;
      
   zmod_poly_normalise(poly);
} 

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
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
         
         zmod_poly_lshift(res, pol1, shift);
         zmod_poly_rshift(res, res, shift);
         
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
         zmod_poly_set_coeff(pol1, num, coeff);
         
         result &= (coeff == zmod_poly_get_coeff(pol1, num));
         
         if (num + 1 > length1) 
         {
            zmod_poly_set_coeff(pol1, num, 0);
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

int test_zmod_poly_mul_naiveKS()
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
         
         unsigned log_length = 0L;
         while ((1L<<log_length) < FLINT_MIN(length1, length2)) log_length++;
         
         if (2*FLINT_BIT_COUNT(modulus) + log_length <= 2*FLINT_BITS)
         {
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            randpoly(pol2, length2, modulus);
         
            zmod_poly_mul_naive(res1, pol1, pol2);
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
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(pol2);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}

int test_zmod_poly_sqr_naiveKS()
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
         
         unsigned log_length = 0L;
         while ((1L<<log_length) < length1) log_length++;
         
         if (2*FLINT_BIT_COUNT(modulus) + log_length <= 2*FLINT_BITS)
         {
         
#if DEBUG
            printf("bits = %ld, length1 = %ld, length2 = %ld, modulus = %ld\n", bits, length1, length2, modulus);
#endif

            randpoly(pol1, length1, modulus);
            
            zmod_poly_sqr_naive(res1, pol1);
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
      
      }
      
      zmod_poly_clear(pol1);
      zmod_poly_clear(res1);  
      zmod_poly_clear(res2);  
   }
   
   return result;
}


void fmpz_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
#endif
   RUN_TEST(zmod_poly_addsub); 
   RUN_TEST(zmod_poly_neg); 
   RUN_TEST(zmod_poly_shift); 
   RUN_TEST(zmod_poly_swap); 
   RUN_TEST(zmod_poly_setequal); 
   RUN_TEST(zmod_poly_getset_coeff); 
   RUN_TEST(zmod_poly_mul_naiveKS); 
   RUN_TEST(zmod_poly_sqr_naiveKS); 
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   fmpz_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


