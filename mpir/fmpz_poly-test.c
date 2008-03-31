/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

fmpz_poly-test.c: Test code for fmpz_poly.c and fmpz_poly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "mpir.h"
#include "test_support.h"
#include "fmpz_poly.h"
#include "mpz_poly.h"

#define VARY_BITS 0
#define SIGNS 1
#define SPARSE 0

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

gmp_randstate_t state;

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
void randpoly(mpz_poly_t pol, long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   pol->length = 0;
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
#if SPARSE
          if (randint(10) == 1) mpz_rrandomb(temp, state, bits);
          else mpz_set_ui(temp, 0);
#else
          mpz_rrandomb(temp, state, bits);
#endif
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

int test_fmpz_poly_to_mpz_poly()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 

   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = randint(5) + 1;
      
      fmpz_poly_init(test_fmpz_poly);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(20);
#if DEBUG
          printf("%ld, %ld\n", length, bits);
#endif
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
           
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly);
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = randint(200) + 1;
      
      fmpz_poly_init(test_fmpz_poly);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(200);
#if DEBUG
          printf("%ld, %ld\n", length, bits);
#endif
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
           
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly);
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_add()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, length3, max_length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   mpz_poly_init(test_poly4); 
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = 100;//randint(1000);
      bits2 = 100;//randint(1000);
      bits3 = 200;//randint(1000);
      
      fmpz_poly_init(test_fmpz_poly);
      fmpz_poly_init(test_fmpz_poly2);
      fmpz_poly_init(test_fmpz_poly3);
      fmpz_poly_init(test_fmpz_poly4);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = 100;//randint(100)+1; 
          length2 = 100;//randint(100)+1;        
          length3 = 100;//randint(100)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          randpoly(test_poly4, length3, bits3); 
          for (ulong i = 0; i < 1; i++)
          {
             mpz_poly_clear(test_poly3);
             mpz_poly_init2(test_poly3, length);
             mpz_poly_add(test_poly3, test_poly, test_poly2);
             //mpz_poly_add(test_poly3, test_poly3, test_poly4);
          }

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly4, test_poly4);
          
          for (ulong i = 0; i < 1000; i++)
          {
             fmpz_poly_clear(test_fmpz_poly3);
             fmpz_poly_init2(test_fmpz_poly3, length);
             fmpz_poly_add(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
             //fmpz_poly_add(test_fmpz_poly3, test_fmpz_poly3, test_fmpz_poly4);
          }

          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);

          result = mpz_poly_equal(test_poly4, test_poly3);

      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }

   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   mpz_poly_clear(test_poly4);
   
   return result; 
}

int test_fmpz_poly_sub()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, length3, max_length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   mpz_poly_init(test_poly4); 
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = randint(1000);
      bits2 = randint(1000);
      bits3 = randint(1000);
      
      fmpz_poly_init(test_fmpz_poly);
      fmpz_poly_init(test_fmpz_poly2);
      fmpz_poly_init(test_fmpz_poly3);
      fmpz_poly_init(test_fmpz_poly4);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100)+1; 
          length2 = randint(100)+1;        
          length3 = randint(100)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          randpoly(test_poly4, length3, bits3); 
          for (ulong i = 0; i < 10; i++)
          {
             mpz_poly_clear(test_poly3);
             mpz_poly_init2(test_poly3, length);
             mpz_poly_sub(test_poly3, test_poly, test_poly2);
             mpz_poly_sub(test_poly3, test_poly3, test_poly4);
          }

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly4, test_poly4);
          
          for (ulong i = 0; i < 10; i++)
          {
             fmpz_poly_clear(test_fmpz_poly3);
             fmpz_poly_init2(test_fmpz_poly3, length);
             fmpz_poly_sub(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
             fmpz_poly_sub(test_fmpz_poly3, test_fmpz_poly3, test_fmpz_poly4);
          }

          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);

          result = mpz_poly_equal(test_poly4, test_poly3);

      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }

   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   mpz_poly_clear(test_poly4);
   
   return result; 
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   //RUN_TEST(fmpz_poly_to_mpz_poly);
   RUN_TEST(fmpz_poly_add);
   //RUN_TEST(fmpz_poly_sub);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(state);
   fmpz_poly_test_all();
   gmp_randclear(state);
   
   mpir_stack_cleanup();

   return 0;
}
