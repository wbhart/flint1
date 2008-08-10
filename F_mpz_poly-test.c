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

#define VARY_BITS 1
#define SIGNS 1
#define SPARSE 1 
#define ITER 1

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

void mpz_randpoly(mpz_poly_t pol, long length, unsigned long maxbits)
{
   unsigned long bits;
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

int test_F_mpz_poly_convert()
{
   mpz_poly_t m_poly1, m_poly2;
   F_mpz_poly_t F_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(m_poly1); 
   mpz_poly_init(m_poly2); 

   for (unsigned long count1 = 1; (count1 < 300000*ITER) && (result == 1) ; count1++)
   {
      bits = random_ulong(200) + 1;
      
      F_mpz_poly_init(F_poly);

      length = random_ulong(20);
#if DEBUG
      printf("%ld, %ld\n",length, bits);
#endif
      mpz_randpoly(m_poly1, length, bits);
           
#if DEBUG
      mpz_poly_print_pretty(m_poly1, "x");
      printf("\n\n");
#endif
      mpz_poly_to_F_mpz_poly(F_poly, m_poly1);
      F_mpz_poly_to_mpz_poly(m_poly2, F_poly);

#if DEBUG
      mpz_poly_print_pretty(m_poly2, "x");
      printf("\n\n");
#endif
          
      result = mpz_poly_equal(m_poly1, m_poly2); 
		if (!result) 
			printf("Error: length = %ld, bits = %ld, length1 = %ld, length2 = %ld\n", length, bits, m_poly1->length, m_poly2->length);
          
      F_mpz_poly_clear(F_poly);
   }
   
   mpz_clear(temp);
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

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


