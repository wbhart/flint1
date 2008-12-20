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

NTL-interface-test.cpp: Test functions for conversion between NTL and FLINT format

Copyright (C) 2007, William Hart

*****************************************************************************/

#include <cstdio>
#include <cstring>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <gmp.h>

#include "NTL-interface.h"
#include "fmpz.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"
#include "fmpz_poly.h"

#include "flint.h"
#include "mpz_poly.h"
#include "long_extras.h"
#include "memory-manager.h"
#include "test-support.h"

#define VARY_BITS 1
#define SIGNS 1
#define SPARSE 1

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

unsigned long randint(unsigned long randsup) 
{
    if (randsup == 0) return 0;
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(mpz_poly_t pol, long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_zero(pol);
   
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
          if (randint(10) == 1) mpz_rrandomb(temp, randstate, bits);
          else mpz_set_ui(temp, 0);
#else
          mpz_rrandomb(temp, randstate, bits);
#endif
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

void randpoly_unsigned(mpz_poly_t pol, long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_zero(pol);
   
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
          mpz_rrandomb(temp, randstate, bits);
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   
   mpz_clear(temp);
} 
 
// generate a random mpz_mat_t with the given number of rows and columns and number of bits per entry
void mpz_randmat(mpz_mat_t mat, ulong r, ulong c, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   for (long i = 0; i < r; i++)
   {
		for (long j = 0; j < c; j++)
		{
#if VARY_BITS
         bits = z_randint(maxbits+1);
#else
         bits = maxbits;
#endif
         if (bits == 0) mpz_set_ui(temp, 0);
         else 
         {
#if SPARSE
            if (z_randint(10) == 1) mpz_rrandomb(temp, randstate, bits);
            else mpz_set_ui(temp, 0);
#else
            mpz_rrandomb(temp, randstate, bits);
#endif
#if SIGNS
            if (z_randint(2)) mpz_neg(temp, temp);
#endif
         }
         mpz_set(mat->entries[i*c+j], temp);
		}
   }
   mpz_clear(temp);
} 

void F_mpz_randmat(F_mpz_mat_t mat, ulong r, ulong c, ulong bits)
{
	mpz_mat_t m_mat;
	mpz_mat_init(m_mat, r, c);
	mpz_randmat(m_mat, r, c, bits);
	mpz_mat_to_F_mpz_mat(mat, m_mat);
	mpz_mat_clear(m_mat);
}

int test_ZZ_to_fmpz()
{
   int result = 1;
   unsigned long limbs, limbs2, randlimbs;
   
   fmpz_t int1, int2;
   
   ZZ z;
   
   for (unsigned long i = 0; (i < 1000000) && (result == 1); i++)
   {
      limbs = randint(100) + 1;
      randlimbs = randint(limbs);
      
      int1 = fmpz_init(limbs);
      
      fmpz_random_limbs2(int1, randlimbs);
      if (z_randint(2)) fmpz_neg(int1, int1);

      fmpz_to_ZZ(z, int1);
      limbs2 = ZZ_limbs(z);
      int2 = fmpz_init(limbs2);
      ZZ_to_fmpz(int2, z);
      
      result = fmpz_equal(int1, int2);
      
      fmpz_clear(int1);
      fmpz_clear(int2);
   }

	return result;
}

int test_ZZ_to_F_mpz()
{
   int result = 1;
   unsigned long bits, limbs, limbs2, randbits;
   
   F_mpz_t int1, int2;
   
   ZZ z;
   
   for (unsigned long i = 0; (i < 1000000) && (result == 1); i++)
   {
      bits = randint(1000)+1;
      randbits = randint(bits);
		limbs = (bits - 1)/FLINT_BITS + 1;
      
      F_mpz_init2(int1, limbs);
      
      F_mpz_random(int1, randbits);
      if (z_randint(2)) F_mpz_neg(int1, int1);

      F_mpz_to_ZZ(z, int1);
      limbs2 = ZZ_limbs(z);
      F_mpz_init2(int2, limbs2);
      ZZ_to_F_mpz(int2, z);
      
      result = F_mpz_equal(int1, int2);
      
      F_mpz_clear(int1);
      F_mpz_clear(int2);
   }

	return result;
}

int test_ZZX_to_fmpz_poly()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   for (unsigned long count1 = 1; (count1 < 3000) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          ZZX ZZX_poly;
          length = random_ulong(1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          randpoly(test_poly, length, bits);
          mpz_poly_normalise(test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_poly_to_ZZX(ZZX_poly, test_fmpz_poly);
          ZZX_to_fmpz_poly(test_fmpz_poly2, ZZX_poly);
          
#if DEBUG
          fmpz_poly_print(test_fmpz_poly); printf("\n"); 
          fmpz_poly_print(test_fmpz_poly2); printf("\n"); 
#endif
          
          result = fmpz_poly_equal(test_fmpz_poly, test_fmpz_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   
   return result;
}

int test_mat_ZZ_to_F_mpz_mat()
{
   int result = 1;
   unsigned long bits, r, c;
   
   F_mpz_mat_t mat1, mat2;
   
   for (unsigned long i = 0; (i < 10000) && (result == 1); i++)
   {
      bits = z_randint(200)+1;

		r = z_randint(50);
		c = z_randint(50);
      
      F_mpz_mat_init(mat1, r, c);
      F_mpz_mat_init(mat2, r, c);

		mat_ZZ zmat;
		zmat.SetDims(r, c);
      
      F_mpz_randmat(mat1, r, c, bits);
      
      F_mpz_mat_to_mat_ZZ(zmat, mat1);
		mat_ZZ_to_F_mpz_mat(mat2, zmat);
		
      result = F_mpz_mat_equal(mat1, mat2);
      
      zmat.kill();
		F_mpz_mat_clear(mat1);
      F_mpz_mat_clear(mat2);
   }

	return result;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(ZZ_to_fmpz); 
   RUN_TEST(ZZ_to_F_mpz); 
   RUN_TEST(ZZX_to_fmpz_poly); 
	RUN_TEST(mat_ZZ_to_F_mpz_mat);
   
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



