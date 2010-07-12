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

F_mpz_LLL_helper-test.c: test module for helper functions for F_mpz_LLL code

Copyright (C) 2010, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz_LLL_helper.h"
#include "F_mpz_mat.h"
#include "d_mat.h"
#include "test-support.h"

#define DEBUG 0    // prints debug information
#define DEBUG2 1 

double random_d()
{
   if (z_randint(2)) return rand()/((double) RAND_MAX + 1);
   else return -rand()/((double) RAND_MAX + 1);
}

void random_d_mat(double ** mat, ulong rows, ulong cols)
{
   ulong i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
	     mat[i][j] = random_d();
}

void F_mpz_test_random(F_mpz_t f, ulong bits)
{
	bits = z_randint(bits);
	
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

// generate a random mpz_mat_t with the given number of rows and columns and number of bits per entry
void mpz_randmat(mpz_mat_t mat, ulong r, ulong c, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);
   
   long i;
   for (i = 0; i < r; i++)
   {
		long j;
		for (j = 0; j < c; j++)
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

// same as for mpz_randmat above, except it creates an F_mpz_mat
// WARNING: do not use for testing of conversion between the two formats
void F_mpz_randmat(F_mpz_mat_t mat, ulong r, ulong c, ulong bits)
{
	mpz_mat_t m_mat;
	mpz_mat_init(m_mat, r, c);
	mpz_randmat(m_mat, r, c, bits);
	mpz_mat_to_F_mpz_mat(mat, m_mat);
	mpz_mat_clear(m_mat);
}

/****************************************************************************

   Test code for F_mpz_LLL helper Routines
   
****************************************************************************/

int test_heuristic_scalar_product()
{
   int result = 1;
   ulong count1, count2, i;

   double ** mat;
   F_mpz_mat_t B;
   int expo[50];

   for (count1 = 0; count1 < 100; count1++)
   {
	  ulong rows = z_randint(50) + 1;
      ulong cols = z_randint(50) + 2;

	  F_mpz_mat_init(B, rows, cols);
      F_mpz_randmat(B, rows, cols, 1000);

      mat = d_mat_init(rows, cols);

	  for (i = 0; i < rows; i++)
	     expo[i] = _F_mpz_vec_to_d_vec_2exp(mat[i], B->rows[i], cols);

      for (count2 = 0; count2 < 1000; count2++)
	  {
	     ulong r1 = z_randint(rows);
		 ulong r2 = z_randint(rows);

		 double d1 = heuristic_scalar_product(mat[r1], mat[r2], cols, 
								B, r1, r2, expo[r1] + expo[r2]);
		 double d2 = heuristic_scalar_product(mat[r1], mat[r1], cols, 
								B, r1, r1, expo[r1] + expo[r1]);
		 double d3 = heuristic_scalar_product(mat[r2], mat[r2], cols, 
								B, r2, r2, expo[r2] + expo[r2]);

		 _F_mpz_vec_add(B->rows[r2], B->rows[r1], B->rows[r2], cols);
         _d_vec_add(mat[r2], mat[r1], mat[r2], cols);

		 double d4 = heuristic_scalar_product(mat[r2], mat[r2], cols, 
								B, r2, r2, expo[r2] + expo[r2]);

		 expo[r2] = _F_mpz_vec_to_d_vec_2exp(mat[r2], B->rows[r2], cols);

		 result = (fabs(d4 - d3 - d2 - 2*d1) < 1.0E-12);

		 if (!result)
		 {
		    printf("Failed d1 = %lf, d2 = %lf, d3 = %lf, d4 = %lf\n", d1, d2, d3, d4);
		 }

	  }

      F_mpz_mat_clear(B);
      d_mat_clear(mat);
   }

   return result;
}

/****************************************************************************

   Main test functions

****************************************************************************/

void F_mpz_LLL_helper_test_all()
{
   int success, all_success = 1;

   RUN_TEST(heuristic_scalar_product);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_LLL_helper_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}

// end of file ****************************************************************
