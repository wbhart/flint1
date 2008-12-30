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

F_zmod_mat-test.c: Test code for F_zmod_mat.c

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "test-support.h"
#include "F_zmod_mat.h"
#include "long_extras.h"

#define VARY_BITS 0
#define SPARSE 0

#define TESTFILE 0 // Set this to test matrix reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

/* 
   Generate a random F_zmod matrix 
   The modulus p and size of the matrix are specified by the parameters
   m->p, m->rows and mat->cols 
*/

void randmat(F_zmod_mat_t mat)
{
   ulong rows = mat->r;
   ulong cols = mat->c;
   ulong p = mat->p;
   ulong * ptr;
              
   pv_iter_s iter;
	
	for (long i = 0; i < rows; i++)
   {
	  PV_ITER_INIT(iter, mat->arr, mat->rows[i]);
	  for (ulong j = 0; j < cols; j++)
		  PV_SET_NEXT(iter, z_randint(p));
   }
} 

int test_F_zmod_mat_add()
{
   int result = 1;
   F_zmod_mat_t mat1, mat2, res;
   unsigned long bits;
   unsigned long modulus;

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      
	   ulong rows = z_randint(200);
	   ulong cols = z_randint(200);
 
#if DEBUG
      printf("rows = %ld, cols = %ld, bits = %ld, modulus = %ld\n", rows, cols, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, rows, cols);
      F_zmod_mat_init(mat2, modulus, rows, cols);
      F_zmod_mat_init(res, modulus, rows, cols);

	   randmat(mat1);
	   randmat(mat2);

		F_zmod_mat_add(res, mat1, mat2);

      ulong i, j, m1, m2, m3;
		for (i = 0; (i < mat1->r) && (result == 1); i++)
		{
			for (j = 0; (j < mat1->c) && (result == 1); j++)
			{
				PV_GET_ENTRY(m1, mat1->arr, mat1->rows[i] + j);
            PV_GET_ENTRY(m2, mat2->arr, mat2->rows[i] + j);
            PV_GET_ENTRY(m3, res->arr, res->rows[i] + j);
            result = (m3 == (m1 + m2) % modulus);
			}
		}

		if (!result)
		{
			printf("i = %ld, j = %ld, m1 = %ld, m2 = %ld, m3 = %ld, modulus = %ld\n", i, j, m1, m2, m3, modulus);
		}
      
		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(mat2);
 		F_zmod_mat_clear(res);
   }

   return result;
}

int test_F_zmod_mat_sub()
{
   int result = 1;
   F_zmod_mat_t mat1, mat2, res;
   unsigned long bits;
   unsigned long modulus;

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      
	   ulong rows = z_randint(200);
	   ulong cols = z_randint(200);
 
#if DEBUG
      printf("rows = %ld, cols = %ld, bits = %ld, modulus = %ld\n", rows, cols, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, rows, cols);
      F_zmod_mat_init(mat2, modulus, rows, cols);
      F_zmod_mat_init(res, modulus, rows, cols);

	   randmat(mat1);
	   randmat(mat2);

		F_zmod_mat_sub(res, mat1, mat2);

      ulong i, j, m1, m2, m3;
		for (i = 0; (i < mat1->r) && (result == 1); i++)
		{
			for (j = 0; (j < mat1->c) && (result == 1); j++)
			{
				PV_GET_ENTRY(m1, mat1->arr, mat1->rows[i] + j);
            PV_GET_ENTRY(m2, mat2->arr, mat2->rows[i] + j);
            PV_GET_ENTRY(m3, res->arr, res->rows[i] + j);
            if (m1 >= m2) result = (m3 == (m1 - m2) % modulus);
				else result = (m3 == (modulus - m2 + m1) % modulus);
			}
		}

		if (!result)
		{
			printf("i = %ld, j = %ld, m1 = %ld, m2 = %ld, m3 = %ld, modulus = %ld\n", i, j, m1, m2, m3, modulus);
		}
      
		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(mat2);
 		F_zmod_mat_clear(res);
   }

   return result;
}

int test_F_zmod_mat_neg()
{
   int result = 1;
   F_zmod_mat_t mat1, res;
   unsigned long bits;
   unsigned long modulus;

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      
	   ulong rows = z_randint(200);
	   ulong cols = z_randint(200);
 
#if DEBUG
      printf("rows = %ld, cols = %ld, bits = %ld, modulus = %ld\n", rows, cols, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, rows, cols);
      F_zmod_mat_init(res, modulus, rows, cols);

	   randmat(mat1);
	   
		F_zmod_mat_neg(res, mat1);

      ulong i, j, m1, m2;
		for (i = 0; (i < mat1->r) && (result == 1); i++)
		{
			for (j = 0; (j < mat1->c) && (result == 1); j++)
			{
				PV_GET_ENTRY(m1, mat1->arr, mat1->rows[i] + j);
            PV_GET_ENTRY(m2, res->arr, res->rows[i] + j);
            result = (m2 == (modulus - m1) % modulus);
			}
		}

		if (!result)
		{
			printf("i = %ld, j = %ld, m1 = %ld, m2 = %ld, modulus = %ld\n", i, j, m1, m2, modulus);
		}
      
		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(res);
   }

   return result;
}

int test_F_zmod_mat_mul_classical()
{
   int result = 1;
   F_zmod_mat_t mat1, mat2, res;
   unsigned long bits;
   unsigned long modulus;

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      
	   ulong r1 = z_randint(200);
	   ulong c1 = z_randint(200) + 1;
		ulong r2 = c1;
		ulong c2 = z_randint(200);
 
#if DEBUG
      printf("r1 = %ld, c1 = %ld, r2 = %ld, c2 = %ld, bits = %ld, modulus = %ld\n", r1, c1, r2, c2, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, r1, c1);
      F_zmod_mat_init(mat2, modulus, r2, c2);
      F_zmod_mat_init(res, modulus, r1, c2);

	   randmat(mat1);
	   randmat(mat2);

		F_zmod_mat_mul_classical(res, mat1, mat2);

		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(mat2);
 		F_zmod_mat_clear(res);
   }

   return result;
}

void zmod_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
#endif
   RUN_TEST(F_zmod_mat_add); 
   RUN_TEST(F_zmod_mat_sub); 
   RUN_TEST(F_zmod_mat_neg); 
   RUN_TEST(F_zmod_mat_mul_classical); 
   
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


