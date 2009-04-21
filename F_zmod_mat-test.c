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
#include "F_mpzmod_mat.h"
#include "packed_vec.h"

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

int test_F_zmod_mat_convert()
{
   int result = 1;
   F_zmod_mat_t mat1, mat2;
	F_mpzmod_mat_t mmat;
   unsigned long bits;
   unsigned long modulus;
	F_mpz mod[1];
	F_mpz_init(mod);

   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1); count1++)
   {
      bits = z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      
	   F_mpz_set_ui(mod, modulus);
		
		ulong rows = z_randint(200);
	   ulong cols = z_randint(200);
 
#if DEBUG
      printf("rows = %ld, cols = %ld, bits = %ld, modulus = %ld\n", rows, cols, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, rows, cols);
      F_zmod_mat_init(mat2, modulus, rows, cols);
      F_mpzmod_mat_init(mmat, mod, rows, cols);

	   randmat(mat1);
	   
		F_zmod_mat_to_F_mpzmod_mat(mmat, mat1);
      F_mpzmod_mat_to_F_zmod_mat(mat2, mmat);

      ulong i, j, m1, m2;
		for (i = 0; (i < mat1->r) && (result == 1); i++)
		{
			for (j = 0; (j < mat1->c) && (result == 1); j++)
			{
				PV_GET_ENTRY(m1, mat1->arr, mat1->rows[i] + j);
            PV_GET_ENTRY(m2, mat2->arr, mat2->rows[i] + j);
            result = (m1 == m2);
			}
		}

		if (!result)
		{
			printf("i = %ld, j = %ld, m1 = %ld, m2 = %ld, modulus = %ld\n", i, j, m1, m2, modulus);
		}
      
		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(mat2);
		F_mpzmod_mat_clear(mmat);
   }

	F_mpz_clear(mod);

   return result;
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
   F_zmod_mat_t mat1, mat2, res, res2;
	F_mpzmod_mat_t mmat1, mmat2, mres;
   unsigned long bits;
   unsigned long modulus;
	F_mpz mod[1];
	F_mpz_init(mod);

   for (unsigned long count1 = 0; (count1 < 1) && (result == 1); count1++)
   {
      bits = 8;//z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      modulus = 251;

		F_mpz_set_ui(mod, modulus);
		
		ulong dim = 256;
		ulong r1 = dim;//z_randint(100);
	   ulong c1 = dim;//z_randint(100) + 1;
		ulong r2 = c1;
		ulong c2 = dim;//z_randint(100);
 
#if DEBUG
      printf("r1 = %ld, c1 = %ld, r2 = %ld, c2 = %ld, bits = %ld, modulus = %ld\n", r1, c1, r2, c2, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, r1, c1);
      F_zmod_mat_init(mat2, modulus, r2, c2);
      F_zmod_mat_init(res, modulus, r1, c2);
      F_zmod_mat_init(res2, modulus, r1, c2);

	   /*F_mpzmod_mat_init(mmat1, mod, r1, c1);
      F_mpzmod_mat_init(mmat2, mod, r2, c2);
      F_mpzmod_mat_init(mres, mod, r1, c2);*/

	   randmat(mat1);
	   randmat(mat2);

		for (ulong i = 0; i < 100; i++) F_zmod_mat_mul_classical(res, mat1, mat2);

		/*F_zmod_mat_to_F_mpzmod_mat(mmat1, mat1);
      F_zmod_mat_to_F_mpzmod_mat(mmat2, mat2);

		F_mpzmod_mat_mul_classical(mres, mmat1, mmat2);
      F_mpzmod_mat_to_F_zmod_mat(res2, mres);

		ulong i, j, m1, m2;
		for (i = 0; (i < res->r) && (result == 1); i++)
		{
			pv_iter_s i1, i2;
			PV_ITER_INIT(i1, res->arr, res->rows[i]);
			PV_ITER_INIT(i2, res2->arr, res2->rows[i]);

			for (j = 0; (j < res->c) && (result == 1); j++)
			{
            PV_GET_NEXT(m1, i1);
            PV_GET_NEXT(m2, i2);
				result &= (m1 == m2);
			}
		}

		if (!result) 
		{
			printf("Error: bits = %ld, i = %ld, j = %ld, rows = %ld, cols = %ld, modulus = %ld, m1 = %ld, m2 = %ld\n", bits, i - 1, j - 1, r1, c2, modulus, m1, m2);
		}

		F_mpzmod_mat_clear(mmat1);
 		F_mpzmod_mat_clear(mmat2);
 		F_mpzmod_mat_clear(mres);*/

		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(mat2);
 		F_zmod_mat_clear(res);
  		F_zmod_mat_clear(res2);
   }

	F_mpz_clear(mod);

   return result;
}

int test_F_zmod_mat_mul_strassen()
{
   int result = 1;
   F_zmod_mat_t mat1, mat2, res, res2;
	F_mpzmod_mat_t mmat1, mmat2, mres;
   unsigned long bits;
   unsigned long modulus;
	F_mpz mod[1];
	F_mpz_init(mod);

   for (unsigned long count1 = 0; (count1 < 1) && (result == 1); count1++)
   {
      bits = 8;//z_randint(FLINT_BITS - 2) + 2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);
      modulus = 251;

		F_mpz_set_ui(mod, modulus);
		
		ulong r1 = 1024;
	   ulong c1 = 1024;
		ulong r2 = c1;
		ulong c2 = 1024;
 
#if DEBUG
      printf("r1 = %ld, c1 = %ld, r2 = %ld, c2 = %ld, bits = %ld, modulus = %ld\n", r1, c1, r2, c2, bits, modulus);
#endif

	   F_zmod_mat_init(mat1, modulus, r1, c1);
      F_zmod_mat_init(mat2, modulus, r2, c2);
      F_zmod_mat_init(res, modulus, r1, c2);
      F_zmod_mat_init(res2, modulus, r1, c2);

	   /*F_mpzmod_mat_init(mmat1, mod, r1, c1);
      F_mpzmod_mat_init(mmat2, mod, r2, c2);
      F_mpzmod_mat_init(mres, mod, r1, c2);*/

	   randmat(mat1);
	   randmat(mat2);

		for (ulong i = 0; i < 100; i++) F_zmod_mat_mul_strassen(res, mat1, mat2);

		/*F_zmod_mat_to_F_mpzmod_mat(mmat1, mat1);
      F_zmod_mat_to_F_mpzmod_mat(mmat2, mat2);

		F_mpzmod_mat_mul_classical(mres, mmat1, mmat2);
      F_mpzmod_mat_to_F_zmod_mat(res2, mres);

		ulong i, j, m1, m2;
		for (i = 0; (i < res->r) && (result == 1); i++)
		{
			pv_iter_s i1, i2;
			PV_ITER_INIT(i1, res->arr, res->rows[i]);
			PV_ITER_INIT(i2, res2->arr, res2->rows[i]);

			for (j = 0; (j < res->c) && (result == 1); j++)
			{
            PV_GET_NEXT(m1, i1);
            PV_GET_NEXT(m2, i2);
				result &= (m1 == m2);
			}
		}

		if (!result) 
		{
			printf("Error: bits = %ld, i = %ld, j = %ld, rows = %ld, cols = %ld, modulus = %ld, m1 = %ld, m2 = %ld\n", bits, i - 1, j - 1, r1, c2, modulus, m1, m2);
		}

		F_mpzmod_mat_clear(mmat1);
 		F_mpzmod_mat_clear(mmat2);
 		F_mpzmod_mat_clear(mres);*/

		F_zmod_mat_clear(mat1);
 		F_zmod_mat_clear(mat2);
 		F_zmod_mat_clear(res);
  		F_zmod_mat_clear(res2);
   }

	F_mpz_clear(mod);

   return result;
}

void zmod_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
#endif
   /*RUN_TEST(F_zmod_mat_convert); 
   RUN_TEST(F_zmod_mat_add); 
   RUN_TEST(F_zmod_mat_sub); 
   RUN_TEST(F_zmod_mat_neg); */
   RUN_TEST(F_zmod_mat_mul_classical); 
   //RUN_TEST(F_zmod_mat_mul_strassen); 
   
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


