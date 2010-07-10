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

d_mat-test.c: test module for matrices with entries that are doubles

Copyright (C) 2010, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "memory-manager.h"
#include "long_extras.h"
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

/****************************************************************************

   Test code for d_mat Routines
   
****************************************************************************/

int test_d_mat_init_clear()
{
   int result = 1;
   ulong count1;

   double ** mat;

   for (count1 = 0; count1 < 10000; count1++)
   {
	  ulong rows = z_randint(50) + 1;
      ulong cols = z_randint(50) + 1;

      mat = d_mat_init(rows, cols);
      random_d_mat(mat, rows, cols);
      d_mat_clear(mat);
   }

   return result;
}

int test_d_mat_row_add_sub()
{
   int result = 1;
   ulong count1, count2;

   double ** mat;
   double * vec;

   for (count1 = 0; count1 < 1000; count1++)
   {
	  ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      mat = d_mat_init(rows, cols);
      random_d_mat(mat, rows, cols);
      vec = malloc(cols*sizeof(double));

	  for (count2 = 0; count2 < 1000; count2++)
	  {
	     ulong r1 = z_randint(rows);
		 ulong r2 = z_randint(rows);

		 d_mat_row_add(vec, mat[r1], mat[r2], 0, cols);
		 d_mat_row_sub(vec, vec, mat[r2], 0, cols);
		 
		 result = (d_mat_row_equal(mat[r1], vec, 0, cols, FLINT_D_BITS - 1));
		 if (!result)
		 {
			 printf("Error: a + b - b != a, rows = %ld, cols = %ld\n", rows, cols);
			 break;
		 }
	  }

      d_mat_clear(mat);
	  free(vec);
   }

   return result;
}

int test_d_vec_scalar_product()
{
   int result = 1;
   ulong count1, count2;

   double ** mat;

   for (count1 = 0; count1 < 1000; count1++)
   {
	  ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      mat = d_mat_init(rows, cols);
      random_d_mat(mat, rows, cols);
     
	  for (count2 = 0; count2 < 1000; count2++)
	  {
	     ulong r1 = z_randint(rows);
		 ulong r2 = z_randint(rows);
		 double s1 = d_vec_scalar_product(mat[r1], mat[r2], cols - 1);
		 double s2 = d_vec_scalar_product(mat[r1] + cols - 1, mat[r2] + cols - 1, 1);
		 double s3 = d_vec_scalar_product(mat[r1], mat[r2], cols);

		 result = ((s1 + s2) == s3);
		 if (!result)
		 {
			 printf("Error: %ld, %ld, %ld\n", s1, s2, s3);
			 break;
		 }
	  }

      d_mat_clear(mat);
   }

   return result;
}

int test_d_vec_norm()
{
   int result = 1;
   ulong count1, count2;

   double ** mat;

   for (count1 = 0; count1 < 1000; count1++)
   {
	  ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      mat = d_mat_init(rows, cols);
      random_d_mat(mat, rows, cols);
     
	  for (count2 = 0; count2 < 1000; count2++)
	  {
	     ulong r1 = z_randint(rows);
		 double s1 = d_vec_norm(mat[r1], cols - 1);
		 double s2 = d_vec_norm(mat[r1] + cols - 1, 1);
		 double s3 = d_vec_norm(mat[r1], cols);

		 result = ((s1 + s2) == s3);
		 if (!result)
		 {
			 printf("Error: %ld, %ld, %ld\n", s1, s2, s3);
			 break;
		 }
	  }

      d_mat_clear(mat);
   }

   return result;
}

/****************************************************************************

   Main test functions

****************************************************************************/

void d_mat_test_all()
{
   int success, all_success = 1;

   RUN_TEST(d_mat_init_clear);
   RUN_TEST(d_mat_row_add_sub);
   RUN_TEST(d_vec_scalar_product);
   RUN_TEST(d_vec_norm);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   d_mat_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}

// end of file ****************************************************************
