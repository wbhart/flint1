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
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include <math.h>
#include "flint.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "mpfr_mat.h"
#include "test-support.h"

#define DEBUG 0    // prints debug information
#define DEBUG2 1 

void random_mpfr(mpfr_t f)
{
   mpfr_urandomb(f, randstate);
   if (z_randint(2)) mpfr_neg(f, f, GMP_RNDN);
}

void random_mpfr_mat(__mpfr_struct ** mat, ulong rows, ulong cols)
{
   ulong i, j;

   for (i = 0; i < rows; i++)
      for (j = 0; j < cols; j++)
         random_mpfr(mat[i] + j);
}

/****************************************************************************

   Test code for d_mat Routines
   
****************************************************************************/

int test_mpfr_mat_init_clear()
{
   int result = 1;
   ulong count1;

   __mpfr_struct ** mat;

   for (count1 = 0; count1 < 10000; count1++)
   {
	  ulong rows = z_randint(50) + 1;
      ulong cols = z_randint(50) + 1;

      mat = mpfr_mat_init(rows, cols);
      random_mpfr_mat(mat, rows, cols);
      mpfr_mat_clear(mat, rows, cols);
   }

   return result;
}

int test__mpfr_vec_scalar_product()
{
   int result = 1;
   ulong count1, count2;
   mpfr_t s1, s2, s3;

   __mpfr_struct ** mat;

   mpfr_init(s1);
   mpfr_init(s2);
   mpfr_init(s3);

   for (count1 = 0; count1 < 200; count1++)
   {
	  ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      mat = mpfr_mat_init(rows, cols);
      random_mpfr_mat(mat, rows, cols);
     
	  for (count2 = 0; count2 < 200; count2++)
	  {
	     ulong r1 = z_randint(rows);
		 ulong r2 = z_randint(rows);
		 _mpfr_vec_scalar_product(s1, mat[r1], mat[r2], cols - 1);
		 _mpfr_vec_scalar_product(s2, mat[r1] + cols - 1, mat[r2] + cols - 1, 1);
		 _mpfr_vec_scalar_product(s3, mat[r1], mat[r2], cols);
       mpfr_add(s1, s1, s2, GMP_RNDN);

		 result = (mpfr_cmp(s1, s3) == 0);
		 if (!result)
		 {
			 printf("Error: %ld, %ld, %ld, %ld\n", rows, cols, r1, r2);
			 break;
		 }
	  }

      mpfr_mat_clear(mat, rows, cols);
   }

   mpfr_clear(s1);
   mpfr_clear(s2);
   mpfr_clear(s3);

   return result;
}

int test__mpfr_vec_scalar_product2()
{
   int result = 1;
   ulong count1, count2, prec;
   mpfr_t s1, s2, s3;

   __mpfr_struct ** mat;

   for (count1 = 0; count1 < 200; count1++)
   {
	  ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      mat = mpfr_mat_init(rows, cols);
      random_mpfr_mat(mat, rows, cols);
      prec = z_randint(200) + MPFR_PREC_MIN;

      mpfr_init2(s1, prec);
      mpfr_init2(s2, prec);
      mpfr_init2(s3, prec);

	  for (count2 = 0; count2 < 200; count2++)
	  {
	     ulong r1 = z_randint(rows);
		 ulong r2 = z_randint(rows);
		 _mpfr_vec_scalar_product2(s1, mat[r1], mat[r2], cols - 1, prec);
		 _mpfr_vec_scalar_product2(s2, mat[r1] + cols - 1, mat[r2] + cols - 1, 1, prec);
		 _mpfr_vec_scalar_product2(s3, mat[r1], mat[r2], cols, prec);
       mpfr_add(s1, s1, s2, GMP_RNDN);

		 result = (mpfr_cmp(s1, s3) == 0);
		 if (!result)
		 {
			 printf("Error: %ld, %ld, %ld, %ld\n", rows, cols, r1, r2);
			 break;
		 }
	  }

      mpfr_clear(s1);
      mpfr_clear(s2);
      mpfr_clear(s3);

      mpfr_mat_clear(mat, rows, cols);
   }

   return result;
}

int test__mpfr_vec_norm()
{
   int result = 1;
   ulong count1, count2;
   mpfr_t s1, s2, s3;

   __mpfr_struct ** mat;

   mpfr_init(s1);
   mpfr_init(s2);
   mpfr_init(s3);

   for (count1 = 0; count1 < 200; count1++)
   {
	   ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      mat = mpfr_mat_init(rows, cols);
      random_mpfr_mat(mat, rows, cols);
     
	   for (count2 = 0; count2 < 200; count2++)
	   {
	      ulong r1 = z_randint(rows);
		   _mpfr_vec_norm(s1, mat[r1], cols - 1);
		   _mpfr_vec_norm(s2, mat[r1] + cols - 1, 1);
  	      _mpfr_vec_norm(s3, mat[r1], cols);
         mpfr_add(s1, s1, s2, GMP_RNDN);

		   result = (mpfr_cmp(s1, s3) == 0);
		   if (!result)
		   {
			   printf("Error: %ld, %ld, %ld\n", rows, cols, r1);
			   break;
		   }
	   }

      mpfr_mat_clear(mat, rows, cols);
   }

   mpfr_clear(s1);
   mpfr_clear(s2);
   mpfr_clear(s3);

   return result;
}

int test__mpfr_vec_norm2()
{
   int result = 1;
   ulong count1, count2, prec;
   mpfr_t s1, s2, s3;

   __mpfr_struct ** mat;

   for (count1 = 0; count1 < 200; count1++)
   {
	   ulong rows = z_randint(100) + 1;
      ulong cols = z_randint(100) + 2;

      prec = z_randint(200) + MPFR_PREC_MIN;

      mat = mpfr_mat_init(rows, cols);
      random_mpfr_mat(mat, rows, cols);
     
      mpfr_init2(s1, prec);
      mpfr_init2(s2, prec);
      mpfr_init2(s3, prec);

	   for (count2 = 0; count2 < 200; count2++)
	   {
	      ulong r1 = z_randint(rows);
		   _mpfr_vec_norm2(s1, mat[r1], cols - 1, prec);
		   _mpfr_vec_norm2(s2, mat[r1] + cols - 1, 1, prec);
  	      _mpfr_vec_norm2(s3, mat[r1], cols, prec);
         mpfr_add(s1, s1, s2, GMP_RNDN);

		   result = (mpfr_cmp(s1, s3) == 0);
		   if (!result)
		   {
			   printf("Error: %ld, %ld, %ld\n", rows, cols, r1);
			   break;
		   }
	   }

      mpfr_clear(s1);
      mpfr_clear(s2);
      mpfr_clear(s3);

      mpfr_mat_clear(mat, rows, cols);
   }

   return result;
}

/****************************************************************************

   Main test functions

****************************************************************************/

void mpfr_mat_test_all()
{
   int success, all_success = 1;

   RUN_TEST(mpfr_mat_init_clear);
   RUN_TEST(_mpfr_vec_scalar_product);
   RUN_TEST(_mpfr_vec_scalar_product2);
   RUN_TEST(_mpfr_vec_norm);
   RUN_TEST(_mpfr_vec_norm2);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   mpfr_mat_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}

// end of file ****************************************************************
