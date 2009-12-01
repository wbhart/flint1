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

zmod_mat-test.c: Test code for zmod_mat.c

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "test-support.h"
#include "zmod_mat.h"
#include "long_extras.h"

#define VARY_BITS 0
#define SPARSE 0

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

/* 
   Generate a random zmod matrix 
   The modulus p and size of the matrix are specified by the parameters
   m->p, m->rows and mat->cols 
*/

void randmat(zmod_mat_t mat)
{
   ulong rows = mat->rows;
   ulong cols = mat->cols;
   ulong p = mat->p;
   ulong * ptr;
              
   long i;
   for (i = 0; i < rows; i++)
   {
	  ptr = mat->arr[i];
	  ulong j;
	  for (j = 0; j < cols; j++)
		  ptr[j] = z_randint(p);
   }
} 

int test_zmod_mat_row_reduce_gauss()
{
   int result = 1;
   zmod_mat_t mat;
   unsigned long bits;
   unsigned long modulus;

   unsigned long count1;
   for (count1 = 0; (count1 < 100) && (result == 1); count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randprime(bits, 0);} while (modulus < 2);
      
	  ulong rows = z_randint(200);
	  ulong cols = z_randint(200);
 
 #if DEBUG
      printf("rows = %ld, cols = %ld, bits = %ld, modulus = %ld\n", rows, cols, bits, modulus);
#endif

	  zmod_mat_init(mat, modulus, rows, cols);

	  randmat(mat);
	  zmod_mat_row_reduce_gauss(mat);
      
	  zmod_mat_clear(mat);
   }

   return result;
}

void zmod_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
#endif
   RUN_TEST(zmod_mat_row_reduce_gauss); 
   
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


