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
#include "mpz_mat.h"
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

// generate a dense random mpz_mat_t with up to the given length and number of bits per entry
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
            mpz_rrandomb(temp, randstate, bits);
#if SIGNS
            if (z_randint(2)) mpz_neg(temp, temp);
#endif
         }
         mpz_set(mat->entries[i*c+j], temp);
		}
   }
   mpz_clear(temp);
} 

void mpz_mat_to_zmod_mat(zmod_mat_t out, mpz_mat_t mat)
{
   mpz_t temp;
   mpz_init(temp);
   
   for (ulong i = 0; i < out->rows; i++)
      for (ulong j = 0; j < out->cols; j++)
         out->arr[i][j] = mpz_mod_ui(temp, mat->entries[i*out->cols + j], out->p);

   mpz_clear(temp);
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

int test_zmod_mat_mul_classical()
{
   int result = 1;

   ulong r1, rc, c2, modulus, bits, bits2;
   mpz_mat_t m_mat1, m_mat2, m_mat3;
   zmod_mat_t z_mat1, z_mat2, z_mat3, z_mat4;

   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r1 = z_randint(50) + 1;
	   rc = z_randint(50) + 1;
      c2 = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r1, rc);
      mpz_mat_init(m_mat2, rc, c2);
      mpz_mat_init(m_mat3, r1, c2);

      zmod_mat_init(z_mat1, modulus, r1, rc);
      zmod_mat_init(z_mat2, modulus, rc, c2);
      zmod_mat_init(z_mat3, modulus, r1, c2);
      zmod_mat_init(z_mat4, modulus, r1, c2);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r1, rc, bits2);
      mpz_randmat(m_mat2, rc, c2, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_mul_classical(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_mul_classical(z_mat3, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat3, z_mat4));

      if (!result)
      {
         printf("Error: r1 = %ld, rc = %ld, c2 = %ld, modulus = %lu\n", r1, rc, c2, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }

   return result;
}

int test_zmod_mat_mul_strassen()
{
   int result = 1;

   ulong r1, rc, c2, modulus, bits, bits2;
   mpz_mat_t m_mat1, m_mat2, m_mat3;
   zmod_mat_t z_mat1, z_mat2, z_mat3, z_mat4;

   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r1 = z_randint(50) + 1;
	   rc = z_randint(50) + 1;
      c2 = z_randint(50) + 1;

      mpz_mat_init(m_mat1, r1, rc);
      mpz_mat_init(m_mat2, rc, c2);
      mpz_mat_init(m_mat3, r1, c2);

      zmod_mat_init(z_mat1, modulus, r1, rc);
      zmod_mat_init(z_mat2, modulus, rc, c2);
      zmod_mat_init(z_mat3, modulus, r1, c2);
      zmod_mat_init(z_mat4, modulus, r1, c2);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r1, rc, bits2);
      mpz_randmat(m_mat2, rc, c2, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_mul_classical(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_mul_strassen(z_mat3, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat3, z_mat4));

      if (!result)
      {
         printf("Error: r1 = %ld, rc = %ld, c2 = %ld, modulus = %lu\n", r1, rc, c2, modulus);
         zmod_mat_print(z_mat3); printf("\n\n");
         zmod_mat_print(z_mat4); printf("\n\n");
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }

   return result;
}

int test_zmod_mat_add()
{
   int result = 1;

   ulong r, c, modulus, bits, bits2;
   mpz_mat_t m_mat1, m_mat2, m_mat3;
   zmod_mat_t z_mat1, z_mat2, z_mat3, z_mat4;

   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat2, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat2, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      mpz_randmat(m_mat2, r, c, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_add(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_add(z_mat3, z_mat1, z_mat2);
      
      result = (zmod_mat_equal(z_mat3, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }
   
   // alias 1 and 2
   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat2, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat2, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      mpz_randmat(m_mat2, r, c, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_add(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_add(z_mat1, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat1, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }
   
   // alias 1 and 3
   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat2, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat2, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      mpz_randmat(m_mat2, r, c, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_add(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_add(z_mat2, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat2, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }
   
   // alias 2 and 3
   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      
      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      
      mpz_mat_add(m_mat3, m_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_add(z_mat3, z_mat1, z_mat1);

      result = (zmod_mat_equal(z_mat3, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat3);
   }

   return result;
}

int test_zmod_mat_sub()
{
   int result = 1;

   ulong r, c, modulus, bits, bits2;
   mpz_mat_t m_mat1, m_mat2, m_mat3;
   zmod_mat_t z_mat1, z_mat2, z_mat3, z_mat4;

   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat2, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat2, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      mpz_randmat(m_mat2, r, c, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_sub(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_sub(z_mat3, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat3, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }

   // alias 1 and 2
   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat2, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat2, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      mpz_randmat(m_mat2, r, c, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_sub(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_sub(z_mat1, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat1, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }

   // alias 1 and 3
   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat2, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat2, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      mpz_randmat(m_mat2, r, c, bits2);

      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat2, m_mat2);

      mpz_mat_sub(m_mat3, m_mat1, m_mat2);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_sub(z_mat2, z_mat1, z_mat2);

      result = (zmod_mat_equal(z_mat2, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat2);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);
   }

   // alias 2 and 3
   for (ulong count1 = 0; count1 < 1000; count1++)
   {
      bits = z_randint(FLINT_BITS-2)+2;
      
      do {modulus = z_randbits(bits);} while (modulus < 2);

	   r = z_randint(50) + 1;
      c = z_randint(50) + 1;
      
      mpz_mat_init(m_mat1, r, c);
      mpz_mat_init(m_mat3, r, c);

      zmod_mat_init(z_mat1, modulus, r, c);
      zmod_mat_init(z_mat3, modulus, r, c);
      zmod_mat_init(z_mat4, modulus, r, c);
      
      bits2 = z_randint(100) + 1;
      mpz_randmat(m_mat1, r, c, bits2);
      
      mpz_mat_to_zmod_mat(z_mat1, m_mat1);
      
      mpz_mat_sub(m_mat3, m_mat1, m_mat1);
      mpz_mat_to_zmod_mat(z_mat4, m_mat3);
      zmod_mat_sub(z_mat3, z_mat1, z_mat1);

      result = (zmod_mat_equal(z_mat3, z_mat4));

      if (!result)
      {
         printf("Error: r = %ld, c = %ld, modulus = %lu\n", r, c, modulus);
         abort();
      }

      zmod_mat_clear(z_mat1);
      zmod_mat_clear(z_mat3);
      zmod_mat_clear(z_mat4);

      mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat3);
   }

   return result;
}

void zmod_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
#endif
   RUN_TEST(zmod_mat_add); 
   RUN_TEST(zmod_mat_sub); 
   RUN_TEST(zmod_mat_row_reduce_gauss); 
   RUN_TEST(zmod_mat_mul_classical); 
   RUN_TEST(zmod_mat_mul_strassen); 
   
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


