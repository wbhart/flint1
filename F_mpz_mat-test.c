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

F_mpz_mat-test.c: Test code for F_mpz_mat.c and F_mpz_mat.h

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "flint.h"
#include "long_extras.h"
#include "mpz_mat.h"
#include "F_mpz_mat.h"
#include "memory-manager.h"
#include "test-support.h"

#define VARY_BITS 1 // random entries have random number of bits up to the limit given
#define SIGNS 1 // random entries will be randomly signed
#define SPARSE 1 // matrices are sparse (triggers more corner cases)
#define ITER 1 // if you want all tests to run longer, increase this

#define TESTFILE 0 // Set this to test matrix reading and writing to a file in the current dir

#define DEBUG 0 // allows easy switching of debugging code on and off when debugging (if inserted)
#define DEBUG2 1 

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

// generate a dense random mpz_mat_t with up to the given length and number of bits per entry
void mpz_randmat_dense(mpz_mat_t mat, ulong r, ulong c, ulong maxbits)
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

int test_F_mpz_mat_convert()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat;
   int result = 1;
   ulong bits, r, c;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      r = z_randint(30);
      c = z_randint(30);
      
		F_mpz_mat_init(F_mat, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 

      bits = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits);
           
      mpz_mat_to_F_mpz_mat(F_mat, m_mat1);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);
          
      result = mpz_mat_equal(m_mat1, m_mat2); 
		if (!result) 
		{
			printf("Error: bits = %ld, r = %ld, c = %ld\n", bits, r, c);
		}
          
      F_mpz_mat_clear(F_mat);
		mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
   }
   
   return result;
}

int test_F_mpz_mat_add()
{
   mpz_mat_t m_mat1, m_mat2, res1, res2;
   F_mpz_mat_t F_mat1, F_mat2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   for (ulong count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30);
      ulong c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(res, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
      mpz_mat_init(res1, r, c); 
      mpz_mat_init(res2, r, c); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      
		F_mpz_mat_add(res, F_mat1, F_mat2);
		F_mpz_mat_to_mpz_mat(res2, res);

      mpz_mat_add(res1, m_mat1, m_mat2);		
		    
      result = mpz_mat_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
		F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(res);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2); 
   }
   
	// test aliasing of res and mat1
   for (ulong count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30);
      ulong c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
      mpz_mat_init(res1, r, c); 
      mpz_mat_init(res2, r, c); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      
		F_mpz_mat_add(F_mat1, F_mat1, F_mat2);
		F_mpz_mat_to_mpz_mat(res2, F_mat1);

      mpz_mat_add(res1, m_mat1, m_mat2);		
		    
      result = mpz_mat_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
		F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2);  
   }
   
   // test aliasing of res and mat2
   for (ulong count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30);
      ulong c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
      mpz_mat_init(res1, r, c); 
      mpz_mat_init(res2, r, c); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      
		F_mpz_mat_add(F_mat2, F_mat1, F_mat2);
		F_mpz_mat_to_mpz_mat(res2, F_mat2);

      mpz_mat_add(res1, m_mat1, m_mat2);		
		    
      result = mpz_mat_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
		F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2);  
   }
   
   return result;
}

int test_F_mpz_mat_sub()
{
   mpz_mat_t m_mat1, m_mat2, res1, res2;
   F_mpz_mat_t F_mat1, F_mat2, res;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   for (ulong count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30);
      ulong c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(res, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
      mpz_mat_init(res1, r, c); 
      mpz_mat_init(res2, r, c); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      
		F_mpz_mat_sub(res, F_mat1, F_mat2);
		F_mpz_mat_to_mpz_mat(res2, res);

      mpz_mat_sub(res1, m_mat1, m_mat2);		
		    
      result = mpz_mat_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
		F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(res);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2); 
   }
   
	// test aliasing of res and mat1
   for (ulong count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30);
      ulong c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
      mpz_mat_init(res1, r, c); 
      mpz_mat_init(res2, r, c); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      
		F_mpz_mat_sub(F_mat1, F_mat1, F_mat2);
		F_mpz_mat_to_mpz_mat(res2, F_mat1);

      mpz_mat_sub(res1, m_mat1, m_mat2);		
		    
      result = mpz_mat_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
		F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2);  
   }
   
   // test aliasing of res and mat2
   for (ulong count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30);
      ulong c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
      mpz_mat_init(res1, r, c); 
      mpz_mat_init(res2, r, c); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      
		mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      
		F_mpz_mat_sub(F_mat2, F_mat1, F_mat2);
		F_mpz_mat_to_mpz_mat(res2, F_mat2);

      mpz_mat_sub(res1, m_mat1, m_mat2);		
		    
      result = mpz_mat_equal(res1, res2); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld, length2 = %ld, bits2 = %ld\n", length1, bits1, length2, bits2);
		}
          
		F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2); 
   }
   
   return result;
}

int test_F_mpz_mat_set()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat1, F_mat2;
   int result = 1;
   ulong bits, r, c;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30);
      c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
   
      mpz_randmat(m_mat1, r, c, bits);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      F_mpz_mat_set(F_mat2, F_mat1);
		F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);
          
      result = mpz_mat_equal(m_mat1, m_mat2); 
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld\n", r, c, bits);
		}
          
      F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
		mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
   }
   
	// aliasing is trivial for set and doesn't need testing

   return result;
}

int test_F_mpz_mat_equal()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat1, F_mat2;
   int result = 1;
   ulong bits1, bits2, length1, length2;
   
   // random mats unlikely to be equal, test against mpz_mat_equal
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
   {
      ulong r = z_randint(30);
		ulong c = z_randint(30);
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      mpz_mat_init(m_mat2, r, c); 
   
		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      mpz_randmat(m_mat1, r, c, bits1);
      mpz_randmat(m_mat2, r, c, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
      result = (mpz_mat_equal(m_mat1, m_mat2) == F_mpz_mat_equal(F_mat1, F_mat2)); 
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits1 = %ld, bits2 = %ld\n", r, c, bits1, bits2);
		}
          
      F_mpz_mat_clear(F_mat1);
		F_mpz_mat_clear(F_mat2);
		mpz_mat_clear(m_mat1);
      mpz_mat_clear(m_mat2);
   }

	// mats are equal
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
   {
      ulong r = z_randint(30);
		ulong c = z_randint(30);
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      
		bits1 = z_randint(200) + 1;
      mpz_randmat(m_mat1, r, c, bits1);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      F_mpz_mat_set(F_mat2, F_mat1);
          
      result = (F_mpz_mat_equal(F_mat1, F_mat2)); 
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld\n", r, c, bits1);
		}
          
      F_mpz_mat_clear(F_mat1);
		F_mpz_mat_clear(F_mat2);
		mpz_mat_clear(m_mat1);
   }

	// mats are aliased
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
   {
      ulong r = z_randint(30);
		ulong c = z_randint(30);
		F_mpz_mat_init(F_mat1, r, c);
      mpz_mat_init(m_mat1, r, c); 
      
		bits1 = z_randint(200) + 1;
      mpz_randmat(m_mat1, r, c, bits1);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
          
      result = (F_mpz_mat_equal(F_mat1, F_mat1)); 
		if (!result) 
		{
			printf("Error: length1 = %ld, bits1 = %ld\n", length1, bits1);
		}
          
      F_mpz_mat_clear(F_mat1);
		mpz_mat_clear(m_mat1);
   }

   return result;
}

int test_F_mpz_mat_resize()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat1, F_mat2;
   int result = 1;
   ulong bits, r, c, r_add, c_add;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30);
      c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      mpz_mat_init(m_mat1, r, c); 
      
      mpz_randmat(m_mat1, r, c, bits);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      F_mpz_mat_set(F_mat2, F_mat1);
		
		r_add = z_randint(30);
      c_add = z_randint(30);
		
		F_mpz_mat_resize(F_mat1, r + r_add, c + c_add);
      F_mpz_mat_resize(F_mat1, r, c);
      
      result = F_mpz_mat_equal(F_mat1, F_mat2); 
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, r_add = %ld, c_add = %ld\n", r, c, bits, r_add, c_add);
		}
          
      F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
		mpz_mat_clear(m_mat1);
   }
   
   return result;
}

int test_F_mpz_mat_neg()
{
   F_mpz_mat_t F_mat1, F_mat2, F_mat3, F_mat4;
   int result = 1;
   ulong bits1, bits2, r, c;
   
   // negate and negate back
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
   {
      r = z_randint(30);
		c = z_randint(30);
      F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
		bits1 = z_randint(200) + 1;
      F_mpz_randmat(F_mat1, r, c, bits1);
      
		F_mpz_mat_neg(F_mat2, F_mat1);
      F_mpz_mat_neg(F_mat2, F_mat2);
          
      result = (F_mpz_mat_equal(F_mat1, F_mat2)); 
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits1 = %ld\n", r, c, bits1);
		}
          
      F_mpz_mat_clear(F_mat1);
		F_mpz_mat_clear(F_mat2);
   }

	// sub equals negate and add, included aliased negation
	for (ulong count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
   {
      r = z_randint(30);
      c = z_randint(30);
      
		F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(F_mat3, r, c);
      F_mpz_mat_init(F_mat4, r, c);
      
		bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      F_mpz_randmat(F_mat1, r, c, bits1);
      F_mpz_randmat(F_mat2, r, c, bits2);
      
		F_mpz_mat_sub(F_mat3, F_mat1, F_mat2);
      F_mpz_mat_neg(F_mat2, F_mat2);
      F_mpz_mat_add(F_mat4, F_mat1, F_mat2);
          
      result = (F_mpz_mat_equal(F_mat3, F_mat4)); 
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits1 = %ld, bits2 = %ld\n", r, c, bits1, bits2);
		}
          
      F_mpz_mat_clear(F_mat1);
		F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
		F_mpz_mat_clear(F_mat4);
   }

   return result;
}

int test_F_mpz_mat_row_mul_ui()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_mat_row_mul_ui(F_mat2, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_mul_ui(temp, m_mat->entries[r1*c+i], mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_mat_row_mul_ui(F_mat, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_mul_ui(temp, m_mat->entries[r1*c+i], mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
    
	mpz_clear(temp);
	
   return result; 
}

int test_F_mpz_mat_row_mul_si()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
		if (z_randint(2)) mult = -mult;
      
		F_mpz_mat_row_mul_si(F_mat2, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);

      for (ulong i = start; i < start + n; i++)
      {
         if (mult > 0L) mpz_mul_ui(temp, m_mat->entries[r1*c+i], mult);
			else
			{
				mpz_mul_ui(temp, m_mat->entries[r1*c+i], -mult);
				mpz_neg(temp, temp);
			}
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
			if (!result) 
		   {
			   printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
			   gmp_printf("%Zd, %Zd\n", temp, m_mat2->entries[r2*c+i]);
			   break;
		   }
      }

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      if (z_randint(2)) mult = -mult;
      
		F_mpz_mat_row_mul_si(F_mat, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         if (mult > 0L) mpz_mul_ui(temp, m_mat->entries[r1*c+i], mult);
			else
			{
				mpz_mul_ui(temp, m_mat->entries[r1*c+i], -mult);
				mpz_neg(temp, temp);
			}
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
         if (!result) 
		   {
			   printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
				gmp_printf("Aliased: %Zd, %Zd\n", temp, m_mat2->entries[r2*c+i]);
			   break;
		   }      
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
      
   mpz_clear(temp);
	
	return result; 
}
      
int test_F_mpz_mat_row_mul_F_mpz()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   
	mpz_t temp, temp2;
	mpz_init(temp);
   mpz_init(temp2);
   F_mpz_t mult;
	F_mpz_init(mult);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(200);
      mpz_rrandomb(temp2, randstate, bits2);
		if (z_randint(2)) mpz_neg(temp2, temp2);
		F_mpz_set_mpz(mult, temp2);

		F_mpz_mat_row_mul_F_mpz(F_mat2, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_mul(temp, m_mat->entries[r1*c+i], temp2);
			
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
			if (!result) 
		   {
			   printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, ", r, c, bits, bits2);
				printf("mult = "); F_mpz_print(mult); printf("\n");
			   gmp_printf("%Zd, %Zd\n", temp, m_mat2->entries[r2*c+i]);
			   break;
		   }
      }

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(200);
      mpz_rrandomb(temp2, randstate, bits2);
		if (z_randint(2)) mpz_neg(temp2, temp2);
      F_mpz_set_mpz(mult, temp2);

		F_mpz_mat_row_mul_F_mpz(F_mat, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_mul(temp, m_mat->entries[r1*c+i], temp2);
			
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
         if (!result) 
		   {
			   printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, ", r, c, bits, bits2);
				printf("mult = "); F_mpz_print(mult); printf("\n");
			   gmp_printf("%Zd, %Zd\n", temp, m_mat2->entries[r2*c+i]);
			   break;
		   }      
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
      
	F_mpz_clear(mult);
	mpz_clear(temp);
	mpz_clear(temp2);
	
   return result; 
}

int test_F_mpz_mat_row_addmul_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      bits2 = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		mpz_mat_init(m_mat3, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      mpz_randmat_dense(m_mat2, r, c, bits2); 
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
      bits3 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits3);
      
		F_mpz_mat_row_addmul_ui(F_mat2, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat2->entries[r2*c+i]);
			mpz_addmul_ui(temp, m_mat->entries[r1*c+i], mult);
         result &= (mpz_cmp(temp, m_mat3->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		mpz_mat_clear(m_mat3); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_mat_row_addmul_ui(F_mat, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat->entries[r2*c+i]);
			mpz_addmul_ui(temp, m_mat->entries[r1*c+i], mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
    
	mpz_clear(temp);
	
   return result; 
}

int test_F_mpz_mat_row_submul_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      bits2 = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		mpz_mat_init(m_mat3, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      mpz_randmat_dense(m_mat2, r, c, bits2); 
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
      bits3 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits3);
      
		F_mpz_mat_row_submul_ui(F_mat2, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat2->entries[r2*c+i]);
			mpz_submul_ui(temp, m_mat->entries[r1*c+i], mult);
         result &= (mpz_cmp(temp, m_mat3->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		mpz_mat_clear(m_mat3); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_mat_row_submul_ui(F_mat, r2, F_mat, r1, start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat->entries[r2*c+i]);
			mpz_submul_ui(temp, m_mat->entries[r1*c+i], mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
    
	mpz_clear(temp);
	
   return result; 
}

int test_F_mpz_mat_row_addmul_2exp_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2, exp;
   ulong mult;
   mpz_t temp, temp2;
   mpz_init(temp);
   mpz_init(temp2);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      bits2 = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		exp = z_randint(100);
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		mpz_mat_init(m_mat3, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      mpz_randmat_dense(m_mat2, r, c, bits2); 
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
      bits3 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits3);
      
		F_mpz_mat_row_addmul_2exp_ui(F_mat2, r2, F_mat, r1, start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);
      
      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat2->entries[r2*c+i]);
			mpz_mul_2exp(temp2, m_mat->entries[r1*c+i], exp);
			mpz_addmul_ui(temp, temp2, mult);
         result &= (mpz_cmp(temp, m_mat3->entries[r2*c+i]) == 0);
      }
      
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}
      
      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		mpz_mat_clear(m_mat3); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      exp = z_randint(100);
		start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_mat_row_addmul_2exp_ui(F_mat, r2, F_mat, r1, start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat->entries[r2*c+i]);
			mpz_mul_2exp(temp2, m_mat->entries[r1*c+i], exp);
			mpz_addmul_ui(temp, temp2, mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
      F_mpz_mat_clear(F_mat);
   }
    
	mpz_clear(temp);
	mpz_clear(temp2);
	
   return result; 
}

int test_F_mpz_mat_row_submul_2exp_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2, exp;
   ulong mult;
   mpz_t temp, temp2;
   mpz_init(temp);
   mpz_init(temp2);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      bits2 = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		exp = z_randint(100);
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		mpz_mat_init(m_mat3, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      mpz_randmat_dense(m_mat2, r, c, bits2); 
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
      bits3 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits3);
      
		F_mpz_mat_row_submul_2exp_ui(F_mat2, r2, F_mat, r1, start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);
      
      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat2->entries[r2*c+i]);
			mpz_mul_2exp(temp2, m_mat->entries[r1*c+i], exp);
			mpz_submul_ui(temp, temp2, mult);
         result &= (mpz_cmp(temp, m_mat3->entries[r2*c+i]) == 0);
      }
      
		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}
      
      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		mpz_mat_clear(m_mat3); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (ulong count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      exp = z_randint(100);
		start = z_randint(c);
		n = z_randint(c-start)+1;

      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits2 = z_randint(FLINT_BITS+1);
		mult = z_randbits(bits2);
      
		F_mpz_mat_row_submul_2exp_ui(F_mat, r2, F_mat, r1, start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      for (ulong i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat->entries[r2*c+i]);
			mpz_mul_2exp(temp2, m_mat->entries[r1*c+i], exp);
			mpz_submul_ui(temp, temp2, mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %ld\n", r, c, bits, bits2, mult);
		}

      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
    
	mpz_clear(temp);
	mpz_clear(temp2);
	
   return result; 
}

int test_F_mpz_mat_row_scalar_mul()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   F_mpz_t F_sp;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   mpz_t m_sp1, m_sp2;
   mpz_init(m_sp1);
   mpz_init(m_sp2);
   
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      bits2 = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      F_mpz_init(F_sp);
      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      mpz_randmat_dense(m_mat2, r, c, bits2); 
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
		F_mpz_mat_row_scalar_product(F_sp, F_mat2, r2, F_mat, r1, start, n);
      F_mpz_get_mpz(m_sp1, F_sp);
      
      mpz_set_ui(m_sp2, 0);
      for (ulong i = start; i < start + n; i++)
      {
         mpz_addmul(m_sp2, m_mat->entries[r1*c+i], m_mat2->entries[r2*c+i]);
      }
      
      result = (mpz_cmp(m_sp1, m_sp2) == 0);

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld\n", r, c, bits, bits2);
		}
      
      F_mpz_clear(F_sp);
      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // alias rows
   for (ulong count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      start = z_randint(c);
		n = z_randint(c-start)+1;

      F_mpz_init(F_sp);
      mpz_mat_init(m_mat, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
		F_mpz_mat_row_scalar_product(F_sp, F_mat, r1, F_mat, r1, start, n);
      F_mpz_get_mpz(m_sp1, F_sp);
      
      mpz_set_ui(m_sp2, 0);
      for (ulong i = start; i < start + n; i++)
      {
         mpz_addmul(m_sp2, m_mat->entries[r1*c+i], m_mat->entries[r1*c+i]);
      }
      
      result = (mpz_cmp(m_sp1, m_sp2) == 0);

		if (!result) 
		{
			printf("Error: r = %ld, c = %ld, bits = %ld\n", r, c, bits);
		}
      
      F_mpz_clear(F_sp);
      mpz_mat_clear(m_mat); 
		F_mpz_mat_clear(F_mat);
   }
       
	mpz_clear(m_sp1);
	mpz_clear(m_sp2);
	
   return result; 
}

void F_mpz_mat_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(F_mpz_mat_convert); 
   RUN_TEST(F_mpz_mat_set); 
   RUN_TEST(F_mpz_mat_equal);  
   RUN_TEST(F_mpz_mat_resize);  
   RUN_TEST(F_mpz_mat_neg); 
   RUN_TEST(F_mpz_mat_add); 
   RUN_TEST(F_mpz_mat_sub); 
   RUN_TEST(F_mpz_mat_row_mul_ui); 
   RUN_TEST(F_mpz_mat_row_mul_si); 
   RUN_TEST(F_mpz_mat_row_mul_F_mpz); 
   RUN_TEST(F_mpz_mat_row_addmul_ui); 
   RUN_TEST(F_mpz_mat_row_submul_ui); 
   RUN_TEST(F_mpz_mat_row_addmul_2exp_ui); 
   RUN_TEST(F_mpz_mat_row_submul_2exp_ui); 
   RUN_TEST(F_mpz_mat_row_scalar_mul); 
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_mat_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


