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
#include <mpfr.h>
#include <time.h>
#include "flint.h"
#include "long_extras.h"
#include "d_mat.h"
#include "mpz_mat.h"
#include "F_mpz_mat.h"
#include "mpfr_mat.h"
#include "memory-manager.h"
#include "test-support.h"

#define VARY_BITS 1 // random entries have random number of bits up to the limit given
#define SIGNS 1 // random entries will be randomly signed
#define SPARSE 1 // matrices are sparse (triggers more corner cases)
#define ITER 1 // if you want all tests to run longer, increase this

#define TESTFILE 0 // Set this to test matrix reading and writing to a file in the current dir

#define DEBUG 0 // allows easy switching of debugging code on and off when debugging (if inserted)
#define DEBUG2 1 

void F_mpz_test_random(F_mpz_t f, ulong bits)
{
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

// generate a dense random mpz_mat_t with up to the given length and number of bits per entry
void mpz_randmat_dense(mpz_mat_t mat, ulong r, ulong c, ulong maxbits)
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

int test__F_mpz_vec_init_clear()
{
   F_mpz * vec;
   int result = 1;
   ulong bits, c;
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      ulong i;
	  bits = z_randint(200)+ 1;
      c = z_randint(30)+1;
	  
      vec = _F_mpz_vec_init(c);
	  for (i = 0; i < c; i++)
		 F_mpz_test_random(vec + i, bits);
	  _F_mpz_vec_clear(vec, c);
   }

   return result; 
}

int test__F_mpz_vec_copy_equal()
{
   F_mpz * vec1, * vec2;
   int result = 1;
   ulong bits, c, c1;
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      ulong i;
	  bits = z_randint(200) + 1;
      c = z_randint(30) + 1;
	  
      vec1 = _F_mpz_vec_init(c);
	  vec2 = _F_mpz_vec_init(c);

	  for (i = 0; i < c; i++)
		 F_mpz_test_random(vec1 + i, bits);

	  _F_mpz_vec_copy(vec2, vec1, c);
	  
	  result = _F_mpz_vec_equal(vec1, vec2, c);

	  _F_mpz_vec_clear(vec1, c);
	  _F_mpz_vec_clear(vec2, c);
   }

   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      ulong i;
	  bits = z_randint(200) + 1;
      c = z_randint(30) + 1;
	  
      vec1 = _F_mpz_vec_init(c);
	  vec2 = _F_mpz_vec_init(c);

	  for (i = 0; i < c; i++)
		 F_mpz_test_random(vec1 + i, bits);

	  _F_mpz_vec_copy(vec2, vec1, c);

	  c1 = z_randint(c);
	  F_mpz_add_ui(vec2 + c1, vec2 + c1, 1);

	  result = (!_F_mpz_vec_equal(vec1, vec2, c));

	  _F_mpz_vec_clear(vec1, c);
	  _F_mpz_vec_clear(vec2, c);
   }

   return result; 
}

int test_F_mpz_mat_swap_rows()
{
   F_mpz_mat_t F_mat;
   F_mpz * vec;
   int result = 1;
   ulong bits, r, c, r1, r2;
   
   F_mpz_t mult;
   F_mpz_init(mult);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      
      F_mpz_randmat(F_mat, r, c, bits); 
          
	  vec = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec, F_mat->rows[r2], c);
	  F_mpz_mat_swap_rows(F_mat, r1, r2);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         result &= (F_mpz_cmp(vec + i, F_mat->rows[r1] + i) == 0);
			
		 if (!result) 
		 {
			 printf("Error: r = %ld, c = %ld, bits = %ld\n", r, c, bits);
			 F_mpz_print(vec + i);
			 F_mpz_print(F_mat->rows[r1] + i);
			 break;
		 }
      }

      _F_mpz_vec_clear(vec, c);

      F_mpz_mat_clear(F_mat);
   }
   
   // alias rows
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      
      F_mpz_randmat(F_mat, r, c, bits); 
          
	  vec = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec, F_mat->rows[r1], c);
	  F_mpz_mat_swap_rows(F_mat, r1, r1);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         result &= (F_mpz_cmp(vec + i, F_mat->rows[r1] + i) == 0);
			
		 if (!result) 
		 {
			 printf("Error: r = %ld, c = %ld, bits = %ld\n", r, c, bits);
			 F_mpz_print(vec + i);
			 F_mpz_print(F_mat->rows[r1] + i);
			 break;
		 }
      }

      _F_mpz_vec_clear(vec, c);

      F_mpz_mat_clear(F_mat);
   }

   return result; 
}

int test_F_mpz_mat_convert()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat;
   int result = 1;
   ulong bits, r, c;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
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
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
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
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
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
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
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
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
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
   
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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
	ulong count1;
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
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

int test_F_mpz_mat_swap()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat1, F_mat2, F_mat3;
   int result = 1;
   ulong bits1, bits2, length1, length2;

   ulong count1;
   for (count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
   {
      ulong r = z_randint(30);
	  ulong c = z_randint(30);
	  F_mpz_mat_init(F_mat1, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(F_mat3, r, c);
      
	  bits1 = z_randint(200) + 1;
      F_mpz_randmat(F_mat1, r, c, bits1);
      F_mpz_randmat(F_mat2, r, c, bits1);
           
      F_mpz_mat_set(F_mat3, F_mat1);
      F_mpz_mat_swap(F_mat2, F_mat1);
          
      result = (F_mpz_mat_equal(F_mat2, F_mat3)); 
	  if (!result) 
	  {
		 printf("Error: r = %ld, c = %ld, bits = %ld\n", r, c, bits1);
	  }
          
      F_mpz_mat_clear(F_mat1);
	  F_mpz_mat_clear(F_mat2);
	  F_mpz_mat_clear(F_mat3);
   }

   return result;
}

int test_F_mpz_mat_resize()
{
   mpz_mat_t m_mat1, m_mat2;
   F_mpz_mat_t F_mat1, F_mat2;
   int result = 1;
   ulong bits, r, c, r_add, c_add;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
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

int test_F_mpz_mat_tofromstring()
{
   mpz_mat_t test_mat;
   F_mpz_mat_t test_F_mpz_mat, test_F_mpz_mat2;
   int result = 1;
   ulong bits;


   ulong count1;
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r1 = z_randint(30);
      ulong c1 = z_randint(30);
      ulong r2 = z_randint(30);
      ulong c2 = z_randint(30);
      
		F_mpz_mat_init(test_F_mpz_mat, r1, c1);
      F_mpz_mat_init(test_F_mpz_mat2, r2, c2);

      bits = z_randint(200) + 1;

      mpz_mat_init(test_mat, r1, c1);
      
		mpz_randmat(test_mat, r1, c1, bits);
           
      mpz_mat_to_F_mpz_mat(test_F_mpz_mat, test_mat);

      char *strbuf = F_mpz_mat_to_string(test_F_mpz_mat);      

      int OK = F_mpz_mat_from_string(test_F_mpz_mat2, strbuf);

      free(strbuf);

      result = F_mpz_mat_equal(test_F_mpz_mat, test_F_mpz_mat2) && OK; 
		if (!result) 
		{
			printf("Error: r1 = %ld, c1 = %ld, r2 = %ld, c2 = %ld\n", r1, c1, r2, c2);
		}
          
		F_mpz_mat_clear(test_F_mpz_mat);
      F_mpz_mat_clear(test_F_mpz_mat2);
      mpz_mat_clear(test_mat); 

   }

   
   return result;
}

int test_F_mpz_mat_tofromstringpretty()
{
   mpz_mat_t test_mat;
   F_mpz_mat_t test_F_mpz_mat, test_F_mpz_mat2;
   int result = 1;
   ulong bits;

   ulong count1;
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r1 = z_randint(30);
      ulong c1 = z_randint(30);
      ulong r2 = z_randint(30);
      ulong c2 = z_randint(30);
      
		F_mpz_mat_init(test_F_mpz_mat, r1, c1);
      F_mpz_mat_init(test_F_mpz_mat2, r2, c2);

      bits = z_randint(200) + 1;

      mpz_mat_init(test_mat, r1, c1);
      
		mpz_randmat(test_mat, r1, c1, bits);
           
      mpz_mat_to_F_mpz_mat(test_F_mpz_mat, test_mat);
      char *strbuf = F_mpz_mat_to_string_pretty(test_F_mpz_mat);      

      int OK = F_mpz_mat_from_string_pretty(test_F_mpz_mat2, strbuf);

      free(strbuf);

      result = F_mpz_mat_equal(test_F_mpz_mat, test_F_mpz_mat2) && OK; 
		if (!result) 
		{
			printf("Error: r1 = %ld, c1 = %ld, r2 = %ld, c2 = %ld\n", r1, c1, r2, c2);
		}
		F_mpz_mat_clear(test_F_mpz_mat);
      F_mpz_mat_clear(test_F_mpz_mat2);
      mpz_mat_clear(test_mat); 
   }  
   return result;
}

int test_F_mpz_mat_neg()
{
   F_mpz_mat_t F_mat1, F_mat2, F_mat3, F_mat4;
   int result = 1;
   ulong bits1, bits2, r, c;
   
   // negate and negate back
	ulong count1;
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
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
	for (count1 = 0; (count1 < 5000*ITER) && (result == 1); count1++)
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

int test__F_mpz_vec_add_sub()
{
   F_mpz_mat_t F_mat, F_mat2, F_mat3;
   int result = 1;
   ulong bits, r, c, r1, r2;
   
   ulong count1;
   // alias sub mat1-mat2
   for (count1 = 0; (count1 < 2000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(F_mat3, r, c);
      
      F_mpz_randmat(F_mat2, r, c, bits); 
      F_mpz_randmat(F_mat3, r, c, bits); 
              
	  _F_mpz_vec_add(F_mat->rows[r2], F_mat2->rows[r2], F_mat3->rows[r1], c);
	  _F_mpz_vec_sub(F_mat->rows[r2], F_mat->rows[r2], F_mat3->rows[r1], c);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         result &= (F_mpz_cmp(F_mat->rows[r2] + i, F_mat2->rows[r2] + i) == 0);
      }

	  if (!result) 
	  {
		 printf("Error: r = %ld, r1 = %ld, r2 = %ld, c = %ld, bits = %ld\n", r, r1, r2, c, bits);
	  }

      F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
   }
   
   // alias sub mat1-mat3
   for (count1 = 0; (count1 < 2000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(F_mat3, r, c);
      
      F_mpz_randmat(F_mat2, r, c, bits); 
      F_mpz_randmat(F_mat3, r, c, bits); 
              
	  _F_mpz_vec_add(F_mat->rows[r2], F_mat2->rows[r2], F_mat3->rows[r1], c);
	  _F_mpz_vec_sub(F_mat3->rows[r1], F_mat->rows[r2], F_mat3->rows[r1], c);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         result &= (F_mpz_cmp(F_mat2->rows[r2] + i, F_mat3->rows[r1] + i) == 0);
      }

	  if (!result) 
	  {
		 printf("Error: r = %ld, r1 = %ld, r2 = %ld, c = %ld, bits = %ld\n", r, r1, r2, c, bits);
	  }

      F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
   }
	
   // alias add mat1-mat2
   for (count1 = 0; (count1 < 2000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(F_mat3, r, c);
      
      F_mpz_randmat(F_mat2, r, c, bits); 
      F_mpz_randmat(F_mat3, r, c, bits); 
              
	  _F_mpz_vec_sub(F_mat->rows[r2], F_mat2->rows[r2], F_mat3->rows[r1], c);
	  _F_mpz_vec_add(F_mat->rows[r2], F_mat->rows[r2], F_mat3->rows[r1], c);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         result &= (F_mpz_cmp(F_mat->rows[r2] + i, F_mat2->rows[r2] + i) == 0);
      }

	  if (!result) 
	  {
		 printf("Error: r = %ld, r1 = %ld, r2 = %ld, c = %ld, bits = %ld\n", r, r1, r2, c, bits);
	  }

      F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
   }
   
   // alias add mat1-mat3
   for (count1 = 0; (count1 < 2000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      F_mpz_mat_init(F_mat3, r, c);
      
      F_mpz_randmat(F_mat2, r, c, bits); 
      F_mpz_randmat(F_mat3, r, c, bits); 
              
	  _F_mpz_vec_sub(F_mat->rows[r2], F_mat2->rows[r2], F_mat3->rows[r1], c);
	  _F_mpz_vec_add(F_mat3->rows[r1], F_mat->rows[r2], F_mat3->rows[r1], c);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         result &= (F_mpz_cmp(F_mat2->rows[r2] + i, F_mat3->rows[r1] + i) == 0);
      }

	  if (!result) 
	  {
		 printf("Error: r = %ld, r1 = %ld, r2 = %ld, c = %ld, bits = %ld\n", r, r1, r2, c, bits);
	  }

      F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
   }
	
   return result; 
}

int test__F_mpz_vec_neg()
{
   F_mpz_mat_t F_mat, F_mat2;
   F_mpz * vec;
   int result = 1;
   ulong bits, r, c, r1, r2;
   
   ulong count1;
   for (count1 = 0; (count1 < 2000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      F_mpz_randmat(F_mat, r, c, bits); 
      F_mpz_randmat(F_mat2, r, c, bits); 
              
	  vec = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec, F_mat2->rows[r2], c);
	  _F_mpz_vec_neg(F_mat->rows[r1], F_mat2->rows[r2], c);
	  _F_mpz_vec_neg(vec, F_mat->rows[r1], c);
	  
      ulong i;
      result = _F_mpz_vec_equal(vec, F_mat2->rows[r2], c);

	  if (!result) 
	  {
		 printf("Error: r = %ld, r1 = %ld, r2 = %ld, c = %ld, bits = %ld\n", r, r1, r2, c, bits);
	  }

      _F_mpz_vec_clear(vec, c);
	  
	  F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }

   // alias rows
   for (count1 = 0; (count1 < 2000) && (result == 1) ; count1++)
   {
      bits = z_randint(200) + 1;
      r = z_randint(30) + 1;  
      c = z_randint(30) + 1;
	  r1 = z_randint(r);
      
	  F_mpz_mat_init(F_mat, r, c);
      
      F_mpz_randmat(F_mat, r, c, bits); 
              
	  vec = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec, F_mat->rows[r1], c);
	  _F_mpz_vec_neg(F_mat->rows[r1], F_mat->rows[r1], c);
	  _F_mpz_vec_neg(F_mat->rows[r1], F_mat->rows[r1], c);
	  
      ulong i;
      result = _F_mpz_vec_equal(vec, F_mat->rows[r1], c);

	  if (!result) 
	  {
		 printf("Error: r = %ld, r1 = %ld, c = %ld, bits = %ld\n", r, r1, c, bits);
	  }

      _F_mpz_vec_clear(vec, c);
	  
	  F_mpz_mat_clear(F_mat);
   }

   return result; 
}

int test__F_mpz_vec_scalar_mul_ui()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_mul_ui(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);

      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_mul_ui(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      ;
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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

int test__F_mpz_vec_scalar_mul_si()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_mul_si(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);

      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_mul_si(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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
      
int test__F_mpz_vec_scalar_mul_F_mpz()
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
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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

		_F_mpz_vec_scalar_mul_F_mpz(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat2);

      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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

		_F_mpz_vec_scalar_mul_F_mpz(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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

int test__F_mpz_vec_scalar_addmul_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_addmul_ui(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);

      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_addmul_ui(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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

int test__F_mpz_vec_scalar_addmul_F_mpz()
{
   F_mpz_mat_t F_mat, F_mat2;
   F_mpz * vec, * vec1, * vec2;
   int result = 1;
   ulong bits, bits2, r, c, r1, r2;
   
   F_mpz_t mult;
   F_mpz_init(mult);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      F_mpz_randmat(F_mat, r, c, bits); 
      F_mpz_randmat(F_mat2, r, c, bits); 
          
      bits2 = z_randint(200);
      F_mpz_test_random(mult, bits2);

	  vec = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec, F_mat2->rows[r2], c);
	  _F_mpz_vec_scalar_addmul_F_mpz(F_mat2->rows[r2], F_mat->rows[r1], c, mult);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         F_mpz_addmul(vec + i, F_mat->rows[r1] + i, mult);
			
         result &= (F_mpz_cmp(vec + i, F_mat2->rows[r2] + i) == 0);
			
		 if (!result) 
		 {
			 printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, ", r, c, bits, bits2);
			 printf("mult = "); F_mpz_print(mult); printf("\n");
			 F_mpz_print(vec + i);
			 F_mpz_print(F_mat2->rows[r2] + i);
			 break;
		 }
      }

      _F_mpz_vec_clear(vec, c);

      F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // alias mats
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat2, r, c);
      
      F_mpz_randmat(F_mat2, r, c, bits); 
          
      bits2 = z_randint(200);
      F_mpz_test_random(mult, bits2);

	  vec1 = _F_mpz_vec_init(c);
	  vec2 = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec2, F_mat2->rows[r2], c);
	  _F_mpz_vec_copy(vec1, F_mat2->rows[r1], c);
	  _F_mpz_vec_scalar_addmul_F_mpz(F_mat2->rows[r2], F_mat2->rows[r1], c, mult);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         F_mpz_addmul(vec2 + i, vec1 + i, mult);
			
         result &= (F_mpz_cmp(vec2 + i, F_mat2->rows[r2] + i) == 0);
			
		 if (!result) 
		 {
			 printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, ", r, c, bits, bits2);
			 printf("mult = "); F_mpz_print(mult); printf("\n");
			 F_mpz_print(vec2 + i);
			 F_mpz_print(F_mat2->rows[r2] + i);
			 break;
		 }
      }

      _F_mpz_vec_clear(vec1, c);
      _F_mpz_vec_clear(vec2, c);

      F_mpz_mat_clear(F_mat2);
   }
	
   F_mpz_clear(mult);
	
   return result; 
}

int test__F_mpz_vec_scalar_submul_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2;
   ulong mult;
   mpz_t temp;
   mpz_init(temp);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_submul_ui(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);

      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_submul_ui(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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

int test__F_mpz_vec_scalar_submul_F_mpz()
{
   F_mpz_mat_t F_mat, F_mat2;
   F_mpz * vec, * vec1, * vec2;
   int result = 1;
   ulong bits, bits2, r, c, r1, r2;
   
   F_mpz_t mult;
   F_mpz_init(mult);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      F_mpz_randmat(F_mat, r, c, bits); 
      F_mpz_randmat(F_mat2, r, c, bits); 
          
      bits2 = z_randint(200);
      F_mpz_test_random(mult, bits2);

	  vec = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec, F_mat2->rows[r2], c);
	  _F_mpz_vec_scalar_submul_F_mpz(F_mat2->rows[r2], F_mat->rows[r1], c, mult);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         F_mpz_submul(vec + i, F_mat->rows[r1] + i, mult);
			
         result &= (F_mpz_cmp(vec + i, F_mat2->rows[r2] + i) == 0);
			
		 if (!result) 
		 {
			 printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, ", r, c, bits, bits2);
			 printf("mult = "); F_mpz_print(mult); printf("\n");
			 F_mpz_print(vec + i);
			 F_mpz_print(F_mat2->rows[r2] + i);
			 break;
		 }
      }

      _F_mpz_vec_clear(vec, c);

      F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // alias mats
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
	  r1 = z_randint(r);
      r2 = z_randint(r);
      
      F_mpz_mat_init(F_mat2, r, c);
      
      F_mpz_randmat(F_mat2, r, c, bits); 
          
      bits2 = z_randint(200);
      F_mpz_test_random(mult, bits2);

	  vec1 = _F_mpz_vec_init(c);
	  vec2 = _F_mpz_vec_init(c);
	  _F_mpz_vec_copy(vec2, F_mat2->rows[r2], c);
	  _F_mpz_vec_copy(vec1, F_mat2->rows[r1], c);
	  _F_mpz_vec_scalar_submul_F_mpz(F_mat2->rows[r2], F_mat2->rows[r1], c, mult);
      
      ulong i;
      for (i = 0; i < c; i++)
      {
         F_mpz_submul(vec2 + i, vec1 + i, mult);
			
         result &= (F_mpz_cmp(vec2 + i, F_mat2->rows[r2] + i) == 0);
			
		 if (!result) 
		 {
			 printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, ", r, c, bits, bits2);
			 printf("mult = "); F_mpz_print(mult); printf("\n");
			 F_mpz_print(vec2 + i);
			 F_mpz_print(F_mat2->rows[r2] + i);
			 break;
		 }
      }

      _F_mpz_vec_clear(vec1, c);
      _F_mpz_vec_clear(vec2, c);

      F_mpz_mat_clear(F_mat2);
   }
	
   F_mpz_clear(mult);
	
   return result; 
}

int test__F_mpz_vec_scalar_addmul_2exp_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2, exp;
   ulong mult;
   mpz_t temp, temp2;
   mpz_init(temp);
   mpz_init(temp2);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_addmul_2exp_ui(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);
      
      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_addmul_2exp_ui(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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

int test__F_mpz_vec_scalar_submul_2exp_ui()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2, exp;
   ulong mult;
   mpz_t temp, temp2;
   mpz_init(temp);
   mpz_init(temp2);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_submul_2exp_ui(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);
      
      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
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
      
		_F_mpz_vec_scalar_submul_2exp_ui(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
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

int test_F_mpz_mat_mul_classical()
{
   mpz_mat_t m_mat1, m_mat2, m_mat3,res1,res2;
   F_mpz_mat_t F_mat1, F_mat2, F_mat3, F_res1, F_res2, F_res3, F_temp;
   int result = 1;
   ulong bits1, bits2, bits3;
   
   ulong count1;
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r1 = z_randint(30)+1;
      ulong c1 = z_randint(30)+1;
      ulong c2 = z_randint(30)+1;

      
      F_mpz_mat_init(F_mat1, r1, c1);
      F_mpz_mat_init(F_mat2, c1, c2);
      F_mpz_mat_init(F_mat3, c1, c2);
      F_mpz_mat_init(F_temp, c1, c2);
      F_mpz_mat_init(F_res1, r1, c2);
      F_mpz_mat_init(F_res2, r1, c2);
      F_mpz_mat_init(F_res3, r1, c2);

      mpz_mat_init(m_mat1, r1, c1); 
      mpz_mat_init(m_mat2, c1, c2); 
      mpz_mat_init(m_mat3, c1, c2); 
      mpz_mat_init(res1, r1, c2); 
      mpz_mat_init(res2, r1, c2); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      bits3 = z_randint(200) + 1;
      
      mpz_randmat(m_mat1, r1, c1, bits1);
      mpz_randmat(m_mat2, c1, c2, bits2);
      mpz_randmat(m_mat3, c1, c2, bits2);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      mpz_mat_to_F_mpz_mat(F_mat3, m_mat3);

      
      F_mpz_mat_add(F_temp,F_mat2,F_mat3);     
      F_mpz_mat_mul_classical(F_res1, F_mat1, F_temp);
      F_mpz_mat_to_mpz_mat(res1, F_res1);


      F_mpz_mat_mul_classical(F_res3,F_mat1,F_mat2);
      F_mpz_mat_mul_classical(F_res2,F_mat1,F_mat3);
      F_mpz_mat_add(F_res2,F_res2,F_res3);
      F_mpz_mat_to_mpz_mat(res2, F_res2);

          
      result = mpz_mat_equal(res1, res2); 
      if (!result) 
      {
         printf("Error: bits1 = %ld, bits2 = %ld, bits3 = %ld, count1 = %ld\n", bits1, bits2, bits3, count1);
      }
          
      F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
      F_mpz_mat_clear(F_temp);
      F_mpz_mat_clear(F_res1);
      F_mpz_mat_clear(F_res2);
      F_mpz_mat_clear(F_res3);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);  
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2); 
   }

//alias testing for the square case
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30)+1;
      
      F_mpz_mat_init(F_mat1, r, r);
      F_mpz_mat_init(F_mat2, r, r);
      F_mpz_mat_init(F_mat3, r, r);
      F_mpz_mat_init(F_temp, r, r);
      F_mpz_mat_init(F_res1, r, r);
      F_mpz_mat_init(F_res2, r, r);
      F_mpz_mat_init(F_res3, r, r);

      mpz_mat_init(m_mat1, r, r); 
      mpz_mat_init(m_mat2, r, r); 
      mpz_mat_init(m_mat3, r, r); 
      mpz_mat_init(res1, r, r); 
      mpz_mat_init(res2, r, r); 

      bits1 = z_randint(200) + 1;
      bits2 = z_randint(200) + 1;
      bits3 = z_randint(200) + 1;
      
      mpz_randmat(m_mat1, r, r, bits1);
      mpz_randmat(m_mat2, r, r, bits2);
      mpz_randmat(m_mat3, r, r, bits2);
           
      mpz_mat_to_F_mpz_mat(F_res2, m_mat1);
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
      mpz_mat_to_F_mpz_mat(F_res3, m_mat3);

      F_mpz_mat_add(F_res1,F_mat2,F_res3);    
      F_mpz_mat_mul_classical(F_res1, F_res2, F_res1);
      F_mpz_mat_to_mpz_mat(res1, F_res1);

      F_mpz_mat_mul_classical(F_res3,F_res2,F_res3);
      F_mpz_mat_mul_classical(F_res2,F_res2,F_mat2);
      F_mpz_mat_add(F_res2,F_res2,F_res3);
      F_mpz_mat_to_mpz_mat(res2, F_res2);
          
      result = mpz_mat_equal(res1, res2); 
      if (!result) 
      {
         printf("Error: bits1 = %ld, bits2 = %ld, bits3 = %ld, count1 = %ld\n", bits1, bits2, bits3, count1);
      }
          
      F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_mat2);
      F_mpz_mat_clear(F_mat3);
      F_mpz_mat_clear(F_temp);
      F_mpz_mat_clear(F_res1);
      F_mpz_mat_clear(F_res2);
      F_mpz_mat_clear(F_res3);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(m_mat2);
      mpz_mat_clear(m_mat3);  
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2); 
   }

//aliasing square case two inputs the same and all inputs same
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(30)+1;
      
      F_mpz_mat_init(F_mat1, r, r);
      F_mpz_mat_init(F_temp, r, r);
      F_mpz_mat_init(F_res1, r, r);
      F_mpz_mat_init(F_res2, r, r);
      mpz_mat_init(m_mat1, r, r); 
      mpz_mat_init(res1, r, r); 
      mpz_mat_init(res2, r, r); 

      bits1 = z_randint(200) + 1;
      
      mpz_randmat(m_mat1, r, r, bits1);
           
      mpz_mat_to_F_mpz_mat(F_mat1, m_mat1);

      F_mpz_mat_add(F_temp,F_mat1,F_mat1);    
      F_mpz_mat_mul_classical(F_res1, F_mat1, F_temp);
      F_mpz_mat_to_mpz_mat(res1, F_res1);

      F_mpz_mat_mul_classical(F_mat1,F_mat1,F_mat1);
      F_mpz_mat_add(F_res2,F_mat1,F_mat1);
      F_mpz_mat_to_mpz_mat(res2, F_res2);
          
      result = mpz_mat_equal(res1, res2); 
      if (!result) 
      {
         printf("Error: bits1 = %ld, bits2 = %ld, bits3 = %ld, count1 = %ld\n", bits1, bits2, bits3, count1);
      }
          
      F_mpz_mat_clear(F_mat1);
      F_mpz_mat_clear(F_temp);
      F_mpz_mat_clear(F_res1);
      F_mpz_mat_clear(F_res2);
      mpz_mat_clear(m_mat1); 
      mpz_mat_clear(res1); 
      mpz_mat_clear(res2); 
   }
   return result;
}

int test__F_mpz_vec_scalar_submul_2exp_F_mpz()
{
   mpz_mat_t m_mat, m_mat2, m_mat3;
   F_mpz_mat_t F_mat, F_mat2;
   F_mpz_t mult;
   int result = 1;
   ulong bits, bits2, bits3, r, c, start, n, r1, r2, exp;
   mpz_t temp, temp2, m_mult;
   mpz_init(temp);
   mpz_init(temp2);
   mpz_init(m_mult);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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

      F_mpz_init(mult);
      
      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		mpz_mat_init(m_mat3, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      F_mpz_mat_init(F_mat2, r, c);
      
      mpz_randmat_dense(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      mpz_randmat_dense(m_mat2, r, c, bits2); 
      mpz_mat_to_F_mpz_mat(F_mat2, m_mat2);
          
      bits3 = z_randint(200);
		F_mpz_test_random(mult, bits3);
      F_mpz_get_mpz(m_mult, mult);

		_F_mpz_vec_scalar_submul_2exp_F_mpz(F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat3, F_mat2);
      
      ulong i;
      for (i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat2->entries[r2*c+i]);
			mpz_mul_2exp(temp2, m_mat->entries[r1*c+i], exp);
			mpz_submul(temp, temp2, m_mult);
         result &= (mpz_cmp(temp, m_mat3->entries[r2*c+i]) == 0);
      }
      
		if (!result) 
		{
			gmp_printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, m_mult = %Zd\n", r, c, bits, bits2, mult);
		}
      
      F_mpz_clear(mult);
      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		mpz_mat_clear(m_mat3); 
		F_mpz_mat_clear(F_mat);
      F_mpz_mat_clear(F_mat2);
   }
   
   // aliased multiply
   for (count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = z_randint(200)+ 1;
      r = z_randint(30)+1;  
      c = z_randint(30)+1;
		r1 = z_randint(r);
      r2 = z_randint(r);
      exp = z_randint(100);
		start = z_randint(c);
		n = z_randint(c-start)+1;

      F_mpz_init(mult);
      mpz_mat_init(m_mat, r, c); 
		mpz_mat_init(m_mat2, r, c); 
		F_mpz_mat_init(F_mat, r, c);
      
      mpz_randmat(m_mat, r, c, bits); 
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
          
      bits3 = z_randint(200);
		F_mpz_test_random(mult, bits3);
      F_mpz_get_mpz(m_mult, mult);

		_F_mpz_vec_scalar_submul_2exp_F_mpz(F_mat->rows[r2] + start, F_mat->rows[r1] + start, n, mult, exp);
      F_mpz_mat_to_mpz_mat(m_mat2, F_mat);

      ulong i;
      for (i = start; i < start + n; i++)
      {
         mpz_set(temp, m_mat->entries[r2*c+i]);
			mpz_mul_2exp(temp2, m_mat->entries[r1*c+i], exp);
			mpz_submul(temp, temp2, m_mult);
         result &= (mpz_cmp(temp, m_mat2->entries[r2*c+i]) == 0);
      }

		if (!result) 
		{
			gmp_printf("Error: r = %ld, c = %ld, bits = %ld, bits2 = %ld, mult = %Zd\n", r, c, bits, bits2, mult);
		}

      F_mpz_clear(mult);
      mpz_mat_clear(m_mat); 
		mpz_mat_clear(m_mat2); 
		F_mpz_mat_clear(F_mat);
   }
    
	mpz_clear(m_mult);
	mpz_clear(temp);
	mpz_clear(temp2);
	
   return result; 
}

int test__F_mpz_vec_scalar_mul()
{
   mpz_mat_t m_mat, m_mat2;
   F_mpz_mat_t F_mat, F_mat2;
   F_mpz_t F_sp;
   int result = 1;
   ulong bits, bits2, r, c, start, n, r1, r2;
   mpz_t m_sp1, m_sp2;
   mpz_init(m_sp1);
   mpz_init(m_sp2);
   
   ulong count1;
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
          
		_F_mpz_vec_scalar_product(F_sp, F_mat2->rows[r2] + start, F_mat->rows[r1] + start, n);
      F_mpz_get_mpz(m_sp1, F_sp);
      
      mpz_set_ui(m_sp2, 0);
      ulong i;
      for (i = start; i < start + n; i++)
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
   for (count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
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
          
		_F_mpz_vec_scalar_product(F_sp, F_mat->rows[r1] + start, F_mat->rows[r1] + start, n);
      F_mpz_get_mpz(m_sp1, F_sp);
      
      mpz_set_ui(m_sp2, 0);
      ulong i;
      for (i = start; i < start + n; i++)
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

int test__F_mpz_vec_to_d_vec_2exp()
{
   double * d1, * d2;
   F_mpz_mat_t F_mat;
   int result = 1;
   ulong bits, r, c, i, j;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      r = z_randint(30) + 1;
      c = z_randint(30) + 2;
      bits = z_randint(200) + 1;
      
	  F_mpz_mat_init(F_mat, r, c);
      d1 = (double *) malloc(c*sizeof(double));
      d2 = (double *) malloc(c*sizeof(double));
      	  
	  F_mpz_randmat(F_mat, r, c, bits);
       
	  i = z_randint(r);
	  ulong l1 = _F_mpz_vec_to_d_vec_2exp(d1, F_mat->rows[i], c/2);
      ulong l2 = _F_mpz_vec_to_d_vec_2exp(d1 + c/2, F_mat->rows[i] + c/2, (c + 1)/2);
      ulong l3 = _F_mpz_vec_to_d_vec_2exp(d2, F_mat->rows[i], c);

	  if (l1 < l2)
	     for (j = 0; j < c/2; j++)
		    d1[j] = ldexp(d1[j], l1 - l2);
      
      if (l2 < l1)
	     for (j = c/2; j < c; j++)
		    d1[j] = ldexp(d1[j], l2 - l1);
      
      result = (l3 == FLINT_MAX(l1, l2) && _d_vec_equal(d1, d2, c, FLINT_D_BITS - 1)); 
		if (!result) 
		{
			printf("Error: bits = %ld, r = %ld, c = %ld\n", bits, r, c);
		}
          
      F_mpz_mat_clear(F_mat);
	  free(d1);
	  free(d2);
   }
   
   return result;
}

int test__F_mpz_vec_to_mpfr_vec()
{
   __mpfr_struct * d1, * d2;
   F_mpz_mat_t F_mat;
   int result = 1;
   ulong bits, r, c, i, j;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      r = z_randint(30) + 1;
      c = z_randint(30) + 2;
      bits = z_randint(200) + 1;
      
	  F_mpz_mat_init(F_mat, r, c);
      d1 = (__mpfr_struct *) malloc(c*sizeof(__mpfr_struct));
      d2 = (__mpfr_struct *) malloc(c*sizeof(__mpfr_struct));
      for (i = 0; i < c; i++)
	  {
		  mpfr_init2(d1 + i, 200);
		  mpfr_init2(d2 + i, 200);
	  }

	  F_mpz_randmat(F_mat, r, c, bits);
       
	  i = z_randint(r);
	  _F_mpz_vec_to_mpfr_vec(d1, F_mat->rows[i], c/2);
      _F_mpz_vec_to_mpfr_vec(d1 + c/2, F_mat->rows[i] + c/2, (c + 1)/2);
      _F_mpz_vec_to_mpfr_vec(d2, F_mat->rows[i], c);

	  result = (_mpfr_vec_equal(d1, d2, c)); 
		if (!result) 
		{
			printf("Error: bits = %ld, r = %ld, c = %ld\n", bits, r, c);
		}
          
      F_mpz_mat_clear(F_mat);
	  for (i = 0; i < c; i++)
	  {
		  mpfr_clear(d1 + i);
		  mpfr_clear(d2 + i);
	  }
      free(d1);
	  free(d2);
   }
   
   return result;
}

int test__F_mpz_vec_2exp_to_mpfr_vec()
{
   __mpfr_struct * d1, * d2;
   F_mpz_mat_t F_mat;
   int * cexp;
   int result = 1;
   ulong bits, r, c, i, j;
   
   ulong count1;
   for (count1 = 0; (count1 < 10000*ITER) && (result == 1) ; count1++)
   {
      r = z_randint(30) + 1;
      c = z_randint(30) + 2;
      bits = z_randint(200) + 1;
      
	  cexp = (int *) malloc(c*sizeof(int));
	  for (i = 0; i < c; i++)
	     cexp[i] = (int) z_randint(200);
	  F_mpz_mat_init(F_mat, r, c);
      d1 = (__mpfr_struct *) malloc(c*sizeof(__mpfr_struct));
      d2 = (__mpfr_struct *) malloc(c*sizeof(__mpfr_struct));
      for (i = 0; i < c; i++)
	  {
		  mpfr_init2(d1 + i, 200);
		  mpfr_init2(d2 + i, 200);
	  }

	  F_mpz_randmat(F_mat, r, c, bits);
       
	  i = z_randint(r);
	  _F_mpz_vec_2exp_to_mpfr_vec(d1, F_mat->rows[i], c/2, cexp);
      _F_mpz_vec_2exp_to_mpfr_vec(d1 + c/2, F_mat->rows[i] + c/2, (c + 1)/2, cexp + c/2);
      _F_mpz_vec_2exp_to_mpfr_vec(d2, F_mat->rows[i], c, cexp);

	  result = (_mpfr_vec_equal(d1, d2, c)); 
		if (!result) 
		{
			printf("Error: bits = %ld, r = %ld, c = %ld\n", bits, r, c);
		}
          
      F_mpz_mat_clear(F_mat);
	  for (i = 0; i < c; i++)
	  {
		  mpfr_clear(d1 + i);
		  mpfr_clear(d2 + i);
	  }
      free(cexp);
	  free(d1);
	  free(d2);
   }
   
   return result;
}

void F_mpz_mat_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(_F_mpz_vec_init_clear); 
   RUN_TEST(_F_mpz_vec_copy_equal); 
   RUN_TEST(_F_mpz_vec_neg); 
   RUN_TEST(_F_mpz_vec_add_sub); 
   RUN_TEST(_F_mpz_vec_scalar_mul_ui); 
   RUN_TEST(_F_mpz_vec_scalar_mul_si); 
   RUN_TEST(_F_mpz_vec_scalar_mul_F_mpz); 
   RUN_TEST(_F_mpz_vec_scalar_addmul_ui); 
   RUN_TEST(_F_mpz_vec_scalar_addmul_F_mpz); 
   RUN_TEST(_F_mpz_vec_scalar_submul_ui); 
   RUN_TEST(_F_mpz_vec_scalar_submul_F_mpz); 
   RUN_TEST(_F_mpz_vec_scalar_addmul_2exp_ui); 
   RUN_TEST(_F_mpz_vec_scalar_submul_2exp_ui); 
   RUN_TEST(_F_mpz_vec_scalar_submul_2exp_F_mpz); 
   RUN_TEST(_F_mpz_vec_scalar_mul); 
   RUN_TEST(_F_mpz_vec_to_d_vec_2exp); 
   RUN_TEST(_F_mpz_vec_to_mpfr_vec); 
   RUN_TEST(_F_mpz_vec_2exp_to_mpfr_vec); 
   RUN_TEST(F_mpz_mat_convert); 
   RUN_TEST(F_mpz_mat_set); 
   RUN_TEST(F_mpz_mat_equal);  
   RUN_TEST(F_mpz_mat_swap);  
   RUN_TEST(F_mpz_mat_swap_rows);  
   RUN_TEST(F_mpz_mat_resize);
   RUN_TEST(F_mpz_mat_tofromstring);
   RUN_TEST(F_mpz_mat_tofromstringpretty);
   RUN_TEST(F_mpz_mat_neg); 
   RUN_TEST(F_mpz_mat_add); 
   RUN_TEST(F_mpz_mat_sub); 
   RUN_TEST(F_mpz_mat_mul_classical);
   
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
