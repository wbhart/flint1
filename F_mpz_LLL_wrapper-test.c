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

F_mpz_LLL-wrapper_test.c: Test code for F_mpz_LLL_wrapper.c and F_mpz_LLL_wrapper.h

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
#include "F_mpz_LLL_wrapper.h"
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

void mpz_mat_randintrel(mpz_mat_t mat, ulong r, ulong c, ulong maxbits)
{
   ulong bits;
   mpz_t temp;
   mpz_init(temp);

   long i;
   for (i = 0; i < r; i++)
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
      mpz_set(mat->entries[i*c + c - 1], temp);
      long j;
      for (j = 0; j < i; j++)
         mpz_set_ui(mat->entries[i*c + j], 0);
      mpz_set_ui(mat->entries[i*c + i], 1);
      for (j = i + 1; j < c - 1; j++)
         mpz_set_ui(mat->entries[i*c + j], 0);
   }

   mpz_clear(temp);

}

void F_mpz_mat_randintrel(F_mpz_mat_t mat, ulong r, ulong c, ulong bits)
{
   mpz_mat_t m_mat;
   mpz_mat_init(m_mat, r, c);
   mpz_mat_randintrel(m_mat, r, c, bits);
   mpz_mat_to_F_mpz_mat(mat, m_mat);
   mpz_mat_clear(m_mat);
}

void _mpfr_vec_clean_scalar_product2(mpfr_t sp, __mpfr_struct * vec1, __mpfr_struct * vec2, int n, mp_prec_t prec)
{
  mpfr_t tmp;
  mpfr_init2(tmp, prec);
  
  mpfr_mul(sp, vec1, vec2, GMP_RNDN);
  
  for (long i = 1; i < n; i++)
  {
     mpfr_mul(tmp, vec1 + i, vec2 + i, GMP_RNDN);
     mpfr_add(sp, sp, tmp, GMP_RNDN);
  }

  mpfr_clear(tmp);

  return;
}

void F_mpz_mat_RQ_factor(F_mpz_mat_t B, __mpfr_struct ** R, __mpfr_struct ** Q, long r, long c, mp_prec_t prec){

//Doing modified GSO will convert Q from B to Q row by row
   mpfr_t tmp;
   mpfr_init2(tmp, prec);

   long i, k, j;
   for (i = 0; i < r; i++)
	   _F_mpz_vec_to_mpfr_vec(Q[i], B->rows[i], c); 

   //iteration k should start with Q = q1, ..., q_(k-1), a_k', ... a_r' then convert a_k to q_k by a subtraction and the other a_j' should be modded by a_k  
   for (k = 0; k < r; k++){
      _mpfr_vec_norm2(R[k]+k, Q[k], c, prec);
      mpfr_sqrt(R[k]+k, R[k]+k, GMP_RNDN);
      for (i = 0; i < c; i++)
         mpfr_div(Q[k]+i, Q[k]+i, R[k]+k, GMP_RNDN);
      for (j = k+1; j < r; j++){
         _mpfr_vec_clean_scalar_product2(R[j]+k, Q[k], Q[j], c, prec);
         for (i = 0; i < c; i++){
            mpfr_mul(tmp, R[j]+k, Q[k]+i, GMP_RNDN);
            mpfr_sub(Q[j]+i, Q[j]+i, tmp, GMP_RNDN);
         }
      }
   }

   mpfr_clear(tmp);
   return;
}


int test_F_mpz_LLL_randintrel()
{
   mpz_mat_t m_mat;
   F_mpz_mat_t F_mat;
   int result = 1;
   ulong bits;
   F_mpz_t fzero;
   F_mpz_init(fzero);
   
   ulong count1;
   for (count1 = 0; (count1 < 10*ITER) && (result == 1) ; count1++)
   {
      ulong r = z_randint(10)+1;
      ulong c = r + 1;

      F_mpz_mat_init(F_mat, r, c);

      mpz_mat_init(m_mat, r, c);

      bits = z_randint(200) + 1;
      
      mpz_mat_randintrel(m_mat, r, c, bits);
           
      mpz_mat_to_F_mpz_mat(F_mat, m_mat);
      
      F_mpz_set_d_2exp(fzero, 2.0, bits);
// good stuff here

      printf("bits = %ld\n", bits); 
      F_mpz_mat_print_pretty(F_mat);

      knapsack_LLL_wrapper_with_removal(F_mat, fzero);

      F_mpz_mat_print_pretty(F_mat);

      mp_prec_t prec;
      prec = 10;

      __mpfr_struct ** Q, ** R;

      Q = mpfr_mat_init2(r, c, prec);
      R = mpfr_mat_init2(r, r, prec);

      F_mpz_mat_RQ_factor(F_mat, R, Q, r, c, prec); 

// should be that RQ = FM_copy
      long j;
      for (j = 0; j < r; j++)
      { 
         mpfr_printf("%.12Rf was R[i][i] for i = %ld\n", R[j] + j, j); 
      }

      mpfr_mat_clear(Q, r, c);
      mpfr_mat_clear(R, r, r);
          
//result here       result = mpz_mat_equal(res1, res2); 
      result = 1; 
      if (!result) 
      {
         printf("Error: bits = %ld, count1 = %ld\n", bits, count1);
      }
          
      F_mpz_mat_clear(F_mat);
      mpz_mat_clear(m_mat);
   }

   F_mpz_clear(fzero);
   return result;
}


void F_mpz_mat_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %d\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(F_mpz_LLL_randintrel); 

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
