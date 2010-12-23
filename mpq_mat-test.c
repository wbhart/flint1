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

mpq_mat-test.c: Test code for mpq_mat.c and mpq_mat.h

Copyright (C) 2010, Andy Novocin, Max Flander, Bill Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "gmp.h"
#include <time.h>
#include "flint.h"
#include "long_extras.h"
#include "mpz_mat.h"
#include "mpq_mat.h"
#include "memory-manager.h"
#include "test-support.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"

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

int test_mpq_mat_inner_product()
{
   mpz_mat_t M;
   mpq_mat_t Mq;
   int result = 1;
   ulong bits, c;
   mpq_t res1, res2;

   int r = 3;   
   ulong count1;
   for (count1 = 0; (count1 < 2000*ITER) && (result == 1) ; count1++)
   {
      bits = z_randint(100);
      c = z_randint(30);
      
      mpz_mat_init(M, r, c); 
		mpq_mat_init(Mq, r, c);
      mpq_init(res1);
      mpq_init(res2);

		mpz_randmat(M, r, c, bits);
      mpz_mat_to_mpq_mat(Mq, M);

      mpq_mat_row_inner_product(res1, Mq, 0, Mq, 1);
      mpq_mat_row_inner_product(res2, Mq, 0, Mq, 2);
      mpq_add(res1, res1, res2);
      mpq_mat_row_add(Mq, 1, Mq, 2);
      mpq_mat_row_inner_product(res2, Mq, 0, Mq, 1);
    
      result = (mpq_cmp(res1, res2) == 0);
		if (!result) 
		{
			printf("Error: bits = %ld, r = %d, c = %ld\n", bits, r, c);
		}

      mpq_clear(res1);
      mpq_clear(res2);
      mpq_mat_clear(Mq);
		mpz_mat_clear(M);
   }
   
   return result;
}


int test_mpq_mat_GS()
{
   mpq_t temp;
   mpq_init(temp);
   mpz_mat_t M;
   mpq_mat_t Mq, mu, GS;
   int result = 1;

   ulong bits, r, c;

   ulong count1;
   for (count1 = 0; (count1 < 30*ITER) && (result == 1) ; count1++)
   {
      r = z_randint(20) + 1;
      c = z_randint(20) + 1;

      mpz_mat_init(M, r, c);
      mpq_mat_init(Mq, r, c);
      mpq_mat_init(mu, r, r);
      mpq_mat_init(GS, r, c);

      bits = z_randint(200) + 1;

      mpz_randmat(M, r, c, bits);

      mpz_mat_to_mpq_mat(Mq, M);

      mpq_mat_GS(mu, GS, Mq);

      int i,j;
      for (i = 1; i < Mq->r ; i++)
      {
         for (j = 0; j < i; j++)
         {


            mpq_mat_row_inner_product(temp, GS, i, GS, j);

            if (mpq_sgn(temp)!=0)
            {
               result=0;
               printf("Error: bits = %ld, r = %ld, c = %ld\n", bits, r, c);
            }
         }
      }

      mpz_mat_clear(M);
      mpq_mat_clear(Mq);
      mpq_mat_clear(mu);
      mpq_mat_clear(GS);
   }
   mpq_clear(temp);
   return result;
}

void mpq_mat_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %d\n", FLINT_BITS);

#if TESTFILE
#endif
   RUN_TEST(mpq_mat_inner_product); 
   RUN_TEST(mpq_mat_GS); 
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   mpq_mat_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}
