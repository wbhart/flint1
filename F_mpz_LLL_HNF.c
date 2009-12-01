/*
    Copyright 2009 William Hart

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
*/

/****************************************************************************

   F_mpz_LLL_HNF.c: Implements the Hermite Normal Form algorithm of 
	                 Havas, Majewski and Matthews. See :
						  http://www.numbertheory.org/pdfs/gcd.pdf
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "gmp.h"
#include "flint.h"
#include "F_mpz_mat.h"
#include "d_mat.h"

/*
   Negate row j before the diagonal
	Negate column j below the diagonal
*/
void F_mpz_mat_HNF_LLL_minus(F_mpz_mat_t L, ulong j)
{
   long i;
   for (i = 0; i < j - 1; i++)
	   F_mpz_neg(L->rows[j] + i, L->rows[j] + i);
	ulong r;
	for (r = j + 1; r < L->r; r++)
		F_mpz_neg(L->rows[r] + j, L->rows[r] + j);
}

void F_mpz_mat_HNF_LLL_swap(F_mpz_mat_t A, F_mpz_mat_t B, F_mpz_mat_t L, F_mpz * D, ulong k)
{
	ulong m = A->r;
	ulong n = A->c;
	
	F_mpz_mat_row_swap(A, k, A, k - 1, 0, n);
   F_mpz_mat_row_swap(B, k, B, k - 1, 0, n);
   
	ulong j;
	for (j = 0; j < k - 1; j++) F_mpz_swap(L->rows[k] + j, L->rows[k - 1] + j);

	F_mpz_t t;
	F_mpz_init(t);

	ulong i;
	for (i = k + 1; i < m; i++)
	{
		F_mpz_mul2(L->rows[i] + k - 1, L->rows[i] + k - 1, L->rows[k] + k - 1);
		F_mpz_addmul(L->rows[i] + k - 1, L->rows[i] + k, D + k - 2);
		F_mpz_divexact(L->rows[i] + k - 1, L->rows[i] + k - 1, D + k - 1);

		F_mpz_mul2(L->rows[i] + k, L->rows[i] + k, L->rows[k] + k - 1);
		F_mpz_submul(L->rows[i] + k, L->rows[i] + k - 1, D + k);
		F_mpz_neg(L->rows[i] + k, L->rows[i] + k);
		F_mpz_divexact(L->rows[i] + k, L->rows[i] + k, D + k - 1);

		F_mpz_mul2(t, D + k - 2, D + k);
		F_mpz_addmul(t, L->rows[k] + k, L->rows[k] + k - 1);
		F_mpz_divexact(D + k - 1, t, D + k - 1);
	}

	F_mpz_clear(t);
}

void F_mpz_mat_HNF_LLL_reduce(F_mpz_mat_t A, F_mpz_mat_t B, F_mpz_mat_t L, 
										  long * col1, long * col2, F_mpz * D, ulong k, ulong i)
{
   ulong n = A->c;
	
	for (*col1 = 0 ; *col1 < n; (*col1)++)
	{
		if (!F_mpz_is_zero(A->rows[i] + (*col1))) break;
	}

	if (F_mpz_sgn(A->rows[i] + (*col1)) < 0)
	{
      F_mpz_mat_HNF_LLL_minus(L, i);
		F_mpz_mat_row_neg(B, i, B, i, 0, n);
	} else
	   (*col1) = n;

	for (*col2 = 0 ; *col2 < n; (*col2)++)
	{
		if (!F_mpz_is_zero(A->rows[k] + (*col2))) break;
	}

   if (F_mpz_sgn(A->rows[k] + (*col2)) < 0)
	{
      F_mpz_mat_HNF_LLL_minus(L, k);
		F_mpz_mat_row_neg(B, k, B, k, 0, n);
	} else
	   (*col2) = n;


	F_mpz_t q;
	F_mpz_init(q);

   if ((*col1) < n)
	{
		F_mpz_fdiv_q(q, A->rows[k] + (*col1), A->rows[i] + (*col1));
	} else
	{
		F_mpz_mul_2exp(q, L->rows[k] + i, 1);
		if (F_mpz_cmpabs(q, D + i) > 0)
			F_mpz_rdiv_q(q, L->rows[k] + i, D + i);
		else
			F_mpz_zero(q);
	}

	if (!F_mpz_is_zero(q))
	{
		F_mpz_mat_row_submul(A, k, A, i, 0, n, q);
      F_mpz_mat_row_submul(B, k, B, i, 0, n, q);
      F_mpz_submul(L->rows[k] + i, q, D + i);
		ulong j;
		for (j = 0; j < i - 1; j++)
		{
         F_mpz_submul(L->rows[k] + j, L->rows[i] + j, q);
		}
	}

	F_mpz_clear(q);
}

/* 
   Computes the Hermite normal form A of G, including the 
	transformation matrix B
*/
void F_mpz_mat_HNF_LLL(F_mpz_mat_t A, F_mpz_mat_t B, F_mpz_mat_t G)
{
   ulong m = G->r;
	ulong n = G->c;
	ulong m1 = 3;
	ulong n1 = 4; /* alpha = m1/n1 */
	ulong k = 1;
	ulong col1, col2;
	F_mpz_mat_t L;

	F_mpz * D = (F_mpz *) flint_heap_alloc(m*sizeof(F_mpz));
   ulong i;
   for (i = 0; i < m; i++) F_mpz_init(D + i);

	F_mpz_mat_init_identity(B, m);
	F_mpz_mat_init(L, m, m);
	F_mpz_mat_set(A, G);

	F_mpz_t temp1, temp2;
	F_mpz_init(temp1);
   F_mpz_init(temp2);

	while (k < m)
	{
      F_mpz_mat_HNF_LLL_reduce(A, B, L, &col1, &col2, D, k, k - 1);
		if ((col1 <= col2) && (col1 < n))
		{
			F_mpz_mat_HNF_LLL_swap(A, B, L, D, k);
			if (k > 1) k--;
			continue;
		} else if ((col1 == n) && (col2 == n))
		{
			F_mpz_mul2(temp1, D + k - 2, D + k);
			F_mpz_addmul(temp1, L->rows[k] + k - 1, L->rows[k] + k - 1);
			F_mpz_mul_ui(temp1, temp1, n1);
			
			F_mpz_mul2(temp2, D + k - 1, D + k - 1);
			F_mpz_mul_ui(temp2, temp2, m1);

			if (F_mpz_cmpabs(temp1, temp2) < 0)
			{
			   F_mpz_mat_HNF_LLL_swap(A, B, L, D, k);
			   if (k > 1) k--;
				continue;
			}
		} 
      
		long i;
		for (i = k - 2; i >= 0; i--)
			F_mpz_mat_HNF_LLL_reduce(A, B, L, &col1, &col2, D, k, i);
		k++;
	}

	F_mpz_clear(temp1);
	F_mpz_clear(temp2);
	ulong i;
	for (i = 0; i < m; i++) F_mpz_clear(D + i);
   flint_heap_free(D);
}
