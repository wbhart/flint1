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
/*****************************************************************************

   zmod_mat.c: Matrices over F_mpz's mod p, for p a multiprecision 
	                prime (FLINT 2.0).
					    Only used for test code at this point - not optimised.

   Copyright (C) 2008, William Hart.
   
*****************************************************************************/

#include "F_mpzmod_mat.h"
#include "flint.h"
#include "F_mpz.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void F_mpzmod_mat_init(F_mpzmod_mat_t mat, F_mpz_t p, ulong rows, ulong cols)
{
   mat->entries = (F_mpz *) flint_heap_alloc(rows*cols);
	mat->rows = (F_mpz **) flint_heap_alloc(rows);
   
   // Set up the rows
   ulong i;
   for (i = 0; i < rows; i++)
      mat->rows[i] = mat->entries + i*cols;

	ulong i;
	for (i = 0; i < rows*cols; i++)
		F_mpz_init(mat->entries + i);

   F_mpz_set(mat->p, p);
	mat->r = rows;
   mat->c = cols;
}

void F_mpzmod_mat_clear(F_mpzmod_mat_t mat)
{
   ulong i;
   for (i = 0; i < mat->r*mat->c; i++)
		F_mpz_clear(mat->entries + i);
	flint_heap_free(mat->rows);
   flint_heap_free(mat->entries);
}

/*******************************************************************************************

   Arithmetic

*******************************************************************************************/

void F_mpzmod_mat_add(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1, F_mpzmod_mat_t mat2)
{
	ulong i;
	for (i = 0; i < mat1->r; i++)
	{
		ulong j;
		for (j = 0; j < mat1->c; j++)
		{
			F_mpz_add(res->rows[i] + j, mat1->rows[i] + j, mat2->rows[i] + j);
			if (F_mpz_cmpabs(res->rows[i] + j, mat1->p) >= 0)
				F_mpz_sub(res->rows[i] + j, res->rows[i] + j, mat1->p);
		}
	}
}

void F_mpzmod_mat_sub(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1, F_mpzmod_mat_t mat2)
{
	ulong i;
	for (i = 0; i < mat1->r; i++)
	{
		ulong j;
		for (j = 0; j < mat1->c; j++)
		{
			F_mpz_sub(res->rows[i] + j, mat1->rows[i] + j, mat2->rows[i] + j);
			if (F_mpz_sgn(res->rows[i] + j) < 0)
				F_mpz_add(res->rows[i] + j, res->rows[i] + j, mat1->p);
		}
	}
}

void F_mpzmod_mat_neg(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1)
{
	ulong i;
	for (i = 0; i < mat1->r; i++)
	{
		ulong j;
		for (j = 0; j < mat1->c; j++)
		{
			if (F_mpz_sgn(mat1->rows[i] + j))
			   F_mpz_sub(res->rows[i] + j, mat1->p, mat1->rows[i] + j);
			else F_mpz_set_ui(res->rows[i] + j, 0L);
		}
	}
}

void F_mpzmod_mat_set(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1)
{
	ulong i;
	for (i = 0; i < mat1->r; i++)
	{
		ulong j;
		for (j = 0; j < mat1->c; j++)
			F_mpz_set(res->rows[i] + j, mat1->rows[i] + j);
	}
}

void F_mpzmod_mat_mul_classical(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1, F_mpzmod_mat_t mat2)
{
   ulong c1 = mat1->c;
	ulong r2 = mat2->r;
   
	if ((c1 != r2) || (c1 == 0))
	{
		printf("FLINT exception : invalid matrix multiplication!\n");
		abort();
	}

	ulong r1 = mat1->r;
   ulong c2 = mat2->c;

	if ((r1 == 0) || (c2 == 0)) return; // no work to do

	F_mpz * temp = (F_mpz *) flint_heap_alloc(c2);

	ulong i;
	for (i = 0; i < c2; i++)
		F_mpz_init(temp + i);
	
   ulong i;
   for (i = 0; i < r1; i++) // for each row of mat1
	{
		F_mpz * c = mat1->rows[i];

		ulong k;
		for (k = 0; k < c2; k++) // do initial scalar product of row 1 of mat2 by c
		   F_mpz_mul2(temp + k, mat2->rows[0] + k, c);

		ulong j;
		for (j = 1; j < c1; j++) // compute scalar product for rows 1 to c1 of mat2
		{
         ulong k;
         for (k = 0; k < c2; k++) // do scalar product of row j of mat2 by c
		   {
            F_mpz_addmul(temp + k, mat2->rows[j] + k, c + j);
         }
		}
			   
		ulong k;
		for (k = 0; k < c2; k++) // do reduction mod p and store in result
	   {
		   F_mpz_mod(res->rows[i] + k, temp + k, mat1->p);
		}
	}

	ulong i;
	for (i = 0; i < c2; i++)
		F_mpz_clear(temp + i);

	flint_heap_free(temp);
}

