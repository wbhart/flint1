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
/******************************************************************************

mpz_mat.c: Matrices over Z, implemented as an array of mpz_t's
           Not intended to be efficient

Copyright (C) 2008, William Hart

******************************************************************************/

#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "mpz_mat.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpz_mat_init(mpz_mat_t mat, ulong r, ulong c)
{
   if ((r) && (c)) mat->entries = (mpz_t*) flint_heap_alloc_bytes(r*c*sizeof(mpz_t));
	else mat->entries = NULL;

   for (long i = 0; i < r*c; i++)
	{
	   mpz_init(mat->entries[i]);
	}
	mat->r = r;
	mat->c = c;
}

void mpz_mat_clear(mpz_mat_t mat)
{
   for (long i = 0; i < mat->r*mat->c; i++)
      mpz_clear(mat->entries[i]);

   if (mat->entries) flint_heap_free(mat->entries);
	mat->entries = NULL;
	mat->r = 0;
	mat->c = 0;
}

int mpz_mat_add(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2)
{
	for (long i = 0; i < mat1->r*mat1->c; i++)
	   mpz_add(res->entries[i], mat1->entries[i], mat2->entries[i]);
}

int mpz_mat_sub(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2)
{
	for (long i = 0; i < mat1->r*mat1->c; i++)
	   mpz_sub(res->entries[i], mat1->entries[i], mat2->entries[i]);
}

// *************** end of file
