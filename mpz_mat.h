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

mpz_poly.h: Matrices over Z, implemented as an array of mpz_t's
            Not intended to be efficient

Copyright (C) 2008, William Hart

******************************************************************************/

#ifndef MPZ_MAT_H
#define MPZ_MAT_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory-manager.h"

typedef struct
{
   mpz_t* entries;
   unsigned long r;
   unsigned long c;
} mpz_mat_struct;


// mpz_poly_t allows reference-like semantics for mpz_poly_struct:
typedef mpz_mat_struct mpz_mat_t[1];

// ------------------------------------------------------
// Initialisation and memory management

void mpz_mat_init(mpz_mat_t mat, ulong r, ulong c);
void mpz_mat_clear(mpz_mat_t mat);

// ------------------------------------------------------
// Comparison

static inline
int mpz_mat_equal(mpz_mat_t mat1, mpz_mat_t mat2)
{
	for (long i = 0; i < mat1->r*mat1->c; i++)
	   if (mpz_cmp(mat1->entries[i], mat2->entries[i])) return 0;

	return 1;
}

// ------------------------------------------------------
// Addition and subtraction

void mpz_mat_add(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2);

void mpz_mat_sub(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2);

// ------------------------------------------------------
// Multiplication

void mpz_mat_mul_classical(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2);

// *************** end of file

#ifdef __cplusplus
 }
#endif
 
#endif
