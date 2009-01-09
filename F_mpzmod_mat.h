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

   F_mpzmod_mat.h: Matrices over F_mpz's mod p, for p a multiprecision 
	                prime (FLINT 2.0).
					    Only used for test code at this point - not optimised.
   
   Copyright (C) 2008, William Hart
   
*****************************************************************************/

#ifndef FLINT_F_MPZMOD_MAT_H
#define FLINT_F_MPZMOD_MAT_H

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "flint.h"
#include "memory-manager.h"
#include "F_mpz.h"
#include "zmod_poly.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
   F_mpz * entries; // array containing the entries
	F_mpz ** rows; // array of pointers to the rows (so rows can be swapped easily)
   ulong r; // number of rows
   ulong c; // number of cols
   F_mpz p[1]; // modulus
} F_mpzmod_mat_struct;

typedef F_mpzmod_mat_struct F_mpzmod_mat_t[1];

/*******************************************************************************************

   Initialisation and memory management

*******************************************************************************************/

void F_mpzmod_mat_init(F_mpzmod_mat_t mat, F_mpz_t p, ulong rows, ulong cols);

void F_mpzmod_mat_clear(F_mpzmod_mat_t mat);

/*******************************************************************************************

   Arithmetic

*******************************************************************************************/

void F_mpzmod_mat_add(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1, F_mpzmod_mat_t mat2);

void F_mpzmod_mat_sub(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1, F_mpzmod_mat_t mat2);

void F_mpzmod_mat_neg(F_mpzmod_mat_t res, F_mpzmod_mat_t mat);

void F_mpzmod_mat_mul_classical(F_mpzmod_mat_t res, F_mpzmod_mat_t mat1, F_mpzmod_mat_t mat2);

/*******************************************************************************************

   Conversions

*******************************************************************************************/

void F_zmod_poly_to_mpzmod_mat_row(F_mpzmod_mat_t mat, ulong row, zmod_poly_t poly);

void F_zmod_poly_to_mpzmod_mat_col(F_mpzmod_mat_t mat, ulong col, zmod_poly_t poly);

void F_mpzmod_mat_col_to_zmod_poly(zmod_poly_t poly, F_mpzmod_mat_t mat, ulong col);

void F_mpzmod_mat_col_to_zmod_poly_shifted(zmod_poly_t poly, F_mpzmod_mat_t mat, ulong col, ulong * shift);

/*******************************************************************************************

   Set/get coefficients

*******************************************************************************************/

/*static inline
void F_mpzmod_mat_set_coeff_ui(F_mpzmod_mat_t mat, ulong row, ulong col, ulong val)
{
   mat->arr[row][col] = val;
}

static inline
ulong F_mpzmod_mat_get_coeff_ui(F_mpzmod_mat_t mat, ulong row, ulong col)
{
   return mat->arr[row][col];
}

static inline
ulong * F_mpzmod_mat_get_coeff_ptr(F_mpzmod_mat_t mat, ulong row, ulong col)
{
   return mat->arr[row] + col;
}
*/

/*******************************************************************************************

   Swap

*******************************************************************************************/

/*static inline
void F_mpzmod_mat_swap_rows(F_mpzmod_mat_t mat, ulong row1, ulong row2)
{
   ulong * temp = mat->arr[row1];
   mat->arr[row1] = mat->arr[row2];
   mat->arr[row2] = temp;
}*/

/*******************************************************************************************

   Elementary row operations

*******************************************************************************************/

void F_mpzmod_mat_row_scalar_addmul_right(F_mpzmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start);

void F_mpzmod_mat_row_scalar_submul_right(F_mpzmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start);

void F_mpzmod_mat_row_scalar_mul_right(F_mpzmod_mat_t mat, ulong row, ulong u, ulong start);

void F_zmod_vec_scalar_mul_range(ulong * r1, ulong * r2, ulong u, ulong p, double p_inv, ulong start, ulong end);

void F_zmod_vec_sub_range(ulong * r1, ulong * r2, ulong p, ulong start, ulong end);

/*******************************************************************************************

   Row reduction

*******************************************************************************************/

ulong F_mpzmod_mat_row_reduce_gauss(F_mpzmod_mat_t mat);

ulong F_mpzmod_mat_row_reduce_gauss_small(F_mpzmod_mat_t mat);

ulong F_mpzmod_mat_row_reduce_gauss_jordan(F_mpzmod_mat_t mat);

/*******************************************************************************************

   Input/output

*******************************************************************************************/

void F_mpzmod_mat_print(F_mpzmod_mat_t mat);

#ifdef __cplusplus
 }
#endif

#endif /* FLINT_F_mpzmod_mat_H */
