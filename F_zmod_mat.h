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

   F_zmod_mat.h: Matrices over (unsigned) long mod p, for p prime with packed
	              representation (using packed_vec) (FLINT 2.0).
   
   Copyright (C) 2008, William Hart
   Copyright (C) 2008, Richard Howell-Peak

*****************************************************************************/

#ifndef FLINT_F_ZMOD_MAT_H
#define FLINT_F_ZMOD_MAT_H

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "flint.h"
#include "memory-manager.h"
#include "mpn_extras.h"
#include "long_extras.h"
#include "zmod_poly.h"
#include "packed_vec.h"
#include "F_mpzmod_mat.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
   pv_s arr; // structure containing the entries
	ulong * rows; // array of offsets of the rows (so rows can be swapped easily)
   ulong r; // number of rows
   ulong c; // number of cols
   ulong p; // modulus
   double p_inv; // precomputed inverse of p
} F_zmod_mat_struct;

typedef F_zmod_mat_struct F_zmod_mat_t[1];

/*******************************************************************************************

   Initialisation and memory management

*******************************************************************************************/

void F_zmod_mat_init(F_zmod_mat_t mat, ulong p, ulong rows, ulong cols);

void F_zmod_mat_init_precomp(F_zmod_mat_t mat, ulong p, double p_inv, 
						                              ulong rows, ulong cols);
void F_zmod_mat_clear(F_zmod_mat_t mat);

/*******************************************************************************************

   Conversions

*******************************************************************************************/

void F_zmod_mat_to_F_mpzmod_mat(F_mpzmod_mat_t res, F_zmod_mat_t mat);

void F_mpzmod_mat_to_F_zmod_mat(F_zmod_mat_t res, F_mpzmod_mat_t mat);

/*******************************************************************************************

   Arithmetic

*******************************************************************************************/

void F_zmod_mat_add(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2);

void F_zmod_mat_sub(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2);

void F_zmod_mat_neg(F_zmod_mat_t res, F_zmod_mat_t mat);

void F_zmod_mat_mul_classical(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2);

void F_zmod_mat_mul_strassen(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2);

/*******************************************************************************************

   Conversions

*******************************************************************************************/

void F_zmod_poly_to_zmod_mat_row(F_zmod_mat_t mat, ulong row, zmod_poly_t poly);

void F_zmod_poly_to_zmod_mat_col(F_zmod_mat_t mat, ulong col, zmod_poly_t poly);

void F_zmod_mat_col_to_zmod_poly(zmod_poly_t poly, F_zmod_mat_t mat, ulong col);

void F_zmod_mat_col_to_zmod_poly_shifted(zmod_poly_t poly, F_zmod_mat_t mat, ulong col, ulong * shift);

/*******************************************************************************************

   Set/get coefficients

*******************************************************************************************/

/*static inline
void F_zmod_mat_set_coeff_ui(F_zmod_mat_t mat, ulong row, ulong col, ulong val)
{
   mat->arr[row][col] = val;
}

static inline
ulong F_zmod_mat_get_coeff_ui(F_zmod_mat_t mat, ulong row, ulong col)
{
   return mat->arr[row][col];
}

static inline
ulong * F_zmod_mat_get_coeff_ptr(F_zmod_mat_t mat, ulong row, ulong col)
{
   return mat->arr[row] + col;
}
*/

/*******************************************************************************************

   Swap

*******************************************************************************************/

/*static inline
void F_zmod_mat_swap_rows(F_zmod_mat_t mat, ulong row1, ulong row2)
{
   ulong * temp = mat->arr[row1];
   mat->arr[row1] = mat->arr[row2];
   mat->arr[row2] = temp;
}*/

/*******************************************************************************************

   Elementary row operations

*******************************************************************************************/

void F_zmod_mat_row_scalar_addmul_right(F_zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start);

void F_zmod_mat_row_scalar_submul_right(F_zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start);

void F_zmod_mat_row_scalar_mul_right(F_zmod_mat_t mat, ulong row, ulong u, ulong start);

void F_zmod_vec_scalar_mul_range(ulong * r1, ulong * r2, ulong u, ulong p, double p_inv, ulong start, ulong end);

void F_zmod_vec_sub_range(ulong * r1, ulong * r2, ulong p, ulong start, ulong end);

/*******************************************************************************************

   Row reduction

*******************************************************************************************/

ulong F_zmod_mat_row_reduce_gauss(F_zmod_mat_t mat);

ulong F_zmod_mat_row_reduce_gauss_small(F_zmod_mat_t mat);

ulong F_zmod_mat_row_reduce_gauss_jordan(F_zmod_mat_t mat);

/*******************************************************************************************

   Input/output

*******************************************************************************************/

void F_zmod_mat_print(F_zmod_mat_t mat);

#ifdef __cplusplus
 }
#endif

#endif /* FLINT_F_ZMOD_MAT_H */
