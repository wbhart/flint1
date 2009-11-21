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

   zmod_mat.h: Matrices over (unsigned) long mod p, for p prime.
   
   Copyright (C) 2008, 2009 William Hart
   Copyright (C) 2008, Richard Howell-Peak

*****************************************************************************/

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "flint.h"
#include "memory-manager.h"
#include "mpn_extras.h"
#include "long_extras.h"
#include "zmod_poly.h"

#ifndef _ZMOD_MAT_H_
#define _ZMOD_MAT_H_

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
   ulong ** arr; // array of pointers to the rows (rows can be swapped easily)
   ulong * ptr; // actual data (only used for allocation and deallocation of memory)
   ulong rows; // number of rows
   ulong cols; // number of cols
   ulong p; // modulus
   double p_inv; // precomputed inverse of p
} zmod_mat_struct;

typedef zmod_mat_struct zmod_mat_t[1];

/*******************************************************************************************

   Initialisation and memory management

*******************************************************************************************/

void zmod_mat_init(zmod_mat_t mat, ulong p, ulong rows, ulong cols);

void zmod_mat_init_precomp(zmod_mat_t mat, ulong p, double p_inv, 
						                              ulong rows, ulong cols);
void zmod_mat_clear(zmod_mat_t mat);

/*******************************************************************************************

   Conversions

*******************************************************************************************/

void zmod_poly_to_zmod_mat_row(zmod_mat_t mat, ulong row, zmod_poly_t poly);

void zmod_poly_to_zmod_mat_col(zmod_mat_t mat, ulong col, zmod_poly_t poly);

void zmod_mat_col_to_zmod_poly(zmod_poly_t poly, zmod_mat_t mat, ulong col);

void zmod_mat_col_to_zmod_poly_shifted(zmod_poly_t poly, zmod_mat_t mat, ulong col, ulong * shift);

/*******************************************************************************************

   Set/get coefficients

*******************************************************************************************/

static inline
void zmod_mat_set_coeff_ui(zmod_mat_t mat, ulong row, ulong col, ulong val)
{
   mat->arr[row][col] = val;
}

static inline
ulong zmod_mat_get_coeff_ui(zmod_mat_t mat, ulong row, ulong col)
{
   return mat->arr[row][col];
}

static inline
ulong * zmod_mat_get_coeff_ptr(zmod_mat_t mat, ulong row, ulong col)
{
   return mat->arr[row] + col;
}


/*******************************************************************************************

   Swap

*******************************************************************************************/

static inline
void zmod_mat_swap_rows(zmod_mat_t mat, ulong row1, ulong row2)
{
   ulong * temp = mat->arr[row1];
   mat->arr[row1] = mat->arr[row2];
   mat->arr[row2] = temp;
}

/*******************************************************************************************

   Elementary row operations

*******************************************************************************************/

void zmod_mat_row_scalar_addmul_right(zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start);

void zmod_mat_row_scalar_submul_right(zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start);

void zmod_mat_row_scalar_mul_right(zmod_mat_t mat, ulong row, ulong u, ulong start);

void zmod_vec_scalar_mul_range(ulong * r1, ulong * r2, ulong u, ulong p, double p_inv, ulong start, ulong end);

void zmod_vec_sub_range(ulong * r1, ulong * r2, ulong p, ulong start, ulong end);

/*******************************************************************************************

   Row reduction

*******************************************************************************************/

ulong zmod_mat_row_reduce_gauss(zmod_mat_t mat);

ulong zmod_mat_row_reduce_gauss_small(zmod_mat_t mat);

ulong zmod_mat_row_reduce_gauss_jordan(zmod_mat_t mat);

/*******************************************************************************************

   Input/output

*******************************************************************************************/

void zmod_mat_print(zmod_mat_t mat);

/*******************************************************************************************

   Input/output

*******************************************************************************************/

ulong zmod_mat_scalar_mul(ulong * r, ulong ** arr, ulong c, ulong n, ulong p, double p_inv)

void zmod_mat_mul_classical(zmod_mat_t prod, zmod_mat_t A, zmod_mat_t B);

#ifdef __cplusplus
 }
#endif

#endif /* _ZMOD_MAT_H_ */
