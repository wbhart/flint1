/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   fmpz_mat.c: Matrices over "flat" multiprecision integers

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "fmpz_mat.h"
#include "fmpz.h"
#include "mpir.h"
#include "memory_manager.h"

/* ==============================================================================

   Memory management

===============================================================================*/

void fmpz_mat_init(fmpz_mat_t mat, ulong r, ulong c)
{
   mat->entries = fmpz_init_array(r*c);
   mat->row_arr = (fmpz_t **) mpir_alloc(r*sizeof(fmpz_t **));
   for (ulong i = 0; i < r; i++) mat->row_arr[i] = mat->entries + c*i;
   mat->rows = r;
   mat->cols = c;
   mat->stride = c;
}

void fmpz_mat_clear(fmpz_mat_t mat)
{
   fmpz_clear_array(mat->entries, mat->rows * mat->cols);
   mpir_free(mat->row_arr);
}

/* ==============================================================================

   Get/set 

===============================================================================*/

/*
   Mat[r][c] = x
*/

void fmpz_mat_set_entry_ui(fmpz_mat_t mat, ulong r, ulong c, ulong x)
{
   fmpz_set_ui(mat->row_arr[r] + c, x);
}

/*
   Return Mat[r][c]
*/

ulong fmpz_mat_get_entry_ui(fmpz_mat_t mat, ulong r, ulong c)
{
   return fmpz_get_ui(mat->row_arr[r] + c);
}

/* ==============================================================================

   Addition/subtraction 

===============================================================================*/

/*
   r1 = r1 + r2
*/

void fmpz_mat_row_add(fmpz_mat_t mat, ulong r1, ulong r2, ulong n)
{
   fmpz_t * row1 = mat->row_arr[r1];
   fmpz_t * row2 = mat->row_arr[r2];
   for (ulong i = 0; i < n; i++)
   {
      fmpz_add(row1 + i, row1 + i, row2 + i);
   } 
}

/*
   r1 = r1 - r2
*/

void fmpz_mat_row_sub(fmpz_mat_t mat, ulong r1, ulong r2, ulong n)
{
   fmpz_t * row1 = mat->row_arr[r1];
   fmpz_t * row2 = mat->row_arr[r2];
   for (ulong i = 0; i < n; i++)
   {
      fmpz_sub(row1 + i, row1 + i, row2 + i);
   }
}

/* ==============================================================================

   Addmul/submul

===============================================================================*/

/*
   r1 = r1 + c*r2
*/

void fmpz_mat_row_addmul_ui(fmpz_mat_t mat, ulong r1, ulong r2, ulong c, ulong n)
{
   fmpz_t * row1 = mat->row_arr[r1];
   fmpz_t * row2 = mat->row_arr[r2];
   for (ulong i = 0; i < n; i++)
   {
      fmpz_addmul_ui(row1 + i, row2 + i, c);
   }
}

/*
   r1 = r1 - c*r2
*/

void fmpz_mat_row_submul_ui(fmpz_mat_t mat, ulong r1, ulong r2, ulong c, ulong n)
{
   fmpz_t * row1 = mat->row_arr[r1];
   fmpz_t * row2 = mat->row_arr[r2];
   for (ulong i = 0; i < n; i++)
   {
      fmpz_submul_ui(row1 + i, row2 + i, c);
   }
}

void fmpz_mat_row_addmul_2exp_ui(fmpz_mat_t mat, ulong r1, ulong r2, ulong c, ulong exp, ulong n, fmpz_t * temp)
{
   fmpz_t * row1 = mat->row_arr[r1];
   fmpz_t * row2 = mat->row_arr[r2];
   for (ulong i = 0; i < n; i++)
   {
      fmpz_mul_2exp(temp, row2 + i, exp);
      fmpz_addmul_ui(row1 + i, temp, c);
   }
}

void fmpz_mat_row_submul_2exp_ui(fmpz_mat_t mat, ulong r1, ulong r2, ulong c, ulong exp, ulong n, fmpz_t * temp)
{
   fmpz_t * row1 = mat->row_arr[r1];
   fmpz_t * row2 = mat->row_arr[r2];
   for (ulong i = 0; i < n; i++)
   {
      fmpz_mul_2exp(temp, row2 + i, exp);
      fmpz_submul_ui(row1 + i, temp, c);
   }
}

/* ==============================================================================

   Print

===============================================================================*/

void fmpz_mat_print(fmpz_mat_t mat, int rows, int cols) 
{
  ulong i, j; 

   printf("[");
   for (i = 0; i < rows; i++) 
   {
      fmpz_t * row = mat->row_arr[i];
      printf("[");
      for (j = 0; j < cols; j++) 
	  { 
	     fmpz_print(row + j); 
	     if (j < cols - 1) printf(" "); 
	  }
      if (i != rows - 1) printf("]\n"); 
   }  
   printf("]]\n"); 
}

// *************** end of file
