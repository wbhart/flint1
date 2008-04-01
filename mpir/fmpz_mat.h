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

   fmpz_mat.h: Matrices over the "flat" multi-precision integer format

   Copyright (C) 2007, William Hart

*****************************************************************************/

#ifndef MPIR_FMPZ_MAT_H
#define MPIR_FMPZ_MAT_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "memory_manager.h"
#include "mpir.h"
#include "fmpz.h"

/*
   The basic matrix type
   If stride = cols, the total number of entries allocated is rows*cols
   The stride parameter allows for block submatrices of another matrix
   to be defined: one sets entries to an entry within an existing matrix,
   rows and cols to the new number of rows and columns of the block 
   starting at that entry, and stride to the old number of columns
*/

typedef struct
{
   fmpz_t * entries;
   ulong rows;
   ulong cols;
   ulong stride;
} fmpz_mat_struct;

// fmpz_mat_t allows reference-like semantics for fmpz_mat_struct:
typedef fmpz_mat_struct fmpz_mat_t[1];

/* ==============================================================================

   Memory management

===============================================================================*/

void fmpz_mat_init(fmpz_mat_t mat, ulong r, ulong c);

void fmpz_mat_clear(fmpz_mat_t mat);

/* ==============================================================================

   Get/set 

===============================================================================*/

void fmpz_mat_set_entry_ui(fmpz_mat_t mat, ulong r, ulong c, ulong x);

ulong fmpz_mat_get_entry_ui(fmpz_mat_t mat, ulong r, ulong c);

#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
