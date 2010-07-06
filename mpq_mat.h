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

mpq_mat.h: Matrices over Q, implemented as an array of mpq_t's
            Not intended to be efficient

Copyright (C) 2010 Andy Novocin, Max Flander

******************************************************************************/

#ifndef MPQ_MAT_H
#define MPQ_MAT_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory-manager.h"

typedef struct
{
   mpq_t* entries;
   unsigned long r;
   unsigned long c;
} mpq_mat_struct;


// mpz_poly_t allows reference-like semantics for mpz_poly_struct:
typedef mpq_mat_struct mpq_mat_t[1];

// ------------------------------------------------------
// Initialisation and memory management

void mpq_mat_init(mpq_mat_t mat, ulong r, ulong c);
void mpq_mat_clear(mpq_mat_t mat);

//------------------------------------------------------
//Conversion

void mpz_mat_to_mpq_mat(mpq_mat_t res, mpz_mat_t mat);

// ----------------------------------------------------
// I/O 

void mpq_mat_print(mpq_mat_t mat);

void mpq_mat_print_pretty(mpq_mat_t mat);

// ----------------------------------------------
//Gram Schmidt Orthogonalization functions

void mpq_mat_row_scalar_product(mpq_mat_t mat1, ulong r1, mpq_mat_t mat2, ulong r2, mpq_t scalar);

void mpq_mat_row_inner_product(mpq_t res, mpq_mat_t mat1, ulong r1, mpq_mat_t mat2, ulong r2);

void mpq_mat_row_add(mpq_mat_t mat1, ulong r1, mpq_mat_t mat2, ulong r2);

// will set mu[i, j] to gso projection for i < j and GS[i] to i'th G-S vector
void mpq_mat_GS(mpq_mat_t mu, mpq_mat_t GS, mpq_mat_t mat);

// returns 1 if delta-eta reduced 0 otherwise
int mpq_mat_is_reduced(mpq_mat_t mu, mpq_mat_t GS, double delta, double eta);
// *************** end of file

#ifdef __cplusplus
 }
#endif
 
#endif
