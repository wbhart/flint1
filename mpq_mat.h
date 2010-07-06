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


// *************** end of file

#ifdef __cplusplus
 }
#endif
 
#endif
