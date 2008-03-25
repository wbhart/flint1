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

   fmpz_poly.h: polynomials over the "flat" multi-precision integer format

   Copyright (C) 2007, William Hart

*****************************************************************************/

#ifndef MPIR_FMPZ_POLY_H
#define MPIR_FMPZ_POLY_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "memory_manager.h"
#include "mpir.h"
#include "fmpz.h"
#include "mpz_poly.h"

typedef struct
{
   fmpz_t * coeffs;
   unsigned long alloc;
   unsigned long length;
} fmpz_poly_struct;

// fmpz_poly_t allows reference-like semantics for fmpz_poly_struct:
typedef fmpz_poly_struct fmpz_poly_t[1];

/* ==============================================================================

   Memory management

===============================================================================*/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, ulong alloc);

void fmpz_poly_realloc(fmpz_poly_t poly, ulong alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, ulong alloc);

void fmpz_poly_clear(fmpz_poly_t poly);

/* ==============================================================================

   Normalisation

===============================================================================*/

static inline
void _fmpz_poly_normalise(fmpz_poly_t poly)
{
   while (poly->length && fmpz_is_zero(poly->coeffs + poly->length - 1))
      poly->length--;
}

/* ==============================================================================

   Conversion

===============================================================================*/

void mpz_poly_to_fmpz_poly(fmpz_poly_t res, mpz_poly_t poly);

void fmpz_poly_to_mpz_poly(mpz_poly_t res, fmpz_poly_t poly);

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

void fmpz_poly_add(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2);


#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
