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

mpz_poly.h: Polynomials over Z, implemented as an array of mpz_t's

Copyright (C) 2007, William Hart and David Harvey

******************************************************************************/

#ifndef MPZ_POLY_H
#define MPZ_POLY_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory_manager.h"

typedef struct
{
   mpz_t* coeffs;
   unsigned long alloc;
   unsigned long length;
} mpz_poly_struct;


// mpz_poly_t allows reference-like semantics for mpz_poly_struct:
typedef mpz_poly_struct mpz_poly_t[1];

// for some functions it's convenient to trick the compiler into letting us
// swap input argument pointers; we use mpz_poly_p for this
typedef mpz_poly_struct* mpz_poly_p;


#define SWAP_MPZ_POLY_PTRS(x, y)    \
do {                                \
   mpz_poly_p zzz_ptr = (x);        \
   (x) = (y);                       \
   (y) = zzz_ptr;                   \
} while (0);


// ------------------------------------------------------
// Initialisation and memory management

void mpz_poly_init(mpz_poly_t poly);

void mpz_poly_clear(mpz_poly_t poly);

void mpz_poly_init2(mpz_poly_t poly, unsigned long alloc);

void mpz_poly_init3(mpz_poly_t poly, unsigned long alloc, unsigned long bits);

void mpz_poly_realloc(mpz_poly_t poly, unsigned long alloc);

// _bits_ only applies to newly allocated coefficients, not existing ones...
void mpz_poly_realloc2(mpz_poly_t poly, unsigned long alloc,
                       unsigned long bits);


// this non-inlined version REQUIRES that alloc > poly->alloc
void __mpz_poly_ensure_alloc(mpz_poly_t poly, unsigned long alloc);

// this is arranged so that the initial comparison (very frequent) is inlined,
// but the actual allocation (infrequent) is not
static inline
void mpz_poly_ensure_alloc(mpz_poly_t poly, unsigned long alloc)
{
   if (alloc > poly->alloc)
      __mpz_poly_ensure_alloc(poly, alloc);
}



// ------------------------------------------------------
// Setting/retrieving coefficients

static inline
mpz_t* mpz_poly_coeff_ptr(mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return NULL;
   return &poly->coeffs[n];
}

static inline
void mpz_poly_get_coeff(mpz_t c, mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(c, 0);
   else
      mpz_set(c, poly->coeffs[n]);
}

static inline
unsigned long mpz_poly_get_coeff_ui(mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return 0;
   return mpz_get_ui(poly->coeffs[n]);
}

static inline
long mpz_poly_get_coeff_si(mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return 0;
   return mpz_get_si(poly->coeffs[n]);
}


void mpz_poly_set_coeff(mpz_poly_t poly, unsigned long n, mpz_t c);

void mpz_poly_set_coeff_ui(mpz_poly_t poly, unsigned long n, unsigned long c);

void mpz_poly_set_coeff_si(mpz_poly_t poly, unsigned long n, long c);

// ------------------------------------------------------
// Length and degree


void mpz_poly_normalise(mpz_poly_t poly);

// ------------------------------------------------------
// Comparison

int mpz_poly_equal(mpz_poly_t poly1, mpz_poly_t poly2);

//--------------------------------------------------------
// Addition/subtraction

void mpz_poly_add(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);

void mpz_poly_sub(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);

#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
