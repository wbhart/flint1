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

   F_mpz_mod_poly.c: Polynomials over F_mpz mod p.
   
   Copyright (C) 2009, William Hart, Andy Novocin

*****************************************************************************/

#include "F_mpz_mod_poly.h"
#include "long_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "flint.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void F_mpz_mod_poly_init(F_mpz_mod_poly_t poly, F_mpz_t P)
{
   poly->coeffs = (unsigned long*) flint_heap_alloc(1);
   poly->P = P;
   poly->alloc = 1;
   poly->length = 0;
}

void F_mpz_mod_poly_init2(F_mpz_mod_poly_t poly, F_mpz_t P, unsigned long alloc)
{
   FLINT_ASSERT(alloc >= 1);
   poly->coeffs = (unsigned long*) flint_heap_alloc(alloc);
   poly->P = P;   
   poly->alloc = alloc;
   poly->length = 0;
}

void F_mpz_mod_poly_clear(F_mpz_mod_poly_t poly)
{
   flint_heap_free(poly->coeffs);
}

