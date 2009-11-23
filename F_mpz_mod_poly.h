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

   F_mpz_mod_poly.h: Polynomials over F_mpz mod p.
   
   Copyright (C) 2009, William Hart, Andy Novocin

*****************************************************************************/

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "flint.h"
#include "F_mpz_poly.h"
#include "memory-manager.h"
#include "mpn_extras.h"
#include "long_extras.h"

#ifndef _F_MPZ_MOD_POLY_H_
#define _F_MPZ_MOD_POLY_H_

#ifdef __cplusplus
 extern "C" {
#endif


typedef struct
{
   F_mpz_t *coeffs;
   unsigned long alloc;
   unsigned long length;
   F_mpz_t P;
} F_mpz_mod_poly_struct;

typedef F_mpz_mod_poly_struct F_mpz_mod_poly_t[1];

void F_mpz_mod_poly_init(F_mpz_mod_poly_t poly, F_mpz_t P);

void F_mpz_mod_poly_init2(F_mpz_mod_poly_t poly, F_mpz_t P, unsigned long alloc);

void F_mpz_mod_poly_realloc(F_mpz_mod_poly_t poly, unsigned long alloc);

void __F_mpz_mod_poly_fit_length(F_mpz_mod_poly_t poly, unsigned long alloc);

void F_mpz_mod_poly_clear(F_mpz_mod_poly_t poly);

void __F_mpz_mod_poly_normalise(F_mpz_mod_poly_t poly);

#ifdef __cplusplus
 }
#endif

#endif /* _F_MPZ_MOD_POLY_H_ */
