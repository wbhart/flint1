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
   F_mpz * coeffs;
   ulong alloc;
   ulong length;
   F_mpz_t P;
} F_mpz_mod_poly_struct;

typedef F_mpz_mod_poly_struct F_mpz_mod_poly_t[1];

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void F_mpz_mod_poly_init(F_mpz_mod_poly_t poly, const F_mpz_t P);

void F_mpz_mod_poly_init2(F_mpz_mod_poly_t poly, const F_mpz_t P, const ulong alloc);

void F_mpz_mod_poly_clear(F_mpz_mod_poly_t poly);

void F_mpz_mod_poly_realloc(F_mpz_mod_poly_t poly, const ulong alloc);

void F_mpz_mod_poly_fit_length(F_mpz_mod_poly_t poly, const ulong length);

/****************************************************************************

   Normalisation/truncation

****************************************************************************/

void _F_mpz_mod_poly_normalise(F_mpz_mod_poly_t poly);

static inline
void _F_mpz_mod_poly_set_length(F_mpz_mod_poly_t poly, const ulong length)
{
	if (poly->length > length) // demote coefficients beyond new length
   {
      for (ulong i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);	
   } 

	poly->length = length;
}

static inline
void F_mpz_mod_poly_truncate(F_mpz_mod_poly_t poly, const ulong length)
{
	if (poly->length > length) // only truncate if necessary
   {
      for (ulong i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);
		poly->length = length;
      _F_mpz_mod_poly_normalise(poly);
   }  
}

/****************************************************************************

   Zero

****************************************************************************/

static inline 
void F_mpz_mod_poly_zero(F_mpz_mod_poly_t poly)
{
   _F_mpz_mod_poly_set_length(poly, 0);
}

/****************************************************************************

   Conversion

****************************************************************************/

void mpz_poly_to_F_mpz_mod_poly(F_mpz_mod_poly_t F_poly, const mpz_poly_t m_poly);

void F_mpz_mod_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_mod_poly_t F_poly);

void zmod_poly_to_F_mpz_mod_poly(F_mpz_mod_poly_t fpol, const zmod_poly_t zpol);

/****************************************************************************

   Attach

****************************************************************************/

/* 
   Attach an F_mpz_mod_poly to an F_mpz_poly. Allows one to temporarily treat
   an F_mpz_mod_poly as though it were an F_mpz_poly.
*/
static inline
void _F_mpz_poly_attach_F_mpz_mod_poly(F_mpz_poly_t out, const F_mpz_mod_poly_t in)
{
   out->coeffs = in->coeffs;
   out->length = in->length;
   out->alloc = in->alloc;
}

/* 
   Attach an F_mpz_poly to an F_mpz_mod_poly. Allows one to return an 
   F_mpz_mod_poly, which is temporarily being treated as an F_mpz_poly,
   to its original F_mpz_mod_poly object.
*/
static inline
void _F_mpz_mod_poly_attach_F_mpz_poly(F_mpz_mod_poly_t out, const F_mpz_poly_t in)
{
   out->coeffs = in->coeffs;
   out->length = in->length;
   out->alloc = in->alloc;
}


/****************************************************************************

   Reduce mod P

****************************************************************************/

void _F_mpz_mod_poly_reduce_coeffs(F_mpz_mod_poly_t poly)
{
   for (ulong i = 0; i < poly->length; i++)
      F_mpz_mod(poly->coeffs + i, poly->coeffs + i, poly->P);
   _F_mpz_mod_poly_normalise(poly);
}

void _F_mpz_poly_reduce_coeffs(F_mpz_poly_t poly, const F_mpz_t P)
{
   for (ulong i = 0; i < poly->length; i++)
      F_mpz_mod(poly->coeffs + i, poly->coeffs + i, P);
}

/****************************************************************************

   Assignment/swap

****************************************************************************/

void F_mpz_mod_poly_swap(F_mpz_mod_poly_t poly1, F_mpz_mod_poly_t poly2);

/****************************************************************************

   Multiplication

****************************************************************************/

void _F_mpz_mod_poly_mul(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2);

void F_mpz_mod_poly_mul(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2);

#ifdef __cplusplus
 }
#endif

#endif /* _F_MPZ_MOD_POLY_H_ */

