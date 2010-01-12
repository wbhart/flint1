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
/****************************************************************************

NTL-interface.h: Header file for NTL-interface.cpp

Copyright (C) 2007, William Hart

*****************************************************************************/

#ifndef FLINT_NTL_INT_H
#define FLINT_NTL_INT_H

#ifdef __cplusplus
 extern "C" {
#endif

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>

#include "flint.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"
#include "fmpz.h"
#include "fmpz_poly.h"

NTL_CLIENT

/*
   Returns the number of limbs taken up by an NTL ZZ
*/

unsigned long ZZ_limbs(const ZZ& z);

/* 
   Convert an NTL ZZ to an fmpz_t
   Assumes the fmpz_t has already been allocated to have sufficient space
*/

void ZZ_to_fmpz(fmpz_t output, const ZZ& z);

/* 
   Convert an NTL ZZ to an F_mpz_t
*/

void ZZ_to_F_mpz(F_mpz_t output, const ZZ& z);

/*
   Convert an fmpz_t to an NTL ZZ
*/

void fmpz_to_ZZ(ZZ& output, const fmpz_t z);

/*
   Convert an F_mpz_t to an NTL ZZ
*/

void F_mpz_to_ZZ(ZZ& output, const F_mpz_t z);

/*
   Convert an fmpz_poly_t to an NTL ZZX
*/

void fmpz_poly_to_ZZX(ZZX& output, const fmpz_poly_t poly);

/*
   Convert an NTL ZZX to an fmpz_poly_t
*/

void ZZX_to_fmpz_poly(fmpz_poly_t output, const ZZX& poly);

/*
   Convert an NTL mat_ZZ to an F_mpz_mat_t
*/


/*
   Convert an F_mpz_mat_t to an NTL mat_ZZ
*/


#ifdef __cplusplus
 }
#endif

#endif
