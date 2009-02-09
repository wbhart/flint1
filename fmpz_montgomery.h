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

   fmpz_montgomery.h: "flat" integer format functions using Montgomery format

   Copyright (C) 2008, William Hart and Gonzalo Tornaria

*****************************************************************************/

#ifndef FMPZ_MONTGOMERY_H
#define FMPZ_MONTGOMERY_H

typedef struct
{
    fmpz_t M; // Holds preconditioned modulus M
    fmpz_t R; // Holds 1/M  mod  2^#
    fmpz_t t; // Holds 2^#
    fmpz_t b; // Holds extra value for mulmod, divmod, mod
} fmpz_montgomery_t[1];


void fmpz_montgomery_init(fmpz_montgomery_t mont, fmpz_t m);

void fmpz_montgomery_clear(fmpz_montgomery_t mont);

void fmpz_montgomery_redc(fmpz_t res, fmpz_t x, fmpz_montgomery_t mont);

void fmpz_montgomery_mulmod_init(fmpz_montgomery_t mont, fmpz_t b, fmpz_t m);

void fmpz_montgomery_divmod_init(fmpz_montgomery_t mont, fmpz_t b, fmpz_t m);

void fmpz_montgomery_mod_init(fmpz_montgomery_t mont, fmpz_t m);

void fmpz_montgomery_mulmod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont);

// montgomery precomputed divmod (same as mulmod)
static inline
void fmpz_montgomery_divmod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont)
{
    fmpz_montgomery_mulmod(res, a, mont);
}

// montgomery precomputed mod (same as mulmod)
static inline
void fmpz_montgomery_mod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont)
{
    fmpz_montgomery_mulmod(res, a, mont);
}

#endif
