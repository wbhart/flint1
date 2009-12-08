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

   fmpz_montgomery.c: "flat" integer format functions using Montgomery format

   Copyright (C) 2008, William Hart and Gonzalo Tornaria

*****************************************************************************/

#include "fmpz.h"
#include "fmpz_montgomery.h"
#include "flint.h"

/// check out
// http://icpc.baylor.edu/Past/icpc2004/RegReport/guan.cse.nsysu.edu.tw/data/montg.pdf

void fmpz_montgomery_clear(fmpz_montgomery_t mont)
{
    fmpz_clear(mont->R);
    fmpz_clear(mont->M);
    fmpz_clear(mont->t);
    if (mont->b)
        fmpz_clear(mont->b);
}

// TODO: add precomputations for the multiplication by M and R
// FIXME: check for m = 0, raise error
// FIXME: check for even m, raise error
void fmpz_montgomery_init(fmpz_montgomery_t mont, fmpz_t m)
{
    if (fmpz_sgn(m) == 0) {
        printf("Error: division by zero!\n");
        abort();
    }

    unsigned long limbs = fmpz_size(m);

    // M := abs(m)
    mont->M = fmpz_init(limbs);
    fmpz_abs(mont->M, m);

    // t := 2^#
    mont->t = fmpz_init(limbs+1);
    F_mpn_clear(mont->t + 1, limbs);
    mont->t[limbs+1] = 1L;
    mont->t[0] = limbs+1;

    // R := 1/M mod t
    mont->R = fmpz_init(limbs);
    fmpz_invert(mont->R, mont->M, mont->t);

    mont->b = 0;
}

// Computes x / 2^# mod M
// Assumes 0 <= x < 2^# * M
void fmpz_montgomery_redc(fmpz_t res, fmpz_t x, fmpz_montgomery_t mont)
{
    unsigned long limbs = fmpz_size(mont->M);

    // q := abs(x) * R  %  2^#
    fmpz_t q = fmpz_init(limbs);
    fmpz_mul_trunc(q, x, mont->R, limbs);

    // a := (abs(x) - q*M)
    fmpz_neg(q, q);
    fmpz_t a = fmpz_init(FLINT_MAX(fmpz_size(q) + limbs, fmpz_size(x)));
    fmpz_abs(a, x);
    fmpz_addmul(a, q, mont->M);
    fmpz_clear(q);

    if (fmpz_is_zero(a)) {
        res[0] = 0;
        fmpz_clear(a);
        return;
    }

    if (fmpz_sgn(a) >= 0) {
        // res := a / 2^#
        F_mpn_copy(res+1, a+limbs+1, fmpz_size(a)-limbs);
        res[0] = a[0] - limbs;
    } else {
        // res := a / 2^# + M
        F_mpn_copy(res+1, a+limbs+1, fmpz_size(a)-limbs);
        res[0] = a[0] + limbs;
        fmpz_add(res, res, mont->M);
    }

    fmpz_clear(a);
}

// initialize montgomery to do multiplication by b mod m
void fmpz_montgomery_mulmod_init(fmpz_montgomery_t mont, fmpz_t b, fmpz_t m)
{
    fmpz_montgomery_init(mont, m);
    unsigned long limbs = fmpz_size(m);
    mont->b = fmpz_init(limbs);
    fmpz_mulmod(mont->b, mont->t, b, m);
}

// initialize montgomery to do division by b mod m (as multiplication by 1/b mod m)
void fmpz_montgomery_divmod_init(fmpz_montgomery_t mont, fmpz_t b, fmpz_t m)
{
    fmpz_montgomery_init(mont, m);
    unsigned long limbs = fmpz_size(m);
    mont->b = fmpz_init(limbs);
    fmpz_divmod(mont->b, mont->t, b, m);
}

// initialize montgomery to do reduction mod m (as multiplication by 1 mod m)
void fmpz_montgomery_mod_init(fmpz_montgomery_t mont, fmpz_t m)
{
    fmpz_montgomery_init(mont, m);
    unsigned long limbs = fmpz_size(m);
    mont->b = fmpz_init(limbs);
    fmpz_mod(mont->b, mont->t, m);
}

// montgomery precomputed mulmod
void fmpz_montgomery_mulmod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont)
{
    fmpz_t x = fmpz_init(fmpz_size(a) + fmpz_size(mont->b));
    fmpz_mul(x, a, mont->b);
    fmpz_montgomery_redc(res, x, mont);
    fmpz_clear(x);
}

