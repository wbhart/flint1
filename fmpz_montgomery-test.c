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

   fmpz_montgomery-test.c: test code for fmpz_montgomery.c/h

   Copyright (C) 2008, William Hart and Gonzalo Tornaria

*****************************************************************************/

#include "fmpz_montgomery.h"

int test_fmpz_montgomery_redc()
{
    mpz_t num1, num2, num3, num4;
    fmpz_t fnum1, fnum2, fnum3;
    unsigned long bits, bits2;
    int result = 1;

    mpz_init(num1);
    mpz_init(num2);
    mpz_init(num3);
    mpz_init(num4);

    unsigned long i;
    for (i = 0; (i < 10000) && (result == 1); i++)
    {
        // m
        bits = random_ulong(1000) + 1;
        do {
            mpz_rrandomb(num1, state, bits);
            const char * s = mpz_get_str(0, 10, num1);
            free((void *)s);
        } while( mpz_sgn(num1) == 0 || mpz_even_p(num1) );

        // x
        bits2 = 2*bits-1;
        mpz_rrandomb(num2, state, bits2);

        fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
        fnum3 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);

        mpz_to_fmpz(fnum1, num1);
        mpz_to_fmpz(fnum2, num2);

        fmpz_montgomery_t mont;

        fmpz_montgomery_init(mont, fnum1);
        fmpz_montgomery_redc(fnum3, fnum2, mont);
        fmpz_check_normalisation(fnum3);
        fmpz_to_mpz(num3, fnum3);
        fmpz_to_mpz(num4, mont->t);
        fmpz_montgomery_clear(mont);

        mpz_invert(num4, num4, num1);
        mpz_mul(num4, num4, num2);
        mpz_mod(num4, num4, num1);

        result = (mpz_cmp(num3, num4) == 0);

        fmpz_clear(fnum1);
        fmpz_clear(fnum2);
        fmpz_clear(fnum3);
    }

    mpz_clear(num1);
    mpz_clear(num2);
    mpz_clear(num3);
    mpz_clear(num4);

    return result;
}

int test_fmpz_montgomery_mulmod()
{
    mpz_t num1, num2, num3, num4, num5;
    fmpz_t fnum1, fnum2, fnum3, fnum4;
    unsigned long bits;
    int result = 1;

    mpz_init(num1);
    mpz_init(num2);
    mpz_init(num3);
    mpz_init(num4);
    mpz_init(num5);

    unsigned long i;
    for (i = 0; (i < 10000) && (result == 1); i++)
    {
        // modulus
        bits = random_ulong(1000) + 1;
        do {
            mpz_rrandomb(num1, state, bits);
        } while( mpz_sgn(num1) == 0 || mpz_even_p(num1) );

        mpz_rrandomb(num2, state, bits);
        mpz_rrandomb(num3, state, bits);

        fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum2 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum3 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum4 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);

        mpz_to_fmpz(fnum1, num1);
        mpz_to_fmpz(fnum2, num2);
        mpz_to_fmpz(fnum3, num3);

        fmpz_montgomery_t mont;

        fmpz_montgomery_mulmod_init(mont, fnum2, fnum1);
        fmpz_montgomery_mulmod(fnum4, fnum3, mont);
        fmpz_check_normalisation(fnum4);
        fmpz_to_mpz(num4, fnum4);
        fmpz_montgomery_clear(mont);

        mpz_mul(num5, num3, num2);
        mpz_mod(num5, num5, num1);

        result = (mpz_cmp(num4, num5) == 0);

        fmpz_clear(fnum1);
        fmpz_clear(fnum2);
        fmpz_clear(fnum3);
        fmpz_clear(fnum4);
    }

    mpz_clear(num1);
    mpz_clear(num2);
    mpz_clear(num3);
    mpz_clear(num4);
    mpz_clear(num5);

    return result;
}

int test_fmpz_montgomery_divmod()
{
    mpz_t num1, num2, num3, num4, num5;
    fmpz_t fnum1, fnum2, fnum3, fnum4;
    unsigned long bits;
    int result = 1;

    mpz_init(num1);
    mpz_init(num2);
    mpz_init(num3);
    mpz_init(num4);
    mpz_init(num5);

    unsigned long i;
    for (i = 0; (i < 10000) && (result == 1); i++)
    {
        // m
        bits = random_ulong(1000) + 1;

        fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum2 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum3 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum4 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);

        do {
            mpz_rrandomb(num1, state, bits);
            const char * s = mpz_get_str(0, 10, num1);
            free((void *)s);
        } while( mpz_sgn(num1) == 0 || mpz_even_p(num1) );

        mpz_to_fmpz(fnum1, num1);

        do {
            mpz_rrandomb(num2, state, bits);
            mpz_to_fmpz(fnum2, num2);
            fmpz_gcd(fnum3, fnum1, fnum2);
        } while (!fmpz_is_one(fnum3));

        mpz_rrandomb(num3, state, bits);
        mpz_to_fmpz(fnum3, num3);

        fmpz_montgomery_t mont;

        fmpz_montgomery_divmod_init(mont, fnum2, fnum1);
        fmpz_montgomery_divmod(fnum4, fnum3, mont);
        fmpz_check_normalisation(fnum4);
        fmpz_to_mpz(num4, fnum4);
        fmpz_montgomery_clear(mont);

        mpz_invert(num5, num2, num1);
        mpz_mul(num5, num5, num3);
        mpz_mod(num5, num5, num1);

        result = (mpz_cmp(num4, num5) == 0);

        fmpz_clear(fnum1);
        fmpz_clear(fnum2);
        fmpz_clear(fnum3);
        fmpz_clear(fnum4);
    }

    mpz_clear(num1);
    mpz_clear(num2);
    mpz_clear(num3);
    mpz_clear(num4);
    mpz_clear(num5);

    return result;
}

int test_fmpz_montgomery_mod()
{
    mpz_t num1, num2, num3, num4;
    fmpz_t fnum1, fnum2, fnum3;
    unsigned long bits;
    int result = 1;

    mpz_init(num1);
    mpz_init(num2);
    mpz_init(num3);
    mpz_init(num4);

    unsigned long i;
    for (i = 0; (i < 10000) && (result == 1); i++)
    {
        // modulus
        bits = random_ulong(1000) + 1;
        do {
            mpz_rrandomb(num1, state, bits);
        } while( mpz_sgn(num1) == 0 || mpz_even_p(num1) );

        mpz_rrandomb(num2, state, bits);

        fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum2 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
        fnum3 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);

        mpz_to_fmpz(fnum1, num1);
        mpz_to_fmpz(fnum2, num2);

        fmpz_montgomery_t mont;

        fmpz_montgomery_mod_init(mont, fnum1);
        fmpz_montgomery_mod(fnum3, fnum2, mont);
        fmpz_check_normalisation(fnum3);
        fmpz_to_mpz(num3, fnum3);
        fmpz_montgomery_clear(mont);

        mpz_mod(num4, num2, num1);

        result = (mpz_cmp(num3, num4) == 0);

        fmpz_clear(fnum1);
        fmpz_clear(fnum2);
        fmpz_clear(fnum3);
    }

    mpz_clear(num1);
    mpz_clear(num2);
    mpz_clear(num3);
    mpz_clear(num4);

    return result;
}
