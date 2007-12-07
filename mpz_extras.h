/*============================================================================

    Copyright (C) 2007, William Hart

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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

#ifndef FLINT_MPZ_EXTRAS_H
#define FLINT_MPZ_EXTRAS_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#define mpz_t mpz_t

/*-----------------------------------------------------------------------------

    Memory Management Functions
    
-----------------------------------------------------------------------------*/

mpz_t* F_mpz_alloc(void);

void F_mpz_release(void);

     
/*-----------------------------------------------------------------------------

    Modular Arithmetic
    
-----------------------------------------------------------------------------*/

void F_mpz_mulmod(mpz_t, mpz_t, mpz_t, mpz_t);

/*
  sets res to a*b modulo p
  assumes a and b are not (much) bigger than p and that res is not p       
*/
static inline void mulmod2(mpz_t res, mpz_t a, mpz_t b, mpz_t p)
{
   mpz_mul(res,a,b);
   mpz_fdiv_r(res,res,p);
}

unsigned long F_mpz_mulmod_ui(mpz_t, mpz_t, mpz_t, unsigned long);

long F_mpz_powm_long(long, long, long);

int F_mpz_sqrtmod(mpz_t, mpz_t, mpz_t);

void F_mpz_sqrtmodpklift(mpz_t, mpz_t, mpz_t, mpz_t);

void F_mpz_sqrtmodptopk(mpz_t, mpz_t, mpz_t, mpz_t, int);

int F_mpz_sqrtmodpk(mpz_t, mpz_t, mpz_t, int);


/*-----------------------------------------------------------------------------

    Number Theoretic
    
-----------------------------------------------------------------------------*/

void F_mpz_CRT(mpz_t, mpz_t, mpz_t, mpz_t, mpz_t, mpz_t);

/*===================================================================================

   Montgomery routines

====================================================================================*/

unsigned long F_mpz_mont_red(mpz_t res, mpz_t a, mpz_t m);

void F_mpz_mont_mul(mpz_t res, mpz_t a, mpz_t b, mpz_t m, mpz_t R, unsigned long n);

void F_mpz_expmod_mont(mpz_t res, mpz_t a, mpz_t exp, mpz_t m);

/*===================================================================================

   Burnikel_Ziegler Division

====================================================================================*/

void F_mpz_divrem_BZ(mpz_t Q, mpz_t R, mpz_t A, mpz_t B);

void F_mpz_rem_BZ(mpz_t R, mpz_t A, mpz_t B);

void F_mpz_mulmod_BZ(mpz_t res, mpz_t a, mpz_t b, mpz_t m);

void F_mpz_expmod_BZ(mpz_t res, mpz_t a, mpz_t exp, mpz_t m);

/*===================================================================================

   Large integer multiplication

====================================================================================*/

void F_mpz_mul(mpz_t res, mpz_t a, mpz_t b);

void __F_mpz_mul(mpz_t res, mpz_t a, mpz_t b, unsigned long twk);

#ifdef __cplusplus
 }
#endif
 
#endif
