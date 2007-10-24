#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <gmp.h>

#ifndef FLINT_Z_H
#define FLINT_Z_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#define Z_t mpz_t

/*-----------------------------------------------------------------------------

    Memory Management Functions
    
-----------------------------------------------------------------------------*/

Z_t* Z_alloc(void);

void Z_release(void);

     
/*-----------------------------------------------------------------------------

    Modular Arithmetic
    
-----------------------------------------------------------------------------*/

void Z_mulmod(Z_t, Z_t, Z_t, Z_t);

/*
  sets res to a*b modulo p
  assumes a and b are not (much) bigger than p and that res is not p       
*/
static inline void mulmod2(mpz_t res, mpz_t a, mpz_t b, mpz_t p)
{
   mpz_mul(res,a,b);
   mpz_fdiv_r(res,res,p);
}

unsigned long Z_mulmod_ui(Z_t, Z_t, Z_t, unsigned long);

long Z_powm_long(long, long, long);

int Z_sqrtmod(Z_t, Z_t, Z_t);

void Z_sqrtmodpklift(Z_t, Z_t, Z_t, Z_t);

void Z_sqrtmodptopk(Z_t, Z_t, Z_t, Z_t, int);

int Z_sqrtmodpk(Z_t, Z_t, Z_t, int);


/*-----------------------------------------------------------------------------

    Number Theoretic
    
-----------------------------------------------------------------------------*/

void Z_CRT(Z_t, Z_t, Z_t, Z_t, Z_t, Z_t);

/*===================================================================================

   Montgomery routines

====================================================================================*/

unsigned long Z_mont_red(mpz_t res, mpz_t a, mpz_t m);

void Z_mont_mul(mpz_t res, mpz_t a, mpz_t b, mpz_t m, mpz_t R, unsigned long n);

void Z_expmod_mont(mpz_t res, mpz_t a, mpz_t exp, mpz_t m);

/*===================================================================================

   Jebelean Division

====================================================================================*/

void Z_divrem_jebelean(mpz_t Q, mpz_t R, mpz_t A, mpz_t B);

void Z_rem_jebelean(mpz_t R, mpz_t A, mpz_t B);

void Z_mulmod_jebelean(mpz_t res, mpz_t a, mpz_t b, mpz_t m);

void Z_expmod_jebelean(mpz_t res, mpz_t a, mpz_t exp, mpz_t m);

/*===================================================================================

   Large integer multiplication

====================================================================================*/

void F_mpz_mul(mpz_t res, mpz_t a, mpz_t b);

void __F_mpz_mul(mpz_t res, mpz_t a, mpz_t b, unsigned long twk);

#ifdef __cplusplus
 }
#endif
 
#endif
