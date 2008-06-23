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

   fmpz.h: "flat" multi-precision integer format

   Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_FMPZ_H
#define FLINT_FMPZ_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "memory-manager.h"
#include "flint.h"
#include "long_extras.h"
#include "zn_poly.h"

typedef mp_limb_t * fmpz_t;

typedef struct
{
   unsigned long* primes;
   unsigned long num_primes;
   unsigned long n; // we have 2^n primes
   unsigned long log_comb; // the number of levels of products we have computed already
   unsigned long log_res; // the number of levels of relative residues computed already
   fmpz_t ** comb; // array of arrays of products
   fmpz_t ** res;  // successive residues r_i^-1 mod r_{i+1} for pairs r_i, r_{i+1}
   fmpz_t ** temp; // temporary space for reductions
   zn_mod_t * mod;
} fmpz_comb_struct;

typedef fmpz_comb_struct fmpz_comb_t[1];

#define ABS(x) (((long) x < 0) ? -x : x)

#define NORM(coeff) \
do { \
   if ((coeff)[0]) \
   { \
      if ((long) (coeff)[0] < 0) \
      { \
         while ((!(coeff)[-(coeff)[0]]) && (coeff)[0]) (coeff)[0]++; \
      } else \
      { \
         while ((!(coeff)[(coeff)[0]]) && (coeff)[0]) (coeff)[0]--; \
      } \
   } \
} while (0);

void fmpz_comb_init(fmpz_comb_t comb, unsigned long * primes, unsigned long num_primes);

void fmpz_comb_clear(fmpz_comb_t comb);

void fmpz_check_normalisation(const fmpz_t x);

static inline
fmpz_t fmpz_init(const unsigned long limbs)
{
   return (fmpz_t) flint_heap_alloc(limbs + 1);
}

static inline
fmpz_t fmpz_realloc(fmpz_t f, const unsigned long limbs)
{
   return (fmpz_t) flint_heap_realloc(f, limbs + 1);
}

static inline
void fmpz_clear(const fmpz_t f)
{
   flint_heap_free(f);
}

void fmpz_print(fmpz_t in);

void fmpz_random_limbs2(fmpz_t x, unsigned long n);

static inline
unsigned long fmpz_size(const fmpz_t x)
{
   long limb = (long) x[0];
   return (unsigned long)  ((limb < 0L) ? -limb : limb);
}

static inline
unsigned long fmpz_bits(const fmpz_t x)
{
   unsigned long limbs = FLINT_ABS(x[0]);
   unsigned long bits = FLINT_BIT_COUNT(x[limbs]);  
   
   if (limbs == 0) return 0;
   return (((limbs-1)<<FLINT_LG_BITS_PER_LIMB) + bits);
}

// returns positive, negative or zero according to sign of x
static inline
long fmpz_sgn(const fmpz_t x)
{
   return (long) x[0];
}

// res := x
// if x == 0, then res needs room only for the control limb
// if x != 0, res needs room for one limb beyond control limb
static inline
void fmpz_set_ui(fmpz_t res, const unsigned long x)
{
   if (x) 
   {
      res[0] = 1UL;
      res[1] = x;
   }
   else
      res[0] = 0UL;
}


// same as fmpz_set_ui
static inline
void fmpz_set_si(fmpz_t res, const long x)
{
   if (x > 0L)
   {
      res[0] = 1L;
      res[1] = x;
   }
   else if (x < 0L)
   {
      res[0] = -1L;
      res[1] = -x;
   }
   else
      res[0] = 0UL;
}


// returns nonzero if op1 == op2
static inline
int fmpz_equal(const fmpz_t op1, const fmpz_t op2)
{
   // if the signs/sizes are different, they can't be equal
   if (op1[0] != op2[0])
      return 0;

   // compare actual limbs
   for (long i = 0; i < fmpz_size(op1); i++)
   {
      if (op1[i+1] != op2[i+1])
         return 0;
   }
     
   return 1;
}

// sets res := op
// doesn't check for aliasing (i.e. if op == res, it will stupidly copy data)
// assumes res has enough room
static inline
void fmpz_set(fmpz_t res, const fmpz_t op)
{
   long i = fmpz_size(op);
   do
   {
      res[i] = op[i];
      i--;
   }
   while (i >= 0);
}

// Absolute value
static inline
void fmpz_abs(fmpz_t res, const fmpz_t x)
{
   if (res != x) // check for aliasing
       fmpz_set(res, x);

   if (fmpz_sgn(res) < 0)
       res[0] = -res[0];
}

// Negative
static inline
void fmpz_neg(fmpz_t res, const fmpz_t x)
{
   if (res != x) // check for aliasing
       fmpz_set(res, x);

   res[0] = -res[0];
}

// res must have enough space for x
void mpz_to_fmpz(fmpz_t res, const mpz_t x);

void fmpz_to_mpz(mpz_t res, const fmpz_t x);

void fmpz_add(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2);

void fmpz_add_ui_inplace(fmpz_t output, const unsigned long x);

void fmpz_add_ui(fmpz_t output, const fmpz_t input, const unsigned long x);

void __fmpz_add_ui_inplace(fmpz_t output, const unsigned long x);

void fmpz_sub(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2);

void fmpz_sub_ui_inplace(fmpz_t output, const unsigned long x);

void fmpz_sub_ui(fmpz_t output, const fmpz_t input, const unsigned long x);

void fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b);

void __fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_mul_ui(fmpz_t output, const fmpz_t input, const unsigned long x);

void fmpz_addmul(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_tdiv(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_fdiv(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_mod(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_tdiv_ui(fmpz_t output, const fmpz_t input, const unsigned long x);

unsigned long fmpz_mod_ui(const fmpz_t input, const unsigned long x);

void fmpz_pow_ui(fmpz_t output, const fmpz_t input, const unsigned long exp);

unsigned long __fmpz_power_of_two(const fmpz_t x);

void fmpz_div_2exp(fmpz_t output, fmpz_t x, unsigned long exp);

void fmpz_mul_2exp(fmpz_t output, fmpz_t x, unsigned long exp);

void fmpz_gcd(fmpz_t output, fmpz_t x1, fmpz_t x2);

void fmpz_invert(fmpz_t res, fmpz_t x, fmpz_t m);

/*
   Computes the binomial coefficient next := bin(n, k) given prev = bin(n, k-1)
   The output is assumed to have enough space for the result, plus one extra limb
   (for efficiency reasons)
   Note: bin(n, k) requires at most n bits to represent it when n and k are positive
   Currently only implemented for positive n and k 
   Todo: implement this for negative n and k
*/

static inline
void __fmpz_binomial_next(fmpz_t next, const fmpz_t prev, const long n, const long k)
{
   fmpz_mul_ui(next, prev, n-k+1);
   fmpz_tdiv_ui(next, next, k);
}

static inline
int fmpz_is_one(const fmpz_t f)
{
   if (f[0] == 1L) return (f[1] == 1L);
   else return 0;
}

static inline
int fmpz_is_zero(const fmpz_t f)
{
   return (f[0] == 0L);
}

static inline
void __fmpz_normalise(const fmpz_t f)
{
   NORM(f);
}

static inline
int fmpz_cmpabs(const fmpz_t f1, const fmpz_t f2)
{
   unsigned long size1 = FLINT_ABS(f1[0]);
   unsigned long size2 = FLINT_ABS(f2[0]);

   if (size1 < size2) return -1;
   if (size1 > size2) return 1;

   return mpn_cmp(f1 + 1, f2 + 1, size1);
   
}

void fmpz_sqrtrem(fmpz_t sqrt, fmpz_t rem, fmpz_t n);

/*
   Computes the unique integer mod m1*m2 which is r1 mod m1 and r2 mod m2
   where m1 is an fmpz_t and m2 is a coprime unsigned long. 

   Assumes both m1 and m2 are reduced modulo their respective moduli.
   c must be set to m1^{-1} mod m2.
   pre must be set to a precomputed inverse of m2
   Assumes the number of bits of m2 is at most FLINT_D_BITS-1
*/

static inline
void fmpz_CRT_ui_precomp(fmpz_t out, fmpz_t r1, fmpz_t m1, unsigned long r2,
                                  unsigned long m2, unsigned long c, double pre)
{
   unsigned long r1mod = fmpz_mod_ui(r1, m2);
   unsigned long s = z_submod(r2, r1mod, m2);
   s = z_mulmod_precomp(s, c, m2, pre); 
   fmpz_t sm1 = fmpz_init(m1[0] + 1);
   fmpz_mul_ui(sm1, m1, s);
   fmpz_add(out, r1, sm1);
   fmpz_clear(sm1);
}

/*
   As for fmpz_CRT_ui_precomp except that it assumes the number of bits of m2 
   is at most FLINT_BITS-1
*/

static inline
void fmpz_CRT_ui2_precomp(fmpz_t out, fmpz_t r1, fmpz_t m1, unsigned long r2,
                                  unsigned long m2, unsigned long c, double pre)
{
   unsigned long r1mod = fmpz_mod_ui(r1, m2);
   unsigned long s = z_submod(r2, r1mod, m2);
   s = z_mulmod2_precomp(s, c, m2, pre); 
   fmpz_t sm1 = fmpz_init(m1[0] + 1);
   fmpz_mul_ui(sm1, m1, s);
   fmpz_add(out, r1, sm1);
   fmpz_clear(sm1);
}

void fmpz_multi_mod_ui(unsigned long * out, fmpz_t in, fmpz_comb_t comb);

void fmpz_multi_crt_ui(fmpz_t output, unsigned long * residues, fmpz_comb_t comb);

void fmpz_multi_crt_sign(fmpz_t output, fmpz_t input, fmpz_comb_t comb);

#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
