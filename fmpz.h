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

typedef mp_limb_t * fmpz_t;

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

static inline
fmpz_t fmpz_init(const unsigned long limbs)
{
   return (fmpz_t) flint_heap_alloc(limbs + 1);
}

static inline
fmpz_t fmpz_stack_init(const unsigned long limbs)
{
   return (fmpz_t) flint_stack_alloc(limbs + 1);
}

static inline
void fmpz_clear(const fmpz_t f)
{
   flint_heap_free(f);
}

static inline
void fmpz_stack_release(void)
{
   flint_stack_release();
}

static inline
unsigned long fmpz_size(const fmpz_t x)
{
   long limb = (long) x[0];
   return (unsigned long)  ((limb < 0) ? -limb : limb);
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
   if (x > 0)
   {
      res[0] = 1L;
      res[1] = x;
   }
   else if (x < 0)
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

// res must have enough space for x
void mpz_to_fmpz(fmpz_t res, const mpz_t x);

void fmpz_to_mpz(mpz_t res, const fmpz_t x);

void fmpz_add(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2);

void fmpz_add_ui_inplace(fmpz_t output, const unsigned long x);

void __fmpz_add_ui_inplace(fmpz_t output, const unsigned long x);

void fmpz_sub(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2);

void fmpz_sub_ui_inplace(fmpz_t output, const unsigned long x);

void fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b);

void __fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_mul_ui(fmpz_t output, const fmpz_t input, const unsigned long x);

void fmpz_addmul(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_div(fmpz_t res, const fmpz_t a, const fmpz_t b);

void fmpz_div_ui(fmpz_t output, const fmpz_t input, const unsigned long x);

void fmpz_pow_ui(fmpz_t output, const fmpz_t input, const unsigned long exp);

/*
   Computes the binomial coefficient next := bin(n, k) given prev = bin(n, k-1)
   The output is assumed to have enough space for the result, plus one extra limb
   (for efficiency reasons)
   Note: bin(n, k) requires at most n bits to represent it when n and k are positive
   Currently only implemented for positive n and k 
   Todo: implement this for negative n and k
*/

static inline
void fmpz_binomial_next(fmpz_t next, const fmpz_t prev, const long n, const long k)
{
   fmpz_mul_ui(next, prev, n-k+1);
   fmpz_div_ui(next, next, k);
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
void fmpz_normalise(const fmpz_t f)
{
   NORM(f);
}

#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
