/******************************************************************************

mpz_poly.h: Polynomials over Z, implemented as an array of mpz_t's

Copyright (C) 2007, William Hart and David Harvey

******************************************************************************/

#ifndef MPZ_POLY_H
#define MPZ_POLY_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "memory-manager.h"
#include "fmpz_poly.h"


typedef struct
{
   mpz_t* coeffs;
   unsigned long alloc;
   unsigned long init;
   unsigned long length;
} mpz_poly_struct;


// mpz_poly_t allows reference-like semantics for mpz_poly_struct:
typedef mpz_poly_struct mpz_poly_t[1];

// for some functions it's convenient to trick the compiler into letting us
// swap input argument pointers; we use mpz_poly_p for this
typedef mpz_poly_struct* mpz_poly_p;


#define SWAP_MPZ_POLY_PTRS(xxx_ptr, yyy_ptr) \
{                                            \
   mpz_poly_p zzz_ptr = xxx_ptr;             \
   xxx_ptr = yyy_ptr;                        \
   yyy_ptr = zzz_ptr;                        \
}


// ------------------------------------------------------
// Initialisation and memory management

void mpz_poly_init(mpz_poly_t poly);
void mpz_poly_clear(mpz_poly_t poly);
void mpz_poly_init2(mpz_poly_t poly, unsigned long alloc);
void mpz_poly_realloc(mpz_poly_t poly, unsigned long alloc);
void mpz_poly_init_upto(mpz_poly_t poly, unsigned long init);


// this non-inlined version REQUIRES that alloc > poly->alloc
void __mpz_poly_ensure_alloc(mpz_poly_t poly, unsigned long alloc);

// this is arranged so that the initial comparison is inlined, but the
// actual allocation is not
static inline
void mpz_poly_ensure_alloc(mpz_poly_t poly, unsigned long alloc)
{
   if (alloc > poly->alloc)
      __mpz_poly_ensure_alloc(poly, alloc);
}



// ------------------------------------------------------
// Setting/retrieving coefficients

mpz_t* mpz_poly_get_coeff_ptr(mpz_poly_t poly, unsigned long n);
void mpz_poly_get_coeff(mpz_t c, mpz_poly_t poly, unsigned long n);
unsigned long mpz_poly_get_coeff_ui(mpz_poly_t poly, unsigned long n);
long mpz_poly_get_coeff_si(mpz_poly_t poly, unsigned long n);

void mpz_poly_set_coeff(mpz_poly_t poly, unsigned long n, mpz_t c);
void mpz_poly_set_coeff_ui(mpz_poly_t poly, unsigned long n, unsigned long c);
void mpz_poly_set_coeff_si(mpz_poly_t poly, unsigned long n, long c);


static inline
mpz_t* _mpz_poly_get_coeff_ptr(mpz_poly_t poly, unsigned long n)
{
   return poly->coeffs + n;
}

static inline
void _mpz_poly_get_coeff(mpz_t c, mpz_poly_t poly, unsigned long n)
{
   mpz_set(c, poly->coeffs[n]);
}

static inline
unsigned long _mpz_poly_get_coeff_ui(mpz_poly_t poly, unsigned long n)
{
   return mpz_get_ui(poly->coeffs[n]);
}

static inline
long _mpz_poly_get_coeff_si(mpz_poly_t poly, unsigned long n)
{
   return mpz_get_si(poly->coeffs[n]);
}

static inline
void _mpz_poly_set_coeff(mpz_poly_t poly, unsigned long n, mpz_t c)
{
   mpz_set(poly->coeffs[n], c);
}

static inline
void _mpz_poly_set_coeff_ui(mpz_poly_t poly, unsigned long n, unsigned long c)
{
   mpz_set_ui(poly->coeffs[n], c);
}

static inline
void _mpz_poly_set_coeff_si(mpz_poly_t poly, unsigned long n, long c)
{
   mpz_set_si(poly->coeffs[n], c);
}


// ------------------------------------------------------
// String conversions and I/O

int mpz_poly_from_string(mpz_poly_t poly, char* s);
char* mpz_poly_to_string(mpz_poly_t poly);
void mpz_poly_print(mpz_poly_t poly);
void mpz_poly_fprint(mpz_poly_t poly, FILE* f);
int mpz_poly_read(mpz_poly_t poly);
int mpz_poly_fread(mpz_poly_t poly, FILE* f);


// ------------------------------------------------------
// Length and degree


void mpz_poly_normalise(mpz_poly_t poly);
int mpz_poly_normalised(mpz_poly_t poly);
void mpz_poly_pad(mpz_poly_t poly, unsigned long length);
void mpz_poly_truncate(mpz_poly_t res, mpz_poly_t poly, unsigned long length);


static inline
unsigned long mpz_poly_length(mpz_poly_t poly)
{
   return poly->length;
}


static inline
long mpz_poly_degree(mpz_poly_t poly)
{
   return (long) poly->length - 1;
}


// ------------------------------------------------------
// Assignment

void mpz_poly_set(mpz_poly_t res, mpz_poly_t poly);


static inline
void mpz_poly_zero(mpz_poly_t poly)
{
   poly->length = 0;
}


static inline
void mpz_poly_swap(mpz_poly_t poly1, mpz_poly_t poly2)
{
    mpz_t* temp1;
    unsigned long temp2;

    temp1 = poly2->coeffs;
    poly2->coeffs = poly1->coeffs;
    poly1->coeffs = temp1;
    
    temp2 = poly1->alloc;
    poly1->alloc = poly2->alloc;
    poly2->alloc = temp2;

    temp2 = poly1->length;
    poly1->length = poly2->length;
    poly2->length = temp2;

    temp2 = poly1->init;
    poly1->init = poly2->init;
    poly2->init = temp2;
}



// ------------------------------------------------------
// Conversions

void mpz_poly_to_fmpz_poly(fmpz_poly_t res, mpz_poly_t poly);
void fmpz_poly_to_mpz_poly(mpz_poly_t res, fmpz_poly_t poly);


// ------------------------------------------------------
// Comparison

int mpz_poly_equal(mpz_poly_t poly1, mpz_poly_t poly2);


// ------------------------------------------------------
// Addition/subtraction

void mpz_poly_add(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_sub(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_neg(mpz_poly_t res, mpz_poly_t poly);


// ------------------------------------------------------
// Shifting

void mpz_poly_lshift(mpz_poly_t res, mpz_poly_t poly, unsigned long k);
void mpz_poly_rshift(mpz_poly_t res, mpz_poly_t poly, unsigned long k);


static inline
void mpz_poly_shift(mpz_poly_t res, mpz_poly_t poly, long k)
{
   if (k >= 0)
      mpz_poly_lshift(res, poly, k);
   else
      mpz_poly_rshift(res, poly, -k);
}


// ------------------------------------------------------
// Scalar multiplication and division

void mpz_poly_scalar_mul(mpz_poly_t res, mpz_poly_t poly, mpz_t c);
void mpz_poly_scalar_mul_ui(mpz_poly_t res, mpz_poly_t poly, unsigned long c);
void mpz_poly_scalar_mul_si(mpz_poly_t res, mpz_poly_t poly, long c);

void mpz_poly_scalar_div(mpz_poly_t res, mpz_poly_t poly, mpz_t c);
void mpz_poly_scalar_div_ui(mpz_poly_t res, mpz_poly_t poly, unsigned long c);
void mpz_poly_scalar_div_si(mpz_poly_t res, mpz_poly_t poly, long c);

void mpz_poly_scalar_div_exact(mpz_poly_t res, mpz_poly_t poly, mpz_t c);
void mpz_poly_scalar_div_exact_ui(mpz_poly_t res, mpz_poly_t poly,
                                  unsigned long c);
void mpz_poly_scalar_div_exact_si(mpz_poly_t res, mpz_poly_t poly, long c);

void mpz_poly_scalar_mod(mpz_poly_t res, mpz_poly_t poly, mpz_t c);
void mpz_poly_scalar_mod_ui(mpz_poly_t res, mpz_poly_t poly, unsigned long c);


// ------------------------------------------------------
// Polynomial multiplication

void mpz_poly_mul(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_mul_naive(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_mul_karatsuba(mpz_poly_t res, mpz_poly_t poly1,
                            mpz_poly_t poly2);
void mpz_poly_mul_SS(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_mul_naive_KS(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_sqr(mpz_poly_t res, mpz_poly_t poly);
void mpz_poly_sqr_naive(mpz_poly_t res, mpz_poly_t poly);
void mpz_poly_sqr_SS(mpz_poly_t res, mpz_poly_t poly);
void mpz_poly_sqr_karatsuba(mpz_poly_t res, mpz_poly_t poly);
void mpz_poly_sqr_naive_KS(mpz_poly_t res, mpz_poly_t poly);


// exported for profiling...
unsigned long _mpz_poly_mul_karatsuba_crossover(unsigned long limbs);


// ------------------------------------------------------
// Polynomial division

void mpz_poly_monic_inverse(mpz_poly_t res, mpz_poly_t poly, unsigned long k);
void mpz_poly_pseudo_inverse(mpz_poly_t res, mpz_poly_t poly, unsigned long k);
void mpz_poly_monic_div(mpz_poly_t quot, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_pseudo_div(mpz_poly_t quot, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_monic_rem(mpz_poly_t rem, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_pseudo_rem(mpz_poly_t rem, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_monic_div_rem(mpz_poly_t quot, mpz_poly_t rem,
                            mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_pseudo_div_rem(mpz_poly_t quot, mpz_poly_t rem, 
                             mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_monic_inverse_naive(mpz_poly_t res, mpz_poly_t poly,
                                  unsigned long k);
void mpz_poly_pseudo_inverse_naive(mpz_poly_t res, mpz_poly_t poly,
                                   unsigned long k);
void mpz_poly_monic_div_naive(mpz_poly_t quot, mpz_poly_t poly1,
                              mpz_poly_t poly2);
void mpz_poly_pseudo_div_naive(mpz_poly_t quot, mpz_poly_t poly1,
                               mpz_poly_t poly2);
void mpz_poly_monic_rem_naive(mpz_poly_t rem, mpz_poly_t poly1,
                              mpz_poly_t poly2);
void mpz_poly_pseudo_rem_naive(mpz_poly_t rem, mpz_poly_t poly1,
                               mpz_poly_t poly2);
void mpz_poly_monic_div_rem_naive(mpz_poly_t quot, mpz_poly_t rem,
                                  mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_pseudo_div_rem_naive(mpz_poly_t quot, mpz_poly_t rem, 
                                   mpz_poly_t poly1, mpz_poly_t poly2);


// ------------------------------------------------------
// GCD and extended GCD

void mpz_poly_content(mpz_t x, mpz_poly_t poly);
unsigned long mpz_poly_content_ui(mpz_poly_t poly);
void mpz_poly_gcd(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2);
void mpz_poly_xgcd(mpz_poly_t res, mpz_poly_t a, mpz_poly_t b,
                   mpz_poly_t poly1, mpz_poly_t poly2);


// ------------------------------------------------------
// Miscellaneous

unsigned long mpz_poly_max_limbs(mpz_poly_t poly);
unsigned long mpz_poly_max_bits(mpz_poly_t poly);
unsigned long mpz_poly_product_max_limbs(mpz_poly_t poly1, mpz_poly_t poly2);
unsigned long mpz_poly_product_max_bits(mpz_poly_t poly1, mpz_poly_t poly2);





// ------------------------------------------------------
// Exported for testing only

void _mpz_poly_mul_kara_recursive(mpz_t* out,
          mpz_t* in1, unsigned long len1, mpz_t* in2, unsigned long len2,
          mpz_t* scratch, unsigned long skip, unsigned long crossover);


// *************** end of file
#endif
