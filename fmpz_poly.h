/****************************************************************************

fmpz_poly.h: Polynomials over Z, implemented as contiguous block of fmpz_t's

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_FMPZ_POLY_H
#define FLINT_FMPZ_POLY_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpn_extras.h"

/****************************************************************************

   fmpz_poly_t
   -----------

fmpz_poly_t represents a dense polynomial in Z[x] using a single block of
memory to hold all the coefficients.

This type is better suited to handling very dense polynomials with relatively
small coefficients, where the memory management overhead of Zpoly_t would
be too expensive.

"coeffs" is an array of limbs of length (alloc * (limbs+1)). Each
coefficient uses limbs+1 limbs. For each coefficient, the first limb is
a sign limb: 0 means positive and 1 means negative. (Zero may be stored as
either positive or negative.) The remaining "limbs" limbs represent the
absolute value of the coefficient, stored in GMP's mpn format.

Only the first "length" coefficients actually represent coefficients of the
polynomial; i.e. it's a polynomial of degree at most length-1. There is no
requirement for coeff[length-1] to be nonzero. If length == 0, this is the
zero polynomial. Obviously always alloc >= length.

There are two classes of functions operating on fmpz_poly_t:

-- The _fmpz_poly_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough space for
   the result. They also NEVER modify the limbs attribute (since this
   would screw up the block size).

-- The fmpz_poly_* functions ASSUME that "coeffs" was allocated via
   flint_malloc, and they MAY free or reallocate "coeffs" using flint_realloc,
   flint_free etc, whenever they feel the need. Furthermore they assume that
   always alloc >= 1.

*/
 
typedef struct
{
   mp_limb_t* coeffs;
   unsigned long alloc;
   unsigned long length;
   unsigned long limbs;
} fmpz_poly_struct;

// fmpz_poly_t allows reference-like semantics for fmpz_poly_struct:
typedef fmpz_poly_struct fmpz_poly_t[1];
typedef fmpz_poly_struct * fmpz_poly_p;

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

#define ABS(x) (((long) x < 0) ? -x : x)

#define SWAP(x_dummy, y_dummy) \
do { \
   fmpz_poly_p swap_temp = x_dummy; \
   x_dummy = y_dummy; \
   y_dummy = swap_temp; \
} while(0);

#define SWAP_PTRS(x_dummy_p, y_dummy_p) \
do { \
   mp_limb_t * swap_temp_p = x_dummy_p; \
   x_dummy_p = y_dummy_p; \
   y_dummy_p = swap_temp_p; \
} while(0);



/*============================================================================
  
    Functions in _fmpz_poly_* layer
    
===============================================================================*/

void _fmpz_poly_stack_init(fmpz_poly_t poly, unsigned long alloc, unsigned long limbs);

void _fmpz_poly_stack_clear(fmpz_poly_t poly);

static inline
mp_limb_t * _fmpz_poly_get_coeff_ptr(fmpz_poly_t poly, unsigned long n)
{
   return poly->coeffs+n*(poly->limbs+1);
}

/* 
   Set "output" to the given coefficient and return the sign
   Assumes length of output is poly->limbs limbs long.  
*/
   
static inline
long _fmpz_poly_get_coeff(mp_limb_t * output, fmpz_poly_t poly,
                          unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) clear_limbs(output, poly->limbs);
   copy_limbs(output, poly->coeffs+n*(poly->limbs+1)+1, poly->limbs);
   return poly->coeffs[n*(poly->limbs+1)];
}

static inline
unsigned long _fmpz_poly_get_coeff_ui(fmpz_poly_t poly, unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   else return poly->coeffs[n*(poly->limbs+1)+1];
}

static inline
long _fmpz_poly_get_coeff_si(fmpz_poly_t poly, unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   if (poly->coeffs[n*(poly->limbs+1)] == 1L) 
                                 return poly->coeffs[n*(poly->limbs+1)+1];
   else return -poly->coeffs[n*(poly->limbs+1)+1];
}

void _fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, unsigned long n);

/* 
   Set a coefficient to the given value having "size" limbs.
   Assumes that the poly->limbs is at least "size".
*/

static inline void _fmpz_poly_set_coeff(fmpz_poly_t poly, unsigned long n, 
                                  mp_limb_t * x, long sign, unsigned long size)
{
   FLINT_ASSERT(poly->limbs >= size);
   copy_limbs(poly->coeffs+n*(poly->limbs+1)+1, x, size);
   poly->coeffs[n*(poly->limbs+1)] = sign;
   if (poly->limbs > size) 
     clear_limbs(poly->coeffs+n*(poly->limbs+1)+size+1, poly->limbs-size);
}

void _fmpz_poly_set_coeff_ui(fmpz_poly_t poly, unsigned long n, unsigned long x);

void _fmpz_poly_set_coeff_si(fmpz_poly_t poly, unsigned long n, long x);

void _fmpz_poly_normalise(fmpz_poly_t poly);

static inline long _fmpz_poly_degree(fmpz_poly_t poly)
{
   return poly->length - 1;
}

static inline unsigned long _fmpz_poly_length(fmpz_poly_t poly)
{
   return poly->length;
}

static inline unsigned long _fmpz_poly_limbs(fmpz_poly_t poly)
{
   return poly->limbs;
}

// These two are the same as above, but normalise the poly first

long fmpz_poly_degree(fmpz_poly_t poly);

unsigned long fmpz_poly_length(fmpz_poly_t poly);


void _fmpz_poly_set(fmpz_poly_t output, fmpz_poly_t input);

/* 
   Zero the polynomial by setting the length to zero.
   Does not set the actual limbs to zero.
*/

static inline void _fmpz_poly_zero(fmpz_poly_t output)
{
   output->length = 0;
}

void _fmpz_poly_swap(fmpz_poly_t x, fmpz_poly_t y);

long _fmpz_poly_bits1(fmpz_poly_t poly_mpn);

long _fmpz_poly_bits(fmpz_poly_t poly_mpn);

unsigned long _fmpz_poly_max_limbs(fmpz_poly_t poly);

int _fmpz_poly_equal(fmpz_poly_p input1, fmpz_poly_p input2);

void _fmpz_poly_neg(fmpz_poly_t output, fmpz_poly_t input);

void _fmpz_poly_truncate(fmpz_poly_t poly, unsigned long trunc);

void _fmpz_poly_add(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2);

void __fmpz_poly_add_coeff_ui(mp_limb_t * output, unsigned long x);

void __fmpz_poly_add_coeff2_ui(mp_limb_t * output, unsigned long x);

void __fmpz_poly_sub_coeff_ui(mp_limb_t * output, unsigned long x);

void __fmpz_poly_add_coeffs(mp_limb_t * coeffs_out, mp_limb_t * coeffs1, mp_limb_t * coeffs2);

void __fmpz_poly_sub_coeffs(mp_limb_t * coeffs_out, mp_limb_t * coeffs1, mp_limb_t * coeffs2);

void __fmpz_poly_addmul_coeffs(mp_limb_t * res, mp_limb_t * a, mp_limb_t * b);

void _fmpz_poly_sub(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_scalar_mul(fmpz_poly_t output, fmpz_poly_t poly, mp_limb_t * x);

void _fmpz_poly_scalar_mul_ui(fmpz_poly_t output, fmpz_poly_t poly, unsigned long x);

void _fmpz_poly_scalar_mul_si(fmpz_poly_t output, fmpz_poly_t poly, long x);

void _fmpz_poly_scalar_div(fmpz_poly_t output, fmpz_poly_t poly, mp_limb_t * x);

void _fmpz_poly_scalar_div_ui(fmpz_poly_t output, fmpz_poly_t poly, unsigned long x);

void _fmpz_poly_scalar_div_si(fmpz_poly_t output, fmpz_poly_t poly, long x);

void _fmpz_poly_scalar_div_exact_ui(fmpz_poly_t output, fmpz_poly_t poly, unsigned long x);

void _fmpz_poly_scalar_div_exact_si(fmpz_poly_t output, fmpz_poly_t poly, long x);

void _fmpz_poly_mul(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_mul_naive(fmpz_poly_t output, fmpz_poly_t input1, 
                                                 fmpz_poly_t input2);
                                                 
void _fmpz_poly_mul_naive_trunc(fmpz_poly_t output, fmpz_poly_t input1, 
                                          fmpz_poly_t input2, unsigned long trunc);

void __fmpz_poly_karamul_recursive(fmpz_poly_t res, fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, unsigned long crossover);

void _fmpz_poly_mul_karatsuba(fmpz_poly_t output, fmpz_poly_t input1, 
                                                 fmpz_poly_t input2);
                                                 
void _fmpz_poly_mul_karatsuba_trunc(fmpz_poly_t output, fmpz_poly_t input1, 
                                           fmpz_poly_t input2, unsigned long trunc);
                                                 
/*
   Multiply two polynomials together using the Kronecker segmentation method.
   Currently assumes that the number of output bits per coefficient is <= 64 and
   is supplied by the parameter "bits"
*/

void _fmpz_poly_mul_KS(fmpz_poly_t output, fmpz_poly_t input1, 
                                       fmpz_poly_t input2);
                                       
void _fmpz_poly_mul_KS_trunc(fmpz_poly_t output, fmpz_poly_p input1, 
                                        fmpz_poly_p input2, unsigned long trunc);

void _fmpz_poly_mul_SS(fmpz_poly_t output, fmpz_poly_p input1, 
                                                            fmpz_poly_p input2);

void _fmpz_poly_sqr(fmpz_poly_t output, fmpz_poly_t input);

void _fmpz_poly_sqr_naive(fmpz_poly_t output, fmpz_poly_t input);

void _fmpz_poly_sqr_karatsuba(fmpz_poly_t output, fmpz_poly_t input);

void _fmpz_poly_left_shift(fmpz_poly_t output, fmpz_poly_t input, 
                                                 unsigned long n);

void _fmpz_poly_right_shift(fmpz_poly_t output, fmpz_poly_t input, unsigned long n);

void _fmpz_poly_div(fmpz_poly_t quotient, fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_rem(fmpz_poly_t remainder, fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_div_rem(fmpz_poly_t quotient, fmpz_poly_t remainder, 
                                     fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_gcd(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_xgcd(fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t output, 
                                      fmpz_poly_t input1, fmpz_poly_t input2);

void _fmpz_poly_content(mp_limb_t * content, fmpz_poly_t a);

/* Zero first n coefficients of poly, regardless of what length is */

void _fmpz_poly_zero_coeffs(fmpz_poly_t poly, unsigned long n);

/*============================================================================
  
    Functions in fmpz_poly_* layer
    
===============================================================================*/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, unsigned long alloc, unsigned long limbs);
                                              
void fmpz_poly_realloc(fmpz_poly_t poly, unsigned long alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, unsigned long alloc);

void fmpz_poly_resize_limbs(fmpz_poly_t poly, unsigned long limbs);

static inline void fmpz_poly_fit_limbs(fmpz_poly_t poly, unsigned long limbs)
{
   if (limbs > poly->limbs) fmpz_poly_resize_limbs(poly, limbs);
}

void fmpz_poly_clear(fmpz_poly_t poly);

void fmpz_poly_set_length(fmpz_poly_t poly, unsigned long length);

static inline void fmpz_poly_truncate(fmpz_poly_t poly, unsigned long length)
{
   FLINT_ASSERT(poly->length >= length);
   fmpz_poly_set_length(poly, length);
}

static inline void fmpz_poly_set_coeff_si(fmpz_poly_t poly, unsigned long n, long x)
{
   fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, 1);
   _fmpz_poly_set_coeff_si(poly, n, x);
   if (n+1 > poly->length) poly->length = n+1;
}

void fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, unsigned long n);

void fmpz_poly_mul(fmpz_poly_t output, fmpz_poly_p input1, fmpz_poly_p input2);

void fmpz_poly_div_naive(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B);

void fmpz_poly_div_karatsuba(fmpz_poly_t Q, fmpz_poly_t DQ, fmpz_poly_t A, fmpz_poly_t B);

// *************** end of file
#endif
