/****************************************************************************

Zpoly_mpn.h: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

There are two entirely separate data formats for polynomials over Z:
  -- Zpoly_t uses an array of mpz_t's (see Zpoly.c and Zpoly.h files)
  -- Zpoly_mpn_t uses a single block of memory with each coefficient occupying
     the same number of limbs 

*****************************************************************************/
#ifndef FLINT_ZPOLY_MPN_H
#define FLINT_ZPOLY_MPN_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint-manager.h"
#include "Zpoly.h"
#include "mpn_extras.h"

/****************************************************************************

   Zpoly_mpn_t
   -----------

Zpoly_mpn_t represents a dense polynomial in Z[x] using a single block of
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

There are two classes of functions operating on Zpoly_mpn_t:

-- The _Zpoly_mpn_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough space for
   the result. They also NEVER modify the limbs attribute (since this
   would screw up the block size).

-- The Zpoly_mpn_* functions ASSUME that "coeffs" was allocated via
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
} Zpoly_mpn_struct;

// Zpoly_mpn_t allows reference-like semantics for Zpoly_mpn_struct:
typedef Zpoly_mpn_struct Zpoly_mpn_t[1];

/*============================================================================
  
    Functions in _Zpoly_mpn_* layer
    
===============================================================================*/

void _Zpoly_mpn_convert_out(Zpoly_t poly_mpz, Zpoly_mpn_t poly_mpn);

void _Zpoly_mpn_convert_in(Zpoly_mpn_t poly_mpn, Zpoly_t poly_mpz);
                     
static inline
mp_limb_t * _Zpoly_mpn_get_coeff_ptr(Zpoly_mpn_t poly, unsigned long n)
{
   return poly->coeffs+n*(poly->limbs+1);
}

/* 
   Set "output" to the given coefficient and return the sign
   Assumes length of output is poly->limbs limbs long.  
*/
   
static inline long _Zpoly_mpn_get_coeff(mp_limb_t * output, Zpoly_mpn_t poly, 
                                                                unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) clear_limbs(output, poly->limbs);
   copy_limbs(output, poly->coeffs+n*(poly->limbs+1)+1, poly->limbs);
   return poly->coeffs[n*(poly->limbs+1)];
}

static inline unsigned long _Zpoly_mpn_get_coeff_ui(Zpoly_mpn_t poly, unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   else return poly->coeffs[n*(poly->limbs+1)+1];
}

static inline long _Zpoly_mpn_get_coeff_si(Zpoly_mpn_t poly, unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   if (poly->coeffs[n*(poly->limbs+1)] == 1L) 
                                 return poly->coeffs[n*(poly->limbs+1)+1];
   else return -poly->coeffs[n*(poly->limbs+1)+1];
}

/* 
   Set a coefficient to the given value having "size" limbs.
   Assumes that the poly->limbs is at least "size".
*/

static inline void _Zpoly_mpn_set_coeff(Zpoly_mpn_t poly, unsigned long n, 
                                  mp_limb_t * x, long sign, unsigned long size)
{
   FLINT_ASSERT(poly->limbs >= size);
   copy_limbs(poly->coeffs+n*(poly->limbs+1)+1, x, size);
   poly->coeffs[n*(poly->limbs+1)] = sign;
   if (poly->limbs > size) 
     clear_limbs(poly->coeffs+n*(poly->limbs+1)+size+1, poly->limbs-size);
}

void _Zpoly_mpn_set_coeff_ui(Zpoly_mpn_t poly, unsigned long n, unsigned long x);

void _Zpoly_mpn_set_coeff_si(Zpoly_mpn_t poly, unsigned long n, long x);

void _Zpoly_mpn_normalise(Zpoly_mpn_t poly);

static inline long _Zpoly_mpn_degree(Zpoly_mpn_t poly)
{
   return poly->length - 1;
}

static inline unsigned long _Zpoly_mpn_length(Zpoly_mpn_t poly)
{
   return poly->length;
}

static inline unsigned long _Zpoly_mpn_limbs(Zpoly_mpn_t poly)
{
   return poly->limbs;
}

// These two are the same as above, but normalise the poly first

long Zpoly_mpn_degree(Zpoly_mpn_t poly);

unsigned long Zpoly_mpn_length(Zpoly_mpn_t poly);


void _Zpoly_mpn_set(Zpoly_mpn_t output, Zpoly_mpn_t input);

/* 
   Zero the polynomial by setting the length to zero.
   Does not set the actual limbs to zero.
*/

static inline void _Zpoly_mpn_zero(Zpoly_mpn_t output)
{
   output->length = 0;
}

void _Zpoly_mpn_swap(Zpoly_mpn_t x, Zpoly_mpn_t y);

int _Zpoly_mpn_equal(Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_add(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_sub(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_negate(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_mpn_scalar_mul(Zpoly_mpn_t poly, mp_limb_t * x);

void _Zpoly_mpn_scalar_mul_ui(Zpoly_mpn_t poly, unsigned long x);

void _Zpoly_mpn_scalar_mul_si(Zpoly_mpn_t poly, long x);

void _Zpoly_mpn_scalar_div(Zpoly_mpn_t poly, mp_limb_t * x);

void _Zpoly_mpn_scalar_div_ui(Zpoly_mpn_t poly, unsigned long x);

void _Zpoly_mpn_scalar_div_si(Zpoly_mpn_t poly, unsigned long x);

void _Zpoly_mpn_mul(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_mul_naive(Zpoly_mpn_t output, Zpoly_mpn_t input1, 
                                                 Zpoly_mpn_t input2);

void _Zpoly_mpn_mul_karatsuba(Zpoly_mpn_t output, Zpoly_mpn_t input1, 
                                                 Zpoly_mpn_t input2);

void _Zpoly_mpn_sqr(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_mpn_sqr_naive(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_mpn_sqr_karatsuba(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_mpn_left_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, 
                                                 unsigned long n);

void _Zpoly_mpn_right_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, unsigned long n);

void _Zpoly_mpn_div(Zpoly_mpn_t quotient, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_rem(Zpoly_mpn_t remainder, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_div_rem(Zpoly_mpn_t quotient, Zpoly_mpn_t remainder, 
                                     Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_gcd(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_xgcd(Zpoly_mpn_t a, Zpoly_mpn_t b, Zpoly_mpn_t output, 
                                      Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mpn_content(mp_limb_t * content, Zpoly_mpn_t a);


/*============================================================================
  
    Functions in Zpoly_mpn_* layer
    
===============================================================================*/

void Zpoly_mpn_init(Zpoly_mpn_t poly, unsigned long alloc,
                                              unsigned long limbs);
                                              
void Zpoly_mpn_realloc(Zpoly_mpn_t poly, unsigned long alloc);

void Zpoly_mpn_clear(Zpoly_mpn_t poly);
                                              

// *************** end of file
#endif
