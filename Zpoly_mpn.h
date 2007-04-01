/****************************************************************************

Zpoly.h: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

There are two entirely separate data formats for polynomials over Z:
  -- Zpoly_t uses an array of mpz_t's (see Zpoly.c and Zpoly.h files)
  -- Zpoly_mpn_t uses a single block of memory with each coefficient occupying
     the same number of limbs 

*****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

/****************************************************************************

   Zpoly_mpn_t
   -----------

Zpoly_mpn_t represents a dense polynomial in Z[x] using a single block of
memory to hold all the coefficients.

This type is better suited to handling very dense polynomials with relatively
small coefficients, where the memory management overhead of Zpoly_t would
be too expensive.

"coeffs" is an array of limbs of length (alloc * (coeff_size+1)). Each
coefficient uses coeff_size+1 limbs. For each coefficient, the first limb is
a sign limb: 0 means positive and 1 means negative. (Zero may be stored as
either positive or negative.) The remaining "coeff_size" limbs represent the
absolute value of the coefficient, stored in GMP's mpn format.

Only the first "length" coefficients actually represent coefficients of the
polynomial; i.e. it's a polynomial of degree at most length-1. There is no
requirement for coeff[length-1] to be nonzero. If length == 0, this is the
zero polynomial. Obviously always alloc >= length.

There are two classes of functions operating on Zpoly_mpn_t:

-- The _Zpoly_mpn_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough space for
   the result. They also NEVER modify the coeff_size attribute (since this
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
   unsigned long coeff_size;
} Zpoly_mpn_struct;

// Zpoly_mpn_t allows reference-like semantics for Zpoly_mpn_struct:
typedef Zpoly_mpn_struct Zpoly_mpn_t[1];

/*============================================================================
  
    Functions in _Zpoly_* layer
    
===============================================================================*/

mp_limb_t * _Zpoly_get_coeff_ptr(Zpoly_mpn_t poly, unsigned long n);

void _Zpoly_get_coeff(mp_limb_t * output, Zpoly_mpn_t poly, unsigned long n);

unsigned long _Zpoly_get_coeff_ui(Zpoly_mpn_t poly, unsigned long n);

long _Zpoly_get_coeff_si(Zpoly_mpn_t poly, unsigned long n);

void _Zpoly_set_coeff(Zpoly_mpn_t poly, unsigned long n, mp_limb_t * x, 
                                                           unsigned long size);

void _Zpoly_set_coeff_ui(Zpoly_mpn_t poly, unsigned long n, unsigned long x);

void _Zpoly_set_coeff_si(Zpoly_mpn_t poly, unsigned long n, long x);

void _Zpoly_normalise(Zpoly_mpn_t poly);

long _Zpoly_get_degree(Zpoly_mpn_t poly);

unsigned long _Zpoly_get_length(Zpoly_mpn_t poly);

unsigned long _Zpoly_get_coeff_size(Zpoly_mpn_t poly);

void _Zpoly_set(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_zero(Zpoly_mpn_t output);

void _Zpoly_swap(Zpoly_mpn_t x, Zpoly_mpn_t y);

int _Zpoly_equal(Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_add(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_sub(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_negate(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_scalar_mul(Zpoly_mpn_t poly, mp_limb_t * x);

void _Zpoly_scalar_mul_ui(Zpoly_mpn_t poly, unsigned long x);

void _Zpoly_scalar_mul_si(Zpoly_mpn_t poly, long x);

void _Zpoly_scalar_div(Zpoly_mpn_t poly, mp_limb_t * x);

void _Zpoly_scalar_div_ui(Zpoly_mpn_t poly, unsigned long x);

void _Zpoly_scalar_div_si(Zpoly_mpn_t poly, unsigned long x);

void _Zpoly_mul(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_mul_naive(Zpoly_mpn_t output, Zpoly_mpn_t input1, 
                                                 Zpoly_mpn_t input2);

void _Zpoly_mul_karatsuba(Zpoly_mpn_t output, Zpoly_mpn_t input1, 
                                                 Zpoly_mpn_t input2);

void _Zpoly_sqr(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_sqr_naive(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_sqr_karatsuba(Zpoly_mpn_t output, Zpoly_mpn_t input);

void _Zpoly_left_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, 
                                                 unsigned long n);

void _Zpoly_right_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, unsigned long n);

void _Zpoly_div(Zpoly_mpn_t quotient, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_rem(Zpoly_mpn_t remainder, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_div_rem(Zpoly_mpn_t quotient, Zpoly_mpn_t remainder, 
                                     Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_gcd(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_xgcd(Zpoly_mpn_t a, Zpoly_mpn_t b, Zpoly_mpn_t output, 
                                      Zpoly_mpn_t input1, Zpoly_mpn_t input2);

void _Zpoly_content(mp_limb_t * content, Zpoly_mpn_t a);

// *************** end of file
