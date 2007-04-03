/****************************************************************************

Zpoly.h: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

There are two entirely separate data formats for polynomials over Z:
  -- Zpoly_t uses an array of mpz_t's
  -- Zpoly_mpn_t uses a single block of memory with each coefficient occupying
     the same number of limbs (see Zpoly_mpn files)

*****************************************************************************/
#ifndef ZPOLY_H
#define ZPOLY_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint-manager.h"

/****************************************************************************

   Zpoly_t
   -----------

Zpoly_t represents a dense polynomial in Z[x] using a single mpz_t for each
coefficient.

"coeffs" is an array of mpz_t's of length "alloc". They are all mpz_init'd.
coeffs[n] is the coefficient of x^n.

Only the first "length" coefficients actually represent coefficients of the
polynomial; i.e. it's a polynomial of degree at most length-1. There is no
requirement for coeff[length-1] to be nonzero. If length == 0, this is the
zero polynomial. Obviously always alloc >= length.

There are two classes of functions operating on Zpoly_t:

-- The _Zpoly_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough mpz_t's for
   the result.

-- The Zpoly_* functions ASSUME that "coeffs" was allocated via
   flint_malloc, and they MAY free or reallocate "coeffs" using flint_realloc,
   flint_free etc, whenever they feel the need. Furthermore they assume that
   always alloc >= 1.

*/

typedef struct
{
   mpz_t* coeffs;
   unsigned long alloc;
   unsigned long length;
} Zpoly_struct;

// Zpoly_t allows reference-like semantics for Zpoly_struct:
typedef Zpoly_struct Zpoly_t[1];


/*============================================================================
  
    Functions in _Zpoly_* layer
    
===============================================================================*/

// returns x^n coefficient with no bounds checking
static inline
mpz_t* _Zpoly_get_coeff_ptr(Zpoly_t poly, unsigned long n)
{
    return &poly->coeffs[n];
}

// copies out x^n coefficient (via mpz_set) with no bounds checking
static inline
void _Zpoly_get_coeff(mpz_t output, Zpoly_t poly,
                             unsigned long n)
{
    mpz_set(output, poly->coeffs[n]);
}

// copies out x^n coefficient (via mpz_get_ui) with no bounds checking
// the coefficient is assumed to fit into an unsigned long
static inline
unsigned long _Zpoly_get_coeff_ui(Zpoly_t poly, unsigned long n)
{
    return mpz_get_ui(poly->coeffs[n]);
}

// copies out x^n coefficient (via mpz_get_si) with no bounds checking
// the coefficient is assumed to fit into a long
static inline
long _Zpoly_get_coeff_si(Zpoly_t poly, unsigned long n)
{
    return mpz_get_si(poly->coeffs[n]);
}


// sets x^n coefficient (via mpz_set) with no bounds checking
// does NOT update poly->length if the index is beyond the end of the poly
static inline
void _Zpoly_set_coeff(Zpoly_t poly, unsigned long n, mpz_t x)
{
    mpz_set(poly->coeffs[n], x);
}

static inline
void _Zpoly_set_coeff_ui(Zpoly_t poly, unsigned long n,
                                unsigned long x)
{
    mpz_set_ui(poly->coeffs[n], x);
}

static inline
void _Zpoly_set_coeff_si(Zpoly_t poly, unsigned long n, long x)
{
    mpz_set_si(poly->coeffs[n], x);
}


// Decreases poly.length to point at the last non-zero coefficient
// (i.e. so that the degree really is length-1)
void _Zpoly_normalise(Zpoly_t poly);

// Normalises the polynomial, and returns its degree
// (the zero polynomial has degree -1 by convention)
long _Zpoly_get_degree(Zpoly_t poly);

// Normalises the polynomial, and returns its length
unsigned long _Zpoly_get_length(Zpoly_t poly);

// output = input
// assumes output.alloc >= input.length
void _Zpoly_set(Zpoly_t output, Zpoly_t input);

// output = 0
static inline
void _Zpoly_zero(Zpoly_t output)
{
   output->length = 0;
}

// swaps coeffs, alloc, length
static inline
void _Zpoly_swap(Zpoly_t x, Zpoly_t y)
{
    mpz_t* temp1;
    unsigned long temp2;

    temp1 = y->coeffs;
    y->coeffs = x->coeffs;
    x->coeffs = temp1;
    
    temp2 = x->alloc;
    x->alloc = y->alloc;
    y->alloc = temp2;

    temp2 = x->length;
    x->length = y->length;
    y->length = temp2;
}


// returns true if input1 is equal to input2. (If the lengths are different, it
// checks that the overhanging coefficients are zero.)
int _Zpoly_equal(Zpoly_t input1, Zpoly_t input2);


// output = input1 + input2
// Assumes output.alloc >= max(input1.length, input2.length)
// all combinations of parameter aliasing are allowed
// output.length is set to the maximum of the two lengths, and no normalisation
// is performed (so the output may have leading zeroes)
void _Zpoly_add(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2);
void _Zpoly_sub(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2);

// output = -input
// assumes output.alloc >= input.length
// output may alias input
void _Zpoly_negate(Zpoly_t output, Zpoly_t input);


// scalar multiplication by x
void _Zpoly_scalar_mul(Zpoly_t poly, mpz_t x);
void _Zpoly_scalar_mul_ui(Zpoly_t poly, unsigned long x);
void _Zpoly_scalar_mul_si(Zpoly_t poly, long x);
// scalar division by x
// todo: what about all the variations... floor, ceiling, truncating,
// exact division, etc?
void _Zpoly_scalar_div(Zpoly_t poly, mpz_t x);
void _Zpoly_scalar_div_ui(Zpoly_t poly, unsigned long x);


// output = input1 * input2
// assumes output.alloc >= input1.length + input2.length - 1
// output may **NOT** alias either input
void _Zpoly_mul(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2);
// as above, but always uses naive multiplication algorithm
void _Zpoly_mul_naive(Zpoly_t output, Zpoly_t input1,
                             Zpoly_t input2);
// as above, but always uses karatsuba
void _Zpoly_mul_karatsuba(Zpoly_t output, Zpoly_t input1,
                                 Zpoly_t input2);

// output = input * input
// output may NOT alias input
void _Zpoly_sqr(Zpoly_t output, Zpoly_t input);
void _Zpoly_sqr_naive(Zpoly_t output, Zpoly_t input);
void _Zpoly_sqr_karatsuba(Zpoly_t output, Zpoly_t input);


// output = x^n * input
// assumes output.alloc >= input.length + n
// output may alias input
void _Zpoly_left_shift(Zpoly_t output, Zpoly_t input,
                              unsigned long n);
// output = input / x^n (truncating division)
// assumes output.alloc >= input.length - n >= 0
// output may alias input
void _Zpoly_right_shift(Zpoly_t output, Zpoly_t input,
                               unsigned long n);


// todo: figure out better division interface once we understand the
// underlying algorithms a bit better.... in particular I'm not even sure to
// what extent these things make sense when the divisor is not monic.
// Perhaps we want specialised versions for the monic case.
// Also perhaps a power series inversion function should go here somewhere;
// i.e. computes inverse modulo a given x^n. And for that matter direct
// power series division, which might be faster than "invert & multiply".

// quotient = input1 / input2 (throw away remainder)
// asumes quotient.alloc >= max(0, input1.length - input2.length)
void _Zpoly_div(Zpoly_t quotient, Zpoly_t input1,
                       Zpoly_t input2);
// remainder = input1 % input2 (throw away quotient)
// assumes remainder.alloc >= input2.length
void _Zpoly_rem(Zpoly_t remainder, Zpoly_t input1,
                       Zpoly_t input2);
void _Zpoly_div_rem(Zpoly_t quotient, Zpoly_t remainder,
                           Zpoly_t input1, Zpoly_t input2);


// output = gcd(input1, input2)
// assumes output.alloc >= max(input1.length, input2.length)
// (or perhaps only requires output.alloc >= length of gcd?)
void _Zpoly_gcd(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2);
// also sets a, b so that a*input1 + b*input2 = output
// (is this even always possible in Z[x]?)
void _Zpoly_xgcd(Zpoly_t a, Zpoly_t b, Zpoly_t output,
                        Zpoly_t input1, Zpoly_t input2);

// Gaussian content of polynomial

void _Zpoly_content(mpz_t content, Zpoly_t a);

// ============================================================================
// Zpoly_* layer ("non-raw" functions)

// note: the implementation of many of these functions will be essentially:
// check for available space in output, reallocate if necessary, then call
// the corresponding raw version, possibly deal with aliasing restrictions
// for raw functions


// allocate coeffs to length 1, sets length = 0 (i.e. zero polynomial)
void Zpoly_init(Zpoly_t poly);

// allocate coeffs to given length, sets length = 0 (i.e. zero polynomial)
void Zpoly_init2(Zpoly_t poly, unsigned long alloc);

// allocate coeffs to given length, with space for coeff_size bits in each
// coefficient, sets length = 0 (i.e. zero polynomial)
void Zpoly_init3(Zpoly_t poly, unsigned long alloc,
                     unsigned long coeff_bits);


// Changes allocated space to be alloc.
// Current value is preserved if it fits, or truncated if it doesn't fit.
void Zpoly_realloc(Zpoly_t poly, unsigned long alloc);

// todo: should have available a version of Zpoly_realloc() that also
// allows the user to specify the space allocated for new mpz's


// This is the actual implementation that's called for Zpoly_ensure_space()
// (see below) if a reallocation is required
void Zpoly_ensure_space_IMPL(Zpoly_t poly, unsigned long alloc);

// Ensures that the polynomial has at least alloc coefficients allocated.
// If the polynomial already has enough space allocated, nothing happens.
// If more space is required, then a realloc occurs, and the space allocated
// will at least *double* from its current amount. This strategy ensures that
// repeated calls have amortised constant cost.
// (NOTE: I'm making only the initial comparison inline, since this will
// happen very frequently; the actual reallocation will be less frequent,
// and will chew up too many bytes of code if I make it inline.)
static inline
void Zpoly_ensure_space(Zpoly_t poly, unsigned long alloc)
{
   if (poly->alloc < alloc)
      Zpoly_ensure_space_IMPL(poly, alloc);
}


// deallocates space (poly becomes uninitialised)
void Zpoly_clear(Zpoly_t poly);


// returns x^n coefficient, or NULL if out of range
mpz_t* Zpoly_get_coeff_ptr(Zpoly_t poly, unsigned long n);

// copies out x^n coefficient (via mpz_set), or retrieves zero if out of range
void Zpoly_get_coeff(mpz_t output, Zpoly_t poly,
                         unsigned long n);
unsigned long Zpoly_get_coeff_ui(Zpoly_t poly, unsigned long n);
long Zpoly_get_coeff_si(Zpoly_t poly, unsigned long n);

// and the rest of them are just like the mpz_raw versions, but reallocate
// on demand

void Zpoly_set_coeff(Zpoly_t poly, unsigned long n, mpz_t x);
void Zpoly_set_coeff_ui(Zpoly_t poly, unsigned long n, unsigned long x);
void Zpoly_set_coeff_si(Zpoly_t poly, unsigned long n, long x);


// sets poly from a string, where the string is just a sequence of
// coefficients in decimal notation, separated by whitespace.
// An example string: "0 -3 45" represents 45x^2 - 3x
// Returns 1 on success, 0 on failure. If failure occurs, then only some
// of the coefficients have been read, and then some garbage was encountered
// that could not be read.
// This is NOT intended to be fast; it's for convenience of debugging only.
int Zpoly_set_from_string(Zpoly_t output, char* s);

// converts poly to a string, in the same format as accepted by
// Zpoly_set_from_string. The output buffer must be large enough for
// the output and the null terminator. See Zpoly_get_string_size().
// This is NOT intended to be fast; it's for convenience of debugging only.
void Zpoly_get_as_string(char* output, Zpoly_t poly);

// Returns the number of bytes necessary to store this poly as a string
// via Zpoly_get_as_string(). The returned value may be an overestimate.
unsigned long Zpoly_get_string_size(Zpoly_t poly);

// Prints polynomial to output (e.g. stdout) in same format as generated by
// Zpoly_get_as_string().
void Zpoly_print(FILE* output, Zpoly_t poly);


static inline
void Zpoly_normalise(Zpoly_t poly)
{
   _Zpoly_normalise(poly);
}

void Zpoly_set(Zpoly_t output, Zpoly_t input);

static inline
void Zpoly_zero(Zpoly_t output)
{
   _Zpoly_zero(output);
}

static inline
int Zpoly_equal(Zpoly_t x, Zpoly_t y)
{
   return _Zpoly_equal(x, y);
}

static inline
void Zpoly_swap(Zpoly_t x, Zpoly_t y)
{
   _Zpoly_swap(x, y);
}

void Zpoly_add(Zpoly_t output, Zpoly_t input1, Zpoly_t input2);
void Zpoly_sub(Zpoly_t output, Zpoly_t input1, Zpoly_t input2);
void Zpoly_negate(Zpoly_t output, Zpoly_t input);
void Zpoly_scalar_mul(Zpoly_t poly, mpz_t x);
void Zpoly_scalar_mul_ui(Zpoly_t poly, unsigned long x);
void Zpoly_scalar_mul_si(Zpoly_t poly, long x);
void Zpoly_scalar_div(Zpoly_t poly, mpz_t x);
void Zpoly_scalar_div_ui(Zpoly_t poly, unsigned long x);
void Zpoly_mul(Zpoly_t output, Zpoly_t input1, Zpoly_t input2);
void Zpoly_mul_naive(Zpoly_t output, Zpoly_t input1,
                         Zpoly_t input2);
void Zpoly_mul_karatsuba(Zpoly_t output, Zpoly_t input1,
                             Zpoly_t input2);
// mul_naive_KS uses a "naive" KS multiplication algorithm (i.e. it doesn't
// try very hard to pack coefficients efficiently, and it uses GMP for the
// underlying multiplication. It's for testing purposes only. The idea is to
// provide a very stable algorithm which can still keep up asymptotically with
// any other multiplication code.)
void Zpoly_mul_naive_KS(Zpoly_t output, Zpoly_t input1,
                            Zpoly_t input2);
void Zpoly_sqr(Zpoly_t output, Zpoly_t input);
void Zpoly_sqr_naive(Zpoly_t output, Zpoly_t input);
void Zpoly_sqr_karatsuba(Zpoly_t output, Zpoly_t input);
void Zpoly_sqr_naive_KS(Zpoly_t output, Zpoly_t input1,
                            Zpoly_t input2);

void Zpoly_left_shift(Zpoly_t output, Zpoly_t input,
                          unsigned long n);
void Zpoly_right_shift(Zpoly_t output, Zpoly_t input,
                           unsigned long n);

void Zpoly_div(Zpoly_t quotient, Zpoly_t input1,
                   Zpoly_t input2);
void Zpoly_rem(Zpoly_t remainder, Zpoly_t input1,
                   Zpoly_t input2);
void Zpoly_div_rem(Zpoly_t quotient, Zpoly_t remainder,
                       Zpoly_t input1, Zpoly_t input2);

void Zpoly_gcd(Zpoly_t output, Zpoly_t input1,
                       Zpoly_t input2);
void Zpoly_xgcd(Zpoly_t a, Zpoly_t b, Zpoly_t output,
                    Zpoly_t input1, Zpoly_t input2);
                    
void Zpoly_content(mpz_t content, Zpoly_t a);


// *************** end of file
#endif
