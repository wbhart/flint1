/*============================================================================

    F_mpz_poly.h: Polynomials over Z (FLINT 2.0 polynomials)

    Copyright (C) 2008, 2009, 2010 William Hart 
    Copyright (C) 2010 Andy Novocin.

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

==============================================================================*/

#ifndef FLINT_F_MPZ_POLY_H
#define FLINT_F_MPZ_POLY_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "mpz_poly.h"
#include "zmod_poly.h"
#include "flint.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL.h"

/*
   WARNING : 
	=======

	When implementing functions for F_mpz_poly, memory leaks can
	occur if the length of a polynomial is shortened without using
	_F_mpz_poly_set_length (which demotes thus releases any mpz_t's 
	that were being used beyong the new length). 

	In some cases it is desirable to use F_mpz_poly_truncate instead. 
	It won't increase the length of a polynomial, but if it is made
	shorter it also normalises it.
*/

/*==============================================================================

   F_mpz_poly_t
   -----------

F_mpz_poly_t represents a dense polynomial in Z[x] 

"coeffs" is an array of F_mpz's (longs), one for each coefficient

There are two things each entry in this array can represent:

1) If the most significant two bits are 01, then the entry represents
an index into an array F_mpz_arr of mpz_t's, defined in the F_mpz module, 
and the mpz_t in that array contains the coefficient.

2) Otherwise, the entry represents a signed coefficient
whose absolute value is no more than FLINT_BIT - 2 bits in length. The
coefficient is stored in twos complement format.

"alloc" is the number of allocated coefficients. Obviously we always have
alloc >= length.

"length" is the length of the polynomial. If length == 0, this is the zero 
polynomial. All functions normalise so that the top coefficient is non-zero.

================================================================================*/
 
typedef struct
{
   F_mpz * coeffs;
   ulong alloc;
   ulong length;
} F_mpz_poly_struct;

// F_mpz_poly_t allows reference-like semantics for F_mpz_poly_struct
typedef F_mpz_poly_struct F_mpz_poly_t[1];

/*****************************************************************************

   F_mpz_poly_factor_t

*****************************************************************************/

/**
 * This is the data type for storing factors for a polynomial
 * It contains an array of polynomials <code>factors</code> that contains the factors of the polynomial.
 * The variable <code>alloc<code> is the number of factors that can be stored in total.
 * <code>num_factors</code> is the number of factors currently stored.
 */
typedef struct
{
	F_mpz_poly_t* factors;
	unsigned long * exponents;
	unsigned long alloc;
	unsigned long num_factors;
} F_mpz_poly_factor_struct;

/**
 * This is the data type actually used allowing us to pass the factor array by reference
 */
typedef F_mpz_poly_factor_struct F_mpz_poly_factor_t[1];


/*===============================================================================

	Memory management

================================================================================*/

/** 
   \fn     void F_mpz_poly_init(F_mpz_poly_t poly)
   \brief  Initialise a polynomial of length zero with zero allocated coefficients
*/
void F_mpz_poly_init(F_mpz_poly_t poly);

/** 
   \fn     void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc)
   \brief  Initialise a polynomial of length zero with the given number of allocated 
	        coefficients
*/
void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc);

/** 
   \fn     void F_mpz_poly_realloc(F_mpz_poly_t poly, const ulong alloc)
   \brief  Reallocates poly to have space for precisely the given number of 
	        coefficients. If alloc = 0, then the polynomial is cleared. If
			  the alloc is smaller than the current length of the polynomial
			  the polynomial is truncated and normalised.
*/
void F_mpz_poly_realloc(F_mpz_poly_t poly, const ulong alloc);

/** 
   \fn     void F_mpz_poly_fit_length(F_mpz_poly_t poly, const ulong length)
   \brief  Expands poly, if necessary, so that it has space for the given number
	        of coefficients. This function never shrinks the polynomial, only 
			  expands it.
*/
void F_mpz_poly_fit_length(F_mpz_poly_t poly, const ulong length);

/** 
   \fn     void F_mpz_poly_clear(F_mpz_poly_t poly)
   \brief  Clear the polynomial, releasing any memory it was using.
*/
void F_mpz_poly_clear(F_mpz_poly_t poly);

/**
   \fn     void F_mpz_poly_factor_init(F_mpz_poly_factor_t fac)
   \brief  Initialises an array of F_mpz_poly's
*/
void F_mpz_poly_factor_init(F_mpz_poly_factor_t fac);

/** 
   \fn     void F_mpz_poly_factor_clear(F_mpz_poly_factor_t fac)
   \brief  Clear the polynomial array, frees any memory being used
*/
void F_mpz_poly_factor_clear(F_mpz_poly_factor_t fac);

/*===============================================================================

   F_mpz_poly_factor_t

================================================================================*/

/**
   \fn     void F_mpz_poly_factor_insert(F_mpz_poly_factor_t fac, 
               F_mpz_poly_t poly, unsigned long exp)
   \brief  Adds an extra element to the array with power exp
 */
void F_mpz_poly_factor_insert(F_mpz_poly_factor_t fac, F_mpz_poly_t poly, unsigned long exp);

/**
   \fn     void F_mpz_poly_factor_concat(F_mpz_poly_factor_t res, F_mpz_poly_factor_t fac)
   \brief  Concatenates array res and array fac and stores in array res
 */
void F_mpz_poly_factor_concat(F_mpz_poly_factor_t res, F_mpz_poly_factor_t fac);

/**
   \fn     void F_mpz_poly_factor_print(F_mpz_poly_factor_t fac)
   \brief  Dumps the array to stdout
 */
void F_mpz_poly_factor_print(F_mpz_poly_factor_t fac);

/*===============================================================================

	Subpolynomials

================================================================================*/

/** 
   \fn     void _F_mpz_poly_attach(F_mpz_poly_t poly1, F_mpz_poly_t poly2)
   \brief  Make poly1 an alias for poly2. Note poly1 must not be reallocated whilst poly2 
           is attached to it.
*/
static inline
void _F_mpz_poly_attach(F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	poly1->coeffs = poly2->coeffs;
	poly1->length = poly2->length;
}

/** 
   \fn     void _F_mpz_poly_attach_shift(F_mpz_poly_t poly1, F_mpz_poly_t poly2)
   \brief  Make poly1 an alias for poly2, but starting at the given coefficient, i.e.
	        as though poly2 had been shifted right by n. 
	        Note poly1 must not be reallocated whilst poly2 is attached to it. If 
			  n > poly2->length then poly1 will have length 0.
*/
static inline
void _F_mpz_poly_attach_shift(F_mpz_poly_t poly1, const F_mpz_poly_t poly2, const ulong n)
{
	poly1->coeffs = poly2->coeffs + n; // set coeffs to start at coeff n

	if (poly2->length >= n) poly1->length = poly2->length - n; // check n is not too large
   else poly1->length = 0;
}

/** 
   \fn     void _F_mpz_poly_attach_truncate(F_mpz_poly_t poly1, F_mpz_poly_t poly2)
   \brief  Make poly1 an alias for poly2, but as though poly2 had been truncated to the 
	        given number of coefficients. 
	        Note poly1 must not be reallocated whilst poly2 is attached to it. 
			  If n > poly2->length then poly1->length is set to poly2->length.
			  The polynomial poly1 is normalised after truncation.

*/
static inline
void _F_mpz_poly_attach_truncate(F_mpz_poly_t poly1, const F_mpz_poly_t poly2, const ulong n)
{
	ulong length;
	
	poly1->coeffs = poly2->coeffs; // set coeffs

	if (poly2->length < n) length = poly2->length; // check that n is not too large
	else length = n;

	while (length && (poly1->coeffs[length - 1] == 0)) length--; // normalise

	poly1->length = length;
}

/*===============================================================================

	Normalisation

================================================================================*/

/** 
   \fn     void _F_mpz_poly_normalise(F_mpz_poly_t poly)
   \brief  Normalise poly so that the leading coefficient is nonzero or the 
	        polynomial has length zero
*/
void _F_mpz_poly_normalise(F_mpz_poly_t poly);

/*===============================================================================

	Coefficient operations

================================================================================*/

/** 
   \fn     void F_mpz_poly_set_coeff_si(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Set coefficient n to the signed long value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff_si(F_mpz_poly_t poly, ulong n, const long x);

/** 
   \fn     void F_mpz_poly_get_coeff_si(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Return coefficient n of poly as a signed long. If n is greater than the degree
	        of poly, then zero is returned.
*/
long F_mpz_poly_get_coeff_si(const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_set_coeff_ui(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Set coefficient n to the unsigned long value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff_ui(F_mpz_poly_t poly, ulong n, const ulong x);

/** 
   \fn     void F_mpz_poly_get_coeff_ui(F_mpz_poly_t poly, ulong n, const long x)
   \brief  Return coefficient n of poly as an unsigned long. If n is greater than the degree
	        of poly, then zero is returned.
*/
ulong F_mpz_poly_get_coeff_ui(const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_set_coeff_mpz(F_mpz_poly_t poly, ulong n, const mpz_t x)
   \brief  Set coefficient n to the mpz_t value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff_mpz(F_mpz_poly_t poly, ulong n, const mpz_t x);

/** 
   \fn     void F_mpz_poly_get_coeff_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong n)
   \brief  Return coefficient n of poly as an mpz_t. If n is greater than the degree
	        of poly, then zero is returned.
*/
void F_mpz_poly_get_coeff_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_set_coeff(F_mpz_poly_t poly, ulong n, const mpz_t x)
   \brief  Set coefficient n to the F_mpz_t value x. Coefficients are numbered
	        from the constant coefficient, starting at zero.
*/
void F_mpz_poly_set_coeff(F_mpz_poly_t poly, ulong n, const F_mpz_t x);

/** 
   \fn     F_mpz * F_mpz_poly_get_coeff_ptr(F_mpz_poly_t poly, ulong n)
   \brief  Return a pointer to coefficient n of poly. Coefficients are numbered
           from the constant coefficient, starting at zero. No check is made
		   to verify that n is in range.
*/
static inline
F_mpz * F_mpz_poly_get_coeff_ptr(F_mpz_poly_t poly, ulong n)
{
   return poly->coeffs + n;
}

/*===============================================================================

	Attributes

================================================================================*/

/** 
   \fn     long F_mpz_poly_degree(const fmpz_poly_t poly)
   \brief  Returns the degree of the polynomial. If the polynomial is zero, then
	        minus one is returned for the degree.
*/
static inline 
long F_mpz_poly_degree(const F_mpz_poly_t poly)
{
   return poly->length - 1;
}

/** 
   \fn     unsigned long F_mpz_poly_length(const fmpz_poly_t poly)
   \brief  Returns the length of the polynomial. The zero polynomial has length zero.
*/
static inline 
ulong F_mpz_poly_length(const F_mpz_poly_t poly)
{
   return poly->length;
}


/*===============================================================================

	Truncation

================================================================================*/

/** 
   \fn     void F_mpz_poly_set_length(F_mpz_poly_t poly, const ulong length)
   \brief  Set the length of the polynomial to the given length. Assumes that
	        all the coefficients are valid coefficients, and does not normalise.
			  If the poly is made shorter, mpz_t's past the end are demoted.
*/
static inline
void _F_mpz_poly_set_length(F_mpz_poly_t poly, const ulong length)
{
	if (poly->length > length) // demote coefficients beyond new length
    {
      ulong i;
      for (i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);	
    } 

	poly->length = length;
}

/** 
   \fn     void F_mpz_poly_truncate(F_mpz_poly_t poly, const ulong length)
   \brief  Truncate the polynomial to the given length. It is permissible for
	        length to be greater than the current length of the polynomial, in 
			  which case nothing will happen.
*/
static inline
void F_mpz_poly_truncate(F_mpz_poly_t poly, const ulong length)
{
	if (poly->length > length) // only truncate if necessary
    {
       ulong i;
       for (i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);
	   poly->length = length;
       _F_mpz_poly_normalise(poly);
    }  
}

/*===============================================================================

	Conversions

================================================================================*/

/**
   Convert from an F_mpz_poly_factor_t to an fmpz_poly_factor_t.
*/
void fmpz_poly_factor_to_F_mpz_poly_factor(F_mpz_poly_factor_t F_fac, fmpz_poly_factor_t f_fac);

/**
   Convert from an fmpz_poly_factor_t to an F_mpz_poly_factor_t.
*/
void F_mpz_poly_factor_to_fmpz_poly_factor(fmpz_poly_factor_t f_fac, F_mpz_poly_factor_t F_fac);

/** 
   \fn     void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly)
   \brief  Convert an mpz_poly_t to an F_mpz_poly_t
*/
void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly);

/** 
   \fn     void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly)
   \brief  Convert an F_mpz_poly_t to an mpz_poly_t
*/
void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly);

/** 
   \fn     void F_mpz_poly_to_fmpz_poly(fmpz_poly_t m_poly, const F_mpz_poly_t F_poly)
   \brief  Convert an F_mpz_poly_t to an fmpz_poly_t
*/
void F_mpz_poly_to_fmpz_poly(fmpz_poly_t m_poly, const F_mpz_poly_t F_poly);

/** 
   \fn     void fmpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const fmpz_poly_t m_poly)
   \brief  Convert an fmpz_poly_t to an F_mpz_poly_t
*/
void fmpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const fmpz_poly_t m_poly);

/** 
   \fn     void F_mpz_poly_to_zmod_poly(zmod_poly_t zpol, const F_mpz_poly_t fpol)
   \brief  Convert an F_mpz_poly_t to a reduced zmod_poly_t
*/
void F_mpz_poly_to_zmod_poly(zmod_poly_t zpol, const F_mpz_poly_t fpol);

/** 
   \fn     void zmod_poly_to_F_mpz_poly(F_mpz_poly_t fpol, const zmod_poly_t zpol)
   \brief  Convert a zmod_poly_t to a F_mpz_poly_t
*/
void zmod_poly_to_F_mpz_poly(F_mpz_poly_t fpol, const zmod_poly_t zpol);

/*===============================================================================

        Input/output 

================================================================================*/

/** 
   \fn     int F_mpz_poly_from_string(F_mpz_poly_t poly, const char* s)
   \brief  Read a polynomial from a string. Format is an integer
           representing the length followed by 2 spaces, followed by a 
           space separated list of coefficients, starting with the
           constant term.
*/
int F_mpz_poly_from_string(F_mpz_poly_t poly, const char * s);

/** 
   \fn     char* F_mpz_poly_to_string(const F_mpz_poly_t poly)
   \brief  Return a char * in standard FLINT format from F_mpz_poly_t
*/
char * F_mpz_poly_to_string(const F_mpz_poly_t poly);

/** 
   \fn     char* F_mpz_poly_to_string_pretty(const F_mpz_poly_t poly, const char * x)
   \brief  Return a formated char * from F_mpz_poly_t with variable named at x
*/
char * F_mpz_poly_to_string_pretty(const F_mpz_poly_t poly, const char * x);

/** 
   \fn     void F_mpz_poly_fprint(const F_mpz_poly_t poly, FILE* f)
   \brief  Prints F_mpz_poly_t to a file stream f in standard FLINT format
*/
void F_mpz_poly_fprint(const F_mpz_poly_t poly, FILE * f);

/** 
   \fn     void F_mpz_poly_fprint_pretty(const F_mpz_poly_t poly, FILE* f, const char * x)
   \brief  Prints F_mpz_poly_t to a file stream f in pretty format with variable names at x
*/
void F_mpz_poly_fprint_pretty(const F_mpz_poly_t poly, FILE * f, const char * x);

/** 
   \fn     void F_mpz_poly_print_pretty(const F_mpz_poly_t poly, const char * x)
   \brief  Prints F_mpz_poly_t to screen in pretty format with variable named at x
*/
void F_mpz_poly_print_pretty(const F_mpz_poly_t poly, const char * x);

/** 
   \fn     int F_mpz_poly_fread(F_mpz_poly_t poly, FILE* f)
   \brief  Reads F_mpz_poly_t from file stream f
*/
int F_mpz_poly_fread(F_mpz_poly_t poly, FILE * f);

/**
   \fn     void F_mpz_poly_print(F_mpz_poly_t poly)
   \brief  Print a polynomial to stdout. Format is an integer
           representing the length followed by 2 spaces, followed by 
           a space separated list of coefficients, starting with the
           constant term.
*/
static inline
void F_mpz_poly_print(F_mpz_poly_t poly)
{
   mpz_poly_t m_poly;
   mpz_poly_init(m_poly);
   F_mpz_poly_to_mpz_poly(m_poly, poly);
   mpz_poly_print(m_poly);
   mpz_poly_clear(m_poly);
}

/*===============================================================================

	Assignment

================================================================================*/

/** 
   \fn     void F_mpz_poly_zero(F_mpz_poly_t poly)
   \brief  Sets poly to the zero polynomial
*/
static inline 
void F_mpz_poly_zero(F_mpz_poly_t poly)
{
   _F_mpz_poly_set_length(poly, 0);
}

/** 
   \fn     void F_mpz_poly_zero_coeffs(F_mpz_poly_t poly, const ulong n)
   \brief  Zero the first n coefficients of poly regardless of its length
*/
static inline
void F_mpz_poly_zero_coeffs(F_mpz_poly_t poly, const ulong n)
{
	ulong i;
   if (n >= poly->length) _F_mpz_poly_set_length(poly, 0);
	else 
	   for (i = 0; i < n; i++)
		   F_mpz_zero(poly->coeffs + i);
}

/** 
   \fn     void F_mpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2)
   \brief  Sets poly1 to equal poly2
*/
void F_mpz_poly_set(F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_swap(F_mpz_poly_t poly1, F_mpz_poly_t poly2)
   \brief  Efficiently swap poly1 and poly2
*/
void F_mpz_poly_swap(F_mpz_poly_t poly1, F_mpz_poly_t poly2);

/*===============================================================================

	Comparison

================================================================================*/

/** 
   \fn     int F_mpz_poly_equal(const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Returns 1 if poly1 and poly2 are equal (arithmetically), otherwise
	        returns 0.
*/
int F_mpz_poly_equal(const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/*===============================================================================

	Coefficient sizes

================================================================================*/

/** 
   \fn     long F_mpz_poly_max_bits1(const F_mpz_poly_t poly)
   \brief  Returns zero if any of the coefficients of poly are mpz_t's. Otherwise,
	        computes the largest number of bits n that any coefficient has and returns
	        -n if a negative coefficient exists in poly, else it returns n. Zero is 
			  returned for the zero polynomial.
*/
long F_mpz_poly_max_bits1(const F_mpz_poly_t poly);

/** 
   \fn     long F_mpz_poly_max_bits(F_mpz_poly_t poly)
   \brief  Computes the largest number of bits n that any coefficient has and returns
	        -n if a negative coefficient exists in poly, else it returns n. Zero is 
			  returned for the zero polynomial.
*/
long F_mpz_poly_max_bits(const F_mpz_poly_t poly);

/** 
   \fn     ulong F_mpz_poly_max_limbs(F_mpz_poly_t poly)
   \brief  Returns the largest number of limbs required to store the absolute value
	        of coefficients of poly. Zero is returned for the zero polynomial.
*/
ulong F_mpz_poly_max_limbs(const F_mpz_poly_t poly);

/*===============================================================================

	Reverse

================================================================================*/

/** 
   \fn     void F_mpz_poly_reverse(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong length)
   \brief  Treats poly as though it were the given length (with leading zeroes if 
	        necessary) and sets res to the reverse polynomial.
*/
void F_mpz_poly_reverse(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong length);

/*===============================================================================

	Negation

================================================================================*/

/** 
   \fn     void F_mpz_poly_neg(F_mpz_poly_t res, const F_mpz_poly_t poly)
   \brief  Set res to the negative of poly.
*/
void F_mpz_poly_neg(F_mpz_poly_t res, const F_mpz_poly_t poly);

/*===============================================================================

	Addition/subtraction

================================================================================*/

/** 
   \fn     void _F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Sets res to the sum of poly1 and poly2. No reallocation is done, i.e. res is assumed to 
	        have enough allocated coefficients for the result.
*/
void _F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Sets res to the sum of poly1 and poly2.
*/
void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void _F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Sets res to the difference of poly1 and poly2. No reallocation is done, i.e. res is assumed to 
	        have enough allocated coefficients for the result.
*/
void _F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Sets res to the difference of poly1 and poly2.
*/
void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/*===============================================================================

	Shifting

================================================================================*/

/** 
   \fn     void F_mpz_poly_left_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n)
   \brief  Multiplies poly by x^n and sets res to the result.
*/
void F_mpz_poly_left_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n);

/** 
   \fn     void F_mpz_poly_right_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n)
   \brief  Divides poly by x^n, discarding any remainder, and sets res to the result.
*/
void F_mpz_poly_right_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n);

/*===============================================================================

	Scalar multiplication

================================================================================*/

/** 
   \fn     void F_mpz_poly_scalar_mul_ui(F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong x)
   \brief  Multiply poly2 by the unsigned long x and set poly1 to the result.
*/
void F_mpz_poly_scalar_mul_ui(F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong x);

/** 
   \fn     void F_mpz_poly_scalar_mul_si(F_mpz_poly_t poly1, F_mpz_poly_t poly2, long x)
   \brief  Multiply poly2 by the signed long x and set poly1 to the result.
*/
void F_mpz_poly_scalar_mul_si(F_mpz_poly_t poly1, F_mpz_poly_t poly2, long x);

/** 
   \fn     void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, F_mpz_poly_t poly2, F_mpz_t x)
   \brief  Multiply poly2 by the F_mpz_t x and set poly1 to the result.
*/
void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, const F_mpz_poly_t poly2, const F_mpz_t x);

/*===============================================================================

	Scalar division

================================================================================*/

/**
   \fn     void F_mpz_poly_scalar_divexact(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t d)
   \brief  Divides polynomial f by the scalar d, assuming division is exact.
*/
void F_mpz_poly_scalar_divexact(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t d);

/**
   \fn     void F_mpz_poly_scalar_smod(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t p)
   \brief  Reduce each coefficient of f modulo p, but normalise to be in the
           range (-p/2, p/2].
*/
void F_mpz_poly_scalar_smod(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t p);


/*===============================================================================

	Bit packing

================================================================================*/

/** 
   \fn     void F_mpz_poly_bit_pack(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                           const ulong bits, const ulong length, const long negate)

   \brief  Pack length coefficients into the given array of limbs, with each coefficient 
	        packed into a bitfield with the given number of bits. If negate is -1L then each 
			  coefficient will be negated before being packed. The array is sign extended to 
			  the end of the array. Each coefficient is stored in bit packed format in 
			  twos complement format and if negative, 1 is borrowed from the next coefficient
			  before it is packed. The number of limbs in the array is n.
*/
void F_mpz_poly_bit_pack(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                         const ulong bits, const ulong length, const long negate);


/** 
   \fn     void F_mpz_poly_bit_pack_unsigned(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                                                    const ulong bits, const ulong length)

   \brief  Pack length coefficients into the given array with each coefficient packed
	        into a bitfield with the given number of bits. The number of limbs in the 
			  array is n.
*/
void F_mpz_poly_bit_pack_unsigned(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                                            const ulong bits, const ulong length);
/** 
   \fn     void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
                                            const ulong bits, const ulong length, const long negate, long negate2)

   \brief  As for F_mpz_poly_bit_pack, except that two arrays are packed. The second
	        is packed from the same data, except the signs are alternated before packing, 
			  starting with the sign given by negative2 for the constant coefficient.

*/
void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
                                 const ulong bits, const ulong length, const long negate, long negate2);

/** 
   \fn     void F_mpz_poly_bit_unpack(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, 
                                                    const ulong bundle, const ulong bits);
   \brief  Coefficients are unpacked from the array from fields of the given number
	        of bits in width. The coefficients are assumed to be signed. If a negative 
			  coefficient is unpacked, the next coefficient gets 1 added to it for the borrow
			  that was effectively made by that previous coefficient. The output coefficients 
			  are written to the coefficients of the F_mpz_poly_t poly_F_mpz, and a total of
			  length coefficients will be written unless there is a final zero coefficient with
			  a borrow, in which case length + 1 coefficients will be written.
*/
void F_mpz_poly_bit_unpack(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array,
                                                    const ulong length, const ulong bits);

/** 
   \fn     void F_mpz_poly_bit_unpack_unsigned(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array,  
                                                                     const ulong length, const ulong bits)

   \brief  Coefficients are unpacked from the array from fields of the given number
	        of bits in width. The coefficients are assumed to be unsigned. The output coefficients 
			  are written to the coefficients of the F_mpz_poly_t poly_F_mpz. A total of length
			  coefficients will be written unless the final coefficient is 0 and there is a borrow
			  in which case length + 1 coefficients will be written.
*/
void F_mpz_poly_bit_unpack_unsigned(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, 
                                                             const ulong length, const ulong bits);

/*===============================================================================

	Byte packing

================================================================================*/

/** 
   \fn     void F_mpz_poly_byte_pack(mp_limb_t * array, const F_mpz_poly_t poly_fmpz,
                   const unsigned long length, const unsigned long coeff_bytes, const long negate)

   \brief  Packs length coefficients of poly_fmpz down to the byte into array, each packed 
	        into a field "bytes" bytes wide.
   
           "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
           coefficients are assumed to be at least a limb wide.

           Assumes 0 < length and 0 < poly_fmpz->length
*/ 
void F_mpz_poly_byte_pack(mp_limb_t * array, const F_mpz_poly_t poly_fmpz,
                   const unsigned long length, const unsigned long coeff_bytes, const long negate);

/** 
   \fn     void F_mpz_poly_byte_pack_unsigned(mp_limb_t * array, const F_mpz_poly_t poly_fmpz,
                                       const unsigned long length, const unsigned long coeff_bytes)

   \brief  Packs length coefficients of poly_fmpz down to the byte into array, each packed 
	        into a field "bytes" bytes wide.
   
           "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
           coefficients are assumed to be at least a limb wide.

           Assumes 0 < length and 0 < poly_fmpz->length and that the coefficients are unsigned
*/ 
void F_mpz_poly_byte_pack_unsigned(mp_limb_t * array, const F_mpz_poly_t poly_fmpz,
                                                       const ulong length, const ulong coeff_bytes);

/** 
   \fn     void F_mpz_poly_byte_unpack_unsigned(F_mpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long coeff_bytes)

   \brief  Unpacks coefficients from array into poly_fmpz. Each coefficient stored is
	        assumed to be packed into a field "bytes" bytes wide. The coefficients are 
			  assumed to be unsigned. 
	
	        It is also assumed that array has one extra (zero) limb beyond what is 
	        required to store the packed coefficients.
   
           The total number of coefficients to be unpacked is given by length.
   
           "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
           coefficients are assumed to be at least a limb wide.

           Assumes 0 < length.
*/ 
void F_mpz_poly_byte_unpack(F_mpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long coeff_bytes);

/** 
   \fn     void F_mpz_poly_byte_unpack(F_mpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long coeff_bytes)

   \brief  Unpacks coefficients from array into poly_fmpz. Each coefficient stored is
	        assumed to be packed into a field "bytes" bytes wide with one bit reserved
	        for a sign bit.
	
	        It is also assumed that array has one extra (zero) limb beyond what is 
	        required to store the packed coefficients.
   
           The total number of coefficients to be unpacked is given by length.
   
           "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
           coefficients are assumed to be at least a limb wide.

           Assumes 0 < length.
*/ 
void F_mpz_poly_byte_unpack_unsigned(F_mpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long coeff_bytes);

/** 
   \fn     void F_mpz_poly_pack_bytes(F_mpz_poly_t res, 
                                             F_mpz_poly_t poly, ulong n, ulong bytes)
   \brief  Pack coefficients of poly into fields with the given number of bytes, in 
           bundles of n, into coefficients of res. It is required that bytes be at 
		   least one limb worth.
*/ 
void F_mpz_poly_pack_bytes(F_mpz_poly_t res, F_mpz_poly_t poly, ulong n, ulong bytes);

/** 
   \fn     void F_mpz_poly_unpack_bytes(F_mpz_poly_t res, 
                                            F_mpz_poly_t poly, ulong n, ulong bytes)
   \brief  Unpack packed coefficients of poly (packed into fields of the given number 
           of bytes width) into res staggering the output by n coefficients for each 
		   large input coefficient (and assuming each large coefficient stores 
		   (2*n - 1) coefficients. It is required that bytes be at least one limb 
		   worth.
*/ 
void F_mpz_poly_unpack_bytes(F_mpz_poly_t res, F_mpz_poly_t poly, ulong n, ulong bytes);

/*===============================================================================

	Multiplication

================================================================================*/

/** 
   \fn     void F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the classical 
	        algorithm.
*/
void F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_mul_classical_trunc_left(F_mpz_poly_t res, 
                       const F_mpz_poly_t poly1, const F_mpz_poly_t poly2, ulong trunc)
   \brief  Multiply poly1 by poly2 and set res to the result but with the bottom trunc
           terms zeroed.
*/
void F_mpz_poly_mul_classical_trunc_left(F_mpz_poly_t res, 
                       const F_mpz_poly_t poly1, const F_mpz_poly_t poly2, ulong trunc);
/** 
   \fn     void _F_mpz_poly_mul_kara_odd_even_recursive(F_mpz * out, F_mpz * in1, ulong len1, 
					               F_mpz * in2, ulong len2, F_mpz * scratch, ulong skip, ulong crossover)
   \brief  Recursive portion of odd/even karatsuba multiplication.

           Input polys are in1 and in2 staggered by skip. We specify a length, 
           len1 and len2 for each of the intputs. Then out will be of length len1 + len2 - 1
			  staggered by skip.

           The scratch buffer should be length len1 + len2 staggered by skip.

           All input/output/scratch polys should be initialised, and shouldn't overlap.

           Must have 1 <= len1 <= len2.

           If len1 * len2 <= crossover, we use the classical multiplication algorithm. 
           The crossover parameter is passed down recursively to subproducts.
*/
void _F_mpz_poly_mul_kara_odd_even_recursive(F_mpz * out, F_mpz * in1, ulong len1, 
					F_mpz * in2, ulong len2, F_mpz * scratch, ulong skip, ulong crossover);

/** 
   \fn     void _F_mpz_poly_mul_kara_recursive(F_mpz_poly_t out, const F_mpz_poly_t in1, const F_mpz_poly_t in2, 
	                                                          F_mpz_poly_t scratch, const ulong crossover)
   \brief  Recursive portion of ordinary karatsuba multiplication. Input lengths are assumed to be 
	        the same. Scratch is assumed to have five times the coefficients of in1 available for scratch
			  space.
*/
void _F_mpz_poly_mul_kara_recursive(F_mpz_poly_t out, const F_mpz_poly_t in1, const F_mpz_poly_t in2, 
	                                                          F_mpz_poly_t scratch, const ulong crossover);

/** 
   \fn     void F_mpz_poly_mul_karatsuba(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the karatsuba method.
*/
void F_mpz_poly_mul_karatsuba(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_mul_karatsuba_trunc_left(F_mpz_poly_t res, 
                                           F_mpz_poly_t poly1, F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result but with the bottom trunc
           terms zeroed.
*/
void F_mpz_poly_mul_karatsuba_trunc_left(F_mpz_poly_t res, 
                                            F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong trunc);
/** 
   \fn     void F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)

   \brief  Multiply poly1 by poly2 using Kronecker segmentation and store the result in res.
*/
void F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_mul_KS2(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)

   \brief  Multiply poly1 by poly2 using David Harvey's KS2 algorithm and store the result in res.
*/
void F_mpz_poly_mul_KS2(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_mul_SS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the Schoenhage-Strassen 
	        algorithm.
*/
void F_mpz_poly_mul_SS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_mul(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result. An attempt is made to choose the 
	        optimal algorithm.
*/
void _F_mpz_poly_mul(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2);
void F_mpz_poly_mul(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     void F_mpz_poly_mul_trunc_left(F_mpz_poly_t res, F_mpz_poly_t poly1, 
                                                 F_mpz_poly_t poly2, ulong trunc)
   \brief  Multiply poly1 by poly2 and set res to the result. An attempt is made to choose the 
	        optimal algorithm. The lower trunc coefficients of res will either be correct or
           set to 0.
*/
void _F_mpz_poly_mul_trunc_left(F_mpz_poly_t res, const F_mpz_poly_t poly1, 
                                              const F_mpz_poly_t poly2, const ulong trunc);
void F_mpz_poly_mul_trunc_left(F_mpz_poly_t res, const F_mpz_poly_t poly1, 
                                              const F_mpz_poly_t poly2, const ulong trunc);

/*===============================================================================

	Powering

================================================================================*/

/** 
   \fn     void F_mpz_poly_pow_ui(F_mpz_poly_t res, 
                                const F_mpz_poly_t poly1, const ulong exp)

   \brief  Set res to poly1^exp.
*/
void F_mpz_poly_pow_ui(F_mpz_poly_t res, const F_mpz_poly_t poly1, const ulong exp);

/*===============================================================================

	Division

================================================================================*/

/** 
   \fn     void F_mpz_poly_divrem_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B and set R to the remainder, i.e. A = B*Q + R.
*/
void F_mpz_poly_divrem_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, 
                                const F_mpz_poly_t A, const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_div_basecase(F_mpz_poly_t Q, const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B computing quotient only, i.e. notionally find A = B*Q + R.
*/
static inline
void F_mpz_poly_div_basecase(F_mpz_poly_t Q, 
                             const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   F_mpz_poly_divrem_basecase(Q, NULL, A, B);
}

/** 
   \fn     void F_mpz_poly_div_divconquer_recursive(F_mpz_poly_t Q, F_mpz_poly_t BQ, 
                                         const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B computing the quotient Q and product of B and Q.
*/
void F_mpz_poly_div_divconquer_recursive(F_mpz_poly_t Q, F_mpz_poly_t BQ, 
                                         const F_mpz_poly_t A, const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_div_divconquer_recursive(F_mpz_poly_t Q, F_mpz_poly_t BQ, 
                                         const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B computing the quotient Q and remainder R such that A = BQ + R.
*/
void F_mpz_poly_divrem_divconquer(F_mpz_poly_t Q, F_mpz_poly_t R, 
                                  const F_mpz_poly_t A, const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_divrem(F_mpz_poly_t Q, F_mpz_poly_t R, 
                                  const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B computing the quotient Q and remainder R such that A = BQ + R.
*/
static inline
void F_mpz_poly_divrem(F_mpz_poly_t Q, F_mpz_poly_t R, 
                                  const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   F_mpz_poly_divrem_divconquer(Q, R, A, B);
}

/*===============================================================================

	Division without remainder

================================================================================*/

/** 
   \fn     void F_mpz_poly_div_divconquer_recursive(F_mpz_poly_t Q, F_mpz_poly_t BQ, 
                                         const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B computing the quotient Q and remainder R such that A = BQ + R,
           but R will be set to contain only the low B->length - 1 terms of the full
           remainder.
*/
void F_mpz_poly_divrem_basecase_low(F_mpz_poly_t Q, F_mpz_poly_t R, 
                                    const F_mpz_poly_t A, const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_div_divconquer_recursive_low(F_mpz_poly_t Q, F_mpz_poly_t BQ, 
                                         const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Divide A by B computing the quotient Q and the low B->length - 1 limbs
           of the product of B and Q.
*/
void F_mpz_poly_div_divconquer_recursive_low(F_mpz_poly_t Q, F_mpz_poly_t BQ, 
                                         const F_mpz_poly_t A, const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_div_divconquer(F_mpz_poly_t Q, const F_mpz_poly_t A, 
                                                               const F_mpz_poly_t B)
   \brief  Divide A by B computing quotient Q only, i.e. notionally find A = B*Q + R.
*/
void F_mpz_poly_div_divconquer(F_mpz_poly_t Q, const F_mpz_poly_t A, 
                                                               const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_div(F_mpz_poly_t Q, const F_mpz_poly_t A, 
                                                               const F_mpz_poly_t B)
   \brief  Divide A by B computing quotient Q only, i.e. notionally find A = B*Q + R.
*/
static inline
void F_mpz_poly_div(F_mpz_poly_t Q, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   F_mpz_poly_div_divconquer(Q, A, B);
}

/*===============================================================================

	Exact division

================================================================================*/

/** 
   \fn     void F_mpz_poly_div(F_mpz_poly_t Q, const F_mpz_poly_t A, const ulong a_len, 
                                              const F_mpz_poly_t B, const ulong b_len)
   \brief  Divide A by B computing quotient Q only, i.e. notionally find A = B*Q + R,
           assuming that the division is exact and treating A as a polynomial of length
           a_len and B as a polynomial of length b_len.
*/
void F_mpz_poly_div_hensel(F_mpz_poly_t Q, const F_mpz_poly_t A, const ulong a_len, 
                                            const F_mpz_poly_t B, const ulong b_len);

/** 
   \fn     void F_mpz_poly_divexact(F_mpz_poly_t Q, const F_mpz_poly_t A, 
                                                                const F_mpz_poly_t B)
   \brief  Divide A by B computing quotient Q only, i.e. notionally find A = B*Q + R,
           assuming that the division is exact, i.e. R = 0.
*/
void F_mpz_poly_divexact(F_mpz_poly_t Q, const F_mpz_poly_t A, const F_mpz_poly_t B);

/*===============================================================================

	Pseudo division

================================================================================*/

/** 
   \fn     void F_mpz_poly_pseudo_divrem_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, 
                            ulong * d, const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Pseudo division of A by B. Returns Q, R and d such that l^d A = QB + R 
           for some R of length less than B, where l is the leading coefficient of 
           B.
*/
void F_mpz_poly_pseudo_divrem_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, 
                            ulong * d, const F_mpz_poly_t A, const F_mpz_poly_t B);

/** 
   \fn     void F_mpz_poly_pseudo_div_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, 
                            ulong * d, const F_mpz_poly_t A, const F_mpz_poly_t B)
   \brief  Pseudo division of A by B. Returns Q and d such that l^d A = QB + R 
           for some R of length less than B, where l is the leading coefficient 
           of B.
*/
static inline
void F_mpz_poly_pseudo_div_basecase(F_mpz_poly_t Q,  
                            ulong * d, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   F_mpz_poly_pseudo_divrem_basecase(Q, NULL, d, A, B);
}

/*===============================================================================

   Derivative

================================================================================*/

/**
   Sets der to the derivative of poly.
*/
void F_mpz_poly_derivative(F_mpz_poly_t der, F_mpz_poly_t poly);

/**
   Computes the Gaussian content c of poly. If all the coefficients are
   negative, the content will be negative also. 
*/
void F_mpz_poly_content(F_mpz_t c, const F_mpz_poly_t poly);

/*===============================================================================

   Evaluation

================================================================================*/

/**
   Evaluate poly at the double val and return the result. No cancellation
   checks are done. The routine is naive and will return an approximate
   result which will not be anywhere near correct unless there has been no 
   cancellation.
*/
double F_mpz_poly_eval_horner_d(F_mpz_poly_t poly, double val);

/**
   Computes the evaluation of poly at the double val, returning the result
   as a normalised mantissa and an exponent. The result should be extremely
   close to the correct value, regardless of cancellations, etc. However, it
   does not guarantee exact rounding (it uses mpfs, not mpfrs internally).
   If d is the return value of the function, the evaluation is given by
   d*2^exp. If d == 0 then exp is undefined.
*/
double F_mpz_poly_eval_horner_d_2exp(long * exp, F_mpz_poly_t poly, double val);

/*===============================================================================

   Miscellaneous Functions

================================================================================*/

/**
   Sets output to the same as input but with each coefficient the absolute
   value of the corresponding coefficient of input.
*/
void F_mpz_poly_scalar_abs(F_mpz_poly_t output, F_mpz_poly_t input);

/*===========================================================================

   stupid F_mpz_poly functions which just wrap fmpz_poly functions

============================================================================*/

/**
   \fn     void F_mpz_poly_gcd(F_mpz_poly_t d, F_mpz_poly_t f, F_mpz_poly_t g)
   \brief  Takes the polynomial gcd of f and g and writes to d
*/
void F_mpz_poly_gcd(F_mpz_poly_t d, F_mpz_poly_t f, F_mpz_poly_t g);


/*===========================================================================

   Computing fast/tight bounds for CLDs
      CLDs:= Coefficients of Logarithmic Derivatives.  f*g'/g

============================================================================*/

/**
   \fn      int _d_2exp_comp(double a, long ap, double b, long bp)
   \brief   An internal comparison function for a*2^ap vs. b*2^bp, outputs
            -2,-1,1, or 2 for the b pair being more than two times larger, 
            >, <=, or less than two times smaller, respectively than the a
            pair. Assumes that a and b are normalised.
*/
int _d_2exp_comp(double a, long ap, double b, long bp);

/**
   \fn      void F_mpz_poly_CLD_bound(F_mpz_t res, F_mpz_poly_t f, ulong n)
   \brief   A new approach for quickly, tightly, finding a bound on the size 
            of the coefficient n of f*g'/g for any integer polynomial factor 
            g of f. This function requires that the polynomial f be 
            squarefree and that f is not divisible by x. 
            If f is irreducible then the bound becomes f*f'/f = f'. So any
            bound smaller than the coeff of x^c of f' can't be correct. 
            (It is not unreasonable to expect the bound from this function to 
            be close to the one we would get by taking the coefficient of 
            f'. We don't test for this in the test code though, as we don't
            have an explicit proven upper bound we can use. An eyeball test
            does confirm that the bound being computed is reasonable though.)
            The algorithm implemented is described in Algorithm 6 of
            "Practical polynomial factoring in polynomial time" by Mark van
            Hoeij, Andy Novocin and William Hart, see:
            http://andy.novocin.com/issac/ISSAC_rewrite_mar_31.pdf      
            Concept is based on the fact that h = sum over all roots of g 
            (alpha_i) of (f / (x-alpha_i))  Every root alpha_i has some 
            abs_val in R.For every possible value of r either 
            B1 = 1/r^(n+1)(a_0 + a_1 r + ... + a_n r^n) or 
            B2 = 1/r^(n+1)(a_{n+1} r^(n+1) + ... + a_(N) r^N) 
            is an upper bound for the coefficient of x^n in f/(x - alpha) 
            where N is the degree of f, a_i is the absolute value of the 
            i'th coeff of f. So we begin with r = 2^0 evaluate B1 and B2 
            and either make r = 2^(-1) or 2^(1) and continue to increase 
            of decrease r until they flip-flop then we refine the severity 
            of our power adjustment (now adjust by 1/2, 1/4, etc) until B1 
            and B2 are relatively close (within a factor 1.5 of each other).  
            This bound must then be multiplied by N the largest possible 
            degree of a factor of f.

*/
void F_mpz_poly_CLD_bound(F_mpz_t res, F_mpz_poly_t f, ulong n);

/*============================================================================

   Naive '_modp' ( := Large moduli ) F_mpz_poly functions

============================================================================*/

/**
   \fn      void F_mpz_poly_div_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, 
                F_mpz_poly_t g, F_mpz_t P, ulong n)
   \brief   A power series modp division.  Using only n coeffs of f and g
               finds the lowest n coeffs of f/g assuming the remainder of f/g
               is 0 mod p. (Designed for padic CLDs) The leading coefficient of 
			   g must be invertible mod P, otherwise the function returns 0. If 
			   the division succeeds, it returns 1.
*/
int F_mpz_poly_div_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, 
                F_mpz_poly_t g, F_mpz_t P, ulong n);

/**
   \fn      void F_mpz_poly_div_upper_trunc_modp( F_mpz_t *res, F_mpz_poly_t f,
                F_mpz_poly_t g, F_mpz_t P, ulong n)
   \brief   Reversed power series division modp.  Using only the top n coeffs of 
               f and g finds the top n coeffs of f/g mod p, assuming that the 
               remainder of f/g is zero mod p.  (Designed for padic CLDS). The 
			   leading coefficient of g must be invertible mod P, otherwise
			   the function returns 0. If the division succeeds, it returns 1.
*/
int F_mpz_poly_div_upper_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, 
               F_mpz_poly_t g, F_mpz_t P, ulong n);

/*****************************************************************************

   Square-Free Factorization

*****************************************************************************/

/**
   \fn      int F_mpz_poly_is_squarefree(F_mpz_poly_t F)
   \brief   Return 1 if F is squarefree, otherwise return 0.
*/
int F_mpz_poly_is_squarefree(F_mpz_poly_t F);

/**
   \fn      void F_mpz_poly_factor_squarefree(F_mpz_poly_factor_t fac, 
               F_mpz_t content, F_mpz_poly_t F)
   \brief   Given poly F, finds the content of F which is stores at content,
               then finds a squarefree factorization stored at fac with 
			   exponents specified in fac. If F has length zero an exception
			   will be raised. If F has length 1 only the content will be set,
			   the contents of fac being undefined.
*/
void F_mpz_poly_factor_squarefree(F_mpz_poly_factor_t fac, 
               F_mpz_t content, F_mpz_poly_t F);

/*****************************************************************************

   NTL based Hensel Lifting procedures using F_mpz_mod_polys

*****************************************************************************/

void F_mpz_poly_build_hensel_tree(long * link, F_mpz_poly_t * v, F_mpz_poly_t * w, zmod_poly_factor_t fac);

void F_mpz_poly_hensel_lift(F_mpz_poly_t Gout, F_mpz_poly_t Hout, F_mpz_poly_t Aout, 
	     F_mpz_poly_t Bout, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t h, F_mpz_poly_t a, 
		                    F_mpz_poly_t b, F_mpz_t p, F_mpz_t p1, F_mpz_t big_P);

void F_mpz_poly_hensel_lift_without_inverse(F_mpz_poly_t Gout, F_mpz_poly_t Hout, 
	     F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t h, F_mpz_poly_t a, F_mpz_poly_t b, 
		                    F_mpz_t p, F_mpz_t p1, F_mpz_t big_P);

void F_mpz_poly_hensel_lift_only_inverse(F_mpz_poly_t Aout, F_mpz_poly_t Bout, F_mpz_poly_t f, 
	     F_mpz_poly_t G, F_mpz_poly_t H, F_mpz_poly_t a, F_mpz_poly_t b, 
		                    F_mpz_t p, F_mpz_t p1, F_mpz_t big_P);

void F_mpz_poly_rec_tree_hensel_lift(long * link, F_mpz_poly_t * v, F_mpz_poly_t * w, 
	        F_mpz_t p, F_mpz_poly_t f, long j, long inv, F_mpz_t p1, F_mpz_t big_P);

void F_mpz_poly_tree_hensel_lift(long * link, F_mpz_poly_t * v, F_mpz_poly_t * w, long e0, 
	            long e1, F_mpz_poly_t monic_f, long inv, long p, long r, F_mpz_t P);

ulong _F_mpz_poly_start_hensel_lift(F_mpz_poly_factor_t lifted_fac, long * link, 
	                   F_mpz_poly_t * v, F_mpz_poly_t * w, F_mpz_poly_t f, 
		                 zmod_poly_factor_t local_fac, ulong target_exp);

ulong _F_mpz_poly_continue_hensel_lift(F_mpz_poly_factor_t lifted_fac, long * link, 
	              F_mpz_poly_t * v, F_mpz_poly_t * w, F_mpz_poly_t f, ulong prev_exp, 
					     ulong current_exp, ulong target_exp, ulong p, ulong r);

void F_mpz_poly_hensel_lift_once(F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, 
	                     zmod_poly_factor_t local_fac, ulong target_exp);

/****************************************************************************

   Naive Zassenhaus

*****************************************************************************/

/**
   This function is unoptimized.  Takes Hensel lifted factors to power 
   P = p^n, the original polynomial F (and it's squarefree exponent),
   and a leading coeff (which might not be needed) ... and inserts the 
   factorisation into final_fac. The end user will not tend to use this 
   function.
*/
void F_mpz_poly_zassenhaus_naive(F_mpz_poly_factor_t final_fac, 
	              F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, 
	                      F_mpz_t P, ulong exp, F_mpz_t lc);

/*
   Given a squarefree polynomial f and an exponent, this function will
   factor f using the Zassenhaus algorithm and merge the factors into
   final_fac with the given exponents. The Zassenhaus implementation is
   not highly optimised and will struggle with more than about 10 actual
   factors, and not too many more local factors.
*/
void F_mpz_poly_factor_zassenhaus(F_mpz_poly_factor_t final_fac, 
								               ulong exp, F_mpz_poly_t f);

/****************************************************************************

   Factoring wrapper after square free part

*****************************************************************************/

/*
   This is a wrapper which makes some choices about primes, Hensel lifting, 
   Zassenhaus, this is the wrapper which does the real stuff.  Call it after 
   squarefree factoring with your square free f and it's eventual exponent.
   The factors of f will be inserted into final_fac then raised to the given 
   exponent. The factors are merged into final_fac rather than replacing what 
   is there already. It is assumed that f has no pure powers of x as factors. 
*/

void F_mpz_poly_factor_sq_fr_prim(F_mpz_poly_factor_t final_fac,
								                  ulong exp, F_mpz_poly_t f);

/*
   An internal version of the above which takes a cutoff (for the number of
   local factors) above which vHN should be used for factoring. Requires
   that f have no power of x factors. A cutoff larger than 10 will likely 
   result in very high times for factoring.
*/
void F_mpz_poly_factor_sq_fr_prim_internal(F_mpz_poly_factor_t final_fac, 
								    ulong exp, F_mpz_poly_t f, ulong cutoff);

/****************************************************************************

   Finally the real F_mpz_poly_factor

*****************************************************************************/

/*
   Returns the factorisation of G in final_fac and the content in cong.
   Does some simple pretests, finds squarefree factors and calls the wrapper.
*/
void F_mpz_poly_factor(F_mpz_poly_factor_t final_fac, 
					                           F_mpz_t cong, F_mpz_poly_t G);

/****************************************************************************

   Hoeij/Novocin approach

*****************************************************************************/

/*
   Has one sub-optimal part, namely we need F_mpz_poly_mul_trunc and left_trunc.  
   This computes bounds for and the CLDs themselves. Stores them in a matrix res 
   which comes out with r + 1 rows and 2*N or length-1 columns (if 2N is smaller 
   then length -1 then just top N and bottom N).  Takes F, Hensel Lifted factors, 
   P and N, in the case we don't choose N large enough we 'could' adapt for 
   restarting... maybe...
*/
void _F_mpz_poly_factor_CLD_mat(F_mpz_mat_t res, F_mpz_poly_t F, 
	                         F_mpz_poly_factor_t lifted_fac, F_mpz_t P, ulong N);

/*
   This function recieves the output of check_if_solved, which means that there is 
   a partition of the factors out there which makes things small. This function does 
   the trial divisions attempting to solve the problem.  Lots of comments in there, 
   and would like to do some specific testing with poorly behaving polynomials.
*/
int _F_mpz_poly_try_to_solve(int num_facs, ulong * part, 
	     F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, 
		             F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc, int safe);

/*
   The complement function to try_to_solve.  Run this first, it is a factorization 
   specific wrapper around F_mpz_mat_col_partition.
*/
int _F_mpz_mat_check_if_solved(F_mpz_mat_t M, ulong r, 
		 F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, 
		             F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc, int safe);

/*
   The actual factoring algorithm.  Set up to accept a prestarted matrix M (use the 
   identity at first) and an array of exponents (0's at first). Attempts to factor the 
   polynomial using all of the data that is available with the current level of Hensel 
   Lifting.  Attempting with a spattering of data at first then using more if it fails.  
   Returns 1 if the problem has been solved and 0 if more Hensel Lifting is needed.  
   Could be improved by not rechecking data and a partial Zassenhaus for some bizarre 
   cases (which I've yet to find or test).
*/
int F_mpz_poly_factor_sq_fr_vHN(F_mpz_poly_factor_t final_fac, 
			 F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, 
		ulong exp, F_mpz_mat_t M, int * cexpo, long U_exp, int hensel_loops);

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
