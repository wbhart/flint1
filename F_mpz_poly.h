/*============================================================================

    F_mpz_poly.h: Polynomials over Z (FLINT 2.0 polynomials)

    Copyright (C) 2008, William Hart 

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
   \brief  Return coefficient n of poly as an usigned long. If n is greater than the degree
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
unsigned long F_mpz_poly_length(const F_mpz_poly_t poly)
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
      for (ulong i = length; i < poly->length; i++)
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
      for (ulong i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);
		poly->length = length;
      _F_mpz_poly_normalise(poly);
   }  
}

/*===============================================================================

	Conversions

================================================================================*/

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
   \brief  Read F_mpz_poly_t from a char *
*/
int F_mpz_poly_from_string(F_mpz_poly_t poly, const char* s);

/** 
   \fn     char* F_mpz_poly_to_string(const F_mpz_poly_t poly)
   \brief  Return a char * in standard FLINT format from F_mpz_poly_t
*/
char* F_mpz_poly_to_string(const F_mpz_poly_t poly);

/** 
   \fn     char* F_mpz_poly_to_string_pretty(const F_mpz_poly_t poly, const char * x)
   \brief  Return a formated char * from F_mpz_poly_t with variable named at x
*/
char* F_mpz_poly_to_string_pretty(const F_mpz_poly_t poly, const char * x);

/** 
   \fn     void F_mpz_poly_fprint(const F_mpz_poly_t poly, FILE* f)
   \brief  Prints F_mpz_poly_t to a file stream f in standard FLINT format
*/
void F_mpz_poly_fprint(const F_mpz_poly_t poly, FILE* f);

/** 
   \fn     void F_mpz_poly_fprint_pretty(const F_mpz_poly_t poly, FILE* f, const char * x)
   \brief  Prints F_mpz_poly_t to a file stream f in pretty format with variable names at x
*/
void F_mpz_poly_fprint_pretty(const F_mpz_poly_t poly, FILE* f, const char * x);

/** 
   \fn     void F_mpz_poly_print(const F_mpz_poly_t poly)
   \brief  Prints F_mpz_poly_t to screen in standard FLINT format
*/
void F_mpz_poly_print(const F_mpz_poly_t poly);

/** 
   \fn     void F_mpz_poly_print_pretty(const F_mpz_poly_t poly, const char * x)
   \brief  Prints F_mpz_poly_t to screen in pretty format with variable named at x
*/
void F_mpz_poly_print_pretty(const F_mpz_poly_t poly, const char * x);

/** 
   \fn     int F_mpz_poly_fread(F_mpz_poly_t poly, FILE* f)
   \brief  Reads F_mpz_poly_t from file stream f
*/
int F_mpz_poly_fread(F_mpz_poly_t poly, FILE* f);

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
	if (n >= poly->length) _F_mpz_poly_set_length(poly, 0);
	else 
	   for (ulong i = 0; i < n; i++)
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
void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, F_mpz_poly_t poly2, F_mpz_t x);

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
   \fn     F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)

   \brief  Multiply poly1 by poly2 using Kronecker segmentation and store the result in res.
*/
void F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     F_mpz_poly_mul_KS2(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)

   \brief  Multiply poly1 by poly2 using David Harvey's KS2 algorithm and store the result in res.
*/
void F_mpz_poly_mul_KS2(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     F_mpz_poly_mul_SS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result, using the Schoenhage-Strassen 
	        algorithm.
*/
void F_mpz_poly_mul_SS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2);

/** 
   \fn     F_mpz_poly_mul(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
   \brief  Multiply poly1 by poly2 and set res to the result. An attempt is made to choose the 
	        optimal algorithm.
*/
void F_mpz_poly_mul(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2);

/*===============================================================================

   New Naive Standard Functions

================================================================================*/

/**
   Should be decently fast, but I'ld like someone else to confirm.
*/
void F_mpz_poly_scalar_div_exact(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t d);

/**
   Might not be what we want in the end, it calls F_mpz_mod on each coeff then 
      compares to p/2 and subtracts if needed, might be faster to expect pre-mod
*/
void F_mpz_poly_smod(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t p);

/**
   Probably fine as is, it's a derivative I mean come on...
*/
void F_mpz_poly_derivative(F_mpz_poly_t der, F_mpz_poly_t poly);

/**
   Copied from fmpz_poly_content, but allows negative contents which might 
      not be math-correct, but I wanted the factors of a poly * content = poly
*/
void F_mpz_poly_content(F_mpz_t c, const F_mpz_poly_t poly);

/**
   No overflow checks, just does horners method with doubles
*/
double F_mpz_poly_eval_horner_d(F_mpz_poly_t poly, double val);

/**
   Uses mpfs in the middle so theres not saved time by the double exp format
      could use dpe's in the middle for arithmetic or mpfr or something, but 
      this should always work for polynomials no matter how large assuming you 
      don't need the exact output...
*/
double F_mpz_poly_eval_horner_d_2exp(long * exp, F_mpz_poly_t poly, double val);

/**
   Needs a test function, it should be fine, but if the input weren't normalized
      or something like that... maybe...
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

/**
   \fn     void F_mpz_poly_div(F_mpz_poly_t d, F_mpz_poly_t f, F_mpz_poly_t g)
   \brief  Takes f/g and writes to d ignoring any remainder
*/
void F_mpz_poly_div(F_mpz_poly_t d, F_mpz_poly_t f, F_mpz_poly_t g);

/**
   \fn     void F_mpz_poly_divrem(F_mpz_poly_t q, F_mpz_poly_t r, 
               F_mpz_poly_t f, F_mpz_poly_t g)
   \brief  Finds polys r,q such that f = qg+r and deg(r) < deg(g)
*/
void F_mpz_poly_divrem(F_mpz_poly_t q, F_mpz_poly_t r,
               F_mpz_poly_t f, F_mpz_poly_t g);

/*============================================================================

   Naive '_modp' ( := Large moduli ) F_mpz_poly functions

============================================================================*/

/**
   \fn     void F_mpz_poly_rem_modp_naive(F_mpz_poly_t R, F_mpz_poly_t A,
               F_mpz_poly_t B, F_mpz_t p)
   \brief  Finds (A modulo B) modulo p and sets this to R.  Uses mod not smod.
               I call it Naive because it doesn't have special code for the 
               fast cases, namely deg(A) == deg(B) or deg(B) + 1 
               lead coeff of B should be invertible mod p
*/
void F_mpz_poly_rem_modp_naive(F_mpz_poly_t R, F_mpz_poly_t A, F_mpz_poly_t B,
               F_mpz_t p);

/**
   \fn      void F_mpz_poly_mulmod_modp_naive(F_mpz_poly_t R, F_mpz_poly_t f,
                F_mpz_poly_t g, F_mpz_poly_t B, F_mpz_t p)
   \brief   Finds (f*g modulo B) modulo p, just does f*g then rem_modp.
*/
void F_mpz_poly_mulmod_modp_naive(F_mpz_poly_t R, F_mpz_poly_t f, 
                F_mpz_poly_t g, F_mpz_poly_t B, F_mpz_t p);

/**
   \fn      void F_mpz_poly_div_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, 
                F_mpz_poly_t g, F_mpz_t P, ulong n)
   \brief   A power series modp division.  Using only n coeffs of f and g
               finds the lowest n coeffs of f/g assuming the remainder of f/g
               is 0 mod p. (Designed for padic CLDs) Right now if the trailing
               coeff of g is not invertible mod p then it resorts to a full
                division.  Want to find a way around that.
*/
void F_mpz_poly_div_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, 
                F_mpz_poly_t g, F_mpz_t P, ulong n);

/**
   \fn      void F_mpz_poly_div_upper_trunc_modp( F_mpz_t *res, F_mpz_poly_t f,
                F_mpz_poly_t g, F_mpz_t P, ulong n)
   \brief   Reversed power series division modp.  Using only the top n coeffs of 
               f and g finds the top n coeffs of f/g mod p, assuming that the 
               remainder of f/g is zero mod p.  (Designed for padic CLDS) Right 
               now, if the leading coeff of g is not invertible mod p then 
               does full division.  In the intended factoring uses g is monic.
*/
void F_mpz_poly_div_upper_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, 
               F_mpz_poly_t g, F_mpz_t P, ulong n);

/*===========================================================================

   New Material for FLINT, computing fast/tight bounds for CLDs
      CLDs:= Coefficients of Logarithmic Derivatives.  f*g'/g

============================================================================*/

/**
   \fn      int _d_2exp_comp(double a, long ap, double b, long bp)
   \brief   A customized comparison function for a*2^ap vs. b*2^bp, outputs
               -2,-1,1, or 2 for the b pair being two times larger, 
               larger, smaller or =, two times smaller than the a pair
*/
int _d_2exp_comp(double a, long ap, double b, long bp);

/**
   \fn      void F_mpz_poly_CLD_bound(F_mpz_t res, F_mpz_poly_t f, ulong N)
   \brief   A new approach for quickly, tightly, finding a bound on the size of
               the Nth coefficient of f*g'/g for any integer factor of f.  
               A good benchmark for comparison is just f'.
*/
void F_mpz_poly_CLD_bound(F_mpz_t res, F_mpz_poly_t f, ulong N);

/*============================================================================

   Square-Free Factorization

============================================================================*/

/**
   \fn      void F_mpz_poly_squarefree(F_mpz_poly_factor_t fac, 
               F_mpz_t content, F_mpz_poly_t F)
   \brief   Given poly F, finds the content of F which is stores at content,
               then finds a squarefree factorization stored at fac with exponents
*/
void F_mpz_poly_squarefree(F_mpz_poly_factor_t fac, 
               F_mpz_t content, F_mpz_poly_t F);

/*****************************************************************************

   NTL Hensel Lifting procedures would like to streamline and FLINT-ize names

*****************************************************************************/

void _Build_Hensel_Tree(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, zmod_poly_factor_t fac);

void _Hensel_Lift(F_mpz_poly_t Gout, F_mpz_poly_t Hout, F_mpz_poly_t Aout, F_mpz_poly_t Bout, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t h, F_mpz_poly_t a, F_mpz_poly_t b, F_mpz_t p, F_mpz_t p1);

void _Rec_Tree_Hensel_Lift(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, F_mpz_t p, F_mpz_poly_t f, long j, long inv, F_mpz_t p1);

void _Tree_Hensel_Lift(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, long e0, long e1, F_mpz_poly_t f, long inv, long p, long r, F_mpz_t P);

/***************************************************

   Naive Zassenhaus

***************************/

/**
   This guy is unoptimized.  Takes Hensel lifted factors to power P, the original polynomial F (and it's squarefree exponent),
    and for some reason a leading coeff which might not be needed... I'll check later, this is devel stuff here.
*/
void F_mpz_poly_zassenhaus_naive(F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc);

/*********

   Factoring wrapper after square free part

*********/

/*
   This is a wrapper which makes some choices about primes, Hensel lifting, Zassenhaus, this is the wrapper which does the real
   stuff.  Call it after squarefree factoring with your square free f and it's eventual exponent.
*/

void F_mpz_poly_factor_sq_fr_prim( F_mpz_poly_factor_t final_fac, ulong exp, F_mpz_poly_t f );

/*==================================
   Finally the real F_mpz_poly_factor
==================================*/

/**
   Doesn't do much, just does some simple pretests squarefree factors and calls the wrapper.
*/
void F_mpz_poly_factor(F_mpz_poly_factor_t final_fac, F_mpz_t cong, F_mpz_poly_t G);

/*===========

   Hoeij/Novocin approach

=========*/

/*
   Has one sub-optimal part, namely we need F_mpz_poly_mul_trunc and left_trunc.  This computes bounds for and the CLDs themselves.
   Stores them in a matrix res which comes out with r + 1 rows and 2*N or length-1 columns
   (if 2N is smaller then length -1 then just top N and bottom N).  Takes F, Hensel Lifted factors, P and N, in the case we don't choose N large enough 
   we 'could' adapt for restarting... maybe...
*/
void _F_mpz_poly_factor_CLD_mat(F_mpz_mat_t res, F_mpz_poly_t F, F_mpz_poly_factor_t lifted_fac, F_mpz_t P, ulong N);

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
