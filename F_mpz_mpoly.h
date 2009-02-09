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

==============================================================================*/

/*
    F_mpz_mpoly.h: Multivariable polynomials over Z (FLINT 2.0 polynomials)

    Copyright (C) 2008, 2009 William Hart 
*/

#ifndef FLINT_F_MPZ_MPOLY_H
#define FLINT_F_MPZ_MPOLY_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "flint.h"
#include "F_mpz.h"
#include "packed_vec.h"

/*==============================================================================

   F_mpz_mpoly_t
   -----------

F_mpz_mpoly_t represents a sparse polynomial in Z[x1, x2, ...., xn] 

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

"vars" is the number of variables the polynomial has.

"exps" is an array of packed vectors (pv_s structs) which represent exponents,
each vector corresponding to a different variable. The structure of each array
is intended to be a black box and accessible only via functions in packed_vec.

"packed" is an array of entries, one ulong for each term of the polynomial. 
Each one is a packed expression which reflects the monomial ordering in that 
if packed[i] < packed[j] then monomial i < monomial j. The information may 
be incomplete, in which case we may have packed[i] == packed[j] even though 
monomial i < monomial j.

"ordering" specifies which polynomial ordering is being used.

"small" is set to true if all the exponent information is stored in "packed".
In this case no monomial exponent information is held in exps.

================================================================================*/
 
typedef enum
{
	LEX,
	REVLEX,
	GRLEX,
	GREVLEX
} ord_t;

typedef struct
{
   F_mpz * coeffs;
   ulong alloc;
   ulong length;
	ulong vars;
	pv_s * exps;
	ulong * packed;
	ord_t ordering;
	int packed_bits;
	int small;
} F_mpz_mpoly_struct;

// F_mpz_poly_t allows reference-like semantics for F_mpz_poly_struct
typedef F_mpz_mpoly_struct F_mpz_mpoly_t[1];

typedef struct F_mpz_mpoly_heap_s
{
   ulong packed;
	ulong j;
	ulong chain;
} F_mpz_mpoly_heap_s;

/*===============================================================================

	Memory management

================================================================================*/

/** 
   \fn     void F_mpz_mpoly_init(F_mpz_mpoly_t poly, ord_t ordering)
   \brief  Initialise a polynomial of length zero with zero allocated coefficients
*/
void F_mpz_mpoly_init(F_mpz_mpoly_t poly, ord_t ordering);

/** 
   \fn     void F_mpz_mpoly_init2(F_mpz_mpoly_t poly, const ulong alloc, 
	                                      const ulong vars, ord_t ordering)
   \brief  Initialise a polynomial of length zero with the given number of allocated 
	        coefficients
*/
void F_mpz_mpoly_init2(F_mpz_mpoly_t poly, const ulong alloc, 
							                    const ulong vars, ord_t ordering);

/** 
   \fn     void F_mpz_mpoly_realloc(F_mpz_mpoly_t poly, const ulong alloc, 
	                                                                const ulong vars)
   \brief  Reallocates poly to have space for precisely the given number of 
	        coefficients. If alloc = 0, then the polynomial is cleared. If
			  the alloc is smaller than the current length of the polynomial
			  the polynomial is truncated and normalised.
*/
void F_mpz_mpoly_realloc(F_mpz_mpoly_t poly, const ulong alloc, const ulong vars);

/** 
   \fn     void F_mpz_mpoly_fit_length(F_mpz_mpoly_t poly, const ulong length)
   \brief  Expands poly, if necessary, so that it has space for the given number
	        of coefficients. This function never shrinks the polynomial, only 
			  expands it.
*/
void F_mpz_mpoly_fit_length(F_mpz_mpoly_t poly, const ulong length);

/** 
   \fn     void F_mpz_mpoly_fit_vars(F_mpz_mpoly_t poly, const ulong vars)
   \brief  Expands poly, if necessary, so that it has space for the given number
	        of variables. This function never shrinks the polynomial, only 
			  expands it.
*/
static inline
void F_mpz_mpoly_fit_vars(F_mpz_mpoly_t poly, const ulong vars)
{
	if (vars > poly->vars) F_mpz_mpoly_realloc(poly, poly->alloc, vars);
}

/** 
   \fn     void F_mpz_mpoly_clear(F_mpz_mpoly_t poly)
   \brief  Clear the polynomial, releasing any memory it was using.
*/
void F_mpz_mpoly_clear(F_mpz_mpoly_t poly);

/*===============================================================================

	Truncation

================================================================================*/

/** 
   \fn     void F_mpz_mpoly_truncate(F_mpz_mpoly_t poly, const ulong length)
   \brief  Truncate the polynomial to the given length. It is permissible for
	        length to be greater than the current length of the polynomial, in 
			  which case nothing will happen.
*/
static inline
void _F_mpz_mpoly_truncate(F_mpz_mpoly_t poly, const ulong length)
{
	if (poly->length > length) // only truncate if necessary
   {
      for (ulong i = length; i < poly->length; i++)
			_F_mpz_demote(poly->coeffs + i);
		poly->length = length;
   }  
}

/*===============================================================================

	Set/get

================================================================================*/

/** 
   \fn     void F_mpz_mpoly_set_coeff_ui(F_mpz_mpoly_t poly, const ulong n, 
	                                                               const ulong x)
   \brief  Set coefficient with index n in poly to the given unsigned long. 
	        Assumes that n is not bigger than the length of poly and that x is not 0.
*/
void F_mpz_mpoly_set_coeff_ui(F_mpz_mpoly_t poly, const ulong n, const ulong x);

/** 
   \fn     void F_mpz_mpoly_get_coeff_ui(F_mpz_mpoly_t poly, ulong n, const long x)
   \brief  Return coefficient n of poly as an usigned long. If n is greater or equal to
	        the length of poly, then zero is returned.
*/
ulong F_mpz_mpoly_get_coeff_ui(const F_mpz_mpoly_t poly, const ulong n);

/** 
   \fn     void F_mpz_mpoly_set_var_exp_grlex(F_mpz_mpoly_t poly, const ulong n, 
	                                                const ulong var, const ulong exp)
   \brief  Assumes the polynomial is small, i.e. all exponent data can
	        be packed into poly->packed. The function sets the exponent of
			  the variable with index var corresponding to monomial with index
			  n in poly to the given exponent.
*/
void F_mpz_mpoly_set_var_exp_packed_grlex(F_mpz_mpoly_t poly, 
								const ulong n, const ulong var, const ulong exp);
/** 
   \fn     void F_mpz_mpoly_set_var_exp(F_mpz_mpoly_t poly, const ulong n, 
	                                                const ulong var, const ulong exp)
   \brief  Set the exponent of the variable with index var corresponding to the 
	        monomial with index n in poly to the given exponent.
*/
void F_mpz_mpoly_set_var_exp(F_mpz_mpoly_t poly, const ulong n, 
									                        const ulong var, const ulong exp);

/** 
   \fn     ulong F_mpz_mpoly_get_var_exp_packed_grlex(F_mpz_mpoly_t poly, const ulong n, 
	                                                                 const ulong var)
   \brief  Get the exponent of the variable with index var corresponding to the 
	        monomial with index n in poly, assuming poly is small with packed
			  monomial exponents and that the ordering is GRLEX.
*/
ulong F_mpz_mpoly_get_var_exp_packed_grlex(F_mpz_mpoly_t poly, const ulong n, const ulong var);

/** 
   \fn     ulong F_mpz_mpoly_get_var_exp(F_mpz_mpoly_t poly, const ulong n, 
	                                                                 const ulong var)
   \brief  Get the exponent of the variable with index var corresponding to the 
	        monomial with index n in poly.
*/
ulong F_mpz_mpoly_get_var_exp(F_mpz_mpoly_t poly, const ulong n, const ulong var);

/*===============================================================================

	Print/read

================================================================================*/

/** 
   \fn     void F_mpz_mpoly_print_pretty(F_mpz_mpoly poly, char ** var_syms)

   \brief  Print the given polynomial taking an array of strings representing
	        the variables of the polynomial.
*/
void F_mpz_mpoly_print_pretty(F_mpz_mpoly_t poly, char ** var_syms);

/*===============================================================================

	Arithmetic

================================================================================*/

/** 
   \fn     void F_mpz_mpoly_add_inplace(F_mpz_mpoly res, ulong n, F_mpz_mpoly poly)

   \brief  Add poly to res, starting at coefficient with index n in res,
	        storing the result in res.
*/
void F_mpz_mpoly_add_inplace(F_mpz_mpoly_t res, ulong n, F_mpz_mpoly_t poly);

/** 
   \fn     void _F_mpz_mpoly_mul_small_mxn(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)

   \brief  Multiply poly1 by poly2, assuming both have length at most 4 and
	        are small. No memory management is done for the result res. The
			  input polynomials are only used from the coefficients with indices
			  n1 and n2 onwards, respectively.
*/
void _F_mpz_mpoly_mul_small_mxn(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2);

/** 
   \fn     void _F_mpz_mpoly_mul_small_Mxn(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								           F_mpz_mpoly_t poly2, ulong n2, ulong n)

   \brief  Multiply poly1 by the n coefficients of poly2 starting at the 
	        coefficient of poly2 with index n2, assuming both polynomials
	        are small. Assumes n is at most 4. No memory management is 
			  done for the result res. 
*/
void _F_mpz_mpoly_mul_small_Mxn(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								           F_mpz_mpoly_t poly2, ulong n2, ulong n);

/** 
   \fn     void F_mpz_mpoly_mul_small_recursive(F_mpz_mpoly_t res, 
				  F_mpz_mpoly_t poly1, F_mpz_mpoly_t poly2, ulong n2, ulong n)

   \brief  Multiply poly1 by the n coefficients of poly2 starting at the 
	        coefficient of poly2 with index n2, assuming both polynomials
	        are small. 
*/
void F_mpz_mpoly_mul_small_recursive(F_mpz_mpoly_t res, 
				  F_mpz_mpoly_t poly1, F_mpz_mpoly_t poly2, ulong n2, ulong n);

/** 
   \fn     void F_mpz_mpoly_mul_small_heap(F_mpz_mpoly_t res, 
				                    F_mpz_mpoly_t poly1, F_mpz_mpoly_t poly2)

   \brief  Multiply poly1 by poly2, assuming both polynomials are small,
	        using a heap algorithm. 
*/
void F_mpz_mpoly_mul_small_heap(F_mpz_mpoly_t res, 
				                    F_mpz_mpoly_t poly1, F_mpz_mpoly_t poly2);


#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
