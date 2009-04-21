/*============================================================================

    F_mpz.h: The FLINT integer format (FLINT 2.0)

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

#ifndef FLINT_F_MPZ_H
#define FLINT_F_MPZ_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"

/* 
   F_mpz_t type
   ============

	The F_mpz_t is a signed integer of FLINT_BITS-2 bits, sign extended to FLINT_BITS bits, unless the most 
	significant bit is zero and the second most significant bit is 1, in which case the bottom FLINT_BITS-2
	bits are an index into the array F_mpz_arr of mpz's.
*/

typedef long F_mpz;
typedef F_mpz F_mpz_t[1];

#define MPZ_BLOCK 16 // number of additional mpz_t's to initialise at a time

// maximum positive value a small coefficient can have
#define COEFF_MAX ((1L<<(FLINT_BITS-2))-1L)

// minimum negative value a small coefficient can have
#define COEFF_MIN (-((1L<<(FLINT_BITS-2))-1L))

// turn an long offset for the F_mpz_arr array into a F_mpz_t style index
#define OFF_TO_COEFF(xxx) ((xxx) | (1L<<(FLINT_BITS - 2))) 

// returns the F_mpz_t style index for the F_mpz_arr array as a long offset
#define COEFF_TO_OFF(xxx) ((xxx) & ((1L<<(FLINT_BITS - 2))-1)) 

#define COEFF_IS_MPZ(xxx) ((xxx>>(FLINT_BITS-2)) == 1L) // is xxx an index into F_mpz_arr?

static gmp_randstate_t F_mpz_state; // Used for random generation in F_mpz_randomm only

/*===============================================================================

	mpz_t memory management

================================================================================*/
 
/** 
   \fn     F_mpz_t _F_mpz_new_mpz(void)
   \brief  Return a new mpz F_mpz_t. The mpz_t's are allocated and initialised
	        in blocks of size MPZ_BLOCK. 
*/
F_mpz _F_mpz_new_mpz(void);

/** 
   \fn     void _F_mpz_clear_mpz(F_mpz_t f)
   \brief  Release the mpz associated to f to the array of unused mpz's. 
	        Assumes f actually represents an mpz.
*/
void _F_mpz_clear_mpz(F_mpz f);

/** 
   \fn     void _F_mpz_cleanup(void)
   \brief  Clear any mpz's still held onto by the F_mpz_t memory management.
*/
void _F_mpz_cleanup(void);

/*===============================================================================

	Promotion/Demotion

================================================================================*/

/** 
   \fn     __mpz_struct * _F_mpz_promote(F_mpz_t f)
   \brief  Promote the F_mpz_t to an mpz_t. No assumption is made about f. The value of f 
	is not preserved. A pointer to an __mpz_struct corresponding to f is returned.
*/
__mpz_struct * _F_mpz_promote(F_mpz_t f);

/** 
   \fn     __mpz_struct * _F_mpz_promote_val(F_mpz_t f)
   \brief  Promote the given F_mpz_t to an mpz_t. No assumption is made about f. The value
	        of f is preserved. A pointer to an __mpz_struct corresponding to f is returned.
*/
__mpz_struct * _F_mpz_promote_val(F_mpz_t f);

/** 
   \fn     void _F_mpz_demote(F_mpz_t f)
   \brief  If f represents an mpz_t then the mpz_t is released. Makes no assumptions about
	        f. Note that f is not set to any value, i.e. it must be set to a small integer 
			  immediately after calling this function!
*/
static inline
void _F_mpz_demote(F_mpz_t f)
{
	if (COEFF_IS_MPZ(*f)) 
	{
		_F_mpz_clear_mpz(*f);
		*f = 0L;
	}
}

/** 
   \fn     void _F_mpz_demote_val(F_mpz_t f)
   \brief  If the F_mpz_t (which is assumed to be an mpz_t) will fit into FLINT_BIT - 2 bits, 
	        it is demoted to a limb instead of an mpz_t, preserving the value, otherwise 
			  nothing happens.
*/
void _F_mpz_demote_val(F_mpz_t f);

/*===============================================================================

	F_mpz_t memory management

================================================================================*/

/** 
   \fn     void F_mpz_init(F_mpz_t f)
   \brief  Allocate an F_mpz_t. A small F_mpz_t is returned (i.e. not 
	        representing an mpz_t).
*/
static inline
void F_mpz_init(F_mpz_t f)
{
	*f = 0L;
}

/** 
   \fn     void F_mpz_init2(F_mpz_t f, ulong limbs)
   \brief  Allocate an F_mpz_t with the given number of limbs. If limbs
	        is zero then a small F_mpz_t is returned (i.e. not representing
			  an mpz_t).
*/
void F_mpz_init2(F_mpz_t f, ulong limbs);

/** 
   \fn     void F_mpz_clear(F_mpz_t f)
   \brief  Clear the given F_mpz_t.
*/
static inline
void F_mpz_clear(F_mpz_t f)
{
   _F_mpz_demote(f);
}

/*===============================================================================

	Random generation

================================================================================*/

/** 
   \fn     void F_mpz_random(F_mpz_t f, const ulong bits)
   \brief  Generate a random F_mpz_t with the given number of bits.
*/
void F_mpz_random(F_mpz_t f, const ulong bits);

/** 
   \fn     void F_mpz_random(F_mpz_t f, const ulong bits)
   \brief  Generate a random F_mpz_t in [0, n) where n is an mpz_t.
*/
void F_mpz_randomm(F_mpz_t f, const mpz_t n);

/*===============================================================================

	Get/set

================================================================================*/

/** 
   \fn     void F_mpz_zero(F_mpz_t f)
   \brief  Set the given F_mpz_t to zero.
*/
static inline
void F_mpz_zero(F_mpz_t f)
{
	_F_mpz_demote(f);
   *f = 0L;
}

/** 
   \fn     void F_mpz_set_si(F_mpz_t f, const long val)
   \brief  Set f to a signed long value val
*/
void F_mpz_set_si(F_mpz_t f, const long val);

/** 
   \fn     void F_mpz_set_ui(F_mpz_t f, const ulong val)
   \brief  Set f to an unsigned long value val
*/
void F_mpz_set_ui(F_mpz_t f, const ulong val);

/** 
   \fn     long F_mpz_get_si(const F_mpz_t f)
   \brief  Return the value of f as a long
*/
long F_mpz_get_si(const F_mpz_t f);

/** 
   \fn     long F_mpz_get_ui(const F_mpz_t f)
   \brief  Return the value of f as an unsigned long
*/
ulong F_mpz_get_ui(const F_mpz_t f);

/** 
   \fn     void F_mpz_get_mpz(mpz_t x, const F_mpz_t f)
   \brief  Returns f as an mpz_t
*/
void F_mpz_get_mpz(mpz_t x, const F_mpz_t f);

/** 
   \fn     double F_mpz_get_d_2exp(long * exp, const F_mpz_t f)
   \brief  Return f as a signed normalised double and a long exponent.
*/
double F_mpz_get_d_2exp(long * exp, const F_mpz_t f);

/** 
   \fn     void F_mpz_set_mpz(F_mpz_t f, const mpz_t x)
   \brief  Sets f to the given mpz_t
*/
void F_mpz_set_mpz(F_mpz_t f, const mpz_t x);

/** 
   \fn     void F_mpz_set_limbs(F_mpz_t f, const mp_limb_t * x, const ulong limbs)
   \brief  Sets f to the array of limbs x which is the given number of 
	        limbs in length and where the least significant limb is 
			  stored first in x.
*/
void F_mpz_set_limbs(F_mpz_t f, const mp_limb_t * x, const ulong limbs);

/** 
   \fn     ulong F_mpz_set_limbs(const mp_limb_t * x, F_mpz_t f)
   \brief  Sets the array of limbs x to the absolute value of f. The 
	        array is assumed to be stored with least significant limb 
			  first. The number of limbs written is returned.
*/
ulong F_mpz_get_limbs(mp_limb_t * x, const F_mpz_t f);

/** 
   \fn     void F_mpz_set(F_mpz_t f, F_mpz_t g)
   \brief  Sets f to the value of g. 
*/
void F_mpz_set(F_mpz_t f, const F_mpz_t g);

/** 
   \fn     void F_mpz_swap(F_mpz_t f, F_mpz_t g)
   \brief  Efficiently swaps the two F_mpz_t's. 
*/
void F_mpz_swap(F_mpz_t f, F_mpz_t g);

/*===============================================================================

	Comparison

================================================================================*/

/** 
   \fn     int F_mpz_equal(const F_mpz_t f, const F_mpz_t g)
   \brief  Returns 1 if the two values are equal, otherwise returns 0.
*/
int F_mpz_equal(const F_mpz_t f, const F_mpz_t g);

/** 
   \fn     int F_mpz_cmpabs(const F_mpz_t f, const F_mpz_t g)
   \brief  Returns a negative int if abs(f) < abs(g), positive if 
	        abs(f) > abs(g) and returns 0 if the two values are equal.
*/
int F_mpz_cmpabs(const F_mpz_t f, const F_mpz_t g);

/*===============================================================================

	Properties

================================================================================*/

/** 
   \fn     ulong F_mpz_entry_size(F_mpz_t f)
   \brief  Returns the number of limbs required to store the absolute value of f.
	        Returns 0 if f is zero.
*/
ulong F_mpz_size(const F_mpz_t f);

/** 
   \fn     int F_mpz_sgn(const F_mpz_t f)
   \brief  Returns 1 if f is positive, -1 if it is negative and 0 if f is zero.
*/
int F_mpz_sgn(const F_mpz_t f);

/** 
   \fn     int F_mpz_is_zero(const F_mpz_t f)
   \brief  Returns 1 if f is zero, 0 otherwise.
*/
static inline
int F_mpz_is_zero(const F_mpz_t f)
{
	return ((*f) == 0L);
}

/** 
   \fn     ulong F_mpz_bits(F_mpz_t f)
   \brief  Returns the number of bits required to store the absolute value of f.
	        Returns 0 if f is zero.
*/
ulong F_mpz_bits(const F_mpz_t f);

/** 
   \fn     __mpz_struct * F_mpz_ptr_mpz(F_mpz f)
   \brief  Returns a pointer to the mpz_t associated with the coefficient f.
	        Assumes f is actually associated with an mpz_t.
*/
__mpz_struct * F_mpz_ptr_mpz(const F_mpz f);

/*===============================================================================

	Input/output

================================================================================*/

/** 
   \fn     void F_mpz_print(F_mpz_t x)
   \brief  Print the given F_mpz_t to stdout.
*/
static inline 
void F_mpz_print(F_mpz_t x)
{
	if (!COEFF_IS_MPZ(*x)) printf("%ld", *x);
	else gmp_printf("%Zd", F_mpz_ptr_mpz(*x));
}

/** 
   \fn     void F_mpz_print(F_mpz_t x)
   \brief  Read an F_mpz_t from stdin. The integer can be a signed multiprecision
	        integer in decimal format.
*/
void F_mpz_read(F_mpz_t f);


/*===============================================================================

	Arithmetic

================================================================================*/

/** 
   \fn     void F_mpz_neg(F_mpz_t f, F_mpz_t g)
   \brief  Sets f to minus g. 
*/
void F_mpz_neg(F_mpz_t f, const F_mpz_t g);

/** 
   \fn     void F_mpz_add(F_mpz_t f, const F_mpz_t g, mpz_t h)
   \brief  Set f to g plus h, where h is an mpz_t. 
*/
void F_mpz_add_mpz(F_mpz_t f, const F_mpz_t g, mpz_t h);

/** 
   \fn     void F_mpz_add(F_mpz_t f, const F_mpz_t g, F_mpz_t h)
   \brief  Set f to g plus h. 
*/
void F_mpz_add(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_sub(F_mpz_t f, const F_mpz_t g, F_mpz_t h)
   \brief  Set f to g minus h. 
*/
void F_mpz_sub(F_mpz_t f, const F_mpz_t g, F_mpz_t h);

/** 
   \fn     void F_mpz_mul_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
   \brief  Multiply g by the unsigned long x and set f to the result.
*/
void F_mpz_mul_ui(F_mpz_t f, const F_mpz_t g, const ulong x);

/** 
   \fn     void F_mpz_mul_si(F_mpz_t f, const F_mpz_t g, const long x)
   \brief  Multiply g by the signed long x and f to the result.
*/
void F_mpz_mul_si(F_mpz_t f, const F_mpz_t g, const long x);

/** 
   \fn     void F_mpz_mul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Multiply g by h and set f to the result.
*/
void F_mpz_mul2(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_mul_2exp(F_mpz_t f, const F_mpz_t g, const ulong exp)
   \brief  Multiply g by 2^exp and set f to the result.
*/
void F_mpz_mul_2exp(F_mpz_t f, const F_mpz_t g, const ulong exp);

/** 
   \fn     void F_mpz_div_2exp(F_mpz_t f, const F_mpz_t g, const ulong exp)
   \brief  Divide g by 2^exp and set f to the result. Rounding is towards zero.
*/
void F_mpz_div_2exp(F_mpz_t f, const F_mpz_t g, const ulong exp);

/** 
   \fn     void F_mpz_add_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
   \brief  Add the unsigned long x to g and set f to the result.
*/
void F_mpz_add_ui(F_mpz_t f, const F_mpz_t g, const ulong x);

/** 
   \fn     void F_mpz_sub_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
   \brief  Subtract the unsigned long x from g and set f to the result.
*/
void F_mpz_sub_ui(F_mpz_t f, const F_mpz_t g, const ulong x);

/** 
   \fn     void F_mpz_addmul_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
   \brief  Multiply g by the unsigned long x and add the result to f, in place.
*/
void F_mpz_addmul_ui(F_mpz_t f, const F_mpz_t g, const ulong x);

/** 
   \fn     void F_mpz_submul_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
   \brief  Multiply g by the unsigned long x and subtract the result from f, in place.
*/
void F_mpz_submul_ui(F_mpz_t f, const F_mpz_t g, const ulong x);

/** 
   \fn     void F_mpz_addmul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Multiply g by h and add the result to f, in place.
*/
void F_mpz_addmul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_submul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Multiply g by h and subtract the result from f, in place.
*/
void F_mpz_submul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_mod(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Set f to g modulo h.
*/
void F_mpz_mod(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_divexact(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Set f to g divided by h, assuming the division is exact.
*/
void F_mpz_divexact(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_fdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Set f to g divided by h, rounded down towards minus infinity.
*/
void F_mpz_fdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_cdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Set f to g divided by h, rounded up towards infinity.
*/
void F_mpz_cdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

/** 
   \fn     void F_mpz_rdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
   \brief  Set f to g divided by h, rounded to nearest, ties rounded towards
	        positive infinity.
*/
void F_mpz_rdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h);

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
