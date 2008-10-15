/*============================================================================

    F_mpz_mat.h: Matrices over Z (FLINT 2.0 matrices)

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

#ifndef FLINT_F_MPZ_MAT_H
#define FLINT_F_MPZ_MAT_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "flint.h"
#include "mpz_mat.h"

/*==============================================================================

   F_mpz_mat_t
   -----------

F_mpz_mat_t represents a dense matrix in M(Z) 

"entries" is an array of limbs, one for each entry

There are two things each entry in this array can represent:

1) If the most significant two bits are 01, then the entry represents
an index into the array "mpz_entries", and the mpz_t in that array
contains the entry.

2) Otherwise, the entry represents a signed entry
whose absolute value is no more than FLINT_BIT - 2 bits in length. The
entry is stored in twos complement format.

"mpz_entries" is an array of entries in mpz_t format (actually an array of
__mpz_struct's). Only entries whose absolute value does not fit into 
FLINT_BITS - 2 bits are stored in an mpz_t.

"alloc" is the number of allocated entries. It is equal to r*c.

"r" is the number of rows of the matrix

"c" is the number of columns

"mpz_alloc" is the number of allocated mpz_entries

"mpz_length" is the number of mpz_t's actually used

"rows" is an array of pointers to the starts of rows of mat

================================================================================*/
 
typedef struct
{
   mp_limb_t * entries;
   __mpz_struct * mpz_entries;
   ulong alloc;
   ulong r;
	ulong c;
	mp_limb_t ** rows;
   ulong mpz_alloc;
   ulong mpz_length;
} F_mpz_mat_struct;

// fmpz_mat_t allows reference-like semantics for F_mpz_mat_struct
typedef F_mpz_mat_struct F_mpz_mat_t[1];

#define MPZ_BLOCK_MAT 16 // number of additional mpz_t's to initialise at a time

// maximum positive value a small entry can have
#define ENTRY_MAX ((1L<<(FLINT_BITS-2))-1L)

// minimum negative value a small entry can have
#define ENTRY_MIN (-((1L<<(FLINT_BITS-2))-1L))

// turn an integer offset for the mpz_entry array into an F_mpz_mat_t entry
#define OFF_TO_ENTRY(xxx) ((xxx) | (1L<<(FLINT_BITS - 2))) 

// returns the index as an integer
#define ENTRY_TO_OFF(xxx) ((xxx) & ((1L<<(FLINT_BITS - 2))-1)) 

#define ENTRY_IS_MPZ(xxx) ((xxx>>(FLINT_BITS-2)) == 1L) // is xxx an index into mpz_entries?

/*===============================================================================

	Memory management

================================================================================*/

/** 
   \fn     void _F_mpz_mat_mpz_entries_new(F_mpz_mat_t mat)
   \brief  Add a new mpz_entry. The mpz_t's are allocated and initialised
	        in blocks size MPZ_BLOCK_MAT. mat->mpz_length is incremented and 
			  mat->mpz_alloc is updated if new mpz_t entries were added.
*/
void _F_mpz_mat_mpz_entries_new(F_mpz_mat_t mat);

/** 
   \fn     void _F_mpz_mat_mpz_entries_clear(F_mpz_mat_t mat)
   \brief  Clear mat's array of mpz_t entries.
*/
void _F_mpz_mat_mpz_entries_clear(F_mpz_mat_t mat);

/** 
   \fn     void F_mpz_mat_init(F_mpz_mat_t mat, ulong r, ulong c)
   \brief  Initialise a matrix with r rows and c columns.
*/
void F_mpz_mat_init(F_mpz_mat_t mat, ulong r, ulong c);

/** 
   \fn     void F_mpz_mat_clear(F_mpz_mat_t mat)
   \brief  Clear the matrix, releasing any memory it was using.
*/
void F_mpz_mat_clear(F_mpz_mat_t mat);

/*===============================================================================

	Coefficient operations

================================================================================*/

/** 
   \fn     __mpz_struct * _F_mpz_entry_promote(F_mpz_mat_t const mat, ulong entry);
   \brief  Promote the given entry of mat to an mpz_t entry. The value
	        of the entry is not preserved. A pointer to an __mpz_struct
			  corresponding to the entry, is returned.
*/
static inline
__mpz_struct * _F_mpz_entry_promote(F_mpz_mat_t mat, const ulong r, const ulong c)
{
   ulong d = mat->rows[r][c];
	if (!ENTRY_IS_MPZ(d)) // entry is small so promote it
	{
	   _F_mpz_mat_mpz_entries_new(mat);
	   mat->rows[r][c] = OFF_TO_ENTRY(mat->mpz_length - 1);
		return mat->mpz_entries + mat->mpz_length - 1;
	} else // entry is large already, just return the pointer
      return mat->mpz_entries + ENTRY_TO_OFF(d);
}

/** 
   \fn     __mpz_struct * _F_mpz_entry_promote_val(F_mpz_mat_t const mat, ulong entry);
   \brief  Promote the given entry of mat to an mpz_t entry. The value
	        of the entry is preserved. A pointer to an __mpz_struct
			  corresponding to the entry, is returned.
*/
static inline
__mpz_struct * _F_mpz_entry_promote_val(F_mpz_mat_t mat, const ulong r, const ulong c)
{
   ulong d = mat->rows[r][c];
	if (!ENTRY_IS_MPZ(d)) // entry is small so promote it
	{
	   _F_mpz_mat_mpz_entries_new(mat);
	   mat->rows[r][c] = OFF_TO_ENTRY(mat->mpz_length - 1);
		__mpz_struct * mpz_ptr = mat->mpz_entries + mat->mpz_length - 1;
		mpz_set_si(mpz_ptr, d);
		return mpz_ptr;
	} else // entry is large already, just return the pointer
      return mat->mpz_entries + ENTRY_TO_OFF(c);
}

/** 
   \fn     _F_mpz_entry_demote_val(F_mpz_mat_t mat, const ulong entry);
   \brief  If the given entry (which is assumed to be an mpz_t) will fit into
	        FLINT_BIT - 2 bits, it is demoted to a limb instead of an mpz_t, preserving
			  the value, otherwise nothing happens.
*/
void _F_mpz_entry_demote_val(F_mpz_mat_t mat, const ulong r, const ulong c);

/** 
   \fn     void _F_mpz_entry_zero(F_mpz_mat_t mat, const ulong entry)
   \brief  Set the given entry to zero
*/
static inline
void _F_mpz_entry_zero(F_mpz_mat_t mat, const ulong r, const ulong c)
{
	mat->rows[r][c] = 0;
}

/** 
   \fn     _F_mpz_entry_set_si(F_mpz_mat_t mat, ulong entry, const long val)
   \brief  Set the given entry of mat to a signed long value val
*/
void _F_mpz_entry_set_si(F_mpz_mat_t mat, const ulong r, const ulong c, const long val);

/** 
   \fn     long _F_mpz_entry_get_si(const F_mpz_mat_t mat, const ulong entry)
   \brief  Return the value of the given entry of mat as a long
*/
long _F_mpz_entry_get_si(const F_mpz_mat_t mat, const ulong r, const ulong c);

/** 
   \fn     _F_mpz_entry_set_ui(F_mpz_mat_t mat, ulong entry, const long val)
   \brief  Set the given entry of mat to an unsigned long value val
*/
void _F_mpz_entry_set_ui(F_mpz_mat_t mat, const ulong r, const ulong c, const ulong val);

/** 
   \fn     long _F_mpz_entry_get_ui(const F_mpz_mat_t mat, const ulong entry)
   \brief  Return the value of the given entry of mat as an unsigned long
*/
ulong _F_mpz_entry_get_ui(const F_mpz_mat_t mat, const ulong r, const ulong c);

/** 
   \fn     _F_mpz_entry_get_mpz(mpz_t x, const F_mpz_mat_t mat, const ulong entry)
   \brief  Returns the given entry of mat as an mpz_t
*/
void _F_mpz_entry_get_mpz(mpz_t x, const F_mpz_mat_t mat, const ulong r, const ulong c);

/** 
   \fn     _F_mpz_entry_set_mpz(F_mpz_mat_t mat, ulong entry, const mpz_t x)
   \brief  Sets the given entry to the given mpz_t
*/
void _F_mpz_entry_set_mpz(F_mpz_mat_t mat, const ulong r, const ulong c, const mpz_t x);

/** 
   \fn     void _F_mpz_entry_set(F_mpz_mat_t mat1, ulong entry1, const F_mpz_mat_t mat2, const ulong entry2)
   \brief  Sets entry1 of mat1 to equal entry2 of mat2. Assumes the entries are distinct.
*/
void _F_mpz_entry_set(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, const ulong r2, const ulong c2);

/** 
   \fn     int _F_mpz_entry_equal(const F_mpz_mat_t mat1, const ulong entry1, const F_mpz_mat_t mat2, const ulong entry2)
   \brief  Returns 1 if entry1 of mat1 is equal to entry2 of mat2, otherwise returns 0.
*/
int _F_mpz_entry_equal(const F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, const ulong r2, const ulong c2);

/** 
   \fn     void _F_mpz_entry_swap(F_mpz_mat_t mat1, ulong entry1, F_mpz_mat_t mat2, ulong entry2)
   \brief  Efficiently swaps entry1 of mat1 and entry2 of mat2. 
*/
void _F_mpz_entry_swap(F_mpz_mat_t mat1, const ulong r1, const ulong c1, F_mpz_mat_t mat2, const ulong r2, const ulong c2);

/** 
   \fn     void _F_mpz_entry_negate(F_mpz_mat_t mat1, ulong entry1, const F_mpz_mat_t mat2, const ulong entry2)
   \brief  Sets entry1 of mat1 to minus entry2 of mat2. Assumes the entries are distinct.
*/
void _F_mpz_entry_neg(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, const ulong r2, const ulong c2);

/** 
   \fn     void _F_mpz_entry_add(F_mpz_mat_t res, const ulong r3, const ulong c3, const F_mpz_mat_t mat1, const ulong r1, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r2, const ulong c2) 

   \brief  Add the given entries of mat1 and mat2 and set the given entry 
	        of res to the result. Assumes the entry of res is distinct from the 
			  other two entries. 
*/
void _F_mpz_entry_add(F_mpz_mat_t res, const ulong r3, const ulong c3, const F_mpz_mat_t mat1, const ulong r1, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r2, const ulong c2); 

/** 
   \fn     void _F_mpz_entry_sub(F_mpz_mat_t res, const ulong r3, const ulong c3, const F_mpz_mat_t mat1, const ulong r1, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r2, const ulong c2) 
   \brief  Subtract the given entries of mat1 and mat2 and set the given entry 
	        of res to the result. Assumes the entry of res is distinct from the 
			  other two entries.
*/
void _F_mpz_entry_sub(F_mpz_mat_t res, const ulong r3, const ulong c3, const F_mpz_mat_t mat1, const ulong r1, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r2, const ulong c2); 

/** 
   \fn     void _F_mpz_entry_mul_ui(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, 
						                                   const ulong r2, const ulong c2, const ulong x)
   \brief  Multiply entry2 of mat2 by the unsigned long x and set entry1 of mat1 to the result.
*/
void _F_mpz_entry_mul_ui(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, 
						                                   const ulong r2, const ulong c2, const ulong x);

/** 
   \fn     void _F_mpz_entry_mul_si(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, 
						                                   const ulong r2, const ulong c2, const long x)
   \brief  Multiply entry2 of mat2 by the signed long x and set entry1 of mat1 to the result.
*/
void _F_mpz_entry_mul_si(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, 
						                                   const ulong r2, const ulong c2, const long x);

/** 
   \fn     void _F_mpz_entry_mul_mpz(F_mpz_mat_t mat1, ulong entry1, F_mpz_mat_t mat2, ulong entry2, mpz_t x)
   \brief  Multiply entry2 of mat2 by the mpz_t x and set entry1 of mat1 to the result.
*/
//void _F_mpz_entry_mul_mpz(F_mpz_mat_t mat1, ulong entry1, F_mpz_mat_t mat2, ulong entry2, mpz_t x);

/** 
   \fn     void _F_mpz_entry_mul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong entry1, 
					                                            const F_mpz_mat_t mat2, const ulong entry2)

   \brief  Multiply entry1 of mat1 by entry2 of mat2 and set entry3 of res to the result.
*/
//void _F_mpz_entry_mul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong entry1, 
//					                                 const F_mpz_mat_t mat2, const ulong entry2);

/** 
   \fn     void _F_mpz_entry_mul_2exp(F_mpz_mat_t mat2, ulong r2, ulong c2, const F_mpz_mat_t mat1, 
									                  const ulong r1, const ulong c1, const ulong exp)
   \brief  Multiply entry1 of mat1 by 2^exp and set entry2 of mat2 to the result.
*/

void _F_mpz_entry_mul_2exp(F_mpz_mat_t mat2, ulong r2, ulong c2, const F_mpz_mat_t mat1, 
									                  const ulong r1, const ulong c1, const ulong exp);
/** 
   \fn     void _F_mpz_entry_add_ui_inplace(F_mpz_mat_t res, ulong entry, const ulong x)
   \brief  Add the unsigned long x to the given entry of res, in place.
*/
//void _F_mpz_entry_add_ui_inplace(F_mpz_mat_t res, ulong entry, const ulong x);

/** 
   \fn     void _F_mpz_entry_sub_ui_inplace(F_mpz_mat_t res, ulong entry, const ulong x)
   \brief  Subtract the unsigned long x from the given entry of res, in place.
*/
//void _F_mpz_entry_sub_ui_inplace(F_mpz_mat_t res, ulong entry, const ulong x);

/** 
   \fn     void _F_mpz_entry_addmul_ui(F_mpz_mat_t res, ulong entry2, const F_mpz_mat_t mat1, 
							                                            const ulong entry1, const ulong x)
   \brief  Multiply entry1 of mat1 by the unsigned long x and add the result to res, in place.
*/
//void _F_mpz_entry_addmul_ui(F_mpz_mat_t res, ulong entry2, const F_mpz_mat_t mat1, 
//							                                 const ulong entry1, const ulong x);

/** 
   \fn     void _F_mpz_entry_submul_ui(F_mpz_mat_t res, ulong entry2, const F_mpz_mat_t mat1, 
							                                            const ulong entry1, const ulong x)
   \brief  Multiply entry1 of mat1 by the unsigned long x and subtract the result from res, in place.
*/
//void _F_mpz_entry_submul_ui(F_mpz_mat_t res, ulong entry2, const F_mpz_mat_t mat1, 
//							                                 const ulong entry1, const ulong x);

/** 
   \fn     _F_mpz_entry_addmul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong entry1, 
						                                       const F_mpz_mat_t mat2, const ulong entry2)
   \brief  Multiply entry1 of mat1 by entry2 of mat2 and add the result to res, in place.
*/
//void _F_mpz_entry_addmul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong entry1, 
//						                                 const F_mpz_mat_t mat2, const ulong entry2);

/** 
   \fn     _F_mpz_entry_submul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong entry1, 
						                                       const F_mpz_mat_t mat2, const ulong entry2)
   \brief  Multiply entry1 of mat1 by entry2 of mat2 and subtract the result from res, inplace.
*/
//void _F_mpz_entry_submul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong entry1, 
//						                                 const F_mpz_mat_t mat2, const ulong entry2);

/** 
   \fn     void F_mpz_mat_set_entry_si(F_mpz_mat_t mat, ulong n, const long x)
   \brief  Set entry n to the signed long value x. Coefficients are numbered
	        from the constant entry, starting at zero.
*/
//void F_mpz_mat_set_entry_si(F_mpz_mat_t mat, ulong n, const long x);

/** 
   \fn     void F_mpz_mat_get_entry_si(F_mpz_mat_t mat, ulong n, const long x)
   \brief  Return entry n of mat as a signed long. If n is greater than the degree
	        of mat, then zero is returned.
*/
//long F_mpz_mat_get_entry_si(const F_mpz_mat_t mat, const ulong n);

/** 
   \fn     void F_mpz_mat_set_entry_ui(F_mpz_mat_t mat, ulong n, const long x)
   \brief  Set entry n to the unsigned long value x. Coefficients are numbered
	        from the constant entry, starting at zero.
*/
//void F_mpz_mat_set_entry_ui(F_mpz_mat_t mat, ulong n, const ulong x);

/** 
   \fn     void F_mpz_mat_get_entry_ui(F_mpz_mat_t mat, ulong n, const long x)
   \brief  Return entry n of mat as an usigned long. If n is greater than the degree
	        of mat, then zero is returned.
*/
//ulong F_mpz_mat_get_entry_ui(const F_mpz_mat_t mat, const ulong n);

/** 
   \fn     void F_mpz_mat_set_entry_mpz(F_mpz_mat_t mat, ulong n, const mpz_t x)
   \brief  Set entry n to the mpz_t value x. Coefficients are numbered
	        from the constant entry, starting at zero.
*/
//void F_mpz_mat_set_entry_mpz(F_mpz_mat_t mat, ulong n, const mpz_t x);

/** 
   \fn     void F_mpz_mat_get_entry_mpz(mpz_t x, const F_mpz_mat_t mat, const ulong n)
   \brief  Return entry n of mat as an mpz_t. If n is greater than the degree
	        of mat, then zero is returned.
*/
//void F_mpz_mat_get_entry_mpz(mpz_t x, const F_mpz_mat_t mat, const ulong n);

/*===============================================================================

	Attributes

================================================================================*/

/*===============================================================================

	Truncation

================================================================================*/

/*===============================================================================

	Conversions

================================================================================*/

/** 
   \fn     void mpz_mat_to_F_mpz_mat(F_mpz_mat_t F_mat, const mpz_mat_t m_mat)
   \brief  Convert an mpz_mat_t to an F_mpz_mat_t
*/
void mpz_mat_to_F_mpz_mat(F_mpz_mat_t F_mat, const mpz_mat_t m_mat);

/** 
   \fn     void F_mpz_mat_to_mpz_mat(mpz_mat_t m_mat, const F_mpz_mat_t F_mat)
   \brief  Convert an F_mpz_mat_t to an mpz_mat_t
*/
void F_mpz_mat_to_mpz_mat(mpz_mat_t m_mat, const F_mpz_mat_t F_mat);


/*===============================================================================

	Assignment

================================================================================*/

/** 
   \fn     void F_mpz_mat_zero(F_mpz_mat_t mat)
   \brief  Sets mat to the zero matrix
*/
/*static inline 
void F_mpz_mat_zero(F_mpz_mat_t mat)
{
   mat->length = 0;
}*/

/** 
   \fn     void F_mpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2)
   \brief  Sets mat1 to equal mat2
*/
void F_mpz_mat_set(F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/*===============================================================================

	Comparison

================================================================================*/

/** 
   \fn     int F_mpz_mat_equal(const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);
   \brief  Returns 1 if mat1 and mat2 are equal (arithmetically), otherwise
	        returns 0.
*/
int F_mpz_mat_equal(const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/*===============================================================================

	Coefficient sizes

================================================================================*/

/** 
   \fn     long F_mpz_mat_max_bits(F_mpz_mat_t mat)
   \brief  Computes the largest number of bits n that any entry has and returns
	        -n if a negative entry exists in mat, else it returns n. Zero is 
			  returned for the zero matrix.
*/
//long F_mpz_mat_max_bits(const F_mpz_mat_t mat);

/** 
   \fn     ulong F_mpz_mat_max_limbs(F_mpz_mat_t mat)
   \brief  Returns the largest number of limbs required to store the absolute value
	        of entries of mat. Zero is returned for the zero matrix.
*/
//ulong F_mpz_mat_max_limbs(const F_mpz_mat_t mat);

/*===============================================================================

	Reverse

================================================================================*/

/*===============================================================================

	Negation

================================================================================*/

/** 
   \fn     void F_mpz_mat_neg(F_mpz_mat_t res, const F_mpz_mat_t mat)
   \brief  Set res to the negative of mat.
*/
void F_mpz_mat_neg(F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/*===============================================================================

	Addition/subtraction

================================================================================*/

/** 
   \fn     void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
   \brief  Sets res to the sum of mat1 and mat2.
*/
void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/** 
   \fn     void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
   \brief  Sets res to the difference of mat1 and mat2.
*/
void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/*===============================================================================

	Shifting

================================================================================*/

/** 
   \fn     void F_mpz_mat_left_shift(F_mpz_mat_t res, const F_mpz_mat_t mat, const ulong n)
   \brief  Multiplies mat by x^n and sets res to the result.
*/
//void F_mpz_mat_left_shift(F_mpz_mat_t res, const F_mpz_mat_t mat, const ulong n);

/** 
   \fn     void F_mpz_mat_right_shift(F_mpz_mat_t res, const F_mpz_mat_t mat, const ulong n)
   \brief  Divides mat by x^n, discarding any remainder, and sets res to the result.
*/
//void F_mpz_mat_right_shift(F_mpz_mat_t res, const F_mpz_mat_t mat, const ulong n);

/*===============================================================================

	Scalar multiplication

================================================================================*/

/** 
   \fn     void F_mpz_mat_row_mul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)

   \brief  Multiply mat2 by the unsigned long x and set mat1 to the result.
*/
void F_mpz_mat_row_mul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x);

/** 
   \fn     void F_mpz_mat_row_mul_si(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, long x)

   \brief  Multiply mat2 by the signed long x and set mat1 to the result.
*/
void F_mpz_mat_row_mul_si(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, long x);
/** 
   \fn     void F_mpz_mat_row_mul_mpz(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, mpz_t x)

   \brief  Multiply mat2 by the mpz_t x and set mat1 to the result.
*/
void F_mpz_mat_row_mul_mpz(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, mpz_t x);

/*===============================================================================

	Scalar addmul/submul

================================================================================*/

/** 
   \fn     void F_mpz_mat_row_addmul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)

   \brief  Multiply mat2 by the unsigned long x and add the result to mat1.
*/
void F_mpz_mat_row_addmul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x);

/** 
   \fn     void F_mpz_mat_row_submul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)

   \brief  Multiply mat2 by the unsigned long x and subtract the result from mat1.
*/
void F_mpz_mat_row_submul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x);

/** 
   \fn     void F_mpz_mat_row_addmul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp, F_mpz_mat_t temp)
   \brief  Multiply entry from row r2 of mat2 by c*2^exp and add to entry from row r1 of mat1.
	        This function takes a row "temp" which must be a 1 x c matrix where c is the column
			  width of mat1 and mat2.
*/

void F_mpz_mat_row_addmul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp, F_mpz_mat_t temp);

/** 
   \fn     void F_mpz_mat_row_submul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp, F_mpz_mat_t temp)
   \brief  Multiply entry from row r2 of mat2 by c*2^exp and subtract from entry from row r1 of mat1.
	        This function takes a row "temp" which must be a 1 x c matrix where c is the column
			  width of mat1 and mat2.
*/

void F_mpz_mat_row_submul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp, F_mpz_mat_t temp);


#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
