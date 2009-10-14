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
#include "F_mpz.h"

/*==============================================================================

   F_mpz_mat_t
   -----------

F_mpz_mat_t represents a dense matrix in M_{r,c}(Z) 

"entries" is an array of limbs, one for each entry

There are two things each entry in this array can represent:

1) If the most significant two bits are 01, then the entry represents
an index into the array "F_mpz_arr" from the F_mpz module, and the 
mpz_t in that array contains the entry.

2) Otherwise, the entry represents a signed entry
whose absolute value is no more than FLINT_BIT - 2 bits in length. The
entry is stored in twos complement format.

"r" is the number of rows of the matrix

"c" is the number of columns

"r_alloc" is the number of rows allocated

"c_alloc" is the number of columns allocated

"rows" is an array of pointers to the starts of rows of mat

================================================================================*/
 
typedef struct
{
   F_mpz * entries;
   ulong r;
	ulong c;
	ulong r_alloc; 
	ulong c_alloc;
	F_mpz ** rows;
} F_mpz_mat_struct;

// fmpz_mat_t allows reference-like semantics for F_mpz_mat_struct
typedef F_mpz_mat_struct F_mpz_mat_t[1];

/*===============================================================================

	Memory management

================================================================================*/

/** 
   \fn     void F_mpz_mat_init(F_mpz_mat_t mat, ulong r, ulong c)
   \brief  Initialise a matrix with r rows and c columns.
*/
void F_mpz_mat_init(F_mpz_mat_t mat, ulong r, ulong c);

/** 
   \fn     void F_mpz_mat_init_identity(F_mpz_mat_t mat, ulong n)
   \brief  Initialise a matrix to the n x n identity matrix.
*/
void F_mpz_mat_init_identity(F_mpz_mat_t mat, ulong n);

/** 
   \fn     void F_mpz_mat_clear(F_mpz_mat_t mat)
   \brief  Clear the matrix, releasing any memory it was using.
*/
void F_mpz_mat_clear(F_mpz_mat_t mat);

/** 
   \fn     void F_mpz_mat_resize(F_mpz_mat_t mat, const ulong r, const ulong c)
   \brief  Resize the matrix to have r rows and c columns. Existing data is preserved
	        where possible (i.e. for old rows and columns upto the new r and c) and
			  zeroes are added elsewhere (if c or r grows).
*/
void F_mpz_mat_resize(F_mpz_mat_t mat, const ulong r, const ulong c);

/*===============================================================================

	Input/Output

================================================================================*/

/** 
   \fn     int mpz_mat_from_string(mpz_mat_t mat, const char *s)
	\brief  Read an mpz_mat_t from a string at s
*/
int mpz_mat_from_string(mpz_mat_t mat, const char *s);

/** 
   \fn     char* mpz_mat_to_string(mpz_mat_t mat)
	\brief  Read a string from an mpz_mat_t
*/
char* mpz_mat_to_string(mpz_mat_t mat);

/** 
   \fn     int mpz_mat_from_string_pretty(mpz_mat_t mat, char *s)
	\brief  Read an mpz_mat_t from a pretty string at s.  A pretty string
                                             starts with [[
*/
int mpz_mat_from_string_pretty(mpz_mat_t mat, char *s);

/** 
   \fn     char* mpz_mat_to_string_pretty(mpz_mat_t mat)
	\brief  Read a pretty string from an mpz_mat_t
*/
char* mpz_mat_to_string_pretty(mpz_mat_t mat);

/** 
   \fn     void mpz_mat_fprint(mpz_mat_t mat, FILE* f)
	\brief  Print an mpz_mat_t to a file stream
*/
void mpz_mat_fprint(mpz_mat_t mat, FILE* f);

/** 
   \fn     void mpz_mat_fprint_pretty(mpz_mat_t mat, FILE* f)
	\brief  Print a pretty format mpz_mat_t to a file stream 
*/
void mpz_mat_fprint_pretty(mpz_mat_t mat, FILE* f);

/** 
   \fn     int mpz_mat_fread(mpz_mat_t mat, FILE* f)
	\brief  Read an mpz_mat_t from a file stream
*/
int mpz_mat_fread(mpz_mat_t mat, FILE* f);

/** 
   \fn     int mpz_mat_fread_pretty(mpz_mat_t mat, FILE* f)
	\brief  Read a an mpz_mat_t from a file stream with pretty formatting
*/
int mpz_mat_fread_pretty(mpz_mat_t mat, FILE* f);

/** 
   \fn     int F_mpz_mat_from_string(F_mpz_mat_t mat, const char *s)
	\brief  Read an F_mpz_mat_t from a string
*/
int F_mpz_mat_from_string(F_mpz_mat_t mat, const char *s);

/** 
   \fn     char* F_mpz_mat_to_string(F_mpz_mat_t mat)
	\brief  Write an F_mpz_mat_t to a string
*/
char* F_mpz_mat_to_string(F_mpz_mat_t mat);

/** 
   \fn     int F_mpz_mat_from_string_pretty(F_mpz_mat_t mat, char *s)
	\brief  Read an F_mpz_mat_t from a string in pretty format
*/
int F_mpz_mat_from_string_pretty(F_mpz_mat_t mat, char *s);

/** 
   \fn     char* F_mpz_mat_to_string_pretty(F_mpz_mat_t mat)
	\brief  Write an F_mpz_mat_t to a string in pretty format
*/
char* F_mpz_mat_to_string_pretty(F_mpz_mat_t mat);

/** 
   \fn     void F_mpz_mat_print(F_mpz_mat_t mat)
   \brief  Print an F_mpz_mat to stdout.
*/
void F_mpz_mat_print(F_mpz_mat_t mat);

/** 
   \fn     void F_mpz_mat_print_pretty(F_mpz_mat_t mat)
	\brief  Prints an F_mpz_mat_t to the screen in pretty format
*/
void F_mpz_mat_print_pretty(F_mpz_mat_t mat);

/** 
   \fn     void F_mpz_mat_fprint(F_mpz_mat_t mat, FILE* f)
	\brief  Prints an F_mpz_mat_t to a file stream
*/
void F_mpz_mat_fprint(F_mpz_mat_t mat, FILE* f);

/** 
   \fn     void F_mpz_mat_fprint_pretty(F_mpz_mat_t mat, FILE* f)
	\brief  Prints an F_mpz_mat_t to a file stream in pretty format
*/
void F_mpz_mat_fprint_pretty(F_mpz_mat_t mat, FILE* f);

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

/** 
   \fn     long F_mpz_mat_set_line_d(double * appv, const F_mpz_mat_t mat, const ulong r, const int n)
   \brief  Sets the entries of appv to the double mantissa of the entries of the given row of mat
	        with respect to a single maximum exponent which is returned.
*/
long F_mpz_mat_set_line_d(double * appv, const F_mpz_mat_t mat, const ulong r, const int n);

/*===============================================================================

	Assignment

================================================================================*/

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
   \fn     void F_mpz_mat_row_add(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, 
	                                 const ulong r1, const F_mpz_mat_t mat2, const ulong r1, 
												                       const ulong start, const ulong n)
   \brief  Sets the given row of res to the sum of the given rows of mat1 and mat2.
*/
void F_mpz_mat_row_add(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, const ulong r1, 
							  const F_mpz_mat_t mat2, const ulong r2, const ulong start, const ulong n);

/** 
   \fn     void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
   \brief  Sets res to the difference of mat1 and mat2.
*/
void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/** 
   \fn     void F_mpz_mat_row_sub(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, 
	                                 const ulong r1, const F_mpz_mat_t mat2, const ulong r1, 
												                       const ulong start, const ulong n)
   \brief  Sets the given row of res to the difference of the given rows of mat1 and mat2.
*/
void F_mpz_mat_row_sub(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, const ulong r1, 
							  const F_mpz_mat_t mat2, const ulong r2, const ulong start, const ulong n);

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
   \fn     void F_mpz_mat_row_mul_F_mpz(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x)

   \brief  Multiply row r2 of mat2 by the mpz_t x and set row r1 of mat1 to the result.
*/
void F_mpz_mat_row_mul_F_mpz(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x);

/*===============================================================================

	Scalar addmul/submul

================================================================================*/

/** 
   \fn     void F_mpz_mat_row_addmul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)

   \brief  Multiply row r2 of mat2 by the unsigned long x and add the result to row r1 of mat1.
*/
void F_mpz_mat_row_addmul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x);

/** 
   \fn     void F_mpz_mat_row_addmul(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x)

   \brief  Multiply row r2 of mat2 by the F_mpz_t x and add the result to row r1 of mat1.
*/
void F_mpz_mat_row_addmul(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x);

/** 
   \fn     void F_mpz_mat_row_submul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)

   \brief  Multiply row r2 of mat2 by the unsigned long x and subtract the result from row r1 of mat1.
*/
void F_mpz_mat_row_submul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x);

/** 
   \fn     void F_mpz_mat_row_submul(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x)

   \brief  Multiply row r2 of mat2 by the F_mpz_t x and subtract the result from row r1 of mat1.
*/
void F_mpz_mat_row_submul(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x);

/** 
   \fn     void F_mpz_mat_row_addmul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp)

	\brief  Multiply entry from row r2 of mat2 by c*2^exp and add to entry from row r1 of mat1.
	        
*/

void F_mpz_mat_row_addmul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp);

/** 
   \fn     void F_mpz_mat_row_submul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp)

	\brief  Multiply entry from row r2 of mat2 by c*2^exp and subtract from entry from row r1 of mat1.
	        
*/

void F_mpz_mat_row_submul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp);


/** 
   \fn     void F_mpz_mat_row_swap(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, 
	                                                ulong r2, ulong start, ulong n)

	\brief  Swap row r1 of mat1 with row r2 of mat2. If mat1 aliases mat2 this is 
	        done efficiently by swapping row pointers.	        
*/
void F_mpz_mat_row_swap(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, 
								                          ulong r2, ulong start, ulong n);

/** 
   \fn     void F_mpz_mat_row_neg(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, 
	                                                ulong r2, ulong start, ulong n)

	\brief  Set row r1 of mat1 to minus row r2 of mat2.
*/
void F_mpz_mat_row_neg(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, 
								                          ulong r2, ulong start, ulong n);


#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
