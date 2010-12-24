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
#include <mpfr.h>

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
   \fn     int F_mpz_mat_from_string(F_mpz_mat_t mat, const char * s)
	\brief  Read an F_mpz_mat_t from a string
*/
int F_mpz_mat_from_string(F_mpz_mat_t mat, const char * s);

/** 
   \fn     char * F_mpz_mat_to_string(F_mpz_mat_t mat)
	\brief  Write an F_mpz_mat_t to a string
*/
char * F_mpz_mat_to_string(F_mpz_mat_t mat);

/** 
   \fn     int F_mpz_mat_from_string_pretty(F_mpz_mat_t mat, char * s)
	\brief  Read an F_mpz_mat_t from a string in pretty format
*/
int F_mpz_mat_from_string_pretty(F_mpz_mat_t mat, char * s);

/** 
   \fn     char * F_mpz_mat_to_string_pretty(F_mpz_mat_t mat)
	\brief  Write an F_mpz_mat_t to a string in pretty format
*/
char * F_mpz_mat_to_string_pretty(F_mpz_mat_t mat);

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
void F_mpz_mat_fprint(F_mpz_mat_t mat, FILE * f);

/** 
   \fn     void F_mpz_mat_fprint_pretty(F_mpz_mat_t mat, FILE* f)
	\brief  Prints an F_mpz_mat_t to a file stream in pretty format
*/
void F_mpz_mat_fprint_pretty(F_mpz_mat_t mat, FILE * f);

/** 
   \fn     int F_mpz_mat_fread(F_mpz_mat_t mat, FILE* f)
	\brief  Read an F_mpz_mat_t from a file stream
*/
int F_mpz_mat_fread(F_mpz_mat_t mat, FILE * f);

/** 
   \fn     int F_mpz_mat_fread_pretty(F_mpz_mat_t mat, FILE* f)
	\brief  Read a F_mpz_mat_t from a file stream in pretty format
                                       useful with fpLLL's generate function
*/
int F_mpz_mat_fread_pretty(F_mpz_mat_t mat, FILE * f);

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
   \fn     long _F_mpz_vec_to_d_vec_2exp(double * appv, const F_mpz * vec, const ulong n)
   \brief  Sets the entries of appv to the double mantissa of the entries of the given vec
	        with respect to a single maximum exponent which is returned.
*/
long _F_mpz_vec_to_d_vec_2exp(double * appv, const F_mpz * vec, const ulong n);

/** 
   \fn     void _F_mpz_vec_to_mpfr_vec(mpfr_t * appv, const F_mpz * vec, const ulong n)
   \brief  Sets the entries of appv to the entries of the given vec.
*/
void _F_mpz_vec_to_mpfr_vec(__mpfr_struct * appv, const F_mpz * vec, const ulong n);

/*===============================================================================

	Assignment

================================================================================*/

/** 
   \fn     void F_mpz_mat_set(fmpz_mat_t mat1, const fmpz_mat_t mat2)
   \brief  Sets mat1 to equal mat2
*/
void F_mpz_mat_set(F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/** 
   \fn     void F_mpz_mat_swap(fmpz_mat_t mat1, const fmpz_mat_t mat2)
   \brief  Efficiently swap mat1 with mat2
*/
static inline
void F_mpz_mat_swap(F_mpz_mat_t mat1, F_mpz_mat_t mat2)
{
   F_mpz * pt = mat1->entries;
   mat1->entries = mat2->entries;
   mat2->entries = pt;
   F_mpz ** ppt = mat1->rows;
   mat1->rows = mat2->rows;
   mat2->rows = ppt;
   ulong t = mat1->r_alloc;
   mat1->r_alloc = mat2->r_alloc;
   mat2->r_alloc = t;
   t = mat1->c_alloc;
   mat1->c_alloc = mat2->c_alloc;
   mat2->c_alloc = t;
   t = mat1->r;
   mat1->r = mat2->r;
   mat2->r = t;
   t = mat1->c;
   mat1->c = mat2->c;
   mat2->c = t;
}


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

	Vector memory management

================================================================================*/

/** 
   \fn     F_mpz * _F_mpz_vec_init(ulong n)
   \brief  Allocate and initialise a vector of length n.
*/
static inline
F_mpz * _F_mpz_vec_init(ulong n)
{
   ulong i;
   F_mpz * vec = (F_mpz *) flint_heap_alloc(n);
   for (i = 0; i < n; i++)
      F_mpz_init(vec + i);
   return vec;
}

/** 
   \fn     void _F_mpz_vec_clear(F_mpz * vec, ulong n)
   \brief  Clear the given vector which is assumed to have been allocated
           at length n.
*/
static inline
void _F_mpz_vec_clear(F_mpz * vec, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
	  F_mpz_clear(vec + i);
   flint_heap_free(vec);
}

/*===============================================================================

	Vector copy and comparison

================================================================================*/

/** 
   \fn     void _F_mpz_vec_copy(F_mpz * vec1, const F_mpz * vec2, ulong n)
   \brief  Copy the contents of vec2 of length n into vec1.
*/
static inline
void _F_mpz_vec_copy(F_mpz * vec1, const F_mpz * vec2, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
      F_mpz_set(vec1 + i, vec2 + i);
}

/** 
   \fn     int _F_mpz_vec_equal(F_mpz * vec1, const F_mpz * vec2, ulong n)
   \brief  Return 1 if vec1 and vec2 of length n are equal, otherwise return 0.
*/
static inline
int _F_mpz_vec_equal(F_mpz * vec1, const F_mpz * vec2, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
      if (F_mpz_cmp(vec1 + i, vec2 + i) != 0) 
		 return 0;
   return 1;
}

/*===============================================================================

	Addition/subtraction

================================================================================*/

/** 
   \fn     void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
   \brief  Sets res to the sum of mat1 and mat2.
*/
void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/** 
   \fn     _F_mpz_vec_add(F_mpz * res, const F_mpz * vec1, const F_mpz * vec2, ulong n)
   \brief  Sets the res to the sum of the given vectors.
*/
void _F_mpz_vec_add(F_mpz * res, const F_mpz * vec1, const F_mpz * vec2, ulong n);

/** 
   \fn     void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
   \brief  Sets res to the difference of mat1 and mat2.
*/
void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2);

/** 
   \fn     _F_mpz_vec_sub(F_mpz * res, const F_mpz * vec1, const F_mpz * vec2, ulong n)
   \brief  Sets tres to the difference of the given vectors.
*/
void _F_mpz_vec_sub(F_mpz * res, const F_mpz * vec1, const F_mpz * vec2, ulong n);

/*===============================================================================

	Scalar multiplication

================================================================================*/

/** 
   \fn     void F_mpz_mat_row_mul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)

   \brief  Multiply vec2 by the unsigned long x and set vec1 to the result..
*/
void _F_mpz_vec_mul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x);

/** 
   \fn     void F_mpz_mat_row_mul_si(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, long x)

   \brief  Multiply vec2 by the signed long x and set vec1 to the result.
*/
void _F_mpz_vec_mul_si(F_mpz * vec1, F_mpz* vec2, ulong n, long x);
/** 
   \fn     _F_mpz_vec_mul_F_mpz(F_mpz * vec1, F_mpz* vec2, ulong n, F_mpz_t x)

   \brief  Multiply vec2 by the mpz_t x and set vec1 to the result.
*/
void _F_mpz_vec_mul_F_mpz(F_mpz * vec1, F_mpz* vec2, ulong n, F_mpz_t x);

/*===============================================================================

	Scalar addmul/submul

================================================================================*/

/** 
   \fn     _F_mpz_vec_addmul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x)

   \brief  Multiply vec2 by the unsigned long x and add the result to vec1.
*/
void _F_mpz_vec_addmul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x);

/** 
   \fn     _F_mpz_vec_addmul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x)

   \brief  Multiply vec2 by the F_mpz_t x and add the result to vec1.
*/
void _F_mpz_vec_addmul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x);

/** 
   \fn     _F_mpz_vec_submul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x)

   \brief  Multiply vec2 by the unsigned long x and subtract the result from vec1.
*/
void _F_mpz_vec_submul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x);

/** 
   \fn     _F_mpz_vec_submul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x)

   \brief  Multiply vec2 by the F_mpz_t x and subtract the result from vec1.
*/
void _F_mpz_vec_submul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x);

/** 
   \fn     _F_mpz_vec_addmul_2exp_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong c, ulong exp)

	\brief  Multiply entry from row r2 of mat2 by c*2^exp and add to entry from row r1 of mat1.	        
*/

void _F_mpz_vec_addmul_2exp_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong c, ulong exp);

/** 
   \fn     _F_mpz_vec_submul_2exp_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong c, ulong exp)

	\brief  Multiply entries from vec2 by c*2^exp and add to entries in vec1.
	        
*/

void _F_mpz_vec_submul_2exp_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong c, ulong exp);

/** 
   \fn     _F_mpz_vec_submul_2exp_F_mpz(F_mpz* vec1, F_mpz * vec2, ulong n, F_mpz_t c, ulong exp)

	\brief  Multiply entries from vec2 by c*2^exp and subtract from entries in vec1.
	        
*/
void _F_mpz_vec_submul_2exp_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t c, ulong exp);

/** 
   \fn     void F_mpz_mat_swap_rows(F_mpz_mat_t mat, ulong r1, ulong r2)

	\brief  Swap rows r1 and r2 of mat. This is done efficiently by swapping 
	        row pointers.	        
*/
static inline
void F_mpz_mat_swap_rows(F_mpz_mat_t mat, ulong r1, ulong r2)
{
	if (r1 == r2) return; // rows are the same, nothing to do

    F_mpz * temp = mat->rows[r1]; // swap rows via temporary
    mat->rows[r1] = mat->rows[r2];
    mat->rows[r2] = temp;
}

/** 
   \fn     void _F_mpz_vec_neg(F_mpz * vec1, F_mpz * vec2, ulong n)

	\brief  Set vec1 of length n to minus vec2.
*/
static inline
void _F_mpz_vec_neg(F_mpz * vec1, F_mpz * vec2, ulong n)
{
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_neg(vec1 + i, vec2 + i);
}

/* ======================================================================================================

   Classical Multiplication for F_mpz_mat.h

=========================================================================================================*/

/** 
   \fn     void _F_mpz_mat_mul_classical(F_mpz_mat_t res, const F_mpz_mat_t mat1,
                                                  const F_mpz_mat_t mat2)

	\brief  Classical multiplication of F_mpz_mat_t's mat1 and mat2 set result to res
                                                  not alias safe	        
*/
void _F_mpz_mat_mul_classical(F_mpz_mat_t res, const F_mpz_mat_t mat1,
                                                  const F_mpz_mat_t mat2);

/** 
   \fn     static inline
           void F_mpz_mat_mul_classical(F_mpz_mat_t P, const F_mpz_mat_t A,
                                                  const F_mpz_mat_t B)

	\brief  Classical multiplication of F_mpz_mat_t's mat1 and mat2 set result to res	        
                                                  alias safe
*/
static inline
void F_mpz_mat_mul_classical(F_mpz_mat_t P, const F_mpz_mat_t A, const F_mpz_mat_t B)
{

	if ((P == A) || (P == B))
	{
		F_mpz_mat_t Pa;
		F_mpz_mat_init(Pa, P->r, P->c);
        _F_mpz_mat_mul_classical(Pa, A, B);
		F_mpz_mat_swap(P, Pa);
		F_mpz_mat_clear(Pa);
		return;
	} else
	   return _F_mpz_mat_mul_classical(P, A, B);
}

/*===========================================================

   Properties

===========================================================*/

/*
   Determine the maximum number of bits b required to store the absolute
   value of the largest coefficient of M. If M has signed coefficients, 
   return -b, else return b. Return the row and column that such an 
   entry was found. WARNING: this need not return the maximum coefficient
   only one of the ones which requires the maximum number of bits.
*/
long F_mpz_mat_max_bits2(ulong * row, ulong * col, const F_mpz_mat_t M);

/*
   Determine the maximum number of bits b required to store the absolute
   value of the largest coefficient of M. If M has signed coefficients, 
   return -b, else return b.
*/
static inline
long F_mpz_mat_max_bits(const F_mpz_mat_t M)
{
   return F_mpz_mat_max_bits2(NULL, NULL, M);
}

/*===========================================================

   Matrix by Scalar Operations

===========================================================*/

/*
   Divides M by the scalar 2^n and stores at res.  Truncates each entry with F_mpz_div_2exp.
*/
void F_mpz_mat_div_2exp(F_mpz_mat_t res, F_mpz_mat_t M, ulong n);

/*
   Scalar multiplication by a power of 2.  res = M * 2^n.
*/
void F_mpz_mat_mul_2exp(F_mpz_mat_t res, F_mpz_mat_t M, ulong n);

/*===========================================================

   Scalar Product

===========================================================*/

/** 
   \fn     void _F_mpz_vec_scalar_product(F_mpz_t sp, F_mpz * vec1, F_mpz * vec2, ulong n)

	\brief  Set sp to the scalar product of vec1 and vec2.
*/
void _F_mpz_vec_scalar_product(F_mpz_t sp, F_mpz * vec1, F_mpz * vec2, ulong n);

/*===========================================================

   Support functions for Z[x] factoring

===========================================================*/

/*
   Compares column a and column b of M returns 1 if they are the same 
   and 0 otherwise.
*/
int F_mpz_mat_col_equal(F_mpz_mat_t M, ulong a, ulong b);

/*
   Set column a of M to be equal to column b
*/
void F_mpz_mat_col_copy(F_mpz_mat_t M, ulong a, ulong b);

/*
   The array _part_ must have room for at least M->c ulongs
   which must initially be set to 0.
   The matrix is then partitioned into equivalence classes of equal 
   columns. The equivalence classes are numbered 1..s say. Upon
   termination, part contains the number of the equivalence class
   for each column. If it is not possible for the RREF of
   the matrix to be a basis, zero is returned.
*/
int F_mpz_mat_col_partition(ulong * part, F_mpz_mat_t M);

/*
   U is set to a readonly window starting at row r0 and col c0 of M
   and consisting of the given number of rows and cols.
*/
void F_mpz_mat_window_init(F_mpz_mat_t U, F_mpz_mat_t M, 
                      ulong r0, ulong c0, ulong rows, ulong cols);

static inline
void F_mpz_mat_window_clear(F_mpz_mat_t U)
{
   if (U->r)
      free(U->rows);
}

/*
   Uses F_mpz_smod to reduce each entry in M mod P taking the representative 
   in (-P/2, P/2]
*/
void F_mpz_mat_smod(F_mpz_mat_t res, F_mpz_mat_t M, F_mpz_t P);

/*
   This is the heart of the Novocin PhD algorithm.  Takes a column of data, col, 
   a large modulus P, and an upperbound for the data in a good vector with exp.  
   Decides if it's worth calling LLL with this data, if not it returns 0, if so 
   it adjusts M in place by truncating the data, augmenting M, and returning the 
   virtual weight of the final column.  May or may not add a row for P, based on 
   a rough estimation.
*/
int _F_mpz_mat_next_col(F_mpz_mat_t M, F_mpz_t P, F_mpz_mat_t col, long exp, long U_exp);


void F_mpz_mat_randintrel(F_mpz_mat_t mat, ulong bits);

void F_mpz_mat_randintrel_little_big(F_mpz_mat_t mat, ulong n, ulong bits1, ulong bits);

void F_mpz_mat_rand_unimodular(F_mpz_mat_t U, ulong bits);

void F_mpz_mat_rand_unimodular_little_big(F_mpz_mat_t U, ulong n, ulong bits1, ulong bits2);

void F_mpz_mat_transpose(F_mpz_mat_t res, F_mpz_mat_t M);

void F_mpz_mat_block_reverse_cols(F_mpz_mat_t res, ulong n, F_mpz_mat_t M);

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
