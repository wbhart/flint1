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

===============================================================================*/
/******************************************************************************

 linear_algebra.h
 Header file for linear_algebra.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef LINALG_H
#define LINALG_H

#include <gmp.h>

#include "common.h"
#include "poly.h"
#include "block_lanczos.h"


#define DUPS 0 // Print info about number of duplicate relations

typedef struct fac_s
{
   unsigned long ind;
   unsigned long exp;
} fac_t;

typedef struct linalg_s
{
   unsigned long * small; // Exponents of small primes in currently evaluated candidate
   fac_t * factor; // An array of factors with exponents corresponding to the current candidate
   unsigned long num_factors; //The length of the factor array for the current candidate
   
   la_col_t * matrix;  // The final sorted F_2 matrix plus possibly some empty columns at the start
   unsigned long columns; // The number of columns in the matrix so far
   
   la_col_t * unmerged; // A new list of unmerged F_2 columns
   unsigned long num_unmerged; // The current number of unmerged F_2 relations
   
   mpz_t * Y_arr; // The Y values corresponding to all relations found
      
   unsigned long * relation; // The list of all relations found 
   unsigned long * curr_rel; // Pointer to where we have got up to in the list of relations found
   unsigned long num_relations; // Total number of relations found so far

   la_col_t ** qsort_arr; // An array of pointers to the unmerged relations for quicksort
} linalg_t;

void linear_algebra_init(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf);
   
void linear_algebra_clear(linalg_t * la_inf, QS_t * qs_inf);

unsigned long merge_sort(linalg_t * la_inf);
      
unsigned long merge_relations(linalg_t * la_inf);

unsigned long insert_relation(linalg_t * la_inf, poly_t * poly_inf, mpz_t Y);

#endif
