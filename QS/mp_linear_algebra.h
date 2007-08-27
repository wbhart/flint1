/******************************************************************************

 mp_linear_algebra.h
 Header file for mp_linear_algebra.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MPLINALG_H
#define MPLINALG_H

#include <gmp.h>
#include <stdio.h>

#include "common.h"
#include "mp_poly.h"
#include "block_lanczos.h"

#define DUPS 0 // Print info about number of duplicate relations

#define TEST3 0 // Checks if X = Y^2 mod N immediately after storing each relation

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
   unsigned long num_lp_unmerged; // The current number of unmerged partial relations
   
   mpz_t * Y_arr; // The Y values corresponding to all relations found
      
   unsigned long * relation; // The list of all relations found 
   unsigned long * curr_rel; // Pointer to where we have got up to in the list of relations found
   unsigned long num_relations; // Total number of relations found so far

   la_col_t ** qsort_arr; // An array of pointers to the unmerged relations for quicksort
   
   char * rel_str;
   FILE * lpnew;
} linalg_t;

void linear_algebra_init(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf);
   
void linear_algebra_clear(linalg_t * la_inf, QS_t * qs_inf);

unsigned long merge_sort(linalg_t * la_inf);
      
unsigned long merge_lp_relations(QS_t * qs_inf, poly_t * poly_inf, linalg_t * la_inf);

unsigned long merge_relations(linalg_t * la_inf);

unsigned long insert_lp_relation(QS_t * qs_inf, linalg_t * la_inf, poly_t * poly_inf, mpz_t Y, mpz_t res);

unsigned long insert_relation(QS_t * qs_inf, linalg_t * la_inf, poly_t * poly_inf, mpz_t Y);

#endif
