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


#define MAX_FACS 25 // Maximum number of different prime factors
                    // a relation can have

typedef struct fac_s
{
   unsigned long ind;
   unsigned long exp;
} fac_t;

typedef struct linalg_s
{
   unsigned long * small;
   fac_t * factor;
   unsigned long num_factors;
   la_col_t * matrix;
   la_col_t * unmerged;
   unsigned long columns;
   unsigned long num_unmerged;
   mpz_t * Y_arr;
   mpz_t * unmerged_Y;
   unsigned long * relation;
   unsigned long * curr_rel;
   unsigned long * unmerged_rels;
   unsigned long * curr_unmerged;
} linalg_t;

void linear_algebra_init(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf);
   
void linear_algebra_clear(linalg_t * la_inf, QS_t * qs_inf);

unsigned long merge_relations(linalg_t * la_inf);

unsigned long insert_relation(linalg_t * la_inf, poly_t * poly_inf, mpz_t Y);

#endif
