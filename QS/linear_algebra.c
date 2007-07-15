/******************************************************************************

 linear_algebra.c
 
 Routines for dealing with building and handling the final F_2 matrix

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../flint.h"
#include "../memory-manager.h"

#include "common.h"
#include "poly.h"
#include "linear_algebra.h"
#include "block_lanczos.h"

/*=========================================================================
   linear_algebra_init:
 
   Function: Allocate space for the various linear algebra structures
 
==========================================================================*/

void linear_algebra_init(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf)
{
   la_col_t * matrix;
   mpz_t * Y_arr;
   
   la_inf->small = (unsigned long *) flint_stack_alloc(SMALL_PRIMES);
   la_inf->factor = (fac_t *) flint_stack_alloc_bytes(sizeof(fac_t)*MAX_FACS);
   matrix = la_inf->matrix = (la_col_t *) flint_stack_alloc_bytes(sizeof(la_col_t)*(qs_inf->num_primes + EXTRA_RELS + 100));
   la_inf->unmerged = la_inf->matrix + qs_inf->num_primes + EXTRA_RELS;
   Y_arr = la_inf->Y_arr = (mpz_t *) flint_stack_alloc_bytes(sizeof(mpz_t)*(qs_inf->num_primes + EXTRA_RELS + 100));
   la_inf->unmerged_Y = Y_arr + qs_inf->num_primes + EXTRA_RELS;
   la_inf->curr_rel = la_inf->relation = (unsigned long *) flint_stack_alloc((qs_inf->num_primes + EXTRA_RELS + 100)*MAX_FACS*2);
   la_inf->curr_unmerged = la_inf->unmerged_rels = la_inf->curr_rel + (qs_inf->num_primes + EXTRA_RELS)*MAX_FACS*2;
   for (unsigned long i = 0; i < qs_inf->num_primes + EXTRA_RELS + 100; i++) 
   {
      clear_col(matrix + i);
      mpz_init2(Y_arr[i], 128);
   }
   
   la_inf->num_unmerged = 0;
   la_inf->columns = 0;
}
   
void linear_algebra_clear(linalg_t * la_inf, QS_t * qs_inf)
{
   la_col_t * matrix = la_inf->matrix;
   la_col_t * unmerged = la_inf->unmerged;
   mpz_t * Y_arr = la_inf->Y_arr;
   
   for (unsigned long i = 0; i < qs_inf->num_primes + EXTRA_RELS + 100; i++) 
   {
      mpz_clear(Y_arr[i]);
   }
   
   for (unsigned long i = 0; i < la_inf->columns; i++) // Clear all used columns
   {
      free_col(matrix + i);
   }
   
   for (unsigned long i = 0; i < la_inf->num_unmerged; i++) // Clear all used columns
   {
      free_col(unmerged + i);
   }
   
   flint_stack_release(); // Clear relation
   flint_stack_release(); // Clear Y_arr
   flint_stack_release(); // Clear matrix and unmerged
   flint_stack_release(); // Clear factor
   flint_stack_release(); // Clear small
}
   
/*==========================================================================
   Merge relations:

   Function: Merge unmerged relations into the matrix
   
===========================================================================*/

unsigned long merge_relations(linalg_t * la_inf)
{
   const unsigned long num_unmerged = la_inf->num_unmerged;
   unsigned long columns = la_inf->columns;
   la_col_t * matrix = la_inf->matrix;
   la_col_t * unmerged = la_inf->unmerged;
   mpz_t * unmerged_Y = la_inf->unmerged_Y;
   mpz_t * Y_arr = la_inf->Y_arr;
   unsigned long * rel_point1 = la_inf->unmerged_rels;
   unsigned long * rel_point2 = la_inf->curr_rel;
    
   for (unsigned long i = 0; i < num_unmerged; i++)
   {
      copy_col(matrix + columns, unmerged + i);
      matrix[columns].orig = columns;
      clear_col(unmerged + i);
      columns++;
      mpz_set(Y_arr[columns], unmerged_Y[i]); 
      for (unsigned long j = 0; j < 2*MAX_FACS; j++)
      {
         rel_point2[j] = rel_point1[j];
      }
      rel_point1 += 2*MAX_FACS;
      rel_point2 += 2*MAX_FACS;
   }
   
   la_inf->curr_rel = rel_point2;
   la_inf->columns = columns;
   la_inf->num_unmerged = 0;
   la_inf->curr_unmerged = la_inf->unmerged_rels;
   
   return num_unmerged;
}

/*==========================================================================
   Insert relation:

   Function: Insert the relation into the matrix and store the Y value
   
===========================================================================*/

unsigned long insert_relation(linalg_t * la_inf, poly_t * poly_inf, mpz_t Y)
{
   la_col_t * unmerged = la_inf->unmerged;
   unsigned long * small = la_inf->small;
   const unsigned long num_factors = la_inf->num_factors; 
   fac_t * factor = la_inf->factor;
   unsigned long num_unmerged = la_inf->num_unmerged;
   unsigned long * curr_unmerged = la_inf->curr_unmerged;
   unsigned long fac_num = 0;
   
   for (unsigned long i = 0; i < SMALL_PRIMES; i++)
   {
       if (small[i] & 1) insert_col_entry(unmerged + num_unmerged, i);
       if (small[i]) 
       {
          curr_unmerged[2*fac_num + 1] = i;
          curr_unmerged[2*fac_num + 2] = small[i];
          fac_num++;
       }
   }
   for (unsigned long i = 0; i < num_factors; i++)
   {
       if (factor[i].exp & 1) insert_col_entry(unmerged + num_unmerged, factor[i].ind);
       curr_unmerged[2*fac_num + 1] = factor[i].ind;
       curr_unmerged[2*fac_num + 2] = factor[i].exp;
       fac_num++;
   }
   curr_unmerged[0] = fac_num;
   la_inf->curr_unmerged += MAX_FACS*2;
   mpz_set(la_inf->unmerged_Y[num_unmerged], Y); 
   la_inf->num_unmerged++;
   
   if (la_inf->num_unmerged == 100)
   {
      return merge_relations(la_inf);
   }
   
   return 0;
}

#include "common.h"

