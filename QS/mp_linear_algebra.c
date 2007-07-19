/******************************************************************************

 mp_linear_algebra.c
 
 Routines for dealing with building and handling the final F_2 matrix

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../flint.h"
#include "../memory-manager.h"

#include "common.h"
#include "mp_poly.h"
#include "mp_linear_algebra.h"
#include "block_lanczos.h"

/*=========================================================================
   linear_algebra_init:
 
   Function: Allocate space for the various linear algebra structures
 
==========================================================================*/

void linear_algebra_init(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf)
{
   la_col_t * matrix;
   mpz_t * Y_arr;
   unsigned long prec = qs_inf->prec+1;
   
   const unsigned long buffer_size = 3*(qs_inf->num_primes + EXTRA_RELS + 100)/2; // Allows for 1/3 of relations to be duplicates
   
   la_inf->small = (unsigned long *) flint_stack_alloc(SMALL_PRIMES);
   la_inf->factor = (fac_t *) flint_stack_alloc_bytes(sizeof(fac_t)*MAX_FACS);
   
   matrix = la_inf->matrix = (la_col_t *) flint_stack_alloc_bytes(sizeof(la_col_t)*(qs_inf->num_primes + EXTRA_RELS + 200));
   la_inf->unmerged = la_inf->matrix + qs_inf->num_primes + EXTRA_RELS + 100;
   Y_arr = la_inf->Y_arr = (mpz_t *) flint_stack_alloc_bytes(sizeof(mpz_t)*buffer_size);
   la_inf->curr_rel = la_inf->relation = (unsigned long *) flint_stack_alloc(buffer_size*MAX_FACS*2);
   la_inf->qsort_arr = (la_col_t **) flint_stack_alloc(100);
   
   for (unsigned long i = 0; i < buffer_size; i++) 
   {
      mpz_init2(Y_arr[i], prec);
   }
   for (unsigned long i = 0; i < qs_inf->num_primes + EXTRA_RELS + 200; i++) 
   {
      matrix[i].weight = 0;
   }
   
   la_inf->num_unmerged = 0;
   la_inf->columns = 0;
   la_inf->num_relations = 0;
}
   
void linear_algebra_clear(linalg_t * la_inf, QS_t * qs_inf)
{
   la_col_t * matrix = la_inf->matrix;
   la_col_t * unmerged = la_inf->unmerged;
   mpz_t * Y_arr = la_inf->Y_arr;
   const unsigned long buffer_size = 3*(qs_inf->num_primes + EXTRA_RELS + 100)/2;
   
   for (unsigned long i = 0; i < buffer_size; i++) 
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
   
   flint_stack_release(); // Clear qsort_array
   flint_stack_release(); // Clear relation
   flint_stack_release(); // Clear Y_arr
   flint_stack_release(); // Clear matrix and unmerged
   flint_stack_release(); // Clear factor
   flint_stack_release(); // Clear small
}

/*=========================================================================
   
   Compare relations:
 
   Function: Compare two relations; used by qsort
 
==========================================================================*/
int relations_cmp(const void *a, const void *b)
{
  la_col_t * ra = *((la_col_t **) a);
  la_col_t * rb = *((la_col_t **) b);
  long point;
  if (ra->weight > rb->weight) return 1;
  else if (ra->weight < rb->weight) return -1;
  
  for (point = ra->weight-1; (ra->data[point] == rb->data[point]) && (point >= 0); point--)
  {
      ;
  }
  
  if (point == -1L) return 0;
  
  if (ra->data[point] > rb->data[point]) return 1;
  else if (ra->data[point] < rb->data[point]) return -1;
}

int relations_cmp2(const void *a, const void *b)
{
  la_col_t * ra = (la_col_t *) a;
  la_col_t * rb = (la_col_t *) b;
  if (!ra->weight) printf("ra error\n");
  long point;
  if (ra->weight > rb->weight) return 1;
  else if (ra->weight < rb->weight) return -1;
  
  for (point = ra->weight-1; (ra->data[point] == rb->data[point]) && (point >= 0); point--)
  {
      ;
  }
  
  if (point == -1L) return 0;

  if (ra->data[point] > rb->data[point]) return 1;
  else if (ra->data[point] < rb->data[point]) return -1;
}
  
/*==========================================================================
   Merge sort:

   Function: Merge a list of sorted new relations into a list of existing
             sorted relations. Sort is done using a merge sort algorithm
             with a short stack.
   
===========================================================================*/

unsigned long merge_sort(linalg_t * la_inf)
{
   la_col_t * matrix = la_inf->matrix;
   long columns = la_inf->columns;
   
   la_col_t ** qsort_arr = la_inf->qsort_arr;
   long num_unmerged = la_inf->num_unmerged;
   
   long dups = 0;
   int comp;
   
   for (long i = columns + num_unmerged - 1L; i >= dups; i--) 
   {
      if (!columns) comp = -1;
      else if (!num_unmerged) comp = 1;
      else 
      {
         comp = relations_cmp2(matrix + columns - 1L, qsort_arr[num_unmerged - 1L]);
      }
      switch (comp)
      {
         case -1: 
         {
            copy_col(matrix + i, qsort_arr[num_unmerged - 1L]);
            clear_col(qsort_arr[num_unmerged - 1L]);
            num_unmerged--;
            break;
         }
         case 1 : 
         {
            copy_col(matrix + i, matrix + columns - 1L);
            columns--;
            break;
         }
         case 0 : 
         {
            free_col(qsort_arr[num_unmerged - 1L]);
            clear_col(qsort_arr[num_unmerged - 1L]);
            num_unmerged--;
            copy_col(matrix + i, matrix + columns - 1L);
            columns--;
            dups++;
            break;
         }
      } 
   }
   
   columns = la_inf->columns + la_inf->num_unmerged - dups;
   
   if (dups)
   {
      for (unsigned long i = 0; i < columns; i++)
      {
         copy_col(matrix + i, matrix + i + dups);
      }
   }
            
   la_inf->columns = columns;
   columns = la_inf->num_unmerged - dups;
   la_inf->num_unmerged = 0;

#if DUPS
   printf("%ld new, %ld dups\n", columns, dups);
#endif   

   return columns;
}

/*==========================================================================
   Merge relations:

   Function: Merge unmerged relations into the matrix
   
===========================================================================*/

unsigned long merge_relations(linalg_t * la_inf)
{
   const unsigned long num_unmerged = la_inf->num_unmerged;
   la_col_t * unmerged = la_inf->unmerged;
   la_col_t ** qsort_arr = la_inf->qsort_arr;
   
   if (num_unmerged)
   {
      for (unsigned long i = 0; i < num_unmerged; i++)
      {
         qsort_arr[i] = unmerged + i;
      }
      qsort(qsort_arr, num_unmerged, sizeof(la_col_t *), relations_cmp);
      
      if ((la_inf->num_relations & 7) == 0) printf("%ld relations found\n", la_inf->num_relations);

      return merge_sort(la_inf);
   }
   
   return 0;
}

/*==========================================================================
   Insert relation:

   Function: Insert the relation into the matrix and store the Y value
   
===========================================================================*/

unsigned long insert_relation(linalg_t * la_inf, poly_t * poly_inf, mpz_t Y)
{
   la_col_t * unmerged = la_inf->unmerged;
   unsigned long num_unmerged = la_inf->num_unmerged;
   
   unsigned long * small = la_inf->small;
   const unsigned long num_factors = la_inf->num_factors; 
   fac_t * factor = la_inf->factor;
   
   unsigned long * curr_rel = la_inf->curr_rel;
   
   unsigned long fac_num = 0; 
   clear_col(unmerged + num_unmerged);
   
   for (unsigned long i = 0; i < SMALL_PRIMES; i++)
   {
       if (small[i] & 1) insert_col_entry(unmerged + num_unmerged, i);
       if (small[i]) 
       {
          curr_rel[2*fac_num + 1] = i;
          curr_rel[2*fac_num + 2] = small[i];
          fac_num++;
       }
   }
   for (unsigned long i = 0; i < num_factors; i++)
   {
       if (factor[i].exp & 1) insert_col_entry(unmerged + num_unmerged, factor[i].ind);
       curr_rel[2*fac_num + 1] = factor[i].ind;
       curr_rel[2*fac_num + 2] = factor[i].exp;
       fac_num++;
   }
   curr_rel[0] = fac_num;
   
   unmerged[num_unmerged].orig = la_inf->num_relations;
   
   mpz_set(la_inf->Y_arr[la_inf->num_relations], Y); 
   
   la_inf->curr_rel += MAX_FACS*2;
   la_inf->num_unmerged++;
   la_inf->num_relations++;
   
   if (la_inf->num_unmerged == 100)
   {
      return merge_relations(la_inf);
   }
   
   return 0;
}

