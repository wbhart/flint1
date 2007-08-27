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
#include "mp_lprels.h"

/*=========================================================================
   linear_algebra_init:
 
   Function: Allocate space for the various linear algebra structures
 
==========================================================================*/

void linear_algebra_init(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf)
{
   la_col_t * matrix;
   mpz_t * Y_arr;
   unsigned long prec = qs_inf->prec+1;
   unsigned long small_primes = qs_inf->small_primes;
   
   const unsigned long buffer_size = 2*(qs_inf->num_primes + EXTRA_RELS + 200); // Allows for 1/2 of relations to be duplicates
   
   la_inf->small = (unsigned long *) flint_stack_alloc(small_primes);
   la_inf->factor = (fac_t *) flint_stack_alloc_bytes(sizeof(fac_t)*MAX_FACS);
   
   matrix = la_inf->matrix = (la_col_t *) flint_stack_alloc_bytes(sizeof(la_col_t)*(qs_inf->num_primes + EXTRA_RELS + 400));
   la_inf->unmerged = la_inf->matrix + qs_inf->num_primes + EXTRA_RELS + 200;
   Y_arr = la_inf->Y_arr = (mpz_t *) flint_stack_alloc_bytes(sizeof(mpz_t)*buffer_size);
   la_inf->curr_rel = la_inf->relation = (unsigned long *) flint_stack_alloc(buffer_size*MAX_FACS*2);
   la_inf->qsort_arr = (la_col_t **) flint_stack_alloc(200);
   la_inf->rel_str = (char *) flint_stack_alloc(MPQS_STRING_LENGTH);
   
   la_inf->lpnew = flint_fopen("lpnew", "w");
   FILE * lprels = flint_fopen("lprels","w");
   fclose(lprels);
    
   for (unsigned long i = 0; i < buffer_size; i++) 
   {
      mpz_init2(Y_arr[i], prec);
   }
   for (unsigned long i = 0; i < qs_inf->num_primes + EXTRA_RELS + 400; i++) 
   {
      matrix[i].weight = 0;
   }
   
   la_inf->num_unmerged = 0;
   la_inf->num_lp_unmerged = 0;
   la_inf->columns = 0;
   la_inf->num_relations = 0;
}
   
void linear_algebra_clear(linalg_t * la_inf, QS_t * qs_inf)
{
   la_col_t * matrix = la_inf->matrix;
   la_col_t * unmerged = la_inf->unmerged;
   mpz_t * Y_arr = la_inf->Y_arr;
   const unsigned long buffer_size = 4*(qs_inf->num_primes + EXTRA_RELS + 200)/2;
   
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
   
   fclose(la_inf->lpnew);
   
   flint_stack_release(); // Clear rel_str
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

unsigned long merge_lp_relations(QS_t * qs_inf, poly_t * poly_inf, linalg_t * la_inf)
{
   FILE * comb;
   unsigned long combined;
   
   fclose(la_inf->lpnew); 
   sort_lp_file("lpnew");
   comb = flint_fopen("comb","w");
   mergesort_lp_file("lprels", "lpnew", "tmp", comb);
   fclose(comb);
   la_inf->lpnew = flint_fopen("lpnew","w");
   
   mpz_t factor;
   mpz_init(factor);
   comb = flint_fopen("comb", "r");
   combined = combine_large_primes(qs_inf, la_inf, poly_inf, comb, factor);
   mpz_clear(factor);
   fclose(comb);
     
   la_inf->num_lp_unmerged = 0;
   
   return combined;         
}

/*==========================================================================
   Insert large prime partial relation:

   Function: Insert the partial relation into the lprels file, return the 
             number of full relations obtained after any sorting and merging
   
===========================================================================*/

unsigned long insert_lp_relation(QS_t * qs_inf, linalg_t * la_inf, poly_t * poly_inf, mpz_t Y, mpz_t res)
{
   char * rel_str = la_inf->rel_str;
   char * rel_ptr = rel_str;
   char Q_str[200];
   char Y_str[200];
   FILE * LPNEW = la_inf->lpnew;

   unsigned long small_primes = qs_inf->small_primes;
   
   unsigned long * small = la_inf->small;
   const unsigned long num_factors = la_inf->num_factors; 
   fac_t * factor = la_inf->factor;
      
   unsigned long fac_num = 0; 
   
   for (unsigned long i = 0; i < small_primes; i++)
   {
       if (small[i]) add_factor(&rel_ptr, (unsigned long) small[i], (unsigned long) i);
   }
   for (unsigned long i = 0; i < num_factors; i++)
   {
       add_factor(&rel_ptr, (unsigned long) factor[i].exp, (unsigned long) factor[i].ind);
   }
   add_0(&rel_ptr);
   
   gmp_sprintf(Y_str, "%Zd\0", Y);
   gmp_sprintf(Q_str, "%Zd\0", res);
   
   fprintf(LPNEW, "%s @ %s :%s\n", Q_str, Y_str, rel_str);
                          
   la_inf->num_lp_unmerged++;
   if ((la_inf->num_lp_unmerged %256) == 0) printf("%ld partials\n", la_inf->num_lp_unmerged);
   
   if (la_inf->num_lp_unmerged == 500)
   {
      return merge_lp_relations(qs_inf, poly_inf, la_inf);
   }
      
   return 0;
}

/*==========================================================================
   Insert relation:

   Function: Insert the relation into the matrix and store the Y value
   
===========================================================================*/

unsigned long insert_relation(QS_t * qs_inf, linalg_t * la_inf, poly_t * poly_inf, mpz_t Y)
{
   la_col_t * unmerged = la_inf->unmerged;
   unsigned long num_unmerged = la_inf->num_unmerged;
   unsigned long small_primes = qs_inf->small_primes;
   
   unsigned long * small = la_inf->small;
   const unsigned long num_factors = la_inf->num_factors; 
   fac_t * factor = la_inf->factor;
   
   unsigned long * curr_rel = la_inf->curr_rel;
   
   unsigned long fac_num = 0; 
   clear_col(unmerged + num_unmerged);
   
   for (unsigned long i = 0; i < small_primes; i++)
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
   
#if TEST3
   mpz_t X, temp, temp2;
   mpz_init(X);
   mpz_init(temp);
   mpz_init(temp2);
   mpz_set_ui(X, 1);
   
   for (unsigned long j = 0; j < curr_rel[0]; j++)
   {
       mpz_set_ui(temp, qs_inf->factor_base[curr_rel[2*j + 1]].p);
       mpz_pow_ui(temp, temp, curr_rel[2*j + 2]); 
       mpz_mul(X, X, temp);
   }
   
   mpz_mod(X, X, qs_inf->mpz_n);
   mpz_mul(temp, Y, Y);
   mpz_mod(temp, temp, qs_inf->mpz_n);
   if (mpz_cmp(X, temp) != 0)
   {
      mpz_add(temp2, temp, X);
      if (mpz_cmp(temp2, qs_inf->mpz_n) != 0) 
      {
         gmp_printf("X = %Zd (mod N) != \nY^2 = %Zd (mod N)\n\n", X, temp);
         gmp_printf("n = %Zd\n", qs_inf->mpz_n);
      }
   }
   mpz_clear(X);
   mpz_clear(temp);
   mpz_clear(temp2);
#endif
   
   la_inf->curr_rel += MAX_FACS*2;
   la_inf->num_unmerged++;
   la_inf->num_relations++;
   
   if (la_inf->num_unmerged == 100)
   {
      return merge_relations(la_inf);
   }
   
   return 0;
}

