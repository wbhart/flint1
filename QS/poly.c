/******************************************************************************

 poly.c
 
 Routines for managing polynomials

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>

#include "../flint.h"
#include "../memory-manager.h"
#include "../long_extras.h"

#include "poly.h"
#include "common.h"

/*=========================================================================
   poly_init:
 
   Function: computes parameters for the polynomials and initialises the 
             various structures required
 
==========================================================================*/

void poly_init(QS_t * qs_inf, poly_t * poly_inf, mpz_t N)
{
   unsigned long num_primes = qs_inf->num_primes;
   unsigned long s = qs_inf->bits/28+1;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long fact_approx, fact, span;
   long min; 
   
   poly_inf->s = s;
     
   poly_inf->A_ind = (unsigned long*) flint_stack_alloc(s);  
   poly_inf->A_modp = (unsigned long*) flint_stack_alloc(s);  
   poly_inf->A_inv2B = (unsigned long**) flint_stack_alloc(s); 
   poly_inf->A_inv = (unsigned long*) flint_stack_alloc(num_primes);  
   poly_inf->soln1 = (unsigned long*) flint_stack_alloc(num_primes); 
   poly_inf->soln2 = (unsigned long*) flint_stack_alloc(num_primes); 
   
   unsigned long ** A_inv2B = poly_inf->A_inv2B;
   
   A_inv2B[0] = (unsigned long *) flint_stack_alloc(num_primes*s);
   
   for (unsigned long i = 1; i < s; i++)
   {
      A_inv2B[i] = A_inv2B[i-1] + num_primes;
   } 
    
   mpz_t temp;
   mpz_init(temp); 
   
   mpz_mul_ui(temp, N, 2*qs_inf->k);
   mpz_sqrt(temp, temp);
   
   mpz_div_ui(temp, temp, 1000);
   poly_inf->target_A = mpz_get_ui(temp);
   
   mpz_root(temp, temp, s);
   fact_approx = mpz_get_ui(temp);
   
   for (fact = 0; fact_approx >= factor_base[fact].p; fact++); 
   
   span = num_primes/s/s/2;
   if (span < 10) span = 10;
   min = fact - span/2;
   if (min < 2) min = 2;
   if (min + span >= qs_inf->num_primes) span = num_primes - min - 1;
   fact = min + span/2;

#if POLY_PARAMS   
   printf("min = FB[%ld], span = %ld, number of factors = %ld\n", min, span, s);
#endif
   
   poly_inf->min = min;
   poly_inf->fact = fact;
   poly_inf->span = span;
   
   mpz_clear(temp); 
}

/*=========================================================================
   compute_A:
 
   Function: Compute a new polynomial A value
             The function attempts to pick A near to an optimal size
 
==========================================================================*/

void compute_A(QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long min = poly_inf->min;
   unsigned long span = poly_inf->span;
   unsigned long s = poly_inf->s;
   unsigned long * A_ind = poly_inf->A_ind;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long factor, i;
   unsigned long diff, best_diff, best1, best2;
   
   unsigned long A;
   
   if (s <= 4) 
   {
       A_ind[0] = long_randint(span) + min;
       do
       {
          A_ind[1] = long_randint(span) + min;
       } while (A_ind[0] == A_ind[1]);
   }
   
   if (s == 2) A = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
   
   if ((s == 3) || (s == 4))
   {
       do
       {
          A_ind[2] = long_randint(span) + min;
       } while ((A_ind[0] == A_ind[2]) || (A_ind[1] == A_ind[2]));
       A = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p;
   }  
   
   if (s == 4)
   {
      factor = (poly_inf->target_A - 1) / A + 1; 
      for (i = min; i < min+span; i++)
      {
         if ((factor_base[i].p > factor) && (i != A_ind[0]) && (i != A_ind[1]) && (i != A_ind[2])) break;
      } 
      if (i == min + span)
      {
         i--;
         while ((i == A_ind[0]) || (i == A_ind[1]) || (i == A_ind[2])) i--;
      }   
      A_ind[3] = i;
      A *= factor_base[A_ind[3]].p;
   }
   
   if (s == 5) 
   {
       A_ind[0] = ((long_randint(span) + min) | 1);
       if (A_ind[0] == min + span) A_ind[0] -= 2;
       
       do
       {
          A_ind[1] = ((long_randint(span) + min) | 1);
          if (A_ind[1] == min + span) A_ind[1] -= 2;
       } while (A_ind[0] == A_ind[1]);
       
       do
       {
          A_ind[2] = ((long_randint(span) + min) | 1);
          if (A_ind[2] == min + span) A_ind[2] -= 2;
       } while ((A_ind[0] == A_ind[2]) || (A_ind[1] == A_ind[2]));
       
       A = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p;
       factor = poly_inf->target_A / A;
       
       for (i = 0; i < 8; i++)
       {
          A_ind[3] = ((long_randint(span) + min) & -2L);
          if (A_ind[3] < min) A_ind[3]+=2;
          
          do
          {
             A_ind[4] = ((long_randint(span) + min) & -2L);
             if (A_ind[4] < min) A_ind[4]+=2;
          } while (A_ind[3] == A_ind[4]);
          
          if (i == 0)
          {
             best_diff = FLINT_ABS(factor_base[A_ind[3]].p * factor_base[A_ind[4]].p - factor);
             best1 = factor_base[A_ind[3]].p;
             best2 = factor_base[A_ind[4]].p;
             continue;
          }
          
          diff = FLINT_ABS(factor_base[A_ind[3]].p * factor_base[A_ind[4]].p - factor);
          
          if (diff < best_diff)
          {
             best_diff = diff;
             best1 = A_ind[3];
             best2 = A_ind[4];
          }
       }
       
       A_ind[3] = best1;
       A_ind[4] = best2;
       A = A * factor_base[A_ind[3]].p * factor_base[A_ind[4]].p;
   }  
   
   poly_inf->A = A;

#if POLY_A
   if ((s == 4) || (s == 5)) printf("A = %ld, target A = %ld\n", A, poly_inf->target_A);
#endif          
}

void poly_clear(void)
{
   flint_stack_release(); // release all A_inv2B[i]
   flint_stack_release(); // release soln1
   flint_stack_release(); // release soln2
   flint_stack_release(); // release A_inv
   flint_stack_release(); // release A_inv2B
   flint_stack_release(); // release A_modp
   flint_stack_release(); // release A_ind
   
}
