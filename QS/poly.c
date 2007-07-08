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
   unsigned long min, fact_approx, fact, span;
   
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
   mpz_div_ui(temp, temp, SIEVE_SIZE/2);
   mpz_root(temp, temp, s);
   fact_approx = mpz_get_ui(temp);
   
   for (fact = 0; fact_approx >= factor_base[fact].p; fact++); 
   
   span = num_primes/s/s/2;
   min = fact - span/2;
   while ((fact*fact)/min - min < span) min--;

   poly_inf->min = min;
   poly_inf->fact = fact;
   poly_inf->span = span;
   
   mpz_clear(temp); 
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
