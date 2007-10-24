/******************************************************************************

 poly.c
 
 Routines for managing polynomials

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../flint.h"
#include "../memory-manager.h"
#include "../long_extras.h"
#include "../longlong_wrapper.h"
#include "../longlong.h"

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
   unsigned long s = (qs_inf->bits-1)/28+1;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long fact_approx, fact, span;
   long min; 
   
   poly_inf->s = s;
     
   poly_inf->B_terms = (unsigned long*) flint_stack_alloc(s);  
   
   poly_inf->A_ind = (unsigned long*) flint_stack_alloc(s);  
   poly_inf->A_modp = (unsigned long*) flint_stack_alloc(s);  
   poly_inf->A_inv2B = (unsigned long**) flint_stack_alloc(s); 
   poly_inf->inv_p2 = (double*) flint_stack_alloc_bytes(s*sizeof(double));  
   poly_inf->A_inv = (unsigned long*) flint_stack_alloc(num_primes);  
   poly_inf->soln1 = (unsigned long*) flint_stack_alloc(num_primes); 
   poly_inf->soln2 = (unsigned long*) flint_stack_alloc(num_primes); 
   
   unsigned long ** A_inv2B = poly_inf->A_inv2B;
   
   A_inv2B[0] = (unsigned long *) flint_stack_alloc(num_primes*s);
   
   mpz_init(poly_inf->C);
   
   for (unsigned long i = 1; i < s; i++)
   {
      A_inv2B[i] = A_inv2B[i-1] + num_primes;
   } 
    
   mpz_t temp;
   mpz_init(temp); 
   
   mpz_mul_ui(temp, N, 2*qs_inf->k);
   mpz_sqrt(temp, temp);
   
   mpz_div_ui(temp, temp, 300);
   poly_inf->target_A = mpz_get_ui(temp);
   
   mpz_root(temp, temp, s);
   fact_approx = mpz_get_ui(temp);
   
   for (fact = 0; fact_approx >= factor_base[fact].p; fact++); 
   
   span = num_primes/s/s/2;
   if (span < 4*s) span = 4*s;
   min = fact - span/2;
   if (min < SMALL_PRIMES) min = SMALL_PRIMES;
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

void poly_clear(poly_t * poly_inf)
{
   mpz_clear(poly_inf->C);
   flint_stack_release(); // release all A_inv2B[i]
   flint_stack_release(); // release soln1
   flint_stack_release(); // release soln2
   flint_stack_release(); // release A_inv
   flint_stack_release(); // release inv_p2
   flint_stack_release(); // release A_inv2B
   flint_stack_release(); // release A_modp
   flint_stack_release(); // release A_ind
   flint_stack_release(); // release B_terms
   
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
   unsigned long factor, i, p;
   unsigned long diff, best_diff, best1, best2;
   
   unsigned long A;
   
   if (s <= 4) 
   {
       A_ind[0] = z_randint(span) + min;
       do
       {
          A_ind[1] = z_randint(span) + min;
       } while (A_ind[0] == A_ind[1]);
   }
   
   if (s == 2) A = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p;
   
   if ((s == 3) || (s == 4))
   {
       do
       {
          A_ind[2] = z_randint(span) + min;
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
       A_ind[0] = ((z_randint(span) + min) | 1);
       if (A_ind[0] == min + span) A_ind[0] -= 2;
       
       do
       {
          A_ind[1] = ((z_randint(span) + min) | 1);
          if (A_ind[1] == min + span) A_ind[1] -= 2;
       } while (A_ind[0] == A_ind[1]);
       
       do
       {
          A_ind[2] = ((z_randint(span) + min) | 1);
          if (A_ind[2] == min + span) A_ind[2] -= 2;
       } while ((A_ind[0] == A_ind[2]) || (A_ind[1] == A_ind[2]));
       
       A = factor_base[A_ind[0]].p * factor_base[A_ind[1]].p * factor_base[A_ind[2]].p;
       factor = poly_inf->target_A / A;
       
       for (i = 0; i < 8; i++)
       {
          A_ind[3] = ((z_randint(span) + min) & -2L);
          if (A_ind[3] < min) A_ind[3]+=2;
          
          do
          {
             A_ind[4] = ((z_randint(span) + min) & -2L);
             if (A_ind[4] < min) A_ind[4]+=2;
          } while (A_ind[3] == A_ind[4]);
          
          if (i == 0)
          {
             best_diff = FLINT_ABS(factor_base[A_ind[3]].p * factor_base[A_ind[4]].p - factor);
             best1 = A_ind[3];
             best2 = A_ind[4];
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
 
   for (i = 0; i < s; i++)
   {
      p = factor_base[A_ind[i]].p;
      poly_inf->inv_p2[i] = z_precompute_inverse(p*p);
   }      
}

/*=========================================================================
   compute B terms:
 
   Function: Compute the terms from which the B values of the polynomials 
             are constructed and compute the starting B coefficient
 
==========================================================================*/

void compute_B_terms(QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long s = poly_inf->s;
   unsigned long * A_ind = poly_inf->A_ind;
   unsigned long * A_modp = poly_inf->A_modp;
   unsigned long * B_terms = poly_inf->B_terms;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long A = poly_inf->A;
   unsigned long B;
   unsigned long p, temp, temp2, i;
   double pinv;
   
   for (i = 0; i < s; i++)
   {
      p = factor_base[A_ind[i]].p;
      pinv = factor_base[A_ind[i]].pinv;
      temp2 = (temp = z_div_64_precomp(A, p, pinv)); 
      A_modp[i] = (temp = z_mod_64_precomp(temp, p, pinv));
      temp = z_invert(temp, p);
      temp = z_mulmod_precomp(temp, qs_inf->sqrts[A_ind[i]], p, pinv);
      if (temp > p/2) temp = p - temp;
      B_terms[i] = temp*temp2;     
   }
   
   B = B_terms[0];
   for (i = 1; i < s; i++)
   {
      B += B_terms[i];
   }
   poly_inf->B = B;
}

/*=========================================================================
   Compute offsets and hypercube polynomial correction factors:
 
   Function: Compute the starting offsets in the sieve for each prime
             and the polynomial correction factors used by the 
             hypercube method
 
==========================================================================*/

void compute_off_adj(QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long num_primes = qs_inf->num_primes;
   unsigned long A = poly_inf->A;
   unsigned long B = poly_inf->B;
   unsigned long * A_inv = poly_inf->A_inv;
   unsigned long ** A_inv2B = poly_inf->A_inv2B;
   unsigned long * B_terms = poly_inf->B_terms;
   unsigned long * soln1 = poly_inf->soln1;
   unsigned long * soln2 = poly_inf->soln2;
   u_int32_t * sqrts = qs_inf->sqrts;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long s = poly_inf->s;
   unsigned long p, temp;
   double pinv;
   
   for (unsigned long i = 2; i < num_primes; i++) // skip k and 2
   {
      p = factor_base[i].p;
      pinv = factor_base[i].pinv;
      
      A_inv[i] = z_invert(z_mod_64_precomp(A, p, pinv), p);
             
      for (unsigned long j = 0; j < s; j++)
      {
         temp = z_mod_64_precomp(B_terms[j], p, pinv);
         temp = z_mulmod_precomp(temp, A_inv[i], p, pinv);
         temp *= 2;
         if (temp >= p) temp -= p;
         A_inv2B[j][i] = temp;
      }
             
      temp = z_mod_64_precomp(B, p, pinv);
      temp = sqrts[i] + p - temp;
      temp *= A_inv[i];
      temp += SIEVE_SIZE/2;
      soln1[i] = z_mod_64_precomp(temp, p, pinv); // Consider using z_mod_precomp
      temp = p - sqrts[i];
      if (temp == p) temp -= p;
      temp = z_mulmod_precomp(temp, A_inv[i], p, pinv);
      temp *= 2;
      if (temp >= p) temp -= p;      
      soln2[i] = temp+soln1[i];
      if (soln2[i] >= p) soln2[i] -= p;
   }  
}

/*=========================================================================
   Compute offsets and hypercube polynomial correction factors:
 
   Function: Compute the starting offsets in the sieve for each prime
             and the polynomial correction factors used by the 
             hypercube method
 
==========================================================================*/
void compute_A_factor_offsets(QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long s = poly_inf->s;
   unsigned long * A_ind = poly_inf->A_ind;
   unsigned long * A_modp = poly_inf->A_modp;
   unsigned long * soln1 = poly_inf->soln1;
   unsigned long * soln2 = poly_inf->soln2;
   unsigned long p, D;
   unsigned long * n = qs_inf->n;
   unsigned long B = poly_inf->B;
   unsigned long temp, temp2, B_modp2, index, p2; 
   prime_t * factor_base = qs_inf->factor_base;
   double * inv_p2 = poly_inf->inv_p2;
   double pinv;
   
   for (unsigned long j = 0; j < s; j++)
   {
      index = A_ind[j];
      p = factor_base[index].p;
      p2 = p*p;
      pinv = factor_base[index].pinv;
      D = z_ll_mod_precomp(n[2], n[1], p*p, inv_p2[j]);    
      if ((long) B < 0) 
      {
         B_modp2 = z_mod_64_precomp(-B, p2, inv_p2[j]);
         B_modp2 = p2 - B_modp2;
         if (B_modp2 == p2) B_modp2 = 0;
      } else
      B_modp2 = z_mod_64_precomp(B, p2, inv_p2[j]);
      temp = B_modp2*A_modp[j];
      temp = z_mod_64_precomp(temp, p, pinv); 
      temp2 = z_invert(temp, p);
      D -= (B_modp2*B_modp2);
      if ((long) D < 0) temp = -z_div_64_precomp(-D, p, pinv);
      else temp = -z_div_64_precomp(-D, p, pinv);
      temp *= temp2;
      temp += SIEVE_SIZE/2;
      if ((long) temp < 0) 
      {
         temp = p - z_mod_64_precomp(-temp, p, pinv);
         if (temp == p) temp = 0;
      }
      else temp = z_mod_64_precomp(temp, p, pinv);
      soln1[index] = temp;
      soln2[index] = -1L;
   }
}          

/*=========================================================================
   Compute C:
 
   Function: Compute the C coefficient of the polynomial with the 
             current A and B values
 
==========================================================================*/

void compute_C(QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long A = poly_inf->A;
   unsigned long B = poly_inf->B;
   mpz_t * C = &poly_inf->C;
   mpz_t * mpz_n = &qs_inf->mpz_n;
   
   if ((long) B < 0L) B = -B;
   mpz_set_ui(*C, B);
   mpz_mul_ui(*C, *C, B);
   mpz_sub(*C, *C, *mpz_n);
   mpz_divexact_ui(*C, *C, A);
} 
