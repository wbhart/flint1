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

 mp_poly.c
 
 Routines for managing polynomials

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "../flint.h"
#include "../mpn_extras.h"
#include "../fmpz.h"
#include "../memory-manager.h"
#include "../long_extras.h"
#include "../longlong_wrapper.h"
#include "../longlong.h"

#include "mp_poly.h"
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
   if (s <= 2) s = 3;
	prime_t * factor_base = qs_inf->factor_base;
   unsigned long fact_approx, fact, span;
   unsigned long sieve_size = qs_inf->sieve_size;
   unsigned long small_primes = qs_inf->small_primes;
   long min; 
   
   poly_inf->s = s;
     
   poly_inf->B = (unsigned long*) flint_stack_alloc(qs_inf->prec+1);  
   poly_inf->B_terms = (unsigned long*) flint_stack_alloc(s*(qs_inf->prec+1));  
   poly_inf->A = (unsigned long*) flint_stack_alloc(qs_inf->prec+1);  
   poly_inf->target_A = (unsigned long*) flint_stack_alloc(qs_inf->prec+1);  
   
   poly_inf->A_ind = (unsigned long*) flint_stack_alloc(s);  
   poly_inf->A_modp = (unsigned long*) flint_stack_alloc(s);  
   poly_inf->A_inv2B = (uint32_t**) flint_stack_alloc(s); 
   poly_inf->inv_p2 = (double*) flint_stack_alloc_bytes(s*sizeof(double));  
   poly_inf->A_inv = (uint32_t *) flint_stack_alloc_bytes(num_primes*sizeof(uint32_t));  
   poly_inf->soln1 = (uint32_t *) flint_stack_alloc_bytes(num_primes*sizeof(uint32_t)); 
   poly_inf->soln2 = (uint32_t *) flint_stack_alloc_bytes(num_primes*sizeof(uint32_t)); 
   poly_inf->posn1 = (uint32_t *) flint_stack_alloc_bytes(num_primes*sizeof(uint32_t)); 
   poly_inf->posn2 = (uint32_t *) flint_stack_alloc_bytes(num_primes*sizeof(uint32_t)); 
   
   uint32_t ** A_inv2B = poly_inf->A_inv2B;
   
   A_inv2B[0] = (uint32_t *) flint_stack_alloc_bytes(num_primes*s*sizeof(uint32_t));
   
   mpz_init(poly_inf->A_mpz);
   mpz_init(poly_inf->B_mpz);
   mpz_init(poly_inf->C);
   
   unsigned long i;
   for (i = 1; i < s; i++)
   {
      A_inv2B[i] = A_inv2B[i-1] + num_primes;
   } 
 
   mpz_t temp;
   mpz_init(temp); 
   
   mpz_mul_ui(temp, N, 2*qs_inf->k);
   mpz_sqrt(temp, temp);
   
   mpz_div_ui(temp, temp, 47*sieve_size/100);
   mpz_to_fmpz(poly_inf->target_A, temp);
   
   mpz_root(temp, temp, s);
   fact_approx = mpz_get_ui(temp);
   
   for (fact = 0; fact_approx >= factor_base[fact].p; fact++); 
   
   span = num_primes/s/s/2;
   if (span < 6*s) span = 6*s;
   min = fact - span/2;
   if (min < (long) small_primes) min = small_primes;
   if (min + span >= qs_inf->num_primes - 1) span = num_primes - min - 1;
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
   mpz_clear(poly_inf->A_mpz);
   mpz_clear(poly_inf->B_mpz);
   mpz_clear(poly_inf->C);
   flint_stack_release(); // release all A_inv2B[i]
   flint_stack_release(); // release posn1
   flint_stack_release(); // release posn2
   flint_stack_release(); // release soln1
   flint_stack_release(); // release soln2
   flint_stack_release(); // release A_inv
   flint_stack_release(); // release inv_p2
   flint_stack_release(); // release A_inv2B
   flint_stack_release(); // release A_modp
   flint_stack_release(); // release A_ind
   flint_stack_release(); // release target_A
   flint_stack_release(); // release A
   flint_stack_release(); // release B_terms
   flint_stack_release(); // release B
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
   unsigned long * A = poly_inf->A;
   unsigned long * target_A = poly_inf->target_A;
   unsigned long * current_A = (unsigned long *) flint_stack_alloc(qs_inf->prec+1);  
   unsigned long * diff = (unsigned long *) flint_stack_alloc(qs_inf->prec+1);  
   unsigned long * best_diff = (unsigned long *) flint_stack_alloc(qs_inf->prec+1);  
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long factor, p;
   unsigned long best1, best2, best3;
   unsigned long odds = s - 3;
   mp_limb_t msl;
   int taken;
   long i, j, k;
   
   A[0] = 1;
   A[1] = 1;
   for (i = 0; i < odds; i++) // Randomly choose the first s-3 prime factors of A with odd indices
   {
      do
      {
         taken = 0;
         A_ind[i] = ((z_randint(span) + min) | 1);
         if (A_ind[i] == min + span) A_ind[i] -= 2;
         for (j = 0; j < i; j++)
         {
            if (A_ind[i] == A_ind[j]) taken = 1;
         }  
      } while (taken);
      msl = mpn_mul_1(A+1, A+1, A[0], factor_base[A_ind[i]].p);
      if (msl) // Compute the product of these s-3 primes
      {
         A[A[0]+1] = msl;
         A[0]++;
      }
   }
 
   for (k = 0; k < 30; k++) // Now try 8 different sets of even index primes as the remaining factors
   {
      F_mpn_copy(current_A, A, A[0] + 1);
      for (i = 0; i < 3; i++) // Randomly choose the last 3 prime factors of A with even indices
      {
         do
         {
            taken = 0;
            A_ind[s-3+i] = ((z_randint(span) + min) & -2L);
            if (A_ind[s-3+i] < min) A_ind[s-3+i] += 2;
            for (j = 0; j < i; j++)
            {
               if (A_ind[s-3+i] == A_ind[s-3+j]) taken = 1;
            }  
         } while (taken);

         msl = mpn_mul_1(current_A+1, current_A+1, current_A[0], factor_base[A_ind[s-3+i]].p);
         if (msl) // Compute the product of these s-3 primes and the odd indexed primes
         {
            current_A[current_A[0]+1] = msl;
            current_A[0]++;
         }
      }
      
      if (k == 0)  // Just store the first difference as the best one
      {
         if (target_A[0] >= current_A[0]) // Compute the difference with the target A
         {
            msl = mpn_sub(best_diff+1, target_A+1, target_A[0], current_A+1, current_A[0]);
            best_diff[0] = target_A[0];
         }
         else 
         {
            msl = mpn_sub(best_diff+1, current_A+1, current_A[0], target_A+1, target_A[0]);
            best_diff[0] = current_A[0];
         }
         if (msl) F_mpn_negate(best_diff+1, best_diff+1, best_diff[0]);
         while ((!best_diff[best_diff[0]]) && (best_diff[0])) best_diff[0]--; // Normalise best_diff
         
         best1 = A_ind[s-3];
         best2 = A_ind[s-2];
         best3 = A_ind[s-1];
   
         continue;
      }

      if (target_A[0] >= current_A[0]) // Compute the difference with the target A
      {
         msl = mpn_sub(diff+1, target_A+1, target_A[0], current_A+1, current_A[0]);
         diff[0] = target_A[0];
      }
      else 
      {
         msl = mpn_sub(diff+1, current_A+1, current_A[0], target_A+1, target_A[0]);
         diff[0] = current_A[0];
      }
      if (msl) F_mpn_negate(diff+1, diff+1, diff[0]);
      while ((!diff[diff[0]]) && (diff[0])) diff[0]--; // Normalise diff

      if ((diff[0] < best_diff[0]) || ((diff[0] == best_diff[0]) && (mpn_cmp(diff+1, best_diff+1, diff[0]) < 0)))  // The new diff is better
      {
         F_mpn_copy(best_diff, diff, diff[0]+1);
         best1 = A_ind[s-3];
         best2 = A_ind[s-2];
         best3 = A_ind[s-1];         
      }
   }    

   A_ind[s-3] = best1; // Multiply A by the product of these 3 primes and store their indices
   A_ind[s-2] = best2;
   A_ind[s-1] = best3;   
   for (i = 0; i < 3; i++) 
   {
      msl = mpn_mul_1(A+1, A+1, A[0], factor_base[A_ind[s+i-3]].p);
      if (msl) 
      {
         A[A[0]+1] = msl;
         A[0]++;
      }
   }
              
#if POLY_A
   mpz_t A_disp, targ_A;
   mpz_init(A_disp);
   mpz_init(targ_A);
   fmpz_to_mpz(A_disp, A);
   fmpz_to_mpz(targ_A, target_A);
   gmp_printf("A = %Zd, target A = %Zd\n", A_disp, targ_A);
   mpz_clear(A_disp);
   mpz_clear(targ_A);
#endif    
 
   /*for (i = 0; i < s; i++)
   {
      p = factor_base[A_ind[i]].p;
      poly_inf->inv_p2[i] = z_precompute_inverse(p*p);
   } */

   fmpz_to_mpz(poly_inf->A_mpz, A);
   
   flint_stack_release(); // release current_A
   flint_stack_release(); // release diff
   flint_stack_release(); // release best_diff     
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
   unsigned long limbs = qs_inf->prec+1;
   unsigned long limbs2;
   unsigned long * A = poly_inf->A;
   unsigned long * B = poly_inf->B;
   unsigned long p, i; 
   unsigned long * temp1 = (unsigned long *) flint_stack_alloc(limbs);
   unsigned long temp;
   mp_limb_t msl;
   double pinv;
   
   for (i = 0; i < s; i++)
   {
      p = factor_base[A_ind[i]].p;
      pinv = z_precompute_inverse(p);
      mpn_divmod_1(temp1 + 1, A + 1, A[0], p);
      temp1[0] = A[0] - (temp1[A[0]] == 0); 
      A_modp[i] = (temp = mpn_mod_1(temp1 + 1, temp1[0], p));
      temp = z_invert(temp, p);
      temp = z_mulmod_precomp(temp, qs_inf->sqrts[A_ind[i]], p, pinv);
      if (temp > p/2) temp = p - temp;
      msl = mpn_mul_1(B_terms + i*limbs + 1, temp1 + 1, temp1[0], temp);
      if (msl) 
      {
         B_terms[i*limbs + temp1[0] + 1] = msl;
         B_terms[i*limbs] = temp1[0] + 1;
      }
      else B_terms[i*limbs] = temp1[0];
#if B_TERMS
      mpz_t temp;
      mpz_init(temp);
      fmpz_to_mpz(temp, B_terms + i*limbs);
      gmp_printf("B_%ld = %Zd\n", i, temp);
      mpz_clear(temp);
#endif
   }
   
   F_mpn_copy(B, B_terms, B_terms[0]+1);  // Set B to the sum of the B terms
   if (limbs > B_terms[0] + 1) F_mpn_clear(B + B_terms[0] + 1, limbs - B_terms[0] - 1);
   for (i = 1; i < s; i++)
   {
      limbs2 = B_terms[i*limbs];
      msl = mpn_add_n(B+1, B+1, B_terms + i*limbs + 1, limbs2);
      if (msl) mpn_add_1(B + limbs2 + 1, B + limbs2 + 1, limbs - limbs2 - 1, msl);
   }
   B[0] = limbs - 1;
   while (!B[B[0]] && B[0]) B[0]--;
#if B_TERMS
   mpz_t temp2;
   mpz_init(temp2);
   fmpz_to_mpz(temp2, B);
   gmp_printf("B = %Zd\n", temp2);
   mpz_clear(temp2);
#endif
   
   flint_stack_release(); // release temp1
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
   unsigned long * A = poly_inf->A;
   unsigned long * B = poly_inf->B;
   uint32_t * A_inv = poly_inf->A_inv;
   uint32_t ** A_inv2B = poly_inf->A_inv2B;
   unsigned long * B_terms = poly_inf->B_terms;
   uint32_t * soln1 = poly_inf->soln1;
   uint32_t* soln2 = poly_inf->soln2;
   uint32_t * sqrts = qs_inf->sqrts;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long sieve_size = qs_inf->sieve_size;
   unsigned long s = poly_inf->s;
   unsigned long p, temp;
   unsigned limbs = qs_inf->prec+1; 
   double pinv;
   
   unsigned long i;
   for (i = 2; i < num_primes; i++) // skip k and 2
   {
      p = factor_base[i].p;
      pinv = factor_base[i].pinv;
      
      A_inv[i] = z_invert(mpn_mod_1(A+1, A[0], p), p);
             
      unsigned long j;
      for (j = 0; j < s; j++)
      {
         temp = mpn_mod_1(B_terms + j*limbs + 1, B_terms[j*limbs], p);
         temp = z_mulmod_precomp(temp, A_inv[i], p, pinv);
         temp *= 2;
         if (temp >= p) temp -= p;
         A_inv2B[j][i] = temp;
      }
             
      temp = mpn_mod_1(B+1, B[0], p);
      temp = sqrts[i] + p - temp;
#if FLINT_BITS == 64
      temp *= A_inv[i];
#else
      temp = z_mulmod2_precomp(temp, A_inv[i], p, pinv);
#endif
      temp += sieve_size/2;
      soln1[i] = z_mod2_precomp(temp, p, pinv); // Consider using long_mod_precomp
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
   uint32_t * soln1 = poly_inf->soln1;
   uint32_t * soln2 = poly_inf->soln2;
   unsigned long p, D;
   unsigned long * n = qs_inf->n;
   unsigned long * B = poly_inf->B;
   unsigned long temp, temp2, B_modp2, index, p2; 
   unsigned long sieve_size = qs_inf->sieve_size;
   prime_t * factor_base = qs_inf->factor_base;
   double * inv_p2 = poly_inf->inv_p2;
   double pinv;
   
   unsigned long j;
   for (j = 0; j < s; j++)
   {
      index = A_ind[j];
/*      p = factor_base[index].p;
      p2 = p*p;
      pinv = factor_base[index].pinv;
      D = z_ll_mod_precomp(n[2], n[1], p*p, inv_p2[j]);    
      if ((long) B < 0) 
      {
         B_modp2 = z_mod2_precomp(-B, p2, inv_p2[j]);
         B_modp2 = p2 - B_modp2;
         if (B_modp2 == p2) B_modp2 = 0;
      } else
      B_modp2 = z_mod2_precomp(B, p2, inv_p2[j]);
      temp = B_modp2*A_modp[j];
      temp = z_mod2_precomp(temp, p, pinv); 
      temp2 = z_invert(temp, p);
      D -= (B_modp2*B_modp2);
      if ((long) D < 0) temp = -z_div2_precomp(-D, p, pinv);
      else temp = -z_div2_precomp(-D, p, pinv);
      temp *= temp2;
      temp += sieve_size/2;
      if ((long) temp < 0) 
      {
         temp = p - z_mod2_precomp(-temp, p, pinv);
         if (temp == p) temp = 0;
      }
      else temp = z_mod2_precomp(temp, p, pinv);
      soln1[index] = temp;*/
      soln2[index] = -1L;
   }
}       

/*=========================================================================
   Compute C:
 
   Function: Compute the C coefficient of the polynomial with the 
             current A and B values
 
==========================================================================*/

void compute_B_C(QS_t * qs_inf, poly_t * poly_inf)
{
   mpz_t * A_mpz = &poly_inf->A_mpz;
   mpz_t * B_mpz = &poly_inf->B_mpz;
   mpz_t * C = &poly_inf->C;
   unsigned long * B = poly_inf->B;
   mpz_t * mpz_n = &qs_inf->mpz_n;

   fmpz_to_mpz(*B_mpz, B);
   mpz_mul(*C, *B_mpz, *B_mpz);
   mpz_sub(*C, *C, *mpz_n);
#if TEST_C
   mpz_t temp;
   mpz_init(temp);
   mpz_mod(temp, *C, *A_mpz);
   if (mpz_cmp_ui(temp, 0) != 0) gmp_printf("B^2 - n = %Zd is not divisible by A = %Zd\n", *C, *A_mpz);
   mpz_clear(temp);
#endif   
   mpz_divexact(*C, *C, *A_mpz);
   qs_inf->sieve_fill = 128-mpz_sizeinbase(*C, 2)+qs_inf->error_bits+13;// 16, 20, 20
} 
