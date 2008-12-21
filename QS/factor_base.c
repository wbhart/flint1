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

 factor_base.c
 
 Routines for generating and maintaining the factor base primes
 including the multiplier

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <math.h>

#include "../flint.h"
#include "../memory-manager.h"
#include "../long_extras.h"

#include "factor_base.h"

/*=========================================================================
   num_FB_primes:
 
   Function: retrieve the number of factor base primes to use from table
             
 
==========================================================================*/

unsigned long tiny_num_FB_primes(unsigned long bits)
{
   unsigned long i;
   
   for (i = 0; i < PTABSIZE_SMALL; i++)
   {
      if (prime_tab_small[i][0] > bits) break;
   }
   
   return prime_tab_small[i-1][1];
}

/*=========================================================================
   sqrts_init:
 
   Function: allocate space for the factorbase primes and associated info
             
 
==========================================================================*/

void tiny_sqrts_init(QS_t * qs_inf)
{
   qs_inf->sqrts = (uint32_t *) flint_stack_alloc_bytes(sizeof(uint32_t)*qs_inf->num_primes);
}

void tiny_sqrts_clear(void)
{
   flint_stack_release();
}

/*=========================================================================
   primes_init:
 
   Function: allocate space for the factorbase primes and associated info
             
 
==========================================================================*/

void tiny_primes_init(QS_t * qs_inf)
{
   unsigned long bits = qs_inf->bits; // set bits to the number of bits of kn
   
   qs_inf->num_primes = tiny_num_FB_primes(bits); 
   
   qs_inf->factor_base = (prime_t *) flint_stack_alloc_bytes(qs_inf->num_primes*sizeof(prime_t));
}

void tiny_primes_clear(void)
{
   flint_stack_release();
}

/*===========================================================================
   Compute Prime Sizes:
 
   Function: Computes the size in bits of each prime in the factor base
 
===========================================================================*/
void tiny_compute_sizes(QS_t * qs_inf)
{
     unsigned long num_primes = qs_inf->num_primes;
     
     qs_inf->sizes = (unsigned char *) flint_stack_alloc_bytes(num_primes);
     unsigned char * sizes = qs_inf->sizes;
     prime_t * factor_base = qs_inf->factor_base;
     
     for (unsigned long i = 0; i < num_primes; i++)
     {
         sizes[i] = (unsigned char) round(log(factor_base[i].p)/log(2.0));
     }
     
     return;
}

void tiny_sizes_clear(void)
{
   flint_stack_release();
}

/*=========================================================================
   Knuth-Schroeppel algorithm:
 
   Function: Find the best multiplier to use (allows 2 as a multiplier).
             The general idea is to find a multiplier k such that kn will
             be faster to factor. This is achieved by making kn a square 
             modulo lots of small primes. These primes will then be factor
             base primes, and the more small factor base primes, the faster
             relations will accumulate, since they hit the sieving interval
             more often. 
             
             Also computes approximate inverses and modular square roots 
             primes that are suitable as factor base primes
 
==========================================================================*/

unsigned long tiny_knuth_schroeppel(QS_t * qs_inf)
{
    float best_factor = -10.0f;
    unsigned long multiplier = 1;
    unsigned long nmod8, mod8, multindex, prime, nmod, mult;
    const unsigned long max_fb_primes = qs_inf->num_primes;
    unsigned long fb_prime = 2; // leave space for the multiplier and 2
    float factors[NUMMULTS];
    float logpdivp;
    double pinv;
    int kron;
    
    uint32_t * sqrts = qs_inf->sqrts;
    
    fmpz_t n = qs_inf->n;
    nmod8 = n[1]%8;
    
    mpz_t r;
        
    for (multindex = 0; multindex < NUMMULTS; multindex++)
    {
       mod8 = ((nmod8*multipliers[multindex])%8);
       factors[multindex] = 0.34657359; // ln2/2 
       if (mod8 == 1) factors[multindex] *= 4.0;   
       if (mod8 == 5) factors[multindex] *= 2.0;   
       factors[multindex] -= (log((float) multipliers[multindex]) / 2.0);
    }
    
    prime = 3;
    while ((prime < KSMAX) && (fb_prime < max_fb_primes))
    {
          pinv = z_precompute_inverse(prime);
          logpdivp = log((float)prime) / (float)prime; // log p / p
          nmod = z_ll_mod_precomp(n[2], n[1], prime, pinv); 
          if (nmod == 0) return prime;
          kron = z_jacobi(nmod, prime); 
          for (multindex = 0; multindex < NUMMULTS; multindex++)
          {
              mult = multipliers[multindex];
              if (mult >= prime) 
              {
                 if (mult >= prime*prime) mult = mult%prime; 
                 else mult = z_mod_precomp(mult, prime, pinv);
              }
              if (mult == 0) factors[multindex] += logpdivp;
              else if (kron*z_jacobi(mult, prime) == 1) 
                 factors[multindex] += 2.0*logpdivp;
          }
          
          prime = z_nextprime(prime, 0);
    }
    
    for (multindex=0; multindex<NUMMULTS; multindex++)
    {
      if (factors[multindex] > best_factor)
      { 
        best_factor = factors[multindex];
        multiplier = multipliers[multindex];
      }
    } 
    
    qs_inf->k = multiplier;
    return 0;
}

/*=========================================================================
   Compute Factor Base:
 
   Function: Compute all the primes p for which n is a quadratic residue 
             mod p. Compute square roots of n modulo each p.
 
==========================================================================*/


unsigned long tiny_compute_factor_base(QS_t * qs_inf)
{
   unsigned long fb_prime = 2;
   unsigned long multiplier = qs_inf->k;
   prime_t * factor_base = qs_inf->factor_base;
   uint32_t * sqrts = qs_inf->sqrts;
   unsigned long num_primes = tiny_num_FB_primes(qs_inf->bits);
   unsigned long prime, nmod;
   double pinv;
   fmpz_t n = qs_inf->n;
   long kron;
    
   factor_base[0].p = multiplier;
   factor_base[0].pinv = z_precompute_inverse(multiplier);
   factor_base[1].p = 2;
   prime = 2;
   
   while (fb_prime < num_primes)
   {
      prime = z_nextprime(prime, 0);
      pinv = z_precompute_inverse(prime);
      nmod = z_ll_mod_precomp(n[2], n[1], prime, pinv); 
      if (nmod == 0) 
      {
         if (z_mod_precomp(multiplier, prime, pinv) != 0) return prime;
      }
      kron = z_jacobi(nmod, prime); 
      if (kron == 1)
      {
         factor_base[fb_prime].p = prime;
         factor_base[fb_prime].pinv = pinv;
         sqrts[fb_prime] = z_sqrtmod(nmod, prime);
         fb_prime++;
      }   
   }
   
   qs_inf->num_primes = fb_prime;
   return 0;
}
