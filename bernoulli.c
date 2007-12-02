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
/****************************************************************************

   bernoulli.c: Finds Bernoulli numbers B_{2k}
                Based on the implementation in SAGE written by David Harvey
                Uses mpz_polys for the calculations.  See bernoulli_fmpz.c
                and bernoulli_zmod.c for use of other polys.
   
   Copyright (C) 2007, David Howden

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "flint.h"
#include "mpz_poly.h"
#include "long_extras.h"

#define TRUE 1;
#define FALSE 0;


/*
   Computes the bernoulli numbers B_0, B_2, ..., B_{p-3}
   for prime p
   
   Requires that res be allocated for (p-1)/2 unsigned longs
   which will hold the result.
   
   If returns 0, then the factoring of p has failed, otherwise
   will always return 1.
*/

int bernoulli_mod_p_mpz(unsigned long *res, unsigned long p)
{
   FLINT_ASSERT(p > 2);
   FLINT_ASSERT(z_isprime(p) == 1);
   
   unsigned long g, g_inv, g_sqr, g_sqr_inv;
   double p_inv = z_precompute_inverse(p);
   g = z_primitive_root(p);
   
   if(!g)
   {
      return FALSE;
   }
   
   g_inv = z_invert(g, p);
   g_sqr = z_mulmod_precomp(g, g, p, p_inv);
   g_sqr_inv = z_mulmod_precomp(g_inv, g_inv, p, p_inv);
   
   unsigned long poly_size = (p-1)/2;
   
   int is_odd = poly_size % 2;
   
   unsigned long g_power, g_power_inv;
   g_power = g_inv;
   g_power_inv = 1;
   
   // constant is (g-1)/2 mod p
   unsigned long constant;
   if(g % 2)
   {
      constant = (g-1)/2;
   }
   else
   {
      constant = (g+p-1)/2;
   }
   
   // fudge holds g^{i^2}, fudge_inv holds g^{-i^2}
   unsigned long fudge, fudge_inv;
   fudge = fudge_inv = 1;
   
   // compute the polynomials F(X) and G(X)
   mpz_poly_t F, G;
   
   mpz_poly_init2(F, poly_size);
   mpz_poly_init2(G, poly_size);
   
   unsigned long i, temp, h;
   
   for(i = 0; i < poly_size; i++)
   {  
      // compute h(g^i)/g^i (h(x) is as in latex notes)
      temp = g * g_power;
            
      h = z_mulmod_precomp(p + constant - (temp / p), g_power_inv, p, p_inv);
      
      g_power = z_mod_precomp(temp, p, p_inv);
      g_power_inv = z_mulmod_precomp(g_power_inv, g_inv, p, p_inv);
      
      // store coefficient g^{i^2} h(g^i)/g^i
      mpz_poly_set_coeff_ui(G, i, z_mulmod_precomp(h, fudge, p, p_inv));
      mpz_poly_set_coeff_ui(F, i, fudge_inv);
      
      // update fudge and fudge_inv
      fudge = z_mulmod_precomp(z_mulmod_precomp(fudge, g_power, p, p_inv), z_mulmod_precomp(g_power, g, p, p_inv), p, p_inv);
      fudge_inv = z_mulmod_precomp(z_mulmod_precomp(fudge_inv, g_power_inv, p, p_inv), z_mulmod_precomp(g_power_inv, g, p, p_inv), p, p_inv);
   }
   
   mpz_poly_set_coeff_ui(F, 0, 0);
   
   // step 2: multiply the polynomials...
   mpz_poly_t product;
   mpz_poly_init(product);
   mpz_poly_mul(product, G, F);
   
   // step 3: assemble the result...   
   unsigned long g_sqr_power, value;
   g_sqr_power = g_sqr;
   fudge = g;

   res[0] = 1;
   
   mpz_t value_coeff;
   mpz_init(value_coeff);
   
   unsigned long value_coeff_ui;

   for(i = 1; i < poly_size; i++)
   {
      mpz_poly_get_coeff(value_coeff, product, i + poly_size);
      value = mpz_fdiv_ui(value_coeff, p);
      
      value = z_mod_precomp(mpz_poly_get_coeff_ui(product, i + poly_size), p, p_inv);
      
      mpz_poly_get_coeff(value_coeff, product, i);
      if(is_odd)
      {
         value = z_mod_precomp(mpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + p - value, p, p_inv);
      }
      else
      {
         value = z_mod_precomp(mpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + value, p, p_inv);
      }
      
      value = z_mulmod_precomp(z_mulmod_precomp(z_mulmod_precomp(4, i, p, p_inv), fudge, p, p_inv), value, p, p_inv);
      value = z_mulmod_precomp(value, z_invert(p+1-g_sqr_power, p), p, p_inv);

      res[i] = value;
      
      g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
      fudge = z_mulmod_precomp(fudge, g_sqr_power, p, p_inv);
      g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
   }
   
   mpz_clear(value_coeff);
   
   mpz_poly_clear(F);
   mpz_poly_clear(G);
   mpz_poly_clear(product);
   
   return TRUE;
}


/*
   Verifies that the ouput of bernoulli_mod_p above is correct.
   
   Takes the result from bernoulli_mod_p (res - an array of (p-1)/2
   unsigned longs), and the prime p.
   
   Returns 0 if res is incorrect, 1 if res is correct.
*/

int verify_bernoulli_mod_p(unsigned long *res, unsigned long p)
{
   unsigned long N, i, product, sum, value, element;
   double p_inv;
   N = (p-1)/2;
   product = 1;
   sum = 0;
   
   p_inv = z_precompute_inverse(p);
   
   for(i = 0; i < N; i++)
   {
      element = res[i];
      value = z_mulmod_precomp(z_mulmod_precomp(product, 2*i+1, p, p_inv), element, p, p_inv);
      sum = z_mod_precomp(sum + value, p, p_inv);
      product = z_mulmod_precomp(product, 4, p, p_inv);
   }
   
   if(z_mod_precomp(sum + 2,  p, p_inv))
   {   
      return FALSE;
   }
   
   return TRUE;
}


/*
   Test function for bernoulli_mod_p
   
   Calculates bernoulli_mod_p for the prime p and verifies the result.
   
   Returs 0 if incorrect, and 1 if correct.
*/

int test_bernoulli_mod_p(unsigned long p)
{
   unsigned long *res = (unsigned long*) flint_stack_alloc((p-1)/2);
   if(!bernoulli_mod_p_mpz(res, p))
   {
      printf("Could not factor p = %d\n", p);
      flint_stack_release();
      return FALSE;
   }
   int result = verify_bernoulli_mod_p(res, p);
   flint_stack_release();
   return result;
}


int main (int argc, char const *argv[])
{
   unsigned long p = 2;
   unsigned long tests = 1000;
   unsigned long fail = 0;
   
   for(unsigned long i = 0; i < tests; i++)
   {
      p = z_nextprime(p);
      if(!test_bernoulli_mod_p(p))
      {
         printf("Fails on p = %d\n", p);
         fail++;
      }
      else
      {
         printf("Works on p = %d\n", p);
      }
   }
   
   printf("\nResults: %d OK, %d FAILED.\n", tests - fail, fail);

   return 0;
}
