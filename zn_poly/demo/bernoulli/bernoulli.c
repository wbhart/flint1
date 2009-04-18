/*
   bernoulli.c:  example program; computes irregular indices for a prime p
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <stdio.h>
#include <gmp.h>
#include "zn_poly.h"


/*
   Finds distinct prime factors of n, in the range k < p <= n.
   Assumes that n does not have any prime factors p <= k.
   Stores them in increasing order starting at res.
   Returns number of factors written.
*/
unsigned
prime_factors_helper (ulong* res, ulong k, ulong n)
{
   if (n == 1)
      return 0;

   ulong i;
   for (i = k + 1; i * i <= n; i++)
   {
      if (n % i == 0)
      {
         // found a factor
         *res = i;
         // remove that factor entirely
         for (n /= i; n % i == 0; n /= i);
         return prime_factors_helper (res + 1, i, n) + 1;
      }
   }
   // no more factors
   *res = n;
   return 1;
}


/*
   Finds distinct prime factors of n.
   Writes them to res in increasing order, returns number of factors found.
*/
unsigned
prime_factors (ulong* res, ulong n)
{
   return prime_factors_helper (res, 1, n);
}


/*
   A really dumb primality test.
*/
int
is_prime (ulong n)
{
   ulong i;
   for (i = 2; i * i <= n; i++)
      if (n % i == 0)
         return 0;
   return 1;
}


/*
   Assuming p >= 3 is prime, find the minimum primitive root mod p.
   Stores it at g, and stores its inverse mod p at g_inv.
   
   Note: this is probably a terrible algorithm, but it's fine compared to the
   running time of bernoulli() for large p.
*/
void
primitive_root (ulong* g, ulong* g_inv, ulong p)
{
   zn_mod_t mod;
   zn_mod_init (mod, p);

   // find prime factors of p-1
   ulong factors[200];
   unsigned num_factors, i;
   num_factors = prime_factors (factors, p - 1);
   
   // loop through candidates
   ulong x;
   for (x = 2; ; x++)
   {
      ZNP_ASSERT (x < p);
      
      // it's a generator if x^((p-1)/q) != 1 for all primes q dividing p-1
      int good = 1;
      for (i = 0; i < num_factors && good; i++)
         good = good && (zn_mod_pow (x, (p - 1) / factors[i], mod) != 1);

      if (good)
      {
         *g = x;
         *g_inv = zn_mod_pow (x, p - 2, mod);
         zn_mod_clear (mod);
         return;
      }
   }
}


/*
   Computes bernoulli numbers B(0), B(2), ..., B((p-3)/2) mod p.
   
   If res != NULL, it stores the result in res, which needs to be of length
   (p-1)/2.
   
   If irregular != NULL, it stores the number of irregular indices at
   irregular[0], follows by the irregular indices. The largest permitted
   index of irregularity is irregular_max; if there are more indices than that
   to store, the program will abort.
   
   p must be a prime >= 3
   p must be < 2^(ULONG_BITS/2)
*/
void
bernoulli (ulong* res, ulong* irregular, unsigned irregular_max, ulong p)
{
   ZNP_ASSERT (p < (1UL << (ULONG_BITS/2)));
   
   ulong g, g_inv;
   primitive_root (&g, &g_inv, p);

   ulong n = (p-1) / 2;

   zn_mod_t mod;
   zn_mod_init (mod, p);
   
   // allocate our own res if the caller didn't
   int own_res = 0;
   if (!res)
   {
      res = (ulong*) malloc (n * sizeof (ulong));
      own_res = 1;
   }
   
   if (irregular)
      irregular[0] = 0;
   
   ulong* G = (ulong*) malloc (n * sizeof (ulong));
   ulong* P = (ulong*) malloc ((2*n - 1) * sizeof (ulong));
   ulong* J = res;
   
   // -------------------------------------------------------------------------
   // Step 1: compute polynomials G(X) and J(X)
   
   // g_pow = g^(i-1), g_pow = g^(-i) at beginning of each iteration
   ulong g_pow = g_inv;
   ulong g_pow_inv = 1;
   
   // bias = (g-1)/2 mod p
   ulong bias = (g - 1 + ((g & 1) ? 0 : p)) / 2;

   // fudge = g^(i^2), fudge_inv = g^(-i^2) at each iteration
   ulong fudge = 1;
   ulong fudge_inv = 1;
   
   ulong i;
   for (i = 0; i < n; i++)
   {
      ulong prod = g * g_pow;
      
      // quo = floor(prod / p)
      // rem = prod % p
      ulong quo = zn_mod_quotient (prod, mod);
      ulong rem = prod - quo * p;

      // h = h(g^i) / g^i mod p
      ulong h = g_pow_inv * zn_mod_sub_slim (bias, quo, mod);
      h = zn_mod_reduce (h, mod);
      
      // update g_pow and g_pow_inv for next iteration
      g_pow = rem;
      g_pow_inv = zn_mod_reduce (g_pow_inv * g_inv, mod);
      
      // X^i coefficient of G(X) is g^(i^2) * h(g^i) / g^i
      G[i] = zn_mod_reduce (h * fudge, mod);
      
      // X^i coefficient of J(X) is g^(-i^2)
      J[i] = fudge_inv;
      
      // update fudge and fudge_inv for next iteration
      fudge = zn_mod_reduce (fudge * g_pow, mod);
      fudge = zn_mod_reduce (fudge * g_pow, mod);
      fudge = zn_mod_reduce (fudge * g, mod);
      fudge_inv = zn_mod_reduce (fudge_inv * g_pow_inv, mod);
      fudge_inv = zn_mod_reduce (fudge_inv * g_pow_inv, mod);
      fudge_inv = zn_mod_reduce (fudge_inv * g, mod);
   }
   
   J[0] = 0;

   // -------------------------------------------------------------------------
   // Step 2: compute product P(X) = G(X) * J(X)
   
   zn_array_mul (P, J, n, G, n, mod);

   // -------------------------------------------------------------------------
   // Step 3: extract output from P(X), and verify result
   
   res[0] = 1;

   // we will verify that \sum_{j=0}^{(p-3)/2} 4^j (2j+1) B(2j) = -2 mod p
   ulong check_accum = 1;
   ulong check_four_pow = 4 % p;

   // g_sqr = g^2
   // g_sqr_inv = g^(-2)
   ulong g_sqr = zn_mod_reduce (g * g, mod);
   ulong g_sqr_inv = zn_mod_reduce (g_inv * g_inv, mod);
   
   // make table with J[i] = (1 - g^(2i+2))(1 - g^(2i+4)) ... (1 - g^(p-3))
   // for 0 <= i < (p-1)/2
   ulong g_sqr_inv_pow = g_sqr_inv;
   J[n-1] = 1;
   for (i = 1; i < n; i++)
   {
      J[n-i-1] = zn_mod_reduce (J[n-i] * (p + 1 - g_sqr_inv_pow), mod);
      g_sqr_inv_pow = zn_mod_reduce (g_sqr_inv_pow * g_sqr_inv, mod);
   }
   
   // fudge = g^(i^2) at each iteration
   fudge = g;
   // g_sqr_pow = g^(2i) at each iteration
   ulong g_sqr_pow = g_sqr;
   
   // prod_inv = [(1 - g^(2i))(1 - g^(2i+2)) ... (1 - g^(p-3))]^(-1)
   // at each iteration (todo: for i == 1, it's experimentally equal to -1/2
   // mod p, need to prove this)
   ulong prod_inv = p - 2;

   for (i = 1; i < n; i++)
   {
      ulong val = (i == (n-1)) ? 0 : P[i + n];
      if (n & 1)
         val = zn_mod_neg (val, mod);
      val = zn_mod_add_slim (val, G[i], mod);
      val = zn_mod_add_slim (val, P[i], mod);

      // multiply by 4 * i * g^(i^2)
      val = zn_mod_reduce (val * fudge, mod);
      val = zn_mod_reduce (val * (2*i), mod);
      val = zn_mod_add_slim (val, val, mod);
      
      // divide by (1 - g^(2i))
      val = zn_mod_reduce (val * prod_inv, mod);
      val = zn_mod_reduce (val * J[i], mod);
      prod_inv = zn_mod_reduce (prod_inv * (1 + p - g_sqr_pow), mod);

      // store output coefficient if requested
      if (!own_res)
         res[i] = val;

      // store irregular index if requested
      if (irregular)
      {
         if (val == 0)
         {
            irregular[0]++;
            if (irregular[0] >= irregular_max)
            {
               printf ("too many irregular indices for p = %lu\n", p);
               abort ();
            }
            irregular[irregular[0]] = 2*i;
         }
      }
         
      // update fudge and g_sqr_pow
      g_sqr_pow = zn_mod_reduce (g_sqr_pow * g, mod);
      fudge = zn_mod_reduce (fudge * g_sqr_pow, mod);
      g_sqr_pow = zn_mod_reduce (g_sqr_pow * g, mod);
      
      // update verification data
      ulong check_term = zn_mod_reduce (check_four_pow * (2*i + 1), mod);
      check_term = zn_mod_reduce (check_term * val, mod);
      check_accum = zn_mod_add_slim (check_accum, check_term, mod);
      check_four_pow = zn_mod_add_slim (check_four_pow, check_four_pow, mod);
      check_four_pow = zn_mod_add_slim (check_four_pow, check_four_pow, mod);
   }

   if (check_accum != p-2)
   {
      printf ("bernoulli failed correctness check for p = %lu\n", p);
      abort ();
   }

   if (own_res)
      free (res);

   free (P);
   free (G);
   zn_mod_clear (mod);
}



/*
   Same as bernoulli(), but handles two primes simultaneously.
   
   p1 and p2 must be distinct primes >= 3.
*/
void
bernoulli_dual(ulong* res1, ulong* irregular1, unsigned irregular1_max, ulong p1,
               ulong* res2, ulong* irregular2, unsigned irregular2_max, ulong p2)
{
   ZNP_ASSERT (p1 < (1UL << (ULONG_BITS/2)));
   ZNP_ASSERT (p2 < (1UL << (ULONG_BITS/2)));
   ZNP_ASSERT (p1 != p2);

   // swap them to make p1 < p2
   if (p1 > p2)
   {
      { ulong temp = p2; p2 = p1; p1 = temp; }
      { ulong* temp = res1; res1 = res2; res2 = temp; }
      { ulong* temp = irregular1; irregular1 = irregular2; irregular2 = temp; }
      { unsigned temp = irregular1_max; irregular1_max = irregular2_max; irregular2_max = temp; }
   }
   
   ulong g1, g_inv1;
   ulong g2, g_inv2;
   primitive_root (&g1, &g_inv1, p1);
   primitive_root (&g2, &g_inv2, p2);
   
   ulong n1 = (p1-1) / 2;
   ulong n2 = (p2-1) / 2;

   zn_mod_t mod1, mod2;
   zn_mod_init (mod1, p1);
   zn_mod_init (mod2, p2);
   
   ulong q = p1 * p2;
   zn_mod_t mod;
   zn_mod_init (mod, q);

   // allocate our own res2 if the caller didn't
   int own_res2 = 0;
   if (!res2)
   {
      res2 = (ulong*) malloc (n2 * sizeof (ulong));
      own_res2 = 1;
   }
   
   if (irregular1)
      irregular1[0] = 0;
   if (irregular2)
      irregular2[0] = 0;
   
   // find idempotents to CRT modulo p1 and p2, i.e.
   // id1 = 1 mod p1, id1 = 0 mod p2
   // id2 = 0 mod p1, id2 = 1 mod p2
   mpz_t p1_mpz, p2_mpz, a1_mpz, a2_mpz, g_mpz;
   mpz_init (p1_mpz);
   mpz_init (p2_mpz);
   mpz_init (a1_mpz);
   mpz_init (a2_mpz);
   mpz_init (g_mpz);

   mpz_set_ui (p1_mpz, p1);
   mpz_set_ui (p2_mpz, p2);
   mpz_gcdext (g_mpz, a1_mpz, a2_mpz, p1_mpz, p2_mpz);

   mpz_mul (a1_mpz, a1_mpz, p1_mpz);
   mpz_mod_ui (a1_mpz, a1_mpz, q);
   ulong id2 = mpz_get_ui (a1_mpz);

   mpz_mul (a2_mpz, a2_mpz, p2_mpz);
   mpz_mod_ui (a2_mpz, a2_mpz, q);
   ulong id1 = mpz_get_ui (a2_mpz);

   mpz_clear (g_mpz);
   mpz_clear (a2_mpz);
   mpz_clear (a1_mpz);
   mpz_clear (p2_mpz);
   mpz_clear (p1_mpz);
   
   ulong* G = (ulong*) malloc (n2 * sizeof (ulong));
   ulong* P = (ulong*) malloc ((2*n2 - 1) * sizeof (ulong));
   ulong* J = res2;
   
   // -------------------------------------------------------------------------
   // Step 1: compute polynomials G(X) and J(X)
   
   // g_pow = g^(i-1), g_pow = g^(-i) at beginning of each iteration
   ulong g_pow1 = g_inv1;
   ulong g_pow2 = g_inv2;
   ulong g_pow_inv1 = 1;
   ulong g_pow_inv2 = 1;
   
   // bias = (g-1)/2 mod p
   ulong bias1 = (g1 - 1 + ((g1 & 1) ? 0 : p1)) / 2;
   ulong bias2 = (g2 - 1 + ((g2 & 1) ? 0 : p2)) / 2;

   // fudge = g^(i^2), fudge_inv = g^(-i^2) at each iteration
   ulong fudge1 = 1;
   ulong fudge2 = 1;
   ulong fudge_inv1 = 1;
   ulong fudge_inv2 = 1;
   
   ulong i;
   for (i = 0; i < n1; i++)
   {
      ulong prod1 = g1 * g_pow1;
      ulong prod2 = g2 * g_pow2;

      // quo = floor(prod / p)
      // rem = prod % p
      ulong quo1 = zn_mod_quotient (prod1, mod1);
      ulong rem1 = prod1 - quo1 * p1;
      ulong quo2 = zn_mod_quotient (prod2, mod2);
      ulong rem2 = prod2 - quo2 * p2;

      // h = h(g^i) / g^i mod p
      ulong h1 = g_pow_inv1 * zn_mod_sub_slim (bias1, quo1, mod1);
      h1 = zn_mod_reduce (h1, mod1);
      ulong h2 = g_pow_inv2 * zn_mod_sub_slim (bias2, quo2, mod2);
      h2 = zn_mod_reduce (h2, mod2);
      
      // update g_pow and g_pow_inv for next iteration
      g_pow1 = rem1;
      g_pow_inv1 = zn_mod_reduce (g_pow_inv1 * g_inv1, mod1);
      g_pow2 = rem2;
      g_pow_inv2 = zn_mod_reduce (g_pow_inv2 * g_inv2, mod2);
      
      // X^i coefficient of G(X) is g^(i^2) * h(g^i) / g^i
      // (combine via CRT)
      ulong Gval1 = zn_mod_reduce (h1 * fudge1, mod1);
      Gval1 = zn_mod_mul (Gval1, id1, mod);
      ulong Gval2 = zn_mod_reduce (h2 * fudge2, mod2);
      Gval2 = zn_mod_mul (Gval2, id2, mod);
      G[i] = zn_mod_add (Gval1, Gval2, mod);
      
      // X^i coefficient of J(X) is g^(-i^2)
      ulong Jval1 = zn_mod_mul (fudge_inv1, id1, mod);
      ulong Jval2 = zn_mod_mul (fudge_inv2, id2, mod);
      J[i] = zn_mod_add (Jval1, Jval2, mod);
      
      // update fudge and fudge_inv for next iteration
      fudge1 = zn_mod_reduce (fudge1 * g_pow1, mod1);
      fudge1 = zn_mod_reduce (fudge1 * g_pow1, mod1);
      fudge1 = zn_mod_reduce (fudge1 * g1, mod1);
      fudge_inv1 = zn_mod_reduce (fudge_inv1 * g_pow_inv1, mod1);
      fudge_inv1 = zn_mod_reduce (fudge_inv1 * g_pow_inv1, mod1);
      fudge_inv1 = zn_mod_reduce (fudge_inv1 * g1, mod1);

      fudge2 = zn_mod_reduce (fudge2 * g_pow2, mod2);
      fudge2 = zn_mod_reduce (fudge2 * g_pow2, mod2);
      fudge2 = zn_mod_reduce (fudge2 * g2, mod2);
      fudge_inv2 = zn_mod_reduce (fudge_inv2 * g_pow_inv2, mod2);
      fudge_inv2 = zn_mod_reduce (fudge_inv2 * g_pow_inv2, mod2);
      fudge_inv2 = zn_mod_reduce (fudge_inv2 * g2, mod2);
   }
   
   // finished with p1, now finish the loop for p2
   for (; i < n2; i++)
   {
      ulong prod2 = g2 * g_pow2;
      
      // quo = floor(prod / p)
      // rem = prod % p
      ulong quo2 = zn_mod_quotient (prod2, mod2);
      ulong rem2 = prod2 - quo2 * p2;

      // h = h(g^i) / g^i mod p
      ulong h2 = g_pow_inv2 * zn_mod_sub_slim (bias2, quo2, mod2);
      h2 = zn_mod_reduce (h2, mod2);
      
      // update g_pow and g_pow_inv for next iteration
      g_pow2 = rem2;
      g_pow_inv2 = zn_mod_reduce (g_pow_inv2 * g_inv2, mod2);
      
      // X^i coefficient of G(X) is g^(i^2) * h(g^i) / g^i
      // (combine via CRT)
      ulong Gval2 = zn_mod_reduce (h2 * fudge2, mod2);
      G[i] = zn_mod_mul (Gval2, id2, mod);
      
      // X^i coefficient of J(X) is g^(-i^2)
      J[i] = zn_mod_mul (fudge_inv2, id2, mod);
      
      // update fudge and fudge_inv for next iteration
      fudge2 = zn_mod_reduce (fudge2 * g_pow2, mod2);
      fudge2 = zn_mod_reduce (fudge2 * g_pow2, mod2);
      fudge2 = zn_mod_reduce (fudge2 * g2, mod2);
      fudge_inv2 = zn_mod_reduce (fudge_inv2 * g_pow_inv2, mod2);
      fudge_inv2 = zn_mod_reduce (fudge_inv2 * g_pow_inv2, mod2);
      fudge_inv2 = zn_mod_reduce (fudge_inv2 * g2, mod2);
   }
   
   J[0] = 0;

   // -------------------------------------------------------------------------
   // Step 2: compute product P(X) = G(X) * J(X)
   
#if 0
   zn_array_mul_fft_dft (P, J, n2, G, n2, 3, mod);
#else
   zn_array_mul (P, J, n2, G, n2, mod);
#endif

   // -------------------------------------------------------------------------
   // Step 3: extract output from P(X), and verify result

   if (res1)
      res1[0] = 1;

   // we will verify that \sum_{j=0}^{(p-3)/2} 4^j (2j+1) B(2j) = -2 mod p
   ulong check_accum1 = 1;
   ulong check_four_pow1 = 4 % p1;

   // g_sqr = g^2
   // g_sqr_inv = g^(-2)
   ulong g_sqr1 = zn_mod_reduce (g1 * g1, mod1);
   ulong g_sqr_inv1 = zn_mod_reduce (g_inv1 * g_inv1, mod1);
   
   // make table with J[i] = (1 - g^(2i+2))(1 - g^(2i+4)) ... (1 - g^(p-3))
   // for 0 <= i < (p-1)/2
   ulong g_sqr_inv_pow1 = g_sqr_inv1;
   J[n1-1] = 1;
   for (i = 1; i < n1; i++)
   {
      J[n1-i-1] = zn_mod_reduce (J[n1-i] * (p1 + 1 - g_sqr_inv_pow1), mod1);
      g_sqr_inv_pow1 = zn_mod_reduce (g_sqr_inv_pow1 * g_sqr_inv1, mod1);
   }
   
   // fudge = g^(i^2) at each iteration
   fudge1 = g1;
   // g_sqr_pow = g^(2i) at each iteration
   ulong g_sqr_pow1 = g_sqr1;
   
   // prod_inv = [(1 - g^(2i))(1 - g^(2i+2)) ... (1 - g^(p-3))]^(-1)
   // at each iteration (todo: for i == 1, it's experimentally equal to -1/2
   // mod p, need to prove this)
   ulong prod_inv1 = p1 - 2;

   for (i = 1; i < n1; i++)
   {
      ulong val = (i == (n1-1)) ? 0 : P[i + n1];
      if (n1 & 1)
         val = zn_mod_neg (val, mod);
      val = zn_mod_add (val, G[i], mod);
      val = zn_mod_add (val, P[i], mod);
      
      // reduce it mod p1
      val = zn_mod_reduce (val, mod1);

      // multiply by 4 * i * g^(i^2)
      val = zn_mod_reduce (val * fudge1, mod1);
      val = zn_mod_reduce (val * (2*i), mod1);
      val = zn_mod_add_slim (val, val, mod1);
      
      // divide by (1 - g^(2i))
      val = zn_mod_reduce (val * prod_inv1, mod1);
      val = zn_mod_reduce (val * J[i], mod1);
      prod_inv1 = zn_mod_reduce (prod_inv1 * (1 + p1 - g_sqr_pow1), mod1);

      // store output coefficient if requested
      if (res1)
         res1[i] = val;

      // store irregular index if requested
      if (irregular1)
      {
         if (val == 0)
         {
            irregular1[0]++;
            if (irregular1[0] >= irregular1_max)
            {
               printf ("too many irregular indices for p = %lu\n", p1);
               abort ();
            }
            irregular1[irregular1[0]] = 2*i;
         }
      }
         
      // update fudge and g_sqr_pow
      g_sqr_pow1 = zn_mod_reduce (g_sqr_pow1 * g1, mod1);
      fudge1 = zn_mod_reduce (fudge1 * g_sqr_pow1, mod1);
      g_sqr_pow1 = zn_mod_reduce (g_sqr_pow1 * g1, mod1);
      
      // update verification data
      ulong check_term1 = zn_mod_reduce (check_four_pow1 * (2*i + 1), mod1);
      check_term1 = zn_mod_reduce (check_term1 * val, mod1);
      check_accum1 = zn_mod_add_slim (check_accum1, check_term1, mod1);
      check_four_pow1 = zn_mod_add_slim (check_four_pow1, check_four_pow1, mod1);
      check_four_pow1 = zn_mod_add_slim (check_four_pow1, check_four_pow1, mod1);
   }

   if (check_accum1 != p1-2)
   {
      printf ("bernoulli_dual failed correctness check for p1 = %lu\n", p1);
      abort ();
   }

   // -------------------------------------------------------------------------
   // Do step 3 again for the second prime

   res2[0] = 1;

   // we will verify that \sum_{j=0}^{(p-3)/2} 4^j (2j+1) B(2j) = -2 mod p
   ulong check_accum2 = 1;
   ulong check_four_pow2 = 4 % p2;

   // g_sqr = g^2
   // g_sqr_inv = g^(-2)
   ulong g_sqr2 = zn_mod_reduce (g2 * g2, mod2);
   ulong g_sqr_inv2 = zn_mod_reduce (g_inv2 * g_inv2, mod2);
   
   // make table with J[i] = (1 - g^(2i+2))(1 - g^(2i+4)) ... (1 - g^(p-3))
   // for 0 <= i < (p-1)/2
   ulong g_sqr_inv_pow2 = g_sqr_inv2;
   J[n2-1] = 1;
   for (i = 1; i < n2; i++)
   {
      J[n2-i-1] = zn_mod_reduce (J[n2-i] * (p2 + 1 - g_sqr_inv_pow2), mod2);
      g_sqr_inv_pow2 = zn_mod_reduce (g_sqr_inv_pow2 * g_sqr_inv2, mod2);
   }
   
   // fudge = g^(i^2) at each iteration
   fudge2 = g2;
   // g_sqr_pow = g^(2i) at each iteration
   ulong g_sqr_pow2 = g_sqr2;
   
   // prod_inv = [(1 - g^(2i))(1 - g^(2i+2)) ... (1 - g^(p-3))]^(-1)
   // at each iteration (todo: for i == 1, it's experimentally equal to -1/2
   // mod p, need to prove this)
   ulong prod_inv2 = p2 - 2;

   for (i = 1; i < n2; i++)
   {
      ulong val = (i == (n2-1)) ? 0 : P[i + n2];
      if (n2 & 1)
         val = zn_mod_neg (val, mod);
      val = zn_mod_add (val, G[i], mod);
      val = zn_mod_add (val, P[i], mod);
      
      // reduce it mod p2
      val = zn_mod_reduce (val, mod2);

      // multiply by 4 * i * g^(i^2)
      val = zn_mod_reduce (val * fudge2, mod2);
      val = zn_mod_reduce (val * (2*i), mod2);
      val = zn_mod_add_slim (val, val, mod2);
      
      // divide by (1 - g^(2i))
      val = zn_mod_reduce (val * prod_inv2, mod2);
      val = zn_mod_reduce (val * J[i], mod2);
      prod_inv2 = zn_mod_reduce (prod_inv2 * (1 + p2 - g_sqr_pow2), mod2);

      // store output coefficient if requested
      if (!own_res2)
         res2[i] = val;

      // store irregular index if requested
      if (irregular2)
      {
         if (val == 0)
         {
            irregular2[0]++;
            if (irregular2[0] >= irregular2_max)
            {
               printf ("too many irregular indices for p = %lu\n", p2);
               abort ();
            }
            irregular2[irregular2[0]] = 2*i;
         }
      }
         
      // update fudge and g_sqr_pow
      g_sqr_pow2 = zn_mod_reduce (g_sqr_pow2 * g2, mod2);
      fudge2 = zn_mod_reduce (fudge2 * g_sqr_pow2, mod2);
      g_sqr_pow2 = zn_mod_reduce (g_sqr_pow2 * g2, mod2);
      
      // update verification data
      ulong check_term2 = zn_mod_reduce (check_four_pow2 * (2*i + 1), mod2);
      check_term2 = zn_mod_reduce (check_term2 * val, mod2);
      check_accum2 = zn_mod_add_slim (check_accum2, check_term2, mod2);
      check_four_pow2 = zn_mod_add_slim (check_four_pow2, check_four_pow2, mod2);
      check_four_pow2 = zn_mod_add_slim (check_four_pow2, check_four_pow2, mod2);
   }

   if (check_accum2 != p2-2)
   {
      printf ("bernoulli_dual failed correctness check for p2 = %lu\n", p2);
      abort ();
   }

   if (own_res2)
      free (res2);
   free (P);
   free (G);

   zn_mod_clear (mod);
   zn_mod_clear (mod2);
   zn_mod_clear (mod1);
}



int
main (int argc, char* argv[])
{
   if (argc == 2)
   {
      ulong i, p = atol (argv[1]);
      ulong irregular[30];
      bernoulli (NULL, irregular, 29, p);

      printf ("irregular indices for p = %lu: ", p);
      for (i = 1; i <= irregular[0]; i++)
         printf ("%lu ", irregular[i]);
      printf ("\n");
   }
   else if (argc == 3)
   {
      ulong i, p1 = atol (argv[1]), p2 = atol (argv[2]);
      ulong irregular1[30];
      ulong irregular2[30];

      bernoulli_dual (NULL, irregular1, 29, p1, NULL, irregular2, 29, p2);

      printf ("irregular indices for p = %lu: ", p1);
      for (i = 1; i <= irregular1[0]; i++)
         printf ("%lu ", irregular1[i]);
      printf ("\n");

      printf ("irregular indices for p = %lu: ", p2);
      for (i = 1; i <= irregular2[0]; i++)
         printf ("%lu ", irregular2[i]);
      printf ("\n");
   }
   else
   {
      printf ("usage:\n");
      printf ("\n");
      printf ("  bernoulli p\n");
      printf ("    prints irregular indices for p\n");
      printf ("\n");
      printf ("  bernoulli p1 p2\n");
      printf ("    prints irregular indices for p1 and p2\n");
   }

   return 0;
}


// end of file ****************************************************************
