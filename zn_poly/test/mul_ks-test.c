/*
   mul_ks-test.c:  test code for functions in mul_ks.c
   
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

#include "support.h"
#include "zn_poly_internal.h"


/*
   Tests zn_array_mul_KSk, for given lengths, reduction algorithm, modulus.
   1 <= k <= 4 indicates which KS variant to call.
   sqr == 1 to test squaring (n2 is ignored).
   Returns 1 on success.
*/
int
testcase_zn_array_mul_KS (int k, size_t n1, size_t n2, int sqr, int redc,
                          const zn_mod_t mod)
{
   // disallow REDC if modulus is even
   if (!(mod->m & 1))
      redc = 0;
      
   if (sqr)
      n2 = n1;

   ulong* buf1 = (ulong*) malloc (sizeof(ulong) * n1);
   ulong* buf2 = sqr ? buf1 : (ulong*) malloc (sizeof(ulong) * n2);
   ulong* ref = (ulong*) malloc (sizeof(ulong) * (n1 + n2 - 1));
   ulong* res = (ulong*) malloc (sizeof(ulong) * (n1 + n2 - 1));
   
   // generate random polys
   size_t i;
   for (i = 0; i < n1; i++)
      buf1[i] = random_ulong (mod->m);
   if (!sqr)
      for (i = 0; i < n2; i++)
         buf2[i] = random_ulong (mod->m);
      
   // compare target implementation against reference implementation
   ref_zn_array_mul (ref, buf1, n1, buf2, n2, mod);

   switch (k)
   {
      case 1: zn_array_mul_KS1 (res, buf1, n1, buf2, n2, redc, mod); break;
      case 2: zn_array_mul_KS2 (res, buf1, n1, buf2, n2, redc, mod); break;
      case 3: zn_array_mul_KS3 (res, buf1, n1, buf2, n2, redc, mod); break;
      case 4: zn_array_mul_KS4 (res, buf1, n1, buf2, n2, redc, mod); break;
      default:
         printf ("oops!\n"); abort ();
   }

   if (redc)
      // correct for REDC reduction
      ref_zn_array_scalar_mul (res, res, n1 + n2 - 1, mod->m - mod->B, mod);

   int success = !zn_array_cmp (ref, res, n1 + n2 - 1);
   
   free (res);
   free (ref);
   if (!sqr)
      free (buf2);
   free (buf1);
   
   return success;
}


/*
   tests zn_array_mul_KSk() on a range of input cases, where 1 <= k <= 4
*/
int
test_zn_array_mul_KSk (unsigned k, int quick)
{
   int success = 1;
   int b, trial, redc;
   size_t n1, n2, t1, t2;
   zn_mod_t mod;

   // first try a dense range of "small" problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (n2 = 1; n2 <= 30 && success; n2 += (quick ? 5 : 1))
   for (n1 = n2; n1 <= 30 && success; n1 += (quick ? 5 : 1))
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 1 : 10) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_mul_KS (k, n1, n2, 0, redc, mod);
      zn_mod_clear (mod);
   }
   
   // now try some random larger problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 3 : 200) && success; trial++)
   {
      size_t t1 = random_ulong (quick ? 250 : 1000) + 1;
      size_t t2 = random_ulong (quick ? 250 : 1000) + 1;
      n1 = ZNP_MAX (t1, t2);
      n2 = ZNP_MIN (t1, t2);
      
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_mul_KS (k, n1, n2, 0, redc, mod);
      zn_mod_clear (mod);
   }

   return success;
}


int test_zn_array_mul_KS1 (int quick)
{
   return test_zn_array_mul_KSk (1, quick);
}

int test_zn_array_mul_KS2 (int quick)
{
   return test_zn_array_mul_KSk (2, quick);
}

int test_zn_array_mul_KS3 (int quick)
{
   return test_zn_array_mul_KSk (3, quick);
}

int test_zn_array_mul_KS4 (int quick)
{
   return test_zn_array_mul_KSk (4, quick);
}




/*
   tests zn_array_mul_KSk() for squaring on a range of input cases,
   where 1 <= k <= 4
*/
int
test_zn_array_sqr_KSk (unsigned k, int quick)
{
   int success = 1;
   int b, trial, redc;
   size_t n;
   zn_mod_t mod;

   // first try a dense range of "small" problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (n = 1; n <= 30 && success; n += (quick ? 5 : 1))
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 1 : 10) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_mul_KS (k, n, n, 1, redc, mod);
      zn_mod_clear(mod);
   }
   
   // now try some random larger problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 3 : 200) && success; trial++)
   {
      n = random_ulong (quick ? 250 : 1000) + 1;
      
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_mul_KS (k, n, n, 1, redc, mod);
      zn_mod_clear (mod);
   }

   return success;
}


int test_zn_array_sqr_KS1 (int quick)
{
   return test_zn_array_sqr_KSk (1, quick);
}

int test_zn_array_sqr_KS2 (int quick)
{
   return test_zn_array_sqr_KSk (2, quick);
}

int test_zn_array_sqr_KS3 (int quick)
{
   return test_zn_array_sqr_KSk (3, quick);
}

int test_zn_array_sqr_KS4 (int quick)
{
   return test_zn_array_sqr_KSk (4, quick);
}



/*
   Tests zn_array_recover_reduce() for given n, b, reduction algorithm
   and modulus.
   
   Doesn't test the s parameter.
   
   Note: running time is quadratic in n, so don't make it too big.
*/
int
testcase_zn_array_recover_reduce (size_t n, unsigned b, int redc,
                                  const zn_mod_t mod)
{
   // disallow REDC if modulus is even
   if (!(mod->m & 1))
      redc = 0;

   ZNP_ASSERT (b >= 1  &&  2 * b <= 3 * ULONG_BITS);

   mpz_t* a;
   size_t i;

   a = (mpz_t*) malloc (sizeof (mpz_t) * n);
   for (i = 0; i < n; i++)
      mpz_init (a[i]);

   // c = 2^b - 1
   mpz_t c;
   mpz_init (c);
   mpz_set_ui (c, 1);
   mpz_mul_2exp (c, c, b);
   mpz_sub_ui (c, c, 1);
   
   mpz_t hi, lo;
   mpz_init (hi);
   mpz_init (lo);
   
   mpz_t temp;
   mpz_init (temp);
   
   // a "small" integer, no more than c
   mpz_t small;
   mpz_init (small);
   mpz_set_ui (small, (b >= 2) ? 3 : 1);
   ZNP_ASSERT (mpz_cmp (small, c) <= 0);

   mpz_t sum1, sum2;
   mpz_init (sum1);
   mpz_init (sum2);

   // make up a list of a[i]'s
   for (i = 0; i < n; i++)
   {
      // make up low digit
      switch (random_ulong (3))
      {
         case 0:
            // some uniform random digit
            mpz_urandomb (lo, randstate, b);
            break;
            
         case 1:
            // a value close to zero
            mpz_urandomm (lo, randstate, small);
            break;
            
         case 2:
            // a value close to the maximum
            // (anything up to and including 2^b - 1)
            mpz_urandomm (lo, randstate, small);
            mpz_sub (lo, c, lo);
            break;
      }
      
      // make up high digit
      switch (random_ulong (3))
      {
         case 0:
            // some uniform random digit
            mpz_urandomm (hi, randstate, c);
            break;
            
         case 1:
            // a value close to zero
            mpz_urandomm (hi, randstate, small);
            break;
            
         case 2:
            // a value close to the maximum
            // (anything up to but NOT including 2^b - 1)
            mpz_urandomm (hi, randstate, small);
            mpz_sub (hi, c, hi);
            mpz_sub_ui (hi, hi, 1);
            break;
      }
      
      // put a[i] = hi*B + lo
      mpz_mul_2exp (a[i], hi, b);
      mpz_add (a[i], a[i], lo);
   }

   // construct the sums in forward and reverse directions
   // i.e. sum1 = a[0] + a[1]*B + ... + a[n-1]*B^(n-1)
   //      sum2 = a[n-1] + a[n-2]*B + ... + a[0]*B^(n-1).
   for (i = 0; i < n; i++)
   {
      mpz_mul_2exp (sum1, sum1, b);
      mpz_add (sum1, sum1, a[n - 1 - i]);
      mpz_mul_2exp (sum2, sum2, b);
      mpz_add (sum2, sum2, a[i]);
   }

   // decompose both sums into sequence of (n+1) base-B digits
   unsigned w = CEIL_DIV (b, ULONG_BITS);
   ZNP_ASSERT (w <= 2);
   ulong* d1 = (ulong*) malloc (sizeof (ulong) * w * (n + 1));
   ulong* d2 = (ulong*) malloc (sizeof (ulong) * w * (n + 1));

   if (w == 1)
   {
      for (i = 0; i <= n; i++)
      {
         mpz_tdiv_r_2exp (temp, sum1, b);
         d1[i] = mpz_get_ui (temp);
         mpz_tdiv_q_2exp (sum1, sum1, b);

         mpz_tdiv_r_2exp (temp, sum2, b);
         d2[i] = mpz_get_ui (temp);
         mpz_tdiv_q_2exp (sum2, sum2, b);
      }
   }
   else
   {
      for (i = 0; i <= n; i++)
      {
         mpz_tdiv_r_2exp (temp, sum1, ULONG_BITS);
         d1[2 * i] = mpz_get_ui (temp);
         mpz_tdiv_q_2exp (sum1, sum1, ULONG_BITS);
         mpz_tdiv_r_2exp (temp, sum1, b - ULONG_BITS);
         d1[2 * i + 1] = mpz_get_ui (temp);
         mpz_tdiv_q_2exp (sum1, sum1, b - ULONG_BITS);

         mpz_tdiv_r_2exp (temp, sum2, ULONG_BITS);
         d2[2 * i] = mpz_get_ui (temp);
         mpz_tdiv_q_2exp (sum2, sum2, ULONG_BITS);
         mpz_tdiv_r_2exp (temp, sum2, b - ULONG_BITS);
         d2[2 * i + 1] = mpz_get_ui (temp);
         mpz_tdiv_q_2exp (sum2, sum2, b - ULONG_BITS);
      }
   }
   
   // shouldn't be any bits left
   ZNP_ASSERT (mpz_cmp_ui (sum1, 0) == 0);
   ZNP_ASSERT (mpz_cmp_ui (sum2, 0) == 0);

   // see if zn_array_recover_reduce() returns the original inputs (mod m)
   ulong* res = (ulong*) malloc (sizeof (ulong) * n);
   zn_array_recover_reduce (res, 1, d1, d2, n, b, redc, mod);

   int success = 1;
   for (i = 0; i < n; i++)
   {
      if (redc)
         // correct for REDC reduction
         mpz_mul_ui (temp, a[i], mod->m - mod->B);
      else
         mpz_set (temp, a[i]);

      mpz_mod_ui (temp, a[i], mod->m);
      success = success && (mpz_get_ui (temp) == res[i]);
   }

   // clean up
   free (res);
   free (d2);
   free (d1);
   mpz_clear (temp);
   mpz_clear (sum2);
   mpz_clear (sum1);
   mpz_clear (lo);
   mpz_clear (hi);
   mpz_clear (small);
   mpz_clear (c);
   for (i = 0; i < n; i++)
      mpz_clear (a[i]);
   free (a);
   
   return 1;
}


/*
   Tests zn_array_recover_reduce() on a range of small problems.
*/
int
test_zn_array_recover_reduce (int quick)
{
   int success = 1;
   int b, trial, redc;
   size_t n;
   zn_mod_t mod;

   for (b = 1; 2 * b <= 3 * ULONG_BITS && success; b++)
   for (n = 1; n <= 15 && success; n++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 10 : 200) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (random_ulong (ULONG_BITS - 1) + 2, 0));
      success = success && testcase_zn_array_recover_reduce (n, b, redc, mod);
      zn_mod_clear (mod);
   }
   
   return success;
}



// end of file ****************************************************************
