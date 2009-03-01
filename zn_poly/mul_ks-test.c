/*
   mul_ks-test.c:  test code for functions in mul_ks.c
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
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
   Tests zn_array_mul_KSk, for given len1, len2, reduction algorithm, modulus.
   1 <= k <= 4 indicates which KS variant to call.
   Returns 1 on success.
*/
int testcase_zn_array_mul_KS(int k, size_t len1, size_t len2, int redc,
                             zn_mod_t mod)
{
   // disallow REDC if modulus is even
   if (!(mod->n & 1))
      redc = 0;

   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len1);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * len2);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * (len1 + len2 - 1));
   ulong* compare = (ulong*) malloc(sizeof(ulong) * (len1 + len2 - 1));
   
   // generate random polys
   size_t i;
   for (i = 0; i < len1; i++)
      buf1[i] = random_ulong(mod->n);
   for (i = 0; i < len2; i++)
      buf2[i] = random_ulong(mod->n);
      
   // compare target implementation against reference implementation
   ref_zn_array_mul(correct, buf1, len1, buf2, len2, mod);

   switch (k)
   {
      case 1:
         zn_array_mul_KS1(compare, buf1, len1, buf2, len2, redc, mod); break;
      case 2:
         zn_array_mul_KS2(compare, buf1, len1, buf2, len2, redc, mod); break;
      case 3:
         zn_array_mul_KS3(compare, buf1, len1, buf2, len2, redc, mod); break;
      case 4:
         zn_array_mul_KS4(compare, buf1, len1, buf2, len2, redc, mod); break;
      default: printf("oops!\n"); abort();
   }

   if (redc)
      // correct for REDC reduction
      ref_zn_array_scalar_mul(compare, compare, len1 + len2 - 1,
                              mod->n - mod->B, mod);

   int success = !zn_array_cmp(correct, compare, len1 + len2 - 1);
   
   free(compare);
   free(correct);
   free(buf2);
   free(buf1);
   
   return success;
}


/*
   tests zn_array_mul_KSk() on a range of input cases,
   where 1 <= k <= 4
*/
int test_zn_array_mul_KSk(unsigned k)
{
   int success = 1;
   int bits, trial, redc;
   size_t len1, len2;

   // first try a dense range of "small" problems

   for (bits = 2; bits <= ULONG_BITS && success; bits++)
   for (len2 = 1; len2 <= 30 && success; len2++)
   for (len1 = len2; len1 <= 30 && success; len1++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < 10 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 0));

      success = success && testcase_zn_array_mul_KS(k, len1, len2, redc, mod);

      zn_mod_clear(mod);
   }
   
   // now try some random larger problems

   for (bits = 2; bits <= ULONG_BITS && success; bits++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < 200 && success; trial++)
   {
      len1 = random_ulong(1000) + 1;
      len2 = random_ulong(1000) + 1;
      
      if (len1 < len2)
      {
         ulong temp = len1;
         len1 = len2;
         len2 = temp;
      }
      
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 0));

      success = success && testcase_zn_array_mul_KS(k, len1, len2, redc, mod);
      
      zn_mod_clear(mod);
   }

   return success;
}


int test_zn_array_mul_KS1() { return test_zn_array_mul_KSk(1); }
int test_zn_array_mul_KS2() { return test_zn_array_mul_KSk(2); }
int test_zn_array_mul_KS3() { return test_zn_array_mul_KSk(3); }
int test_zn_array_mul_KS4() { return test_zn_array_mul_KSk(4); }



/*
   Tests zn_array_mul_KSk for squaring, for given len, reduction algorithm,
   modulus.
   1 <= k <= 4 indicates which KS variant to call.
   Returns 1 on success.
*/
int testcase_zn_array_sqr_KS(int k, size_t len, int redc, zn_mod_t mod)
{
   // disallow REDC if modulus is even
   if (!(mod->n & 1))
      redc = 0;

   ulong* buf = (ulong*) malloc(sizeof(ulong) * len);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * (2 * len - 1));
   ulong* compare = (ulong*) malloc(sizeof(ulong) * (2 * len - 1));

   // generate random poly
   size_t i;
   for (i = 0; i < len; i++)
      buf[i] = random_ulong(mod->n);
      
   // compare target implementation against reference implementation
   ref_zn_array_mul(correct, buf, len, buf, len, mod);

   switch (k)
   {
      case 1: zn_array_mul_KS1(compare, buf, len, buf, len, redc, mod); break;
      case 2: zn_array_mul_KS2(compare, buf, len, buf, len, redc, mod); break;
      case 3: zn_array_mul_KS3(compare, buf, len, buf, len, redc, mod); break;
      case 4: zn_array_mul_KS4(compare, buf, len, buf, len, redc, mod); break;
      default: printf("oops!\n"); abort();
   }

   if (redc)
      // correct for REDC reduction
      ref_zn_array_scalar_mul(compare, compare, 2 * len - 1,
                              mod->n - mod->B, mod);

   int success = !zn_array_cmp(correct, compare, 2 * len - 1);
   
   free(compare);
   free(correct);
   free(buf);
   
   return success;
}



/*
   tests zn_array_mul_KSk() for squaring on a range of input cases,
   where 1 <= k <= 4
*/
int test_zn_array_sqr_KSk(unsigned k)
{
   int success = 1;
   int bits, trial, redc;
   size_t len;

   // first try a dense range of "small" problems

   for (bits = 2; bits <= ULONG_BITS && success; bits++)
   for (len = 1; len <= 30 && success; len++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < 10 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 0));
      
      success = success && testcase_zn_array_sqr_KS(k, len, redc, mod);
      
      zn_mod_clear(mod);
   }
   
   // now try some random larger problems

   for (bits = 2; bits <= ULONG_BITS && success; bits++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < 200 && success; trial++)
   {
      len = random_ulong(1000) + 1;
      
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 0));

      success = success && testcase_zn_array_sqr_KS(k, len, redc, mod);
      
      zn_mod_clear(mod);
   }

   return success;
}


int test_zn_array_sqr_KS1() { return test_zn_array_sqr_KSk(1); }
int test_zn_array_sqr_KS2() { return test_zn_array_sqr_KSk(2); }
int test_zn_array_sqr_KS3() { return test_zn_array_sqr_KSk(3); }
int test_zn_array_sqr_KS4() { return test_zn_array_sqr_KSk(4); }


/*
   Tests zn_array_recip_fix_reduce() for given len, bits, reduction algorithm
   and modulus.
   
   Doesn't test the _skip_ parameter.
   
   Note: running time is quadratic in len, so don't make it too big.
*/
int testcase_zn_array_recip_fix_reduce(size_t len, unsigned bits,
                                       int redc, const zn_mod_t mod)
{
   // disallow REDC if modulus is even
   if (!(mod->n & 1))
      redc = 0;

   ZNP_ASSERT(bits >= 1 && 2*bits <= 3*ULONG_BITS);

   mpz_t* a;
   size_t i;

   a = (mpz_t*) malloc(sizeof(mpz_t) * len);
   for (i = 0; i < len; i++)
      mpz_init(a[i]);

   // c = 2^bits - 1
   mpz_t c;
   mpz_init(c);
   mpz_set_ui(c, 1);
   mpz_mul_2exp(c, c, bits);
   mpz_sub_ui(c, c, 1);
   
   mpz_t hi, lo;
   mpz_init(hi);
   mpz_init(lo);
   
   mpz_t temp;
   mpz_init(temp);
   
   // a "small" integer, no more than c
   mpz_t small;
   mpz_init(small);
   mpz_set_ui(small, (bits >= 2) ? 3 : 1);
   ZNP_ASSERT(mpz_cmp(small, c) <= 0);

   mpz_t sum1, sum2;
   mpz_init(sum1);
   mpz_init(sum2);

   // make up a list of a[i]'s
   for (i = 0; i < len; i++)
   {
      // make up low digit
      switch (random_ulong(3))
      {
         case 0:
            // some uniform random digit
            mpz_urandomb(lo, randstate, bits);
            break;
            
         case 1:
            // a value close to zero
            mpz_urandomm(lo, randstate, small);
            break;
            
         case 2:
            // a value close to the maximum
            // (anything up to and including 2^bits - 1)
            mpz_urandomm(lo, randstate, small);
            mpz_sub(lo, c, lo);
            break;
      }
      
      // make up high digit
      switch (random_ulong(3))
      {
         case 0:
            // some uniform random digit
            mpz_urandomm(hi, randstate, c);
            break;
            
         case 1:
            // a value close to zero
            mpz_urandomm(hi, randstate, small);
            break;
            
         case 2:
            // a value close to the maximum
            // (anything up to but NOT including 2^bits - 1)
            mpz_urandomm(hi, randstate, small);
            mpz_sub(hi, c, hi);
            mpz_sub_ui(hi, hi, 1);
            break;
      }
      
      // put a[i] = hi*B + lo
      mpz_mul_2exp(a[i], hi, bits);
      mpz_add(a[i], a[i], lo);
   }

   // construct the sums in forward and reverse directions
   // i.e. sum1 = a[0] + a[1]*B + ... + a[len-1]*B^(len-1)
   //      sum2 = a[len-1] + a[len-2]*B + ... + a[0]*B^(len-1).
   for (i = 0; i < len; i++)
   {
      mpz_mul_2exp(sum1, sum1, bits);
      mpz_add(sum1, sum1, a[len-1-i]);
      mpz_mul_2exp(sum2, sum2, bits);
      mpz_add(sum2, sum2, a[i]);
   }

   // decompose both sums into sequence of (len+1) base-B digits
   unsigned words = CEIL_DIV(bits, ULONG_BITS);
   ZNP_ASSERT(words <= 2);
   ulong* digits1 = (ulong*) malloc(sizeof(ulong) * words * (len + 1));
   ulong* digits2 = (ulong*) malloc(sizeof(ulong) * words * (len + 1));

   if (words == 1)
   {
      for (i = 0; i <= len; i++)
      {
         mpz_tdiv_r_2exp(temp, sum1, bits);
         digits1[i] = mpz_get_ui(temp);
         mpz_tdiv_q_2exp(sum1, sum1, bits);

         mpz_tdiv_r_2exp(temp, sum2, bits);
         digits2[i] = mpz_get_ui(temp);
         mpz_tdiv_q_2exp(sum2, sum2, bits);
      }
   }
   else
   {
      for (i = 0; i <= len; i++)
      {
         mpz_tdiv_r_2exp(temp, sum1, ULONG_BITS);
         digits1[2*i] = mpz_get_ui(temp);
         mpz_tdiv_q_2exp(sum1, sum1, ULONG_BITS);
         mpz_tdiv_r_2exp(temp, sum1, bits - ULONG_BITS);
         digits1[2*i+1] = mpz_get_ui(temp);
         mpz_tdiv_q_2exp(sum1, sum1, bits - ULONG_BITS);

         mpz_tdiv_r_2exp(temp, sum2, ULONG_BITS);
         digits2[2*i] = mpz_get_ui(temp);
         mpz_tdiv_q_2exp(sum2, sum2, ULONG_BITS);
         mpz_tdiv_r_2exp(temp, sum2, bits - ULONG_BITS);
         digits2[2*i+1] = mpz_get_ui(temp);
         mpz_tdiv_q_2exp(sum2, sum2, bits - ULONG_BITS);
      }
   }
   
   // shouldn't be any bits left
   ZNP_ASSERT(mpz_cmp_ui(sum1, 0) == 0);
   ZNP_ASSERT(mpz_cmp_ui(sum2, 0) == 0);

   // see if zn_array_recip_fix_reduce() returns the original inputs (mod n)
   ulong* b = (ulong*) malloc(sizeof(ulong) * len);
   zn_array_recip_fix_reduce(b, 1, digits1, digits2, len, bits, redc, mod);
   int success = 1;
   for (i = 0; i < len; i++)
   {
      if (redc)
         // correct for REDC reduction
         mpz_mul_ui(temp, a[i], mod->n - mod->B);
      else
         mpz_set(temp, a[i]);

      mpz_mod_ui(temp, a[i], mod->n);
      success = success && (mpz_get_ui(temp) == b[i]);
   }

   // clean up
   free(b);
   free(digits2);
   free(digits1);
   mpz_clear(temp);
   mpz_clear(sum2);
   mpz_clear(sum1);
   mpz_clear(lo);
   mpz_clear(hi);
   mpz_clear(small);
   mpz_clear(c);
   for (i = 0; i < len; i++)
      mpz_clear(a[i]);
   free(a);
   
   return 1;
}


/*
   Tests zn_array_recip_fix_reduce() on a range of small problems.
*/
int test_zn_array_recip_fix_reduce()
{
   int success = 1;
   int bits, trial, redc;
   size_t len;

   for (bits = 1; 2*bits <= 3*ULONG_BITS && success; bits++)
   for (len = 1; len <= 15 && success; len++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < 200 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(random_ulong(ULONG_BITS - 1) + 2, 0));

      success = success &&
                     testcase_zn_array_recip_fix_reduce(len, bits, redc, mod);

      zn_mod_clear(mod);
   }
   
   return success;
}



// end of file ****************************************************************
