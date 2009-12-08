/*
   pack-test.c:  test code for functions in pack.c
   
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
   Helper function for ref_zn_array_pack().

   Sets x = 2^k * (op[0] + op[1]*2^b + ... + op[n-1]*2^((n-1)*b)).
   
   Running time is soft-linear in output length.
*/
void
ref_zn_array_pack_helper (mpz_t x, const ulong* op, size_t n, unsigned b,
                          unsigned k)
{
   ZNP_ASSERT (n >= 1);

   if (n == 1)
   {
      // base case
      mpz_set_ui (x, op[0]);
      mpz_mul_2exp (x, x, k);
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_t y;
      mpz_init (y);
      ref_zn_array_pack_helper (x, op, n / 2, b, k);
      ref_zn_array_pack_helper (y, op + n / 2, n - n / 2, b, 0);
      mpz_mul_2exp (y, y, (n / 2) * b + k);
      mpz_add (x, x, y);
      mpz_clear (y);
   }
}


/*
   Reference implementation of zn_array_pack().
   
   (doesn't take into account the s or r parameters)
*/
void
ref_zn_array_pack (mp_limb_t* res, const ulong* op, size_t n, unsigned b,
                   unsigned k)
{
   mpz_t x;
   mpz_init (x);
   ref_zn_array_pack_helper (x, op, n, b, k);
   mpz_to_mpn (res, CEIL_DIV (n * b + k, GMP_NUMB_BITS), x);
   mpz_clear (x);
}


/*
   Helper function for ref_zn_array_unpack().

   Inverse operation of ref_zn_array_pack_helper(); each output coefficient
   occupies ceil(b / ULONG_BITS) ulongs.
   
   Running time is soft-linear in output length.
*/
void
ref_zn_array_unpack_helper (ulong* res, const mpz_t op, size_t n, unsigned b,
                            unsigned k)
{
   ZNP_ASSERT (n >= 1);
   ZNP_ASSERT (mpz_sizeinbase (op, 2) <= n * b + k);

   unsigned w = CEIL_DIV (b, ULONG_BITS);

   mpz_t y;
   mpz_init (y);
   
   if (n == 1)
   {
      // base case
      unsigned i;
      mpz_tdiv_q_2exp (y, op, k);
      for (i = 0; i < w; i++)
      {
         res[i] = mpz_get_ui (y);
         mpz_tdiv_q_2exp (y, y, ULONG_BITS);
      }
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_tdiv_q_2exp (y, op, (n / 2) * b + k);
      ref_zn_array_unpack_helper (res + w * (n / 2), y, n - n / 2, b, 0);
      mpz_tdiv_r_2exp (y, op, (n / 2) * b + k);
      ref_zn_array_unpack_helper (res, y, n / 2, b, k);
   }
   
   mpz_clear (y);
}


/*
   Reference implementation of zn_array_unpack().
*/
void
ref_zn_array_unpack (ulong* res, const mp_limb_t* op, size_t n, unsigned b,
                     unsigned k)
{
   mpz_t x;
   mpz_init (x);
   mpn_to_mpz (x, op, CEIL_DIV (n * b + k, GMP_NUMB_BITS));
   ref_zn_array_unpack_helper (res, x, n, b, k);
   mpz_clear (x);
}



/*
   tests zn_array_pack() once for given n, b, k
*/
int
testcase_zn_array_pack (size_t n, unsigned b, unsigned k)
{
   ZNP_ASSERT (b >= 1);
   ZNP_ASSERT (n >= 1);

   int success = 1;
   ulong* in = (ulong*) malloc (sizeof (ulong) * n);

   size_t size = CEIL_DIV (n * b + k, GMP_NUMB_BITS);
   mp_limb_t* res = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (size + 2));
   mp_limb_t* ref = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (size + 2));

   // sentries to check buffer overflow
   res[0] = res[size + 1] = ref[0] = ref[size + 1] = 0x1234;

   // generate random data: at most b bits per input coefficient, possibly less
   unsigned rand_bits = (b >= ULONG_BITS) ? ULONG_BITS : b;
   rand_bits = random_ulong (rand_bits) + 1;
   ulong max = (rand_bits == ULONG_BITS)
                        ? ((ulong)(-1)) : ((1UL << rand_bits) - 1);
   size_t i;
   for (i = 0; i < n; i++)
      in[i] = random_ulong (max);

   // run target and reference implementation
   zn_array_pack (res + 1, in, n, 1, b, k, 0);
   ref_zn_array_pack (ref + 1, in, n, b, k);
   
   // check sentries
   success = success && (res[0] == 0x1234);
   success = success && (ref[0] == 0x1234);
   success = success && (res[size + 1] == 0x1234);
   success = success && (ref[size + 1] == 0x1234);
   // check correct result
   success = success && (mpn_cmp (res + 1, ref + 1, size) == 0);

   free (ref);
   free (res);
   free (in);
   
   return success;
}



/*
   tests zn_array_pack() on a range of input cases
*/
int
test_zn_array_pack (int quick)
{
   int success = 1;
   unsigned b, k;
   size_t n;

   for (b = 1; b < 3 * ULONG_BITS && success; b++)
   for (n = 1; n < (quick ? 100 : 200) && success;
        n += (quick ? (n < 5 ? 1 : 13) : 1))
   for (k = 0; k < 160; k += 20)
      success = success && testcase_zn_array_pack (n, b, k);
   
   return success;
}



/*
   tests zn_array_unpack() once for given n, b, k
*/
int
testcase_zn_array_unpack (size_t n, unsigned b, unsigned k)
{
   size_t buf_size = CEIL_DIV (n * b + k, GMP_NUMB_BITS);
   size_t size = n * CEIL_DIV (b, ULONG_BITS);

   mp_limb_t* buf = (mp_limb_t*) malloc (sizeof (mp_limb_t) * buf_size);
   ulong* res = (ulong*) malloc (sizeof (ulong) * (size + 2));
   ulong* ref = (ulong*) malloc (sizeof (ulong) * (size + 2));

   // sentries to check buffer overflow
   res[0] = res[size + 1] = ref[0] = ref[size + 1] = 0x1234;

   // generate random data
   mpz_t x;
   mpz_init (x);
   mpz_urandomb (x, randstate, n * b);
   mpz_mul_2exp (x, x, k);
   mpz_to_mpn (buf, buf_size, x);
   mpz_clear (x);

   // run target and reference implementation
   zn_array_unpack (res + 1, buf, n, b, k);
   ref_zn_array_unpack (ref + 1, buf, n, b, k);
   
   int success = 1;

   // check sentries
   success = success && (res[0] == 0x1234);
   success = success && (ref[0] == 0x1234);
   success = success && (res[size + 1] == 0x1234);
   success = success && (ref[size + 1] == 0x1234);
   // check correct result
   success = success && (zn_array_cmp (res + 1, ref + 1, size) == 0);
   
   free (ref);
   free (res);
   free (buf);

   return success;
}


/*
   tests zn_array_unpack() on a range of input cases
*/
int
test_zn_array_unpack (int quick)
{
   int success = 1;
   unsigned b, k;
   size_t n;

   for (b = 1; b < 3 * ULONG_BITS && success; b++)
   for (n = 1; n < (quick ? 100 : 200) && success;
        n += (quick ? (n < 5 ? 1 : 13) : 1))
   for (k = 0; k < 160; k += 19)
      success = success && testcase_zn_array_unpack (n, b, k);
   
   return success;
}


// end of file ****************************************************************
