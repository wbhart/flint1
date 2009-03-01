/*
   pack-test.c:  test code for functions in pack.c
   
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
   Helper function for ref_zn_array_pack().

   Sets x = 2^lead * (op[0] + op[1]*2^bits + ... + op[len-1]*2^((len-1)*bits)).
   
   Running time is soft-linear in output length.
*/
void ref_zn_array_pack_helper(mpz_t x, const ulong* op, size_t len,
                              unsigned bits, unsigned lead)
{
   ZNP_ASSERT(len >= 1);

   if (len == 1)
   {
      // base case
      mpz_set_ui(x, op[0]);
      mpz_mul_2exp(x, x, lead);
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_t y;
      mpz_init(y);
      ref_zn_array_pack_helper(x, op, len/2, bits, lead);
      ref_zn_array_pack_helper(y, op + len/2, len - len/2, bits, 0);
      mpz_mul_2exp(y, y, (len/2) * bits + lead);
      mpz_add(x, x, y);
      mpz_clear(y);
   }
}


/*
   Reference implementation of zn_array_pack().
   
   (doesn't take into account the _skip_ or _requested_ parameters)
*/
void ref_zn_array_pack(mp_limb_t* res, const ulong* op, size_t len,
                       unsigned bits, unsigned lead)
{
   mpz_t x;
   mpz_init(x);
   ref_zn_array_pack_helper(x, op, len, bits, lead);
   mpz_to_mpn(res, CEIL_DIV(len * bits + lead, GMP_NUMB_BITS), x);
   mpz_clear(x);
}


/*
   Helper function for ref_zn_array_unpack().

   Inverse operation of ref_zn_array_pack_helper(); each output coefficient
   occupies ceil(b / ULONG_BITS) ulongs.
   
   Running time is soft-linear in output length.
*/
void ref_zn_array_unpack_helper(ulong* res, const mpz_t op, size_t len,
                                unsigned bits, unsigned lead)
{
   ZNP_ASSERT(len >= 1);
   ZNP_ASSERT(mpz_sizeinbase(op, 2) <= len * bits + lead);

   unsigned words = CEIL_DIV(bits, ULONG_BITS);

   mpz_t y;
   mpz_init(y);
   
   if (len == 1)
   {
      // base case
      unsigned i;
      mpz_tdiv_q_2exp(y, op, lead);
      for (i = 0; i < words; i++)
      {
         res[i] = mpz_get_ui(y);
         mpz_tdiv_q_2exp(y, y, ULONG_BITS);
      }
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_tdiv_q_2exp(y, op, (len/2) * bits + lead);
      ref_zn_array_unpack_helper(res + words * (len/2), y, len - len/2,
                                 bits, 0);
      mpz_tdiv_r_2exp(y, op, (len/2) * bits + lead);
      ref_zn_array_unpack_helper(res, y, len/2, bits, lead);
   }
   
   mpz_clear(y);
}


/*
   Reference implementation of zn_array_unpack().
*/
void ref_zn_array_unpack(ulong* res, const mp_limb_t* op, size_t len,
                         unsigned bits, unsigned lead)
{
   mpz_t x;
   mpz_init(x);
   mpn_to_mpz(x, op, CEIL_DIV(len * bits + lead, GMP_NUMB_BITS));
   ref_zn_array_unpack_helper(res, x, len, bits, lead);
   mpz_clear(x);
}



/*
   tests zn_array_pack() once for given _len_ and _bits_ and _lead_
*/
int testcase_zn_array_pack(size_t len, unsigned bits, unsigned lead)
{
   ZNP_ASSERT(bits >= 1);
   ZNP_ASSERT(len >= 1);

   int success = 1;
   ulong* in = (ulong*) malloc(sizeof(ulong) * len);

   size_t out_size = CEIL_DIV(len * bits + lead, GMP_NUMB_BITS);
   mp_limb_t* out1 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * (out_size + 2));
   mp_limb_t* out2 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * (out_size + 2));

   // sentries to check buffer overflow
   out1[0] = out1[out_size + 1] = out2[0] = out2[out_size + 1] = 0x1234;

   // generate random data: at most b bits per input coefficient, possibly less
   unsigned rand_bits = (bits >= ULONG_BITS) ? ULONG_BITS : bits;
   rand_bits = random_ulong(rand_bits) + 1;
   ulong max = (rand_bits == ULONG_BITS)
                        ? ((ulong)(-1)) : ((1UL << rand_bits) - 1);
   size_t i;
   for (i = 0; i < len; i++)
      in[i] = random_ulong(max);

   // run target and reference implementation
   zn_array_pack(out1 + 1, in, len, 1, bits, lead, 0);
   ref_zn_array_pack(out2 + 1, in, len, bits, lead);
   
   // check sentries
   success = success && (out1[0] == 0x1234);
   success = success && (out2[0] == 0x1234);
   success = success && (out1[out_size + 1] == 0x1234);
   success = success && (out2[out_size + 1] == 0x1234);
   // check correct result
   success = success && (mpn_cmp(out1 + 1, out2 + 1, out_size) == 0);

   free(out2);
   free(out1);
   free(in);
   
   return success;
}



/*
   tests zn_array_pack() on a range of input cases
*/
int test_zn_array_pack()
{
   int success = 1;
   unsigned bits, lead;
   size_t len;

   for (bits = 1; bits < 3*ULONG_BITS && success; bits++)
   for (len = 1; len < 200 && success; len++)
   for (lead = 0; lead < 160; lead += 20)
      success = success && testcase_zn_array_pack(len, bits, lead);
   
   return success;
}



/*
   tests zn_array_unpack() once for given len and bits and lead
*/
int testcase_zn_array_unpack(size_t len, unsigned bits, unsigned lead)
{
   size_t buf_size = CEIL_DIV(len * bits + lead, GMP_NUMB_BITS);
   unsigned words = CEIL_DIV(bits, ULONG_BITS);
   size_t out_size = words * len;

   mp_limb_t* buf = (mp_limb_t*) malloc(sizeof(mp_limb_t) * buf_size);
   ulong* out1 = (ulong*) malloc(sizeof(ulong) * (out_size + 2));
   ulong* out2 = (ulong*) malloc(sizeof(ulong) * (out_size + 2));

   // sentries to check buffer overflow
   out1[0] = out1[out_size + 1] = out2[0] = out2[out_size + 1] = 0x1234;

   // generate random data
   mpz_t x;
   mpz_init(x);
   mpz_urandomb(x, randstate, len * bits);
   mpz_mul_2exp(x, x, lead);
   mpz_to_mpn(buf, buf_size, x);
   mpz_clear(x);

   // run target and reference implementation
   zn_array_unpack(out1 + 1, buf, len, bits, lead);
   ref_zn_array_unpack(out2 + 1, buf, len, bits, lead);
   
   int success = 1;

   // check sentries
   success = success && (out1[0] == 0x1234);
   success = success && (out2[0] == 0x1234);
   success = success && (out1[out_size + 1] == 0x1234);
   success = success && (out2[out_size + 1] == 0x1234);
   // check correct result
   success = success && (zn_array_cmp(out1 + 1, out2 + 1, out_size) == 0);
   
   free(out2);
   free(out1);
   free(buf);

   return success;
}


/*
   tests zn_array_unpack() on a range of input cases
*/
int test_zn_array_unpack()
{
   int success = 1;
   unsigned bits, lead;
   size_t len;

   for (bits = 1; bits < 3*ULONG_BITS && success; bits++)
   for (len = 1; len < 200 && success; len++)
   for (lead = 0; lead < 160; lead += 19)
      success = success && testcase_zn_array_unpack(len, bits, lead);
   
   return success;
}


// end of file ****************************************************************
