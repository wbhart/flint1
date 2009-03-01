/*
   ref_mul.c:  reference zn_array_mul implementation
   
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
   Sets x = op[0] + op[1]*B + ... + op[len-1]*B^(len-1), where B = 2^bits.
   
   Running time is soft-linear in output length.
*/
void pack(mpz_t x, const ulong* op, size_t len, unsigned bits)
{
   ZNP_ASSERT(len >= 1);

   if (len == 1)
   {
      // base case
      mpz_set_ui(x, op[0]);
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_t y;
      mpz_init(y);
      pack(x, op, len/2, bits);
      pack(y, op + len/2, len - len/2, bits);
      mpz_mul_2exp(y, y, (len/2) * bits);
      mpz_add(x, x, y);
      mpz_clear(y);
   }
}



/*
   Inverse operation of pack(), with output coefficients reduced mod n.
   
   Running time is soft-linear in output length.
*/
void unpack(ulong* res, const mpz_t op, size_t len, unsigned bits, ulong n)
{
   ZNP_ASSERT(len >= 1);
   ZNP_ASSERT(mpz_sizeinbase(op, 2) <= len * bits);

   mpz_t y;
   mpz_init(y);
   
   if (len == 1)
   {
      // base case
      mpz_set(y, op);
      mpz_fdiv_r_ui(y, y, n);
      *res = mpz_get_ui(y);
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_tdiv_q_2exp(y, op, (len/2) * bits);
      unpack(res + len/2, y, len - len/2, bits, n);
      mpz_tdiv_r_2exp(y, op, (len/2) * bits);
      unpack(res, y, len/2, bits, n);
   }
   
   mpz_clear(y);
}



/*
   Reference implementation of zn_array_mul().
   Very simple kronecker substitution, uses GMP for multiplication.
*/
void ref_zn_array_mul(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   
   mpz_t x, y;
   mpz_init(x);
   mpz_init(y);
   
   unsigned bits = 2 * mod->bits + ceil_lg(len2);
   unsigned words = CEIL_DIV(bits, ULONG_BITS);
   
   pack(x, op1, len1, bits);
   pack(y, op2, len2, bits);
   mpz_mul(x, x, y);
   unpack(res, x, len1 + len2 - 1, bits, mod->n);
   
   mpz_clear(y);
   mpz_clear(x);
}



/*
   Reference implementation of zn_array_midmul().
   Just calls ref_zn_array_mul() and extracts relevant part of output.
*/
void ref_zn_array_midmul(ulong* res, const ulong* op1, size_t len1,
                         const ulong* op2, size_t len2, const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   
   ulong* temp = (ulong*) malloc((len1 + len2 - 1) * sizeof(ulong));
   ref_zn_array_mul(temp, op1, len1, op2, len2, mod);
   ulong i;
   for (i = 0; i < len1 - len2 + 1; i++)
      res[i] = temp[i + len2 - 1];
   free(temp);
}



/*
   Reference implementation of negacyclic multiplication.
   
   Multiplies op1[0, len) by op2[0, len) negacyclically, puts result
   into res[0, len).
*/
void ref_zn_array_negamul(ulong* res, const ulong* op1, const ulong* op2,
                          size_t len, const zn_mod_t mod)
{
   ulong* temp = (ulong*) malloc(sizeof(ulong) * 2 * len);
   
   ref_zn_array_mul(temp, op1, len, op2, len, mod);
   temp[2*len - 1] = 0;
   
   mpz_t x;
   mpz_init(x);

   size_t i;
   for (i = 0; i < len; i++)
   {
      mpz_set_ui(x, temp[i]);
      mpz_sub_ui(x, x, temp[i + len]);
      mpz_mod_ui(x, x, mod->n);
      res[i] = mpz_get_ui(x);
   }
   
   mpz_clear(x);
   free(temp);
}


/*
   Reference implementation of scalar multiplication.
   
   Multiplies op[0, len) by x, puts result in res[0, len).
*/
void ref_zn_array_scalar_mul(ulong* res, const ulong* op, size_t len,
                             ulong x, const zn_mod_t mod)
{
   mpz_t y;
   mpz_init(y);
   
   size_t i;
   for (i = 0; i < len; i++)
   {
      mpz_set_ui(y, op[i]);
      mpz_mul_ui(y, y, x);
      mpz_mod_ui(y, y, mod->n);
      res[i] = mpz_get_ui(y);
   }
   
   mpz_clear(y);
}


// end of file ****************************************************************
