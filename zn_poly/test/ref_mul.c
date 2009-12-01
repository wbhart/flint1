/*
   ref_mul.c:  reference implementations for polynomial multiplication,
               middle product, scalar multiplication, integer middle product
   
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
#include <string.h>


/*
   Sets x = op[0] + op[1]*B + ... + op[n-1]*B^(n-1), where B = 2^b.
   
   Running time is soft-linear in output length.
*/
void
pack (mpz_t x, const ulong* op, size_t n, unsigned b)
{
   ZNP_ASSERT (n >= 1);

   if (n == 1)
   {
      // base case
      mpz_set_ui (x, op[0]);
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_t y;
      mpz_init (y);
      pack (x, op, n / 2, b);
      pack (y, op + n / 2, n - n / 2, b);
      mpz_mul_2exp (y, y, (n / 2) * b);
      mpz_add (x, x, y);
      mpz_clear (y);
   }
}



/*
   Inverse operation of pack(), with output coefficients reduced mod m.
   
   Running time is soft-linear in output length.
*/
void
unpack (ulong* res, const mpz_t op, size_t n, unsigned b, ulong m)
{
   ZNP_ASSERT (n >= 1);
   ZNP_ASSERT (mpz_sizeinbase (op, 2) <= n * b);

   mpz_t y;
   mpz_init(y);
   
   if (n == 1)
   {
      // base case
      mpz_set (y, op);
      mpz_fdiv_r_ui (y, y, m);
      *res = mpz_get_ui (y);
   }
   else
   {
      // recursively split into top and bottom halves
      mpz_tdiv_q_2exp (y, op, (n / 2) * b);
      unpack (res + n / 2, y, n - n / 2, b, m);
      mpz_tdiv_r_2exp (y, op, (n / 2) * b);
      unpack (res, y, n / 2, b, m);
   }
   
   mpz_clear (y);
}



/*
   Reference implementation of zn_array_mul().
   Very simple Kronecker substitution, uses GMP for multiplication.
*/
void
ref_zn_array_mul (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   
   mpz_t x, y;
   mpz_init (x);
   mpz_init (y);
   
   unsigned b = 2 * mod->bits + ceil_lg (n2);
   unsigned words = CEIL_DIV (b, ULONG_BITS);
   
   pack (x, op1, n1, b);
   pack (y, op2, n2, b);
   mpz_mul (x, x, y);
   unpack (res, x, n1 + n2 - 1, b, mod->m);
   
   mpz_clear (y);
   mpz_clear (x);
}



/*
   Reference implementation of zn_array_mulmid().
   Just calls ref_zn_array_mul() and extracts relevant part of output.
*/
void
ref_zn_array_mulmid (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   
   ulong* temp = (ulong*) malloc ((n1 + n2 - 1) * sizeof (ulong));
   ref_zn_array_mul (temp, op1, n1, op2, n2, mod);
   ulong i;
   for (i = 0; i < n1 - n2 + 1; i++)
      res[i] = temp[i + n2 - 1];
   free (temp);
}



/*
   Reference implementation of negacyclic multiplication.
   
   Multiplies op1[0, n) by op2[0, n) negacyclically, puts result into
   res[0, n).
*/
void
ref_zn_array_negamul (ulong* res, const ulong* op1, const ulong* op2,
                      size_t n, const zn_mod_t mod)
{
   ulong* temp = (ulong*) malloc (sizeof (ulong) * 2 * n);
   
   ref_zn_array_mul (temp, op1, n, op2, n, mod);
   temp[2 * n - 1] = 0;
   
   mpz_t x;
   mpz_init (x);

   size_t i;
   for (i = 0; i < n; i++)
   {
      mpz_set_ui (x, temp[i]);
      mpz_sub_ui (x, x, temp[i + n]);
      mpz_mod_ui (x, x, mod->m);
      res[i] = mpz_get_ui (x);
   }
   
   mpz_clear (x);
   free (temp);
}



/*
   Reference implementation of scalar multiplication.
   
   Multiplies op[0, n) by x, puts result in res[0, n).
*/
void
ref_zn_array_scalar_mul (ulong* res, const ulong* op, size_t n,
                         ulong x, const zn_mod_t mod)
{
   mpz_t y;
   mpz_init (y);
   
   size_t i;
   for (i = 0; i < n; i++)
   {
      mpz_set_ui (y, op[i]);
      mpz_mul_ui (y, y, x);
      mpz_mod_ui (y, y, mod->m);
      res[i] = mpz_get_ui (y);
   }
   
   mpz_clear (y);
}



/*
   Reference implementation of mpn_smp.
   
   Computes SMP(op1[0, n1), op2[0, n2)), stores result at res[0, n1 - n2 + 3).
*/
void
ref_mpn_smp (mp_limb_t* res,
             const mp_limb_t* op1, size_t n1,
             const mp_limb_t* op2, size_t n2)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   mp_limb_t* prod = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (n1 + n2));
   
   // first compute the ordinary product
   mpn_mul (prod, op1, n1, op2, n2);

   // now we want to remove the cross-terms that could possibly interfere
   // with the result we want, i.e. in the following diagram, we want only
   // contributions from O's, but mpn_mul has given us all of O, A and X,
   // and we will remove the A's.
   
   // OOOOAAXX
   // AOOOOAAX
   // AAOOOOAA
   // XAAOOOOA
   // XXAAOOOO
   
   int which;      // 0 == bottom-left corner, 1 == top-right corner
   size_t diag;    // 0 == closest to diagonal, 1 == next diagonal
   size_t i, x, y, off;
   mp_limb_t lo, hi;
   
   for (which = 0; which <= 1; which++)
   for (diag = 0; diag < ZNP_MIN (n2 - 1, 2); diag++)
   for (i = 0; i < n2 - 1 - diag; i++)
   {
      x = n2 - 2 - i - diag;
      y = i;
      if (which)
      {
         x = n1 - 1 - x;
         y = n2 - 1 - y;
      }
      off = x + y;
      hi = mpn_mul_1 (&lo, op1 + x, 1, op2[y]);
      mpn_sub_1 (prod + off, prod + off, n1 + n2 - off, lo);
      mpn_sub_1 (prod + off + 1, prod + off + 1, n1 + n2 - off - 1, hi);
   }

   // copy the result to the output array
   memcpy (res, prod + n2 - 1, sizeof (mp_limb_t) * (n1 - n2 + 2));
   res[n1 - n2 + 2] = (n2 > 1) ? prod[n1 + 1] : 0;
   free (prod);
}


/*
   Reference implementation of mpn_mulmid.

   Let P = product op1 * op2. Computes P[n2 + 1, n1), stores result at
   res[2, n1 - n2 + 1).
*/
void
ref_mpn_mulmid (mp_limb_t* res, const mp_limb_t* op1, size_t n1,
                const mp_limb_t* op2, size_t n2)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   mp_limb_t* prod = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (n1 + n2));
   
   // compute the ordinary product
   mpn_mul (prod, op1, n1, op2, n2);
   
   // copy relevant segment to output
   if (n1 > n2)
      memcpy (res + 2, prod + n2 + 1, sizeof (mp_limb_t) * (n1 - n2 - 1));
   
   free (prod);
}


// end of file ****************************************************************
