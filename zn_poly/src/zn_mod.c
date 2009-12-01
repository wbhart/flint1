/*
   zn_mod.c:  functions operating on zn_mod_t objects
   
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

#include "zn_poly_internal.h"


void
zn_mod_init (zn_mod_t mod, ulong m)
{
   ZNP_ASSERT (m >= 2);
   
   mod->m = m;
   mod->bits = ceil_lg (m);
   
   mpz_t x, y;
   mpz_init (x);
   mpz_init (y);

   // compute B and B^2 mod m
   mpz_set_ui (x, 1);
   mpz_mul_2exp (x, x, ULONG_BITS);
   mpz_mod_ui (x, x, m);
   mod->B = mpz_get_ui (x);
   
   mpz_set_ui (x, 1);
   mpz_mul_2exp (x, x, 2*ULONG_BITS);
   mpz_mod_ui (x, x, m);
   mod->B2 = mpz_get_ui (x);
   
   // compute sh1 and inv1
   mod->sh1 = ceil_lg (m) - 1;
   mpz_set_ui (x, 1);
   mpz_mul_2exp (x, x, mod->sh1 + 1);
   mpz_sub_ui (x, x, m);
   mpz_mul_2exp (x, x, ULONG_BITS);
   mpz_fdiv_q_ui (x, x, m);
   mpz_add_ui (x, x, 1);
   mod->inv1 = mpz_get_ui (x);

   // compute sh2, sh3, inv2, m_norm
   unsigned ell = floor_lg (m) + 1;
   mod->sh2 = ULONG_BITS - ell;
   mod->sh3 = ell - 1;
   mod->m_norm = m << mod->sh2;
   mpz_set_ui (x, 1);
   mpz_mul_2exp (x, x, ell);
   mpz_sub_ui (x, x, m);
   mpz_mul_2exp (x, x, ULONG_BITS);
   mpz_sub_ui (x, x, 1);
   mpz_fdiv_q_ui (x, x, m);
   mod->inv2 = mpz_get_ui (x);
   
   // compute inv3, if m is odd
   if (m & 1)
   {
      // m^(-1) = m mod 8
      ulong minv = m;
      // lift 2-adically
      int i;
      for (i = 3; i < ULONG_BITS; i <<= 1)
         minv = 2 * minv -  m * minv * minv;
      mod->inv3 = minv;
   }

   mpz_clear (y);
   mpz_clear (x);
}


void
zn_mod_clear (zn_mod_t mod)
{
   // nothing to do yet, but maybe one day there will be
}



ulong
zn_mod_pow2 (int k, const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (k > -ULONG_BITS && k < ULONG_BITS);
   
   if (k == 0)
      return 1;
   
   if (k > 0)
      return zn_mod_reduce (1UL << k, mod);
   
   return zn_mod_pow (zn_mod_divby2 (1, mod), -k, mod);
}



ulong
zn_mod_pow (ulong x, long k, const zn_mod_t mod)
{
   ZNP_ASSERT (k >= 0);

   // repeated squaring
   ulong prod = 1;
   ulong x_pow = x;
   for (; k; k >>= 1)
   {
      if (k & 1)
         prod = zn_mod_mul (prod, x_pow, mod);
      x_pow = zn_mod_mul (x_pow, x_pow, mod);
   }
   return prod;
}



ulong
zn_mod_invert (ulong x, const zn_mod_t mod)
{
   ZNP_ASSERT (x < mod->m);

   // for now just use GMP
   mpz_t a, m;
   mpz_init (a);
   mpz_set_ui (a, x);
   mpz_init (m);
   mpz_set_ui (m, mod->m);
   int success = mpz_invert (a, a, m);
   x = mpz_get_ui (a);
   mpz_clear (m);
   mpz_clear (a);

   return success ? x : 0;
}


// end of file ****************************************************************
