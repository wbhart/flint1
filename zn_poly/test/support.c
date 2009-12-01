/*
   support.c:  various support routines for test, profiling and tuning code
   
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


unsigned test_bitsizes[] = {2, 3, 8,
                            ULONG_BITS/2 - 1, ULONG_BITS/2, ULONG_BITS/2 + 1,
                            ULONG_BITS - 2, ULONG_BITS - 1, ULONG_BITS};
                  
unsigned num_test_bitsizes = sizeof (test_bitsizes) / sizeof (unsigned);


gmp_randstate_t randstate;


void
mpz_to_mpn (mp_limb_t* res, size_t n, const mpz_t op)
{
   ZNP_ASSERT (mpz_size (op) <= n);
   size_t count_p;
   
   mpz_export (res, &count_p, -1, sizeof (mp_limb_t), 0, GMP_NAIL_BITS, op);
   // zero-pad remaining buffer
   for (; count_p < n; count_p++)
      res[count_p] = 0;
}


void
mpn_to_mpz (mpz_t res, const mp_limb_t* op, size_t n)
{
   ZNP_ASSERT (n >= 1);
   mpz_import (res, n, -1, sizeof (mp_limb_t), 0, GMP_NAIL_BITS, op);
}


ulong
random_ulong (ulong max)
{
   return gmp_urandomm_ui (randstate, max);
}


ulong
random_ulong_bits (unsigned b)
{
   return gmp_urandomb_ui (randstate, b);
}


ulong
random_modulus (unsigned b, int require_odd)
{
   ZNP_ASSERT(b >= 2 && b <= ULONG_BITS);

   if (require_odd)
   {
      if (b == 2)
         return 3;
      return (1UL << (b - 1)) + 2 * random_ulong_bits (b - 2) + 1;
   }
   else
      return (1UL << (b - 1)) + random_ulong_bits (b - 1);
}


void
zn_array_print (const ulong* x, size_t n)
{
   size_t i;
   printf ("[");
   for (i = 0; i < n; i++)
   {
      if (i)
         printf (", ");
      printf ("%lu", x[i]);
   }
   printf ("]");
}


void
ZNP_mpn_random2 (mp_limb_t* res, size_t n)
{
   size_t i;
   
   mpn_random2 (res, n);
   
   if (random_ulong (2))
      for (i = 0; i < n; i++)
         res[i] ^= GMP_NUMB_MASK;
}


// end of file ****************************************************************
