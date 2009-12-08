/*
   invert-test.c:  test code for functions in invert.c
   
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
   Tests zn_array_invert() for a given series length and modulus.
*/
int
testcase_zn_array_invert (size_t n, const zn_mod_t mod)
{
   ulong* op = (ulong*) malloc (sizeof (ulong) * n);
   ulong* res = (ulong*) malloc (sizeof (ulong) * n);
   ulong* check = (ulong*) malloc (sizeof (ulong) * (2 * n - 1));
   
   // make up random input poly
   size_t i;
   op[0] = 1;
   for (i = 1; i < n; i++)
      op[i] = random_ulong (mod->m);
      
   // compute inverse
   zn_array_invert (res, op, n, mod);

   // multiply by original series and check we get 1
   ref_zn_array_mul (check, op, n, res, n, mod);
   
   int success = (check[0] == 1);
   for (i = 1; i < n; i++)
      success = success && (check[i] == 0);

   free (check);
   free (res);
   free (op);
   
   return success;
}


/*
   Tests zn_array_invert() on a range of problems.
*/
int
test_zn_array_invert (int quick)
{
   int success = 1;
   int b, trial;
   size_t n;
   zn_mod_t mod;

   // first try a dense range of "small" problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (n = 1; n <= 60 && success; n++)
   for (trial = 0; trial < (quick ? 1 : 10) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_invert (n, mod);
      zn_mod_clear (mod);
   }
   
   // now try a few larger random problems

   for (b = 2; b <= ULONG_BITS && success;
        b += (quick ? random_ulong (3) + 1 : 1))
   for (trial = 0; trial < (quick ? 1 : 5) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (b, 0));
      n = random_ulong (quick ? 2000 : 10000) + 1;
      success = success && testcase_zn_array_invert (n, mod);
      zn_mod_clear (mod);
   }

   return success;
}


// end of file ****************************************************************
