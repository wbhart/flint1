/*
   invert-test.c:  test code for functions in invert.c
   
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
   Tests zn_array_invert() for a given series length and modulus.
*/
int testcase_zn_array_invert(size_t len, const zn_mod_t mod)
{
   ulong* op = (ulong*) malloc(sizeof(ulong) * len);
   ulong* res = (ulong*) malloc(sizeof(ulong) * len);
   ulong* check = (ulong*) malloc(sizeof(ulong) * (2 * len - 1));
   
   // make up random input poly
   size_t i;
   op[0] = 1;    // todo: generalise this to any invertible element
   for (i = 1; i < len; i++)
      op[i] = random_ulong(mod->n);
      
   // compute inverse
   zn_array_invert(res, op, len, mod);

   // multiply by original series and check we get 1
   ref_zn_array_mul(check, op, len, res, len, mod);
   
   int success = (check[0] == 1);
   for (i = 1; i < len; i++)
      success = success && (check[i] == 0);

   free(check);
   free(res);
   free(op);
   
   return success;
}


/*
   Tests zn_array_invert() on a range of problems.
*/
int test_zn_array_invert()
{
   int success = 1;
   int bits, trial;
   size_t len;

   // first try a dense range of "small" problems

   for (bits = 2; bits <= ULONG_BITS && success; bits++)
   for (len = 1; len <= 100 && success; len++)
   for (trial = 0; trial < 10 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 0));
      
      success = success && testcase_zn_array_invert(len, mod);
      
      zn_mod_clear(mod);
   }
   
   // now try a few larger random problems

   for (bits = 2; bits <= ULONG_BITS && success; bits++)
   for (trial = 0; trial < 5 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 0));
      
      len = random_ulong(10000) + 1;
      success = success && testcase_zn_array_invert(len, mod);
      
      zn_mod_clear(mod);
   }

   return success;
}


// end of file ****************************************************************
