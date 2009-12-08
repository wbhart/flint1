/*
   nuss-test.c:  test code for functions in nussbaumer.c
   
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
   Tests nuss_mul, for given lgL and modulus.
   sqr == 1 means to test squaring.
   Returns 1 on success.
*/
int
testcase_nuss_mul (unsigned lgL, int sqr, const zn_mod_t mod)
{
   ulong n = 1UL << lgL;

   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n);
   ulong* buf2 = sqr ? buf1 : (ulong*) malloc (sizeof (ulong) * n);
   ulong* ref = (ulong*) malloc (sizeof (ulong) * n);
   ulong* res = (ulong*) malloc (sizeof (ulong) * n);

   // generate random polys
   ulong i;
   for (i = 0; i < n; i++)
      buf1[i] = random_ulong (mod->m);
   if (!sqr)
      for (i = 0; i < n; i++)
         buf2[i] = random_ulong (mod->m);

   // allocate scratch space for nuss_mul
   pmfvec_t vec1, vec2;
   pmfvec_init_nuss (vec1, lgL, mod);
   pmfvec_init_nuss (vec2, lgL, mod);
      
   // compare target implementation against reference implementation
   ref_zn_array_negamul (ref, buf1, buf2, n, mod);
   nuss_mul (res, buf1, buf2, vec1, vec2);
   ulong x = nuss_mul_fudge (lgL, 0, mod);
   ref_zn_array_scalar_mul (res, res, n, x, mod);
   int success = !zn_array_cmp (ref, res, n);
   
   pmfvec_clear (vec2);
   pmfvec_clear (vec1);
   
   free (res);
   free (ref);
   if (!sqr)
      free (buf2);
   free (buf1);
   
   return success;
}


/*
   tests nuss_mul() on a range of input cases (multiplication and squaring)
*/
int
test_nuss_mul (int quick)
{
   int success = 1;
   int i, trial;
   unsigned lgL;
   zn_mod_t mod;

   for (i = 0; i < num_test_bitsizes; i++)
   for (lgL = 2; lgL <= (quick ? 11 : 13) && success; lgL++)
   for (trial = 0; trial < (quick ? 1 : 5) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
      success = success && testcase_nuss_mul (lgL, 0, mod);
      success = success && testcase_nuss_mul (lgL, 1, mod);
      zn_mod_clear (mod);
   }
   
   return success;
}


// end of file ****************************************************************
