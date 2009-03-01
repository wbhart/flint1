/*
   nussbaumer-test.c:  test code for functions in nussbaumer.c
   
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
   Tests nussbaumer_mul, for given lgL and modulus.
   Returns 1 on success.
*/
int testcase_nussbaumer_mul(unsigned lgL, const zn_mod_t mod)
{
   ulong len = 1UL << lgL;

   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * len);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * len);
   ulong* compare = (ulong*) malloc(sizeof(ulong) * len);

   // generate random polys
   ulong i;
   for (i = 0; i < len; i++)
      buf1[i] = random_ulong(mod->n);
   for (i = 0; i < len; i++)
      buf2[i] = random_ulong(mod->n);

   // allocate scratch space for nussbaumer_mul
   zn_pmf_vec_t vec1, vec2;
   zn_pmf_vec_init_nussbaumer(vec1, lgL, mod);
   zn_pmf_vec_init_nussbaumer(vec2, lgL, mod);
      
   // compare target implementation against reference implementation
   ref_zn_array_negamul(correct, buf1, buf2, len, mod);
   nussbaumer_mul(compare, buf1, buf2, vec1, vec2);
   ulong scale = nussbaumer_mul_get_fudge(lgL, 0, mod);
   ref_zn_array_scalar_mul(compare, compare, len, scale, mod);
   int success = !zn_array_cmp(correct, compare, len);
   
   zn_pmf_vec_clear(vec2);
   zn_pmf_vec_clear(vec1);
   
   free(compare);
   free(correct);
   free(buf2);
   free(buf1);
   
   return success;
}


/*
   Tests nussbaumer_mul for squaring, for given lgL and modulus.
   Returns 1 on success.
*/
int testcase_nussbaumer_sqr(unsigned lgL, const zn_mod_t mod)
{
   ulong len = 1UL << lgL;

   ulong* buf = (ulong*) malloc(sizeof(ulong) * len);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * len);
   ulong* compare = (ulong*) malloc(sizeof(ulong) * len);

   // generate random poly
   ulong i;
   for (i = 0; i < len; i++)
      buf[i] = random_ulong(mod->n);

   // allocate scratch space for nussbaumer_mul
   zn_pmf_vec_t vec;
   zn_pmf_vec_init_nussbaumer(vec, lgL, mod);
      
   // compare target implementation against reference implementation
   ref_zn_array_negamul(correct, buf, buf, len, mod);
   nussbaumer_mul(compare, buf, buf, vec, NULL);
   ulong scale = nussbaumer_mul_get_fudge(lgL, 1, mod);
   ref_zn_array_scalar_mul(compare, compare, len, scale, mod);
   int success = !zn_array_cmp(correct, compare, len);
   
   zn_pmf_vec_clear(vec);
   
   free(compare);
   free(correct);
   free(buf);
   
   return success;
}


/*
   tests nussbaumer_mul() on a range of input cases (multiplication and
   squaring)
*/
int test_nussbaumer_mul()
{
   int success = 1;
   int i, trial;
   unsigned lgL;

   for (i = 0; i < num_test_bitsizes; i++)
   for (lgL = 2; lgL <= 13 && success; lgL++)
   for (trial = 0; trial < 5 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(test_bitsizes[i], 1));
      
      success = success && testcase_nussbaumer_mul(lgL, mod);
      success = success && testcase_nussbaumer_sqr(lgL, mod);
      
      zn_mod_clear(mod);
   }
   
   return success;
}


// end of file ****************************************************************
