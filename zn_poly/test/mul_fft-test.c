/*
   mul_fft-test.c:  test code for functions in mul_fft.c and mul_fft_dft.c
   
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
   Tests zn_array_mul_fft, for given lengths and modulus.
   
   If use_scale is set, zn_array_mul_fft() gets called with a random x
   (post-scaling factor), otherwise gets called with x == 1.
   
   If sqr == 1, tests squaring (n2 is ignored), otherwise ordinary
   multiplication.
   
   Returns 1 on success.
*/
int
testcase_zn_array_mul_fft (size_t n1, size_t n2, int sqr, int use_scale,
                           const zn_mod_t mod)
{
   if (sqr)
      n2 = n1;

   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n1);
   ulong* buf2 = sqr ? buf1 : (ulong*) malloc (sizeof (ulong) * n2);
   ulong* ref = (ulong*) malloc (sizeof (ulong) * (n1 + n2 - 1));
   ulong* res = (ulong*) malloc (sizeof (ulong) * (n1 + n2 - 1));

   // generate random polys
   size_t i;
   for (i = 0; i < n1; i++)
      buf1[i] = random_ulong (mod->m);
   if (!sqr)
      for (i = 0; i < n2; i++)
         buf2[i] = random_ulong (mod->m);

   ulong x = use_scale ? random_ulong (mod->m) : 1;

   // compare target implementation against reference implementation
   ref_zn_array_mul (ref, buf1, n1, buf2, n2, mod);
   ref_zn_array_scalar_mul (ref, ref, n1 + n2 - 1, x, mod);
      
   zn_array_mul_fft (res, buf1, n1, buf2, n2, x, mod);
   ulong y = zn_array_mul_fft_fudge (n1, n2, sqr, mod);
   ref_zn_array_scalar_mul (res, res, n1 + n2 - 1, y, mod);
   
   int success = !zn_array_cmp (ref, res, n1 + n2 - 1);
   
   free (res);
   free (ref);
   if (!sqr)
      free (buf2);
   free (buf1);
   
   return success;
}



/*
   tests zn_array_mul_fft() on a range of input cases
*/
int
test_zn_array_mul_or_sqr_fft (int sqr, int quick)
{
   int success = 1;
   int i, trial, use_scale;
   size_t n1, n2;
   zn_mod_t mod;

   // first try a dense range of "small" problems

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (n2 = 1; n2 <= 50 && success; n2 += (quick ? 3 : 1))
   for (n1 = n2; n1 <= 50 && (!sqr || n1 <= n2) && success;
        n1 += (quick ? 3 : 1))
   for (use_scale = 0; use_scale <= 1 && success; use_scale++)
   for (trial = 0; trial < (quick ? 1 : 3) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
      success = success &&
                     testcase_zn_array_mul_fft (n1, n2, sqr, use_scale, mod);
      zn_mod_clear (mod);
   }
   
   // now try some random larger problems
   // and temporarily change the nussbaumer thresholds so we use that
   // code sometimes
   unsigned thresh;

   for (i = 0; i < num_test_bitsizes && success; i++)
   {
      unsigned b = test_bitsizes[i];
      unsigned* c = sqr ? &(tuning_info[b].nuss_sqr_thresh)
                        : &(tuning_info[b].nuss_mul_thresh);

      for (use_scale = 0; use_scale <= 1 && success; use_scale++)
      for (thresh = 2; thresh <= 8 && success; thresh += (quick ? 4 : 1))
      {
         unsigned save_thresh = *c;
         *c = thresh;
      
         size_t t1 = random_ulong (quick ? 3000 : 10000) + 1;
         size_t t2 = sqr ? t1 : (random_ulong (quick ? 3000 : 10000) + 1);
         n1 = ZNP_MAX (t1, t2);
         n2 = ZNP_MIN (t1, t2);
         
         zn_mod_init (mod, random_modulus (b, 1));
         success = success &&
                      testcase_zn_array_mul_fft (n1, n2, sqr, use_scale, mod);
         zn_mod_clear (mod);

         *c = save_thresh;
      }
   }

   return success;
}


int
test_zn_array_mul_fft (int quick)
{
   return test_zn_array_mul_or_sqr_fft (0, quick);
}


int
test_zn_array_sqr_fft (int quick)
{
   return test_zn_array_mul_or_sqr_fft (1, quick);
}



/*
   Tests zn_array_mul_fft_dft, for given lengths, lgT, modulus.
   Returns 1 on success.
*/
int
testcase_zn_array_mul_fft_dft (size_t n1, size_t n2, unsigned lgT,
                               const zn_mod_t mod)
{
   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (ulong) * n2);
   ulong* ref = (ulong*) malloc (sizeof (ulong) * (n1 + n2 - 1));
   ulong* res = (ulong*) malloc (sizeof (ulong) * (n1 + n2 - 1));

   // generate random polys
   size_t i;
   for (i = 0; i < n1; i++)
      buf1[i] = random_ulong (mod->m);
   for (i = 0; i < n2; i++)
      buf2[i] = random_ulong (mod->m);
      
   // compare target implementation against reference implementation
   ref_zn_array_mul (ref, buf1, n1, buf2, n2, mod);
   zn_array_mul_fft_dft (res, buf1, n1, buf2, n2, lgT, mod);
   int success = !zn_array_cmp (ref, res, n1 + n2 - 1);
   
   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}



/*
   Tests zn_array_mulmid_fft, for given n1, n2, modulus.

   If use_scale is set, zn_array_mulmid_fft() gets called with a random x
   (post-scaling factor), otherwise gets called with x == 1.

   Returns 1 on success.
*/
int
testcase_zn_array_mulmid_fft (size_t n1, size_t n2, int use_scale,
                              const zn_mod_t mod)
{
   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (ulong) * n2);
   ulong* ref = (ulong*) malloc (sizeof (ulong) * (n1 - n2 + 1));
   ulong* res = (ulong*) malloc (sizeof (ulong) * (n1 - n2 + 1));

   // generate random polys
   size_t i;
   for (i = 0; i < n1; i++)
      buf1[i] = random_ulong (mod->m);
   for (i = 0; i < n2; i++)
      buf2[i] = random_ulong (mod->m);
      
   ulong x = use_scale ? random_ulong (mod->m) : 1;

   // compare target implementation against reference implementation
   ref_zn_array_mulmid (ref, buf1, n1, buf2, n2, mod);
   ref_zn_array_scalar_mul (ref, ref, n1 - n2 + 1, x, mod);

   zn_array_mulmid_fft (res, buf1, n1, buf2, n2, x, mod);
   ulong y = zn_array_mulmid_fft_fudge (n1, n2, mod);
   ref_zn_array_scalar_mul (res, res, n1 - n2 + 1, y, mod);
   
   int success = !zn_array_cmp (ref, res, n1 - n2 + 1);

   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}


/*
   tests zn_array_mulmid_fft() on a range of input cases
*/
int
test_zn_array_mulmid_fft (int quick)
{
   int success = 1;
   int i, trial, use_scale;
   size_t n1, n2;
   zn_mod_t mod;

   // first try a dense range of "small" problems

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (n2 = 1; n2 <= 50 && success; n2 += (quick ? 3 : 1))
   for (n1 = n2; n1 <= 50 && success; n1 += (quick ? 3 : 1))
   for (use_scale = 0; use_scale <= 1 && success; use_scale++)
   for (trial = 0; trial < (quick ? 1 : 3) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
      success = success &&
                   testcase_zn_array_mulmid_fft (n1, n2, use_scale, mod);
      zn_mod_clear (mod);
   }

   // now try some random larger problems
   // and temporarily change the nussbaumer thresholds so we use that
   // code sometimes
   
   ulong thresh;
   
   for (i = 0; i < num_test_bitsizes && success; i++)
   for (use_scale = 0; use_scale <= 1 && success; use_scale++)
   for (thresh = 2; thresh <= 8; thresh += (quick ? 4 : 1))
   {
      unsigned b = test_bitsizes[i];

      ulong save_thresh = tuning_info[b].nuss_mul_thresh;
      tuning_info[b].nuss_mul_thresh = thresh;
   
      size_t t1 = random_ulong (quick ? 3000 : 10000) + 1;
      size_t t2 = random_ulong (quick ? 3000 : 10000) + 1;
      n1 = ZNP_MAX (t1, t2);
      n2 = ZNP_MIN (t1, t2);
      
      zn_mod_init (mod, random_modulus (b, 1));
      success = success &&
                   testcase_zn_array_mulmid_fft (n1, n2, use_scale, mod);
      zn_mod_clear (mod);

      tuning_info[b].nuss_mul_thresh = save_thresh;
   }

   return success;
}



/*
   tests zn_array_mul_fft_dft() on a range of input cases
*/
int
test_zn_array_mul_fft_dft (int quick)
{
   int success = 1;
   int i, trial;
   unsigned lgT;
   size_t n1, n2;
   zn_mod_t mod;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (n2 = 1; n2 <= 30 && success; n2 += (quick ? random_ulong (2) + 1 : 1))
   for (n1 = n2; n1 <= 30 && success; n1 += (quick ? random_ulong (2) + 1 : 1))
   for (lgT = 0; lgT < 5 && success; lgT++)
   for (trial = 0; trial < (quick ? 1 : 3) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
      success = success && testcase_zn_array_mul_fft_dft (n1, n2, lgT, mod);
      zn_mod_clear (mod);
   }
   
   return success;
}


// end of file ****************************************************************
