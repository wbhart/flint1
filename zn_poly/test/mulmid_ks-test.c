/*
   mulmid_ks-test.c:  test code for functions in mulmid_ks.c
   
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
   Tests zn_array_mulmid_KSk, for given lengths, reduction algorithm, modulus.
   1 <= k <= 4 indicates which KS variant to call.
   Returns 1 on success.
*/
int
testcase_zn_array_mulmid_KS (int k, size_t n1, size_t n2,
                             int redc, const zn_mod_t mod)
{
   // disallow REDC if modulus is even
   if (!(mod->m & 1))
      redc = 0;

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
      
   // compare target implementation against reference implementation
   ref_zn_array_mulmid (ref, buf1, n1, buf2, n2, mod);

   switch (k)
   {
      case 1:  zn_array_mulmid_KS1 (res, buf1, n1, buf2, n2, redc, mod); break;
      case 2:  zn_array_mulmid_KS2 (res, buf1, n1, buf2, n2, redc, mod); break;
      case 3:  zn_array_mulmid_KS3 (res, buf1, n1, buf2, n2, redc, mod); break;
      case 4:  zn_array_mulmid_KS4 (res, buf1, n1, buf2, n2, redc, mod); break;
      default:
         printf ("oops!\n"); abort ();
   }

   if (redc)
      // correct for REDC reduction
      ref_zn_array_scalar_mul (res, res, n1 - n2 + 1, mod->m - mod->B, mod);

   int success = !zn_array_cmp (ref, res, n1 - n2 + 1);

   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}


/*
   tests zn_array_mulmid_KSk() on a range of input cases,
   where 1 <= k <= 4
*/
int
test_zn_array_mulmid_KSk (unsigned k, int quick)
{
   int success = 1;
   int b, trial, redc;
   size_t n1, n2;
   zn_mod_t mod;

   // first try a dense range of "small" problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (n2 = 1; n2 <= 30 && success; n2 += (quick ? 5 : 1))
   for (n1 = n2; n1 <= 30 && success; n1 += (quick ? 5 : 1))
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 1 : 10) && success; trial++)
   {
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_mulmid_KS (k, n1, n2, redc, mod);
      zn_mod_clear(mod);
   }
   
   // now try some random larger problems

   for (b = 2; b <= ULONG_BITS && success; b++)
   for (redc = 0; redc < 2 && success; redc++)
   for (trial = 0; trial < (quick ? 3 : 200) && success; trial++)
   {
      size_t t1 = random_ulong (quick ? 250 : 1000) + 1;
      size_t t2 = random_ulong (quick ? 250 : 1000) + 1;
      n1 = ZNP_MAX (t1, t2);
      n2 = ZNP_MIN (t1, t2);
      
      zn_mod_init (mod, random_modulus (b, 0));
      success = success && testcase_zn_array_mulmid_KS (k, n1, n2, redc, mod);
      zn_mod_clear(mod);
   }

   return success;
}


int test_zn_array_mulmid_KS1 (int quick)
{
   return test_zn_array_mulmid_KSk (1, quick);
}

int test_zn_array_mulmid_KS2 (int quick)
{
   return test_zn_array_mulmid_KSk (2, quick);
}

int test_zn_array_mulmid_KS3 (int quick)
{
   return test_zn_array_mulmid_KSk (3, quick);
}

int test_zn_array_mulmid_KS4 (int quick)
{
   return test_zn_array_mulmid_KSk (4, quick);
}


// end of file ****************************************************************
