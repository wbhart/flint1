/*
   mpn_mulmid-test.c:  test code for functions in mpn_mulmid.c
   
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
   Tests mpn_smp for given n1 >= n2 >= 1.
*/
int 
testcase_mpn_smp_basecase (size_t n1, size_t n2)
{
   size_t n3 = n1 - n2 + 3;
   mp_limb_t* buf1 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n1);
   mp_limb_t* buf2 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n2);
   mp_limb_t* ref = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n3);
   mp_limb_t* res = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n3);

   // generate random inputs
   ZNP_mpn_random2 (buf1, n1);
   ZNP_mpn_random2 (buf2, n2);
   
   // compare target against reference implementation
   ZNP_mpn_smp_basecase (res, buf1, n1, buf2, n2);
   ref_mpn_smp (ref, buf1, n1, buf2, n2);
   int success = !mpn_cmp (ref, res, n3);
   
   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}


/*
   Tests mpn_smp for a range of n1, n2.
*/
int
test_mpn_smp_basecase (int quick)
{
   int success = 1;
   size_t n1, n2;
   ulong trial;

   for (n2 = 1; n2 <= 30 && success; n2++)
   for (n1 = n2; n1 <= 30 && success; n1++)
   for (trial = 0; trial < (quick ? 30 : 3000) && success; trial++)
      success = success && testcase_mpn_smp_basecase (n1, n2);

   return success;
}


/*
   Tests mpn_smp_kara for given n.
*/
int
testcase_mpn_smp_kara (size_t n)
{
   mp_limb_t* buf1 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (2 * n - 1));
   mp_limb_t* buf2 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n);
   mp_limb_t* ref = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (n + 2));
   mp_limb_t* res = (mp_limb_t*) malloc (sizeof (mp_limb_t) * (n + 2));

   // generate random inputs
   ZNP_mpn_random2 (buf1, 2 * n - 1);
   ZNP_mpn_random2 (buf2, n);
   
   // compare target against reference implementation
   ZNP_mpn_smp_kara (res, buf1, buf2, n);
   ref_mpn_smp (ref, buf1, 2 * n - 1, buf2, n);
   int success = !zn_array_cmp (ref, res, n + 2);
   
   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}


/*
   Tests mpn_smp_kara for a range of n.
*/
int
test_mpn_smp_kara (int quick)
{
   int success = 1;
   size_t n;
   ulong trial;

   // first a dense range of small problems
   for (n = 2; n <= 30 && success; n++)
   for (trial = 0; trial < (quick ? 300 : 30000) && success; trial++)
      success = success && testcase_mpn_smp_kara (n);
   
   // now a few larger problems too
   for (trial = 0; trial < (quick ? 100 : 3000) && success; trial++)
   {
      n = random_ulong (3 * ZNP_mpn_smp_kara_thresh) + 2;
      success = success && testcase_mpn_smp_kara (n);
   }

   return success;
}


int
testcase_mpn_smp (size_t n1, size_t n2)
{
   size_t n3 = n1 - n2 + 3;
   mp_limb_t* buf1 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n1);
   mp_limb_t* buf2 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n2);
   mp_limb_t* ref = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n3);
   mp_limb_t* res = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n3);

   // generate random inputs
   ZNP_mpn_random2 (buf1, n1);
   ZNP_mpn_random2 (buf2, n2);
   
   // temporarily lower the karatsuba threshold for more stringent testing
   unsigned long temp_thresh = ZNP_mpn_smp_kara_thresh;
   ZNP_mpn_smp_kara_thresh = 5;
   
   // compare target against reference implementation
   ZNP_mpn_smp (res, buf1, n1, buf2, n2);
   ref_mpn_smp (ref, buf1, n1, buf2, n2);
   int success = !zn_array_cmp (ref, res, n3);

   ZNP_mpn_smp_kara_thresh = temp_thresh;
   
   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}


/*
   Tests mpn_smp for a range of n1, n2.
*/
int
test_mpn_smp (int quick)
{
   int success = 1;
   size_t n1, n2;
   ulong trial;

   for (n2 = 1; n2 <= 30 && success; n2++)
   for (n1 = n2; n1 <= 30 && success; n1++)
   for (trial = 0; trial < (quick ? 30 : 3000) && success; trial++)
      success = success && testcase_mpn_smp (n1, n2);

   return success;
}



int
testcase_mpn_mulmid (size_t n1, size_t n2)
{
   size_t n3 = n1 - n2 + 3;
   mp_limb_t* buf1 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n1);
   mp_limb_t* buf2 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n2);
   mp_limb_t* ref = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n3);
   mp_limb_t* res = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n3);

   // generate random inputs
   ZNP_mpn_random2 (buf1, n1);
   ZNP_mpn_random2 (buf2, n2);
   
   // compare target against reference implementation
   ZNP_mpn_mulmid (res, buf1, n1, buf2, n2);
   ref_mpn_mulmid (ref, buf1, n1, buf2, n2);
   int success = (n3 <= 4) || !mpn_cmp (ref + 2, res + 2, n3 - 4);

   free (res);
   free (ref);
   free (buf2);
   free (buf1);
   
   return success;
}



/*
   Tests mpn_mulmid for a range of n1, n2.
*/
int
test_mpn_mulmid (int quick)
{
   int success = 1;
   size_t n1, n2;
   ulong trial;

   for (n2 = 1; n2 <= 30 && success; n2++)
   for (n1 = n2; n1 <= 30 && success; n1++)
   for (trial = 0; trial < (quick ? 30 : 3000) && success; trial++)
      success = success && testcase_mpn_mulmid (n1, n2);

   return success;
}


// end of file ****************************************************************
