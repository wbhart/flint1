/*
   mpn_mulmid-profile.c:  routines for profiling mpn middle products
   
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


#include <math.h>
#include "support.h"
#include "profiler.h"
#include "zn_poly_internal.h"
#include "zn_poly.h"


double
profile_mpn_mul (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   size_t n1 = info->n1;
   size_t n2 = info->n2;

   ulong* buf1 = (ulong*) malloc (sizeof (mp_limb_t) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (mp_limb_t) * n2);
   ulong* buf3 = (ulong*) malloc (sizeof (mp_limb_t) * (n1 + n2));

   // generate random inputs
   mpn_random (buf1, n1);
   mpn_random (buf2, n2);

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      ZNP_mpn_mul (buf3, buf1, n1, buf2, n2);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      ZNP_mpn_mul (buf3, buf1, n1, buf2, n2);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf3);
   free (buf2);
   free (buf1);

   return cycle_diff (t0, t1);
}


double
profile_mpn_mulmid_fallback (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   size_t n1 = info->n1;
   size_t n2 = info->n2;

   ulong* buf1 = (ulong*) malloc (sizeof (mp_limb_t) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (mp_limb_t) * n2);
   ulong* buf3 = (ulong*) malloc (sizeof (mp_limb_t) * (n1 - n2 + 3));

   // generate random inputs
   mpn_random (buf1, n1);
   mpn_random (buf2, n2);

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      ZNP_mpn_mulmid_fallback (buf3, buf1, n1, buf2, n2);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      ZNP_mpn_mulmid_fallback (buf3, buf1, n1, buf2, n2);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf3);
   free (buf2);
   free (buf1);

   return cycle_diff (t0, t1);
}


double
profile_mpn_smp (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   size_t n1 = info->n1;
   size_t n2 = info->n2;

   ulong* buf1 = (ulong*) malloc (sizeof (mp_limb_t) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (mp_limb_t) * n2);
   ulong* buf3 = (ulong*) malloc (sizeof (mp_limb_t) * (n1 - n2 + 3));

   // generate random inputs
   mpn_random (buf1, n1);
   mpn_random (buf2, n2);

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      ZNP_mpn_smp (buf3, buf1, n1, buf2, n2);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      ZNP_mpn_smp (buf3, buf1, n1, buf2, n2);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf3);
   free (buf2);
   free (buf1);

   return cycle_diff (t0, t1);
}


double
profile_mpn_smp_basecase (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   size_t n1 = info->n1;
   size_t n2 = info->n2;

   ulong* buf1 = (ulong*) malloc (sizeof (mp_limb_t) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (mp_limb_t) * n2);
   ulong* buf3 = (ulong*) malloc (sizeof (mp_limb_t) * (n1 - n2 + 3));

   // generate random inputs
   mpn_random (buf1, n1);
   mpn_random (buf2, n2);

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      ZNP_mpn_smp_basecase (buf3, buf1, n1, buf2, n2);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      ZNP_mpn_smp_basecase (buf3, buf1, n1, buf2, n2);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf3);
   free (buf2);
   free (buf1);

   return cycle_diff (t0, t1);
}


double
profile_mpn_smp_kara (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   size_t n = info->n;

   ulong* buf1 = (ulong*) malloc (sizeof (mp_limb_t) * (2 * n - 1));
   ulong* buf2 = (ulong*) malloc (sizeof (mp_limb_t) * n);
   ulong* buf3 = (ulong*) malloc (sizeof (mp_limb_t) * (n + 2));

   // generate random inputs
   mpn_random (buf1, 2 * n - 1);
   mpn_random (buf2, n);

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      ZNP_mpn_smp_kara (buf3, buf1, buf2, n);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      ZNP_mpn_smp_kara (buf3, buf1, buf2, n);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf3);
   free (buf2);
   free (buf1);

   return cycle_diff (t0, t1);
}



// end of file ****************************************************************
