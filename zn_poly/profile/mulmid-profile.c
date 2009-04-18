/*
   mulmid-profile.c:  routines for profiling middle products
   
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


/*
   Wrapper functions to make zn_array_mulmid, zn_array_mulmid_fallback and
   zn_array_mulmid_fft look as if they have a redc flag.
*/

void
zn_array_mulmid_wrapper (ulong* res,
                         const ulong* op1, size_t n1,
                         const ulong* op2, size_t n2,
                         int redc, const zn_mod_t mod)
{
   zn_array_mulmid (res, op1, n1, op2, n2, mod);
}


void
zn_array_mulmid_fallback_wrapper (ulong* res,
                                  const ulong* op1, size_t n1,
                                  const ulong* op2, size_t n2,
                                  int redc, const zn_mod_t mod)
{
   zn_array_mulmid_fallback (res, op1, n1, op2, n2, 0, mod);
}


void
zn_array_mulmid_fft_wrapper (ulong* res,
                             const ulong* op1, size_t n1,
                             const ulong* op2, size_t n2,
                             int redc, const zn_mod_t mod)
{
   // call the FFT code with the correct scaling factor
   ulong x = zn_array_mulmid_fft_fudge (n1, n2, mod);
   zn_array_mulmid_fft (res, op1, n1, op2, n2, x, mod);
}


double
profile_mulmid (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   size_t n1 = info->n1;
   size_t n2 = info->n2;

   zn_mod_t mod;
   zn_mod_init (mod, info->m);
   
   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n1);
   ulong* buf2 = (ulong*) malloc (sizeof (ulong) * n2);
   ulong* buf3 = (ulong*) malloc (sizeof (ulong) * (n1 - n2 + 1));

   // generate random inputs
   size_t i;
   for (i = 0; i < n1; i++)
      buf1[i] = random_ulong (info->m);
   for (i = 0; i < n2; i++)
      buf2[i] = random_ulong (info->m);

   void (*target)(ulong*, const ulong*, size_t, const ulong*, size_t,
                  int, const zn_mod_t);

   int redc;
   
   switch (info->algo)
   {
      case ALGO_MULMID_BEST:      target = zn_array_mulmid_wrapper; break;
      case ALGO_MULMID_FALLBACK:  target = zn_array_mulmid_fallback_wrapper;
                                  break;
      case ALGO_MULMID_KS1:       target = zn_array_mulmid_KS1; redc = 0;
                                  break;
      case ALGO_MULMID_KS1_REDC:  target = zn_array_mulmid_KS1; redc = 1;
                                  break;
      case ALGO_MULMID_KS2:       target = zn_array_mulmid_KS2; redc = 0;
                                  break;
      case ALGO_MULMID_KS2_REDC:  target = zn_array_mulmid_KS2; redc = 1;
                                  break;
      case ALGO_MULMID_KS3:       target = zn_array_mulmid_KS3; redc = 0;
                                  break;
      case ALGO_MULMID_KS3_REDC:  target = zn_array_mulmid_KS3; redc = 1;
                                  break;
      case ALGO_MULMID_KS4:       target = zn_array_mulmid_KS4; redc = 0;
                                  break;
      case ALGO_MULMID_KS4_REDC:  target = zn_array_mulmid_KS4; redc = 1;
                                  break;
      case ALGO_MULMID_FFT:       target = zn_array_mulmid_fft_wrapper; break;
      default: abort ();
   }
   
   // warm up
   ulong j;
   for (j = 0; j < count/4; j++)
      target (buf3, buf1, n1, buf2, n2, redc, mod);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      target (buf3, buf1, n1, buf2, n2, redc, mod);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf3);
   free (buf2);
   free (buf1);
   
   zn_mod_clear (mod);

   return cycle_diff (t0, t1);
}


// end of file ****************************************************************
