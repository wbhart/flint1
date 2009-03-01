/*
   midmul-profile.c:  routines for profiling middle products
   
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


#include <math.h>
#include "support.h"
#include "profiler.h"
#include "zn_poly_internal.h"


/*
   Wrapper functions to make zn_array_midmul_fallback and
   zn_array_midmul_fft look as if they have the same signature as
   zn_array_midmul.
*/

void zn_array_midmul_fallback_wrapper(
               ulong* res, const ulong* op1, size_t len1,
               const ulong* op2, size_t len2, const zn_mod_t mod)
{
   zn_array_midmul_fallback(res, op1, len1, op2, len2, 0, mod);
}


void zn_array_midmul_fft_wrapper(
               ulong* res, const ulong* op1, size_t len1,
               const ulong* op2, size_t len2, const zn_mod_t mod)
{
   // call the FFT code with the correct scaling factor
   ulong scale = zn_array_midmul_fft_get_fudge(len1, len2, mod);
   zn_array_midmul_fft(res, op1, len1, op2, len2, scale, mod);
}


double profile_midmul(void* arg, unsigned long count)
{
   profile_midmul_info_struct* info = (profile_midmul_info_struct*) arg;

   zn_mod_t mod;
   zn_mod_init(mod, info->n);
   
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * 2 * info->len);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * info->len);
   ulong* buf3 = (ulong*) malloc(sizeof(ulong) * (info->len + 1));

   // generate random inputs
   size_t i;
   for (i = 0; i < 2 * info->len; i++)
      buf1[i] = random_ulong(info->n);
   for (i = 0; i < info->len; i++)
      buf2[i] = random_ulong(info->n);

   void (*target)(ulong*, const ulong*, size_t, const ulong*,
                  size_t, const zn_mod_t);
   
   switch (info->algo)
   {
      case ALGO_MIDMUL_BEST:      target = zn_array_midmul; break;
      case ALGO_MIDMUL_FALLBACK:  target = zn_array_midmul_fallback_wrapper;
                                  break;
/*
      case ALGO_MIDMUL_KS1:       target = zn_array_midmul_KS1; break;
      case ALGO_MIDMUL_KS2:       target = zn_array_midmul_KS2; break;
      case ALGO_MIDMUL_KS3:       target = zn_array_midmul_KS3; break;
      case ALGO_MIDMUL_KS4:       target = zn_array_midmul_KS4; break;
*/
      case ALGO_MIDMUL_FFT:       target = zn_array_midmul_fft_wrapper; break;
      default: abort();
   }
   
   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      target(buf3, buf1, 2 * info->len, buf2, info->len, mod);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter();

   for (j = 0; j < count; j++)
      target(buf3, buf1, 2 * info->len, buf2, info->len, mod);

   cycle_count_t t1 = get_cycle_counter();
   
   free(buf3);
   free(buf2);
   free(buf1);
   
   zn_mod_clear(mod);

   return cycle_diff(t0, t1);
}


// end of file ****************************************************************
