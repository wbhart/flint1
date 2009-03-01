/*
   invert-profile.c:  routines for profiling power series inversion
   
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


double profile_invert(void* arg, unsigned long count)
{
   profile_invert_info_struct* info = (profile_invert_info_struct*) arg;
   
   if (info->algo == ALGO_INVERT_NTL)
      return profile_invert_ntl(arg, count);

   zn_mod_t mod;
   zn_mod_init(mod, info->n);
   
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * info->len);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * info->len);

   // generate random inputs
   size_t i;
   for (i = 0; i < info->len; i++)
      buf1[i] = random_ulong(info->n);
   // todo: change this next line when we can handle non-monic input
   buf1[0] = 1;

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      zn_array_invert(buf2, buf1, info->len, mod);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter();

   for (j = 0; j < count; j++)
      zn_array_invert(buf2, buf1, info->len, mod);

   cycle_count_t t1 = get_cycle_counter();
   
   free(buf2);
   free(buf1);
   
   zn_mod_clear(mod);

   return cycle_diff(t0, t1);
}


// end of file ****************************************************************
