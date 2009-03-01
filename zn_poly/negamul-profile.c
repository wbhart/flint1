/*
   negamul-profile.c:  routines for profiling negacyclic multiplication
                       and squaring
   
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
#include "zn_poly.h"


double profile_negamul(void* arg, unsigned long count)
{
   profile_negamul_info_struct* info = (profile_negamul_info_struct*) arg;

   ulong j;
   
   zn_mod_t mod;
   zn_mod_init(mod, info->n);
   
   size_t len = 1UL << info->lg_len;
   
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len);
   ulong* buf2 = info->squaring ? buf1
                                : ((ulong*) malloc(sizeof(ulong) * len));
   ulong* buf3 = (ulong*) malloc(sizeof(ulong) * len);

   cycle_count_t t0, t1;

   // generate random inputs
   size_t i;
   for (i = 0; i < len; i++)
      buf1[i] = random_ulong(info->n);
   for (i = 0; i < len; i++)
      buf2[i] = random_ulong(info->n);

   if (info->algo == ALGO_NEGAMUL_FALLBACK)
   {
      // KS version
      ulong* temp = (ulong*) malloc(sizeof(ulong) * 2 * len);
      
      // warm up
      for (j = 0; j < count/4; j++)
      {
         zn_array_mul(temp, buf1, len, buf2, len, mod);
         zn_array_sub(buf3, temp, temp + len, len, mod);
      }

      // do the actual profile
      t0 = get_cycle_counter();

      for (j = 0; j < count; j++)
      {
         zn_array_mul(temp, buf1, len, buf2, len, mod);
         zn_array_sub(buf3, temp, temp + len, len, mod);
      }
      
      t1 = get_cycle_counter();
      
      free(temp);
   }
   else if (info->algo == ALGO_NEGAMUL_NUSSBAUMER)
   {
      // nussbaumer version
      
      zn_pmf_vec_t vec1, vec2;
      zn_pmf_vec_init_nussbaumer(vec1, info->lg_len, mod);
      zn_pmf_vec_init_nussbaumer(vec2, info->lg_len, mod);

      // warm up
      for (j = 0; j < count/4; j++)
         nussbaumer_mul(buf3, buf1, buf2, vec1, vec2);

      // do the actual profile
      t0 = get_cycle_counter();

      for (j = 0; j < count; j++)
         nussbaumer_mul(buf3, buf1, buf2, vec1, vec2);
      
      t1 = get_cycle_counter();
      
      zn_pmf_vec_clear(vec2);
      zn_pmf_vec_clear(vec1);
   }
   else
      abort();

   free(buf3);
   if (!info->squaring)
      free(buf2);
   free(buf1);
   
   zn_mod_clear(mod);
   
   return cycle_diff(t0, t1);
}



// end of file ****************************************************************
