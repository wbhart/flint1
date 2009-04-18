/*
   negamul-profile.c:  routines for profiling negacyclic multiplication
                       and squaring
   
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
profile_negamul (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   ulong j;
   
   zn_mod_t mod;
   zn_mod_init (mod, info->m);
   
   size_t n = 1UL << info->lgL;
   
   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n);
   ulong* buf2 = info->sqr ? buf1 : ((ulong*) malloc (sizeof (ulong) * n));
   ulong* buf3 = (ulong*) malloc (sizeof (ulong) * n);

   cycle_count_t t0, t1;

   // generate random inputs
   size_t i;
   for (i = 0; i < n; i++)
      buf1[i] = random_ulong (info->m);
   for (i = 0; i < n; i++)
      buf2[i] = random_ulong (info->m);

   if (info->algo == ALGO_NEGAMUL_FALLBACK)
   {
      // KS version
      ulong* temp = (ulong*) malloc (sizeof (ulong) * 2 * n);
      
      // warm up
      for (j = 0; j < count/4; j++)
      {
         zn_array_mul (temp, buf1, n, buf2, n, mod);
         zn_array_sub (buf3, temp, temp + n, n, mod);
      }

      // do the actual profile
      t0 = get_cycle_counter ();

      for (j = 0; j < count; j++)
      {
         zn_array_mul (temp, buf1, n, buf2, n, mod);
         zn_array_sub (buf3, temp, temp + n, n, mod);
      }
      
      t1 = get_cycle_counter ();
      
      free(temp);
   }
   else if (info->algo == ALGO_NEGAMUL_NUSS)
   {
      // nussbaumer version
      
      pmfvec_t vec1, vec2;
      pmfvec_init_nuss (vec1, info->lgL, mod);
      pmfvec_init_nuss (vec2, info->lgL, mod);

      // warm up
      for (j = 0; j < count/4; j++)
         nuss_mul (buf3, buf1, buf2, vec1, vec2);

      // do the actual profile
      t0 = get_cycle_counter ();

      for (j = 0; j < count; j++)
         nuss_mul (buf3, buf1, buf2, vec1, vec2);
      
      t1 = get_cycle_counter ();
      
      pmfvec_clear (vec2);
      pmfvec_clear (vec1);
   }
   else
      abort ();

   free (buf3);
   if (!info->sqr)
      free (buf2);
   free (buf1);
   
   zn_mod_clear (mod);
   
   return cycle_diff (t0, t1);
}



// end of file ****************************************************************
