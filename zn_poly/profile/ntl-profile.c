/*
   ntl-profile.c:  routines for profiling NTL
   
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
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>


extern "C" double
profile_mul_ntl (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   cycle_count_t t0, t1;
   
   if (info->m < (ulong) NTL_SP_BOUND)
   {
      // zz_pX version
   
      NTL::zz_pX f1, f2, g;
      NTL::zz_p::init (info->m);
      
      size_t i;
      for (i = 0; i < info->n1; i++)
         SetCoeff (f1, i, random_ulong (info->m));
      for (i = 0; i < info->n2; i++)
         SetCoeff (f2, i, random_ulong (info->m));

      if (info->sqr)
      {
         // warm up
         ulong j;
         for (j = 0; j < count; j++)
            sqr (g, f1);
            
         t0 = get_cycle_counter ();

         for (j = 0; j < count; j++)
            sqr (g, f1);

         t1 = get_cycle_counter ();
      }
      else
      {
         // warm up
         ulong j;
         for (j = 0; j < count; j++)
            mul (g, f1, f2);
            
         t0 = get_cycle_counter ();

         for (j = 0; j < count; j++)
            mul (g, f1, f2);

         t1 = get_cycle_counter ();
      }
   }
   else
   {
      // ZZ_pX version

      NTL::ZZ_pX f1, f2, g;
      NTL::ZZ_p::init (NTL::to_ZZ (info->m));
      
      size_t i;
      for (i = 0; i < info->n1; i++)
         SetCoeff (f1, i, random_ulong (info->m));
      for (i = 0; i < info->n2; i++)
         SetCoeff (f2, i, random_ulong (info->m));

      if (info->sqr)
      {
         // warm up
         ulong j;
         for (j = 0; j < count; j++)
            sqr (g, f1);
            
         t0 = get_cycle_counter ();

         for (j = 0; j < count; j++)
            sqr (g, f1);

         t1 = get_cycle_counter ();
      }
      else
      {
         // warm up
         ulong j;
         for (j = 0; j < count; j++)
            mul (g, f1, f2);
            
         t0 = get_cycle_counter ();

         for (j = 0; j < count; j++)
            mul (g, f1, f2);

         t1 = get_cycle_counter ();
      }
   }
   
   return cycle_diff (t0, t1);
}



extern "C" double
profile_invert_ntl (void* arg, unsigned long count)
{
   profile_info_struct* info = (profile_info_struct*) arg;

   cycle_count_t t0, t1;
   
   if (info->m < (ulong) NTL_SP_BOUND)
   {
      // zz_pX version
   
      NTL::zz_pX f1, f2, g;
      NTL::zz_p::init (info->m);
      
      size_t i;
      SetCoeff (f1, 0, 1);
      for (i = 1; i < info->n; i++)
         SetCoeff (f1, i, random_ulong (info->m));

      // warm up
      ulong j;
      for (j = 0; j < count; j++)
         InvTrunc (g, f1, info->n);
         
      t0 = get_cycle_counter ();

      for (j = 0; j < count; j++)
         InvTrunc (g, f1, info->n);

      t1 = get_cycle_counter ();
   }
   else
   {
      // ZZ_pX version

      NTL::ZZ_pX f1, f2, g;
      NTL::ZZ_p::init (NTL::to_ZZ (info->m));

      size_t i;
      SetCoeff (f1, 0, 1);
      for (i = 1; i < info->n; i++)
         SetCoeff (f1, i, random_ulong (info->m));

      // warm up
      ulong j;
      for (j = 0; j < count; j++)
         InvTrunc (g, f1, info->n);
         
      t0 = get_cycle_counter ();

      for (j = 0; j < count; j++)
         InvTrunc (g, f1, info->n);

      t1 = get_cycle_counter ();
   }

   return cycle_diff (t0, t1);
}


// end of file ****************************************************************
