/*
   array-profile-main.c:  program for profiling simple array operations
   
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


char* type_str[3] = {"add", "sub", "bfly"};
char* speed_str[2] = {"safe", "slim"};


void prof_main(int argc, char* argv[])
{
   ulong type, speed;
   double result, spread;
   
   printf("\n");

   // profile various butterfly loops
   for (type = 0; type < 3; type++)
   for (speed = 0; speed < 2; speed++)
   {
      ulong arg[2];
      arg[0] = type;
      arg[1] = speed;
      
      result = profile(&spread, NULL, profile_bfly, arg, 1.0) / 1000;
      
      printf(" %4s %s, cycles/coeff = %6.2lf (%.1lf%%)\n",
             type_str[type], speed_str[speed], result, 100 * spread);
   }
   
   // profile mpn_add_n and mpn_sub_n
   for (type = 0; type < 2; type++)
   {
      ulong arg[1];
      arg[0] = type;

      result = profile(&spread, NULL, profile_mpn_aors, arg, 1.0) / 1000;
   
      printf(" mpn_%s_n, cycles/limb  = %6.2lf (%.1lf%%)\n",
             type_str[type], result, 100 * spread);
   }

   // profile zn_array_scalar_mul
   {
      ulong arg[2];
      arg[1] = 0;
   
      arg[0] = ULONG_BITS - 1;
      result = profile(&spread, NULL, profile_scalar_mul, arg, 1.0) / 1000;
      printf("scalar_mul (> half-word), cycles/coeff = %6.2lf (%.1lf%%)\n",
             result, 100 * spread);

      arg[0] = ULONG_BITS/2 - 1;
      result = profile(&spread, NULL, profile_scalar_mul, arg, 1.0) / 1000;
      printf("scalar_mul (< half-word), cycles/coeff = %6.2lf (%.1lf%%)\n",
             result, 100 * spread);
   }

   // profile zn_array_scalar_mul with REDC
   {
      ulong arg[2];
      arg[1] = 1;

      arg[0] = ULONG_BITS;
      result = profile(&spread, NULL, profile_scalar_mul, arg, 1.0) / 1000;
      printf("scalar_mul (> half-word, non-slim, REDC), "
             "cycles/coeff = %6.2lf (%.1lf%%)\n",
             result, 100 * spread);

      arg[0] = ULONG_BITS - 1;
      result = profile(&spread, NULL, profile_scalar_mul, arg, 1.0) / 1000;
      printf("scalar_mul (> half-word, REDC), "
             "cycles/coeff = %6.2lf (%.1lf%%)\n",
             result, 100 * spread);

      arg[0] = ULONG_BITS/2 - 1;
      result = profile(&spread, NULL, profile_scalar_mul, arg, 1.0) / 1000;
      printf("scalar_mul (< half-word, REDC), "
             "cycles/coeff = %6.2lf (%.1lf%%)\n",
             result, 100 * spread);
   }
}


// end of file ****************************************************************
