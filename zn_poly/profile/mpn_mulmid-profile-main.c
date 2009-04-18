/*
   mpn_mulmid-profile-main.c:  program for profiling mpn middle products
   
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


#include <string.h>
#include <math.h>
#include "support.h"
#include "profiler.h"
#include "zn_poly_internal.h"
#include "zn_poly.h"


/*
   Performs and prints one line of profiling output, for the given
   integer length
*/
void
do_line (size_t n)
{
   profile_info_t info;

   printf ("n = %5lu", n);
   fflush (stdout);

   double spread, result;

   info->n1 = 2 * n - 1;
   info->n2 = n;
   result = profile (&spread, NULL, profile_mpn_smp_basecase,
                     &info, 1.0);
   printf (", %.3le (%.1lf%%)", result, 100 * spread);

   if (n >= 2)
   {
      info->n = n;
      result = profile (&spread, NULL, profile_mpn_smp_kara,
                        &info, 1.0);
      printf (", %.3le (%.1lf%%)", result, 100 * spread);
   }
   else
      printf (", N/A             ");

   result = profile (&spread, NULL, profile_mpn_smp, &info, 1.0);
   printf (", %.3le (%.1lf%%)", result, 100 * spread);

   info->n1 = info->n2 = n;
   result = profile (&spread, NULL, profile_mpn_mul, &info, 1.0);
   printf (", %.3le (%.1lf%%)", result, 100 * spread);

   info->n1 = 2 * n - 1;
   info->n2 = n;
   result = profile (&spread, NULL, profile_mpn_mulmid_fallback, &info, 1.0);
   printf (", %.3le (%.1lf%%)", result, 100 * spread);

   printf ("\n");
}


#if __cplusplus
extern "C"
#endif
void
prof_main (int argc, char* argv[])
{
   // read command line arguments
   
   // if you do "length <nnn>" then only that length will be profiled
   // otherwise it ranges over various lengths

   int do_one_length = 0;
   ulong chosen_length = 0;
   
   int i;
   for (i = 1; i < argc; i++)
   {
      if (!strcmp (argv[i], "length"))
      {
         do_one_length = 1;
         chosen_length = atol (argv[++i]);
      }
      else
      {
         printf ("unknown option %s\n", argv[i]);
         exit (1);
      }
   }
   
   int j;
   size_t n;

   printf ("fields: smp_basecase, smp_kara, "
           "smp, mpn_mul, mulmid_fallback\n");

   if (do_one_length)
   {
      do_line (chosen_length);
   }
   else
   {
      // loop over lengths, spaced out logarithmically
      for (n = 1; n <= 100; n++)
         do_line (n);
   }
}

// end of file ****************************************************************
