/*
   negamul-profile-main.c:  program for profiling negacyclic multiplication
   
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
   Performs and prints one line of profiling output, for the given number
   of bits and polynomial length.
*/
void
do_line (unsigned b, unsigned lgL, int sqr)
{
   profile_info_t info;
   
   info->lgL = lgL;

   // choose an odd modulus exactly _bits_ bits long
   info->m = (1UL << (b - 1)) + 2 * random_ulong (1UL << (b - 2)) + 1;
   info->sqr = sqr;

   size_t n = 1UL << lgL;

   printf ("len = %6lu, bits = %2u", n, b);
   fflush (stdout);

   double result, spread;
   
   info->algo = ALGO_NEGAMUL_FALLBACK;
   result = profile (&spread, NULL, profile_negamul, &info, 1.0);
   printf (", fallback = %.3le (%.1lf%%)", result, 100 * spread);
   fflush (stdout);

   info->algo = ALGO_NEGAMUL_NUSS;
   result = profile (&spread, NULL, profile_negamul, &info, 1.0);
   printf (", nuss = %.3le (%.1lf%%)", result, 100 * spread);
   fflush (stdout);
   
   printf ("\n");
}


void
prof_main (int argc, char* argv[])
{
   // read command line arguments
   
   // if you do "bits <nnn>" then only that number of bits will be profiled
   // otherwise it ranges over various bitsizes
   
   // if you do "lgL <nnn>" then only that length will be profiled
   // otherwise it ranges over various lengths

   // can also include "sqr" anywhere, which means to profile squaring

   int do_one_bits = 0;
   int chosen_bits = 0;
   int sqr = 0;

   int do_one_length = 0;
   unsigned chosen_lgL = 0;
   
   int i;
   for (i = 1; i < argc; i++)
   {
      if (!strcmp (argv[i], "bits"))
      {
         do_one_bits = 1;
         chosen_bits = atoi (argv[++i]);
      }
      else if (!strcmp (argv[i], "lgL"))
      {
         do_one_length = 1;
         chosen_lgL = atol (argv[++i]);
      }
      else if (!strcmp (argv[i], "sqr"))
      {
         sqr = 1;
      }
      else
      {
         printf ("unknown option %s\n", argv[i]);
         exit (1);
      }
   }
   
   // bitsizes to use by default if none are selected
   unsigned bitsizes[9] = {4, 8, 16, 24, 32, 40, 48, 56, 64};

   unsigned lgL;

   if (do_one_bits && do_one_length)
   {
      do_line (chosen_bits, chosen_lgL, sqr);
   }
   else if (do_one_bits && !do_one_length)
   {
      // loop over lengths
      for (lgL = 4; lgL <= 16; lgL++)
         do_line (chosen_bits, lgL, sqr);
   }
   else if (!do_one_bits && do_one_length)
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof (bitsizes) / sizeof (bitsizes[0]); i++)
         do_line (bitsizes[i], chosen_lgL, sqr);
   }
   else   // neither bits nor length is fixed
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof (bitsizes) / sizeof (bitsizes[0]); i++)
      {
         // loop over lengths, spaced out logarithmically
         for (lgL = 4; lgL <= 16; lgL++)
            do_line (bitsizes[i], lgL, sqr);
         
         printf ("-------------------------------------------\n");
      }
   }
}


// end of file ****************************************************************
