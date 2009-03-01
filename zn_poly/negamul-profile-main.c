/*
   negamul-profile-main.c:  program for profiling negacyclic multiplication
   
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
void do_line(unsigned bits, unsigned lg_len, int squaring)
{
   profile_negamul_info_t info;
   
   info->lg_len = lg_len;

   // choose an odd modulus exactly _bits_ bits long
   info->n = (1UL << (bits - 1))
                + 2 * random_ulong(1UL << (bits - 2)) + 1;
   info->squaring = squaring;

   size_t len = 1UL << lg_len;

   printf("len = %6lu, bits = %2u", len, bits);
   fflush(stdout);

   double result, spread;
   
   info->algo = ALGO_NEGAMUL_FALLBACK;
   result = profile(&spread, NULL, profile_negamul, &info, 1.0);
   printf(", fallback = %.3le (%.1lf%%)", result, 100 * spread);
   fflush(stdout);

   info->algo = ALGO_NEGAMUL_NUSSBAUMER;
   result = profile(&spread, NULL, profile_negamul, &info, 1.0);
   printf(", nussbaumer = %.3le (%.1lf%%)", result, 100 * spread);
   fflush(stdout);
   
   printf("\n");
}


void prof_main(int argc, char* argv[])
{
   // read command line arguments
   
   // if you do "bits <nnn>" then only that number of bits will be profiled
   // otherwise it ranges over various bitsizes
   
   // if you do "lg_len <nnn>" then only that length will be profiled
   // otherwise it ranges over various lengths

   // can also include "sqr" anywhere, which means to profile squaring

   int do_one_bits = 0;
   int chosen_bits = 0;
   int squaring = 0;

   int do_one_length = 0;
   unsigned chosen_lg_len = 0;
   
   int i;
   for (i = 1; i < argc; i++)
   {
      if (!strcmp(argv[i], "bits"))
      {
         do_one_bits = 1;
         chosen_bits = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i], "lg_len"))
      {
         do_one_length = 1;
         chosen_lg_len = atol(argv[++i]);
      }
      else if (!strcmp(argv[i], "sqr"))
      {
         squaring = 1;
      }
      else
      {
         printf("unknown option %s\n", argv[i]);
         exit(1);
      }
   }
   
   // bitsizes to use by default if none are selected
   unsigned bitsizes[9] = {4, 8, 16, 24, 32, 40, 48, 56, 64};

   unsigned lg_len;
   unsigned bits;

   if (do_one_bits && do_one_length)
   {
      do_line(chosen_bits, chosen_lg_len, squaring);
   }
   else if (do_one_bits && !do_one_length)
   {
      // loop over lengths
      for (lg_len = 4; lg_len <= 16; lg_len++)
         do_line(chosen_bits, lg_len, squaring);
   }
   else if (!do_one_bits && do_one_length)
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof(bitsizes) / sizeof(bitsizes[0]); i++)
         do_line(bitsizes[i], chosen_lg_len, squaring);
   }
   else   // neither bits nor length is fixed
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof(bitsizes) / sizeof(bitsizes[0]); i++)
      {
         // loop over lengths, spaced out logarithmically
         for (lg_len = 4; lg_len <= 16; lg_len++)
            do_line(bitsizes[i], lg_len, squaring);
         
         printf("-------------------------------------------\n");
      }
   }
}


// end of file ****************************************************************
