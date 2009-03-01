/*
   invert-profile-main.c:  program for profiling power series inversion
                           (zn_poly algorithms, and optionally NTL)
   
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


/*
   Note: if this file is compiled with the constant PROFILE_NTL defined, then
   it will include support for NTL profiling, and needs to be compiled as C++.
   Otherwise it can be compiled as a C program.
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
   
   active[i] is a flag indicating whether the algorithm indexed by i
   should be profiled (see ALGO_INVERT_xyz constants).
*/
void do_line(int* active, unsigned bits, size_t len)
{
   char* names[] = {"best", "ntl"};

   profile_invert_info_t info;
   info->len = len;

   // choose an odd modulus exactly _bits_ bits long
   info->n = (1UL << (bits - 1))
                + 2 * random_ulong(1UL << (bits - 2)) + 1;

   printf("len = %5lu, bits = %2u", len, bits);
   fflush(stdout);

   int algo;
   for (algo = 0; algo < 2; algo++)
   {
      if (active[algo])
      {
         info->algo = algo;
         
         double spread, result;
         result = profile(&spread, NULL, profile_invert, &info, 1.0);
         printf(", %s = %.3le (%.1lf%%)", names[algo], result, 100 * spread);
         fflush(stdout);
      }
   }
   
   printf("\n");
}


#if __cplusplus
extern "C"
#endif
void prof_main(int argc, char* argv[])
{
   // read command line arguments
   
   // can include the strings "best", "ntl"
   // to select various algorithms
   
   // if you do "bits <nnn>" then only that number of bits will be profiled
   // otherwise it ranges over various bitsizes
   
   // if you do "length <nnn>" then only that length will be profiled
   // otherwise it ranges over various lengths

   int active[2] = {0, 0};
   int any_active = 0;
   
   int do_one_bits = 0;
   int chosen_bits = 0;

   int do_one_length = 0;
   ulong chosen_length = 0;
   
   int i;
   for (i = 1; i < argc; i++)
   {
      if (!strcmp(argv[i], "best"))
         active[ALGO_INVERT_BEST] = any_active = 1;
      else if (!strcmp(argv[i], "ntl"))
         active[ALGO_INVERT_NTL] = any_active = 1;

      else if (!strcmp(argv[i], "bits"))
      {
         do_one_bits = 1;
         chosen_bits = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i], "length"))
      {
         do_one_length = 1;
         chosen_length = atol(argv[++i]);
      }
      else
      {
         printf("unknown option %s\n", argv[i]);
         exit(1);
      }
   }
   
   if (!any_active)
      active[0] = 1;     // profile plain multiplication if nothing selected

   // bitsizes to use by default if none are selected
   unsigned bitsizes[9] = {4, 8, 16, 24, 32, 40, 48, 56, 64};

   int j;
   size_t len;
   unsigned bits;

   if (do_one_bits && do_one_length)
   {
      do_line(active, chosen_bits, chosen_length);
   }
   else if (do_one_bits && !do_one_length)
   {
      // loop over lengths, spaced out logarithmically
      for (j = 0; j < 120; j++)
      {
         size_t new_len = (size_t) floor(pow(1.1, (double) j));
         if (new_len == len)
            continue;
         len = new_len;
         
         do_line(active, chosen_bits, len);
      }
   }
   else if (!do_one_bits && do_one_length)
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof(bitsizes) / sizeof(bitsizes[0]); i++)
         do_line(active, bitsizes[i], chosen_length);
   }
   else   // neither bits nor length is fixed
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof(bitsizes) / sizeof(bitsizes[0]); i++)
      {
         // loop over lengths, spaced out logarithmically
         for (j = 0; j < 120; j++)
         {
            size_t new_len = (size_t) floor(pow(1.1, (double) j));
            if (new_len == len)
               continue;
            len = new_len;

            do_line(active, bitsizes[i], len);
         }
         
         printf("-------------------------------------------\n");
      }
   }
}

// end of file ****************************************************************
