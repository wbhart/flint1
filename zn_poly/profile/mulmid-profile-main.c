/*
   mulmid-profile-main.c:  program for profiling middle products
   
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
   of bits and polynomial length (it does a (2*n) x n middle product).
   
   active[i] is a flag indicating whether the algorithm indexed by i
   should be profiled (see ALGO_MULMID_xyz constants).
*/
void
do_line (int* active, unsigned b, size_t n)
{
   char* names[] = {"best", "fallback", "ks1", "ks1_redc", "ks2", "ks2_redc",
                    "ks3", "ks3_redc", "ks4", "ks4_redc", "fft"};

   profile_info_t info;
   
   info->n1 = 2 * n;
   info->n2 = n;

   // choose an odd modulus exactly b bits long
   info->m = (1UL << (b - 1)) + 2 * random_ulong (1UL << (b - 2)) + 1;

   printf ("len = %5lu, bits = %2u", n, b);
   fflush (stdout);

   int algo;
   for (algo = 0; algo < 11; algo++)
   {
      if (active[algo])
      {
         info->algo = algo;
         double result, spread;
         result = profile (&spread, NULL, profile_mulmid, &info, 1.0);
         printf (", %s = %.3le (%.1lf%%)", names[algo], result, 100 * spread);
         fflush (stdout);
      }
   }
   
   printf ("\n");
}


void
prof_main (int argc, char* argv[])
{
   // read command line arguments
   
   // can include the strings "best", "ks1", "ks1_redc", "ks2", "ks2_redc",
   //                         "ks3", "ks3_redc", "ks4", "ks4_redc",
   //                         "fallback", "fft"
   // to select various algorithms
   
   // if you do "bits <nnn>" then only that number of bits will be profiled
   // otherwise it ranges over various bitsizes
   
   // if you do "length <nnn>" then only that length will be profiled
   // otherwise it ranges over various lengths

   int active[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
   int any_active = 0;
   
   int do_one_bits = 0;
   int chosen_bits = 0;

   int do_one_length = 0;
   ulong chosen_length = 0;
   
   int i;
   for (i = 1; i < argc; i++)
   {
      if (!strcmp (argv[i], "best"))
         active[ALGO_MULMID_BEST] = any_active = 1;
      else if (!strcmp (argv[i], "fallback"))
         active[ALGO_MULMID_FALLBACK] = any_active = 1;
      else if (!strcmp (argv[i], "ks1"))
         active[ALGO_MULMID_KS1] = any_active = 1;
      else if (!strcmp (argv[i], "ks1_redc"))
         active[ALGO_MULMID_KS1_REDC] = any_active = 1;
      else if (!strcmp (argv[i], "ks2"))
         active[ALGO_MULMID_KS2] = any_active = 1;
      else if (!strcmp (argv[i], "ks2_redc"))
         active[ALGO_MULMID_KS2_REDC] = any_active = 1;
      else if (!strcmp (argv[i], "ks3"))
         active[ALGO_MULMID_KS3] = any_active = 1;
      else if (!strcmp (argv[i], "ks3_redc"))
         active[ALGO_MULMID_KS3_REDC] = any_active = 1;
      else if (!strcmp (argv[i], "ks4"))
         active[ALGO_MULMID_KS4] = any_active = 1;
      else if (!strcmp (argv[i], "ks4_redc"))
         active[ALGO_MULMID_KS4_REDC] = any_active = 1;
      else if (!strcmp (argv[i], "fft"))
         active[ALGO_MULMID_FFT] = any_active = 1;

      else if (!strcmp (argv[i], "bits"))
      {
         do_one_bits = 1;
         chosen_bits = atoi (argv[++i]);
      }
      else if (!strcmp (argv[i], "length"))
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
   
   if (!any_active)
      active[0] = 1;     // profile plain multiplication if nothing selected

   // bitsizes to use by default if none are selected
   unsigned bitsizes[9] = {4, 8, 16, 24, 32, 40, 48, 56, 64};

   int j;
   size_t len;

   if (do_one_bits && do_one_length)
   {
      do_line (active, chosen_bits, chosen_length);
   }
   else if (do_one_bits && !do_one_length)
   {
      // loop over lengths, spaced out logarithmically
      for (j = 0; j < 120; j++)
      {
         size_t new_len = (size_t) floor (pow (1.1, (double) j));
         if (new_len == len)
            continue;
         len = new_len;
         
         do_line (active, chosen_bits, len);
      }
   }
   else if (!do_one_bits && do_one_length)
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof (bitsizes) / sizeof (bitsizes[0]); i++)
         do_line (active, bitsizes[i], chosen_length);
   }
   else   // neither bits nor length is fixed
   {
      // loop over bitsizes in above table
      for (i = 0; i < sizeof (bitsizes) / sizeof (bitsizes[0]); i++)
      {
         // loop over lengths, spaced out logarithmically
         for (j = 0; j < 120; j++)
         {
            size_t new_len = (size_t) floor (pow (1.1, (double) j));
            if (new_len == len)
               continue;
            len = new_len;

            do_line (active, bitsizes[i], len);
         }
         
         printf("-------------------------------------------\n");
      }
   }
}

// end of file ****************************************************************
