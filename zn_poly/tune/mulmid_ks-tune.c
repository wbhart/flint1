/*
   mulmid_ks-tune.c:  tuning program for mulmid_ks.c module
   
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
#include <stdint.h>
#include "support.h"
#include "zn_poly_internal.h"
#include "profiler.h"


/*
   For each modulus bitsize, finds approximate threshold between KS1, KS2, and
   KS4 middle product algorithms.

   Store these in the global tuning table, and writes some logging information
   to _flog_.
*/
void
tune_mulmid_KS (FILE* flog, int verbose)
{
   unsigned b;

   fprintf (flog, "KS1/2/4 mulmid: ");
   fflush (flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.0001;
   
   // run tuning process for each modulus size
   for (b = 2; b <= ULONG_BITS; b++)
   {
      // thresholds for KS1->KS2 and KS2->KS4
      unsigned long thresh[2];
   
      // which = 0 means comparing KS1 vs KS2
      // which = 1 means comparing KS2 vs KS4
      unsigned which;
      for (which = 0; which < 2; which++)
      {
         double result[2];

         profile_info_t info[2];
         info[0]->m = info[1]->m = (1UL << (b - 1)) + 1;

         info[0]->algo = which ? ALGO_MULMID_KS2_REDC : ALGO_MULMID_KS1_REDC;
         info[1]->algo = which ? ALGO_MULMID_KS4_REDC : ALGO_MULMID_KS2_REDC;
         
         // find an upper bound, where 2nd algorithm appears to be safely
         // ahead of 1st algorithm
         size_t upper;
         int found = 0;
         for (upper = 45; upper <= 16384 && !found; upper = 2 * upper)
         {
            info[0]->n1 = info[1]->n1 = 2 * upper;
            info[0]->n2 = info[1]->n2 = upper;
            
            result[0] = profile (NULL, NULL, profile_mulmid, info[0], speed);
            result[1] = profile (NULL, NULL, profile_mulmid, info[1], speed);
            
            if (result[1] < 0.95 * result[0])
               found = 1;
         }

         if (!found)
         {
            // couldn't find a reasonable upper bound
            thresh[which] = SIZE_MAX;
            continue;
         }
         
         // find a lower bound, where 1st algorithm appears to be safely
         // ahead of 2nd algorithm
         size_t lower;
         found = 0;
         for (lower = upper / 2; lower >= 2 && !found; lower = lower / 2)
         {
            info[0]->n1 = info[1]->n1 = 2 * lower;
            info[0]->n2 = info[1]->n2 = lower;
            
            result[0] = profile (NULL, NULL, profile_mulmid, info[0], speed);
            result[1] = profile (NULL, NULL, profile_mulmid, info[1], speed);
            
            if (result[1] > 1.05 * result[0])
               found = 1;
         }

         if (!found)
         {
            // couldn't find a reasonable lower bound
            thresh[which] = 0;
            continue;
         }

         // subdivide [lower, upper] into intervals and sample at each endpoint
         double ratio = (double) upper / (double) lower;
         const int max_intervals = 30;
         size_t points[max_intervals + 1];
         double score[max_intervals + 1];
         unsigned i;
         for (i = 0; i <= max_intervals; i++)
         {
            points[i] = ceil (lower * pow (ratio, (double) i / max_intervals));
            info[0]->n1 = info[1]->n1 = 2 * points[i];
            info[0]->n2 = info[1]->n2 = points[i];
            result[0] = profile (NULL, NULL, profile_mulmid, info[0], speed);
            result[1] = profile (NULL, NULL, profile_mulmid, info[1], speed);
            score[i] = result[1] / result[0];
         }
         
         // estimate threshold
         unsigned count = 0;
         for (i = 0; i <= max_intervals; i++)
            if (score[i] > 1.0)
               count++;
         thresh[which] = (size_t)
              ceil (lower * pow (ratio, (double) count / (max_intervals + 1)));
      }

      if (verbose)
      {
         fprintf (flog, "\nbits = %u, cross to KS2 at ", b);

         if (thresh[0] == SIZE_MAX)
            fprintf (flog, "infinity");
         else
            fprintf (flog, "%lu", thresh[0]);

         fprintf (flog, ", cross to KS4 at ");

         if (thresh[1] == SIZE_MAX)
            fprintf (flog, "infinity");
         else
            fprintf (flog, "%lu", thresh[1]);
      }
      else
         fprintf (flog, ".");

      fflush (flog);

      tuning_info[b].mulmid_KS2_thresh = thresh[0];
      tuning_info[b].mulmid_KS4_thresh = thresh[1];
   }

   fprintf (flog, "\n");
}


// end of file ****************************************************************
