/*
   nuss-tune.c:  tuning program for negacyclic nussbaumer multiplication
   
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
#include "zn_poly_internal.h"
#include "profiler.h"



// the maximum threshold we allow for switching from KS to nussbaumer
#define MAX_NUSS_THRESH 14


/*
   For each modulus size, finds threshold between KS and nussbaumer
   multiplication (or squaring, if the squaring flag is set). Store these in
   the global threshold table, and writes some logging information to flog.
*/
void
tune_nuss (FILE* flog, int sqr, int verbose)
{
   unsigned b;

   fprintf (flog, "      nuss %s: ", sqr ? "sqr" : "mul");
   fflush (flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.001;

   // run tuning process for each modulus size
   for (b = 2; b <= ULONG_BITS; b++)
   {
      unsigned thresh;
   
      profile_info_t info[2];
      info[0]->sqr = info[1]->sqr = sqr;
      info[0]->m = info[1]->m = (1UL << (b - 1)) + 1;

      info[0]->algo = ALGO_NEGAMUL_FALLBACK;
      info[1]->algo = ALGO_NEGAMUL_NUSS;

      for (thresh = 3; thresh < MAX_NUSS_THRESH; thresh++)
      {
         double result[2];
      
         info[0]->lgL = info[1]->lgL = thresh;
         
         // need nussbaumer to win three times
         unsigned trial;
         for (trial = 0; trial < 3; trial++)
         {
            result[0] = profile (NULL, NULL, profile_negamul, info[0], speed);
            result[1] = profile (NULL, NULL, profile_negamul, info[1], speed);

            if (result[0] < result[1])
               break;
         }
         
         if (trial == 3)
            break;    // found threshold
      }
      
      // If it looks like KS is always going to win, just always use KS
      // (for instance this might happen if GMP's huge integer multiplication
      // improves stupendously)
      if (thresh == MAX_NUSS_THRESH)
         thresh = 1000;
      
      if (sqr)
         tuning_info[b].nuss_sqr_thresh = thresh;
      else
         tuning_info[b].nuss_mul_thresh = thresh;

      if (verbose)
      {
         fprintf (flog, "\nbits = %u, cross to Nussbaumer at length ", b);
         if (thresh == 1000)
            fprintf (flog, "infinity");
         else
            fprintf (flog, "%lu", 1UL << thresh);
      }
      else
         fprintf (flog, ".");

      fflush (flog);
   }
   
   fprintf (flog, "\n");
}


// end of file ****************************************************************
