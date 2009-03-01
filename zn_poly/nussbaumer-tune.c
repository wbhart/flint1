/*
   nussbaumer-tune.c:  tuning program for negacyclic nussbaumer multiplication
   
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
#include "zn_poly_internal.h"
#include "profiler.h"



// the maximum crossover we allow for switching from KS to nussbaumer
#define MAX_NUSSBAUMER_CROSSOVER 14


/*
   For each modulus size, finds crossover points between KS and nussbaumer
   multiplication (or squaring, if the squaring flag is set). Store these in
   the global crossover table, and writes some logging information to flog.
*/
void tune_nussbaumer(FILE* flog, int squaring, int verbose)
{
   unsigned bits;

   fprintf(flog, "Nussbaumer %s: ", squaring ? "sqr" : "mul");
   fflush(flog);

   // how long we are willing to wait for each profile run
   const double speed = 0.001;

   // run tuning process for each modulus size
   for (bits = 2; bits <= ULONG_BITS; bits++)
   {
      unsigned crossover;
   
      profile_negamul_info_t info[2];
      info[0]->squaring = info[1]->squaring = squaring;
      info[0]->n = info[1]->n = (1UL << (bits - 1)) + 1;

      info[0]->algo = ALGO_NEGAMUL_FALLBACK;
      info[1]->algo = ALGO_NEGAMUL_NUSSBAUMER;

      for (crossover = 3; crossover < MAX_NUSSBAUMER_CROSSOVER; crossover++)
      {
         double result[2];
      
         info[0]->lg_len = info[1]->lg_len = crossover;
         
         // need nussbaumer to win three times
         unsigned trial;
         for (trial = 0; trial < 3; trial++)
         {
            result[0] = profile(NULL, NULL, profile_negamul, info[0], speed);
            result[1] = profile(NULL, NULL, profile_negamul, info[1], speed);

            if (result[0] < result[1])
               break;
         }
         
         if (trial == 3)
            break;    // found crossover
      }
      
      // If it looks like KS is always going to win, just always use KS
      // (for instance this might happen if GMP's huge integer multiplication
      // improves stupendously)
      if (crossover == MAX_NUSSBAUMER_CROSSOVER)
         crossover = 1000;
      
      if (squaring)
         tuning_info[bits].nuss_sqr_crossover = crossover;
      else
         tuning_info[bits].nuss_mul_crossover = crossover;

      if (verbose)
      {
         fprintf(flog, "\nbits = %u, cross to Nussbaumer at length ", bits);
         if (crossover == 1000)
            fprintf(flog, "infinity");
         else
            fprintf(flog, "%lu", 1UL << crossover);
      }
      else
         fprintf(flog, ".");

      fflush(flog);
   }
   
   fprintf(flog, "\n");
}


// end of file ****************************************************************
