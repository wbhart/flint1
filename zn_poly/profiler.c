/*
   profiler.c:  some profiling routines
   
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
#include "profiler.h"


#include <sys/time.h>
#include <sys/resource.h>


double cycle_scale_factor;


/*
   This function eats some CPU cycles. The number of cycles eaten is roughly
   proportional to the count parameter. This function exists only to ensure
   that the compiler is not smart enough to optimise away our cycle-eating.
*/
void use_up_cycles(unsigned long count)
{
   for (; count; count--)
   {
      unsigned long x[3] = {0, 1, 2};
      unsigned long y[3] = {0, 1, 2};
      unsigned long z[5];
      zn_mod_t mod;
      zn_mod_init(mod, 3);
      zn_array_mul(z, x, 3, y, 3, mod);
      zn_mod_clear(mod);
   }
}


double estimate_cycle_scale_factor()
{
   unsigned long count;
   
   for (count = 1; count < ULONG_MAX/4; )
   {
      struct rusage usage1, usage2;
      cycle_count_t t1, t2;

      // try using up some cycles and time how long (in microseconds) it takes
      getrusage(RUSAGE_SELF, &usage1);
      t1 = get_cycle_counter();
      use_up_cycles(count);
      t2 = get_cycle_counter();
      getrusage(RUSAGE_SELF, &usage2);
      
      long u1 = usage1.ru_utime.tv_usec;
      long u2 = usage2.ru_utime.tv_usec;
      
      if (u2 < u1)
         continue;

      // if we've used up at least 0.1s, estimate number of cycles per second
      if ((u2 - u1) > 100000)
         return (t2 - t1) / (u2 - u1) * 1000000;
         
      count *= 2;
   }
   
   abort();
}


void calibrate_cycle_scale_factor()
{
   fprintf(stderr, "Calibrating cycle counter... ");
   fflush(stderr);
   cycle_scale_factor = estimate_cycle_scale_factor();
   fprintf(stderr, "ok (%.2le)\n", cycle_scale_factor);
}



// borrowed from gcc documentation:
int compare_doubles(const void *a, const void *b)
{
   const double *da = (const double*) a;
   const double *db = (const double*) b;
     
   return (*da > *db) - (*da < *db);
}


double profile(double* spread, unsigned* samples,
               double (*target)(void* arg, unsigned long count),
               void* arg, double limit)
{
   const unsigned max_times = 50;
   double times[max_times];

   // convert limit from seconds to cycles
   limit *= cycle_scale_factor;

   unsigned n;
   unsigned long count = 1;
   double elapsed = 0.0;

   for (n = 0; n < max_times  &&  elapsed < limit; n++)
   {
      double time = target(arg, count);
      elapsed += time;
      times[n] = time / count;

      if (time < limit/100)
         count++;
      if (count > 10)
         count = 10;
   }

   // get median, lower quartile, upper quartile of measured times
   qsort(times, n, sizeof(double), compare_doubles);
   double median = 0.5 * (times[(n-1)/2] + times[n/2]);
   double q1 = 0.25 * (times[(n-1)/4] + times[n/4] +
                       times[(n+1)/4] + times[(n+2)/4]);
   double q3 = 0.25 * (times[3*n/4] + times[(3*n-1)/4] +
                       times[(3*n-2)/4] + times[(3*n-3)/4]);
   
   if (samples)
      *samples = n;

   if (spread)
      *spread = (q3 - q1) / median;

   return median;
}



// end of file ****************************************************************
