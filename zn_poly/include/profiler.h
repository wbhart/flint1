/*
   profiler.h:  header file for profiling routines
   
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

   ============================================================================
   
   NOTE: this file includes code adapted from GMP 4.2.1 and FLINT 1.0.
   Please see the file COPYING for further information about licensing of
   this code.

*/

#ifndef ZNP_PROFILER_H
#define ZNP_PROFILER_H

#include <stdlib.h>
#include <stdio.h>
#include "zn_poly.h"


#ifdef __cplusplus
extern "C" {
#endif


/*
   Estimates the number of "cycles" per second (i.e. as reported by
   get_cycle_counter(), which might not actually be "cycles", but anyway).
   
   This is not supposed to be an accurate measurement; it's just a ballpark
   estimate, for the purpose of calibrating profiling/tuning runs to take
   a reasonable amount of time (w.r.t. typical human patience).
*/
double
estimate_cycle_scale_factor ();


/*
   Global scaling factor, set at startup via calibrate_cycle_scale_factor().
*/
extern double cycle_scale_factor;


/*
   Sets cycle_scale_factor according to estimate_cycle_scale_factor(),
   and prints some logging information to stderr.
*/
void
calibrate_cycle_scale_factor ();



/*
   ZNP_HAVE_CYCLE_COUNTER flag is set if a cycle counter is available on this
   platform. If not, we can't do any profiling or tuning.
   
   ZNP_PROFILE_CLOCK_MULTIPLIER is a multiplier used to decide how long to
   let profiling/tuning run for, in terms of the cycle count. This is a pretty
   rough estimate at the moment.

*/

#if defined(__GNUC__)


#if defined(__i386__) || defined(__x86_64__)

#define ZNP_HAVE_CYCLE_COUNTER 1

#define ZNP_PROFILE_CLOCK_MULTIPLIER 100

typedef unsigned long long  cycle_count_t;

typedef unsigned long long  cycle_diff_t;


ZNP_INLINE cycle_count_t
get_cycle_counter ()
{
   // hmmm... we're assuming "unsigned" is always 32 bits?
   unsigned hi;
   unsigned lo;

   __asm __volatile__ (
      "\t"
      "rdtsc            \n\t"
      "movl %%edx, %0   \n\t"
      "movl %%eax, %1   \n\t" 
       : "=r" (hi), "=r" (lo)
       : 
       : "%edx", "%eax");
       
   return (((unsigned long long)(hi)) << 32) + ((unsigned long long) lo);
}

#endif



#if defined (_ARCH_PPC) || defined (__powerpc__) || defined (__POWERPC__)     \
  || defined (__ppc__) || defined(__ppc64__)

#define ZNP_HAVE_CYCLE_COUNTER 1

#define ZNP_PROFILE_CLOCK_MULTIPLIER 1

typedef unsigned long long  cycle_count_t;

typedef unsigned long long  cycle_diff_t;


ZNP_INLINE cycle_count_t
get_cycle_counter ()
{
   // hmmm... we're assuming "unsigned" is always 32 bits?
   unsigned hi1, hi2;
   unsigned lo;

   do {
      __asm __volatile__ (
         "\t"
         "mftbu   %0 \n\t"
         "mftb    %1 \n\t"
         "mftbu   %2 \n\t"
       : "=r" (hi1), "=r" (lo), "=r" (hi2));
   }
   while (hi1 != hi2);

   return (((unsigned long long)(hi2)) << 32) + ((unsigned long long) lo);
}

#endif


#endif



#ifndef ZNP_HAVE_CYCLE_COUNTER

#define ZNP_HAVE_CYCLE_COUNTER 0
#define ZNP_PROFILE_CLOCK_MULTIPLIER 0

typedef int  cycle_count_t;
typedef int  cycle_diff_t;

ZNP_INLINE cycle_count_t
get_cycle_counter ()
{
   printf ("called get_cycle_counter() without a cycle counter!\n");
   abort ();
   return 0;
}

#endif



ZNP_INLINE cycle_diff_t
cycle_diff (cycle_count_t start, cycle_count_t stop)
{
   return stop - start;
}



/*
   Repeatedly runs target with parameter arg, with increasing values of count
   (count not allowed to get bigger than 10).
   
   Expects target to be a function that runs a certain profile count times
   and returns a total cycle count.
   
   Collects statistics regarding return values of target (i.e. cycle counts),
   in particular the median, upper and lower quartiles.
   
   Runs for at most limit seconds (as measured by the target!).
   Calls target at least once, and at most 50 times.
   
   If spread != NULL, it stores there the ratio
       (interquartile range) / (median).
       
   If samples != NULL, it stores there the number of times target was called.
*/
double
profile (double* spread, unsigned* samples,
         double (*target)(void* arg, unsigned long count),
         void* arg, double limit);


#ifdef __cplusplus
}
#endif

#endif

// end of file ****************************************************************
