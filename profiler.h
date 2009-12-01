/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/******************************************************************************

 Timing/profiling

 (C) 2007 William Hart and David Harvey

******************************************************************************/


#undef ulong // namespace pollution
#include <time.h>
#if ! defined (__TINYC__)
#include <sys/time.h>
#endif
#undef ulong 
#define ulong unsigned long

#ifndef FLINT_PROFILER_H
#define FLINT_PROFILER_H

#ifdef __cplusplus
 extern "C" {
#endif

/* 
   Quick and dirty timing (gives wall and user time - useful for parallel programming)

   Example usage:

   timeit_t t0;

   ///....

   timeit_start(t0);

   /// do stuff, take some time

   timeit_stop(t0);

   printf("cpu = %ld ms  wall = %ld ms\n", t0->cpu, t0->wall);
*/

typedef struct
{
   long cpu;
   long wall;
} timeit_struct;

typedef timeit_struct timeit_t[1];

static inline
void timeit_start(timeit_t t)
{
#ifndef __TINYC__
   struct timeval tv;
   gettimeofday(&tv, 0);
   t->wall = - tv.tv_sec * 1000 - tv.tv_usec / 1000;
#else
   t->wall = 0;
#endif
   t->cpu = - clock() * 1000 / CLOCKS_PER_SEC;
}


static inline
void timeit_stop(timeit_t t)
{
#ifndef __TINYC__
   struct timeval tv;
   gettimeofday(&tv, 0);
   t->wall += tv.tv_sec * 1000 + tv.tv_usec / 1000;
#else
   t->wall = 0;
#endif
   t->cpu += clock() * 1000 / CLOCKS_PER_SEC;
}

// More refined timing routines
// number of independent global clocks

#define FLINT_NUM_CLOCKS 20


// If this flag is set, profiling will use a cycle counter *if one is
// available* (otherwise this flag is ignored)
#define FLINT_USE_CYCLE_COUNTER 1

// cycles/second
#define FLINT_CLOCKSPEED 2800121000.0


extern double clock_last[FLINT_NUM_CLOCKS];
extern double clock_accum[FLINT_NUM_CLOCKS];



#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))

// Relative timings on X86 machines running gcc

#define FLINT_HAVE_CYCLE_COUNTER 1

static inline double get_cycle_counter()
{
   // dirty: do we need to ensure these are 32-bit types?
   unsigned hi;
   unsigned lo;

   __asm("rdtsc; movl %%edx,%0; movl %%eax,%1" 
       : "=r" (hi), "=r" (lo)
       : 
       : "%edx", "%eax");

   return (double) hi * (1 << 30) * 4 + lo;
}

#else

#define FLINT_HAVE_CYCLE_COUNTER 0

#endif



/*
Here we define FLINT_CLOCK_SCALE_FACTOR, which converts the output of
get_current_time() into microseconds
*/
#if FLINT_HAVE_CYCLE_COUNTER && FLINT_USE_CYCLE_COUNTER
// microseconds per cycle
#define FLINT_CLOCK_SCALE_FACTOR (1000000.0 / FLINT_CLOCKSPEED)
#else
// we'll use getrusage, which is already in microseconds
#define FLINT_CLOCK_SCALE_FACTOR 1.0
#endif


static inline double get_current_time()
{
#if FLINT_HAVE_CYCLE_COUNTER && FLINT_USE_CYCLE_COUNTER
   return get_cycle_counter();
#else
   // user time in microseconds
   struct rusage x;
   getrusage(RUSAGE_SELF, &x);
   return x.ru_utime.tv_sec * 1000000.0 +  x.ru_utime.tv_usec;
#endif
}


static inline void init_clock(unsigned long n)
{
   clock_accum[n] = 0.0;
}

static inline void init_all_clocks()
{
   unsigned long i;
   for (i = 0; i < FLINT_NUM_CLOCKS; i++)
      clock_accum[i] = 0.0;
}

static inline double get_clock(unsigned long n)
{
   return clock_accum[n] * FLINT_CLOCK_SCALE_FACTOR;
}

static inline void start_clock(unsigned long n)
{
   clock_last[n] = get_current_time();
}

static inline void stop_clock(unsigned long n)
{
   double now = get_current_time();
   clock_accum[n] += (now - clock_last[n]);
}


/******************************************************************************

   Framework for repeatedly sampling a single target
   
******************************************************************************/

// A profiling target (a function called with one argument and an
// iteration count.)
typedef void (*profile_target_t)(void* arg, unsigned long count);

/*
   Calls target(arg) repeatedly, adjusting the iteration count based on the
   observed running times.
   
   The target function should use clock #0 (i.e. with start_clock() and
   stop_clock()) to mark which code should be timed.
   
   Stores minimum/maximum time per iteration (in microseconds) in min and max
   (either may be NULL, in which case the value is not stored).
*/
void prof_repeat(double* min, double* max, profile_target_t target, void* arg);


// Timing runs need to last at least this many microseconds to be counted:
#define DURATION_THRESHOLD 200000
// Microseconds per timing run that the prof_repeat function aims for:
#define DURATION_TARGET 300000

#ifdef __cplusplus
 }
#endif
 

#endif // #ifndef FLINT_PROFILER_H
