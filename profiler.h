/******************************************************************************

 Timing/profiling

 (C) 2006 William Hart and David Harvey

******************************************************************************/


#include <sys/time.h>
#include <sys/resource.h>


#ifndef FLINT_PROFILER_H
#define FLINT_PROFILER_H

#define FLINT_NUM_CLOCKS 20
#define FLINT_USE_CYCLE_COUNTER 1


extern double clock_last[FLINT_NUM_CLOCKS];
extern double clock_accum[FLINT_NUM_CLOCKS];


// Relative timings on X86 machines running gcc

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))

#define FLINT_HAVE_CYCLE_COUNTER 1

// Make it around order of magnitude of microseconds, assuming about 1GHz cpu:
#define FLINT_CLOCK_SCALE_FACTOR 0.001

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
   for (unsigned long i = 0; i < FLINT_NUM_CLOCKS; i++)
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


// ============================================================================
// Support for repeatedly running a given routine

typedef void (*profile_target_t)(void* arg);

// Timing runs need to last at least this many microseconds to be counted:
#define DURATION_THRESHOLD 200000
// Number of seconds per timing run that the run_profile function aims for:
#define DURATION_TARGET 300000


/*
This function repeatedly calls the function "target" with the parameter
"arg". See profiler.c for the strategy for selecting the number of loops.
It returns the minimum time (per call). If range is non-NULL, it stores the
difference between the maximum and minimum times (per call).
 */
double run_profile(double* range, profile_target_t target, void* arg);


#endif // #ifndef FLINT_PROFILER_H
