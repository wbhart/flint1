/******************************************************************************

 Timing/profiling

 (C) 2007 William Hart and David Harvey

******************************************************************************/

#include "profiler.h"


/*
clock_last[i] is the last read clock value for clock #i.
clock_accum[i] is the total time attributed to clock #i so far.
These should not be read directly; use get_clock(i) instead.
 */
double clock_last[FLINT_NUM_CLOCKS];
double clock_accum[FLINT_NUM_CLOCKS];


// end of file ****************************************************************
