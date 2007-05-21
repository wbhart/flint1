/******************************************************************************

 Command-line profiling utility

 (C) 2007 William Hart and David Harvey

******************************************************************************/

#include "profiler.h"


// A function that takes two arguments and an iteration count
// (i.e. the type of function that is getting profiled).
typedef void (*prof2d_Sampler_t)(unsigned long x, unsigned long y,
                                 unsigned long count);


// A function that runs a bunch of profiles
typedef void (*prof2d_Driver_t)(char* params);


// A function that returns a string (the description of the target)
typedef char* (*prof2d_DriverString_t)(char* params);


void prof2d_set_sampler(prof2d_Sampler_t sampler);
void prof2d_sample(unsigned long x, unsigned long y);
void prof2d_start();
void prof2d_stop();


/* ============================================================================

   Imported data from the auto-generated table file
   
   (See make-profile-tables.py.)
   
=============================================================================*/

// name of module being profiled.
extern char* prof_module_name;


extern int prof2d_target_count;

extern char* prof2d_target_name[];
extern prof2d_Driver_t prof2d_Driver_list[];
extern prof2d_DriverString_t prof2d_DriverString_list[];


// end of file ****************************************************************
