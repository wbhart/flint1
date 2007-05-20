/******************************************************************************

 Command-line profiling utility

 (C) 2007 William Hart and David Harvey

******************************************************************************/

#include "profiler.h"


// A function that takes two arguments and an iteration count
// (i.e. the type of function that is getting profiled).
typedef void (*prof2d_exec_t)(unsigned long x, unsigned long y,
                              unsigned long count);


typedef void (*prof2d_main_t)();


void prof2d_exec(unsigned long x, unsigned long y);


/* ============================================================================

   Imported data from the auto-generated table file
   
   (See make-profile-tables.py.)
   
=============================================================================*/

// name of module being profiled.
extern char* prof_module_name;


extern int prof2d_target_count;

extern char* prof2d_target_name[];
extern char* prof2d_target_string[];
extern prof2d_exec_t prof2d_target_exec[];
extern prof2d_main_t prof2d_target_main[];


// end of file ****************************************************************
