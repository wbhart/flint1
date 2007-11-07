/****************************************************************************

   zmod_poly-profile.c : Profiling code for zmod_poly

   Copyright (C) 2007, David Howden

*****************************************************************************/

#include "profiler-main.h"
#include "zmod_poly.h"
#include "long_extras.h"
#include "flint.h"
#include <string.h>
#include <math.h>

#define PRIME 7 // Prime p to use for Z/pZ

// ============================================================================


void sample_zmod_poly_KS(unsigned long n, void* arg, unsigned long count)
{
   zmod_poly_t poly;
   zmod_poly_init(poly, PRIME);
   zmod_poly_clear(poly);
}


char* profDriverString_zmod_poly_KS(char* params)
{
   return
   "zmod_poly KS mul.\n"
   "Parameters: n_min, n_max, n_skip.\n";
}


char* profDriverDefaultParams_zmod_poly_KS()
{
   return "1 1000 1";
}


void profDriver_zmod_poly_KS(char* params)
{
   unsigned long n_min, n_max, n_skip;
   sscanf(params, "%ld %ld %ld", &n_min, &n_max, &n_skip);

   prof1d_set_sampler(sample_zmod_poly_KS);
   
   for (unsigned long n = n_min; n <= n_max; n += n_skip)
      prof1d_sample(n, NULL);
}

// end of file ****************************************************************
