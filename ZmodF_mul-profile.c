/****************************************************************************

ZmodF_mul-profile.c

Profiling for ZmodF_mul

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "ZmodF_mul.h"
#include "flint.h"
#include <string.h>
#include <math.h>


// ============================================================================


// yuck, need to make this a global, something wrong with the design...
unsigned long ZmodF_mul_depth = 0;


void sample_ZmodF_mul(unsigned long n, void* arg, unsigned long count)
{
   ZmodF_mul_info_t info;
   
   // this function assumes n is legal for the requested algorithm
   
   if (ZmodF_mul_depth == 0)
      ZmodF_mul_info_init(info, n, 0);
   else if (ZmodF_mul_depth == 1)
      ZmodF_mul_info_init_plain(info, n, 0);
   else if (ZmodF_mul_depth == 2)
      ZmodF_mul_info_init_threeway(info, n, 0);
   else
      ZmodF_mul_info_init_negacyclic(info, n, ZmodF_mul_depth, 0);

   mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x2 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));

   profiler_random_limbs(x1, n);
   x1[n] = 0;
   profiler_random_limbs(x2, n);
   x2[n] = 0;
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_mul_info_mul(info, x3, x1, x2);

   prof_stop();

   ZmodF_mul_info_clear(info);
   free(x3);
   free(x2);
   free(x1);
}


char* profDriverString_ZmodF_mul(char* params)
{
   return
   "ZmodF_mul over various n and multiplication algorithms.\n"
   "Parameters: n_min, n_max, n_skip, depth. Depth == 0 means\n"
   "select algorithm automatically. Depth == 1 means use plain algorithm.\n"
   "Depth == 2 means use threeway algorithm. Otherwise negacyclic algorithm\n"
   "is used with indicated depth.\n";
}


char* profDriverDefaultParams_ZmodF_mul()
{
   return "20 1000 1 4";
}


void profDriver_ZmodF_mul(char* params)
{
   unsigned long n_min, n_max, n_skip, depth;

   sscanf(params, "%ld %ld %ld %ld", &n_min, &n_max, &n_skip, &depth);

   prof1d_set_sampler(sample_ZmodF_mul);
   
   ZmodF_mul_depth = depth;
   for (unsigned long n = n_min; n <= n_max; n += n_skip)
   {
      // ensure n is legal for given algorithm
      if ((depth == 2) && (n % 3))
         continue;
      if (depth > 2 && ((n * FLINT_BITS) & ((1 << depth) - 1)))
         continue;
      
      prof1d_sample(n, NULL);
   }
}




// ============================================================================


void sample_ZmodF_sqr(unsigned long n, void* arg, unsigned long count)
{
   ZmodF_mul_info_t info;
   
   // this function assumes n is legal for the requested algorithm

   if (ZmodF_mul_depth == 0)
      ZmodF_mul_info_init(info, n, 1);
   else if (ZmodF_mul_depth == 1)
      ZmodF_mul_info_init_plain(info, n, 1);
   else if (ZmodF_mul_depth == 2)
      ZmodF_mul_info_init_threeway(info, n, 1);
   else
      ZmodF_mul_info_init_negacyclic(info, n, ZmodF_mul_depth, 1);

   mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));

   profiler_random_limbs(x1, n);
   x1[n] = 0;
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_mul_info_sqr(info, x3, x1);

   prof_stop();

   ZmodF_mul_info_clear(info);
   free(x3);
   free(x1);
}


char* profDriverString_ZmodF_sqr(char* params)
{
   return
   "ZmodF_sqr over various n and multiplication algorithms.\n"
   "Parameters: n_min, n_max, n_skip, depth. Depth == 0 means\n"
   "select algorithm automatically. Depth == 1 means use plain algorithm.\n"
   "Depth == 2 means use threeway algorithm. Otherwise negacyclic algorithm\n"
   "is used with indicated depth.\n";
}


char* profDriverDefaultParams_ZmodF_sqr()
{
   return "20 1000 1 4";
}


void profDriver_ZmodF_sqr(char* params)
{
   unsigned long n_min, n_max, n_skip, depth;

   sscanf(params, "%ld %ld %ld %ld", &n_min, &n_max, &n_skip, &depth);

   prof1d_set_sampler(sample_ZmodF_sqr);
   
   ZmodF_mul_depth = depth;
   for (unsigned long n = n_min; n <= n_max; n += n_skip)
   {
      // ensure n is legal for given algorithm
      if ((depth == 2) && (n % 3))
         continue;
      if (depth > 2 && ((n * FLINT_BITS) & ((1 << depth) - 1)))
         continue;
      
      prof1d_sample(n, NULL);
   }
}


// end of file ****************************************************************
