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


void sample_ZmodF_mul(unsigned long n, unsigned long depth,
                      unsigned long count)
{
   ZmodF_mul_info_t info;
   
   if (depth == 0)
      ZmodF_mul_info_init(info, n, 0);
   else if (depth == 1)
      ZmodF_mul_info_init_plain(info, n);
   else if (depth == 2)
      ZmodF_mul_info_init_threeway(info, n);
   else
      ZmodF_mul_info_init_negacyclic(info, n, depth);

   mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x2 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));

   profiler_random_limbs(x1, n);
   x1[n] = 0;
   profiler_random_limbs(x2, n);
   x2[n] = 0;
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_mul_info_mul(info, x3, x1, x2);

   prof2d_stop();

   ZmodF_mul_info_clear(info);
   free(x3);
   free(x2);
   free(x1);
}


char* prof2dDriverString_ZmodF_mul(char* params)
{
   return
   "ZmodF_sqr over various n and multiplication algorithms.\n"
   "Input parameters are n_min, n_max, n_skip, depth. Depth == 0 means\n"
   "select algorithm automatically. Depth == 1 means use plain algorithm.\n"
   "Depth == 2 means use threeway algorithm. Otherwise negacyclic algorithm\n"
   "is used with indicated depth.\n"
   "Output fields are n and depth.\n";
}


void prof2dDriver_ZmodF_mul(char* params)
{
   int n_min, n_max, n_skip, depth;

   if (strlen(params) == 0)
   {
      // default parameters:
      n_min = 20;
      n_max = 1000;
      n_skip = 1;
      depth = 4;
   }
   else
   {
      sscanf(params, "%d %d %d %d", &n_min, &n_max, &n_skip, &depth);
   }

   prof2d_set_sampler(sample_ZmodF_mul);

   for (unsigned long n = n_min; n <= n_max; n += n_skip)
      prof2d_sample(n, depth);
}




// ============================================================================


void sample_ZmodF_sqr(unsigned long n, unsigned long depth,
                      unsigned long count)
{
   ZmodF_mul_info_t info;
   
   if (depth == 0)
      ZmodF_mul_info_init(info, n, 1);
   else if (depth == 1)
      ZmodF_mul_info_init_plain(info, n);
   else if (depth == 2)
      ZmodF_mul_info_init_threeway(info, n);
   else
      ZmodF_mul_info_init_negacyclic(info, n, depth);

   mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));

   profiler_random_limbs(x1, n);
   x1[n] = 0;
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_mul_info_sqr(info, x3, x1);

   prof2d_stop();

   ZmodF_mul_info_clear(info);
   free(x3);
   free(x1);
}


char* prof2dDriverString_ZmodF_sqr(char* params)
{
   return
   "ZmodF_sqr over various n and multiplication algorithms.\n"
   "Input parameters are n_min, n_max, n_skip, depth. Depth == 0 means\n"
   "select algorithm automatically. Depth == 1 means use plain algorithm.\n"
   "Depth == 2 means use threeway algorithm. Otherwise negacyclic algorithm\n"
   "is used with indicated depth.\n"
   "Output fields are n and depth.\n";
}


void prof2dDriver_ZmodF_sqr(char* params)
{
   int n_min, n_max, n_skip, depth;

   if (strlen(params) == 0)
   {
      // default parameters:
      n_min = 20;
      n_max = 1000;
      n_skip = 1;
      depth = 4;
   }
   else
   {
      sscanf(params, "%d %d %d %d", &n_min, &n_max, &n_skip, &depth);
   }

   prof2d_set_sampler(sample_ZmodF_sqr);

   for (unsigned long n = n_min; n <= n_max; n += n_skip)
      prof2d_sample(n, depth);
}


// end of file ****************************************************************
