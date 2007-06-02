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


extern unsigned long min_negacyclic_depth;
extern unsigned long negacyclic_threshold_table[];


void sample_ZmodF_mul(unsigned long n, unsigned long depth,
                      unsigned long count)
{
   // force the ZmodF_mul_precomp machinery to use requested transform depth
   // (depth == 0 means use default, depth == 1 means use mpn_mul_n)
   unsigned long save0 = negacyclic_threshold_table[0];
   unsigned long save1 = negacyclic_threshold_table[1];
   unsigned long save_min = min_negacyclic_depth;

   // todo: the above thing is really hackish, needs to be fixed

   if (depth == 1)
   {
      negacyclic_threshold_table[0] = n+1;
   }
   else if (depth >= 2)
   {
      negacyclic_threshold_table[0] = n;
      negacyclic_threshold_table[1] = n+1;
      min_negacyclic_depth = depth;
   }

   ZmodF_mul_precomp_t info;
   ZmodF_mul_precomp_init(info, n, 0);
   
   mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x2 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
   mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));

   profiler_random_limbs(x1, n);
   x1[n] = 0;
   profiler_random_limbs(x2, n);
   x2[n] = 0;
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_mul_precomp(info, x3, x1, x2);

   prof2d_stop();

   ZmodF_mul_precomp_clear(info);
   free(x3);
   free(x2);
   free(x1);

   negacyclic_threshold_table[0] = save0;
   negacyclic_threshold_table[1] = save1;
   min_negacyclic_depth = save_min;
}


char* prof2dDriverString_ZmodF_mul(char* params)
{
   return "ZmodF_mul over various n and negacyclic transform depths";
}


/*
Parameters for this target are:
   n_min, n_max: minimum and maximum n
   n_skip
   depth: transform depth (0 means use default, 1 means use mpn_mul_n)
*/
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


// end of file ****************************************************************
