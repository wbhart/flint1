/****************************************************************************

ZmodFpoly-profile.c

Profiling for ZmodFpoly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "ZmodFpoly.h"
#include "flint.h"
#include <string.h>
#include <math.h>


// ============================================================================


void sample_ZmodFpoly_FFT(unsigned long length, unsigned long n,
                          unsigned long count)
{
   unsigned long m = ceil_log2(2*length);
   
   ZmodFpoly_t poly;
   ZmodFpoly_init(poly, m, n, 1);
   
   // todo: need to generate random data here
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
   {
      poly->length = length;
      ZmodFpoly_FFT(poly, 2*length);
   }

   prof2d_stop();
   
   ZmodFpoly_clear(poly);
}


char* prof2dDriverString_ZmodFpoly_FFT(char* params)
{
   return "ZmodFpoly_FFT over various transform lengths and coefficient sizes";
}


/*
Parameters for this target are:
   length_min: minimum truncation length
   length_max: maximum truncation length
   length_ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_ZmodFpoly_FFT(char* params)
{
   int length_min, length_max, n_count;
   double length_ratio;

   if (strlen(params) == 0)
   {
      // default parameters:
      length_min = 100;
      length_max = 200;
      length_ratio = 1.1;
      n_count = 6;
   }
   else
   {
      sscanf(params, "%d %d %lf %d", &length_min, &length_max,
                                     &length_ratio, &n_count);
   }

   prof2d_set_sampler(sample_ZmodFpoly_FFT);

   for (unsigned long length = length_min; length < length_max;
        length = (int)(ceil(length_ratio * length)))
   {
      unsigned long m = ceil_log2(2*length);

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n <= n_count * n_skip; n += n_skip)
         prof2d_sample(length, n);
   }
}


// ============================================================================

void sample_ZmodFpoly_IFFT(unsigned long length, unsigned long n,
                           unsigned long count)
{
   unsigned long m = ceil_log2(length);
   
   ZmodFpoly_t poly;
   ZmodFpoly_init(poly, m, n, 1);
   poly->length = length;
   
   // todo: need to generate random data here
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodFpoly_IFFT(poly);

   prof2d_stop();
   
   ZmodFpoly_clear(poly);
}


char* prof2dDriverString_ZmodFpoly_IFFT(char* params)
{
   return "ZmodFpoly_IFFT over various transform lengths and coefficient sizes";
}


/*
Parameters for this target are:
   length_min: minimum truncation length
   length_max: maximum truncation length
   length_ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_ZmodFpoly_IFFT(char* params)
{
   int length_min, length_max, n_count;
   double length_ratio;

   if (strlen(params) == 0)
   {
      // default parameters:
      length_min = 100;
      length_max = 200;
      length_ratio = 1.1;
      n_count = 6;
   }
   else
   {
      sscanf(params, "%d %d %lf %d", &length_min, &length_max,
                                     &length_ratio, &n_count);
   }

   prof2d_set_sampler(sample_ZmodFpoly_IFFT);

   for (unsigned long length = length_min; length < length_max;
        length = (int)(ceil(length_ratio * length)))
   {
      unsigned long m = ceil_log2(length);

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n <= n_count * n_skip; n += n_skip)
         prof2d_sample(length, n);
   }
}



// ============================================================================

void sample_ZmodFpoly_negacyclic_convolution(
      unsigned long depth, unsigned long n, unsigned long count)
{
   ZmodFpoly_t poly1, poly2, poly3;
   ZmodFpoly_init(poly1, depth, n, 1);
   ZmodFpoly_init(poly2, depth, n, 1);
   ZmodFpoly_init(poly3, depth, n, 1);

   unsigned long size = 1 << depth;
   for (unsigned long i = 0; i < size; i++)
   {
      profiler_random_limbs(poly1->coeffs[i], n+1);
      profiler_random_limbs(poly2->coeffs[i], n+1);
   }

   unsigned long twist = (2*n*FLINT_BITS) >> depth;
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodFpoly_negacyclic_convolution(poly3, poly1, poly2);

   prof2d_stop();
   
   ZmodFpoly_clear(poly3);
   ZmodFpoly_clear(poly2);
   ZmodFpoly_clear(poly1);
}


char* prof2dDriverString_ZmodFpoly_negacyclic_convolution(char* params)
{
   return "ZmodFpoly_negacyclic_convolution over various depths and coefficient lengths";
}


/*
Parameters for this target are:
   depth_min: minimum depth
   depth_max: maximum depth
   n_min: minimum n to try
   n_max: maximum n to try
*/
void prof2dDriver_ZmodFpoly_negacyclic_convolution(char* params)
{
   int depth_min, depth_max, n_min, n_max;

   if (strlen(params) == 0)
   {
      // default parameters:
      depth_min = 3;
      depth_max = 8;
      n_min = 1;
      n_max = 8;
   }
   else
   {
      sscanf(params, "%d %d %d %d", &depth_min, &depth_max, &n_min, &n_max);
   }

   prof2d_set_sampler(sample_ZmodFpoly_negacyclic_convolution);

   for (unsigned long depth = depth_min; depth <= depth_max; depth++)
   {
      for (unsigned long n = n_min; n <= n_max; n++)
      {
         if ((2*n*FLINT_BITS) % (1 << depth))
            continue;
         
         prof2d_sample(depth, n);
      }
   }
}


// ============================================================================


void sample_mpn_mul_n(
      unsigned long n, unsigned long dummy, unsigned long count)
{
   mp_limb_t* buf;
   buf = malloc(4*n*sizeof(mp_limb_t));

   profiler_random_limbs(buf, 2*n);
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      mpn_mul_n(buf + 2*n, buf, buf+n, n);

   prof2d_stop();
   
   free(buf);
}


char* prof2dDriverString_mpn_mul_n(char* params)
{
   return "mpn_mul_n for various n";
}


/*
Parameters for this target are:
   n_min
   n_max
*/
void prof2dDriver_mpn_mul_n(char* params)
{
   int n_min, n_max;

   if (strlen(params) == 0)
   {
      // default parameters:
      n_min = 30;
      n_max = 200;
   }
   else
   {
      sscanf(params, "%d %d", &n_min, &n_max);
   }

   prof2d_set_sampler(sample_mpn_mul_n);

   for (unsigned long n = n_min; n <= n_max; n++)
      prof2d_sample(n, 0);
}


// end of file ****************************************************************
