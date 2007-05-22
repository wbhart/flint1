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
   poly->length = length;
   
   // todo: need to generate random data here
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodFpoly_FFT(poly, 2*length);

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

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS_PER_LIMB);
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
      unsigned long m = ceil_log2(2*length);

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS_PER_LIMB);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n <= n_count * n_skip; n += n_skip)
         prof2d_sample(length, n);
   }
}


// ============================================================================


void sample_ssfft(unsigned long length, unsigned long n, unsigned long count)
{
   unsigned long m = ceil_log2(2*length);
   
   ZmodFpoly_t poly;
   ZmodFpoly_init(poly, m, n, 1);
   poly->length = length;
   
   // todo: need to generate random data here
   
   prof2d_start();

   for (unsigned long i = 0; i < count; i++)
      ssfft_fft(poly->coeffs, 1, m, length, 2*length, 0,
                (4*n*FLINT_BITS_PER_LIMB) >> m, n, poly->scratch);

   prof2d_stop();
   
   ZmodFpoly_clear(poly);
}


char* prof2dDriverString_ssfft(char* params)
{
   return "legacy ssfft code";
}


/*
Parameters for this target are:
   length_min: minimum truncation length
   length_max: maximum truncation length
   length_ratio: e.g. 1.03 means increase by 3% at a time
   n_count: number of coefficient lengths to test
*/
void prof2dDriver_ssfft(char* params)
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

   prof2d_set_sampler(sample_ssfft);

   for (unsigned long length = length_min; length < length_max;
        length = (int)(ceil(length_ratio * length)))
   {
      unsigned long m = ceil_log2(2*length);

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS_PER_LIMB);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n <= n_count * n_skip; n += n_skip)
         prof2d_sample(length, n);
   }
}


// end of file ****************************************************************
