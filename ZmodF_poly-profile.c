/****************************************************************************

ZmodF_poly-profile.c

Profiling for ZmodF_poly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "ZmodF_poly.h"
#include "flint.h"
#include <string.h>
#include <math.h>


// ============================================================================


void sample_ZmodF_poly_FFT(unsigned long length, unsigned long n,
                          unsigned long count)
{
   unsigned long m = ceil_log2(2*length);
   
   ZmodF_poly_t poly;
   ZmodF_poly_init(poly, m, n, 1);
   
   // todo: need to generate random data here
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
   {
      poly->length = length;
      ZmodF_poly_FFT(poly, 2*length);
   }

   prof_stop();
   
   ZmodF_poly_clear(poly);
}


char* profDriverString_ZmodF_poly_FFT(char* params)
{
   return "ZmodF_poly_FFT over various truncation lengths and coefficient sizes.\n"
   "Parameters are: min truncation length; max truncation length; ratio between\n"
   "consecutive truncation lengths; number of coefficient lengths to try.";
}

char* profDriverDefaultParams_ZmodF_poly_FFT()
{
   return "100 200 1.1 6";
}


void profDriver_ZmodF_poly_FFT(char* params)
{
   unsigned long length_min, length_max, n_count;
   double length_ratio;

   sscanf(params, "%ld %ld %lf %ld", &length_min, &length_max,
                                     &length_ratio, &n_count);

   prof2d_set_sampler(sample_ZmodF_poly_FFT);

   for (unsigned long length = length_min; length < length_max;
        length = (int)(ceil(length_ratio * length)))
   {
      unsigned long m = ceil_log2(2*length);

      // restrict coefficient lengths so that appropriate roots of unity
      // are available
      unsigned long n_skip = (1 << m) / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n <= n_count * n_skip; n += n_skip)
         prof2d_sample(length, n);
   }
}


// ============================================================================

void sample_ZmodF_poly_IFFT(unsigned long length, unsigned long n,
                           unsigned long count)
{
   unsigned long m = ceil_log2(length);
   
   ZmodF_poly_t poly;
   ZmodF_poly_init(poly, m, n, 1);
   poly->length = length;
   
   // todo: need to generate random data here
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_poly_IFFT(poly);

   prof_stop();
   
   ZmodF_poly_clear(poly);
}


char* profDriverString_ZmodF_poly_IFFT(char* params)
{
   return "ZmodF_poly_IFFT over various truncation lengths and coefficient sizes.\n"
   "Parameters are: min truncation length; max truncation length; ratio between\n"
   "consecutive truncation lengths; number of coefficient lengths to try.";
}

char* profDriverDefaultParams_ZmodF_poly_IFFT()
{
   return "100 200 1.1 6";
}


void profDriver_ZmodF_poly_IFFT(char* params)
{
   unsigned long length_min, length_max, n_count;
   double length_ratio;

   sscanf(params, "%ld %ld %lf %ld", &length_min, &length_max,
                                     &length_ratio, &n_count);

   prof2d_set_sampler(sample_ZmodF_poly_IFFT);

   for (unsigned long length = length_min; length < length_max;
        length = (int)(ceil(length_ratio * length)))
   {
      unsigned long m = ceil_log2(length);

      // restrict coefficient lengths so that appropriate roots of unity
      // are available
      unsigned long n_skip = (1 << m) / (4*FLINT_BITS);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n <= n_count * n_skip; n += n_skip)
         prof2d_sample(length, n);
   }
}



// ============================================================================

void sample_ZmodF_poly_negacyclic_convolution(
      unsigned long depth, unsigned long n, unsigned long count)
{
   ZmodF_poly_t poly1, poly2, poly3;
   ZmodF_poly_init(poly1, depth, n, 1);
   ZmodF_poly_init(poly2, depth, n, 1);
   ZmodF_poly_init(poly3, depth, n, 1);

   unsigned long size = 1 << depth;
   for (unsigned long i = 0; i < size; i++)
   {
      profiler_random_limbs(poly1->coeffs[i], n+1);
      profiler_random_limbs(poly2->coeffs[i], n+1);
   }

   unsigned long twist = (2*n*FLINT_BITS) >> depth;
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      ZmodF_poly_negacyclic_convolution(poly3, poly1, poly2);

   prof_stop();
   
   ZmodF_poly_clear(poly3);
   ZmodF_poly_clear(poly2);
   ZmodF_poly_clear(poly1);
}


char* profDriverString_ZmodF_poly_negacyclic_convolution(char* params)
{
   return "ZmodF_poly_negacyclic_convolution over various depths and coefficient sizes.\n"
   "Parameters are: min depth; max depth; min coeff length; max coeff length.";
}

char* profDriverDefaultParams_ZmodF_poly_negacyclic_convolution()
{
   return "3 8 1 8";
}


void profDriver_ZmodF_poly_negacyclic_convolution(char* params)
{
   unsigned long depth_min, depth_max, n_min, n_max;

   sscanf(params, "%ld %ld %ld %ld", &depth_min, &depth_max,
                                     &n_min, &n_max);

   prof2d_set_sampler(sample_ZmodF_poly_negacyclic_convolution);

   for (unsigned long depth = depth_min; depth <= depth_max; depth++)
   {
      for (unsigned long n = n_min; n <= n_max; n++)
      {
         // restrict coefficient lengths so that appropriate roots of unity
         // are available
         if ((2*n*FLINT_BITS) % (1 << depth))
            continue;
         
         prof2d_sample(depth, n);
      }
   }
}


// end of file ****************************************************************
