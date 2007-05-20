/****************************************************************************

ZmodFpoly-profile.c

Profiling for ZmodFpoly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "ZmodFpoly.h"
#include "flint.h"


// ============================================================================


char prof2dString_ZmodFpoly_FFT[] =
   "ZmodFpoly_FFT over various transform lengths and coefficient sizes";


void prof2dMain_ZmodFpoly_FFT()
{
   for (unsigned long length = 30; length < 50; length++)
   {
      unsigned long m = ceil_log2(2*length);

      unsigned long n_skip = (1 << m) / (4*FLINT_BITS_PER_LIMB);
      if (n_skip == 0)
         n_skip = 1;

      for (unsigned long n = n_skip; n < 10; n += n_skip)
      {
         prof2d_exec(length, n);
      }
   }
}


void prof2dExec_ZmodFpoly_FFT(unsigned long length, unsigned long n,
                              unsigned long count)
{
   unsigned long m = ceil_log2(2*length);
   
   ZmodFpoly_t poly;
   ZmodFpoly_init(poly, m, n, 1);
   poly->length = length;
   
   // todo: need to generate random data here
   
   for (unsigned long i = 0; i < count; i++)
      ZmodFpoly_FFT(poly, 2*length);
   
   ZmodFpoly_clear(poly);
}


// ============================================================================


char prof2dString_stupidfunc[] =
   "stupidfunc for various x and y";


void prof2dExec_stupidfunc(unsigned long x, unsigned long y, unsigned long count)
{
   for (unsigned long i = 0; i < count; i++)
      for (unsigned long j = 0; j < x; j++)
         for (unsigned long k = 0; k < y; k++)
         { }
}


// end of file ****************************************************************
