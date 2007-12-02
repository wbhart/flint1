/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/*
   ZmodF_mul-tune
   
   Program for tuning the ZmodF_mul module.
   
   This program writes to standard output an automatically tuned version of
   ZmodF_mul-tuning.c.
   
   (If DEBUG is set, it also writes logging info to standard error.)
   
   (C) 2007 David Harvey and William Hart
*/

#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "test-support.h"
#include "profiler.h"
#include "ZmodF_mul.h"
#include "ZmodF_mul-tuning.h"


#define DEBUG 1


typedef struct {
   int algo;
   int squaring;
   unsigned long n;
   unsigned long depth, m, k;
} sample_info_t;



// arg should point to a sample_info_t
void sample_mul(void* arg, unsigned long count)
{
   ZmodF_mul_info_t info;
   sample_info_t* z = (sample_info_t*) arg;
   
   switch (z->algo)
   {
      case ZMODF_MUL_ALGO_PLAIN:
         ZmodF_mul_info_init_plain(info, z->n, z->squaring);
         break;
      
      case ZMODF_MUL_ALGO_THREEWAY:
         ZmodF_mul_info_init_threeway(info, z->n, z->squaring);
         break;

      case ZMODF_MUL_ALGO_FFT:
         ZmodF_mul_info_init_fft(info, z->n, z->depth, z->m, z->k, z->squaring);
         break;
   }

   mp_limb_t* in1 = (mp_limb_t*) flint_stack_alloc(z->n + 1);
   mp_limb_t* in2 = (mp_limb_t*) flint_stack_alloc(z->n + 1);
   mp_limb_t* out = (mp_limb_t*) flint_stack_alloc(z->n + 1);

   urandom_limbs(in1, z->n + 1);
   urandom_limbs(in2, z->n + 1);
   
   if (z->squaring)
      in2 = in1;
   
   // warm up
   for (unsigned long i = 0; i < count/4; i++)
      ZmodF_mul_info_mul(info, out, in1, in2);
      
   // time it
   start_clock(0);
   for (unsigned long i = 0; i < count; i++)
      ZmodF_mul_info_mul(info, out, in1, in2);
   stop_clock(0);
   
   flint_stack_release();
   flint_stack_release();
   flint_stack_release();

   ZmodF_mul_info_clear(info);
}



/*
   Compares two ZmodF_mul algorithms for a specific n.

   algo1/algo2 are:
   0 for plain algorithm
   1 for threeway algorithm
   > 2 indicates FFT of given depth

   Returns nonzero if algo2 is more efficient than algo1 for given n.

   (n will be rounded up automatically to satisfy whatever divisibility
   conditions are required by the requested algorithms.)
*/
int algo_compare(unsigned long n, unsigned long squaring,
                 unsigned long algo1, unsigned long algo2,
                 FILE* f)
{
   sample_info_t info1, info2;

   info1.squaring = info2.squaring = squaring;
   info1.n = info2.n = n;

   if (algo1 == 0)
      info1.algo = ZMODF_MUL_ALGO_PLAIN;
   else if (algo1 == 1)
      info1.algo = ZMODF_MUL_ALGO_THREEWAY;
   else
   {
      info1.algo = ZMODF_MUL_ALGO_FFT;
      info1.depth = algo1;
      info1.m = info1.k = 0;
   }
   
   if (algo2 == 0)
      info2.algo = ZMODF_MUL_ALGO_PLAIN;
   else if (algo2 == 1)
      info2.algo = ZMODF_MUL_ALGO_THREEWAY;
   else
   {
      info2.algo = ZMODF_MUL_ALGO_FFT;
      info2.depth = algo2;
      info2.m = info1.k = 0;
   }

   // round up n appropriately
   unsigned long round = 1;
   if (algo1 == 1 || algo2 == 1)
      round = 3;
   if (algo1 > FLINT_LG_BITS_PER_LIMB || algo2 > FLINT_LG_BITS_PER_LIMB)
      round <<= (FLINT_MAX(algo1, algo2) - FLINT_LG_BITS_PER_LIMB);
   n = (((n-1) / round) + 1) * round;
   
   double time1, time2;

   prof_repeat(&time1, NULL, sample_mul, &info1);
   prof_repeat(&time2, NULL, sample_mul, &info2);
   
#if DEBUG
   fprintf(f, "n = %ld, ", n);

   if (algo1 == 0)
      fprintf(f, "plain");
   else if (algo1 == 1)
      fprintf(f, "threeway");
   else
      fprintf(f, "FFT %ld", algo1);
      
   fprintf(f, " vs ");

   if (algo2 == 0)
      fprintf(f, "plain");
   else if (algo2 == 1)
      fprintf(f, "threeway");
   else
      fprintf(f, "FFT %ld", algo2);
   
   if (time2 < time1)
      fprintf(f, ", 2nd wins");
   else
      fprintf(f, ", 1st wins");

   fprintf(f, " (%lf vs %lf)\n", time1, time2);
#endif

   return time2 < time1;
}



/*
Finds crossover value of n to get from algo1 to algo2.
If start != 0, then it's a starting estimate.
*/
unsigned long algo_threshold(unsigned long algo1, unsigned long algo2,
                             unsigned long squaring, unsigned long start,
                             FILE* f)
{
   // find upper bound
   unsigned long hi = start ? start : 100;
   while (!algo_compare(hi, squaring, algo1, algo2, f))
      hi *= 2;

   hi *= 2;

#if DEBUG
   fprintf(f, "upper bound is %ld\n\n", hi);
#endif

   // find lower bound
   unsigned long lo = hi / 2;
   while (algo_compare(lo, squaring, algo1, algo2, f))
      lo /= 2;

   lo /= 2;

#if DEBUG
   fprintf(f, "lower bound is %ld\n\n", lo);
#endif

   // shrink interval until we reach tolerance of 10%
   while (hi > 1.1 * lo)
   {
      unsigned long mid = (unsigned long) sqrt(1.0 * hi * lo);
      double range = 1.0 * hi / lo;
      if (algo_compare(mid, squaring, algo1, algo2, f))
      {
         lo = (unsigned long) (pow(range, -0.15) * lo);
         hi = (unsigned long) (pow(range, -0.3) * hi);
      }
      else
      {
         lo = (unsigned long) (pow(range, 0.3) * lo);
         hi = (unsigned long) (pow(range, 0.15) * hi);
      }
#if DEBUG
      fprintf(f, "interval is [%ld, %ld], ratio = %lf\n",
              lo, hi, 1.0 * hi / lo);
#endif
   }
   
   return (unsigned long) sqrt(1.0 * hi * lo);
}


int main(int argc, char* argv[])
{
   FILE* fout = stdout;
   FILE* flog = stderr;

   test_support_init();

   fprintf(fout, "/*\n");
   fprintf(fout, "   Tuning values for ZmodF_mul module\n");
   fprintf(fout, "\n");
   fprintf(fout, "   Automatically generated by ZmodF_mul-tune program\n");
   fprintf(fout, "*/\n\n");
   fprintf(fout, "#include \"ZmodF_mul-tuning.h\"\n");
   fprintf(fout, "#include \"ZmodF_mul.h\"\n");
   fprintf(fout, "\n");
   fflush(fout);

   for (int squaring = 0; squaring <= 1; squaring++)
   {
      char* type = squaring ? "sqr" : "mul";
      
      // plain/threeway threshold
      unsigned long n;
      for (n = 3; algo_compare(n, squaring, 1, 0, flog); n += 3);
      fprintf(fout, "unsigned long ZmodF_%s_plain_threeway_threshold = %ld;\n",
              type, n);
      fflush(fout);

      if (!squaring)
         ZmodF_mul_plain_threeway_threshold = n;
      else
         ZmodF_sqr_plain_threeway_threshold = n;

      // plain/fft threshold
      n = algo_threshold(0, 3, squaring, 0, flog);
      fprintf(fout, "unsigned long ZmodF_%s_plain_fft_threshold = %ld;\n",
              type, n);
      fflush(fout);

      if (!squaring)
         ZmodF_mul_plain_fft_threshold = n;
      else
         ZmodF_sqr_plain_fft_threshold = n;

      // threeway/fft threshold
      n = algo_threshold(1, 4, squaring, 0, flog);
      fprintf(fout, "unsigned long ZmodF_%s_threeway_fft_threshold = %ld;\n",
              type, n);
      fflush(fout);

      if (!squaring)
         ZmodF_mul_threeway_fft_threshold = n;
      else
         ZmodF_sqr_threeway_fft_threshold = n;

      // fft thresholds between different depths
      fprintf(fout, "unsigned long ZmodF_%s_fft_table[20] =\n   {", type);
      unsigned long depth;
      for (depth = 3; depth < 10; depth++)
      {
         n = algo_threshold(depth, depth+1, squaring, 0, flog);
         if (!squaring)
            ZmodF_mul_fft_table[depth - 3] = n;
         else
            ZmodF_sqr_fft_table[depth - 3] = n;

         fprintf(fout, "%ld, ", n);
         fflush(fout);
      }

      fprintf(fout, "0};\n\n");
      if (!squaring)
         ZmodF_mul_fft_table[depth - 3] = 0;
      else
         ZmodF_sqr_fft_table[depth - 3] = 0;
   }

   fprintf(fout, "\n");
   fprintf(fout, "// end of file *********************************\n");

   test_support_cleanup();
   return 0;
}



// end of file ****************************************************************
