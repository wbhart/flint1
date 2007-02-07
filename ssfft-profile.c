/*****************************************************************************

 ssfft-profile.cpp
 profiling for ssfft.c

 Copyright (C) 2006, David Harvey

****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "ssfft.h"
#include "profiler.h"


gmp_randstate_t ssfft_profile_randstate;


typedef struct fft_arg_t
{
   unsigned long m, z, g;
   int e;
   unsigned long ru, rU;
   unsigned long n;
   mp_limb_t* buf;
   mp_limb_t** buf_ptrs;
   mp_limb_t** scratch;
} fft_arg_t;


/*
Initialises m, n, buf, buf_ptrs, scratch.
Users needs to initialise z, g, e, ru, rU.
 */
void fft_arg_init(fft_arg_t* arg, unsigned long m, unsigned long n)
{
   arg->m = m;
   unsigned long M = 1 << m;
   arg->n = n;
   arg->buf = (mp_limb_t*) malloc((n+1) * (M+1) * sizeof(mp_limb_t));
   for (unsigned long i = 0; i < (n+1)*(M+1); i++)
      arg->buf[i] = gmp_urandomb_ui(ssfft_profile_randstate,
                                    FLINT_BITS_PER_LIMB);
   arg->buf_ptrs = (mp_limb_t**) malloc((M+1) * sizeof(mp_limb_t*));
   for (unsigned long i = 0; i <= M; i++)
      arg->buf_ptrs[i] = arg->buf + (n+1) * i;
   arg->scratch = arg->buf_ptrs + M;
}

void fft_arg_clear(fft_arg_t* arg)
{
   free(arg->buf_ptrs);
   free(arg->buf);
}

void run_ifft_size4_bits(void* x)
{
   fft_arg_t* arg = (fft_arg_t*) x;
   ssfft_ifft_size4_bits(arg->buf_ptrs, 1, arg->z, arg->g, arg->e,
                         arg->ru, arg->rU, arg->n, arg->scratch);
}


void run_threaded_fft()
{
   unsigned long m = 6;
   unsigned long M = 1 << m;
   unsigned long i;
   
   fft_arg_t arg;
   fft_arg_init(&arg, m, (1 << 20) / FLINT_BITS_PER_LIMB);
   arg.z = 3*M/8;
   arg.g = 3*M/4;
   arg.e = 0;
   arg.ru = 0;
   arg.rU = 4 * arg.n * FLINT_BITS_PER_LIMB / M;

   for (i = 0; i < 400; i++)
#if 1
      // threaded version
      ssfft_fft_threaded(arg.buf_ptrs, 1, m, arg.z, arg.g, arg.ru, arg.rU,
                         arg.n);
#else
      // non-threaded version
      ssfft_fft(arg.buf_ptrs, 1, m, arg.z, arg.g, arg.ru, arg.rU,
                         arg.n, arg.scratch);
#endif
}



int main()
{
   gmp_randinit_default(ssfft_profile_randstate);

   run_threaded_fft();

   return 0;


   fft_arg_t arg;
   fft_arg_init(&arg, 2, 8);
   arg.z = 4;
   arg.g = 4;
   arg.e = 0;
   arg.rU = 2 * arg.n * FLINT_BITS_PER_LIMB / 4;

   for (arg.ru = 0; arg.ru < 2; arg.ru++)
   {
      printf("\n------------------------------\n");
      printf("M = 4, z = %ld, g = %ld, ru = %ld\n\n", arg.z, arg.g, arg.ru);
      
      for (unsigned long i = 0; i < 10; i++)
         printf("%lf\n", run_profile(NULL, run_ifft_size4_bits, &arg));
   }

   fft_arg_clear(&arg);
   return 0;
}



// ------------------------------------------------------------
// old profiling code, disabled for now

#if 0

double do_trials(unsigned long m, unsigned long n, unsigned long z,
                 unsigned long g, unsigned long num_trials)
{
   unsigned long r = (4 * n * FLINT_BITS_PER_LIMB) >> m;

   gmp_randstate_t state;
   gmp_randinit_default(state);

   // alloc space and make up random data
   unsigned long M = 1 << m;
   mp_limb_t* data = (mp_limb_t*) malloc((n+1) * (M+1) * sizeof(mp_limb_t));
   mp_limb_t** data_ptrs = (mp_limb_t**) malloc((M+1) * sizeof(mp_limb_t*));
   for (unsigned long i = 0; i <= M; i++)
      data_ptrs[i] = data + i*(n+1);
   
   for (unsigned long i = 0; i < (M+1)*(n+1); i++)
      data[i] = gmp_urandomb_ui(state, FLINT_BITS_PER_LIMB);
   
   init_clock(0);
   start_clock(0);

   for (unsigned long j = 0; j < num_trials; j++)
      ssfft_fft(data_ptrs, 1, m, z, g, 0, r, n, data_ptrs + M);

   stop_clock(0);

   // clean up
   free(data);
   free(data_ptrs);
   
   return get_clock(0) / 1000000;
}


// Performs several runs of do_trials, with varying num_trials.
// It tries to quickly find a good num_trials that will give accurate results.
// It returns the average time (per multiplication) and the difference between
// the shortest and longest times.
pair<double, double> test_case(unsigned long m, unsigned long n,
                               unsigned long z, unsigned long g)
{
   vector<pair<unsigned long, double> > times;
   // the number of timings that were at least 0.02 seconds:
   unsigned long good_count = 0;
   
   // first try one loop
   unsigned long num_trials = 1;
   double result = do_trials(m, n, z, g, 1);
   times.push_back(make_pair(num_trials, result));
   if (result > 0.8)
      good_count++;
   
   // now keep doubling until we get at least 0.02 seconds
   while (times.back().second < 0.8)
   {
      num_trials <<= 1;
      result = do_trials(m, n, z, g, num_trials);
      times.push_back(make_pair(num_trials, result));
      if (result > 0.8)
         good_count++;
   }
   
   // now we've got enough data to make a decent guess at a good value
   // for num_trials
   while (good_count < 5)
   {
      num_trials = (unsigned long)
                         ceil(1.0 * num_trials / times.back().second);
      result = do_trials(m, n, z, g, num_trials);
      times.push_back(make_pair(num_trials, result));
      if (result > 0.8)
         good_count++;
   }
   
   // extract the "good" times
   vector<double> good_times;
   for (vector<pair<unsigned long, double> >::iterator ptr = times.begin();
        ptr != times.end(); ptr++)
   {
      if (ptr->second > 0.8)
         good_times.push_back(ptr->second / ptr->first);
   }
   
   // get max/min/average
   double avg_time = 0.0;
   double max_time = good_times.front();
   double min_time = good_times.back();
   for (vector<double>::iterator ptr = good_times.begin();
        ptr != good_times.end(); ptr++)
   {
      max_time = max(max_time, *ptr);
      min_time = min(min_time, *ptr);
      avg_time += *ptr;
   }
   avg_time /= good_count;
   
   return make_pair(avg_time, max_time - min_time);
}


int main()
{
   double max_total_bits = 40000000;
   double coeff_bits_ratio = 1.05;
   double degree_ratio = 1.05;

#if 1
   // speed data generation up a bit by skipping lots of points
   coeff_bits_ratio *= coeff_bits_ratio;
   degree_ratio *= degree_ratio;
#endif
#if 1
   // skip even more points
   coeff_bits_ratio *= coeff_bits_ratio;
   degree_ratio *= degree_ratio;
#endif

#if 0
   // this version generates nicely spread out degrees
   vector<unsigned long> degree_list;
   degree_list.push_back(1);
   for (unsigned long n = 1; n < log(max_total_bits) / log(degree_ratio); n++)
   {
      unsigned long degree = (unsigned long) pow(degree_ratio, (double)n);
      if (degree_list.back() != degree)
         degree_list.push_back(degree);
   }
#else
   // this version generates degrees that imply transforms whose length
   // is just under a power of two
   vector<unsigned long> degree_list;
   for (unsigned long n = 1; n < log(max_total_bits) / log(2.0); n++)
   {
      degree_list.push_back((1 << n) - 1);
//      degree_list.push_back((1 << n) - 2);
//      if ((1 << n) >= 3)
//         degree_list.push_back((1 << n) - 3);
   }
#endif

   vector<unsigned long> coeff_bits_list;
   coeff_bits_list.push_back(1);
   for (unsigned long n = 1; n < log(max_total_bits) / log(coeff_bits_ratio); n++)
   {
      unsigned long coeff_bits = (unsigned long) pow(coeff_bits_ratio, (double)n);
      if (coeff_bits_list.back() != coeff_bits)
         coeff_bits_list.push_back(coeff_bits);
   }

   for (unsigned long i = 0; i < degree_list.size(); i++)
   {
      unsigned long degree = degree_list[i];
      if (degree < 3) continue;
      
      for (unsigned long j = 0; j < coeff_bits_list.size(); j++)
      {
         unsigned long coeff_bits = coeff_bits_list[j];
         if (degree * coeff_bits >= max_total_bits)
            continue;

         unsigned long m = (unsigned long) ceil(log((double)(2*degree+1)) / log(2.0));
         unsigned long g = 2*degree+1;
         unsigned long z = degree+1;
         unsigned long n = (unsigned long) ceil(1.0 * coeff_bits / FLINT_BITS_PER_LIMB);
         unsigned long M = 1 << m;

#if 1
         // don't allow sqrt2
         unsigned long min_nB = M / 2;
#else
         // allow sqrt2
         unsigned long min_nB = M / 4;
#endif

         // round up min_nB to a multiple of B
         unsigned long min_n = (unsigned long)
                                   ceil(1.0 * min_nB / FLINT_BITS_PER_LIMB);

         // stay above the diagonal
         if (n < min_n)
            continue;
         
         // round n up to a multiple of min_n
         n = ((unsigned long) ceil(1.0 * n / min_n)) * min_n;

         pair<double, double> result = test_case(m, n, z, g);
         printf("%ld %ld %lf %lf\n", degree, coeff_bits,
                log(result.first) / log(10.0), result.second);
         fflush(stdout);
      }
   }
}
#endif

// end of file ***************************************************************
