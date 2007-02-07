/*****************************************************************************

 ssmul-profiling.cpp
 profiling for ssmul.c

 Copyright (C) 2006, David Harvey

****************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include <vector>
#include <cmath>
#include <utility>
#include "flint.h"
#include "ssmul.h"
#include "Zvec.h"
#include "profiler.h"
#include "time.h"

#define REAL_TIME 0

void timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
}

// Timing runs need to last at least this many seconds to be counted:
#define DURATION_THRESHOLD 0.2
// Number of seconds per timing run that the test_case function aims for:
#define DURATION_TARGET 0.3

using namespace std;

// *************************************************************************
// profiling SSMul

double ssmul_do_trials(unsigned long degree, unsigned long coeff_bits,
                       unsigned long num_trials)
{
   struct timeval result;
   struct timeval start;
   struct timeval stop;
   
   gmp_randstate_t state;
   gmp_randinit_default(state);

   // alloc space and make up random polys
   mpz_t* data1 = (mpz_t*) malloc((degree + 1) * sizeof(mpz_t));
   mpz_t* data2 = (mpz_t*) malloc((degree + 1) * sizeof(mpz_t));
   mpz_t* data3 = (mpz_t*) malloc((2*degree + 1) * sizeof(mpz_t));
   
   unsigned long log_length = 0;
   while ((1<<log_length) <= degree) log_length++;
   
   double elapsed;
   
   if ((2*coeff_bits+log_length < FLINT_BITS_PER_LIMB) && (degree > 9))
   {   
       data1[0]->_mp_d = (mp_limb_t*) malloc((degree + 1) * sizeof(mp_limb_t));
       data2[0]->_mp_d = (mp_limb_t*) malloc((degree + 1) * sizeof(mp_limb_t));
       data3[0]->_mp_d = (mp_limb_t*) malloc((2*degree + 1) * sizeof(mp_limb_t));
       data1[0]->_mp_alloc=1;
       data1[0]->_mp_size=0;
       data2[0]->_mp_alloc=1;
       data2[0]->_mp_size=0;
       data3[0]->_mp_alloc=1;
       data3[0]->_mp_size=0;
       for (unsigned long i = 1; i <= degree; i++)
       {
           data1[i]->_mp_d = data1[0]->_mp_d+i;
           data1[i]->_mp_alloc=1;
           data1[i]->_mp_size=0;
           data2[i]->_mp_d = data2[0]->_mp_d+i;
           data2[i]->_mp_alloc=1;
           data2[i]->_mp_size=0;
       }
       for (unsigned long i = 1; i <= 2*degree; i++)
       {
           data3[i]->_mp_d = data3[0]->_mp_d+i;
           data3[i]->_mp_alloc=1;
           data3[i]->_mp_size=0;
       }
   } else
   {
      for (unsigned long i = 0; i <= degree; i++)
      { 
          mpz_init(data1[i]);
          mpz_init(data2[i]);
      }
      for (unsigned long i = 0; i <= 2*degree; i++)
      {
         mpz_init(data3[i]);
      }
   }

   for (unsigned long j = 0; j <= degree; j++)
   {
      mpz_urandomb(data1[j], state, coeff_bits);
      /*if (gmp_urandomm_ui(state, 2))
         mpz_neg(data1[j], data1[j]);*/
         
      mpz_urandomb(data2[j], state, coeff_bits);
      /*if (gmp_urandomb_ui(state, 2))
         mpz_neg(data2[j], data2[j]);*/
   }
   Zvec a, b, c;
   a.length = degree+1;
   b.length = degree+1;
   c.length = 2*degree+1;
   a.coords = data1;
   b.coords = data2;
   c.coords = data3;

#if REAL_TIME
   gettimeofday(&start,NULL);
#else
   init_clock(0);
#endif
   
   for (unsigned long j = 0; j < num_trials; j++)
   {
#if !REAL_TIME
      start_clock(0);
#endif
      //SSMul(c, a, b, coeff_bits);
      Zvec_mul(c,a,b);
      //Zvec_karamul(c, a, b, 2*coeff_bits+log_length);
#if !REAL_TIME
      stop_clock(0);
#endif
   }
#if REAL_TIME
   gettimeofday(&stop,NULL);
   
   timeval_subtract(&result,&stop,&start);
   
   elapsed = (double) (result.tv_sec*1000000UL + result.tv_usec)/1000000.0f;
#endif
   //elapsed = ((double) (end - start)) / CLOCKS_PER_SEC;
   
   // clean up
   if ((2*coeff_bits+log_length >= FLINT_BITS_PER_LIMB) || (degree < 9))
   {
      for (unsigned long i = 0; i <= degree; i++)
      {
         mpz_clear(data1[i]);
         mpz_clear(data2[i]);
      }
      for (unsigned long i = 0; i <= 2*degree; i++)
         mpz_clear(data3[i]);
   } else
   {
      free(data1[0]->_mp_d);
      free(data2[0]->_mp_d);
      free(data3[0]->_mp_d);
   }

   free(data1);
   free(data2);
   free(data3);
   
   gmp_randclear(state);

#if REAL_TIME
   return elapsed;
#else
   return get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0;
#endif
}


// Performs several runs of ssmul_do_trials, with varying num_trials.
// It tries to quickly find a good num_trials that will give accurate results.
// It returns the shortest time (per multiplication) and the difference between
// the shortest and longest times.
pair<double, double> ssmul_test_case(unsigned long degree, unsigned long coeff_bits)
{
    // the number of timings that were at least DURATION_THRESHOLD seconds:
    unsigned long good_count = 0;
    double max_time, min_time;

    // first try one loop
    unsigned long num_trials = 1;
    double last_time = ssmul_do_trials(degree, coeff_bits, 1);

    // loop until we have enough good times
    while (1)
    {
       double per_trial = last_time / num_trials;

       // if the last recorded time was long enough, record it
       if (last_time > DURATION_THRESHOLD)
       {
          if (good_count)
          {
             max_time = max(max_time, per_trial);
             min_time = min(min_time, per_trial);
          }
          else
             max_time = min_time = per_trial;

          if (++good_count == 5)
             return make_pair(min_time, max_time - min_time);
       }

       // adjust num_trials so that the elapsed time gravitates towards
       // DURATION_TARGET; num_trials can be changed by a factor of
       // at most 25%, and must be at least 1
       if (last_time < 0.0001)
          last_time = 0.0001;
       double adjust_ratio = DURATION_TARGET / last_time;
       if (adjust_ratio > 1.25)
          adjust_ratio = 1.25;
       if (adjust_ratio < 0.75)
          adjust_ratio = 0.75;
       num_trials = (unsigned long) ceil(adjust_ratio * num_trials);
       // just to be safe:
       if (num_trials == 0)
          num_trials = 1;

       // run another trial
       last_time = ssmul_do_trials(degree, coeff_bits, num_trials);
    }
}

void ssmul_run_profile()
{
   double max_total_bits = 16000000;
   double degree_ratio = 1.03*1.03;
   double coeff_bits_ratio = 1.03*1.03;

   vector<unsigned long> degree_list;
   degree_list.push_back(1);
   for (unsigned long n = 1; n < log(max_total_bits) / log(degree_ratio); n++)
   {
      unsigned long degree = (unsigned long) pow(degree_ratio, (double)n);
      if (degree_list.back() != degree)
         degree_list.push_back(degree);
   }

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
      if (degree < 55) continue;
      //if (degree < 4000) continue;
      for (unsigned long j = 0; j < coeff_bits_list.size(); j++)
      {
         unsigned long coeff_bits = coeff_bits_list[j];
         //if (coeff_bits < 1042000) continue;  // skip these for now
         if (degree * coeff_bits < max_total_bits)
         {
            pair<double, double> result = ssmul_test_case(degree, coeff_bits);
            printf("%ld %ld %lf %lf\n", degree, coeff_bits,
                   log(result.first) / log(10.0), result.second);
            fflush(stdout);
         }
      }
   }
}


// *************************************************************************
// profiling ffts

double fft_do_trials(unsigned long m, unsigned long n, unsigned long num_trials)
{
   unsigned long r = (n * FLINT_BITS_PER_LIMB) >> (m-1);

   gmp_randstate_t state;
   gmp_randinit_default(state);

   // alloc space and make up random polys
   unsigned long M = 1 << m;
   mp_limb_t* data = (mp_limb_t*) malloc((n+1) * (M+1) * sizeof(mp_limb_t));

   for (unsigned long i = 0; i < M*(n+1); i++)
      data[i] = gmp_urandomb_ui(state, FLINT_BITS_PER_LIMB);
   
   mp_limb_t** pointarr = (mp_limb_t**) malloc((M+1)*sizeof(mp_limb_t*));
   for (unsigned long i = 0; i <= M; i++)
      pointarr[i] = data + i*(n+1);

   init_clock(0);
   start_clock(0);

   for (unsigned long j = 0; j < num_trials; j++)
      fft_main(pointarr, 1, 0, r, m, pointarr + M, n, 1, -1);

   stop_clock(0);

   free(data);
   free(pointarr);

   return get_clock(0) / 1000000;
}


// Performs several runs of fft_do_trials, with varying num_trials.
// It tries to quickly find a good num_trials that will give accurate results.
// It returns the average time (per multiplication) and the difference between
// the shortest and longest times.
pair<double, double> fft_test_case(unsigned long m, unsigned long n)
{
   vector<pair<unsigned long, double> > times;
   // the number of timings that were at least 0.02 seconds:
   unsigned long good_count = 0;
   
   // first try one loop
   unsigned long num_trials = 1;
   double result = fft_do_trials(m, n, 1);
   times.push_back(make_pair(num_trials, result));
   if (result > 0.02)
      good_count++;
   
   // now keep doubling until we get at least 0.02 seconds
   while (times.back().second < 0.02)
   {
      num_trials <<= 1;
      result = fft_do_trials(m, n, num_trials);
      times.push_back(make_pair(num_trials, result));
      if (result > 0.02)
         good_count++;
   }
   
   // now we've got enough data to make a decent guess at a good value
   // for num_trials
   while (good_count < 5)
   {
      num_trials = (unsigned long)
                         ceil(0.025 * num_trials / times.back().second);
      result = fft_do_trials(m, n, num_trials);
      times.push_back(make_pair(num_trials, result));
      if (result > 0.02)
         good_count++;
   }
   
   // extract the "good" times
   vector<double> good_times;
   for (vector<pair<unsigned long, double> >::iterator ptr = times.begin();
        ptr != times.end(); ptr++)
   {
      if (ptr->second > 0.02)
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


void fft_run_profile()
{
   double max_total_bits = 40000000;
   double coeff_bits_ratio = 1.05;
   double degree_ratio = 1.05;

#if 1
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
      degree_list.push_back((1 << n) - 1);
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
         unsigned long n = (unsigned long) ceil(1.0 * coeff_bits / FLINT_BITS_PER_LIMB);
         unsigned long M = 1 << m;

         unsigned long min_nB = M / 2;
         // round up min_nB to a multiple of B
         unsigned long min_n = (unsigned long)
                                   ceil(1.0 * min_nB / FLINT_BITS_PER_LIMB);
         // stay above the diagonal
         if (n < min_n)
            continue;
         
         // round n up to a multiple of min_n
         n = ((unsigned long) ceil(1.0 * n / min_n)) * min_n;
         
         pair<double, double> result = fft_test_case(m, n);
         printf("%ld %ld %lf %lf\n", degree, coeff_bits,
                log(result.first) / log(10.0), result.second);
         fflush(stdout);
      }
   }
}


int main()
{
 //  ssmul_do_trials(16384, 5, 1);
   ssmul_run_profile();
//   fft_run_profile();
}


// end of file ***************************************************************
