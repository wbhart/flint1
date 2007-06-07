#define TRIALS 1
#define BITS 1000000UL
//#define SSMUL
//#define GMP
//#define TEST

#ifdef TEST
#define SSMUL
#define GMP
#endif

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "Z_mpn.h"
#include "profiler.h"

/* Runs Z_fast_mul through random data.

TRIALS = number of trials to perform.
LOGLENGTH = log of degree+1.
BITS = number of bits to use in each coefficient. */

unsigned long run_Z_Mul(unsigned long num_trials, unsigned coeff_bits, int fast, unsigned long tweak)
{
   unsigned long i;
   
   // alloc some space
   mpz_t data1;
   mpz_t data2;
   mpz_t data3;
   mpz_t data4;
   
   mpz_init(data1);
   mpz_init(data2);
   mpz_init(data3);
   mpz_init(data4);
   
   gmp_randstate_t state;
   gmp_randinit_default(state);
    
   for (i = 0; i < num_trials; i++)
   {
      // make up random polys
      if (i%20==0)
      {
         mpz_urandomb(data1, state, coeff_bits);
         //mpz_urandomb(data2, state, coeff_bits);
            /*if (gmp_urandomb_ui(state, 1))
               mpz_neg(data1, data1);
            if (gmp_urandomb_ui(state, 1))
               mpz_neg(data2, data2);*/
      }
     

#ifdef SSMUL
      // compute product using Z_mpn_mul
      Z_mul(data3, data1, data1);
#endif
#ifdef GMP
      // compute product using GMP
      mpz_mul(data4, data1, data1);
#endif
#ifdef TEST
      if (mpz_cmp(data3,data4)!=0) printf("Failure!!\n");
#endif

      start_clock(0);
      if (fast)
      {
         Z_mul(data3, data1, data1, tweak);
      } else
      {
         mpz_mul(data4, data1, data1);
      }
      stop_clock(0);
          
   }   
   // clean up
   mpz_clear(data1);
   mpz_clear(data2);
   mpz_clear(data3);
   mpz_clear(data4);
}


int main (int argc, const char * argv[])
{
   double time1, time2;
   unsigned long besti, imax;
   double best;
   unsigned long tweak;
   unsigned long bits;
   
   for (unsigned long words = 1000UL; words < 17100000; words=floor(words*pow(2.0,1.0/64.0))) 
   {
       bits = 64*words;
       printf("%ld words\n",words);
       besti = 0;
       best = 1000.0;
       if (bits/64>200000)
       {
          imax = 3;
          tweak = 1;
       } else
       {
          tweak = 1;
          imax = 5;
       }
       for (unsigned long i = 0; i < imax; i++)
       {
          init_clock(0);
          run_Z_Mul(TRIALS, bits, 1, tweak);
          time1 = get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0;
          if (time1 < best) 
          {
             best = time1;
             besti = tweak;
          }
          tweak <<= 2;
       }
       printf("FLINT = %lf ", besti, best); 
       init_clock(0);
       run_Z_Mul(TRIALS, bits, 0, 0);
       time2 = get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0;
       printf("GMP = %lf ",time2);   
       printf("ratio = %lf, best = %ld\n",time2/best,besti);
       
       //printf("%ld\n",bits);
       //run_Z_Mul(TRIALS, bits, 1);
       
   }

   return 0;
}
