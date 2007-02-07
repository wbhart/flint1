#define TRIALS 5
#define BITS 1000000UL
//#define SSMUL
//#define GMP
#define TEST

#ifdef TEST
#define SSMUL
#define GMP
#endif

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "Z-ssmul.h"
#include "profiler.h"

/* Runs Z_fast_mul through random data.

TRIALS = number of trials to perform.
LOGLENGTH = log of degree+1.
BITS = number of bits to use in each coefficient. */

void run_Z_Mul(unsigned long num_trials, unsigned coeff_bits, int fast)
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
         mpz_urandomb(data2, state, coeff_bits);
            /*if (gmp_urandomb_ui(state, 1))
               mpz_neg(data1, data1);
            if (gmp_urandomb_ui(state, 1))
               mpz_neg(data2, data2);*/
      }
     
/*
#ifdef SSMUL
      // compute product using SSMul
      Z_fast_mul(data3, data1, data2);
#endif
#ifdef GMP
      // compute product using GMP
      mpz_mul(data4, data1, data2);
#endif
#ifdef TEST
      if (mpz_cmp(data3,data4)!=0) printf("Failure!!\n");
#endif
*/
      if (fast)
      {
         Z_fast_mul(data3, data1, data2);
      } else
      {
         mpz_mul(data4, data1, data2);
      }

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
   for (unsigned long bits = 100000UL; bits < 1000000; bits+=1234) 
   {
       printf("%ld\n",bits);
       init_clock(0);
       start_clock(0);
       run_Z_Mul(TRIALS, bits, 1);
       stop_clock(0);
       time1 = get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0;
       printf("Time fast = %lf\n",time1); 
       init_clock(0);
       start_clock(0);
       run_Z_Mul(TRIALS, bits, 0);
       stop_clock(0);
       time2 = get_clock(0) / FLINT_CLOCK_SCALE_FACTOR / 1800000000.0;
       printf("Time mpz_mul = %lf\n",time2);   
       if (time1 > time2) printf("***********************************\n");
       
       /*printf("%ld\n",bits);
       run_Z_Mul(TRIALS, bits, 1);*/
       
   }

   return 0;
}
