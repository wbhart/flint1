#define TRIALS 1

#define TEST
//#define TIMING


#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "flint.h"
#include "Z_mpn.h"
#include "profiler.h"

/* Runs Z_fast_mul through random data.

num_trials = number of trials to perform.
coeff_bits = number of bits to use in each coefficient. */

unsigned long run_Z_Mul(unsigned long num_trials, unsigned long coeff_bits, int fast, unsigned long tweak)
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
     

#ifdef TEST
      // compute product using Z_mpn_mul
      Z_mul(data3, data1, data1);
      // compute product using GMP
      mpz_mul(data4, data1, data1);
      
      if (mpz_size(data3) != mpz_size(data4)) printf("Huh! %ld, %ld\n", mpz_size(data3), mpz_size(data4));

      if (mpz_cmp(data3,data4)!=0) printf("Failure!!\n");
#endif

#ifdef TIMING
      start_clock(0);
      if (fast)
      {
         __Z_mul(data3, data1, data1, tweak);
      } else
      {
         mpz_mul(data4, data1, data1);
      }
      stop_clock(0);
#endif          
   }   
   
   // clean up
   gmp_randclear(state);
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
   
   for (unsigned long words = 1000UL; words < 500000000; words=floor(words*pow(2.0,1.0/64.0))) 
   {
       bits = 64*words;
#ifdef TIMING
       printf("%ld words\n",words);
       besti = 0;
       best = 1000.0;
       if (bits/64>20000)
       {
          imax = 2;
          tweak = 1;
       } else
       {
          tweak = 1;
          imax = 8;
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
#else       
       printf("%ld\n",bits);
       run_Z_Mul(TRIALS, bits, 1, 1);
#endif       
   }

   return 0;
}
