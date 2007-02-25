#define TRIALS 10
#define LENGTH 2560
#define BITS 6000

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "ssmul.h"
#include "Zvec.h"

/* Runs SSMul through random data.

TRIALS = number of trials to perform.
LOGLENGTH = log of degree+1.
BITS = number of bits to use in each coefficient. */

void run_SSMul(unsigned long num_trials, unsigned long length,
               unsigned coeff_bits)
{
   // alloc some space
   mpz_t* data1 = (mpz_t*) malloc(length * sizeof(mpz_t));
   mpz_t* data2 = (mpz_t*) malloc(length * sizeof(mpz_t));
   mpz_t* data3 = (mpz_t*) malloc((2*length-1) * sizeof(mpz_t));
   for (unsigned i = 0; i < length; i++)
   {
      mpz_init(data1[i]);
      mpz_init(data2[i]);
   }
   for (unsigned i = 0; i < 2*length-1; i++)
   {
      mpz_init(data3[i]);
   }
   
   gmp_randstate_t state;
   gmp_randinit_default(state);
   
   Zvec a,b,c;
   a->coords = data1;
   b->coords = data2;
   c->coords = data3;
   a->length = length;
   b->length = length;
   c->length = 2*length-1;
   
   for (unsigned i = 0; i < num_trials; i++)
   {
      // make up random polys
      if (i%2==0)
      {
         for (unsigned j = 0; j < length; j++)
         {
            mpz_urandomb(data1[j], state, coeff_bits);
            mpz_urandomb(data2[j], state, coeff_bits);
            /*if (gmp_urandomb_ui(state, 1))
               mpz_neg(data1[j], data1[j]);
            if (gmp_urandomb_ui(state, 1))
               mpz_neg(data2[j], data2[j]);*/
         }
      }
      
      // compute product using SSMul
      Zvec_mul(c,a,b);
      printf("Polynomial multiplication %ld completed.\n",i+1);
   }
   
   // clean up
   for (unsigned i = 0; i < length; i++)
   {
      mpz_clear(data1[i]);
      mpz_clear(data2[i]);
   }
   for (unsigned i = 0; i < 2*length-1; i++)
   {
      mpz_clear(data3[i]);
   }
   
   free(data1);
   free(data2);
   free(data3);
}


int main (int argc, const char * argv[])
{
   run_SSMul(TRIALS, LENGTH, BITS);
   
   return 0;
}
