#define TRIALS 10000
#define LOGLENGTH 10 //log input poly length
#define BITS 80 //input coefficient bits

#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "radixmul.h"

/* Runs SSMul through random data.

TRIALS = number of trials to perform.
LOGLENGTH = log of degree+1.
BITS = number of bits to use in each coefficient. */

void run_RadixMul(unsigned long num_trials, unsigned long log_size,
               unsigned coeff_bits)
{
   mp_limb_t primes[3] = {9223372036854775783U,9223372036854775643U,
                          9223372036854775549U};
   
   unsigned long size = 1 << (log_size + 1);

   // alloc some space
   mpz_t* data1 = (mpz_t*) malloc(size * sizeof(mpz_t));
   mpz_t* data2 = (mpz_t*) malloc(size * sizeof(mpz_t));
   mpz_t* out = (mpz_t*) malloc((2*size-1) * sizeof(mpz_t));
   for (unsigned i = 0; i < size; i++)
   {
      mpz_init(data1[i]);
      mpz_init(data2[i]);
   }
   for (unsigned i = 0; i < 2*size-1; i++)
   {
      mpz_init(out[i]);
   }
   
   gmp_randstate_t state;
   gmp_randinit_default(state);
    
   for (unsigned i = 0; i < num_trials; i++)
   {
      if (i%20==0)
      {
         // make up random polys
         for (unsigned j = 0; j < size; j++)
         {
            mpz_urandomb(data1[j], state, coeff_bits);
            mpz_urandomb(data2[j], state, coeff_bits);
         }
      }
      
      // compute product using SSMul
      RadixMul(out, data1, data2, log_size, 3, primes);
   }
   
   // clean up
   for (unsigned i = 0; i < size; i++)
   {
      mpz_clear(data1[i]);
      mpz_clear(data2[i]);
   }
   for (unsigned i = 0; i < 2*size-1; i++)
   {
      mpz_clear(out[i]);
   }
   
   free(data1);
   free(data2);
   free(out);
}


int main (int argc, const char * argv[])
{
   run_RadixMul(TRIALS, LOGLENGTH, BITS);
   
   return 0;
}
