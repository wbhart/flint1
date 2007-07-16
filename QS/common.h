/******************************************************************************

 common.h
 Common header file for tinyQS

 (C) 2006 William Hart

******************************************************************************/

#ifndef TINYQS_COMMON_H
#define TINYQS_COMMON_H

#include "../fmpz.h"

#define SIEVE_SIZE 9000

#define SMALL_PRIMES 4

#define EXTRA_RELS 64 // number of additional relations to find above the number of primes

#define MAX_FACS 25 // Maximum number of different prime factors
                    // a relation can have

typedef struct prime_s
{
   unsigned long p; // prime
   double pinv; // precomputed inverse
} prime_t;

typedef struct QS_s
{
   fmpz_t n; // Number to factor = kn when multiplier is found
   mpz_t mpz_n;
   unsigned long k; // Multiplier
   unsigned long bits;
   unsigned long num_primes;
   prime_t * factor_base; 
   unsigned long * sqrts;
   unsigned char * sizes;
   unsigned long * prime_count;
} QS_t;

#endif
