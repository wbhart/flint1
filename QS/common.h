/******************************************************************************

 common.h
 Common header file for tinyQS

 (C) 2006 William Hart

******************************************************************************/

#ifndef TINYQS_COMMON_H
#define TINYQS_COMMON_H

#include "../fmpz.h"

#define SIEVE_SIZE 16384

typedef struct prime_s
{
   unsigned long p; // prime
   double pinv; // precomputed inverse
} prime_t;

typedef struct QS_s
{
   fmpz_t n; // Number to factor = kn when multiplier is found
   unsigned long k; // Multiplier
   unsigned long bits;
   unsigned long num_primes;
   prime_t * factor_base; 
   unsigned long * sqrts;
   unsigned char * sizes;
} QS_t;

#endif
