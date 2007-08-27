/******************************************************************************

 common.h
 Common header file for tinyQS

 (C) 2006 William Hart

******************************************************************************/

#ifndef QS_COMMON_H
#define QS_COMMON_H

#include <sys/types.h>

#include "../fmpz.h"

// For each bitsize, this table stores, in order:
// bitsize, num_primes, sieve_size, small_primes and large_prime/factor_base[num_primes-1]

static const unsigned long prime_tab[][5] =
{
   {32, 30, 2500, 4, 1},
   {40, 50, 3000, 4, 1},
   {50, 80, 3500, 5, 1},
   {60, 100, 4000, 5, 1},
   {70, 300, 6000, 6, 1},
   {80, 400, 8000, 6, 1},
   {90, 500, 10000, 7, 1},
   {100, 650, 13000, 7, 1},
   {110, 800, 15000, 7, 1}, // 31 digits
   {120, 1000, 20000, 7, 2},
   {130, 800, 32000, 9, 1}, // 41 digits
   {140, 1200, 28000, 8, 10}, 
   {150, 1800, 32000, 8, 10},
   {160, 2000, 40000, 8, 30}, 
   {170, 2200, 64000, 9, 1}, // 50 digits 
   {180, 2400, 64000, 9, 35},
   {190, 2700, 64000, 10, 40}, 
   {200, 3600, 64000, 10, 60}, // 60 digits  5200
   {210, 6000, 64000, 12, 60},
   {220, 11000, 64000, 15, 70},
   {230, 8500, 64000, 17, 80}, // 70 digits
   {240, 24000, 64000, 19, 80},
   {250, 24000, 64000, 19, 80},
   {260, 55000, 128000, 25, 100},
   {270, 55000, 128000, 27, 100}
};

#define PTABSIZE (sizeof(prime_tab)/(5*sizeof(unsigned long)))

#define SIEVE_BLOCK 64000

#define SECOND_PRIME 3000 // 3000 6400

#define EXTRA_RELS 64L // number of additional relations to find above the number of primes
                         
#define MAX_FACS 60 // Maximum number of different prime factors
                    // a relation can have 25, 30
                    
typedef struct prime_s
{
   u_int32_t p; // prime
   double pinv; // precomputed inverse
} prime_t;

typedef struct QS_s
{
   fmpz_t n; // Number to factor = kn when multiplier is found
   mpz_t mpz_n;
   unsigned long k; // Multiplier
   unsigned long bits;
   unsigned long prec; // Number of limbs required to hold B
   unsigned long num_primes;
   unsigned long sieve_size;
   unsigned long error_bits;
   unsigned long sieve_fill;
   unsigned long small_primes;
   unsigned long large_prime;
   prime_t * factor_base; 
   u_int32_t * sqrts;
   unsigned char * sizes;
   unsigned long * prime_count;
} QS_t;

#endif
