/******************************************************************************

 common.h
 Common header file for tinyQS

 (C) 2006 William Hart

******************************************************************************/

#ifndef QS_COMMON_H
#define QS_COMMON_H

#include "../fmpz.h"

// For each bitsize, this table stores, in order:
// bitsize, num_primes, sieve_size, small_primes and large_prime/factor_base[num_primes-1]

static const unsigned long prime_tab[][5] =
{
   {32, 30, 2500, 4, 10},
   {40, 50, 3000, 4, 10},
   {50, 80, 3500, 5, 10},
   {60, 100, 4000, 5, 10},
   {70, 150, 6000, 6, 15},
   {80, 200, 8000, 6, 15},
   {90, 200, 10000, 7, 15},
   {100, 250, 13000, 7, 20},
   {110, 300, 17000, 7, 20},
   {120, 500, 20000, 7, 20},
   {130, 550, 24000, 8, 25}, // 37 digits
   {140, 800, 28000, 8, 25}, // 40 digits
   {150, 1000, 32000, 8, 30},
   {160, 2200, 40000, 8, 30}, 
   {170, 2200, 64000, 9, 35}, // 50 digits 
   {180, 2400, 64000, 9, 35},
   {190, 2700, 64000, 10, 40}, 
   {200, 3600, 64000, 10, 60}, // 60 digits  5200
   {210, 6000, 64000, 12, 60},
   {220, 11000, 64000, 15, 70},
   {230, 8500, 64000, 17, 80}, // 70 digits
   {240, 24000, 3*64000, 19, 80},
   {250, 24000, 3*64000, 19, 80},
   {260, 55000, 128000, 25, 100},
   {270, 55000, 128000, 27, 100}
};

#define PTABSIZE (sizeof(prime_tab)/(5*sizeof(unsigned long)))

#define SIEVE_BLOCK 64000

#define SECOND_PRIME 3500

#define EXTRA_RELS 64L // number of additional relations to find above the number of primes
                         // -150L, -400L
                         
#define MAX_FACS 60 // Maximum number of different prime factors
                    // a relation can have 25, 30
                    
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
   unsigned long prec; // Number of limbs required to hold B
   unsigned long num_primes;
   unsigned long sieve_size;
   unsigned long error_bits;
   unsigned long sieve_fill;
   unsigned long small_primes;
   unsigned long large_prime;
   prime_t * factor_base; 
   unsigned long * sqrts;
   unsigned char * sizes;
   unsigned long * prime_count;
} QS_t;

#endif
