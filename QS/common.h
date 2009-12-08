/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/******************************************************************************

 common.h
 Common header file for tinyQS

 (C) 2006 William Hart

******************************************************************************/

#ifndef QS_COMMON_H
#define QS_COMMON_H

#include <stdint.h>

#include "../fmpz.h"
#include "../flint.h"

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
   {200, 3600, 64000, 10, 60}, // 60 digits
   {210, 6000, 64000, 12, 60},
   {220, 7500, 64000, 15, 70},
   {230, 8500, 64000, 17, 80}, // 70 digits
   {240, 18000, 64000, 19, 80}, 
   {250, 24000, 64000, 19, 80}, // 75 digits
   {260, 55000, 128000, 25, 100}, // 80 digits
   {270, 64000, 128000, 27, 100}
};

#if FLINT_BITS == 64
#define TINY_BITS 74
#else
#define TINY_BITS 54
#endif

#define PTABSIZE (sizeof(prime_tab)/(5*sizeof(unsigned long)))

#define SIEVE_BLOCK 64000

#define SECOND_PRIME 3000 // 3000 6400

#define EXTRA_RELS 64L // number of additional relations to find above the number of primes
                         
#define MAX_FACS 60 // Maximum number of different prime factors
                    // a relation can have 25, 30
                    
#define SMALL_PRIMES 8 // Used by tinyQS Todo: make this a variable

#define SIEVE_SIZE 64000 // Used by tinyQS Todo: make this a variable

#define QS_INFO 0 // Print some info about what is being factored, etc

#define TEST 0 // Test that squares have come out of the linear algebra phase

#define TEST2 0 // When the large prime variant is not used we can check X^2 = Y^2 mod N

#define CURVES 1

#define PRINT_FACTORS 1

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;
                    
typedef struct hash_s
{
	int offset;
	char size;
} hash_entry;

typedef struct prime_s
{
   uint32_t p; // prime
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
   uint32_t * sqrts;
   unsigned char * sizes;
   unsigned long * prime_count;
} QS_t;

#endif
