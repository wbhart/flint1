/******************************************************************************

 factor_base.h
 Header file for factor_base.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef FACTORBASE_H
#define FACTORBASE_H

#include <gmp.h>

#include "tinyQS.h"
#include "common.h"

#define KSMAX 1000

static const unsigned long prime_tab[][2] =
{
   {32, 30},
   {40, 50},
   {50, 80},
   {60, 100},
   {70, 150},
   {80, 200},
   {90, 300},
   {100, 350},
   {110, 400},
   {120, 450},
   {130, 550}
};

#define PTABSIZE (sizeof(prime_tab)/(2*sizeof(unsigned long)))

unsigned long num_FB_primes(unsigned long bits);

void sqrts_init(QS_t * qs_inf);

void sqrts_clear(void);

void compute_sizes(QS_t * qs_inf);

void sizes_clear(void);

//Knuth-Schroeppel multipliers and a macro to count them

static const unsigned long multipliers[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 
                                                23, 29, 31, 37, 41, 43};

#define NUMMULTS (sizeof(multipliers)/sizeof(unsigned long))
#define max_mult_size 6 // number of bits of maximum multiplier

void primes_clear(void);

void primes_init(QS_t * qs_inf);

unsigned long knuth_schroeppel(QS_t * qs_inf);

unsigned long compute_factor_base(QS_t * qs_inf);

#endif
