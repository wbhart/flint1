/******************************************************************************

 tinyQS.h
 Header file for mpQS.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MPQS_H
#define MPQS_H

#include <gmp.h>

#include "mp_factor_base.h"
#include "common.h"

#define QS_INFO 1 // Print some info about what is being factored, etc

#define MINBITS 64 // Smallest bits including multiplier that can be factored

#define TEST 0 // Test that squares have come out of the linear algebra phase

#define TEST2 0 // When the large prime variant is not used we can X^2 = Y^2 mod N

#define CURVES 1

#define PRINT_FACTORS 1

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;

#endif
