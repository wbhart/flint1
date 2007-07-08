/******************************************************************************

 tinyQS.h
 Header file for tinyQS.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef TINYQS_H
#define TINYQS_H

#include <gmp.h>

#include "factor_base.h"
#include "common.h"

#define MAXBITS 115
#define EXTRA_RELS 64 // number of additional relations to find above the number of primes

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;

#endif
