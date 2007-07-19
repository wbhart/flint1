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

#define MINBITS 128 // Smallest bits including multiplier that can be factored

#define TEST 0

#define PRINT_FACTORS 1

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;

#endif
