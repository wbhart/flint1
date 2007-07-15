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

#define QS_INFO 0 // Print some info about what is being factored, etc

#define MAXBITS 128 // Largest bits including multiplier that can be factored

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;

#endif
