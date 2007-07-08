/******************************************************************************

 poly.h
 Header file for poly.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef POLY_H
#define POLY_H

#include <gmp.h>

#include "common.h"

typedef struct poly_s
{
    unsigned long s;
    unsigned long fact, span, min;
     
    unsigned long * A_ind;
    unsigned long * A_modp;
    unsigned long * A_inv;
    unsigned long * soln1;
    unsigned long * soln2;
    unsigned long ** A_inv2B; 
} poly_t;

void poly_init(QS_t * qs_inf, poly_t * poly_inf, mpz_t N);

void poly_clear(void);

#endif
