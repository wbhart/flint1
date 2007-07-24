/******************************************************************************

 mp_poly.h
 Header file for poly.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MPPOLY_H
#define MPPOLY_H

#include <gmp.h>

#include "common.h"

#define POLY_PARAMS 1 // Print the parameters being used to choose polynomials

#define POLY_A 0 // Print target A and actual A

#define TEST_C 0 // Test polynomial coefficients 

#define B_TERMS 0 // Print out the B_terms

typedef struct poly_s
{
    unsigned long s;
    unsigned long fact, span, min;
    unsigned long * target_A;
    unsigned long * A;
    unsigned long * B;
    mpz_t A_mpz;
    mpz_t B_mpz;
    mpz_t C;
     
    unsigned long * A_ind;
    unsigned long * A_modp;
    unsigned long * A_inv;
    unsigned long * soln1;
    unsigned long * soln2;
    char * * posn1;
    char * * posn2;
    unsigned long ** A_inv2B; 
    double * inv_p2;
    
    unsigned long * B_terms;
} poly_t;

void poly_init(QS_t * qs_inf, poly_t * poly_inf, mpz_t N);

void poly_clear(poly_t * poly_inf);

void compute_A(QS_t * qs_inf, poly_t * poly_inf);

void compute_B_terms(QS_t * qs_inf, poly_t * poly_inf);

void compute_off_adj(QS_t * qs_inf, poly_t * poly_inf);

void compute_A_factor_offsets(QS_t * qs_inf, poly_t * poly_inf);

void compute_B_C(QS_t * qs_inf, poly_t * poly_inf);

#endif
