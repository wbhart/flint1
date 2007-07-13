/******************************************************************************

 sieve.h
 Header file for sieve.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef SIEVE_H
#define SIEVE_H

#include <gmp.h>

#include "common.h"

#define RELATIONS 1 // Print out relations as they are generated
#define POLYS 0 // Print out polynomials and offsets in candidate evaluation

void do_sieving(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);
                
void update_offsets(unsigned long poly_add, unsigned long * poly_corr, 
                                        QS_t * qs_inf, poly_t * poly_inf);

unsigned long evaluate_sieve(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);

int evaluate_candidate(QS_t * qs_inf, poly_t * poly_inf, 
                                     unsigned long i, unsigned char * sieve);

#endif
