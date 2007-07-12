/******************************************************************************

 sieve.h
 Header file for sieve.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef SIEVE_H
#define SIEVE_H

#include <gmp.h>

#include "common.h"

void do_sieving(unsigned long poly_add, unsigned long * poly_corr, 
                QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);

#endif
