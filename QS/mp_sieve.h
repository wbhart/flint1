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

 mp_sieve.h
 Header file for sieve.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MPSIEVE_H
#define MPSIEVE_H

#include <gmp.h>

#include "mp_linear_algebra.h"
#include "common.h"

#define RELATIONS 0 // Print out relations as they are generated

#define POLYS 0 // Print out polynomials and offsets in candidate evaluation

void get_sieve_params(QS_t * qs_inf);

void do_sieving(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve, 
                  unsigned long first_prime, unsigned long second_prime, 
                  unsigned long M, int first, int last);
                                                  
void do_sieving2(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);

void do_sieving3(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve, 
                  unsigned long first_prime, unsigned long second_prime, 
                  unsigned long M);

void do_sieving4(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve, 
                  unsigned long first_prime, unsigned long second_prime, 
                  unsigned long M);
                
void update_offsets(unsigned long poly_add, uint32_t * poly_corr, 
                                        QS_t * qs_inf, poly_t * poly_inf);

unsigned long evaluate_sieve(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);

unsigned long evaluate_candidate(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf, 
                                     unsigned long i, unsigned char * sieve);
                                     
#endif
