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

 sieve.h
 Header file for sieve.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef SIEVE_H
#define SIEVE_H

#include <gmp.h>

#include "linear_algebra.h"
#include "common.h"

#define RELATIONS 1 // Print out relations as they are generated
#define POLYS 0 // Print out polynomials and offsets in candidate evaluation

void do_sieving(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);
                
void update_offsets(unsigned long poly_add, unsigned long * poly_corr, 
                                        QS_t * qs_inf, poly_t * poly_inf);

unsigned long evaluate_sieve(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve);

unsigned long evaluate_candidate(QS_t * qs_inf, poly_t * poly_inf, 
                                     unsigned long i, unsigned char * sieve);
                                     
#endif
