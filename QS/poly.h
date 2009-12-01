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

 poly.h
 Header file for poly.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef POLY_H
#define POLY_H

#include <gmp.h>

#include "common.h"

#define POLY_PARAMS 0 // Print the parameters being used to choose polynomials
#define POLY_A 0 // Print target A and actual A

typedef struct poly_s
{
    unsigned long s;
    unsigned long fact, span, min;
    unsigned long target_A;
    unsigned long A;
    unsigned long B;
    mpz_t C;
     
    unsigned long * A_ind;
    unsigned long * A_modp;
    unsigned long * A_inv;
    unsigned long * soln1;
    unsigned long * soln2;
    unsigned long ** A_inv2B; 
    double * inv_p2;
    
    unsigned long * B_terms;
} poly_t;

void tiny_poly_init(QS_t * qs_inf, poly_t * poly_inf, mpz_t N);

void tiny_poly_clear(poly_t * poly_inf);

void tiny_compute_A(QS_t * qs_inf, poly_t * poly_inf);

void tiny_compute_B_terms(QS_t * qs_inf, poly_t * poly_inf);

void tiny_compute_C(QS_t * qs_inf, poly_t * poly_inf);

void tiny_compute_off_adj(QS_t * qs_inf, poly_t * poly_inf);

void tiny_compute_A_factor_offsets(QS_t * qs_inf, poly_t * poly_inf);

#endif
