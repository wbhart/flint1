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

 factor_base.h
 Header file for factor_base.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef FACTORBASE_H
#define FACTORBASE_H

#include <gmp.h>

#include "tinyQS.h"
#include "common.h"

#define KSMAX 100

static const unsigned long prime_tab_small[][2] =
{
   {0, 50},
	{30, 60},
   {40, 70},
   {50, 80},
   {60, 100},
   {70, 140},
   {80, 180},
   {90, 200},
   {100, 250},
   {110, 300},
   {120, 500},
   {130, 550}
};

#define PTABSIZE_SMALL (sizeof(prime_tab_small)/(2*sizeof(unsigned long)))

unsigned long tiny_num_FB_primes(unsigned long bits);

void tiny_sqrts_init(QS_t * qs_inf);

void tiny_sqrts_clear(void);

void tiny_compute_sizes(QS_t * qs_inf);

void tiny_sizes_clear(void);

//Knuth-Schroeppel multipliers and a macro to count them

static const unsigned long multipliers[] = {1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 
                                            15, 17, 19, 21, 22, 23, 26, 29, 
                                            30, 31, 33, 34, 35, 37, 38, 41, 
                                            42, 43, 47};

#define NUMMULTS (sizeof(multipliers)/sizeof(unsigned long))
#define max_mult_size 6 // number of bits of maximum multiplier

void tiny_primes_clear(void);

void tiny_primes_init(QS_t * qs_inf);

unsigned long tiny_knuth_schroeppel(QS_t * qs_inf);

unsigned long tiny_compute_factor_base(QS_t * qs_inf);

#endif
