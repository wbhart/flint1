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

#define MINBITS 40 // Smallest bits including multiplier that can be factored

#define TEST 0 // Test that squares have come out of the linear algebra phase

#define TEST2 0 // When the large prime variant is not used we can X^2 = Y^2 mod N

#define CURVES 1

#define PRINT_FACTORS 1

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;

#endif
