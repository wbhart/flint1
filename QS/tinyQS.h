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
 Header file for tinyQS.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef TINYQS_H
#define TINYQS_H

#include <gmp.h>

#include "common.h"

#define QS_INFO 0 // Print some info about what is being factored, etc

#define MAXBITS 81 // Largest bits including multiplier that can be factored

#define TEST 0

#define PRINT_FACTORS 0

typedef struct F_mpz_fact_s
{
   mpz_t * fact;
   unsigned long num;         
} F_mpz_factor_t;

#endif
