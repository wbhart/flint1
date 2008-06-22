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
/****************************************************************************

test-support.h: Support code for test modules

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_TEST_SUPPORT_H
#define FLINT_TEST_SUPPORT_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "profiler.h"

#define RUN_TEST(targetfunc) \
{ \
   timeit_t t0; \
   timeit_start(t0); \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();  \
   timeit_stop(t0); \
   printf("Cpu = %ld ms  Wall = %ld ms ", t0->cpu, t0->wall); \
   printf(success ? "... ok\n" : "... FAIL!\n"); \
   all_success = all_success && success;               \
} 

// set up and clean up the test module
void test_support_init();
void test_support_cleanup();


// a single GMP random state object that test modules may use
extern gmp_randstate_t randstate;


// returns random unsigned long in [0, max)
unsigned long random_ulong(unsigned long max);

unsigned long random_ulong2(unsigned long max);

// returns random limb
mp_limb_t random_limb();

// writes "limbs" random limbs to dest, using mpz_rrandomb
void random_limbs(mp_limb_t* dest, unsigned long limbs);

// writes "limbs" random limbs to dest, uniform distribution
void urandom_limbs(mp_limb_t* dest, unsigned long limbs);


// converts mpz to mpn, writes "limbs" limbs (zero-padded), ignores sign of src
void mpz_to_mpn(mp_limb_t* dest, unsigned long limbs, mpz_t src);

// convert mpn to mpz, reads exactly "limbs" limbs
void mpn_to_mpz(mpz_t dest, mp_limb_t* src, unsigned long limbs);


#define TEST_MPZ_COUNT 5

#ifdef __cplusplus
 }
#endif

#endif

// *************** end of file
