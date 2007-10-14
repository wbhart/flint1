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


// set up and clean up the test module
void test_support_init();
void test_support_cleanup();


// a single GMP random state object that test modules may use
extern gmp_randstate_t randstate;


// returns random unsigned long in [0, max)
unsigned long random_ulong(unsigned long max);

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
