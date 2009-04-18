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

test-support.c: Support code for test modules

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "test-support.h"


THREAD gmp_randstate_t randstate;
THREAD int rand_initialised = 0;


// a bunch of global mpz's, guaranteed to be init'd
THREAD mpz_t test_mpz[TEST_MPZ_COUNT];


void test_support_init()
{
   gmp_randinit_default(randstate);
   
   for (unsigned long i = 0; i < TEST_MPZ_COUNT; i++)
      mpz_init(test_mpz[i]);
}


void test_support_cleanup()
{
   gmp_randclear(randstate);

   for (unsigned long i = 0; i < TEST_MPZ_COUNT; i++)
      mpz_clear(test_mpz[i]);
}


unsigned long random_ulong(unsigned long max)
{
   if (!rand_initialised) 
	{
		test_support_init();
		rand_initialised = 1;
	}
	return gmp_urandomm_ui(randstate, max);
}

unsigned long random_ulong2(unsigned long max)
{
   if (!rand_initialised) 
	{
		test_support_init();
		rand_initialised = 1;
	}
	unsigned long bits = FLINT_BIT_COUNT(max);
   mpz_t rand;
   mpz_init(rand);
   mpz_rrandomb(rand, randstate, bits);
   unsigned long randi = mpz_get_ui(rand) % max; 
   mpz_clear(rand);
   return randi;
}

mp_limb_t random_limb()
{
   if (!rand_initialised) 
	{
		test_support_init();
		rand_initialised = 1;
	}
	return gmp_urandomb_ui(randstate, FLINT_BITS);
}


void urandom_limbs(mp_limb_t* dest, unsigned long limbs)
{
   if (!rand_initialised) 
	{
		test_support_init();
		rand_initialised = 1;
	}
	for (unsigned long i = 0; i < limbs; i++)
      dest[i] = gmp_urandomb_ui(randstate, FLINT_BITS);
}


void random_limbs(mp_limb_t* dest, unsigned long limbs)
{
   if (!rand_initialised) 
	{
		test_support_init();
		rand_initialised = 1;
	}
	mpz_rrandomb(test_mpz[0], randstate, limbs*FLINT_BITS);

   if (random_ulong(2))
   {
      // GMP always sets the high bit equal to 1,
      // so with probability 1/2 we flip all the bits
      mpz_set_ui(test_mpz[1], 1);
      mpz_mul_2exp(test_mpz[1], test_mpz[1], limbs*FLINT_BITS);
      mpz_sub_ui(test_mpz[1], test_mpz[1], 1);
      mpz_sub(test_mpz[0], test_mpz[1], test_mpz[0]);
   }

   memset(dest, 0, limbs * sizeof(mp_limb_t));
   mpz_export(dest, NULL, -1, sizeof(mp_limb_t), 0, 0, test_mpz[0]);
}


void mpz_to_mpn(mp_limb_t* dest, unsigned long limbs, mpz_t src)
{
   memset(dest, 0, limbs * sizeof(mp_limb_t));
   mpz_export(dest, NULL, -1, sizeof(mp_limb_t), 0, 0, src);
}


void mpn_to_mpz(mpz_t dest, mp_limb_t* src, unsigned long limbs)
{
   mpz_import(dest, limbs, -1, sizeof(mp_limb_t), 0, 0, src);
}


// end of file ****************************************************************
