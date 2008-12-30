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

packed_vec-test.c: Test code for packed_vec.c and packed_vec.h

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "memory-manager.h"
#include "test-support.h"
#include "long_extras.h"
#include "packed_vec.h"

#define VARY_BITS 1
#define SPARSE 1

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

#define ITER 1 // To make the tests longer multiply this by an appropriate amount

/* Generate a random unpacked vector with the given bits per entry. */

void randarray(mp_limb_t * vec, long length, unsigned long bits)
{
   for (ulong i = 0; i < length; i++)
      vec[i] = z_randbits(bits);
} 

int test_PV_GET_SET_NEXT()
{
   int result = 1;
	pv_s pv;
   mp_limb_t * vec;
	ulong temp, i;

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
	{
      ulong bits = z_randint(FLINT_BITS) + 1;
		ulong length = z_randint(200);
		
		vec = (mp_limb_t *) flint_heap_alloc(length);
		randarray(vec, length, bits);

      pv_init(&pv, length, pv_bit_fit(bits));

		pv_iter_s iter;
		PV_ITER_INIT(iter, pv, 0);

		for (i = 0; i < length; i++)
		{
			PV_SET_NEXT(iter, vec[i]);
		}

		PV_ITER_INIT(iter, pv, 0);

		for (i = 0; (i < length) && (result == 1); i++)
		{
			PV_GET_NEXT(temp, iter);
			result = (temp == vec[i]);
		}
      
		if (!result)
		{
			i--;
			printf("bits = %ld, length = %ld, i = %ld, temp = %ld, vec[i] = %ld\n", bits, length, i, temp, vec[i]);
		}
	}

	return result;
}

int test_PV_GET_SET_PREV()
{
   int result = 1;
	pv_s pv;
   mp_limb_t * vec;
	ulong temp;
	long i;

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
	{
      ulong bits = z_randint(FLINT_BITS) + 1;
		ulong length = z_randint(200);
		
		vec = (mp_limb_t *) flint_heap_alloc(length);
		randarray(vec, length, bits);

      pv_init(&pv, length, pv_bit_fit(bits));

		pv_iter_s iter;
		PV_ITER_INIT(iter, pv, length - 1);

		for (i = length - 1; i >= 0; i--)
		{
			PV_SET_PREV(iter, vec[i]);
		}

		PV_ITER_INIT(iter, pv, length - 1);

		for (i = length - 1; (i >= 0) && (result == 1); i--)
		{
			PV_GET_PREV(temp, iter);
			result = (temp == vec[i]);
		}
      
		if (!result)
		{
			i++;
			printf("bits = %ld, length = %ld, i = %ld, temp = %ld, vec[i] = %ld\n", bits, length, i, temp, vec[i]);
		}
	}

	return result;
}

int test_PV_GET_SET_ENTRY()
{
   int result = 1;
	pv_s pv;
   mp_limb_t * vec;
	ulong temp;
	long i;

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
	{
      ulong bits = z_randint(FLINT_BITS) + 1;
		ulong length = z_randint(200)+1;
		
		vec = (mp_limb_t *) flint_heap_alloc(length);
		randarray(vec, length, bits);

      pv_init(&pv, length, pv_bit_fit(bits));

		pv_iter_s iter;
		PV_ITER_INIT(iter, pv, length - 1);

		for (i = length - 1; i >= 0; i--)
		{
			PV_SET_PREV(iter, vec[i]);
		}

		ulong entry;
		for (i = 0; (i < 100) && (result == 1); i++)
		{
			entry = z_randint(length);
			PV_GET_ENTRY(temp, pv, entry);
			result = (temp == vec[entry]);
		}
      
		if (!result)
		{
			printf("bits = %ld, length = %ld, entry = %ld, temp = %ld, vec[entry] = %ld, vec[entry-1] = %ld, vec[entry+1] = %ld\n", bits, length, entry, temp, vec[entry], vec[entry-1], vec[entry+1]);
		}
	}

	return result;
}

int test_pv_set_bits()
{
   int result = 1;
	pv_s pv;
   mp_limb_t * vec;
	ulong temp, i;

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
	{
      ulong bits = z_randint(FLINT_BITS) + 1;
		ulong length = z_randint(200) + 1;
		
		vec = (mp_limb_t *) flint_heap_alloc(length);
		randarray(vec, length, bits);

      pv_init(&pv, length, FLINT_BITS);

		pv_iter_s iter;
		PV_ITER_INIT(iter, pv, 0);

		for (i = 0; i < length; i++)
		{
			PV_SET_NEXT(iter, vec[i]);
		}

		pv.length = length;

		pv_set_bits(&pv, pv_bit_fit(bits));
		pv_set_bits(&pv, pv_bit_fit(bits + z_randint(FLINT_BITS - bits + 1)));
		pv_set_bits(&pv, pv_bit_fit(bits));
			
		PV_ITER_INIT(iter, pv, 0);

		for (i = 0; (i < length) && (result == 1); i++)
		{
			PV_GET_NEXT(temp, iter);
			result = (temp == vec[i]);
		}
      
		if (!result)
		{
			i--;
			printf("bits = %ld, length = %ld, i = %ld, temp = %ld, vec[i] = %ld\n", bits, length, i, temp, vec[i]);
		}
	}

	return result;
}

void zmod_poly_test_all()
{
   int success, all_success = 1;

	RUN_TEST(PV_GET_SET_NEXT); 
   RUN_TEST(PV_GET_SET_PREV); 
   RUN_TEST(PV_GET_SET_ENTRY); 
   RUN_TEST(pv_set_bits); 

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   zmod_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}
