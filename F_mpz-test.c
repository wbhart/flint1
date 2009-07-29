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

F_mpz-test.c: Test code for F_mpz.c and F_mpz.h

Copyright (C) 2008, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include <pthread.h>
#include "flint.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "test-support.h"
#include "F_mpz.h"

#define VARY_BITS 1 // random coefficients have random number of bits up to the limit given
#define SIGNS 1 // random coefficients will be randomly signed
#define ITER 1 // if you want all tests to run longer, increase this

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // allows easy switching of debugging code on and off when debugging (if inserted)
#define DEBUG2 1 

/*
   Generate a random F_mpz_t with the given number of bits.
	If SIGNS is 1 then the value will be randomly signed, else it will be positive.
	Warning : do not use this function to test F_mpz_set_mpz!!
*/

void F_mpz_test_random(F_mpz_t f, ulong bits)
{
	if (bits == 0)
	{
		F_mpz_zero(f);
      return;
	}
	
	mpz_t temp;
	mpz_init(temp);
	
	pthread_mutex_lock(&F_mpz_random_mutex);
	mpz_rrandomb(temp, randstate, bits);
	pthread_mutex_unlock(&F_mpz_random_mutex);
#if SIGNS
	if (z_randint(2)) mpz_neg(temp, temp);
#endif
   
	F_mpz_set_mpz(f, temp);

   mpz_clear(temp);
}

int test_F_mpz_getset_si()
{
   F_mpz_t f;
   int result = 1;
   ulong bits, val_bits;
	long val, val2;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      
		for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random long
         val_bits = z_randint(FLINT_BITS - 1);
         val = z_randbits(val_bits);
         if (z_randint(2)) val = -val;
              
	      F_mpz_set_si(f, val);
         val2 = F_mpz_get_si(f);

		   result = (val2 == val);
		   if (!result)
	      {
			   printf("Error: val = %ld, val2 = %ld\n", val, val2);
		   }
		}

      F_mpz_clear(f);
   }
   
   return result; 
}

int test_F_mpz_getset_ui()
{
   F_mpz_t f;
   int result = 1;
   ulong bits, val_bits, val, val2;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random long
         val_bits = z_randint(FLINT_BITS);
         val = z_randbits(val_bits);
              
	      F_mpz_set_ui(f, val);
         val2 = F_mpz_get_ui(f);

		   result = (val2 == val);
		   if (!result)
	      {
			   printf("Error: val = %u, val2 = %u\n", val, val2);
		   }
		}

      F_mpz_clear(f);
   }
   
   return result; 
}

int test_F_mpz_getset_mpz()
{
   F_mpz_t f;
   int result = 1;
   ulong bits;
   mpz_t val, val2;
   ulong val_bits;
	mpz_init(val);
	mpz_init(val2);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;      
      
		F_mpz_init2(f, z_randint(10));

      F_mpz_test_random(f, bits); 
		    
      for (ulong count2 = 0; (count2 < 100) && result == 1; count2++)
      {
         val_bits = z_randint(200);
         pthread_mutex_lock(&F_mpz_random_mutex);
	      mpz_rrandomb(val, randstate, val_bits);
         pthread_mutex_unlock(&F_mpz_random_mutex);
	     
			if (z_randint(2)) mpz_neg(val, val);
              
			F_mpz_set_mpz(f, val);
         F_mpz_get_mpz(val2, f);
				  
			result = (mpz_cmp(val2, val) == 0);
			if (!result)
			{
			   gmp_printf("Error: val = %Zd, val2 = %Zd\n", val, val2);
			}
      }

      F_mpz_clear(f);
   }
   
   mpz_clear(val);
	mpz_clear(val2);
	return result; 
}

int test_F_mpz_get_d_2exp()
{
   mpz_t m1;
	double fd, md;
   F_mpz_t f1;
   int result = 1;
   ulong bits, expf, expm;
   
   mpz_init(m1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      
		fd = F_mpz_get_d_2exp(&expf, f1);
		md = mpz_get_d_2exp(&expm, m1);

      result = ((fd == md) && (expf == expm)); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, md = %f, fd = %f, expf = %ld, expm = %ld\n", bits, m1, md, fd, expf, expm);
		}
          
      F_mpz_clear(f1);
   }
   
   mpz_clear(m1);
   
   return result;
}

int test_F_mpz_getset_limbs()
{
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits, limbs, limbs2;
	mp_limb_t * arr;
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(1000)+ 1; 
      limbs = (bits - 1)/FLINT_BITS + 1;
		arr = flint_heap_alloc(limbs);

		F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f1, bits); 
		    
         limbs2 = F_mpz_get_limbs(arr, f1);
         F_mpz_set_limbs(f2, arr, limbs2);
			if (F_mpz_sgn(f1) < 0) F_mpz_neg(f2, f2);
				  
		   result = (F_mpz_equal(f2, f1));
		   if (!result)
		   {
		  	   printf("Error: f1 = "); F_mpz_print(f1); 
			   printf(" f2 = "); F_mpz_print(f2); printf("\n");
		   }
		}

      flint_heap_free(arr);
		F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
   return result; 
}

int test_F_mpz_set()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits;
   
   mpz_init(m1); 
   mpz_init(m2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_set(f2, f1);
		F_mpz_get_mpz(m2, f2);
          
      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// check aliasing
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_set(f1, f1);
		F_mpz_get_mpz(m2, f1);
          
      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

int test_F_mpz_equal()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits, bits2;
   
   mpz_init(m1); 
   mpz_init(m2); 

   // Check in case when operands are equal
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      F_mpz_get_mpz(m1, f1);
           
      F_mpz_set(f2, f1);
          
      result = (F_mpz_equal(f1, f2)); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd\n", bits, m1);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// Check in case when operands are not likely equal
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
          
      result = ((F_mpz_equal(f1, f2) && (mpz_cmp(m1, m2) == 0))
			    || (!F_mpz_equal(f1, f2) && (mpz_cmp(m1, m2) != 0))); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd\n", bits, bits2, m1, m2);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// check aliasing
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
          
      result = (F_mpz_equal(f1, f1)); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd\n", bits, m1);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

int test_F_mpz_swap()
{
   mpz_t m1, m2, m3, m4;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits, bits2;
   
   mpz_init(m1); 
   mpz_init(m2); 
   mpz_init(m3); 
   mpz_init(m4); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      F_mpz_swap(f2, f1);
		F_mpz_get_mpz(m3, f1);
      F_mpz_get_mpz(m4, f2);
          
      result = ((mpz_cmp(m1, m4) == 0) && (mpz_cmp(m2, m3) == 0)); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// check aliasing
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_swap(f1, f1);
		F_mpz_get_mpz(m2, f1);
          
      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   mpz_clear(m3);
   mpz_clear(m4);
   
   return result;
}

int test_F_mpz_neg()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits;
   
   mpz_init(m1); 
   mpz_init(m2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_neg(f2, f1);
		F_mpz_get_mpz(m2, f2);
      mpz_neg(m2, m2);

      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// check aliasing
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_neg(f1, f1);
		F_mpz_get_mpz(m2, f1);
      mpz_neg(m2, m2);

      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

int test_F_mpz_abs()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits;
   
   mpz_init(m1); 
   mpz_init(m2); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_abs(f2, f1);
		F_mpz_get_mpz(m2, f2);
      mpz_abs(m1, m1);

      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
	// check aliasing
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      
      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_abs(f1, f1);
		F_mpz_get_mpz(m2, f1);
      mpz_abs(m1, m1);

      result = (mpz_cmp(m1, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd\n", bits, m1, m2);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

int test_F_mpz_add()
{
   mpz_t m1, m2, m3, m4;
   F_mpz_t f1, f2, f3;
   int result = 1;
   ulong bits, bits2;
   
   mpz_init(m1); 
   mpz_init(m2); 
   mpz_init(m3); 
   mpz_init(m4); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      F_mpz_init(f3);

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      F_mpz_add(f3, f2, f1);
		F_mpz_get_mpz(m3, f3);
      mpz_add(m4, m1, m2);
          
      result = (mpz_cmp(m4, m3) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
      F_mpz_clear(f3);
   }
   
	// check aliasing of operands 1 and 2
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_add(f1, f1, f2);
		F_mpz_get_mpz(m3, f1);
      mpz_add(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of operands 1 and 3
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_add(f2, f1, f2);
		F_mpz_get_mpz(m3, f2);
      mpz_add(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of all operands
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      
		F_mpz_add(f1, f1, f1);
		F_mpz_get_mpz(m2, f1);
      mpz_add(m3, m1, m1);

      result = (mpz_cmp(m3, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd\n", bits, m1, m2, m3);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   mpz_clear(m3);
   mpz_clear(m4);
   
   return result;
}

int test_F_mpz_sub()
{
   mpz_t m1, m2, m3, m4;
   F_mpz_t f1, f2, f3;
   int result = 1;
   ulong bits, bits2;
   
   mpz_init(m1); 
   mpz_init(m2); 
   mpz_init(m3); 
   mpz_init(m4); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      F_mpz_init(f3);

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      F_mpz_sub(f3, f1, f2);
		F_mpz_get_mpz(m3, f3);
      mpz_sub(m4, m1, m2);
          
      result = (mpz_cmp(m4, m3) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
      F_mpz_clear(f3);
   }
   
	// check aliasing of operands 1 and 2
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_sub(f1, f1, f2);
		F_mpz_get_mpz(m3, f1);
      mpz_sub(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of operands 1 and 3
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_sub(f2, f1, f2);
		F_mpz_get_mpz(m3, f2);
      mpz_sub(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of all operands
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      
		F_mpz_sub(f1, f1, f1);
		F_mpz_get_mpz(m2, f1);
      mpz_sub(m3, m1, m1);

      result = (mpz_cmp(m3, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd\n", bits, m1, m2, m3);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   mpz_clear(m3);
   mpz_clear(m4);
   
   return result;
}

int test_F_mpz_mul_ui()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, val_bits, val;
	mpz_t m1, m2;

	mpz_init(m1);
   mpz_init(m2);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(g, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, g);
			
			F_mpz_mul_ui(f, g, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_mul_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_mul_ui(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_mul_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	
	return result; 
}

int test_F_mpz_mul_si()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, val_bits;
	long val;
	mpz_t m1, m2;

	mpz_init(m1);
   mpz_init(m2);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(g, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS);
         val = z_randbits(val_bits);
			if (z_randint(2)) val = -val;
              
	      F_mpz_get_mpz(m1, g);
			
			F_mpz_mul_si(f, g, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_mul_si(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS);
         val = z_randbits(val_bits);
         if (z_randint(2)) val = -val;

	      F_mpz_get_mpz(m1, f);
			
			F_mpz_mul_si(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_mul_si(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	
	return result; 
}

int test_F_mpz_mul()
{
   mpz_t m1, m2, m3, m4;
   F_mpz_t f1, f2, f3;
   int result = 1;
   ulong bits, bits2;
   
   mpz_init(m1); 
   mpz_init(m2); 
   mpz_init(m3); 
   mpz_init(m4); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      F_mpz_init(f3);

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      F_mpz_mul2(f3, f2, f1);
		F_mpz_get_mpz(m3, f3);
      mpz_mul(m4, m1, m2);
          
      result = (mpz_cmp(m4, m3) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
      F_mpz_clear(f3);
   }
   
	// check aliasing of operands 1 and 2
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_mul2(f1, f1, f2);
		F_mpz_get_mpz(m3, f1);
      mpz_mul(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of operands 1 and 3
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_mul2(f2, f1, f2);
		F_mpz_get_mpz(m3, f2);
      mpz_mul(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of all operands
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
           
      F_mpz_get_mpz(m1, f1);
      
		F_mpz_mul2(f1, f1, f1);
		F_mpz_get_mpz(m2, f1);
      mpz_mul(m3, m1, m1);

      result = (mpz_cmp(m3, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd\n", bits, m1, m2, m3);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   mpz_clear(m3);
   mpz_clear(m4);
   
   return result;
}

int test_F_mpz_mul_2exp()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, exp;
	mpz_t m1, m2;

	mpz_init(m1);
   mpz_init(m2);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(g, bits); 

		   // Generate a random unsigned long
         exp = z_randint(200);
              
	      F_mpz_get_mpz(m1, g);
			
			F_mpz_mul_2exp(f, g, exp);
         
			F_mpz_get_mpz(m2, f);

			mpz_mul_2exp(m1, m1, exp);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd, exp = %ld\n", m1, m2, exp);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         exp = z_randint(200);
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_mul_2exp(f, f, exp);
         
			F_mpz_get_mpz(m2, f);

			mpz_mul_2exp(m1, m1, exp);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	
	return result; 
}

int test_F_mpz_add_ui()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, val_bits, val;
	mpz_t m1, m2, m3;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(g, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_add_ui(f, g, val);
         
			F_mpz_get_mpz(m3, f);

			mpz_add_ui(m1, m2, val);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_add_ui(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_add_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	
	return result; 
}

int test_F_mpz_sub_ui()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, val_bits, val;
	mpz_t m1, m2, m3;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(g, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_sub_ui(f, g, val);
         
			F_mpz_get_mpz(m3, f);

			mpz_sub_ui(m1, m2, val);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_sub_ui(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_sub_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	
	return result; 
}

int test_F_mpz_addmul_ui()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, bits2, val_bits, val;
	mpz_t m1, m2, m3;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_addmul_ui(f, g, val);
         
			F_mpz_get_mpz(m3, f);

			mpz_addmul_ui(m1, m2, val);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m3 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_addmul_ui(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_addmul_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	
	return result; 
}

int test_F_mpz_submul_ui()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, bits2, val_bits, val;
	mpz_t m1, m2, m3;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_submul_ui(f, g, val);
         
			F_mpz_get_mpz(m3, f);

			mpz_submul_ui(m1, m2, val);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m3 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_submul_ui(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_submul_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	
	return result; 
}

int test_F_mpz_addmul()
{
   F_mpz_t f, g, h;
   int result = 1;
   ulong bits, bits2, bits3;
	mpz_t m1, m2, m3, m4;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   mpz_init(m4);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
      F_mpz_init2(h, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 
         bits3 = z_randint(200)+ 1;
         F_mpz_test_random(h, bits3); 
     
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			F_mpz_get_mpz(m3, h);
			
			F_mpz_addmul(f, g, h);
         
			F_mpz_get_mpz(m4, f);

			mpz_addmul(m1, m2, m3);

		   result = (mpz_cmp(m1, m4) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m4 = %Zd\n", m1, m4);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
      F_mpz_clear(h);
   }
   
	// Check aliasing of arguments 1 and 2
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 
         
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_addmul(f, f, g);
         
			F_mpz_get_mpz(m3, f);

			mpz_addmul(m1, m1, m2);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m3 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing of arguments 1 and 3
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 
         
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_addmul(g, f, g);
         
			F_mpz_get_mpz(m3, g);

			mpz_addmul(m2, m1, m2);

		   result = (mpz_cmp(m2, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m2 = %Zd, m3 = %Zd\n", m2, m3);
		   }
		}

      F_mpz_clear(f);
   }
   
	// Check aliasing of all arguments
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         
			F_mpz_get_mpz(m1, f);
			
			F_mpz_addmul(f, f, f);
         
			F_mpz_get_mpz(m2, f);

			mpz_addmul(m1, m1, m1);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	mpz_clear(m4);
	
	return result; 
}

int test_F_mpz_submul()
{
   F_mpz_t f, g, h;
   int result = 1;
   ulong bits, bits2, bits3;
	mpz_t m1, m2, m3, m4;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   mpz_init(m4);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
      F_mpz_init2(h, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 
         bits3 = z_randint(200)+ 1;
         F_mpz_test_random(h, bits3); 
     
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			F_mpz_get_mpz(m3, h);
			
			F_mpz_submul(f, g, h);
         
			F_mpz_get_mpz(m4, f);

			mpz_submul(m1, m2, m3);

		   result = (mpz_cmp(m1, m4) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m4 = %Zd\n", m1, m4);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
      F_mpz_clear(h);
   }
   
   // Check aliasing of arguments 1 and 2
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 
         
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_submul(f, f, g);
         
			F_mpz_get_mpz(m3, f);

			mpz_submul(m1, m1, m2);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m3 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }

   // Check aliasing of arguments 1 and 3
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         bits2 = z_randint(200)+ 1;
         F_mpz_test_random(g, bits2); 
         
	      F_mpz_get_mpz(m1, f);
			F_mpz_get_mpz(m2, g);
			
			F_mpz_submul(g, f, g);
         
			F_mpz_get_mpz(m3, g);

			mpz_submul(m2, m1, m2);

		   result = (mpz_cmp(m2, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m2 = %Zd, m3 = %Zd\n", m2, m3);
		   }
		}

      F_mpz_clear(f);
   }

   // Check aliasing of all arguments
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {

		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			bits = z_randint(200)+ 1;
         F_mpz_test_random(f, bits); 
         
			F_mpz_get_mpz(m1, f);
			
			F_mpz_submul(f, f, f);
         
			F_mpz_get_mpz(m2, f);

			mpz_submul(m1, m1, m1);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	mpz_clear(m4);
	
	return result; 
}

int test_F_mpz_mod_ui()
{
   F_mpz_t f, g;
   int result = 1;
   ulong bits, val_bits, val;
	mpz_t m1, m2, m3;

	mpz_init(m1);
   mpz_init(m2);
   mpz_init(m3);
   
   for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
      F_mpz_init2(g, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(g, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
			if (val == 0L) val = 1L;
              
	      F_mpz_get_mpz(m2, g);
			
			F_mpz_mod_ui(f, g, val);
         
			F_mpz_get_mpz(m3, f);

			mpz_mod_ui(m1, m2, val);

		   result = (mpz_cmp(m1, m3) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m3);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   // Check aliasing
	for (ulong count1 = 0; (count1 < 10000*ITER) && (result == 1); count1++)
   {
      bits = z_randint(200)+ 1;
      
		// Make f a random number of limbs in size to start with
		F_mpz_init2(f, z_randint(10));
          
      for (ulong count2 = 0; (count2 < 100) && (result == 1); count2++)
		{
			F_mpz_test_random(f, bits); 

		   // Generate a random unsigned long
         val_bits = z_randint(FLINT_BITS + 1);
         val = z_randbits(val_bits);
			if (val == 0L) val = 1L;
              
	      F_mpz_get_mpz(m1, f);
			
			F_mpz_mod_ui(f, f, val);
         
			F_mpz_get_mpz(m2, f);

			mpz_mod_ui(m1, m1, val);

		   result = (mpz_cmp(m1, m2) == 0);
		   if (!result)
	      {
			   gmp_printf("Error: m1 = %Zd, m2 = %Zd\n", m1, m2);
		   }
		}

      F_mpz_clear(f);
      F_mpz_clear(g);
   }
   
   mpz_clear(m1);
	mpz_clear(m2);
	mpz_clear(m3);
	
	return result; 
}

int test_F_mpz_mod()
{
   mpz_t m1, m2, m3, m4;
   F_mpz_t f1, f2, f3;
   int result = 1;
   ulong bits, bits2;
   
   mpz_init(m1); 
   mpz_init(m2); 
   mpz_init(m3); 
   mpz_init(m4); 

   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      F_mpz_init(f3);

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
      if (F_mpz_is_zero(f2)) 
			F_mpz_set_ui(f2, 1L);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
		
		F_mpz_mod(f3, f1, f2);
		F_mpz_get_mpz(m3, f3);
      mpz_mod(m4, m1, m2);
          
      result = (mpz_cmp(m4, m3) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
      F_mpz_clear(f3);
   }
   
	// check aliasing of operands 1 and 2
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
      if (F_mpz_is_zero(f2)) 
			F_mpz_set_ui(f2, 1L);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_mod(f1, f1, f2);
		F_mpz_get_mpz(m3, f1);
      mpz_mod(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of operands 1 and 3
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      bits2 = z_randint(200) + 1;
      F_mpz_test_random(f2, bits2);
      if (F_mpz_is_zero(f2)) 
			F_mpz_set_ui(f2, 1L);
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		F_mpz_mod(f2, f1, f2);
		F_mpz_get_mpz(m3, f2);
      mpz_mod(m4, m1, m2);

      result = (mpz_cmp(m3, m4) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }

   // check aliasing of all operands
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      if (F_mpz_is_zero(f1)) 
			F_mpz_set_ui(f1, 1L);
          
      F_mpz_get_mpz(m1, f1);
      
		F_mpz_mod(f1, f1, f1);
		F_mpz_get_mpz(m2, f1);
      mpz_mod(m3, m1, m1);

      result = (mpz_cmp(m3, m2) == 0); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd\n", bits, m1, m2, m3);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   mpz_clear(m3);
   mpz_clear(m4);
   
   return result;
}

int test_F_mpz_invert()
{
   mpz_t m1, m2, m3, m4;
   F_mpz_t f1, f2, f3;
   int result = 1;
   ulong bits, bits2;
	int val1, val2;
   
   mpz_init(m1); 
   mpz_init(m2); 
   mpz_init(m3); 
   mpz_init(m4); 

   for (ulong count1 = 0; (count1 < 100*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      F_mpz_init(f3);

      bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      do
		{
			bits2 = z_randint(200) + 1;
         F_mpz_test_random(f2, bits2);
		} while (F_mpz_is_zero(f2));
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
		
		val1 = F_mpz_invert(f3, f1, f2);
		F_mpz_get_mpz(m3, f3);
      val2 = mpz_invert(m4, m1, m2);
          
      result = (((val1 == 0) && (val2 == 0)) || (mpz_cmp(m4, m3) == 0)); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
      F_mpz_clear(f3);
   }
   // check aliasing of operands 1 and 2
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      do
		{
			bits2 = z_randint(200) + 1;
         F_mpz_test_random(f2, bits2);
		} while (F_mpz_is_zero(f2));
           
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		val1 = F_mpz_invert(f1, f1, f2);
		F_mpz_get_mpz(m3, f1);
      val2 = mpz_invert(m4, m1, m2);

      result = (((val1 == 0) && (val2 == 0)) || (mpz_cmp(m4, m3) == 0)); 
		if (!result) 
		{
			printf("val1 = %ld, val2 = %ld\n", val1, val2);
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
   // check aliasing of operands 1 and 3
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      F_mpz_init(f2);
      
		bits = z_randint(200) + 1;
      F_mpz_test_random(f1, bits);
      do
		{
			bits2 = z_randint(200) + 1;
         F_mpz_test_random(f2, bits2);
		} while (F_mpz_is_zero(f2));
           
           
      F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);
      
		val1 = F_mpz_invert(f2, f1, f2);
		F_mpz_get_mpz(m3, f2);
      val2 = mpz_invert(m4, m1, m2);

      result = (((val1 == 0) && (val2 == 0)) || (mpz_cmp(m4, m3) == 0)); 
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, bits2 = %ld, m1 = %Zd, m2 = %Zd, m3 = %Zd, m4 = %Zd\n", bits, bits2, m1, m2, m3, m4);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
   // check aliasing of all operands
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init(f1);
      
		do
		{
			bits = z_randint(200) + 1;
         F_mpz_test_random(f1, bits);
		} while (F_mpz_is_zero(f1));
           
          
      F_mpz_get_mpz(m1, f1);
      
		val1 = F_mpz_invert(f1, f1, f1);
		F_mpz_get_mpz(m2, f1);
      val2 = mpz_invert(m3, m1, m1);

      result = ((val1 == 0) && (val2 == 0));
		if (!result) 
		{
			gmp_printf("Error: bits = %ld, val1 = %Zd, val2 = %Zd\n", bits, val1, val2);
		}
          
      F_mpz_clear(f1);
   }

   mpz_clear(m1);
   mpz_clear(m2);
   mpz_clear(m3);
   mpz_clear(m4);
   
   return result;
}

int test_F_mpz_size()
{
   mpz_t m1;
   F_mpz_t f1;
   int result = 1;
   ulong bits;
	ulong m_limbs, f_limbs;
   
   mpz_init(m1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));

      bits = z_randint(500) + 1;
      F_mpz_test_random(f1, bits);
      
		F_mpz_get_mpz(m1, f1);

		m_limbs = mpz_size(m1);
     
		f_limbs = F_mpz_size(f1);
          
      result = (m_limbs == f_limbs); 
		if (!result) 
		{
			printf("Error: bits = %ld, m_limbs = %ld, f_limbs = %ld\n", bits, m_limbs, f_limbs);
		}
          
      F_mpz_clear(f1);
   }
   
   mpz_clear(m1);
   
   return result;
}

int test_F_mpz_sgn()
{
   mpz_t m1;
   F_mpz_t f1;
   int result = 1;
   ulong bits;
	int m_sign, f_sign;
   
   mpz_init(m1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));

      bits = z_randint(500) + 1;
      F_mpz_test_random(f1, bits);
      
		F_mpz_get_mpz(m1, f1);

		m_sign = mpz_sgn(m1);
     
		f_sign = F_mpz_sgn(f1);
          
      result = (m_sign == f_sign); 
		if (!result) 
		{
			printf("Error: bits = %ld, m_sign = %d, f_sign = %d\n", bits, m_sign, f_sign);
		}
          
      F_mpz_clear(f1);
   }
   
   mpz_clear(m1);
   
   return result;
}

int test_F_mpz_is_one()
{
   mpz_t m1;
   F_mpz_t f1;
   int result = 1;
   ulong bits;
	
   mpz_init(m1); 

   // Test case where f1 is 1
	F_mpz_init2(f1, z_randint(10));

	F_mpz_set_ui(f1, 1L);
	result = (F_mpz_is_one(f1));
	if (!result)
	{
		printf("Error: f1 is %ld\n", F_mpz_get_ui(f1));
	}
   
	F_mpz_clear(f1);
   
   // Case where f1 is not necessarily 1
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));

      bits = z_randint(500) + 1;
      F_mpz_test_random(f1, bits);
      
		F_mpz_get_mpz(m1, f1);

		if (mpz_cmp_ui(m1, 1L)) // if f1 is not 1
		{
         result = (!F_mpz_is_one(f1)); 
		   if (!result) 
		   {
			   printf("Error: bits = %ld = %d\n", bits);
		   }
		}

      F_mpz_clear(f1);
   }
   
   mpz_clear(m1);
   
   return result;
}

int test_F_mpz_is_m1()
{
   mpz_t m1;
   F_mpz_t f1;
   int result = 1;
   ulong bits;
	
   mpz_init(m1); 

   // Test case where f1 is 1
	F_mpz_init2(f1, z_randint(10));

	F_mpz_set_si(f1, -1L);
	result = (F_mpz_is_m1(f1));
	if (!result)
	{
		printf("Error: f1 is %ld\n", F_mpz_get_ui(f1));
	}

	F_mpz_clear(f1);
   
   // Case where f1 is not necessarily 1
	for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));

      bits = z_randint(500) + 1;
      F_mpz_test_random(f1, bits);
      
		F_mpz_get_mpz(m1, f1);

		if (mpz_cmp_si(m1, -1L)) // if f1 is not 1
		{
         result = (!F_mpz_is_m1(f1)); 
		   if (!result) 
		   {
			   printf("Error: bits = %ld = %d\n", bits);
		   }
		}

      F_mpz_clear(f1);
   }
   
   mpz_clear(m1);
   
   return result;
}

int test_F_mpz_cmp()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits1, bits2;
	int m_cmp, f_cmp;
   
   mpz_init(m1); 
   mpz_init(m2); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits1 = z_randint(500) + 1;
      bits2 = z_randint(500) + 1;
      F_mpz_test_random(f1, bits1);
      F_mpz_test_random(f2, bits2);
      
		F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);

		m_cmp = mpz_cmp(m1, m2);
     
		f_cmp = F_mpz_cmp(f1, f2);
          
      result = (((m_cmp == 0) && (f_cmp == 0)) || ((m_cmp > 0) && (f_cmp > 0)) || ((m_cmp < 0) && (f_cmp < 0))); 
		if (!result) 
		{
			printf("Error: bits1 = %ld, bits2 = %ld, m_cmp = %d, f_cmp = %d\n", bits1, bits2, m_cmp, f_cmp);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

int test_F_mpz_cmpabs()
{
   mpz_t m1, m2;
   F_mpz_t f1, f2;
   int result = 1;
   ulong bits1, bits2;
	int m_cmp, f_cmp;
   
   mpz_init(m1); 
   mpz_init(m2); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));
      F_mpz_init2(f2, z_randint(10));

      bits1 = z_randint(500) + 1;
      bits2 = z_randint(500) + 1;
      F_mpz_test_random(f1, bits1);
      F_mpz_test_random(f2, bits2);
      
		F_mpz_get_mpz(m1, f1);
      F_mpz_get_mpz(m2, f2);

		m_cmp = mpz_cmpabs(m1, m2);
     
		f_cmp = F_mpz_cmpabs(f1, f2);
          
      result = (((m_cmp == 0) && (f_cmp == 0)) || ((m_cmp > 0) && (f_cmp > 0)) || ((m_cmp < 0) && (f_cmp < 0))); 
		if (!result) 
		{
			printf("Error: bits1 = %ld, bits2 = %ld, m_cmp = %d, f_cmp = %d\n", bits1, bits2, m_cmp, f_cmp);
		}
          
      F_mpz_clear(f1);
      F_mpz_clear(f2);
   }
   
   mpz_clear(m1);
   mpz_clear(m2);
   
   return result;
}

int test_F_mpz_bits()
{
   mpz_t m1;
   F_mpz_t f1;
   int result = 1;
   ulong bits;
	ulong m_bits, f_bits;
   
   mpz_init(m1); 
   
   for (ulong count1 = 0; (count1 < 100000*ITER) && (result == 1); count1++)
   {
      F_mpz_init2(f1, z_randint(10));

      bits = z_randint(500) + 1;
      F_mpz_test_random(f1, bits);
      
		F_mpz_get_mpz(m1, f1);

		m_bits = mpz_sizeinbase(m1, 2);
     
		f_bits = F_mpz_bits(f1);
          
      result = ((m_bits == f_bits) || ((m_bits == 1) && (f_bits == 0))); 
		if (!result) 
		{
			printf("Error: bits = %ld, m_bits = %ld, f_bits = %ld\n", bits, m_bits, f_bits);
		}
          
      F_mpz_clear(f1);
   }
   
   mpz_clear(m1);
   
   return result;
}

int test_F_mpz_comb_init_clear()
{
   int result = 1;
      
   for (ulong i = 0; (i < 100*ITER) && (result == 1); i++)
   {
      ulong n = random_ulong(10);
      ulong num_primes = (1L << n);
      ulong * primes = (ulong *) flint_heap_alloc(num_primes);
      ulong p = z_nextprime((1UL << (FLINT_BITS-1)) - 10000000L, 0);
      
		for (ulong i = 0; i < num_primes; i++)
      {
	      primes[i] = p;
	      p = z_nextprime(p, 0);
      }

#if DEBUG
      printf("n = %ld, num_primes = %ld\n", n, num_primes);
#endif

      F_mpz_comb_t comb;
      F_mpz_comb_init(comb, primes, num_primes);
      F_mpz_comb_clear(comb);
      flint_heap_free(primes);
   }
      
   return result;
}

int test_F_mpz_multi_CRT_ui_unsigned()
{
   int result = 1;
   F_mpz_t input;
   mpz_t num1;
   ulong * output, * output2;
      
   mpz_init(num1);

   for (ulong i = 0; (i < 10000*ITER) && (result == 1); i++)
   {
      ulong bits = z_randint(300)+1;

#if FLINT_BITS == 32
      double primes_per_limb = 1.0325;
#elif FLINT_BITS == 64
      double primes_per_limb = 1.016;
#endif

	  ulong num_primes = (bits*primes_per_limb)/FLINT_BITS + 1;

     ulong * primes = (ulong *) flint_heap_alloc(num_primes);
     ulong prime = z_nextprime((1UL << (FLINT_BITS-1)) - 10000000L, 0);
      
	  for (ulong j = 0; j < num_primes; j++)
     {
        primes[j] = prime;
        prime = z_nextprime(prime, 0);
     }

	  F_mpz_init(input);

     F_mpz_test_random(input, bits);
     F_mpz_abs(input, input);
	  F_mpz_get_mpz(num1, input);

     output = (ulong *) flint_heap_alloc(num_primes);
     output2 = (ulong *) flint_heap_alloc(num_primes);

     F_mpz_comb_t comb;
     F_mpz_comb_init(comb, primes, num_primes);
     F_mpz_multi_mod_ui(output, input, comb);
      
     F_mpz_t temp;
	  F_mpz_init(temp);

     F_mpz_multi_CRT_ui_unsigned(temp, output, comb);
     result &= F_mpz_equal(temp, input);
      
     for (ulong k = 0; k < num_primes; k++)
     {
        output2[k] = F_mpz_mod_ui(temp, input, primes[k]);
		  result &= (output[k] == output2[k]);
     }

     if (!result)
	  {
		  printf("Error: bits = %ld, num_primes = %ld\n", bits, num_primes);
		  F_mpz_print(temp); printf("\n");
		  F_mpz_print(input); printf("\n");
	  }

	  F_mpz_clear(temp);

     F_mpz_comb_clear(comb);
     F_mpz_clear(input);
     flint_heap_free(output);
     flint_heap_free(output2);
     flint_heap_free(primes);
  }

  mpz_clear(num1);
         
  return result;
}

int test_F_mpz_multi_CRT_ui()
{
   int result = 1;
   F_mpz_t input;
   mpz_t num1;
   ulong * output, * output2;
      
   mpz_init(num1);

   for (ulong i = 0; (i < 10000*ITER) && (result == 1); i++)
   {
      ulong bits = z_randint(300)+1;

#if FLINT_BITS == 32
      double primes_per_limb = 1.0325;
#elif FLINT_BITS == 64
      double primes_per_limb = 1.016;
#endif

	  ulong num_primes = ((bits + 1)*primes_per_limb)/FLINT_BITS + 1;

     ulong * primes = (ulong *) flint_heap_alloc(num_primes);
     ulong prime = z_nextprime((1UL << (FLINT_BITS-1)) - 10000000L, 0);
      
	  for (ulong j = 0; j < num_primes; j++)
     {
        primes[j] = prime;
        prime = z_nextprime(prime, 0);
     }

	  F_mpz_init(input);

     F_mpz_test_random(input, bits);
     F_mpz_get_mpz(num1, input);

     output = (ulong *) flint_heap_alloc(num_primes);
     output2 = (ulong *) flint_heap_alloc(num_primes);

     F_mpz_comb_t comb;
     F_mpz_comb_init(comb, primes, num_primes);
     F_mpz_multi_mod_ui(output, input, comb);
      
     F_mpz_t temp;
	  F_mpz_init(temp);

     F_mpz_multi_CRT_ui(temp, output, comb);
     result &= F_mpz_equal(temp, input);
      
	  if (!result)
	  {
		  printf("Error: bits = %ld, num_primes = %ld\n", bits, num_primes);
		  F_mpz_print(temp); printf("\n");
		  F_mpz_print(input); printf("\n");
	  }

     for (ulong k = 0; (k < num_primes) && (result == 1); k++)
     {
        output2[k] = F_mpz_mod_ui(temp, input, primes[k]);
		  result &= (output[k] == output2[k]);
		  
		  if (!result)
	     {
		     printf("Error: bits = %ld, num_primes = %ld\n", bits, num_primes);
		     printf("Error: k = %ld, output[k] = %ld, output2[k] = %ld\n", k, output[k], output2[k]);
	     }
     }

	  F_mpz_clear(temp);

     F_mpz_comb_clear(comb);
     F_mpz_clear(input);
     flint_heap_free(output);
     flint_heap_free(output2);
     flint_heap_free(primes);
  }

  mpz_clear(num1);
         
  return result;
}

void F_mpz_poly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);
   semaphore_init();

#if TESTFILE
#endif
#pragma omp parallel sections
	{
#pragma omp section
	{
	RUN_TEST(F_mpz_getset_ui); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_getset_si); 
	}
#pragma omp section
	{
	RUN_TEST(F_mpz_getset_mpz); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_getset_limbs); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_get_d_2exp); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_set); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_equal); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_swap); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_neg); 
	}
#pragma omp section
	{
   RUN_TEST(F_mpz_abs); 
   RUN_TEST(F_mpz_add); 
   RUN_TEST(F_mpz_sub); 
   RUN_TEST(F_mpz_mul_ui); 
   RUN_TEST(F_mpz_mul_si); 
   RUN_TEST(F_mpz_mul); 
   RUN_TEST(F_mpz_mul_2exp); 
   RUN_TEST(F_mpz_add_ui); 
   RUN_TEST(F_mpz_sub_ui); 
   RUN_TEST(F_mpz_addmul_ui); 
   RUN_TEST(F_mpz_submul_ui); 
   RUN_TEST(F_mpz_addmul); 
   RUN_TEST(F_mpz_submul); 
   RUN_TEST(F_mpz_mod_ui); 
   RUN_TEST(F_mpz_mod); 
   RUN_TEST(F_mpz_invert); 
   RUN_TEST(F_mpz_size);
	RUN_TEST(F_mpz_sgn);
	RUN_TEST(F_mpz_cmp);
	RUN_TEST(F_mpz_cmpabs);
	RUN_TEST(F_mpz_is_one);
	RUN_TEST(F_mpz_is_m1);
	RUN_TEST(F_mpz_bits); 
   RUN_TEST(F_mpz_comb_init_clear); 
	RUN_TEST(F_mpz_multi_CRT_ui_unsigned);
	RUN_TEST(F_mpz_multi_CRT_ui);
	}
	}
	
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_poly_test_all();
   test_support_cleanup();
	_F_mpz_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


