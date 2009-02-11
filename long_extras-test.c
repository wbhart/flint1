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

long_extras-test.c: Test code for long_extras.c and long_extras.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "test-support.h"
#include "memory-manager.h"

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

#if FLINT_BITS == 64
int test_z_mulmod32_precomp()
{
   uint32_t ninv;
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 10000) && (result == 1); count++)
   { 
      bits = z_randint(32)+1;
      n = z_randbits(bits)+1;
      
      ninv = z_precompute_inverse32(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(n);
         b = random_ulong2(n);
         
         res1 = z_mulmod32_precomp(a, b, n, ninv);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_mul(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG2              
         if (res1 != res2)
         {
            printf("a = %ld, b = %ld, n = %ld, ninv = %ld, res1 = %ld, res2 = %ld\n", a, b, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_b);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}
#endif

int test_z_intsqrt()
{
   unsigned long n;
   unsigned long a, res1, bits;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 20000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS/2-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      a = n*n;
         
      for (unsigned long i = 0; i < 100; i++)
		 res1 = z_intsqrt(a);
         
#if DEBUG2            
      if (res1 != n)
      {
         printf("a = %ld, n = %ld, res1 = %ld\n", a, n, res1);
      }
#endif
         
      result = (res1 == n);
   } 
     
   return result;
}

int test_z_intcuberoot()
{
   unsigned long n;
   unsigned long a, res1, bits;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 5000000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS/3-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      a = n*n*n;

		res1 = z_intcuberoot(a);
         
      if (res1 != n)
      {
         printf("a = %ld, n = %ld, res1 = %ld\n", a, n, res1);
      }
         
      result = (res1 == n);
   } 
     
   for (unsigned long count = 0; (count < 5000000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS/3-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      a = n*n*n;
      ulong a1 = (n + 1)*(n + 1)*(n + 1);
		a = a + z_randint(a1 - a);

		res1 = z_intcuberoot(a);
         
      if (res1 != n)
      {
         printf("a = %ld, n = %ld, res1 = %ld\n", a, n, res1);
      }
         
      result = (res1 == n);
   } 
     
   return result;
}

int test_z_pow()
{
   unsigned long n, res, res2;
   unsigned long bits;

	mpz_t n_mpz;
	mpz_init(n_mpz);
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 1000000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS - 1) + 1;
      n = random_ulong((1UL<<bits) - 1) + 1;
      
      ulong max_exp = 0;
		mpz_set_ui(n_mpz, n);
		while (mpz_size(n_mpz) == 1)
		{
         max_exp++;
			mpz_set_ui(n_mpz, n);
		   mpz_pow_ui(n_mpz, n_mpz, max_exp + 1);
			if (n == 1) break;
		}
      
		ulong exp = z_randint(max_exp + 1);
		res = z_pow(n, exp);
		mpz_set_ui(n_mpz, n);
		mpz_pow_ui(n_mpz, n_mpz, exp);
		res2 = mpz_get_ui(n_mpz);

		result = (res == res2);

      if (!result)
      {
         printf("n, exp, res, res2\n", n, exp, res, res2);
      }
   }
   
	mpz_clear(n_mpz);

   return result;
}

int test_z_gcd()
{
   long a, b, c, res;
	
	ulong bits1, bits2, bits3;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      bits1 = z_randint(FLINT_BITS - 2) + 1;
      bits2 = z_randint(FLINT_BITS - 2) + 1;
      
		do
		{
			a = z_randbits(bits1);
			b = z_randbits(bits2);
		} while (z_gcd(a, b) != 1);

		if (z_randint(2)) a = -a;
      if (z_randint(2)) b = -b;  

		bits3 = FLINT_BITS - 1 - FLINT_MAX(bits1, bits2);
		c = z_randbits(bits3) + 1; 

      res = z_gcd(a*c, b*c);
		result = (res == c);

      if (!result)
      {
         printf("a = %ld, b = %ld, c = %ld, res = %ld\n", a, b, c, res);
      }
   }
   
   return result;
}

int test_z_invert()
{
   unsigned long a, b, res, res2, bits;
   
   int result = 1;

	mpz_t prod;
	mpz_init(prod);
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS - 2) + 2;
      
		do
		{
			a = z_randbits(bits);
         b = z_randprime(bits, 0);
		} while (z_gcd(a, b) != 1);

		a = a % b;

      res = z_invert(a, b);
		mpz_set_ui(prod, a);
		mpz_mul_ui(prod, prod, res);
      
		res2 = mpz_mod_ui(prod, prod, b);
		
		result = (res2 == 1);

      if (!result)
      {
         printf("a = %ld, b = %ld, res = %ld, res2 = %ld\n", a, b, res, res2);
      }
   }
   
	mpz_clear(prod);

   return result;
}

int test_z_gcd_invert()
{
   unsigned long a, b, c, u, res, res2, bits, bits3;
   
   int result = 1;

	mpz_t prod;
	mpz_init(prod);
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS - 3) + 2;
      
		do
		{
			a = z_randbits(bits);
         b = z_randprime(bits, 0);
		} while (z_gcd(a, b) != 1);

		a = a % b;

		bits3 = FLINT_BITS - 1 - bits;
		c = z_randbits(bits3) + 1; 

      res = z_gcd_invert(&u, a*c, b*c);

		mpz_set_ui(prod, a);
		mpz_mul_ui(prod, prod, u);
      
		res2 = mpz_mod_ui(prod, prod, b);
		
		result = ((res == c) && (res2 == 1));

      if (!result)
      {
         printf("a = %ld, b = %ld, c = %ld, res = %ld\n", a, b, c, u, res, res2);
      }
   }

	mpz_clear(prod);
   
   return result;
}

int test_z_xgcd()
{
   long a, b, c, s, t, res, res2;
	ulong bits1, bits2, bits3;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      bits1 = z_randint(FLINT_BITS - 2) + 1;
      bits2 = z_randint(FLINT_BITS - 2) + 1;
      
		do
		{
			a = z_randbits(bits1);
         b = z_randbits(bits2);
		} while (z_gcd(a, b) != 1);
      
		if (z_randint(2)) a = -a;
		if (z_randint(2)) b = -b;

		bits3 = FLINT_BITS - 1 - FLINT_MAX(bits1, bits2);
		c = z_randbits(bits3) + 1; 

      res = z_xgcd(&s, &t, a*c, b*c);
		res2 = s*a*c + t*b*c;
		result = ((res == c) && (res2 == c));

      if (!result)
      {
         printf("a = %ld, b = %ld, c = %ld, res = %ld\n", a, b, c, res, res2);
      }
   }
   
   return result;
}

int test_z_mod_precomp()
{
   double ninv;
   unsigned long n;
   unsigned long a, res1, res2, bits;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 20000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS/2)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(n*n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_mod_precomp(a, n, ninv);
                  
         res2 = a%n;
         
#if DEBUG               
         if (res1 != res2)
         {
            printf("a = %ld, n = %ld, ninv = %lf, res1 = %ld, res2 = %ld\n", a, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   } 
     
   return result;
}

int test_z_div2_precomp()
{
   double ninv;
   unsigned long n, bits;
   unsigned long a, res1, res2;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 20000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
         
      ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         bits = z_randint(FLINT_BITS-1)+1;
         a = random_ulong2((1UL<<bits)-1)+1;
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_div2_precomp(a, n, ninv);
                  
         res2 = a / n;
         
#if DEBUG           
         if (res1 != res2)
         {
            printf("a = %ld, n = %ld, ninv = %lf, res1 = %ld, res2 = %ld\n", a, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   } 
     
   return result;
}

int test_z_mod2_precomp()
{
   double ninv;
   unsigned long n, bits;
   unsigned long a, res1, res2;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 20000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
         
      ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         bits = z_randint(FLINT_BITS-1)+1;
         a = random_ulong2((1UL<<bits)-1)+1;
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_mod2_precomp(a, n, ninv);
                  
         res2 = a % n;
         
#if DEBUG             
         if (res1 != res2)
         {
            printf("a = %ld, n = %ld, ninv = %lf, res1 = %ld, res2 = %ld\n", a, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   } 
     
   return result;
}

int test_z_ll_mod_precomp()
{
   double ninv;
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         bits = z_randint(FLINT_D_BITS-1)+1;
         a = random_ulong2((1UL<<bits)-1)+1;
         b = random_ulong2(-1L);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_ll_mod_precomp(a, b, n, ninv);
                  
         mpz_set_ui(mpz_a, a);
         mpz_mul_2exp(mpz_res, mpz_a, FLINT_BITS);
         mpz_add_ui(mpz_res, mpz_res, b);
         mpz_set_ui(mpz_n, n);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG          
         if (res1 != res2)
         {
            printf("a = %ld, b = %ld, n = %ld, ninv = %e, res1 = %ld, res2 = %ld\n", a, b, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_b);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}

int test_z_addmod()
{
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong2(n);
         b = random_ulong2(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_addmod(a, b, n);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_add(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG              
         if (res1 != res2)
         {
            printf("a = %ld, b = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, b, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_b);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}

int test_z_submod()
{
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong2(n);
         b = random_ulong2(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_submod(a, b, n);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_sub(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG              
         if (res1 != res2)
         {
            printf("a = %ld, b = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, b, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_b);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}

int test_z_negmod()
{
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong2(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_negmod(a, n);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_neg(mpz_res, mpz_a);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG              
         if (res1 != res2)
         {
            printf("a = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, b, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}

int test_z_mulmod_precomp()
{
   double ninv;
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong2(n);
         b = random_ulong2(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_mulmod_precomp(a, b, n, ninv);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_mul(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG              
         if (res1 != res2)
         {
            printf("a = %ld, b = %ld, n = %ld, ninv = %lf, res1 = %ld, res2 = %ld\n", a, b, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_b);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}

int test_z_mulmod2_precomp()
{
   double ninv;
   unsigned long n;
   unsigned long a, b, res1, res2, bits;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1)+1;
      
      ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong2(n);
         b = random_ulong2(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_mulmod2_precomp(a, b, n, ninv);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_mul(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG              
         if (res1 != res2)
         {
            printf("a = %ld, b = %ld, n = %ld, ninv = %lf, res1 = %ld, res2 = %ld\n", a, b, n, ninv, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }
   
   mpz_clear(mpz_a);
   mpz_clear(mpz_b);
   mpz_clear(mpz_n);
   mpz_clear(mpz_res); 
   
   return result;
}

int test_z_powmod()
{
   unsigned long n;
   unsigned long a, exp, res1, res2, bits;
   
   mpz_t mpz_res, mpz_a, mpz_n;
   mpz_init(mpz_res);
   mpz_init(mpz_a);
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(n); 
         bits = z_randint(FLINT_BITS-1)+1;
         exp = random_ulong2((1UL<<bits)-1)+1;
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_powmod(a, exp, n);
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_powm_ui(mpz_res, mpz_a, exp, mpz_n);
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG            
         if (res1 != res2)
         {
            printf("a = %ld, exp = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, exp, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_a);
   mpz_clear(mpz_n); 

   return result;
}

int test_z_powmod2()
{
   unsigned long n;
   unsigned long a, exp, res1, res2, bits;
   
   mpz_t mpz_res, mpz_a, mpz_n;
   mpz_init(mpz_res);
   mpz_init(mpz_a);
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(n); 
         bits = z_randint(FLINT_BITS-1)+1;
         exp = random_ulong2((1UL<<bits)-1)+1;
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_powmod2(a, exp, n);
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_powm_ui(mpz_res, mpz_a, exp, mpz_n);
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG            
         if (res1 != res2)
         {
            printf("a = %ld, exp = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, exp, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_a);
   mpz_clear(mpz_n); 

   return result;
}

int test_z_powmod_precomp()
{
   unsigned long n;
   unsigned long a, exp, res1, res2, bits;
   double ninv;

   mpz_t mpz_res, mpz_a, mpz_n;
   mpz_init(mpz_res);
   mpz_init(mpz_a);
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 

		ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(n); 
         bits = z_randint(FLINT_BITS-1)+1;
         exp = random_ulong2((1UL<<bits)-1)+1;
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_powmod_precomp(a, exp, n, ninv);
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_powm_ui(mpz_res, mpz_a, exp, mpz_n);
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG            
         if (res1 != res2)
         {
            printf("a = %ld, exp = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, exp, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_a);
   mpz_clear(mpz_n); 

   return result;
}

int test_z_powmod2_precomp()
{
   unsigned long n;
   unsigned long a, exp, res1, res2, bits;
   double ninv;

   mpz_t mpz_res, mpz_a, mpz_n;
   mpz_init(mpz_res);
   mpz_init(mpz_a);
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 

		ninv = z_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(n); 
         bits = z_randint(FLINT_BITS-1)+1;
         exp = random_ulong2((1UL<<bits)-1)+1;
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = z_powmod2_precomp(a, exp, n, ninv);
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_powm_ui(mpz_res, mpz_a, exp, mpz_n);
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG            
         if (res1 != res2)
         {
            printf("a = %ld, exp = %ld, n = %ld, res1 = %ld, res2 = %ld\n", a, exp, n, res1, res2);
         }
#endif
         
         result = (res1 == res2);
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_a);
   mpz_clear(mpz_n); 

   return result;
}

int test_z_sqrtmod()
{
   unsigned long p = 0;
   unsigned long a, res1, bits;
   
   mpz_t mpz_res, mpz_p, mpz_temp;
   mpz_init(mpz_res);
   mpz_init(mpz_p);
   mpz_init(mpz_temp);
          
   int result = 1;
   
   for (unsigned long count = 0; (count < 10000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_D_BITS-1)+1;
      p = random_ulong((1UL<<bits)-1UL)+1; 
      p = z_nextprime(p, 0);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(p); 
         
         for (unsigned long count3 = 0; count3 < 1; count3++)   
            res1 = z_sqrtmod(a, p);
            
         if (res1)
         {
            mpz_set_ui(mpz_temp, res1);
            mpz_set_ui(mpz_p, p);
            mpz_powm_ui(mpz_res, mpz_temp, 2UL, mpz_p);
            result &= (mpz_cmp_ui(mpz_res, a) == 0);
         
#if DEBUG
            if (mpz_cmp_ui(mpz_res,a))
            {
               gmp_printf("res1 = %ld, p = %ld, a = %ld, sqrt^2 = %Zd\n", res1, p, a, mpz_res);
            }
#endif
         }
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_p);
   mpz_clear(mpz_temp);
   
   return result;
}

int test_z_cuberootmod()
{
   unsigned long p = 0;
   unsigned long a, res1, res2;
   unsigned long cuberoot1, bits;
   
   mpz_t mpz_res, mpz_p, mpz_temp;
   mpz_init(mpz_res);
   mpz_init(mpz_p);
   mpz_init(mpz_temp);
          
   int result = 1;
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
#if FLINT_BITS == 64
      bits = z_randint(38)+2;
#else 
      bits = z_randint(29)+2;
#endif
      p = random_ulong((1UL<<bits)-1)+3;
      p = z_nextprime(p, 0);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong2(p-1)+1; 
#if DEBUG 
         printf("bits = %ld, p = %ld, a = %ld\n", bits, p, a);
#endif
         
         for (unsigned long count = 0; count < 10; count++)   
            res1 = z_cuberootmod(&cuberoot1, a, p);
         if ((res1 == 0) && (p % 3 == 2) && (a != 0)) result == 0;
         else if (res1)
         {
            mpz_set_ui(mpz_temp, res1);
            mpz_set_ui(mpz_p, p);
            mpz_powm_ui(mpz_res, mpz_temp, 3UL, mpz_p);
            result &= (mpz_cmp_ui(mpz_res,a) == 0);
            mpz_mul_ui(mpz_temp, mpz_temp, cuberoot1);
            mpz_mod(mpz_temp, mpz_temp, mpz_p);
            mpz_powm_ui(mpz_res, mpz_temp, 3UL, mpz_p);
            result &= (mpz_cmp_ui(mpz_res,a) == 0);
            mpz_mul_ui(mpz_temp, mpz_temp, cuberoot1);
            mpz_powm_ui(mpz_res, mpz_temp, 3UL, mpz_p);
            result &= (mpz_cmp_ui(mpz_res,a) == 0);
         
#if DEBUG
            if (mpz_cmp_ui(mpz_res,a))
            {
               gmp_printf("res1 = %ld, p = %ld, a = %ld, cuberoot^3 = %Zd\n", res1, p, a, mpz_res);
            }
#endif
         }
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_p);
   mpz_clear(mpz_temp);
   
   return result;
}

int test_z_jacobi()
{
	unsigned long m, i;
	long a;
	int res1, res2;
	mpz_t a2, m2;
	mpz_init(a2);
	mpz_init(m2);
   int result = 1;

	for (i = 0; (i < 1000000) && (result == 1); i++)
	{
	   ulong bits1 = z_randint(FLINT_BITS);
		ulong bits2 = z_randint(FLINT_BITS-1);

		a = z_randbits(bits1); 
		do {
			m = 2*z_randbits(bits2) + 1;
			if ((long) m < 0L) m = -m;
		} while ((a!= 0) && (z_gcd(a,m) > 1));

		if (z_randint(2)) a = -a;

#if DEBUG
        printf("(%ld/%ld), bits1 = %ld, bits2 = %ld\n", a, m, bits1, bits2);
#endif

		mpz_set_si(a2, a);
		mpz_set_ui(m2, m);

		res1 = z_jacobi(a, m);
		res2 = mpz_jacobi(a2, m2);

		result = (res1 == res2);

#if DEBUG2
		if (!result)
		{
			printf("FAIL (%ld/%ld) FLINT:%d GMP:%d gcd:%ld\n", a, m, res1, res2, z_gcd(a,m));
		}
#endif
	}

	mpz_clear(m2);
	mpz_clear(a2);
	
	return result;
}

int test_z_legendre_precomp()
{
	unsigned long m, i;
	long a;
	int res1, res2;
	mpz_t a2, m2;
	mpz_init(a2);
	mpz_init(m2);
   int result = 1;
	double m_inv;

	for (i = 0; (i < 100000) && (result == 1); i++)
	{
	   ulong bits1 = z_randint(FLINT_BITS);
		ulong bits2 = z_randint(FLINT_BITS-2)+2;

		a = z_randbits(bits1); 
		m = z_randprime(bits2, 0);
	   if (m == 2) m++;
		a = a % m;

#if DEBUG
        printf("(%ld/%ld), bits1 = %ld, bits2 = %ld\n", a, m, bits1, bits2);
#endif

		mpz_set_si(a2, a);
		mpz_set_ui(m2, m);

		m_inv = z_precompute_inverse(m);
		res1 = z_legendre_precomp(a, m, m_inv);
		res2 = mpz_legendre(a2, m2);

		result = (res1 == res2);

#if DEBUG2
		if (!result)
		{
			printf("FAIL (%ld/%ld) FLINT:%d GMP:%d gcd:%ld\n", a, m, res1, res2, z_gcd(a,m));
		}
#endif
	}
	
	mpz_clear(m2);
	mpz_clear(a2);
	
	return result;
}

int test_z_nextprime()
{
   unsigned long n;
   unsigned long res1, res2;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      unsigned long bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      for (unsigned long i = 0; i < 1; i++)
      {
         mpz_nextprime(mpz_n, mpz_n);
         n = z_nextprime(n, 0);
      }
      res1 = n;
      res2 = mpz_get_ui(mpz_n);
#if DEBUG
      if (res1 != res2) printf("res1 = %ld, res2 = %ld\n", res1, res2);
#endif
      result = (res1 == res2);
   }  
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_ispseudoprime_fermat()
{
   unsigned long n;
   unsigned long res;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_D_BITS-1)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      ulong i = z_randint(1000000)+2;
	   if ((i % res) == 0) i++;
	   result = (z_ispseudoprime_fermat(res, i) == 1);

#if DEBUG
      if (!result) printf("i = %ld, n = %ld\n", i, n);
#endif

   }  
   
#define FERMAT_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; count < FERMAT_COUNT; count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		n *= z_randint(bits) + 2;
      ulong i = z_randint(1000000)+2;
	   if ((i % n) == 0) i++;
	   if (z_ispseudoprime_fermat(n, i) != 1) comp++; 
	}
   double frac = (double) comp / (double) FERMAT_COUNT;
	result = (frac > 0.9); 
	if (!result) printf("Error: only %lf%% of composites declared composite\n", frac*100.0);

   mpz_clear(mpz_n); 

   return result;
}

int test_z_ispseudoprime_lucas()
{
   ulong n;
   ulong res;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (ulong count = 0; (count < 100000) && (result == 1); count++)
   { 
	   ulong bits = z_randint(FLINT_BITS - 1) + 1;
      n = random_ulong((1UL<<bits) - 1UL) + 256; 
      mpz_set_ui(mpz_n, n);

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_ispseudoprime_lucas(res) == 1);
		if (!result) printf("Error: prime not pseudoprime, n = %ld\n", n);
   } 

#define LUCAS_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; count < LUCAS_COUNT; count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		n *= z_randint(bits) + 2;
      if (z_ispseudoprime_lucas(n) != 1) comp++; 
	}
   double frac = (double) comp / (double) LUCAS_COUNT;
	result = (frac > 0.99); 
	if (!result) printf("Error: only %lf%% of composites declared composite\n", frac*100.0);

   mpz_clear(mpz_n); 

   return result;
}

int test_z_ispseudoprime_lucas_ab()
{
   ulong n;
   ulong res;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (ulong count = 0; (count < 100000) && (result == 1); count++)
   { 
	   ulong bits = z_randint(FLINT_BITS - 1) + 1;
      n = random_ulong((1UL<<bits) - 1UL) + 256; 
      mpz_set_ui(mpz_n, n);

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);
		int a = 1;
      int b = 1;

      result = (z_ispseudoprime_lucas_ab(res, a, b) == 1);
		if (!result) printf("Error: prime not pseudoprime, n = %ld, a = %d, b = %d\n", n, a, b);
   } 

#define LUCAS2_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; count < LUCAS2_COUNT; count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		n *= z_randint(bits) + 2;
      int a = 1;
      int b = 1;
      if (z_ispseudoprime_lucas_ab(n, a, b) != 1) comp++; 
	}
   double frac = (double) comp / (double) LUCAS2_COUNT;
	result = (frac > 0.50); 
	if (!result) printf("Error: only %lf%% of composites declared composite\n", frac*100.0);

   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprobab_prime()
{
   unsigned long n;
   unsigned long res;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	   unsigned long bits = z_randint(FLINT_BITS-2)+2;
	   n = random_ulong((1UL<<bits)-2UL)+2; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprobab_prime(res));
   }  

#define PRIME_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; (count < PRIME_COUNT) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		n *= z_randint(bits) + 2;
      if (z_isprobab_prime(n)) result = 0; 
	}
   
	if (!result) printf("Error : n = %ld\n", n);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprobab_prime_precomp()
{
   unsigned long n;
   unsigned long res;
   double ninv;

   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-2)+2;
	  n = z_randint((1UL<<bits)-2UL)+2; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      ninv = z_precompute_inverse(res);
		result = (z_isprobab_prime_precomp(res, ninv));
		if (!result) printf("Error: %ld is reported composite!\n", n, res);
   }  

#define PRIME_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; (count < PRIME_COUNT) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		if (n == 2) n++;
		ulong n1 = z_randint(bits) + 2;
		if ((n1 & 1L) == 0) n1++;
		n *= n1;
      ninv = z_precompute_inverse(n);
	   if (z_isprobab_prime_precomp(n, ninv)) result = 0; 
	}
   if (!result) printf("Error: n = %ld\n", n);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_ispseudoprime_fibonacci_precomp()
{
   unsigned long n;
   unsigned long res;
   double ninv;

   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-2)+2;
	  n = z_randint((1UL<<bits)-3UL)+3; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      ninv = z_precompute_inverse(res);
		result = ((res == 5) || (z_ispseudoprime_fibonacci_precomp(res, ninv)));
		if (!result) printf("Error: %ld is reported composite!\n", n, res);
   }  

#define PRIME_COUNT 100000
	
	ulong pseudo = 0;
	for (ulong count = 0; (count < PRIME_COUNT) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		if (n == 2) n++;
		ulong n1 = z_randint(bits) + 2;
		if ((n1 & 1L) == 0) n1++;
		n *= n1;
		ninv = z_precompute_inverse(n);
	   if ((n % 5) != 0) 
			if (z_ispseudoprime_fibonacci_precomp(n, ninv)) pseudo++; 
	}

	result = (pseudo < 50);

   if (!result) printf("Error: %ld pseudoprimes in PRIME_COUNT\n", pseudo);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprobab_prime_BPSW()
{
   unsigned long n;
   unsigned long res;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	   unsigned long bits = z_randint(FLINT_BITS-2)+2;
	   n = random_ulong((1UL<<bits)-2UL)+2; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = ((res < 13) || (z_isprobab_prime_BPSW(res)));
   }  
	
	if (!result) printf("Error : n = %ld\n", res);

#define PRIME_COUNT2 100000
	
	ulong comp = 0;
	for (ulong count = 0; (count < PRIME_COUNT2) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		n *= z_randint(bits) + 2;
      if (z_isprobab_prime_BPSW(n)) result = 0; 
	}
   
	if (!result) printf("Error : n = %ld\n", n);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_miller_rabin_precomp()
{
   unsigned long n;
   unsigned long res;
   double ninv;
	ulong pseudo = 0;

   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-2)+2;
	  n = z_randint((1UL<<bits)-2UL)+2; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      ninv = z_precompute_inverse(res);
		result = (z_miller_rabin_precomp(res, ninv, 5));
		if (!result) printf("Error: %ld is reported composite!\n", n, res);
   }  

#define PRIME_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; (count < PRIME_COUNT) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		if (n == 2) n++;
		ulong n1 = z_randint(bits) + 2;
		if ((n1 & 1L) == 0) n1++;
		n *= n1;
      ninv = z_precompute_inverse(n);
	   if (z_miller_rabin_precomp(n, ninv, 5)) pseudo++; 
	}
	result = (pseudo < 10);
   if (!result) printf("Error: %ld composites declared prime!\n", pseudo);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprime()
{
   unsigned long n;
   unsigned long res;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 10000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-1)+1;
	  n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprime(res));
   }  

	for (unsigned long count = 0; (count < 10000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-2)+2;
	  n = random_ulong((1UL<<bits)-2UL)+2; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprime(res));
   }  

#define PRIME_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; (count < PRIME_COUNT) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		n *= z_randint(bits) + 2;
      if (z_isprime(n)) result = 0; 
	}
   
	if (!result) printf("Error : n = %ld\n", n);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprime_precomp()
{
   unsigned long n;
   unsigned long res;
   double ninv;

   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 10000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-2)+2;
	  n = z_randint((1UL<<bits)-2UL)+2; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      ninv = z_precompute_inverse(res);
		result = (z_isprime_precomp(res, ninv));
	   if (!result) printf("Error: %ld is reported composite!\n", res);
   }

   for (unsigned long count = 0; (count < 10000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_BITS-2)+2;
	  n = z_randint((1UL<<bits)-2UL)+2; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      ninv = z_precompute_inverse(res);
		result = (z_isprime_precomp(res, ninv));
	   if (!result) printf("Error: %ld is reported composite!\n", res);
   }  

#define PRIME_COUNT 100000
	
	ulong comp = 0;
	for (ulong count = 0; (count < PRIME_COUNT) && (result == 1); count++)
   { 
		ulong bits = z_randint(FLINT_BITS/2 - 2) + 2;
      n = z_randprime(bits, 0); 
		if (n == 2) n++;
		ulong n1 = z_randint(bits) + 2;
		if ((n1 & 1L) == 0L) n1++;
		n *= n1;
      ninv = z_precompute_inverse(n);
	   if (z_isprime_precomp(n, ninv)) result = 0; 
	}
   if (!result) printf("Error: n = %ld\n", n);
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprime_pocklington()
{
   unsigned long n;
   unsigned long res, res2;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 10000L) && (result == 1); count++)
   { 
	   unsigned long bits = z_randint(FLINT_BITS-2)+2;
	   n = random_ulong((1UL<<bits)-2UL)+2; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprime_pocklington(res, 100));
		if (result != 1) 
		{
			printf("%d\n", result);
			result = 0;
		}
   }  
   
	for (unsigned long count = 0; (count < 100000L) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint((FLINT_BITS-1)/2 - 1)+1;
	  n = random_ulong((1UL<<bits)-1UL)+1; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

	   n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res2 = mpz_get_ui(mpz_n);

      int resval = z_isprime_pocklington(res*res2, 100);

#if DEBUG2
	  if (resval < 0L) printf("Failure %ld, %ld, %d\n", res, res2, resval);
#endif

	  result = (!resval);

#if DEBUG2
	  if (!result) printf("res = %ld, res2 = %ld\n", res, res2);
#endif

   }  
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_isprime_nm1()
{
   unsigned long n;
   unsigned long res, res2;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 10000L) && (result == 1); count++)
   { 
	   unsigned long bits = z_randint(FLINT_BITS-2)+2;
	   n = random_ulong((1UL<<bits)-2UL)+2; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprime_nm1(res, 200));
		if (result != 1) 
		{
			printf("%d, %ld\n", result, res);
			result = 0;
		}
   }  
   
	for (unsigned long count = 0; (count < 100000L) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint((FLINT_BITS-1)/2 - 1)+1;
	  n = random_ulong((1UL<<bits)-1UL)+1; 
     mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

	   n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res2 = mpz_get_ui(mpz_n);

      int resval = z_isprime_nm1(res*res2, 200);

#if DEBUG2
	  if (resval < 0L) printf("Failure %ld, %ld, %d\n", res, res2, resval);
#endif

	  result = (!resval);

#if DEBUG2
	  if (!result) printf("res = %ld, res2 = %ld\n", res, res2);
#endif

   }  
   
   mpz_clear(mpz_n); 

   return result;
}

int test_z_remove()
{
   unsigned long n;
   unsigned long bits;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 5000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS - 1) + 1;
      n = random_ulong((1UL<<bits) - 1) + 1;
      
      ulong p = 2;
		for (unsigned long count = 0; (count < 200) && (result == 1); count++)
		{
         ulong oldn = n;
		   int exp = z_remove(&n, p);
			result &= (oldn == n*z_pow(p, exp));
			p = z_nextprime(p, 0);
		}
                  
      if (!result)
      {
         printf("n = %ld, p = %ld\n", n, p);
      }

		p = 2;
	   for (unsigned long count = 0; (count < 200) && (result == 1); count++)
		{
         result &= ((n % p) != 0);
		   p = z_nextprime(p, 0);
		}
                                         
      if (!result)
      {
         printf("n = %ld, p = %ld\n", n, p);
      }
   }
   
   return result;
}

int test_z_remove_precomp()
{
   unsigned long n;
   unsigned long bits;
	double pinv;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 5000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS - 1) + 1;
      n = random_ulong((1UL<<bits) - 1) + 1;
      
      ulong p = 2;
	   for (unsigned long count = 0; (count < 200) && (result == 1); count++)
		{
         ulong oldn = n;
			pinv = z_precompute_inverse(p);
			int exp = z_remove_precomp(&n, p, pinv);
			result &= (oldn == n*z_pow(p, exp));
			p = z_nextprime(p, 0);
		}
                  
      if (!result)
      {
         printf("n = %ld, p = %ld\n", n, p);
      }

		p = 2;
	   for (unsigned long count = 0; (count < 200) && (result == 1); count++)
		{
         result &= ((n % p) != 0);
			p = z_nextprime(p, 0);
		}
                                         
      if (!result)
      {
         printf("n = %ld, p = %ld\n", n, p);
      }
   }
   
   return result;
}

int test_z_CRT()
{
   unsigned long x1, x2, n1, n2;
   unsigned long res;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 500000) && (result == 1); count++)
   { 
      unsigned long bits = z_randint(FLINT_BITS-2) + 2;
      unsigned long bits1 = z_randint(bits-1) + 1;
      unsigned long bits2 = bits - bits1;
      
      n1 = random_ulong((1UL<<bits1)-1UL)+1; 
      do n2 = random_ulong((1UL<<bits2)-1UL)+1;
      while (z_gcd(n1, n2) != 1); 
      
      x1 = random_ulong2(n1);
      x2 = random_ulong2(n2);
      
#if DEBUG
      printf("x1 = %ld, n1 = %ld, x2 = %ld, n2 = %ld\n", x1, n1, x2, n2);
#endif

      for (unsigned long i = 0; i < 10; i++)
      {
         res = z_CRT(x1, n1, x2, n2);
      }
      result = (((res % n1) == x1) && ((res % n2) == x2));
      
#if DEBUG
      if (!result) printf("res = %ld\n", res);
#endif
   }  
   
   return result;
}

int test_z_issquarefree()
{
   unsigned long n, n1, n2, n3;

   int result = 1;
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      do
      {
         n1 = random_ulong(1L<<(FLINT_BITS/3));
         n2 = random_ulong(1L<<(FLINT_BITS/3));
         n2 = n2*n2;
      } while (n2 == 1);
      
      n = n1*n2;
      
#if DEBUG
      printf("n1 = %ld, n2 = %ld\n", n1, n2);
#endif

      result = !z_issquarefree(n, 1);
   }  

   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      n = 1;
      n1 = 1;
      n2 = random_ulong(999998)+2;
      
      do
      {
         n = n*n1;
         for (unsigned long i = 0; i < random_ulong(3)+1; i++)
         {
            n1 = z_nextprime(n1, 0);
         } 
      } while (n*n1 < n2);
      
#if DEBUG
      printf("%ld\n", n);
#endif

      result = z_issquarefree(n, 1);
   }  
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      n1 = z_randprime(FLINT_BITS/3, 0);
		do {n2 = z_randprime(FLINT_BITS/3, 0);} while (n2 == n1);
      do {n3 = z_randprime(FLINT_BITS/3, 0);} while ((n3 == n1) || (n3 == n2));
      
      n = n1*n2*n3;
      
#if DEBUG
      printf("%ld\n", n);
#endif

      result = z_issquarefree(n, 1);
   }  
   
   return result;
}

int test_z_factor_trial()
{
   unsigned long n, prod, orig_n;
   factor_t factors;
   int i;

   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      orig_n = random_ulong(1000000)+1;
           
      n = z_factor_trial(&factors, orig_n);
      
      prod = n;
      for (i = 0; i < factors.num; i++)
      {
          prod *= z_pow(factors.p[i], factors.exp[i]);
      }
      
      result = (prod == orig_n);

#if DEBUG
      if (!result)
      {
         printf("n = %ld: [", orig_n);
         for (i = 0; i < factors.num - 1; i++)
         {
            printf("%ld, %ld; ", factors.p[i], factors.exp[i]);
         }
         printf("%ld, %ld", factors.p[i], factors.exp[i]);
         if (n != 1) printf("; %ld, 1]\n", n);
         else printf("]\n");
      }
#endif

   }  
   
   return result;
}

int test_z_factor_SQUFOF()
{
   unsigned long n, factor, bits;

   int result = 1;
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      do 
      {
         bits = z_randint(FLINT_BITS - 1)+1;
         n = random_ulong((1UL<<bits)-1)+3;
         n|=1;
      } while (z_isprobab_prime(n));     
      
#if DEBUG
      printf("n = %ld\n");
#endif

      factor = z_factor_SQUFOF(n);
      
      result = (n == factor*(n/factor));

      if (!result)
      {
         printf("n = %ld\n", n);
         printf("factors = %ld, %ld\n", factor, n/factor);
      }
   }  
   
   return result;
}

int test_z_factor()
{
   unsigned long n, prod, orig_n, bits;
   factor_t factors;
   int i, factored;

   int result = 1;
   
   for (unsigned long count = 0; (count < 5000) && (result == 1); count++)
   { 
      bits = z_randint(FLINT_BITS-1)+1;
      orig_n = random_ulong((1UL<<bits)-1)+2;


#if DEBUG
      printf("Factoring n = %ld....\n", orig_n);
#endif
           
      z_factor(&factors, orig_n, 1);
      
      prod = 1;
      for (i = 0; i < factors.num; i++)
      {
         prod *= z_pow(factors.p[i], factors.exp[i]);
      }
      
      result = (prod == orig_n);
 
      if (!result)
      {
         printf("n = %ld: [", orig_n);
         for (i = 0; i < factors.num - 1; i++)
         {
            printf("%ld, %ld; ", factors.p[i], factors.exp[i]);
         }
         printf("%ld, %ld]\n", factors.p[i], factors.exp[i]);
      }
   }  
   
   return result;
}

int test_z_factor_partial()
{
   unsigned long n, prod, cofactor, out, orig_n, limit;
   factor_t factors;
   int i;

   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); )
   { 
      orig_n = random_ulong(1000000)+2;
      limit = z_intsqrt(orig_n);

      n = z_factor_partial(&factors, orig_n, limit, 1);
      
      prod = 1;
      for (i = 0; i < factors.num; i++)
      {
          prod *= z_pow(factors.p[i], factors.exp[i]);
      }
      
	  if (n*prod != orig_n) result = 0;
	  if (prod <= limit) result = 0;
	  count++;

#if DEBUG2
      if (!result)
      {
         printf("n = %ld: [", orig_n);
         for (i = 0; i < factors.num - 1; i++)
         {
            printf("%ld, %ld; ", factors.p[i], factors.exp[i]);
         }
         printf("%ld, %ld", factors.p[i], factors.exp[i]);
         printf("]\nprod = %ld, cofactor = %ld, limit = %ld\n", prod, n, limit);
      }
#endif

   }  
   
   for (int i = 0; (i < 10000) && (result == 1); i++)
   {
	  do {
#if FLINT_BITS == 64
		 n = z_randint(z_pow(10, 17))+1;
#else
		 n = z_randint(z_pow(2, 28))+1;
#endif
	  } while (z_isprobab_prime(n));

	  limit = z_randint(n);
	  cofactor = z_factor_partial(&factors, n, limit, 1);
	  out = 1;

	  for (int j = 0; j < factors.num; j++)
	  {
		 out*=z_pow(factors.p[j], factors.exp[j]);
	  }

	  if (out*cofactor != n || out <= limit)
	  {
		 result = 0;
#if DEBUG2
		 printf("failed to factor %ld got to %ld, limit = %ld\n", n, out, limit);
		 for (int j = 0; j < factors.num; j++)
		 {
			printf("%ld^%ld ", factors.p[j], factors.exp[j]);
		 }
#endif
	  }
   }

   return result;
}

int test_z_primitive_root()
{
   unsigned long p = 2;
   unsigned long r;
   
   for(int i = 0; i < 100000; i++) {
      p = z_nextprime(p, 0);
      r = z_primitive_root(p);
      if(r == 0) {
         printf("Fails on p=%d\n", p);
         return 0;
      }
#if DEBUG
      printf("primitive root of p=%d is %d\n", p, r);
#endif
   }
   return 1;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;


#if FLINT_BITS == 64
	RUN_TEST(z_mulmod32_precomp); 
#endif
	RUN_TEST(z_intsqrt);
   RUN_TEST(z_intcuberoot);
   RUN_TEST(z_pow);
   RUN_TEST(z_gcd);
   RUN_TEST(z_invert);
   RUN_TEST(z_gcd_invert);
   RUN_TEST(z_xgcd);
   RUN_TEST(z_primitive_root);
   RUN_TEST(z_mod_precomp);
   RUN_TEST(z_div2_precomp);
   RUN_TEST(z_mod2_precomp);
   RUN_TEST(z_ll_mod_precomp);
   RUN_TEST(z_addmod);
   RUN_TEST(z_submod);
   RUN_TEST(z_negmod);
   RUN_TEST(z_mulmod_precomp);
   RUN_TEST(z_mulmod2_precomp);
   RUN_TEST(z_powmod);
   RUN_TEST(z_powmod2);
   RUN_TEST(z_powmod_precomp);
   RUN_TEST(z_powmod2_precomp);
   RUN_TEST(z_legendre_precomp);
   RUN_TEST(z_jacobi);
   RUN_TEST(z_sqrtmod);
   RUN_TEST(z_cuberootmod);
   RUN_TEST(z_ispseudoprime_fermat);
   RUN_TEST(z_ispseudoprime_lucas);
   RUN_TEST(z_ispseudoprime_lucas_ab);
	RUN_TEST(z_ispseudoprime_fibonacci_precomp);
   RUN_TEST(z_isprobab_prime);
   RUN_TEST(z_isprobab_prime_precomp);
   RUN_TEST(z_isprobab_prime_BPSW);
   RUN_TEST(z_isprime);
   RUN_TEST(z_isprime_precomp);
   RUN_TEST(z_isprime_pocklington);
   RUN_TEST(z_isprime_nm1);
   RUN_TEST(z_miller_rabin_precomp);
   RUN_TEST(z_nextprime);
	RUN_TEST(z_remove);
   RUN_TEST(z_remove_precomp);
   RUN_TEST(z_CRT);
   RUN_TEST(z_issquarefree);
   RUN_TEST(z_factor_trial);
   RUN_TEST(z_factor_SQUFOF);
   RUN_TEST(z_factor);
   RUN_TEST(z_factor_partial);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   fmpz_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


