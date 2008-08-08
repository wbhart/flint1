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

/*int test_z_mulmod32_precomp()
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
      n = random_ulong2((1UL<<bits)-1)+1;
      
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
}*/

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
   unsigned long n, ninv_hi, ninv_lo;
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
   unsigned long n, ninv_hi, ninv_lo;
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
      p = z_nextprime(p);
      
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
      p = z_nextprime(p);
      
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
		m = z_randprime(bits2);
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
         n = z_nextprime(n);
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
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_D_BITS-1)+1;
	  n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprime(res));
   }  
   
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
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_D_BITS-1)+1;
	  n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      mpz_nextprime(mpz_n, mpz_n);
      res = mpz_get_ui(mpz_n);

      result = (z_isprime_pocklington(res, 100));
   }  
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
	  unsigned long bits = z_randint(FLINT_D_BITS/2 - 1)+1;
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
   unsigned long n, n1, n2;

   int result = 1;
   
   for (unsigned long count = 0; (count < 5000000) && (result == 1); count++)
   { 
      do
      {
         n1 = random_ulong(100);
         n2 = random_ulong(100);
         n2 = n2*n2;
      } while ((n1*n2 > 65535) || (n2 == 1));
      
      n = n1*n2;
      
#if DEBUG
      printf("n1 = %ld, n2 = %ld\n", n1, n2);
#endif

      result = !z_issquarefree(n);
   }  

   for (unsigned long count = 0; (count < 500000) && (result == 1); count++)
   { 
      n = 1;
      n1 = 1;
      n2 = random_ulong(999998)+2;
      
      do
      {
         n = n*n1;
         for (unsigned long i = 0; i < random_ulong(3)+1; i++)
         {
            n1 = z_nextprime(n1);
         } 
      } while (n*n1 < n2);
      
#if DEBUG
      printf("%ld\n", n);
#endif

      result = z_issquarefree(n);
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
           
      for (unsigned long j = 0; j < 10; j++)
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
      } while (z_isprime(n));     
      
#if DEBUG
      printf("n = %ld\n");
#endif

      for (unsigned long j = 0; j < 10; j++)
      {
         factor = z_factor_SQUFOF(n);
      }
      
      if (factor) result = (n == factor*(n/factor));

#if DEBUG
      if (!factor) printf("%ld failed to factor\n", n);
      if (factor) printf("%ld factored\n", n);
      if (!result)
      {
         printf("n = %ld\n", n);
         printf("factors = %ld, %ld\n", factor, n/factor);
      }
#endif

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
           
      for (unsigned long j = 0; j < 1; j++)
         factored = z_factor(&factors, orig_n);
      
      if (factored)
      {
         prod = 1;
         for (i = 0; i < factors.num; i++)
         {
            prod *= z_pow(factors.p[i], factors.exp[i]);
         }
      
         result = (prod == orig_n);
      } 
#if DEBUG
      else printf("%ld didn't factor\n", orig_n);
#endif

#if DEBUG      
      if (!result)
      {
         printf("n = %ld: [", orig_n);
         for (i = 0; i < factors.num - 1; i++)
         {
            printf("%ld, %ld; ", factors.p[i], factors.exp[i]);
         }
         printf("%ld, %ld]\n", factors.p[i], factors.exp[i]);
      }
#endif

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
      orig_n = random_ulong(1000000)+1;
      limit = z_intsqrt(orig_n);

      n = z_factor_partial(&factors, orig_n, limit);
      
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
		 n = z_randint(-1L)+1;
#endif
	  } while (z_isprime(n));

	  limit = z_randint(n);
	  cofactor = z_factor_partial(&factors, n, limit);
	  out = 1;

	  for (int j = 0; j < factors.num; j++)
	  {
		 out*=z_pow(factors.p[j], factors.exp[j]);
	  }

	  if (out*cofactor != n || out <= limit) 
	  {
		 result = 0;
#if DEBUG2
		 printf("failed to factor %ld got to %ld\n", n, out);
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
      p = z_nextprime(p);
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


   //RUN_TEST(z_mulmod32_precomp); // Not current available due to lack of proof of code
   RUN_TEST(z_intsqrt);
   RUN_TEST(z_primitive_root);
   RUN_TEST(z_mod_precomp);
   RUN_TEST(z_div2_precomp);
   RUN_TEST(z_mod2_precomp);
   RUN_TEST(z_ll_mod_precomp);
   RUN_TEST(z_mulmod_precomp);
   RUN_TEST(z_mulmod2_precomp);
   RUN_TEST(z_powmod);
   RUN_TEST(z_powmod2);
   RUN_TEST(z_legendre_precomp);
   RUN_TEST(z_jacobi);
   RUN_TEST(z_sqrtmod);
   RUN_TEST(z_cuberootmod);
   RUN_TEST(z_ispseudoprime_fermat);
   RUN_TEST(z_isprime);
   RUN_TEST(z_isprime_pocklington);
   RUN_TEST(z_nextprime);
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


