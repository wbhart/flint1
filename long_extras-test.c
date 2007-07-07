/****************************************************************************

fmpz_poly-test.c: Test code for fmpz_poly.c and fmpz_poly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "test-support.h"

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");

int test_long_mod_precomp()
{
   double ninv;
   unsigned long n;
   unsigned long a, res1, res2;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 20000) && (result == 1); count++)
   { 
      n = random_ulong(4503599627370494UL)+1;
      
      ninv = long_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(18446744073709551615UL);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = long_mod_precomp(a, n, ninv);
                  
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

int test_long_mulmod_precomp()
{
   double ninv;
   unsigned long n;
   unsigned long a, b, res1, res2;
   
   int result = 1;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 

   for (unsigned long count = 0; (count < 2000) && (result == 1); count++)
   { 
      n = random_ulong(4503599627370494UL)+1;
      
      ninv = long_precompute_inverse(n);
      
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong(n);
         b = random_ulong(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = long_mulmod_precomp(a, b, n, ninv);
                  
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_mul(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG2              
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

int test_long_powmod0()
{
   unsigned long n, ninv_hi, ninv_lo;
   unsigned long a, exp, res1, res2;
   
   mpz_t mpz_res, mpz_a, mpz_n;
   mpz_init(mpz_res);
   mpz_init(mpz_a);
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 200) && (result == 1); count++)
   { 
      n = random_ulong(4503599627370494UL)+1; 
      if (!n) n++;
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(n); 
         exp = random_ulong(9223372036854775808UL); 
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = long_powmod0(a, exp, n);
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_powm_ui(mpz_res, mpz_a, exp, mpz_n);
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG2               
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
 
int test_long_mulmod_precomp2()
{
   unsigned long n, ninv_hi, ninv_lo;
   unsigned long a, b, res1, res2;
   
   mpz_t mpz_a, mpz_b, mpz_n, mpz_res;
   mpz_init(mpz_a);
   mpz_init(mpz_b);
   mpz_init(mpz_n);
   mpz_init(mpz_res); 
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      n = random_ulong(18446744073709551614UL)+1;
      
      long_precompute_inverse2(&ninv_hi, &ninv_lo, n);
      for (unsigned long count2 = 0; (count2 < 1000) && (result == 1); count2++)
      {
         a = random_ulong(n);
         b = random_ulong(n);
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = long_mulmod_precomp2(a, b, n, ninv_hi, ninv_lo);
            
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_b, b);
         mpz_set_ui(mpz_n, n);
         mpz_mul(mpz_res, mpz_a, mpz_b);
         mpz_mod(mpz_res, mpz_res, mpz_n);       
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG              
         if (res1 != res2)
         {
            printf("ninv_hi = %ld, ninv_lo = %ld\n", ninv_hi, ninv_lo);
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

int test_long_powmod()
{
   unsigned long n, ninv_hi, ninv_lo;
   unsigned long a, exp, res1, res2;
   
   mpz_t mpz_res, mpz_a, mpz_n;
   mpz_init(mpz_res);
   mpz_init(mpz_a);
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 200) && (result == 1); count++)
   { 
      n = random_ulong(18446744073709551614UL)+1; 
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(n); 
         exp = random_ulong(9223372036854775807UL); 
         
         for (unsigned long count = 0; count < 100; count++)   
            res1 = long_powmod(a, exp, n);
         mpz_set_ui(mpz_a, a);
         mpz_set_ui(mpz_n, n);
         mpz_powm_ui(mpz_res, mpz_a, exp, mpz_n);
         res2 = mpz_get_ui(mpz_res);
         
#if DEBUG2               
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

int test_long_cuberootmod()
{
   unsigned long p = 0;
   unsigned long a, res1, res2;
   unsigned long cuberoot1;
   
   mpz_t mpz_res, mpz_p, mpz_temp;
   mpz_init(mpz_res);
   mpz_init(mpz_p);
   mpz_init(mpz_temp);
          
   int result = 1;
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      p = random_ulong(1099511627776); 
      p = long_nextprime(p);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(p); 
         
         for (unsigned long count = 0; count < 10; count++)   
            res1 = long_cuberootmod(&cuberoot1, a, p);
         
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

int test_long_sqrtmod()
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
      bits = long_randint(63)+1;
      p = random_ulong((1UL<<bits)-1UL)+1; 
      p = long_nextprime(p);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(p); 
         
         for (unsigned long count = 0; count < 1; count++)   
            res1 = long_sqrtmod(a, p);
            
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

int test_long_sqrtmod0()
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
      bits = long_randint(52)+1;
      p = random_ulong((1UL<<bits)-1UL)+1; 
      p = long_nextprime(p);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(p); 
         
         for (unsigned long count = 0; count < 1; count++)   
            res1 = long_sqrtmod0(a, p);
            
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

int test_long_nextprime()
{
   unsigned long n;
   unsigned long res1, res2;
   
   mpz_t mpz_n;
   mpz_init(mpz_n);
       
   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      unsigned long bits = long_randint(63)+1;
      n = random_ulong((1UL<<bits)-1UL)+1; 
      mpz_set_ui(mpz_n, n);

#if DEBUG
      printf("n = %ld\n", n);
#endif

      for (unsigned long i = 0; i < 1; i++)
      {
         mpz_nextprime(mpz_n, mpz_n);
         n = long_nextprime(n);
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

int test_long_CRT()
{
   unsigned long x1, x2, n1, n2;
   unsigned long res;
   
   int result = 1;
   
   for (unsigned long count = 0; (count < 500000) && (result == 1); count++)
   { 
      unsigned long bits = long_randint(62)+2;
      unsigned long bits1 = long_randint(bits-1) + 1;
      unsigned long bits2 = bits - bits1;
      
      n1 = random_ulong((1UL<<bits1)-1UL)+1; 
      do n2 = random_ulong((1UL<<bits2)-1UL)+1;
      while (long_gcd(n1, n2) != 1); 
      
      x1 = random_ulong(n1);
      x2 = random_ulong(n2);
      
#if DEBUG
      printf("x1 = %ld, n1 = %ld, x2 = %ld, n2 = %ld\n", x1, n1, x2, n2);
#endif

      for (unsigned long i = 0; i < 10; i++)
      {
         res = long_CRT(x1, x2, n1, n2);
      }
      result = (((res % n1) == x1) && ((res % n2) == x2));
      
#if DEBUG
      if (!result) printf("res = %ld\n", res);
#endif
   }  
   
   return result;
}

int test_long_issquarefree()
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

      result = !long_issquarefree(n);
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
            n1 = long_nextprime(n1);
         } 
      } while (n*n1 < n2);
      
#if DEBUG
      printf("%ld\n", n);
#endif

      result = long_issquarefree(n);
   }  
   
   return result;
}

int test_long_factor_trial()
{
   unsigned long n, prod, orig_n;
   factor_t factors;
   int i;

   int result = 1;
   
   for (unsigned long count = 0; (count < 100000) && (result == 1); count++)
   { 
      orig_n = random_ulong(1000000);
           
      for (unsigned long j = 0; j < 10; j++)
         n = long_factor_trial(&factors, orig_n);
      
      prod = n;
      for (i = 0; i < factors.num; i++)
      {
          prod *= long_pow(factors.p[i], factors.exp[i]);
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

int test_long_factor_SQUFOF()
{
   unsigned long n, factor;

   int result = 1;
   
   for (unsigned long count = 0; (count < 5000) && (result == 1); count++)
   { 
      do 
      {
         n = random_ulong(1000000000000)+1;
         n|=1;
      } while (long_isprime(n));     
      
      for (unsigned long j = 0; j < 10; j++)
      {
         factor = long_factor_SQUFOF(n);
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

int test_long_factor()
{
   unsigned long n, prod, orig_n;
   factor_t factors;
   int i, factored;

   int result = 1;
   
   for (unsigned long count = 0; (count < 5000) && (result == 1); count++)
   { 
      orig_n = random_ulong(1000000000000000000);

#if DEBUG
         printf("Factoring n = %ld....\n", orig_n);
#endif
           
      for (unsigned long j = 0; j < 1; j++)
         factored = long_factor(&factors, orig_n);
      
      if (factored)
      {
         prod = 1;
         for (i = 0; i < factors.num; i++)
         {
            prod *= long_pow(factors.p[i], factors.exp[i]);
         }
      
         result = (prod == orig_n);
      } else printf("%ld didn't factor\n", orig_n);

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

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(long_mod_precomp);
   RUN_TEST(long_mulmod_precomp);
   RUN_TEST(long_powmod0);
   RUN_TEST(long_sqrtmod0);
   RUN_TEST(long_mulmod_precomp2);
   RUN_TEST(long_powmod);
   RUN_TEST(long_sqrtmod);
   RUN_TEST(long_cuberootmod);
   RUN_TEST(long_nextprime);
   RUN_TEST(long_CRT);
   RUN_TEST(long_issquarefree);
   RUN_TEST(long_factor_trial);
   RUN_TEST(long_factor_SQUFOF);
   RUN_TEST(long_factor);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   fmpz_poly_test_all();
   test_support_cleanup();

   return 0;
}


