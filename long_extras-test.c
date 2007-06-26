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
      n = random_ulong(7742837468472453352UL)+1;
      
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
         
#if DEBUG2               
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
   
   for (unsigned long count = 0; (count < 100) && (result == 1); count++)
   { 
      n = random_ulong(8428374684722343453); 
      if (!n) n++;
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(n); 
         exp = random_ulong(n); 
         
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
   
   mpz_t mpz_res, mpz_p;
   mpz_init(mpz_res);
   mpz_init(mpz_p);
          
   int result = 1;
   
   for (unsigned long count = 0; (count < 1000) && (result == 1); count++)
   { 
      p = random_ulong(1099511627776); 
      p = long_nextprime(p);
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      {
         a = random_ulong(p); 
         
         for (unsigned long count = 0; count < 10; count++)   
            res1 = long_cuberootmod(a, p);
         
         
         if ((res1 == 0) && (p % 3 == 2) && (a != 0)) result == 0;
         else if (res1)
         {
            mpz_set_ui(mpz_res, res1);
            mpz_set_ui(mpz_p, p);
            mpz_powm_ui(mpz_res, mpz_res, 3UL, mpz_p);
         
#if DEBUG
            if (mpz_cmp_ui(mpz_res,a))
            {
               gmp_printf("res1 = %ld, p = %ld, a = %ld, cuberoot^3 = %Zd\n", res1, p, a, mpz_res);
            }
#endif
            result = (mpz_cmp_ui(mpz_res,a) == 0);
         }
      }
   }  
   
   mpz_clear(mpz_res);
   mpz_clear(mpz_p);
   
   return result;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(long_mulmod_precomp2);
   RUN_TEST(long_powmod);
   RUN_TEST(long_cuberootmod);
   
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


