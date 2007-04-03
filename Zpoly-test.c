/****************************************************************************

Zpoly-test.c: Test code for Zpoly.c and Zpoly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "Zpoly.h"
#include "flint-manager.h"
#include <stdio.h>
#include <string.h>


gmp_randstate_t Zpoly_test_randstate;


// tests whether the given polynomial is equal to the one given by the string
// (only for testing purposes in this file)
int Zpoly_equal_str(Zpoly_t poly, char* s)
{
   Zpoly_t poly2;
   Zpoly_init(poly2);
   Zpoly_set_from_string(poly2, s);
   int result = Zpoly_equal(poly, poly2);
   Zpoly_clear(poly2);
   return result;
}



// all test functions return 1 on success, 0 on failure


int test_Zpoly_get_coeff()
{
   int success = 1;
   Zpoly_t poly;
   mpz_t x;

   Zpoly_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_ui(poly->coeffs[0], 47);
   mpz_set_ui(poly->coeffs[1], 48);
   mpz_set_ui(poly->coeffs[2], 49);
   
   Zpoly_get_coeff(x, poly, 0);
   success = success && !mpz_cmp_ui(x, 47);
   Zpoly_get_coeff(x, poly, 1);
   success = success && !mpz_cmp_ui(x, 48);
   Zpoly_get_coeff(x, poly, 2);
   success = success && !mpz_cmp_ui(x, 0);
   
   Zpoly_clear(poly);
   mpz_clear(x);
   return success;
}


int test_Zpoly_get_coeff_ui()
{
   int success = 1;
   Zpoly_t poly;
   mpz_t x;

   Zpoly_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_ui(poly->coeffs[0], 47);
   mpz_set_ui(poly->coeffs[1], 48);
   mpz_set_ui(poly->coeffs[2], 49);
   
   success = success && (Zpoly_get_coeff_ui(poly, 0) == 47);
   success = success && (Zpoly_get_coeff_ui(poly, 1) == 48);
   success = success && (Zpoly_get_coeff_ui(poly, 2) == 0);
   
   Zpoly_clear(poly);
   mpz_clear(x);
   return success;
}


int test_Zpoly_get_coeff_si()
{
   int success = 1;
   Zpoly_t poly;
   mpz_t x;

   Zpoly_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_si(poly->coeffs[0], 47);
   mpz_set_si(poly->coeffs[1], -48);
   mpz_set_si(poly->coeffs[2], 49);
   
   success = success && (Zpoly_get_coeff_si(poly, 0) == 47);
   success = success && (Zpoly_get_coeff_si(poly, 1) == -48);
   success = success && (Zpoly_get_coeff_si(poly, 2) == 0);
   
   Zpoly_clear(poly);
   mpz_clear(x);
   return success;
}


int test_Zpoly_set_from_string()
{
   int success = 1;
   Zpoly_t poly;
   Zpoly_init(poly);

   Zpoly_set_from_string(poly, "");
   success = success && (poly->length == 0);

   Zpoly_set_from_string(poly, "   \t\n\r   ");
   success = success && (poly->length == 0);
   
   Zpoly_set_from_string(poly, "47");
   success = success && (poly->length == 1);
   success = success && (Zpoly_get_coeff_ui(poly, 0) == 47);
   
   Zpoly_set_from_string(poly, "   47   ");
   success = success && (poly->length == 1);
   success = success && (Zpoly_get_coeff_ui(poly, 0) == 47);
   
   Zpoly_set_from_string(poly, "   47 0   -23  ");
   success = success && (poly->length == 3);
   success = success && (Zpoly_get_coeff_ui(poly, 0) == 47);
   success = success && (Zpoly_get_coeff_ui(poly, 1) == 0);
   success = success && (Zpoly_get_coeff_si(poly, 2) == -23);
   
   // todo: also test a few cases where Zpoly_set_from_string()
   // should return 0
   
   Zpoly_clear(poly);
   return success;
}



int test_Zpoly_get_as_string()
{
   int success = 1;
   char buf[1000];
   
   Zpoly_t poly;
   Zpoly_init2(poly, 10);

   Zpoly_get_as_string(buf, poly);
   success = success && !strcmp(buf, "");
   
   poly->length = 2;
   mpz_set_si(poly->coeffs[1], -57);
   Zpoly_get_as_string(buf, poly);
   success = success && !strcmp(buf, "0 -57");
   success = success &&
             (Zpoly_get_string_size(poly) >= strlen("0 -57") + 1);
   
   Zpoly_clear(poly);
   return success;
}


int test__Zpoly_normalise()
{
   int success = 1;
   Zpoly_t poly;
   Zpoly_init2(poly, 10);
   
   poly->length = 3;
   _Zpoly_normalise(poly);
   success = success && (poly->length == 0);
   
   poly->length = 3;
   mpz_set_ui(poly->coeffs[1], 5);
   _Zpoly_normalise(poly);
   success = success && (poly->length == 2);

   Zpoly_clear(poly);
   return success;
}


int test__Zpoly_set()
{
   int success = 1;
   Zpoly_t poly1, poly2;
   Zpoly_init2(poly1, 10);
   Zpoly_init2(poly2, 10);

   Zpoly_set_from_string(poly1, "42 -5 0 3");
   _Zpoly_set(poly2, poly1);
   success = success && Zpoly_equal_str(poly2, "42 -5 0 3");

   Zpoly_clear(poly1);
   Zpoly_clear(poly2);
   return success;
}


int test__Zpoly_equal()
{
   int success = 1;
   Zpoly_t poly1, poly2;
   Zpoly_init2(poly1, 10);
   Zpoly_init2(poly2, 10);

   Zpoly_set_from_string(poly1, "42 -5 0 3");
   Zpoly_set_from_string(poly2, "42 -5 0 3");
   success = success && _Zpoly_equal(poly1, poly2);

   Zpoly_set_from_string(poly1, "42 -5 0 3");
   Zpoly_set_from_string(poly2, "42 -5 0 3 0");
   success = success && _Zpoly_equal(poly1, poly2);

   Zpoly_set_from_string(poly1, "42 -5 0 3");
   Zpoly_set_from_string(poly2, "42 -5 0 3 1");
   success = success && !_Zpoly_equal(poly1, poly2);

   Zpoly_set_from_string(poly1, "42 -5 0 3 4");
   Zpoly_set_from_string(poly2, "42 -5 0 3");
   success = success && !_Zpoly_equal(poly1, poly2);

   Zpoly_set_from_string(poly1, "42 -6 0 3");
   Zpoly_set_from_string(poly2, "42 -5 0 3");
   success = success && !_Zpoly_equal(poly1, poly2);

   Zpoly_set_from_string(poly1, "");
   Zpoly_set_from_string(poly2, "42 -5 0 3");
   success = success && !_Zpoly_equal(poly1, poly2);

   Zpoly_set_from_string(poly1, "");
   Zpoly_set_from_string(poly2, "0");
   success = success && _Zpoly_equal(poly1, poly2);

   Zpoly_clear(poly1);
   Zpoly_clear(poly2);
   return success;
}


int test__Zpoly_add()
{
   int success = 1;
   Zpoly_t poly[3];
   Zpoly_init2(poly[0], 10);
   Zpoly_init2(poly[1], 10);
   Zpoly_init2(poly[2], 10);

   // try various combinations of argument aliasing
   for (int i = 0; i < 3; i++)
   {
      Zpoly_t* target = &poly[i];

      Zpoly_set_from_string(poly[0], "");
      Zpoly_set_from_string(poly[1], "");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "");

      Zpoly_set_from_string(poly[0], "1");
      Zpoly_set_from_string(poly[1], "");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "1");

      Zpoly_set_from_string(poly[0], "");
      Zpoly_set_from_string(poly[1], "1");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "1");

      Zpoly_set_from_string(poly[0], "-1");
      Zpoly_set_from_string(poly[1], "1");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "");

      Zpoly_set_from_string(poly[0], "-1 0 47");
      Zpoly_set_from_string(poly[1], "0 0 0 0 8");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "-1 0 47 0 8");

      Zpoly_set_from_string(poly[0], "0 0 0 0 8");
      Zpoly_set_from_string(poly[1], "-1 0 47");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "-1 0 47 0 8");
   }

   Zpoly_set_from_string(poly[0], "2 3 4");
   Zpoly_set_from_string(poly[1], "123 456 789 123 456");
   _Zpoly_add(poly[1], poly[0], poly[0]);
   success = success && Zpoly_equal_str(poly[1], "4 6 8");

   Zpoly_set_from_string(poly[0], "2 3 4");
   _Zpoly_add(poly[0], poly[0], poly[0]);
   success = success && Zpoly_equal_str(poly[0], "4 6 8");

   Zpoly_clear(poly[0]);
   Zpoly_clear(poly[1]);
   Zpoly_clear(poly[2]);
   return success;
}


int test__Zpoly_sub()
{
   int success = 1;
   Zpoly_t poly[3];
   Zpoly_init2(poly[0], 10);
   Zpoly_init2(poly[1], 10);
   Zpoly_init2(poly[2], 10);

   // try various combinations of argument aliasing
   for (int i = 0; i < 3; i++)
   {
      Zpoly_t* target = &poly[i];

      Zpoly_set_from_string(poly[0], "");
      Zpoly_set_from_string(poly[1], "");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "");

      Zpoly_set_from_string(poly[0], "1");
      Zpoly_set_from_string(poly[1], "");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "1");

      Zpoly_set_from_string(poly[0], "");
      Zpoly_set_from_string(poly[1], "1");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "-1");

      Zpoly_set_from_string(poly[0], "-1");
      Zpoly_set_from_string(poly[1], "1");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "-2");

      Zpoly_set_from_string(poly[0], "-1 0 47");
      Zpoly_set_from_string(poly[1], "0 0 0 0 8");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "-1 0 47 0 -8");

      Zpoly_set_from_string(poly[0], "0 0 0 0 8");
      Zpoly_set_from_string(poly[1], "-1 0 47");
      Zpoly_set_from_string(poly[2], "123 456 789 123 456");
      _Zpoly_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal_str(*target, "1 0 -47 0 8");
   }

   Zpoly_set_from_string(poly[0], "2 3 4");
   Zpoly_set_from_string(poly[1], "123 456 789 123 456");
   _Zpoly_sub(poly[1], poly[0], poly[0]);
   success = success && Zpoly_equal_str(poly[1], "");

   Zpoly_set_from_string(poly[0], "2 3 4");
   _Zpoly_sub(poly[0], poly[0], poly[0]);
   success = success && Zpoly_equal_str(poly[0], "");

   Zpoly_clear(poly[0]);
   Zpoly_clear(poly[1]);
   Zpoly_clear(poly[2]);
   return success;
}


int test__Zpoly_negate()
{
   int success = 1;
   Zpoly_t poly1, poly2;
   Zpoly_init2(poly1, 10);
   Zpoly_init2(poly2, 10);

   // out-of-place

   Zpoly_set_from_string(poly1, "");
   Zpoly_set_from_string(poly2, "123 456 789 123 456");
   _Zpoly_negate(poly2, poly1);
   success = success && Zpoly_equal_str(poly2, "");

   Zpoly_set_from_string(poly1, "0 2 -5 6");
   Zpoly_set_from_string(poly2, "123 456 789 123 456");
   _Zpoly_negate(poly2, poly1);
   success = success && Zpoly_equal_str(poly2, "0 -2 5 -6");

   // in-place

   Zpoly_set_from_string(poly1, "");
   _Zpoly_negate(poly1, poly1);
   success = success && Zpoly_equal_str(poly1, "");

   Zpoly_set_from_string(poly1, "0 2 -5 6");
   _Zpoly_negate(poly1, poly1);
   success = success && Zpoly_equal_str(poly1, "0 -2 5 -6");

   Zpoly_clear(poly1);
   Zpoly_clear(poly2);
   return success;
}


#if 0
int test__Zpoly_mul()
{
   return 0;
}
#endif


int test__Zpoly_mul_naive()
{
   int success = 1;
   Zpoly_t poly1, poly2, poly3;
   Zpoly_init(poly1);
   Zpoly_init(poly2);
   Zpoly_init2(poly3, 20);

   Zpoly_set_from_string(poly1, "");
   Zpoly_set_from_string(poly2, "");
   _Zpoly_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal_str(poly3, "");
   
   Zpoly_set_from_string(poly1, "1 2 3");
   Zpoly_set_from_string(poly2, "");
   _Zpoly_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal_str(poly3, "");
   
   Zpoly_set_from_string(poly1, "1 2 3");
   Zpoly_set_from_string(poly2, "2");
   _Zpoly_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal_str(poly3, "2 4 6");
   
   Zpoly_set_from_string(poly1, "-3 4 0 2 56");
   Zpoly_set_from_string(poly2, "48 -2 3");
   _Zpoly_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal_str(poly3, "-144 198 -17 108 2684 -106 168");

   return success;
}


int test_Zpoly_set_coeff()
{
   return 0;
}


int test_Zpoly_set_coeff_ui()
{
   return 0;
}


int test_Zpoly_set_coeff_si()
{
   return 0;
}


int test_Zpoly_set()
{
   return 0;
}


int test_Zpoly_add()
{
   return 0;
}


int test_Zpoly_sub()
{
   return 0;
}


int test_Zpoly_negate()
{
   return 0;
}


int test_Zpoly_mul()
{
   return 0;
}


int test_Zpoly_mul_naive()
{
   return 0;
}


int test_Zpoly_mul_naive_KS()
{
   // todo: test inplace multiplication too

   int success = 1;
   
   unsigned long max_degree = 10;
   unsigned long max_bitsize = 10;
   Zpoly_t poly[4];
   for (unsigned long i = 0; i < 4; i++)
      Zpoly_init2(poly[i], max_degree*2 + 1);
   mpz_t temp;
   mpz_init(temp);

   unsigned long degree[2];
   unsigned long bitsize[2];

   for (degree[0] = 1; degree[0] <= max_degree; degree[0]++)
      for (degree[1] = 1; degree[1] <= max_degree; degree[1]++)
         for (bitsize[0] = 1; bitsize[0] <= max_bitsize; bitsize[0]++)
            for (bitsize[1] = 1; bitsize[1] <= max_bitsize; bitsize[1]++)
               for (unsigned long trial = 0; trial < 10; trial++)
               {
                  // generate random polys
                  for (unsigned long j = 0; j < 2; j++)
                  {
                     Zpoly_zero(poly[j]);
                     for (unsigned long i = 0; i < degree[j]; i++)
                     {
                        unsigned long bits = gmp_urandomm_ui(
                                       Zpoly_test_randstate, bitsize[j]+1);
                        mpz_rrandomb(temp, Zpoly_test_randstate, bits);
                        if (gmp_urandomb_ui(Zpoly_test_randstate, 1))
                           mpz_neg(temp, temp);
                        Zpoly_set_coeff(poly[j], i, temp);
                     }
                  }
                  
                  // compute product using naive multiplication and by
                  // naive KS, and compare answers
                  Zpoly_mul_naive(poly[2], poly[0], poly[1]);
                  Zpoly_mul_naive_KS(poly[3], poly[0], poly[1]);
                  success = success && Zpoly_equal(poly[2], poly[3]);
               }

   for (unsigned long i = 0; i < 4; i++)
      Zpoly_clear(poly[i]);
   mpz_clear(temp);

   return success;
}


#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void Zpoly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(Zpoly_get_coeff);
   RUN_TEST(Zpoly_get_coeff_ui);
   RUN_TEST(Zpoly_get_coeff_si);
   RUN_TEST(Zpoly_set_from_string);
   RUN_TEST(Zpoly_get_as_string);
   RUN_TEST(_Zpoly_normalise);
   RUN_TEST(_Zpoly_set);
   RUN_TEST(_Zpoly_equal);
   RUN_TEST(_Zpoly_add);
   RUN_TEST(_Zpoly_sub);
   RUN_TEST(_Zpoly_negate);
#if 0
   // disabled for the moment, since implementation is currently the same
   // as _Zpoly_mul_naive
   RUN_TEST(_Zpoly_mul);
#endif
   RUN_TEST(_Zpoly_mul_naive);
   RUN_TEST(Zpoly_set_coeff);
   RUN_TEST(Zpoly_set_coeff_ui);
   RUN_TEST(Zpoly_set_coeff_si);
   RUN_TEST(Zpoly_set);
   RUN_TEST(Zpoly_add);
   RUN_TEST(Zpoly_sub);
   RUN_TEST(Zpoly_negate);
   RUN_TEST(Zpoly_mul);
   RUN_TEST(Zpoly_mul_naive);
   RUN_TEST(Zpoly_mul_naive_KS);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   Zpoly_test_all();

   return 0;
}


// *************** end of file
