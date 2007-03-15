/****************************************************************************

Zpoly-test.c: Test code for Zpoly.c and Zpoly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "Zpoly.h"
#include <stdio.h>
#include <string.h>


// tests whether the given polynomial is equal to the one given by the string
// (only for testing purposes in this file)
int Zpoly_equal(Zpoly_mpz_t poly, char* s)
{
   Zpoly_mpz_t poly2;
   Zpoly_mpz_init(poly2);
   Zpoly_mpz_set_from_string(poly2, s);
   int result = Zpoly_mpz_equal(poly, poly2);
   Zpoly_mpz_clear(poly2);
   return result;
}



// all test functions return 1 on success, 0 on failure


int test_Zpoly_mpz_get_coeff()
{
   int success = 1;
   Zpoly_mpz_t poly;
   mpz_t x;

   Zpoly_mpz_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_ui(poly->coeffs[0], 47);
   mpz_set_ui(poly->coeffs[1], 48);
   mpz_set_ui(poly->coeffs[2], 49);
   
   Zpoly_mpz_get_coeff(x, poly, 0);
   success = success && !mpz_cmp_ui(x, 47);
   Zpoly_mpz_get_coeff(x, poly, 1);
   success = success && !mpz_cmp_ui(x, 48);
   Zpoly_mpz_get_coeff(x, poly, 2);
   success = success && !mpz_cmp_ui(x, 0);
   
   Zpoly_mpz_clear(poly);
   mpz_clear(x);
   return success;
}


int test_Zpoly_mpz_get_coeff_ui()
{
   int success = 1;
   Zpoly_mpz_t poly;
   mpz_t x;

   Zpoly_mpz_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_ui(poly->coeffs[0], 47);
   mpz_set_ui(poly->coeffs[1], 48);
   mpz_set_ui(poly->coeffs[2], 49);
   
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 0) == 47);
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 1) == 48);
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 2) == 0);
   
   Zpoly_mpz_clear(poly);
   mpz_clear(x);
   return success;
}


int test_Zpoly_mpz_get_coeff_si()
{
   int success = 1;
   Zpoly_mpz_t poly;
   mpz_t x;

   Zpoly_mpz_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_si(poly->coeffs[0], 47);
   mpz_set_si(poly->coeffs[1], -48);
   mpz_set_si(poly->coeffs[2], 49);
   
   success = success && (Zpoly_mpz_get_coeff_si(poly, 0) == 47);
   success = success && (Zpoly_mpz_get_coeff_si(poly, 1) == -48);
   success = success && (Zpoly_mpz_get_coeff_si(poly, 2) == 0);
   
   Zpoly_mpz_clear(poly);
   mpz_clear(x);
   return success;
}


int test_Zpoly_mpz_set_from_string()
{
   int success = 1;
   Zpoly_mpz_t poly;
   Zpoly_mpz_init(poly);

   Zpoly_mpz_set_from_string(poly, "");
   success = success && (poly->length == 0);

   Zpoly_mpz_set_from_string(poly, "   \t\n\r   ");
   success = success && (poly->length == 0);
   
   Zpoly_mpz_set_from_string(poly, "47");
   success = success && (poly->length == 1);
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 0) == 47);
   
   Zpoly_mpz_set_from_string(poly, "   47   ");
   success = success && (poly->length == 1);
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 0) == 47);
   
   Zpoly_mpz_set_from_string(poly, "   47 0   -23  ");
   success = success && (poly->length == 3);
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 0) == 47);
   success = success && (Zpoly_mpz_get_coeff_ui(poly, 1) == 0);
   success = success && (Zpoly_mpz_get_coeff_si(poly, 2) == -23);
   
   Zpoly_mpz_clear(poly);
   return success;
}



int test_Zpoly_mpz_get_as_string()
{
   int success = 1;
   char buf[1000];
   
   Zpoly_mpz_t poly;
   Zpoly_mpz_init2(poly, 10);

   Zpoly_mpz_get_as_string(buf, poly);
   success = success && !strcmp(buf, "");
   
   poly->length = 2;
   mpz_set_si(poly->coeffs[1], -57);
   Zpoly_mpz_get_as_string(buf, poly);
   success = success && !strcmp(buf, "0 -57");
   success = success &&
             (Zpoly_mpz_get_string_size(poly) >= strlen("0 -57") + 1);
   
   Zpoly_mpz_clear(poly);
   return success;
}


int test_Zpoly_mpz_raw_normalise()
{
   int success = 1;
   Zpoly_mpz_t poly;
   Zpoly_mpz_init2(poly, 10);
   
   poly->length = 3;
   Zpoly_mpz_raw_normalise(poly);
   success = success && (poly->length == 0);
   
   poly->length = 3;
   mpz_set_ui(poly->coeffs[1], 5);
   Zpoly_mpz_raw_normalise(poly);
   success = success && (poly->length == 2);

   Zpoly_mpz_clear(poly);
   return success;
}


int test_Zpoly_mpz_raw_set()
{
   int success = 1;
   Zpoly_mpz_t poly1, poly2;
   Zpoly_mpz_init2(poly1, 10);
   Zpoly_mpz_init2(poly2, 10);

   Zpoly_mpz_set_from_string(poly1, "42 -5 0 3");
   Zpoly_mpz_raw_set(poly2, poly1);
   success = success && Zpoly_equal(poly2, "42 -5 0 3");

   Zpoly_mpz_clear(poly1);
   Zpoly_mpz_clear(poly2);
   return success;
}


int test_Zpoly_mpz_raw_equal()
{
   int success = 1;
   Zpoly_mpz_t poly1, poly2;
   Zpoly_mpz_init2(poly1, 10);
   Zpoly_mpz_init2(poly2, 10);

   Zpoly_mpz_set_from_string(poly1, "42 -5 0 3");
   Zpoly_mpz_set_from_string(poly2, "42 -5 0 3");
   success = success && Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_set_from_string(poly1, "42 -5 0 3");
   Zpoly_mpz_set_from_string(poly2, "42 -5 0 3 0");
   success = success && Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_set_from_string(poly1, "42 -5 0 3");
   Zpoly_mpz_set_from_string(poly2, "42 -5 0 3 1");
   success = success && !Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_set_from_string(poly1, "42 -5 0 3 4");
   Zpoly_mpz_set_from_string(poly2, "42 -5 0 3");
   success = success && !Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_set_from_string(poly1, "42 -6 0 3");
   Zpoly_mpz_set_from_string(poly2, "42 -5 0 3");
   success = success && !Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_set_from_string(poly1, "");
   Zpoly_mpz_set_from_string(poly2, "42 -5 0 3");
   success = success && !Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_set_from_string(poly1, "");
   Zpoly_mpz_set_from_string(poly2, "0");
   success = success && Zpoly_mpz_raw_equal(poly1, poly2);

   Zpoly_mpz_clear(poly1);
   Zpoly_mpz_clear(poly2);
   return success;
}


int test_Zpoly_mpz_raw_add()
{
   int success = 1;
   Zpoly_mpz_t poly[3];
   Zpoly_mpz_init2(poly[0], 10);
   Zpoly_mpz_init2(poly[1], 10);
   Zpoly_mpz_init2(poly[2], 10);

   // try various combinations of argument aliasing
   for (int i = 0; i < 3; i++)
   {
      Zpoly_mpz_t* target = &poly[i];

      Zpoly_mpz_set_from_string(poly[0], "");
      Zpoly_mpz_set_from_string(poly[1], "");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "");

      Zpoly_mpz_set_from_string(poly[0], "1");
      Zpoly_mpz_set_from_string(poly[1], "");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "1");

      Zpoly_mpz_set_from_string(poly[0], "");
      Zpoly_mpz_set_from_string(poly[1], "1");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "1");

      Zpoly_mpz_set_from_string(poly[0], "-1");
      Zpoly_mpz_set_from_string(poly[1], "1");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "");

      Zpoly_mpz_set_from_string(poly[0], "-1 0 47");
      Zpoly_mpz_set_from_string(poly[1], "0 0 0 0 8");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "-1 0 47 0 8");

      Zpoly_mpz_set_from_string(poly[0], "0 0 0 0 8");
      Zpoly_mpz_set_from_string(poly[1], "-1 0 47");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_add(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "-1 0 47 0 8");
   }

   Zpoly_mpz_set_from_string(poly[0], "2 3 4");
   Zpoly_mpz_set_from_string(poly[1], "123 456 789 123 456");
   Zpoly_mpz_raw_add(poly[1], poly[0], poly[0]);
   success = success && Zpoly_equal(poly[1], "4 6 8");

   Zpoly_mpz_set_from_string(poly[0], "2 3 4");
   Zpoly_mpz_raw_add(poly[0], poly[0], poly[0]);
   success = success && Zpoly_equal(poly[0], "4 6 8");

   Zpoly_mpz_clear(poly[0]);
   Zpoly_mpz_clear(poly[1]);
   Zpoly_mpz_clear(poly[2]);
   return success;
}


int test_Zpoly_mpz_raw_sub()
{
   int success = 1;
   Zpoly_mpz_t poly[3];
   Zpoly_mpz_init2(poly[0], 10);
   Zpoly_mpz_init2(poly[1], 10);
   Zpoly_mpz_init2(poly[2], 10);

   // try various combinations of argument aliasing
   for (int i = 0; i < 3; i++)
   {
      Zpoly_mpz_t* target = &poly[i];

      Zpoly_mpz_set_from_string(poly[0], "");
      Zpoly_mpz_set_from_string(poly[1], "");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "");

      Zpoly_mpz_set_from_string(poly[0], "1");
      Zpoly_mpz_set_from_string(poly[1], "");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "1");

      Zpoly_mpz_set_from_string(poly[0], "");
      Zpoly_mpz_set_from_string(poly[1], "1");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "-1");

      Zpoly_mpz_set_from_string(poly[0], "-1");
      Zpoly_mpz_set_from_string(poly[1], "1");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "-2");

      Zpoly_mpz_set_from_string(poly[0], "-1 0 47");
      Zpoly_mpz_set_from_string(poly[1], "0 0 0 0 8");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "-1 0 47 0 -8");

      Zpoly_mpz_set_from_string(poly[0], "0 0 0 0 8");
      Zpoly_mpz_set_from_string(poly[1], "-1 0 47");
      Zpoly_mpz_set_from_string(poly[2], "123 456 789 123 456");
      Zpoly_mpz_raw_sub(*target, poly[0], poly[1]);
      success = success && Zpoly_equal(*target, "1 0 -47 0 8");
   }

   Zpoly_mpz_set_from_string(poly[0], "2 3 4");
   Zpoly_mpz_set_from_string(poly[1], "123 456 789 123 456");
   Zpoly_mpz_raw_sub(poly[1], poly[0], poly[0]);
   success = success && Zpoly_equal(poly[1], "");

   Zpoly_mpz_set_from_string(poly[0], "2 3 4");
   Zpoly_mpz_raw_sub(poly[0], poly[0], poly[0]);
   success = success && Zpoly_equal(poly[0], "");

   Zpoly_mpz_clear(poly[0]);
   Zpoly_mpz_clear(poly[1]);
   Zpoly_mpz_clear(poly[2]);
   return success;
}


int test_Zpoly_mpz_raw_negate()
{
   int success = 1;
   Zpoly_mpz_t poly1, poly2;
   Zpoly_mpz_init2(poly1, 10);
   Zpoly_mpz_init2(poly2, 10);

   // out-of-place

   Zpoly_mpz_set_from_string(poly1, "");
   Zpoly_mpz_set_from_string(poly2, "123 456 789 123 456");
   Zpoly_mpz_raw_negate(poly2, poly1);
   success = success && Zpoly_equal(poly2, "");

   Zpoly_mpz_set_from_string(poly1, "0 2 -5 6");
   Zpoly_mpz_set_from_string(poly2, "123 456 789 123 456");
   Zpoly_mpz_raw_negate(poly2, poly1);
   success = success && Zpoly_equal(poly2, "0 -2 5 -6");

   // in-place

   Zpoly_mpz_set_from_string(poly1, "");
   Zpoly_mpz_raw_negate(poly1, poly1);
   success = success && Zpoly_equal(poly1, "");

   Zpoly_mpz_set_from_string(poly1, "0 2 -5 6");
   Zpoly_mpz_raw_negate(poly1, poly1);
   success = success && Zpoly_equal(poly1, "0 -2 5 -6");

   Zpoly_mpz_clear(poly1);
   Zpoly_mpz_clear(poly2);
   return success;
}


#if 0
int test_Zpoly_mpz_raw_mul()
{
   return 0;
}
#endif


int test_Zpoly_mpz_raw_mul_naive()
{
   int success = 1;
   Zpoly_mpz_t poly1, poly2, poly3;
   Zpoly_mpz_init(poly1);
   Zpoly_mpz_init(poly2);
   Zpoly_mpz_init2(poly3, 20);

   Zpoly_mpz_set_from_string(poly1, "");
   Zpoly_mpz_set_from_string(poly2, "");
   Zpoly_mpz_raw_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal(poly3, "");
   
   Zpoly_mpz_set_from_string(poly1, "1 2 3");
   Zpoly_mpz_set_from_string(poly2, "");
   Zpoly_mpz_raw_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal(poly3, "");
   
   Zpoly_mpz_set_from_string(poly1, "1 2 3");
   Zpoly_mpz_set_from_string(poly2, "2");
   Zpoly_mpz_raw_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal(poly3, "2 4 6");
   
   Zpoly_mpz_set_from_string(poly1, "-3 4 0 2 56");
   Zpoly_mpz_set_from_string(poly2, "48 -2 3");
   Zpoly_mpz_raw_mul_naive(poly3, poly1, poly2);
   success = success && Zpoly_equal(poly3, "-144 198 -17 108 2684 -106 168");

   return success;
}


int test_Zpoly_mpz_set_coeff()
{
   return 0;
}


int test_Zpoly_mpz_set_coeff_ui()
{
   return 0;
}


int test_Zpoly_mpz_set_coeff_si()
{
   return 0;
}


int test_Zpoly_mpz_set()
{
   return 0;
}


int test_Zpoly_mpz_add()
{
   return 0;
}


int test_Zpoly_mpz_sub()
{
   return 0;
}


int test_Zpoly_mpz_negate()
{
   return 0;
}


int test_Zpoly_mpz_mul()
{
   return 0;
}


int test_Zpoly_mpz_mul_naive()
{
   return 0;
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

   RUN_TEST(Zpoly_mpz_get_coeff);
   RUN_TEST(Zpoly_mpz_get_coeff_ui);
   RUN_TEST(Zpoly_mpz_get_coeff_si);
   RUN_TEST(Zpoly_mpz_set_from_string);
   RUN_TEST(Zpoly_mpz_get_as_string);
   RUN_TEST(Zpoly_mpz_raw_normalise);
   RUN_TEST(Zpoly_mpz_raw_set);
   RUN_TEST(Zpoly_mpz_raw_equal);
   RUN_TEST(Zpoly_mpz_raw_add);
   RUN_TEST(Zpoly_mpz_raw_sub);
   RUN_TEST(Zpoly_mpz_raw_negate);
#if 0
   // disabled for the moment, since implementation is currently the same
   // as Zpoly_mpz_raw_mul_naive
   RUN_TEST(Zpoly_mpz_raw_mul);
#endif
   RUN_TEST(Zpoly_mpz_raw_mul_naive);
   RUN_TEST(Zpoly_mpz_set_coeff);
   RUN_TEST(Zpoly_mpz_set_coeff_ui);
   RUN_TEST(Zpoly_mpz_set_coeff_si);
   RUN_TEST(Zpoly_mpz_set);
   RUN_TEST(Zpoly_mpz_add);
   RUN_TEST(Zpoly_mpz_sub);
   RUN_TEST(Zpoly_mpz_negate);
   RUN_TEST(Zpoly_mpz_mul);
   RUN_TEST(Zpoly_mpz_mul_naive);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   Zpoly_test_all();

   return 0;
}


// *************** end of file
