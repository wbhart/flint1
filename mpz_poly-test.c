/****************************************************************************

mpz_poly-test.c: Test code for mpz_poly.c and mpz_poly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "flint.h"
#include "mpz_poly.h"


#if 0    // old test code, disabled for the moment


gmp_randstate_t mpz_poly_test_randstate;


// tests whether the given polynomial is equal to the one given by the string
// (only for testing purposes in this file)
int mpz_poly_equal_str(mpz_poly_t poly, char* s)
{
   mpz_poly_t poly2;
   mpz_poly_init(poly2);
   mpz_poly_set_from_string(poly2, s);
   int result = mpz_poly_equal(poly, poly2);
   mpz_poly_clear(poly2);
   return result;
}



// all test functions return 1 on success, 0 on failure


int test_mpz_poly_get_coeff()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_t x;

   mpz_poly_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_ui(poly->coeffs[0], 47);
   mpz_set_ui(poly->coeffs[1], 48);
   mpz_set_ui(poly->coeffs[2], 49);
   
   mpz_poly_get_coeff(x, poly, 0);
   success = success && !mpz_cmp_ui(x, 47);
   mpz_poly_get_coeff(x, poly, 1);
   success = success && !mpz_cmp_ui(x, 48);
   mpz_poly_get_coeff(x, poly, 2);
   success = success && !mpz_cmp_ui(x, 0);
   
   mpz_poly_clear(poly);
   mpz_clear(x);
   return success;
}


int test_mpz_poly_get_coeff_ui()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_t x;

   mpz_poly_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_ui(poly->coeffs[0], 47);
   mpz_set_ui(poly->coeffs[1], 48);
   mpz_set_ui(poly->coeffs[2], 49);
   
   success = success && (mpz_poly_get_coeff_ui(poly, 0) == 47);
   success = success && (mpz_poly_get_coeff_ui(poly, 1) == 48);
   success = success && (mpz_poly_get_coeff_ui(poly, 2) == 0);
   
   mpz_poly_clear(poly);
   mpz_clear(x);
   return success;
}


int test_mpz_poly_get_coeff_si()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_t x;

   mpz_poly_init2(poly, 3);
   mpz_init(x);
   
   poly->length = 2;
   mpz_set_si(poly->coeffs[0], 47);
   mpz_set_si(poly->coeffs[1], -48);
   mpz_set_si(poly->coeffs[2], 49);
   
   success = success && (mpz_poly_get_coeff_si(poly, 0) == 47);
   success = success && (mpz_poly_get_coeff_si(poly, 1) == -48);
   success = success && (mpz_poly_get_coeff_si(poly, 2) == 0);
   
   mpz_poly_clear(poly);
   mpz_clear(x);
   return success;
}


int test_mpz_poly_set_from_string()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_poly_init(poly);

   mpz_poly_set_from_string(poly, "");
   success = success && (poly->length == 0);

   mpz_poly_set_from_string(poly, "   \t\n\r   ");
   success = success && (poly->length == 0);
   
   mpz_poly_set_from_string(poly, "47");
   success = success && (poly->length == 1);
   success = success && (mpz_poly_get_coeff_ui(poly, 0) == 47);
   
   mpz_poly_set_from_string(poly, "   47   ");
   success = success && (poly->length == 1);
   success = success && (mpz_poly_get_coeff_ui(poly, 0) == 47);
   
   mpz_poly_set_from_string(poly, "   47 0   -23  ");
   success = success && (poly->length == 3);
   success = success && (mpz_poly_get_coeff_ui(poly, 0) == 47);
   success = success && (mpz_poly_get_coeff_ui(poly, 1) == 0);
   success = success && (mpz_poly_get_coeff_si(poly, 2) == -23);
   
   // todo: also test a few cases where mpz_poly_set_from_string()
   // should return 0
   
   mpz_poly_clear(poly);
   return success;
}



int test_mpz_poly_get_as_string()
{
   int success = 1;
   char buf[1000];
   
   mpz_poly_t poly;
   mpz_poly_init2(poly, 10);

   mpz_poly_get_as_string(buf, poly);
   success = success && !strcmp(buf, "");
   
   poly->length = 2;
   mpz_set_si(poly->coeffs[1], -57);
   mpz_poly_get_as_string(buf, poly);
   success = success && !strcmp(buf, "0 -57");
   success = success &&
             (mpz_poly_get_string_size(poly) >= strlen("0 -57") + 1);
   
   mpz_poly_clear(poly);
   return success;
}


int test__mpz_poly_normalise()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_poly_init2(poly, 10);
   
   poly->length = 3;
   _mpz_poly_normalise(poly);
   success = success && (poly->length == 0);
   
   poly->length = 3;
   mpz_set_ui(poly->coeffs[1], 5);
   _mpz_poly_normalise(poly);
   success = success && (poly->length == 2);

   mpz_poly_clear(poly);
   return success;
}


int test__mpz_poly_set()
{
   int success = 1;
   mpz_poly_t poly1, poly2;
   mpz_poly_init2(poly1, 10);
   mpz_poly_init2(poly2, 10);

   mpz_poly_set_from_string(poly1, "42 -5 0 3");
   _mpz_poly_set(poly2, poly1);
   success = success && mpz_poly_equal_str(poly2, "42 -5 0 3");

   mpz_poly_clear(poly1);
   mpz_poly_clear(poly2);
   return success;
}


int test__mpz_poly_equal()
{
   int success = 1;
   mpz_poly_t poly1, poly2;
   mpz_poly_init2(poly1, 10);
   mpz_poly_init2(poly2, 10);

   mpz_poly_set_from_string(poly1, "42 -5 0 3");
   mpz_poly_set_from_string(poly2, "42 -5 0 3");
   success = success && _mpz_poly_equal(poly1, poly2);

   mpz_poly_set_from_string(poly1, "42 -5 0 3");
   mpz_poly_set_from_string(poly2, "42 -5 0 3 0");
   success = success && _mpz_poly_equal(poly1, poly2);

   mpz_poly_set_from_string(poly1, "42 -5 0 3");
   mpz_poly_set_from_string(poly2, "42 -5 0 3 1");
   success = success && !_mpz_poly_equal(poly1, poly2);

   mpz_poly_set_from_string(poly1, "42 -5 0 3 4");
   mpz_poly_set_from_string(poly2, "42 -5 0 3");
   success = success && !_mpz_poly_equal(poly1, poly2);

   mpz_poly_set_from_string(poly1, "42 -6 0 3");
   mpz_poly_set_from_string(poly2, "42 -5 0 3");
   success = success && !_mpz_poly_equal(poly1, poly2);

   mpz_poly_set_from_string(poly1, "");
   mpz_poly_set_from_string(poly2, "42 -5 0 3");
   success = success && !_mpz_poly_equal(poly1, poly2);

   mpz_poly_set_from_string(poly1, "");
   mpz_poly_set_from_string(poly2, "0");
   success = success && _mpz_poly_equal(poly1, poly2);

   mpz_poly_clear(poly1);
   mpz_poly_clear(poly2);
   return success;
}


int test__mpz_poly_add()
{
   int success = 1;
   mpz_poly_t poly[3];
   mpz_poly_init2(poly[0], 10);
   mpz_poly_init2(poly[1], 10);
   mpz_poly_init2(poly[2], 10);

   // try various combinations of argument aliasing
   for (int i = 0; i < 3; i++)
   {
      mpz_poly_t* target = &poly[i];

      mpz_poly_set_from_string(poly[0], "");
      mpz_poly_set_from_string(poly[1], "");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_add(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "");

      mpz_poly_set_from_string(poly[0], "1");
      mpz_poly_set_from_string(poly[1], "");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_add(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "1");

      mpz_poly_set_from_string(poly[0], "");
      mpz_poly_set_from_string(poly[1], "1");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_add(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "1");

      mpz_poly_set_from_string(poly[0], "-1");
      mpz_poly_set_from_string(poly[1], "1");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_add(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "");

      mpz_poly_set_from_string(poly[0], "-1 0 47");
      mpz_poly_set_from_string(poly[1], "0 0 0 0 8");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_add(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "-1 0 47 0 8");

      mpz_poly_set_from_string(poly[0], "0 0 0 0 8");
      mpz_poly_set_from_string(poly[1], "-1 0 47");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_add(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "-1 0 47 0 8");
   }

   mpz_poly_set_from_string(poly[0], "2 3 4");
   mpz_poly_set_from_string(poly[1], "123 456 789 123 456");
   _mpz_poly_add(poly[1], poly[0], poly[0]);
   success = success && mpz_poly_equal_str(poly[1], "4 6 8");

   mpz_poly_set_from_string(poly[0], "2 3 4");
   _mpz_poly_add(poly[0], poly[0], poly[0]);
   success = success && mpz_poly_equal_str(poly[0], "4 6 8");

   mpz_poly_clear(poly[0]);
   mpz_poly_clear(poly[1]);
   mpz_poly_clear(poly[2]);
   return success;
}


int test__mpz_poly_sub()
{
   int success = 1;
   mpz_poly_t poly[3];
   mpz_poly_init2(poly[0], 10);
   mpz_poly_init2(poly[1], 10);
   mpz_poly_init2(poly[2], 10);

   // try various combinations of argument aliasing
   for (int i = 0; i < 3; i++)
   {
      mpz_poly_t* target = &poly[i];

      mpz_poly_set_from_string(poly[0], "");
      mpz_poly_set_from_string(poly[1], "");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_sub(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "");

      mpz_poly_set_from_string(poly[0], "1");
      mpz_poly_set_from_string(poly[1], "");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_sub(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "1");

      mpz_poly_set_from_string(poly[0], "");
      mpz_poly_set_from_string(poly[1], "1");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_sub(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "-1");

      mpz_poly_set_from_string(poly[0], "-1");
      mpz_poly_set_from_string(poly[1], "1");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_sub(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "-2");

      mpz_poly_set_from_string(poly[0], "-1 0 47");
      mpz_poly_set_from_string(poly[1], "0 0 0 0 8");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_sub(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "-1 0 47 0 -8");

      mpz_poly_set_from_string(poly[0], "0 0 0 0 8");
      mpz_poly_set_from_string(poly[1], "-1 0 47");
      mpz_poly_set_from_string(poly[2], "123 456 789 123 456");
      _mpz_poly_sub(*target, poly[0], poly[1]);
      success = success && mpz_poly_equal_str(*target, "1 0 -47 0 8");
   }

   mpz_poly_set_from_string(poly[0], "2 3 4");
   mpz_poly_set_from_string(poly[1], "123 456 789 123 456");
   _mpz_poly_sub(poly[1], poly[0], poly[0]);
   success = success && mpz_poly_equal_str(poly[1], "");

   mpz_poly_set_from_string(poly[0], "2 3 4");
   _mpz_poly_sub(poly[0], poly[0], poly[0]);
   success = success && mpz_poly_equal_str(poly[0], "");

   mpz_poly_clear(poly[0]);
   mpz_poly_clear(poly[1]);
   mpz_poly_clear(poly[2]);
   return success;
}


int test__mpz_poly_negate()
{
   int success = 1;
   mpz_poly_t poly1, poly2;
   mpz_poly_init2(poly1, 10);
   mpz_poly_init2(poly2, 10);

   // out-of-place

   mpz_poly_set_from_string(poly1, "");
   mpz_poly_set_from_string(poly2, "123 456 789 123 456");
   _mpz_poly_negate(poly2, poly1);
   success = success && mpz_poly_equal_str(poly2, "");

   mpz_poly_set_from_string(poly1, "0 2 -5 6");
   mpz_poly_set_from_string(poly2, "123 456 789 123 456");
   _mpz_poly_negate(poly2, poly1);
   success = success && mpz_poly_equal_str(poly2, "0 -2 5 -6");

   // in-place

   mpz_poly_set_from_string(poly1, "");
   _mpz_poly_negate(poly1, poly1);
   success = success && mpz_poly_equal_str(poly1, "");

   mpz_poly_set_from_string(poly1, "0 2 -5 6");
   _mpz_poly_negate(poly1, poly1);
   success = success && mpz_poly_equal_str(poly1, "0 -2 5 -6");

   mpz_poly_clear(poly1);
   mpz_poly_clear(poly2);
   return success;
}


#if 0
int test__mpz_poly_mul()
{
   return 0;
}
#endif


int test__mpz_poly_mul_naive()
{
   int success = 1;
   mpz_poly_t poly1, poly2, poly3;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init2(poly3, 20);

   mpz_poly_set_from_string(poly1, "");
   mpz_poly_set_from_string(poly2, "");
   _mpz_poly_mul_naive(poly3, poly1, poly2);
   success = success && mpz_poly_equal_str(poly3, "");
   
   mpz_poly_set_from_string(poly1, "1 2 3");
   mpz_poly_set_from_string(poly2, "");
   _mpz_poly_mul_naive(poly3, poly1, poly2);
   success = success && mpz_poly_equal_str(poly3, "");
   
   mpz_poly_set_from_string(poly1, "1 2 3");
   mpz_poly_set_from_string(poly2, "2");
   _mpz_poly_mul_naive(poly3, poly1, poly2);
   success = success && mpz_poly_equal_str(poly3, "2 4 6");
   
   mpz_poly_set_from_string(poly1, "-3 4 0 2 56");
   mpz_poly_set_from_string(poly2, "48 -2 3");
   _mpz_poly_mul_naive(poly3, poly1, poly2);
   success = success && mpz_poly_equal_str(poly3, "-144 198 -17 108 2684 -106 168");

   mpz_poly_clear(poly1);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly3);
   return success;
}


int test_mpz_poly_monic_inverse()
{
   int success = 1;
   mpz_poly_t poly1, poly2;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);

   for (unsigned long deg1 = 2; deg1 <= 10; deg1++)
   {
      for (unsigned long trial = 0; trial < 20; trial++)
      {
         // generate random input poly
         mpz_poly_set_coeff_ui(poly1, deg1, 1);
         for (unsigned long i = 0; i < deg1; i++)
            mpz_poly_set_coeff_si(poly1, i, gmp_urandomb_ui(mpz_poly_test_randstate, 10) - 512);
      
         // try computing inverses to various lengths
         for (unsigned long deg2 = deg1; deg2 <= 50; deg2++)
         {
            mpz_poly_t poly3;
            mpz_poly_init(poly3);
            
            mpz_poly_monic_inverse(poly3, poly1, deg2);
            if (poly3->length-1 != deg2)
               success = 0;
            else
            {
               // check correctness by multiplying back together
               mpz_poly_mul(poly2, poly1, poly3);
               success = success && !mpz_cmp_ui(poly2->coeffs[deg1+deg2], 1);
               for (unsigned long i = 0; i < deg2; i++)
                  success = success && !mpz_sgn(poly2->coeffs[deg1+i]);
            }
            
            mpz_poly_clear(poly3);
         }
      }
   }
   
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
   return success;
}


int test_mpz_poly_set_coeff()
{
   return 0;
}


int test_mpz_poly_set_coeff_ui()
{
   return 0;
}


int test_mpz_poly_set_coeff_si()
{
   return 0;
}


int test_mpz_poly_set()
{
   return 0;
}


int test_mpz_poly_add()
{
   return 0;
}


int test_mpz_poly_sub()
{
   return 0;
}


int test_mpz_poly_negate()
{
   return 0;
}


int test_mpz_poly_mul()
{
   return 0;
}


int test_mpz_poly_mul_naive()
{
   return 0;
}


int test_mpz_poly_mul_naive_KS()
{
   // todo: test inplace multiplication too

   int success = 1;
   
   unsigned long max_degree = 10;
   unsigned long max_bitsize = 10;
   mpz_poly_t poly[4];
   for (unsigned long i = 0; i < 4; i++)
      mpz_poly_init2(poly[i], max_degree*2 + 1);
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
                     mpz_poly_zero(poly[j]);
                     for (unsigned long i = 0; i < degree[j]; i++)
                     {
                        unsigned long bits = gmp_urandomm_ui(
                                       mpz_poly_test_randstate, bitsize[j]+1);
                        mpz_rrandomb(temp, mpz_poly_test_randstate, bits);
                        if (gmp_urandomb_ui(mpz_poly_test_randstate, 1))
                           mpz_neg(temp, temp);
                        mpz_poly_set_coeff(poly[j], i, temp);
                     }
                  }
                  
                  // compute product using naive multiplication and by
                  // naive KS, and compare answers
                  mpz_poly_mul_naive(poly[2], poly[0], poly[1]);
                  mpz_poly_mul_naive_KS(poly[3], poly[0], poly[1]);
                  success = success && mpz_poly_equal(poly[2], poly[3]);
               }

   for (unsigned long i = 0; i < 4; i++)
      mpz_poly_clear(poly[i]);
   mpz_clear(temp);

   return success;
}


int test_mpz_poly_sqr_naive_KS()
{
   // todo: test inplace multiplication too

   int success = 1;
   
   unsigned long max_degree = 10;
   unsigned long max_bitsize = 10;
   mpz_poly_t poly[3];
   for (unsigned long i = 0; i < 3; i++)
      mpz_poly_init2(poly[i], max_degree*2 + 1);
   mpz_t temp;
   mpz_init(temp);

   unsigned long degree;
   unsigned long bitsize;

   for (degree = 1; degree <= max_degree; degree++)
      for (bitsize = 1; bitsize <= max_bitsize; bitsize++)
         for (unsigned long trial = 0; trial < 10; trial++)
         {
            // generate random polys
            mpz_poly_zero(poly[0]);
            for (unsigned long i = 0; i < degree; i++)
            {
               unsigned long bits = gmp_urandomm_ui(
                                          mpz_poly_test_randstate, bitsize+1);
               mpz_rrandomb(temp, mpz_poly_test_randstate, bits);
               if (gmp_urandomb_ui(mpz_poly_test_randstate, 1))
                  mpz_neg(temp, temp);
               mpz_poly_set_coeff(poly[0], i, temp);
            }
            
            // compute product using naive multiplication and by
            // naive KS, and compare answers
            mpz_poly_sqr_naive(poly[1], poly[0]);
            mpz_poly_sqr_naive_KS(poly[2], poly[0]);
            success = success && mpz_poly_equal(poly[1], poly[2]);
         }

   for (unsigned long i = 0; i < 3; i++)
      mpz_poly_clear(poly[i]);
   mpz_clear(temp);

   return success;
}


#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void mpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(mpz_poly_get_coeff);
   RUN_TEST(mpz_poly_get_coeff_ui);
   RUN_TEST(mpz_poly_get_coeff_si);
   RUN_TEST(mpz_poly_set_from_string);
   RUN_TEST(mpz_poly_get_as_string);
   RUN_TEST(_mpz_poly_normalise);
   RUN_TEST(_mpz_poly_set);
   RUN_TEST(_mpz_poly_equal);
   RUN_TEST(_mpz_poly_add);
   RUN_TEST(_mpz_poly_sub);
   RUN_TEST(_mpz_poly_negate);
#if 0
   // disabled for the moment, since implementation is currently the same
   // as _mpz_poly_mul_naive
   RUN_TEST(_mpz_poly_mul);
#endif
   RUN_TEST(_mpz_poly_mul_naive);
   RUN_TEST(mpz_poly_set_coeff);
   RUN_TEST(mpz_poly_set_coeff_ui);
   RUN_TEST(mpz_poly_set_coeff_si);
   RUN_TEST(mpz_poly_set);
   RUN_TEST(mpz_poly_add);
   RUN_TEST(mpz_poly_sub);
   RUN_TEST(mpz_poly_negate);
   RUN_TEST(mpz_poly_mul);
   RUN_TEST(mpz_poly_mul_naive);
   RUN_TEST(mpz_poly_mul_naive_KS);
#if 0
   // disabled for the moment, because we haven't written the naive
   // squaring routine yet
   RUN_TEST(mpz_poly_sqr_naive_KS);
#endif
   RUN_TEST(mpz_poly_monic_inverse);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(mpz_poly_test_randstate);
   mpz_poly_test_all();

   return 0;
}

#endif


int main()
{
   mpz_poly_t x, y;
   mpz_poly_init(x);
   mpz_poly_init(y);

   mpz_poly_set_coeff_ui(x, 0, 2);
   mpz_poly_set_coeff_si(x, 2, 1);
   mpz_poly_monic_inverse(y, x, 20);
   mpz_poly_print(x); printf("\n");
   mpz_poly_print(y); printf("\n");

   mpz_poly_clear(x);
   mpz_poly_clear(y);
}


// *************** end of file
