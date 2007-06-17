/****************************************************************************

mpz_poly-test.c: Test code for mpz_poly.c and mpz_poly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "flint.h"
#include "mpz_poly.h"
#include "test-support.h"



// tests whether the given polynomial is equal to the one given by the string
// (only for testing purposes in this file)
int mpz_poly_equal_str(mpz_poly_t poly, char* s)
{
   mpz_poly_t poly2;
   mpz_poly_init(poly2);
   mpz_poly_from_string(poly2, s);
   int result = mpz_poly_equal(poly, poly2);
   mpz_poly_clear(poly2);
   return result;
}


/****************************************************************************

   Setting/retrieving coefficients

****************************************************************************/


int test_mpz_poly_get_coeff_ptr()
{
   int success = 1;
   mpz_poly_t poly;

   mpz_poly_init(poly);
   mpz_poly_init_upto(poly, 3);
   poly->length = 2;
   
   success = success && (mpz_poly_get_coeff_ptr(poly, 0) == poly->coeffs);
   success = success && (mpz_poly_get_coeff_ptr(poly, 1) == poly->coeffs + 1);
   success = success && (mpz_poly_get_coeff_ptr(poly, 2) == NULL);
   
   mpz_poly_clear(poly);
   return success;
}


int test_mpz_poly_get_coeff()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_t x;

   mpz_poly_init(poly);
   mpz_poly_init_upto(poly, 3);
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

   mpz_poly_init(poly);
   mpz_poly_init_upto(poly, 3);
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

   mpz_poly_init(poly);
   mpz_poly_init_upto(poly, 3);
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


int test_mpz_poly_set_coeff()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_t x, y, zero;

   mpz_poly_init(poly);
   mpz_init(x);
   mpz_init(y);
   mpz_init(zero);
   mpz_set_ui(x, 42);
   mpz_set_ui(y, 37);

   mpz_poly_set_coeff(poly, 2, x);
   success = success && mpz_poly_equal_str(poly, "3  0 0 42");

   mpz_poly_set_coeff(poly, 5, y);
   success = success && mpz_poly_equal_str(poly, "6  0 0 42 0 0 37");

   mpz_poly_set_coeff(poly, 1, y);
   success = success && mpz_poly_equal_str(poly, "6  0 37 42 0 0 37");

   mpz_poly_set_coeff(poly, 5, x);
   success = success && mpz_poly_equal_str(poly, "6  0 37 42 0 0 42");

   mpz_poly_set_coeff(poly, 2, zero);
   success = success && mpz_poly_equal_str(poly, "6  0 37 0 0 0 42");

   mpz_poly_set_coeff(poly, 8, zero);
   success = success && mpz_poly_equal_str(poly, "6  0 37 0 0 0 42");

   mpz_poly_set_coeff(poly, 5, zero);
   success = success && mpz_poly_equal_str(poly, "2  0 37");
   
   mpz_poly_truncate(poly, poly, 1);
   success = success && mpz_poly_equal_str(poly, "0");

   mpz_poly_set_coeff(poly, 3, x);
   success = success && mpz_poly_equal_str(poly, "4  0 0 0 42");

   mpz_poly_set_coeff(poly, 3, zero);
   success = success && mpz_poly_equal_str(poly, "0");

   mpz_poly_clear(poly);
   mpz_clear(y);
   mpz_clear(x);
   mpz_clear(zero);

   return success;
}


int test_mpz_poly_set_coeff_ui()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_poly_init(poly);

   mpz_poly_set_coeff_ui(poly, 2, 42);
   success = success && mpz_poly_equal_str(poly, "3  0 0 42");

   mpz_poly_set_coeff_ui(poly, 5, 37);
   success = success && mpz_poly_equal_str(poly, "6  0 0 42 0 0 37");

   mpz_poly_set_coeff_ui(poly, 1, 37);
   success = success && mpz_poly_equal_str(poly, "6  0 37 42 0 0 37");

   mpz_poly_set_coeff_ui(poly, 5, 42);
   success = success && mpz_poly_equal_str(poly, "6  0 37 42 0 0 42");

   mpz_poly_set_coeff_ui(poly, 2, 0);
   success = success && mpz_poly_equal_str(poly, "6  0 37 0 0 0 42");

   mpz_poly_set_coeff_ui(poly, 8, 0);
   success = success && mpz_poly_equal_str(poly, "6  0 37 0 0 0 42");

   mpz_poly_set_coeff_ui(poly, 5, 0);
   success = success && mpz_poly_equal_str(poly, "2  0 37");
   
   mpz_poly_truncate(poly, poly, 1);
   success = success && mpz_poly_equal_str(poly, "0");

   mpz_poly_set_coeff_ui(poly, 3, 42);
   success = success && mpz_poly_equal_str(poly, "4  0 0 0 42");

   mpz_poly_set_coeff_ui(poly, 3, 0);
   success = success && mpz_poly_equal_str(poly, "0");

   mpz_poly_clear(poly);
   return success;
}


int test_mpz_poly_set_coeff_si()
{
   int success = 1;
   mpz_poly_t poly;
   mpz_poly_init(poly);

   mpz_poly_set_coeff_si(poly, 2, 42);
   success = success && mpz_poly_equal_str(poly, "3  0 0 42");

   mpz_poly_set_coeff_si(poly, 5, -37);
   success = success && mpz_poly_equal_str(poly, "6  0 0 42 0 0 -37");

   mpz_poly_set_coeff_si(poly, 1, -37);
   success = success && mpz_poly_equal_str(poly, "6  0 -37 42 0 0 -37");

   mpz_poly_set_coeff_si(poly, 5, 42);
   success = success && mpz_poly_equal_str(poly, "6  0 -37 42 0 0 42");

   mpz_poly_set_coeff_si(poly, 2, 0);
   success = success && mpz_poly_equal_str(poly, "6  0 -37 0 0 0 42");

   mpz_poly_set_coeff_si(poly, 8, 0);
   success = success && mpz_poly_equal_str(poly, "6  0 -37 0 0 0 42");

   mpz_poly_set_coeff_si(poly, 5, 0);
   success = success && mpz_poly_equal_str(poly, "2  0 -37");
   
   mpz_poly_truncate(poly, poly, 1);
   success = success && mpz_poly_equal_str(poly, "0");

   mpz_poly_set_coeff_si(poly, 3, 42);
   success = success && mpz_poly_equal_str(poly, "4  0 0 0 42");

   mpz_poly_set_coeff_si(poly, 3, 0);
   success = success && mpz_poly_equal_str(poly, "0");

   mpz_poly_clear(poly);
   return success;
}


/****************************************************************************

   Conversions

****************************************************************************/


int test_mpz_poly_to_fmpz_poly()
{
   return 0;
}


int test_fmpz_poly_to_mpz_poly()
{
   return 0;
}


/****************************************************************************

   String conversions and I/O

****************************************************************************/


int test_mpz_poly_from_string()
{
   return 0;
/*
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
*/
}



int test_mpz_poly_to_string()
{
   return 0;
/*
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
*/
}


int test_mpz_poly_fprint()
{
   return 0;
}


int test_mpz_poly_fread()
{
   return 0;
}


/****************************************************************************

   Length and degree

****************************************************************************/


int test_mpz_poly_normalise()
{
   return 0;
/*
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
*/
}


int test_mpz_poly_normalised()
{
   return 0;
}


int test_mpz_poly_pad()
{
   return 0;
}


int test_mpz_poly_truncate()
{
   return 0;
}


int test_mpz_poly_length()
{
   return 0;
}


int test_mpz_poly_degree()
{
   return 0;
}


/****************************************************************************

   Assignment

****************************************************************************/

int test_mpz_poly_set()
{
   return 0;
/*
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
*/
}


int test_mpz_poly_swap()
{
   return 0;
}


/****************************************************************************

   Comparison

****************************************************************************/


int test_mpz_poly_equal()
{
   return 0;
/*
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
*/
}


/****************************************************************************

   Addition/subtraction

****************************************************************************/

int test_mpz_poly_add()
{
   return 0;
/*
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
*/
}


int test_mpz_poly_sub()
{
   return 0;
/*
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
*/
}


int test_mpz_poly_neg()
{
   return 0;
/*
   int success = 1;
   mpz_poly_t poly1, poly2;
   mpz_poly_init2(poly1, 10);
   mpz_poly_init2(poly2, 10);

   // out-of-place

   mpz_poly_set_from_string(poly1, "");
   mpz_poly_set_from_string(poly2, "123 456 789 123 456");
   _mpz_poly_neg(poly2, poly1);
   success = success && mpz_poly_equal_str(poly2, "");

   mpz_poly_set_from_string(poly1, "0 2 -5 6");
   mpz_poly_set_from_string(poly2, "123 456 789 123 456");
   _mpz_poly_neg(poly2, poly1);
   success = success && mpz_poly_equal_str(poly2, "0 -2 5 -6");

   // in-place

   mpz_poly_set_from_string(poly1, "");
   _mpz_poly_neg(poly1, poly1);
   success = success && mpz_poly_equal_str(poly1, "");

   mpz_poly_set_from_string(poly1, "0 2 -5 6");
   _mpz_poly_neg(poly1, poly1);
   success = success && mpz_poly_equal_str(poly1, "0 -2 5 -6");

   mpz_poly_clear(poly1);
   mpz_poly_clear(poly2);
   return success;
*/
}



/****************************************************************************

   Shifting

****************************************************************************/


int test_mpz_poly_lshift()
{
   return 0;
}


int test_mpz_poly_rshift()
{
   return 0;
}


int test_mpz_poly_shift()
{
   return 0;
}


/****************************************************************************

   Scalar multiplication and division

****************************************************************************/


int test_mpz_poly_scalar_mul()
{
   return 0;
}


int test_mpz_poly_scalar_mul_ui()
{
   return 0;
}


int test_mpz_poly_scalar_mul_si()
{
   return 0;
}


int test_mpz_poly_scalar_div()
{
   return 0;
}


int test_mpz_poly_scalar_div_ui()
{
   return 0;
}


int test_mpz_poly_scalar_div_si()
{
   return 0;
}


int test_mpz_poly_scalar_div_exact()
{
   return 0;
}


int test_mpz_poly_scalar_div_exact_ui()
{
   return 0;
}


int test_mpz_poly_scalar_div_exact_si()
{
   return 0;
}


int test_mpz_poly_mod()
{
   return 0;
}


int test_mpz_poly_mod_ui()
{
   return 0;
}


/****************************************************************************

   Polynomial multiplication

****************************************************************************/


int test_mpz_poly_mul()
{
   // also should test squaring
   return 0;
}


int test_mpz_poly_mul_naive()
{
   // also should test squaring
   return 0;
/*
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
*/
}


int test_mpz_poly_mul_karatsuba()
{
   // also should test squaring
   return 0;
}


int test_mpz_poly_mul_SS()
{
   // also should test squaring
   return 0;
}


int test_mpz_poly_mul_naive_KS()
{
   // also should test squaring
   return 0;
/*
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
                                       randstate, bitsize[j]+1);
                        mpz_rrandomb(temp, randstate, bits);
                        if (gmp_urandomb_ui(randstate, 1))
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
*/
}





/****************************************************************************

   Polynomial division

****************************************************************************/

int test_mpz_poly_monic_inverse()
{
   return 0;
/*
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
            mpz_poly_set_coeff_si(poly1, i, gmp_urandomb_ui(randstate, 10) - 512);
      
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
*/
}


int test_mpz_poly_pseudo_inverse()
{
   return 0;
}


int test_mpz_poly_monic_div()
{
   return 0;
}


int test_mpz_poly_pseudo_div()
{
   return 0;
}


int test_mpz_poly_monic_rem()
{
   return 0;
}


int test_mpz_poly_pseudo_rem()
{
   return 0;
}


int test_mpz_poly_monic_div_rem()
{
   return 0;
}


int test_mpz_poly_pseudo_div_rem()
{
   return 0;
}


int test_mpz_poly_monic_inverse_naive()
{
   return 0;
}


int test_mpz_poly_pseudo_inverse_naive()
{
   return 0;
}


int test_mpz_poly_monic_div_naive()
{
   return 0;
}


int test_mpz_poly_pseudo_div_naive()
{
   return 0;
}


int test_mpz_poly_monic_rem_naive()
{
   return 0;
}


int test_mpz_poly_pseudo_rem_naive()
{
   return 0;
}


int test_mpz_poly_monic_div_rem_naive()
{
   return 0;
}


int test_mpz_poly_pseudo_div_rem_naive()
{
   return 0;
}


/****************************************************************************

   GCD and extended GCD

****************************************************************************/


int test_mpz_poly_content()
{
   return 0;
}


int test_mpz_poly_content_ui()
{
   return 0;
}


int test_mpz_poly_gcd()
{
   return 0;
}


int test_mpz_poly_xgcd()
{
   return 0;
}


/****************************************************************************

   Miscellaneous

****************************************************************************/

int test_mpz_poly_max_limbs()
{
   return 0;
}


int test_mpz_poly_max_bits()
{
   return 0;
}



/****************************************************************************

   Main test code

****************************************************************************/

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void mpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(mpz_poly_get_coeff_ptr);
   RUN_TEST(mpz_poly_get_coeff);
   RUN_TEST(mpz_poly_get_coeff_ui);
   RUN_TEST(mpz_poly_get_coeff_si);
   RUN_TEST(mpz_poly_set_coeff);
   RUN_TEST(mpz_poly_set_coeff_ui);
   RUN_TEST(mpz_poly_set_coeff_si);
/*
   RUN_TEST(mpz_poly_to_fmpz_poly);
   RUN_TEST(fmpz_poly_to_mpz_poly);
   RUN_TEST(mpz_poly_from_string);
   RUN_TEST(mpz_poly_to_string);
   RUN_TEST(mpz_poly_fprint);
   RUN_TEST(mpz_poly_fread);
   RUN_TEST(mpz_poly_normalise);
   RUN_TEST(mpz_poly_normalised);
   RUN_TEST(mpz_poly_pad);
   RUN_TEST(mpz_poly_length);
   RUN_TEST(mpz_poly_degree);
   RUN_TEST(mpz_poly_set);
   RUN_TEST(mpz_poly_swap);
   RUN_TEST(mpz_poly_equal);
   RUN_TEST(mpz_poly_add);
   RUN_TEST(mpz_poly_sub);
   RUN_TEST(mpz_poly_neg);
   RUN_TEST(mpz_poly_lshift);
   RUN_TEST(mpz_poly_rshift);
   RUN_TEST(mpz_poly_shift);
   RUN_TEST(mpz_poly_scalar_mul);
   RUN_TEST(mpz_poly_scalar_mul_ui);
   RUN_TEST(mpz_poly_scalar_mul_si);
   RUN_TEST(mpz_poly_scalar_div);
   RUN_TEST(mpz_poly_scalar_div_ui);
   RUN_TEST(mpz_poly_scalar_div_si);
   RUN_TEST(mpz_poly_scalar_div_exact);
   RUN_TEST(mpz_poly_scalar_div_exact_ui);
   RUN_TEST(mpz_poly_scalar_div_exact_si);
   RUN_TEST(mpz_poly_mod);
   RUN_TEST(mpz_poly_mod_ui);
   RUN_TEST(mpz_poly_mul);
   RUN_TEST(mpz_poly_mul_naive);
   RUN_TEST(mpz_poly_mul_karatsuba);
   RUN_TEST(mpz_poly_mul_SS);
   RUN_TEST(mpz_poly_mul_naive_KS);
   RUN_TEST(mpz_poly_monic_inverse);
   RUN_TEST(mpz_poly_monic_inverse);
   RUN_TEST(mpz_poly_pseudo_inverse);
   RUN_TEST(mpz_poly_monic_div);
   RUN_TEST(mpz_poly_pseudo_div);
   RUN_TEST(mpz_poly_monic_rem);
   RUN_TEST(mpz_poly_pseudo_rem);
   RUN_TEST(mpz_poly_monic_div_rem);
   RUN_TEST(mpz_poly_pseudo_div_rem);
   RUN_TEST(mpz_poly_monic_inverse_naive);
   RUN_TEST(mpz_poly_pseudo_inverse_naive);
   RUN_TEST(mpz_poly_monic_div_naive);
   RUN_TEST(mpz_poly_pseudo_div_naive);
   RUN_TEST(mpz_poly_monic_rem_naive);
   RUN_TEST(mpz_poly_pseudo_rem_naive);
   RUN_TEST(mpz_poly_monic_div_rem_naive);
   RUN_TEST(mpz_poly_pseudo_div_rem_naive);
   RUN_TEST(mpz_poly_content);
   RUN_TEST(mpz_poly_content_ui);
   RUN_TEST(mpz_poly_gcd);
   RUN_TEST(mpz_poly_xgcd);
   RUN_TEST(mpz_poly_max_limbs);
   RUN_TEST(mpz_poly_max_bits);
*/
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   test_support_init();
   mpz_poly_test_all();
   test_support_cleanup();

   return 0;
}


// *************** end of file
