/****************************************************************************

Zpoly-test.c: Test code for Zpoly.c and Zpoly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "Zpoly.h"
#include <stdio.h>



// all test functions return 1 on success, 0 on failure


int test_Zpoly_mpz_get_coeff()
{
   Zpoly_mpz_t poly;
   mpz_t x;

   Zpoly_mpz_init2(poly, 3);
   mpz_init(x);
   int success = 1;
   
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
   Zpoly_mpz_t poly;
   mpz_t x;

   Zpoly_mpz_init2(poly, 3);
   mpz_init(x);
   int success = 1;
   
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
   Zpoly_mpz_t poly;
   mpz_t x;

   Zpoly_mpz_init2(poly, 3);
   mpz_init(x);
   int success = 1;
   
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



#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!!\n");


void Zpoly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(Zpoly_mpz_get_coeff);
   RUN_TEST(Zpoly_mpz_get_coeff_ui);
   RUN_TEST(Zpoly_mpz_get_coeff_si);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!!\n");
}



int main()
{

   Zpoly_test_all();

/*
   char *s;

   Zpoly_mpz_t poly1, poly2, poly3;
   Zpoly_mpz_init(poly1);
   Zpoly_mpz_init(poly2);
   Zpoly_mpz_init(poly3);

   Zpoly_mpz_set_from_string(poly1, " -1 2 0  34   -123123123123123123 34 ");
   Zpoly_mpz_set_from_string(poly2, "0 49 28");
   Zpoly_mpz_add(poly3, poly1, poly2);

   s = malloc(Zpoly_mpz_get_string_size(poly1));
   Zpoly_mpz_get_as_string(s, poly1);
   printf(" first poly is: %s\n", s);
   free(s);
   
   s = malloc(Zpoly_mpz_get_string_size(poly2));
   Zpoly_mpz_get_as_string(s, poly2);
   printf("second poly is: %s\n", s);
   free(s);
   
   s = malloc(Zpoly_mpz_get_string_size(poly3));
   Zpoly_mpz_get_as_string(s, poly3);
   printf("  their sum is: %s\n", s);
   free(s);
   
   Zpoly_mpz_clear(poly1);
   Zpoly_mpz_clear(poly2);
   Zpoly_mpz_clear(poly3);
*/

   return 0;
}


// *************** end of file
