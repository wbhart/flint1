/****************************************************************************

Zpoly-test.c: Test code for Zpoly.c and Zpoly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"
#include "flint-manager.h"

gmp_randstate_t Zpoly_test_randstate;

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
int test_Zpoly_mpn_convert()
{
   Zpoly_t test_poly, test_poly2;
   Zpoly_mpn_t test_mpn_poly;
   mpz_t temp;
   mpz_init(temp);
   int result;
   
   Zpoly_init3(test_poly, 100, 100); //alloc, coeff_bits
   Zpoly_init3(test_poly2, 100, 100); //alloc, coeff_bits
   _Zpoly_mpn_init(test_mpn_poly, 100, 100/FLINT_BITS_PER_LIMB+1);
   test_poly->length = 100;
   for (unsigned long i = 0; i < 100; i++)
   {
      mpz_rrandomb(temp, Zpoly_test_randstate, 100);
      Zpoly_set_coeff(test_poly, i, temp);
   } 
   _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
   _Zpoly_mpn_convert_out(test_poly2, test_mpn_poly);
   
   result = Zpoly_equal(test_poly, test_poly2);
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   _Zpoly_mpn_clear(test_mpn_poly);
   return result;
}

void Zpoly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(Zpoly_mpn_convert);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   Zpoly_test_all();

   return 0;
}


