/****************************************************************************

Zpoly_mpn-test.c: Test code for Zpoly_mpn.c and Zpoly_mpn.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"
#include "flint-manager.h"

#define DEBUG 0 // prints debug information

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
   int result = 1;
   unsigned long bits, length;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      _Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          _Zpoly_mpn_realloc(test_mpn_poly, length);
          Zpoly_realloc(test_poly, length);
          Zpoly_realloc(test_poly2, length);
          test_poly->length = length;
          for (unsigned long i = 0; i < length; i++)
          {
             mpz_rrandomb(temp, Zpoly_test_randstate, bits);
             Zpoly_set_coeff(test_poly, i, temp);
          } 
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");*/
#endif
          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          _Zpoly_mpn_convert_out(test_poly2, test_mpn_poly);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = Zpoly_equal(test_poly, test_poly2);
      }   
      _Zpoly_mpn_clear(test_mpn_poly);
   }
   
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   
   return result;
}

void Zpoly_mpn_test_all()
{
   int success, all_success = 1;

   RUN_TEST(Zpoly_mpn_convert);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   Zpoly_mpn_test_all();

   return 0;
}


