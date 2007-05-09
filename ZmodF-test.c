/******************************************************************************

 ZmodF-test.c: test module for ZmodF module

 Copyright (C) 2007, David Harvey
 
******************************************************************************/

#include <stdio.h>
#include "ZmodF.h"

gmp_randstate_t ZmodF_test_randstate;



int test_ZmodF_normalise()
{
   return 0;
}


int test_ZmodF_fast_reduce()
{
   return 0;
}


int test_ZmodF_neg()
{
   return 0;
}


int test_ZmodF_mul()
{
   return 0;
}


int test_ZmodF_sqr()
{
   return 0;
}


int test_ZmodF_short_div_2exp()
{
   return 0;
}


int test_ZmodF_mul_Bexp()
{
   return 0;
}


int test_ZmodF_div_Bexp_sub()
{
   return 0;
}


int test_ZmodF_div_Bexp_add()
{
   return 0;
}


int test_ZmodF_sub_mul_Bexp()
{
   return 0;
}


int test_ZmodF_mul_pseudosqrt2_n_odd()
{
   return 0;
}


int test_ZmodF_mul_pseudosqrt2_n_even()
{
   return 0;
}


int test_ZmodF_mul_2exp()
{
   return 0;
}


int test_ZmodF_mul_sqrt2exp()
{
   return 0;
}


int test_ZmodF_sub_mul_2exp()
{
   return 0;
}


int test_ZmodF_forward_butterfly_2exp()
{
   return 0;
}


int test_ZmodF_forward_butterfly_sqrt2exp()
{
   return 0;
}


int test_ZmodF_inverse_butterfly_sqrt2exp()
{
   return 0;
}


int test_ZmodF_simple_butterfly()
{
   return 0;
}



#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void ZmodF_test_all()
{
   int success, all_success = 1;

   RUN_TEST(ZmodF_normalise);
   RUN_TEST(ZmodF_fast_reduce);
   RUN_TEST(ZmodF_neg);
   RUN_TEST(ZmodF_mul);
   RUN_TEST(ZmodF_sqr);
   RUN_TEST(ZmodF_short_div_2exp);
   RUN_TEST(ZmodF_mul_Bexp);
   RUN_TEST(ZmodF_div_Bexp_sub);
   RUN_TEST(ZmodF_div_Bexp_add);
   RUN_TEST(ZmodF_sub_mul_Bexp);
   RUN_TEST(ZmodF_mul_pseudosqrt2_n_odd);
   RUN_TEST(ZmodF_mul_pseudosqrt2_n_even);
   RUN_TEST(ZmodF_mul_2exp);
   RUN_TEST(ZmodF_mul_sqrt2exp);
   RUN_TEST(ZmodF_sub_mul_2exp);
   RUN_TEST(ZmodF_forward_butterfly_2exp);
   RUN_TEST(ZmodF_forward_butterfly_sqrt2exp);
   RUN_TEST(ZmodF_inverse_butterfly_sqrt2exp);
   RUN_TEST(ZmodF_simple_butterfly);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(ZmodF_test_randstate);
   ZmodF_test_all();

   return 0;
}



// end of file ****************************************************************
