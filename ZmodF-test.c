/******************************************************************************

 ZmodF-test.c: test module for ZmodF module

 Copyright (C) 2007, David Harvey
 
******************************************************************************/

#include <stdio.h>
#include "ZmodF.h"

gmp_randstate_t ZmodF_test_randstate;
mpz_t global_mpz;   // to avoid frequent mpz_init calls


/*
Prints the ZmodF_t to the given stream in hex, each limb in a separate block,
most significant limb (i.e. the overflow limb) first.
*/
void ZmodF_print(FILE* stream, ZmodF_t x, unsigned long n)
{
   for (long i = n; i >= 0; i--)
#if FLINT_BITS_PER_LIMB == 64
      fprintf(stream, "%016lx ", x[i]);
#else
      fprintf(stream, "%08lx ", x[i]);
#endif
}


/*
Generates a random ZmodF_t with at most overflow_bits used in the
overflow limb. We use mpz_rrandomb to get long strings of 0's and 1's.
*/
void ZmodF_random(ZmodF_t x, unsigned long n, unsigned long overflow_bits)
{
   ZmodF_zero(x, n);

   mpz_rrandomb(global_mpz, ZmodF_test_randstate, (n+1)*FLINT_BITS_PER_LIMB);
   mpz_export(x, NULL, -1, sizeof(mp_limb_t), 0, 0, global_mpz);

   // GMP has a "bug" where the top bit of the output of mpz_rrandomb
   // is always set. So we flip everything with probability 1/2.
   if (gmp_urandomb_ui(ZmodF_test_randstate, 1))
      for (unsigned long i = 0; i <= n; i++)
         x[i] = ~x[i];

   // Now copy the sign bit downwards so that only overflow_bits bits are used.
   if ((mp_limb_signed_t) x[n] >= 0)
      x[n] &= (1UL << overflow_bits) - 1;
   else
      x[n] |= ~((1UL << overflow_bits) - 1);
}



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
   mpz_init(global_mpz);
   
   ZmodF_test_all();
   
   return 0;
}



// end of file ****************************************************************
