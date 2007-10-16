/****************************************************************************

fmpz-test.c: Test code for fmpz.c and fmpz.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "test-support.h"
#include "fmpz.h"

#define SIGNS 1

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

gmp_randstate_t state;

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
int test_fmpz_convert()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
       fmpz_to_mpz(num2, fnum1);
       
       fmpz_clear(fnum1);
       
       result = (mpz_cmp(num1, num2) == 0);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_size()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
              
       result = (mpz_size(num1) == fmpz_size(fnum1));

       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_bits()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);

       result = (mpz_sizeinbase(num1, 2) == fmpz_bits(fnum1)) 
                  || ((mpz_cmp_ui(num1, 0) == 0) && (fmpz_bits(fnum1) == 0));
       
#if DEBUG2
       if (!result)
       {
          printf("bits = %ld, bits2 = %ld\n", mpz_sizeinbase(num1, 2), fmpz_bits(fnum1));
          gmp_printf("%Zd\n", num1);
       }
#endif

       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_sgn()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
              
       result = (((long) mpz_sgn(num1) > 0) && ((long) fmpz_sgn(fnum1) > 0))
             || (((long) mpz_sgn(num1) < 0) && ((long) fmpz_sgn(fnum1) < 0))
             || (((long) mpz_sgn(num1) == 0) && ((long) fmpz_sgn(fnum1) == 0));
#if DEBUG2
       if (!result)
       {
          printf("sign = %ld, sign2 = %ld\n", mpz_sgn(num1), fmpz_sgn(fnum1));
          gmp_printf("%Zd\n", num1);
       }
#endif

       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_set_si()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits;
   long x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   fnum1 = fmpz_init(0);
   fmpz_set_si(fnum1, 0);
   mpz_set_si(num1, 0);
   fmpz_to_mpz(num2, fnum1);
   
   result = (mpz_cmp(num1, num2) == 0);
   
   fmpz_clear(fnum1);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(63)+1;
       x = random_ulong(1L<<bits);
       
#if SIGNS
       if (random_ulong(2)) x = -x;
#endif
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       
       mpz_set_si(num1, x);
       fmpz_set_si(fnum1, x);
       fmpz_to_mpz(num2, fnum1);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_set_ui()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   fnum1 = fmpz_init(0);
   fmpz_set_ui(fnum1, 0);
   mpz_set_ui(num1, 0);
   fmpz_to_mpz(num2, fnum1);
   
   result = (mpz_cmp(num1, num2) == 0);
   
   fmpz_clear(fnum1);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(63)+1;
       x = random_ulong(1L<<bits);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       
       mpz_set_ui(num1, x);
       fmpz_set_ui(fnum1, x);
       fmpz_to_mpz(num2, fnum1);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_set_equal()
{
   mpz_t num1, num2;
   fmpz_t fnum1, fnum2;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
       
       fmpz_set(fnum2, fnum1);
              
       result = (fmpz_equal(fnum1, fnum2));
       
#if DEBUG2
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear(fnum1);
   }
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       do
       {
          bits2 = random_ulong(1000);

          mpz_rrandomb(num2, state, bits2);
#if SIGNS
          if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       } while (mpz_cmp(num1, num2) == 0);
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
              
       result = (!fmpz_equal(fnum1, fnum2));
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_add()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t fnum1, fnum2, fnum3;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       bits2 = random_ulong(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
       fnum3 = fmpz_init(FLINT_MAX((long)bits2/FLINT_BITS, (long)bits/FLINT_BITS)+1);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       fmpz_add(fnum3, fnum1, fnum2);
       mpz_add(num4, num1, num2);
       
       fmpz_to_mpz(num3, fnum3);
              
       result = (mpz_cmp(num3, num4) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
       fmpz_clear(fnum3);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   
   return result;
}

int test_fmpz_sub()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t fnum1, fnum2, fnum3;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       bits2 = random_ulong(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
       fnum3 = fmpz_init(FLINT_MAX((long)bits2/FLINT_BITS, (long)bits/FLINT_BITS)+1);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       fmpz_sub(fnum3, fnum1, fnum2);
       mpz_sub(num4, num1, num2);
       
       fmpz_to_mpz(num3, fnum3);
              
       result = (mpz_cmp(num3, num4) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
       fmpz_clear(fnum3);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   
   return result;
}

int test_fmpz_mul()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t fnum1, fnum2, fnum3;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       bits2 = random_ulong(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
       fnum3 = fmpz_init(FLINT_MAX((long)(bits+bits2-1)/FLINT_BITS,0)+1);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       fmpz_mul(fnum3, fnum1, fnum2);
       mpz_mul(num4, num1, num2);
       
       fmpz_to_mpz(num3, fnum3);
              
       result = (mpz_cmp(num3, num4) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
       fmpz_clear(fnum3);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   
   return result;
}

int test___fmpz_mul()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t fnum1, fnum2, fnum3;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       bits2 = random_ulong(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
       fnum3 = fmpz_init(FLINT_MAX((long)(bits+bits2-1)/FLINT_BITS,0)+2);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       __fmpz_mul(fnum3, fnum1, fnum2);
       mpz_mul(num4, num1, num2);
       
       fmpz_to_mpz(num3, fnum3);
              
       result = (mpz_cmp(num3, num4) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
       fmpz_clear(fnum3);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   
   return result;
}

int test_fmpz_addmul()
{
   mpz_t num1, num2, num3, num4, num5;
   fmpz_t fnum1, fnum2, fnum3, fnum4;
   unsigned long bits, bits2, bits3;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
   mpz_init(num5);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       bits2 = random_ulong(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       
       bits3 = random_ulong(1000);

       mpz_rrandomb(num3, state, bits3);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num3, num3);
#endif
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
       fnum3 = fmpz_init(FLINT_MAX((long)(bits3-1)/FLINT_BITS,0)+1);
       fnum4 = fmpz_init(FLINT_MAX((long)bits3/FLINT_BITS, (long)(bits+bits2)/FLINT_BITS)+1);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       mpz_to_fmpz(fnum3, num3);
       
       fmpz_set(fnum4, fnum3);
       fmpz_addmul(fnum4, fnum1, fnum2);
       mpz_set(num4, num3);
       mpz_addmul(num4, num1, num2);
       
       fmpz_to_mpz(num5, fnum4);
              
       result = (mpz_cmp(num5, num4) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
       fmpz_clear(fnum3);
       fmpz_clear(fnum4);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   mpz_clear(num5);
   
   return result;
}

int test_fmpz_div()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t fnum1, fnum2, fnum3;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       do
       {
          bits2 = random_ulong(1000);

          mpz_rrandomb(num2, state, bits2);
#if SIGNS
          if (random_ulong(2)) mpz_neg(num2, num2);
#endif
       } while (mpz_cmp_ui(num2, 0) == 0);
       
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1)/FLINT_BITS,0)+1);
       fnum2 = fmpz_init(FLINT_MAX((long)(bits2-1)/FLINT_BITS,0)+1);
       fnum3 = fmpz_init(FLINT_MAX((long)(bits-bits2)/FLINT_BITS,0)+2);

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       fmpz_div(fnum3, fnum1, fnum2);
       mpz_tdiv_q(num4, num1, num2);
       
       fmpz_to_mpz(num3, fnum3);
              
       result = (mpz_cmp(num3, num4) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
       fmpz_clear(fnum3);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   
   return result;
}

int test_fmpz_add_ui_inplace()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits, bits2, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(63)+1;
       x = random_ulong(1L<<bits);
       
       bits2 = random_ulong(1000);
       fnum1 = fmpz_init((long) FLINT_MAX(bits, bits2)/FLINT_BITS+1);
       
       mpz_rrandomb(num1, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       mpz_to_fmpz(fnum1, num1);
       fmpz_add_ui_inplace(fnum1, x);
       mpz_add_ui(num1, num1, x);
       fmpz_to_mpz(num2, fnum1);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test___fmpz_add_ui_inplace()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits, bits2, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(63)+1;
       x = random_ulong(1L<<bits);
       
       bits2 = random_ulong(1000);
       fnum1 = fmpz_init((long) FLINT_MAX(bits, bits2)/FLINT_BITS+1);
       
       mpz_rrandomb(num1, state, bits2);

       mpz_to_fmpz(fnum1, num1);
       __fmpz_add_ui_inplace(fnum1, x);
       mpz_add_ui(num1, num1, x);
       fmpz_to_mpz(num2, fnum1);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_sub_ui_inplace()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long bits, bits2, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(63)+1;
       x = random_ulong(1L<<bits);
       
       bits2 = random_ulong(1000);
       fnum1 = fmpz_init((long) FLINT_MAX(bits, bits2)/FLINT_BITS+1);
       
       mpz_rrandomb(num1, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       mpz_to_fmpz(fnum1, num1);
       fmpz_sub_ui_inplace(fnum1, x);
       mpz_sub_ui(num1, num1, x);
       fmpz_to_mpz(num2, fnum1);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_mul_ui()
{
   mpz_t num1, num2;
   fmpz_t fnum1, fnum2;
   unsigned long bits, bits2, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(63)+1;
       x = random_ulong(1L<<bits);
       
       bits2 = random_ulong(1000);
       fnum1 = fmpz_init((long) FLINT_MAX(bits2-1, 0)/FLINT_BITS+1);
       fnum2 = fmpz_init((long) FLINT_MAX(bits2-1, 0)/FLINT_BITS+2);
       
       mpz_rrandomb(num1, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       mpz_to_fmpz(fnum1, num1);
       fmpz_mul_ui(fnum2, fnum1, x);
       mpz_mul_ui(num1, num1, x);
       fmpz_to_mpz(num2, fnum2);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_div_ui()
{
   mpz_t num1, num2;
   fmpz_t fnum1, fnum2;
   unsigned long bits, bits2, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       do
       {
          bits = random_ulong(63)+1;
          x = random_ulong(1L<<bits);
       } while (x == 0);
       
       bits2 = random_ulong(1000);
       fnum1 = fmpz_init((long) FLINT_MAX(bits2-1, 0)/FLINT_BITS+1);
       fnum2 = fmpz_init((long) FLINT_MAX(bits2-1, 0)/FLINT_BITS+1);
       
       mpz_rrandomb(num1, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       mpz_to_fmpz(fnum1, num1);
       fmpz_div_ui(fnum2, fnum1, x);
       mpz_tdiv_q_ui(num1, num1, x);
       fmpz_to_mpz(num2, fnum2);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_pow_ui()
{
   mpz_t num1, num2;
   fmpz_t fnum1, fnum2;
   unsigned long bits, bits2, x;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 10000) && (result == 1); i++)
   {
       bits = random_ulong(10)+1;
       x = random_ulong(1L<<bits);
       
       bits2 = random_ulong(150);
       fnum1 = fmpz_init((long) FLINT_MAX(bits2-1, 0)/FLINT_BITS+1);
       fnum2 = fmpz_init((long) FLINT_MAX(bits2*x-1, 0)/FLINT_BITS+1);
       
       mpz_rrandomb(num1, state, bits2);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif

       mpz_to_fmpz(fnum1, num1);
       fmpz_pow_ui(fnum2, fnum1, x);
       mpz_pow_ui(num1, num1, x);
       fmpz_to_mpz(num2, fnum2);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_is_one()
{
   mpz_t num1;
   fmpz_t fnum1;
   unsigned long bits;
   int result = 1;
   int test1, test2;
   
   mpz_init(num1);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(random_ulong(1000)+1);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1), 0)/FLINT_BITS+1);
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
       test1 = (mpz_cmp_ui(num1, 1L) == 0);
       test2 = fmpz_is_one(fnum1);
               
       result = (test1 == test2);
#if DEBUG
       if (!result) gmp_printf("%Zd, %d\n", num1, test1);
#endif
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   
   return result;
}

int test_fmpz_is_zero()
{
   mpz_t num1;
   fmpz_t fnum1;
   unsigned long bits;
   int result = 1;
   int test1, test2;
   
   mpz_init(num1);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = random_ulong(random_ulong(1000)+1);
       fnum1 = fmpz_init(FLINT_MAX((long)(bits-1), 0)/FLINT_BITS+1);
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
       test1 = (mpz_cmp_ui(num1, 0L) == 0);
       test2 = fmpz_is_zero(fnum1);
               
       result = (test1 == test2);
#if DEBUG
       if (!result) gmp_printf("%Zd, %d\n", num1, test1);
#endif
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   
   return result;
}

int test_fmpz_normalise()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long limbs, limbs2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 50000) && (result == 1); i++)
   {
       limbs = random_ulong(1000);
       limbs2 = random_ulong(1000);
       fnum1 = fmpz_init(limbs+limbs2);
       mpz_random(num1, limbs);
#if SIGNS
       if (random_ulong(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
       
       for (unsigned long j = FLINT_ABS(fnum1[0])+1; j < limbs + limbs2 + 1; j++)
       {
          fnum1[j] = 0L;
       }
       fnum1[0] = limbs + limbs2;
       fmpz_normalise(fnum1);
              
       result = (mpz_size(num1) == fmpz_size(fnum1));

       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_binomial_next()
{
   mpz_t num1, num2;
   fmpz_t fnum1;
   unsigned long n, m;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 10000) && (result == 1); i++)
   {
       n = random_ulong(1000);
       m = random_ulong(n+1);
       
       fnum1 = fmpz_init((long) FLINT_MAX(n-1, 0)/FLINT_BITS+2);
       
       fmpz_set_ui(fnum1, 1L);
       for (long j = 1; j <= m; j++)
          fmpz_binomial_next(fnum1, fnum1, n, j);
       mpz_bin_uiui(num1, n, m);
       fmpz_to_mpz(num2, fnum1);
             
       result = (mpz_cmp(num1, num2) == 0);
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(fmpz_convert);
   RUN_TEST(fmpz_size);
   RUN_TEST(fmpz_bits);
   RUN_TEST(fmpz_sgn);
   RUN_TEST(fmpz_set_si);
   RUN_TEST(fmpz_set_ui);
   RUN_TEST(fmpz_set_equal);
   RUN_TEST(fmpz_add);
   RUN_TEST(fmpz_add_ui_inplace);
   RUN_TEST(__fmpz_add_ui_inplace);
   RUN_TEST(fmpz_sub);
   RUN_TEST(fmpz_sub_ui_inplace);
   RUN_TEST(fmpz_mul);
   RUN_TEST(fmpz_mul_ui);
   RUN_TEST(__fmpz_mul);
   RUN_TEST(fmpz_addmul);
   RUN_TEST(fmpz_div);
   RUN_TEST(fmpz_div_ui);
   RUN_TEST(fmpz_pow_ui);
   RUN_TEST(fmpz_is_one);
   RUN_TEST(fmpz_is_zero);
   RUN_TEST(fmpz_normalise);
   RUN_TEST(fmpz_binomial_next);
      
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   gmp_randinit_default(state);
   fmpz_poly_test_all();
   gmp_randclear(state);
   test_support_cleanup();

   return 0;
}
