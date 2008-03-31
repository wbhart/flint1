/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

fmpz-test.c: Test code for fmpz.c and fmpz.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "mpir.h"
#include "test_support.h"
#include "memory_manager.h"
#include "fmpz.h"
#include "mpn_extras.h"

#define SIGNS 1

#define DEBUG 0 // prints debug information
#define DEBUG2 1 
#define SIGNS 1

gmp_randstate_t state;

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
int test_fmpz_block_init_realloc_clear()
{
   int result = 1;
   table_entry entry[1];
   
   for (ulong i = 0; (i < 10000) && (result == 1); i++)
   {
       ulong n = randint(100)+1;
       fmpz_block_init2(entry, n);
       mpn_random2(entry->block_ptr, (n+1)*MPIR_BLOCK);
       mp_limb_t * temp = (mp_limb_t *) mpir_alloc((n+1)*MPIR_BLOCK*sizeof(mp_limb_t));
       F_mpn_copy(temp, entry->block_ptr, (n+1)*MPIR_BLOCK);
       fmpz_block_realloc(entry, n+randint(100)+1);
       fmpz_block_realloc(entry, n);
       result = (mpn_cmp(entry->block_ptr, temp, (n+1)*MPIR_BLOCK) == 0);
       fmpz_block_clear(entry);
       mpir_free(temp);
   }
   
   for (ulong i = 0; (i < 10000) && (result == 1); i++)
   {
       ulong n = randint(100)+1;
       fmpz_block_init(entry);
       fmpz_block_realloc(entry, n);
       mpn_random2(entry->block_ptr, (n+1)*MPIR_BLOCK);
       mp_limb_t * temp = (mp_limb_t *) mpir_alloc((n+1)*MPIR_BLOCK*sizeof(mp_limb_t));
       F_mpn_copy(temp, entry->block_ptr, (n+1)*MPIR_BLOCK);
       fmpz_block_realloc(entry, n+randint(100)+1);
       fmpz_block_realloc(entry, n);
       result = (mpn_cmp(entry->block_ptr, temp, (n+1)*MPIR_BLOCK) == 0);
       fmpz_block_clear(entry);
       mpir_free(temp);
   }
   
   for (ulong i = 0; (i < 10000) && (result == 1); i++)
   {
       ulong n = randint(100)+1;
       ulong m = randint(MPIR_BLOCK)+1;
       fmpz_block_init2_small(entry, m, n);
       mpn_random2(entry->block_ptr, m*(n+1));
       mp_limb_t * temp = (mp_limb_t *) mpir_alloc(m*(n+1)*sizeof(mp_limb_t));
       F_mpn_copy(temp, entry->block_ptr, m*(n+1));
       fmpz_block_realloc(entry, n+randint(100)+1);
       fmpz_block_realloc(entry, n);
       result = (mpn_cmp(entry->block_ptr, temp, m*(n+1)) == 0);
       fmpz_block_clear(entry);
       mpir_free(temp);
   }
   
   for (ulong i = 0; (i < 10000) && (result == 1); i++)
   {
       ulong n = randint(100)+1;
       ulong m = randint(MPIR_BLOCK)+1;
       fmpz_block_init_small(entry, m);
       fmpz_block_realloc(entry, n);
       mpn_random2(entry->block_ptr, m*(n+1));
       mp_limb_t * temp = (mp_limb_t *) mpir_alloc(m*(n+1)*sizeof(mp_limb_t));
       F_mpn_copy(temp, entry->block_ptr, m*(n+1));
       fmpz_block_realloc(entry, n+randint(100)+1);
       fmpz_block_realloc(entry, n);
       result = (mpn_cmp(entry->block_ptr, temp, m*(n+1)) == 0);
       fmpz_block_clear(entry);
       mpir_free(temp);
   }
   
   return result;
}


int test_fmpz_init_fit_limbs_clear()
{
   int result = 1;
   fmpz_t * int1;

   for (ulong i = 0; (i < 100) && (result == 1); i++)
   {
       int1 = fmpz_init();
       for (ulong j = 0; j < 100; j++)
       {
          fmpz_fit_limbs(int1, randint(100)+1);
       }
       fmpz_clear(int1);
   }
   
   return result;
}

int test_fmpz_init_fit_limbs_clear_array()
{
   int result = 1;
   fmpz_t * int_arr;

   for (ulong i = 0; (i < 100) && (result == 1); i++)
   {
       ulong num_ints = randint(200)+1;
       int_arr = fmpz_init_array(num_ints);
       for (ulong j = 0; j < 1000; j++)
       {
          fmpz_fit_limbs(int_arr + randint(num_ints), randint(100)+1);
       }
       fmpz_clear_array(int_arr, num_ints);
   }
   
   return result;
}

int test_fmpz_to_mpz()
{
   mpz_t num1, num2;
   fmpz_t * fnum1;
   ulong bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);
#if DEBUG
       printf("Bits = %ld\n", bits);
#endif
       fnum1 = fmpz_init();
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
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
   fmpz_t * fnum1;
   ulong bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);
       fnum1 = fmpz_init();
       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif
       mpz_to_fmpz(fnum1, num1);
              
       result = (mpz_size(num1) == fmpz_size(fnum1));

       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_get_set_ui()
{
   int result = 1;
   fmpz_t * int_arr, * int1;

   for (ulong i = 0; (i < 100) && (result == 1); i++)
   {
       ulong num_ints = randint(200)+1;
       int_arr = fmpz_init_array(num_ints);
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          ulong coord = randint(num_ints);
          ulong val = randbits(MPIR_BITS);
          fmpz_set_ui(int_arr + coord, val);
          fmpz_fit_limbs(int_arr + coord, randint(100)+1);
          ulong val2 = fmpz_get_ui(int_arr + coord);
          result &= (val == val2);
#if DEBUG 
          if (!result) printf("%ld, %ld\n", val, val2);
#endif
       }
       fmpz_clear_array(int_arr, num_ints);
   }
   
   for (ulong i = 0; (i < 100) && (result == 1); i++)
   {
       int1 = fmpz_init();
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          ulong val = randbits(MPIR_BITS);
          fmpz_set_ui(int1, val);
          fmpz_fit_limbs(int1, randint(100)+1);
          ulong val2 = fmpz_get_ui(int1);
          result &= (val == val2);
#if DEBUG 
          if (!result) printf("%ld, %ld\n", val, val2);
#endif
       }
       fmpz_clear(int1);
   }

   return result;
}

int test_fmpz_neg()
{
   mpz_t num1, num2, num3;
   fmpz_t * fnum1, * fnum2;
   ulong bits;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif

       fnum1 = fmpz_init();
       fnum2 = fmpz_init();

       mpz_to_fmpz(fnum1, num1);
       
       fmpz_neg(fnum2, fnum1);
       mpz_neg(num2, num1);
       
       fmpz_to_mpz(num3, fnum2);
              
       result = (mpz_cmp(num2, num3) == 0);
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   
   return result;
}

int test_fmpz_add()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t * fnum1, * fnum2, * fnum3;
   ulong bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif

       bits2 = randint(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (randint(2)) mpz_neg(num2, num2);
#endif
       
       fnum1 = fmpz_init();
       fnum2 = fmpz_init();
       fnum3 = fmpz_init();

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
   fmpz_t * fnum1, * fnum2, * fnum3;
   ulong bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif

       bits2 = randint(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (randint(2)) mpz_neg(num2, num2);
#endif
       
       fnum1 = fmpz_init();
       fnum2 = fmpz_init();
       fnum3 = fmpz_init();

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

int test_fmpz_addmul_ui()
{
   mpz_t num1, num2;
   fmpz_t * fnum1, * fnum2;
   ulong bits, bits2, bits3;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif

       bits2 = randint(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (randint(2)) mpz_neg(num2, num2);
#endif
       
       bits3 = randint(MPIR_BITS);
       ulong c = randbits(bits3);
#if DEBUG
       printf("Bits = %ld, bits2 = %ld, bits3 = %ld\n", bits, bits2, bits3);
#endif
       
       fnum1 = fmpz_init();
       fnum2 = fmpz_init();

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       fmpz_addmul_ui(fnum1, fnum2, c);
       mpz_addmul_ui(num1, num2, c);
       
       fmpz_to_mpz(num2, fnum1);
              
       result = (mpz_cmp(num1, num2) == 0);
#if DEBUG
       if (!result) gmp_printf("%Zd\n\n%Zd\n\n\n", num1, num2);
#endif
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_submul_ui()
{
   mpz_t num1, num2;
   fmpz_t * fnum1, * fnum2;
   ulong bits, bits2, bits3;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif

       bits2 = randint(1000);

       mpz_rrandomb(num2, state, bits2);
#if SIGNS
       if (randint(2)) mpz_neg(num2, num2);
#endif
       
       bits3 = randint(MPIR_BITS);
       ulong c = randbits(bits3);
#if DEBUG
       printf("Bits = %ld, bits2 = %ld, bits3 = %ld\n", bits, bits2, bits3);
#endif
       
       fnum1 = fmpz_init();
       fnum2 = fmpz_init();

       mpz_to_fmpz(fnum1, num1);
       mpz_to_fmpz(fnum2, num2);
       
       fmpz_submul_ui(fnum1, fnum2, c);
       mpz_submul_ui(num1, num2, c);
       
       fmpz_to_mpz(num2, fnum1);
              
       result = (mpz_cmp(num1, num2) == 0);
#if DEBUG
       if (!result) gmp_printf("%Zd\n\n%Zd\n\n\n", num1, num2);
#endif
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_mul_2exp()
{
   mpz_t num1, num2;
   fmpz_t * fnum1, * fnum2;
   ulong bits, exp;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000);

       mpz_rrandomb(num1, state, bits);
#if SIGNS
       if (randint(2)) mpz_neg(num1, num1);
#endif

       exp = randint(1000);
       
#if DEBUG
       printf("Bits = %ld, bits2 = %ld\n", bits, bits2);
#endif
       
       fnum1 = fmpz_init();
       fnum2 = fmpz_init();

       mpz_to_fmpz(fnum1, num1);
       
       fmpz_mul_2exp(fnum2, fnum1, exp);
       mpz_mul_2exp(num2, num1, exp);
       
       fmpz_to_mpz(num1, fnum2);
              
       result = (mpz_cmp(num1, num2) == 0);
#if DEBUG
       if (!result) gmp_printf("%Zd\n\n%Zd\n\n\n", num1, num2);
#endif
       
       fmpz_clear(fnum1);
       fmpz_clear(fnum2);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(fmpz_block_init_realloc_clear);
   RUN_TEST(fmpz_init_fit_limbs_clear);
   RUN_TEST(fmpz_init_fit_limbs_clear_array);
   RUN_TEST(fmpz_to_mpz);
   RUN_TEST(fmpz_size);
   RUN_TEST(fmpz_get_set_ui);
   RUN_TEST(fmpz_neg);
   RUN_TEST(fmpz_add);
   RUN_TEST(fmpz_sub);
   RUN_TEST(fmpz_addmul_ui);
   RUN_TEST(fmpz_submul_ui);
   RUN_TEST(fmpz_mul_2exp);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(state);
   fmpz_poly_test_all();
   gmp_randclear(state);
   
   mpir_stack_cleanup();

   return 0;
}
