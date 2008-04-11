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
   fmpz_t entry[MPIR_BLOCK];
   long temp_arr[MPIR_BLOCK];
   mp_limb_t * temp;
       
   for (ulong i = 0; (i < 100000) && (result == 1); i++)
   {
       ulong n = randint(100)+1;
       ulong m = randint(MPIR_BLOCK)+1;
#if DEBUG
       printf("n = %ld, m = %ld", n, m);
#endif
       fmpz_block_init2(entry, m, n);
       if (n == 1)
       {
          fmpz_t * ptr = entry;
          for (ulong i = 0; i < m; i++, ptr++)
          {
             temp_arr[i] = randbits(MPIR_BITS-2);
             ptr->_mp_d = (mp_limb_t *) temp_arr[i];
          }
       } else
       {
          mpn_random2(entry->_mp_d, n*m);
          temp = (mp_limb_t *) mpir_alloc(m*n*sizeof(mp_limb_t));
          F_mpn_copy(temp, entry->_mp_d, m*n);
       }
       ulong newn = n+randint(100);
#if DEBUG
       printf(", newn = %ld\n", newn);
#endif
       fmpz_block_realloc(entry, newn);
       fmpz_block_realloc(entry, n);
       if ((n == 1) && (newn == 1))
       {
          fmpz_t * ptr = entry;
          for (ulong i = 0; i < m; i++, ptr++)
          {
             result &= (temp_arr[i] == (long) ptr->_mp_d);
          }        
       } else if (n == 1)
       {
          fmpz_t * ptr = entry;
          for (ulong i = 0; i < m; i++, ptr++)
          {
             if ((long) ptr->_mp_size < 0L) result &= (temp_arr[i] == -ptr->_mp_d[0]);
             else if ((long) ptr->_mp_size > 0L) result &= (temp_arr[i] == ptr->_mp_d[0]);
             else result &= (temp_arr[i] == 0L);
#if DEBUG
             if (!result) printf("%ld, %ld\n", temp_arr[i], ptr->_mp_d[0]);
#endif
          }                  
       } else result &= (mpn_cmp(entry->_mp_d, temp, m*n) == 0);
       fmpz_block_clear(entry);
       if (n != 1) mpir_free(temp);
   }
   
   for (ulong i = 0; (i < 10000) && (result == 1); i++)
   {
       ulong n = randint(100)+1;
       ulong m = randint(MPIR_BLOCK)+1;
#if DEBUG
       printf("n = %ld, m = %ld", n, m);
#endif
       fmpz_block_init(entry, m);
       fmpz_block_realloc(entry, n);
       if (n == 1)
       {
          fmpz_t * ptr = entry;
          for (ulong i = 0; i < m; i++, ptr++)
          {
             temp_arr[i] = randbits(MPIR_BITS-2);
             ptr->_mp_d = (mp_limb_t *) temp_arr[i];
          }
       } else
       {
          mpn_random2(entry->_mp_d, n*m);
          temp = (mp_limb_t *) mpir_alloc(m*n*sizeof(mp_limb_t));
          F_mpn_copy(temp, entry->_mp_d, m*n);
       }
       ulong newn = n+randint(100);
#if DEBUG
       printf(", newn = %ld\n", newn);
#endif
       fmpz_block_realloc(entry, newn);
       fmpz_block_realloc(entry, n);
       if ((n == 1) && (newn == 1))
       {
          fmpz_t * ptr = entry;
          for (ulong i = 0; i < m; i++, ptr++)
          {
             result &= (temp_arr[i] == (long) ptr->_mp_d);
          }        
       } else if (n == 1)
       {
          fmpz_t * ptr = entry;
          for (ulong i = 0; i < m; i++, ptr++)
          {
             if ((long) ptr->_mp_size < 0L) result &= (temp_arr[i] == -ptr->_mp_d[0]);
             else if ((long) ptr->_mp_size > 0L) result &= (temp_arr[i] == ptr->_mp_d[0]);
             else result &= (temp_arr[i] == 0L);
#if DEBUG
             if (!result) printf("%ld, %ld\n", temp_arr[i], ptr->_mp_d[0]);
#endif
          }                  
       } else result &= (mpn_cmp(entry->_mp_d, temp, m*n) == 0);
       fmpz_block_clear(entry);
       if (n != 1) mpir_free(temp);
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

int test_mpz_to_fmpz()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(1000);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
       
          fmpz_to_mpz(num2, fnum_arr + coeff);
          
          result &= (mpz_cmp(num1, num2) == 0);
#if DEBUG
       if (!result) gmp_printf("%Zd, %Zd\n", num1, num2);
#endif
       }
              
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_realloc_array()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 3000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 100) && (result == 1); j++)
       {
          bits = randint(1000);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
          
          ulong new_length = length + randint(30);
#if DEBUG
          printf("%ld, %ld\n", length, new_length);
#endif
          fnum_arr = fmpz_realloc_array(fnum_arr, length, new_length);
          length = new_length;
          
          fmpz_to_mpz(num2, fnum_arr + coeff);
          
          result &= (mpz_cmp(num1, num2) == 0);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}


int test_fmpz_size()
{
   mpz_t num1;
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(1000);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
       
          result &= (fmpz_size(fnum_arr + coeff) == mpz_size(num1));
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   
   return result;
}

int test_fmpz_get_set_ui()
{
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(MPIR_BITS)+1;
          ulong c = randbits(bits);
          ulong coeff = randint(length);
          fmpz_set_ui(fnum_arr + coeff, c);
       
          ulong c1 = fmpz_get_ui(fnum_arr + coeff);
          
          result &= (c == c1);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   return result;
}

int test_fmpz_get_set_si()
{
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(MPIR_BITS)+1;
          long c = randbits(bits);
          if (randint(2)) c = -c;
          ulong coeff = randint(length);
          fmpz_set_si(fnum_arr + coeff, c);
       
          long c1 = fmpz_get_si(fnum_arr + coeff);
          
          result &= (c == c1);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   return result;
}

int test_fmpz_neg()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(200);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
          fmpz_neg(fnum_arr + coeff, fnum_arr + coeff);
          fmpz_to_mpz(num2, fnum_arr + coeff);
          mpz_neg(num2, num2);          
          result &= (mpz_cmp(num1, num2) == 0);
#if DEBUG2
          if (!result) gmp_printf("%Zd, %Zd\n", num1, num2);
#endif      
       }
       
       fmpz_clear_array(fnum_arr, length);
   }
      
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_set_equal()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr, * fnum;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(200);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          fnum = fmpz_init();
          mpz_to_fmpz(fnum, num1);
          ulong coeff = randint(length);
          fmpz_set(fnum_arr + coeff, fnum);
          fmpz_to_mpz(num2, fnum_arr + coeff);
          
          result &= (mpz_cmp(num1, num2) == 0);
          fmpz_clear(fnum);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
      
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test__fmpz_add_IMM()
{
   mpz_t num1, num2, num3;
   fmpz_t * fnum_arr, * fnum1;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(200);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
          long c = randbits(MPIR_BITS - 2);
          if (randint(2)) c = -c;
          _fmpz_add_IMM(fnum_arr + coeff, fnum_arr + coeff, c);
          fmpz_to_mpz(num2, fnum_arr + coeff);
          if (c < 0L) mpz_sub_ui(num3, num1, -c);
          else mpz_add_ui(num3, num1, c);
          
          result &= (mpz_cmp(num3, num2) == 0);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
      
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(200);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          
          fnum1 = fmpz_init();
          mpz_to_fmpz(fnum1, num1);
          long c = randbits(MPIR_BITS - 2);
          if (randint(2)) c = -c;
          _fmpz_add_IMM(fnum_arr + coeff, fnum1, c);
          fmpz_to_mpz(num2, fnum_arr + coeff);
          if (c < 0L) mpz_sub_ui(num3, num1, -c);
          else mpz_add_ui(num3, num1, c);
          
          result &= (mpz_cmp(num3, num2) == 0);
          fmpz_clear(fnum1);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   
   return result;
}

int test__fmpz_sub_IMM()
{
   mpz_t num1, num2, num3;
   fmpz_t * fnum_arr, * fnum1;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(200);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
          long c = randbits(MPIR_BITS - 2);
          if (randint(2)) c = -c;
          _fmpz_sub_IMM(fnum_arr + coeff, fnum_arr + coeff, c);
          fmpz_to_mpz(num2, fnum_arr + coeff);
          if (c < 0L) mpz_add_ui(num3, num1, -c);
          else mpz_sub_ui(num3, num1, c);
          
          result &= (mpz_cmp(num3, num2) == 0);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
      
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(200);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          
          fnum1 = fmpz_init();
          mpz_to_fmpz(fnum1, num1);
          long c = randbits(MPIR_BITS - 2);
          if (randint(2)) c = -c;
          _fmpz_sub_IMM(fnum_arr + coeff, fnum1, c);
          fmpz_to_mpz(num2, fnum_arr + coeff);
          if (c < 0L) mpz_add_ui(num3, num1, -c);
          else mpz_sub_ui(num3, num1, c);
          
          result &= (mpz_cmp(num3, num2) == 0);
          fmpz_clear(fnum1);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   
   return result;
}

int test_fmpz_add()
{
   mpz_t num1, num2, num3, num4;
   fmpz_t * fnum_arr, * fnum1, * fnum2;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(500);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          bits2 = randint(500);
          mpz_rrandomb(num2, state, bits2);
#if SIGNS
          if (randint(2)) mpz_neg(num2, num2);
#endif

#if DEBUG
       printf("%ld, %ld\n", bits, bits2);
#endif

          ulong coeff = randint(length);
          fnum1 = fmpz_init();
          fnum2 = fmpz_init();
          mpz_to_fmpz(fnum1, num1);
          mpz_to_fmpz(fnum2, num2);
          
          fmpz_add(fnum_arr + coeff, fnum1, fnum2);
          fmpz_to_mpz(num4, fnum_arr + coeff);
          
          mpz_add(num3, num1, num2);
          
          result &= (mpz_cmp(num4, num3) == 0);
          fmpz_clear(fnum1);
          fmpz_clear(fnum2);
          
#if DEBUG
       if (!result) gmp_printf("%Zd, %Zd, %Zd, %Zd\n", num1, num2, num3, num4);
#endif
       }
       
       fmpz_clear_array(fnum_arr, length);
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
   fmpz_t * fnum_arr, * fnum1, * fnum2;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   mpz_init(num4);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(500);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          bits2 = randint(500);
          mpz_rrandomb(num2, state, bits2);
#if SIGNS
          if (randint(2)) mpz_neg(num2, num2);
#endif

#if DEBUG
       printf("%ld, %ld\n", bits, bits2);
#endif

          ulong coeff = randint(length);
          fnum1 = fmpz_init();
          fnum2 = fmpz_init();
          mpz_to_fmpz(fnum1, num1);
          mpz_to_fmpz(fnum2, num2);
          
          fmpz_sub(fnum_arr + coeff, fnum1, fnum2);
          fmpz_to_mpz(num4, fnum_arr + coeff);
          
          mpz_sub(num3, num1, num2);
          
          result &= (mpz_cmp(num4, num3) == 0);
          fmpz_clear(fnum1);
          fmpz_clear(fnum2);
          
#if DEBUG
       if (!result) gmp_printf("%Zd, %Zd, %Zd, %Zd\n", num1, num2, num3, num4);
#endif
       }
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   mpz_clear(num4);
   
   return result;
}

int test_fmpz_addmul_ui()
{
   mpz_t num1, num2, num3;
   fmpz_t * fnum_arr, * fnum;
   unsigned long bits, bits2, bits3, c;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
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

#if DEBUG
          printf("bits = %ld, bits2 = %ld\n", bits, bits2);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
          
          fnum = fmpz_init();
          mpz_to_fmpz(fnum, num2);
          
          bits3 = randint(MPIR_BITS)+1;
          c = randbits(bits3);
          
          fmpz_addmul_ui(fnum_arr + coeff, fnum, c);
          mpz_addmul_ui(num1, num2, c);
          fmpz_to_mpz(num3, fnum_arr + coeff);
          
          result &= (mpz_cmp(num1, num3) == 0);
          
          fmpz_clear(fnum);
       }
       
#if DEBUG
       if (!result) gmp_printf("%ld, %Zd, %Zd\n", c, num1, num3);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   
   return result;
}

int test_fmpz_submul_ui()
{
   mpz_t num1, num2, num3;
   fmpz_t * fnum_arr, * fnum;
   unsigned long bits, bits2, bits3, c;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   mpz_init(num3);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
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

#if DEBUG
          printf("bits = %ld, bits2 = %ld\n", bits, bits2);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
          
          fnum = fmpz_init();
          mpz_to_fmpz(fnum, num2);
          
          bits3 = randint(MPIR_BITS)+1;
          c = randbits(bits3);
          
          fmpz_submul_ui(fnum_arr + coeff, fnum, c);
          mpz_submul_ui(num1, num2, c);
          fmpz_to_mpz(num3, fnum_arr + coeff);
          
          result &= (mpz_cmp(num1, num3) == 0);
          
          fmpz_clear(fnum);
       }
       
#if DEBUG
       if (!result) gmp_printf("%ld, %Zd, %Zd\n", c, num1, num3);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   mpz_clear(num3);
   
   return result;
}

int test_fmpz_mul_2exp()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr, * fnum;
   unsigned long bits, exp;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(1000);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          
          exp = randint(300);
          
          fnum = fmpz_init();
          mpz_to_fmpz(fnum, num1);
          fmpz_mul_2exp(fnum_arr + coeff, fnum, exp);
          mpz_mul_2exp(num2, num1, exp);
          fmpz_to_mpz(num1, fnum_arr + coeff);
          
          result &= (mpz_cmp(num2, num1) == 0);
          
          fmpz_clear(fnum);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_random_bits()
{
   mpz_t num1, num2;
   fmpz_t * fnum1, * fnum2;
   ulong bits, exp;
   int result = 1;
   
   mpz_init(num1);
      
   for (unsigned long i = 0; (i < 100000) && (result == 1); i++)
   {
       bits = randint(1000)+1;
       
#if DEBUG
       printf("Bits = %ld\n", bits);
#endif
       
       fnum1 = fmpz_init();

       fmpz_random(fnum1, bits);
       fmpz_to_mpz(num1, fnum1);
       
       ulong bits1 = fmpz_bits(fnum1);
       ulong bits2 = mpz_sizeinbase(num1, 2);  
         
       result = (((bits1 == bits2) || ((bits1 == 0) && (bits2 == 1))) && (bits1 <= bits));
#if DEBUG2
       if (!result) printf("Bits = %ld, bits1 = %ld, bits2 = %ld\n", bits, bits1, bits2);
#endif
       
       fmpz_clear(fnum1);
   }
   
   mpz_clear(num1);
   
   return result;
}

int test_fmpz_get_d()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(1000);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          mpz_to_fmpz(fnum_arr + coeff, num1);
       
          double d1 = fmpz_get_d(fnum_arr + coeff);
          double d2 = mpz_get_d(num1);
          
          result &= (d1 == d2);
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

int test_fmpz_get_d_2exp()
{
   mpz_t num1, num2;
   fmpz_t * fnum_arr;
   unsigned long bits, bits2;
   int result = 1;
   
   mpz_init(num1);
   mpz_init(num2);
   
   for (unsigned long i = 0; (i < 1000) && (result == 1); i++)
   {
       ulong length = randint(100)+1;
       fnum_arr = fmpz_init_array(length);
       
       for (ulong j = 0; (j < 1000) && (result == 1); j++)
       {
          bits = randint(1000);
          mpz_rrandomb(num1, state, bits);
#if SIGNS
          if (randint(2)) mpz_neg(num1, num1);
#endif
          ulong coeff = randint(length);
          fmpz_set(fnum_arr + coeff, num1);
       
          ulong exp1, exp2;
          double d1 = fmpz_get_d_2exp(&exp1, fnum_arr + coeff);
          double d2 = mpz_get_d_2exp(&exp2, num1);
          
          result &= ((d1 == d2) && (exp1 == exp2));
       }
       
#if DEBUG
       if (!result) gmp_printf("%Zd\n", num1);
#endif
       
       fmpz_clear_array(fnum_arr, length);
   }
   
   mpz_clear(num1);
   mpz_clear(num2);
   
   return result;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(fmpz_block_init_realloc_clear);
   RUN_TEST(fmpz_init_fit_limbs_clear_array);
   RUN_TEST(mpz_to_fmpz);
   RUN_TEST(fmpz_set_equal);
   RUN_TEST(fmpz_realloc_array);
   RUN_TEST(fmpz_size);
   RUN_TEST(fmpz_get_set_ui);
   RUN_TEST(fmpz_get_set_si);
   RUN_TEST(fmpz_neg);
   RUN_TEST(_fmpz_add_IMM);
   RUN_TEST(_fmpz_sub_IMM);
   RUN_TEST(fmpz_add);
   RUN_TEST(fmpz_sub);
   RUN_TEST(fmpz_addmul_ui);
   RUN_TEST(fmpz_submul_ui);
   RUN_TEST(fmpz_mul_2exp);
   RUN_TEST(fmpz_random_bits);
   RUN_TEST(fmpz_get_d);
   RUN_TEST(fmpz_get_d_2exp);
   
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
