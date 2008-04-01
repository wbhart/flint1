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

fmpz_mat-test.c: Test code for fmpz_mat.c and fmpz_mat.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "mpir.h"
#include "test_support.h"
#include "fmpz_mat.h"

#define VARY_BITS 0
#define SIGNS 1
#define SPARSE 0

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

gmp_randstate_t state;

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
int test_fmpz_mat_init_clear()
{
   fmpz_mat_t mat;
   int result = 1;
   ulong r, c;
   
   for (ulong count1 = 1; (count1 < 10000) && (result == 1) ; count1++)
   {
      r = randint(100) + 1;
      c = randint(100) + 1;
      
      fmpz_mat_init(mat, r, c);
      fmpz_mat_clear(mat);
   }
         
   return result;
}

int test_fmpz_mat_get_set_entry_ui()
{
   fmpz_mat_t mat;
   int result = 1;
   ulong r, c, r1, c1, val1, val2, bits;
   
   for (ulong count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      r = randint(100) + 1;
      c = randint(100) + 1;
      
      fmpz_mat_init(mat, r, c);
      for (ulong i = 0; i < 10000; i++)
      {
         r1 = randint(r);
         c1 = randint(c);
         bits = randint(MPIR_BITS) + 1;
         val1 = randbits(bits);
         fmpz_mat_set_entry_ui(mat, r1, c1, val1);
         val2 = fmpz_mat_get_entry_ui(mat, r1, c1);
         result &= (val1 == val2);
      }
      fmpz_mat_clear(mat);
   }
         
   return result;
}

void fmpz_mat_test_all()
{
   int success, all_success = 1;

   RUN_TEST(fmpz_mat_init_clear);
   RUN_TEST(fmpz_mat_get_set_entry_ui);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(state);
   fmpz_mat_test_all();
   gmp_randclear(state);
   
   mpir_stack_cleanup();

   return 0;
}
