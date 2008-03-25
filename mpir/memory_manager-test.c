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

memory_manager-test.c: test module for memory_manager

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "mpir.h"
#include "memory_manager.h"
#include "test_support.h"

#define DEBUG 0    // prints debug information
#define DEBUG2 1 

int test_mpir_alloc_realloc_free()
{    
    for (ulong count = 0; count < 30000; count++)
    {
        ulong bytes = randint(100) + 1;
        ulong bytes2 = randint(100) + bytes + 1;

        char * ptr = (char *) mpir_alloc(bytes);
        for (ulong i = 0; i < bytes; i++)
           ptr[i] = randbyte();
        ptr = mpir_realloc(ptr, bytes2);
        for (ulong i = 0; i < bytes2; i++)
           ptr[i] = randbyte();
        mpir_free(ptr);
    }
    
    return 1;
}

int test_mpir_stack_alloc_release()
{    
    for (ulong count = 0; count < 3000; count++)
    {
        ulong num = randint(100) + 1;

        for (ulong i = 0; i < num; i++)
        {
           ulong bytes = randint(100) + 1;

           char * ptr = (char *) mpir_stack_alloc(bytes);
           for (ulong i = 0; i < bytes; i++)
              ptr[i] = randbyte();
        }

        for (ulong i = 0; i < num; i++)
        {
           mpir_stack_release();
        }
    }
    
    return 1;
}

int test_mpir_stack_alloc_release_small()
{    
    for (ulong count = 0; count < 3000; count++)
    {
        ulong num = randint(100) + 1;

        for (ulong i = 0; i < num; i++)
        {
           ulong bytes = randint(100) + 1;

           char * ptr = (char *) mpir_stack_alloc_small(bytes);
           for (ulong i = 0; i < bytes; i++)
              ptr[i] = randbyte();
        }

        for (ulong i = 0; i < num; i++)
        {
           mpir_stack_release_small();
        }
    }
    
    return 1;
}

int test_mpir_aligned_alloc_realloc_free()
{    
    for (ulong count = 0; count < 30000; count++)
    {
        ulong bytes = randint(100) + 1;
        ulong bytes2 = randint(100) + bytes + 1;
        
        char * ptr = (char *) mpir_aligned_alloc(bytes);
        for (ulong i = 0; i < bytes; i++)
           ptr[i] = randbyte();
        ptr = mpir_aligned_realloc(ptr, bytes2);
        for (ulong i = 0; i < bytes2; i++)
           ptr[i] = randbyte();
        mpir_aligned_free(ptr);
    }
    
    return 1;
}


/****************************************************************************

   Main test functions

****************************************************************************/

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");


void mpir_test_all()
{
   int success, all_success = 1;

   RUN_TEST(mpir_alloc_realloc_free);
   RUN_TEST(mpir_stack_alloc_release);
   RUN_TEST(mpir_stack_alloc_release_small);
   RUN_TEST(mpir_aligned_alloc_realloc_free);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   mpir_test_all();
   
   mpir_stack_cleanup();

   return 0;
}



// end of file ****************************************************************
