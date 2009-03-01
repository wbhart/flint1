/*
   test.c:  command line program for running test routines
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <string.h>
#include "support.h"
#include "zn_poly_internal.h"


typedef struct
{
   // command line option string for this target
   char* option;
   
   // function name to print
   char* name;
   
   // pointer to actual test function; should return 1 if test passes
   int (*func)();

} test_target_t;


extern int test_zn_array_mul_KS1();
extern int test_zn_array_mul_KS2();
extern int test_zn_array_mul_KS3();
extern int test_zn_array_mul_KS4();
extern int test_zn_array_sqr_KS1();
extern int test_zn_array_sqr_KS2();
extern int test_zn_array_sqr_KS3();
extern int test_zn_array_sqr_KS4();
extern int test_zn_array_recip_fix_reduce();
extern int test_zn_array_pack();
extern int test_zn_array_unpack();
extern int test_zn_array_mul_fft();
extern int test_zn_array_sqr_fft();
extern int test_zn_array_mul_fft_dft();
extern int test_zn_array_midmul_fft();
extern int test_nussbaumer_mul();
extern int test_zn_pmf_vec_fft_small();
extern int test_zn_pmf_vec_fft_factor();
extern int test_zn_pmf_vec_ifft_small();
extern int test_zn_pmf_vec_ifft_factor();
extern int test_zn_pmf_vec_fft_transposed_small();
extern int test_zn_pmf_vec_fft_transposed_factor();
extern int test_zn_pmf_vec_ifft_transposed_small();
extern int test_zn_pmf_vec_ifft_transposed_factor();
extern int test_zn_array_invert();


test_target_t targets[] = {
   {"mul_KS1",
    "zn_array_mul_KS1",
    test_zn_array_mul_KS1},
    
   {"mul_KS2",
    "zn_array_mul_KS2",
    test_zn_array_mul_KS2},
    
   {"mul_KS3",
    "zn_array_mul_KS3",
    test_zn_array_mul_KS3},
    
   {"mul_KS4",
    "zn_array_mul_KS4",
    test_zn_array_mul_KS4},
    
   {"sqr_KS1",
    "zn_array_sqr_KS1",
    test_zn_array_sqr_KS1},
    
   {"sqr_KS2",
    "zn_array_sqr_KS2",
    test_zn_array_sqr_KS2},
    
   {"sqr_KS3",
    "zn_array_sqr_KS3",
    test_zn_array_sqr_KS3},
    
   {"sqr_KS4",
    "zn_array_sqr_KS4",
    test_zn_array_sqr_KS4},
    
   {"recip",
    "zn_array_recip_fix_reduce",
    test_zn_array_recip_fix_reduce},
    
   {"pack",
    "zn_array_pack",
    test_zn_array_pack},
    
   {"unpack",
    "zn_array_unpack",
    test_zn_array_unpack},
    
   {"mul_fft",
    "zn_array_mul_fft",
    test_zn_array_mul_fft},
    
   {"sqr_fft",
    "zn_array_sqr_fft",
    test_zn_array_sqr_fft},

   {"mul_fft_dft",
    "zn_array_mul_fft_dft",
    test_zn_array_mul_fft_dft},
    
   {"midmul_fft",
    "zn_array_midmul_fft",
    test_zn_array_midmul_fft},
    
   {"nussbaumer",
    "nussbaumer_mul",
    test_nussbaumer_mul},

   {"fft_small",
    "zn_pmf_vec_fft_small",
    test_zn_pmf_vec_fft_small},
    
   {"fft_factor",
    "zn_pmf_vec_fft_factor",
    test_zn_pmf_vec_fft_factor},
    
   {"ifft_small",
    "zn_pmf_vec_ifft_small",
    test_zn_pmf_vec_ifft_small},
    
   {"ifft_factor",
    "zn_pmf_vec_ifft_factor",
    test_zn_pmf_vec_ifft_factor},

   {"fft_transposed_small",
    "zn_pmf_vec_fft_transposed_small",
    test_zn_pmf_vec_fft_transposed_small},
    
   {"fft_transposed_factor",
    "zn_pmf_vec_fft_transposed_factor",
    test_zn_pmf_vec_fft_transposed_factor},
    
   {"ifft_transposed_small",
    "zn_pmf_vec_ifft_transposed_small",
    test_zn_pmf_vec_ifft_transposed_small},

   {"ifft_transposed_factor",
    "zn_pmf_vec_ifft_transposed_factor",
    test_zn_pmf_vec_ifft_transposed_factor},

   {"invert",
    "zn_array_invert",
    test_zn_array_invert},

};

const unsigned num_targets = sizeof(targets) / sizeof(targets[0]);


int run_test(test_target_t* target)
{
   printf("%s()... ", target->name);
   fflush(stdout);
   int success = target->func();
   printf("%s\n", success ? "ok" : "FAIL!");
   return success;
}


int main(int argc, char* argv[])
{
   gmp_randinit_default(randstate);
   
   if (argc == 1)
   {
      printf("usage: test <target1> <target2> ...\n\n");
      printf("Available targets:\n\n");
      printf("   %25s: runs all tests\n", "all");
      int i;
      for (i = 0; i < num_targets; i++)
         printf("   %25s: %s\n", targets[i].option, targets[i].name);
      return 0;
   }

   int all_success = 1, success, i, j;

   for (j = 1; j < argc; j++)
   {
      if (!strcmp(argv[j], "all"))
      {
         for (i = 0; i < num_targets; i++)
            all_success = all_success && run_test(&targets[i]);
      }
      else
      {
         int found = -1;
         for (i = 0; i < num_targets; i++)
            if (!strcmp(argv[j], targets[i].option))
               found = i;
         if (found == -1)
         {
            printf("unknown target string \"%s\"\n", argv[j]);
            return 0;
         }

         all_success = all_success && run_test(&targets[found]);
      }
   }

   printf("\n");
   if (all_success)
      printf("All tests passed.\n");
   else
      printf("At least one test FAILED!\n");
   
   gmp_randclear(randstate);

   return !all_success;
}


// end of file ****************************************************************
