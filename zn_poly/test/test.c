/*
   test.c:  command line program for running test routines
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
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
   // function name to print
   char* name;
   
   // pointer to actual test function; should return 1 if test passes
   // parameter quick is 1 if only a quick test is wanted
   int (*func)(int quick);

}
test_target_t;


extern int test_zn_array_mul_KS1 (int quick);
extern int test_zn_array_mul_KS2 (int quick);
extern int test_zn_array_mul_KS3 (int quick);
extern int test_zn_array_mul_KS4 (int quick);
extern int test_zn_array_sqr_KS1 (int quick);
extern int test_zn_array_sqr_KS2 (int quick);
extern int test_zn_array_sqr_KS3 (int quick);
extern int test_zn_array_sqr_KS4 (int quick);
extern int test_zn_array_mulmid_KS1 (int quick);
extern int test_zn_array_mulmid_KS2 (int quick);
extern int test_zn_array_mulmid_KS3 (int quick);
extern int test_zn_array_mulmid_KS4 (int quick);
extern int test_zn_array_recover_reduce (int quick);
extern int test_zn_array_pack (int quick);
extern int test_zn_array_unpack (int quick);
extern int test_zn_array_mul_fft (int quick);
extern int test_zn_array_sqr_fft (int quick);
extern int test_zn_array_mul_fft_dft (int quick);
extern int test_zn_array_mulmid_fft (int quick);
extern int test_nuss_mul (int quick);
extern int test_pmfvec_fft_dc (int quick);
extern int test_pmfvec_fft_huge (int quick);
extern int test_pmfvec_ifft_dc (int quick);
extern int test_pmfvec_ifft_huge (int quick);
extern int test_pmfvec_tpfft_dc (int quick);
extern int test_pmfvec_tpfft_huge (int quick);
extern int test_pmfvec_tpifft_dc (int quick);
extern int test_pmfvec_tpifft_huge (int quick);
extern int test_zn_array_invert (int quick);
extern int test_mpn_smp_basecase (int quick);
extern int test_mpn_smp_kara (int quick);
extern int test_mpn_smp (int quick);
extern int test_mpn_mulmid (int quick);


test_target_t targets[] = {
   {"mpn_smp_basecase",
    test_mpn_smp_basecase},

   {"mpn_smp_kara",
    test_mpn_smp_kara},

   {"mpn_smp",
    test_mpn_smp},

   {"mpn_mulmid",
    test_mpn_mulmid},

   {"zn_array_recover_reduce",
    test_zn_array_recover_reduce},
    
   {"zn_array_pack",
    test_zn_array_pack},
    
   {"zn_array_unpack",
    test_zn_array_unpack},
    
   {"zn_array_mul_KS1",
    test_zn_array_mul_KS1},
    
   {"zn_array_mul_KS2",
    test_zn_array_mul_KS2},
    
   {"zn_array_mul_KS3",
    test_zn_array_mul_KS3},
    
   {"zn_array_mul_KS4",
    test_zn_array_mul_KS4},
    
   {"zn_array_sqr_KS1",
    test_zn_array_sqr_KS1},
    
   {"zn_array_sqr_KS2",
    test_zn_array_sqr_KS2},
    
   {"zn_array_sqr_KS3",
    test_zn_array_sqr_KS3},
    
   {"zn_array_sqr_KS4",
    test_zn_array_sqr_KS4},
    
   {"zn_array_mulmid_KS1",
    test_zn_array_mulmid_KS1},
    
   {"zn_array_mulmid_KS2",
    test_zn_array_mulmid_KS2},
    
   {"zn_array_mulmid_KS3",
    test_zn_array_mulmid_KS3},
    
   {"zn_array_mulmid_KS4",
    test_zn_array_mulmid_KS4},
    
   {"nuss_mul",
    test_nuss_mul},

   {"pmfvec_fft_dc",
    test_pmfvec_fft_dc},
    
   {"pmfvec_fft_huge",
    test_pmfvec_fft_huge},
    
   {"pmfvec_ifft_dc",
    test_pmfvec_ifft_dc},
    
   {"pmfvec_ifft_huge",
    test_pmfvec_ifft_huge},

   {"pmfvec_tpfft_dc",
    test_pmfvec_tpfft_dc},
    
   {"pmfvec_tpfft_huge",
    test_pmfvec_tpfft_huge},
    
   {"pmfvec_tpifft_dc",
    test_pmfvec_tpifft_dc},

   {"pmfvec_tpifft_huge",
    test_pmfvec_tpifft_huge},

   {"zn_array_mul_fft",
    test_zn_array_mul_fft},
    
   {"zn_array_sqr_fft",
    test_zn_array_sqr_fft},

   {"zn_array_mulmid_fft",
    test_zn_array_mulmid_fft},
    
   {"zn_array_mul_fft_dft",
    test_zn_array_mul_fft_dft},
    
   {"zn_array_invert",
    test_zn_array_invert},

};

const unsigned num_targets = sizeof (targets) / sizeof (targets[0]);


int
run_test (test_target_t* target, int quick)
{
   printf ("%s()... ", target->name);
   fflush (stdout);
   int success = target->func (quick);
   printf ("%s\n", success ? "ok" : "FAIL!");
   return success;
}


void
usage ()
{
   printf ("usage: test [ -quick ] <target1> <target2> ...\n\n");
   printf ("Available targets:\n\n");
   printf ("   all (runs all tests)\n");
   int i;
   for (i = 0; i < num_targets; i++)
      printf ("   %s\n", targets[i].name);
}


int
main (int argc, char* argv[])
{
   gmp_randinit_default (randstate);
   
   int all_success = 1, any_targets = 0, quick = 0, success, i, j;

   for (j = 1; j < argc; j++)
   {
      if (!strcmp (argv[j], "-quick"))
         quick = 1;
      else if (!strcmp (argv[j], "all"))
      {
         any_targets = 1;
         for (i = 0; i < num_targets; i++)
            all_success = all_success && run_test (&targets[i], quick);
      }
      else
      {
         int found = -1;
         for (i = 0; i < num_targets; i++)
            if (!strcmp (argv[j], targets[i].name))
               found = i;
         if (found == -1)
         {
            printf ("unknown target string \"%s\"\n", argv[j]);
            return 0;
         }

         any_targets = 1;
         all_success = all_success && run_test (&targets[found], quick);
      }
   }
   
   if (!any_targets)
   {
      usage ();
      return 0;
   }

   printf("\n");
   if (all_success)
      printf ("All tests passed.\n");
   else
      printf ("At least one test FAILED!\n");
   
   gmp_randclear (randstate);

   return !all_success;
}


// end of file ****************************************************************
