/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

kara-profile.c

Comparative profiling for karatsuba multiplication

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler.h"
#include "flint.h"
#include "test-support.h"
#include "mpz_poly.h"
#include "fmpz_poly.h"
#include <stdio.h>


typedef struct
{
   unsigned long length1, length2, bits1, bits2;
   int which;  // 0 for mpz_poly, 1 for fmpz_poly
} arg_t;


void target(void* y, unsigned long count)
{
   arg_t* arg = (arg_t*) y;

   mpz_t x;
   mpz_init(x);
   
   mpz_poly_t in1, in2, out;
   mpz_poly_init(in1);
   mpz_poly_init(in2);
   mpz_poly_init(out);

   fmpz_poly_t in1f, in2f, outf;
   fmpz_poly_init2(in1f, arg->length1, (arg->bits1-1)/FLINT_BITS+1);
   fmpz_poly_init2(in2f, arg->length2, (arg->bits2-1)/FLINT_BITS+1);
   _fmpz_poly_stack_init(outf, arg->length1 + arg->length2 - 1,
                         in1f->limbs + in2f->limbs + 1);

   for (unsigned long i = 0; i < arg->length1; i++)
   {
      mpz_urandomb(x, randstate, arg->bits1);
      mpz_poly_set_coeff(in1, i, x);
   }
   mpz_poly_to_fmpz_poly(in1f, in1);

   for (unsigned long i = 0; i < arg->length2; i++)
   {
      mpz_urandomb(x, randstate, arg->bits2);
      mpz_poly_set_coeff(in2, i, x);
   }
   mpz_poly_to_fmpz_poly(in2f, in2);

   start_clock(0);
   
   if (arg->which)
   {
      for (unsigned long i = 0; i < count; i++)
         _fmpz_poly_mul_karatsuba(outf, in1f, in2f);
   }
   else
   {
      for (unsigned long i = 0; i < count; i++)
         mpz_poly_mul_karatsuba(out, in1, in2);
   }
   
   stop_clock(0);
   
   _fmpz_poly_stack_clear(outf);
   fmpz_poly_clear(in2f);
   fmpz_poly_clear(in1f);

   mpz_poly_clear(out);
   mpz_poly_clear(in2);
   mpz_poly_clear(in1);
   mpz_clear(x);
}


/*
command line arguments are:

length1, length2, bits1, bits2
*/

int main(int argc, char* argv[])
{
   if (argc != 5)
   {
      printf("expected four arguments (see source)\n");
      return 0;
   }
   
   test_support_init();
   
   arg_t arg;
   arg.length1 = atoi(argv[1]);
   arg.length2 = atoi(argv[2]);
   arg.bits1 = atoi(argv[3]);
   arg.bits2 = atoi(argv[4]);
   
   double min_time, max_time;

   for (unsigned long i = 0; i < 3; i++)
   {
      arg.which = 0;
      prof_repeat(&min_time, &max_time, target, &arg);
      printf(" mpz_poly: min = %.3le, \tmax = %.3le\n", min_time, max_time);
      fflush(stdout);

      arg.which = 1;
      prof_repeat(&min_time, &max_time, target, &arg);
      printf("fmpz_poly: min = %.3le, \tmax = %.3le\n", min_time, max_time);
      fflush(stdout);
   }

   test_support_cleanup();
   
   return 0;
}


// end of file ****************************************************************
