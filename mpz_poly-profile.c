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

mpz_poly-profile.c

Profiling for mpz_poly

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "profiler-main.h"
#include "mpz_poly.h"
#include "flint.h"
#include "test-support.h"
#include <string.h>
#include <math.h>


// ============================================================================


/*
this function samples multiplying polynomials of lengths len1 and len2
using mpz_poly_mul_karatsuba

arg should point to an unsigned long, giving the coefficient bitlengths
*/
void sample_mpz_poly_mul_karatsuba_mixlengths(
     unsigned long len1, unsigned long len2, void* arg, unsigned long count)
{
   unsigned long bits = *(unsigned long*) arg;
   
   mpz_poly_t poly1, poly2, poly3;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init(poly3);

   mpz_t x;
   mpz_init(x);
   for (unsigned long i = 0; i < len1; i++)
   {
      mpz_urandomb(x, randstate, bits);
      if (random_ulong(2)) mpz_neg(x, x);
      mpz_poly_set_coeff(poly1, i, x);
   }
   for (unsigned long i = 0; i < len2; i++)
   {
      mpz_urandomb(x, randstate, bits);
      if (random_ulong(2)) mpz_neg(x, x);
      mpz_poly_set_coeff(poly2, i, x);
   }
   mpz_clear(x);
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      mpz_poly_mul_karatsuba(poly3, poly1, poly2);

   prof_stop();
   
   mpz_poly_clear(poly3);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
}


char* profDriverString_mpz_poly_mul_karatsuba_mixlengths(char* params)
{
   return "mpz_poly_mul_karatubsa for distinct input lengths and fixed\n"
   "coefficient size. Parameters are: max length; length skip; coefficient size (in bits)\n";
}

char* profDriverDefaultParams_mpz_poly_mul_karatsuba_mixlengths()
{
   return "50 1 100";
}


void profDriver_mpz_poly_mul_karatsuba_mixlengths(char* params)
{
   unsigned long max_length, skip, bits;

   sscanf(params, "%ld %ld %ld", &max_length, &skip, &bits);

   prof2d_set_sampler(sample_mpz_poly_mul_karatsuba_mixlengths);

   test_support_init();

   for (unsigned long len1 = skip; len1 <= max_length; len1 += skip)
      for (unsigned long len2 = skip; len2 <= len1; len2 += skip)
         prof2d_sample(len1, len2, &bits);

   test_support_cleanup();
}

// ============================================================================


/*
this function samples multiplying polynomials of lengths len1 and len2
using mpz_poly_mul_karatsuba

arg should point to an unsigned long, giving the coefficient bitlengths
*/
void sample__mpz_poly_mul_kara_recursive_mixlengths(
     unsigned long len1, unsigned long len2, void* arg, unsigned long count)
{
   unsigned long bits = *(unsigned long*) arg;
   
   mpz_poly_t poly1, poly2, poly3;
   mpz_poly_init(poly1);
   mpz_poly_init(poly2);
   mpz_poly_init(poly3);
   mpz_poly_ensure_alloc(poly3, len1 + len2);

   // allocate scratch space
   unsigned long scratch_len = len1 + len2;
   mpz_t* scratch = (mpz_t*) malloc(scratch_len * sizeof(mpz_t));
   for (unsigned long i = 0; i < scratch_len; i++)
      mpz_init2(scratch[i], 2*bits + 100);

   mpz_t x;
   mpz_init(x);
   for (unsigned long i = 0; i < len1; i++)
   {
      mpz_urandomb(x, randstate, bits);
      if (random_ulong(2)) mpz_neg(x, x);
      mpz_poly_set_coeff(poly1, i, x);
   }
   for (unsigned long i = 0; i < len2; i++)
   {
      mpz_urandomb(x, randstate, bits);
      if (random_ulong(2)) mpz_neg(x, x);
      mpz_poly_set_coeff(poly2, i, x);
   }
   mpz_clear(x);

   unsigned long crossover =
                      _mpz_poly_mul_karatsuba_crossover(bits / FLINT_BITS);
   
   prof_start();

   for (unsigned long i = 0; i < count; i++)
      _mpz_poly_mul_kara_recursive(poly3->coeffs, poly1->coeffs, len1,
                                   poly2->coeffs, len2, scratch, 1, crossover);

   prof_stop();

   for (unsigned long i = 0; i < scratch_len; i++)
      mpz_clear(scratch[i]);
   free(scratch);
   
   mpz_poly_clear(poly3);
   mpz_poly_clear(poly2);
   mpz_poly_clear(poly1);
}


char* profDriverString__mpz_poly_mul_kara_recursive_mixlengths(char* params)
{
   return "_mpz_poly_mul_kara_recursive for distinct input lengths and fixed\n"
   "coefficient size. Parameters are: max length; length skip; coefficient size (in bits)\n";
}

char* profDriverDefaultParams__mpz_poly_mul_kara_recursive_mixlengths()
{
   return "50 1 100";
}


void profDriver__mpz_poly_mul_kara_recursive_mixlengths(char* params)
{
   unsigned long max_length, skip, bits;

   sscanf(params, "%ld %ld %ld", &max_length, &skip, &bits);

   prof2d_set_sampler(sample__mpz_poly_mul_kara_recursive_mixlengths);

   test_support_init();

   for (unsigned long len2 = skip; len2 <= max_length; len2 += skip)
      for (unsigned long len1 = skip; len1 <= len2; len1 += skip)
         prof2d_sample(len1, len2, &bits);

   test_support_cleanup();
}


// end of file ****************************************************************
