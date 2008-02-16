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
/*
   Demo FLINT program for computing the q-expansion of the delta function.
   
   (C) 2007 David Harvey and William Hart
*/

#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "fmpz_poly.h"


int main(int argc, char* argv[])
{
   if (argc != 2)
   {
      printf("Syntax: delta_qexp <integer>\n");
      printf("where <integer> is the number of terms to compute\n");
      return 0;
   }

   // number of terms to compute
   long N = atoi(argv[1]);

   // compute coefficients of F(q)^2
   long* values = malloc(sizeof(long) * N);
   for (long i = 0; i < N; i++)
      values[i] = 0;

   long stop = (long) ceil((-1.0 + sqrt(1.0 + 8.0*N)) / 2.0);
   
   for (long i = 0; i <= stop; i++)
   {
      long index1 = i*(i+1)/2;
      long value1 = (i & 1) ? (-2*i-1) : (2*i+1);
      for (long j = 0; j <= stop; j++)
      {
         long index2 = j*(j+1)/2;
         if (index1 + index2 >= N)
            break;
         long value2 = (j & 1) ? (-2*j-1) : (2*j+1);
         values[index1 + index2] += value1 * value2;
      }
   }

   // Create some polynomial objects
   fmpz_poly_t F2, F4, F8;
   fmpz_poly_init(F2);
   fmpz_poly_init(F4);
   fmpz_poly_init(F8);
   
   fmpz_poly_fit_length(F2, N);
   
   for (long i = 0; i < N; i++)
      fmpz_poly_set_coeff_si(F2, i, values[i]);

   free(values);
   
   // compute F^4, truncated to length N
   fmpz_poly_mul_trunc_n(F4, F2, F2, N);
   
   // compute F^8, truncated to length N
   fmpz_poly_mul_trunc_n(F8, F4, F4, N);
   
   // print out last coefficient
   fmpz_t coeff = fmpz_poly_get_coeff_ptr(F8, N-1);
   printf("coefficient of q^%d is ", N);
   fmpz_print(coeff); 
   printf("\n");
   
   // clean up
   fmpz_poly_clear(F8);
   fmpz_poly_clear(F4);
   fmpz_poly_clear(F2);

   return 0;
}
