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

   unsigned long limbs = (FLINT_BITS_PER_LIMB == 32) ? 6 : 3;

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
   fmpz_poly_init(F2, 2*N, limbs);
   fmpz_poly_init(F4, 2*N, limbs);
   fmpz_poly_init(F8, 2*N, limbs);

   // Initialise F2 with coefficients of F(q)^2
   fmpz_poly_ensure_space(F2, N);
   for (long i = 0; i < N; i++)
      fmpz_poly_set_coeff_si(F2, i, values[i]);

   free(values);
   printf("Done F2\n");
   // compute F^4, truncated to length N
   _fmpz_poly_mul_KS(F4, F2, F2);
   fmpz_poly_set_length(F4, N);
   printf("Done F4\n");
   
   // compute F^8, truncated to length N
   _fmpz_poly_mul_KS(F8, F4, F4);
   fmpz_poly_set_length(F8, N);
   printf("Done F8\n");
   
   // print out last coefficient
   mpz_t x;
   mpz_init(x);
   fmpz_poly_get_coeff_mpz(x, F8, N-1);
   gmp_printf("coefficient of q^%d is %Zd\n", N, x);
   mpz_clear(x);
   
   // clean up
   fmpz_poly_clear(F8);
   fmpz_poly_clear(F4);
   fmpz_poly_clear(F2);

   return 0;
}
