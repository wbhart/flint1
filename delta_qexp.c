/*
   Demo FLINT program for computing the q-expansion of the delta function.
   
   (C) 2007 David Harvey and William Hart
*/

#include <stdio.h>
#include <gmp.h>
#include <math.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"


void print_poly(Zpoly_mpn_t x_mpn)
{
   Zpoly_t x;
   Zpoly_init2(x, x_mpn->length);
   _Zpoly_mpn_convert_out(x, x_mpn);
   Zpoly_print(stdout, x);
   Zpoly_clear(x);
}


int main(int argc, char* argv[])
{
   if (argc != 2)
   {
      printf("Syntax: delta_qexp <integer>\n");
      printf("where <integer> is the number of terms to compute\n");
      return 0;
   }
   
   // number of terms to compute
   long n = atoi(argv[1]);

   // initialise polynomial objects
   // (play it safe with allocated space for the moment)
   Zpoly_mpn_t F2, F4, F8;
   unsigned long limbs = (FLINT_BITS_PER_LIMB == 32) ? 6 : 3;
   Zpoly_mpn_init(F2, 2*n, limbs);
   Zpoly_mpn_init(F4, 2*n, limbs);
   Zpoly_mpn_init(F8, 2*n, limbs);
   
   // set F2 := whatever it's supposed to be
   long* values = malloc(sizeof(long) * n);
   for (long i = 0; i < n; i++)
      values[i] = 0;

   long stop = (long) ceil((-1.0 + sqrt(1.0 + 8.0*n))/2.0);
   for (long i = 0; i <= stop; i++)
   {
      long index1 = i*(i+1)/2;
      long value1 = (i & 1) ? (-2*i-1) : (2*i+1);
      for (long j = 0; j <= stop; j++)
      {
         long index2 = j*(j+1)/2;
         if (index1 + index2 >= n)
            break;
         long value2 = (j & 1) ? (-2*j-1) : (2*j+1);
         values[index1 + index2] += value1 * value2;
      }
   }

   // NB: could avoid reallocations in the next loop by calling
   //    Zpoly_mpn_ensure_space(F2, n);
   for (long i = 0; i < n; i++)
      Zpoly_mpn_set_coeff_si(F2, i, values[i]);

   free(values);

   // compute F^4, truncated to length n
   _Zpoly_mpn_mul_KS(F4, F2, F2);
   Zpoly_mpn_set_length(F4, n);

   // compute F^8, truncated to length n
   _Zpoly_mpn_mul_KS(F8, F4, F4);
   Zpoly_mpn_set_length(F8, n);
   
   // print out last coefficient
   Zpoly_t output;
   Zpoly_init2(output, n);
   _Zpoly_mpn_convert_out(output, F8);
   gmp_printf("coefficient of q^%d is %Zd\n", n, output->coeffs[n-1]);
   Zpoly_clear(output);
   
   // clean up
   Zpoly_mpn_clear(F8);
   Zpoly_mpn_clear(F4);
   Zpoly_mpn_clear(F2);

   return 0;
}


// Here's what it *should* look like:
/*
int main()
{
   if (argc != 2)
   {
      printf("Syntax: delta_qexp <integer>\n");
      printf("where <integer> is the number of terms to compute\n");
      return 0;
   }
   
   // number of terms to compute
   long n = atoi(argv[1]);

   // initialise polynomial objects
   Zpoly_mpn_t F, Fpow;
   Zpoly_mpn_init(F);
   Zpoly_mpn_init(Fpow);

   // set F := whatever it's supposed to be
   for (long i = 0; 1; i++)
   {
      long index = i * (i+1) / 2;
      if (index >= n)
         break;
      Zpoly_mpn_set_coeff_si(F, index, (i & 1) ? -(2*i+1) : (2*i+1));
   }

   // compute F^8, truncated to length n
   Zpoly_mpn_mul(Fpow, F, F);
   Zpoly_mpn_truncate(Fpow, n);

   Zpoly_mpn_mul(Fpow, Fpow, Fpow);
   Zpoly_mpn_truncate(Fpow, n);

   Zpoly_mpn_mul(Fpow, Fpow, Fpow);
   Zpoly_mpn_truncate(Fpow, n);
   
   // print out last coefficient
   mpz_t x;
   mpz_init(x);
   Zpoly_mpn_get_coeff_mpz(Fpow, x, n-1);
   
   gmp_printf("coefficient of q^%d is %Z\n", n-1, x);
   
   // clean up
   mpz_clear(x);
   Zpoly_mpn_clear(Fpow);
   Zpoly_mpn_clear(F);
   return 0;
}
*/
