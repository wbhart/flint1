/*
   Demo FLINT program for computing the q-expansion of the delta function.
   
   (C) 2007 David Harvey
*/

#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "Zpoly_mpn.h"


void print_poly(Zpoly_mpn_t x)
{
   for (int i = 0; i < x->length; i++)
   {
      if ((signed long) x->coeffs[i*(x->limbs+1)] > 0)
         printf("+%d ", x->coeffs[i*(x->limbs+1)+1]);
      else if ((signed long) x->coeffs[i*(x->limbs+1)] < 0)
         printf("-%d ", x->coeffs[i*(x->limbs+1)+1]);
      else
         printf("0 ");
   }
   printf("\n");
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
   Zpoly_mpn_t F, Fpow;
   Zpoly_mpn_init(F, 2*n, 1);
   Zpoly_mpn_init(Fpow, 2*n, 3);
   
   F->length = n;
   for (long i = 0; i < n; i++)
      _Zpoly_mpn_set_coeff_ui(F, i, 0);
      
   // set F := whatever it's supposed to be
   for (long i = 0; 1; i++)
   {
      long index = i * (i+1) / 2;
      if (index >= n)
         break;
      _Zpoly_mpn_set_coeff_si(F, index, (i & 1) ? -(2*i+1) : (2*i+1));
   }

   print_poly(F);

   // compute F^8, truncated to length n
   _Zpoly_mpn_mul_KS(Fpow, F, F);
   Fpow->length = n;
   print_poly(Fpow);

   _Zpoly_mpn_mul_KS(Fpow, Fpow, Fpow);
   Fpow->length = n;
   print_poly(Fpow);   // this one is wrong

   _Zpoly_mpn_mul_KS(Fpow, Fpow, Fpow);
   Fpow->length = n;
   
   // print out last coefficient (or at least one word of it)
//   printf("coefficient of q^%d is %d\n", n-1, _Zpoly_mpn_get_coeff_ui(Fpow, n-1));
   
   // clean up
   Zpoly_mpn_clear(Fpow);
   Zpoly_mpn_clear(F);
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
