/******************************************************************************

 tinyQS.c
 
 Implementation of a tiny hypercube MPQS

 (C) 2006 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <string.h>

#include "../fmpz.h"
#include "../long_extras.h"
#include "../memory-manager.h"

#include "tinyQS.h"
#include "factor_base.h"
#include "poly.h"

/*===========================================================================
   Collect relations:

   Function: Sets up batches of polynomials
             Do the sieving
             Evaluate candidates
             
   Returns:  The number of relations found with this batch of polynomials

===========================================================================*/

unsigned long collect_relations(QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long s = poly_inf->s;
   unsigned long * poly_corr;
   unsigned long relations = 0;
   unsigned long ** A_inv2B = poly_inf->A_inv2B;
   unsigned long i, j;
   
   compute_A(qs_inf, poly_inf);
   compute_B_terms(qs_inf, poly_inf);
   compute_off_adj(qs_inf, poly_inf);
   
   for (i = 0; i < (1<<(s-1)); i++)
   {
      for (j = 0; j < s; j++)
      {
         if (((i>>j) & 1) != 0) break;
      }
      poly_corr = A_inv2B[j];
           
      compute_A_factor_offsets(qs_inf, poly_inf);
      
      relations += 3;    
   }
   
   return relations;
}

/*===========================================================================
   Main Quadratic Sieve Factoring Routine:

   Function: Finds the factors of a number using the quadratic sieve
             Assume n is odd, not a perfect power and not a prime 
             Returns 0 if factorisation was unsuccessful
             Returns 1 if factorisation was successful
             If a small factor is found, it is returned and the QS is not run

===========================================================================*/

int F_mpz_factor_tinyQS(F_mpz_factor_t factors, mpz_t N)
{
   unsigned long small_factor;
   unsigned long rels_found = 0;
   
   mpz_t temp;
   QS_t qs_inf; 
   poly_t poly_inf;
   
   qs_inf.bits = mpz_sizeinbase(N,2);
   if (qs_inf.bits > MAXBITS) return 0; // Number too big for tinyQS 
   
   qs_inf.n = (fmpz_t) flint_stack_alloc(3);
   qs_inf.n[2] = 0;
   mpz_to_fmpz(qs_inf.n, N); // set n to the number to be factored
   
   small_factor = knuth_schroeppel(&qs_inf); // Compute multiplier and some FB primes
   if (small_factor) goto cleanup_2;
   
   primes_init(&qs_inf);
   sqrts_init(&qs_inf);
   
   mpz_init(temp);
   mpz_set(temp, N);
   mpz_mul_ui(temp, temp, qs_inf.k);
   qs_inf.bits = mpz_sizeinbase(temp,2);
   mpz_to_fmpz(qs_inf.n, temp); // set n to the number to be factored times k
   mpz_clear(temp);
      
   if (qs_inf.bits > MAXBITS) 
   {
      small_factor = 0; // kn too big for tinyQS
      goto cleanup_1;
   }
   
   small_factor = compute_factor_base(&qs_inf); // Computes the factor base primes and modular square roots
   if (small_factor) goto cleanup_1;
   compute_sizes(&qs_inf);
   
   poly_init(&qs_inf, &poly_inf, N);
   
   while (rels_found < qs_inf.num_primes + EXTRA_RELS)
   {
      rels_found += collect_relations(&qs_inf, &poly_inf);
   }
   
   small_factor = 1; // sieve was successful
   poly_clear();
   sizes_clear();
cleanup_1:
   sqrts_clear(); // release modular square roots
   primes_clear(); // release factor_base
cleanup_2:
   flint_stack_release(); // release n
   
   return small_factor;    
}

/*===========================================================================
   Main Program:

   Function: Factors a user specified number using the quadratic sieve


===========================================================================*/
/*int main(int argc, unsigned char *argv[])
{
    mpz_t N;
    mpz_init(N); 
    
    F_mpz_factor_t factors;
    
    printf("Input number to factor [ <= 35 decimal digits ] : "); 
    gmp_scanf("%Zd", N); getchar();
    
    F_mpz_factor_tinyQS(factors, N);
    
    mpz_clear(N);
}*/

int main(int argc, unsigned char *argv[])
{
    mpz_t N;
    mpz_init(N); 
    
    F_mpz_factor_t factors;
    unsigned long factor;
    unsigned long failed = 0;
    unsigned long small_factors = 0;
    unsigned long succeed = 0;
    unsigned long bits1, bits2, i;
    
    for (i = 0; i < 1000; i++)
    {
       mpz_set_ui(N, long_nextprime(long_randint(100000000000000000UL)+1UL));
       mpz_mul_ui(N, N, long_nextprime(long_randint(40000000000000000UL)+1UL));
       //bits1 = long_randint(41UL)+13UL;
       //bits2 = long_randint(22UL)+13UL;
       //mpz_mul_ui(N, N, long_nextprime(long_randint((1UL<<bits1)-1UL)+1UL));
       //mpz_mul_ui(N, N, long_nextprime(long_randint((1UL<<bits2)-1UL)+1UL));

#if QS_INFO
       gmp_printf("Factoring %Zd\n", N);
#endif

       factor = F_mpz_factor_tinyQS(factors, N);
       if (!factor) failed++;
       if (factor > 1) small_factors++;
       if (factor == 1) succeed++; 
    }
    
    printf("TinyQS succeeded %ld times, found a small factor %ld times\n", succeed, small_factors);
    printf("and failed %ld times\n", failed);
    
    mpz_clear(N);
}
