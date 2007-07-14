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
#include "sieve.h"

/*===========================================================================
   Collect relations:

   Function: Sets up batches of polynomials
             Do the sieving
             Evaluate candidates
             
   Returns:  The number of relations found with this batch of polynomials

===========================================================================*/

unsigned long collect_relations(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
{
   unsigned long s = poly_inf->s;
   unsigned long * poly_corr;
   unsigned long relations = 0;
   unsigned long ** A_inv2B = poly_inf->A_inv2B;
   unsigned long poly_index, j;
   unsigned long poly_add;
   
   compute_A(qs_inf, poly_inf);
   compute_B_terms(qs_inf, poly_inf);
   compute_off_adj(qs_inf, poly_inf);
   compute_A_factor_offsets(qs_inf, poly_inf);
   compute_C(qs_inf, poly_inf);          
      
   for (poly_index = 1; poly_index < (1<<(s-1)); poly_index++)
   {
      for (j = 0; j < s; j++)
      {
         if (((poly_index >> j) & 1UL) != 0UL) break;
      }
      
      poly_add = (((poly_index >> j) & 2UL) != 0UL);
      
      poly_corr = A_inv2B[j];
           
      do_sieving(qs_inf, poly_inf, sieve);
      
      relations += evaluate_sieve(qs_inf, poly_inf, sieve);
      
      update_offsets(poly_add, poly_corr, qs_inf, poly_inf);
      
      if (poly_add) poly_inf->B += (2*poly_inf->B_terms[j]); 
      else poly_inf->B -= (2*poly_inf->B_terms[j]); 
      
      compute_C(qs_inf, poly_inf);          
      
      compute_A_factor_offsets(qs_inf, poly_inf);    
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
   
   QS_t qs_inf; 
   poly_t poly_inf;
   
   qs_inf.bits = mpz_sizeinbase(N,2);
   if (qs_inf.bits > MAXBITS) return 0; // Number too big for tinyQS 
   
   qs_inf.n = (fmpz_t) flint_stack_alloc(3);
   qs_inf.n[2] = 0;
   mpz_to_fmpz(qs_inf.n, N); // set n to the number to be factored
   
   small_factor = knuth_schroeppel(&qs_inf); // Compute multiplier and some FB primes
   if (small_factor) goto cleanup_2;

#if QS_INFO
   printf("Multiplier = %ld\n", qs_inf.k);
#endif
   
   mpz_init(qs_inf.mpz_n);
   mpz_set(qs_inf.mpz_n, N);
   mpz_mul_ui(qs_inf.mpz_n, qs_inf.mpz_n, qs_inf.k);
   qs_inf.bits = mpz_sizeinbase(qs_inf.mpz_n,2);
   mpz_to_fmpz(qs_inf.n, qs_inf.mpz_n); // set n to the number to be factored times k
   
   primes_init(&qs_inf);
   sqrts_init(&qs_inf);
      
   if (qs_inf.bits > MAXBITS) 
   {
      small_factor = 0; // kn too big for tinyQS
      goto cleanup_1;
   }
   
   small_factor = compute_factor_base(&qs_inf); // Computes the factor base primes and modular square roots
   if (small_factor) goto cleanup_1;
   compute_sizes(&qs_inf);
   
   poly_init(&qs_inf, &poly_inf, N);
   
   unsigned char * sieve = (unsigned char *) flint_stack_alloc_bytes(SIEVE_SIZE+1);
   while (rels_found < qs_inf.num_primes + EXTRA_RELS)
   {
      rels_found += collect_relations(&qs_inf, &poly_inf, sieve);
   }
   flint_stack_release(); // release sieve
   
   small_factor = 1; // sieve was successful
   poly_clear(&poly_inf);
   sizes_clear();
cleanup_1:
   sqrts_clear(); // release modular square roots
   primes_clear(); // release factor_base
cleanup_2:
   flint_stack_release(); // release n
   mpz_clear(qs_inf.mpz_n);
   
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
    
    for (i = 0; i < 100; i++)
    {
       mpz_set_ui(N, long_nextprime(long_randint(4000000000000000000UL)+1UL));
       mpz_mul_ui(N, N, long_nextprime(long_randint(4000000000000000000UL)+1UL));
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
