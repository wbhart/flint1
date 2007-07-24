/******************************************************************************

 mpQS.c
 
 A full implementation of the Quadratic Sieve using multi-precision arithmetic

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
#include "../flint.h"
#include "../mpn_extras.h"

#include "mpQS.h"
#include "mp_factor_base.h"
#include "mp_poly.h"
#include "mp_sieve.h"
#include "mp_linear_algebra.h"
#include "block_lanczos.h"

/*===========================================================================
   Square Root:

   Function: Compute the square root of the product of all the partial 
             relations and take it mod p

===========================================================================*/

static inline void square_root(mpz_t X, mpz_t Y, QS_t * qs_inf, linalg_t * la_inf, 
   u_int64_t * nullrows, unsigned long ncols, unsigned long l, mpz_t N)
{
   unsigned long position;
   unsigned long * relation = la_inf->relation;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long * prime_count = qs_inf->prime_count;
   unsigned long num_primes = qs_inf->num_primes;
   mpz_t * Y_arr = la_inf->Y_arr;
   
   mpz_t pow;
   mpz_init(pow);
      
   memset(prime_count, 0, num_primes*sizeof(unsigned long));
      
   mpz_set_ui(X, 1);
   mpz_set_ui(Y, 1);
   
   for (unsigned long i = 0; i < ncols; i++)
   {
      if (get_null_entry(nullrows, i, l)) 
      {
         position = la_inf->matrix[i].orig*2*MAX_FACS;
         for (unsigned long j = 0; j < relation[position]; j++)
         {
            prime_count[relation[position+2*j+1]] +=
               (relation[position+2*j+2]);
         }
         mpz_mul(Y, Y, Y_arr[la_inf->matrix[i].orig]);
         if (i % 10 == 0) mpz_mod(Y, Y, N);
      }
   }

   for (unsigned long i = 0; i < num_primes; i++)
   {
      if (prime_count[i]) 
      {
         mpz_set_ui(pow, factor_base[i].p);
         mpz_powm_ui(pow, pow, prime_count[i]/2, N);
         mpz_mul(X, X, pow);
      } 
      if (i%10 == 0) mpz_mod(X, X, N);
   }

#if TEST
   for (unsigned long i = 0; i < num_primes; i++)
   {
      if ((prime_count[i] %2) != 0) printf("Error %ld, %ld, %ld\n", l, i, prime_count[i]);
   }
#endif

   mpz_clear(pow);
}

/*===========================================================================
   Collect relations:

   Function: Sets up batches of polynomials
             Do the sieving
             Evaluate candidates
             
   Returns:  The number of relations found with this batch of polynomials

===========================================================================*/

unsigned long collect_relations(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
{
   unsigned long s = poly_inf->s;
   unsigned long * poly_corr;
   unsigned long relations = 0;
   unsigned long ** A_inv2B = poly_inf->A_inv2B;
   unsigned long poly_index, j;
   unsigned long poly_add;
   unsigned long * B = poly_inf->B;
   unsigned long * B_terms = poly_inf->B_terms;
   unsigned long sieve_size = qs_inf->sieve_size;
   unsigned long small_primes = qs_inf->small_primes;
   unsigned long limbs = qs_inf->prec + 1;
   unsigned long limbs2;
   mp_limb_t msl;
   
#if CURVES
   static unsigned long count = 0;
   count++;
   if ((count & 7) == 0) printf("%ld curves\n", count*((1<<(s-1))-1));
#endif
   
   compute_A(qs_inf, poly_inf);
   compute_B_terms(qs_inf, poly_inf);
   compute_off_adj(qs_inf, poly_inf);
   compute_A_factor_offsets(qs_inf, poly_inf);
   compute_B_C(qs_inf, poly_inf);          
      
   for (poly_index = 1; poly_index < (1<<(s-1)); poly_index++)
   {
      for (j = 0; j < s; j++)
      {
         if (((poly_index >> j) & 1UL) != 0UL) break;
      }
      
      poly_add = (((poly_index >> j) & 2UL) != 0UL);
      
      poly_corr = A_inv2B[j];
           
      if (sieve_size <= SIEVE_BLOCK)
         do_sieving2(qs_inf, poly_inf, sieve);
      else
      {
         unsigned long blocks = sieve_size/SIEVE_BLOCK;
         unsigned long offset = SIEVE_BLOCK;
         unsigned long sieve_fill = qs_inf->sieve_fill;
        
         memset(sieve, sieve_fill, sieve_size);
         *(sieve+sieve_size) = 255;
         
         do_sieving(qs_inf, poly_inf, sieve, small_primes, SECOND_PRIME, SIEVE_BLOCK, 1, 0);
         for (long i = 1; i < blocks - 1; i++, offset += SIEVE_BLOCK)
            do_sieving(qs_inf, poly_inf, sieve + offset, small_primes, SECOND_PRIME, SIEVE_BLOCK, 0, 0);
         do_sieving(qs_inf, poly_inf, sieve + offset, small_primes, SECOND_PRIME, sieve_size - offset, 0, 1);
         
         do_sieving3(qs_inf, poly_inf, sieve, SECOND_PRIME, qs_inf->num_primes, sieve_size);
      }
         
      relations += evaluate_sieve(la_inf, qs_inf, poly_inf, sieve);
      
      update_offsets(poly_add, poly_corr, qs_inf, poly_inf);
      
      limbs2 = B_terms[j*limbs];
      if (((poly_add) && ((long)(limbs2 ^ B[0]) >= 0L)) || ((!poly_add) && ((long)(limbs2 ^ B[0]) < 0L)))
      {
         msl = mpn_add_n(B+1, B+1, B_terms + j*limbs + 1, limbs2);
         msl += mpn_add_n(B+1, B+1, B_terms + j*limbs + 1, limbs2);
         if (msl) mpn_add_1(B + limbs2 + 1, B + limbs2 + 1, limbs - limbs2 - 1, msl);
      } else 
      {
         msl = mpn_sub_n(B+1, B+1, B_terms + j*limbs + 1, limbs2);
         msl += mpn_sub_n(B+1, B+1, B_terms + j*limbs + 1, limbs2);
         if ((msl) && (limbs2 < limbs - 1)) msl = mpn_sub_1(B + limbs2 + 1, B + limbs2 + 1, limbs - limbs2 - 1, msl);
         if (msl)
         {
            negate_limbs(B+1, B+1, limbs - 1);
            B[0] = -B[0];
         }
      }
      if ((long) B[0] < 0)
      {
         B[0] = 1L - limbs;
         while (!B[FLINT_ABS(B[0])] && B[0]) B[0]++;
      } else
      {
         B[0] = limbs - 1L;
         while (!B[B[0]] && B[0]) B[0]--;
      }
      
      compute_B_C(qs_inf, poly_inf);          
      
      //compute_A_factor_offsets(qs_inf, poly_inf);    
   }
   
   relations += merge_relations(la_inf);
   
   return relations;
}

/*===========================================================================
   Main Quadratic Sieve Factoring Routine:

   Function: Finds the factors of a number using the quadratic sieve
             Assume n is odd, not a perfect power, positive and not a prime 
             Returns 0 if factorisation was unsuccessful
             Returns 1 if factorisation was successful
             If a small factor is found, it is returned and the QS is not run

===========================================================================*/

int F_mpz_factor_mpQS(F_mpz_factor_t factors, mpz_t N)
{
   unsigned long small_factor;
   unsigned long rels_found = 0;
   unsigned long prec;
   
   QS_t qs_inf; 
   poly_t poly_inf;
   linalg_t la_inf;
   
   qs_inf.bits = mpz_sizeinbase(N,2);
   if (qs_inf.bits < MINBITS) return 0; // Number too small for mpQS 
   
   prec = (qs_inf.bits + max_mult_size - 1)/FLINT_BITS + 1;
   
   qs_inf.n = (fmpz_t) flint_stack_alloc(prec+1);
   qs_inf.n[prec] = 0;
   mpz_to_fmpz(qs_inf.n, N); // set n to the number to be factored
   
   qs_inf.prec = (prec + 1)/2;
   
   small_factor = knuth_schroeppel(&qs_inf); // Compute multiplier and some FB primes
   if (small_factor) goto cleanup_2;

#if QS_INFO
   printf("Multiplier = %ld\n", qs_inf.k);
#endif
   
   mpz_init(qs_inf.mpz_n);
   mpz_set(qs_inf.mpz_n, N);
   mpz_mul_ui(qs_inf.mpz_n, qs_inf.mpz_n, qs_inf.k);
   qs_inf.bits = mpz_sizeinbase(qs_inf.mpz_n, 2);
   if (qs_inf.bits < MINBITS) 
   {
      small_factor = 0; // Number too small for mpQS 
      goto cleanup_1a;
   }
   mpz_to_fmpz(qs_inf.n, qs_inf.mpz_n); // set n to the number to be factored times k
   
   primes_init(&qs_inf);
   sqrts_init(&qs_inf);
      
   if (qs_inf.bits < MINBITS) 
   {
      small_factor = 0; // kn too big for tinyQS
      goto cleanup_1;
   }
   
   small_factor = compute_factor_base(&qs_inf); // Computes the factor base primes and modular square roots
   if (small_factor) goto cleanup_1;
   compute_sizes(&qs_inf);
   get_sieve_params(&qs_inf);
   
   poly_init(&qs_inf, &poly_inf, N);
   linear_algebra_init(&la_inf, &qs_inf, &poly_inf);
   
   unsigned char * sieve = (unsigned char *) flint_stack_alloc_bytes(qs_inf.sieve_size+1);
   while (rels_found < qs_inf.num_primes + EXTRA_RELS)
   {
      rels_found += collect_relations(&la_inf, &qs_inf, &poly_inf, sieve);
   }
   flint_stack_release(); // release sieve
   
   la_col_t * matrix = la_inf.matrix;
   unsigned long ncols = qs_inf.num_primes + EXTRA_RELS;
   unsigned long nrows = qs_inf.num_primes;

   reduce_matrix(&nrows, &ncols, matrix); // Do some filtering on the matrix
   
   u_int64_t* nullrows;
   do {
      nullrows = block_lanczos(nrows, 0, ncols, matrix); // Linear algebra (block Lanczos)
   } while (nullrows == NULL); 
   
   unsigned long i, j, mask;
     
   for (i = 0, mask = 0; i < ncols; i++)
      mask |= nullrows[i];

   for (i = j = 0; i < 64; i++) {
		if (mask & ((u_int64_t)(1) << i))
			j++;
   }

#if QS_INFO
   printf("%ld nullspace vectors found\n", j);
#endif

   qs_inf.prime_count = (unsigned long *) flint_stack_alloc(qs_inf.num_primes);
   
   mpz_t X, Y, F, Q, R;
   mpz_init(X);
   mpz_init(Y);
   mpz_init(F);
   mpz_init(Q);
   mpz_init(R);
   
   mpz_set(F, N);
    
#if PRINT_FACTORS
   gmp_printf("Factors of %Zd:\n", N);
#endif

   for (unsigned long l = 0; l < 64; l++)
   {
      if (mask & ((u_int64_t)(1) << l))
      {
         square_root(X, Y, &qs_inf, &la_inf, nullrows, ncols, l, N); 
         mpz_sub(X, X, Y);
         mpz_gcd(X, X, N);
         if ((mpz_cmp(X, N) != 0) && (mpz_cmp_ui(X, 1) != 0))
         {
#if PRINT_FACTORS
            gmp_printf("%Zd\n", X);
#endif
            if (mpz_probab_prime_p(X, 10)) 
            {
               mpz_fdiv_qr(Q, R, F, X);
               if (mpz_cmp_ui(R, 0) == 0) mpz_set(F, Q);
            }
            if (mpz_cmp_ui(F, 1) == 0) break; 
         }
      }
   }
   
   small_factor = 1; // sieve was successful
   mpz_clear(Q);
   mpz_clear(R);
   mpz_clear(F);
   mpz_clear(X);
   mpz_clear(Y);
   flint_stack_release(); // release prime_count
   linear_algebra_clear(&la_inf, &qs_inf);
   poly_clear(&poly_inf);
   sizes_clear();
cleanup_1:
   sqrts_clear(); // release modular square roots
   primes_clear(); // release factor_base
cleanup_1a:
   mpz_clear(qs_inf.mpz_n);
cleanup_2:
   flint_stack_release(); // release n

   return small_factor;    
}

/*===========================================================================
   Main Program:

   Function: Factors a user specified number using the quadratic sieve


===========================================================================*/
int main(int argc, unsigned char *argv[])
{
    mpz_t N;
    mpz_init(N); 
    
    F_mpz_factor_t factors;
    
    printf("Input number to factor [ >= 40 decimal digits ] : "); 
    gmp_scanf("%Zd", N); getchar();
    
    F_mpz_factor_mpQS(factors, N);
    
    mpz_clear(N);
}

/*int main(int argc, unsigned char *argv[])
{
    mpz_t N;
    mpz_init(N); 
    
    F_mpz_factor_t factors;
    
    unsigned long factor;
    unsigned long failed = 0;
    unsigned long small_factors = 0;
    unsigned long succeed = 0;
    unsigned long bits1, bits2, bits3, i;
    
    for (i = 0; i < 1; i++)
    {
       mpz_set_ui(N, long_nextprime(long_randint(4000000000000000000UL)));
       mpz_mul_ui(N, N, long_nextprime(long_randint(4000000000000000000UL)));
       mpz_mul_ui(N, N, long_nextprime(long_randint(4000000000000000UL)));
       mpz_mul_ui(N, N, long_nextprime(long_randint(1000000000UL)));
       //bits1 = long_randint(41UL)+13UL;
       //bits2 = long_randint(22UL)+13UL;
       //bits3 = long_randint(22UL)+13UL;
       //mpz_set_ui(N, long_nextprime(long_randint((1UL<<bits1)-1UL)+1UL));
       //mpz_mul_ui(N, N, long_nextprime(long_randint((1UL<<bits2)-1UL)+1UL));
       //mpz_mul_ui(N, N, long_nextprime(long_randint((1UL<<bits3)-1UL)+1UL));

#if QS_INFO
       gmp_printf("Factoring %Zd\n", N);
#endif

       factor = F_mpz_factor_mpQS(factors, N);
       if (!factor) failed++;
       if (factor > 1) small_factors++;
       if (factor == 1) succeed++; 
    }
    
    printf("mpQS succeeded %ld times, found a small factor %ld times\n", succeed, small_factors);
    printf("and failed %ld times\n", failed);
    
    mpz_clear(N);
}*/
