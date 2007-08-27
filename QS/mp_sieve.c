/******************************************************************************

 mp_sieve.c
 
 Routines for doing and managing sieving

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../flint.h"
#include "../long_extras.h"

#include "common.h"
#include "mp_poly.h"
#include "mp_linear_algebra.h"
#include "mp_sieve.h"

void get_sieve_params(QS_t * qs_inf)
{
   unsigned long bits = qs_inf->bits;
   unsigned long i;
   
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long num_primes = qs_inf->num_primes;
   unsigned long mult = qs_inf->num_primes;
   
   for (i = 0; i < PTABSIZE; i++)
   {
      if (prime_tab[i][0] > bits) break;
   }
   
   qs_inf->sieve_size = prime_tab[i-1][2]; 
   qs_inf->small_primes = prime_tab[i-1][3]; 
   qs_inf->large_prime = prime_tab[i-1][4]*factor_base[num_primes-1].p;
   qs_inf->error_bits = round(log(qs_inf->large_prime)/log(2.0))+3; // 2, 5, 6 
   printf("Error bits = %ld\n", qs_inf->error_bits);
}

void do_sieving2(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
{
   unsigned long num_primes = qs_inf->num_primes;
   u_int32_t * soln1 = poly_inf->soln1;
   u_int32_t * soln2 = poly_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long sieve_size = qs_inf->sieve_size;
   unsigned char * end = sieve + sieve_size;
   unsigned char * sizes = qs_inf->sizes;
   unsigned long sieve_fill = qs_inf->sieve_fill;
   unsigned long small_primes = qs_inf->small_primes;
   unsigned char * bound;
   unsigned char * pos1;
   unsigned char * pos2;
   unsigned long size;
   unsigned long p;
   
   memset(sieve, sieve_fill, sieve_size);
   *end = 255;
   
   const unsigned long second_prime = FLINT_MIN(SECOND_PRIME, num_primes);
   
   for (unsigned long prime = small_primes; prime < second_prime; prime++) 
   {
      if (soln2[prime] == -1) continue;
      
      p = factor_base[prime].p;
      size = sizes[prime];
      pos1 = sieve + soln1[prime];
      pos2 = sieve + soln2[prime];
      bound = end - p;
        
      while (bound - pos1 > 0)  
      {  
         (*pos1)+=size, pos1+=p, (*pos2)+=size, pos2+=p;
      }
      if ((end - pos1 > 0) && (end - pos2 > 0))
      { 
         (*pos1)+=size, pos1+=p, (*pos2)+=size, pos2+=p;
      }
      if (end - pos2 > 0)
      { 
         (*pos2)+=size;
      }
      if (end - pos1 > 0)
      { 
         (*pos1)+=size;
      }
   }
   
   for (unsigned long prime = second_prime; prime < num_primes; prime++) 
   {
      p = factor_base[prime].p;
      size = sizes[prime];
      pos1 = sieve + soln1[prime];
      pos2 = sieve + soln2[prime];
        
      if (end - pos2 > 0)
      { 
         (*pos2)+=size;
      }
      if (end - pos1 > 0)
      { 
         (*pos1)+=size;
      }
   }
}

void update_offsets(unsigned long poly_add, u_int32_t * poly_corr, 
             QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long num_primes = qs_inf->num_primes;
   u_int32_t * soln1 = poly_inf->soln1;
   u_int32_t * soln2 = poly_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long p, correction;
   
   for (unsigned long prime = 2; prime < num_primes; prime++) 
   {
      if (soln2[prime] == -1) continue;
      p = factor_base[prime].p;
      correction = (poly_add ? p - poly_corr[prime] : poly_corr[prime]);
      soln1[prime] += correction;
      if (soln1[prime] >= p) soln1[prime] -= p;
      soln2[prime] += correction;
      if (soln2[prime] >= p) soln2[prime] -= p; 
   }
   
}  


void do_sieving(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve, 
                  unsigned long first_prime, unsigned long second_prime, 
                  unsigned long M, int first, int last)
{
   unsigned long num_primes = qs_inf->num_primes;
   u_int32_t * soln1 = poly_inf->soln1;
   u_int32_t * soln2 = poly_inf->soln2;
   u_int32_t * posn1 = poly_inf->posn1;
   u_int32_t * posn2 = poly_inf->posn2;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long small_primes = qs_inf->small_primes;
   unsigned long p, correction;
   register unsigned char * position;
   unsigned char * end = sieve + M;
   unsigned char * sizes = qs_inf->sizes;
   register unsigned char * pos1;
   unsigned char * pos2;
   register unsigned char * bound;  
   unsigned long size;
   long diff;
   
   for (unsigned long prime = first_prime; prime < second_prime; prime++) 
   {
      if (soln2[prime] == -1) continue;
      
      p = factor_base[prime].p;
      size = sizes[prime];
      if (first)
      {
         pos1 = sieve + soln1[prime];
         pos2 = sieve + soln2[prime];
      } else 
      {
         pos1 = sieve + posn1[prime];
         pos2 = sieve + posn2[prime];         
      }
      
      bound = end - p;
        
      while (bound - pos1 > 0)  
      {  
         (*pos1)+=size, pos1+=p, (*pos2)+=size, pos2+=p;
      }
      if ((end - pos1 > 0) && (end - pos2 > 0))
      { 
         (*pos1)+=size, pos1+=p, (*pos2)+=size, pos2+=p;
      }
      if (end - pos2 > 0)
      { 
         (*pos2)+=size, pos2+=p;
      }
      if (end - pos1 > 0)
      { 
         (*pos1)+=size, pos1+=p;
      }

      if (!last)
      { 
         posn1[prime] = pos1 - sieve;
         posn2[prime] = pos2 - sieve;
      }
   }
   
}

void do_sieving3(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve, 
                  unsigned long first_prime, unsigned long second_prime, 
                  unsigned long M)
{
   u_int32_t * soln1 = poly_inf->soln1;
   u_int32_t * soln2 = poly_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long p;
   unsigned char * end = sieve + M;
   unsigned char * sizes = qs_inf->sizes;
   register unsigned char * pos1;
   unsigned char * pos2;
   unsigned long size;
   
   for (unsigned long prime = first_prime; prime < second_prime; prime++) 
   {
      if (soln2[prime] == -1) continue;
      
      p = factor_base[prime].p;
      size = sizes[prime];
      pos1 = sieve + soln1[prime];
      pos2 = sieve + soln2[prime];
              
      while (end - pos2 > 0)
      { 
         (*pos2)+=size, pos2+=p;
      }
      while (end - pos1 > 0)
      { 
         (*pos1)+=size, pos1+=p;
      } 
   }  
}

/*==========================================================================
   evaluate_candidate:

   Function: determine whether a given sieve entry is a relation 

===========================================================================*/

unsigned long evaluate_candidate(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf, 
                          unsigned long i, unsigned char * sieve)
{
   unsigned long bits, exp, extra_bits, modp, prime;
   unsigned long num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   u_int32_t * soln1 = poly_inf->soln1;
   u_int32_t * soln2 = poly_inf->soln2;
   unsigned long * small = la_inf->small;
   unsigned long sieve_fill = qs_inf->sieve_fill;
   unsigned long sieve_size = qs_inf->sieve_size;
   fac_t * factor = la_inf->factor;
   mpz_t * A = &poly_inf->A_mpz;
   mpz_t * B = &poly_inf->B_mpz;
   unsigned long error_bits = qs_inf->error_bits;
   unsigned long small_primes = qs_inf->small_primes;
   unsigned long large_prime = qs_inf->large_prime;
   unsigned long num_factors = 0;
   unsigned long j;
   mpz_t * C = &poly_inf->C;
   unsigned long relations = 0;
   double pinv;
   const unsigned long second_prime = FLINT_MIN(SECOND_PRIME, num_primes);
   
   mpz_t X, Y, res, p;
   mpz_init(X); 
   mpz_init(Y); 
   mpz_init(res); 
   mpz_init(p); 

#if POLYS
   printf("X = %ld\n", i);    
   gmp_printf("%ZdX^2+2*%ZdX+%Zd\n", A, B, C);
#endif

   mpz_set_ui(X, i);
   mpz_sub_ui(X, X, sieve_size/2); //X
              
   mpz_mul(Y, X, *A);
   mpz_add(Y, Y, *B);  // Y = AX+B
   mpz_add(res, Y, *B);  
   
   mpz_mul(res, res, X);  
   mpz_add(res, res, *C); // res = AX^2+2BX+C
           
   bits = mpz_sizeinbase(res, 2);
   bits -= error_bits; 
   extra_bits = 0;
   
   mpz_set_ui(p, 2); // divide out by powers of 2
   exp = mpz_remove(res, res, p);

#if RELATIONS
   if (exp) printf("2^%ld ", exp);
#endif
   extra_bits += exp;
   small[1] = exp;
     
   if (factor_base[0].p != 1) // divide out powers of the multiplier
   {
      mpz_set_ui(p, factor_base[0].p);
      exp = mpz_remove(res, res, p);
      if (exp) extra_bits += qs_inf->sizes[0];
      small[0] = exp;
#if RELATIONS
      if (exp) printf("%ld^%ld ", factor_base[0].p, exp); 
#endif
   } else small[0] = 0;
   
   for (unsigned long j = 2; j < small_primes; j++) // pull out small primes
   {
      prime = factor_base[j].p;
      pinv = factor_base[j].pinv;
      modp = z_mod2_precomp(i, prime, pinv);
      if ((modp == soln1[j]) || (modp == soln2[j]))
      {
         mpz_set_ui(p, prime);
         exp = mpz_remove(res, res, p);
         extra_bits += qs_inf->sizes[j];
         small[j] = exp;
#if RELATIONS
         gmp_printf("%Zd^%ld ", p, exp); 
#endif
      } else small[j] = 0;
   }
   
   if (extra_bits + sieve[i] > bits+sieve_fill)
   {
      sieve[i] += extra_bits - sieve_fill;
      for (j = small_primes; (j < second_prime) && (extra_bits < sieve[i]); j++) // pull out remaining primes
      {
         prime = factor_base[j].p;
         pinv = factor_base[j].pinv;
         modp = z_mod2_precomp(i, prime, pinv);
         if (soln2[j] != -1)
         {
            if ((modp == soln1[j]) || (modp == soln2[j]))
            {
               mpz_set_ui(p, prime);
               exp = mpz_remove(res, res, p);          
#if RELATIONS
               gmp_printf("%Zd^%ld ", p, exp);
#endif
               extra_bits += qs_inf->sizes[j];
               factor[num_factors].ind = j;
               factor[num_factors++].exp = exp; 
            }
         } else
         {
            mpz_set_ui(p, prime);
            exp = mpz_remove(res, res, p);
            factor[num_factors].ind = j;
            factor[num_factors++].exp = exp+1; 
#if RELATIONS
            if (exp) gmp_printf("%Zd^%ld ", p, exp);
#endif
         }    
      }
      for ( ; (j < num_primes) && (extra_bits < sieve[i]); j++) // pull out remaining primes
      {
         if ((i == soln1[j]) || (i == soln2[j]))
         {
            prime = factor_base[j].p;
            mpz_set_ui(p, prime);
            exp = mpz_remove(res, res, p);          
#if RELATIONS
            gmp_printf("%Zd^%ld ", p, exp);
#endif
            extra_bits += qs_inf->sizes[j];
            factor[num_factors].ind = j;
            factor[num_factors++].exp = exp; 
         }     
      }
      if (mpz_cmpabs_ui(res, 1) == 0) // We've found a relation
      {
         unsigned long * A_ind = poly_inf->A_ind;
         for (unsigned long i = 0; i < poly_inf->s; i++) // Commit any outstanding A factors
         {
            if (A_ind[i] >= j)
            {
               factor[num_factors].ind = A_ind[i];
               factor[num_factors++].exp = 1; 
            }
         }
         la_inf->num_factors = num_factors;
         relations += insert_relation(qs_inf, la_inf, poly_inf, Y);  // Insert the relation in the matrix
         if (la_inf->num_relations >= qs_inf->num_primes + EXTRA_RELS + 200)
         {
            printf("Error: too many duplicate relations!\n");
            abort();
         }
         goto cleanup;
      } else if(mpz_cmpabs_ui(res, large_prime) < 0) 
      {
         if (mpz_sgn(res) < 0) mpz_neg(res, res);
         unsigned long * A_ind = poly_inf->A_ind;
         for (unsigned long i = 0; i < poly_inf->s; i++) // Commit any outstanding A factors
         {
            if (A_ind[i] >= j)
            {
               factor[num_factors].ind = A_ind[i];
               factor[num_factors++].exp = 1; 
            }
         }
         la_inf->num_factors = num_factors;
         relations += insert_lp_relation(qs_inf, la_inf, poly_inf, Y, res);  // Insert the relation in the matrix                    
         goto cleanup;
      }
   }
  
cleanup:
#if RELATIONS
   printf("\n");
#endif
   mpz_clear(X);
   mpz_clear(Y);
   mpz_clear(res);
   mpz_clear(p);
      
   return relations;
}

/*==========================================================================
   evaluateSieve:

   Function: searches sieve for relations and sticks them into a matrix

===========================================================================*/
unsigned long evaluate_sieve(linalg_t * la_inf, QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
{
   unsigned long i = 0;
   unsigned long j = 0;
   unsigned long * sieve2 = (unsigned long *) sieve;
   unsigned long sieve_size = qs_inf->sieve_size;
   unsigned long rels = 0;
     
   while (j < sieve_size/sizeof(unsigned long))
   {
#if FLINT_BITS == 64
      while (!(sieve2[j] & 0x8080808080808080U)) j++;
#else
      while (!(sieve2[j] & 0x80808080U)) j++;
#endif
      i = j*sizeof(unsigned long);
      while ((i < (j+1)*sizeof(unsigned long)) && (i < sieve_size))
      {
         if (sieve[i] > 128) 
         {
             rels += evaluate_candidate(la_inf, qs_inf, poly_inf, i, sieve);
         }
         i++;
      }
      j++;
   }
   return rels;
}
