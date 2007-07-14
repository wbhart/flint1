/******************************************************************************

 sieve.c
 
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
#include "poly.h"

void do_sieving(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
{
   unsigned long num_primes = qs_inf->num_primes;
   unsigned long * soln1 = poly_inf->soln1;
   unsigned long * soln2 = poly_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long p, correction;
   register unsigned char * position;
   unsigned char * end = sieve + SIEVE_SIZE;
   unsigned char * sizes = qs_inf->sizes;
   register unsigned char * pos1;
   unsigned char * pos2;
   register unsigned char * bound;  
   unsigned long size;
   long diff;
   
   memset(sieve, 0, SIEVE_SIZE);
   *end = 255;
   
   for (unsigned long prime = 7; prime < num_primes; prime++) 
   {
      if (soln2[prime] == -1L) continue;
      
      p = factor_base[prime].p;
      size = sizes[prime];
      pos1 = sieve + soln1[prime];
      pos2 = sieve + soln2[prime];
      diff = pos2 - pos1;
      bound = end - p*4;
        
      while (bound - pos1 > 0)  
      {  
         (*pos1)+=size, (*(pos1+diff))+=size, pos1+=p;
         (*pos1)+=size, (*(pos1+diff))+=size, pos1+=p;
         (*pos1)+=size, (*(pos1+diff))+=size, pos1+=p;
         (*pos1)+=size, (*(pos1+diff))+=size, pos1+=p;
      }
      while ((end - pos1 > 0) && (end - pos1 - diff > 0))
      { 
         (*pos1)+=size, (*(pos1+diff))+=size, pos1+=p;
      }
      pos2 = pos1+diff;
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

void update_offsets(unsigned long poly_add, unsigned long * poly_corr, 
             QS_t * qs_inf, poly_t * poly_inf)
{
   unsigned long num_primes = qs_inf->num_primes;
   unsigned long * soln1 = poly_inf->soln1;
   unsigned long * soln2 = poly_inf->soln2;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long p, correction;

   for (unsigned long prime = 2; prime < num_primes; prime++) 
   {
      p = factor_base[prime].p;
      correction = (poly_add ? p - poly_corr[prime] : poly_corr[prime]);
      soln1[prime] += correction;
      if (soln1[prime] >= p) soln1[prime] -= p;
      if (soln2[prime] == -1L) continue;
      soln2[prime] += correction;
      if (soln2[prime] >= p) soln2[prime] -= p; 
   }
}  

/*==========================================================================
   evaluate_candidate:

   Function: determine whether a given sieve entry is a relation or 
             partial relation

===========================================================================*/
int evaluate_candidate(QS_t * qs_inf, poly_t * poly_inf, 
                          unsigned long i, unsigned char * sieve)
{
   unsigned long bits, exp, extra_bits, modp, prime;
   unsigned long num_primes = qs_inf->num_primes;
   prime_t * factor_base = qs_inf->factor_base;
   unsigned long * soln1 = poly_inf->soln1;
   unsigned long * soln2 = poly_inf->soln2;
   unsigned long A = poly_inf->A;
   unsigned long B = poly_inf->B;
   mpz_t * C = &poly_inf->C;
   int relation = 0;
   double pinv;
   
   mpz_t X, Y, res, p;
   mpz_init(X); 
   mpz_init(Y); 
   mpz_init(res); 
   mpz_init(p); 
    
#if POLYS
   printf("X = %ld\n", i);
   printf("%ldX^2+2*%ldX+%ld\n", A, B, C);
#endif

   mpz_set_ui(X, i);
   mpz_sub_ui(X, X, SIEVE_SIZE/2); //X
              
   mpz_mul_ui(Y, X, A);
   if ((long) B < 0) 
   {
      mpz_sub_ui(Y, Y, -B);  // Y = AX+B
      mpz_sub_ui(res, Y, -B);  
   }
   else 
   {
      mpz_add_ui(Y, Y, B);
      mpz_add_ui(res, Y, B);
   }
   mpz_mul(res, res, X);  
   mpz_add(res, res, *C); // res = AX^2+2BX+C
           
   bits = mpz_sizeinbase(res, 2);
   bits -= 10; 
   extra_bits = 0;
   
   mpz_set_ui(p, 2); // divide out by powers of 2
   exp = mpz_remove(res, res, p);

#if RELATIONS
   if (exp) printf("2^%ld ", exp);
#endif
   extra_bits += exp;
     
   if (factor_base[0].p != 1) // divide out powers of the multiplier
   {
      mpz_set_ui(p, factor_base[0].p);
      exp = mpz_remove(res, res, p);
      if (exp) extra_bits += qs_inf->sizes[0];
#if RELATIONS
      if (exp) printf("%ld^%ld ", factor_base[0].p, exp); 
#endif
   }
   
   for (unsigned long j = 2; j < 7; j++) // pull out small primes
   {
      prime = factor_base[j].p;
      pinv = factor_base[j].pinv;
      modp = long_mod63_precomp(i, prime, pinv);
      if ((modp == soln1[j]) || (modp == soln2[j]))
      {
         mpz_set_ui(p, prime);
         exp = mpz_remove(res, res, p);
         if (exp) extra_bits += qs_inf->sizes[j];
#if RELATIONS
         if (exp) gmp_printf("%Zd^%ld ", p, exp); 
#endif
      }
   }
   
   if (extra_bits + sieve[i] > bits)
   {
      sieve[i] += extra_bits;
      for (unsigned long j = 7; (j < num_primes) && (extra_bits < sieve[i]); j++) // pull out remaining primes
      {
         prime = factor_base[j].p;
         pinv = factor_base[j].pinv;
         modp = long_mod63_precomp(i, prime, pinv);
         if (soln2[j] != -1L)
         {
            if ((modp == soln1[j]) || (modp == soln2[j]))
            {
               mpz_set_ui(p, prime);
               exp = mpz_remove(res, res, p);
#if RELATIONS
               if (exp) gmp_printf("%Zd^%ld ", p, exp);
#endif
               if (exp) extra_bits += qs_inf->sizes[j]; 
            }
         } else
         {
            mpz_set_ui(p, prime);
            exp = mpz_remove(res, res, p);
#if RELATIONS
            if (exp) gmp_printf("%Zd^%ld ", p, exp);
#endif
            //if (exp) extra_bits += qs_inf->sizes[j];
         }    
      }
      if (mpz_cmpabs_ui(res, 1) == 0) 
      {
         relation = 1;
         goto cleanup;
      }
   }
#if RELATIONS
   printf("\n");
#endif
  
cleanup:
   mpz_clear(X);
   mpz_clear(Y);
   mpz_clear(res);
   mpz_clear(p);
      
   return relation;
}

/*==========================================================================
   evaluateSieve:

   Function: searches sieve for relations and sticks them into a matrix

===========================================================================*/
unsigned long evaluate_sieve(QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
{
   unsigned long i = 0;
   unsigned long j=0;
   unsigned long * sieve2 = (unsigned long *) sieve;
   unsigned long rels = 0;
     
   while (j < SIEVE_SIZE/sizeof(unsigned long))
   {
      while (!(sieve2[j] & 0xE0E0E0E0E0E0E0E0U)) j++;
      i = j*sizeof(unsigned long);
      while ((i < (j+1)*sizeof(unsigned long)) && (i < SIEVE_SIZE))
      {
         if (sieve[i] > 51) 
         {
            if (evaluate_candidate(qs_inf, poly_inf, i, sieve)) rels++;
         }
         i++;
      }
      j++;
   }
   
   return rels;
}
