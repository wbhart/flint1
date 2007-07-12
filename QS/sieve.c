/******************************************************************************

 sieve.c
 
 Routines for doing and managing sieving

 (C) 2006 William Hart

******************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <string.h>

#include "../flint.h"

#include "common.h"
#include "poly.h"

void do_sieving(unsigned long poly_add, unsigned long * poly_corr, 
             QS_t * qs_inf, poly_t * poly_inf, unsigned char * sieve)
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
   
   memset(sieve, SIEVE_SIZE, 0);
   
   for (unsigned long prime = 2; prime < num_primes; prime++) 
   {
      if (soln2[prime] == -1L) continue;
      
      p = factor_base[prime].p;
      size = sizes[prime];
      pos1 = sieve+soln1[prime];
      pos2 = sieve+soln2[prime];
      diff = pos2-pos1;
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
      
      correction = poly_add ? p - poly_corr[prime] : poly_corr[prime];
      soln1[prime] += correction;
      while (soln1[prime] >= p) soln1[prime] -= p;
      soln2[prime] += correction;
      while (soln2[prime] >= p) soln2[prime] -= p;   
   }
}
