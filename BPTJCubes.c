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
/****************************************************************************

   BPTJ_cubes.c: Finds solutions to x^3 + y^3 + z^3 = k
                Based on the algorithm of Beck, Pine, Tarrant and Jensen
                
   Simultaneously searches for solutions for k = k1, k2, k3
   Searches from T = START to STOP
   
   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "long_extras.h"

#define START 1000000 
#define STOP 200000000
#define k1 1965
#define k2 1986
#define k3 1991

#define NUMPRIMES 4000
#define TABLESIZE 1800000 // Must be a multiple of CACHEBLOCK
#define CACHEBLOCK 60000
#define RABIN 6

mpz_t temp, temp2, D;

gmp_randstate_t randstate;

static inline
int test_root(unsigned long root, unsigned long T, unsigned long k)
{
   long d;
   
   if (root > T/2)
   {
      mpz_set_ui(D, root);
      mpz_pow_ui(D, D, 3);
      mpz_sub_ui(D, D, k);
   } else
   {
      mpz_set_ui(D, T-root);
      mpz_pow_ui(D, D, 3);
      mpz_add_ui(D, D, k);   
   }
   
   mpz_mul_2exp(D, D, 2);
   mpz_divexact_ui(D, D, T);
   mpz_set_ui(temp2, T);
   mpz_submul_ui(D, temp2, T);
   mpz_mul_ui(D, D, 3); 
   
   if (mpz_sgn(D) >= 0)
   {
      mpz_sqrtrem(D, temp2, D);
      if (!mpz_sgn(temp2))
      {
         if (!mpz_fdiv_r_ui(temp, D, 3))
         {
            d = mpz_fdiv_r_ui(temp, D, 6);
            if ((((3*T)%6) == d) || (((3*T)%6) == 6-d)) return 1;
         }
      }
   }
   
   return 0;

}

int main()
{
   gmp_randinit_default(randstate);
   
   FILE * file1 = fopen("output.log","w");
   
   unsigned long T = z_nextprime(START, 0);
   double Tinv;
   unsigned long cuberoot1;
   unsigned long root1, root2, root3;
   unsigned long s = 0, t, p;
   unsigned char * current;
   
   unsigned char * table = (unsigned char *) malloc(TABLESIZE);
   unsigned long * mod = (unsigned long *) malloc(NUMPRIMES*sizeof(unsigned long));
   unsigned long * prime = (unsigned long *) malloc(NUMPRIMES*sizeof(unsigned long));
   
   s = 3;
   unsigned long i;
   for (i = 0; i < NUMPRIMES; i++)
   {
      prime[i] = s;
      s = z_nextprime(s, 0);
   }
   unsigned long i;
   for (i = 0; i < NUMPRIMES; i++)
   {
      s = (T%prime[i]);
      if (s == 0) s = prime[i];
      mod[i] = prime[i]-s;
   }
   
   
   mpz_init(temp);
   mpz_init(temp2);
   mpz_init(D);
   
   unsigned long Ttab = 0;
   
   while (T < STOP)
   {
      memset(table, 0, TABLESIZE);
      unsigned long offset;
      for (offset = 0; offset < TABLESIZE; offset+=CACHEBLOCK)
      {
        unsigned long i;
        for (i = 0; i < NUMPRIMES; i++)
        {
          s = mod[i];
          p = prime[i];
          current = table + offset;
          for ( ; s < CACHEBLOCK; s += p) current[s] = 1;
          s -= CACHEBLOCK;
          mod[i] = s;
        }
      }
      
      mpz_set_ui(temp, T);
      while (!mpz_probab_prime_p(temp,3))
      {
         T += 2;
         Ttab += 2;
         while (table[Ttab])
         {
            T += 2;
            Ttab += 2;
         }  
         mpz_set_ui(temp, T);
      } 
      Tinv = z_precompute_inverse(T);
      
      while (Ttab < TABLESIZE)
      {
         root1 = z_cuberootmod(&cuberoot1, k1, T);
         if (root1)
         {
            if (cuberoot1 != 1)
            {
               root2 = z_mulmod_precomp(root1, cuberoot1, Tinv, T);
               root3 = z_mulmod_precomp(root2, cuberoot1, Tinv, T);
            }
            if (test_root(root1, T, k1)) 
            {
               fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k1, T, root1);
               fflush(file1);
               abort();
            }
         
            if (cuberoot1 != 1)
            {
               if (test_root(root2, T, k1)) 
               {
                  fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k1, T, root2);
                  fflush(file1);
                  abort();
               }
               if (test_root(root3, T, k1)) 
               {
                  fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k1, T, root3);
                  fflush(file1);
                  abort();
               }
            }
         } 
         root1 = z_cuberootmod(&cuberoot1, k2, T);
         if (root1)
         {
            if (cuberoot1 != 1)
            {
               root2 = z_mulmod_precomp(root1, cuberoot1, Tinv, T);
               root3 = z_mulmod_precomp(root2, cuberoot1, Tinv, T);
            }
            if (test_root(root1, T, k2)) 
            {
               fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k2, T, root1);
               fflush(file1);
               abort();
            }
         
            if (cuberoot1 != 1)
            {
               if (test_root(root2, T, k2)) 
               {
                  fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k2, T, root2);
                  fflush(file1);
                  abort();
               }
               if (test_root(root3, T, k2)) 
               {
                  fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k2, T, root3);
                  fflush(file1);
                  abort();
               }
            }
         } 
         root1 = z_cuberootmod(&cuberoot1, k3, T);
         if (root1)
         {
            if (cuberoot1 != 1)
            {
               root2 = z_mulmod_precomp(root1, cuberoot1, Tinv, T);
               root3 = z_mulmod_precomp(root2, cuberoot1, Tinv, T);
            }
            if (test_root(root1, T, k3)) 
            {
               fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k3, T, root1);
               fflush(file1);
               abort();
            }
         
            if (cuberoot1 != 1)
            {
               if (test_root(root2, T, k3)) 
               {
                  fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k3, T, root2);
                  fflush(file1);
                  abort();
               }
               if (test_root(root3, T, k3)) 
               {
                  fprintf(file1,"k = %ld, T = %ld, root = %ld\n", k3, T, root3);
                  fflush(file1);
                  abort();
               }
            }
         } 

      do
      {
         do 
         {
            T += 2;
            Ttab += 2;
         } while (table[Ttab] && (Ttab < TABLESIZE));
         Tinv = z_precompute_inverse(T);
      } while (!z_isprime_precomp(T, Tinv) && (Ttab < TABLESIZE));
            
      t++;
      if ((t%1000000UL) == 0) 
      {
         fprintf(file1,"Checkpoint T = %ld\n", T);
         fflush(file1);
      }
      
      }
      Ttab -= TABLESIZE;
   }
   
   mpz_clear(temp);
   mpz_clear(temp2);
   mpz_clear(D);
   gmp_randclear(randstate);
   
   return 0;
}
