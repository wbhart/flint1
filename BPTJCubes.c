/****************************************************************************

   BPTJ_cubes.c: Finds solutions to x^3 + y^3 + z^3 = k
                Based on the algorithm of Beck, Pine, Tarrant and Jensen
   
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

// Checked 42 up to 1.35x10^12 (from 6.5x10^11 - Siksek's machine
// Checked 114 up to 10^11 (from 100000)
// Check 195, 290, 452 up to 10^11 (from 100000)

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
   /*int cubes[72];
   int OK[72];
   for (unsigned long i = 0; i < 72; i++) 
   {
      OK[i] = 0;
      cubes[i] = 0;
   }
   for (unsigned long i = 0; i < 72; i++) 
   {
      cubes[(i*i*i)%72] = 1;
   }  
   for (unsigned long i = 0; i < 72; i++)
   {
       for (unsigned long j = 0; j < 72; j++)
       {
           if (cubes[(i*i*i+j*j*j-k)%72]) OK[(i+j)%72] = 1;
       }
   }*/
    
   unsigned long T = long_nextprime(START);
   unsigned long Tinv_hi, Tinv_lo;
   unsigned long cuberoot1;
   unsigned long root1, root2, root3;
   unsigned long s = 0, t, p;
   //unsigned long Tmod72, oldT;
   unsigned char * current;
   
   unsigned char * table = (unsigned char *) malloc(TABLESIZE);
   unsigned long * mod = (unsigned long *) malloc(NUMPRIMES*sizeof(unsigned long));
   unsigned long * prime = (unsigned long *) malloc(NUMPRIMES*sizeof(unsigned long));
   
   s = 3;
   for (unsigned long i = 0; i < NUMPRIMES; i++)
   {
      prime[i] = s;
      s = long_nextprime(s);
   }
   for (unsigned long i = 0; i < NUMPRIMES; i++)
   {
      s = (T%prime[i]);
      if (s == 0) s = prime[i];
      mod[i] = prime[i]-s;
   }
   
   
   mpz_init(temp);
   mpz_init(temp2);
   mpz_init(D);
   
   //Tmod72 = T%72;
   
   unsigned long Ttab = 0;
   
   while (T < STOP)
   {
      memset(table, 0, TABLESIZE);
      for (unsigned long offset = 0; offset < TABLESIZE; offset+=CACHEBLOCK)
      {
        for (unsigned long i = 0; i < NUMPRIMES; i++)
        {
          s = mod[i];
          p = prime[i];
          current = table + offset;
          for ( ; s < CACHEBLOCK; s += p) current[s] = 1;
          s -= CACHEBLOCK;
          mod[i] = s;
        }
      }
      
      //oldT = T;
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
      long_precompute_inverse2(&Tinv_hi, &Tinv_lo, T);
      
      //Tmod72 += (T-oldT);
      //while (Tmod72 >= 72) Tmod72-=72;

      while (Ttab < TABLESIZE)
      {
            
      //if (OK[Tmod72])
      //{
         root1 = long_cuberootmod(&cuberoot1, k1, T);
         if (root1)
         {
            if (cuberoot1 != 1)
            {
               root2 = long_mulmod_precomp2(root1, cuberoot1, Tinv_hi, Tinv_lo, T);
               root3 = long_mulmod_precomp2(root2, cuberoot1, Tinv_hi, Tinv_lo, T);
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
         root1 = long_cuberootmod(&cuberoot1, k2, T);
         if (root1)
         {
            if (cuberoot1 != 1)
            {
               root2 = long_mulmod_precomp2(root1, cuberoot1, Tinv_hi, Tinv_lo, T);
               root3 = long_mulmod_precomp2(root2, cuberoot1, Tinv_hi, Tinv_lo, T);
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
         root1 = long_cuberootmod(&cuberoot1, k3, T);
         if (root1)
         {
            if (cuberoot1 != 1)
            {
               root2 = long_mulmod_precomp2(root1, cuberoot1, Tinv_hi, Tinv_lo, T);
               root3 = long_mulmod_precomp2(root2, cuberoot1, Tinv_hi, Tinv_lo, T);
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
      //}
                    
      //oldT = T;
      do
      {
         do 
         {
            T += 2;
            Ttab += 2;
         } while (table[Ttab] && (Ttab < TABLESIZE));
         long_precompute_inverse2(&Tinv_hi, &Tinv_lo, T);
      } while (!long_isprime_precomp2(T, Tinv_hi, Tinv_lo) && (Ttab < TABLESIZE));
      
      //Tmod72 += (T-oldT);
      //while (Tmod72 >= 72) Tmod72-=72;
      
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
