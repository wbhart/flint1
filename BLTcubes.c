/****************************************************************************

   BLT_cubes.c: Finds solutions to x^3 + y^3 + z^3 = 42
                Based on the algorithm of Brown, Lioen, Te Riele
   
   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "long_extras.h"

#define K 400000
#define k 42
#define eps 0.000015746
#define theta 3.47602664

#define r 1 //not used at present

#define LO 650000000000
#define HI 1000000000000
#define SKIP 50

static inline int inbounds(unsigned long n)
{
   if (((n < HI) && (n > LO)) || ((n > -HI) && (n < -LO)))
   return 1;         
   else return 0;
}

int main()
{
   FILE * file1 = fopen("output.log","w");
   int cubes[72];
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
           if (cubes[(i*i*i+j*j*j-42)%72]) OK[(i+j)%72] = 1;
       }
   }
   
   const long clim = (long) (K/theta/theta);
   const long blim = (long) (K/theta);
   const long alim = K;
   
   mpz_t temp2, nsq, rem;
   mpz_init(temp2);
   mpz_init(nsq);
   mpz_init(rem);
   
   long n, w, v;
   long asq, kbc3;
   
   long winv, gcd, temp;
   
   long z, zsq;
   
   long D;
   
   long nmod72;
   
   for (long c = -clim; c < clim; c++)
   {
      for (long b = -blim; b < blim; b++)
      {
          n = -alim*alim*alim + k*b*b*b + k*k*c*c*c + 3*k*alim*b*c;
          w = b*b+alim*c;
          v = k*c*c+alim*b;
          
          asq = alim*alim;
          kbc3 = 3*k*b*c;
                    
          for (long a = -alim; a < alim; a++)
          {
            if (inbounds(n))
            {
              nmod72 = (n%72);
              if (nmod72 < 0) nmod72+=72;
              
              if ((a%3 == 1) && OK[nmod72]) 
              {
                 // gcd(x, y) = a*x + b*y
                 if ((w & 1 == 1) || (n & 1 == 1))
                 {
                    gcd = long_gcd_invert(&winv, w, n);
                    if ((gcd == 1) || (gcd == -1)) 
                    {
                       z = ((winv*v)%n);
                       if (!z) z+=n;
                       if ((z ^ n) >= 0) z = -z;
              
                       //D = (4*(k-z*z*z)/n - n*n)/3;
                       if (z >= 0) 
                       {
                          mpz_set_ui(temp2, z);
                          mpz_mul_ui(temp2, temp2, z);
                          mpz_mul_ui(temp2, temp2, z);
                          mpz_neg(temp2, temp2);
                       }
                       else 
                       {
                          mpz_set_ui(temp2, -z);
                          mpz_mul_ui(temp2, temp2, -z);
                          mpz_mul_ui(temp2, temp2, -z);
                       }
                       mpz_add_ui(temp2, temp2, k);
                       if (n > 0) mpz_fdiv_qr_ui(temp2, rem, temp2, n);
                       else
                       {
                          mpz_cdiv_qr_ui(temp2, rem, temp2, -n);
                          mpz_neg(temp2, temp2);
                       }
                       if (mpz_sgn(rem) == 0)
                       {
                          mpz_mul_ui(temp2, temp2, 4);
                          if (n > 0) 
                          {
                             mpz_set_ui(nsq, n);
                             mpz_mul_ui(nsq, nsq, n);
                          }                            
                          else
                          {
                             mpz_set_ui(nsq, -n);
                             mpz_mul_ui(nsq, nsq, -n);
                          }
                          mpz_sub(temp2, temp2, nsq);
                          mpz_fdiv_qr_ui(temp2, rem, temp2, 3);
                          if (!mpz_sgn(rem) && mpz_perfect_square_p(temp2))
                          {
                             fprintf(file1,"%ld, %ld, %ld\n", a, b, c);
                             fflush(file1);
                             mpz_clear(temp2);
                             mpz_clear(nsq);
                             mpz_clear(rem);
                             fclose(file1);
                             return 0;
                          }
                       }
                    } 
                 }
              }
            } else
            {
               while (!inbounds(n) && (a < alim))
               {
                  a += SKIP;
                  n = a*a*a + k*b*b*b + k*k*c*c*c - 3*k*a*b*c;
               }
               a-=(SKIP);
               n = a*a*a + k*b*b*b + k*k*c*c*c - 3*k*a*b*c;
               while (!inbounds(n) && (a < alim))
               {
                  a++;
                  n = a*a*a + k*b*b*b + k*k*c*c*c - 3*k*a*b*c;
               }
               a--;
               n = a*a*a + k*b*b*b + k*k*c*c*c - 3*k*a*b*c;
               w = b*b-a*c;
               v = k*c*c-a*b;
               asq = a*a;
            }
              
            n += (3*(asq+a)-kbc3+1);
            w -= c;
            v -= b;
            asq += (2*a+1);                             
          }
          if ((b&1023) == 0) fprintf(file1,"Checkpoint %ld, %ld\n", b, c);
          fflush(file1);
      }
   }
   
   mpz_clear(temp2);
   mpz_clear(nsq);
   mpz_clear(rem);
   fclose(file1);
   return 0;
}

