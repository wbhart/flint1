/******************************************************************************

 long_extras.h
 Header file for long_extras.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef LONGEXTRAS_H
#define LONGEXTRAS_H

#include <math.h>

typedef struct factor_s
{
   int num;
   unsigned long p[15];
   unsigned long exp[15];
} factor_t;

unsigned long long_randint(unsigned long limit);

double long_precompute_inverse(unsigned long n);

unsigned long long_mod_precomp(unsigned long a, unsigned long n, double ninv);

unsigned long long_div_1_precomp(unsigned long a, unsigned long n, double ninv);

unsigned long long_mod_1_precomp(unsigned long a, unsigned long n, double ninv);

unsigned long long_mod_2_precomp(unsigned long a_hi, unsigned long a_lo, 
                                             unsigned long n, double ninv);

unsigned long long_mulmod_precomp(unsigned long a, unsigned long b, 
                                         unsigned long n, double ninv);
                                         
unsigned long long_mulmod_1_precomp(unsigned long a, unsigned long b, unsigned long n,
                        double ninv);
                                         
unsigned long long_powmod(unsigned long a, long exp, unsigned long n);

unsigned long long_powmod_precomp(unsigned long a, long exp, 
                                     unsigned long n, double ninv);
                                     
int long_jacobi_precomp(unsigned long a, unsigned long p, double pinv);

int long_isprime(unsigned long n);

unsigned long long_nextprime(unsigned long n);

unsigned long long_pow(unsigned long a, unsigned long exp);

unsigned long long_powmod(unsigned long a, long exp, unsigned long n);
                                                    
unsigned long long_sqrtmod(unsigned long a, unsigned long p); 

unsigned long long_cuberootmod(unsigned long * cuberoot1, unsigned long a, 
       unsigned long p);

unsigned long long_cuberootmod(unsigned long * cuberoot1, 
                               unsigned long a, unsigned long p);

unsigned long long_invert(unsigned long a, unsigned long p);

long long_gcd_invert(long* a, long x, long y);

long long_extgcd(long* a, long* b, long x, long y);

unsigned long long_gcd(long x, long y);

static inline unsigned long long_intsqrt(unsigned long n)
{
   return (unsigned long) floor(sqrt((double)n));
}

static inline int long_issquare(long x)
{
   static int mod64[64] = {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0}; 
   static int mod65[65] = {1,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,1};
   static int mod_ui[63] = {1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0};
   
   if (x < 0) return 0;
   if (!mod64[x%64]) return 0;
   if (!mod_ui[x%63]) return 0;
   if (!mod65[x%65]) return 0;
   unsigned long sqroot = (unsigned long) sqrt((double)x);
   return (x == sqroot*sqroot);
}

unsigned long long_CRT(unsigned long x1, unsigned long x2, 
                       unsigned long n1, unsigned long n2);
                       
int long_issquarefree(unsigned long n);

int long_remove_1_precomp(unsigned long * n, unsigned long p, double pinv);

int long_remove(unsigned long * n, unsigned long p);

unsigned long long_factor_trial(factor_t * factors, unsigned long n);

unsigned long long_factor_SQUFOF(unsigned long n);

int long_factor(factor_t * factors, unsigned long n);

#endif






