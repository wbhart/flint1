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

unsigned long z_randint(unsigned long limit);

double z_precompute_inverse(unsigned long n);

unsigned long z_mod_precomp(unsigned long a, unsigned long n, double ninv);

unsigned long z_div_64_precomp(unsigned long a, unsigned long n, double ninv);

unsigned long z_mod_64_precomp(unsigned long a, unsigned long n, double ninv);

unsigned long z_ll_mod_precomp(unsigned long a_hi, unsigned long a_lo, 
                                             unsigned long n, double ninv);

unsigned long z_mulmod_precomp(unsigned long a, unsigned long b, 
                                         unsigned long n, double ninv);
                                         
unsigned long z_mulmod_64_precomp(unsigned long a, unsigned long b, unsigned long n,
                        double ninv);
                                         
unsigned long z_powmod(unsigned long a, long exp, unsigned long n);

unsigned long z_powmod_64(unsigned long a, long exp, unsigned long n);

unsigned long z_powmod_precomp(unsigned long a, long exp, 
                                     unsigned long n, double ninv);
                                     
unsigned long z_powmod_64_precomp(unsigned long a, long exp, 
                                     unsigned long n, double ninv);
                                     
#if FLINT_BITS == 64
#define z_div2_precomp z_div_64_precomp
#define z_mod2_precomp z_mod_64_precomp
#define z_mulmod2_precomp z_mulmod_64_precomp
#define z_powmod2 z_powmod_64
#define z_powmod2_precomp z_powmod_64_precomp
#else
#define z_div2_precomp z_div_64_precomp
#define z_mod2_precomp z_mod_precomp
#define z_mulmod2_precomp z_mulmod_precomp
#define z_powmod2 z_powmod
#define z_powmod2_precomp z_powmod_precomp
#endif


int z_jacobi_precomp(unsigned long a, unsigned long p, double pinv);

int z_isprime(unsigned long n);

unsigned long z_nextprime(unsigned long n);

unsigned long z_pow(unsigned long a, unsigned long exp);

unsigned long z_powmod(unsigned long a, long exp, unsigned long n);
                                                    
unsigned long z_sqrtmod(unsigned long a, unsigned long p); 

unsigned long z_cuberootmod(unsigned long * cuberoot1, unsigned long a, 
       unsigned long p);

unsigned long z_cuberootmod(unsigned long * cuberoot1, 
                               unsigned long a, unsigned long p);

unsigned long z_invert(unsigned long a, unsigned long p);

long z_gcd_invert(long* a, long x, long y);

long z_extgcd(long* a, long* b, long x, long y);

unsigned long z_gcd(long x, long y);

static inline unsigned long z_intsqrt(unsigned long n)
{
   return (unsigned long) floor(sqrt((double)n));
}

static inline int z_issquare(long x)
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

unsigned long z_CRT(unsigned long x1, unsigned long x2, 
                       unsigned long n1, unsigned long n2);
                       
int z_issquarefree(unsigned long n);

int z_remove_precomp(unsigned long * n, unsigned long p, double pinv);

int z_remove(unsigned long * n, unsigned long p);

unsigned long z_factor_trial(factor_t * factors, unsigned long n);

unsigned long z_factor_SQUFOF(unsigned long n);

int z_factor(factor_t * factors, unsigned long n);

unsigned long z_primitive_root(unsigned long p);

#endif






