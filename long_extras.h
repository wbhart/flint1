/******************************************************************************

 long_extras.h
 Header file for long_extras.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef LONGEXTRAS_H
#define LONGEXTRAS_H

#include <math.h>

unsigned long long_randint(unsigned long limit);

void long_precompute_inverse2(unsigned long * ninv_hi, 
                            unsigned long * ninv_lo, unsigned long n);
                            
unsigned long long_mod_precomp2(unsigned long a_hi, unsigned long a_lo, 
         unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo);
         
unsigned long long_mulmod_precomp2(unsigned long a, unsigned long b, unsigned long n,
                        unsigned long ninv_hi, unsigned long ninv_lo);

int long_miller_rabin_precomp2(unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo, unsigned long reps);

int long_isprime_precomp2(unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo);

unsigned long long_pow(unsigned long a, unsigned long exp);

unsigned long long_powmod(unsigned long a, long exp, unsigned long n);

unsigned long long_powmod_precomp2(unsigned long a, long exp, unsigned long n,
                            unsigned long ninv_hi, unsigned long ninv_lo);
                            
unsigned long long_cuberootmod_precomp2(unsigned long * cuberoot1, 
                               unsigned long a, unsigned long pinv_hi, 
                               unsigned long pinv_lo, unsigned long p);

unsigned long long_nextprime(unsigned long n);

unsigned long long_invert(unsigned long a, unsigned long p);

long long_gcd_invert(long* a, long x, long y);

long long_extgcd(long* a, long* b, long x, long y);

unsigned long long_gcd(long x, long y);

static inline int long_issquare(long x)
{
   static int mod64[64] = {1,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0}; 
   static int mod65[65] = {1,1,0,0,1,0,0,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,1,0,0,1};
   static int mod63[63] = {1,1,0,0,1,0,0,1,0,1,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0};
   
   if (x < 0) return 0;
   if (!mod64[x%64]) return 0;
   if (!mod63[x%63]) return 0;
   if (!mod65[x%65]) return 0;
   unsigned long sqroot = (unsigned long) sqrt(x);
   return (x == sqroot*sqroot);
}


#endif






