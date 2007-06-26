/****************************************************************************

   long_extras.c: extra functions for longs and unsigned longs

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"

/*  
   Returns a pseudorandom integer in the range [0, limit)
   limit must be no more than 4294967291
   todo: get rid of divisions, allow limit = 2^32 or higher
*/

unsigned long long_randint(unsigned long limit) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%limit;
}

/* 
   Computes a 2 limb approximate inverse, i.e. 2 limbs of 2^B / n
   Requires that n be no more than 63 bits
*/

void long_precompute_inverse2(unsigned long * ninv_hi, 
                            unsigned long * ninv_lo, unsigned long n)
{
   unsigned long norm;
   unsigned long rem, temp;
   
#if UDIV_NEEDS_NORMALIZATION
   count_lead_zeros(norm, n);

#if DEBUG2
   printf("n = %lu, lead zeroes = %ld, normalised n = %lu\n", n, norm, n<<norm);
#endif
   udiv_qrnnd(*ninv_hi, rem, 1UL<<norm, 0UL, n<<norm);
   udiv_qrnnd(*ndiv_lo, temp, rem<<norm, 0UL, n<<norm);  
#else
   udiv_qrnnd(*ninv_hi, rem, 1UL, 0UL, n);
   udiv_qrnnd(*ninv_lo, temp, rem, 0UL, n);
#endif
}

/* 
    Returns a_hi a_lo % n given precomputed approx inverse ninv_quot ninv_rem
    Assumes the result fits in a single limb, e.g. a_hi a_lo < n^2
    Operation is *unsigned*
    Requires that n be no more than 63 bits
*/

unsigned long long_mod_precomp2(unsigned long a_hi, unsigned long a_lo, 
         unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo)
{
   unsigned long p2, p3, p4;
   unsigned long t1, t2, t3, t4;
   
   // p3 = top limb of (a_hi a_lo) * (ninv_quot ninv_rem) 
   // approx (may be too small by 1)
   p2 = 0UL;
   p3 = a_hi * ninv_hi;
   umul_ppmm(t2, t1, a_lo, ninv_hi);
   umul_ppmm(t4, t3, a_hi, ninv_lo);
   add_ssaaaa(p3, p2, p3, p2, t2, t1);
   add_ssaaaa(p3, p2, p3, p2, t4, t3);
   
   t4 = a_hi;
   t3 = a_lo;
   umul_ppmm(t2, t1, p3, n);
   sub_ddmmss(t4, t3, t4, t3, t2, t1);
   if (t4 || (!t4 && (t3 >= n))) return t3-n;
   else return t3;   
}
 
/* 
   Computes a*b mod n, given a precomputed inverse ninv_quot ninv_rem
   Assumes a an b are both in [0,n). There is no restriction on a*b, 
   i.e. it can be two limbs
   Requires that n be no more than 63 bits
*/
unsigned long long_mulmod_precomp2(unsigned long a, unsigned long b, unsigned long n,
                        unsigned long ninv_hi, unsigned long ninv_lo)
{
   unsigned long p1, p2;
   
   umul_ppmm(p2, p1, a, b);
   return long_mod_precomp2(p2, p1, n, ninv_hi, ninv_lo);
}                       

/*
    returns a^exp
*/

unsigned long long_pow(unsigned long a, unsigned long exp)
{
   if (exp == 0) return 1;
   if (a == 1) return 1;
   
   unsigned long power = a;

   for (unsigned long i = 1; i < exp; i++)
      power *= a;

   return power;
}

/*
   Returns a^exp modulo n
   Assumes a is reduced mod n
   Requires that n be no more than 63 bits
*/

unsigned long long_powmod(unsigned long a, long exp, unsigned long n)
{
   unsigned long x, y;
   unsigned long ninv_hi, ninv_lo;

   unsigned long e;

   if (exp < 0)
      e = (unsigned long) -exp;
   else
      e = exp;

   long_precompute_inverse2(&ninv_hi, &ninv_lo, n);
      
   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = long_mulmod_precomp2(x, y, n, ninv_hi, ninv_lo);
      y = long_mulmod_precomp2(y, y, n, ninv_hi, ninv_lo);
      e = e >> 1;
   }

   if (exp < 0) x = long_invert(x, n);

   return x;
}

/*
   Returns a^exp modulo n given a precomputed inverse
   Assumes a is reduced mod n
   Requires that n be no more than 63 bits
*/

unsigned long long_powmod_precomp2(unsigned long a, long exp, unsigned long n,
                            unsigned long ninv_hi, unsigned long ninv_lo)
{
   unsigned long x, y;
   
   unsigned long e;

   if (exp < 0)
      e = (unsigned long) -exp;
   else
      e = exp;
   
   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = long_mulmod_precomp2(x, y, n, ninv_hi, ninv_lo);
      y = long_mulmod_precomp2(y, y, n, ninv_hi, ninv_lo);
      e = e >> 1;
   }

   if (exp < 0) x = long_invert(x, n);

   return x;
}

/* 
   Currently assumes p = 1 mod 3
*/

unsigned long long_cuberootmod_precomp2(unsigned long * cuberoot1, unsigned long a, 
       unsigned long pinv_hi, unsigned long pinv_lo, unsigned long p)
{
   unsigned long x;
   
   if ((p % 3) == 2)
   {
      *cuberoot1 = 1;
      return long_powmod_precomp2(a, 2*((p+1)/3)-1, p, pinv_hi, pinv_lo);
   }
   
   unsigned long e=0;
   unsigned long q = p-1;
   unsigned long l;
   unsigned long n = 2;
   unsigned long z, y, r, temp, temp2, b, m, s, t;
   
   r = 1;
   
   while ((q%3) == 0)
   {
      q = q/3;
      e++;
   }
   l = q%3;
   
   x = long_powmod_precomp2(a, (q-l)/3, p, pinv_hi, pinv_lo);
   temp = long_powmod_precomp2(a, l, p, pinv_hi, pinv_lo);
   temp2 = long_powmod_precomp2(x, 3UL, p, pinv_hi, pinv_lo);
   b = long_mulmod_precomp2(temp, temp2, p, pinv_hi, pinv_lo);
   if (l == 2) x = long_mulmod_precomp2(a, x, p, pinv_hi, pinv_lo);
      
   while(long_powmod_precomp2(n, (p-1)/3, p, pinv_hi, pinv_lo)==1) n++;
   
   z = long_powmod_precomp2(n, q, p, pinv_hi, pinv_lo);
   y = z;
   r = e;
   
   while (b!=1)
   {
      s = long_powmod_precomp2(b, 3UL, p, pinv_hi, pinv_lo);
      m = 1;
      while(s!=1) 
      {
         s = long_powmod_precomp2(s, 3UL, p, pinv_hi, pinv_lo);
         m++;
      }
      if(m==r) return(0);
      t = long_powmod_precomp2(y, long_pow(3UL, r-m-1UL), p, pinv_hi, pinv_lo);
      y = long_powmod_precomp2(t, 3UL, p, pinv_hi, pinv_lo);
      r = m;
      x = long_mulmod_precomp2(t, x, p, pinv_hi, pinv_lo);
      b = long_mulmod_precomp2(y, b, p, pinv_hi, pinv_lo);
   }
   
   if (r==1) *cuberoot1 = y;
   else *cuberoot1 = long_powmod_precomp2(y, long_pow(3UL, r-1), p, pinv_hi, pinv_lo);
   if (l==2) return(x);
   else return(long_invert(x, p));
}

/* 
    returns the next prime after n (does not check if the result is too big)
*/

unsigned long long_nextprime(unsigned long n)
{
   mpz_t temp;
   mpz_init(temp);
     
   mpz_set_ui(temp,n);
   mpz_nextprime(temp,temp);
   
   unsigned long prime = mpz_get_ui(temp);
   
   mpz_clear(temp);
   
   return prime;
}

/* 
    returns the inverse of a modulo p
*/

unsigned long long_invert(unsigned long a, unsigned long p)
{
   long u1=1, u3=a;
   long v1=0, v3=p;
   long t1=0, t3=0;
   long quot;
   
   while (v3)
   {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               t1 = u1; u1 = v1; v1 = t1;
               t3 = u3; u3 = v3; v3 = t3;
            } else 
            {
               t1 = u1 - v1; u1 = v1; v1 = t1;
               t3 = u3 - v3; u3 = v3; v3 = t3;
            }
         } else if (quot < (v3<<1))
         {  
            t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
            t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
         } else
         {
            t1 = u1 - v1*3; u1 = v1; v1 = t1;
            t3 = u3 - v3*3; u3 = v3; v3 = t3;
         }
      } else
      {
         quot=u3/v3;
         t1 = u1 - v1*quot; u1 = v1; v1 = t1;
         t3 = u3 - v3*quot; u3 = v3; v3 = t3;
      }
   } 
   
   if (u1<0) u1+=p;
   
   return u1;
}

/* 
     returns gcd(x, y) = a*x + b*y. If gcd = 1 then a = x^-1 mod y
*/

long long_gcd_invert(long* a, long x, long y)
{
   long u1=1; 
   long u2=0; 
   long t1; 
   long u3, v3;
   long quot, rem;
   long xsign = 0;
   
   if (x < 0) {
      x = -x;
      xsign = 1;
   }

   if (y < 0) {
      y = -y;
   }
   
   u3 = x, v3 = y;

   while (v3) {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               rem = u3;
               u2 = u1; u1 = t1; u3 = v3;
               v3 = rem;
            } else 
            {
               t1 = u2; u2 = u1 - u2; u1 = t1; u3 = v3;
               v3 = quot;
            }
         } else if (quot < (v3<<1))
         {  
            t1 = u2; u2 = u1 - (u2<<1); u1 = t1; u3 = v3;
            v3 = quot-u3;
         } else
         {
            t1 = u2; u2 = u1 - 3*u2; u1 = t1; u3 = v3;
            v3 = quot-(u3<<1);
         }
      } else
      {
         quot=u3/v3;
         rem = u3 - v3*quot;
         t1 = u2; u2 = u1 - quot*u2; u1 = t1; u3 = v3;
         v3 = rem;
      }
   }
   if (xsign)
      u1 = -u1;

   *a = u1;
   
   return u3;
}

/* 
     returns gcd(x, y) = a*x + b*y.
*/

long long_extgcd(long* a, long* b, long x, long y)
{
   long u1=1, v1=0;
   long u2=0, v2=1;
   long t1, t2;
   long u3, v3;
   long quot, rem;
   long xsign = 0;
   long ysign = 0;
   
   if (x < 0) {
      x = -x;
      xsign = 1;
   }

   if (y < 0) {
      y = -y;
      ysign = 1;
   }
   
   u3 = x, v3 = y;

   while (v3) {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            t2 = v2; 
            if (quot < 0)
            { 
               rem = u3;
               u2 = u1; u1 = t1; u3 = v3;
               v2 = v1; v1 = t2; v3 = rem;
            } else 
            {
               t1 = u2; u2 = u1 - u2; u1 = t1; u3 = v3;
               v2 = v1 - v2; v1 = t2; v3 = quot;
            }
         } else if (quot < (v3<<1))
         {  
            t1 = u2; u2 = u1 - (u2<<1); u1 = t1; u3 = v3;
            t2 = v2; v2 = v1 - (v2<<1); v1 = t2; v3 = quot-u3;
         } else
         {
            t1 = u2; u2 = u1 - 3*u2; u1 = t1; u3 = v3;
            t2 = v2; v2 = v1 - 3*v2; v1 = t2; v3 = quot-(u3<<1);
         }
      } else
      {
         quot=u3/v3;
         rem = u3 - v3*quot;
         t1 = u2; u2 = u1 - quot*u2; u1 = t1; u3 = v3;
         t2 = v2; v2 = v1 - quot*v2; v1 = t2; v3 = rem;
      }
   }
   if (xsign)
      u1 = -u1;

   if (ysign)
      v1 = -v1;

   *a = u1;
   *b = v1;
   
   return u3;
}

/* 
     returns gcd(x, y)
*/

unsigned long long_gcd(long x, long y)
{
   if (x < 0) {
      x = -x;
   }

   if (y < 0) {
      y = -y;
   }

   long u3 = x, v3 = y;
   long quot, rem;

   while (v3) {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               rem = u3;
               u3 = v3;
               v3 = rem;
            } else 
            {
               u3 = v3;
               v3 = quot;
            }
         } else if (quot < (v3<<1))
         {  
            u3 = v3;
            v3 = quot-u3;
         } else
         {
            u3 = v3;
            v3 = quot-(u3<<1);
         }
      } else
      {
         rem = u3 % v3;
         u3 = v3;
         v3 = rem;
      }
   }

   return u3;
}





