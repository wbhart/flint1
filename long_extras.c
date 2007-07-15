/****************************************************************************

   long_extras.c: extra functions for longs and unsigned longs

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/types.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"

/*  
   Returns a pseudorandom integer in the range [0, limit)
   todo: get rid of divisions
*/

unsigned long long_randint(unsigned long limit) 
{
#if FLINTBITS == 32
    static u_int64_t randval = 4035456057U;
    randval = ((u_int64_t)randval*(u_int64_t)1025416097U+(u_int64_t)286824430U)%(u_int64_t)4294967311U;
    
    return (unsigned long)randval%limit;
#else
    static unsigned long randval = 4035456057U;
    static unsigned long randval2 = 6748392731U;
    randval = ((unsigned long)randval*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    randval2 = ((unsigned long)randval2*(unsigned long)1647637699U+(unsigned long)286824428U)%(unsigned long)4294967357U;
    
    return (unsigned long)(randval+(randval2<<32))%limit;
#endif
}

/* 
   Computes a double floating point approximate inverse, 
   i.e. 53 bits of 1 / n 
   Requires that n be no more than 53 bits
*/

double long_precompute_inverse(unsigned long n)
{
   return (double) 1 / (double) n;
}

/* 
    Returns a % n given a precomputed approx inverse ninv
    Operation is *unsigned*
    Requires that n be no more than 53 bits and _a_ be less than n^2
*/

unsigned long long_mod_precomp(unsigned long a, unsigned long n, double ninv)
{
   if (a < n) return a;
   unsigned long quot = (unsigned long) ((double) a * ninv);
   unsigned long rem = a - quot*n;
   if (rem >= n) return rem - n;
   else return rem;
}

/* 
    Returns a / n given a precomputed approx inverse ninv
    Operation is *unsigned*
    Requires that n be no more than 63 bits but there are no
    restrictions on _a_
*/

unsigned long long_div63_precomp(unsigned long a, unsigned long n, double ninv)
{
   if (a < n) return 0;
   unsigned long quot = (unsigned long) ((double) a * ninv);
   long rem = a - quot*n;
   if (rem < (long)(-n)) quot -= (unsigned long) ((double) (-rem) * ninv);
   else if (rem >= (long) n) quot += (unsigned long) ((double) rem * ninv);
   else if (rem < 0L) return quot - 1;
   else return quot;
   rem = a - quot*n;
   if (rem >= (long) n) return quot + 1;
   else if (rem < 0L) return quot - 1;
   else return quot;
}

/* 
    Returns a % n given a precomputed approx inverse ninv
    Operation is *unsigned*
    Requires that n be no more than 63 bits but there are no
    restrictions on _a_
*/

unsigned long long_mod63_precomp(unsigned long a, unsigned long n, double ninv)
{
   if (a < n) return a;
   unsigned long quot = (unsigned long) ((double) a * ninv);
   long rem = a - quot*n;
   if (rem < (long)(-n)) quot -= (unsigned long) ((double) (-rem) * ninv);
   else if (rem >= (long) n) quot += (unsigned long) ((double) rem * ninv);
   else if (rem < 0L) return rem + n;
   else return rem;
   rem = a - quot*n;
   if (rem >= (long) n) return rem - n;
   else if (rem < 0L) return rem + n;
   else return rem;
}

/*
   Computes a_hi a_lo mod n given an approximate inverse
   Assumes n is no more than 63 bits, but there are no
   restrictions on a_hi or a_lo
   Operation is unsigned.
*/

unsigned long long_mod2_precomp(unsigned long a_hi, unsigned long a_lo, 
                                             unsigned long n, double ninv)
{
   unsigned long t1;
   unsigned long norm, q, r, orig_n;
   
   if (a_hi >= n) 
   {
      if (((n>>32) == 0) && (a_hi >= n*n)) a_hi = a_hi%n;
      else a_hi = long_mod_precomp(a_hi, n, ninv);
   }
   
#if UDIV_NEEDS_NORMALIZATION
   count_lead_zeros(norm, n);
   udiv_qrnnd(q, r, (a_hi<<norm) + (a_lo>>(FLINT_BITS-norm)), a_lo<<norm, n<<norm);
#else
   udiv_qrnnd(q, r, a_hi, a_lo, n);
#endif
   
   return r;
}

/* 
   Computes a*b mod n, given a precomputed inverse ninv
   Assumes a an b are both in [0,n). 
   Requires a*b to be no more than n^2
   Requires that n be no more than 53 bits
*/

unsigned long long_mulmod_precomp(unsigned long a, unsigned long b, 
                                         unsigned long n, double ninv)
{
   unsigned long quot = (unsigned long) ((double) a * (double) b * ninv);
   long rem = a*b - quot*n;
   if (rem < 0) 
   {
      rem += n;
      if (rem < 0) return rem + n;
   } else if (rem >= n) return rem - n;
   return rem;
}

/*
   Returns a^exp modulo n
   Assumes a is reduced mod n
   Requires that n be no more than 53 bits
*/

unsigned long long_powmod0(unsigned long a, long exp, unsigned long n)
{
   double ninv = long_precompute_inverse(n);
   
   unsigned long x, y;
   
   unsigned long e;

   if (exp < 0)
      e = (unsigned long) -exp;
   else
      e = exp;
   
   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = long_mulmod_precomp(x, y, n, ninv);
      y = long_mulmod_precomp(y, y, n, ninv);
      e = e >> 1;
   }

   if (exp < 0) x = long_invert(x, n);

   return x;
} 

/*
   Returns a^exp modulo n given a precomputed inverse
   Assumes a is reduced mod n
   Requires that n be no more than 53 bits
*/

unsigned long long_powmod_precomp0(unsigned long a, long exp, 
                                     unsigned long n, double ninv)
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
      if (e & 1) x = long_mulmod_precomp(x, y, n, ninv);
      y = long_mulmod_precomp(y, y, n, ninv);
      e = e >> 1;
   }

   if (exp < 0) x = long_invert(x, n);

   return x;
} 

/* 
   Computes the Jacobi symbol of _a_ modulo p
   Assumes p is a prime of no more than 53 bits and that _a_
   is reduced modulo p
*/

int long_jacobi_precomp(unsigned long a, unsigned long p, double pinv)
{
   if (a == 0) return 0;  
   if (long_powmod_precomp0(a, (p-1)/2, p, pinv) == p-1) return -1;
   else return 1;                                            
}
                      
/* 
   Computes a square root of _a_ modulo p.
   Assumes p is a prime of no more than 53 bits,
   that _a_ is reduced modulo p. 
   Returns 0 if _a_ is a quadratic non-residue modulo p.
*/
unsigned long long_sqrtmod0(unsigned long a, unsigned long p) 
{
     unsigned int r, k, m;
     unsigned long p1, b, g, bpow, gpow, res;
     double pinv;
         
     if ((a==0) || (a==1)) 
     {
        return a;
     }
     
     pinv = long_precompute_inverse(p);
     
     if (long_jacobi_precomp(a, p, pinv) == -1) return 0;
     
     if ((p&3)==3)
     {
        return long_powmod_precomp0(a, (p+1)/4, p, pinv);
     }
     
     r = 0;
     p1 = p-1;
     
     do {
        p1>>=1UL; 
        r++;
     } while ((p1&1UL) == 0);
 
     b = long_powmod_precomp0(a, p1, p, pinv);
     
     for (k=2UL; ;k++)
     {
         if (long_jacobi_precomp(k, p, pinv) == -1) break;
     }
     
     g = long_powmod_precomp0(k, p1, p, pinv);
     res = long_powmod_precomp0(a, (p1+1)/2, p, pinv);
     if (b == 1UL) 
     {
        return res;
     }
        
     while (b != 1)
     {
           bpow = b;
           for (m = 1; (m <= r-1) && (bpow != 1); m++)
           {
               bpow = long_mulmod_precomp(bpow, bpow, p, pinv);
           }
           gpow = g;
           for (int i = 1; i < r-m; i++)
           {
               gpow = long_mulmod_precomp(gpow, gpow, p, pinv);
           }
           res = long_mulmod_precomp(res, gpow, p, pinv);
           gpow = long_mulmod_precomp(gpow, gpow, p, pinv);
           b = long_mulmod_precomp(b, gpow, p, pinv);
           gpow = g;
           r = m;
     }
     
     return res;
}

/* 
   Computes a 2 limb approximate inverse, i.e. 2 limbs of 2^B / n
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
   if (t4 || (!t4 && (t3 >= n))) 
   {
      t3-=n;
      if (t3 >= n) return t3-n;
   } 
   return t3;    
}
 
/* 
   Computes a*b mod n, given a precomputed inverse ninv_quot ninv_rem
   Assumes a an b are both in [0,n). There is no restriction on a*b, 
   i.e. it can be two limbs
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
   Requires that n be no more than 63 bits if exp < 0
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
   Requires that n be no more than 63 bits if exp < 0
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
   Computes the Jacobi symbol of _a_ modulo p
   Assumes p is a prime and that _a_ is reduced modulo p
*/

int long_jacobi_precomp2(unsigned long a, unsigned long p, 
                            unsigned long pinv_hi, unsigned long pinv_lo)
{
   if (a == 0) return 0;  
   if (long_powmod_precomp2(a, (p-1)/2, p, pinv_hi, pinv_lo) == p-1) return -1;
   else return 1;                                            
}

/* 
   Computes a square root of _a_ modulo p.
   Assumes p is a prime and that _a_ is reduced modulo p. 
   Returns 0 if _a_ is a quadratic non-residue modulo p.
*/
unsigned long long_sqrtmod(unsigned long a, unsigned long p) 
{
     unsigned int r, k, m;
     unsigned long p1, b, g, bpow, gpow, res;
     unsigned long pinv_hi; 
     unsigned long pinv_lo;
         
     if ((a==0) || (a==1)) 
     {
        return a;
     }
     
     long_precompute_inverse2(&pinv_hi, &pinv_lo, p);
     
     if (long_jacobi_precomp2(a, p, pinv_hi, pinv_lo) == -1) return 0;
     
     if ((p&3)==3)
     {
        return long_powmod_precomp2(a, (p+1)/4, p, pinv_hi, pinv_lo);
     }
     
     r = 0;
     p1 = p-1;
     
     do {
        p1>>=1UL; 
        r++;
     } while ((p1&1UL) == 0);
 
     b = long_powmod_precomp2(a, p1, p, pinv_hi, pinv_lo);
     
     for (k=2UL; ;k++)
     {
         if (long_jacobi_precomp2(k, p, pinv_hi, pinv_lo) == -1) break;
     }
     
     g = long_powmod_precomp2(k, p1, p, pinv_hi, pinv_lo);
     res = long_powmod_precomp2(a, (p1+1)/2, p, pinv_hi, pinv_lo);
     if (b == 1UL) 
     {
        return res;
     }
        
     while (b != 1)
     {
           bpow = b;
           for (m = 1; (m <= r-1) && (bpow != 1); m++)
           {
               bpow = long_mulmod_precomp2(bpow, bpow, p, pinv_hi, pinv_lo);
           }
           gpow = g;
           for (int i = 1; i < r-m; i++)
           {
               gpow = long_mulmod_precomp2(gpow, gpow, p, pinv_hi, pinv_lo);
           }
           res = long_mulmod_precomp2(res, gpow, p, pinv_hi, pinv_lo);
           gpow = long_mulmod_precomp2(gpow, gpow, p, pinv_hi, pinv_lo);
           b = long_mulmod_precomp2(b, gpow, p, pinv_hi, pinv_lo);
           gpow = g;
           r = m;
     }
     
     return res;
}

/* 
   Computes a cube root of _a_ mod p for a prime p and returns a cube 
   root of unity if the cube roots of _a_ are distinct else the cube 
   root is set to 1
   If _a_ is not a cube modulo p then 0 is returned
   This function assumes _a_ is not 0 and that _a_ is reduced modulo p
*/

unsigned long long_cuberootmod(unsigned long * cuberoot1, unsigned long a, 
       unsigned long p)
{
   unsigned long x;
   unsigned long pinv_hi; 
   unsigned long pinv_lo;
    
   long_precompute_inverse2(&pinv_hi, &pinv_lo, p);
   
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
      if(m>=r) return(0);
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
   Tests whether n is an a-Strong Pseudo Prime
   Assumes d is set to the largest odd factor of n-1
*/
static inline
int SPRP(unsigned long a, unsigned long d, unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo)
{
      unsigned long t = d;
      unsigned long y;
      
      y = long_powmod_precomp2(a, t , n, ninv_hi, ninv_lo);
      while ((t != n-1) && (y != 1) && (y != n-1))
      {
         y = long_mulmod_precomp2(y, y, n, ninv_hi, ninv_lo);
         t <<= 1;
      }
      if ((y != n-1) && ((t&1) == 0)) return 0;
      return 1;
}

/* 
    Miller-Rabin primality test. If reps is set to 5 a couple of 
    pseudoprimes on average will pass the test out of each 10^11 tests. 
    Every increase of reps by 1 decreases the chance or composites 
    passing by a factor of 4. 
*/
     
int long_miller_rabin_precomp2(unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo, unsigned long reps)
{
   unsigned long d = n-1, a, t, y;
   
   do {
      d>>=1UL; 
   } while ((d&1UL) == 0);
      
   for (unsigned long i = 0; i < reps; i++)
   {
      a = long_randint(n-2)+1UL;
      t = d;
      y = long_powmod_precomp2(a, t , n, ninv_hi, ninv_lo);
      while ((t != n-1) && (y != 1UL) && (y != n-1))
      {
         y = long_mulmod_precomp2(y, y, n, ninv_hi, ninv_lo);
         t <<= 1UL;
      }
      if ((y != n-1) && ((t&1UL) == 0UL)) return 0;
   }
   return 1;
}

/* 
   For n < 10^16 this is a deterministic prime test. 
   For n > 10^16 then it is a probabalistic test at present.
   Todo: use the table here: http://oldweb.cecm.sfu.ca/pseudoprime/
   to make this into an unconditional primality test.
   This test is intended to be run after checking for divisibility by
   primes up to 257 say.
*/
int long_isprime_precomp2(unsigned long n, unsigned long ninv_hi, unsigned long ninv_lo)
{
   unsigned long d = n-1;
   
   do {
      d>>=1; 
   } while ((d&1) == 0);
   
   if (n < 9080191UL) 
   { 
      if (SPRP(31UL, d, n, ninv_hi, ninv_lo) && SPRP(73UL, d, n, ninv_hi, ninv_lo)) return 1;
      else return 0;
   }
   if (n < 4759123141UL)
   {
      if (SPRP(2UL, d, n, ninv_hi, ninv_lo) && SPRP(7UL, d, n, ninv_hi, ninv_lo) && SPRP(61UL, d, n, ninv_hi, ninv_lo)) return 1;
      else return 0;
   }
   if (n < 1122004669633UL)
   {
      if (SPRP(2UL, d, n, ninv_hi, ninv_lo) && SPRP(13UL, d, n, ninv_hi, ninv_lo) && SPRP(23UL, d, n, ninv_hi, ninv_lo) && SPRP(1662803UL, d, n, ninv_hi, ninv_lo)) 
         if (n != 46856248255981UL) return 1;
      return 0;
   }
   if (n < 10000000000000000UL)
   {
      if (SPRP(2UL, d, n, ninv_hi, ninv_lo) && SPRP(3UL, d, n, ninv_hi, ninv_lo) && SPRP(7UL, d, n, ninv_hi, ninv_lo) && SPRP(61UL, d, n, ninv_hi, ninv_lo) && SPRP(24251UL, d, n, ninv_hi, ninv_lo)) 
         if (n != 46856248255981UL) return 1;
      return 0;
   }
   return long_miller_rabin_precomp2(n, ninv_hi, ninv_lo, 6);  
}

/* 
   For n < 10^16 this is a deterministic prime test. 
   For n > 10^16 then it is a probabalistic test at present.
   Todo: use the table here: http://oldweb.cecm.sfu.ca/pseudoprime/
   to make this into an unconditional primality test.
   This test is intended to be run after checking for divisibility by
   primes up to 257 say.
*/
int long_isprime(unsigned long n)
{
   unsigned long ninv_hi, ninv_lo;
   unsigned long d = n-1;
   
   long_precompute_inverse2(&ninv_hi, &ninv_lo, n);

   do {
      d>>=1; 
   } while ((d&1) == 0);
   
   if (n < 9080191UL) 
   { 
      if (SPRP(31UL, d, n, ninv_hi, ninv_lo) && SPRP(73UL, d, n, ninv_hi, ninv_lo)) return 1;
      else return 0;
   }
   if (n < 4759123141UL)
   {
      if (SPRP(2UL, d, n, ninv_hi, ninv_lo) && SPRP(7UL, d, n, ninv_hi, ninv_lo) && SPRP(61UL, d, n, ninv_hi, ninv_lo)) return 1;
      else return 0;
   }
   if (n < 1122004669633UL)
   {
      if (SPRP(2UL, d, n, ninv_hi, ninv_lo) && SPRP(13UL, d, n, ninv_hi, ninv_lo) && SPRP(23UL, d, n, ninv_hi, ninv_lo) && SPRP(1662803UL, d, n, ninv_hi, ninv_lo)) 
         if (n != 46856248255981UL) return 1;
      return 0;
   }
   if (n < 10000000000000000UL)
   {
      if (SPRP(2UL, d, n, ninv_hi, ninv_lo) && SPRP(3UL, d, n, ninv_hi, ninv_lo) && SPRP(7UL, d, n, ninv_hi, ninv_lo) && SPRP(61UL, d, n, ninv_hi, ninv_lo) && SPRP(24251UL, d, n, ninv_hi, ninv_lo)) 
         if (n != 46856248255981UL) return 1;
      return 0;
   }
   return long_miller_rabin_precomp2(n, ninv_hi, ninv_lo, 6);  
}

unsigned int nextmod30[] = 
{
   1, 6, 5, 4, 3, 2, 1, 4, 3, 2, 1, 2, 1, 4, 3, 2, 1, 2, 1,
   4, 3, 2, 1, 6, 5, 4, 3, 2, 1, 2
};

unsigned int nextindex[] = 
{
   1, 7, 7, 7, 7, 7, 7, 11, 11, 11, 11, 13, 13, 17, 17, 17, 17, 19, 19,
   23, 23, 23, 23, 29, 29, 29, 29, 29, 29, 1
};

unsigned int primes[] =
{
   2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
   191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,
   281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,
   389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,
   491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,
   607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,
   719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,
   829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,
   953,967,971,977,983,991,997
};

#define NEXTPRIME_PRIMES 54
#define NUMBER_OF_PRIMES 168

/* 
    Returns the next prime after n 
    Assumes the result will fit in an unsigned long
*/

unsigned long long_nextprime(unsigned long n)
{
   if (n < 7) 
   {
      if (n<2) return 2;
      n++;
      n|=1;
      return n;  
   }
   
   unsigned long index = n%30;
   n+=nextmod30[index];
   index = nextindex[index];
         
   if (n <= primes[NEXTPRIME_PRIMES-1])
   {
      if (n == 7) return 7;
      if (n == 11) return 11;
      if (n == 13) return 13;
      
      while (((n%7)==0)||((n%11)==0)||((n%13)==0))
      {
         n += nextmod30[index];
         index = nextindex[index];
      }
      return n;
   }
    
   unsigned int * moduli = (unsigned int *) flint_stack_alloc_bytes(NEXTPRIME_PRIMES * sizeof(unsigned int));

   for (unsigned int i = 3; i < NEXTPRIME_PRIMES; i++)
      moduli[i] = (n % primes[i]);
      
   while (1) 
   {
      unsigned int composite = 0;

      unsigned int diff, acc, pr;;
      
      diff = nextmod30[index];
      
      /* First check residues */
      for (unsigned int i = 3; i < NEXTPRIME_PRIMES; i++)
	  {
	     composite |= (moduli[i] == 0);
	     acc = moduli[i] + diff;
	     pr = primes[i];
	     moduli[i] = acc >= pr ? acc - pr : acc;
	   }
       if (composite)
       {
	      n += diff;
          index = nextindex[index];
          continue;
       }
       
       /* Miller-Rabin test */
      if (long_isprime(n)) break;
      else
      {
         n += diff;
         index = nextindex[index];
      }   
   }
   
   flint_stack_release(); 
   
   return n;
}

/* 
    returns the inverse of a modulo p
*/

unsigned long long_invert(unsigned long a, unsigned long p)
{
   if (a == 0) return 0;
   
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

/*
   Return 0 <= a < n1*n2 such that a mod n1 = x1 and a mod n2 = x2
   Assumes gcd(n1, n2) = 1 and that n1*n2 is at most 62 bits
   Assumes x1 is reduced modulo n1 and x2 is reduced modulo n2
*/

unsigned long long_CRT(unsigned long x1, unsigned long x2, 
                       unsigned long n1, unsigned long n2)
{
     unsigned long n, res, ch;
     unsigned long ninv_hi, ninv_lo;
     
     n = n1*n2;
     if (n == 1) return 0;
     long_precompute_inverse2(&ninv_hi, &ninv_lo, n);
     
     res = long_invert(n2,n1);
     res = long_mulmod_precomp2(res, n2, n, ninv_hi, ninv_lo);
     res = long_mulmod_precomp2(res, x1, n, ninv_hi, ninv_lo);
     
     ch = long_invert(n1,n2);
     ch = long_mulmod_precomp2(ch, n1, n, ninv_hi, ninv_lo);
     ch = long_mulmod_precomp2(ch, x2, n, ninv_hi, ninv_lo);
     
     res = res+ch;
     if (res >= n) return res - n;
     else return res;
}

#define SQFREE_TF_PRIMES_LIMIT 168
#define SQFREE_TF_CUTOFF 1000000

int long_issquarefree_trial(unsigned long n)
{
   unsigned long quot, rem;
   
   if ((n&1) == 0)
   {
      if ((n&3) == 0) return 0;
      else n = (n>>1);
   }
   for (unsigned long i = 1; (i < SQFREE_TF_PRIMES_LIMIT) && (primes[i]*primes[i] <= n); i++)
   {
      quot = n/primes[i];
      rem = n - quot*primes[i];
      if (rem == 0) 
      { 
         if ((quot % primes[i]) == 0) return 0;
         else n = quot;
      }
   }
   return 1;
}

/*
   Tests if n is squarefree or not
   Currently only works for numbers up to 65535
*/

int long_issquarefree(unsigned long n)
{
   if (n < SQFREE_TF_CUTOFF) return long_issquarefree_trial(n);
   else 
   {
      printf("Not implemented yet!\n");
      abort();
   }
}

/*
   Removes the highest power of p possible from n and 
   returns the exponent to which it appeared in n
   n can be up to 63 bits
*/

int long_remove63_precomp(unsigned long * n, unsigned long p, double pinv)
{
   unsigned long quot, rem;
   int exp = 0;
   
   quot = long_div63_precomp(*n, p, pinv);
   rem = (*n) - quot*p;
   printf("n = %ld, quot = %ld, rem = %ld\n", *n, quot, rem);
   while (rem == 0); 
   {
      exp++;
      (*n) = quot;
      quot = long_div63_precomp(*n, p, pinv);
      rem = (*n) - quot*p;
   } 
   return exp;
}
   
/*
   Removes the highest power of p possible from n and 
   returns the exponent to which it appeared in n
*/

int long_remove(unsigned long * n, unsigned long p)
{
   unsigned long exp; 
   int i;
   unsigned long powp[7]; // One more than I calculate is necessary for 64 bits
   unsigned long quot, rem;
   
   if (p == 2)
   {
      count_trail_zeros(exp, *n);
      if (exp)
      {
         *n = ((*n)>>exp);
         return exp;      
      }
   }
   
   powp[0] = p;
   
   for (i = 0; ; i++)
   {
      quot = *n/powp[i];
      rem = *n - quot*powp[i];
      if (rem != 0) break;
      powp[i + 1] = powp[i] * powp[i];
      *n = quot;
   }
   
   exp = (1<<i) - 1;
   
   while (i > 0)
   {
      i--;
      quot = *n/powp[i];
      rem = *n - quot*powp[i];
      if (rem == 0) 
      {
         exp += (1<<i);
         *n = quot;
      }     
   } 
   
   return exp;
}

#define TF_CUTOFF 168
#define TF_FACTORS_IN_LIMB 8 // how many factors above the TF cutoff could
                             // an integer one limb wide have

/*
   Finds all the factors of n by trial factoring up to some limit
   Returns the cofactor after removing these factors
*/

unsigned long long_factor_trial(factor_t * factors, unsigned long n)
{
   int num_factors = 0;
   int exp;
   
   for (unsigned long i = 0; (i < TF_CUTOFF) && (primes[i]*primes[i] <= n); i++)
   {
      exp = long_remove(&n, primes[i]);
      if (exp)
      {
         factors->p[num_factors] = primes[i];
         factors->exp[num_factors] = exp;
         num_factors++;
      }      
   }
       
   factors->num = num_factors;
   
   return n;
}

/*
   Square forms factoring algorithm of Shanks
   Adapted from the (simplified) algorithm as described by 
   Gower and Wagstaff Math of Comp. (Preprint May 2007)
*/

#define SQUFOF_ITERS 50000

unsigned long _long_factor_SQUFOF(unsigned long n)
{
   unsigned long sqroot = long_intsqrt(n);
   unsigned long p = sqroot;
   unsigned long q = n - sqroot*sqroot;
   
   if (q == 0) 
   {
      return sqroot;
   }
   
   unsigned long l = 1 + 2*long_intsqrt(2*p);
   unsigned long l2 = l/2;
   unsigned long iq, pnext;
   unsigned long qarr[50];
   unsigned long qupto = 0;
   unsigned long qlast = 1;
   unsigned long i, j, t, r;
   
   for (i = 0; i < SQUFOF_ITERS; i++)
   {
      iq = (sqroot + p)/q;
      pnext = iq*q - p;
      if (q <= l) 
      {
         if ((q & 1) == 0) 
         {
            qarr[qupto] = q/2;
            qupto++;
            if (qupto >= 50) return 0;
         } else if (q <= l2)
         {
            qarr[qupto] = q;
            qupto++;
            if (qupto >= 50) return 0;
         }
      }
      t = qlast + iq*(p - pnext);
	  qlast = q;
	  q = t;
	  p = pnext;
	  if ((i&1) == 1) continue;
	  if (!long_issquare(q)) continue;
	  r = long_intsqrt(q);
	  if (qupto == 0) break;
	  for (j = 0; j < qupto; j++)	
         if (r == qarr[j]) goto cont;
      break;
      cont: ; if (r == 1) return 0;
   }
   
   if (i == SQUFOF_ITERS) return 0; // taken too long, give up
   
   qlast = r;
   p = p + r*((sqroot - p)/r);
   q = (n - p*p)/qlast;
   
   for (j = 0; j < SQUFOF_ITERS; j++)
   {	
	  iq = (sqroot + p)/q;
	  pnext = iq*q - p;
	  if (p == pnext) break;
	  t = qlast + iq*(p - pnext);
	  qlast = q;
	  q = t;
	  p = pnext;
   }
   
   if ((q & 1) == 0) q /= 2;
   
   return q;
}

/* 
   Factor n using as many rounds of SQUFOF as it takes
   Assumes trial factoring of n has already been done and that
   n is not a prime
*/

unsigned long long_factor_SQUFOF(unsigned long n)
{
   unsigned long factor = _long_factor_SQUFOF(n);
   unsigned long multiplier;
   unsigned long quot, rem, kn;
   unsigned long s1, s2;
   
   if (factor) return factor;
   
   for (unsigned long i = 1; (i < NUMBER_OF_PRIMES) && !factor; i++)
   {
      multiplier = primes[i];
      count_lead_zeros(s1, multiplier);
      s1 = FLINT_BITS - s1;
      count_lead_zeros(s2, n);
      if (s1 > s2) return 0; // kn is more than one limb 
      kn = multiplier*n;
      factor = _long_factor_SQUFOF(kn);
      if (factor) 
      {
         quot = factor/multiplier;
         rem = factor - quot*multiplier;
         if (!rem) factor = quot;
         if ((factor == 1) || (factor == n)) factor = 0;
      }
   }
   return factor; 
}

void insert_factor(factor_t * factors, unsigned long p)
{
   int i = 0;
   
   for (i = 0; i < factors->num; i++)
   {
      if (factors->p[i] == p)
      {
         factors->exp[i]++;
         break;
      }
   }
   if (i == factors->num)
   {
      factors->p[i] = p;
      factors->exp[i] = 1;
      factors->num++;
   }
}

/*
   Find the factors of n
   This function may fail if n is so large that SQUFOF multipliers
   (rarely up to 1000) times n push it over a limb in size
   If factoring fails (very rare), this function returns 0, 
   else it returns 1
*/

int long_factor(factor_t * factors, unsigned long n)
{
   unsigned long cofactor;
   unsigned long factor_arr[TF_FACTORS_IN_LIMB];
   unsigned long cutoff = primes[TF_CUTOFF-1]*primes[TF_CUTOFF-1];
   unsigned long factors_left = 1;
   unsigned long factor;
   
   cofactor = long_factor_trial(factors, n);
      
   if (cofactor != 1)
   {
      factor = factor_arr[0] = cofactor;
      
      while (factors_left > 0)
      {
         factor = factor_arr[factors_left-1]; 
         if ((factor < cutoff) || long_isprime(factor))
         {
            insert_factor(factors, factor);
            factors_left--;
         } else
         {
            factor = factor_arr[factors_left] = long_factor_SQUFOF(factor);
            if (!factor_arr[factors_left]) return 0;
            factor_arr[factors_left-1] /= factor;
            factors_left++;
         }
      }
      return 1;
   } 
} 
