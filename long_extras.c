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
    randval = ((u_int64_t)randval*(u_int64_t)1025416097U+(u_int64_t)286824428U)%(u_int64_t)4294967311U;
    
    return (unsigned long)randval%limit;
#else
    static unsigned long randval = 4035456057U;
    static unsigned long randval2 = 6748392731U;
    randval = ((unsigned long)randval*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    randval2 = ((unsigned long)randval2*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    
    return (unsigned long)(randval+(randval2<<32))%limit;
#endif
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
   Computes the Jacobi symbol of _a_ modulo p
   Assumes p is a prime of no more than 63 bits and that _a_
   is reduced modulo p
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
   Assumes p is a prime of no more than 63 bits,
   that _a_ is reduced modulo p and is a quadratic
   residue modulo p. Returns 0 if a is a quadratic
   non-residue modulo p.
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
   Assumes p is no more than 63 bits 
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
   Assumes n is no more than 63 bits
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
    passing by a factor of 4. Assumes n is no more than 63 bits
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
   7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
   191,193,197,199,211,223,227,229,233,239,241,251
};

#define NUMBER_OF_PRIMES 51

/* 
    Returns the next prime after n 
    Assumes the result will fit in an unsigned long
    Assumes n is at most 63 bits
    N.B: Test currently fails for 63 bits.
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
         
   if (n <= primes[NUMBER_OF_PRIMES-1])
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
    
   unsigned int * moduli = (unsigned int *) flint_stack_alloc_bytes(NUMBER_OF_PRIMES * sizeof(unsigned int));

   for (unsigned int i = 0; i < NUMBER_OF_PRIMES; i++)
      moduli[i] = (n % primes[i]);
      
   while (1) 
   {
      unsigned int composite = 0;

      unsigned int diff, acc, pr;;
      
      diff = nextmod30[index];
      
      /* First check residues */
      for (unsigned int i = 0; i < NUMBER_OF_PRIMES; i++)
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
   Assumes gcd(n1, n2) = 1 and that n1*n2 is at most 63 bits
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




