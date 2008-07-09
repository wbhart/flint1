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

   long_extras.c: extra functions for longs and unsigned longs

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"

/* 
   Generate a random integer in the range [0, limit) 
   If limit == 0, return a random limb
*/

unsigned long z_randint(unsigned long limit) 
{
#if FLINT_BITS == 32
    static uint64_t randval = 4035456057U;
    randval = ((uint64_t)randval*(uint64_t)1025416097U+(uint64_t)286824430U)%(uint64_t)4294967311U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)randval%limit;
#else
    static unsigned long randval = 4035456057U;
    static unsigned long randval2 = 6748392731U;
    randval = ((unsigned long)randval*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    randval2 = ((unsigned long)randval2*(unsigned long)1647637699U+(unsigned long)286824428U)%(unsigned long)4294967357U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)(randval+(randval2<<32))%limit;
#endif
}

/*
   Generate a random integer with up to the given number of bits [0, FLINT_BITS]
*/

unsigned long z_randbits(unsigned long bits)
{
   return z_randint(l_shift(1L, bits));
}

/* 
   Computes a double floating point approximate inverse, 
   i.e. 53 bits of 1 / n 
*/

double z_precompute_inverse(unsigned long n)
{
   return (double) 1 / (double) n;
}

/* 
    Returns a % n given a precomputed approx inverse ninv
    Operation is *unsigned*
    Requires that n be no more than FLINT_D_BITS bits and _a_ be less than n^2
*/

unsigned long z_mod_precomp(unsigned long a, unsigned long n, double ninv)
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
    Requires that n be no more than FLINT_BITS-1 bits but there are no
    restrictions on _a_
*/

unsigned long z_div_64_precomp(unsigned long a, unsigned long n, double ninv)
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
    Requires that n be no more than FLINT_BITS-1 bits but there are no
    restrictions on _a_
*/

unsigned long z_mod_64_precomp(unsigned long a, unsigned long n, double ninv)
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
   Assumes n is no more than FLINT_BITS-1 bits, but there are no
   restrictions on a_hi or a_lo
   Operation is unsigned.
*/

unsigned long z_ll_mod_precomp(unsigned long a_hi, unsigned long a_lo, 
                                             unsigned long n, double ninv)
{
   unsigned long t1;
   unsigned long norm, q, r, orig_n;
   
   if (a_hi >= n) 
   {
      if (((n>>(FLINT_BITS/2)) == 0) && (a_hi >= n*n)) a_hi = a_hi%n;
      else a_hi = z_mod2_precomp(a_hi, n, ninv);
   }
   
#if UDIV_NEEDS_NORMALIZATION
   count_lead_zeros(norm, n);
   udiv_qrnnd(q, r, (a_hi<<norm) + (a_lo>>(FLINT_BITS-norm)), a_lo<<norm, n<<norm);
   r >>= norm;
#else
   udiv_qrnnd(q, r, a_hi, a_lo, n);
#endif
   
   return r;
}

/*
   I don't trust the code below. It is for precomputed inverses of 32 bits.
   It should work without modification up to 32 bits, but does not.

   The code is a brutalisation of code found in GMP. */

#if PREINV32

#define invert_limb32(invxl, xl)                  \
  do {                                          \
    invxl = (~(((unsigned long)xl)<<32))/xl; \
  } while (0)

#define LIMB_HIGHBIT_TO_MASK32(n)                                 \
  ((n) & (1<<31)) ? (~0) : (0))

uint32_t z_mod32_precomp(unsigned long n64, uint32_t d, uint32_t di)				
{									
    uint32_t xh, xl, nh, nl, nmask, nadj, q1;			
    unsigned long x;							

    nh = ((uint32_t) (n64 >> 32));								
    nl = ((uint32_t) n64);							
    nmask = LIMB_HIGHBIT_TO_MASK (nl);				
    nadj = nl + (nmask & d);				
    x = (unsigned long) di * (unsigned long) (nh - nmask);			
    x += (((unsigned long)nh)<<32) + (unsigned long) nadj;
    q1 = ~(x>>32);								
    x = (unsigned long) q1 * (unsigned long) d; 
    x += n64;
    xh = (uint32_t) (x >> 32);								
    xl = (uint32_t) x;							
    xh -= d;
    return xl + ((xh) & d);						
} 

uint32_t z_precompute_inverse32(unsigned long n)
{
   unsigned long norm;
   count_lead_zeros(norm, n);
   uint32_t ninv;
   invert_limb32(ninv, (n<<(norm-32)));
   return ninv;
}

unsigned long z_mulmod32_precomp(unsigned long a, unsigned long b, 
                                         unsigned long n, uint32_t ninv)
{
   unsigned long norm;
   unsigned long prod = a*b;
   count_lead_zeros(norm, n);
   norm -= 32;
   unsigned long res = (unsigned long) z_mod32_precomp(prod<<norm, n<<norm, ninv);
   res >>= norm;
   if (res >= n) res -= n;
   return res;
}

#endif

/* 
   Computes a*b mod n, given a precomputed inverse ninv
   Assumes a an b are both in [0,n). 
   Requires that n be no more than FLINT_D_BITS bits
*/

unsigned long z_mulmod_precomp(unsigned long a, unsigned long b, 
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
   Computes a*b mod n, given a precomputed inverse ninv
   Assumes a an b are both in [0,n). There is no restriction on a*b, 
   i.e. it can be two limbs
*/
unsigned long z_mulmod_64_precomp(unsigned long a, unsigned long b, unsigned long n,
                        double ninv)
{
   unsigned long p1, p2;
   
   umul_ppmm(p2, p1, a, b);
   return z_ll_mod_precomp(p2, p1, n, ninv);
}                       

/*
   Returns a^exp modulo n
   Assumes a is reduced mod n
   Requires that n be no more than FLINT_D_BITS bits
   There are no restrictions on exp, which can also be negative.
*/

unsigned long z_powmod(unsigned long a, long exp, unsigned long n)
{
   double ninv = z_precompute_inverse(n);
   
   unsigned long x, y;
   
   unsigned long e;

   if (exp < 0)
      e = (unsigned long) -exp;
   else
      e = exp;
   
   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = z_mulmod_precomp(x, y, n, ninv);
      y = z_mulmod_precomp(y, y, n, ninv);
      e = e >> 1;
   }

   if (exp < 0) x = z_invert(x, n);

   return x;
} 

/*
   Returns a^exp modulo n
   Assumes a is reduced mod n
   Requires that n be no more than FLINT_BITS-1 bits
   There are no restrictions on exp, which can also be negative.
*/

unsigned long z_powmod_64(unsigned long a, long exp, unsigned long n)
{
   double ninv = z_precompute_inverse(n);
   
   unsigned long x, y;
   
   unsigned long e;

   if (exp < 0)
      e = (unsigned long) -exp;
   else
      e = exp;
   
   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = z_mulmod_64_precomp(x, y, n, ninv);
      y = z_mulmod_64_precomp(y, y, n, ninv);
      e = e >> 1;
   }

   if (exp < 0) x = z_invert(x, n);

   return x;
} 

/*
   Returns a^exp modulo n given a precomputed inverse
   Assumes a is reduced mod n
   Requires that n be no more than FLINT_D_BITS bits
   There are no restrictions on exp, which may also be negative
*/

unsigned long z_powmod_precomp(unsigned long a, long exp, 
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
      if (e & 1) x = z_mulmod_precomp(x, y, n, ninv);
      y = z_mulmod_precomp(y, y, n, ninv);
      e = e >> 1;
   }

   if (exp < 0) x = z_invert(x, n);

   return x;
} 

/*
   Returns a^exp modulo n given a precomputed inverse
   Assumes a is reduced mod n
   Requires that n be no more than FLINT_BITS-1 bits
   There are no restrictions on exp, which may also be negative
*/

unsigned long z_powmod_64_precomp(unsigned long a, long exp, 
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
      if (e & 1) x = z_mulmod_64_precomp(x, y, n, ninv);
      y = z_mulmod_64_precomp(y, y, n, ninv);
      e = e >> 1;
   }

   if (exp < 0) x = z_invert(x, n);

   return x;
} 

/* 
   Computes the Jacobi symbol of _a_ modulo p
   Assumes p is a prime of no more than FLINT_BITS-1 bits and that _a_
   is reduced modulo p
*/

int z_jacobi_precomp(unsigned long a, unsigned long p, double pinv)
{
   if (a == 0) return 0;  
   if (z_powmod2_precomp(a, (p-1)/2, p, pinv) == p-1) return -1;
   else return 1;                                            
}
                      
/* 
   Computes a square root of _a_ modulo p.
   Assumes p is a prime of no more than FLINT_BITS-1 bits,
   that _a_ is reduced modulo p. 
   Returns 0 if _a_ is a quadratic non-residue modulo p.
*/
unsigned long z_sqrtmod(unsigned long a, unsigned long p) 
{
     unsigned int r, k, m;
     unsigned long p1, b, g, bpow, gpow, res;
     double pinv;
         
     if ((a==0) || (a==1)) 
     {
        return a;
     }
     
     pinv = z_precompute_inverse(p);
     
     if (z_jacobi_precomp(a, p, pinv) == -1) return 0;
     
     if ((p&3)==3)
     {
        return z_powmod2_precomp(a, (p+1)/4, p, pinv);
     }
     
     r = 0;
     p1 = p-1;
     
     do {
        p1>>=1UL; 
        r++;
     } while ((p1&1UL) == 0);
 
     b = z_powmod2_precomp(a, p1, p, pinv);
     
     for (k=2UL; ;k++)
     {
         if (z_jacobi_precomp(k, p, pinv) == -1) break;
     }
     
     g = z_powmod2_precomp(k, p1, p, pinv);
     res = z_powmod2_precomp(a, (p1+1)/2, p, pinv);
     if (b == 1UL) 
     {
        return res;
     }
        
     while (b != 1)
     {
           bpow = b;
           for (m = 1; (m <= r-1) && (bpow != 1); m++)
           {
               bpow = z_mulmod2_precomp(bpow, bpow, p, pinv);
           }
           gpow = g;
           for (int i = 1; i < r-m; i++)
           {
               gpow = z_mulmod2_precomp(gpow, gpow, p, pinv);
           }
           res = z_mulmod2_precomp(res, gpow, p, pinv);
           gpow = z_mulmod2_precomp(gpow, gpow, p, pinv);
           b = z_mulmod2_precomp(b, gpow, p, pinv);
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
   This function assumes _a_ is reduced modulo p
   Requires p be no more than FLINT_BITS-1 bits
*/

unsigned long z_cuberootmod(unsigned long * cuberoot1, unsigned long a, 
       unsigned long p)
{
   unsigned long x;
   double pinv; 
    
   if (a == 0) return 0;

   pinv = z_precompute_inverse(p);
   
   if ((p % 3) == 2)
   {
      *cuberoot1 = 1;
      return z_powmod2_precomp(a, 2*((p+1)/3)-1, p, pinv);
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
   
   x = z_powmod2_precomp(a, (q-l)/3, p, pinv);
   temp = z_powmod2_precomp(a, l, p, pinv);
   temp2 = z_powmod2_precomp(x, 3UL, p, pinv);
   b = z_mulmod2_precomp(temp, temp2, p, pinv);
   if (l == 2) x = z_mulmod2_precomp(a, x, p, pinv);
      
   while(z_powmod2_precomp(n, (p-1)/3, p, pinv)==1) n++;
   
   z = z_powmod2_precomp(n, q, p, pinv);
   y = z;
   r = e;
   
   while (b!=1)
   {
      s = z_powmod2_precomp(b, 3UL, p, pinv);
      m = 1;
      while(s!=1) 
      {
         s = z_powmod2_precomp(s, 3UL, p, pinv);
         m++;
      }
      if(m>=r) return(0);
      t = z_powmod2_precomp(y, z_pow(3UL, r-m-1UL), p, pinv);
      y = z_powmod2_precomp(t, 3UL, p, pinv);
      r = m;
      x = z_mulmod2_precomp(t, x, p, pinv);
      b = z_mulmod2_precomp(y, b, p, pinv);
   }
   
   if (r==1) *cuberoot1 = y;
   else *cuberoot1 = z_powmod2_precomp(y, z_pow(3UL, r-1), p, pinv);
   if (l==2) return(x);
   else return(z_invert(x, p));
}

/*
    returns a^exp
*/

unsigned long z_pow(unsigned long a, unsigned long exp)
{
   if (exp == 0) return 1;
   if (a == 1) return 1;
   
   unsigned long power = a;

   for (unsigned long i = 1; i < exp; i++)
      power *= a;

   return power;
}

/* 
   Tests whether n is an a-Strong Pseudo Prime
   Assumes d is set to the largest odd factor of n-1
   Assumes n is at most FLINT_D_BITS bits
   Requires _a_ to be reduced mod n
*/
static inline
int SPRP(unsigned long a, unsigned long d, unsigned long n, double ninv)
{
      unsigned long t = d;
      unsigned long y;
      
      y = z_powmod_precomp(a, t , n, ninv);
      while ((t != n-1) && (y != 1) && (y != n-1))
      {
         y = z_mulmod_precomp(y, y, n, ninv);
         t <<= 1;
      }
      if ((y != n-1) && ((t&1) == 0)) return 0;
      return 1;
}

/* 
   Tests whether n is an a-Strong Pseudo Prime
   Assumes d is set to the largest odd factor of n-1
   Assumes n is at most FLINT_BITS-1 bits
   Requires _a_ to be reduced mod n
*/
static inline
int SPRP_64(unsigned long a, unsigned long d, unsigned long n, double ninv)
{
      unsigned long t = d;
      unsigned long y;
      
      y = z_powmod_64_precomp(a, t , n, ninv);
      while ((t != n-1) && (y != 1) && (y != n-1))
      {
         y = z_mulmod_64_precomp(y, y, n, ninv);
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
    
    Requires n be no more than FLINT_BITS-1 bits
*/
     
int z_miller_rabin_precomp(unsigned long n, double ninv, unsigned long reps)
{
   unsigned long d = n-1, a, t, y;
   
   do {
      d>>=1UL; 
   } while ((d&1UL) == 0);
      
   for (unsigned long i = 0; i < reps; i++)
   {
      a = z_randint(n-2)+1UL;
      t = d;
      y = z_powmod2_precomp(a, t , n, ninv);
      while ((t != n-1) && (y != 1UL) && (y != n-1))
      {
         y = z_mulmod2_precomp(y, y, n, ninv);
         t <<= 1UL;
      }
      if ((y != n-1) && ((t&1UL) == 0UL)) return 0;
   }
   return 1;
}

/* 
   This is a deterministic prime test up to 10^16. 
   Todo: use the table here: http://oldweb.cecm.sfu.ca/pseudoprime/
   to make this into an unconditional primality test for larger n 
   This test is intended to be run after checking for divisibility by
   primes up to 257 say.
   Requires n is no more than FLINT_BITS-1 bits
*/
int z_isprime_precomp(unsigned long n, double ninv)
{
   unsigned long d = n-1;
   
   do {
      d>>=1; 
   } while ((d&1) == 0);
   
   if (n < 9080191UL) 
   { 
      if (SPRP(31UL, d, n, ninv) && SPRP(73UL, d, n, ninv)) return 1;
      else return 0;
   }
#if FLINT_BITS == 64
   if (n < 4759123141UL)
   {
      if (SPRP(2UL, d, n, ninv) && SPRP(7UL, d, n, ninv) && SPRP(61UL, d, n, ninv)) return 1;
      else return 0;
   }
   if (n < 1122004669633UL)
   {
      if (SPRP(2UL, d, n, ninv) && SPRP(13UL, d, n, ninv) && SPRP(23UL, d, n, ninv) && SPRP(1662803UL, d, n, ninv)) 
         if (n != 46856248255981UL) return 1;
      return 0;
   }
   if (n < 10000000000000000UL)
   {
      if (SPRP_64(2UL, d, n, ninv) && SPRP_64(3UL, d, n, ninv) && SPRP_64(7UL, d, n, ninv) && SPRP_64(61UL, d, n, ninv) && SPRP_64(24251UL, d, n, ninv)) 
         if (n != 46856248255981UL) return 1;
      return 0;
   }
   return z_miller_rabin_precomp(n, ninv, 6); 
#else
      if (SPRP(2UL, d, n, ninv) && SPRP(7UL, d, n, ninv) && SPRP(61UL, d, n, ninv)) return 1;
      else return 0;
#endif 
}

/* 
   This is a deterministic prime test up to 10^16. 
   Todo: use the table here: http://oldweb.cecm.sfu.ca/pseudoprime/
   to make this into an unconditional primality test for larger n 
   This test is intended to be run after checking for divisibility by
   primes up to 257 say.
   Requires n to be at most FLINT_BITS-1 bits
*/

int z_isprime(unsigned long n)
{
   double ninv;

   ninv = z_precompute_inverse(n);

   return z_isprime_precomp(n, ninv);
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

unsigned long z_nextprime(unsigned long n)
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
      if (z_isprime(n)) break;
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

unsigned long z_invert(unsigned long a, unsigned long p)
{
   if (a == 0) return 0;
   if (a == 1) return 1; // Important to optimise for Newton inversion
   
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
     We ensure a is reduced mod y
*/

long z_gcd_invert(long* a, long x, long y)
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
               t1 = u2; u2 = u1; u1 = t1; u3 = v3;
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

   if (u1 < 0L) u1 += y;
   *a = u1;
   
   return u3;
}

/* 
     returns gcd(x, y) = a*x + b*y.
*/

long z_extgcd(long* a, long* b, long x, long y)
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
               t1 = u2; u2 = u1; u1 = t1; u3 = v3;
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

unsigned long z_gcd(long x, long y)
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
   Assumes gcd(n1, n2) = 1 and that n1*n2 is at most FLINT_BITS-1 bits
   Assumes x1 is reduced modulo n1 and x2 is reduced modulo n2
   Requires n1*n2 to be at most FLINT_BITS-1 bits
*/

unsigned long z_CRT(unsigned long x1, unsigned long n1,  
                        unsigned long x2, unsigned long n2)
{
     unsigned long n, res, ch;
     double ninv;
     
     n = n1*n2;
     if (n == 1) return 0;
     ninv = z_precompute_inverse(n);
     
     res = z_invert(n2,n1);
     res = z_mulmod2_precomp(res, n2, n, ninv);
     res = z_mulmod2_precomp(res, x1, n, ninv);
     
     ch = z_invert(n1,n2);
     ch = z_mulmod2_precomp(ch, n1, n, ninv);
     ch = z_mulmod2_precomp(ch, x2, n, ninv);
     
     res = res+ch;
     if (res >= n) return res - n;
     else return res;
}

#define SQFREE_TF_PRIMES_LIMIT 168
#define SQFREE_TF_CUTOFF 1000000

int z_issquarefree_trial(unsigned long n)
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

int z_issquarefree(unsigned long n)
{
   if (n < SQFREE_TF_CUTOFF) return z_issquarefree_trial(n);
   else 
   {
      printf("Not implemented yet!\n");
      abort();
   }
}

/*
   Removes the highest power of p possible from n and 
   returns the exponent to which it appeared in n
   n can be up to FLINT_BITS-1 bits
*/

int z_remove_precomp(unsigned long * n, unsigned long p, double pinv)
{
   unsigned long quot, rem;
   int exp = 0;
   
   quot = z_div2_precomp(*n, p, pinv);
   rem = (*n) - quot*p;
   while (rem == 0); 
   {
      exp++;
      (*n) = quot;
      quot = z_div2_precomp(*n, p, pinv);
      rem = (*n) - quot*p;
   } 
   return exp;
}

/*
   Removes the highest power of p possible from n and 
   returns the exponent to which it appeared in n
*/

int z_remove(unsigned long * n, unsigned long p)
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

unsigned long z_factor_trial(factor_t * factors, unsigned long n)
{
   int num_factors = 0;
   int exp;
   
   for (unsigned long i = 0; (i < TF_CUTOFF) && (primes[i]*primes[i] <= n); i++)
   {
      exp = z_remove(&n, primes[i]);
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
   Finds prime factors of n by trial factoring wrt some list of precomputed primes
   or until the product of the factors found is greater than the provided limit
   Returns the cofactor after removing these factors and sets prod to the product
   of the removed factors
   Assumes n is not prime
*/

unsigned long z_factor_partial_trial(factor_t * factors, unsigned long * prod, unsigned long n, unsigned long limit)
{
   int num_factors = 0;
   int exp;
   *prod = 1;
   
   for (unsigned long i = 0; (i < TF_CUTOFF) && ((primes[i]-1)*(primes[i]-1) <= n); i++)
   {
      exp = z_remove(&n, primes[i]);
      if (exp)
      {
         factors->p[num_factors] = primes[i];
         factors->exp[num_factors] = exp;
         num_factors++;
		 *prod *= z_pow(primes[i], exp);
		 if (*prod > limit) 
		 {
            factors->num = num_factors;
	        return n;
		 }
      }  
   }
       
   factors->num = num_factors;
   return 0;
}

/*
   Partially factors n into primes. Keeps finding prime factors until the product
   of all the factors found surpasses limit.
   Currently just calls z_factor_partial_trial and returns 0 if this was unable to
   factorise the number sufficiently well, otherwise it returns the cofactor after 
   removing the highest powers of all the prime factors found that it can remove
*/

unsigned long z_factor_partial(factor_t * factors, unsigned long n, unsigned long limit)
{
   unsigned long prod;
   unsigned long cofactor = z_factor_partial_trial(factors, &prod, n, limit);
   if (prod <= limit) return 0;
   else return cofactor;
}

/*
   Square forms factoring algorithm of Shanks
   Adapted from the (simplified) algorithm as described by 
   Gower and Wagstaff Math of Comp. (Preprint May 2007)
*/

#define SQUFOF_ITERS 50000

unsigned long _z_factor_SQUFOF(unsigned long n)
{
   unsigned long sqroot = z_intsqrt(n);
   unsigned long p = sqroot;
   unsigned long q = n - sqroot*sqroot;
   
   if (q == 0) 
   {
      return sqroot;
   }
   
   unsigned long l = 1 + 2*z_intsqrt(2*p);
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
	  if (!z_issquare(q)) continue;
	  r = z_intsqrt(q);
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

unsigned long z_factor_SQUFOF(unsigned long n)
{
   unsigned long factor = _z_factor_SQUFOF(n);
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
      factor = _z_factor_SQUFOF(kn);
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

int z_factor(factor_t * factors, unsigned long n)
{
   unsigned long cofactor;
   unsigned long factor_arr[TF_FACTORS_IN_LIMB];
   unsigned long cutoff = primes[TF_CUTOFF-1]*primes[TF_CUTOFF-1];
   unsigned long factors_left = 1;
   unsigned long factor;
   
   cofactor = z_factor_trial(factors, n);
      
   if (cofactor != 1)
   {
      factor = factor_arr[0] = cofactor;
      
      while (factors_left > 0)
      {
         factor = factor_arr[factors_left-1]; 
         if ((factor < cutoff) || z_isprime(factor))
         {
            insert_factor(factors, factor);
            factors_left--;
         } else
         {
            factor = factor_arr[factors_left] = z_factor_SQUFOF(factor);
            if (!factor_arr[factors_left]) return 0;
            factor_arr[factors_left-1] /= factor;
            factors_left++;
         }
      }
      return 1;
   } 
   return 1;
} 

/*
   Finds the smallest primitive root of the prime p
   
   Initial attempt - could be extremely sub-optimal!
*/

unsigned long z_primitive_root(unsigned long p)
{
   FLINT_ASSERT(p > 2);
   FLINT_ASSERT(z_isprime(p) == 1);
   
   unsigned long res;
   factor_t factors;
   
   if(z_factor(&factors, (p - 1)) == 0)
   {
      return 0;
   }
   
   res = 2;
   
   int i = 0;
   do {
      if(z_powmod(res, (p-1) / factors.p[i], p) == 1) {
         res++;
         i = 0;
      } else {
         i++;
      }
   } while(i != factors.num);
   
   return res;
}

unsigned long z_primitive_root_precomp(unsigned long p, double p_inv)
{
   FLINT_ASSERT(p > 2);
   FLINT_ASSERT(z_isprime(p) == 1);
   
   unsigned long res;
   factor_t factors;
   
   if(z_factor(&factors, (p - 1)) == 0)
   {
      return 0;
   }
   
   res = 2;
   
   int i = 0;
   do {
      if(z_powmod_precomp(res, (p-1) / factors.p[i], p, p_inv) == 1) {
         res++;
         i = 0;
      } else {
         i++;
      }
   } while(i != factors.num);
   
   return res;
}

