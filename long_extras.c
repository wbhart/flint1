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

   Copyright (C) 2007, 2008 William Hart
   Copyright (C) 2008, Peter Shrimpton

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "QS/tinyQS.h"


#define MIN_HOLF 0xFFFFFUL
#define MAX_HOLF 0x1FFFFFFFFFUL
#define HOLF_MULTIPLIER 480
#define HOLF_ITERS 50000



/* 
   Generate a random integer in the range [0, limit) 
   If limit == 0, return a random limb
*/

const unsigned int z_primes[] =
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
   Generates a random prime of _bits_ bits.
	If proved is o (false) the prime is not proven prime, otherwise it is.
*/

unsigned long z_randprime(unsigned long bits, int proved)
{
	unsigned long limit, rand;
   
	if (bits < 2)
	{
		printf("FLINT Exception: attempt to generate prime < 2!\n");
		abort();
	}
   
	if (bits == FLINT_BITS)
	{
		do
		{
			rand = z_randbits(bits);
      
#if FLINT_BITS == 32
		}  while (rand > 4294967290UL);
#else
	    }  while (rand > 18446744073709551556UL);
#endif
	    rand = z_nextprime(rand, proved);

    } else
    {
	   do
	   {
	      rand = z_randbits(bits);
		  rand = z_nextprime(rand, proved);
	   } while ((rand >> bits) > 0L);
   }
   
   return rand;
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

#if FLINT_BITS == 64

/*
   Computes a 32 bit approximation to 1/xl - 1
	Here xl is thought of as being between 0.5 and (2^32 - 1)/2^32
	This is basically a 33 bit approximation to 1/xl without the leading 1
	The product of this and a 64 bit integer will give the quotient which
	will at worst be 1 too small
*/

#define invert_limb32(invxl, xl)                  \
  do {                                          \
    invxl = ((~(((unsigned long)xl)<<32))/xl); \
  } while (0)

/*
   -1 if high bit set
	0 otherwise
*/

#define LIMB_HIGHBIT_TO_MASK32(n)                                 \
  (((n) & (1<<31)) ? (~0) : (0))


/*uint32_t z_mod32_precomp(unsigned long n64, uint32_t d, uint32_t di)				
{									
    uint32_t xh, xl, nh, nl, nmask, nadj, q1;			
    unsigned long x;							

    nh = ((uint32_t) (n64 >> 32));								
    nl = ((uint32_t) n64);							
    nmask = LIMB_HIGHBIT_TO_MASK32(nl);	// nmask = -1 if high bit of nl is set, else nmask = 0			
    nadj = nl + (nmask & d); // nadj = nl if high bit of nl was 0, else it is nl + d - 2^32	

	 // di may be too small, so that if we compute q = floor( (2^32+di)*nh + (2^32+di)*nl/2^32 )
	 // we will not be more than 1 out (after shifting by 32) even if we only take notice of 
	 // the top bit of nl
	 // thus if that top bit is 0 we can use q = (2^32+di)*nh + nl
	 // if the top bit is 1 then q = (2^32+di)*nh + nl may be 2 too small after shifting by 32
	 // we need to add di/2 to q to ensure that we are no more than one out after shifting by 32
	                          
    x = (unsigned long) di * (unsigned long) (nh - nmask); // x = di * nh if high bit of nl was 0
	                                                        // else x = di * (nh + 1)
    x += (((unsigned long)nh)<<32) + (unsigned long) nadj; // x = di * nh + n64 if high bit of nl was 0
	                                                // else x = di * (nh + 1) + 2^32 * nh - 2^32 + nl + d
    q1 = ~(x>>32); // if high bit of nl is 0 then we want q = (2^32+di)*nh >> 32 
	                // else we get q = (2^32+di)*nh + nl + (di + d) - 2^32
	                // and q1 = 2^32 - q - 1							
    x = (unsigned long) q1 * (unsigned long) d; 
    x += n64; // now x = n64 + 2^32*d - q*d - d
    xh = (uint32_t) (x >> 32);								
    xl = (uint32_t) x;							
    xh -= d; // makes up for the 2^32, and now xh is -1 or 0
    return xl + ((xh) & d); // if xh is -1 we have gone too far, so add d, else we are done
} */

uint32_t z_mod32_precomp(unsigned long n64, uint32_t d, uint32_t di)				
{									
   uint32_t xl, nh;
	unsigned long q, x;
    
   nh = ((uint32_t) (n64 >> 32));								
    
	q = (unsigned long) nh * (unsigned long) di;					
   q += (((unsigned long)nh)<<32);	// Compensate, di is 2^FLINT_BITS too small 	

	// As we only use nh which is only 32 bits, the quotient may be incorrect in 
	// the last 2 bits, i.e. it may be up to 3 too small
	// we adjust for this below

   x = (q>>32)*(unsigned long)d;
	x = n64 - x;

	if (x>>32)							
   {									
	   x -= (unsigned long) d;
		if (x>>32) x -= (unsigned long) d;														
   }									
   
	xl = (uint32_t) x;
	if (xl >= d) xl -= d;	
   
	return xl;							
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
   Computes the Legendre symbol of _a_ modulo p
   Assumes p is a prime of no more than FLINT_BITS-1 bits and that _a_
   is reduced modulo p
*/

int z_legendre_precomp(unsigned long a, unsigned long p, double pinv)
{
   if (a == 0) return 0;  
   if (z_powmod2_precomp(a, (p-1)/2, p, pinv) == p-1) return -1;
   else return 1;                                            
}
   
/*
    Calculates the jacobi symbol (x/y)
    Assumes that gcd(x,y) = 1 and y is odd
*/

int z_jacobi(long x, unsigned long y)
{
	unsigned long b, temp;
	long a;
	int s, exp;
	pre_inv2_t inv;
	a = x;
	b = y;
	s = 1;
	if (a < 0)
	{
		if (((b-1)/2)%2 == 1)
		{
			s = -s;
		}
		a = -a;
	}
	a = a % b;
	if (a == 0)
	{
		if (b == 1) return 1;
		else return 0;
	}
	while (b != 1)
	{
		a = a % b;
	    exp = z_remove(&a, 2);
		if (((exp*(b*b-1))/8)%2 == 1) // we are only interested in values mod 8, 
		{	                          //so overflows don't matter here
			s = -s;
		}
		temp = a;
		a = b;
		b = temp;
		if (((a-1)*(b-1)/4)%2 == 1) // we are only interested in values mod 4, 
		{	                        //so overflows don't matter here
			s = -s;
		}
	}

	return s;
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
     
     if (z_legendre_precomp(a, p, pinv) == -1) return 0;
     
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
         if (z_legendre_precomp(k, p, pinv) == -1) break;
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
   Checks to see if n is  fermat pseudoprime with base i
   Assumes n doesn't divide i
*/

int z_ispseudoprime_fermat(unsigned long const n, unsigned long const i)
{
    /* 
	   By Fermats little thm if i^n-1 is not congruent to 1 then n is not prime 
    */
	return (z_powmod2(i, n-1, n) == 1);
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
   This test is intended to be run after checking for divisibility by
   primes up to 257 say, but returns a correct result if n is odd and n > 2.
   Requires n is no more than FLINT_BITS-1 bits
*/
int z_isprobab_prime_precomp(unsigned long n, double ninv)
{
#if FLINT_BITS == 64
   if (n >= 10000000000000000UL) return z_isprobab_prime_BPSW(n);
#endif

	unsigned long d = n-1;
   
   do {
      d>>=1; 
   } while ((d&1) == 0);
   

    if (n < 2047)
   {
      if (SPRP(2UL, d, n, ninv)) return 1;
      else return 0;
   }
   
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
   
	if (SPRP_64(2UL, d, n, ninv) && SPRP_64(3UL, d, n, ninv) && SPRP_64(7UL, d, n, ninv) && SPRP_64(61UL, d, n, ninv) && SPRP_64(24251UL, d, n, ninv)) 
      if (n != 46856248255981UL) return 1;
   return 0;
   
#else
   if (SPRP(2UL, d, n, ninv) && SPRP(7UL, d, n, ninv) && SPRP(61UL, d, n, ninv)) return 1;
      else return 0;
#endif 
}

/* 
   This is currently only a deterministic primality test if FLINT can factor n - 1.
   Note that this is currently the case if n - 1 has less than FLINT_BITS - 10 bits
	after removal of factors below 1000.
	It is assumed that n is greater than 2 and odd.
*/

int z_isprime_precomp(unsigned long n, double ninv)
{

#if FLINT_BITS == 64
   if (n >= 10000000000000000UL)
	{
        return z_isprime_pocklington(n, -1L);
	}
#endif

	unsigned long d = n-1;
   
   do {
      d>>=1; 
   } while ((d&1) == 0);
   
   if (n < 2047)
   {
      if (SPRP(2UL, d, n, ninv)) return 1;
      else return 0;
   }
   
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
   
	if (SPRP_64(2UL, d, n, ninv) && SPRP_64(3UL, d, n, ninv) && SPRP_64(7UL, d, n, ninv) && SPRP_64(61UL, d, n, ninv) && SPRP_64(24251UL, d, n, ninv)) 
      if (n != 46856248255981UL) return 1;
   return 0;
   
#else
   if (SPRP(2UL, d, n, ninv) && SPRP(7UL, d, n, ninv) && SPRP(61UL, d, n, ninv)) return 1;
      else return 0;
#endif 
}


/*
    This is a very dense representation of the prime indicator
    function for odd primes < 4096.  In that range,
    oddprime_lt_4096 is probably about the fastest primality
    test we can get.

    One might be tempted to do something clever, like pack everything
    according to its value mod 30.  In practice, that either requires
    branching, or some very dirty tricks.  In the end, this is 2-5
    times faster than everything I tried.
*/

static inline int z_oddprime_lt_4096(unsigned long n) {
    static unsigned long oddprime_indicator[] =
    {
        0x816d129a64b4cb6eUL,0x2196820d864a4c32UL,0xa48961205a0434c9UL,
        0x4a2882d129861144UL,0x834992132424030UL, 0x148a48844225064bUL,
        0xb40b4086c304205UL, 0x65048928125108a0UL,0x80124496804c3098UL,
        0xc02104c941124221UL,0x804490000982d32UL, 0x220825b082689681UL,
        0x9004265940a28948UL,0x6900924430434006UL,0x12410da408088210UL,
        0x86122d22400c060UL, 0x110d301821b0484UL, 0x14916022c044a002UL,
        0x92094d204a6400cUL, 0x4ca2100800522094UL,0xa48b081051018200UL,
        0x34c108144309a25UL, 0x2084490880522502UL,0x241140a218003250UL,
        0xa41a00101840128UL, 0x2926000836004512UL,0x10100480c0618283UL,
        0xc20c26584822006dUL,0x4520582024894810UL,0x10c0250219002488UL,
        0x802832ca01140868UL,0x60901300264b0400UL
    };
    unsigned long q = n/2;
    unsigned long x = (q&63);
    return (oddprime_indicator[q/64]&(1UL<<x))>>x;
}



/* 
   This is a deterministic prime test up to 10^16. 
   Requires n to be at most FLINT_BITS-1 bits
	For numbers greater than 10^16 there are no known
	counterexamples to the conjecture that a composite 
	will never be declared prime.
*/

int z_isprobab_prime(unsigned long n)
{
   if (n <= 1UL) return 0;
   if (n == 2UL) return 1;
   if ((n & 1UL) == 0) return 0;
   if (n < 4096UL) return z_oddprime_lt_4096(n);
  
   double ninv;

   ninv = z_precompute_inverse(n);
   
   return z_isprobab_prime_precomp(n, ninv);
}

/* 
   This is currently only a deterministic primality test if FLINT can factor n - 1.
   Note that this is currently the case if n - 1 has less than FLINT_BITS - 10 bits
	after removal of factors below 1000.
*/

int z_isprime(unsigned long n)
{
   if (n <= 1UL) return 0;
   if (n == 2UL) return 1;
   if ((n & 1UL) == 0) return 0;
   if (n < 4096UL) return z_oddprime_lt_4096(n);

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

#define NEXTPRIME_PRIMES 54
#define NUMBER_OF_PRIMES 168

/* 
    Returns the next prime after n 
    Assumes the result will fit in an unsigned long
	 If proved is o (false) the prime is not proven prime, 
	 otherwise it is.
*/

unsigned long z_nextprime(unsigned long n, int proved)
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
         
   if (n <= z_primes[NEXTPRIME_PRIMES-1])
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
      moduli[i] = (n % z_primes[i]);
      
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
	      pr = z_primes[i];
	      moduli[i] = acc >= pr ? acc - pr : acc;
	   }
      if (composite)
      {
	      n += diff;
         index = nextindex[index];
         continue;
      }
       
      if ((!proved && z_isprobab_prime(n)) || (proved && z_isprime(n))) break;
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
   Proves that n is prime using a Pocklington-Lehmer test.
   Returns 0 if composite, 1 if prime and -1 if it failed 
	to prove either way.
*/

int z_isprime_pocklington(unsigned long const n, unsigned long const iterations)
{
	int i, j, k, pass, test, exp;
	unsigned long sqrt, n1, f, factor, prod, temp, b, c;
	factor_t factors;
	pre_inv_t inv;
	
	/* 2 is prime, all other even numbers are not */
	
	if (n%2 == 0){
		if (n == 2) return 1;
		else return 0;
	}
	
	n1 = n - 1;
	f = 1;
	factors.num = 0;
   sqrt = z_intsqrt(n1);

   z_factor_partial(&factors, n1, sqrt, 1);
   inv = z_precompute_inverse(n);

	//for (i = factors.num-1; i >= 0 ; i--)		//Quicker?? keeps exponents down so maybe.
	for (i = 0; i < factors.num; i++)
	{
		/*pass = 0;
		unsigned long exp = n1/factors.p[i];
		
		for (j = 2; j < iterations && pass == 0; j++)
		{

			if (z_fermat_pseudoprime(n, j)){
				if ((z_gcd(z_powmod2(j, exp, n) - 1, n) == 1)){
					pass = 1;
				}
			}
		}
		
		if (pass == 0)
		{
			if (j == iterations){
				//printf("%ld: failed on factor:%ld exp:%ld n1:%ld\n", n, factors.p[i], exp, n1);
				return -2;
			} else {
				return 0;
			}
		}*/
		
		/*
			Another method found in http://www.jstor.org/stable/2005583?seq=1
		*/
		
		pass = 0;
		c = 1;
		unsigned long exp = n1/factors.p[i];
		
		for (j = 2; j < iterations && pass == 0; j++)
		{
			if (z_ispseudoprime_fermat(n, j))
			{
				b = z_submod(z_powmod2(j, exp, n), 1, n);
				if (b != 0)
				{
					c = z_mulmod2_precomp(c, b, n, inv);	
					pass = 1;
				}
			} else return 0;
			if (c == 0)
			{
				return 0;
			}
		}
		if(j == iterations){
			return -1;
		}
	}

	if (z_gcd(c, n) != 1) return 0;
	return 1;
}

/*
   The "n-1 primality test from Crandall and Pomerance, 
	"Primes: a computational perspective 2nd Ed.", pp 178.
*/

int z_isprime_nm1(unsigned long const n, unsigned long const iterations)
{
	if (n < 214) return z_isprime(n);
	
	int i, j, k, pass, test, exp;
	unsigned long cuberoot, n1, f, factor, prod, temp, b, c1, c2, cofactor;
	factor_t factors;
	pre_inv_t inv;
	
	/* 2 is prime, all other even numbers are not */
	
	if (n % 2 == 0)
	{
		if (n == 2) return 1;
		else return 0;
	}
	
	n1 = n - 1;
	f = 1;
	factors.num = 0;
	cuberoot = z_intcuberoot(n);
	
	cofactor = z_factor_partial(&factors, n1, cuberoot, 1);

	inv = z_precompute_inverse(n);

	pass = 0;
	
	for (j = 2; (j < iterations) && (pass == 0); j++)
	{
		pass = 1;

		if (!z_ispseudoprime_fermat(n, j)) return 0;
			
		for (i = 0; (i < factors.num) && (pass == 1); i++)
	   {
		   unsigned long exp = n1/factors.p[i];
		
			ulong g = z_gcd(z_powmod2(j, exp, n) - 1, n);
			if ((g > 1) && (g < n)) return 0;	
			if (g == n) pass = 0;
		}
	}

	if (j == iterations) return -1;

	ulong F = (n-1)/cofactor;
	
	if (F >= z_intsqrt(n)) return 1;
	
	c2 = (n - 1)/(F*F);
	c1 = cofactor - c2*F;
	if (2*c2 <= c1) return 1;
	if (z_issquare(c1*c1 - 4*c2)) return 0;
	else return 1;

}

/*
    Computes the pair (V_m,V_m+1) using a Lucas chain with the rules:
    V_2j=V_j*V_j-2 mod n
    V_(2j+1)=V_j*V_(j+1)-a mod n
    
*/

pair_t z_lchain_mod_precomp(ulong m, ulong a, ulong n, pre_inv2_t ninv)
{
	pair_t current, old;
	int length, i;
	ulong power, xy, xx, yy;
	
	old.x = 2;
	old.y = a;
	
	length = FLINT_BIT_COUNT(m);
	power = (1L<<(length-1));
	for (i = 0; i < length; i++)
	{
		xy = z_mulmod2_precomp(old.x, old.y, n, ninv);
		xy = z_mod2_precomp(xy - a + n, n, ninv);
		if (((m & power) > 0))
		{
			yy = z_mulmod2_precomp(old.y, old.y, n, ninv);
			yy = z_mod2_precomp(yy - 2L + n, n, ninv);
			current.x = xy;
			current.y = yy;
		} else 
		{
			xx = z_mulmod2_precomp(old.x, old.x, n, ninv);
			xx = z_mod2_precomp(xx -2L + n, n, ninv);
			current.x = xx;
			current.y = xy;
		}
		power = (power>>1);
		old = current;
	}

	return current;
}

/*
    Test to see if _n_ is a Lucas Pseudoprime. Returns 0 if 
    _n_ is composite or fails gcd(n, 2*a*b(a*a - 4*b)) = 1, -1 if delta
    is square and 1 if it is a Lucas Pseudoprime w.r.t x^2 - ax + b  
	 Assumes n has been checked for primality using trial factoring up to 256.
	 The absolute values of a and b should be < 128.
*/

int z_ispseudoprime_lucas_ab(ulong n, int a, int b)
{
	int b2, delta;
	ulong m, A, left, right, binv;
	pre_inv2_t inv;
	pair_t V;
	
	delta = a*a - 4*b;
	if (z_issquare(delta)) return -1;
	
	if (z_gcd(n, (long) (2*a*b*delta)) != 1) return 0;
	
	m = (n - z_jacobi(delta, n))/2;
	
	if (b < 0) A = a*a*z_invert(n + b, n) - 2;
	else A = a*a*z_invert(b, n) - 2;
	
	inv = z_precompute_inverse(n);
	if ((long) A < 0L) A = A + n;
	
	A = z_mod2_precomp(A, n, inv);
	
	V = z_lchain_mod_precomp(m, A, n, inv);
	
	left = z_mulmod2_precomp(A, V.x, n, inv);
	right = z_mulmod2_precomp(2, V.y, n, inv);
	
	if (left == right) return 1;
	else return 0;
}

/*
   A second version of the Lucas pseudoprime test. This uses the algorithm described in 
   http://www.jstor.org/sici?sici=0025-5718(198010)35%3A152%3C1391%3ALP%3E2.0.CO%3B2-N
   where it picks D, P, Q using algorithm A.
	Assumes n has been checked for primality using trial factoring up to 256.
*/

int z_ispseudoprime_lucas(ulong n)
{
	int i, D, Q;
	long A;
	ulong left, right;
	pair_t V;
	D = 0;
	Q = 0;
	D = 0;
	
	if (n % 2 == 0)
	{
		if (n == 2) return 1;
		else return 0;
	}
	
	for (i = 0; i < 100; i++)
	{
		D = 5 + 2*i;
		if (z_gcd(D, n) != 1) return 0;
		if (i % 2 == 1) D = -D;
		if (z_jacobi(D, n) == -1) break;
	}

	if (i == 100)
	{
		if (z_issquare(n)) return -1;
		else return 1;
	}
   
	Q = (1 - D)/4;
	A = z_invert(Q + n, n) - 2;
	
	pre_inv_t inv = z_precompute_inverse(n);
	V = z_lchain_mod_precomp(n+1, A, n, inv);
	
	left = z_mulmod2_precomp(A, V.x, n, inv);
	right = z_mulmod2_precomp(2, V.y, n, inv);
	
	if (left == right) return 1;
	else return 0;
}

pair_t z_fchain_mod_precomp(ulong m, ulong n, pre_inv2_t ninv)
{
	pair_t current, old;
	int length;
	ulong power, xy, xx, yy;
	
	old.x = 2;
	old.y = n - 3;
	
	length = FLINT_BIT_COUNT(m);
	power = (1L<<(length-1));

	for ( ; length > 0; length--)
	{
		xy = z_mulmod2_precomp(old.x, old.y, n, ninv);
		
		xy = xy + 3L;
		if (xy >= n) xy = xy - n;
		
		if (m & power)
		{
			current.y = z_mulmod2_precomp(old.y, old.y, n, ninv) - 2L;
			if ((long) current.y < 0L) current.y += n;
			current.x = xy;
		} else 
		{
			current.x = z_mulmod2_precomp(old.x, old.x, n, ninv);
			if (current.x >= 2L) current.x -= 2;
			else current.x = current.x - 2L + n;
			current.y = xy;
		}
		power = (power>>1);
		old = current;
	}

	return current;
}

/*
	Checks to see if n is a fibonacci pseudoprime. Assumes that 
	gcd(10, n) == 1
*/

int z_ispseudoprime_fibonacci_precomp(unsigned long n, pre_inv2_t inv)
{
	unsigned long m, left, right;
	pair_t V;
	m = (n - z_jacobi(5L, n))/2;
	V = z_fchain_mod_precomp(m, n, inv);
		
	return (z_mulmod2_precomp(n - 3L, V.x, n, inv) == z_mulmod2_precomp(2L, V.y, n, inv));
}

#define INV10 ((double) 1 / (double) 10)

/* 
   There are no known exceptions to the conjecture that this is a primality test
	if n > 13.
*/

int z_isprobab_prime_BPSW(unsigned long n)
{
	int nmod10;
	pre_inv2_t inv = z_precompute_inverse(n);
	nmod10 = z_mod2_precomp(n, 10, INV10);

	if (nmod10 == 3 || nmod10 == 7)
	{
		if (z_ispseudoprime_fermat(n, 2) == 0) return 0;
		return z_ispseudoprime_fibonacci_precomp(n, inv);
		
	} else
	{
		unsigned long d;
		d = n - 1;
		while ((d & 1) == 0) d >>= 1;

		if (SPRP_64(2L, d, n, inv) == 0) return 0;
		if (z_ispseudoprime_lucas(n) == 0) return 0;
		return 1;
	}
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

long z_xgcd(long* a, long* b, long x, long y)
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
   for (unsigned long i = 1; (i < SQFREE_TF_PRIMES_LIMIT) && (z_primes[i]*z_primes[i] <= n); i++)
   {
      quot = n/z_primes[i];
      rem = n - quot*z_primes[i];
      if (rem == 0) 
      { 
         if ((quot % z_primes[i]) == 0) return 0;
         else n = quot;
      }
   }
   return 1;
}

/*
   Tests if n is squarefree or not
*/

int z_issquarefree(unsigned long n, int proved)
{
   if (n < SQFREE_TF_CUTOFF) return z_issquarefree_trial(n);
   else 
   {
		if (!z_issquarefree_trial(n)) return 0;
			
		factor_t factors;
		
		z_factor(&factors, n, proved);
		for (ulong i = 0; i < factors.num; i++)
		{
			if ((factors.exp[i] & 1) == 0) return 0;
		}

		return 1;
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
   
   if (p == 2)
   {
      count_trail_zeros(exp, *n);
      if (exp)
      {
         *n = ((*n)>>exp);
      }
      return exp;      
   }

	quot = z_div2_precomp(*n, p, pinv);
   rem = (*n) - quot*p;
   while (rem == 0) 
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
      }
      return exp;      
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

void insert_factorpower(factor_t * factors, unsigned long p, unsigned long e)
{
   int i = 0;
   
   for (i = 0; i < factors->num; i++)
   {
      if (factors->p[i] == p)
      {
         factors->exp[i]+= e;
         break;
      }
   }
   if (i == factors->num)
   {
      factors->p[i] = p;
      factors->exp[i] = e;
      factors->num++;
   }
}

/*
   Finds all the factors of n by trial factoring up to some limit
   Returns the cofactor after removing these factors
*/

unsigned long _z_factor_trial(factor_t * factors, unsigned long n, unsigned long cutoff)
{
    int num_factors = 0;
    int exp;

    for (unsigned long i = 0; (i < cutoff) && (z_primes[i]*z_primes[i] <= n); i++)
    {
        exp = z_remove(&n, z_primes[i]);
        if (exp)
        {
/*            printf("found factor by trial division, %ld^%ld\n",z_primes[i],exp);*/
            factors->p[num_factors] = z_primes[i];
            factors->exp[num_factors] = exp;
            num_factors++;
        }

    }
       
    factors->num = num_factors;

    return n;
}

unsigned long z_factor_trial(factor_t * factors, unsigned long n) {
    return _z_factor_trial(factors, n, TF_CUTOFF);
}


#define ETF_CUTOFF 1000000
#define ETF_NUM_PRIMES 78330
#define ETF_SIEVE_SIZE 499500

unsigned long z_extended_primes[ETF_NUM_PRIMES];

/*
    Computes primes up to 1 million.  Primes up to 1000 are stored in _primes_, and
    primes between 1000 and 1 million are stored in extended_primes.
*/
static inline void z_initialize_extended_primes()
{
    static int already_computed = 0;
    if (!already_computed) { 
        already_computed = 1;
        z_compute_extended_primes();
/*        printf("have sieved\n");*/
    }
}

/*
    DO NOT call this, except for benchmarking purposes.  This should only ever get
    called once during a running process, and z_initialize_extended_primes does
    exactly that.
*/
void z_compute_extended_primes()
{
    unsigned int sieve[ETF_SIEVE_SIZE];
    unsigned int p,q,oldq = 0,n_found=0;
    unsigned int i;

    memset(sieve,1,ETF_SIEVE_SIZE);
    for (i = 1; z_primes[i]*z_primes[i] < 1000; i++) /* skip 2 */
    {
        p = z_primes[i];
        q = 1000-1000%p;
        q = q+(1+q%2)*p;
        q = q/2 - 500;   /* index of first multiple of p left in the sieve: first odd multiple after 1000 */
        while (q < ETF_SIEVE_SIZE)
        {
            sieve[q] = 0;
            q+=p;
        }
    }

    for (; i < TF_CUTOFF; i++)
    {
        p = z_primes[i];
        q = (p*p)/2 - 500;   /* index of first multiple of p left in the sieve: p^2 */
        for (unsigned int j = oldq; j < q; j++) 
        {
            if (sieve[j]) 
            {
                z_extended_primes[n_found] = 2*(j+500)+1;
                n_found++;
            }
        }
        oldq = q;
        while (q < ETF_SIEVE_SIZE)
        {
            sieve[q] = 0;
            q+=p;
        }
    }
    for (unsigned int j = oldq; j < ETF_NUM_PRIMES; j++) 
    {
        if (sieve[j]) 
        {
            z_extended_primes[n_found] = 2*(j+500)+1;
            n_found++;
        }
    }
}


/*
    Perform trial division on primes 1000 < p < 10^6.  Called in the middle of z_factor, so
    we need to make sure that we don't duplicate primes in the factor list.  This is really
    slow, and we typically expect other factor methods to work -- so dump out the first
    factor we find.
*/

unsigned long z_factor_trial_extended3(factor_t *factors, unsigned long n)
{
    int exp;
    unsigned long i = 0;
    unsigned long p = z_extended_primes[i];
    
    for (; (i < TF_CUTOFF) && (p*p*p <= 8*n); p = z_extended_primes[++i])
    {
        exp = z_remove(&n, z_primes[i]);
        if (exp)
        {
            insert_factorpower(factors,p,exp);
        }
        p = z_extended_primes[i];
    }

    return n;
}



unsigned long z_factor_trial_extended(unsigned long n)
{
    z_initialize_extended_primes();

    for (unsigned long i = 0; (i < ETF_NUM_PRIMES) && (z_extended_primes[i]*z_extended_primes[i] <= n); i++)
    {
        if (n%z_extended_primes[i] == 0UL)
            return z_extended_primes[i];
    }

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
   
   for (unsigned long i = 0; (i < TF_CUTOFF) && ((z_primes[i]-1)*(z_primes[i]-1) <= n); i++)
   {
      exp = z_remove(&n, z_primes[i]);
      if (exp)
      {
         factors->p[num_factors] = z_primes[i];
         factors->exp[num_factors] = exp;
         num_factors++;
		   *prod *= z_pow(z_primes[i], exp);
		 if (*prod > limit) break;
      }  
   }
       
   factors->num = num_factors;
   return n;
}

/*
   Factors _n_ until the product of the factor found is > _limit_. It puts the factors
   in _factors_ and returns the cofactor.
	If proved is 0 (false) the factors are not proved prime, otherwise the result is
	proved.
*/

unsigned long z_factor_partial(factor_t * factors, unsigned long n, unsigned long limit, int proved)
{
	unsigned long prod, cofactor;
	unsigned long factor_arr[TF_FACTORS_IN_LIMB];
	unsigned long exp_arr[TF_FACTORS_IN_LIMB];
	unsigned long cutoff = z_primes[TF_CUTOFF-1]*z_primes[TF_CUTOFF-1];
	unsigned long factors_left = 1;
	unsigned long factor;
   unsigned long exp;

   cofactor = z_factor_partial_trial(factors, &prod, n, limit);
   if (prod != n && prod <= limit)
	{
	   factor = factor_arr[0] = cofactor;
      exp_arr[0] = 1;
      		
		while (factors_left > 0 && prod <= limit)
		{
         factor = factor_arr[factors_left-1];

         if (factor >= cutoff)
         {
            if (cofactor = z_factor_235power(factor, &exp) )
            {
               exp_arr[factors_left-1] *= exp;
               factor_arr[factors_left-1] = factor = cofactor;
            } 
				if ((factor >= cutoff) && !(proved && z_isprime(factor)) && !(!proved && z_isprobab_prime(factor)))
            {
               if (
                  (
                    (factor < MAX_HOLF) && 
                    (factor > MIN_HOLF) &&
                    (cofactor = z_factor_HOLF(factor,HOLF_ITERS))
                  ) ||
                  (
                    (factor > MAX_HOLF) &&
                    (cofactor = z_factor_HOLF(factor,100))
                  ) ||
                  ( cofactor = z_factor_SQUFOF(factor) ) ||
                  ( cofactor = z_factor_trial_extended(factor) )
               )
               {
                  exp_arr[factors_left] = exp_arr[factors_left-1];
                  factor_arr[factors_left] = cofactor;
                  factor_arr[factors_left-1] /= cofactor;
                  factors_left++;
               } else
               {
                  printf("Error : failed to factor %ld\n", n);
                  abort();
               }
            } else
            {
               insert_factorpower(factors, factor, exp_arr[factors_left-1]);
               prod *= z_pow(factor, exp_arr[factors_left-1]);
               factors_left--;
            }
			} else
         {
            insert_factorpower(factors, factor, exp_arr[factors_left-1]);
            prod *= z_pow(factor, exp_arr[factors_left-1]);
            factors_left--;
			}
		}
		return n/prod;
	} 
	return cofactor;
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

unsigned long _z_ll_factor_SQUFOF(ulong n_hi, ulong n_lo)
{
   mp_limb_t n[2];
	mp_limb_t sqrt[2];
	mp_limb_t rem[2];
	
	n[0] = n_lo;
	n[1] = n_hi;
	
	int num = mpn_sqrtrem(sqrt, rem, n, 2);
	
	unsigned long sqroot = sqrt[0];
   unsigned long p = sqroot;
   unsigned long q = rem[0];
   
   if ((q == 0) || (num == 0))
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

	umul_ppmm(rem[1], rem[0], p, p);
   sub_ddmmss(sqrt[1], sqrt[0], n[1], n[0], rem[1], rem[0]);
	if (sqrt[1])
	{
		int norm;
	   count_lead_zeros(norm, sqrt[1]);
	   udiv_qrnnd(q, rem[0], sqrt[1]<<norm + r_shift(sqrt[0], FLINT_BITS-norm), sqrt[0]<<norm, qlast); 
      q >>= norm;
	} else
	{
		q = sqrt[0]/qlast;
	}

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
    unsigned long s1, s2, i;
   
    for (i = 1; (i < NUMBER_OF_PRIMES) && !factor; i++)
    {
        multiplier = z_primes[i];
        s1 = FLINT_BIT_COUNT(multiplier);
        count_lead_zeros(s2, n);

        if (s1 > s2) // kn is possibly more than one limb 
        {
            mp_limb_t multn[2];
            umul_ppmm(multn[1], multn[0], multiplier, n);
            if (multn[1] == 0)
            {
                kn = multiplier*n;
                factor = _z_factor_SQUFOF(kn);
            } else
            {
                factor = _z_ll_factor_SQUFOF(multn[1], multn[0]);
            }
        } else
        {
            kn = multiplier*n;
            factor = _z_factor_SQUFOF(kn);
        }

        if (factor) 
        {
            quot = factor/multiplier;
            rem = factor - quot*multiplier;
            if (!rem) factor = quot;
            if ((factor == 1) || (factor == n)) factor = 0;
        }
    }

    if (i == NUMBER_OF_PRIMES) return 0;

    return factor; 
}

/*
   This is a reasonable implementation of Bill Hart's "one line" factor algorithm. (HOLF)
   WARNING: this is more than one line.
*/

unsigned long z_factor_HOLF(unsigned long n,unsigned long iters)
{
    unsigned long orig_n=n, in, square, sqrti, mod, factor, factoring = iters, iin;
    n*=HOLF_MULTIPLIER;

    iin = 0;
    in = n;
    while (factoring && (iin < in))
    {
        sqrti = z_intsqrt(in);
        sqrti++;
        square = sqrti*sqrti;
        mod = square-in;
        if (z_issquare2(mod,&factor)) 
        {
            sqrti -= factor;
            factor = z_gcd(orig_n,sqrti);  
            if (factor != 1) 
            { 
                return factor;
            }
        }     
        factoring--;    
        iin = in;
        in += n;
    }

    return 0;
}


unsigned long z_factor_tinyQS(unsigned long n) {
    F_mpz_factor_t factors;
    mpz_t N;
    unsigned long factor;
    mpz_set_ui(N,n);
    factors.fact = malloc(64*sizeof(mpz_t *));
    factor = (unsigned long) F_mpz_factor_tinyQS_silent(factors,N);
    if (factor > 1UL)
        factor = mpz_get_ui(factors.fact[0]);

    for(int i=0;i<factors.num;i++)
        mpz_clear(factors.fact[i]);
    free(factors.fact);

    return factor;
}



/*
   Find the factors of n.
   If proved is 0 the factors are not proved prime, otherwise the
	result is proved.
*/

#define MAX_SQUFOF 0xFFFFFFFFFFFFF

void z_factor(factor_t * factors, unsigned long n, int proved)
{
    unsigned long cofactor;
    unsigned long factor_arr[TF_FACTORS_IN_LIMB];
    unsigned long exp_arr[TF_FACTORS_IN_LIMB];
    unsigned long cutoff;
    unsigned long factors_left = 1, k;
    unsigned long factor;
    unsigned long exp;

    if  (n < 0x7FFFFFFFUL)
        if (n < 0xFFFFFUL) cutoff = TF_CUTOFF;     /* can't beat trial division under 20 bits */
        else cutoff = 31;                          /* 7-bit primes for 31-bit n */
    else 
        if (n < 0x3FFFFFFFFUL) cutoff = 97;        /* 9-bit primes for 35-bit n */
        else cutoff = TF_CUTOFF;                   /* 10-bit primes for the rest*/

    cofactor = _z_factor_trial(factors, n, cutoff);
    if (cofactor == 1UL) return;
    factor = factor_arr[0] = cofactor;
    exp_arr[0] = 1;
    cutoff = z_primes[cutoff-1]*z_primes[cutoff-1];

    while (factors_left > 0)
    {
        factor = factor_arr[factors_left-1];

        if (factor >= cutoff)
		  {
			  if( cofactor = z_factor_235power(factor, &exp) ) 
           {
              exp_arr[factors_left-1] *= exp;
              factor_arr[factors_left-1] = factor = cofactor;
           }
           if (
              (factor >= cutoff) &&
              !(proved && z_isprime(factor)) && 
              !(!proved && z_isprobab_prime(factor))
           )
			  {
				  if (
                 (
                    (factor < MAX_HOLF) && 
                    (factor > MIN_HOLF) &&
                    (cofactor = z_factor_HOLF(factor,HOLF_ITERS))
                 ) ||
                 (cofactor = z_factor_SQUFOF(factor)) ||
                 ( cofactor = z_factor_trial_extended(factor) ) ||
                 ( cofactor = z_factor_tinyQS(factor) ) 
              )
				  {
					   exp_arr[factors_left] = exp_arr[factors_left-1];
                  factor_arr[factors_left] = cofactor;
                  factor_arr[factors_left-1] /= cofactor;
                  factors_left++;
				  } else
				  {
                 printf("Error : failed to factor %ld\n", n);
                 abort();
				  }
			  } else
			  {
				  insert_factorpower(factors, factor, exp_arr[factors_left-1]);
              factors_left--;
			  }
		  } else
		  {
			  insert_factorpower(factors, factor, exp_arr[factors_left-1]);
           factors_left--;
		  }
    }
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
   
   z_factor(&factors, (p - 1), 1);
   
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
   
   z_factor(&factors, (p - 1), 1);
   
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

/*unsigned long z_intsqrt(unsigned long r)
{
#if FLINT_BITS == 32
    float x, z;
	union {
	  float f;
      unsigned long l;
	} temp;

	temp.f = (float) r;
	temp.l = (0xbe6ec85eUL - temp.l)>>1; // estimate of 1/sqrt(y) 
	x =  temp.f;
	z =  (float) r*0.5;                        
    x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
    unsigned long res =  is + ((is+1)*(is+1) <= r);
    return res - (res*res > r);
#else
    unsigned long is;
	
	double x, z;
	union {
	  double f;
      unsigned long l;
	} temp;

	temp.f = (double) r;
	temp.l = (0xbfcdd90a00000000UL - temp.l)>>1; // estimate of 1/sqrt(y) 
	x =  temp.f;
	z =  (double) r*0.5;                        
    x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
	x = (1.5*x) - (x*x)*(x*z);
    is = (unsigned long) (x*(double) r);
    unsigned long res =  is + ((is+1)*(is+1) <= r);
    return res - (res*res > r);
#endif
}*/

/*
   Returns the integer part of the square root of r
   i.e. floor(sqrt(r))
*/
unsigned long z_intsqrt(unsigned long r)
{
   unsigned long is;
#if FLINT_BITS == 64
   if (r >= 70368744177664UL)
   {
	  is = (unsigned long) floor(sqrt((double) r));
   } else
   {
#endif
	   is = (unsigned long) floorf(sqrtf((float) r));
#if FLINT_BITS == 64
   }
#endif
   unsigned long res =  is + ((is+1)*(is+1) <= r);
   return res - (res*res > r);
}

unsigned long z_intcuberoot(unsigned long n)
{
    unsigned long newg, g;
    unsigned long oldg = n;

    if (n < 8) // avoid infinite loops
    {
       if (n == 0) return 0;
       return 1;
    }

	g = (unsigned long) ceil(pow(n, 0.333333333));
	
    do 
	{
		newg = 2L*g + n/(g*g);
		newg = newg/3L;
	
        if (g == newg || newg == oldg) break;
        
        oldg = g;
		g = newg;
	} while (1);

    if (newg*newg*newg > n) return newg - 1;
    else return newg;
}

unsigned long z_intfifthroot(unsigned long n)
{
	unsigned long newg, g;
    unsigned long oldg = n;
    
    if (n < 243) // special cases can lead to infinite loops
    {
       if (n == 0) return 0;
       if (n < 32) return 1;
       return 2;
    }

	g = (unsigned long)ceil(pow(n, 0.2)); 
	
	do 
	{
        newg = g*g;
		newg = 4L*g + n/(newg*newg);
		newg = newg/5L;
	
        if (newg == g || newg == oldg) break;
		oldg = g;
		g = newg;
	} while (1);

	if (newg*newg*newg*newg*newg > n) return newg - 1;
    else return newg;
}

/*
    Factors perfect squares, cubes, and fifth-powers.

    We do trial division up to 2^10 > (2^64)^(1/7).  Hence, if an unsigned word is a perfect
    power, it must be a square, cube, or fifth power.  If a number if a number is *not* a power
    modulo some n, it cannot be a power (of course, the converse does not hold).  We have found
    four moduli such that this test rejects 98.7% of non-powers.

    By keeping track of which possible powers (of 2,3,5) a power may be, we can now make a very 
    small number of tests to determine if the number is indeed a perfect power.

    We have not optimized the order the tests come in.  That will require some thought (or, brute
    force).
*/
unsigned long z_factor_235power(unsigned long n, unsigned long *exp)
{
    static char mod63[63] = {7,7,4,0,5,4,0,5,6,5,4,4,0,4,4,0,5,4,5,4,4,0,5,4,0,5,4,6,7,4,0,4,4,0,4,6,7,5,4,0,4,4,0,5,4,4,5,4,0,5,4,0,4,4,4,6,4,0,5,4,0,4,6};
    static char mod61[61] = {7,7,0,3,1,1,0,0,2,3,0,6,1,5,5,1,1,0,0,1,3,4,1,2,2,1,0,3,2,4,0,0,4,2,3,0,1,2,2,1,4,3,1,0,0,1,1,5,5,1,6,0,3,2,0,0,1,1,3,0,7};
    static char mod44[44] = {7,7,0,2,3,3,0,2,2,3,0,6,7,2,0,2,3,2,0,2,3,6,0,6,2,3,0,2,2,2,0,2,6,7,0,2,3,3,0,2,2,2,0,6};
    static char mod31[31] = {7,7,3,0,3,5,4,1,3,1,1,0,0,0,1,2,3,0,1,1,1,0,0,2,0,5,4,2,1,2,6};
    char t;
    
    t = mod31[n%31];
    if (!t) return 0;
    t&= mod44[n%44];
    if (!t) return 0;
    t&= mod61[n%61];
    if (!t) return 0;
    t&= mod63[n%63];
    if (t&1) {
        unsigned long y = (unsigned long) sqrt((double)n);
        if (n == y*y)
        {
            *exp = 2;
            return y;
        }
    }
    if (t&2) {
        unsigned long y = (unsigned long) round(pow((double)n,1/3.));
        if (n == y*y*y)
        {
            *exp = 3;
            return y;
        }
    }
    if (t&4) {
        unsigned long y = (unsigned long) round(pow((double)n,1/5.));
        if (n == y*y*y*y*y)
        {
            *exp = 5;
            return y;
        }
    }
    return 0;
}
