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
/******************************************************************************

 long_extras.h
 Header file for long_extras.c.

 (C) 2006 William Hart
 
 Some of the macros in this file were borrowed from GMP, (C) Free Software Foundation

******************************************************************************/

#ifndef LONGEXTRAS_H
#define LONGEXTRAS_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <math.h>

#include "longlong_wrapper.h"
#include "longlong.h"

//======================================================================================
//
//   The code in this section is borrowed from the GMP library v 4.2.1
//   See gmp-impl.h, (C) Free Software Foundation
//
//======================================================================================

#define invert_limb(invxl, xl)                  \
  do {                                          \
    mp_limb_t dummy;                            \
    udiv_qrnnd (invxl, dummy, ~(xl), ~(0L), xl);  \
  } while (0)
  
#define LIMB_HIGHBIT_TO_MASK(n)                                 \
  (((mp_limb_signed_t) -1 >> 1) < 0                             \
   ? (mp_limb_signed_t) (n) >> (FLINT_BITS - 1)              \
   : (n) & (1L<<(FLINT_BITS-1)) ? (~ (mp_limb_t) 0L) : (0L))

#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)				\
  do {									\
    mp_limb_t _n2, _n10, _nmask, _nadj, _q1;				\
    mp_limb_t _xh, _xl;							\
    _n2 = (nh);								\
    _n10 = (nl);							\
    _nmask = LIMB_HIGHBIT_TO_MASK (_n10);				\
    _nadj = _n10 + (_nmask & (d));					\
    umul_ppmm (_xh, _xl, di, _n2 - _nmask);				\
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);			\
    _q1 = ~_xh;								\
    umul_ppmm (_xh, _xl, _q1, d);					\
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);				\
    _xh -= (d);					/* xh = 0 or -1 */	\
    (r) = _xl + ((d) & _xh);						\
    (q) = _xh - _q1;							\
  } while (0)

//=====================================================================================

typedef struct factor_s
{
   int num;
   unsigned long p[15];
   unsigned long exp[15];
} factor_t;

#define pre_inv_t double
#define pre_inv2_t double

unsigned long z_randint(unsigned long limit);

double z_precompute_inverse(unsigned long n);

static inline
unsigned long z_addmod(unsigned long a, unsigned long b, unsigned long p)
{
   unsigned long neg1 = p - a;
   if (neg1 > b)
      return a + b;
   else 
      return b - neg1;
}

static inline
unsigned long z_submod(unsigned long a, unsigned long b, unsigned long p)
{
   if (a < b)
      return p + a - b;
   else
      return a - b;
}

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

int z_isprime_precomp(unsigned long n, double ninv);

unsigned long z_nextprime(unsigned long n);

unsigned long z_pow(unsigned long a, unsigned long exp);

unsigned long z_powmod(unsigned long a, long exp, unsigned long n);
                                                    
unsigned long z_sqrtmod(unsigned long a, unsigned long p); 

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

unsigned long z_CRT(unsigned long x1, unsigned long n1, 
                        unsigned long x2, unsigned long n2);
                       
int z_issquarefree(unsigned long n);

int z_remove_precomp(unsigned long * n, unsigned long p, double pinv);

int z_remove(unsigned long * n, unsigned long p);

unsigned long z_factor_trial(factor_t * factors, unsigned long n);

unsigned long z_factor_SQUFOF(unsigned long n);

int z_factor(factor_t * factors, unsigned long n);

unsigned long z_primitive_root(unsigned long p);
unsigned long z_primitive_root_precomp(unsigned long p, double p_inv);

#ifdef __cplusplus
 }
#endif
 
#endif






