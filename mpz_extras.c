/*============================================================================

    Copyright (C) 2007, William Hart

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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "mpz_extras.h"
#include "flint.h"
#include "mpn_extras.h"
#include "F_mpn_mul-tuning.h"
#include "memory-manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"

#define DEBUG2 1
#define DEBUG 0

/*
   Memory manager to allocate a single mpz_t. It returns a pointer to the mpz_t. 
   mpz_t's should be released in the order they were allocated.
*/

#define RESALLOC 100 //allocate this many mpz_t's at once to save on overheads

mpz_t** reservoir; // Array of pointers to mpz_t's in the reservoir
unsigned long rescount=0; //Next available mpz_t in reservoir
unsigned long currentalloc=0; //total number of mpz_t's in reservoir

mpz_t* F_mpz_alloc(void)
{
   static int initialised = 0;
   static mpz_t** tempres;
   mpz_t* alloc_d;
   
   //allocate another block of mpz_t's if none are currently allocated, or the reservoir is depleted
   if (rescount==currentalloc) // need more limb_memp_t's
   {
      if (!initialised) 
      {
         reservoir = (mpz_t**)malloc(RESALLOC*sizeof(mpz_t*)); //allocate space for the array of pointers
         reservoir[0] = (mpz_t*)malloc(RESALLOC*sizeof(mpz_t)); //allocate space for the mpz_t's
         for (unsigned long i=0; i<RESALLOC-1; i++)
         {
             reservoir[i+1]=reservoir[i]+1; //initialise the array
             mpz_init(*reservoir[i]); //initialise the mpz_t's
         }
         mpz_init(*reservoir[RESALLOC-1]);
         rescount=0;
         initialised = 1;
         currentalloc = RESALLOC;
      } else
      {
         //copy old reservoir into larger one
         tempres = reservoir;
         reservoir = (mpz_t**)malloc((currentalloc+RESALLOC)*sizeof(mpz_t*));
         reservoir[currentalloc] = (mpz_t*)malloc(RESALLOC*sizeof(mpz_t));  
         memcpy(reservoir,tempres,currentalloc*sizeof(mpz_t*)); 
         for (unsigned long i=currentalloc; i<RESALLOC+currentalloc-1; i++)
         {
             reservoir[i+1]=reservoir[i]+1; //initialise the array
             mpz_init(*reservoir[i]); //initialise the mpz_t's
         }
         mpz_init(*reservoir[currentalloc+RESALLOC-1]);
         
         currentalloc+=RESALLOC;
         //free old reservoir
         free(tempres);  
      }       
   }
   
   alloc_d = reservoir[rescount];
   rescount++;
   return alloc_d;
}

/*
    Release a mpz_t back into the reservoir
*/

void F_mpz_release(void)
{   
    rescount--;
}

/* 
    sets res to a*b modulo p
    assumes res is not p
    Does not assume a and b are reduced mod p
*/

void F_mpz_mulmod(mpz_t res, mpz_t a, mpz_t b, mpz_t p)
{
     mpz_t* temp = F_mpz_alloc();
     
     mpz_fdiv_r(*temp,a,p);
     mpz_fdiv_r(res,b,p);
     mpz_mul(res,*temp,res);
     mpz_fdiv_r(res,res,p);
     
     F_mpz_release();
     
     return;
}

/* 
     Sets res to a*b modulo p
     Does not assume a and b are reduced mod p
*/

unsigned long F_mpz_mulmod_ui(mpz_t res, mpz_t a, mpz_t b, unsigned long p)
{
   unsigned long p1, p2, al, bl;
   
   al = mpz_fdiv_r_ui(res, a, p);
   bl = mpz_fdiv_r_ui(res, b, p);
   
   umul_ppmm(p2, p1, al, bl);
   
   unsigned long norm, q, r;
   
   if (p2 >= p) p2 %= p;
   
#if UDIV_NEEDS_NORMALIZATION
   count_lead_zeros(norm, p);
   udiv_qrnnd(q, r, (p2<<norm) + (p1>>(FLINT_BITS-norm)), p1<<norm, p<<norm);
#else
   udiv_qrnnd(q, r, p2, p1, p);
#endif
   
   mpz_set_ui(res, r);
   
   return r;
}

/* 
     Sets res to the square root of a modulo p for a prime p
     Returns 0 if a is not a square modulo p
*/

int F_mpz_sqrtmod(mpz_t res, mpz_t a, mpz_t p) 
{
     int r,k,m;
     mpz_t* p1 = F_mpz_alloc();
     mpz_t* two = F_mpz_alloc();
     mpz_t* b = F_mpz_alloc();
     mpz_t* g = F_mpz_alloc();
     mpz_t* bpow = F_mpz_alloc();
     mpz_t* gpow = F_mpz_alloc();
     mpz_t* mk = F_mpz_alloc();
     
     if (mpz_kronecker(a,p)!=1) 
     {
         mpz_set_ui(res,0);
         return 0;   //return 0 if a is not a square mod p
     }
     
     if ((mpz_cmp_ui(a,0)==0)||(mpz_cmp_ui(a,1)==0)) 
     {
        mpz_set(res,a);
        return 1;
     }
     
     if ((mpz_tstbit(p,0)==1)||(mpz_tstbit(p,1)==1))
     {
        mpz_add_ui(*p1,p,1);
        mpz_fdiv_q_2exp(*p1,*p1,2);
        mpz_powm(res,a,*p1,p);
        return 1;
     }
     
     mpz_set_ui(*two,2);
     
     mpz_sub_ui(*p1,p,1);
     r = mpz_remove(*p1,*p1,*two);
     mpz_powm(*b,a,*p1,p);
     for (k=2; ;k++)
     {
         if (mpz_ui_kronecker(k,p) == -1) break;
     }
     mpz_set_ui(*mk,k);
     mpz_powm(*g,*mk,*p1,p);
     mpz_add_ui(*p1,*p1,1);
     mpz_divexact_ui(*p1,*p1,2);
     mpz_powm(res,a,*p1,p);
     if (!mpz_cmp_ui(*b,1)) return 1;
     
     while (mpz_cmp_ui(*b,1))
     {
           mpz_set(*bpow,*b);
           for (m=1; (m<=r-1) && (mpz_cmp_ui(*bpow,1));m++)
           {
               mpz_powm_ui(*bpow,*bpow,2,p);
           }
           mpz_set(*gpow,*g);
           for (int i = 1;i<= r-m-1;i++)
           {
               mpz_powm_ui(*gpow,*gpow,2,p);
           };
           F_mpz_mulmod(res,res,*gpow,p);
           mpz_powm_ui(*gpow,*gpow,2,p);
           F_mpz_mulmod(*b,*b,*gpow,p);
           mpz_set(*gpow,*g);
           r = m;
     }
          
     F_mpz_release();F_mpz_release();F_mpz_release();F_mpz_release();
     F_mpz_release();F_mpz_release();F_mpz_release();
     
     return 1;
}

/* 
     Computes the square root of a modulo p^k when given z, the square root mod p^(k-1)
*/

void __sqrtmodpow(mpz_t res, mpz_t z, mpz_t a, mpz_t pk, mpz_t tempsqpow, mpz_t inv)
{
     mpz_mul_ui(inv,z,2);
     mpz_invert(inv,inv,pk);
     mpz_set(tempsqpow,a);
     mpz_submul(tempsqpow,z,z);
     mpz_fdiv_r(tempsqpow,tempsqpow,pk);
     F_mpz_mulmod(tempsqpow,tempsqpow,inv,pk);
     mpz_add(tempsqpow,tempsqpow,z);
     mpz_set(res,tempsqpow);
     
     return;
} 

/* 
     Computes the square root of a modulo p^k when given z, the square root mod p^{k-1}
*/

void F_mpz_sqrtmodpklift(mpz_t res, mpz_t z, mpz_t a, mpz_t pk)
{
     mpz_t* tempsqpow = F_mpz_alloc();
     mpz_t* inv = F_mpz_alloc();
     
     __sqrtmodpow(res, z, a, pk, *tempsqpow, *inv);
     
     F_mpz_release();F_mpz_release();
}

/* 
     computes the square root of a modulo p^k when given the square root modulo p
*/

void F_mpz_sqrtmodptopk(mpz_t res, mpz_t sqrt, mpz_t a, mpz_t p, int k)
{
     mpz_t* tempsqpow = F_mpz_alloc();
     mpz_t* inv = F_mpz_alloc();
     mpz_t* pk = F_mpz_alloc();
     
     mpz_set(res,sqrt);
     mpz_set(*pk,p);
     for (int i = 2; i<=k; i++)
     {
            mpz_mul(*pk,*pk,p);
            __sqrtmodpow(res,res,a,*pk, *tempsqpow, *inv);
     }
     
     F_mpz_release();F_mpz_release();F_mpz_release();
}

/* 
     Computes the square root of a modulo p^k 
     Returns 0 if the square root of a does not exist mod p
*/

int F_mpz_sqrtmodpk(mpz_t res, mpz_t a, mpz_t p, int k)
{
     int exists;
     mpz_t* sqrtmodp = F_mpz_alloc();
     
     exists = F_mpz_sqrtmod(*sqrtmodp,a,p);
     F_mpz_sqrtmodptopk(res,*sqrtmodp,a,p,k);
     
     F_mpz_release();
     
     return exists;
}


/* 
     Find res mod n=n1*n2 such that res = x1 mod n1 and res = x2 mod n2
*/
 
void F_mpz_CRT(mpz_t res, mpz_t n, mpz_t x1, mpz_t x2, mpz_t n1, mpz_t n2)
{
     mpz_t* ntemp = F_mpz_alloc();
     mpz_t* restemp = F_mpz_alloc();
     mpz_t* chtemp = F_mpz_alloc();
     
     mpz_mul(*ntemp,n1,n2);
     mpz_invert(*restemp,n2,n1);
     F_mpz_mulmod(*restemp,res,n2,*ntemp);
     F_mpz_mulmod(*restemp,*restemp,x1,*ntemp);
     
     mpz_invert(*chtemp,n1,n2);
     F_mpz_mulmod(*chtemp,*chtemp,n1,*ntemp);
     F_mpz_mulmod(*chtemp,*chtemp,x2,*ntemp);
     
     mpz_add(res,*restemp,*chtemp);
     mpz_fdiv_r(res,res,*ntemp);
     mpz_set(n,*ntemp);
     
     F_mpz_release();F_mpz_release();F_mpz_release();
     
     return;
}

/*
    Compute the Montgomery reduced form of a mod m
    Returns n such that m < 2^n (n will be divisible by FLINT_BITS)
    Assumes a is already reduced mod m
*/

unsigned long F_mpz_mont_red(mpz_t res, mpz_t a, mpz_t m)
{
   unsigned long n = mpz_size(m)*FLINT_BITS;
   
   mpz_mul_2exp(res, a, n);
   mpz_mod(res, res, m);
   
   return n;
}

/* 
    Compute the Montgomery multiplication r = a*b mod m assuming a and b are in 
    Montgomery form with respect to 2^n where m < 2^n and R is -m mod 2^n
*/

void F_mpz_mont_mul(mpz_t res, mpz_t a, mpz_t b, mpz_t m, mpz_t R, unsigned long n)
{
   mpz_t x, s;
   mpz_init(x);
   mpz_init(s);
   
   mpz_mul(x, a, b);
   mpz_fdiv_r_2exp(s, x, n);
   mpz_mul(s, s, R);
   mpz_fdiv_r_2exp(s, s, n);
   mpz_mul(res, s, m);
   mpz_add(res, res, x);
   mpz_fdiv_q_2exp(res, res, n);
   
   if (mpz_cmp(res, m) >= 0) mpz_sub(res, res, m);
    
   mpz_clear(x);
   mpz_clear(s);
}

/* 
    Compute a^exp mod m using Montgomery reduction
    Requires that m is odd and positive and that exp is positive
*/

void F_mpz_expmod_mont(mpz_t res, mpz_t a, mpz_t exp, mpz_t m)
{
   unsigned long n;
   unsigned long bits = mpz_sizeinbase(exp, 2);
   mpz_t aRED;
   mpz_t powRED;
   mpz_t R;
   mpz_t temp;
   int flag = 0;
   
   mpz_init(aRED);
   mpz_init(powRED);
   mpz_init(R);
   mpz_init(temp);
   
   n = F_mpz_mont_red(aRED, a, m);
   
   mpz_set_ui(temp, 1);
   mpz_mul_2exp(temp, temp, n);
   mpz_invert(R, m, temp);
   mpz_sub(R, temp, R);
   if (mpz_cmp(R, temp) == 0) mpz_sub(R, R, temp);
   
   mpz_set(powRED, aRED);
#ifdef DEBUG
   gmp_printf("powRED = %Zd\n", powRED);
#endif
   
   for (unsigned long i = 0; i < bits - 1; i++)
   {
      if (mpz_tstbit(exp, i))
      {
         if (flag) F_mpz_mont_mul(res, res, powRED, m, R, n);
         else 
         {
            mpz_set(res, powRED);
            flag = 1;
         }
      }
      F_mpz_mont_mul(powRED, powRED, powRED, m, R, n);
#ifdef DEBUG
      gmp_printf("powRED = %Zd\n", powRED);
#endif
   }
   
   if (flag) F_mpz_mont_mul(res, res, powRED, m, R, n);
   else mpz_set(res, powRED);
   
   mpz_set_ui(temp, 1);
   F_mpz_mont_mul(res, res, temp, m, R, n);
   
   mpz_clear(temp);
   mpz_clear(R);
   mpz_clear(powRED);
   mpz_clear(aRED);
}

void F_mpz_divrem_BZ(mpz_t Q, mpz_t R, mpz_t A, mpz_t B)
{
   unsigned long n = mpz_size(B);
   unsigned long m = mpz_size(A) - n;
   
   if ((long) m < 0)
   {
      mpz_set_ui(Q, 0);
      mpz_set(R, A);
      return;
   }
   
   if (m < 64) 
   {
      mpz_fdiv_qr(Q, R, A, B);
      return;
   }
   
   unsigned long k = m/2;
   
   mpz_t * B0 = F_mpz_alloc();
   mpz_t * B1 = F_mpz_alloc();
   mpz_t * A0 = F_mpz_alloc();
   mpz_t * A1 = F_mpz_alloc();
   mpz_t * Q0 = F_mpz_alloc();
   mpz_t * Q1 = F_mpz_alloc();
   mpz_t * R0 = F_mpz_alloc();
   mpz_t * R1 = F_mpz_alloc();
   mpz_t * temp = F_mpz_alloc();
   mpz_t * temp2 = F_mpz_alloc();
   mpz_t * temp3 = F_mpz_alloc();
   
   
   mpz_fdiv_q_2exp(*B1, B, FLINT_BITS*k);
   mpz_fdiv_q_2exp(*A1, A, FLINT_BITS*2*k);
   
   F_mpz_divrem_BZ(*Q1, *R1, *A1, *B1);
   mpz_fdiv_r_2exp(*B0, B, FLINT_BITS*k);
   mpz_fdiv_r_2exp(*A0, A, FLINT_BITS*2*k);
   mpz_mul_2exp(*temp, *R1, FLINT_BITS*2*k);
   mpz_add(*temp, *temp, *A0);
   mpz_mul_2exp(*temp2, *Q1, FLINT_BITS*k);
   mpz_mul(*temp2, *temp2, *B0);
   mpz_sub(*temp, *temp, *temp2);
   mpz_mul_2exp(*temp2, B, FLINT_BITS*k);
   
   while (mpz_cmp_ui(*temp, 0) < 0) 
   {
      mpz_sub_ui(*Q1, *Q1, 1);
      mpz_add(*temp, *temp, *temp2);
   }
   mpz_fdiv_q_2exp(*temp2, *temp, FLINT_BITS*k); 
   F_mpz_divrem_BZ(*Q0, *R0, *temp2, *B1);
   
   mpz_fdiv_r_2exp(*temp2, *temp, FLINT_BITS*k);
   mpz_mul_2exp(R, *R0, FLINT_BITS*k);
   mpz_add(R, R, *temp2);
   mpz_submul(R, *Q0, *B0);
   while (mpz_cmp_ui(R, 0) < 0) 
   {
      mpz_sub_ui(*Q0, *Q0, 1);
      mpz_add(R, R, B);
   }
   mpz_mul_2exp(Q, *Q1, FLINT_BITS*k);
   mpz_add(Q, Q, *Q0);
   
   F_mpz_release(); F_mpz_release(); F_mpz_release(); F_mpz_release();
   F_mpz_release(); F_mpz_release(); F_mpz_release(); F_mpz_release();
   F_mpz_release(); F_mpz_release(); F_mpz_release();
}

void F_mpz_rem_BZ(mpz_t R, mpz_t A, mpz_t B)
{
   unsigned long n = mpz_size(B);
   unsigned long m = mpz_size(A) - n;
   
   if ((long) m < 0)
   {
      mpz_set(R, A);
      return;
   }
   
   if (m < 64) 
   {
      mpz_fdiv_r(R, A, B);
      return;
   }
   
   unsigned long k = m/2;
   
   mpz_t * B0 = F_mpz_alloc();
   mpz_t * B1 = F_mpz_alloc();
   mpz_t * A0 = F_mpz_alloc();
   mpz_t * A1 = F_mpz_alloc();
   mpz_t * Q0 = F_mpz_alloc();
   mpz_t * Q1 = F_mpz_alloc();
   mpz_t * R0 = F_mpz_alloc();
   mpz_t * R1 = F_mpz_alloc();
   mpz_t * temp = F_mpz_alloc();
   mpz_t * temp2 = F_mpz_alloc();
   mpz_t * temp3 = F_mpz_alloc();
   
   
   mpz_fdiv_q_2exp(*B1, B, FLINT_BITS*k);
   mpz_fdiv_q_2exp(*A1, A, FLINT_BITS*2*k);
   
   F_mpz_divrem_BZ(*Q1, *R1, *A1, *B1);
   mpz_fdiv_r_2exp(*B0, B, FLINT_BITS*k);
   mpz_fdiv_r_2exp(*A0, A, FLINT_BITS*2*k);
   mpz_mul_2exp(*temp, *R1, FLINT_BITS*2*k);
   mpz_add(*temp, *temp, *A0);
   mpz_mul_2exp(*temp2, *Q1, FLINT_BITS*k);
   mpz_mul(*temp2, *temp2, *B0);
   mpz_sub(*temp, *temp, *temp2);
   mpz_mul_2exp(*temp2, B, FLINT_BITS*k);
   
   while (mpz_cmp_ui(*temp, 0) < 0) 
   {
      mpz_sub_ui(*Q1, *Q1, 1);
      mpz_add(*temp, *temp, *temp2);
   }
   mpz_fdiv_q_2exp(*temp2, *temp, FLINT_BITS*k); 
   F_mpz_divrem_BZ(*Q0, *R0, *temp2, *B1);
   
   mpz_fdiv_r_2exp(*temp2, *temp, FLINT_BITS*k);
   mpz_mul_2exp(R, *R0, FLINT_BITS*k);
   mpz_add(R, R, *temp2);
   mpz_submul(R, *Q0, *B0);
   while (mpz_cmp_ui(R, 0) < 0) 
   {
      mpz_add(R, R, B);
   }
   
   F_mpz_release(); F_mpz_release(); F_mpz_release(); F_mpz_release();
   F_mpz_release(); F_mpz_release(); F_mpz_release(); F_mpz_release();
   F_mpz_release(); F_mpz_release(); F_mpz_release();
}

void F_mpz_mulmod_BZ(mpz_t res, mpz_t a, mpz_t b, mpz_t m)
{
   mpz_t * temp = F_mpz_alloc();
   
   mpz_mul(*temp, a, b);
   F_mpz_rem_BZ(res, *temp, m);
   
   F_mpz_release();
}


void F_mpz_expmod_BZ(mpz_t res, mpz_t a, mpz_t exp, mpz_t m)
{
   unsigned long n;
   unsigned long bits = mpz_sizeinbase(exp, 2);
   mpz_t aRED;
   mpz_t powRED;
   mpz_t temp;
   int flag = 0;
   
   mpz_init(aRED);
   mpz_init(powRED);
   mpz_init(temp);
   
   mpz_set(powRED, a);
#if DEBUG
   gmp_printf("powRED = %Zd\n", powRED);
#endif
   
   for (unsigned long i = 0; i < bits - 1; i++)
   {
      if (mpz_tstbit(exp, i))
      {
         if (flag) F_mpz_mulmod_BZ(res, res, powRED, m);
         else 
         {
            mpz_set(res, powRED);
            flag = 1;
         }
      }
      F_mpz_mulmod_BZ(powRED, powRED, powRED, m);
#if DEBUG
      gmp_printf("powRED = %Zd\n", powRED);
#endif
   }
   
   if (flag) F_mpz_mulmod_BZ(res, res, powRED, m);
   else mpz_set(res, powRED);
   
   mpz_clear(temp);
   mpz_clear(powRED);
   mpz_clear(aRED);
}

/*
   Large integer multiplication code
*/

void __F_mpz_mul(mpz_t res, mpz_t a, mpz_t b, unsigned long twk)
{
   unsigned long sa = mpz_size(a);
   unsigned long sb = mpz_size(b);

   if (sa+sb > FLINT_FFT_LIMBS_CROSSOVER) 
   {
      unsigned long s1 = (FLINT_BIT_COUNT(a->_mp_d[sa-1]) + FLINT_BIT_COUNT(b->_mp_d[sb-1]) <= FLINT_BITS);
   
      mp_limb_t* output = 
         (mp_limb_t*) flint_stack_alloc(sa + sb);
      __F_mpn_mul(output, a->_mp_d, sa, b->_mp_d, sb, twk);
      mpz_import(res, sa+sb-s1, -1, sizeof(mp_limb_t), 0, 0, output);
      if (mpz_sgn(res) != mpz_sgn(a)*mpz_sgn(b)) mpz_neg(res,res);
      flint_stack_release();
   } else mpz_mul(res, a, b);
}

void F_mpz_mul(mpz_t res, mpz_t a, mpz_t b)
{   
   unsigned long sa = mpz_size(a);
   unsigned long sb = mpz_size(b);

   if (sa+sb > FLINT_FFT_LIMBS_CROSSOVER) 
   {
      unsigned long s1 = (FLINT_BIT_COUNT(a->_mp_d[sa-1]) + FLINT_BIT_COUNT(b->_mp_d[sb-1]) <= FLINT_BITS);
      mp_limb_t* output = 
         (mp_limb_t*) flint_stack_alloc(sa + sb);
      F_mpn_mul(output, a->_mp_d, sa, b->_mp_d, sb);
      mpz_import(res, sa+sb-s1, -1, sizeof(mp_limb_t), 0, 0, output);
      if (mpz_sgn(res) != mpz_sgn(a)*mpz_sgn(b)) mpz_neg(res,res);
      flint_stack_release();
   } else mpz_mul(res, a, b);
}

