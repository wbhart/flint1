/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   fmpz.h: "flat" multi-precision integer format

   Copyright (C) 2007, William Hart

*****************************************************************************/

#ifndef MPIR_FMPZ_H
#define MPIR_FMPZ_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <string.h>
#include <stdio.h>
#include <gmp.h>
#include "memory_manager.h"
#include "mpir.h"

#define fmpz_t __mpz_struct

#if MPIR_BITS == 64
#define TAB_START(xxx) ((fmpz_t *) ((uint64_t) xxx & (((uint64_t)-1)<<MPIR_LG_ALIGN)))
#else
#define TAB_START(xxx) ((fmpz_t *) ((uint32_t) xxx & (((uint32_t)-1)<<MPIR_LG_ALIGN)))
#endif

#define IMM_MAX (1UL<<(MPIR_BITS-2))-1UL

/* ==============================================================================

   Block memory management

===============================================================================*/

void fmpz_block_init(fmpz_t * entry, ulong m);

void fmpz_block_init2(fmpz_t * entry, ulong m, ulong n);

void fmpz_block_realloc(fmpz_t * entry, ulong n);

/*
   If n is greater than the current value of n for the block, reallocate
   the block so that the integers in the block have n limbs plus a 
   sign/size limb
*/

static inline 
int fmpz_block_fit_limbs(fmpz_t * entry, ulong n)
{
   if UNLIKELY(n > entry->_mp_alloc) 
   {
      fmpz_block_realloc(entry, n);
      return 1;
   } else return 0;
}

void fmpz_block_clear(fmpz_t * entry);

/* ==============================================================================

   fmpz_t memory management

===============================================================================*/

fmpz_t * fmpz_init_array(ulong count);

static inline
fmpz_t * fmpz_init(void)
{
   return fmpz_init_array(1L);
}

void fmpz_clear_array(fmpz_t * f, ulong count);

static inline
void fmpz_clear(fmpz_t * f)
{
   fmpz_clear_array(f, 1L);
}

/* 
   Make the integer f (and all others in the same block) have space
   for n limbs
*/

static inline
void fmpz_realloc(fmpz_t * f, ulong n)
{
   fmpz_block_realloc(TAB_START(f), n);
}

fmpz_t * fmpz_realloc_array(fmpz_t * arr, ulong old_count, ulong count);

/* 
   Make the integer f (and all others in the same block) have space 
   for at least n limbs
   Integers will not be shrunk if the new n is smaller than the current 
   value of n for the block
*/

static inline
int fmpz_fit_limbs(fmpz_t * f, ulong n)
{
   return fmpz_block_fit_limbs(TAB_START(f), n);
}

/* ==============================================================================

   Properties

===============================================================================*/

static inline 
mp_limb_t * fmpz_data(fmpz_t * in)
{
   return in->_mp_d;
}

static inline
ulong fmpz_size(fmpz_t * in)
{
   ulong alloc = in->_mp_alloc;
   if (alloc > 1L) return MPIR_ABS(in->_mp_size);
   if ((long) in->_mp_d == 0L) return 0L;
   return 1L;
}

static inline
ulong fmpz_bits(fmpz_t * in)
{
   if (in->_mp_alloc == 1L)
   {
      long i_int = (long) in->_mp_d;
      if (i_int) return MPIR_BIT_COUNT(MPIR_ABS(i_int));
      else return 0L;
   }
   ulong limbs = fmpz_size(in);
   mp_limb_t * dp = in->_mp_d;
   if (limbs)
   {
      return ((limbs-1)<<MPIR_LG_BITS) + MPIR_BIT_COUNT(dp[limbs-1]);
   } else return 0;
}

/* ==============================================================================

   Conversion

===============================================================================*/

void mpz_to_fmpz(fmpz_t * fnum, mpz_t num);

void fmpz_to_mpz(mpz_t num, fmpz_t * fnum);

extern double __gmpn_get_d(mp_limb_t * fp, ulong limbs, ulong size, long exp);

void gmpz_set(mpz_ptr w, mpz_srcptr u);

double fmpz_get_d(fmpz_t * f);

double fmpz_get_d_2exp(long * exp, fmpz_t * f);

/* ==============================================================================

   Random generation

===============================================================================*/

void fmpz_random(fmpz_t * f, ulong bits);

void fmpz_randomm(fmpz_t * out, fmpz_t * in);

/* ==============================================================================

   Read/print

===============================================================================*/

void fmpz_print(fmpz_t * in);

void fmpz_fread(fmpz_t * in, FILE * f);

static inline
void fmpz_read(fmpz_t * in)
{
   fmpz_fread(in, stdin);
}

/* ==============================================================================

   Number theoretical

===============================================================================*/

int fmpz_probab_prime_p(fmpz_t * p, ulong n);

/* ==============================================================================

   Set/get

===============================================================================*/

static inline
void fmpz_zero(fmpz_t * f)
{
   if (f->_mp_alloc == 1L) f->_mp_d = (mp_limb_t *) 0L;
   else f->_mp_size = 0L;
}

void fmpz_set(fmpz_t * out, fmpz_t * f);

/*
   Set res to the unsigned long x
*/

static inline
void fmpz_set_ui(fmpz_t * res, ulong x)
{
   if (res->_mp_alloc > 1L)
   { 
      if (x == 0L) res->_mp_size = 0L;
      else
      {
         res->_mp_size = 1L;
         res->_mp_d[0] = x;
      }
   } else if (x > IMM_MAX)
   {
      fmpz_fit_limbs(res, 2L);
      res->_mp_size = 1L;
      res->_mp_d[0] = x;      
   } else res->_mp_d = (mp_limb_t *) x;
}

/*
   Return the least significant limb of the absolute value of x
*/

static inline
unsigned long fmpz_get_ui(fmpz_t * res)
{
   ulong alloc = res->_mp_alloc;
   if (alloc > 1L) return mpz_get_ui(res);

   long r_int = (long) res->_mp_d;
   if (r_int > 0L) return r_int;
   else return -r_int;
}

/*
   Set res to the signed long x
*/

static inline
void fmpz_set_si(fmpz_t * res, long x)
{
   if (res->_mp_alloc > 1L)
   {  
      if (x < 0L)
      {
         res->_mp_size = -1L;
         res->_mp_d[0] = -x;
      } else if (x > 0L)
      {
         res->_mp_size = 1L;
         res->_mp_d[0] = x;
      } else res->_mp_size = 0L;
   } else if (MPIR_ABS(x) > IMM_MAX)
   {
      fmpz_fit_limbs(res, 2L);
      if (x < 0L)
      {
         res->_mp_size = -1L;
         res->_mp_d[0] = -x;
      } else if (x > 0L)
      {
         res->_mp_size = 1L;
         res->_mp_d[0] = x;
      } else res->_mp_size = 0L;      
   } else res->_mp_d = (mp_limb_t *) x;
}

/*
   Return the least significant limb of x as a signed quantity
*/

static inline
long fmpz_get_si(fmpz_t * res)
{
   if (res->_mp_alloc > 1L) return mpz_get_si(res);

   return (long) res->_mp_d;
}

/* ==============================================================================

   Negation

===============================================================================*/

void fmpz_neg(fmpz_t * out, fmpz_t * f);

/* ==============================================================================

   Comparison

===============================================================================*/

static inline
int fmpz_is_zero(fmpz_t * x)
{
   if (x->_mp_alloc == 1L) return (((long) x->_mp_d) == 0L);
   else return (x->_mp_size == 0L);
}

int fmpz_equal(fmpz_t * f2, fmpz_t * f1);

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

void _fmpz_add_IMM(fmpz_t * out, fmpz_t * f1, long c);

void fmpz_add(fmpz_t * out, fmpz_t * f1, fmpz_t * f2);

void _fmpz_sub_IMM(fmpz_t * out, fmpz_t * f1, long c);

void fmpz_sub(fmpz_t * out, fmpz_t * f1, fmpz_t * f2);

/* ==============================================================================

   Addmul/submul

===============================================================================*/

void fmpz_addmul_ui (fmpz_t * w, fmpz_t * x, ulong y);

void fmpz_submul_ui(fmpz_t * w, fmpz_t * x, ulong y);

/* ==============================================================================

   Multiplication

===============================================================================*/

void fmpz_mul_2exp(fmpz_t * w, fmpz_t * u, ulong exp);


#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
