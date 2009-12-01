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

 flint.h
 Main header file for FLINT.

 (C) 2006 William Hart and David Harvey

******************************************************************************/

#ifndef FLINT_H
#define FLINT_H

#include <limits.h>
#include <assert.h>
#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include "longlong_wrapper.h"

#ifdef __cplusplus
 extern "C" {
#endif

#if 0
#define FLINT_ASSERT assert
#else
#define FLINT_ASSERT(zzz_dummy)
#endif

#ifndef __USE_ISOC99
#define round(x) floor(x + 0.5)
#endif


#define FLINT_MAX(zzz1, zzz2) ((zzz1) > (zzz2) ? (zzz1) : (zzz2))
#define FLINT_MIN(zzz1, zzz2) ((zzz1) > (zzz2) ? (zzz2) : (zzz1))
#define FLINT_ABS(zzz) ((long)(zzz) < 0 ? (-zzz) : (zzz))


/*
 FLINT_BITS is the number of bits per limb.
*/
#if ULONG_MAX == 4294967295U
#define FLINT_BITS 32
#define FLINT_D_BITS 32
#define FLINT_LG_BITS_PER_LIMB 5
#define FLINT_BYTES_PER_LIMB 4
#define FLINT_LG_BYTES_PER_LIMB 2

#elif ULONG_MAX == 18446744073709551615U
#define FLINT_BITS 64
#define FLINT_D_BITS 53
#define FLINT_LG_BITS_PER_LIMB 6
#define FLINT_BYTES_PER_LIMB 8
#define FLINT_LG_BYTES_PER_LIMB 3

#else
// only 32 and 64 bits are supported
#error FLINT requires that unsigned long is 32 bits or 64 bits
#endif

/* 
Cache/branch hints to speed up reads from data in memory/branches
*/
#if defined(__GNUC__) 
#define FLINT_PREFETCH(addr,n) __builtin_prefetch((unsigned long*)addr+n,1,0) 
#define LIKELY(cond)    (__builtin_expect ((cond) != 0, 1))
#define UNLIKELY(cond)  (__builtin_expect ((cond) != 0, 0))
#elif defined(_MSC_VER) && _MSC_VER >= 1400
#define FLINT_PREFETCH(addr,n) PreFetchCacheLine(PF_TEMPORAL_LEVEL_1, (unsigned long*)addr+n)
#define LIKELY(cond)
#define UNLIKELY(cond)
#else
#define FLINT_PREFETCH(addr,n)
#define LIKELY(cond) (cond)
#define UNLIKELY(cond) (cond)
#endif

/*
Thread stuff
*/

#define THREAD

#ifdef FLINT_TEST_SUPPORT_H 
#define FLINT_THREAD_CLEANUP \
	do { \
      if (omp_get_thread_num() != 0) \
	   { \
		   flint_stack_cleanup(); \
	      flint_test_support_cleanup(); \
		} \
	} while (0);
#else
#define FLINT_THREAD_CLEANUP \
	do { \
      if (omp_get_thread_num() != 0) \
	   { \
		   flint_stack_cleanup(); \
      } \
	} while (0);
#endif 

/*
Cache size in bytes.
*/
#define FLINT_CACHE_SIZE 65536

#define FLINT_POL_DIV_1_LENGTH 10

#define ulong unsigned long

#if FLINT_BITS == 32
#define half_ulong uint16_t
#define half_long int16_t
#define HALF_FLINT_BITS 16
#else
#define half_ulong uint32_t
#define half_long int32_t
#define HALF_FLINT_BITS 32
#endif

#if defined(__GNUC__) 
#if FLINT_BITS == 64
#define count_lead_zeros(a,b) \
   a = __builtin_clzll(b);
#define count_trail_zeros(a,b) \
   a = __builtin_ctzll(b);
#else
#define count_lead_zeros(a,b) \
   a = __builtin_clzl(b);
#define count_trail_zeros(a,b) \
   a = __builtin_ctzl(b);
#endif
#else
#error Currently FLINT only compiles with GCC
#endif

/* 
   On some platforms arithmetic shifts by FLINT_BITS don't yield all zeros 
   So we define these macros for use in situations where this would be a problem
*/

static inline unsigned long r_shift(unsigned long in, unsigned long shift)
{
   if (shift == FLINT_BITS) return 0L;
   return (in>>shift);
}

static inline unsigned long l_shift(unsigned long in, unsigned long shift)
{
   if (shift == FLINT_BITS) return 0L;
   return (in<<shift);
}

static inline unsigned long FLINT_BIT_COUNT(unsigned long x)
{
   unsigned long zeros = FLINT_BITS;
   if (x) count_lead_zeros(zeros, x);
   return FLINT_BITS - zeros;
}

/*
Returns ceil(log2(x)).
If x == 0, returns 0.
*/
static inline unsigned long ceil_log2(unsigned long x)
{
   unsigned long result = 0;
   if (x > 1)
   {
      x--;
      result = FLINT_BIT_COUNT(x);
   }
   return result;
}



#ifdef __cplusplus
 }
#endif

#endif

// end of file ****************************************************************
