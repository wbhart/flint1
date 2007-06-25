/******************************************************************************

 flint.h
 Main header file for FLINT.

 (C) 2006 William Hart and David Harvey

******************************************************************************/

#include <limits.h>
#include <assert.h>
#include <gmp.h>
#include "longlong_wrapper.h"


#ifndef FLINT_H
#define FLINT_H


#define FLINT_INLINE inline


#if 0
#define FLINT_ASSERT assert
#else
#define FLINT_ASSERT(zzz_dummy)
#endif



#define FLINT_MAX(zzz1, zzz2) ((zzz1) > (zzz2) ? (zzz1) : (zzz2))
#define FLINT_MIN(zzz1, zzz2) ((zzz1) > (zzz2) ? (zzz2) : (zzz1))


/*
 FLINT_BITS is the number of bits per limb.
*/
#if ULONG_MAX == 4294967295U
#define FLINT_BITS 32
#define FLINT_LG_BITS_PER_LIMB 5
#define FLINT_BYTES_PER_LIMB 4
#define FLINT_LG_BYTES_PER_LIMB 2

#elif ULONG_MAX == 18446744073709551615U
#define FLINT_BITS 64
#define FLINT_LG_BITS_PER_LIMB 6
#define FLINT_BYTES_PER_LIMB 8
#define FLINT_LG_BYTES_PER_LIMB 3

#else
// only 32 and 64 bits are supported
#error FLINT requires that unsigned long is 32 bits or 64 bits
#endif

/* 
Cache hints to speed up reads from data in memory
*/
#if defined(__GNUC__) && __GNUC__ >= 3
#define FLINT_PREFETCH(addr,n) __builtin_prefetch((unsigned long*)addr+n,1,0) 
#elif defined(_MSC_VER) && _MSC_VER >= 1400
#define FLINT_PREFETCH(addr,n) PreFetchCacheLine(PF_TEMPORAL_LEVEL_1, (unsigned long*)addr+n)
#else
#define FLINT_PREFETCH(addr,n) /* nothing */
#endif

/*
Cache size in bytes.
*/
#define FLINT_CACHE_SIZE 65536

#define FLINT_POL_DIV_1_LENGTH 10

/*
If this flag is set, we are strictly only using GMP's documented interface.

If not set, we sometimes use hacks to speed things up (for example accessing
the mpz_t::_mp_d member directly).
 */
#define FLINT_GMP_COMPLIANT 0

#if FLINT_BITS == 32
#define half_ulong u_int16_t
#define half_long int16_t
#define HALF_FLINT_BITS 16
#else
#define half_ulong u_int32_t
#define half_long int32_t
#define HALF_FLINT_BITS 32
#endif

#define count_lead_zeros(a,b) \
   a = __builtin_clzll(b);

/*
Returns ceil(log2(x)).
If x == 0, returns 0.
*/
static inline unsigned long ceil_log2(unsigned long x)
{
   unsigned long result = 0;
   if (x)
   {
      x--;
      while (x)
      {
         x >>= 1;
         result++;
      }
   }
   return result;
}

/* Horrendously slow function for counting the number of bits a value takes up */

static inline unsigned long FLINT_BIT_COUNT(unsigned long x)
{
   unsigned long bits = FLINT_BITS;
   unsigned long mask = 1L<<(FLINT_BITS-1);
   while ((bits > 1) && !(mask & x)) 
   {
      bits--;
      mask>>=1;
   }
   return bits;
}

#endif


// end of file ****************************************************************
