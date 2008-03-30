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
/******************************************************************************

 mpir.h
 Main header file for MPIR.

 (C) 2006 William Hart and David Harvey

******************************************************************************/

#ifndef MPIR_H
#define MPIR_H

#include <limits.h>
#include <assert.h>
#include <gmp.h>
#include <math.h>
#include <stdint.h>
#include "longlong_wrapper.h"

#ifdef __cplusplus
 extern "C" {
#endif

#define ulong unsigned long

#if 0
#define MPIR_ASSERT assert
#else
#define MPIR_ASSERT(zzz_dummy)
#endif

#define MPIR_MAX(zzz1, zzz2) ((zzz1) > (zzz2) ? (zzz1) : (zzz2))
#define MPIR_MIN(zzz1, zzz2) ((zzz1) > (zzz2) ? (zzz2) : (zzz1))
#define MPIR_ABS(zzz) ((long)(zzz) < 0 ? (-zzz) : (zzz))

/*
 MPIR_BITS is the number of bits per limb.
*/
#if ULONG_MAX == 4294967295U
#define MPIR_BITS 32
#define MPIR_D_BITS 32
#define MPIR_LG_BITS 5
#define MPIR_BYTES 4
#define MPIR_LG_BYTES 2
#define MPIR_ALIGN 8
#define MPIR_LG_ALIGN 3

#elif ULONG_MAX == 18446744073709551615U
#define MPIR_BITS 64
#define MPIR_D_BITS 53
#define MPIR_LG_BITS 6
#define MPIR_BYTES 8
#define MPIR_LG_BYTES 3
#define MPIR_ALIGN 64
#define MPIR_LG_ALIGN 6

#else
// only 32 and 64 bits are supported
#error MPIR requires that ulong is 32 bits or 64 bits
#endif

#define MPIR_BLOCK MPIR_ALIGN
#define MPIR_LG_BLOCK MPIR_LG_ALIGN

#if defined(__GNUC__) 
#if MPIR_BITS == 64
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
#error Currently MPIR only compiles with GCC
#endif

static inline ulong MPIR_BIT_COUNT(ulong x)
{
   ulong zeros = 0L;
   if (x) count_lead_zeros(zeros, x);
   return MPIR_BITS-zeros;
}

/* 
   On some platforms arithmetic shifts by MPIR_BITS don't yield all zeros 
   So we define these macros for use in situations where this would be a problem
*/

static inline ulong MPIR_RSHIFT(ulong in, ulong shift)
{
   if (shift == MPIR_BITS) return 0L;
   return (in>>shift);
}

static inline ulong MPIR_LSHIFT(ulong in, ulong shift)
{
   if (shift == MPIR_BITS) return 0L;
   return (in<<shift);
}

#ifdef __cplusplus
 }
#endif

#endif

// end of file ****************************************************************
