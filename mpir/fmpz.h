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
 
#include <gmp.h>
#include "memory_manager.h"
#include "mpir.h"

#define TOP_MASK (1L<<(MPIR_BITS-1))
#if ULONG_MAX != ((size_t)-1) // For LLP64 machines
#define ENTRY_PTR(xxx) ((table_entry*) ((uint64_t) xxx & (((uint64_t)-1)<<MPIR_LG_BLOCK)))
#else
#define ENTRY_PTR(xxx) ((table_entry*) ((ulong) xxx & (((ulong)-1)<<MPIR_LG_BLOCK)))
#endif

#define OFFSET(xxx) ((ulong) xxx & ((1L<<MPIR_LG_BLOCK)-1L))

#define FMPZ_DATA(xxx) ENTRY_PTR(xxx)->n

/*
   Returns a pointer to the actual data in the given fmpz_t
   The first limb is a sign/size limb, the remainder is the 
   absolute value of the integer
*/

typedef char fmpz_t; // fmpz_t * will have char * arithmetic, but will never be
                    // dereferenced directly

typedef mp_limb_t block; // A block is a set of MPIR_BLOCK integers contiguously in memory

#if ULONG_MAX != ((size_t)-1) // For LLP64 machines
typedef struct
{
   block * block_ptr;
   long n;
   ulong pad; 
} table_entry;
#else // For LP64 and ILP32 machines
typedef struct
{
   block * block_ptr;
   long n;
} table_entry;
#endif

/* ==============================================================================

   Block memory management

===============================================================================*/

void fmpz_block_init(table_entry * entry);

void fmpz_block_init2(table_entry * entry, long n);

void fmpz_block_init_single(table_entry * entry);

void fmpz_block_init2_single(table_entry * entry, long n);

void fmpz_block_realloc(table_entry * entry, long n);

void fmpz_block_fit_limbs(table_entry * entry, long n);

fmpz_t * fmpz_realloc_array(fmpz_t * arr, ulong old_count, ulong count);

void fmpz_block_clear(table_entry * entry);

/* ==============================================================================

   fmpz_t memory management

===============================================================================*/

fmpz_t * fmpz_init_array(ulong count);

fmpz_t * fmpz_init(void);

void fmpz_clear(fmpz_t * f);

void fmpz_clear_array(fmpz_t * f, ulong count);

void fmpz_realloc(fmpz_t* f, ulong n);

void fmpz_fit_limbs(fmpz_t* f, ulong n);

/* ==============================================================================

   Properties

===============================================================================*/

static inline 
mp_limb_t * fmpz_data(fmpz_t * int_in)
{
   table_entry * entry = ENTRY_PTR(int_in);
   return entry->block_ptr + (entry->n+1) * OFFSET(int_in);
}

static inline
ulong fmpz_size(fmpz_t * in)
{
   return MPIR_ABS(fmpz_data(in)[0]);
}

/* ==============================================================================

   Conversion

===============================================================================*/

void mpz_to_fmpz(fmpz_t * res, mpz_t x);

void fmpz_to_mpz(mpz_t res, fmpz_t * x);

/* ==============================================================================

   Set/get

===============================================================================*/

void fmpz_set(fmpz_t * out, fmpz_t * f);

/*
   Set res to the unsigned long x
*/

static inline
void fmpz_set_ui(fmpz_t * res, unsigned long x)
{
   fmpz_fit_limbs(res, 1);
   
   mp_limb_t * data = fmpz_data(res);

   if (x) 
   {
      data[0] = 1UL;
      data[1] = x;
   }
   else
      data[0] = 0UL;
}

/*
   Return the least significant limb of the absolute value of x
*/

static inline
unsigned long fmpz_get_ui(fmpz_t * res)
{
   mp_limb_t * data = fmpz_data(res);

   if (data[0]) return data[1];
   else return 0UL;
}

/* ==============================================================================

   Comparison

===============================================================================*/

static inline
int fmpz_is_zero(fmpz_t * x)
{
   return (fmpz_data(x)[0] == 0L);
}

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

void fmpz_add(fmpz_t * out, fmpz_t * f1, fmpz_t * f2);

void fmpz_sub(fmpz_t * out, fmpz_t * f1, fmpz_t * f2);


#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
