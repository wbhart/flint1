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

#define TOP_MASK (1L<<(MPIR_BITS-1))
#if ULONG_MAX != ((size_t)-1) // For LLP64 machines
#define ENTRY_PTR(xxx) ((table_entry*) ((uint64_t) xxx & (((uint64_t)-1)<<MPIR_LG_BLOCK)))
#else
#define ENTRY_PTR(xxx) ((table_entry*) ((ulong) xxx & (((ulong)-1)<<MPIR_LG_BLOCK)))
#endif

#define OFFSET(xxx) ((ulong) xxx & ((1L<<MPIR_LG_BLOCK)-1L))

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
   ulong n;
   ulong count; 
   block * pad1;
   block * pad2;
   block * pad3;
   block * pad4;
   block * pad5;
   block * pad6;
} table_entry;
#else // For LP64 and ILP32 machines
typedef struct
{
   block * block_ptr;
   ulong n;
   ulong count; 
   ulong pad1; 
   ulong pad2; 
   ulong pad3; 
   ulong pad4; 
   ulong pad5; 
} table_entry;
#endif

/* ==============================================================================

   Block memory management

===============================================================================*/

void fmpz_block_init(table_entry * entry);

void fmpz_block_init2(table_entry * entry, ulong n);

void fmpz_block_init_small(table_entry * entry, ulong m);

void fmpz_block_init2_small(table_entry * entry, ulong m, ulong n);

void fmpz_block_realloc(table_entry * entry, ulong n);

/*
   If n is greater than the current value of n for the block, reallocate
   the block so that the integers in the block have n limbs plus a 
   sign/size limb
*/

static inline 
int fmpz_block_fit_limbs(table_entry * entry, ulong n)
{
   if (n > entry->n) 
   {
      fmpz_block_realloc(entry, n);
      return 1;
   } else return 0;
}

fmpz_t * fmpz_realloc_array(fmpz_t * arr, ulong old_count, ulong count);

void fmpz_block_clear(table_entry * entry);

/* ==============================================================================

   fmpz_t memory management

===============================================================================*/

fmpz_t * fmpz_init_array(ulong count);

fmpz_t * fmpz_init(void);

void fmpz_clear(fmpz_t * f);

void fmpz_clear_array(fmpz_t * f, ulong count);

/* 
   Make the integer f (and all others in the same block) have space
   for n limbs
*/

static inline
void fmpz_realloc(fmpz_t* f, ulong n)
{
   fmpz_block_realloc(ENTRY_PTR(f), n);
}

/* 
   Make the integer f (and all others in the same block) have space 
   for at least n limbs
   Integers will not be shrunk if the new n is smaller than the current 
   value of n for the block
*/

static inline
int fmpz_fit_limbs(fmpz_t* f, ulong n)
{
   return fmpz_block_fit_limbs(ENTRY_PTR(f), n);
}

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
   if (ENTRY_PTR(in)->n) return MPIR_ABS(fmpz_data(in)[0]);
   else return 0;
}

/* ==============================================================================

   Conversion

===============================================================================*/

void mpz_to_fmpz(fmpz_t * res, mpz_t x);

void fmpz_to_mpz(mpz_t res, fmpz_t * x);

/* ==============================================================================

   Read/print

===============================================================================*/

void fmpz_print(fmpz_t * in);

void fmpz_fread(fmpz_t * in, FILE * f);

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

   Negation

===============================================================================*/

void fmpz_neg(fmpz_t * out, fmpz_t * f);

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
