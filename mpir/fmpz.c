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

   fmpz.c: "flat" integer format

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "fmpz.h"
#include "mpir.h"
#include "memory_manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "mpn_extras.h"

#define SWAP_PTRS(xxx, yyy) \
do { \
   fmpz_t * temp_ptr_zzz = yyy; \
   yyy = xxx; \
   xxx = temp_ptr_zzz; \
} while (0);

#define SWAP_SIZES(xxx, yyy) \
do { \
   ulong temp_size_zzz = yyy; \
   yyy = xxx; \
   xxx = temp_size_zzz; \
} while (0);

#define NORM(coeff) \
do { \
   if ((coeff)[0]) \
   { \
      if ((long) (coeff)[0] < 0) \
      { \
         while ((!(coeff)[-(coeff)[0]]) && (coeff)[0]) (coeff)[0]++; \
      } else \
      { \
         while ((!(coeff)[(coeff)[0]]) && (coeff)[0]) (coeff)[0]--; \
      } \
   } \
} while (0);

/* ==============================================================================

   Block memory management

===============================================================================*/

/*
   Initialises a block, but doesn't actually allocate any memory
   fmpz_block_realloc and fmpz_block_clear can be called on the block
*/

void fmpz_block_init(table_entry * entry)
{
   entry->n = 0L;
}

/*
   Initialises a block, allocating space for MPIR_BLOCK integers of size n
*/

void fmpz_block_init2(table_entry * entry, long n)
{
   entry->n = n;
  /* 
     Block has MPIR_BLOCK integers, n+1 limbs per integer, MPIR_BYTES bytes per limb
  */
   if (n) entry->block_ptr = (block *) mpir_alloc(((n+1)<<MPIR_LG_BLOCK)<<MPIR_LG_BYTES);
}

/*
   Initialises a block for storing a single integer
*/

void fmpz_block_init_single(table_entry * entry)
{
   entry->n = (1L | TOP_MASK);
   entry->block_ptr = (block *) mpir_alloc(2<<MPIR_LG_BYTES);
   entry->block_ptr[0] = 0L;
}

/*
   Allocates a block with a single integer of n limbs plus a sign/size limb
*/

void fmpz_block_init2_single(table_entry * entry, long n)
{
   if (!n) n++;
   entry->n = (n | TOP_MASK);
   entry->block_ptr = (block *) mpir_alloc((n+1)<<MPIR_LG_BYTES);
   entry->block_ptr[0] = 0L;
}

/*
   Reallocate the integers in a block to have n limbs plus a sign/size limb
*/

void fmpz_block_realloc(table_entry * entry, long n)
{
   long oldn = entry->n;
   
   if (n == 0L)
   {
      fmpz_block_clear(entry);
      if (oldn < 0L) fmpz_block_init_single(entry);
      else fmpz_block_init(entry);
      return;
   }
   
   if (oldn < 0L)
   {
      long oldn = (entry->n ^ TOP_MASK);
      if (oldn <= n)
      {
         entry->block_ptr = (block *) mpir_realloc(entry->block_ptr, (n+1)<<MPIR_LG_BYTES);
      } else
      {
         block * new_block = (block *) mpir_alloc((n+1)<<MPIR_LG_BYTES);
         F_mpn_copy(new_block, entry->block_ptr, n+1);
         free(entry->block_ptr);
         entry->block_ptr = new_block;
      }
      entry->n = (n | TOP_MASK);
   } else
   {
      if (oldn == 0L)
      {
         entry->block_ptr = (block *) mpir_alloc(((n+1)<<MPIR_LG_BLOCK)<<MPIR_LG_BYTES);
         // Now zero all the integers in the block
         mp_limb_t * ptr = (mp_limb_t*) entry->block_ptr;
         for (ulong i = 0; i < MPIR_BLOCK; i++)
         {
            ptr[i*(n+1)] = 0L;
         }
      } else if (oldn <= n)
      {
         entry->block_ptr = (block *) mpir_realloc(entry->block_ptr, ((n+1)<<MPIR_LG_BLOCK)<<MPIR_LG_BYTES);
         mp_limb_t * ptr = entry->block_ptr + (MPIR_BLOCK - 1)*(oldn+1);
         mp_limb_t * ptr2 = entry->block_ptr + (MPIR_BLOCK - 1)*(n+1);
         for (int i = MPIR_BLOCK - 1; i > 0; i--, ptr -= (oldn+1), ptr2 -= (n+1))
         {
            F_mpn_copy(ptr2, ptr, oldn+1);
         }
      } else
      {
         mp_limb_t * ptr = entry->block_ptr + (MPIR_BLOCK - 1)*(oldn+1);
         block * new_block = (block *) mpir_alloc(((n+1)<<MPIR_LG_BLOCK)<<MPIR_LG_BYTES);
         mp_limb_t * ptr2 = new_block + (MPIR_BLOCK - 1)*(n+1);
         for (int i = MPIR_BLOCK - 1; i >= 0; i--, ptr -= (oldn+1), ptr2 -= (n+1))
         {
            F_mpn_copy(ptr2, ptr, n+1);
         }
         free(entry->block_ptr);
         entry->block_ptr = new_block;
      }

      entry->n = n;
   }
}

/*
   If n is greater than the current value of n for the block, reallocate
   the block so that the integers in the block have n limbs plus a 
   sign/size limb
*/

void fmpz_block_fit_limbs(table_entry * entry, long n)
{
   long oldn = (entry->n & ~TOP_MASK);
   if (n > oldn) fmpz_block_realloc(entry, n);
}

/*
   Clear a block (whether single or not) of integers
*/

void fmpz_block_clear(table_entry * entry)
{
   if (entry->n & ~TOP_MASK) free(entry->block_ptr);
}

/* ==============================================================================

   fmpz_t memory management

===============================================================================*/

/*
   Initialise an array of count integers and return a pointer to the first
*/

fmpz_t * fmpz_init_array(ulong count)
{
   ulong blocks = ((count-1)>>MPIR_LG_BLOCK) + 1;
   table_entry* tab = (table_entry*) mpir_aligned_alloc(blocks<<MPIR_LG_BLOCK);
   for (ulong i = 0; i < blocks; i++)
   {
      fmpz_block_init(tab + i);  
   }
   return (fmpz_t*) tab;
}

/*
   Initialise a single integer
*/

fmpz_t * fmpz_init(void)
{
   table_entry * entry = (table_entry*) mpir_aligned_alloc(MPIR_BLOCK);
   fmpz_block_init_single(entry);
   return (fmpz_t*) entry;
}

/*
   Clear a single integer
   This must not be called on individual integers in an array of integers
*/

void fmpz_clear(fmpz_t * f)
{
   fmpz_block_clear((table_entry*) f);
   mpir_aligned_free((void*) f);
}

/*
   Clear the array of count integers where f is a pointer to the first
   Partial arrays cannot be cleared. The full array must be cleared
   with count set to the current number of integers in the array
*/

void fmpz_clear_array(fmpz_t * f, ulong count)
{
   ulong blocks = ((count-1)>>MPIR_LG_BLOCK) + 1;
   table_entry * ptr = (table_entry*) f;
   for (ulong i = 0; i < blocks; i++)
   {
      fmpz_block_clear(ptr + i);  
   }  
   mpir_aligned_free((void*) f);
}

/* 
   Make the integer f (and all others in the same block) have space
   for n limbs
*/

void fmpz_realloc(fmpz_t* f, ulong n)
{
   fmpz_block_realloc(ENTRY_PTR(f), n);
}

/*
   Realloc an array of fmpz_t's from length old_count to count
*/

fmpz_t * fmpz_realloc_array(fmpz_t * arr, ulong old_count, ulong count)
{
   ulong blocks = ((count-1)>>MPIR_LG_BLOCK) + 1;
   ulong old_blocks = ((old_count-1)>>MPIR_LG_BLOCK) + 1;
   if (blocks > old_blocks)
   {
      table_entry * tab = (table_entry *) mpir_aligned_realloc(arr, blocks<<MPIR_LG_BLOCK); 
      for (ulong i = old_blocks; i < blocks; i++)
      {
         fmpz_block_init(tab + i);  
      }
      return (fmpz_t *) tab;
   } else if (blocks < old_blocks)
   {
      table_entry * tab = (table_entry *) arr; 
      for (ulong i = blocks; i < old_blocks; i++)
      {
         fmpz_block_clear(tab + i);  
      }
      return (fmpz_t *) mpir_aligned_realloc(arr, blocks<<MPIR_LG_BLOCK);
   }

   return arr;
}

/* 
   Make the integer f (and all others in the same block) have space 
   for at least n limbs
   Integers will not be shrunk if the new n is smaller than the current 
   value of n for the block
*/

void fmpz_fit_limbs(fmpz_t* f, ulong n)
{
   fmpz_block_fit_limbs(ENTRY_PTR(f), n);
}

/* ==============================================================================

   Conversion

===============================================================================*/

void mpz_to_fmpz(fmpz_t * res, mpz_t x)
{
   if (mpz_sgn(x))
   {
      fmpz_fit_limbs(res, mpz_size(x));
      mp_limb_t * rp = fmpz_data(res);
      size_t countp;
      mpz_export(rp + 1, &countp, -1, sizeof(mp_limb_t), 0, 0, x);
      rp[0] = ((long) mpz_sgn(x) > 0L) ? (long) countp : (long) -countp;
   }
   else
   {
      fmpz_fit_limbs(res, 1);
      mp_limb_t * rp = fmpz_data(res);
      rp[0] = 0L;
   }
}

void fmpz_to_mpz(mpz_t res, fmpz_t * x)
{
   mp_limb_t * xp = fmpz_data(x);
   long ss = xp[0];
   
   if (ss)
   {
      mpz_import(res, MPIR_ABS(ss), -1, sizeof(mp_limb_t), 0, 0, xp + 1);

      if (ss < 0L)
         mpz_neg(res, res);
   } else
      mpz_set_ui(res, 0L);
   
}

/* ==============================================================================

   Set/get

===============================================================================*/

void fmpz_set(fmpz_t * out, fmpz_t * f)
{
   unsigned long d1n;

   if (out != f)
   {
      ulong dn = MPIR_ABS(fmpz_data(f)[0]);
      fmpz_fit_limbs(out, MPIR_MAX(dn, 1L));
      mp_limb_t * rp = fmpz_data(out);
      mp_limb_t * dp = fmpz_data(f);
      F_mpn_copy(rp, dp, dn + 1); 
   }
}

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

void fmpz_add(fmpz_t * out, fmpz_t * f1, fmpz_t * f2)
{
   long cy;
   unsigned long d1n = MPIR_ABS(fmpz_data(f1)[0]);
   unsigned long d2n = MPIR_ABS(fmpz_data(f2)[0]);
   
   if (d1n < d2n) 
   {
      SWAP_PTRS(f1, f2);
      SWAP_SIZES(d1n, d2n);
   } 

   fmpz_fit_limbs(out, MPIR_MAX(d1n, 1L));
   mp_limb_t * rp = fmpz_data(out);
   mp_limb_t * d1p = fmpz_data(f1);
   mp_limb_t * d2p = fmpz_data(f2); 

   if (!d2n)
   {
      if (!d1n) rp[0] = 0L;
      else
      {
         if (rp != d1p) F_mpn_copy(rp, d1p, d1n+1);
      }
   } else if ((long) (d1p[0] ^ d2p[0]) >= 0L)
   {
      rp[0] = d1p[0];
      cy = mpn_add(rp+1, d1p+1, d1n, d2p+1, d2n);
      if (cy) 
      {
         fmpz_fit_limbs(out, d1n + 1);
         rp = fmpz_data(out);
         rp[d1n+1] = cy;
         if ((long) rp[0] < 0L) rp[0]--;
         else rp[0]++;
      }
   } else
   {
      cy = 0;
      if (d1n != d2n) cy = 1;
      else cy = mpn_cmp(d1p+1, d2p+1, d1n); 
          
      if (cy == 0) rp[0] = 0L;
      else if (cy > 0) 
      {
         mpn_sub(rp+1, d1p+1, d1n, d2p+1, d2n);
         rp[0] = d1p[0];
         NORM(rp);
      }
      else
      {
         mpn_sub_n(rp+1, d2p+1, d1p+1, d1n);
         rp[0] = -d1p[0];
         NORM(rp);
      }
   }
}

void fmpz_sub(fmpz_t * out, fmpz_t * f1, fmpz_t * f2)
{
   long cy;
   unsigned long d1n = MPIR_ABS(fmpz_data(f1)[0]);
   unsigned long d2n = MPIR_ABS(fmpz_data(f2)[0]);
   
   int in_order = 1;
   if (d1n < d2n) 
   {
      SWAP_PTRS(f1, f2);
      SWAP_SIZES(d1n, d2n);
      in_order = 0;
   } 

   fmpz_fit_limbs(out, MPIR_MAX(d1n, 1L));
   mp_limb_t * rp = fmpz_data(out);
   mp_limb_t * d1p = fmpz_data(f1);
   mp_limb_t * d2p = fmpz_data(f2); 
   
   if (!d2n)
   {
      if (!d1n) rp[0] = 0L;
      else
      {
         if (d1p != rp) F_mpn_copy(rp, d1p, d1n+1);
         if (!in_order) rp[0] = -rp[0];
      }
   } else if ((long) (d1p[0] ^ d2p[0]) < 0)
   {
      if (in_order) rp[0] = d1p[0];
      else rp[0] = -d1p[0];
      cy = mpn_add(rp+1, d1p+1, d1n, d2p+1, d2n);
      if (cy) 
      {
         fmpz_fit_limbs(out, d1n+1);
         rp = fmpz_data(out);
         rp[d1n+1] = cy;
         if ((long) rp[0] < 0) rp[0]--;
         else rp[0]++;
      }
   } else
   {
      cy = 0;
      if (d1n != d2n) cy = 1;
      else cy = mpn_cmp(d1p+1, d2p+1, d1n); 
          
      if (cy == 0) rp[0] = 0L;
      else if (cy > 0) 
      {
         mpn_sub(rp+1, d1p+1, d1n, d2p+1, d2n);
         if (in_order) rp[0] = d1p[0];
         else rp[0] = -d1p[0];
         NORM(rp);
      }
      else
      {
         mpn_sub_n(rp+1, d2p+1, d1p+1, d1n);
         if (in_order) rp[0] = -d1p[0];
         else rp[0] = d1p[0];
         NORM(rp);
      }
   }
}

// *************** end of file
