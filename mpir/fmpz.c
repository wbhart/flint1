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

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>

#include "fmpz.h"
#include "mpir.h"
#include "memory_manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "mpn_extras.h"

#define SWAP_FMPZ_PTRS(xxx, yyy) \
do { \
   fmpz_t * temp_ptr_zzz = yyy; \
   yyy = xxx; \
   xxx = temp_ptr_zzz; \
} while (0);

#define SWAP_LIMB_PTRS(xxx, yyy) \
do { \
   mp_limb_t * temp_ptr_zzz = yyy; \
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
   Initialises a block with MPIR_BLOCK integers in it, but 
   doesn't actually allocate any memory fmpz_block_realloc and 
   fmpz_block_clear can be called on the block
*/

void fmpz_block_init(table_entry * entry)
{
   entry->n = 0L;
   entry->count = MPIR_BLOCK;
}

/*
   Initialises a block, allocating space for MPIR_BLOCK integers 
   of size 0 <= n. Memory is actually allocated if n != 0
*/

void fmpz_block_init2(table_entry * entry, ulong n)
{
   entry->n = n;
   entry->count = MPIR_BLOCK;
  /* 
     Block has m integers, n+1 limbs per integer, MPIR_BYTES bytes per limb
  */
   if (n) entry->block_ptr = (block *) mpir_alloc(((n+1)<<MPIR_LG_BLOCK)<<MPIR_LG_BYTES);
}

/*
   Initialises a block with 1 <= m <= MPIR_BLOCK integers in it, but 
   doesn't actually allocate any memory fmpz_block_realloc and 
   fmpz_block_clear can be called on the block
*/

void fmpz_block_init_small(table_entry * entry, ulong m)
{
   entry->n = 0L;
   entry->count = m;
}

/*
   Initialises a block, allocating space for 1 <= m <= MPIR_BLOCK integers 
   of size 0 <= n. Memory is actually allocated if n != 0
*/

void fmpz_block_init2_small(table_entry * entry, ulong m, ulong n)
{
   entry->n = n;
   entry->count = m;
   /* 
      Block has m integers, n+1 limbs per integer, MPIR_BYTES bytes per limb
   */
   if (n) entry->block_ptr = (block *) mpir_alloc(((n+1)*m)<<MPIR_LG_BYTES);
}

/*
   Reallocate the integers in a block to have n limbs plus a sign/size limb
*/

void fmpz_block_realloc(table_entry * entry, ulong n)
{
   long oldn = entry->n;
   
   if (n == 0L)
   {
      fmpz_block_clear(entry);
      fmpz_block_init_small(entry, entry->count);
      return;
   }
   
   ulong m = entry->count;

   if (oldn == 0L)
   {
      entry->block_ptr = (block *) mpir_alloc(((n+1)*m)<<MPIR_LG_BYTES);
      // Now zero all the integers in the block
      mp_limb_t * ptr = (mp_limb_t*) entry->block_ptr;
      for (ulong i = 0; i < m; i++)
      {
         ptr[i*(n+1)] = 0L;
      }
   } else if (oldn < n)
   {
      entry->block_ptr = (block *) mpir_realloc(entry->block_ptr, ((n+1)*m)<<MPIR_LG_BYTES);
      mp_limb_t * ptr = entry->block_ptr + (m - 1)*(oldn+1);
      mp_limb_t * ptr2 = entry->block_ptr + (m - 1)*(n+1);
      for (int i = m - 1; i > 0; i--, ptr -= (oldn+1), ptr2 -= (n+1))
      {
         F_mpn_copy(ptr2, ptr, oldn+1);
      }
   } else if (oldn > n)
   {
      block * new_block = (block *) mpir_alloc(((n+1)*m)<<MPIR_LG_BYTES);
      mp_limb_t * ptr = entry->block_ptr + (m - 1)*(oldn+1);
      mp_limb_t * ptr2 = new_block + (m - 1)*(n+1);
      for (int i = m - 1; i >= 0; i--, ptr -= (oldn+1), ptr2 -= (n+1))
      {
         F_mpn_copy(ptr2, ptr, n+1);
      }
      free(entry->block_ptr);
      entry->block_ptr = new_block;
   }

   entry->n = n;
}

/*
   Clear a block (whether small or not) of integers
*/

void fmpz_block_clear(table_entry * entry)
{
   if (entry->n) free(entry->block_ptr);
}

/* ==============================================================================

   fmpz_t memory management

===============================================================================*/

/*
   Initialise an array of count != 0 integers and return a pointer to the first
*/

fmpz_t * fmpz_init_array(ulong count)
{
   ulong blocks = ((count-1)>>MPIR_LG_BLOCK) + 1;
   ulong blocks2 = (count>>MPIR_LG_BLOCK);
   table_entry* tab = (table_entry*) mpir_aligned_alloc(blocks<<MPIR_LG_BLOCK);
   ulong i;
   for (i = 0; i < blocks2; i++)
   {
      fmpz_block_init(tab + i);  
   }
   if (blocks != blocks2) fmpz_block_init_small(tab + i, count - (blocks2<<MPIR_LG_BLOCK));
   return (fmpz_t*) tab;
}

/*
   Initialise a single integer
*/

fmpz_t * fmpz_init(void)
{
   table_entry * entry = (table_entry*) mpir_aligned_alloc(MPIR_BLOCK);
   fmpz_block_init_small(entry, 1L);
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
   Realloc an array of fmpz_t's from length old_count to count
*/

fmpz_t * fmpz_realloc_array(fmpz_t * arr, ulong old_count, ulong count)
{
   ulong blocks = ((count-1)>>MPIR_LG_BLOCK) + 1;
   ulong blocks2 = (count>>MPIR_LG_BLOCK);
   ulong old_blocks = ((old_count-1)>>MPIR_LG_BLOCK) + 1;
   ulong old_blocks2 = (old_count>>MPIR_LG_BLOCK);
   if (count > old_count)
   {
      if (blocks > old_blocks)
      {
         table_entry * tab = (table_entry *) mpir_aligned_realloc(arr, blocks<<MPIR_LG_BLOCK); 
         if (old_blocks != old_blocks2)
         {
            table_entry * entry = (table_entry *) tab + old_blocks2;
            entry->count = MPIR_BLOCK;
            if (entry->n) entry->block_ptr = (block *) mpir_realloc(entry->block_ptr, ((entry->n + 1)<<MPIR_LG_BLOCK)<<MPIR_LG_BYTES);
         }
         ulong i;
         for (i = old_blocks; i < blocks2; i++)
         {
            fmpz_block_init(tab + i);  
         }
         if (blocks != blocks2) fmpz_block_init_small(tab + i, count - (blocks2<<MPIR_LG_BLOCK));
         return (fmpz_t *) tab;
      } else if (blocks == old_blocks)
      {
         table_entry * entry = (table_entry *) arr + old_blocks2;
         entry->count = count - (old_blocks2<<MPIR_LG_BLOCK);
         if (entry->n) entry->block_ptr = (block *) mpir_realloc(entry->block_ptr, (entry->n + 1)*(entry->count<<MPIR_LG_BYTES));
         return (fmpz_t *) arr;
      }
   } else if (count == old_count)
   {
   } else
   {
      printf("Error: attempted to shrink array!\n");
      abort();
   }

   return arr;
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

   String functions

===============================================================================*/

void fmpz_print(fmpz_t * in)
{
   mpz_t x;
   mpz_init(x);
   fmpz_to_mpz(x, in);
   gmp_printf("%Zd", x); 
   mpz_clear(x);
}

void fmpz_fread(fmpz_t * in, FILE * f)
{
   mpz_t x;
   mpz_init(x);
   mpz_inp_str(x, f, 10);
   mpz_to_fmpz(in, x);
   mpz_clear(x);
}

/* ==============================================================================

   Random generation

===============================================================================*/

/*
  This is not a serious random generator, it is just here for testing 
  purposes at this stage
  Bits must be non-zero
*/

void fmpz_random(fmpz_t * f, ulong bits)
{
   ulong limbs = ((bits-1)>>MPIR_LG_BITS)+1;
   ulong rem = (bits & (MPIR_BITS - 1));
   
   fmpz_fit_limbs(f, limbs);
   mp_limb_t * fp = fmpz_data(f);
   fp[0] = limbs;
   mpn_random(fp + 1, limbs);
   if (rem)
   {
      ulong mask = ((1L<<rem)-1L);
      fp[limbs] &= mask;
   }
   NORM(fp);
}

/* ==============================================================================

   Set/get

===============================================================================*/

void fmpz_set(fmpz_t * out, fmpz_t * f)
{
   if (out != f)
   {
      unsigned long d1n;
      mp_limb_t * dp = fmpz_data(f);
      ulong dn = MPIR_ABS(dp[0]);
      int re = fmpz_fit_limbs(out, MPIR_MAX(dn, 1L));
      mp_limb_t * rp = fmpz_data(out);
      if (re) dp = fmpz_data(f);
      F_mpn_copy(rp, dp, dn + 1); 
   }
}

/* ==============================================================================

   Negation

===============================================================================*/

void fmpz_neg(fmpz_t * out, fmpz_t * f)
{
   if (out != f)
   {
      unsigned long d1n;
      mp_limb_t * dp = fmpz_data(f);
      ulong dn = MPIR_ABS(dp[0]);
      int re = fmpz_fit_limbs(out, MPIR_MAX(dn, 1L));
      mp_limb_t * rp = fmpz_data(out);
      if (re) dp = fmpz_data(f);
      F_mpn_copy(rp + 1, dp + 1, dn); 
      rp[0] = -dp[0];
   } else
   { 
      mp_limb_t * rp = fmpz_data(out);
      rp[0] = -rp[0];
   }
}

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

void fmpz_add(fmpz_t * out, fmpz_t * f1, fmpz_t * f2)
{
   long cy;
   mp_limb_t * d1p = fmpz_data(f1);
   mp_limb_t * d2p = fmpz_data(f2);
   unsigned long d1n = MPIR_ABS(d1p[0]);
   unsigned long d2n = MPIR_ABS(d2p[0]);
   
   if (d1n < d2n) 
   {
      SWAP_FMPZ_PTRS(f1, f2);
      SWAP_SIZES(d1n, d2n);
      SWAP_LIMB_PTRS(d1p, d2p);
   } 

   int d11 = 0;
   int re;
   if ((OFFSET(out) == 0) & (d1n))
   {
      if (MPIR_BIT_COUNT(d1p[d1n]) == MPIR_BITS) 
      {
         re = fmpz_fit_limbs(out, d1n+1);
         d11 = 1;
      } else re = fmpz_fit_limbs(out, d1n);
   } else re = fmpz_fit_limbs(out, MPIR_MAX(d1n, 1L));

   mp_limb_t * rp = fmpz_data(out);
   if (re)
   {
      d1p = fmpz_data(f1);
      d2p = fmpz_data(f2); 
   }

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
      if (d1n == 1)
      {
         add_ssaaaa(cy, rp[1], 0L, d1p[1], 0L, d2p[1]);
      } else 
         cy = mpn_add(rp+1, d1p+1, d1n, d2p+1, d2n);
      if (cy) 
      {
         if (!d11)
         {
            re = fmpz_fit_limbs(out, d1n + 1);
            if (re) rp = fmpz_data(out);
         }
         rp[d1n+1] = cy;
         if ((long) rp[0] < 0L) rp[0]--;
         else rp[0]++;
      }
   } else
   {
      cy = 0;
      if (d1n == 1)
      {
         mp_limb_t diff = d1p[1] - d2p[1];
         if (d1p[1] > d2p[1]) 
         {
            rp[1] = d1p[1] - d2p[1];
            rp[0] = d1p[0];
         } else if (d1p[1] < d2p[1])
         {
            rp[1] = d2p[1] - d1p[1];
            rp[0] = -d1p[0];
         } else rp[0] = 0L;
      } else
      {
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
}

void fmpz_sub(fmpz_t * out, fmpz_t * f1, fmpz_t * f2)
{
   long cy;
   mp_limb_t * d1p = fmpz_data(f1);
   mp_limb_t * d2p = fmpz_data(f2);
   unsigned long d1n = MPIR_ABS(d1p[0]);
   unsigned long d2n = MPIR_ABS(d2p[0]);
   
   int in_order = 1;

   if (d1n < d2n) 
   {
      SWAP_FMPZ_PTRS(f1, f2);
      SWAP_SIZES(d1n, d2n);
      SWAP_LIMB_PTRS(d1p, d2p);
      in_order = 0;
   } 

   int d11 = 0;
   int re;
   if ((OFFSET(out) == 0) & (d1n))
   {
      if ((MPIR_BIT_COUNT(d1p[d1n]) == MPIR_BITS) 
         || ((d2n == d1n) && (MPIR_BIT_COUNT(d2p[d2n]) == MPIR_BITS)))
         {
            re = fmpz_fit_limbs(out, d1n+1);
            d11 = 1;
         } else re = fmpz_fit_limbs(out, d1n);
   } else re = fmpz_fit_limbs(out, MPIR_MAX(d1n, 1L));

   mp_limb_t * rp = fmpz_data(out);
   if (re)
   {
      d1p = fmpz_data(f1);
      d2p = fmpz_data(f2); 
   }
   
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
      if (d2n == 1)
      {
         cy = mpn_add_1(rp+1, d1p+1, d1n, d2p[1]);
      } else cy = mpn_add(rp+1, d1p+1, d1n, d2p+1, d2n);
      if (cy) 
      {
         if (!d11)
         {
            re = fmpz_fit_limbs(out, d1n+1);
            if (re) rp = fmpz_data(out);
         }
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
         if (d2n == 1) mpn_sub_1(rp+1, d1p+1, d1n, d2p[1]);
         else mpn_sub(rp+1, d1p+1, d1n, d2p+1, d2n);
         if (in_order) rp[0] = d1p[0];
         else rp[0] = -d1p[0];
         NORM(rp);
      }
      else
      {
         if (d1n == 1) rp[1] = d2p[1] - d1p[1];
         else mpn_sub_n(rp+1, d2p+1, d1p+1, d1n);
         if (in_order) rp[0] = -d1p[0];
         else rp[0] = d1p[0];
         NORM(rp);
      }
   }
}

/* ==============================================================================

   Addmul/submul

===============================================================================*/

/*
   w = w + x*y
*/

void fmpz_addmul_ui(fmpz_t * w, fmpz_t * x, ulong y)
{
   ulong  xn, wn, wss, new_wn, min_n, dn;
   mp_limb_t * xp = fmpz_data(x);
   mp_limb_t * wp = fmpz_data(w);
   mp_limb_t cy;
   int re;

   xn = xp[0];
   if (xn == 0 || y == 0)
     return;

   mp_limb_t sub = xn;
   xn = MPIR_ABS(xn);

   wss = wp[0];
   if (wss == 0L)
   {
      re = fmpz_fit_limbs(w, xn+1);
      if (re) 
      {
         wp = fmpz_data(w);
         xp = fmpz_data(x);
      }
      cy = mpn_mul_1(wp+1, xp+1, xn, y);
      wp[xn+1] = cy;
      xn += (cy != 0L);
      wp[0] = ((long) sub >= 0L ? xn : -xn);
      return;
   }

   sub ^= wss;
   wn = MPIR_ABS(wss);
   new_wn = MPIR_MAX(wn, xn);
   re = fmpz_fit_limbs(w, new_wn+1);
   if (re)
   {
      wp = fmpz_data(w);
      xp = fmpz_data(x);
   }
   min_n = MPIR_MIN(wn, xn);

   if ((long) sub >= 0L)
   {
      /* addmul of absolute values */

      cy = mpn_addmul_1 (wp+1, xp+1, min_n, y);
      wp += min_n;
      xp += min_n;

      dn = xn - wn;
      if (dn)
      {
         mp_limb_t cy2;
         if ((long) dn > 0L)
            cy2 = mpn_mul_1(wp+1, xp+1, dn, y);
         else
         {
            dn = -dn;
            cy2 = 0L;
         }
         cy = cy2 + mpn_add_1(wp+1, wp+1, dn, cy);
      }
      wp[dn+1] = cy;
      new_wn += (cy != 0L);
      wp -= min_n;
   } else
   {
      /* submul of absolute values */

      cy = mpn_submul_1(wp+1, xp+1, min_n, y);
      if (wn >= xn)
      {
         /* if w bigger than x, then propagate borrow through it */
         if (wn != xn)
            cy = mpn_sub_1(wp+xn+1, wp+xn+1, wn-xn, cy);

         if (cy != 0)
         {
            /* Borrow out of w, take twos complement negative to get
               absolute value, flip sign of w.  */
            wp[new_wn+1] = ~-cy;  /* extra limb is 0-cy */
            F_mpn_com(wp+1, wp+1, new_wn);
            new_wn++;
            mpn_add_1(wp+1, wp+1, new_wn, 1L); // No carry
            wss = -wss;
         }
      }
      else /* wsize < xsize */
      {
         /* x bigger than w, so want x*y-w.  Submul has given w-x*y, so
            take twos complement and use an mpn_mul_1 for the rest.  */

         mp_limb_t cy2;

         /* -(-cy*b^n + w-x*y) = (cy-1)*b^n + ~(w-x*y) + 1 */
         F_mpn_com(wp+1, wp+1, wn);
         cy += mpn_add_1(wp+1, wp+1, wn, 1L);
         cy -= 1;

         /* If cy-1 == -1 then hold that -1 for latter.  mpn_submul_1 never
            returns cy==UINT_MAX so that value always indicates a -1. */
         cy2 = (cy == -1L);
         cy += cy2;
         mp_limb_t cy3;
         cy3 = mpn_mul_1(wp+wn+1, xp+wn+1, xn-wn, y);               
         cy = cy3 + mpn_add_1(wp+wn+1, wp+wn+1, xn-wn, cy);    
         wp[new_wn+1] = cy;
         new_wn += (cy != 0);

          /* Apply any -1 from above.  The value at wp+wsize is non-zero
             because y!=0 and the high limb of x will be non-zero.  */
         if (cy2)
            mpn_sub_1(wp+wn+1, wp+wn+1, new_wn-wn, 1L);

          wss = -wss;
        }

        /* submul can produce high zero limbs due to cancellation, both when w
           has more limbs or x has more  */
        wp[0] = new_wn;
        NORM(wp);
    }

    wp[0] = ((long) wss >= 0L ? new_wn : -new_wn);
    NORM(wp);
}

/*
   w = w - x*y
*/

void fmpz_submul_ui(fmpz_t * w, fmpz_t * x, ulong y)
{
   ulong  xn, wn, wss, new_wn, min_n, dn;
   mp_limb_t * xp = fmpz_data(x);
   mp_limb_t * wp = fmpz_data(w);
   mp_limb_t cy;
   int re;

   xn = xp[0];
   if (xn == 0 || y == 0)
     return;

   mp_limb_t sub = ~xn;
   xn = MPIR_ABS(xn);

   wss = wp[0];
   if (wss == 0L)
   {
      re = fmpz_fit_limbs(w, xn+1);
      if (re) 
      {
         wp = fmpz_data(w);
         xp = fmpz_data(x);
      }
      cy = mpn_mul_1(wp+1, xp+1, xn, y);
      wp[xn+1] = cy;
      xn += (cy != 0L);
      wp[0] = ((long) sub >= 0L ? xn : -xn);
      return;
   }

   sub ^= wss;
   wn = MPIR_ABS(wss);
   new_wn = MPIR_MAX(wn, xn);
   re = fmpz_fit_limbs(w, new_wn+1);
   if (re)
   {
      wp = fmpz_data(w);
      xp = fmpz_data(x);
   }
   min_n = MPIR_MIN(wn, xn);

   if ((long) sub >= 0L)
   {
      /* addmul of absolute values */

      cy = mpn_addmul_1 (wp+1, xp+1, min_n, y);
      wp += min_n;
      xp += min_n;

      dn = xn - wn;
      if (dn)
      {
         mp_limb_t cy2;
         if ((long) dn > 0L)
            cy2 = mpn_mul_1(wp+1, xp+1, dn, y);
         else
         {
            dn = -dn;
            cy2 = 0L;
         }
         cy = cy2 + mpn_add_1(wp+1, wp+1, dn, cy);
      }
      wp[dn+1] = cy;
      new_wn += (cy != 0L);
      wp -= min_n;
   } else
   {
      /* submul of absolute values */

      cy = mpn_submul_1(wp+1, xp+1, min_n, y);
      if (wn >= xn)
      {
         /* if w bigger than x, then propagate borrow through it */
         if (wn != xn)
            cy = mpn_sub_1(wp+xn+1, wp+xn+1, wn-xn, cy);

         if (cy != 0)
         {
            /* Borrow out of w, take twos complement negative to get
               absolute value, flip sign of w.  */
            wp[new_wn+1] = ~-cy;  /* extra limb is 0-cy */
            F_mpn_com(wp+1, wp+1, new_wn);
            new_wn++;
            mpn_add_1(wp+1, wp+1, new_wn, 1L); // No carry
            wss = -wss;
         }
      }
      else /* wsize < xsize */
      {
         /* x bigger than w, so want x*y-w.  Submul has given w-x*y, so
            take twos complement and use an mpn_mul_1 for the rest.  */

         mp_limb_t cy2;

         /* -(-cy*b^n + w-x*y) = (cy-1)*b^n + ~(w-x*y) + 1 */
         F_mpn_com(wp+1, wp+1, wn);
         cy += mpn_add_1(wp+1, wp+1, wn, 1L);
         cy -= 1;

         /* If cy-1 == -1 then hold that -1 for latter.  mpn_submul_1 never
            returns cy==UINT_MAX so that value always indicates a -1. */
         cy2 = (cy == -1L);
         cy += cy2;
         mp_limb_t cy3;
         cy3 = mpn_mul_1(wp+wn+1, xp+wn+1, xn-wn, y);               
         cy = cy3 + mpn_add_1(wp+wn+1, wp+wn+1, xn-wn, cy);    
         wp[new_wn+1] = cy;
         new_wn += (cy != 0);

          /* Apply any -1 from above.  The value at wp+wsize is non-zero
             because y!=0 and the high limb of x will be non-zero.  */
         if (cy2)
            mpn_sub_1(wp+wn+1, wp+wn+1, new_wn-wn, 1L);

          wss = -wss;
        }

        /* submul can produce high zero limbs due to cancellation, both when w
           has more limbs or x has more  */
        wp[0] = new_wn;
        NORM(wp);
    }

    wp[0] = ((long) wss >= 0L ? new_wn : -new_wn);
    NORM(wp);
}

/* ==============================================================================

   Multiplication

===============================================================================*/

/*
   w = u*2^exp for 0 <= exp
*/

void fmpz_mul_2exp(fmpz_t * w, fmpz_t * u, ulong exp)
{
   mp_limb_t * up = fmpz_data(u);
   long uss = up[0];
   ulong un = MPIR_ABS(uss);
   ulong wn;
   ulong limb_cnt;
   mp_limb_t * wp;
   mp_limb_t wlimb;

   if (uss == 0L)
   {
      fmpz_fit_limbs(w, 1L);
      fmpz_data(w)[0] = 0L;
      return;
   }

   limb_cnt = (exp>>MPIR_LG_BITS);
   exp &= (MPIR_BITS-1);
  
   int re = fmpz_fit_limbs(w, un + limb_cnt + (exp != 0L));
   if (re) up = fmpz_data(u);
   
   wp = fmpz_data(w);
   wn = un + limb_cnt;

   if (exp)
   {
      wlimb = mpn_lshift(wp + limb_cnt + 1, up + 1, un, exp);
      if (wlimb)
	  {
	     wp[wn + 1] = wlimb;
	     wn++;
      }
   } else
   {
      F_mpn_copy(wp + limb_cnt + 1, up + 1, un);
   }

   /* Zero all whole limbs at low end.  Do it here and not before calling
      mpn_lshift, so as not to lose data when U == W.  */
   F_mpn_clear(wp + 1, limb_cnt);

   wp[0] = (uss >= 0L ? wn : -wn);
}


// *************** end of file
