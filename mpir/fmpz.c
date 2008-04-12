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
#include <math.h>

#include "fmpz.h"
#include "mpir.h"
#include "memory_manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "test_support.h"

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

#define NORM(xxx, coeffxxx) \
do { \
   if ((long) xxx->_mp_size < 0L) \
   { \
      while ((xxx->_mp_size) && (!(coeffxxx)[-xxx->_mp_size - 1])) xxx->_mp_size++; \
   } else if ((long) xxx->_mp_size > 0L) \
   { \
      while ((xxx->_mp_size) && (!(coeffxxx)[xxx->_mp_size - 1])) xxx->_mp_size--; \
   } \
} while (0);

#define MPN_COPY_INCR(dst, src, n)                      \
  do {                                                  \
    ASSERT ((n) >= 0);                                  \
    ASSERT (MPN_SAME_OR_INCR_P (dst, src, n));          \
    if ((n) != 0)                                       \
      {                                                 \
	mp_size_t __n = (n) - 1;                        \
	mp_ptr __dst = (dst);                           \
	mp_srcptr __src = (src);                        \
	mp_limb_t __x;                                  \
	__x = *__src++;                                 \
	if (__n != 0)                                   \
	  {                                             \
	    do                                          \
	      {                                         \
		*__dst++ = __x;                         \
		__x = *__src++;                         \
	      }                                         \
	    while (--__n);                              \
	  }                                             \
	*__dst++ = __x;                                 \
      }                                                 \
  } while (0)

/* ==============================================================================

   Block memory management

===============================================================================*/

/*
   Initialises a block with 0 < m < MPIR_BLOCK  integers in it, but 
   doesn't actually allocate any memory 
   fmpz_block_realloc and fmpz_block_clear can be called on the block
*/

void fmpz_block_init(fmpz_t * entry, ulong m)
{
   entry->_mp_size = m;
   entry->_mp_alloc = 0L;
   entry->_mp_d = (mp_limb_t *) 0L;
   entry++;
   for (ulong i = 1; i < m; i++, entry++)
   {
      entry->_mp_alloc = 0L;
      entry->_mp_d = (mp_limb_t *) 0L;
   }
}

/*
   Initialises a block, allocating space for MPIR_BLOCK integers 
   of size 0 < n. 
*/

void fmpz_block_init2(fmpz_t * entry, ulong m, ulong n)
{
   if (n == 0L)
   {
      fmpz_block_init(entry, m);
      return;
   }

   ulong intsize = (n<<MPIR_LG_BYTES);
   mp_limb_t * block_ptr = (mp_limb_t *) mpir_alloc(intsize*m+sizeof(mp_limb_t));
   block_ptr[0] = m;
   block_ptr++;
   for (ulong i = 0; i < m; i++, entry++, block_ptr += n)
   {
      entry->_mp_alloc = n;
      entry->_mp_size = 0L;
      entry->_mp_d = block_ptr;
   }
}

/*
   Reallocate the integers in a block to have n limbs plus a sign/size limb
*/

void fmpz_block_realloc(fmpz_t * entry, ulong n)
{
   long oldn = entry->_mp_alloc;
   
   if (oldn < n)
   {
      if (oldn == 0L)
      {
         ulong intsize = (n<<MPIR_LG_BYTES);
         ulong m = entry->_mp_size;
         mp_limb_t * new_block = (mp_limb_t *) mpir_alloc(intsize*m+sizeof(mp_limb_t));
         new_block[0] = m;
         new_block++;
         mp_limb_t * ptr2 = new_block;
         long int_i;
         for (ulong i = 0; i < m; i++, ptr2 += n, entry ++)
         {
            int_i = ((long) entry->_mp_d);
            if (int_i > 0L) 
            {
               entry->_mp_size = 1L;
               ptr2[0] = int_i;
            } else if (int_i < 0L) 
            {
               entry->_mp_size = -1L;
               ptr2[0] = -int_i;
            } else entry->_mp_size = 0L;
            entry->_mp_alloc = n;
            entry->_mp_d = ptr2;
         }
         return;
      }
      ulong intsize = (n<<MPIR_LG_BYTES);
      mp_limb_t * block_ptr = entry->_mp_d - 1;
      ulong m = block_ptr[0];
      block_ptr++;
      mp_limb_t * new_block = (mp_limb_t *) mpir_alloc(intsize*m+sizeof(mp_limb_t));
      new_block[0] = m;
      new_block++;
      mp_limb_t * ptr = block_ptr + (m - 1)*oldn;
      mp_limb_t * ptr2 = new_block + (m - 1)*n;
      entry += (m-1);
      for (int i = m - 1; i >= 0; i--, ptr -= oldn, ptr2 -= n, entry--)
      {
         entry->_mp_d = ptr2;
         entry->_mp_alloc = n;
         F_mpn_copy(ptr2, ptr, oldn);
      }
      free(block_ptr - 1);
   } else if (oldn > n)
   {
      mp_limb_t * block_ptr = entry->_mp_d - 1;
      ulong m = block_ptr[0];
      if (n == 0L)
      {
         for (ulong i = 0; i < m; i++, entry++)
         {
             if (entry->_mp_size > 0L) entry->_mp_d = (mp_limb_t *) entry->_mp_d[0];
             else if (entry->_mp_size < 0L) entry->_mp_d = (mp_limb_t *) -entry->_mp_d[0];
             else entry->_mp_d = (mp_limb_t *) 0L;
             entry->_mp_alloc = 0L;
         }
         free(block_ptr);
         return;
      }
      block_ptr++;
      ulong intsize = (n<<MPIR_LG_BYTES);
      mp_limb_t * new_block = (mp_limb_t *) mpir_alloc(intsize*m+sizeof(mp_limb_t));
      new_block[0] = m;
      new_block++;
      mp_limb_t * ptr = block_ptr + (m - 1)*oldn;
      mp_limb_t * ptr2 = new_block + (m - 1)*n;
      entry += (m-1);
      for (int i = m - 1; i >= 0; i--, ptr -= oldn, ptr2 -= n, entry--)
      {
         entry->_mp_d = ptr2;
         entry->_mp_alloc = n;
         F_mpn_copy(ptr2, ptr, n);
      }
      free(block_ptr - 1);
   }
}

/*
   Clear a block (whether small or not) of integers
*/

void fmpz_block_clear(fmpz_t * entry)
{
   if (entry->_mp_alloc > 0L) 
   {
      mp_limb_t * data = entry->_mp_d - 1;
      free(data);
   }
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
   fmpz_t * tab = (fmpz_t *) mpir_aligned_alloc(count*sizeof(fmpz_t));
   ulong i;
   fmpz_t * ptr = tab;
   for (i = 0; i < blocks2; i++, ptr+= MPIR_BLOCK)
   {
      fmpz_block_init(ptr, MPIR_BLOCK);  
   }
   if (blocks != blocks2) fmpz_block_init(ptr, count - (blocks2<<MPIR_LG_BLOCK));
   return (fmpz_t *) tab;
}

/*
   Clear the array of count integers where f is a pointer to the first
   Partial arrays cannot be cleared. The full array must be cleared
   with count set to the current number of integers in the array
*/

void fmpz_clear_array(fmpz_t * f, ulong count)
{
   ulong blocks = ((count-1)>>MPIR_LG_BLOCK) + 1;
   fmpz_t * ptr = (fmpz_t *) f;
   for (ulong i = 0; i < blocks; i++, ptr += MPIR_BLOCK)
   {
      fmpz_block_clear(ptr);  
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
      fmpz_t * tab = (fmpz_t *) mpir_aligned_realloc(arr, count*sizeof(fmpz_t)); 
      fmpz_t * entry = (fmpz_t *) tab + (old_blocks2<<MPIR_LG_BLOCK);
      if (blocks > old_blocks)
      {
         if (old_blocks != old_blocks2)
         {
            ulong n = entry->_mp_alloc;
            if (n == 0L)
            {
               ulong m = entry->_mp_size;
               entry->_mp_size = MPIR_BLOCK;
               entry += m;
               for (ulong i = m; i < MPIR_BLOCK; i++, entry++)
               {
                  entry->_mp_alloc = 0L;
                  entry->_mp_d = (mp_limb_t *) 0L;
               }              
            } else  
            {
               mp_limb_t * block_ptr = entry->_mp_d - 1;
               ulong m = block_ptr[0];
               ulong intsize = (n<<MPIR_LG_BYTES);
               block_ptr = (mp_limb_t *) mpir_realloc(block_ptr, (intsize<<MPIR_LG_BLOCK)+sizeof(mp_limb_t));
               block_ptr[0] = MPIR_BLOCK;
               block_ptr++;
               ulong i;
               for (i = 0; i < m; i++, entry++, block_ptr += n)
               {
                  entry->_mp_d = block_ptr;
               }
               for (; i < MPIR_BLOCK; i++, entry++, block_ptr += n)
               {
                  entry->_mp_alloc = n;
                  entry->_mp_size = 0L;
                  entry->_mp_d = block_ptr;
               }
            } 
         }
         ulong i;
         for (i = old_blocks; i < blocks2; i++, entry += MPIR_BLOCK)
         {
            fmpz_block_init(entry, MPIR_BLOCK);  
         }
         if (blocks != blocks2) fmpz_block_init(entry, count - (blocks2<<MPIR_LG_BLOCK));
         return (fmpz_t *) tab;
      } else if (blocks == old_blocks)
      {
         ulong n = entry->_mp_alloc;
         ulong newm = count - (old_blocks2<<MPIR_LG_BLOCK);
         if (n == 0L)
         {
            ulong m = entry->_mp_size;
            entry->_mp_size = newm;
            entry += m;
            for (ulong i = m; i < newm; i++, entry++)
            {
               entry->_mp_alloc = 0L;
               entry->_mp_d = (mp_limb_t *) 0L;
            }              
         } else 
         {
            mp_limb_t * block_ptr = entry->_mp_d - 1;
            ulong m = block_ptr[0];
            ulong intsize = (n<<MPIR_LG_BYTES);
            block_ptr = (mp_limb_t *) mpir_realloc(block_ptr, intsize*newm+sizeof(mp_limb_t));
            block_ptr[0] = newm;
            block_ptr++;
            ulong i;
            for (i = 0; i < m; i++, entry++, block_ptr += n)
            {
               entry->_mp_d = block_ptr;
            }
            for (; i < newm; i++, entry++, block_ptr += n)
            {
               entry->_mp_alloc = n;
               entry->_mp_size = 0L;
               entry->_mp_d = block_ptr;
            }
         }
         return (fmpz_t *) tab;
      }
   } else if (count == old_count)
   {
      return arr;
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

void mpz_to_fmpz(fmpz_t * fnum, mpz_t num)
{
   ulong size = mpz_size(num);
   
   if ((size > 1L) || (fnum->_mp_alloc > 0L))
   {
      fmpz_fit_limbs(fnum, size);
      mpz_set(fnum, num);
   } else
   {
      ulong m_int = mpz_get_ui(num);
      if (m_int > IMM_MAX) 
      {
         fmpz_fit_limbs(fnum, 1L);
         mpz_set(fnum, num);   
      } else
      {
         fnum->_mp_d = (mp_limb_t *) mpz_get_si(num);
      }
   }
}

void fmpz_to_mpz(mpz_t num, fmpz_t * fnum)
{
   if (fnum->_mp_alloc > 0L) mpz_set(num, fnum);
   else
   {
      mpz_set_si(num, (long) fnum->_mp_d);
   }
}

double fmpz_get_d(fmpz_t * f)
{
   if (f->_mp_alloc == 0L)
   {
      long int_f = (long) f->_mp_d;
      if (int_f == 0) return 0.0;
      ulong int_abs = MPIR_ABS(int_f);
      if (int_f < 0L) return __gmpn_get_d(&int_abs, 1L, -1L, 0L);
      else return __gmpn_get_d(&int_abs, 1L, 1L, 0L);
   } else return mpz_get_d(f);
}

double fmpz_get_d_2exp(long * exp, fmpz_t * f)
{
   if (f->_mp_alloc == 0L)
   {
      double d;
      long int_f = (long) f->_mp_d;
      if (int_f == 0L) 
      {
         (*exp) = 0L;
         return 0.0;
      }
      ulong int_abs = MPIR_ABS(int_f);
      (*exp) = MPIR_BIT_COUNT(int_abs);
      if (int_f < 0L) return __gmpn_get_d(&int_abs, 1L, -1L, -*exp);
      else return __gmpn_get_d(&int_f, 1L, 1L, -*exp);
   } else return mpz_get_d_2exp(exp, f);
}

/* ==============================================================================

   String functions

===============================================================================*/

/*
   These are not meant to be efficient at this point
*/

void fmpz_print(fmpz_t * in)
{
   if (in->_mp_alloc == 0L) printf("%ld", (long) in->_mp_d);
   else gmp_printf("%Zd", in);
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
  We require bits to be non-zero
*/

void fmpz_random(fmpz_t * f, ulong bits)
{
   if ((f->_mp_alloc == 0L) && (bits <= MPIR_BITS - 2))
   {
      ulong temp;
      mpn_random(&temp, 1L);
      ulong mask = ((1L<<bits)-1L);
      f->_mp_d = (mp_limb_t *) (temp & mask);
      return;
   }
   ulong limbs = ((bits-1)>>MPIR_LG_BITS)+1;
   ulong rem = (bits & (MPIR_BITS - 1));
   
   fmpz_fit_limbs(f, limbs);
   mp_limb_t * fp = f->_mp_d;
   f->_mp_size = limbs;
   mpn_random(fp, limbs);
   if (rem)
   {
      ulong mask = ((1L<<rem)-1L);
      fp[limbs-1] &= mask;
   }
   NORM(f, fp);
}

void fmpz_randomm(fmpz_t * out, fmpz_t * in)
{
   if ((in->_mp_alloc > 0L) && (out->_mp_alloc > 0L)) 
   {
      fmpz_fit_limbs(out, MPIR_ABS(in->_mp_size));
      mpz_urandomm(out, state, in);
   } else
   {
      mpz_t temp1, temp2;
      mpz_init(temp1);
      mpz_init(temp2);
      fmpz_to_mpz(temp1, in);
      mpz_urandomm(temp2, state, temp1);
      mpz_to_fmpz(out, temp2);
      mpz_clear(temp1);
      mpz_clear(temp2);
   }
}

/* ==============================================================================

   Number theoretical

===============================================================================*/

int fmpz_probab_prime_p(fmpz_t * p, ulong n)
{
   if (p->_mp_alloc > 0L) return mpz_probab_prime_p(p, n);
   else
   {
      mpz_t temp1;
      mpz_init(temp1);
      fmpz_to_mpz(temp1, p);
      int prime = mpz_probab_prime_p(temp1, n);
      mpz_clear(temp1); 
      return prime;     
   }
}

/* ==============================================================================

   Set/get

===============================================================================*/

void fmpz_set(fmpz_t * w, fmpz_t * u)
{
  if (u != w) 
  {
     if (u->_mp_alloc == 0L)
     {
        fmpz_set_si(w, (long) u->_mp_d);
     } else
     {
        mp_limb_t * wp, * up;
        mp_size_t usize, size;

        usize = u->_mp_size;
        if (usize == 0L) 
        {
           fmpz_set_ui(w, 0L);
           return;
        }
        size = MPIR_ABS (usize);

        fmpz_fit_limbs(w, size); 

        wp = w->_mp_d;
        up = u->_mp_d;
        
        MPN_COPY_INCR (wp, up, size);
        w->_mp_size = usize;
      }
   }
}

/* ==============================================================================

   Negation

===============================================================================*/

void fmpz_neg(fmpz_t * out, fmpz_t * f)
{
   if (out != f) fmpz_set(out, f);
   if (out->_mp_alloc == 0L) out->_mp_d = (mp_limb_t *) -((long) out->_mp_d);
   else out->_mp_size = -out->_mp_size;
}

/* ==============================================================================

   Comparison

===============================================================================*/

int fmpz_equal(fmpz_t * f2, fmpz_t * f1)
{
   if UNLIKELY(f1 == f2) return 1;
   ulong a1 = f1->_mp_alloc;
   ulong a2 = f1->_mp_alloc;
   if ((a1 > 0L) && (a2 > 0L)) return (mpz_cmp(f1, f2) == 0);
   if (a1 == 0L)
   {
      if (a2 == 0L) return (f1->_mp_d == f2->_mp_d);
      return (mpz_cmp_si(f2, (long) f1->_mp_d) == 0);
   }
   return (mpz_cmp_si(f1, (long) f2->_mp_d) == 0);
}

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

/* 
   Adds a signed quantity no bigger than MPIR_BITS - 2 bits
*/

void _fmpz_add_IMM(fmpz_t * out, fmpz_t * f1, long c)
{
   if (f1->_mp_alloc == 0L)
   {
      fmpz_set_si(out, c + (long) f1->_mp_d); 
   } else
   {
      fmpz_fit_limbs(out, MPIR_ABS(f1->_mp_size) + 1);
      if (c >= 0L) mpz_add_ui(out, f1, c);
      else mpz_sub_ui(out, f1, -c);
   }
}

void fmpz_add(fmpz_t * out, fmpz_t * f1, fmpz_t * f2)
{
   if (f1->_mp_alloc == 0L)
   {
      _fmpz_add_IMM(out, f2, (long) f1->_mp_d);
      return;
   } else if (f2->_mp_alloc == 0L)
   {
      _fmpz_add_IMM(out, f1, (long) f2->_mp_d);
      return;
   } 
   
   long cy;
   unsigned long d1n = MPIR_ABS(f1->_mp_size);
   unsigned long d2n = MPIR_ABS(f2->_mp_size);
   
   if (d1n < d2n) 
   {
      SWAP_FMPZ_PTRS(f1, f2);
      SWAP_SIZES(d1n, d2n);
   } 

   int re;
   re = fmpz_fit_limbs(out, d1n + 1);

   mp_limb_t * rp = fmpz_data(out);
   mp_limb_t * d1p = fmpz_data(f1);
   mp_limb_t * d2p = fmpz_data(f2); 
   
   if UNLIKELY(!d2n)
   {
      if UNLIKELY(!d1n) out->_mp_size = 0L;
      else
      {
         if (rp != d1p) 
         {
            MPN_COPY_INCR(rp, d1p, d1n);
            out->_mp_size =  f1->_mp_size;
         }
      }
   } else if ((long) (f1->_mp_size ^ f2->_mp_size) >= 0L)
   {
      out->_mp_size = f1->_mp_size;
      cy = mpn_add(rp, d1p, d1n, d2p, d2n);
      if (cy) 
      {
         rp[d1n] = cy;
         if ((long) out->_mp_size < 0L) out->_mp_size--;
         else out->_mp_size++;
      }
   } else
   {
      cy = 0;
      if (d1n != d2n) cy = 1;
      else cy = mpn_cmp(d1p, d2p, d1n); 
          
      if (cy == 0) out->_mp_size = 0L;
      else if (cy > 0) 
      {
         mpn_sub(rp, d1p, d1n, d2p, d2n);
         out->_mp_size = f1->_mp_size;
         NORM(out, rp);
      }
      else
      {
         mpn_sub_n(rp, d2p, d1p, d1n);
         out->_mp_size = -f1->_mp_size;
         NORM(out, rp);
      }
   }
}

/* 
   Subtracts a signed quantity no bigger than MPIR_BITS - 2 bits
*/

void _fmpz_sub_IMM(fmpz_t * out, fmpz_t * f1, long c)
{
   if (f1->_mp_alloc == 0L)
   {
      fmpz_set_si(out, ((long) f1->_mp_d) - c); 
   } else
   {
      fmpz_fit_limbs(out, MPIR_ABS(f1->_mp_size) + 1);
      if (c >= 0L) mpz_sub_ui(out, f1, c);
      else mpz_add_ui(out, f1, -c);
   }
}

void fmpz_sub(fmpz_t * out, fmpz_t * f1, fmpz_t * f2)
{
   if (f1->_mp_alloc == 0L)
   {
      _fmpz_sub_IMM(out, f2, (long) f1->_mp_d);
      fmpz_neg(out, out);
      return;
   } else if (f2->_mp_alloc == 0L)
   {
      _fmpz_sub_IMM(out, f1, (long) f2->_mp_d);
      return;
   } 

   long cy;
   unsigned long d1n = MPIR_ABS(f1->_mp_size);
   unsigned long d2n = MPIR_ABS(f2->_mp_size);
   int in_order = 1;
   
   if (d1n < d2n) 
   {
      SWAP_FMPZ_PTRS(f1, f2);
      SWAP_SIZES(d1n, d2n);
      in_order = 0;
   } 

   int re;
   re = fmpz_fit_limbs(out, d1n + 1);

   mp_limb_t * rp = fmpz_data(out);
   mp_limb_t * d1p = fmpz_data(f1);
   mp_limb_t * d2p = fmpz_data(f2); 
   
   if UNLIKELY(!d2n)
   {
      if UNLIKELY(!d1n) out->_mp_size = 0L;
      else
      {
         if (rp != d1p) 
         {
            MPN_COPY_INCR(rp, d1p, d1n);
            out->_mp_size =  f1->_mp_size;
            if (!in_order) out->_mp_size = -out->_mp_size;
         }
      }
   } else if ((long) (f1->_mp_size ^ f2->_mp_size) < 0L)
   {
      if (in_order) out->_mp_size = f1->_mp_size;
      else out->_mp_size = -f1->_mp_size;
      if (d1n == 1)
      {
         add_ssaaaa(cy, rp[0], 0L, d1p[0], 0L, d2p[0]);
      } else 
         cy = mpn_add(rp, d1p, d1n, d2p, d2n);
      if (cy) 
      {
         rp[d1n] = cy;
         if ((long) out->_mp_size < 0L) out->_mp_size--;
         else out->_mp_size++;
      }
   } else
   {
      if (d1n == 1)
      {
         if (d1p[0] > d2p[0]) 
         {
            rp[0] = d1p[0] - d2p[0];
            if (in_order) out->_mp_size = f1->_mp_size;
            else out->_mp_size = -f1->_mp_size;
         } else if (d1p[0] < d2p[0])
         {
            rp[0] = d2p[0] - d1p[0];
            if (in_order) out->_mp_size = -f1->_mp_size;
            else out->_mp_size = f1->_mp_size;
         } else out->_mp_size = 0L;
      } else
      {
         cy = 0;
         if (d1n != d2n) cy = 1;
         else cy = mpn_cmp(d1p, d2p, d1n); 
          
         if (cy == 0) out->_mp_size = 0L;
         else if (cy > 0) 
         {
            mpn_sub(rp, d1p, d1n, d2p, d2n);
            if (in_order) out->_mp_size = f1->_mp_size;
            else out->_mp_size = -f1->_mp_size;
            NORM(out, rp);
         }
         else
         {
            mpn_sub_n(rp, d2p, d1p, d1n);
            if (in_order) out->_mp_size = -f1->_mp_size;
            else out->_mp_size = f1->_mp_size;
            NORM(out, rp);
         }
      }
   }
}

/* ==============================================================================

   Addmul/submul

===============================================================================*/
#define mpn_incr_u(p,incr)                              \
  do {                                                  \
    mp_limb_t __x;                                      \
    mp_ptr __p = (p);                                   \
    if (__builtin_constant_p (incr) && (incr) == 1)     \
      {                                                 \
        while (++(*(__p++)) == 0)                       \
          ;                                             \
      }                                                 \
    else                                                \
      {                                                 \
        __x = *__p + (incr);                            \
        *__p = __x;                                     \
        if (__x < (incr))                               \
          while (++(*(++__p)) == 0)                     \
            ;                                           \
      }                                                 \
  } while (0)

#define mpn_decr_u(p,incr)                              \
  do {                                                  \
    mp_limb_t __x;                                      \
    mp_ptr __p = (p);                                   \
    if (__builtin_constant_p (incr) && (incr) == 1)     \
      {                                                 \
        while ((*(__p++))-- == 0)                       \
          ;                                             \
      }                                                 \
    else                                                \
      {                                                 \
        __x = *__p;                                     \
        *__p = __x - (incr);                            \
        if (__x < (incr))                               \
          while ((*(++__p))-- == 0)                     \
            ;                                           \
      }                                                 \
  } while (0)

#define MPN_INCR_U(ptr, size, n)   mpn_incr_u (ptr, n)
#define MPN_DECR_U(ptr, size, n)   mpn_decr_u (ptr, n)

#define mpn_com_n(d,s,n)                                \
  do {                                                  \
    mp_ptr     __d = (d);                               \
    mp_srcptr  __s = (s);                               \
    mp_size_t  __n = (n);                               \
    ASSERT (__n >= 1);                                  \
    ASSERT (MPN_SAME_OR_SEPARATE_P (__d, __s, __n));    \
    do                                                  \
      *__d++ = (~ *__s++) & GMP_NUMB_MASK;              \
    while (--__n);                                      \
  } while (0)

/*
   w = w + x*y
*/
 
void _fmpz_addmul_ui(fmpz_t * w, fmpz_t * x, ulong y)
{
   ulong  xn, wn, wss, new_wn, min_n, dn;
   mp_limb_t * xp;
   mp_limb_t * wp;
   mp_limb_t cy;
   
   xn = x->_mp_size;
   if UNLIKELY(xn == 0 || y == 0L)
     return;

   mp_limb_t sub = xn;
   xn = MPIR_ABS(xn);

   if (w->_mp_alloc == 0L)
   {
      if ((long) w->_mp_d < 0L) wss = -1L;
      else if ((long) w->_mp_d > 0L) wss = 1L;
      else wss = 0L;
   } else wss = w->_mp_size;
   
   if UNLIKELY(wss == 0L)
   {
      fmpz_fit_limbs(w, xn+1);
      wp = fmpz_data(w);
      xp = fmpz_data(x);
      cy = mpn_mul_1(wp, xp, xn, y);
      wp[xn] = cy;
      xn += (cy != 0L);
      w->_mp_size = ((long) sub >= 0L ? xn : -xn);
      return;
   }

   sub ^= wss;
   wn = MPIR_ABS(wss);
   new_wn = MPIR_MAX(wn, xn);
   fmpz_fit_limbs(w, new_wn+1);
   wp = fmpz_data(w);
   xp = fmpz_data(x);
   
   min_n = MPIR_MIN(wn, xn);

   if ((long) sub >= 0L)
   {
      // addmul of absolute values 

      cy = mpn_addmul_1 (wp, xp, min_n, y);
      wp += min_n;
      xp += min_n;

      dn = xn - wn;
      if (dn)
      {
         mp_limb_t cy2;
         if ((long) dn > 0L)
            cy2 = mpn_mul_1(wp, xp, dn, y);
         else
         {
            dn = -dn;
            cy2 = 0L;
         }
         cy = cy2 + mpn_add_1(wp, wp, dn, cy);
      }
      wp[dn] = cy;
      new_wn += (cy != 0L);
      wp -= min_n;
   } else
   {
      // submul of absolute values 

      cy = mpn_submul_1(wp, xp, min_n, y);
      if (wn >= xn)
      {
         // if w bigger than x, then propagate borrow through it 
         if (wn != xn)
            cy = mpn_sub_1(wp+xn, wp+xn, wn-xn, cy);

         if (cy != 0)
         {
            // Borrow out of w, take twos complement negative to get
            //   absolute value, flip sign of w.  
            wp[new_wn] = ~-cy;  // extra limb is 0-cy 
            mpn_com_n(wp, wp, new_wn);
            new_wn++;
            MPN_INCR_U(wp, new_wn, (1L));
            wss = -wss;
         }
      }
      else // wsize < xsize 
      {
         // x bigger than w, so want x*y-w.  Submul has given w-x*y, so
         //   take twos complement and use an mpn_mul_1 for the rest.  

         mp_limb_t cy2;

         // -(-cy*b^n + w-x*y) = (cy-1)*b^n + ~(w-x*y) + 1 
         mpn_com_n(wp, wp, wn);
         cy += mpn_add_1(wp, wp, wn, 1L);
         cy -= 1;

         // If cy-1 == -1 then hold that -1 for latter.  mpn_submul_1 never
         //   returns cy==UINT_MAX so that value always indicates a -1. 
         cy2 = (cy == -1L);
         cy += cy2;
         mp_limb_t cy3;
         cy3 = mpn_mul_1(wp+wn, xp+wn, xn-wn, y);               
         cy = cy3 + mpn_add_1(wp+wn, wp+wn, xn-wn, cy);    
         wp[new_wn] = cy;
         new_wn += (cy != 0);

          // Apply any -1 from above.  The value at wp+wsize is non-zero
          //   because y!=0 and the high limb of x will be non-zero.  
         if (cy2)
            MPN_DECR_U (wp+wn, new_wn-wn, (1L));

          wss = -wss;
        }

        // we can have high zero limbs due to cancellation, both when w
        // has more limbs or x has more  
        w->_mp_size = new_wn;
        NORM(w, wp);
    }

    w->_mp_size = ((long) wss >= 0L ? new_wn : -new_wn);
    NORM(w, wp);
}

void fmpz_addmul_ui(fmpz_t * w, fmpz_t * x, ulong y)
{
   if (x->_mp_alloc == 0L)
   {
      long x_int = (long) x->_mp_d;
      ulong xabs = MPIR_ABS(x_int);
      mp_limb_t prod[2];
      umul_ppmm(prod[1], prod[0], xabs, y);
      if ((prod[1] == 0L) && (prod[0] <= IMM_MAX))
      {
         if (x_int < 0L) _fmpz_sub_IMM(w, w, prod[0]);
         else _fmpz_add_IMM(w, w, prod[0]);
      } else
      {
         if (prod[1] == 0L)
         {
            fmpz_fit_limbs(w, MPIR_MAX(fmpz_size(w)+1, 2L));
            if (x_int < 0L) mpz_sub_ui(w, w, prod[0]);
            else mpz_add_ui(w, w, prod[0]);
         } else
         {
            fmpz_fit_limbs(w, MPIR_MAX(fmpz_size(w)+1, 3L));
            __mpz_struct temp;
            temp._mp_d = prod;
            temp._mp_alloc = 2L;
            if (x_int < 0L) temp._mp_size = -2L;
            else temp._mp_size = 2L;
            fmpz_add(w, w, &temp);
         }      
      }
   } else _fmpz_addmul_ui(w, x, y);
}

/*
   w = w - x*y
*/

void _fmpz_submul_ui(fmpz_t * w, fmpz_t * x, ulong y)
{
   ulong  xn, wn, wss, new_wn, min_n, dn;
   mp_limb_t * xp;
   mp_limb_t * wp;
   mp_limb_t cy;
   
   xn = x->_mp_size;
   if UNLIKELY(xn == 0 || y == 0L)
      return;

   mp_limb_t sub = ~xn;
   xn = MPIR_ABS(xn);

   if (w->_mp_alloc == 0L)
   {
      if ((long) w->_mp_d < 0L) wss = -1L;
      else if ((long) w->_mp_d > 0L) wss = 1L;
      else wss = 0L;
   } else wss = w->_mp_size;
   
   if UNLIKELY(wss == 0L)
   {
      fmpz_fit_limbs(w, xn+1);
      wp = fmpz_data(w);
      xp = fmpz_data(x);
      cy = mpn_mul_1(wp, xp, xn, y);
      wp[xn] = cy;
      xn += (cy != 0L);
      w->_mp_size = ((long) sub >= 0L ? xn : -xn);
      return;
   }

   sub ^= wss;
   wn = MPIR_ABS(wss);
   new_wn = MPIR_MAX(wn, xn);
   fmpz_fit_limbs(w, new_wn+1);
   wp = fmpz_data(w);
   xp = fmpz_data(x);
   
   min_n = MPIR_MIN(wn, xn);

   if ((long) sub >= 0L)
   {
      // addmul of absolute values 

      cy = mpn_addmul_1 (wp, xp, min_n, y);
      wp += min_n;
      xp += min_n;

      dn = xn - wn;
      if (dn)
      {
         mp_limb_t cy2;
         if ((long) dn > 0L)
            cy2 = mpn_mul_1(wp, xp, dn, y);
         else
         {
            dn = -dn;
            cy2 = 0L;
         }
         cy = cy2 + mpn_add_1(wp, wp, dn, cy);
      }
      wp[dn] = cy;
      new_wn += (cy != 0L);
      wp -= min_n;
   } else
   {
      // submul of absolute values 

      cy = mpn_submul_1(wp, xp, min_n, y);
      if (wn >= xn)
      {
         // if w bigger than x, then propagate borrow through it 
         if (wn != xn)
            cy = mpn_sub_1(wp+xn, wp+xn, wn-xn, cy);

         if (cy != 0)
         {
            // Borrow out of w, take twos complement negative to get
            //   absolute value, flip sign of w.  
            wp[new_wn] = ~-cy;  // extra limb is 0-cy 
            mpn_com_n(wp, wp, new_wn);
            new_wn++;
            MPN_INCR_U(wp, new_wn, (1L));
            wss = -wss;
         }
      }
      else // wsize < xsize 
      {
         // x bigger than w, so want x*y-w.  Submul has given w-x*y, so
         //   take twos complement and use an mpn_mul_1 for the rest.  

         mp_limb_t cy2;

         // -(-cy*b^n + w-x*y) = (cy-1)*b^n + ~(w-x*y) + 1 
         mpn_com_n(wp, wp, wn);
         cy += mpn_add_1(wp, wp, wn, 1L);
         cy -= 1;

         // If cy-1 == -1 then hold that -1 for latter.  mpn_submul_1 never
         //   returns cy==UINT_MAX so that value always indicates a -1. 
         cy2 = (cy == -1L);
         cy += cy2;
         mp_limb_t cy3;
         cy3 = mpn_mul_1(wp+wn, xp+wn, xn-wn, y);               
         cy = cy3 + mpn_add_1(wp+wn, wp+wn, xn-wn, cy);    
         wp[new_wn] = cy;
         new_wn += (cy != 0);

          // Apply any -1 from above.  The value at wp+wsize is non-zero
          //   because y!=0 and the high limb of x will be non-zero.  
         if (cy2)
            MPN_DECR_U (wp+wn, new_wn-wn, (1L));

          wss = -wss;
        }

        // we can have high zero limbs due to cancellation, both when w
        // has more limbs or x has more  
        w->_mp_size = new_wn;
        NORM(w, wp);
    }

    w->_mp_size = ((long) wss >= 0L ? new_wn : -new_wn);
    NORM(w, wp);
}

void fmpz_submul_ui(fmpz_t * w, fmpz_t * x, ulong y)
{
   if (x->_mp_alloc == 0L)
   {
      long x_int = (long) x->_mp_d;
      ulong xabs = MPIR_ABS(x_int);
      mp_limb_t prod[2];
      umul_ppmm(prod[1], prod[0], xabs, y);
      if ((prod[1] == 0L) && (prod[0] <= IMM_MAX))
      {
         if (x_int < 0L) _fmpz_add_IMM(w, w, prod[0]);
         else _fmpz_sub_IMM(w, w, prod[0]);
      } else
      {
         if (prod[1] == 0L)
         {
            fmpz_fit_limbs(w, MPIR_MAX(fmpz_size(w)+1, 2L));
            if (x_int < 0L) mpz_add_ui(w, w, prod[0]);
            else mpz_sub_ui(w, w, prod[0]);
         } else
         {
            fmpz_fit_limbs(w, MPIR_MAX(fmpz_size(w)+1, 3L));
            __mpz_struct temp;
            temp._mp_d = prod;
            temp._mp_alloc = 2L;
            if (x_int < 0L) temp._mp_size = -2L;
            else temp._mp_size = 2L;
            fmpz_sub(w, w, &temp);
         }      
      }
   } else _fmpz_submul_ui(w, x, y);
}

/* ==============================================================================

   Multiplication

===============================================================================*/

/*
   w = u*2^exp for 0 <= exp
*/

void _fmpz_mul_2exp(fmpz_t * w, fmpz_t * u, ulong exp)
{
   long uss = u->_mp_size;
   ulong un = MPIR_ABS(uss);
   ulong wn;
   ulong limb_cnt;
   mp_limb_t * up;
   mp_limb_t * wp;
   mp_limb_t wlimb;

   if UNLIKELY(uss == 0L)
   {
      if (w->_mp_alloc == 0L) w->_mp_d = (mp_limb_t *) 0L;
      else w->_mp_size = 0L;
      return;
   }

   limb_cnt = (exp>>MPIR_LG_BITS);
   exp &= (MPIR_BITS-1);
  
   fmpz_fit_limbs(w, un + limb_cnt + (exp != 0L));
   
   up = fmpz_data(u);
   wp = fmpz_data(w);
   wn = un + limb_cnt;

   if (exp)
   {
      wlimb = mpn_lshift(wp + limb_cnt, up, un, exp);
      if (wlimb)
	  {
	     wp[wn] = wlimb;
	     wn++;
      }
   } else
   {
      MPN_COPY_INCR(wp + limb_cnt, up, un);
   }

   // Zero all whole limbs at low end.  Do it here and not before calling
   // mpn_lshift, so as not to lose data when U == W.  
   F_mpn_clear(wp, limb_cnt);

   w->_mp_size = (uss >= 0L ? wn : -wn);
}

void fmpz_mul_2exp(fmpz_t * w, fmpz_t * u, ulong exp)
{
   if (u->_mp_alloc == 0L)
   {
      long u_int = (long) u->_mp_d; 
      ulong u_abs = MPIR_ABS(u_int);
      ulong bits = MPIR_BIT_COUNT(u_abs);
      if (bits + exp <= MPIR_BITS - 2) 
      {
         fmpz_set_si(w, u_int<<exp);
      } else 
      {
         __mpz_struct temp;
         temp._mp_d = &u_abs;
         if (u_int < 0L) temp._mp_size = -1L;
         else if (u_int > 0L) temp._mp_size = 1L;
         else temp._mp_size = 0L;
         temp._mp_alloc = 2L;
         _fmpz_mul_2exp(w, &temp, exp);
      }
   } else 
   {
      _fmpz_mul_2exp(w, u, exp);
   }
}

// *************** end of file
