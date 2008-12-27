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
/****************************************************************************

   packed_vec.h: vectors with packed representation of either 8, 16, 32 
	              or 64 bits per entry (FLINT 2.0)

   Copyright (C) 2008, William Hart 

*****************************************************************************/

#ifndef FLINT_PV_H
#define FLINT_PV_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "memory-manager.h"
#include "flint.h"
#include "mpn_extras.h"

#define BIT_FIDDLE 1

#if BIT_FIDDLE
	 
/*
   A structure representing a packed vector with either 8, 16, 32
	or 64 bits per entry.
	"bits" is the number of bits per entry.
	"log_bits" is log_2 of the number of bits per entry.
	"alloc" is the number of entries currently allocated.
	"length" is the number of entries actually used.
	"entries" is a pointer to the actual limbs storing the entries.
	"pack" is the number of entries packed into each limb.
	"log_pack" is log_2 of the number of entries per limb.
*/

typedef struct
{
   mp_limb_t * entries;
	ulong alloc;
	ulong length;
	int bits;
	int log_bits;
	int pack; 
	int log_pack; 
} pv_s;

/*
  A structure for iterating through entries in a packed vector.
  "limb" is the current limb of the array representing the vector,
  "mask" is a bit mask which masks the current entry
  "shift" is the shift required to get the entry to the least 
  significant bits of a limb.
  The other parameters are simply copies of the corresponding 
  parameters from the pv_s for this iterator.
 */

typedef struct
{
   mp_limb_t * entries;
	ulong limb;
   ulong mask;
	int shift;
	int bits;
	int log_bits;
	int pack;
	int log_pack;
} pv_iter_s;

/*
   Initialise an iterator "iter_xx" for a packed vector "pv_xxx" to 
	correspond to entry "entry_xxx"
*/

#define PV_ITER_INIT(pv_xxx, iter_xxx, entry_xxx) \
do { \
  iter_xxx.limb = (entry_xxx >> (pv_xxx).log_pack); \
  iter_xxx.shift = ((entry_xxx & ((pv_xxx).pack - 1)) << (pv_xxx).log_bits); \
  if ((pv_xxx).bits == FLINT_BITS) iter_xxx.mask = -1L; \
  else iter_xxx.mask = ((1UL << (pv_xxx).bits) - 1UL) << iter_xxx.shift; \
  iter_xxx.entries = (pv_xxx).entries; \
  iter_xxx.bits = (pv_xxx).bits; \
  iter_xxx.log_bits = (pv_xxx).log_bits; \
  iter_xxx.pack = (pv_xxx).pack; \
  iter_xxx.log_pack = (pv_xxx).log_pack; \
} while (0)

/*
   Given an iterator iter_xxx, set xxx to the current entry 
	of the vector and iterate to the next entry.
*/

#define PV_GET_NEXT(xxx, iter_xxx) \
do { \
   if (iter_xxx.bits == FLINT_BITS) \
   { \
	   xxx = iter_xxx.entries[iter_xxx.limb]; \
		iter_xxx.limb++; \
   } else \
	{ \
		xxx = ((iter_xxx.entries[iter_xxx.limb] & iter_xxx.mask) >> iter_xxx.shift); \
      iter_xxx.mask <<= iter_xxx.bits; \
	   if (iter_xxx.mask) \
	      iter_xxx.shift += iter_xxx.bits; \
	   else \
	   { \
         iter_xxx.mask = ((1UL << iter_xxx.bits) - 1UL); \
		   iter_xxx.shift = 0; \
		   iter_xxx.limb++; \
	   } \
	} \
} while (0)

/*
   Given an iterator iter_xxx, set the current entry of the 
	vector to xxx and iterate to the next entry.
*/

#define PV_SET_NEXT(iter_xxx, xxx) \
do { \
   if (iter_xxx.bits == FLINT_BITS) \
   { \
	   iter_xxx.entries[iter_xxx.limb] = xxx; \
		iter_xxx.limb++; \
   } else \
	{ \
		iter_xxx.entries[iter_xxx.limb] = (iter_xxx.entries[iter_xxx.limb] & ~iter_xxx.mask) + (xxx << iter_xxx.shift); \
      iter_xxx.mask <<= iter_xxx.bits; \
	   if (iter_xxx.mask) \
	      iter_xxx.shift += iter_xxx.bits; \
	   else \
	   { \
         iter_xxx.mask = ((1UL << iter_xxx.bits) - 1UL); \
		   iter_xxx.shift = 0; \
		   iter_xxx.limb++; \
	   } \
	} \
} while (0)

/*
   Given an iterator iter_xxx, set xxx to the current entry of the 
	vector and iterate to the previous entry.
*/

#define PV_GET_PREV(xxx, iter_xxx) \
do { \
   if (iter_xxx.bits == FLINT_BITS) \
   { \
	   xxx = iter_xxx.entries[iter_xxx.limb]; \
		iter_xxx.limb--; \
   } else \
	{ \
		xxx = ((iter_xxx.entries[iter_xxx.limb] & iter_xxx.mask) >> iter_xxx.shift); \
      iter_xxx.mask >>= iter_xxx.bits; \
	   if (iter_xxx.mask) \
	      iter_xxx.shift -= iter_xxx.bits; \
	   else \
	   { \
         iter_xxx.mask = (((1UL << iter_xxx.bits) - 1UL) << (FLINT_BITS - iter_xxx.bits)); \
		   iter_xxx.shift = FLINT_BITS - iter_xxx.bits; \
		   iter_xxx.limb--; \
	   } \
	} \
} while (0)

/*
   Given a packed vector pv_xxx and an iterator iter_xxx, set
   the current entry of the vector to xxx and iterate to the
	previous entry.
*/

#define PV_SET_PREV(iter_xxx, xxx) \
do { \
   if (iter_xxx.bits == FLINT_BITS) \
   { \
	   iter_xxx.entries[iter_xxx.limb] = xxx; \
		iter_xxx.limb--; \
   } else \
	{ \
		iter_xxx.entries[iter_xxx.limb] = (iter_xxx.entries[iter_xxx.limb] & ~iter_xxx.mask) + (xxx << iter_xxx.shift); \
      iter_xxx.mask >>= iter_xxx.bits; \
	   if (iter_xxx.mask) \
	      iter_xxx.shift -= iter_xxx.bits; \
	   else \
	   { \
         iter_xxx.mask = (((1UL << iter_xxx.bits) - 1UL) << (FLINT_BITS - iter_xxx.bits)); \
		   iter_xxx.shift = FLINT_BITS - iter_xxx.bits; \
		   iter_xxx.limb--; \
	   } \
	} \
} while (0)

/*
   Given a packed vector "pv_xxx" and an entry "entry_xxx", set
	xxx to the given entry of the vector.
*/

#define PV_GET_ENTRY(xxx, pv_xxx, entry_xxx) \
do { \
  int xxx_limb = (entry_xxx >> pv_xxx.log_pack); \
  int xxx_shift = ((entry_xxx & (pv_xxx.pack - 1)) << pv_xxx.log_bits); \
  ulong xxx_mask; \
  if (pv_xxx.bits == FLINT_BITS) xxx_mask = -1L; \
  else xxx_mask = (1UL << pv_xxx.bits) - 1UL; \
  xxx = ((pv_xxx.entries[xxx_limb] >> xxx_shift) & xxx_mask); \
} while (0)

/*
   Given a packed vector "pv_xxx" set the given entry "entry_xxx" to 
	xxx.
*/

#define PV_SET_ENTRY(pv_xxx, entry_xxx, xxx) \
do { \
  int xxx_limb = (entry_xxx >> pv_xxx.log_pack); \
  int xxx_shift = ((entry_xxx & (pv_xxx.pack - 1)) << pv_xxx.log_bits); \
  ulong xxx_mask; \
  if (pv_xxx.bits == FLINT_BITS) xxx_mask = -1L; \
  else xxx_mask = ((1UL << pv_xxx.bits) - 1UL) << xxx_shift); \
  pv_xxx.entries[xxx_limb] = (pv_xxx.entries[xxx_limb] & ~xxx_mask) + (xxx << xxx_shift); \
} while (0)

#else

/*
   A structure representing a packed vector with either 8, 16, 32
	or 64 bits per entry.
	"bits" is the number of bits per entry.
	"alloc" is the number of entries currently allocated.
	"length" is the number of entries actually used.
	"entries" is a pointer to the actual limbs storing the entries.
*/

typedef struct
{
   mp_limb_t * entries;
	ulong alloc;
	ulong length;
	int bits;
} pv_s;

/*
  A structure for iterating through entries in a packed vector.
  "limb" is the current limb of the array representing the vector,
  "mask" is a bit mask which masks the current entry
  "shift" is the shift required to get the entry to the least 
  significant bits of a limb.
  The other parameters are simply copies of the corresponding 
  parameters from the pv_s for this iterator.
 */

typedef struct
{
   mp_limb_t * entries;
	ulong i;
	int bits;
} pv_iter_s;

/*
   Initialise an iterator "iter_xx" for a packed vector "pv_xxx" to 
	correspond to entry "entry_xxx"
*/

#define PV_ITER_INIT(pv_xxx, iter_xxx, entry_xxx) \
do { \
  iter_xxx.entries = (pv_xxx).entries; \
  iter_xxx.bits = (pv_xxx).bits; \
  iter_xxx.i = entry_xxx; \
} while (0)

/*
   Given an iterator iter_xxx, set xxx to the current entry 
	of the vector and iterate to the next entry.
*/

#define PV_GET_NEXT(xxx, iter_xxx) \
do { \
   int bits = iter_xxx.bits; \
	if (bits == 8) (xxx) = (ulong) ((uint8_t *) iter_xxx.entries)[iter_xxx.i]; \
	else if (bits == 16) (xxx) = (ulong) ((uint16_t *) iter_xxx.entries)[iter_xxx.i]; \
	else if (bits == 32) (xxx) = (ulong) ((uint32_t *) iter_xxx.entries)[iter_xxx.i]; \
   else if (bits == 64) (xxx) = (ulong) ((uint64_t *) iter_xxx.entries)[iter_xxx.i]; \
	iter_xxx.i++; \
} while (0)

/*
   Given an iterator iter_xxx, set the current entry of the 
	vector to xxx and iterate to the next entry.
*/

#define PV_SET_NEXT(iter_xxx, xxx) \
do { \
   int bits = iter_xxx.bits; \
	if (bits == 8) ((uint8_t *) iter_xxx.entries)[iter_xxx.i] = (uint8_t) (xxx); \
	else if (bits == 16) ((uint16_t *) iter_xxx.entries)[iter_xxx.i] = (uint16_t) (xxx); \
	else if (bits == 32) ((uint32_t *) iter_xxx.entries)[iter_xxx.i] = (uint32_t) (xxx); \
	else if (bits == 64) ((uint64_t *) iter_xxx.entries)[iter_xxx.i] = (uint64_t) (xxx); \
	iter_xxx.i++; \
} while (0)

/*
   Given an iterator iter_xxx, set xxx to the current entry of the 
	vector and iterate to the previous entry.
*/

#define PV_GET_PREV(xxx, iter_xxx) \
do { \
   int bits = iter_xxx.bits; \
	if (bits == 8) (xxx) = (ulong) ((uint8_t *) iter_xxx.entries)[iter_xxx.i]; \
	else if (bits == 16) (xxx) = (ulong) ((uint16_t *) iter_xxx.entries)[iter_xxx.i]; \
	else if (bits == 32) (xxx) = (ulong) ((uint32_t *) iter_xxx.entries)[iter_xxx.i]; \
	else if (bits == 64) (xxx) = (ulong) ((uint64_t *) iter_xxx.entries)[iter_xxx.i]; \
	iter_xxx.i--; \
} while (0)

/*
   Given a packed vector pv_xxx and an iterator iter_xxx, set
   the current entry of the vector to xxx and iterate to the
	previous entry.
*/

#define PV_SET_PREV(iter_xxx, xxx) \
do { \
   int bits = iter_xxx.bits; \
	if (bits == 8) ((uint8_t *) iter_xxx.entries)[iter_xxx.i] = (uint8_t) (xxx); \
	else if (bits == 16) ((uint16_t *) iter_xxx.entries)[iter_xxx.i] = (uint16_t) (xxx); \
	else if (bits == 32) ((uint32_t *) iter_xxx.entries)[iter_xxx.i] = (uint32_t) (xxx); \
	else if (bits == 64) ((uint64_t *) iter_xxx.entries)[iter_xxx.i] = (uint64_t) (xxx); \
	iter_xxx.i--; \
} while (0)

/*
   Given a packed vector "pv_xxx" and an entry "entry_xxx", set
	xxx to the given entry of the vector.
*/

#define PV_GET_ENTRY(xxx, pv_xxx, entry_xxx) \
do { \
   int bits = pv_xxx.bits; \
	if (bits == 8) (xxx) = (ulong) ((uint8_t *) pv_xxx.entries)[entry_xxx]; \
	else if (bits == 16) (xxx) = (ulong) ((uint16_t *) pv_xxx.entries)[entry_xxx]; \
	else if (bits == 32) (xxx) = (ulong) ((uint32_t *) pv_xxx.entries)[entry_xxx]; \
	else if (bits == 64) (xxx) = (ulong) ((uint64_t *) pv_xxx.entries)[entry_xxx]; \
} while (0)

/*
   Given a packed vector "pv_xxx" set the given entry "entry_xxx" to 
	xxx.
*/

#define PV_SET_ENTRY(pv_xxx, entry_xxx, xxx) \
do { \
   int bits = pv_xxx.bits; \
	if (bits == 8) ((uint8_t *) pv_xxx.entries)[entry_xxx] = (uint8_t) (xxx); \
	else if (bits == 16) ((uint16_t *) pv_xxx.entries)[entry_xxx] = (uint16_t) (xxx); \
	else if (bits == 32) ((uint32_t *) pv_xxx.entries)[entry_xxx] = (uint32_t) (xxx); \
	else if (bits == 64) ((uint64_t *) pv_xxx.entries)[entry_xxx] = (uint64_t) (xxx); \
} while (0)

#endif

static inline 
int pv_bit_fit(int bits)
{
	int b = FLINT_BITS;
#if FLINT_BITS == 64
	if (bits <= 32) b = 32;
#endif
	if (bits <= 16) b = 16;
   if (bits <= 8) b = 8;
	return b;
}

/*
   Initialise a packed vector pointed to by vec with space for at least the
	given number of entries and with the given number of bits per entry.
*/

void pv_init(pv_s * vec, ulong entries, int bits);

/*
   Reallocate the packed vector pointed to by vec to have space for at least the
	given number of entries.
*/

void pv_realloc(pv_s * vec, ulong entries);

/*
   If entries is greater than the number of allocated entries this function will
	realloc the given vector, at least doubling the amount of allocated entries
	if reallocation occurs and allowing space for the given number of entries.
*/

static inline
void pv_fit_length(pv_s * vec, ulong entries)
{
	if (entries == 0) return;

	if (entries > vec->alloc)
	{
		if (entries < 2*vec->alloc) pv_realloc(vec, 2*vec->alloc);
		else pv_realloc(vec, entries);
	}
}

/* 
   Clear the given packed vector, free'ing any memory used by it.
*/

static inline
void pv_clear(pv_s * vec)
{
	if (vec->entries) flint_heap_free(vec->entries);
	vec->entries = NULL;
	vec->alloc = 0;
   vec->length = 0;
}

/* 
   Set the number of bits per entry to the given number of bits,
	preserving all the data. The packed vector is reallocated if
	necessary.
*/

void pv_set_bits(pv_s * vec, int bits);

#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
