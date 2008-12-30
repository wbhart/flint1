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

   packed_vec.c: vectors with packed representation of either 8, 16, 32 
	              or 64 bits per entry (FLINT 2.0)

   Copyright (C) 2008, William Hart

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "packed_vec.h"
#include "flint.h"
#include "memory-manager.h"

void pv_init(pv_s * vec, ulong entries, int bits)
{
   if (entries && bits)
	{
		ulong limbs = (entries*bits - 1)/FLINT_BITS + 1;
		vec->alloc = (limbs*FLINT_BITS)/bits;
		vec->entries = (mp_limb_t *) flint_heap_alloc(limbs);
		vec->bits = bits;
#if BIT_FIDDLE
		vec->log_bits = FLINT_BIT_COUNT(bits) - 1;
		vec->pack = FLINT_BITS/bits;
		vec->log_pack = FLINT_BIT_COUNT(vec->pack) - 1;
		vec->length = 0;
#endif
	} else
	{
		vec->entries = NULL;
		vec->alloc = entries;
		vec->bits = bits;
		vec->length = 0;
	}
}

void pv_realloc(pv_s * vec, ulong entries)
{
	if (entries)
	{
		ulong limbs = (entries*vec->bits - 1)/FLINT_BITS + 1;
		   
		if (vec->entries)
		{
			if (entries == vec->alloc) return;
			
			vec->entries = (mp_limb_t *) flint_heap_realloc(vec->entries, limbs);
		} else
         vec->entries = (mp_limb_t *) flint_heap_alloc(limbs);

		vec->alloc = (limbs*FLINT_BITS)/vec->bits;
	} else
	{
		if (vec->entries) flint_heap_free(vec->entries);
		vec->entries = NULL;
		vec->alloc = 0;
	}
}

void pv_set_bits(pv_s * vec, int bits)
{
	if (bits == vec->bits) return; // nothing to do

	if (bits == FLINT_BITS)
	{
		pv_iter_s iter;
      PV_ITER_INIT(iter, *vec, 0);

		mp_limb_t * temp = (mp_limb_t *) flint_heap_alloc(vec->alloc);

		for (ulong i = 0; i < vec->length; i++)
		   PV_GET_NEXT(temp[i], iter);

		if (vec->entries) flint_heap_free(vec->entries);
		vec->entries = temp;
		vec->alloc = vec->alloc;

		vec->bits = bits;
#if BIT_FIDDLE
		vec->log_bits = FLINT_LG_BITS_PER_LIMB;
		vec->pack = 1;
		vec->log_pack = 0;
#endif
	} else if (vec->bits == FLINT_BITS)
	{ 
		mp_limb_t * temp = vec->entries;
			
		ulong limbs = (vec->alloc*bits - 1)/FLINT_BITS + 1;
		vec->entries = (mp_limb_t *) flint_heap_alloc(limbs);
      vec->alloc = (limbs*FLINT_BITS)/bits;

		vec->bits = bits;
#if BIT_FIDDLE
		vec->log_bits = FLINT_BIT_COUNT(bits) - 1;
		vec->pack = FLINT_BITS/bits;
		vec->log_pack = FLINT_BIT_COUNT(vec->pack) - 1;
#endif
      
		pv_iter_s iter;
      PV_ITER_INIT(iter, *vec, 0);

		for (ulong i = 0; i < vec->length; i++)
		{
			PV_SET_NEXT(iter, temp[i]);
		}
		
		flint_heap_free(temp);
	} else
	{
		pv_iter_s iter1, iter2;
		PV_ITER_INIT(iter1, *vec, 0);

		ulong limbs = (vec->alloc*bits - 1)/FLINT_BITS + 1;
		vec->entries = (mp_limb_t *) flint_heap_alloc(limbs);
      vec->alloc = (limbs*FLINT_BITS)/bits;

		vec->bits = bits;
#if BIT_FIDDLE
		vec->log_bits = FLINT_BIT_COUNT(bits) - 1;
		vec->pack = FLINT_BITS/bits;
		vec->log_pack = FLINT_BIT_COUNT(vec->pack) - 1;
#endif
		PV_ITER_INIT(iter2, *vec, 0);
     
      ulong temp;
		
		for (ulong i = 0; i < vec->length; i++)
		{   
			PV_GET_NEXT(temp, iter1);
			PV_SET_NEXT(iter2, temp);
		}
      
		if (vec->entries) flint_heap_free(iter1.entries);
	}
}


// *************** end of file
