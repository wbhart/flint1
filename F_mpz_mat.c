/*============================================================================

    F_mpz_mat.c: Matrices over Z (FLINT 2.0 matrices)

    Copyright (C) 2008, William Hart 

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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "flint.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz_mat.h"
#include "mpz_mat.h"

/*===============================================================================

	Memory management

================================================================================*/

void _F_mpz_mat_mpz_entries_new(F_mpz_mat_t mat)
{
	if (mat->mpz_length == mat->mpz_alloc) // time to allocate MPZ_BLOCK_MAT more mpz_t's
	{
	   if (mat->mpz_alloc) // realloc mpz_t's
			mat->mpz_entries = (__mpz_struct*) flint_heap_realloc_bytes(mat->mpz_entries, (mat->mpz_alloc + MPZ_BLOCK_MAT)*sizeof(__mpz_struct));
		else // first time alloc of mpz_t'
			mat->mpz_entries = (__mpz_struct*) flint_heap_alloc_bytes(MPZ_BLOCK_MAT*sizeof(__mpz_struct));	
		
		// initialise the new mpz_t's
		for (ulong i = mat->mpz_alloc; i < mat->mpz_alloc + MPZ_BLOCK_MAT; i++)
			mpz_init(mat->mpz_entries + i);
		mat->mpz_alloc += MPZ_BLOCK_MAT;
	}

	mat->mpz_length++;
}

void _F_mpz_mat_mpz_entries_clear(F_mpz_mat_t mat)
{
   for (ulong i = 0; i < mat->mpz_alloc; i++) // clear any initialised mpz_t's
	   mpz_clear(mat->mpz_entries + i);

	flint_heap_free(mat->mpz_entries); // clear mpz_t array itself
	mat->mpz_entries = NULL;
	mat->mpz_alloc = 0;
	mat->mpz_length = 0;
}

void F_mpz_mat_init(F_mpz_mat_t mat, const ulong r, const ulong c)
{
   if ((r) && (c)) // allocate space for r*c small entries
   {
      mat->entries = (mp_limb_t *) flint_heap_alloc(r*c);
		F_mpn_clear(mat->entries, r*c);
		mat->rows = (mp_limb_t **) flint_heap_alloc(r);
		for (ulong i = 0; i < r; i++)
		   mat->rows[i] = mat->entries + i*c;
   } else mat->entries = NULL;
   
	mat->mpz_entries = NULL;

   mat->r = r;
	mat->c = c;
   mat->mpz_alloc = 0;
   mat->mpz_length = 0;
}

void F_mpz_mat_clear(F_mpz_mat_t mat)
{
   if (mat->entries) 
	{
		flint_heap_free(mat->entries); // clean up ordinary entries
		flint_heap_free(mat->rows);
	}
   if (mat->mpz_entries) _F_mpz_mat_mpz_entries_clear(mat); // clean up mpz_t entries
   mat->entries = NULL;
	mat->r = 0;
	mat->c = 0;
}

/*===============================================================================

	Coefficient operations

================================================================================*/

void _F_mpz_entry_demote_val(F_mpz_mat_t mat, const ulong r, const ulong c)
{
   mp_limb_t * row = mat->rows[r];
	
	__mpz_struct * mpz_ptr = mat->mpz_entries + ENTRY_TO_OFF(row[c]);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L) // entry is zero
	{
		_F_mpz_entry_zero(mat, r, c);
	} else if (size == 1L) // entry is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= ENTRY_MAX) row[c] = uval;
	} else if (size == -1L) // entry is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= ENTRY_MAX) row[c] = -uval;
	}
	// don't do anything if entry has to be multi precision
}

void _F_mpz_entry_set_si(F_mpz_mat_t mat, const ulong r, const ulong c, const long val)
{
   ulong uval = FLINT_ABS(val);
	
	if (uval > ENTRY_MAX) // val is large
	{
		__mpz_struct * mpz_entry = _F_mpz_entry_promote(mat, r, c);
		mpz_set_si(mpz_entry, val);
	} else mat->rows[r][c] = val; // val is small
}

void _F_mpz_entry_set_ui(F_mpz_mat_t mat, const ulong r, const ulong c, const ulong val)
{
   if (val > ENTRY_MAX) // val is large
	{
		__mpz_struct * mpz_entry = _F_mpz_entry_promote(mat, r, c);
		mpz_set_ui(mpz_entry, val);
	} else mat->rows[r][c] = val; // val is small
}

long _F_mpz_entry_get_si(const F_mpz_mat_t mat, const ulong r, const ulong c)
{
   ulong d = mat->rows[r][c];

	if (!ENTRY_IS_MPZ(d)) return d; // entry is small
	return mpz_get_si(mat->mpz_entries + ENTRY_TO_OFF(d)); // entry is large
}

ulong _F_mpz_entry_get_ui(const F_mpz_mat_t mat, const ulong r, const ulong c)
{
   ulong d = mat->rows[r][c];

	if (!ENTRY_IS_MPZ(d)) return d; // entry is small
	return mpz_get_ui(mat->mpz_entries + ENTRY_TO_OFF(d)); //entry is large
}

void _F_mpz_entry_get_mpz(mpz_t x, const F_mpz_mat_t mat, const ulong r, const ulong c)
{
	mp_limb_t d = mat->rows[r][c];

	if (!ENTRY_IS_MPZ(d)) mpz_set_si(x, d); // set x to small entry
	else mpz_set(x, mat->mpz_entries + ENTRY_TO_OFF(d)); // set x to large entry
}

extern double __gmpn_get_d(mp_limb_t *, size_t, size_t, long);

double _F_mpz_entry_get_d_2exp(long * exp, const F_mpz_mat_t mat, const ulong r, const ulong c)
{
   mp_limb_t d = mat->rows[r][c];

	if (!ENTRY_IS_MPZ(d))
   {
      if (d == 0L) 
      {
         (*exp) = 0L;
         return 0.0;
      }
      ulong d_abs = FLINT_ABS(d);
      (*exp) = FLINT_BIT_COUNT(d_abs);
      if ((long) d < 0L) return __gmpn_get_d(&d_abs, 1L, -1L, -*exp);
      else return __gmpn_get_d(&d, 1L, 1L, -*exp);
   } else 
	   return mpz_get_d_2exp(exp, mat->mpz_entries + ENTRY_TO_OFF(d));
}

void _F_mpz_entry_set_mpz(F_mpz_mat_t mat, const ulong r, const ulong c, const mpz_t x)
{
   long size = x->_mp_size;
	
	if (size == 0L) // x is zero
	{
		_F_mpz_entry_zero(mat, r, c);
	} else if (size == 1L) // x is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(x);
		if (uval <= ENTRY_MAX) mat->rows[r][c] = uval;
		else 
		{
			__mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat, r, c);
			mpz_set_ui(mpz_ptr, uval);
		}
	} else if (size == -1L) // x is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(x);
		if (uval <= ENTRY_MAX) mat->rows[r][c] = -uval;
		else 
		{
			__mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat, r, c);
			mpz_set_ui(mpz_ptr, uval);
			mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // x is more than one limb
	{
		__mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat, r, c);
		mpz_set(mpz_ptr, x);
	}			
}

void _F_mpz_entry_set(F_mpz_mat_t mat1, ulong r1, ulong c1, const F_mpz_mat_t mat2, const ulong r2, const ulong c2)
{
   ulong d = mat2->rows[r2][c2];
   
	if (!ENTRY_IS_MPZ(d)) // entry is small
	{
		mat1->rows[r1][c1] = d;
	} else // entry is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
		mpz_set(mpz_ptr, mat2->mpz_entries + ENTRY_TO_OFF(d));
	}
}

int _F_mpz_entry_equal(const F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, const ulong r2, const ulong c2)
{
	ulong d1 = mat1->rows[r1][c1];
   ulong d2 = mat2->rows[r2][c2];

	if (!ENTRY_IS_MPZ(d1)) return (d1 == d2); // if c2 is large it can't be equal to c1
	else if (!ENTRY_IS_MPZ(d2)) return 0; // c1 is large, so if c2 isn't....
	else return (!mpz_cmp(mat1->mpz_entries + ENTRY_TO_OFF(d1), mat2->mpz_entries + ENTRY_TO_OFF(d2))); 
}

void _F_mpz_entry_swap(F_mpz_mat_t mat1, const ulong r1, const ulong c1, F_mpz_mat_t mat2, const ulong r2, const ulong c2)
{
	ulong d1 = mat1->rows[r1][c1];
   ulong d2 = mat2->rows[r2][c2];

   if (!ENTRY_IS_MPZ(d1))
	{
		if (!ENTRY_IS_MPZ(d2)) // both entries are small
		{
	      mat1->rows[r1][c1] = d2;
	      mat2->rows[r2][c2] = d1;
		} else // c1 is small, c2 is large
		{
			__mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
			mpz_swap(mpz_ptr, mat2->mpz_entries + ENTRY_TO_OFF(d2));
			mat2->rows[r2][c2] = d1;
		}
	} else
	{
      if (!ENTRY_IS_MPZ(d2)) // c2 is small, c1 is large
		{
         __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat2, r2, c2);
			mpz_swap(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(d1));
			mat1->rows[r1][c1] = d2;
		} else // both entries are large
		{
			mpz_swap(mat1->mpz_entries + ENTRY_TO_OFF(d1), mat2->mpz_entries + ENTRY_TO_OFF(d2));
		}
	}
}

void _F_mpz_entry_neg(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, const ulong r2, const ulong c2)
{
   ulong d = mat2->rows[r2][c2];

	if (!ENTRY_IS_MPZ(d)) // entry is small
	{
		mat1->rows[r1][c1] = -d;
	} else // entry is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
		mpz_neg(mpz_ptr, mat2->mpz_entries + ENTRY_TO_OFF(d));
	}
}

void _F_mpz_entry_add(F_mpz_mat_t res, const ulong r3, const ulong c3, const F_mpz_mat_t mat1, const ulong r1, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r2, const ulong c2) 
{
	mp_limb_t d1 = mat1->rows[r1][c1];
	mp_limb_t d2 = mat2->rows[r2][c2];
	
	if (!ENTRY_IS_MPZ(d1))
	{
	   if (!ENTRY_IS_MPZ(d2)) // both entries are small
		{
			_F_mpz_entry_set_si(res, r3, c3, d1 + d2);
		} else // d1 is small, d2 is large
		{
         __mpz_struct * mpz3 = _F_mpz_entry_promote(res, r3, c3);
			__mpz_struct * mpz2 = mat2->mpz_entries + ENTRY_TO_OFF(d2);
			if ((long) d1 < 0L) mpz_sub_ui(mpz3, mpz2, -d1);	
		   else mpz_add_ui(mpz3, mpz2, d1);
			_F_mpz_entry_demote_val(res, r3, c3); // entries may have cancelled
		}
	} else
	{
		if (!ENTRY_IS_MPZ(d2)) // c1 is large, c2 is small
		{
         __mpz_struct * mpz3 = _F_mpz_entry_promote(res, r3, c3);
			__mpz_struct * mpz1 = mat1->mpz_entries + ENTRY_TO_OFF(d1);
			if ((long) d2 < 0L) mpz_sub_ui(mpz3, mpz1, -d2);	
			else mpz_add_ui(mpz3, mpz1, d2);
			_F_mpz_entry_demote_val(res, r3, c3); // entries may have cancelled
		} else // c1 and c2 are large
		{
         __mpz_struct * mpz3 = _F_mpz_entry_promote(res, r3, c3);
			__mpz_struct * mpz1 = mat1->mpz_entries + ENTRY_TO_OFF(d1);
			__mpz_struct * mpz2 = mat2->mpz_entries + ENTRY_TO_OFF(d2);
			mpz_add(mpz3, mpz1, mpz2);
			_F_mpz_entry_demote_val(res, r3, c3); // entries may have cancelled
		}
	}
}

void _F_mpz_entry_sub(F_mpz_mat_t res, const ulong r3, const ulong c3, const F_mpz_mat_t mat1, const ulong r1, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r2, const ulong c2) 
{
	mp_limb_t d1 = mat1->rows[r1][c1];
	mp_limb_t d2 = mat2->rows[r2][c2];
	
	if (!ENTRY_IS_MPZ(d1))
	{
	   if (!ENTRY_IS_MPZ(d2)) // both entries are small
		{
			_F_mpz_entry_set_si(res, r3, c3, d1 - d2);
		} else // d1 is small, d2 is large
		{
         __mpz_struct * mpz3 = _F_mpz_entry_promote(res, r3, c3);
			__mpz_struct * mpz2 = mat2->mpz_entries + ENTRY_TO_OFF(d2);
			if ((long) d1 < 0L) 
			{
				mpz_add_ui(mpz3, mpz2, -d1);	
				mpz_neg(mpz3, mpz3);
			} else mpz_ui_sub(mpz3, d1, mpz2);
			_F_mpz_entry_demote_val(res, r3, c3); // entries may have cancelled
		}
	} else
	{
		if (!ENTRY_IS_MPZ(d2)) // d1 is large, d2 is small
		{
         __mpz_struct * mpz3 = _F_mpz_entry_promote(res, r3, c3);
			__mpz_struct * mpz1 = mat1->mpz_entries + ENTRY_TO_OFF(d1);
			if ((long) d2 < 0L) mpz_add_ui(mpz3, mpz1, -d2);	
			else mpz_sub_ui(mpz3, mpz1, d2);
			_F_mpz_entry_demote_val(res, r3, c3); // entries may have cancelled
		} else // c1 and c2 are large
		{
         __mpz_struct * mpz3 = _F_mpz_entry_promote(res, r3, c3);
			__mpz_struct * mpz1 = mat1->mpz_entries + ENTRY_TO_OFF(d1);
			__mpz_struct * mpz2 = mat2->mpz_entries + ENTRY_TO_OFF(d2);
			mpz_sub(mpz3, mpz1, mpz2);
			_F_mpz_entry_demote_val(res, r3, c3); // entries may have cancelled
		}
	}
}

void _F_mpz_entry_mul_ui(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, 
						                                   const ulong r2, const ulong c2, const ulong x)
{
	ulong d2 = mat2->rows[r2][c2];

	if (!ENTRY_IS_MPZ(d2)) // entry2 is small
	{
		mp_limb_t prod[2];
		ulong ud2 = FLINT_ABS(d2);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], ud2, x); 
		if (!prod[1]) // result fits in one limb
		{
			_F_mpz_entry_set_ui(mat1, r1, c1, prod[0]);
			if ((long) d2 < 0L) _F_mpz_entry_neg(mat1, r1, c1, mat1, r1, c1);
		} else // result takes two limbs
		{
		   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
			// two limbs, least significant first, native endian, no nails, stored in prod
         mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if ((long) d2 < 0L) mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // entry2 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
      mpz_mul_ui(mpz_ptr, mat2->mpz_entries + ENTRY_TO_OFF(d2), x);
	}
}

void _F_mpz_entry_mul_si(F_mpz_mat_t mat1, const ulong r1, const ulong c1, const F_mpz_mat_t mat2, 
						                                   const ulong r2, const ulong c2, const long x)
{
	ulong d2 = mat2->rows[r2][c2];

	if (!ENTRY_IS_MPZ(d2)) // entry2 is small
	{
		mp_limb_t prod[2];
		ulong ud2 = FLINT_ABS(d2);
		ulong ux = FLINT_ABS(x);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], ud2, ux); 
		if (!prod[1]) // result fits in one limb
		{
			_F_mpz_entry_set_ui(mat1, r1, c1, prod[0]);
			if ((long) (d2 ^ x) < 0L) _F_mpz_entry_neg(mat1, r1, c1, mat1, r1, c1);
		} else // result takes two limbs
		{
		   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
         // two limbs, least significant first, native endian, no nails, stored in prod
			mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if ((long) (d2 ^ x) < 0L) mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // entry2 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);
      mpz_mul_si(mpz_ptr, mat2->mpz_entries + ENTRY_TO_OFF(d2), x);
	}
}

void _F_mpz_entry_mul_mpz(F_mpz_mat_t mat1, ulong r1, ulong c1, F_mpz_mat_t mat2, ulong r2, ulong c2, mpz_t x)
{
	ulong d2 = mat2->rows[r2][c2];
   
	if (mpz_size(x) <= 1) // 1 limb to multiply by
	{
		long x_limb = mpz_get_ui(x);
      _F_mpz_entry_mul_ui(mat1, r1, c1, mat2, r2, c2, x_limb); 
		if (mpz_sgn(x) < 0) _F_mpz_entry_neg(mat1, r1, c1, mat1, r1, c1);
		return;
	}

   // more than one limb to multiply by
	__mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat1, r1, c1);

	if (!ENTRY_IS_MPZ(d2)) // entry2 is small
	   mpz_mul_si(mpz_ptr, x, d2);
	else // entry2 is large
	   mpz_mul(mpz_ptr, mat2->mpz_entries + ENTRY_TO_OFF(d2), x);
}

/*void _F_mpz_entry_mul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong r, const ulong c1, 
					                                 const F_mpz_mat_t mat2, const ulong r, const ulong c2)
{
	ulong c1 = mat1->entries[entry1];
   
	if (!ENTRY_IS_MPZ(c1)) // c1 is small
	{
		_F_mpz_entry_mul_si(res, entry3, mat2, entry2, c1);
      return;
	}
	
	__mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, entry3);
	ulong c2 = mat2->entries[entry2];
   	
	if (!ENTRY_IS_MPZ(c2)) // c1 is large, c2 is small
		mpz_mul_si(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(c1), c2);
   else // c1 and c2 are large
	   mpz_mul(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(c1), mat2->mpz_entries + ENTRY_TO_OFF(c2));
}*/

void _F_mpz_entry_mul_2exp(F_mpz_mat_t mat2, ulong r2, ulong c2, const F_mpz_mat_t mat1, 
									                  const ulong r1, const ulong c1, const ulong exp)
{
	long d = mat1->rows[r1][c1];

	if (!ENTRY_IS_MPZ(d)) // entry is small
	{
		ulong dabs = FLINT_ABS(d);
		ulong bits = FLINT_BIT_COUNT(dabs);
		if (bits + exp <= FLINT_BITS - 2) // result will fit in small
		{
			_F_mpz_entry_set_si(mat2, r2, c2, d<<exp);
		} else // result is large
		{
			__mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat2, r2, c2);
         mpz_set_si(mpz_ptr, d); 
	      mpz_mul_2exp(mpz_ptr, mpz_ptr, exp);
		}
	} else // entry is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat2, r2, c2);
      mpz_mul_2exp(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(d), exp);   
	}
}

void _F_mpz_entry_add_ui_inplace(F_mpz_mat_t res, ulong r, ulong c, const ulong x)
{
	ulong d = res->rows[r][c];

	if (!ENTRY_IS_MPZ(d)) // entry is small
	{
      mp_limb_t sum[2];
		if ((long) d >= 0L) // both operands non-negative
		{
			add_ssaaaa(sum[1], sum[0], 0, d, 0, x);
			if (sum[1] == 0) _F_mpz_entry_set_ui(res, r, c, sum[0]); // result fits in 1 limb
			else // result takes two limbs
			{
				mpz_t temp;
				temp->_mp_d = sum;
				temp->_mp_size = 2; // result is sum of two non-negative numbers and is hence non-negative
				__mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, r, c);
				mpz_set(mpz_ptr, temp);
			}
		} else // entry is negative, x positive
		{
			if (-d > x) _F_mpz_entry_set_si(res, r, c, x + d); // can't overflow as entry is small and x smaller
			else _F_mpz_entry_set_ui(res, r, c, x + d); // won't be negative and has to be less than x
		}
	} else
	{
		__mpz_struct * mpz_ptr = res->mpz_entries + ENTRY_TO_OFF(d);
		mpz_add_ui(mpz_ptr, mpz_ptr, x);
		_F_mpz_entry_demote_val(res, r, c); // cancellation may have occurred
	}
}

void _F_mpz_entry_sub_ui_inplace(F_mpz_mat_t res, ulong r, ulong c, const ulong x)
{
	ulong d = res->rows[r][c];

	if (!ENTRY_IS_MPZ(d)) // entry is small
	{
      mp_limb_t sum[2];
		if ((long) d < 0L) // entry negative, x positive, so difference is negative
		{
			add_ssaaaa(sum[1], sum[0], 0, -d, 0, x);
			if (sum[1] == 0) 
			{   
				_F_mpz_entry_set_ui(res, r, c, sum[0]); // result fits in 1 limb
				_F_mpz_entry_neg(res, r, c, res, r, c);
			} else // result takes two limbs
			{
				mpz_t temp;
				temp->_mp_d = sum;
				temp->_mp_size = -2; // result is negative number minus negative number, hence negative
				__mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, r, c);
				mpz_set(mpz_ptr, temp);
			}
		} else // entry is non-negative, x non-negative
		{
			if (x < d) _F_mpz_entry_set_ui(res, r, c, d - x); // won't be negative and is smaller than c
			else 
			{
				_F_mpz_entry_set_ui(res, r, c, x - d); // positive or zero
				_F_mpz_entry_neg(res, r, c, res, r, c);
			}
		}
	} else
	{
		__mpz_struct * mpz_ptr = res->mpz_entries + ENTRY_TO_OFF(d);
		mpz_sub_ui(mpz_ptr, mpz_ptr, x);
		_F_mpz_entry_demote_val(res, r, c); // cancellation may have occurred
	}
}

void _F_mpz_entry_addmul_ui(F_mpz_mat_t res, ulong r2, ulong c2, const F_mpz_mat_t mat1, 
							                                 const ulong r1, const ulong c1, const ulong x)
{
	ulong d1 = mat1->rows[r1][c1];
   if ((x == 0) || (d1 == 0)) return; // product is zero
   
	ulong r = res->rows[r2][c2];
	if (r == 0) 
	{
		_F_mpz_entry_mul_ui(res, r2, c2, mat1, r1, c1, x); // we are adding product to 0
		return;
	}

	if (!ENTRY_IS_MPZ(d1)) // c1 is small
	{
      mp_limb_t prod[2];
	   ulong ud1 = FLINT_ABS(d1);
      
		umul_ppmm(prod[1], prod[0], ud1, x); // compute product

		if (prod[1] == 0) // product fits in one limb
		{
			if ((long) d1 < 0L) _F_mpz_entry_sub_ui_inplace(res, r2, c2, prod[0]);
			else _F_mpz_entry_add_ui_inplace(res, r2, c2, prod[0]);
			return;
		} else if ((prod[1] == 1) && (!ENTRY_IS_MPZ(r)) && ((long)(r ^ d1) < 0L))
		{
			// only chance at cancellation is if product is one bit past a limb
			// and res is small and opposite sign to this product
			
			ulong ur = FLINT_ABS(r);
			if (ur > prod[0]) // cancellation will occur
			{
				_F_mpz_entry_set_ui(res, r2, c2, prod[0] - ur);
				if ((long) r > 0L) _F_mpz_entry_neg(res, r2, c2, res, r2, c2);
				return;
			} 
		}
		
		// in all remaining cases res is either big already, or will be big in the end
	   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, r2, c2);
		if (!ENTRY_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);
      
		mpz_t temp; // set up a temporary, cheap mpz_t to contain prod
	   temp->_mp_d = prod;
	   temp->_mp_size = ((long) d1 < 0L ? -2 : 2);
	   mpz_add(mpz_ptr, mpz_ptr, temp);
		_F_mpz_entry_demote_val(res, r2, c2); // cancellation may have occurred
	} else // c1 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, r2, c2);
		if (!ENTRY_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);

      mpz_addmul_ui(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(d1), x);
		_F_mpz_entry_demote_val(res, r2, c2); // cancellation may have occurred
	}
}

void _F_mpz_entry_submul_ui(F_mpz_mat_t res, ulong r2, ulong c2, const F_mpz_mat_t mat1, 
							                                 const ulong r1, const ulong c1, const ulong x)
{
	ulong d1 = mat1->rows[r1][c1];
   if ((x == 0) || (d1 == 0)) return; // product is zero
   
	ulong r = res->rows[r2][c2];
	if (r == 0) 
	{
		_F_mpz_entry_mul_ui(res, r2, c2, mat1, r1, c1, x); // we are adding product to 0
		_F_mpz_entry_neg(res, r2, c2, res, r2, c2);
		return;
	}

	mp_limb_t prod[2];
	
	if (!ENTRY_IS_MPZ(d1)) // c1 is small
	{
      ulong ud1 = FLINT_ABS(d1);
      
		umul_ppmm(prod[1], prod[0], ud1, x); // compute product

		if (prod[1] == 0) // product fits in one limb
		{
			if ((long) d1 < 0L) _F_mpz_entry_add_ui_inplace(res, r2, c2, prod[0]);
			else _F_mpz_entry_sub_ui_inplace(res, r2, c2, prod[0]);
			return;
		} else if ((prod[1] == 1) && (!ENTRY_IS_MPZ(r)) && ((long)(r ^ d1) >= 0L))
		{
			// only chance at cancellation is if product is one bit past a limb
			// and res is small and same sign as this product
			
			ulong ur = FLINT_ABS(r);
			if (ur > prod[0]) // cancellation will occur
			{
				_F_mpz_entry_set_ui(res, r2, c2, prod[0] - ur);
				if ((long) r > 0L) _F_mpz_entry_neg(res, r2, c2, res, r2, c2);
				return;
			} 
		}
		
		// in all remaining cases res is either big already, or will be big in the end
	   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, r2, c2);
		if (!ENTRY_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);

		mpz_t temp; // set up a temporary, cheap mpz_t to contain prod
	   temp->_mp_d = prod;
	   temp->_mp_size = ((long) d1 < 0L ? -2 : 2);
	   mpz_sub(mpz_ptr, mpz_ptr, temp);
		_F_mpz_entry_demote_val(res, r2, c2); // cancellation may have occurred
	} else // c1 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_entry_promote(res, r2, c2);
		if (!ENTRY_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);

      mpz_submul_ui(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(d1), x);
		_F_mpz_entry_demote_val(res, r2, c2); // cancellation may have occurred
	}
}

/*void _F_mpz_entry_addmul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong r, const ulong c1, 
						                                 const F_mpz_mat_t mat2, const ulong r, const ulong c2)
{
	ulong c1 = mat1->entries[entry1];
	
	if (!ENTRY_IS_MPZ(c1)) // c1 is small
	{
		if ((long) c1 < 0) _F_mpz_entry_submul_ui(res, entry3, mat2, entry2, -c1);
		else _F_mpz_entry_addmul_ui(res, entry3, mat2, entry2, c1);
		return;
	} 

	ulong c2 = mat2->entries[entry2];
   
	if (!ENTRY_IS_MPZ(c2)) // c2 is small
	{
		if ((long) c2 < 0) _F_mpz_entry_submul_ui(res, entry3, mat1, entry1, -c2);
		else _F_mpz_entry_addmul_ui(res, entry3, mat1, entry1, c2);
		return;
	} 

	// both c1 and c2 are large
   __mpz_struct * mpz_ptr = _F_mpz_entry_promote_val(res, entry3);
	
   mpz_addmul(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(c1), mat2->mpz_entries + ENTRY_TO_OFF(c2));
	_F_mpz_entry_demote_val(res, entry3); // cancellation may have occurred
}

void _F_mpz_entry_submul(F_mpz_mat_t res, ulong entry3, const F_mpz_mat_t mat1, const ulong r, const ulong c1, 
						                                 const F_mpz_mat_t mat2, const ulong r, const ulong c2)
{
	ulong c1 = mat1->entries[entry1];
	
	if (!ENTRY_IS_MPZ(c1)) // c1 is small
	{
		if ((long) c1 < 0) _F_mpz_entry_addmul_ui(res, entry3, mat2, entry2, -c1);
		else _F_mpz_entry_submul_ui(res, entry3, mat2, entry2, c1);
		return;
	} 

	ulong c2 = mat2->entries[entry2];
   
	if (!ENTRY_IS_MPZ(c2)) // c2 is small
	{
		if ((long) c2 < 0) _F_mpz_entry_addmul_ui(res, entry3, mat1, entry1, -c2);
		else _F_mpz_entry_submul_ui(res, entry3, mat1, entry1, c2);
		return;
	} 

	// both c1 and c2 are large
  __mpz_struct * mpz_ptr = _F_mpz_entry_promote_val(res, entry3);
	
   mpz_submul(mpz_ptr, mat1->mpz_entries + ENTRY_TO_OFF(c1), mat2->mpz_entries + ENTRY_TO_OFF(c2));
	_F_mpz_entry_demote_val(res, entry3); // cancellation may have occurred
}*/

ulong _F_mpz_entry_size(F_mpz_mat_t mat, const ulong r, const ulong c)
{
	ulong d = mat->rows[r][c];

	if (d == 0) return 0;
   if (!ENTRY_IS_MPZ(d)) // c1 is small
	{
		return 1;
	}

	return mpz_size(mat->mpz_entries + ENTRY_TO_OFF(d));
}

ulong _F_mpz_entry_print(F_mpz_mat_t mat, const ulong r, const ulong c)
{
   ulong d = mat->rows[r][c];

	if (!ENTRY_IS_MPZ(d)) 
	   printf("%ld", d);
   else
	   gmp_printf("%Zd", mat->mpz_entries + ENTRY_TO_OFF(d));
}

/* ==============================================================================

   Random generation

===============================================================================*/

/*
  These are not serious random generators, they are just here for testing 
  purposes at this stage

  We require bits to be non-zero
*/

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

void F_mpz_entry_random(F_mpz_mat_t mat, const ulong r, const ulong c, const ulong bits)
{
   mp_limb_t d = mat->rows[r][c];
   
	if (bits <= FLINT_BITS - 2)
   {
      ulong temp;
      mpn_random(&temp, 1L);
      ulong mask = ((1L<<bits)-1L);
      mat->rows[r][c] = temp & mask;
      return;
   }
   
	ulong limbs = ((bits-1)>>FLINT_LG_BITS_PER_LIMB)+1;
   ulong rem = (bits & (FLINT_BITS - 1));
   
   __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat, r, c);
   mpz_realloc2(mpz_ptr, bits);
	
	mp_limb_t * fp = mpz_ptr->_mp_d;
   mpz_ptr->_mp_size = limbs;
   mpn_random(fp, limbs);
   if (rem)
   {
      ulong mask = ((1L<<rem)-1L);
      fp[limbs-1] &= mask;
   }
   NORM(mpz_ptr, fp);
	_F_mpz_entry_demote_val(mat, r, c);
}

void F_mpz_entry_randomm(F_mpz_mat_t mat, const ulong r, const ulong c, const mpz_t in)
{
   if (mpz_size(in) > 1) 
   {
      __mpz_struct * mpz_ptr = _F_mpz_entry_promote(mat, r, c);
		mpz_urandomm(mpz_ptr, state, in);
		_F_mpz_entry_demote_val(mat, r, c);
   } else
   {
      ulong val = mpz_get_ui(in);
		ulong rnd = (val == 0 ? 0L : z_randint(val));
		_F_mpz_entry_set_ui(mat, r, c, rnd);
   }
}

/* ==============================================================================

   Input/Output

===============================================================================*/

void F_mpz_entry_read(F_mpz_mat_t mat, const ulong r, const ulong c)
{
	mpz_t temp;
	mpz_init(temp);

	mpz_inp_str(temp, stdin, 10);
	_F_mpz_entry_set_mpz(mat, r, c, temp);

	mpz_clear(temp);
}

void F_mpz_mat_print(F_mpz_mat_t mat) 
{
   ulong i, j; 
   ulong r = mat->r;
   ulong c = mat->c;
	
   printf("[");
   for (i = 0; i < r; i++) 
   {
      printf("[");
      for (j = 0; j < c; j++) 
	   { 
	      _F_mpz_entry_print(mat, i, j); 
	      if (j < c - 1) printf(" "); 
	   }
      if (i != r - 1) printf("]\n"); 
   }  
   printf("]]\n"); 
}

/*===============================================================================

	Conversions

================================================================================*/

void mpz_mat_to_F_mpz_mat(F_mpz_mat_t F_mat, const mpz_mat_t m_mat)
{
	for (ulong r = 0; r < m_mat->r; r++)
	{
		ulong row = r*m_mat->c;
		for (ulong c = 0; c < m_mat->c; c++)
		_F_mpz_entry_set_mpz(F_mat, r, c, m_mat->entries[row+c]);
	}
}

void F_mpz_mat_to_mpz_mat(mpz_mat_t m_mat, const F_mpz_mat_t F_mat)
{
	for (ulong r = 0; r < m_mat->r; r++)
	{
		ulong row = r*m_mat->c;
		for (ulong c = 0; c < m_mat->c; c++)
	   _F_mpz_entry_get_mpz(m_mat->entries[row+c], F_mat, r, c);
	}
}

int F_mpz_mat_set_line_d(double * appv, const F_mpz_mat_t mat, const ulong r, const int n)
{
   long * exp, i, maxexp = 0L;
   exp = (long *) malloc(n * sizeof(long)); 
  
   for (i = 0; i < n; i++)
   {
      appv[i] = _F_mpz_entry_get_d_2exp(&exp[i], mat, r, i);
      if (exp[i] > maxexp) maxexp = exp[i];
   }

   for (i = 0; i < n; i++) appv[i] = ldexp(appv[i], exp[i] - maxexp);

   free(exp);
   return maxexp;
}

/*===============================================================================

	Assignment/swap

================================================================================*/

void F_mpz_mat_set(F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
	ulong r = mat1->r;
	ulong c = mat1->c;
	
	if (mat1 != mat2)
	{
		for (ulong i = 0; i < r; i++)
			for (ulong j = 0; j < c; j++)
			   _F_mpz_entry_set(mat1, i, j, mat2, i, j);
	}
}

/*void F_mpz_mat_swap(F_mpz_mat_t mat1, F_mpz_mat_t mat2)
{
	if (mat1 == mat2) return;

	ulong temp = mat1->length;
	mat1->length = mat2->length;
	mat2->length = temp;
	
	temp = mat1->alloc;
	mat1->alloc = mat2->alloc;
	mat2->alloc = temp;
	
	temp = mat1->mpz_alloc;
	mat1->mpz_alloc = mat2->mpz_alloc;
	mat2->mpz_alloc = temp;
	
	temp = mat1->mpz_length;
	mat1->mpz_length = mat2->mpz_length;
	mat2->mpz_length = temp;

	ulong * temp_c = mat1->entries;
	mat1->entries = mat2->entries;
	mat2->entries = temp_c;

   __mpz_struct * temp_m = mat1->mpz_entries;
	mat1->mpz_entries = mat2->mpz_entries;
	mat2->mpz_entries = temp_m;

   return;
}*/

/*===============================================================================

	Comparison

================================================================================*/

int F_mpz_mat_equal(const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   if (mat1 == mat2) return 1; // same matrix
	
	ulong r = mat1->r;
	ulong c = mat1->c;

	for (ulong i = 0; i < r; i++) // check if entries the same
		for (ulong j = 0; j < c; j++)
			if (!_F_mpz_entry_equal(mat1, i, j, mat2, i, j)) return 0;

	return 1;
}

/*===============================================================================

	Coefficient sizes

================================================================================*/

/*long F_mpz_mat_max_bits(const F_mpz_mat_t mat)
{
	int sign = 0;
	ulong max = 0;
   ulong bits = 0;
   ulong i;
	ulong c;

	// search until we find an mpz_t entry or one of at least FLINT_BITS - 2 bits
	for (i = 0; i < mat->length; i++) 
	{
		c = mat->entries[i];
		if (ENTRY_IS_MPZ(c)) break; // found an mpz_t entry
      if ((long) c < 0L) 
		{
			sign = 1;
         bits = FLINT_BIT_COUNT(-c);
		} else bits = FLINT_BIT_COUNT(c);
		if (bits > max) max = bits;
		if (max == FLINT_BITS - 2) break; // entry is at least FLINT_BITS - 2 bits
	}

   // search through mpz entries for largest size in bits
	for ( ; i < mat->length; i++)
	{
		c = mat->entries[i];
      if (ENTRY_IS_MPZ(c))
		{
			__mpz_struct * mpz_ptr = mat->mpz_entries + ENTRY_TO_OFF(c);
			if (mpz_sgn(mpz_ptr) < 0) sign = 1;
			bits = mpz_sizeinbase(mpz_ptr, 2);
			if (bits > max) max = bits;
		} else if ((long) c < 0L) sign = 1; // still need to check the sign of small entries
	}

	if (sign) return -max;
	else return max;
}

ulong F_mpz_mat_max_limbs(const F_mpz_mat_t mat)
{
	if (mat->length == 0) return 0; // matrix is zero

	ulong max = 1; // all other entries have at least one limb
   ulong limbs;
	ulong c;

   // search through mpz entries for one of largest size
	for (ulong i = 0; i < mat->length; i++)
	{
		c = mat->entries[i];
      if (ENTRY_IS_MPZ(c))
		{
			limbs = mpz_size(mat->mpz_entries + ENTRY_TO_OFF(c));
			if (limbs > max) max = limbs;
		} 
	}

	return max;
}*/


/*===============================================================================

	Negation

================================================================================*/

void F_mpz_mat_neg(F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
	ulong r = mat1->r;
	ulong c = mat1->c;
	
	for (ulong i = 0; i < r; i++)
		for (ulong j = 0; j < c; j++)
		   _F_mpz_entry_neg(mat1, i, j, mat2, i, j);
}
/*===============================================================================

	Addition/subtraction

================================================================================*/

void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   ulong r = mat1->r;
	ulong c = mat1->c;
		
	for (ulong i = 0; i < r; i++) // add up to the length of the shorter mat
		for (ulong j = 0; j < c; j++)
			_F_mpz_entry_add(res, i, j, mat1, i, j, mat2, i, j);   
}

void F_mpz_mat_row_add(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, 
							  const ulong r1, const F_mpz_mat_t mat2, const ulong r2, 
							                         const ulong start, const ulong n)
{
   for (ulong i = start; i < start + n; i++) 
		_F_mpz_entry_add(res, r3, i, mat1, r1, i, mat2, r2, i);   
}

void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   ulong r = mat1->r;
	ulong c = mat1->c;
		
	for (ulong i = 0; i < r; i++) // add up to the length of the shorter mat
		for (ulong j = 0; j < c; j++)
			_F_mpz_entry_sub(res, i, j, mat1, i, j, mat2, i, j);   
}

void F_mpz_mat_row_sub(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, 
							  const ulong r1, const F_mpz_mat_t mat2, const ulong r2, 
							                          const ulong start, const ulong n)
{
   for (ulong i = start; i < start + n; i++) 
		_F_mpz_entry_sub(res, r3, i, mat1, r1, i, mat2, r2, i);   
}


/*===============================================================================

	Scalar multiplication

================================================================================*/

void F_mpz_mat_row_mul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)
{
	// either scalar of input mat is zero
	if (x == 0L)
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_zero(mat1, r1, i);

		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_set(mat1, r1, i, mat2, r2, i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		_F_mpz_entry_mul_ui(mat1, r1, i, mat2, r2, i, x);
}

void F_mpz_mat_row_mul_si(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, long x)
{
	// either scalar of input mat is zero
	if (x == 0L)
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_zero(mat1, r1, i);

		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_set(mat1, r1, i, mat2, r2, i);

		return;
	}
	
	if (x == -1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_neg(mat1, r1, i, mat2, r2, i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		_F_mpz_entry_mul_si(mat1, r1, i, mat2, r2, i, x);
}

void F_mpz_mat_row_mul_mpz(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, mpz_t x)
{
	// either scalar or input mat is zero
	if (mpz_cmpabs_ui(x, 0L) == 0)
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_zero(mat1, r1, i);

		return;
	}
	
	// special cases, muliply by +/- 1
	if (mpz_cmpabs_ui(x, 1L) == 0)
	{
	   if (mpz_sgn(x) < 0) 
		   for (ulong i = start; i < start + n; i++)
		      _F_mpz_entry_neg(mat1, r1, i, mat2, r2, i);
		else 
			for (ulong i = start; i < start + n; i++)
		      _F_mpz_entry_set(mat1, r1, i, mat2, r2, i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		_F_mpz_entry_mul_mpz(mat1, r1, i, mat2, r2, i, x);
}

/*===============================================================================

	Scalar addmul/submul

================================================================================*/

void F_mpz_mat_row_addmul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)
{
	// scalar is zero
	if (x == 0L)
		return;
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_add(mat1, r1, i, mat1, r1, i, mat2, r2, i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		_F_mpz_entry_addmul_ui(mat1, r1, i, mat2, r2, i, x);
}

void F_mpz_mat_row_submul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)
{
	// scalar is zero
	if (x == 0L)
		return;
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_sub(mat1, r1, i, mat1, r1, i, mat2, r2, i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		_F_mpz_entry_submul_ui(mat1, r1, i, mat2, r2, i, x);
}

void F_mpz_mat_row_addmul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp, F_mpz_mat_t temp)
{
	// scalar is zero, nothing to add
	if (c == 0)
	{
	   return;
	}
	
	// scalar is 1, just add 2^exp times the entry
	if (c == 1)
	{
	   for (ulong i = start; i < start + n; i++)
		{
			_F_mpz_entry_mul_2exp(temp, 0, i, mat2, r2, i, exp);
         _F_mpz_entry_add(mat1, r1, i, mat1, r1, i, temp, 0, i);
		}

		return;
	}
	
	// exp is 1, just do addmul
	if (exp == 0)
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_addmul_ui(mat1, r1, i, mat2, r2, i, c);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
	{
		_F_mpz_entry_mul_2exp(temp, 0, i, mat2, r2, i, exp);
	   _F_mpz_entry_addmul_ui(mat1, r1, i, temp, 0, i, c);
	}
}

void F_mpz_mat_row_submul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp, F_mpz_mat_t temp)
{
	// scalar is zero, nothing to add
	if (c == 0)
	{
	   return;
	}
	
	// scalar is 1, just add 2^exp times the entry
	if (c == 1)
	{
	   for (ulong i = start; i < start + n; i++)
		{
			_F_mpz_entry_mul_2exp(temp, 0, i, mat2, r2, i, exp);
         _F_mpz_entry_sub(mat1, r1, i, mat1, r1, i, temp, 0, i);
		}

		return;
	}
	
	// exp is 1, just do addmul
	if (exp == 0)
	{
	   for (ulong i = start; i < start + n; i++)
		   _F_mpz_entry_submul_ui(mat1, r1, i, mat2, r2, i, c);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
	{
		_F_mpz_entry_mul_2exp(temp, 0, i, mat2, r2, i, exp);
	   _F_mpz_entry_submul_ui(mat1, r1, i, temp, 0, i, c);
	}
}

