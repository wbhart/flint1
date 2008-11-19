/*============================================================================

    F_mpz.c: The FLINT integer format (FLINT 2.0)

    Copyright (C) 2008, William Hart 

	 This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/

#include <stdint.h>
#include <string.h>
#include <math.h>

#include "flint.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz.h"

/*===============================================================================

	mpz_t memory management

================================================================================*/

// The array of mpz's used by the F_mpz type
__mpz_struct * F_mpz_arr;

// Total number of mpz's initialised in F_mpz_arr;
ulong F_mpz_allocated;

// An array of indices of mpz's which are not being used presently. These are stored with the second
// most significant bit set so that they are a valid index as per the F_mpz_t type below.
F_mpz * F_mpz_unused_arr;

// The number of mpz's not being used presently
ulong F_mpz_num_unused = 0;

F_mpz _F_mpz_new_mpz(void)
{
	if (!F_mpz_num_unused) // time to allocate MPZ_BLOCK more mpz_t's
	{
	   if (F_mpz_allocated) // realloc mpz_t's and unused array
		{
			F_mpz_arr = (__mpz_struct *) flint_heap_realloc_bytes(F_mpz_arr, (F_mpz_allocated + MPZ_BLOCK)*sizeof(__mpz_struct));
			F_mpz_unused_arr = (F_mpz *) flint_heap_realloc_bytes(F_mpz_unused_arr, (F_mpz_allocated + MPZ_BLOCK)*sizeof(F_mpz));
		} else // first time alloc of mpz_t's and unused array
		{
			F_mpz_arr = (__mpz_struct*) flint_heap_alloc_bytes(MPZ_BLOCK*sizeof(__mpz_struct));	
			F_mpz_unused_arr = (F_mpz *) flint_heap_alloc_bytes(MPZ_BLOCK*sizeof(F_mpz));
		}
		
		// initialise the new mpz_t's and unused array
		for (ulong i = 0; i < MPZ_BLOCK; i++)
		{
			mpz_init(F_mpz_arr + F_mpz_allocated + i);
			F_mpz_unused_arr[F_mpz_num_unused] = OFF_TO_COEFF(F_mpz_allocated + i);
         F_mpz_num_unused++;
		}

		F_mpz_num_unused--;
		F_mpz_allocated += MPZ_BLOCK;

		return F_mpz_unused_arr[F_mpz_num_unused];
	} else // unused mpz's are available
	{
		F_mpz_num_unused--;
		return F_mpz_unused_arr[F_mpz_num_unused];
	}
}

void _F_mpz_clear_mpz(F_mpz f)
{
   F_mpz_unused_arr[F_mpz_num_unused] = f;
   F_mpz_num_unused++;	
}

void _F_mpz_cleanup(void)
{
	for (ulong i = 0; i < F_mpz_num_unused; i++)
	{
		mpz_clear(F_mpz_arr + COEFF_TO_OFF(F_mpz_unused_arr[i]));
	}
	
	if (F_mpz_num_unused) free(F_mpz_unused_arr);
	if (F_mpz_allocated) free(F_mpz_arr);
}

/*===============================================================================

	Promotion/Demotion

================================================================================*/

__mpz_struct * _F_mpz_promote(F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) *f = _F_mpz_new_mpz(); // f is small so promote it first
	// if f is large already, just return the pointer
      
   return F_mpz_arr + COEFF_TO_OFF(*f);
}

__mpz_struct * _F_mpz_promote_val(F_mpz_t f)
{
   F_mpz c = *f;
	if (!COEFF_IS_MPZ(c)) // f is small so promote it
	{
	   *f = _F_mpz_new_mpz();
	   __mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(*f);
		mpz_set_si(mpz_ptr, c);
		return mpz_ptr;
	} else // f is large already, just return the pointer
      return F_mpz_arr + COEFF_TO_OFF(*f);
}

void _F_mpz_demote_val(F_mpz_t f)
{
   __mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(*f);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L) // value is zero
	{
		*f = 0;
	} else if (size == 1L) // value is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= COEFF_MAX) *f = (F_mpz) uval;
	} else if (size == -1L) // value is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= COEFF_MAX) *f = (F_mpz) -uval;
	}
	// don't do anything if value has to be multi precision
}

/*===============================================================================

	F_mpz_t memory management

================================================================================*/

void F_mpz_init2(F_mpz_t f, ulong limbs)
{
	if (limbs) 
	{
		*f = _F_mpz_new_mpz();
		_mpz_realloc(F_mpz_arr + COEFF_TO_OFF(*f), limbs);
		return;
	} else 
	{
		*f = 0;
		return;
	}
}

/*===============================================================================

	Get/set

================================================================================*/

void F_mpz_set_si(F_mpz_t f, const long val)
{
   if (FLINT_ABS(val) > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(f);
		mpz_set_si(mpz_coeff, val);
	} else 
	{
		_F_mpz_demote(f);
		*f = val; // val is small
	}
}

void F_mpz_set_ui(F_mpz_t f, const ulong val)
{
   if (val > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(f);
		mpz_set_ui(mpz_coeff, val);
	} else 
	{
		_F_mpz_demote(f);
		*f = val; // val is small
	}
}

long F_mpz_get_si(const F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) return *f; // value is small
	return mpz_get_si(F_mpz_arr + COEFF_TO_OFF(*f)); // value is large
}

long F_mpz_get_ui(const F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) return *f; // value is small
	return mpz_get_ui(F_mpz_arr + COEFF_TO_OFF(*f)); // value is large
}

void F_mpz_get_mpz(mpz_t x, const F_mpz_t f)
{
	if (!COEFF_IS_MPZ(*f)) mpz_set_si(x, *f); // set x to small value
	else mpz_set(x, F_mpz_arr + COEFF_TO_OFF(*f)); // set x to large value
}

void F_mpz_set_mpz(F_mpz_t f, const mpz_t x)
{
   long size = x->_mp_size;
	
	if (size == 0L) // x is zero
	{
		F_mpz_zero(f);
	} else if (size == 1L) // x is positive and 1 limb
	{
	   F_mpz_set_ui(f, mpz_get_ui(x));
	} else if (size == -1L) // x is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(x);
		if (uval <= COEFF_MAX) // x is small
		{
		   _F_mpz_demote(f);
		   *f = -uval;
		} else // x is large but one limb
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			mpz_set_ui(mpz_ptr, uval);
			mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // x is more than one limb
	{
		__mpz_struct * mpz_ptr = _F_mpz_promote(f);
		mpz_set(mpz_ptr, x);
	}			
}

void F_mpz_set(F_mpz_t f, const F_mpz_t g)
{
	if (f == g) return; // aliased inputs
	
	if (!COEFF_IS_MPZ(*g)) // g is small
	{
		_F_mpz_demote(f);
		*f = *g;
	} else // g is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_promote(f);
		mpz_set(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*g));
	}
}

void F_mpz_swap(F_mpz_t f, F_mpz_t g)
{
	if (f == g) return; // swapping not required
	
	if (!COEFF_IS_MPZ(*f))
	{
		if (!COEFF_IS_MPZ(*g)) // both values are small
		{
	      F_mpz t = *f;
			*f = *g;
	      *g = t;
		} else // f is small, g is large
		{
			F_mpz t = *f;
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			mpz_set(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*g));
			_F_mpz_demote(g);
         *g = t;
		}
	} else
	{
      if (!COEFF_IS_MPZ(*g)) // g is small, f is large
		{
         F_mpz t = *g;
			__mpz_struct * mpz_ptr = _F_mpz_promote(g);
			mpz_set(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*f));
			_F_mpz_demote(f);
         *f = t;
		} else // both values are large
		{
			mpz_swap(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g));
		}
	}
}

/*===============================================================================

	Comparison

================================================================================*/

int F_mpz_equal(const F_mpz_t f, const F_mpz_t g)
{
	if (f == g) return 1; // aliased inputs
	
	if (!COEFF_IS_MPZ(*g)) return (*f == *g); // if f is large it can't be equal to g
	else if (!COEFF_IS_MPZ(*g)) return 0; // c1 is large, so if c2 isn't....
	else return (mpz_cmp(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)) == 0); 
}

/*===============================================================================

	Arithmetic

================================================================================*/

void F_mpz_neg(F_mpz_t f1, const F_mpz_t f2)
{
   if (!COEFF_IS_MPZ(*f2)) // coeff is small
	{
		F_mpz t = -*f2; // Need to save value in case of aliasing
		_F_mpz_demote(f1);
		*f1 = t;
	} else // coeff is large
	{
	   // No need to retain value in promotion, as if aliased, both already large
		__mpz_struct * mpz_ptr = _F_mpz_promote(f1);
		mpz_neg(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*f2));
	}
}

void F_mpz_add(F_mpz_t f, const F_mpz_t g, F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // both inputs are small
		{
			F_mpz_set_si(f, c1 + c2);
		} else // g is small, h is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // g is saved and h is large
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			if (c1 < 0L) mpz_sub_ui(mpz3, mpz2, -c1);	
		   else mpz_add_ui(mpz3, mpz2, c1);
			_F_mpz_demote_val(f); // may have cancelled
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // g is large, h is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // h is saved and g is large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			if (c2 < 0L) mpz_sub_ui(mpz3, mpz1, -c2);	
			else mpz_add_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(f); // may have cancelled
		} else // g and h are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // aliasing means f is already large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			mpz_add(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(f); // may have cancelled
		}
	}
}

void F_mpz_sub(F_mpz_t f, const F_mpz_t g, F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // both inputs are small
		{
			F_mpz_set_si(f, c1 - c2);
		} else // g is small, h is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // g is saved and h is large
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			if (c1 < 0L) 
			{
				mpz_add_ui(mpz3, mpz2, -c1);
				mpz_neg(mpz3, mpz3);
			} else mpz_ui_sub(mpz3, c1, mpz2);
			_F_mpz_demote_val(f); // may have cancelled
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // g is large, h is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // h is saved and g is large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			if (c2 < 0L) mpz_add_ui(mpz3, mpz1, -c2);	
			else mpz_sub_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(f); // may have cancelled
		} else // g and h are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // aliasing means f is already large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			mpz_sub(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(f); // may have cancelled
		}
	}
}

/*
void _F_mpz_mul_ui(F_mpz_poly_t poly1, ulong coeff1, const F_mpz_poly_t poly2, 
						                                   const ulong coeff2, const ulong x)
{
	ulong c2 = poly2->coeffs[coeff2];

	if (!COEFF_IS_MPZ(c2)) // coeff2 is small
	{
		mp_limb_t prod[2];
		ulong uc2 = FLINT_ABS(c2);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], uc2, x); 
		if (!prod[1]) // result fits in one limb
		{
			_F_mpz_set_ui(poly1, coeff1, prod[0]);
			if ((long) c2 < 0L) _F_mpz_neg(poly1, coeff1, poly1, coeff1);
		} else // result takes two limbs
		{
		   __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
			// two limbs, least significant first, native endian, no nails, stored in prod
         mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if ((long) c2 < 0L) mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // coeff2 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
      mpz_mul_ui(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c2), x);
	}
}

void _F_mpz_mul_si(F_mpz_poly_t poly1, ulong coeff1, const F_mpz_poly_t poly2, 
						                                   const ulong coeff2, const long x)
{
	ulong c2 = poly2->coeffs[coeff2];

	if (!COEFF_IS_MPZ(c2)) // coeff2 is small
	{
		mp_limb_t prod[2];
		ulong uc2 = FLINT_ABS(c2);
		ulong ux = FLINT_ABS(x);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], uc2, ux); 
		if (!prod[1]) // result fits in one limb
		{
			_F_mpz_set_ui(poly1, coeff1, prod[0]);
			if ((long) (c2 ^ x) < 0L) _F_mpz_neg(poly1, coeff1, poly1, coeff1);
		} else // result takes two limbs
		{
		   __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
         // two limbs, least significant first, native endian, no nails, stored in prod
			mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if ((long) (c2 ^ x) < 0L) mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // coeff2 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
      mpz_mul_si(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c2), x);
	}
}

void _F_mpz_mul_mpz(F_mpz_poly_t poly1, ulong coeff1, F_mpz_poly_t poly2, ulong coeff2, mpz_t x)
{
	ulong c2 = poly2->coeffs[coeff2];
   
	if (mpz_size(x) <= 1) // 1 limb to multiply by
	{
		long x_limb = mpz_get_ui(x);
      _F_mpz_mul_ui(poly1, coeff1, poly2, coeff2, x_limb); 
		if (mpz_sgn(x) < 0) _F_mpz_neg(poly1, coeff1, poly1, coeff1);
		return;
	}

   // more than one limb to multiply by
	__mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);

	if (!COEFF_IS_MPZ(c2)) // coeff2 is small
	   mpz_mul_si(mpz_ptr, x, c2);
	else // coeff2 is large
	   mpz_mul(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c2), x);
}

void _F_mpz_mul(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                 const F_mpz_poly_t poly2, const ulong coeff2)
{
	ulong c1 = poly1->coeffs[coeff1];
   
	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
		_F_mpz_mul_si(res, coeff3, poly2, coeff2, c1);
      return;
	}
	
	__mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff3);
	ulong c2 = poly2->coeffs[coeff2];
   	
	if (!COEFF_IS_MPZ(c2)) // c1 is large, c2 is small
		mpz_mul_si(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1), c2);
   else // c1 and c2 are large
	   mpz_mul(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1), poly2->mpz_coeffs + COEFF_TO_OFF(c2));
}

void _F_mpz_add_ui_inplace(F_mpz_poly_t res, ulong coeff, const ulong x)
{
	ulong c = res->coeffs[coeff];

	if (!COEFF_IS_MPZ(c)) // coeff is small
	{
      mp_limb_t sum[2];
		if ((long) c >= 0L) // both operands non-negative
		{
			add_ssaaaa(sum[1], sum[0], 0, c, 0, x);
			if (sum[1] == 0) _F_mpz_set_ui(res, coeff, sum[0]); // result fits in 1 limb
			else // result takes two limbs
			{
				mpz_t temp;
				temp->_mp_d = sum;
				temp->_mp_size = 2; // result is sum of two non-negative numbers and is hence non-negative
				__mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff);
				mpz_set(mpz_ptr, temp);
			}
		} else // coeff is negative, x positive
		{
			if (-c > x) _F_mpz_set_si(res, coeff, x + c); // can't overflow as coeff is small and x smaller
			else _F_mpz_set_ui(res, coeff, x + c); // won't be negative and has to be less than x
		}
	} else
	{
		__mpz_struct * mpz_ptr = res->mpz_coeffs + COEFF_TO_OFF(c);
		mpz_add_ui(mpz_ptr, mpz_ptr, x);
		_F_mpz_demote_val(res, coeff); // cancellation may have occurred
	}
}

void _F_mpz_sub_ui_inplace(F_mpz_poly_t res, ulong coeff, const ulong x)
{
	ulong c = res->coeffs[coeff];

	if (!COEFF_IS_MPZ(c)) // coeff is small
	{
      mp_limb_t sum[2];
		if ((long) c < 0L) // coeff negative, x positive, so difference is negative
		{
			add_ssaaaa(sum[1], sum[0], 0, -c, 0, x);
			if (sum[1] == 0) 
			{   
				_F_mpz_set_ui(res, coeff, sum[0]); // result fits in 1 limb
				_F_mpz_neg(res, coeff, res, coeff);
			} else // result takes two limbs
			{
				mpz_t temp;
				temp->_mp_d = sum;
				temp->_mp_size = -2; // result is negative number minus negative number, hence negative
				__mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff);
				mpz_set(mpz_ptr, temp);
			}
		} else // coeff is non-negative, x non-negative
		{
			if (x < c) _F_mpz_set_ui(res, coeff, c - x); // won't be negative and is smaller than c
			else 
			{
				_F_mpz_set_ui(res, coeff, x - c); // positive or zero
				_F_mpz_neg(res, coeff, res, coeff);
			}
		}
	} else
	{
		__mpz_struct * mpz_ptr = res->mpz_coeffs + COEFF_TO_OFF(c);
		mpz_sub_ui(mpz_ptr, mpz_ptr, x);
		_F_mpz_demote_val(res, coeff); // cancellation may have occurred
	}
}

void _F_mpz_addmul_ui(F_mpz_poly_t res, ulong coeff2, const F_mpz_poly_t poly1, 
							                                 const ulong coeff1, const ulong x)
{
	ulong c1 = poly1->coeffs[coeff1];
   if ((x == 0) || (c1 == 0)) return; // product is zero
   
	ulong r = res->coeffs[coeff2];
	if (r == 0) 
	{
		_F_mpz_mul_ui(res, coeff2, poly1, coeff1, x); // we are adding product to 0
		return;
	}

	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
      mp_limb_t prod[2];
	   ulong uc1 = FLINT_ABS(c1);
      
		umul_ppmm(prod[1], prod[0], uc1, x); // compute product

		if (prod[1] == 0) // product fits in one limb
		{
			if ((long) c1 < 0L) _F_mpz_sub_ui_inplace(res, coeff2, prod[0]);
			else _F_mpz_add_ui_inplace(res, coeff2, prod[0]);
			return;
		} else if ((prod[1] == 1) && (!COEFF_IS_MPZ(r)) && ((long)(r ^ c1) < 0L))
		{
			// only chance at cancellation is if product is one bit past a limb
			// and res is small and opposite sign to this product
			
			ulong ur = FLINT_ABS(r);
			if (ur > prod[0]) // cancellation will occur
			{
				_F_mpz_set_ui(res, coeff2, prod[0] - ur);
				if ((long) r > 0L) _F_mpz_neg(res, coeff2, res, coeff2);
				return;
			} 
		}
		
		// in all remaining cases res is either big already, or will be big in the end
	   __mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff2);
		if (!COEFF_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);
      
		mpz_t temp; // set up a temporary, cheap mpz_t to contain prod
	   temp->_mp_d = prod;
	   temp->_mp_size = ((long) c1 < 0L ? -2 : 2);
	   mpz_add(mpz_ptr, mpz_ptr, temp);
		_F_mpz_demote_val(res, coeff2); // cancellation may have occurred
	} else // c1 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff2);
		if (!COEFF_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);

      mpz_addmul_ui(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1), x);
		_F_mpz_demote_val(res, coeff2); // cancellation may have occurred
	}
}

void _F_mpz_submul_ui(F_mpz_poly_t res, ulong coeff2, const F_mpz_poly_t poly1, 
							                                 const ulong coeff1, const ulong x)
{
	ulong c1 = poly1->coeffs[coeff1];
   if ((x == 0) || (c1 == 0)) return; // product is zero
   
	ulong r = res->coeffs[coeff2];
	if (r == 0) 
	{
		_F_mpz_mul_ui(res, coeff2, poly1, coeff1, x); // we are adding product to 0
		_F_mpz_neg(res, coeff2, res, coeff2);
		return;
	}

	mp_limb_t prod[2];
	
	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
      ulong uc1 = FLINT_ABS(c1);
      
		umul_ppmm(prod[1], prod[0], uc1, x); // compute product

		if (prod[1] == 0) // product fits in one limb
		{
			if ((long) c1 < 0L) _F_mpz_add_ui_inplace(res, coeff2, prod[0]);
			else _F_mpz_sub_ui_inplace(res, coeff2, prod[0]);
			return;
		} else if ((prod[1] == 1) && (!COEFF_IS_MPZ(r)) && ((long)(r ^ c1) >= 0L))
		{
			// only chance at cancellation is if product is one bit past a limb
			// and res is small and same sign as this product
			
			ulong ur = FLINT_ABS(r);
			if (ur > prod[0]) // cancellation will occur
			{
				_F_mpz_set_ui(res, coeff2, prod[0] - ur);
				if ((long) r > 0L) _F_mpz_neg(res, coeff2, res, coeff2);
				return;
			} 
		}
		
		// in all remaining cases res is either big already, or will be big in the end
	   __mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff2);
		if (!COEFF_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);

		mpz_t temp; // set up a temporary, cheap mpz_t to contain prod
	   temp->_mp_d = prod;
	   temp->_mp_size = ((long) c1 < 0L ? -2 : 2);
	   mpz_sub(mpz_ptr, mpz_ptr, temp);
		_F_mpz_demote_val(res, coeff2); // cancellation may have occurred
	} else // c1 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff2);
		if (!COEFF_IS_MPZ(r)) // res is small
		   mpz_set_si(mpz_ptr, r);

      mpz_submul_ui(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1), x);
		_F_mpz_demote_val(res, coeff2); // cancellation may have occurred
	}
}

void _F_mpz_addmul(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
						                                 const F_mpz_poly_t poly2, const ulong coeff2)
{
	ulong c1 = poly1->coeffs[coeff1];
	
	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
		if ((long) c1 < 0) _F_mpz_submul_ui(res, coeff3, poly2, coeff2, -c1);
		else _F_mpz_addmul_ui(res, coeff3, poly2, coeff2, c1);
		return;
	} 

	ulong c2 = poly2->coeffs[coeff2];
   
	if (!COEFF_IS_MPZ(c2)) // c2 is small
	{
		if ((long) c2 < 0) _F_mpz_submul_ui(res, coeff3, poly1, coeff1, -c2);
		else _F_mpz_addmul_ui(res, coeff3, poly1, coeff1, c2);
		return;
	} 

	// both c1 and c2 are large
   __mpz_struct * mpz_ptr = _F_mpz_promote_val(res, coeff3);
	
   mpz_addmul(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1), poly2->mpz_coeffs + COEFF_TO_OFF(c2));
	_F_mpz_demote_val(res, coeff3); // cancellation may have occurred
}

void _F_mpz_submul(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
						                                 const F_mpz_poly_t poly2, const ulong coeff2)
{
	ulong c1 = poly1->coeffs[coeff1];
	
	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
		if ((long) c1 < 0) _F_mpz_addmul_ui(res, coeff3, poly2, coeff2, -c1);
		else _F_mpz_submul_ui(res, coeff3, poly2, coeff2, c1);
		return;
	} 

	ulong c2 = poly2->coeffs[coeff2];
   
	if (!COEFF_IS_MPZ(c2)) // c2 is small
	{
		if ((long) c2 < 0) _F_mpz_addmul_ui(res, coeff3, poly1, coeff1, -c2);
		else _F_mpz_submul_ui(res, coeff3, poly1, coeff1, c2);
		return;
	} 

	// both c1 and c2 are large
  __mpz_struct * mpz_ptr = _F_mpz_promote_val(res, coeff3);
	
   mpz_submul(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1), poly2->mpz_coeffs + COEFF_TO_OFF(c2));
	_F_mpz_demote_val(res, coeff3); // cancellation may have occurred
}*/
