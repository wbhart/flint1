/*============================================================================

    F_mpz_poly.c: Polynomials over Z (FLINT 2.0 polynomials)

    Copyright (C) 2007, David Harvey (Odd/even Karatsuba) 
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

#include <stdint.h>
#include <string.h>
#include <math.h>

#include "mpz_poly.h"
#include "flint.h"
#include "F_mpz_poly.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "long_extras.h"
#include "zmod_poly.h"

/*===============================================================================

	Memory management

================================================================================*/

void _F_mpz_poly_mpz_coeffs_new(F_mpz_poly_t poly)
{
	if (poly->mpz_length == poly->mpz_alloc) // time to allocate MPZ_BLOCK more mpz_t's
	{
	   if (poly->mpz_alloc) // realloc mpz_t's
			poly->mpz_coeffs = (__mpz_struct*) flint_heap_realloc_bytes(poly->mpz_coeffs, (poly->mpz_alloc + MPZ_BLOCK)*sizeof(__mpz_struct));
		else // first time alloc of mpz_t'
			poly->mpz_coeffs = (__mpz_struct*) flint_heap_alloc_bytes(MPZ_BLOCK*sizeof(__mpz_struct));	
		
		// initialise the new mpz_t's
		for (ulong i = poly->mpz_alloc; i < poly->mpz_alloc + MPZ_BLOCK; i++)
			mpz_init(poly->mpz_coeffs + i);
		poly->mpz_alloc += MPZ_BLOCK;
	}

	poly->mpz_length++;
}

void _F_mpz_poly_mpz_coeffs_clear(F_mpz_poly_t poly)
{
   for (ulong i = 0; i < poly->mpz_alloc; i++) // clear any initialised mpz_t's
	   mpz_clear(poly->mpz_coeffs + i);

	flint_heap_free(poly->mpz_coeffs); // clear mpz_t array itself
	poly->mpz_coeffs = NULL;
	poly->mpz_alloc = 0;
	poly->mpz_length = 0;
}

void F_mpz_poly_init(F_mpz_poly_t poly)
{
   poly->coeffs = NULL;
   poly->mpz_coeffs = NULL;
   
   poly->alloc = 0;
   poly->length = 0;
   poly->mpz_alloc = 0;
   poly->mpz_length = 0;
}

void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc)
{
   if (alloc) // allocate space for alloc small coeffs
   {
      poly->coeffs = (mp_limb_t *) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
   }
   else poly->coeffs = NULL;
   poly->mpz_coeffs = NULL;

   poly->alloc = alloc;
   poly->length = 0;
   poly->mpz_alloc = 0;
   poly->mpz_length = 0;
}

void F_mpz_poly_realloc(F_mpz_poly_t poly, const ulong alloc)
{
   if (!alloc) // alloc == 0, clear up
   {
         F_mpz_poly_clear(poly);
			return;
   }  
   
	if (poly->alloc) // realloc
	{
		poly->coeffs = (mp_limb_t*) flint_heap_realloc(poly->coeffs, alloc);
		if (alloc > poly->alloc)
		   F_mpn_clear(poly->coeffs + poly->alloc, alloc - poly->alloc);
	} else // nothing allocated already so do it now
	{
		poly->coeffs = (mp_limb_t*) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
	}
   
   poly->alloc = alloc;
   
   F_mpz_poly_truncate(poly, alloc); // truncate actual data if necessary   
}

void F_mpz_poly_fit_length(F_mpz_poly_t poly, const ulong length)
{
   ulong alloc = length;
   
	if (alloc <= poly->alloc) return;

   // at least double number of allocated coeffs
	if (alloc < 2*poly->alloc) alloc = 2*poly->alloc; 
   
   F_mpz_poly_realloc(poly, alloc);
}

void F_mpz_poly_clear(F_mpz_poly_t poly)
{
   if (poly->coeffs) flint_heap_free(poly->coeffs); // clean up ordinary coeffs
   if (poly->mpz_coeffs) _F_mpz_poly_mpz_coeffs_clear(poly); // clean up mpz_t coeffs
   poly->coeffs = NULL;
	poly->alloc = 0;
	poly->length = 0;
}

/*===============================================================================

	Normalisation

================================================================================*/

void _F_mpz_poly_normalise(F_mpz_poly_t poly)
{
   ulong length = poly->length;
	
	while ((length) && (!poly->coeffs[length - 1])) length--;

	poly->length = length;
}

/*===============================================================================

	Coefficient operations

================================================================================*/

void _F_mpz_demote_val(F_mpz_poly_t poly, const ulong coeff)
{
   __mpz_struct * mpz_ptr = poly->mpz_coeffs + COEFF_TO_OFF(poly->coeffs[coeff]);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L) // coefficient is zero
	{
		_F_mpz_zero(poly, coeff);
	} else if (size == 1L) // coefficient is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = uval;
	} else if (size == -1L) // coefficient is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = -uval;
	}
	// don't do anything if coeff has to be multi precision
}

void _F_mpz_set_si(F_mpz_poly_t poly, ulong coeff, const long val)
{
   ulong uval = FLINT_ABS(val);
	
	if (uval > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(poly, coeff);
		mpz_set_si(mpz_coeff, val);
	} else poly->coeffs[coeff] = val; // val is small
}

void _F_mpz_set_ui(F_mpz_poly_t poly, ulong coeff, const ulong val)
{
   if (val > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(poly, coeff);
		mpz_set_ui(mpz_coeff, val);
	} else poly->coeffs[coeff] = val; // val is small
}

long _F_mpz_get_si(const F_mpz_poly_t poly, const ulong coeff)
{
   ulong c = poly->coeffs[coeff];

	if (!COEFF_IS_MPZ(c)) return c; // coeff is small
	return mpz_get_si(poly->mpz_coeffs + COEFF_TO_OFF(c)); // coeff is large
}

ulong _F_mpz_get_ui(const F_mpz_poly_t poly, const ulong coeff)
{
   ulong c = poly->coeffs[coeff];

	if (!COEFF_IS_MPZ(c)) return c; // coeff is small
	return mpz_get_ui(poly->mpz_coeffs + COEFF_TO_OFF(c)); //coeff is large
}

void _F_mpz_get_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong coeff)
{
	mp_limb_t c = poly->coeffs[coeff];

	if (!COEFF_IS_MPZ(c)) mpz_set_si(x, c); // set x to small coeff
	else mpz_set(x, poly->mpz_coeffs + COEFF_TO_OFF(c)); // set x to large coeff
}

void _F_mpz_set_mpz(F_mpz_poly_t poly, ulong coeff, const mpz_t x)
{
   long size = x->_mp_size;
	
	if (size == 0L) // x is zero
	{
		_F_mpz_zero(poly, coeff);
	} else if (size == 1L) // x is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(x);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = uval;
		else 
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(poly, coeff);
			mpz_set_ui(mpz_ptr, uval);
		}
	} else if (size == -1L) // x is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(x);
		if (uval <= COEFF_MAX) poly->coeffs[coeff] = -uval;
		else 
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(poly, coeff);
			mpz_set_ui(mpz_ptr, uval);
			mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // x is more than one limb
	{
		__mpz_struct * mpz_ptr = _F_mpz_promote(poly, coeff);
		mpz_set(mpz_ptr, x);
	}			
}

void _F_mpz_set(F_mpz_poly_t poly1, ulong coeff1, const F_mpz_poly_t poly2, const ulong coeff2)
{
   ulong c = poly2->coeffs[coeff2];
   
	if (!COEFF_IS_MPZ(c)) // coeff is small
	{
		poly1->coeffs[coeff1] = c;
	} else // coeff is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
		mpz_set(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c));
	}
}

int _F_mpz_equal(const F_mpz_poly_t poly1, const ulong coeff1, const F_mpz_poly_t poly2, const ulong coeff2)
{
	ulong c1 = poly1->coeffs[coeff1];
   ulong c2 = poly2->coeffs[coeff2];

	if (!COEFF_IS_MPZ(c1)) return (c1 == c2); // if c2 is large it can't be equal to c1
	else if (!COEFF_IS_MPZ(c2)) return 0; // c1 is large, so if c2 isn't....
	else return (!mpz_cmp(poly1->mpz_coeffs + COEFF_TO_OFF(c1), poly2->mpz_coeffs + COEFF_TO_OFF(c2))); 
}

void _F_mpz_swap(F_mpz_poly_t poly1, ulong coeff1, F_mpz_poly_t poly2, ulong coeff2)
{
	ulong c1 = poly1->coeffs[coeff1];
   ulong c2 = poly2->coeffs[coeff2];

   if (!COEFF_IS_MPZ(c1))
	{
		if (!COEFF_IS_MPZ(c2)) // both coefficients are small
		{
	      poly1->coeffs[coeff1] = c2;
	      poly2->coeffs[coeff2] = c1;
		} else // c1 is small, c2 is large
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
			mpz_swap(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c2));
			poly2->coeffs[coeff2] = c1;
		}
	} else
	{
      if (!COEFF_IS_MPZ(c2)) // c2 is small, c1 is large
		{
         __mpz_struct * mpz_ptr = _F_mpz_promote(poly2, coeff2);
			mpz_swap(mpz_ptr, poly1->mpz_coeffs + COEFF_TO_OFF(c1));
			poly1->coeffs[coeff1] = c2;
		} else // both coefficients are large
		{
			mpz_swap(poly1->mpz_coeffs + COEFF_TO_OFF(c1), poly2->mpz_coeffs + COEFF_TO_OFF(c2));
		}
	}
}

void _F_mpz_neg(F_mpz_poly_t poly1, ulong coeff1, const F_mpz_poly_t poly2, const ulong coeff2)
{
   ulong c = poly2->coeffs[coeff2];

	if (!COEFF_IS_MPZ(c)) // coeff is small
	{
		poly1->coeffs[coeff1] = -c;
	} else // coeff is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_promote(poly1, coeff1);
		mpz_neg(mpz_ptr, poly2->mpz_coeffs + COEFF_TO_OFF(c));
	}
}

void _F_mpz_add(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                 const F_mpz_poly_t poly2, const ulong coeff2) 
{
	mp_limb_t c1 = poly1->coeffs[coeff1];
	mp_limb_t c2 = poly2->coeffs[coeff2];
	
	if (!COEFF_IS_MPZ(c1))
	{
	   if (!COEFF_IS_MPZ(c2)) // both coefficients are small
		{
			_F_mpz_set_si(res, coeff3, c1 + c2);
		} else // c1 is small, c2 is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + COEFF_TO_OFF(c2);
			if ((long) c1 < 0L) mpz_sub_ui(mpz3, mpz2, -c1);	
		   else mpz_add_ui(mpz3, mpz2, c1);
			_F_mpz_demote_val(res, coeff3); // coefficients may have cancelled
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // c1 is large, c2 is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + COEFF_TO_OFF(c1);
			if ((long) c2 < 0L) mpz_sub_ui(mpz3, mpz1, -c2);	
			else mpz_add_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(res, coeff3); // coefficients may have cancelled
		} else // c1 and c2 are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + COEFF_TO_OFF(c2);
			mpz_add(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(res, coeff3); // coefficients may have cancelled
		}
	}
}

void _F_mpz_sub(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                 const F_mpz_poly_t poly2, const ulong coeff2) 
{
	mp_limb_t c1 = poly1->coeffs[coeff1];
	mp_limb_t c2 = poly2->coeffs[coeff2];
	
	if (!COEFF_IS_MPZ(c1))
	{
	   if (!COEFF_IS_MPZ(c2)) // both coefficients are small
		{
			_F_mpz_set_si(res, coeff3, c1 - c2);
		} else // c1 is small, c2 is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + COEFF_TO_OFF(c2);
			if ((long) c1 < 0L) 
			{
				mpz_add_ui(mpz3, mpz2, -c1);	
				mpz_neg(mpz3, mpz3);
			} else mpz_ui_sub(mpz3, c1, mpz2);
			_F_mpz_demote_val(res, coeff3); // coefficients may have cancelled
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // c1 is large, c2 is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + COEFF_TO_OFF(c1);
			if ((long) c2 < 0L) mpz_add_ui(mpz3, mpz1, -c2);	
			else mpz_sub_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(res, coeff3); // coefficients may have cancelled
		} else // c1 and c2 are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(res, coeff3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + COEFF_TO_OFF(c2);
			mpz_sub(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(res, coeff3); // coefficients may have cancelled
		}
	}
}

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
}

void F_mpz_poly_set_coeff_si(F_mpz_poly_t poly, ulong n, const long x)
{
   F_mpz_poly_fit_length(poly, n + 1);
   
	if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      for (ulong i = poly->length; i + 1 < n; i++)
         _F_mpz_zero(poly, i);
      poly->length = n+1;
   }
   
	_F_mpz_set_si(poly, n, x);
   _F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

void F_mpz_poly_set_coeff_ui(F_mpz_poly_t poly, ulong n, const ulong x)
{
   F_mpz_poly_fit_length(poly, n+1);

   if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      for (long i = poly->length; i + 1 < n; i++)
         _F_mpz_zero(poly, i); 
      poly->length = n+1;
   }

   _F_mpz_set_ui(poly, n, x);
   _F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

void F_mpz_poly_set_coeff_mpz(F_mpz_poly_t poly, ulong n, const mpz_t x)
{
   F_mpz_poly_fit_length(poly, n+1);

   if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      for (long i = poly->length; i + 1 < n; i++)
         _F_mpz_zero(poly, i); 
      poly->length = n+1;
   }

   _F_mpz_set_mpz(poly, n, x);
	_F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

long F_mpz_poly_get_coeff_si(const F_mpz_poly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
      return 0;
   
	return _F_mpz_get_si(poly, n);
}

ulong F_mpz_poly_get_coeff_ui(const F_mpz_poly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
      return 0;
   
	return _F_mpz_get_ui(poly, n);
}

void F_mpz_poly_get_coeff_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
	{
		mpz_set_ui(x, 0);
		return;
   }
   
	_F_mpz_get_mpz(x, poly, n);
	return;
}


/*===============================================================================

	Conversions

================================================================================*/

void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly)
{
	F_mpz_poly_fit_length(F_poly, m_poly->length);

	F_poly->length = m_poly->length;
   
	for (ulong i = 0; i < m_poly->length; i++)
		_F_mpz_set_mpz(F_poly, i, m_poly->coeffs[i]);
}

void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly)
{
	mpz_poly_ensure_alloc(m_poly, F_poly->length);

   m_poly->length = F_poly->length;
   
   for (ulong i = 0; i < F_poly->length; i++)
	   _F_mpz_get_mpz(m_poly->coeffs[i], F_poly, i);
}

/*===============================================================================

	Assignment/swap

================================================================================*/

void F_mpz_poly_set(F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	ulong length = poly2->length;
	
	if (poly1 != poly2)
	{
		F_mpz_poly_fit_length(poly1, poly2->length);

		for (ulong i = 0; i < poly2->length; i++)
			_F_mpz_set(poly1, i, poly2, i);
		
		poly1->length = poly2->length;
	}
}

void F_mpz_poly_swap(F_mpz_poly_t poly1, F_mpz_poly_t poly2)
{
	if (poly1 == poly2) return;

	ulong temp = poly1->length;
	poly1->length = poly2->length;
	poly2->length = temp;
	
	temp = poly1->alloc;
	poly1->alloc = poly2->alloc;
	poly2->alloc = temp;
	
	temp = poly1->mpz_alloc;
	poly1->mpz_alloc = poly2->mpz_alloc;
	poly2->mpz_alloc = temp;
	
	temp = poly1->mpz_length;
	poly1->mpz_length = poly2->mpz_length;
	poly2->mpz_length = temp;

	ulong * temp_c = poly1->coeffs;
	poly1->coeffs = poly2->coeffs;
	poly2->coeffs = temp_c;

   __mpz_struct * temp_m = poly1->mpz_coeffs;
	poly1->mpz_coeffs = poly2->mpz_coeffs;
	poly2->mpz_coeffs = temp_m;

   return;
}

/*===============================================================================

	Comparison

================================================================================*/

int F_mpz_poly_equal(const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   if (poly1 == poly2) return 1; // same polynomial
	if (poly1->length != poly2->length) return 0; // check if lengths the same

	for (ulong i = 0; i < poly1->length; i++) // check if coefficients the same
		if (!_F_mpz_equal(poly1, i, poly2, i)) return 0;

	return 1;
}

/*===============================================================================

	Coefficient sizes

================================================================================*/

long F_mpz_poly_max_bits(const F_mpz_poly_t poly)
{
	int sign = 0;
	ulong max = 0;
   ulong bits = 0;
   ulong i;
	ulong c;

	// search until we find an mpz_t coefficient or one of at least FLINT_BITS - 2 bits
	for (i = 0; i < poly->length; i++) 
	{
		c = poly->coeffs[i];
		if (COEFF_IS_MPZ(c)) break; // found an mpz_t coeff
      if ((long) c < 0L) 
		{
			sign = 1;
         bits = FLINT_BIT_COUNT(-c);
		} else bits = FLINT_BIT_COUNT(c);
		if (bits > max) max = bits;
		if (max == FLINT_BITS - 2) break; // coeff is at least FLINT_BITS - 2 bits
	}

   // search through mpz coefficients for largest size in bits
	for ( ; i < poly->length; i++)
	{
		c = poly->coeffs[i];
      if (COEFF_IS_MPZ(c))
		{
			__mpz_struct * mpz_ptr = poly->mpz_coeffs + COEFF_TO_OFF(c);
			if (mpz_sgn(mpz_ptr) < 0) sign = 1;
			bits = mpz_sizeinbase(mpz_ptr, 2);
			if (bits > max) max = bits;
		} else if ((long) c < 0L) sign = 1; // still need to check the sign of small coefficients
	}

	if (sign) return -max;
	else return max;
}

ulong F_mpz_poly_max_limbs(const F_mpz_poly_t poly)
{
	if (poly->length == 0) return 0; // polynomial is zero

	ulong max = 1; // all other coefficients have at least one limb
   ulong limbs;
	ulong c;

   // search through mpz coefficients for one of largest size
	for (ulong i = 0; i < poly->length; i++)
	{
		c = poly->coeffs[i];
      if (COEFF_IS_MPZ(c))
		{
			limbs = mpz_size(poly->mpz_coeffs + COEFF_TO_OFF(c));
			if (limbs > max) max = limbs;
		} 
	}

	return max;
}

/*===============================================================================

	Reverse

================================================================================*/

void F_mpz_poly_reverse(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong length)
{
   long i;
   
   F_mpz_poly_fit_length(res, length);
	
	if (poly != res) // not the same polynomial
   {
      for (i = 0; i < FLINT_MIN(length, poly->length); i++)
         _F_mpz_set(res, length - i - 1, poly, i); // copy over extant coefficients in reverse

      for ( ; i < length; i++) // set other coefficients to zero
         _F_mpz_zero(res, length - i - 1);

   } else // same polynomial
   {
      for (i = 0; i < length/2; i++)
      {
         // swap extant coefficients
			if (length - i - 1 < res->length) _F_mpz_swap(res, i, res, length - i - 1); 
			else
			{
				_F_mpz_set(res, length - i - 1, res, i); // for other coefficients "swap" with zero
			   _F_mpz_zero(res, i);
		   }
		}
      // if length is odd we missed a coefficient in swapping pairs, it may need to be set to zero
		if ((length & 1) && (i >= poly->length)) _F_mpz_zero(res, i); 
   }
	
	res->length = length;
   _F_mpz_poly_normalise(res); // new leading coeff, which was trailing coeff, may now be zero
}

/*===============================================================================

	Negation

================================================================================*/

void F_mpz_poly_neg(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
	F_mpz_poly_fit_length(res, poly->length);
	
	for (ulong i = 0; i < poly->length; i++)
		_F_mpz_neg(res, i, poly, i);

	res->length = poly->length;
}

/*===============================================================================

	Addition/subtraction

================================================================================*/

void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);
	ulong shorter = FLINT_MIN(poly1->length, poly2->length);

   F_mpz_poly_fit_length(res, longer);
	
   for (ulong i = 0; i < shorter; i++) // add up to the length of the shorter poly
      _F_mpz_add(res, i, poly1, i, poly2, i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (ulong i = shorter; i < poly1->length; i++)
         _F_mpz_set(res, i, poly1, i);

   if (poly2 != res) // copy any remaining coefficients from poly2
      for (ulong i = shorter; i < poly2->length; i++)
         _F_mpz_set(res, i, poly2, i);
   
   if (poly1->length == poly2->length)
   {
      res->length = poly1->length;
      _F_mpz_poly_normalise(res); // there may have been cancellation
   } else
      res->length = longer;
}

void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);
	ulong shorter = FLINT_MIN(poly1->length, poly2->length);

   F_mpz_poly_fit_length(res, longer);
	
   for (ulong i = 0; i < shorter; i++) // add up to the length of the shorter poly
      _F_mpz_sub(res, i, poly1, i, poly2, i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (ulong i = shorter; i < poly1->length; i++)
         _F_mpz_set(res, i, poly1, i);

   // careful, it is *always* necessary to negate coeffs from poly2, even if this is already res
	for (ulong i = shorter; i < poly2->length; i++) 
      _F_mpz_neg(res, i, poly2, i);

   if (poly1->length == poly2->length)
   {
      res->length = poly1->length;
      _F_mpz_poly_normalise(res); // there may have been cancellation
   } else
      res->length = longer;
}

/*===============================================================================

	Shifting

================================================================================*/

void F_mpz_poly_left_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n)
{
   if (n == 0) // special case, no shift
	{
		if (res != poly) F_mpz_poly_set(res, poly);
		return;
	}
	
	if (poly->length == 0) // nothing to shift
	{
		res->length = 0;
		return;
	}
	
	F_mpz_poly_fit_length(res, poly->length + n);
	
	// copy in reverse order to avoid writing over unshifted coeffs
	for (long i = poly->length - 1; i >= 0; i--) _F_mpz_set(res, i + n, poly, i);

   // insert n zeroes
	for (ulong i = 0; i < n; i++) _F_mpz_zero(res, i);
   
   res->length = poly->length + n;
}

void F_mpz_poly_right_shift(F_mpz_poly_t res, const F_mpz_poly_t poly, const ulong n)
{
   if (poly->length <= n) 
   {
      F_mpz_poly_zero(res);
      return;
   }

   F_mpz_poly_fit_length(res, poly->length - n);
	
	// copy in forward order to avoid writing over unshifted coeffs
	for (ulong i = 0; i < poly->length - n; i++) _F_mpz_set(res, i, poly, i + n);
	
	res->length = poly->length - n;
}

/*===============================================================================

	Interleave

================================================================================*/

void F_mpz_poly_interleave(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2)
{
	ulong length = FLINT_MAX(2*poly1->length - 1, 2*poly2->length);
	F_mpz_poly_fit_length(res, length + 1); // extra one in case length is odd

	ulong shorter = FLINT_MIN(poly1->length, poly2->length);
   ulong i;

	for (i = 0; i < shorter; i++)
	{
		_F_mpz_set(res, 2*i, poly1, i);
		_F_mpz_set(res, 2*i+1, poly2, i);
	}

   for ( ; i < poly1->length; i++)
	{
		_F_mpz_set(res, 2*i, poly1, i);
		_F_mpz_zero(res, 2*i+1);
	}

	for ( ; i < poly2->length; i++)
	{
		_F_mpz_zero(res, 2*i);
		_F_mpz_set(res, 2*i+1, poly2, i);
	}

	res->length = length;
}

void F_mpz_poly_interleave_small(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2)
{
	ulong length = FLINT_MAX(2*poly1->length - 1, 2*poly2->length);
	F_mpz_poly_fit_length(res, length + 1); // extra one in case length is odd

	ulong shorter = FLINT_MIN(poly1->length, poly2->length);
   ulong i;

	for (i = 0; i < shorter; i++)
	{
		res->coeffs[2*i] =  poly1->coeffs[i];
		res->coeffs[2*i+1] =  poly2->coeffs[i];
	}

   for ( ; i < poly1->length; i++)
	{
		res->coeffs[2*i] =  poly1->coeffs[i];
		res->coeffs[2*i+1] = 0;
	}

	for ( ; i < poly2->length; i++)
	{
		res->coeffs[2*i] = 0;
		res->coeffs[2*i+1] =  poly2->coeffs[i];
	}

	res->length = length;
}

/*===============================================================================

	Scalar multiplication

================================================================================*/

void F_mpz_poly_scalar_mul_ui(F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong x)
{
	// either scalar of input poly is zero
	if ((x == 0L) || (poly2->length == 0))  
	{
	   F_mpz_poly_zero(poly1);
		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   F_mpz_poly_set(poly1, poly2);
		return;
	}
	
	F_mpz_poly_fit_length(poly1, poly2->length);
	
	for (ulong i = 0; i < poly2->length; i++) _F_mpz_mul_ui(poly1, i, poly2, i, x);

	poly1->length = poly2->length;
}

void F_mpz_poly_scalar_mul_si(F_mpz_poly_t poly1, F_mpz_poly_t poly2, long x)
{
	// either scalar or input poly is zero
	if ((x == 0L) || (poly2->length == 0)) 
	{
	   F_mpz_poly_zero(poly1);
		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   F_mpz_poly_set(poly1, poly2);
		return;
	}
	
	// special case, multiply by -1
	if (x == -1L) 
	{
	   F_mpz_poly_neg(poly1, poly2);
		return;
	}
	
	F_mpz_poly_fit_length(poly1, poly2->length);
	
	for (ulong i = 0; i < poly2->length; i++) _F_mpz_mul_si(poly1, i, poly2, i, x);

	poly1->length = poly2->length;
}

void F_mpz_poly_scalar_mul_mpz(F_mpz_poly_t poly1, F_mpz_poly_t poly2, mpz_t x)
{
	// either scalar or input poly is zero
	if ((mpz_cmpabs_ui(x, 0L) == 0) || (poly2->length == 0)) 
	{
	   F_mpz_poly_zero(poly1);
		return;
	}
	
	// special cases, muliply by +/- 1
	if (mpz_cmpabs_ui(x, 1L) == 0)
	{
	   if (mpz_sgn(x) < 0) F_mpz_poly_neg(poly1, poly2);
		else F_mpz_poly_set(poly1, poly2);
		return;
	}
	
	F_mpz_poly_fit_length(poly1, poly2->length);
	
	for (ulong i = 0; i < poly2->length; i++) _F_mpz_mul_mpz(poly1, i, poly2, i, x);

	poly1->length = poly2->length;
}

/*===============================================================================

	Classical multiplication

================================================================================*/

void _F_mpz_poly_sqr_classical(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
   F_mpz_poly_fit_length(res, 2*poly->length - 1);
	res->length = 2*poly->length - 1;

   for (ulong i = 0; i < res->length; i++)
      _F_mpz_zero(res, i);
   
   // off-diagonal products
   for (ulong i = 1; i < poly->length; i++)
	{
		ulong c = poly->coeffs[i];
	   if (c)
		{
			if (!COEFF_IS_MPZ(c))
			{
				if ((long) c < 0L) 
					for (ulong j = 0; j < i; j++)
                  _F_mpz_submul_ui(res, i + j, poly, j, -c);
				else
               for (ulong j = 0; j < i; j++)
                  _F_mpz_addmul_ui(res, i + j, poly, j, c);
			} else
		      for (ulong j = 0; j < i; j++)
               _F_mpz_addmul(res, i+j, poly, i, poly, j);
		}
	}

   // double the off-diagonal products
   for (ulong i = 1; i < res->length - 1; i++)
      _F_mpz_add(res, i, res, i, res, i);
      
   // add in diagonal products
   for (ulong i = 0; i < poly->length; i++)
      _F_mpz_addmul(res, 2*i, poly, i, poly, i);
}

void F_mpz_poly_sqr_classical(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
   if (!poly->length)
   {
      // input is zero
      F_mpz_poly_zero(res);
      return;
   }
   
   ulong length = 2*poly->length - 1;
   
   if (res == poly)
   {
      // output is in place, so need a temporary
      F_mpz_poly_t temp;
      F_mpz_poly_init(temp);
      _F_mpz_poly_sqr_classical(temp, poly);
      F_mpz_poly_swap(temp, res);
      F_mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      _F_mpz_poly_sqr_classical(res, poly);
   }
}

void _F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;

	F_mpz_poly_fit_length(res, len1 + len2 - 1);
      
   if ((len1 == 1) && (len2 == 1)) // Special case if the length of both inputs is 1
   {
      _F_mpz_mul(res, 0, poly1, 0, poly2, 0);      
   } else // Ordinary case
   {
      long i, j;
      
      // Set res[i] = poly1[i]*poly2[0] 
      if (poly2->coeffs[0])
			for (i = 0; i < len1; i++)
            _F_mpz_mul(res, i, poly1, i, poly2, 0);
		else 
			for (i = 0; i < len1; i++)
            _F_mpz_zero(res, i);

      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      if (poly1->coeffs[len1 - 1])
		   for (i = 1; i < len2; i++)
            _F_mpz_mul(res, i + len1 - 1, poly1, len1 - 1, poly2, i);  
		else 
         for (i = 1; i < len2; i++)
            _F_mpz_zero(res, i + len1 - 1);
      
      // out[i+j] += in1[i]*in2[j] 
      for (i = 0; i < len1 - 1; i++)
      {      
         ulong c = poly1->coeffs[i];
			if (c)
			{
				if (!COEFF_IS_MPZ(c))
				{
					if ((long) c < 0L) 
						for (j = 1; j < len2; j++)
                     _F_mpz_submul_ui(res, i + j, poly2, j, -c);
					else
                  for (j = 1; j < len2; j++)
                     _F_mpz_addmul_ui(res, i + j, poly2, j, c);
				} else
					for (j = 1; j < len2; j++)
                  _F_mpz_addmul(res, i + j, poly1, i, poly2, j);
			}
      }
   } 
   
   res->length = len1 + len2 - 1;
}

void F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if (poly1 == poly2) F_mpz_poly_sqr_classical(res, poly1);
 
	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init(output);
		_F_mpz_poly_mul_classical(output, poly1, poly2);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
		_F_mpz_poly_mul_classical(res, poly1, poly2);
}

/*===============================================================================

	Karatsuba multiplication

================================================================================*/

void _F_mpz_poly_mul_kara_recursive(F_mpz_poly_t out, ulong ostart, F_mpz_poly_t in1, ulong istart1, 
											ulong len1, F_mpz_poly_t in2, ulong istart2, ulong len2, 
											F_mpz_poly_t scratch, ulong sstart, ulong skip, ulong crossover)
{
   // ==================== base cases
   
   if (len1 == 1)
   {
      // special case, just scalar multiplication
      for (ulong i = 0; i < len2; i++)
         _F_mpz_mul(out, ostart + i*skip, in1, istart1, in2, istart2 + i*skip);
      return;
   }
   
   if (len1 * len2 < crossover)
   {
      // switch to naive multiplication

		if (in2->coeffs[istart2])
			for (ulong i = 0; i < len1; i++)
            _F_mpz_mul(out, ostart + i*skip, in1, istart1 + i*skip, in2, istart2);
		else 
			for (ulong i = 0; i < len1; i++)
            _F_mpz_zero(out, ostart + i*skip);

      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      const ulong term = istart1 + (len1 - 1)*skip;
		if (in1->coeffs[istart1+(len1 - 1)*skip])
		   for (ulong i = 1; i < len2; i++)
            _F_mpz_mul(out, ostart + (i + len1 - 1)*skip, in1, term, in2, istart2 + i*skip);  
	 	else 
         for (ulong i = 1; i < len2; i++)
            _F_mpz_zero(out, ostart + (i + len1 - 1)*skip);
      
      for (ulong i = 0; i < len1 - 1; i++)
		{
			ulong c = in1->coeffs[istart1 + i*skip];
			const ulong term = istart1 + i*skip;
			if (!COEFF_IS_MPZ(c))
			{
				if ((long) c < 0L) 
					for (ulong j = 1; j < len2; j++)
                  _F_mpz_submul_ui(out, ostart + (i+j)*skip, in2, istart2 + j*skip, -c);
				else
               for (ulong j = 1; j < len2; j++)
                  _F_mpz_addmul_ui(out, ostart + (i+j)*skip, in2, istart2 + j*skip, c);
			} else
				for (ulong j = 1; j < len2; j++)
               _F_mpz_addmul(out, ostart + (i+j)*skip, in1, term, in2, istart2 + j*skip);
		}    

		return;
   }

   ulong i, j;
	
	// ==================== recursive case

   // Let in1 = A1(x^2) + x*B1(x^2) + x^(2*floor(len1/2))*C1,
   // where A1, B1 have length floor(len1/2),
   // and C1 is the leading term of in1 if len1 is odd

   // Similarly for in2 = A2(x^2) + x*B2(x^2) + x^(2*floor(len2/2))*C2.
   
   // Put A1 + B1 into even slots of scratch space
   // (uses len1/2 scratch slots)
   for (i = 0; i < len1/2; i++)
      _F_mpz_add(scratch, sstart + 2*i*skip, in1, istart1 + 2*i*skip, in1, istart1 + 2*i*skip + skip);

   // Put A2 + B2 into remaining even slots of scratch space
   // (uses len2/2 slots of scratch)
   for (j = 0; j < len2/2; j++)
      _F_mpz_add(scratch, sstart + 2*(i+j)*skip, in2, istart2 + 2*j*skip, in2, istart2 + 2*j*skip + skip);

   // The following three recursive calls all use the odd slots of the current
   // scratch array as the next layer's scratch space
   
   // Put product (A1+B1)*(A2+B2) into odd slots of output array
   _F_mpz_poly_mul_kara_recursive(out, ostart + skip, scratch, sstart, len1/2, scratch, sstart + 2*i*skip, len2/2,
                                scratch, sstart + skip, 2*skip, crossover);

   // Put product x^2*(B1*B2) into even slots of output array
   // (except first slot, which is an implied zero)
   _F_mpz_poly_mul_kara_recursive(out, ostart + 2*skip, in1, istart1 + skip, len1/2, in2, istart2 + skip,
                                len2/2, scratch, sstart + skip, 2*skip, crossover);

   // Put product A1*A2 into even slots of scratch space
   _F_mpz_poly_mul_kara_recursive(scratch, sstart, in1, istart1, len1/2, in2, istart2, len2/2,
                                scratch, sstart + skip, 2*skip, crossover);
                            
   // Subtract A1*A2 and B1*B2 from (A1+B1)*(A2+B2) to get (A1*B2 + A2*B1)
   // in odd slots of output
   for (ulong i = 0; i < len1/2 + len2/2 - 1; i++)
   {
      _F_mpz_sub(out, ostart + 2*i*skip + skip, out, ostart + 2*i*skip + skip, out, ostart + 2*(i+1)*skip);
      _F_mpz_sub(out, ostart + 2*i*skip + skip, out, ostart + 2*i*skip + skip, scratch, sstart + 2*i*skip);
   }
      
   // Add A1*A2 to x^2*(B1*B2) into even slots of output
   _F_mpz_set(out, ostart, scratch, sstart);
   for (ulong i = 1; i < len1/2 + len2/2 - 1; i++)
      _F_mpz_add(out, ostart + 2*i*skip, out, ostart + 2*i*skip, scratch, sstart + 2*i*skip);
   
   // Now we have the product (A1(x^2) + x*B1(x^2)) * (A2(x^2) + x*B2(x^2))
   // in the output array. Still need to handle C1 and C2 terms.
   
   if (len1 & 1)
   {
      const ulong term1 = istart1 + skip*(len1-1);
	   const ulong term2 = istart2 + skip*(len2-1);
	   if (len2 & 1)
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         for (ulong i = 0; i < len2-2; i++)
            _F_mpz_addmul(out, ostart + (i+len1-1)*skip, in1, term1, in2, istart2 + i*skip);
         _F_mpz_mul(out, ostart + (len1+len2-3)*skip, in1, term1, in2, istart2 + (len2-2)*skip);

         // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
         for (ulong i = 0; i < len1-1; i++)
            _F_mpz_addmul(out, ostart + (i+len2-1)*skip, in2, term2, in1, istart1 + i*skip);
            
         // final C1*C2 term
         _F_mpz_mul(out, ostart + (len1+len2-2)*skip, in1, term1, in2,term2);
      }
      else
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         for (ulong i = 0; i < len2-1; i++)
            _F_mpz_addmul(out, ostart + (i+len1-1)*skip, in1, term1, in2, istart2 + i*skip);
         _F_mpz_mul(out, ostart + (len1+len2-2)*skip, in1, term1, in2, term2);
      }
   }
   else if (len2 & 1)
   {
      const ulong term1 = istart1 + skip*(len1-1);
	   const ulong term2 = istart2 + skip*(len2-1);
	   // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
      for (ulong i = 0; i < len1-1; i++)
         _F_mpz_addmul(out, ostart + (i+len2-1)*skip, in2, term2, in1, istart1 + i*skip);
      _F_mpz_mul(out, ostart + (len1+len2-2)*skip, in2, term2, in1, term1);
   }
}

ulong F_mpz_poly_kara_table[] = {8, 27, 10, 10, 8, 6, 5, 5, 5, 4, 4, 2};
ulong F_mpz_poly_kara_table_size = 12;

ulong _F_mpz_poly_mul_karatsuba_crossover(ulong limbs)
{
   ulong crossover;

   if (limbs >= F_mpz_poly_kara_table_size)
      crossover = 0;
   else
      crossover = F_mpz_poly_kara_table[limbs - 1];

   return crossover * crossover;
}

void _F_mpz_poly_mul_karatsuba(F_mpz_poly_t res, F_mpz_poly_t poly1,
                            F_mpz_poly_t poly2)
{
	// number of output coefficients, and a rough upper bound on the number
   // of limbs needed for each one
   ulong length = poly1->length + poly2->length - 1;
   
   // allocate scratch space for lower-level karatsuba routine
   F_mpz_poly_t scratch;
	F_mpz_poly_init(scratch);
	F_mpz_poly_fit_length(scratch, length + 1);
	scratch->length = length + 1;
	
	// look up crossover parameter (i.e. when to switch from classical to
   // karatsuba multiplication) based on coefficient size
   ulong bits1 = FLINT_ABS(F_mpz_poly_max_bits(poly1));
	ulong bits2 = FLINT_ABS(F_mpz_poly_max_bits(poly2));
	ulong log_length = 0;
	while ((1L<<log_length) < poly1->length) log_length++;
	ulong limbs = (bits1 + bits2 + log_length - 1)/FLINT_BITS + 1;

	ulong crossover = _F_mpz_poly_mul_karatsuba_crossover(limbs);
   
   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      F_mpz_poly_t temp;
      F_mpz_poly_init2(temp, length);
		
      _F_mpz_poly_mul_kara_recursive(
            temp, 0, poly1, 0, poly1->length,
            poly2, 0, poly2->length, scratch, 0, 1, crossover);

      temp->length = length;
		
		F_mpz_poly_swap(temp, res);
      F_mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      _F_mpz_poly_mul_kara_recursive(
            res, 0, poly1, 0, poly1->length,
            poly2, 0, poly2->length, scratch, 0, 1, crossover);

		res->length = length;
   }
   
   F_mpz_poly_clear(scratch);
}

void F_mpz_poly_mul_karatsuba(F_mpz_poly_t res, F_mpz_poly_t poly1,
                            F_mpz_poly_t poly2)
{
   if (!poly1->length || !poly2->length)
   {
      // one of the polys is zero
      F_mpz_poly_zero(res);
      return;
   }
   
   /*if (poly1 == poly2)
   {
      // polys are identical, so call specialised squaring routine
      F_mpz_poly_sqr_karatsuba(res, poly1);
      return;
   }*/

   F_mpz_poly_fit_length(res, poly1->length + poly2->length - 1);
	
	// rearrange parameters to make poly1 no longer than poly2
   if (poly1->length <= poly2->length)
      _F_mpz_poly_mul_karatsuba(res, poly1, poly2);
	else 
      _F_mpz_poly_mul_karatsuba(res, poly2, poly1);
}

/*===============================================================================

	Bit/byte/limb packing

================================================================================*/

/*
   Read the next coefficient, xxxarr, from the polynomial, negate it if xxxneg is -1,
	subtract borrow. If the coefficient is negative, set borrow to be 1. Mask the coefficient
	against xxxmask.
*/
#define NEXT_COEFF(xxxarr, xxxborr, xxxcoeff, xxxmask, xxxneg) \
	do { \
	   xxxcoeff = (xxxarr ^ xxxneg) - xxxneg - xxxborr; \
		xxxborr = 0UL; \
		if (xxxcoeff < 0) xxxborr = 1UL; \
		xxxcoeff &= xxxmask; \
	} while (0)

void F_mpz_poly_bit_pack(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                         const long bitwidth, const ulong length, const long negate)
{   
   ulong i, k, skip;
   ulong * coeff_m = poly_F_mpz->coeffs;
   mp_limb_t * end_point;
    
   ulong temp;
   half_ulong lower;
   long coeff;
   long borrow;
   mp_limb_t extend;
   
   long bits = bitwidth;
   int sign = (bits < 0);
   if (sign) bits = ABS(bits);
   
   ulong coeffs_per_limb = FLINT_BITS/bits;

   const ulong mask = (1UL<<bits)-1;
      
      k = 0; skip = 0;
      coeff = 0; borrow = 0L; temp = 0;
      
      end_point = poly_F_mpz->coeffs + length;
         
      while (coeff_m < end_point)
      {
         if ((ulong) coeff_m & 7 == 0) FLINT_PREFETCH(coeff_m, 64);

         // k is guaranteed to be less than FLINT_BITS at this point
         while ((k < HALF_FLINT_BITS) && (coeff_m < end_point))
         {
            NEXT_COEFF(*coeff_m, borrow, coeff, mask, negate);
				//if (sign) 
			   temp += (coeff << k);
            //else temp += (__get_next_coeff_unsigned(coeff_m, &coeff) << k);
            coeff_m++; k += bits;
         }
         
			// k may exceed FLINT_BITS at this point but is less than 3*HALF_FLINT_BITS
			if (k > FLINT_BITS)
         {
            // if k > FLINT_BITS write out a whole limb and read in remaining bits of coeff
            array[skip] = temp;
            skip++;
            temp = (coeff >> (bits+FLINT_BITS-k));
            k = (k - FLINT_BITS);
            // k < HALF_FLINT_BITS
         } else
         {
            // k <= FLINT_BITS
            if (k >= HALF_FLINT_BITS)
            {
               // if k >= HALF_FLINT_BITS store bottom HALF_FLINT_BITS bits
               lower = (half_ulong) temp;
               k -= HALF_FLINT_BITS;
               temp >>= HALF_FLINT_BITS;
               
					// k is now <= HALF_FLINT_BITS
               while ((k < HALF_FLINT_BITS) && (coeff_m < end_point))
               {
                  NEXT_COEFF(*coeff_m, borrow, coeff, mask, negate);
				      //if (sign) 
						temp+=(coeff << k);
                  //else temp+=(__get_next_coeff_unsigned(coeff_m, &coeff) << k);
                  coeff_m++; k += bits;
               }

               // k may again exceed FLINT_BITS bits but is less than 3*HALF_FLINT_BITS
               if (k > FLINT_BITS)
               {
                  // if k > FLINT_BITS, write out bottom HALF_FLINT_BITS bits (along with HALF_FLINT_BITS bits from lower)
                  // read remaining bits from coeff and reduce k by HALF_FLINT_BITS
                  array[skip] = (temp << HALF_FLINT_BITS) + (ulong) lower;
                  skip++;
                  temp >>= HALF_FLINT_BITS;
                  temp += ((coeff >> (bits+FLINT_BITS-k)) << HALF_FLINT_BITS);
                  k = (k - HALF_FLINT_BITS);
                  // k < FLINT_BITS and we are ready to read next coefficient if there is one
               } else if (k >= HALF_FLINT_BITS) 
               {
                  // k <= FLINT_BITS
                  // if k >= HALF_FLINT_BITS write out bottom HALF_FLINT_BITS bits (along with lower)
                  // and reduce k by HALF_FLINT_BITS
                  k -= HALF_FLINT_BITS;
                  array[skip] = (temp << HALF_FLINT_BITS) + lower;
                  temp >>= HALF_FLINT_BITS;
                  skip++;
                  // k is now less than or equal to HALF_FLINT_BITS and we are now ready to read 
                  // the next coefficient if there is one
               } else
               {
                  // k < HALF_FLINT_BITS
                  // there isn't enough to write out a whole FLINT_BITS bits, so put it all 
                  // together in temp
                  temp = (temp << HALF_FLINT_BITS) + lower;
                  k += HALF_FLINT_BITS;
                  // k is now guaranteed to be less than FLINT_BITS and we are ready for the
                  // next coefficient if there is one
               }
            } // if
         } // else
      } // while

      // sign extend the last FLINT_BITS bits we write out
      if (skip < n)
      {
        if (borrow) temp += (-1UL << k);
        array[skip] = temp;
        skip++;
      } 
}

void F_mpz_poly_bit_pack_unsigned(mp_limb_t * array, ulong n, const F_mpz_poly_t poly_F_mpz,
                                                            const ulong bits, const ulong length)
{   
   ulong i, k, skip;
   ulong * coeff_m = poly_F_mpz->coeffs;
   mp_limb_t * end_point;
    
   ulong temp;
   half_ulong lower;
   long coeff;
   
   ulong coeffs_per_limb = FLINT_BITS/bits;

      k = 0; skip = 0;
      coeff = 0; temp = 0;
      
      end_point = poly_F_mpz->coeffs + length;
         
      while (coeff_m < end_point)
      {
         if ((ulong) coeff_m & 7 == 0) FLINT_PREFETCH(coeff_m, 64);

         // k is guaranteed to be less than FLINT_BITS at this point
         while ((k < HALF_FLINT_BITS) && (coeff_m < end_point))
         {
            coeff = *coeff_m;
				temp += (coeff << k);
            coeff_m++; k += bits;
         }
         
			// k may exceed FLINT_BITS at this point but is less than 3*HALF_FLINT_BITS
			if (k > FLINT_BITS)
         {
            // if k > FLINT_BITS write out a whole limb and read in remaining bits of coeff
            array[skip] = temp;
            skip++;
            temp = (coeff >> (bits+FLINT_BITS-k));
            k = (k - FLINT_BITS);
            // k < HALF_FLINT_BITS
         } else
         {
            // k <= FLINT_BITS
            if (k >= HALF_FLINT_BITS)
            {
               // if k >= HALF_FLINT_BITS store bottom HALF_FLINT_BITS bits
               lower = (half_ulong) temp;
               k -= HALF_FLINT_BITS;
               temp >>= HALF_FLINT_BITS;
               
					// k is now <= HALF_FLINT_BITS
               while ((k < HALF_FLINT_BITS) && (coeff_m < end_point))
               {
                  coeff = *coeff_m;
				      temp += (coeff << k);
                  coeff_m++; k += bits;
               }

               // k may again exceed FLINT_BITS bits but is less than 3*HALF_FLINT_BITS
               if (k > FLINT_BITS)
               {
                  // if k > FLINT_BITS, write out bottom HALF_FLINT_BITS bits (along with HALF_FLINT_BITS bits from lower)
                  // read remaining bits from coeff and reduce k by HALF_FLINT_BITS
                  array[skip] = (temp << HALF_FLINT_BITS) + (ulong) lower;
                  skip++;
                  temp >>= HALF_FLINT_BITS;
                  temp += ((coeff >> (bits+FLINT_BITS-k)) << HALF_FLINT_BITS);
                  k = (k - HALF_FLINT_BITS);
                  // k < FLINT_BITS and we are ready to read next coefficient if there is one
               } else if (k >= HALF_FLINT_BITS) 
               {
                  // k <= FLINT_BITS
                  // if k >= HALF_FLINT_BITS write out bottom HALF_FLINT_BITS bits (along with lower)
                  // and reduce k by HALF_FLINT_BITS
                  k -= HALF_FLINT_BITS;
                  array[skip] = (temp << HALF_FLINT_BITS) + lower;
                  temp >>= HALF_FLINT_BITS;
                  skip++;
                  // k is now less than or equal to HALF_FLINT_BITS and we are now ready to read 
                  // the next coefficient if there is one
               } else
               {
                  // k < HALF_FLINT_BITS
                  // there isn't enough to write out a whole FLINT_BITS bits, so put it all 
                  // together in temp
                  temp = (temp << HALF_FLINT_BITS) + lower;
                  k += HALF_FLINT_BITS;
                  // k is now guaranteed to be less than FLINT_BITS and we are ready for the
                  // next coefficient if there is one
               }
            } // if
         } // else
      } // while

      // sign extend the last FLINT_BITS bits we write out
      if (skip < n)
      {
        array[skip] = temp;
        skip++;
      } 
}

void F_mpz_poly_bit_unpack(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, ulong n,
                                                    const ulong length, const ulong bits)
{
   ulong k, skip;

   ulong temp;
   ulong full_limb;
   ulong carry;
     
   const ulong mask = (1UL << bits) - 1;
   const ulong sign_mask = (1UL << (bits-1));

   ulong s;
   mp_limb_t * coeff_m = poly_F_mpz->coeffs;
   mp_limb_t * end_point;

      k = 0; skip = 0; carry = 0UL; temp = 0L;
      end_point = poly_F_mpz->coeffs + length;
      
      while (coeff_m < end_point)
      {
         // read in a full limb
         full_limb = array[skip];
         temp += l_shift(full_limb, k);
         s = FLINT_BITS - k;
         k += s;
         while ((k >= bits) && (coeff_m < end_point))
         {
            if (!(temp & sign_mask)) 
            {
               *coeff_m += ((temp & mask) + carry);
               carry = 0UL;
            }  
            else
            {
               *coeff_m -= (((-temp) & mask) - carry);
               carry = 1UL;
            }
            coeff_m ++;
            temp >>= bits;
            k -= bits;
         }

         // k is now less than bits
         // read in remainder of full_limb
         temp += l_shift(r_shift(full_limb, s), k);
         k += (FLINT_BITS - s);
       
         while ((k >= bits) && (coeff_m < end_point))
         {
            if (!(temp & sign_mask)) 
            {
               *coeff_m += ((temp & mask) + carry);
               carry = 0UL;
            }
            else
            {
               *coeff_m -= (((-temp) & mask) - carry);
					carry = 1UL;
            }
            
				coeff_m ++;
            temp >>= bits;
            k -= bits;
         }

         // k is now less than bits
         skip++;
      }

   _F_mpz_poly_normalise(poly_F_mpz);
}

void F_mpz_poly_bit_unpack_unsigned(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, ulong n, 
                                                             const ulong length, const ulong bits)
{
   ulong k, l, skip;

   ulong temp;
   ulong full_limb;
    
   const ulong mask = (1UL<<bits)-1;

   ulong s;
   ulong * coeff_m = poly_F_mpz->coeffs;
   ulong * end_point;
       
      k = 0; skip = 0; temp = 0;
      end_point = poly_F_mpz->coeffs + poly_F_mpz->length;
      
      while (coeff_m < end_point)
      {
         if (skip & 7 == 0) FLINT_PREFETCH(array + skip, 64);
         // read in a full limb
         full_limb = array[skip];
         temp += l_shift(full_limb, k);
         s = FLINT_BITS - k;
         k += s;
         
			while ((k >= bits) && (coeff_m < end_point))
         {
            *coeff_m += (temp & mask);
            coeff_m++;
            temp >>= bits;
            k -= bits;
         }

         // k is now less than bits
         // read in remainder of full_limb
         temp += l_shift(r_shift(full_limb, s), k);
         k += (FLINT_BITS - s);
       
         while ((k >= bits) && (coeff_m < end_point))
         {
            *coeff_m += (temp & mask);
            coeff_m++;
            temp >>= bits;
            l++;
            k -= bits;
         }
         // k is now less than bits
         skip++;
      }

   _F_mpz_poly_normalise(poly_F_mpz);
}
 
/*
  Bit packing used with David Harvey's KS2 algorithm

*/
void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
                              const long bitwidth, const ulong length, const long negate, long negate2)

{   
   ulong k, skip;
   ulong * coeff_m = poly_F_mpz->coeffs;
   
   ulong temp, temp2;
   half_ulong lower, lower2;
   long coeff, coeff2;
   long borrow, borrow2;
   
   long bits = bitwidth;
   int sign = (bits < 0);
   if (sign) bits = ABS(bits);
   
   ulong coeffs_per_limb = FLINT_BITS/bits;

   const ulong mask = (1UL<<bits) - 1;
      
	   k = 0; skip = 0;
      coeff = 0; coeff2 = 0; borrow = 0L; borrow2 = 0L, temp = 0, temp2 = 0;
		
      mp_limb_t * end_point = poly_F_mpz->coeffs + length;
         
      while (coeff_m < end_point)
      {
         if ((ulong) coeff_m & 7 == 0) FLINT_PREFETCH(coeff_m, 64);

         // k is guaranteed to be less than FLINT_BITS at this point
         while ((k < HALF_FLINT_BITS) && (coeff_m < end_point))
         {
            NEXT_COEFF(*coeff_m, borrow, coeff, mask, negate);
				NEXT_COEFF(*coeff_m, borrow2, coeff2, mask, negate2);
				negate2 = -negate2 - 1L;
				//if (sign) 
			   temp += (coeff << k);
            temp2 += (coeff2 << k);
            //else temp += (__get_next_coeff_unsigned(coeff_m, &coeff) << k);
            coeff_m++; k += bits;
         }
         
			// k may exceed FLINT_BITS at this point but is less than 3*HALF_FLINT_BITS
			if (k > FLINT_BITS)
         {
            // if k > FLINT_BITS write out a whole limb and read in remaining bits of coeff
            array[skip] = temp;
            array2[skip] = temp2;
            skip++;
            temp = (coeff >> (bits+FLINT_BITS-k));
            temp2 = (coeff2 >> (bits+FLINT_BITS-k));
            k = (k - FLINT_BITS);
            // k < HALF_FLINT_BITS
         } else
         {
            // k <= FLINT_BITS
            if (k >= HALF_FLINT_BITS)
            {
               // if k >= HALF_FLINT_BITS store bottom HALF_FLINT_BITS bits
               lower = (half_ulong) temp;
               lower2 = (half_ulong) temp2;
               k -= HALF_FLINT_BITS;
               temp >>= HALF_FLINT_BITS;
               temp2 >>= HALF_FLINT_BITS;
               
					// k is now <= HALF_FLINT_BITS
               while ((k < HALF_FLINT_BITS) && (coeff_m < end_point))
               {
                  NEXT_COEFF(*coeff_m, borrow, coeff, mask, negate);
				      NEXT_COEFF(*coeff_m, borrow2, coeff2, mask, negate2);
						negate2 = -negate2 - 1L;
				      //if (sign) 
						temp+=(coeff << k);
                  temp2+=(coeff2 << k);
                  //else temp+=(__get_next_coeff_unsigned(coeff_m, &coeff) << k);
                  coeff_m++; k += bits;
               }

               // k may again exceed FLINT_BITS bits but is less than 3*HALF_FLINT_BITS
               if (k > FLINT_BITS)
               {
                  // if k > FLINT_BITS, write out bottom HALF_FLINT_BITS bits (along with HALF_FLINT_BITS bits from lower)
                  // read remaining bits from coeff and reduce k by HALF_FLINT_BITS
                  array[skip] = (temp << HALF_FLINT_BITS) + (ulong) lower;
                  array2[skip] = (temp2 << HALF_FLINT_BITS) + (ulong) lower2;
                  skip++;
                  temp >>= HALF_FLINT_BITS;
                  temp2 >>= HALF_FLINT_BITS;
                  temp += ((coeff >> (bits+FLINT_BITS-k)) << HALF_FLINT_BITS);
                  temp2 += ((coeff2 >> (bits+FLINT_BITS-k)) << HALF_FLINT_BITS);
                  k = (k - HALF_FLINT_BITS);
                  // k < FLINT_BITS and we are ready to read next coefficient if there is one
               } else if (k >= HALF_FLINT_BITS) 
               {
                  // k <= FLINT_BITS
                  // if k >= HALF_FLINT_BITS write out bottom HALF_FLINT_BITS bits (along with lower)
                  // and reduce k by HALF_FLINT_BITS
                  k -= HALF_FLINT_BITS;
                  array[skip] = (temp << HALF_FLINT_BITS) + lower;
                  array2[skip] = (temp2 << HALF_FLINT_BITS) + lower2;
                  temp >>= HALF_FLINT_BITS;
                  temp2 >>= HALF_FLINT_BITS;
                  skip++;
                  // k is now less than or equal to HALF_FLINT_BITS and we are now ready to read 
                  // the next coefficient if there is one
               } else
               {
                  // k < HALF_FLINT_BITS
                  // there isn't enough to write out a whole FLINT_BITS bits, so put it all 
                  // together in temp
                  temp = (temp << HALF_FLINT_BITS) + lower;
                  temp2 = (temp2 << HALF_FLINT_BITS) + lower2;
                  k += HALF_FLINT_BITS;
                  // k is now guaranteed to be less than FLINT_BITS and we are ready for the
                  // next coefficient if there is one
               }
            } // if
         } // else
      } // while

      // sign extend the last FLINT_BITS bits we write out
      if (skip < n)
      {
        //if (borrow) temp+= (-1UL << k);
        array[skip] = temp;
        array2[skip] = temp2;
        skip++;
      } 
}


/*===============================================================================

	Kronecker Segmentation multiplication

================================================================================*/

void _F_mpz_poly_mul_KS(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2)
{
   long sign1 = 0L;
   long sign2 = 0L;
   
   ulong length1 = input1->length;
   ulong length2 = input2->length;
   
   ulong final_length = length1 + length2 - 1;
   
	F_mpz_poly_fit_length(output, final_length);

   long bits1, bits2;
   int bitpack = 0;
   
   bits1 = F_mpz_poly_max_bits(input1);
   bits2 = (input1 == input2) ? bits1 : F_mpz_poly_max_bits(input2);
      
   ulong sign = ((bits1 < 0) || (bits2 < 0));
   ulong length = length2;
   ulong log_length = 0L;
   while ((1<<log_length) < length) log_length++;
   ulong bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   if (bits - sign <= FLINT_BITS - 2) bitpack = 1;
   
   if ((long) input1->coeffs[length1 - 1] < 0L) sign1 = -1L;
   
   if (input1 != input2)
   {
      if ((long) input2->coeffs[length2 - 1] < 0L) sign2 = -1L;
   } else sign2 = sign1;
   
   mp_limb_t * int1, * int2, * int3;
	ulong n1, n2;

   if (bitpack)
   {
      n1 = (bits*length1 - 1)/FLINT_BITS + 1;
		n2 = (bits*length2 - 1)/FLINT_BITS + 1;

		int1 = (mp_limb_t*) flint_stack_alloc(n1);
      if (input1 != input2)
         int2 = (mp_limb_t*) flint_stack_alloc(n2);
        
      if (sign) 
		{
			bits = -1L*bits;
         if (input1 != input2)
            F_mpz_poly_bit_pack(int2, n2, input2, bits, length2, sign2);
         F_mpz_poly_bit_pack(int1, n1, input1, bits, length1, sign1);
         bits=ABS(bits);
		} else
		{
         if (input1 != input2)
            F_mpz_poly_bit_pack_unsigned(int2, n2, input2, bits, length2);
         F_mpz_poly_bit_pack_unsigned(int1, n1, input1, bits, length1);
		}
   } 
	
   if (input1 == input2)
      int2 = int1;
   
   int3 = (mp_limb_t*) flint_stack_alloc(n1 + n2);
           
   mp_limb_t msl = F_mpn_mul(int3, int1, n1, int2, n2);
   
   int3[n1 + n2 - 1] = msl;
   
   output->length = length1 + length2 - 1;
  
   for (unsigned long i = 0; i < output->length; i++)
      _F_mpz_zero(output, i);
      
   if (bitpack)
   {
      if (sign) F_mpz_poly_bit_unpack(output, int3, n1 + n2, length1 + length2 - 1, bits);  
      else F_mpz_poly_bit_unpack_unsigned(output, int3, n1 + n2, length1 + length2 - 1, bits);  
   } 
	
   flint_stack_release(); // release int3
   if (input1 != input2)
      flint_stack_release(); // release int2
   flint_stack_release(); // release int1
     
   if ((long) (sign1 ^ sign2) < 0L) F_mpz_poly_neg(output, output);
}

/*
  David Harvey's KS2 algorithm.
  This is faster than ordinary KS for about lengths 300-4000. (The overhead prevents it
  from being faster for smaller multiplications, the FFT's quasilinear run time prevents
  it from being faster for larger lengths).
  It uses the identity f(x)*g(x) = (f(x)*g(x) + f(-x)*g(-x))/2 + (f(x)*g(x) - f(-x)*g(-x))/2
  to reduce multiplication to two half sized multiplications which is faster when the 
  multiplication is worse than quasi-linear time. The multiplications are half the size
  because the first bracketed term only has the even exponent coefficients and the other
  bracketed term has the odd coefficients. In both cases each of the output coefficients 
  of the multiplications can overflow right through the following (zero) coefficient.
*/
void _F_mpz_poly_mul_KS2(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2)
{
   long sign1 = 0L;
	long sign2 = 0L;
   long sign1a, sign2a;
   
   ulong length1 = input1->length;
   ulong length2 = input2->length;
   
   ulong final_length = length1 + length2 - 1;
   
	F_mpz_poly_fit_length(output, final_length);

   if ((long) input1->coeffs[length1 - 1] < 0L) sign1 = -1L;
   if (length1 & 1) sign1a = sign1;
	else sign1a = -sign1 - 1L;

   if (input1 != input2)
   {
      if ((long) input2->coeffs[length2 - 1] < 0L) sign2 = -1L;
   } else sign2 = sign1;
   if (length2 & 1) sign2a = sign2;
	else sign2a = -sign2 - 1L;

	long bits1, bits2;
   int bitpack = 0;
   
   bits1 = F_mpz_poly_max_bits(input1);
   bits2 = (input1 == input2) ? bits1 : F_mpz_poly_max_bits(input2);
      
   ulong sign = ((bits1 < 0) || (bits2 < 0));
   ulong length = length2;
   ulong log_length = 0L;
   while ((1<<log_length) < length) log_length++;
	ulong bits_max = FLINT_MAX(ABS(bits1), ABS(bits2));
   ulong bits = FLINT_MAX(ABS(bits1) + ABS(bits2) + log_length + 1, 2*bits_max); 
	if ((bits%2) == 1) bits++;
   
	if (bits - sign <= FLINT_BITS - 2) bitpack = 1;
   
   mp_limb_t * int1, * int2, * int1b, * int2b, * int3, * int3b, * int4;
	ulong n1, n2;
   
	if (bitpack)
   {
      bits /= 2;
	   n1 = (bits*length1 - 1)/FLINT_BITS + 1;
		n2 = (bits*length2 - 1)/FLINT_BITS + 1;
		int1 = (mp_limb_t*) flint_stack_alloc(n1);
      int1b = (mp_limb_t*) flint_stack_alloc(n1);
      if (input1 != input2)
		{
			int2 = (mp_limb_t*) flint_stack_alloc(n2);
         int2b = (mp_limb_t*) flint_stack_alloc(n2);
		}

      bits = -1L*bits;
      if (input1 != input2)
         F_mpz_poly_bit_pack2(int2, int2b, n2, input2, bits, length2, sign2, sign2a);
      F_mpz_poly_bit_pack2(int1, int1b, n1, input1, bits, length1, sign1, sign1a);  
      bits=ABS(bits);
		bits *= 2;
   } 
	
   if (input1 == input2)
	{
		int2 = int1;
		int2b = int1b;
	}
   
   int3 = (mp_limb_t *) flint_stack_alloc(n1 + n2);
   int3b = (mp_limb_t *) flint_stack_alloc(n1 + n2);
   int4 = (mp_limb_t *) flint_stack_alloc(n1 + n2);
          
   mp_limb_t msl = F_mpn_mul(int3, int1, n1, int2, n2);   	
	int3[n1 + n2 - 1] = msl;
   
	msl = F_mpn_mul(int3b, int1b, n1, int2b, n2);	
	int3b[n1 + n2 - 1] = msl;
   
	F_mpn_clear(int4, n1 + n2);
	if ((length1 ^ length2) & 1)
	{
	   mpn_add_n(int4, int3, int3b, n1 + n2);
      mpn_rshift(int4, int4, n1 + n2, 1 + bits/2);
	   mpn_sub_n(int3, int3, int3b, n1 + n2);
      mpn_rshift(int3, int3, n1 + n2, 1);
	} else
   {
	   mpn_sub_n(int4, int3, int3b, n1 + n2);
      mpn_rshift(int4, int4, n1 + n2, 1 + bits/2);
	   mpn_add_n(int3, int3, int3b, n1 + n2);
      mpn_rshift(int3, int3, n1 + n2, 1);
	}
   
   F_mpz_poly_t prod1, prod2;
	F_mpz_poly_init2(prod1, (final_length + 1)/2);
	F_mpz_poly_init2(prod2, (final_length + 1)/2);

	for (ulong i = 0; i < (final_length + 1)/2; i++)
	{   
		_F_mpz_zero(prod1, i);
      _F_mpz_zero(prod2, i);
	}

	for (ulong i = 0; i < final_length; i++)
	{   
		_F_mpz_zero(output, i);
	}

   prod1->length = (final_length + 1)/2;
	prod2->length = final_length/2;
   
	if (bitpack)
   {
      F_mpz_poly_bit_unpack(prod1, int3, n1 + n2, (final_length + 1)/2, bits);  
      F_mpz_poly_bit_unpack(prod2, int4, n1 + n2, final_length/2, bits);  
   } 
   	
	F_mpz_poly_interleave_small(output, prod1, prod2);
   output->length = final_length;

	F_mpz_poly_clear(prod1);
   F_mpz_poly_clear(prod2);
	
   flint_stack_release(); // release int4
   flint_stack_release(); // release int3b
   flint_stack_release(); // release int3
   if (input1 != input2)
	{
		flint_stack_release(); // release int2b
      flint_stack_release(); // release int2
	}
   flint_stack_release(); // release int1b
   flint_stack_release(); // release int1
   
	if ((long) (sign1 ^ sign2) < 0L) F_mpz_poly_neg(output, output);
}

void F_mpz_poly_mul_KS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init(output);
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_KS(output, poly1, poly2);
		else _F_mpz_poly_mul_KS(output, poly2, poly1);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
	{
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_KS(res, poly1, poly2);
		else _F_mpz_poly_mul_KS(res, poly2, poly1);
	}		
}

void F_mpz_poly_mul_KS2(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init(output);
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_KS2(output, poly1, poly2);
		else _F_mpz_poly_mul_KS2(output, poly2, poly1);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
	{
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_KS2(res, poly1, poly2);
		else _F_mpz_poly_mul_KS2(res, poly2, poly1);
	}		
}

