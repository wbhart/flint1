/*============================================================================

    F_mpz_poly.c: Polynomials over Z (FLINT 2.0 polynomials)

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
   ulong r = res->coeffs[coeff3];
	__mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff3);
	if (!COEFF_IS_MPZ(r)) // res is small
		mpz_set_si(mpz_ptr, r);

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
   ulong r = res->coeffs[coeff3];
	__mpz_struct * mpz_ptr = _F_mpz_promote(res, coeff3);
	if (!COEFF_IS_MPZ(r)) // res is small
		mpz_set_si(mpz_ptr, r);

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

long F_mpz_poly_max_bits(F_mpz_poly_t poly)
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

ulong F_mpz_poly_max_limbs(F_mpz_poly_t poly)
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

	Multiplication

================================================================================*/

void _F_mpz_poly_sqr_classical(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
   F_mpz_poly_fit_length(res, 2*poly->length - 1);
	res->length = 2*poly->length - 1;

   for (ulong i = 0; i < res->length; i++)
      _F_mpz_zero(res, i);
   
   // off-diagonal products
   for (ulong i = 1; i < poly->length; i++)
      for (ulong j = 0; j < i; j++)
         _F_mpz_addmul(res, i+j, poly, i, poly, j);
         
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
