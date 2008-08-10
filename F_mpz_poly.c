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
	if (poly->mpz_length == poly->mpz_alloc) // time to allocate some more mpz_t's
	{
	   if (poly->mpz_alloc) 
		   poly->mpz_coeffs = (__mpz_struct*) flint_heap_realloc_bytes(poly->mpz_coeffs, (poly->mpz_alloc + MPZ_BLOCK)*sizeof(__mpz_struct));
		else
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
   for (ulong i = 0; i < poly->mpz_alloc; i++)
	   mpz_clear(poly->mpz_coeffs + i);

	flint_heap_free(poly->mpz_coeffs);
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
   if (alloc)
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
   if (!alloc) // clear up
   {
         F_mpz_poly_clear(poly);
			return;
   }  
   
	if (poly->alloc) 
	{
		poly->coeffs = (mp_limb_t*) flint_heap_realloc(poly->coeffs, alloc);
		if (alloc > poly->alloc)
		   F_mpn_clear(poly->coeffs + poly->alloc, alloc - poly->alloc);
	}
   else 
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

   if (alloc < 2*poly->alloc) alloc = 2*poly->alloc;
   
   F_mpz_poly_realloc(poly, alloc);
}

void F_mpz_poly_clear(F_mpz_poly_t poly)
{
   if (poly->coeffs) flint_heap_free(poly->coeffs);
   if (poly->mpz_coeffs) _F_mpz_poly_mpz_coeffs_clear(poly);
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
   __mpz_struct * mpz_ptr = poly->mpz_coeffs + FMPZ_MPZ(poly->coeffs[coeff]);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L)
	{
		poly->coeffs[coeff] = 0;
	} else if (size == 1L)
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval < FMPZ_MPZ_MASK) poly->coeffs[coeff] = uval;
	} else if (size == -1L)
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval < FMPZ_MPZ_MASK) poly->coeffs[coeff] = uval + FMPZ_SIGN_MASK;
	}
}

void _F_mpz_set_si(F_mpz_poly_t poly, ulong coeff, const long val)
{
   mp_limb_t c = poly->coeffs[coeff];
	ulong uval = FLINT_ABS(val);
	
	if (uval > FMPZ_UVAL_MASK) // val is large
	{
		__mpz_struct * mpz_coeff;
		if (!FMPZ_IS_MPZ(c)) // coeff is small
		{
			_F_mpz_promote(poly, coeff);
		   mpz_coeff = poly->mpz_coeffs + poly->mpz_length - 1;
		} else mpz_coeff = poly->mpz_coeffs + FMPZ_MPZ(c);
		mpz_set_si(mpz_coeff, val);
	} else // val is small
	{
		if (val < 0L) poly->coeffs[coeff] = FMPZ_SIGN_MASK | uval;
		else poly->coeffs[coeff] = uval;
	}
}

void _F_mpz_get_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong coeff)
{
	mp_limb_t c = poly->coeffs[coeff];

	if (!FMPZ_IS_MPZ(c)) // coeff is small
	{
		if ((long) c < 0L) mpz_set_si(x, -FMPZ_UVAL(c));
		else mpz_set_ui(x, c);
	} else mpz_set(x, poly->mpz_coeffs + FMPZ_MPZ(c));
}

void _F_mpz_set_mpz(F_mpz_poly_t poly, ulong coeff, const mpz_t x)
{
   long size = x->_mp_size;
	
	if (size == 0L)
	{
		poly->coeffs[coeff] = 0;
	} else if (size == 1L)
	{
	   ulong uval = mpz_get_ui(x);
		if (uval < FMPZ_MPZ_MASK) poly->coeffs[coeff] = uval;
		else 
		{
			mp_limb_t c = poly->coeffs[coeff];
			__mpz_struct * mpz_ptr;
			if (!FMPZ_IS_MPZ(c)) // coeff is small
			{
				_F_mpz_promote(poly, coeff);
            mpz_ptr = poly->mpz_coeffs + poly->mpz_length - 1;
			} else // coeff is large
			   mpz_ptr = poly->mpz_coeffs + FMPZ_MPZ(c);
			mpz_set_ui(mpz_ptr, uval);
		}
	} else if (size == -1L)
   {
	   ulong uval = mpz_get_ui(x);
		if (uval < FMPZ_MPZ_MASK) poly->coeffs[coeff] = uval + FMPZ_SIGN_MASK;
		else 
		{
			mp_limb_t c = poly->coeffs[coeff];
			__mpz_struct * mpz_ptr;
			if (!FMPZ_IS_MPZ(c)) // coeff is small
			{
				_F_mpz_promote(poly, coeff);
            mpz_ptr = poly->mpz_coeffs + poly->mpz_length - 1;
			} else // coeff is large
			   mpz_ptr = poly->mpz_coeffs + FMPZ_MPZ(c);
			mpz_set_ui(mpz_ptr, uval);
			mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else 
	{
		mp_limb_t c = poly->coeffs[coeff];
		__mpz_struct * mpz_ptr;
		if (!FMPZ_IS_MPZ(c)) // coeff is small
		{
			_F_mpz_promote(poly, coeff);
         mpz_ptr = poly->mpz_coeffs + poly->mpz_length - 1;
		} else // coeff is large
			mpz_ptr = poly->mpz_coeffs + FMPZ_MPZ(c);
		mpz_set(mpz_ptr, x);
		_F_mpz_demote_val(poly, coeff);
	}			
}

void _F_mpz_add(F_mpz_poly_t res, ulong coeff3, const F_mpz_poly_t poly1, const ulong coeff1, 
					                                 const F_mpz_poly_t poly2, const ulong coeff2) 
{
	mp_limb_t c1 = poly1->coeffs[coeff1];
	mp_limb_t c2 = poly2->coeffs[coeff2];
	
	if (!FMPZ_IS_MPZ(c1))
	{
	   if (!FMPZ_IS_MPZ(c2)) // both coefficients are small
		{
		   long sum; 
	      if ((long) (c1 ^ c2) < 0L) // coefficients have opposite signs
			{
				if ((long) c1 < 0L) sum = c2 - FMPZ_UVAL(c1);
				else sum = c1 - FMPZ_UVAL(c2);
			} else // coefficients have same sign
			{
            if ((long) c1 < 0L) sum = -FMPZ_UVAL(c1) - FMPZ_UVAL(c2);
				else sum = c1 + c2;
			}
			_F_mpz_set_si(res, coeff3, sum);
		} else // c1 is small, c2 is large
		{
         mp_limb_t c3 = res->coeffs[coeff3];
	      __mpz_struct * mpz3;
			if (!FMPZ_IS_MPZ(c3)) // c3 is small
			{
				_F_mpz_promote(res, coeff3);
				mpz3 = res->mpz_coeffs + res->mpz_length - 1;
			} else // c3 is large
			   mpz3 = res->mpz_coeffs + FMPZ_MPZ(c3);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + FMPZ_MPZ(c2);
			if ((long) c1 < 0L) // c1 is negative
			{
				mpz_sub_ui(mpz3, mpz2, FMPZ_UVAL(c1));	
				mpz_neg(mpz3, mpz3);
			} else mpz_add_ui(mpz3, mpz2, c1); // c1 is positive
			_F_mpz_demote_val(res, coeff3);
		}
	} else
	{
		if (!FMPZ_MPZ(c2)) // c1 is large, c2 is small
		{
         mp_limb_t c3 = res->coeffs[coeff3];
	      __mpz_struct * mpz3;
			if (!FMPZ_IS_MPZ(c3)) // c3 is small
			{
				_F_mpz_promote(res, coeff3);
				mpz3 = res->mpz_coeffs + res->mpz_length - 1;
			} else // c3 is large
			   mpz3 = res->mpz_coeffs + FMPZ_MPZ(c3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + FMPZ_MPZ(c1);
			if ((long) c2 < 0L) // c2 is negative
			{
				mpz_sub_ui(mpz3, mpz1, FMPZ_UVAL(c2));	
			} else mpz_add_ui(mpz3, mpz1, c2); // c2 is positive
			_F_mpz_demote_val(res, coeff3);
		} else // c1 and c2 are large
		{
         mp_limb_t c3 = res->coeffs[coeff3];
	      __mpz_struct * mpz3;
			if (!FMPZ_IS_MPZ(c3)) // c3 is small
			{
				_F_mpz_promote(res, coeff3);
				mpz3 = res->mpz_coeffs + res->mpz_length - 1;
			} else // c3 is large
			   mpz3 = res->mpz_coeffs + FMPZ_MPZ(c3);
			__mpz_struct * mpz1 = poly1->mpz_coeffs + FMPZ_MPZ(c1);
			__mpz_struct * mpz2 = poly2->mpz_coeffs + FMPZ_MPZ(c2);
			mpz_add(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(res, coeff3);
		}
	}
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

