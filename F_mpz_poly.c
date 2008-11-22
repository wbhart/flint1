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
#include "F_mpz.h"
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

void F_mpz_poly_init(F_mpz_poly_t poly)
{
   poly->coeffs = NULL;
   
   poly->alloc = 0;
   poly->length = 0;
}

void F_mpz_poly_init2(F_mpz_poly_t poly, const ulong alloc)
{
   if (alloc) // allocate space for alloc small coeffs
   {
      poly->coeffs = (F_mpz *) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
   }
   else poly->coeffs = NULL;
   
   poly->alloc = alloc;
   poly->length = 0;
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
		F_mpz_poly_truncate(poly, alloc);

		poly->coeffs = (F_mpz *) flint_heap_realloc(poly->coeffs, alloc);
		if (alloc > poly->alloc)
		   F_mpn_clear(poly->coeffs + poly->alloc, alloc - poly->alloc);
	} else // nothing allocated already so do it now
	{
		poly->coeffs = (mp_limb_t*) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
	}
   
   poly->alloc = alloc;  
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
   for (ulong i = 0; i < poly->alloc; i++) // Clean up any mpz_t's
		_F_mpz_demote(poly->coeffs + i);
	if (poly->coeffs) flint_heap_free(poly->coeffs); // clean up ordinary coeffs
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

void F_mpz_poly_set_coeff_si(F_mpz_poly_t poly, ulong n, const long x)
{
   F_mpz_poly_fit_length(poly, n + 1);
   
	if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      for (ulong i = poly->length; i + 1 < n; i++)
         F_mpz_zero(poly->coeffs + i);
      poly->length = n+1;
   }
   
	F_mpz_set_si(poly->coeffs + n, x);
   _F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

void F_mpz_poly_set_coeff_ui(F_mpz_poly_t poly, ulong n, const ulong x)
{
   F_mpz_poly_fit_length(poly, n+1);

   if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      for (long i = poly->length; i + 1 < n; i++)
         F_mpz_zero(poly->coeffs + i); 
      poly->length = n+1;
   }

   F_mpz_set_ui(poly->coeffs + n, x);
   _F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

void F_mpz_poly_set_coeff_mpz(F_mpz_poly_t poly, ulong n, const mpz_t x)
{
   F_mpz_poly_fit_length(poly, n+1);

   if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      for (long i = poly->length; i + 1 < n; i++)
         F_mpz_zero(poly->coeffs + i); 
      poly->length = n+1;
   }

   F_mpz_set_mpz(poly->coeffs + n, x);
	_F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

long F_mpz_poly_get_coeff_si(const F_mpz_poly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
      return 0;
   
	return F_mpz_get_si(poly->coeffs + n);
}

ulong F_mpz_poly_get_coeff_ui(const F_mpz_poly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
      return 0;
   
	return F_mpz_get_ui(poly->coeffs + n);
}

void F_mpz_poly_get_coeff_mpz(mpz_t x, const F_mpz_poly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
	{
		mpz_set_ui(x, 0);
		return;
   }
   
	F_mpz_get_mpz(x, poly->coeffs + n);
	return;
}


/*===============================================================================

	Conversions

================================================================================*/

void mpz_poly_to_F_mpz_poly(F_mpz_poly_t F_poly, const mpz_poly_t m_poly)
{
	F_mpz_poly_fit_length(F_poly, m_poly->length);

	_F_mpz_poly_set_length(F_poly, m_poly->length);
   
	for (ulong i = 0; i < m_poly->length; i++)
		F_mpz_set_mpz(F_poly->coeffs + i, m_poly->coeffs[i]);
}

void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly)
{
	mpz_poly_ensure_alloc(m_poly, F_poly->length);

   m_poly->length = F_poly->length;
   
   for (ulong i = 0; i < F_poly->length; i++)
	   F_mpz_get_mpz(m_poly->coeffs[i], F_poly->coeffs + i);
}

/*===============================================================================

	Assignment/swap

================================================================================*/

void F_mpz_poly_set(F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	if (poly1 != poly2) // aliasing is trivial
	{
		ulong length = poly2->length;
	
	   F_mpz_poly_fit_length(poly1, poly2->length);

		for (ulong i = 0; i < poly2->length; i++)
			F_mpz_set(poly1->coeffs + i, poly2->coeffs + i);
		
		_F_mpz_poly_set_length(poly1, poly2->length);
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
	
	F_mpz * temp_c = poly1->coeffs;
	poly1->coeffs = poly2->coeffs;
	poly2->coeffs = temp_c;

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
		if (!F_mpz_equal(poly1->coeffs + i, poly2->coeffs + i)) return 0;

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
	F_mpz c;

	// search until we find an mpz_t coefficient or one of at least FLINT_BITS - 2 bits
	for (i = 0; i < poly->length; i++) 
	{
		c = poly->coeffs[i];
		if (COEFF_IS_MPZ(c)) break; // found an mpz_t coeff
      if (c < 0L) 
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
			__mpz_struct * mpz_ptr = F_mpz_ptr_mpz(c);
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

	ulong max = 1; // all coefficients have at least one limb
   ulong limbs;
	ulong c;

   // search through mpz coefficients for one of largest size
	for (ulong i = 0; i < poly->length; i++)
	{
		c = poly->coeffs[i];
      if (COEFF_IS_MPZ(c))
		{
			limbs = mpz_size(F_mpz_ptr_mpz(c));
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
         F_mpz_set(res->coeffs + length - i - 1, poly->coeffs + i); // copy over extant coefficients in reverse

      for ( ; i < length; i++) // set other coefficients to zero
         F_mpz_zero(res->coeffs + length - i - 1);

   } else // same polynomial
   {
      for (i = 0; i < length/2; i++)
      {
         // swap extant coefficients
			if (length - i - 1 < res->length) F_mpz_swap(res->coeffs + i, res->coeffs + length - i - 1); 
			else
			{
				F_mpz_set(res->coeffs + length - i - 1, res->coeffs + i); // for other coefficients "swap" with zero
			   F_mpz_zero(res->coeffs + i);
		   }
		}
      // if length is odd we missed a coefficient in swapping pairs, it may need to be set to zero
		if ((length & 1) && (i >= poly->length)) F_mpz_zero(res->coeffs + i); 
   }
	
	_F_mpz_poly_set_length(res, length);
   _F_mpz_poly_normalise(res); // new leading coeff, which was trailing coeff, may now be zero
}

/*===============================================================================

	Negation

================================================================================*/

void F_mpz_poly_neg(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
	F_mpz_poly_fit_length(res, poly->length);
	
	for (ulong i = 0; i < poly->length; i++)
		F_mpz_neg(res->coeffs + i, poly->coeffs + i);

	_F_mpz_poly_set_length(res, poly->length);
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
      F_mpz_add(res->coeffs + i, poly1->coeffs + i, poly2->coeffs + i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (ulong i = shorter; i < poly1->length; i++)
         F_mpz_set(res->coeffs + i, poly1->coeffs + i);

   if (poly2 != res) // copy any remaining coefficients from poly2
      for (ulong i = shorter; i < poly2->length; i++)
         F_mpz_set(res->coeffs + i, poly2->coeffs + i);
   
   if (poly1->length == poly2->length)
   {
      _F_mpz_poly_set_length(res, poly1->length);
      _F_mpz_poly_normalise(res); // there may have been cancellation
   } else
      _F_mpz_poly_set_length(res, longer);
}

void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);
	ulong shorter = FLINT_MIN(poly1->length, poly2->length);

   F_mpz_poly_fit_length(res, longer);
	
   for (ulong i = 0; i < shorter; i++) // add up to the length of the shorter poly
      F_mpz_sub(res->coeffs + i, poly1->coeffs + i, poly2->coeffs + i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (ulong i = shorter; i < poly1->length; i++)
         F_mpz_set(res->coeffs + i, poly1->coeffs + i);

   // careful, it is *always* necessary to negate coeffs from poly2, even if this is already res
	for (ulong i = shorter; i < poly2->length; i++) 
      F_mpz_neg(res->coeffs + i, poly2->coeffs + i);

   if (poly1->length == poly2->length)
   {
      _F_mpz_poly_set_length(res, poly1->length);
      _F_mpz_poly_normalise(res); // there may have been cancellation
   } else
      _F_mpz_poly_set_length(res, longer);
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
		_F_mpz_poly_set_length(res, 0);
		return;
	}
	
	F_mpz_poly_fit_length(res, poly->length + n);
	
	// copy in reverse order to avoid writing over unshifted coeffs
	for (long i = poly->length - 1; i >= 0; i--) 
		F_mpz_set(res->coeffs + i + n, poly->coeffs + i);

   // insert n zeroes
	for (ulong i = 0; i < n; i++) F_mpz_zero(res->coeffs + i);
   
   _F_mpz_poly_set_length(res, poly->length + n);
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
	for (ulong i = 0; i < poly->length - n; i++) 
		F_mpz_set(res->coeffs + i, poly->coeffs + i + n);
	
	_F_mpz_poly_set_length(res, poly->length - n);
}

/*===============================================================================

	Interleave

================================================================================*/

/*void F_mpz_poly_interleave(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2)
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
}*/

/*===============================================================================

	Scalar multiplication

================================================================================*/

void F_mpz_poly_scalar_mul_ui(F_mpz_poly_t poly1, F_mpz_poly_t poly2, ulong x)
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
	
	F_mpz_poly_fit_length(poly1, poly2->length);
	
	for (ulong i = 0; i < poly2->length; i++) 
		F_mpz_mul_ui(poly1->coeffs + i, poly2->coeffs + i, x);

	_F_mpz_poly_set_length(poly1, poly2->length);
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
	
	for (ulong i = 0; i < poly2->length; i++) 
		F_mpz_mul_si(poly1->coeffs + i, poly2->coeffs + i, x);

	_F_mpz_poly_set_length(poly1, poly2->length);
}

void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, F_mpz_poly_t poly2, F_mpz_t x)
{
	// either scalar or input poly is zero
	if ((*x == 0) || (poly2->length == 0)) 
	{
	   F_mpz_poly_zero(poly1);
		return;
	}
	
	// special case, muliply by 1
	if (*x == 1L)
	{
	   F_mpz_poly_set(poly1, poly2);
		return;
	}
	
	// special case, muliply by -1
	if (*x == -1L)
	{
	   F_mpz_poly_neg(poly1, poly2);
		return;
	}
	
	F_mpz_poly_fit_length(poly1, poly2->length);
	
	for (ulong i = 0; i < poly2->length; i++) 
		F_mpz_mul2(poly1->coeffs + i, poly2->coeffs + i, x);

	_F_mpz_poly_set_length(poly1, poly2->length);
}

/*===============================================================================

	Classical multiplication

================================================================================*/

void _F_mpz_poly_sqr_classical(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
   F_mpz_poly_fit_length(res, 2*poly->length - 1);
	_F_mpz_poly_set_length(res, 2*poly->length - 1);

   for (ulong i = 0; i < res->length; i++)
      F_mpz_zero(res->coeffs + i);
   
   // off-diagonal products
   for (ulong i = 1; i < poly->length; i++)
	{
		F_mpz c = poly->coeffs[i];
	   if (c)
		{
			if (!COEFF_IS_MPZ(c))
			{
				if (c < 0L) 
					for (ulong j = 0; j < i; j++)
                  F_mpz_submul_ui(res->coeffs + i + j, poly->coeffs + j, -c);
				else
               for (ulong j = 0; j < i; j++)
                  F_mpz_addmul_ui(res->coeffs + i + j, poly->coeffs + j, c);
			} else
		      for (ulong j = 0; j < i; j++)
               F_mpz_addmul(res->coeffs + i+j, poly->coeffs + i, poly->coeffs + j);
		}
	}

   // double the off-diagonal products
   for (ulong i = 1; i < res->length - 1; i++)
      F_mpz_add(res->coeffs + i, res->coeffs + i, res->coeffs + i);
      
   // add in diagonal products
   for (ulong i = 0; i < poly->length; i++)
      F_mpz_addmul(res->coeffs + 2*i, poly->coeffs + i, poly->coeffs + i);
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
      F_mpz_poly_init2(temp, length);
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
      F_mpz_mul2(res->coeffs, poly1->coeffs, poly2->coeffs);      
   } else // Ordinary case
   {
      long i, j;
      
      // Set res[i] = poly1[i]*poly2[0] 
      if (poly2->coeffs[0])
			for (i = 0; i < len1; i++)
            F_mpz_mul2(res->coeffs + i, poly1->coeffs + i, poly2->coeffs);
		else 
			for (i = 0; i < len1; i++)
            F_mpz_zero(res->coeffs + i);

      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      if (poly1->coeffs[len1 - 1])
		   for (i = 1; i < len2; i++)
            F_mpz_mul2(res->coeffs + i + len1 - 1, poly1->coeffs + len1 - 1, poly2->coeffs + i);  
		else 
         for (i = 1; i < len2; i++)
            F_mpz_zero(res->coeffs + i + len1 - 1);
      
      // out[i+j] += in1[i]*in2[j] 
      for (i = 0; i < len1 - 1; i++)
      {      
         F_mpz c = poly1->coeffs[i];
			if (c)
			{
				if (!COEFF_IS_MPZ(c))
				{
					if (c < 0L) 
						for (j = 1; j < len2; j++)
                     F_mpz_submul_ui(res->coeffs + i + j, poly2->coeffs + j, -c);
					else
                  for (j = 1; j < len2; j++)
						   F_mpz_addmul_ui(res->coeffs + i + j, poly2->coeffs + j, c);
				} else
					for (j = 1; j < len2; j++)
                  F_mpz_addmul(res->coeffs + i + j, poly1->coeffs + i, poly2->coeffs + j);
			}
      }
   } 
   
   _F_mpz_poly_set_length(res, len1 + len2 - 1);
}

void F_mpz_poly_mul_classical(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if (poly1 == poly2) F_mpz_poly_sqr_classical(res, poly1); // Aliased inputs
 
	if ((poly1 == res) || (poly2 == res)) // aliased input and output
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init2(output, poly1->length + poly2->length - 1);
		_F_mpz_poly_mul_classical(output, poly1, poly2);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
		_F_mpz_poly_mul_classical(res, poly1, poly2);
}

/*===============================================================================

	Karatsuba multiplication

================================================================================*/

void _F_mpz_poly_mul_kara_recursive(F_mpz * out, F_mpz * in1, 
											ulong len1, F_mpz * in2, ulong len2, 
											F_mpz * scratch, ulong skip, ulong crossover)
{
   // Base cases:
   
   if (len1 == 1)
   {
      // special case, just scalar multiplication
      for (ulong i = 0; i < len2; i++)
         F_mpz_mul2(out + i*skip, in1, in2 + i*skip);
      return;
   }
   
   if (len1 * len2 < crossover)
   {
      // switch to naive multiplication

		if (in2)
			for (ulong i = 0; i < len1; i++)
            F_mpz_mul2(out + i*skip, in1 + i*skip, in2);
		else 
			for (ulong i = 0; i < len1; i++)
            F_mpz_zero(out + i*skip);

      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      const ulong term = (len1 - 1)*skip;
		if (in1[(len1 - 1)*skip])
		   for (ulong i = 1; i < len2; i++)
            F_mpz_mul2(out + (i + len1 - 1)*skip, in1 + term, in2 + i*skip);  
	 	else 
         for (ulong i = 1; i < len2; i++)
            F_mpz_zero(out + (i + len1 - 1)*skip);
      
      for (ulong i = 0; i < len1 - 1; i++)
		{
			F_mpz c = in1[i*skip];
			const ulong term = i*skip;
			if (!COEFF_IS_MPZ(c))
			{
				if (c < 0L) 
					for (ulong j = 1; j < len2; j++)
                  F_mpz_submul_ui(out + (i+j)*skip, in2 + j*skip, -c);
				else
               for (ulong j = 1; j < len2; j++)
                  F_mpz_addmul_ui(out + (i+j)*skip, in2 + j*skip, c);
			} else
				for (ulong j = 1; j < len2; j++)
               F_mpz_addmul(out + (i+j)*skip, in1 + term, in2 + j*skip);
		}    

		return;
   }

   ulong i, j;
	
	// Recursive case:

   // Let in1 = A1(x^2) + x*B1(x^2) + x^(2*floor(len1/2))*C1,
   // where A1, B1 have length floor(len1/2),
   // and C1 is the leading term of in1 if len1 is odd

   // Similarly for in2 = A2(x^2) + x*B2(x^2) + x^(2*floor(len2/2))*C2.
   
   // Put A1 + B1 into even slots of scratch space
   // (uses len1/2 scratch slots)
   for (i = 0; i < len1/2; i++)
      F_mpz_add(scratch + 2*i*skip, in1 + 2*i*skip, in1 + 2*i*skip + skip);

   // Put A2 + B2 into remaining even slots of scratch space
   // (uses len2/2 slots of scratch)
   for (j = 0; j < len2/2; j++)
      F_mpz_add(scratch + 2*(i+j)*skip, in2 + 2*j*skip, in2 + 2*j*skip + skip);

   // The following three recursive calls all use the odd slots of the current
   // scratch array as the next layer's scratch space
   
   // Put product (A1+B1)*(A2+B2) into odd slots of output array
   _F_mpz_poly_mul_kara_recursive(out + skip, scratch, len1/2, scratch + 2*i*skip, len2/2,
                                scratch + skip, 2*skip, crossover);

   // Put product x^2*(B1*B2) into even slots of output array
   // (except first slot, which is an implied zero)
   _F_mpz_poly_mul_kara_recursive(out + 2*skip, in1 + skip, len1/2, in2 + skip,
                                len2/2, scratch + skip, 2*skip, crossover);

   // Put product A1*A2 into even slots of scratch space
   _F_mpz_poly_mul_kara_recursive(scratch, in1, len1/2, in2, len2/2,
                                scratch + skip, 2*skip, crossover);
                            
   // Subtract A1*A2 and B1*B2 from (A1+B1)*(A2+B2) to get (A1*B2 + A2*B1)
   // in odd slots of output
   for (ulong i = 0; i < len1/2 + len2/2 - 1; i++)
   {
      F_mpz_sub(out + 2*i*skip + skip, out + 2*i*skip + skip, out + 2*(i+1)*skip);
      F_mpz_sub(out + 2*i*skip + skip, out + 2*i*skip + skip, scratch + 2*i*skip);
   }
      
   // Add A1*A2 to x^2*(B1*B2) into even slots of output
   F_mpz_set(out, scratch);
   for (ulong i = 1; i < len1/2 + len2/2 - 1; i++)
      F_mpz_add(out + 2*i*skip, out + 2*i*skip, scratch + 2*i*skip);
   
   // Now we have the product (A1(x^2) + x*B1(x^2)) * (A2(x^2) + x*B2(x^2))
   // in the output array. Still need to handle C1 and C2 terms.
   
   if (len1 & 1)
   {
      const ulong term1 = skip*(len1-1);
	   const ulong term2 = skip*(len2-1);
	   if (len2 & 1)
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         for (ulong i = 0; i < len2-2; i++)
            F_mpz_addmul(out + (i+len1-1)*skip, in1 + term1, in2 + i*skip);
         F_mpz_mul2(out + (len1+len2-3)*skip, in1 + term1, in2 + (len2-2)*skip);

         // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
         for (ulong i = 0; i < len1-1; i++)
            F_mpz_addmul(out + (i+len2-1)*skip, in2 + term2, in1 + i*skip);
            
         // final C1*C2 term
         F_mpz_mul2(out + (len1+len2-2)*skip, in1 + term1, in2 + term2);
      }
      else
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         for (ulong i = 0; i < len2-1; i++)
            F_mpz_addmul(out + (i+len1-1)*skip, in1 + term1, in2 + i*skip);
         F_mpz_mul2(out + (len1+len2-2)*skip, in1 + term1, in2 + term2);
      }
   }
   else if (len2 & 1)
   {
      const ulong term1 = skip*(len1-1);
	   const ulong term2 = skip*(len2-1);
	   // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
      for (ulong i = 0; i < len1-1; i++)
         F_mpz_addmul(out + (i+len2-1)*skip, in2 + term2, in1 + i*skip);
      F_mpz_mul2(out + (len1+len2-2)*skip, in2 + term2, in1 + term1);
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
	F_mpz_poly_init2(scratch, length + 1);
	scratch->length = length + 1;
	
	// look up crossover parameter (i.e. when to switch from classical to
   // karatsuba multiplication) based on coefficient size
   ulong bits1 = 64;//FLINT_ABS(F_mpz_poly_max_bits(poly1));
	ulong bits2 = 64;//FLINT_ABS(F_mpz_poly_max_bits(poly2));
	ulong log_length = 0;
	while ((1L<<log_length) < poly1->length) log_length++;
	ulong limbs = (bits1 + bits2 + log_length - 1)/FLINT_BITS + 1;

	ulong crossover = _F_mpz_poly_mul_karatsuba_crossover(limbs);
   
   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      F_mpz_poly_t temp;
      F_mpz_poly_init2(temp, length);
		
      _F_mpz_poly_mul_kara_recursive(temp->coeffs, poly1->coeffs, poly1->length,
            poly2->coeffs, poly2->length, scratch->coeffs, 1, crossover);

      temp->length = length;
		
		F_mpz_poly_swap(temp, res);
      F_mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      _F_mpz_poly_mul_kara_recursive(res->coeffs, poly1->coeffs, poly1->length,
            poly2->coeffs, poly2->length, scratch->coeffs, 1, crossover);

		_F_mpz_poly_set_length(res, length);
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
      // polys are identical, call specialised squaring routine
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
/*#define NEXT_COEFF(xxxarr, xxxborr, xxxcoeff, xxxmask, xxxneg) \
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
}*/
 
/*
  Bit packing used with David Harvey's KS2 algorithm

*/
/*void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
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
}*/


/*===============================================================================

	Kronecker Segmentation multiplication

================================================================================*/

/*void _F_mpz_poly_mul_KS(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2)
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
}*/

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
/*void _F_mpz_poly_mul_KS2(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2)
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
}*/
