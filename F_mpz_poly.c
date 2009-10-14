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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>

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

       Input/output 

================================================================================*/

int F_mpz_poly_from_string(F_mpz_poly_t poly, const char* s)
{
   int ok;
   
   mpz_poly_t p;
   mpz_poly_init(p);
   ok = mpz_poly_from_string(p, s);
   if (ok)
   {
      mpz_poly_to_F_mpz_poly(poly, p);
   }
   mpz_poly_clear(p);
   
   return ok;
}

char* F_mpz_poly_to_string(const F_mpz_poly_t poly)
{
   char* buf;
   mpz_poly_t m_poly;
   mpz_poly_init(m_poly);
   F_mpz_poly_to_mpz_poly(m_poly, poly);
   buf = mpz_poly_to_string(m_poly);
   mpz_poly_clear(m_poly);
   return buf;
}

char* F_mpz_poly_to_string_pretty(const F_mpz_poly_t poly, const char * x)
{
   char* buf;
   mpz_poly_t m_poly;
   mpz_poly_init(m_poly);
   F_mpz_poly_to_mpz_poly(m_poly, poly);
   buf = mpz_poly_to_string_pretty(m_poly, x);
   mpz_poly_clear(m_poly);
   return buf;
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
		if (!F_mpz_equal(poly1->coeffs + i, poly2->coeffs + i)) 
		   return 0;

	return 1;
}

/*===============================================================================

	Coefficient sizes

================================================================================*/

long F_mpz_poly_max_bits1(const F_mpz_poly_t poly)
{
	int sign = 0;
	F_mpz or = 0;
   F_mpz c;
	F_mpz mask = (1L<<(FLINT_BITS-3)); // see if the largest bit possible in a small is set
   ulong i;

	// search until we find an mpz_t coefficient or one of at least FLINT_BITS - 2 bits
	for (i = 0; i < poly->length; i++) 
	{
		c = poly->coeffs[i];
		if (c < 0L) 
		{
			sign = 1;
			or = or | -c;
		} else
		{
			if (COEFF_IS_MPZ(c)) return 0; // found an mpz_t coeff
         or = or | c;
		}
		if (or & mask) break; // found a coeff with FLINT_BIT - 2 bits
	}

	if (!sign) // if no negative coefficient yet keep searching
	   for ( ; i < poly->length; i++)  
	      if (poly->coeffs[i] < 0L) return -(FLINT_BITS - 2); // only happens if we hit
	                                                       // break in previous loop

	if (sign) return -FLINT_BIT_COUNT(or); // return bits n or -n if negative found
	else return FLINT_BIT_COUNT(or);
}

long F_mpz_poly_max_bits(const F_mpz_poly_t poly)
{
	int sign = 0;
	ulong max = 0;
   ulong bits = 0;
   ulong max_limbs = 1;
	ulong size;
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
		if (bits > max) 
		{
			max = bits;
		   if (max == FLINT_BITS - 2) break; // coeff is at least FLINT_BITS - 2 bits
		}
	}

   // search through mpz coefficients for largest size in bits
	
	for ( ; i < poly->length; i++)
	{
		c = poly->coeffs[i];
      if (COEFF_IS_MPZ(c))
		{
			__mpz_struct * mpz_ptr = F_mpz_ptr_mpz(c);
			if (mpz_sgn(mpz_ptr) < 0) sign = 1;
			size = mpz_size(mpz_ptr);
			if (size > max_limbs)
			{
			   max_limbs = size;
				mp_limb_t * data = mpz_ptr->_mp_d;
			   bits = FLINT_BIT_COUNT(data[max_limbs - 1]);
				max = bits;
			} else if (size == max_limbs)
			{
				mp_limb_t * data = mpz_ptr->_mp_d;
			   bits = FLINT_BIT_COUNT(data[max_limbs - 1]);
			   if (bits > max) max = bits;
			}
		} else if ((long) c < 0L) sign = 1; // still need to check the sign of small coefficients
	}
	

	if (sign) return -(max + FLINT_BITS*(max_limbs - 1));
	else return max + FLINT_BITS*(max_limbs - 1);
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

void _F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	ulong longer = FLINT_MAX(poly1->length, poly2->length);
	ulong shorter = FLINT_MIN(poly1->length, poly2->length);

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

void F_mpz_poly_add(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);

	F_mpz_poly_fit_length(res, longer);
	
	_F_mpz_poly_add(res, poly1, poly2);
}

void _F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);
	ulong shorter = FLINT_MIN(poly1->length, poly2->length);

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

void F_mpz_poly_sub(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);

	F_mpz_poly_fit_length(res, longer);
   
	_F_mpz_poly_sub(res, poly1, poly2);
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

void F_mpz_poly_interleave(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2)
{
	ulong length = FLINT_MAX(2*poly1->length - 1, 2*poly2->length);
	F_mpz_poly_fit_length(res, length + 1); // extra one in case length is odd
	                                        // it will be set to zero

	ulong shorter = FLINT_MIN(poly1->length, poly2->length);
   ulong i;

	for (i = 0; i < shorter; i++)
	{
		F_mpz_set(res->coeffs + 2*i, poly1->coeffs + i);
		F_mpz_set(res->coeffs + 2*i+1, poly2->coeffs + i);
	}

   for ( ; i < poly1->length; i++)
	{
		F_mpz_set(res->coeffs + 2*i, poly1->coeffs + i);
		F_mpz_zero(res->coeffs + 2*i+1);
	}

	for ( ; i < poly2->length; i++)
	{
		F_mpz_zero(res->coeffs + 2*i);
		F_mpz_set(res->coeffs + 2*i+1, poly2->coeffs + i);
	}

	_F_mpz_poly_set_length(res, length);
}

void F_mpz_poly_interleave_small(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2)
{
	ulong length = FLINT_MAX(2*poly1->length - 1, 2*poly2->length);
	F_mpz_poly_fit_length(res, length + 1); // extra one in case length is odd
	                                        // it will be set to zero

	ulong shorter = FLINT_MIN(poly1->length, poly2->length);
   ulong i;

	for (i = 0; i < shorter; i++)
	{
		_F_mpz_demote(res->coeffs + 2*i);
		_F_mpz_demote(res->coeffs + 2*i + 1);
		res->coeffs[2*i] =  poly1->coeffs[i];
		res->coeffs[2*i+1] =  poly2->coeffs[i];
	}

   for ( ; i < poly1->length; i++)
	{
		_F_mpz_demote(res->coeffs + 2*i);
		_F_mpz_demote(res->coeffs + 2*i + 1);
		res->coeffs[2*i] =  poly1->coeffs[i];
		res->coeffs[2*i+1] = 0;
	}

	for ( ; i < poly2->length; i++)
	{
		_F_mpz_demote(res->coeffs + 2*i);
		_F_mpz_demote(res->coeffs + 2*i + 1);
		res->coeffs[2*i] = 0;
		res->coeffs[2*i+1] =  poly2->coeffs[i];
	}

	_F_mpz_poly_set_length(res, length);
}

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

   if ((len1 == 0) || (len2 == 0))
   {
      res->length = 0;
      return;
   }
   
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

   F_mpz_poly_fit_length(res, poly1->length + poly2->length - 1);
   
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

void _F_mpz_poly_mul_kara_odd_even_recursive(F_mpz * out, F_mpz * in1, 
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
   _F_mpz_poly_mul_kara_odd_even_recursive(out + skip, scratch, len1/2, scratch + 2*i*skip, len2/2,
                                scratch + skip, 2*skip, crossover);

   // Put product x^2*(B1*B2) into even slots of output array
   // (except first slot, which is an implied zero)
   _F_mpz_poly_mul_kara_odd_even_recursive(out + 2*skip, in1 + skip, len1/2, in2 + skip,
                                len2/2, scratch + skip, 2*skip, crossover);

   // Put product A1*A2 into even slots of scratch space
   _F_mpz_poly_mul_kara_odd_even_recursive(scratch, in1, len1/2, in2, len2/2,
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

void _F_mpz_poly_mul_kara_recursive(F_mpz_poly_t res, const F_mpz_poly_t a, 
										  const F_mpz_poly_t b, F_mpz_poly_t scratch, const ulong crossover)
{
   if ((crossover < 4) && (a->length == 2 && b->length == 2)) // special case for 2 x 2 multiplication
	{
      F_mpz_mul2(res->coeffs, a->coeffs, b->coeffs); // r[0] = a[0] * b[0]
      F_mpz_add(scratch->coeffs, a->coeffs, a->coeffs + 1); // s[0] = a[0] + a[1];
      F_mpz_mul2(res->coeffs + 2, a->coeffs + 1, b->coeffs + 1); // r[2] = a[1] * b[1];
      F_mpz_add(scratch->coeffs + 1, b->coeffs, b->coeffs + 1); // s[1] = b[0] + b[1];
      F_mpz_mul2(res->coeffs + 1, scratch->coeffs, scratch->coeffs + 1); // r[1] = s[0] * s[1] 
      F_mpz_sub(res->coeffs + 1, res->coeffs + 1, res->coeffs); // r[1] -= r[0]
      F_mpz_sub(res->coeffs + 1, res->coeffs + 1, res->coeffs + 2); // r[1] -= r[2]
      _F_mpz_poly_set_length(res, 3);
      
      return;
   }
   
   if ((a->length + b->length <= crossover) ||  (a->length <= 1) || (b->length <= 1)) // too short for karatsuba
   {
      _F_mpz_poly_mul_classical(res, a, b);
      
      return;
   }  
    
	/*
      (a1 + a2*x)*(b1 + b2*x) = a1*b1 + a2*b2*x^2 + (a1 + a2)*(b1 + b2)*x - a1*b1*x - a2*b2*x;
   */

	F_mpz_poly_t temp, a1, a2, b1, b2;
        
	ulong n = a->length; // common length of a and b
	ulong n1 = (n + 1)/2; // length of bottom half of a (top half may be shorter by 1)

   _F_mpz_poly_attach_truncate(a1, a, n1); // attach polys to top and bottom halves of a and b
	_F_mpz_poly_attach_truncate(b1, b, n1); 
	_F_mpz_poly_attach_shift(a2, a, n1); 
   _F_mpz_poly_attach_shift(b2, b, n1); 
     
   /* 
       From 0 for a1->length + b1->length - 1 and from 2 * n1 for 2 * a2->length - 1
       will be written directly to, so we only need to clean the coefficients in between
   */
   
	ulong start = a1->length + b1->length - 1;
	if (!a1->length || !b1->length) start = 0;
	for (ulong i = start; i < 2*n1; i++)
      F_mpz_zero(res->coeffs + i);
		
   F_mpz_poly_t asum, bsum, prodsum, scratch2;
     
	asum->coeffs = scratch->coeffs;
	asum->length = 0;
   bsum->coeffs = scratch->coeffs + n1;
   bsum->length = 0;
   prodsum->coeffs = scratch->coeffs + 2*n1;
   prodsum->length = 0;
     
   // res_lo = a1 * b1
	if (a1->length && b1->length) 
	{
		if (a1->length < b1->length) 
		   _F_mpz_poly_mul_kara_odd_even_recursive(res->coeffs, a1->coeffs, a1->length, b1->coeffs, b1->length, scratch->coeffs, 1, crossover);
      else if (a1->length > b1->length)
		   _F_mpz_poly_mul_kara_odd_even_recursive(res->coeffs, b1->coeffs, b1->length, a1->coeffs, a1->length, scratch->coeffs, 1, crossover);
	   else
		   _F_mpz_poly_mul_kara_recursive(res, a1, b1, scratch, crossover);
	   _F_mpz_poly_set_length(res, a1->length + b1->length - 1); // need to set length because odd/even kara can't do this
	} else _F_mpz_poly_set_length(res, 0);

	// res_hi = a2 * b2 (lengths are guaranteed to be the same)
	temp->coeffs = res->coeffs + 2*n1;
	temp->length = 0;
   _F_mpz_poly_mul_kara_recursive(temp, a2, b2, scratch, crossover); 

	_F_mpz_poly_add(asum, a1, a2); // asum = a1 + a2
   _F_mpz_poly_add(bsum, b1, b2); // bsum = b1 + b2
  
   // prodsum = asum * bsum   
	scratch2->coeffs = scratch->coeffs + 4*n1 - 1;  
	scratch2->length = 0;
	
   if (asum->length && bsum->length) 
	{
		if (asum->length < bsum->length) 
		   _F_mpz_poly_mul_kara_odd_even_recursive(prodsum->coeffs, asum->coeffs, asum->length, bsum->coeffs, bsum->length, scratch2->coeffs, 1, crossover);
      else if (asum->length > bsum->length)
		   _F_mpz_poly_mul_kara_odd_even_recursive(prodsum->coeffs, bsum->coeffs, bsum->length, asum->coeffs, asum->length, scratch2->coeffs, 1, crossover);
	   else
		   _F_mpz_poly_mul_kara_recursive(prodsum, asum, bsum, scratch2, crossover);
	   _F_mpz_poly_set_length(prodsum, asum->length + bsum->length - 1);  // need to set length because odd/even kara can't do this
	} else _F_mpz_poly_set_length(prodsum, 0);
   
	// prodsum = prodsum - res_lo
   _F_mpz_poly_attach(temp, res); // res will already have the length of res_lo
	
	_F_mpz_poly_sub(prodsum, prodsum, temp);
    
   // set the final length now
	_F_mpz_poly_set_length(res, a->length + b->length - 1);

   // prodsum = prodsum - res_hi
   _F_mpz_poly_attach_shift(temp, res, 2*n1); // temp will get right length because we just set res->length
	
	_F_mpz_poly_sub(prodsum, prodsum, temp);

	// res_mid += prodsum
   temp->coeffs = res->coeffs + n1;
   temp->length = prodsum->length;
   _F_mpz_poly_add(temp, temp, prodsum);
}

ulong F_mpz_poly_kara_table[] = {20, 21, 20, 20, 16, 16, 14, 16, 14, 11, 10, 10, 10, 10};
ulong F_mpz_poly_kara_table_size = 14;

ulong _F_mpz_poly_mul_karatsuba_crossover(ulong limbs)
{
   ulong crossover;

   if (limbs >= F_mpz_poly_kara_table_size)
      crossover = 0;
   else
      crossover = F_mpz_poly_kara_table[limbs - 1];

   return crossover * crossover;
}

void _F_mpz_poly_mul_karatsuba_odd_even(F_mpz_poly_t res, F_mpz_poly_t poly1,
                            F_mpz_poly_t poly2)
{
	// number of output coefficients
   ulong length = poly1->length + poly2->length - 1;
   
   F_mpz_poly_t scratch;
	ulong crossover;

	if (poly1->length > 1)
	{
		// allocate scratch space for lower-level karatsuba routine
      F_mpz_poly_init2(scratch, length + 1);
	   scratch->length = length + 1;
	
	   // look up crossover parameter (i.e. when to switch from classical to
      // karatsuba multiplication) based on coefficient size
      ulong bits1 = FLINT_ABS(F_mpz_poly_max_bits(poly1));
	   ulong bits2 = FLINT_ABS(F_mpz_poly_max_bits(poly2));
	   ulong log_length = 0;
	   while ((1L<<log_length) < poly1->length) log_length++;
	   ulong limbs = (bits1 + bits2 + log_length - 1)/FLINT_BITS + 1;

	   crossover = _F_mpz_poly_mul_karatsuba_crossover(limbs);
	} else
	{
		crossover = 0;
		F_mpz_poly_init(scratch);
	}

   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      F_mpz_poly_t temp;
      F_mpz_poly_init2(temp, length);
		
      _F_mpz_poly_mul_kara_odd_even_recursive(temp->coeffs, poly1->coeffs, poly1->length,
            poly2->coeffs, poly2->length, scratch->coeffs, 1, crossover);

      temp->length = length;
		
		F_mpz_poly_swap(temp, res);
      F_mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      _F_mpz_poly_mul_kara_odd_even_recursive(res->coeffs, poly1->coeffs, poly1->length,
            poly2->coeffs, poly2->length, scratch->coeffs, 1, crossover);

		_F_mpz_poly_set_length(res, length);
   }
   
   F_mpz_poly_clear(scratch);
}

void _F_mpz_poly_mul_karatsuba(F_mpz_poly_t output, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      F_mpz_poly_zero(output);
      return;
   }
   
   if ((poly1->length <= 1) || (poly2->length <= 1)) // too short for karatsuba
   {
      F_mpz_poly_mul_classical(output, poly1, poly2);
      
      return;
   }
	
	ulong crossover;

   ulong limbs1 = F_mpz_poly_max_limbs(poly1);
   ulong limbs2 = F_mpz_poly_max_limbs(poly2);

   if (limbs1 + limbs2 >= 19) crossover = 0;
   else crossover = 19 - limbs1 - limbs2;

	if (poly1->length + poly2->length <= crossover) // too short for karatsuba
   {
      F_mpz_poly_mul_classical(output, poly1, poly2);
      
      return;
   }
   
	F_mpz_poly_t scratch;
   F_mpz_poly_init2(scratch, 5*poly1->length);
   
	if (output == poly1 || output == poly2)
   {
      // output is inplace, so need a temporary
      F_mpz_poly_t temp;
      F_mpz_poly_init2(temp, poly1->length + poly2->length - 1);
      _F_mpz_poly_mul_kara_recursive(temp, poly1, poly2, scratch, crossover);

      F_mpz_poly_swap(temp, output);
      F_mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      _F_mpz_poly_mul_kara_recursive(output, poly1, poly2, scratch, crossover);
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
   if (poly1->length < poly2->length)
      _F_mpz_poly_mul_karatsuba_odd_even(res, poly1, poly2);
	else if (poly1->length > poly2->length)
      _F_mpz_poly_mul_karatsuba_odd_even(res, poly2, poly1);
	else
		_F_mpz_poly_mul_karatsuba(res, poly1, poly2);
}

/*===============================================================================

	Bit packing

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
                         const ulong bits, const ulong length, const long negate)
{   
   ulong i, k, skip;
   ulong * coeff_m = poly_F_mpz->coeffs;
   mp_limb_t * end_point;
    
   ulong temp;
   half_ulong lower;
   long coeff;
   long borrow;
   mp_limb_t extend;
   
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

void F_mpz_poly_bit_unpack(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array,
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
            _F_mpz_demote(coeff_m);
				*coeff_m = ((temp & mask) + carry);
            carry = 0UL;
         }  
         else
         {
            _F_mpz_demote(coeff_m);
				*coeff_m = -(((-temp) & mask) - carry);
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
            _F_mpz_demote(coeff_m);
				*coeff_m = ((temp & mask) + carry);
            carry = 0UL;
         }
         else
         {
            _F_mpz_demote(coeff_m);
				*coeff_m = -(((-temp) & mask) - carry);
			   carry = 1UL;
         }
            
		   coeff_m ++;
         temp >>= bits;
         k -= bits;
      }

      // k is now less than bits
      skip++;
   }

   if ((carry) && (poly_F_mpz->coeffs[length - 1] == 0))
	{
		_F_mpz_demote(coeff_m);
		*coeff_m = 1L;
      poly_F_mpz->length = length + 1;
	} else
	{
		poly_F_mpz->length = length;
		_F_mpz_poly_normalise(poly_F_mpz);
	}
}

void F_mpz_poly_bit_unpack_unsigned(F_mpz_poly_t poly_F_mpz, const mp_limb_t * array, 
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
   end_point = poly_F_mpz->coeffs + length;
      
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
         _F_mpz_demote(coeff_m);
			*coeff_m = (temp & mask);
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
         _F_mpz_demote(coeff_m);
			*coeff_m = (temp & mask);
         coeff_m++;
         temp >>= bits;
         l++;
         k -= bits;
      }
      // k is now less than bits
      skip++;
   }

   poly_F_mpz->length = length;
	_F_mpz_poly_normalise(poly_F_mpz);
}
 
/*
  Bit packing used with David Harvey's KS2 algorithm

*/
void F_mpz_poly_bit_pack2(mp_limb_t * array, mp_limb_t * array2, ulong n, const F_mpz_poly_t poly_F_mpz, 
                              const ulong bits, const ulong length, const long negate, long negate2)

{   
   ulong k, skip;
   ulong * coeff_m = poly_F_mpz->coeffs;
   
   ulong temp, temp2;
   half_ulong lower, lower2;
   long coeff, coeff2;
   long borrow, borrow2;
   
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
		   temp += (coeff << k);
         temp2 += (coeff2 << k);
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
				   temp+=(coeff << k);
               temp2+=(coeff2 << k);
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
      array[skip] = temp;
      array2[skip] = temp2;
      skip++;
   } 
}

void F_mpz_poly_bit_unpack2(F_mpz_poly_t poly_F_mpz,
									 const mp_limb_t * array, const mp_limb_t * array2,
                             const ulong length, const ulong bits)
{
   ulong k, skip, s;

   ulong temp, temp2;
   ulong full_limb, full_limb2;
   ulong carry, carry2;
     
   const ulong mask = (1UL << bits) - 1;
   const ulong sign_mask = (1UL << (bits-1));

   mp_limb_t * coeff_m = poly_F_mpz->coeffs;
   mp_limb_t * end_point;

   s = 1; k = 0; skip = 0; 
	carry = 0UL; temp = 0L;
   carry2 = 0UL; temp2 = 0L;

   end_point = poly_F_mpz->coeffs + length;

   // k is zero
   // read in remainder of full_limb
   full_limb = array[skip];
   full_limb2 = array2[skip];
   temp += r_shift(full_limb, s);
   temp2 += r_shift(full_limb2, s);
   k = FLINT_BITS - s;

   while ((k >= bits) && (coeff_m < end_point))
   {
      if (!(temp2 & sign_mask)) 
      {
         _F_mpz_demote(coeff_m);
			*coeff_m = ((temp2 & mask) + carry2);
          carry2 = 0UL;
      } else
      {
         _F_mpz_demote(coeff_m);
		   *coeff_m = -(((-temp2) & mask) - carry2);
			carry2 = 1UL;
      }
      
		coeff_m ++;
      
		if (coeff_m < end_point)
		{
			if (!(temp & sign_mask)) 
         {
            _F_mpz_demote(coeff_m);
			   *coeff_m = ((temp & mask) + carry);
             carry = 0UL;
         } else
         {
            _F_mpz_demote(coeff_m);
		      *coeff_m = -(((-temp) & mask) - carry);
			   carry = 1UL;
         }
		}
            
		coeff_m ++;

      temp >>= bits;
      temp2 >>= bits;
      k -= bits;
   }
      
   skip++;
      
   while (coeff_m < end_point)
   {
      // read in a full limb
      full_limb = array[skip];
      full_limb2 = array2[skip];
      temp += l_shift(full_limb, k);
      temp2 += l_shift(full_limb2, k);
      s = FLINT_BITS - k;
      k += s;
      while ((k >= bits) && (coeff_m < end_point))
      {
         if (!(temp2 & sign_mask)) 
         {
            _F_mpz_demote(coeff_m);
			   *coeff_m = ((temp2 & mask) + carry2);
             carry2 = 0UL;
         } else
         {
            _F_mpz_demote(coeff_m);
		      *coeff_m = -(((-temp2) & mask) - carry2);
			   carry2 = 1UL;
         }
       
	   	coeff_m ++;
      
	   	if (coeff_m < end_point)
		   {
			   if (!(temp & sign_mask)) 
            {
               _F_mpz_demote(coeff_m);
			      *coeff_m = ((temp & mask) + carry);
                carry = 0UL;
            } else
            {
               _F_mpz_demote(coeff_m);
		         *coeff_m = -(((-temp) & mask) - carry);
			      carry = 1UL;
            }
	   	}

         coeff_m ++;
         temp >>= bits;
         temp2 >>= bits;
         k -= bits;
      }

      // k is now less than bits
      // read in remainder of full_limb
      temp += l_shift(r_shift(full_limb, s), k);
      temp2 += l_shift(r_shift(full_limb2, s), k);
      k += (FLINT_BITS - s);
       
      while ((k >= bits) && (coeff_m < end_point))
      {
         if (!(temp2 & sign_mask)) 
         {
            _F_mpz_demote(coeff_m);
			   *coeff_m = ((temp2 & mask) + carry2);
             carry2 = 0UL;
         } else
         {
            _F_mpz_demote(coeff_m);
		      *coeff_m = -(((-temp2) & mask) - carry2);
			   carry2 = 1UL;
         }
       
	   	coeff_m ++;
      
	   	if (coeff_m < end_point)
		   {
			   if (!(temp & sign_mask)) 
            {
               _F_mpz_demote(coeff_m);
			      *coeff_m = ((temp & mask) + carry);
                carry = 0UL;
            } else
            {
               _F_mpz_demote(coeff_m);
		         *coeff_m = -(((-temp) & mask) - carry);
			      carry = 1UL;
            }
	   	}
            
		   coeff_m ++;
         temp >>= bits;
         temp2 >>= bits;
         k -= bits;
      }
      
      skip++;
   }

   poly_F_mpz->length = length;
   _F_mpz_poly_normalise(poly_F_mpz);
}

/*===============================================================================

	Byte packing

================================================================================*/

/*
   Write a single limb to array where shift_1 and shift_2 are the amounts to shift next_limb left
	to put it into the right place in array (offset_limb being the limb of array to write it to) and
	the amount we need to shift the limb right once we are done with the part we are writing
	temp is merely a temporary limb used in the computation
*/
void __F_mpz_poly_write_next_limb(mp_limb_t * array, unsigned long * temp, unsigned long * offset_limb, 
             const unsigned long next_limb, const unsigned long shift_1, const unsigned long shift_2)
{
   *temp += l_shift(next_limb, shift_1);
   array[*offset_limb] = *temp + ((l_shift(1UL, shift_1)-1)&array[*offset_limb]);
   (*offset_limb)++;
   *temp = r_shift(next_limb, shift_2);
}

/*
   Write a single limb to array where shift_1 and shift_2 are the amounts to shift next_limb left
	to put it into the right place in array (offset_limb being the limb of array to write it to) and
	the amount we need to shift the limb right once we are done with the part we are writing
	For this function we assume that there is nothing else to be written to the output limb (except 
	zeroes)
	temp is merely a temporary limb used in the computation
*/
void __F_mpz_poly_write_whole_limb(mp_limb_t * array, unsigned long * temp, unsigned long * offset_limb, 
             const unsigned long next_limb, const unsigned long shift_1, const unsigned long shift_2)
{
   *temp += l_shift(next_limb, shift_1);
   array[*offset_limb] = *temp;
   (*offset_limb)++;
   *temp = r_shift(next_limb, shift_2);
}

void F_mpz_poly_byte_pack(mp_limb_t * array, const F_mpz_poly_t poly_fmpz,
                   const ulong length, const ulong coeff_bytes, const long negate)
{
   if (coeff_bytes == 2L)
	{
		F_mpz * coeffs = poly_fmpz->coeffs;
		short int borrow = 0;
		short int * array2 = (short int *) array;
		ulong i;
		short int temp = 0;
		if (negate > 0L)
		{
			for (i = 0; i < length; i++)
		   {
            temp = ((short int) coeffs[i]) - borrow;
				array2[i] = temp;
				if (temp < 0) borrow = 1;
				else borrow = 0;
		   }
		} else
		{
			for (i = 0; i < length; i++)
		   {
            temp = ((short int) -coeffs[i]) - borrow;
				array2[i] = temp;
				if (temp < 0) borrow = 1;
				else borrow = 0;
		   }
		}

		return;
	}
	
	F_mpz * coeff_m = poly_fmpz->coeffs;
    
   const ulong limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const ulong extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   // Start limb of the current coefficient within array
   ulong coeff_limb;
   // Additional offset in bytes of current coefficient within array
   ulong coeff_byte;
    
   // Where we are up to in the current coefficient: limbs + bytes
   ulong offset_limb;
    
   // Next limb to be written to bytes
   ulong next_limb;
   ulong temp;
   
	// the notional sign/size limb of the F_mpz
	long ssize;
   
   ulong shift_1, shift_2;
   
   fmpz_t scratch = (fmpz_t) flint_stack_alloc(limbs_per_coeff + 2); // 1 extra for sign
   
   fmpz_t co;

	mp_limb_t * coeff;
	mp_limb_t coeff_small;
    
   // when a coefficient is negative, we need to borrow from the next coefficient
   int borrow;
    
   ulong j;
    
   coeff_limb = 0; // initialise
   coeff_byte = 0;
   offset_limb = 0;
   temp = 0;
   borrow = 0;
            
   F_mpz * next_point = coeff_m + length;
         
   while (coeff_m < next_point)
   {
       // compute shifts to be used
       shift_1 = coeff_byte << 3;
       shift_2 = FLINT_BITS - shift_1;
                     
       // determine sign/size of next coefficient
		 if (!COEFF_IS_MPZ(*coeff_m)) // coeff is small
		 {
			 ssize = (long) (*coeff_m != 0L);
			 if (*coeff_m < 0L) // coefficient is negative
			 {
				 ssize = -1L; 
				 coeff_small = -*coeff_m;
				 coeff = &coeff_small;
			 } else // coefficient is positive
			    coeff = coeff_m;
		 } else // coeff is an mpz_t
		 {
			 __mpz_struct * ptr = F_mpz_ptr_mpz(*coeff_m);
			 coeff = ptr->_mp_d;
			 ssize = (long) ptr->_mp_size;
		 }
       
		 /* Coefficient is negative after borrow */
		 if (((negate > 0L) && (ssize - (long) borrow < 0L)) || ((negate < 0L) && (-ssize - (long) borrow < 0L)))
       {
          // mpz_t's store the absolute value only, so add 1 after borrow, then complement (i.e. negate after borrow)
          if (borrow) // add 1 then subtract borrow yields no change
          {
             if (ssize == 0L) next_limb = ~0L;
             else next_limb = ~coeff[0];
             co = coeff;
          } else 
          {
             if (negate > 0L) // don't negate, thus just add 1
				 {
				    if (ssize == 0L) // add 1 to zero gives 1
					 {
						 scratch[0] = 1L;
						 ssize = 1L;
					 } else if (ssize < 0L) // add 1 to a negative
					 {
						 mpn_sub_1(scratch, coeff, -ssize, 1L);
						 ssize += (long) (scratch[-ssize - 1] == 0L);
					 } else // add 1 to a positive
					 {   
						 mp_limb_t cry = mpn_add_1(scratch, coeff, ssize, 1L);
						 if (cry)
						 {
							 scratch[ssize] = cry;
							 ssize += 1;
						 }
					 }
				 } else // consider coeff to be negated, and add 1
				 {
				    if (ssize == 0L) // add 1 to zero
					 {
						 scratch[0] = 1L;
						 ssize = 1L;
					 } else if (ssize < 0L) // add 1 to a negative considered to be negated
					 {
						 mp_limb_t cry = mpn_add_1(scratch, coeff, -ssize, 1L);
						 if (cry)
						 {
							 scratch[-ssize] = cry;
							 ssize -= 1;
						 }
					 } else // add 1 to a positive considered to be negated
					 {   
						 mp_limb_t cry = mpn_sub_1(scratch, coeff, ssize, 1L);
						 ssize -= (long) (scratch[ssize - 1] == 0L);
					 }
				 }
				 if (ssize == 0L) next_limb = ~0L;
             else next_limb = ~scratch[0];
             co = scratch;
          }

          // deal with first limb of coefficient
          if (limbs_per_coeff == 0) // careful if we don't have even a full limb
          {
             if (coeff_m == next_point - 1) // this is the last coefficient, extend sign
             {
                __F_mpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                temp += l_shift(-1UL, shift_1);
                array[offset_limb] = temp;
                offset_limb++;
             } else // the coeff needs to be masked before writing
             {
                next_limb &= ((1UL<<(extra_bytes_per_coeff<<3)) - 1);
                __F_mpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                array[offset_limb] = temp;
             }
          } else
          {
             __F_mpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             // deal with remaining limbs
             for (j = 1; j < ABS(ssize); j++)
             {
                next_limb = ~co[j];
                __F_mpz_poly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             }
             // write remaining part of coefficient and fill 
             // remainder of coeff_bytes with binary 1's
             if ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte) 
             {
                temp += l_shift(-1UL, shift_1);
                array[offset_limb] = temp;
                offset_limb++;
             }
             for ( ; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = -1UL;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = -1UL;
                offset_limb++;
             }
          }
          temp = 0;
             
          borrow = 1;
       }        
       else // non-negative after borrow, and mpz_t's store absolute value only, so simply subtract 1 if borrow
       {
          if (borrow) 
          {
             if (negate > 0L) // do not negate, simply subtract 1
				 {
				    if (ssize == 0L) // subtract 1 from 0 gives -1
					 {
						 scratch[0] = 1L;
						 ssize = -1L;
					 } else if (ssize > 0L)
					 {
						 mpn_sub_1(scratch, coeff, ssize, 1L); // subtract 1 from positive
						 ssize -= (long) (scratch[ssize - 1] == 0L);
					 } else // subtract 1 from negative 
					 {   
						 mp_limb_t cry = mpn_add_1(scratch, coeff, -ssize, 1L);
						 if (cry)
						 {
							 scratch[-ssize] = cry;
							 ssize -= 1;
						 }
					 }
				 } else // consider coefficient to be negated, then subtract 1, i.e. add 1
				 {
				    if (ssize == 0L) // add 1 to 0 gives 1
					 {
						 scratch[0] = 1L;
						 ssize = 1L;
					 } else if (ssize > 0L) // add 1 to positive
					 {
						 mp_limb_t cry = mpn_add_1(scratch, coeff, ssize, 1L);
						 if (cry)
						 {
							 scratch[ssize] = cry;
							 ssize += 1;
						 }
					 } else // add 1 to negative
					 {   
						 mpn_sub_1(scratch, coeff, -ssize, 1L);
						 ssize += (long) (scratch[-ssize - 1] == 0L);
					 }
				 }
				 co = scratch;
          } else // no borrow so nothing to do
			 {
             co = coeff;
          } 
    
          
          /* Coefficient is positive after borrow (have to allow negative ssize, as we didn't bother making this positive*/
          if (ssize != 0L)
          {
             // deal with first limb of coefficient
             next_limb = co[0];
             __F_mpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             if (shift_2 == FLINT_BITS) temp = 0;
             // deal with remaining limbs
             for (j = 1; j < ABS(ssize); j++)
             {
                next_limb = co[j];
                __F_mpz_poly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             }
             // write remaining part of coefficient and extend with binary zeros
             for (; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = temp;
					 temp = 0;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB) + extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = temp;
					 temp = 0;
                offset_limb++;
             }
             temp = 0;
             borrow = 0;
          }
          /* Coefficient is zero after borrow */
          else 
          {
             temp = ((l_shift(1UL, shift_1) - 1) & array[offset_limb]); // write any bits still hanging around
             for ( ; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++) // write out zeroes for this coefficient
             {
                array[offset_limb] = temp;
					 temp = 0;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB) + extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = temp;
					 temp = 0;
                offset_limb++;
             }

             temp = 0;
             borrow = 0;
          }
       }
       
       // update information for next coefficient
       coeff_limb += limbs_per_coeff;
       coeff_byte += extra_bytes_per_coeff;
       if (coeff_byte > FLINT_BYTES_PER_LIMB) 
       {
          coeff_byte -= FLINT_BYTES_PER_LIMB;
          coeff_limb++;
       }
       offset_limb = coeff_limb;
          
       coeff_m++; // next coefficient
   }

	flint_stack_release();
}

 void F_mpz_poly_byte_pack_unsigned(mp_limb_t * array, const F_mpz_poly_t poly_fmpz,
                                                       const ulong length, const ulong coeff_bytes)
{
   F_mpz * coeff_m = poly_fmpz->coeffs;
    
   const ulong limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const ulong extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   // Start limb of the current coefficient within array
   ulong coeff_limb;
   // Additional offset in bytes of current coefficient within array
   ulong coeff_byte;
    
   // Where we are up to in the current coefficient: limbs + bytes
   ulong offset_limb;
    
   // Next limb to be written to bytes
   ulong next_limb;
   ulong temp;
   
	// the notional sign/size limb of the F_mpz
	long ssize;
   
   ulong shift_1, shift_2;
   
   mp_limb_t * coeff;
	mp_limb_t coeff_small;
     
   ulong j;
    
   coeff_limb = 0; // initialise
   coeff_byte = 0;
   offset_limb = 0;
   temp = 0;
            
   F_mpz * next_point = coeff_m + length;
         
   while (coeff_m < next_point)
   {
       // compute shifts to be used
       shift_1 = coeff_byte << 3;
       shift_2 = FLINT_BITS - shift_1;
                     
       // determine size of next coefficient
		 if (!COEFF_IS_MPZ(*coeff_m)) // coeff is small
		 {
			 ssize = (long) ((*coeff_m) != 0L);
			 coeff = coeff_m;
		 } else // coeff is an mpz_t
		 {
			 __mpz_struct * ptr = F_mpz_ptr_mpz(*coeff_m);
			 coeff = ptr->_mp_d;
			 ssize = (long) ptr->_mp_size;
		 }
          
       if (ssize != 0L) // Coefficient is non-zero
       {
          // deal with first limb of coefficient
          next_limb = coeff[0];
          __F_mpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
          if (shift_2 == FLINT_BITS) temp = 0;
          // deal with remaining limbs
          for (j = 1; j < ssize; j++)
          {
             next_limb = coeff[j];
             __F_mpz_poly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
          }
          // write remaining part of coefficient and extend with binary zeros
          for ( ; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
          {
             array[offset_limb] = temp; // may have something left over to write out
             temp = 0;
          }
          while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
              limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB) + extra_bytes_per_coeff + coeff_byte)
          {
             array[offset_limb] = temp; // may have something left over to write out
             temp = 0;
				 offset_limb++;
          }
          temp = 0UL;
		 }
       /* Coefficient is zero */
       else 
       {
          temp = ((l_shift(1UL, shift_1) - 1) & array[offset_limb]); // write any bits still hanging around
          for ( ; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++) // write out zeroes for this coefficient
          {
             array[offset_limb] = temp; // may have something left over to write out
				 temp = 0;
          }
          while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
              limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB) + extra_bytes_per_coeff + coeff_byte)
          {
             array[offset_limb] = temp; // may have something left over to write out
				 temp = 0;
             offset_limb++;
          }

          temp = 0;
       }
       
       // update information for next coefficient
       coeff_limb += limbs_per_coeff;
       coeff_byte += extra_bytes_per_coeff;
       if (coeff_byte > FLINT_BYTES_PER_LIMB) 
       {
          coeff_byte -= FLINT_BYTES_PER_LIMB;
          coeff_limb++;
       }
       offset_limb = coeff_limb;
          
       coeff_m++; // next coefficient
   }
}
 
/*
   Unpack a single coefficient of the given number of bytes from the given limb and byte of array 
	into the output array of limbs, assuming coefficients are unsigned
*/
static inline void __F_mpz_poly_unpack_bytes(mp_limb_t * output, const mp_limb_t * array, 
            const ulong limb_start, const ulong byte_start, const ulong num_bytes)
{
    const ulong limbs_to_extract = (num_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const ulong extra_bytes_to_extract = num_bytes 
                                    - (limbs_to_extract<<FLINT_LG_BYTES_PER_LIMB);
                                    
    ulong next_limb;
    ulong temp = 0;
    
    // the limb we are up to in the array and output respectively
    ulong coeff_limb = limb_start;
    ulong output_limb = 0;

    ulong shift_1, shift_2;
    
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS - shift_1;
    
    temp = array[coeff_limb];
    coeff_limb++;
    while (output_limb < limbs_to_extract) // extract limbs
    {
       next_limb = r_shift(temp,shift_1);
       temp = array[coeff_limb];
       coeff_limb++;
       next_limb += l_shift(temp,shift_2);
       output[output_limb] = next_limb;
       output_limb++;
    }
    if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start) // extract any remaining bytes
    {
       next_limb = r_shift(temp,shift_1);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    } else // there will be overhang of final bytes
    {
       next_limb = r_shift(temp,shift_1);
       temp = array[coeff_limb];
       next_limb += l_shift(temp,shift_2);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    }
}

/*
   Unpack a single coefficient of the given number of bytes from the given limb and byte of array 
	into the output array of limbs, assuming coefficients are signed
*/
static inline ulong __F_mpz_poly_unpack_signed_bytes(mp_limb_t * output, 
        const mp_limb_t * array, const ulong limb_start, const ulong byte_start, const ulong num_bytes)
{
    const ulong limbs_to_extract = (num_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const ulong extra_bytes_to_extract = num_bytes 
                                    - (limbs_to_extract<<FLINT_LG_BYTES_PER_LIMB);
                                    
    ulong next_limb;
    ulong temp = 0;
    
    // the limb we are up to in the array and output respectively
    ulong coeff_limb = limb_start;
    ulong output_limb = 0;

    ulong shift_1, shift_2;
    
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS - shift_1;
    
    ulong sign;
    
	 /* determine sign of coefficient */
    if (byte_start + extra_bytes_to_extract > FLINT_BYTES_PER_LIMB) // last few bytes will hang over last full limb
    {
       sign = array[limb_start+limbs_to_extract+1]&(1UL<<(((byte_start 
            + extra_bytes_to_extract - FLINT_BYTES_PER_LIMB)<<3)-1));
    } else if (byte_start + extra_bytes_to_extract == FLINT_BYTES_PER_LIMB) // last few bytes will end at limb boundary
    {
       sign = array[limb_start+limbs_to_extract]&(1UL<<(FLINT_BITS-1));
    } else if (byte_start + extra_bytes_to_extract == 0) // there will be no bytes to extract after last limb 
    {
       sign = array[limb_start+limbs_to_extract-1]&(1UL<<(FLINT_BITS-1));
    } else // generic case where last few bytes will not push to the end of a limb
    {
       sign = array[limb_start+limbs_to_extract]&(1UL<<(((byte_start 
            + extra_bytes_to_extract)<<3)-1));
    }
    
    if (sign) // coefficient is negative
    {
        temp = ~array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract) // write full limbs
        {
           next_limb = r_shift(temp,shift_1);
           temp = ~array[coeff_limb];
           coeff_limb++;
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start) // write final bytes
        {
           next_limb = r_shift(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else // there will be overhang of final bytes
        {
           next_limb = r_shift(temp,shift_1);
           temp = ~array[coeff_limb];
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    else // coefficient is positive
    {
        temp = array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract) // write full limbs
        {
           next_limb = r_shift(temp,shift_1);
           temp = array[coeff_limb];
           coeff_limb++;
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start) // write final bytes
        {
           next_limb = r_shift(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else // there will be overhang of the last few bytes
        {
           next_limb = r_shift(temp,shift_1);
           temp = array[coeff_limb];
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    return sign;
}   

void F_mpz_poly_byte_unpack_unsigned(F_mpz_poly_t poly_m, const mp_limb_t * array,
                               const ulong length, const ulong coeff_bytes)
{
   const ulong limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const ulong extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   const ulong limbs = ((coeff_bytes - 1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
   
   mp_limb_t * temp = (mp_limb_t*) flint_stack_alloc(limbs + 2); // allocate temporary space to unpack into
   
   ulong limb_upto = 0;
   ulong byte_offset = 0;

	mpz_t t;
   
   ulong old_length = poly_m->length; // ensure output poly is big enough to hold result
	F_mpz_poly_fit_length(poly_m, length);
	_F_mpz_poly_set_length(poly_m, length);
	F_mpz * coeff_m = poly_m->coeffs;
	if (length > old_length) F_mpn_clear(coeff_m + old_length, length - old_length); // clear any new coefficients

   for (ulong i = 0; i < length; i++)
   {
       if ((*coeff_m == 0L) && (limbs != 1)) // we are adding to a zero coeff, so just set
		 {
			 __mpz_struct * mpz_ptr = _F_mpz_promote(coeff_m);
          if (limbs + 1 > mpz_ptr->_mp_alloc) mpz_realloc(mpz_ptr, limbs + 1);
			 mp_limb_t * data = mpz_ptr->_mp_d;
			 __F_mpz_poly_unpack_bytes(data, array, limb_upto, // unpack directly into output mpz_t
                                             byte_offset, coeff_bytes);
			 ulong size = limbs;
			 while ((size) && (data[size - 1] == 0L)) size--; // normalise
			 mpz_ptr->_mp_size = size;
			 if (size <= 1) _F_mpz_demote_val(coeff_m); // coefficient may be less than FLINT_BITS - 2 bits
		 } else // add to existing coefficient
		 {
			 F_mpn_clear(temp, limbs + 2);
          __F_mpz_poly_unpack_bytes(temp + 1, array, limb_upto, 
                                             byte_offset, coeff_bytes); // unpack into temporary
          temp[0] = limbs;
          NORM(temp); // normalise
       
          if (temp[0] == 1) // single limb
		       F_mpz_add_ui(coeff_m, coeff_m, temp[1]); 
		    else if (temp[0]) // large non-zero coefficient
		    {
             t->_mp_size = temp[0]; // set up temporary mpz_t
			    t->_mp_d = temp + 1;
			    F_mpz_add_mpz(coeff_m, coeff_m, t);
		    }
		 }
      
       limb_upto += limbs_per_coeff; // get ready for next coefficient
       
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
       coeff_m++;
   }         

   flint_stack_release();
   _F_mpz_poly_normalise(poly_m);
}

void F_mpz_poly_byte_unpack(F_mpz_poly_t poly_m, const mp_limb_t * array,
                               const ulong length, const ulong coeff_bytes)
{	
	const ulong limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const ulong extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   const ulong limbs = ((coeff_bytes - 1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
   
   mp_limb_t * temp = (mp_limb_t*) flint_stack_alloc(limbs + 2); // allocate temporary to unpack into
   
   ulong limb_upto = 0;
   ulong byte_offset = 0;
   
   ulong sign;
   ulong borrow = 0;

	mpz_t t;
   
   ulong old_length = poly_m->length; // ensure output poly is big enough for result
	F_mpz_poly_fit_length(poly_m, length);
	_F_mpz_poly_set_length(poly_m, length);
	F_mpz * coeff_m = poly_m->coeffs;
	if (length > old_length) F_mpn_clear(coeff_m + old_length, length - old_length); // clear any new coeffs

   for (ulong i = 0; i < length; i++)
   {
       F_mpn_clear(temp, limbs + 2);
       sign = __F_mpz_poly_unpack_signed_bytes(temp + 1, array, limb_upto, 
                                             byte_offset, coeff_bytes); // unpack into temp
       if (sign) temp[0] = -limbs;
       else temp[0] = limbs;
       NORM(temp); // normalise

		 if (sign && !borrow) // if sign is negative, subtract 1 as otherwise we'll just have the complement 
       {
          if (temp[0] == 0L) // 0 minus 1 is -1
			 {
             temp[0] = -1L;
				 temp[1] = 1L;
			 } else
			 {
				 mp_limb_t cry = mpn_add_1(temp + 1, temp + 1, -temp[0], 1L); // minus 1 from negative
			    if (cry)
			    {
				    temp[-temp[0] + 1] = cry;
				    temp[0]--;
			    }
			 } // positive is not possible as we have (sign)
       }

       if (borrow && !sign)  // unsigned coeff which we have to add the borrow
		 {
			 if (temp[0] == 0L) // add 1 to 0 gives 1
			 {
             temp[0] = 1L;
				 temp[1] = 1L;
			 } else
			 {
				 mp_limb_t cry = mpn_add_1(temp + 1, temp + 1, temp[0], 1L); // add 1 to positive
			    if (cry)
			    {
				    temp[temp[0] + 1] = cry;
				    temp[0]++;
			    }
		    } // negative is not possible as we do not have (sign)
		 }
       
       if (temp[0] == 1) // we have a single limb, add it to output coefficient
		    F_mpz_add_ui(coeff_m, coeff_m, temp[1]);
		 else if (temp[0] == -1L) // we have a single negative limb
 		    F_mpz_sub_ui(coeff_m, coeff_m, temp[1]);
		 else if (temp[0]) // non-zero larger coefficient
		 {
          t->_mp_size = temp[0]; // set up temporary mpz_t
			 t->_mp_d = temp + 1;
			 F_mpz_add_mpz(coeff_m, coeff_m, t);
		 }
      
       borrow = 0;
       if (sign) borrow = 1;
       
       limb_upto += limbs_per_coeff; // get ready for next coefficient
       
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
       coeff_m++;
   } 

   flint_stack_release();
   _F_mpz_poly_normalise(poly_m);
}

void F_mpz_poly_pack_bytes(F_mpz_poly_t res, F_mpz_poly_t poly, ulong n, ulong bytes)
{
    // special case, poly is zero
	if (poly->length == 0)
	{
		F_mpz_poly_zero(res);
		return;
	}
	
	long i, j;
	
	ulong max_bits = FLINT_ABS(F_mpz_poly_max_bits(poly));
	
	// length after packing
	ulong short_length = (poly->length - 1)/n + 1;
	ulong limbs = ((n*bytes - 1) >> FLINT_LG_BYTES_PER_LIMB) + 1;
	
	F_mpz_poly_fit_length(res, short_length);
	
	// pack coefficients

	for (long i = 0; i < short_length - 1; i++)
	{
		long j = i*n;
		
		F_mpz_poly_t poly_p;
	   __mpz_struct * coeff;
      mp_limb_t * coeff_r;

	   _F_mpz_poly_attach_shift(poly_p, poly, j);
	   poly_p->length = FLINT_MIN(poly_p->length, n);
	   _F_mpz_poly_normalise(poly_p);
	   
	   long negate = 1L;
	   if ((poly_p->length) && (F_mpz_sgn(poly_p->coeffs + poly_p->length - 1) < 0)) 
		   negate = -1L;

	   
	   coeff = _F_mpz_promote(res->coeffs + i);
	   // one extra limb for byte pack
	   _mpz_realloc(coeff, limbs);
	   coeff_r = coeff->_mp_d;
	   F_mpn_clear(coeff_r, limbs);

		if (poly_p->length) F_mpz_poly_byte_pack(coeff_r, poly_p, n, bytes, negate);

	   // normalise number of limbs
	   coeff->_mp_size = limbs;
	   while ((coeff->_mp_size) && !(coeff_r[coeff->_mp_size - 1])) coeff->_mp_size--;
	   if (negate < 0L) coeff->_mp_size = -coeff->_mp_size;	
	   _F_mpz_demote_val(res->coeffs + i); // coeff may end up small
	}
	
	i = short_length - 1;
	j = i*n;

	F_mpz_poly_t poly_p;
	__mpz_struct * coeff;
   mp_limb_t * coeff_r;

	_F_mpz_poly_attach_shift(poly_p, poly, j);
	
	long negate = 1L;
	if (poly_p->length)
	{	
	   if (F_mpz_sgn(poly_p->coeffs + poly_p->length - 1) < 0) negate = -1L;
	}
	
	
	coeff = _F_mpz_promote(res->coeffs + i);
	_mpz_realloc(coeff, limbs);
	coeff_r = coeff->_mp_d;
   F_mpn_clear(coeff_r, limbs);
	
	if (poly_p->length)
	{	
	     F_mpz_poly_byte_pack(coeff_r, poly_p, poly_p->length, bytes, negate);
	}
	
	// normalise number of limbs
	coeff->_mp_size = limbs;
	while ((coeff->_mp_size) && !(coeff_r[coeff->_mp_size - 1])) coeff->_mp_size--;
	if (negate < 0L) coeff->_mp_size = -coeff->_mp_size;	
	_F_mpz_demote_val(res->coeffs + i); // coeff may end up small
	

	res->length = short_length;
	_F_mpz_poly_normalise(res);
}

void F_mpz_poly_unpack_bytes(F_mpz_poly_t res, F_mpz_poly_t poly, ulong n, ulong bytes)
{
    // Special case, length zero
	if (poly->length == 0)
	{
		F_mpz_poly_zero(res);
		return;
	}
	
	long i, j;
	F_mpz_poly_t poly_r;

	// one extra limb for byte_unpack
	ulong limbs = ((2*n)*bytes*8 - 1)/FLINT_BITS + 2; // number of limbs of each large coeff
	ulong length_max = n*poly->length + n - 1;

	F_mpz_poly_fit_length(res, length_max);
	
	// zero coeffs, as we will be adding to them
	res->length = length_max;
    for (i = 0; i < length_max; i++)
		F_mpz_zero(res->coeffs + i);

	mp_limb_t * arr = flint_heap_alloc(limbs);
	
	for (long i = 0; i < poly->length; i+=2)
	{
		long j;
		F_mpz_poly_t poly_r;
      int negate = 0;
		if (F_mpz_sgn(poly->coeffs + i) < 0) negate = 1; 
	    
		_F_mpz_poly_attach_shift(poly_r, res, i*n);
		poly_r->alloc = 2*n;
		poly_r->length = 2*n;

		if (negate) // negate if necessary
		   for (j = 0; j < 2*n; j++) F_mpz_neg(poly_r->coeffs + j, poly_r->coeffs + j);
		 
		F_mpn_clear(arr, limbs);
		
		F_mpz_get_limbs(arr, poly->coeffs + i);
		
		
		F_mpz_poly_byte_unpack(poly_r, arr, 2*n, bytes);
		

		if (negate) // negate if necessary
		   for (j = 0; j < 2*n; j++) F_mpz_neg(poly_r->coeffs + j, poly_r->coeffs + j);
	}

	for (long i = 1; i < poly->length; i+=2)
	{
		long j;
		F_mpz_poly_t poly_r;
      int negate = 0;
		if (F_mpz_sgn(poly->coeffs + i) < 0) negate = 1; 
	    
		_F_mpz_poly_attach_shift(poly_r, res, i*n);
		poly_r->alloc = 2*n;
		poly_r->length = 2*n;
		
		if (negate) // negate existing coefficients then add to them
			for (j = 0; j < 2*n; j++) F_mpz_neg(poly_r->coeffs + j, poly_r->coeffs + j);
	    
		F_mpn_clear(arr, limbs);
		
		F_mpz_get_limbs(arr, poly->coeffs + i);
		
		
		F_mpz_poly_byte_unpack(poly_r, arr, 2*n, bytes);
		

		if (negate) // then negate back if necessary
		   for (j = 0; j < 2*n; j++) F_mpz_neg(poly_r->coeffs + j, poly_r->coeffs + j);
	}

	flint_heap_free(arr);

	_F_mpz_poly_normalise(res);
}

/*===============================================================================

	Kronecker Segmentation multiplication

================================================================================*/

void _F_mpz_poly_mul_KS(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2, const long bits_in)
{
   long sign1;
   long sign2;
   
   ulong length1 = input1->length;
   ulong length2 = input2->length;
   
   ulong final_length = length1 + length2 - 1;
   
	F_mpz_poly_fit_length(output, final_length);

   long bits1, bits2;
   int bitpack = 0;
   ulong sign, bits;
	
   if (bits_in)
	{
      sign = (bits_in < 0);
		bits = FLINT_ABS(bits_in);
      if (bits - sign <= FLINT_BITS - 2) bitpack = 1; // we want output bits to fit into a small F_mpz for bitpacking
	} else 
	{
		bits1 = F_mpz_poly_max_bits(input1); // compute the maximum number of coeff bits for each poly
      bits2 = (input1 == input2) ? bits1 : F_mpz_poly_max_bits(input2); 

	   sign = ((bits1 < 0) || (bits2 < 0)); // an extra bit if any coefficients are signed
      ulong length = length2;
      ulong log_length = 0L;
      while ((1<<log_length) < length) log_length++;
      bits = FLINT_ABS(bits1) + FLINT_ABS(bits2) + log_length + sign; // total number of output bits
      if (bits - sign <= FLINT_BITS - 2) bitpack = 1; // we want output bits to fit into a small F_mpz for bitpacking
	}

	ulong bytes = ((bits - 1)>>3) + 1; // otherwise we byte pack with this number of bytes per output coefficient
   
   mp_limb_t * int1, * int2, * int3;
	ulong n1, n2;

   if (bitpack)
   {
      sign1 = 0L;
		sign2 = 0L;
		
		if (input1->coeffs[length1 - 1] < 0L) sign1 = -1L; // leading coefficient is negative
   
      if (input1 != input2)
      {
         if (input2->coeffs[length2 - 1] < 0L) sign2 = -1L;
      } else sign2 = sign1;
   
      n1 = (bits*length1 - 1)/FLINT_BITS + 1; // number of limbs for large integers
		n2 = (bits*length2 - 1)/FLINT_BITS + 1;

		int1 = (mp_limb_t *) flint_stack_alloc(n1);
      if (input1 != input2)
         int2 = (mp_limb_t *) flint_stack_alloc(n2);
        
      if (sign) // coeffs are signed
		{
			if (input1 != input2)
            F_mpz_poly_bit_pack(int2, n2, input2, bits, length2, sign2);
         F_mpz_poly_bit_pack(int1, n1, input1, bits, length1, sign1);
		} else // coeffs are unsigned
		{
         if (input1 != input2)
            F_mpz_poly_bit_pack_unsigned(int2, n2, input2, bits, length2);
         F_mpz_poly_bit_pack_unsigned(int1, n1, input1, bits, length1);
		}
   } else
	{
		sign1 = 1L;
		sign2 = 1L;
		
		if (F_mpz_sgn(input1->coeffs + length1 - 1) < 0) sign1 = -1L; // leading coeff is negative
   
      if (input1 != input2)
      {
         if (F_mpz_sgn(input2->coeffs + length2 - 1) < 0) sign2 = -1L;
      } else sign2 = sign1;
   
      n1 = ((bytes*length1 - 1)>>FLINT_LG_BYTES_PER_LIMB) + 1; // number of limbs for large integers
		int1 = flint_stack_alloc(n1 + 1); // extra limb required
      
		if (input1 != input2)
		{
         n2 = ((bytes*length2 - 1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
			int2 = flint_stack_alloc(n2 + 1); // extra limb required
		}

		if (sign) // coefficients are signed
		{
			
			F_mpz_poly_byte_pack(int1, input1, length1, bytes, sign1);
         
         if (input1 != input2)
            F_mpz_poly_byte_pack(int2, input2, length2, bytes, sign2);
			
		} else
		{
			
			F_mpz_poly_byte_pack_unsigned(int1, input1, length1, bytes);

         if (input1 != input2)
            F_mpz_poly_byte_pack_unsigned(int2, input2, length2, bytes);
			
		}
	}
	
   if (input1 == input2) // aliased inputs
      int2 = int1;
   
   int3 = (mp_limb_t *) flint_stack_alloc(n1 + n2); // allocate space for product large integer
	        
   mp_limb_t msl = F_mpn_mul(int3, int1, n1, int2, n2); // multiply large integers
   
   int3[n1 + n2 - 1] = msl;
   
   if (bitpack)
   {
      if (sign) F_mpz_poly_bit_unpack(output, int3, length1 + length2 - 1, bits);  // signed coeffs
      else F_mpz_poly_bit_unpack_unsigned(output, int3, length1 + length2 - 1, bits);  // unsigned coeffs
   } else
   {
		
		if (sign) 
			F_mpz_poly_byte_unpack(output, int3, length1 + length2 - 1, bytes); // signed coeffs 
      else F_mpz_poly_byte_unpack_unsigned(output, int3, length1 + length2 - 1, bytes); // unsigned coeffs
		
   }
	
   flint_stack_release(); // release int3
   if (input1 != input2)
      flint_stack_release(); // release int2
   flint_stack_release(); // release int1
     
   if ((sign1 ^ sign2) < 0L) F_mpz_poly_neg(output, output); // one of the leading coeffs was negative
}

/*
  David Harvey's KS2 algorithm.
  This is faster than ordinary KS for about lengths 100-4000. (The overhead prevents it
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

   if (input1->coeffs[length1 - 1] < 0L) sign1 = -1L;
   if (length1 & 1) sign1a = sign1;
	else sign1a = -sign1 - 1L;

   if (input1 != input2)
   {
      if (input2->coeffs[length2 - 1] < 0L) sign2 = -1L;
   } else sign2 = sign1;
   if (length2 & 1) sign2a = sign2;
	else sign2a = -sign2 - 1L;

	long bits1, bits2;
   int bitpack = 0;
   
   bits1 = F_mpz_poly_max_bits1(input1);
   bits2 = (input1 == input2) ? bits1 : F_mpz_poly_max_bits1(input2);
      
   ulong sign = ((bits1 < 0) || (bits2 < 0));
   ulong length = length2;
   ulong log_length = 0L;
   while ((1<<log_length) < length) log_length++;
	ulong bits_max = FLINT_MAX(FLINT_ABS(bits1), FLINT_ABS(bits2));
   ulong bits = FLINT_MAX(FLINT_ABS(bits1) + FLINT_ABS(bits2) + log_length + 1, 2*bits_max); 
	if ((bits%2) == 1) bits++;
   
	if ((bits - sign <= FLINT_BITS - 2) && (bits1) && (bits2)) bitpack = 1;
   
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

      if (input1 != input2)
         F_mpz_poly_bit_pack2(int2, int2b, n2, input2, bits, length2, sign2, sign2a);
      F_mpz_poly_bit_pack2(int1, int1b, n1, input1, bits, length1, sign1, sign1a);  
      bits *= 2;
   } 
	
   if (input1 == input2)
	{
		int2 = int1;
		int2b = int1b;
	}
   
   int3 = (mp_limb_t *) flint_stack_alloc(n1 + n2 + 1);
   int3b = (mp_limb_t *) flint_stack_alloc(n1 + n2);
   int4 = (mp_limb_t *) flint_stack_alloc(n1 + n2 + 1);

	int3[n1 + n2] = 0L; // unpacking code may read past the end
   int4[n1 + n2] = 0L;
          
   mp_limb_t msl = F_mpn_mul(int3, int1, n1, int2, n2);   	
	int3[n1 + n2 - 1] = msl;
   
	msl = F_mpn_mul(int3b, int1b, n1, int2b, n2);	
	int3b[n1 + n2 - 1] = msl;
   
	if ((length1 ^ length2) & 1)
	{
	   mpn_add_n(int4, int3, int3b, n1 + n2);
		mpn_rshift(int4, int4, n1 + n2, bits/2);
      mpn_sub_n(int3, int3, int3b, n1 + n2);
	} else
   {
	   mpn_sub_n(int4, int3, int3b, n1 + n2);
      mpn_rshift(int4, int4, n1 + n2, bits/2);
      mpn_add_n(int3, int3, int3b, n1 + n2);
	}
   
	if (bitpack)
   {
      F_mpz_poly_bit_unpack2(output, int4, int3, final_length, bits);  
   } 
   
   flint_stack_release(); // release int4
   flint_stack_release(); // release int3b
   flint_stack_release(); // release int3

   if (bitpack)
	{
		if (input1 != input2)
	   {
		   flint_stack_release(); // release int2b
         flint_stack_release(); // release int2
	   }
      flint_stack_release(); // release int1b
      flint_stack_release(); // release int1
	}

	if ((sign1 ^ sign2) < 0L) F_mpz_poly_neg(output, output);
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
		F_mpz_poly_init2(output, poly1->length + poly2->length - 1);
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_KS(output, poly1, poly2, 0);
		else _F_mpz_poly_mul_KS(output, poly2, poly1, 0);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
	{
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_KS(res, poly1, poly2, 0);
		else _F_mpz_poly_mul_KS(res, poly2, poly1, 0);
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

/*===============================================================================

	Conversion to and from FFT polynomials

================================================================================*/

/* 
   Convert length coefficients of an F_mpz_poly_t to an
   already initialised ZmodF_poly_t. Each coefficient will
   be represented mod p = 2^Bn+1 where n is given by the field
   n of the ZmodF_poly_t. Coefficients will be assumed to 
   be in the range [-p/2, p/2].
   Assumes 0 < length <= poly_fmpz->length 
	The maximum number of bits per coefficient is returned (negative
	if any of the coeffs were negative) if bits_known is set to zero.
*/
   
long F_mpz_poly_to_ZmodF_poly(ZmodF_poly_t poly_f, const F_mpz_poly_t poly_fmpz, 
                                                                const ulong length, const ulong bits_known)
{
   ulong size_f = poly_f->n + 1;
   F_mpz * coeffs_m = poly_fmpz->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;
   mp_limb_t * coeff;
	mp_limb_t temp;

   ulong mask = -1L;
   long bits = 0;
   long limbs = 0;
   long sign = 1;
   
   long size_j;
	int signed_c;
   
   for (ulong i = 0; i < length; i++)
   {
      signed_c = 0;
		long c = *coeffs_m;
		
		if (!COEFF_IS_MPZ(c)) // coeff is small
		{
			size_j = 1L;
			if (c < 0L) 
			{
				signed_c = 1;
				temp = -c;
				coeff = &temp;
			} else
				coeff = coeffs_m;
		} else // coeff is an mpz_t
		{
			__mpz_struct * mpz_ptr = F_mpz_ptr_mpz(*coeffs_m);
		   size_j = mpz_ptr->_mp_size;
			if (size_j < 0L) 
			{
				signed_c = 1;
				size_j = -size_j;
			}
			coeff = mpz_ptr->_mp_d;
		}

		if (!bits_known)
		{
			if (signed_c) sign = -1L;
      
		   if (size_j > limbs + 1) // coeff is at least 1 limb bigger than biggest so far, reset limbs and bits
         {
            limbs = size_j - 1;
            bits = FLINT_BIT_COUNT(coeff[size_j - 1]); 
            if (bits == FLINT_BITS) mask = 0L;
            else mask = -1L - ((1L<<bits) - 1);  
         } else if (size_j == limbs + 1) // coeff is same size as previous biggest in limbs
         {
            if (coeff[size_j - 1] & mask) // see if we have more bits than before
            {
               bits = FLINT_BIT_COUNT(coeff[size_j - 1]);   
               if (bits == FLINT_BITS) mask = 0L;
               else mask = -1L - ((1L<<bits) - 1);
            }
         }
		}

      if (signed_c) // write out FFT coefficient, ensuring sign is correct
      {
         F_mpn_negate(coeffs_f[i], coeff, size_j); 
         F_mpn_set(coeffs_f[i] + size_j, size_f - size_j); 
      } else
      {
         F_mpn_copy(coeffs_f[i], coeff, size_j); 
         F_mpn_clear(coeffs_f[i] + size_j, size_f - size_j); 
      }
      

		coeffs_m++;
   }

   poly_f->length = length; 
   
   return sign*(FLINT_BITS*limbs+bits);  
}

/* 
   Convert a ZmodF_poly_t to an F_mpz_poly_t. Coefficients will
   be taken to be in the range [-p/2, p/2] where p = 2^nB+1.
   Assumes 0 < poly_f->length 
*/

void ZmodF_poly_to_F_mpz_poly(F_mpz_poly_t poly_fmpz, const ZmodF_poly_t poly_f, const long sign)
{
   ulong n = poly_f->n;
   
   F_mpz * coeffs_m = poly_fmpz->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;

   _F_mpz_poly_set_length(poly_fmpz, poly_f->length);
   
	if (sign)
   {
      
		for (ulong i = 0; i < poly_f->length; i++)
      {
         ZmodF_normalise(coeffs_f[i], n);
         
			__mpz_struct * mpz_ptr = _F_mpz_promote(coeffs_m);
         if (mpz_ptr->_mp_alloc < n) _mpz_realloc(mpz_ptr, n);
			mp_limb_t * data = mpz_ptr->_mp_d;
			
			if (coeffs_f[i][n - 1] >> (FLINT_BITS - 1) || coeffs_f[i][n])
         {
            F_mpn_negate(data, coeffs_f[i], n);
            mpn_add_1(data, data, n, 1L);
            long size = n;
			   while ((size) && (data[size - 1] == 0L)) size--; // normalise
			   mpz_ptr->_mp_size = -size;
			   if (size >= -1L) _F_mpz_demote_val(coeffs_m); // coefficient may be less than FLINT_BITS - 2 bits
         } else
         {
            F_mpn_copy(data, coeffs_f[i], n); 
			   ulong size = n;
			   while ((size) && (data[size - 1] == 0L)) size--; // normalise
			   mpz_ptr->_mp_size = size;
			   if (size <= 1) _F_mpz_demote_val(coeffs_m); // coefficient may be less than FLINT_BITS - 2 bits
         }
         
			coeffs_m++;
      }
		
   } else 
   {
      
		for (ulong i = 0; i < poly_f->length; i++)
      {
         ZmodF_normalise(coeffs_f[i], n);
         
			__mpz_struct * mpz_ptr = _F_mpz_promote(coeffs_m);
         if (mpz_ptr->_mp_alloc < n) _mpz_realloc(mpz_ptr, n);
			mp_limb_t * data = mpz_ptr->_mp_d;
			F_mpn_copy(data, coeffs_f[i], n); 
			ulong size = n;
			while ((size) && (data[size - 1] == 0L)) size--; // normalise
			mpz_ptr->_mp_size = size;
			if (size <= 1) _F_mpz_demote_val(coeffs_m); // coefficient may be less than FLINT_BITS - 2 bits

			coeffs_m++;
      }
		
   }
   
   _F_mpz_poly_normalise(poly_fmpz);   
}

/*===============================================================================

	Schoenhage-Strassen multiplication

================================================================================*/

void _F_mpz_poly_mul_SS(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2, const long bits_in)
{
   ulong length1 = input1->length;
   
   unsigned long length2;
   if (input1 != input2)
      length2 = input2->length;
   else length2 = length1;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      F_mpz_poly_zero(output);
      return;
   }
   
   ulong output_bits;
	ulong log_length2 = 0;
   
	ulong log_length = 0; // length of largest polynomial
   while ((1<<log_length) < length1) log_length++;
      
	if (!bits_in)
	{
		ulong size1 = F_mpz_poly_max_limbs(input1); // determine size of input coefficients
      ulong size2 = F_mpz_poly_max_limbs(input2);
   
      if (input1 != input2) while ((1<<log_length2) < length2) log_length2++; // length of shortest
      else (log_length2 = log_length);
   
      /* Start with an upper bound on the number of bits needed */
   
      output_bits = FLINT_BITS * (size1 + size2) + log_length2 + 1; // guess at size of output coeffs
	} else
	{
		output_bits = FLINT_ABS(bits_in);
	}
   
	if (output_bits <= length1) // round up so that sqrt2 trick can be used
      output_bits = (((output_bits - 1) >> (log_length - 1)) + 1) << (log_length - 1);
   else // simple round up for FFT length
      output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
      
   ulong n = (output_bits - 1) / FLINT_BITS + 1; // size of FFT coeffs
   
   ZmodF_poly_t poly1, poly2, res;
   long bits1, bits2;
   
   ZmodF_poly_stack_init(poly1, log_length + 1, n, 1); // initialise FFT polynomials
   if (input1 != input2) ZmodF_poly_stack_init(poly2, log_length + 1, n, 1);
   ZmodF_poly_stack_init(res, log_length + 1, n, 1);
   
   bits1 = F_mpz_poly_to_ZmodF_poly(poly1, input1, length1, bits_in); // put coefficients into FFT polys
   if (input1 != input2) bits2 = F_mpz_poly_to_ZmodF_poly(poly2, input2, length2, bits_in);
   else bits2 = bits1;
   
	ulong sign = 0;
   
	if (!bits_in)
	{
		if ((bits1 < 0L) || (bits2 < 0L)) 
      {
         sign = 1;  
         bits1 = ABS(bits1);
         bits2 = ABS(bits2);
      }
   
      /* Recompute the length of n now that we know how large everything really is */
   
      output_bits = bits1 + bits2 + log_length2 + sign;
   
      if (output_bits <= length1) // round output bits for sqrt2 trick
         output_bits = (((output_bits - 1) >> (log_length - 1)) + 1) << (log_length - 1);
      else // ordinary round up
         output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
      
      n = (output_bits - 1) / FLINT_BITS + 1;
   
      ZmodF_poly_decrease_n(poly1, n); // decrease the number of limbs for the FFT coefficients to new n
      if (input1 != input2) ZmodF_poly_decrease_n(poly2, n);
      ZmodF_poly_decrease_n(res, n);
	} else
	{
		if (bits_in < 0L) sign = 1;
	}
                    
   if (input1 != input2) ZmodF_poly_convolution(res, poly1, poly2); // do convolution
   else ZmodF_poly_convolution(res, poly1, poly1); // aliased inputs
   ZmodF_poly_normalise(res); // normalise convolution outputs
         
   F_mpz_poly_fit_length(output, length1 + length2 - 1);
   output->length = length1 + length2 - 1; // set output length
   
	ZmodF_poly_to_F_mpz_poly(output, res, sign); // write output
   
   ZmodF_poly_stack_clear(res); // clean up
   if (input1 != input2) ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
}

void F_mpz_poly_mul_SS(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init2(output, poly1->length + poly2->length - 1);
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_SS(output, poly1, poly2, 0);
		else _F_mpz_poly_mul_SS(output, poly2, poly1, 0);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
	{
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_SS(res, poly1, poly2, 0);
		else _F_mpz_poly_mul_SS(res, poly2, poly1, 0);
	}		
}

/*===============================================================================

	Multiplication

================================================================================*/

void _F_mpz_poly_mul(F_mpz_poly_t output, F_mpz_poly_t input1, F_mpz_poly_t input2)
{
   if ((input1->length == 0) || (input2->length == 0)) // special case, length == 0
   {
      F_mpz_poly_zero(output);
      return;
   }

   if ((input1->length <= 2) && (input2->length <= 2)) // karatsuba seems to be unconditionally faster
   {
      _F_mpz_poly_mul_karatsuba(output, input1, input2);
      return;
   }
   
	if (input1->length + input2->length <= 128) // we don't want to check max limbs or max bits more than once
	{
	   ulong limbs1 = F_mpz_poly_max_limbs(input1);
      ulong limbs2 = F_mpz_poly_max_limbs(input2);

      if (limbs1 + limbs2 <= 512/FLINT_BITS)
      {
         _F_mpz_poly_mul_KS(output, input1, input2, 0);
         return;
      }
   
      if (input1->length + input2->length <= 32) 
      {
         _F_mpz_poly_mul_karatsuba(output, input1, input2);
         return;
      }
	}
   
   long bits2 = F_mpz_poly_max_bits(input2);
   long bits1 = (input1 == input2) ? bits2 : F_mpz_poly_max_bits(input1);
   
   ulong sign = ((bits1 < 0) || (bits2 < 0)); // an extra bit if any coefficients are signed
   ulong length = input2->length; // length of shortest poly
   ulong log_length = 0L;
   while ((1<<log_length) < length) log_length++;
   ulong bits = FLINT_ABS(bits1) + FLINT_ABS(bits2) + log_length + sign; // total number of output bits
   if (sign) bits = -bits; // coefficients are signed

   bits1 = FLINT_ABS(bits1);
   bits2 = FLINT_ABS(bits2);

	if (bits1 + bits2 <= 470)
   {
      _F_mpz_poly_mul_KS(output, input1, input2, bits);
      return;
   }
   
   if (3*(bits1 + bits2) >= input1->length + input2->length) // the "diagonal"
   {
      _F_mpz_poly_mul_SS(output, input1, input2, bits);
      return;
   } 
   
   _F_mpz_poly_mul_KS(output, input1, input2, bits);     
}

void F_mpz_poly_mul(F_mpz_poly_t res, F_mpz_poly_t poly1, F_mpz_poly_t poly2)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init2(output, poly1->length + poly2->length - 1);
		if (poly1->length >= poly2->length) _F_mpz_poly_mul(output, poly1, poly2);
		else _F_mpz_poly_mul(output, poly2, poly1);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
	{
		F_mpz_poly_fit_length(res, poly1->length + poly2->length - 1);
      if (poly1->length >= poly2->length) _F_mpz_poly_mul(res, poly1, poly2);
		else _F_mpz_poly_mul(res, poly2, poly1);
	}		
}

