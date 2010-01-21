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
#include <time.h>
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
#include "F_mpz_mod_poly.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL_fast_d.h"
#include "F_mpz_LLL_wrapper.h"

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
   ulong i;
   for (i = 0; i < poly->alloc; i++) // Clean up any mpz_t's
		_F_mpz_demote(poly->coeffs + i);
	if (poly->coeffs) flint_heap_free(poly->coeffs); // clean up ordinary coeffs
}

void F_mpz_poly_factor_init(F_mpz_poly_factor_t fac)
{

   fac->alloc = 5;
   fac->num_factors = 0;
   fac->factors = (F_mpz_poly_t *) flint_heap_alloc_bytes(sizeof(F_mpz_poly_t)*5);
   fac->exponents = (unsigned long *) flint_heap_alloc(5);
   unsigned long i;
   for (i = 0; i < 5; i++)
	   F_mpz_poly_init(fac->factors[i]);
}

void F_mpz_poly_factor_clear(F_mpz_poly_factor_t fac)
{
	unsigned long i;
	for (i = 0; i < fac->alloc; i++)
	   F_mpz_poly_clear(fac->factors[i]);
	free(fac->factors);
	free(fac->exponents);
}

/*===============================================================================

   F_mpz_poly_factor_t

================================================================================*/

void F_mpz_poly_factor_insert(F_mpz_poly_factor_t fac, F_mpz_poly_t poly, unsigned long exp)
{
   if (poly->length <= 1) return;   
// how much space left in the array?, 
// if none make a new one twice as big (for efficiency) and copy contents across
   if(fac->alloc == fac->num_factors)
   {
      fac->factors = (F_mpz_poly_t *) flint_heap_realloc_bytes(fac->factors, sizeof(F_mpz_poly_t)*2*fac->alloc);
      fac->exponents = (unsigned long *) flint_heap_realloc(fac->exponents, 2*fac->alloc);
      unsigned long i;
      for (i = fac->alloc; i < 2*fac->alloc; i++)
         F_mpz_poly_init(fac->factors[i]);
      fac->alloc = 2*fac->alloc;
   } 

   F_mpz_poly_set(fac->factors[fac->num_factors], poly);
   fac->exponents[fac->num_factors] = exp;
   fac->num_factors++;

}

void F_mpz_poly_factor_concat(F_mpz_poly_factor_t res, F_mpz_poly_factor_t fac)
{
   for(unsigned long i = 0; i < fac->num_factors; i++)
      F_mpz_poly_factor_insert(res, fac->factors[i], fac->exponents[i]);
}

void F_mpz_poly_factor_print(F_mpz_poly_factor_t fac)
{
   for(unsigned long i = 0; i < fac->num_factors; i++)
   {
      F_mpz_poly_print(fac->factors[i]); 
      printf(" ^ %ld\n", fac->exponents[i]);
   }	
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
      ulong i;
      for (i = poly->length; i < n; i++)
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
      long i;
      for (i = poly->length; i < n; i++)
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
      long i;
      for (i = poly->length; i < n; i++)
         F_mpz_zero(poly->coeffs + i); 
      poly->length = n+1;
   }

   F_mpz_set_mpz(poly->coeffs + n, x);
	_F_mpz_poly_normalise(poly); // we may have set leading coefficient to zero
}

void F_mpz_poly_set_coeff_F_mpz(F_mpz_poly_t poly, ulong n, const F_mpz_t x)
{
   F_mpz_poly_fit_length(poly, n+1);

   if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      long i;
      for (i = poly->length; i < n; i++)
         F_mpz_zero(poly->coeffs + i); 
      poly->length = n+1;
   }

   F_mpz_set(poly->coeffs + n, x);
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
   
	ulong i;
	for (i = 0; i < m_poly->length; i++)
		F_mpz_set_mpz(F_poly->coeffs + i, m_poly->coeffs[i]);
}

void F_mpz_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_poly_t F_poly)
{
	mpz_poly_ensure_alloc(m_poly, F_poly->length);

   m_poly->length = F_poly->length;
   
   ulong i;
   for (i = 0; i < F_poly->length; i++)
	   F_mpz_get_mpz(m_poly->coeffs[i], F_poly->coeffs + i);
}

void F_mpz_poly_to_zmod_poly(zmod_poly_t zpol, const F_mpz_poly_t fpol)
{
   unsigned long p = zpol->p;

   if (fpol->length == 0) 
   {
      zmod_poly_zero(zpol);
      return;
   } 

   F_mpz_t temp;
   F_mpz_init(temp);
   
   zmod_poly_fit_length(zpol, fpol->length);
   zpol->length = fpol->length;
   
   unsigned long i;
   for (i = 0; i < fpol->length; i++)
   {
      *(zpol->coeffs + i) = F_mpz_mod_ui(temp, fpol->coeffs + i, p);
   }

   __zmod_poly_normalise(zpol);

   F_mpz_clear(temp);
}

void zmod_poly_to_F_mpz_poly(F_mpz_poly_t fpol, const zmod_poly_t zpol)
{
   unsigned long p = zpol->p;

   if (zpol->length == 0) 
   {
      F_mpz_poly_zero(fpol);
      return;
   } 
   
   F_mpz_poly_fit_length(fpol, zpol->length);
   fpol->length = zpol->length;
   
   unsigned long i;
   for (i = 0; i < zpol->length; i++)
   {
      F_mpz_set_si(fpol->coeffs+i, zpol->coeffs[i]);
   }

   _F_mpz_poly_normalise(fpol);

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

void F_mpz_poly_fprint(const F_mpz_poly_t poly, FILE* f)
{
   char* s = F_mpz_poly_to_string(poly);
   fputs(s, f);
   free(s);
}

void F_mpz_poly_fprint_pretty(const F_mpz_poly_t poly, FILE* f, const char * x)
{
   char* s = F_mpz_poly_to_string_pretty(poly, x);
   fputs(s, f);
   free(s);
}

void F_mpz_poly_print_pretty(const F_mpz_poly_t poly, const char * x)
{
   F_mpz_poly_fprint_pretty(poly, stdout, x);
}

int F_mpz_poly_fread(F_mpz_poly_t poly, FILE* f)
{
   int ok;
   
   mpz_poly_t p;
   mpz_poly_init(p);
   ok = mpz_poly_fread(p, f);
   if (ok)
   {
      mpz_poly_to_F_mpz_poly(poly, p);
   }
   mpz_poly_clear(p);
   
   return ok;
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

		ulong i;
		for (i = 0; i < poly2->length; i++)
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

	ulong i;
	for (i = 0; i < poly1->length; i++) // check if coefficients the same
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
	ulong i;
	for (i = 0; i < poly->length; i++)
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
	
	ulong i;
	for (i = 0; i < poly->length; i++)
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

   ulong i;
   for (i = 0; i < shorter; i++) // add up to the length of the shorter poly
      F_mpz_add(res->coeffs + i, poly1->coeffs + i, poly2->coeffs + i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (i = shorter; i < poly1->length; i++)
         F_mpz_set(res->coeffs + i, poly1->coeffs + i);

   if (poly2 != res) // copy any remaining coefficients from poly2
      for (i = shorter; i < poly2->length; i++)
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

   ulong i;
   for (i = 0; i < shorter; i++) // add up to the length of the shorter poly
      F_mpz_sub(res->coeffs + i, poly1->coeffs + i, poly2->coeffs + i);   
   
   if (poly1 != res) // copy any remaining coefficients from poly1
      for (i = shorter; i < poly1->length; i++)
         F_mpz_set(res->coeffs + i, poly1->coeffs + i);

   // careful, it is *always* necessary to negate coeffs from poly2, even if this is already res
	for (i = shorter; i < poly2->length; i++) 
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
	long i;
	for (i = poly->length - 1; i >= 0; i--) 
		F_mpz_set(res->coeffs + i + n, poly->coeffs + i);

   // insert n zeroes
	for (i = 0; i < n; i++) F_mpz_zero(res->coeffs + i);
   
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
	ulong i;
	for (i = 0; i < poly->length - n; i++) 
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
	
	ulong i;
	for (i = 0; i < poly2->length; i++) 
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
	
	ulong i;
	for (i = 0; i < poly2->length; i++) 
		F_mpz_mul_si(poly1->coeffs + i, poly2->coeffs + i, x);

	_F_mpz_poly_set_length(poly1, poly2->length);
}

void F_mpz_poly_scalar_mul(F_mpz_poly_t poly1, const F_mpz_poly_t poly2, const F_mpz_t x)
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
	
	ulong i;
	for (i = 0; i < poly2->length; i++) 
   {
		F_mpz_mul2(poly1->coeffs + i, poly2->coeffs + i, x);
   }

	_F_mpz_poly_set_length(poly1, poly2->length);
}

/*===============================================================================

	Classical multiplication

================================================================================*/

void _F_mpz_poly_sqr_classical(F_mpz_poly_t res, const F_mpz_poly_t poly)
{
   F_mpz_poly_fit_length(res, 2*poly->length - 1);
	_F_mpz_poly_set_length(res, 2*poly->length - 1);

   ulong i;
   for (i = 0; i < res->length; i++)
      F_mpz_zero(res->coeffs + i);
   
   // off-diagonal products
   for (i = 1; i < poly->length; i++)
	{
		F_mpz c = poly->coeffs[i];
	   if (c)
		{
			if (!COEFF_IS_MPZ(c))
			{
				ulong j;
			   if (c < 0L) 
					for (j = 0; j < i; j++)
                  F_mpz_submul_ui(res->coeffs + i + j, poly->coeffs + j, -c);
				else
               for (j = 0; j < i; j++)
                  F_mpz_addmul_ui(res->coeffs + i + j, poly->coeffs + j, c);
			} else
         {
            ulong j;
		      for (j = 0; j < i; j++)
               F_mpz_addmul(res->coeffs + i+j, poly->coeffs + i, poly->coeffs + j);
         }
		}
	}

   // double the off-diagonal products
   for (i = 1; i < res->length - 1; i++)
      F_mpz_add(res->coeffs + i, res->coeffs + i, res->coeffs + i);
      
   // add in diagonal products
   for (i = 0; i < poly->length; i++)
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

void _F_mpz_poly_mul_classical_trunc_left(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2, ulong trunc)
{
   ulong len1 = poly1->length;
   ulong len2 = poly2->length;
   ulong start;
   long i;

   if ((len1 == 0) || (len2 == 0) || (trunc >= len1 + len2 - 1))
   {
      res->length = 0;
      return;
   }
   
   if ((len1 == 1) && (len2 == 1)) // Special case if the length of both inputs is 1
   {
      F_mpz_mul2(res->coeffs, poly1->coeffs, poly2->coeffs);      
   } else // Ordinary case
   {
      long j;
      
      // Set res[i] = poly1[i]*poly2[0] 
      if (poly2->coeffs[0])
			for (i = trunc; i < len1; i++)
            F_mpz_mul2(res->coeffs + i, poly1->coeffs + i, poly2->coeffs);
		else 
			for (i = trunc; i < len1; i++)
            F_mpz_zero(res->coeffs + i);

      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      if (len1 > trunc) start = 1;
      else start = trunc - len1 + 1;

      if (poly1->coeffs[len1 - 1])
		   for (i = start; i < len2; i++)
            F_mpz_mul2(res->coeffs + i + len1 - 1, poly1->coeffs + len1 - 1, poly2->coeffs + i);  
		else 
         for (i = start; i < len2; i++)
            F_mpz_zero(res->coeffs + i + len1 - 1);
      
      // out[i+j] += in1[i]*in2[j] 
      for (i = 0; i < len1 - 1; i++)
      {      
         F_mpz c = poly1->coeffs[i];
			if (c)
			{
				if (trunc > i) start = trunc - i;
            else start = 1;
            
            if (!COEFF_IS_MPZ(c))
				{
					if (c < 0L) 
						for (j = start; j < len2; j++)
                     F_mpz_submul_ui(res->coeffs + i + j, poly2->coeffs + j, -c);
					else
                  for (j = start; j < len2; j++)
						   F_mpz_addmul_ui(res->coeffs + i + j, poly2->coeffs + j, c);
				} else
					for (j = start; j < len2; j++)
                  F_mpz_addmul(res->coeffs + i + j, poly1->coeffs + i, poly2->coeffs + j);
			}
      }
   } 
      
   for (i = 0; i < FLINT_MIN(trunc, len1 + len2 - 1); i++)
      F_mpz_zero(res->coeffs + i);
      
   res->length = len1 + len2 - 1;
   if (trunc >= len1 + len2 - 1) _F_mpz_poly_normalise(res);   
}

void F_mpz_poly_mul_classical_trunc_left(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2, ulong trunc)
{
	if ((poly1->length == 0) || (poly2->length == 0)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

   F_mpz_poly_fit_length(res, poly1->length + poly2->length - 1);
   
	if ((poly1 == res) || (poly2 == res)) // aliased input and output
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init2(output, poly1->length + poly2->length - 1);
		_F_mpz_poly_mul_classical_trunc_left(output, poly1, poly2, trunc);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
		_F_mpz_poly_mul_classical_trunc_left(res, poly1, poly2, trunc);
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
      ulong i;
      for (i = 0; i < len2; i++)
         F_mpz_mul2(out + i*skip, in1, in2 + i*skip);
      return;
   }
   
   if (len1 * len2 < crossover)
   {
      // switch to naive multiplication

		ulong i;
		if (in2)
			for (i = 0; i < len1; i++)
            F_mpz_mul2(out + i*skip, in1 + i*skip, in2);
		else 
			for (i = 0; i < len1; i++)
            F_mpz_zero(out + i*skip);

      // Set res[i+len1-1] = in1[len1-1]*in2[i]
      const ulong term = (len1 - 1)*skip;
		if (in1[(len1 - 1)*skip])
		   for (i = 1; i < len2; i++)
            F_mpz_mul2(out + (i + len1 - 1)*skip, in1 + term, in2 + i*skip);  
	 	else 
         for (i = 1; i < len2; i++)
            F_mpz_zero(out + (i + len1 - 1)*skip);
      
      for (i = 0; i < len1 - 1; i++)
		{
			F_mpz c = in1[i*skip];
			const ulong term = i*skip;
			ulong j;
			if (!COEFF_IS_MPZ(c))
			{
				if (c < 0L) 
					for (j = 1; j < len2; j++)
                  F_mpz_submul_ui(out + (i+j)*skip, in2 + j*skip, -c);
				else
               for (j = 1; j < len2; j++)
                  F_mpz_addmul_ui(out + (i+j)*skip, in2 + j*skip, c);
			} else
				for (j = 1; j < len2; j++)
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
   for (i = 0; i < len1/2 + len2/2 - 1; i++)
   {
      F_mpz_sub(out + 2*i*skip + skip, out + 2*i*skip + skip, out + 2*(i+1)*skip);
      F_mpz_sub(out + 2*i*skip + skip, out + 2*i*skip + skip, scratch + 2*i*skip);
   }
      
   // Add A1*A2 to x^2*(B1*B2) into even slots of output
   F_mpz_set(out, scratch);
   for (i = 1; i < len1/2 + len2/2 - 1; i++)
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
         ulong i;
         for (i = 0; i < len2-2; i++)
            F_mpz_addmul(out + (i+len1-1)*skip, in1 + term1, in2 + i*skip);
         F_mpz_mul2(out + (len1+len2-3)*skip, in1 + term1, in2 + (len2-2)*skip);

         // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
         for (i = 0; i < len1-1; i++)
            F_mpz_addmul(out + (i+len2-1)*skip, in2 + term2, in1 + i*skip);
            
         // final C1*C2 term
         F_mpz_mul2(out + (len1+len2-2)*skip, in1 + term1, in2 + term2);
      }
      else
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         ulong i;
         for (i = 0; i < len2-1; i++)
            F_mpz_addmul(out + (i+len1-1)*skip, in1 + term1, in2 + i*skip);
         F_mpz_mul2(out + (len1+len2-2)*skip, in1 + term1, in2 + term2);
      }
   }
   else if (len2 & 1)
   {
      const ulong term1 = skip*(len1-1);
	   const ulong term2 = skip*(len2-1);
	   // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
      ulong i;
      for (i = 0; i < len1-1; i++)
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
	ulong i;
	for (i = start; i < 2*n1; i++)
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

		res->length = length;
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

void _F_mpz_poly_karatrunc_left_recursive(F_mpz_poly_t res, const F_mpz_poly_t a, const F_mpz_poly_t b, F_mpz_poly_t scratch, const ulong crossover, const ulong trunc)
{
   F_mpz_poly_t temp, temp2;
   
   long non_zero = a->length + b->length - trunc - 1;
   
   if (non_zero <= 0)
   {
      long i;
      for (i = 0; i < (long) (a->length + b->length - 1); i++)
         F_mpz_zero(res->coeffs + i);
      
      res->length = 0;

      return;
   }
   
   if ((a->length <= 1) || (b->length <= 1) || (non_zero == 1)) 
   {
      _F_mpz_poly_mul_classical_trunc_left(res, a, b, trunc);
      return;
   }
   
   if ((a->length == 2 && b->length == 2) && (crossover < 4) && (!trunc)) 
   {
      F_mpz_mul2(res->coeffs, a->coeffs, b->coeffs); 
         
      F_mpz_add(scratch->coeffs, a->coeffs, a->coeffs + 1);
      F_mpz_add(scratch->coeffs + 1, b->coeffs, b->coeffs + 1);
         
      F_mpz_mul2(res->coeffs + 2, a->coeffs + 1, b->coeffs + 1); 
         
      F_mpz_mul2(res->coeffs + 1, scratch->coeffs, scratch->coeffs + 1); 
         
      F_mpz_sub(res->coeffs + 1, res->coeffs + 1, res->coeffs);
      F_mpz_sub(res->coeffs + 1, res->coeffs + 1, res->coeffs + 2);
      
            
      res->length = a->length + b->length - 1;
      
      return;
   }
   
   if ((a->length + b->length <= crossover) || ((a->length == 2) && (b->length == 2)))
   {
      _F_mpz_poly_mul_classical_trunc_left(res, a, b, trunc);
      
      return;
   }   
        
   F_mpz_poly_t a1, a2, b1, b2;
      
   ulong l2 = 0;
   ulong old_length;
     
   ulong m = (a->length + 1)/2;
   _F_mpz_poly_attach_truncate(a1, a, m);
   _F_mpz_poly_attach_shift(a2, a, m);
   
   if (m < b->length) //ordinary case
   {
      /*
         (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
      */
      
      _F_mpz_poly_attach_truncate(b1, b, m);
      _F_mpz_poly_attach_shift(b2, b, m);
   
      /* 
         from 0 for a1->length + b1->length - 1 will be directly written to, 
         as will all coeffs from 2*m onwards.
      */
      ulong start = a1->length + b1->length - 1;
      if (!a1->length || !b1->length) start = 0;
      ulong i;
      for (i = start; i < 2*m; i++)
         F_mpz_zero(res->coeffs + i);
  
      F_mpz_poly_t asum, bsum, prodsum, scratch2;
     
      asum->length = m;
      asum->coeffs = scratch->coeffs;
      
      bsum->length = m;
      bsum->coeffs = scratch->coeffs + m;
      
      prodsum->length = 2*m - 1;
      prodsum->coeffs = scratch->coeffs + 2*m;    

      /*
         (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
      */
      
      // res_lo = a1*b1
      if (trunc > m) 
      {
         if (a1->length >= b1->length) _F_mpz_poly_karatrunc_left_recursive(res, a1, b1, scratch, crossover, trunc - m);
         else _F_mpz_poly_karatrunc_left_recursive(res, b1, a1, scratch, crossover, trunc - m);
      }
      else 
      {
         if (a1->length >= b1->length) _F_mpz_poly_karatrunc_left_recursive(res, a1, b1, scratch, crossover, 0);
         else _F_mpz_poly_karatrunc_left_recursive(res, b1, a1, scratch, crossover, 0);
      }

      // res_hi = a2*b2
      temp->coeffs = res->coeffs + 2*m;
      
      if (trunc > 2*m) _F_mpz_poly_karatrunc_left_recursive(temp, a2, b2, scratch, crossover, trunc - 2*m);
      else _F_mpz_poly_karatrunc_left_recursive(temp, a2, b2, scratch, crossover, 0);
      
      if (trunc < 3*m - 1)
      {
         // asum = a1+a2
         _F_mpz_poly_add(asum, a1, a2);
         // bsum = b1+b2
         _F_mpz_poly_add(bsum, b1, b2);
         // prodsum = asum*bsum
         
         scratch2->coeffs = scratch->coeffs + 4*m - 1;
         
         if (trunc > m) 
         {
            if (asum->length > bsum->length) _F_mpz_poly_karatrunc_left_recursive(prodsum, asum, bsum, scratch2, crossover, trunc - m);
            else _F_mpz_poly_karatrunc_left_recursive(prodsum, bsum, asum, scratch2, crossover, trunc - m);
         } else 
         {
            if (asum->length > bsum->length) _F_mpz_poly_karatrunc_left_recursive(prodsum, asum, bsum, scratch2, crossover, 0);
            else _F_mpz_poly_karatrunc_left_recursive(prodsum, bsum, asum, scratch2, crossover, 0);
         }

         long i;
         for (i = prodsum->length; i < 2*m - 1; i++)
            F_mpz_zero(prodsum->coeffs + i);
               
         // prodsum = prodsum - res_lo
         temp->coeffs = res->coeffs;
         temp->length = 2*m - 1;
         _F_mpz_poly_sub(prodsum, prodsum, temp);
       
         // prodsum = prodsum - res_hi
         temp->coeffs = res->coeffs + 2*m;
         temp->length = a2->length + b2->length - 1;
         _F_mpz_poly_sub(prodsum, prodsum, temp);
      
         // res_mid += prodsum
         temp->coeffs = res->coeffs + m;
         temp->length = prodsum->length;
         _F_mpz_poly_add(temp, temp, prodsum);     
      }
      
      res->length = a->length + b->length - 1;
      
   } else 
   {
      F_mpz_poly_t scratch2, temp1; 

      while ((1<<l2) < m) l2++;
      if ((1<<l2) < a->length) m = (1<<l2);
      
      _F_mpz_poly_attach_truncate(a1, a, m);
      _F_mpz_poly_attach_shift(a2, a, m);
      
      /* 
         from 0 for a1->length + b->length - 1 will be directly written to.
      */
      ulong start = a1->length + b->length - 1;
      if (!a1->length) start = 0;
      ulong i;
      for (i = start; i < a->length + b->length - 1; i++)
         F_mpz_zero(res->coeffs + i);
  
      // res_lo = a1*b
      if (trunc < a1->length + b->length - 1) 
      {
         if (a1->length >= b->length) _F_mpz_poly_karatrunc_left_recursive(res, a1, b, scratch, crossover, trunc);
         else _F_mpz_poly_karatrunc_left_recursive(res, b, a1, scratch, crossover, trunc);
      }

      //temp = a2*b
      temp->coeffs = scratch->coeffs;
      temp->length = a2->length + b->length - 1;
      
      scratch2->coeffs = scratch->coeffs + temp->length;
      
      if (trunc > m)
      {
         if (b->length <= a2->length) _F_mpz_poly_karatrunc_left_recursive(temp, a2, b, scratch2, crossover, trunc - m);
         else _F_mpz_poly_karatrunc_left_recursive(temp, b, a2, scratch2, crossover, trunc - m);
      } else
      {
         if (b->length <= a2->length) _F_mpz_poly_karatrunc_left_recursive(temp, a2, b, scratch2, crossover, 0);
         else _F_mpz_poly_karatrunc_left_recursive(temp, b, a2, scratch2, crossover, 0);
      }
      
      // res_mid += temp
      temp1->coeffs = res->coeffs + m;
      temp1->length = temp->length;
      
      _F_mpz_poly_add(temp1, temp1, temp); 
  
      res->length = a->length + b->length - 1;
   } 
   
   ulong i;
   for (i = 0; i < trunc; i++)
      F_mpz_zero(res->coeffs + i);  
}

void _F_mpz_poly_mul_karatsuba_trunc_left(F_mpz_poly_t output, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2, ulong trunc)
{
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      F_mpz_poly_zero(output);
      return;
   }
   
   if ((poly1->length <= 1) || (poly2->length <= 1)) // too short for karatsuba
   {
      F_mpz_poly_mul_classical_trunc_left(output, poly1, poly2, trunc);
      
      return;
   }
	
	ulong crossover;

   ulong limbs1 = F_mpz_poly_max_limbs(poly1);
   ulong limbs2 = F_mpz_poly_max_limbs(poly2);

   if (limbs1 + limbs2 >= 19) crossover = 0;
   else crossover = 19 - limbs1 - limbs2;

	if (poly1->length + poly2->length <= crossover) // too short for karatsuba
   {
      F_mpz_poly_mul_classical_trunc_left(output, poly1, poly2, trunc);
      
      return;
   }
   
	F_mpz_poly_t scratch;
   F_mpz_poly_init2(scratch, 5*poly1->length);
   
	if (output == poly1 || output == poly2)
   {
      // output is inplace, so need a temporary
      F_mpz_poly_t temp;
      F_mpz_poly_init2(temp, poly1->length + poly2->length - 1);
      _F_mpz_poly_karatrunc_left_recursive(temp, poly1, poly2, scratch, crossover, trunc);

      F_mpz_poly_swap(temp, output);
      F_mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      _F_mpz_poly_karatrunc_left_recursive(output, poly1, poly2, scratch, crossover, trunc);
   }

	F_mpz_poly_clear(scratch);
}

void F_mpz_poly_mul_karatsuba_trunc_left(F_mpz_poly_t res, F_mpz_poly_t poly1,
                              F_mpz_poly_t poly2, ulong trunc)
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
      F_mpz_poly_sqr_karatsuba_trunc_left(res, poly1, trunc);
      return;
   }*/

   F_mpz_poly_fit_length(res, poly1->length + poly2->length - 1);
	// rearrange parameters to make poly1 no longer than poly2
   if (poly1->length >= poly2->length)
      _F_mpz_poly_mul_karatsuba_trunc_left(res, poly1, poly2, trunc);
	else 
      _F_mpz_poly_mul_karatsuba_trunc_left(res, poly2, poly1, trunc);
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

   ulong i;
   for (i = 0; i < length; i++)
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

   ulong i;
   for (i = 0; i < length; i++)
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

	for (i = 0; i < short_length - 1; i++)
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
	
	for (i = 0; i < poly->length; i+=2)
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

	for (i = 1; i < poly->length; i+=2)
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
   {
      int2 = int1;
      n2 = n1;
   }
   
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
   
   ulong i;
   for (i = 0; i < length; i++)
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
      
		ulong i;
		for (i = 0; i < poly_f->length; i++)
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
      
		ulong i;
		for (i = 0; i < poly_f->length; i++)
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

void _F_mpz_poly_mul(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2)
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

void F_mpz_poly_mul(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2)
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

void _F_mpz_poly_mul_trunc_left(F_mpz_poly_t output, const F_mpz_poly_t input1, const F_mpz_poly_t input2, ulong trunc)
{
   if ((input1->length == 0) || (input2->length == 0) || (input1->length + input2->length <= trunc + 1)) // special case, length == 0
   {
      F_mpz_poly_zero(output);
      return;
   }

   if ((input1->length <= 2) && (input2->length <= 2)) // karatsuba seems to be unconditionally faster
   {
      _F_mpz_poly_mul_karatsuba_trunc_left(output, input1, input2, trunc);
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
         _F_mpz_poly_mul_karatsuba_trunc_left(output, input1, input2, trunc);
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

void F_mpz_poly_mul_trunc_left(F_mpz_poly_t res, const F_mpz_poly_t poly1, const F_mpz_poly_t poly2, ulong trunc)
{
   if ((poly1->length == 0) || (poly2->length == 0) || (poly1->length + poly2->length <= trunc + 1)) // special case if either poly is zero
   {
      F_mpz_poly_zero(res);
      return;
   }

	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_poly_t output; // create temporary
		F_mpz_poly_init2(output, poly1->length + poly2->length - 1);
		if (poly1->length >= poly2->length) _F_mpz_poly_mul_trunc_left(output, poly1, poly2, trunc);
		else _F_mpz_poly_mul_trunc_left(output, poly2, poly1, trunc);
		F_mpz_poly_swap(output, res); // swap temporary with real output
		F_mpz_poly_clear(output);
	} else // ordinary case
	{
		F_mpz_poly_fit_length(res, poly1->length + poly2->length - 1);
      if (poly1->length >= poly2->length) _F_mpz_poly_mul_trunc_left(res, poly1, poly2, trunc);
		else _F_mpz_poly_mul_trunc_left(res, poly2, poly1, trunc);
	}		
}
/*===============================================================================

	Division with remainder

================================================================================*/

void F_mpz_poly_divrem_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   if (B->length == 0)
   {
      printf("Exception : Divide by zero in F_mpz_poly_divrem_basecase.\n");
      abort();
   }
   
   ulong coeff = A->length;
   ulong B_length = B->length;

   F_mpz_poly_t qB;
   
   F_mpz * coeffs_A = A->coeffs;
   F_mpz * coeffs_B = B->coeffs;
   F_mpz * B_lead = coeffs_B + B_length - 1; 
   F_mpz * coeff_Q;
   F_mpz * coeffs_R;

   int want_rem = 1;
   F_mpz_poly_struct Rs;
   if (R == NULL) 
   {
      want_rem = 0;
      R = &Rs;
      F_mpz_poly_init(R);
   }

   while (coeff >= B_length)
   {
      if (F_mpz_cmpabs(coeffs_A + coeff - 1, B_lead) >= 0) break;
      else coeff--;   
   }
   
   if (want_rem) F_mpz_poly_set(R, A);
   
   if (coeff >= B_length)
   {
      F_mpz_poly_fit_length(Q, coeff - B_length + 1);    
      Q->length = coeff - B_length + 1;
   } else 
   {
      F_mpz_poly_zero(Q);
      return;
   }
    
   if (!want_rem) 
   {
      F_mpz_poly_fit_length(R, coeff);
      R->length = coeff;
      ulong i;
      for (i = B_length - 1; i < coeff; i++)
         F_mpz_set(R->coeffs + i, A->coeffs + i);
   }
   
   coeffs_R = R->coeffs; 

   F_mpz_poly_t Bsub;
   if (want_rem) _F_mpz_poly_attach(Bsub, B);
   else _F_mpz_poly_attach_truncate(Bsub, B, B_length - 1);
   ulong Bsub_length = B_length;
   
   F_mpz_poly_init2(qB, B_length);
   while (coeff >= B_length)
   {
      coeff_Q = Q->coeffs + coeff - B_length;

      if (F_mpz_cmpabs(coeffs_R + coeff - 1, B_lead) < 0) F_mpz_zero(coeff_Q);
      else
      {
         F_mpz_fdiv_q(coeff_Q, coeffs_R + coeff - 1, B_lead);
         
         F_mpz_poly_scalar_mul(qB, Bsub, coeff_Q); 
         
         coeffs_R = R->coeffs;
         
         F_mpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R + coeff - Bsub_length;
         R_sub->length = Bsub->length;
         _F_mpz_poly_sub(R_sub, R_sub, qB);
      }
      
      if ((!want_rem) && (Bsub->length >= coeff - B_length + 1))
      {
         Bsub->coeffs++;
         Bsub->length--;
         Bsub_length--;
      }
      
      coeff--;
   }
         
   F_mpz_poly_clear(qB);
   
   if (want_rem) _F_mpz_poly_normalise(R);
}

void F_mpz_poly_div_divconquer_recursive(F_mpz_poly_t Q, F_mpz_poly_t BQ, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   if (A->length < B->length)
   {
      F_mpz_poly_zero(Q);
      F_mpz_poly_zero(BQ);

      return;
   }
   
   // A->length is now >= B->length
   
   ulong crossover = 16;
   
   if (A->length - B->length + 1 <= crossover) 
   {
      /*
         Use the classical algorithm to compute the
         quotient and remainder, then use A - R to compute BQ
      */
      
      F_mpz_poly_t Rb;
      F_mpz_poly_init(Rb);
      F_mpz_poly_divrem_basecase(Q, Rb, A, B);
      F_mpz_poly_fit_length(BQ, A->length);
      F_mpz_poly_sub(BQ, A, Rb);
      F_mpz_poly_clear(Rb);
      
      return;
   }
   
   F_mpz_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   ulong n1 = (B->length + 1)/2;
   ulong n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _F_mpz_poly_attach_shift(d1, B, n2);
   _F_mpz_poly_attach_truncate(d2, B, n2);
   _F_mpz_poly_attach_shift(d3, B, n1);
   _F_mpz_poly_attach_truncate(d4, B, n1);
   
   if (A->length < 2*B->length - 1)
   {
      /* Convert unbalanced division into a 2*q - 1 by q division */
      F_mpz_poly_t t_A, t_B, t_B2;
      
      ulong q = A->length - B->length + 1;
      ulong q2 = B->length - q;

      _F_mpz_poly_attach_shift(t_A, A, A->length - 2*q + 1);
      _F_mpz_poly_attach_shift(t_B, B, q2);
      _F_mpz_poly_attach_truncate(t_B2, B, q2);
      
      F_mpz_poly_init(d1q1);
      F_mpz_poly_div_divconquer_recursive(Q, d1q1, t_A, t_B); 
      
      /*
         Compute d2q1 = Q*t_B2
         It is of length q2*q terms
      */
      
      F_mpz_poly_init(d2q1);
      F_mpz_poly_mul(d2q1, Q, t_B2);
      
      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length at most n1+n2-1
      */
      
      F_mpz_poly_fit_length(BQ, FLINT_MAX(d1q1->length + B->length - q, d2q1->length));
      F_mpz_poly_left_shift(BQ, d1q1, B->length - q);
      F_mpz_poly_clear(d1q1);
      _F_mpz_poly_add(BQ, BQ, d2q1);
      F_mpz_poly_clear(d2q1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      ulong shift = A->length - 2*B->length + 1;
      _F_mpz_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = d1*q1 is length at most 2*B->length-1
      */
      
      F_mpz_poly_init(d1q1);
      F_mpz_poly_init(q1);
      
      F_mpz_poly_div_divconquer_recursive(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      F_mpz_poly_init(dq1);
      
      F_mpz_poly_fit_length(dq1, d1q1->length + shift);
      F_mpz_poly_left_shift(dq1, d1q1, shift);
      F_mpz_poly_clear(d1q1);
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      F_mpz_poly_init(t);
      F_mpz_poly_sub(t, A, dq1);
      F_mpz_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      F_mpz_poly_init(q2);
      F_mpz_poly_init(dq2);
      F_mpz_poly_div_divconquer_recursive(q2, dq2, t, B); 
      F_mpz_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      F_mpz_poly_fit_length(Q, FLINT_MAX(q1->length + shift, q2->length));
      
      F_mpz_poly_left_shift(Q, q1, shift);
      F_mpz_poly_clear(q1);
      F_mpz_poly_add(Q, Q, q2);
      F_mpz_poly_clear(q2);
      
      /*
         Write out BQ = dq1 + dq2
      */
      
      F_mpz_poly_fit_length(BQ, FLINT_MAX(dq1->length, dq2->length));
      
      F_mpz_poly_add(BQ, dq1, dq2);
      F_mpz_poly_clear(dq1);
      F_mpz_poly_clear(dq2);
      
      return;
   } 
   
   // n2 + B->length - 1 < A->length <= n1 + n2 + B->length - 1
    
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length at most n1 (and at least 1), 
      a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n1-1)+ other terms, so it has 
      length at most 2*n1-1 
   */
      
   _F_mpz_poly_attach_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length n1
      d1q1 = d1*q1 is length at most 2*n1-1
   */
      
   F_mpz_poly_init(d1q1);
   F_mpz_poly_init(q1);
   F_mpz_poly_div_divconquer_recursive(q1, d1q1, p1, d1); 
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being at most length n1+n2-1
   */  
   
   F_mpz_poly_init(d2q1);
   F_mpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most 2*n1+n2-1
   */
   
   F_mpz_poly_init2(dq1, FLINT_MAX(d1q1->length + n2, d2q1->length));
   F_mpz_poly_left_shift(dq1, d1q1, n2);
   F_mpz_poly_clear(d1q1);
   _F_mpz_poly_add(dq1, dq1, d2q1);
   F_mpz_poly_clear(d2q1);
   
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1
      which has length at most 2*n1+n2-1, but we are not interested 
      in up to the first n1 coefficients, so it has 
      effective length at most n1+n2-1
   */
   
   F_mpz_poly_init2(t, FLINT_MAX(A->length-n2, dq1->length));
   F_mpz_poly_right_shift(t, A, n2);
   _F_mpz_poly_sub(t, t, dq1);
   F_mpz_poly_truncate(t, B->length - 1);
   
   /*
      Compute q2 = t div d1
      It is at most an n1+n2-1 by n1 division, so
      the length of q2 will be at most n2
      Also compute d1q2 of length at most n1+n2-1
   */
   
   F_mpz_poly_init(d1q2);
   F_mpz_poly_init(q2);
   F_mpz_poly_div_divconquer_recursive(q2, d1q2, t, d1); 
   F_mpz_poly_clear(t);
      
   /*
      Compute d2q2 = d2*q2 which is of length 
      at most n1+n2-1
   */
   
   F_mpz_poly_init(d2q2);
   F_mpz_poly_mul(d2q2, d2, q2);
   
   /*
      Compute dq2 = d1*q2*x^n2 + d2q2
      which is of length at most n1+2*n2-1
   */
   
   F_mpz_poly_init2(dq2, FLINT_MAX(d1q2->length+n2, d2q2->length));
   F_mpz_poly_left_shift(dq2, d1q2, n2);
   F_mpz_poly_clear(d1q2);
   _F_mpz_poly_add(dq2, dq2, d2q2);
   F_mpz_poly_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length at most n1+n2
   */
   
   F_mpz_poly_fit_length(Q, FLINT_MAX(q1->length+n2, q2->length));
   F_mpz_poly_left_shift(Q, q1, n2);
   F_mpz_poly_clear(q1);
   _F_mpz_poly_add(Q, Q, q2);
   F_mpz_poly_clear(q2);
   
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length at most 2*(n1+n2)-1
   */
   
   F_mpz_poly_fit_length(BQ, FLINT_MAX(n2 + dq1->length, dq2->length));
   F_mpz_poly_left_shift(BQ, dq1, n2);
   _F_mpz_poly_add(BQ, BQ, dq2);
   
   F_mpz_poly_clear(dq2);
   F_mpz_poly_clear(dq1);
}

void F_mpz_poly_divrem_divconquer(F_mpz_poly_t Q, F_mpz_poly_t R, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   F_mpz_poly_t QB;
   
   F_mpz_poly_init(QB);
   
   F_mpz_poly_div_divconquer_recursive(Q, QB, A, B);
   
   F_mpz_poly_fit_length(R, A->length);
   _F_mpz_poly_sub(R, A, QB);
   _F_mpz_poly_normalise(R);
   
   F_mpz_poly_clear(QB);
}

/*===============================================================================

	Division without remainder

================================================================================*/

void F_mpz_poly_divrem_basecase_low(F_mpz_poly_t Q, F_mpz_poly_t R, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   if (B->length == 0)
   {
      printf("Exception : Divide by zero in F_mpz_poly_divrem_basecase_low.\n");
      abort();
   }
   
   ulong coeff = A->length;
   ulong B_length = B->length;

   F_mpz_poly_t qB;
   
   F_mpz * coeffs_A = A->coeffs;
   F_mpz * coeffs_B = B->coeffs;
   F_mpz * B_lead = coeffs_B + B_length - 1; 
   F_mpz * coeff_Q;
   F_mpz * coeffs_R;

   while (coeff >= B_length)
   {
      if (F_mpz_cmpabs(coeffs_A + coeff - 1, B_lead) >= 0) break;
      else coeff--;   
   }
   
   F_mpz_poly_set(R, A);
   
   if (coeff >= B_length)
   {
      F_mpz_poly_fit_length(Q, coeff - B_length + 1);    
      Q->length = coeff - B_length + 1;
   } else 
   {
      F_mpz_poly_zero(Q);
      F_mpz_poly_truncate(R, B_length - 1);
      return;
   }   
   
   coeffs_R = R->coeffs; 

   F_mpz_poly_t Bsub;
   _F_mpz_poly_attach_truncate(Bsub, B, B_length - 1);
   
   F_mpz_poly_init2(qB, Bsub->length);
   while (coeff >= B_length)
   {
      coeff_Q = Q->coeffs + coeff - B_length;

      if (F_mpz_cmpabs(coeffs_R + coeff - 1, B_lead) < 0) F_mpz_zero(coeff_Q);
      else
      {
         F_mpz_fdiv_q(coeff_Q, coeffs_R + coeff - 1, B_lead);
         
         F_mpz_poly_scalar_mul(qB, Bsub, coeff_Q); 
         
         coeffs_R = R->coeffs;
         
         F_mpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R + coeff - B_length;
         R_sub->length = Bsub->length;
         _F_mpz_poly_sub(R_sub, R_sub, qB);
      }
      
      coeff--;
   }
         
   F_mpz_poly_clear(qB);
   
   F_mpz_poly_truncate(R, B_length - 1);
}

void F_mpz_poly_div_divconquer_recursive_low(F_mpz_poly_t Q, F_mpz_poly_t BQ, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   if (A->length < B->length)
   {
      F_mpz_poly_zero(Q);
      F_mpz_poly_zero(BQ);
      
      return;
   }
   
   // A->length is now >= B->length
   
   ulong crossover = 16;
   
   if (A->length - B->length + 1 <= crossover) 
   {
      /*
         Use the classical algorithm to compute the
         quotient and low half of the remainder, then 
         truncate A-R to compute BQ
      */
      
      F_mpz_poly_t Rb;
      F_mpz_poly_init(Rb);
      F_mpz_poly_divrem_basecase_low(Q, Rb, A, B);
      F_mpz_poly_fit_length(BQ, A->length);
      _F_mpz_poly_sub(BQ, A, Rb);
      F_mpz_poly_clear(Rb);
      F_mpz_poly_truncate(BQ, B->length - 1);
      
      return;
   }
   
   F_mpz_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   ulong n1 = (B->length + 1)/2;
   ulong n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _F_mpz_poly_attach_shift(d1, B, n2);
   _F_mpz_poly_attach_truncate(d2, B, n2);
   _F_mpz_poly_attach_shift(d3, B, n1);
   _F_mpz_poly_attach_truncate(d4, B, n1);
   
   if (A->length < 2*B->length - 1)
   {
      /* Convert unbalanced division into a 2*q - 1 by q division */
      F_mpz_poly_t t_A, t_B, t_B2;
      
      ulong q = A->length - B->length + 1;
      ulong q2 = B->length - q;

      _F_mpz_poly_attach_shift(t_A, A, A->length - 2*q + 1);
      _F_mpz_poly_attach_shift(t_B, B, q2);
      _F_mpz_poly_attach_truncate(t_B2, B, q2);
      
      F_mpz_poly_init(d1q1);
      F_mpz_poly_div_divconquer_recursive_low(Q, d1q1, t_A, t_B); 
      
      /*
         Compute d2q1 = Q*t_B2
         It is of length q2*q terms
      */
      
      F_mpz_poly_init(d2q1);
      F_mpz_poly_mul(d2q1, Q, t_B2);
      
      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length at most n1+n2-1
      */
      
      F_mpz_poly_fit_length(BQ, FLINT_MAX(d1q1->length + B->length - q, d2q1->length));
      F_mpz_poly_left_shift(BQ, d1q1, B->length - q);
      F_mpz_poly_clear(d1q1);
      _F_mpz_poly_add(BQ, BQ, d2q1);
      F_mpz_poly_clear(d2q1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length - 1
      // We call this polynomial p1
      
      ulong shift = A->length - 2*B->length + 1;
      _F_mpz_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = d1*q1 is truncated to length at most B->length-1
      */
      
      F_mpz_poly_init(d1q1);
      F_mpz_poly_init(q1);
      
      F_mpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length - B->length
         dq1 is normalised since d1q1 was
      */
   
      F_mpz_poly_init(dq1);
      
      F_mpz_poly_fit_length(dq1, d1q1->length + shift);
      F_mpz_poly_left_shift(dq1, d1q1, shift);
      F_mpz_poly_clear(d1q1);
      
      /*
         Compute t = A - dq1 
         We truncate, leaving at most A->length - B->length 
         significant terms
      */
   
      F_mpz_poly_init(t);
      F_mpz_poly_sub(t, A, dq1);
      F_mpz_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
         dq2 has length at most B->length - 1
      */
   
      F_mpz_poly_init(q2);
      F_mpz_poly_init(dq2);
      F_mpz_poly_div_divconquer_recursive_low(q2, dq2, t, B); 
      F_mpz_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      F_mpz_poly_fit_length(Q, FLINT_MAX(q1->length+shift, q2->length));
      
      F_mpz_poly_left_shift(Q, q1, shift);
      F_mpz_poly_clear(q1);
      _F_mpz_poly_add(Q, Q, q2);
      F_mpz_poly_clear(q2);
      
      /*
         Write out BQ = dq1 + dq2
      */
      
      F_mpz_poly_fit_length(BQ, FLINT_MAX(dq1->length, dq2->length));
      
      _F_mpz_poly_add(BQ, dq1, dq2);
      F_mpz_poly_truncate(BQ, B->length - 1);
      F_mpz_poly_clear(dq1);
      F_mpz_poly_clear(dq2);
      
      return;
   } 
   
   // A->length == 2*B->length - 1
    
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length n1, 
      a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n1-1)+ other terms, so it has 
      length 2*n1-1 
   */
      
   _F_mpz_poly_attach_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is a 2*n1-1 by n1 division so 
      q1 ends up being length n1
      d1q1 = d1*q1 is truncated to length n1-1
   */
      
   F_mpz_poly_init(d1q1);
   F_mpz_poly_init(q1);
   F_mpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, d1); 
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length n1+n2-1
   */  
   
   F_mpz_poly_init(d2q1);
   F_mpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length n1+n2-1
   */
   
   F_mpz_poly_init2(dq1, FLINT_MAX(d1q1->length + n2, d2q1->length));
   F_mpz_poly_left_shift(dq1, d1q1, n2);
   F_mpz_poly_clear(d1q1);
   _F_mpz_poly_add(dq1, dq1, d2q1);
   F_mpz_poly_clear(d2q1);
   
   /*
      Compute t = a1*x^(n1+n2-1) + a2*x^(n1-1) - dq1
      which has length 2*n1+n2-1, but we are not interested 
      in up to the first n1 coefficients, so it has 
      effective length n1+n2-1
   */
   
   F_mpz_poly_init2(t, FLINT_MAX(A->length - n2, dq1->length));
   F_mpz_poly_right_shift(t, A, n2);
   _F_mpz_poly_sub(t, t, dq1);
   F_mpz_poly_truncate(t, B->length - 1);
   
   /*
      Compute q2 = t div d1
      It is an n1+n2-1 by n1 division, so
      the length of q2 will be n2
      Also compute d1q2 truncated to length n1-1
   */
   
   F_mpz_poly_init(d1q2);
   F_mpz_poly_init(q2);
   F_mpz_poly_div_divconquer_recursive_low(q2, d1q2, t, d1); 
   F_mpz_poly_clear(t);
      
   /*
      Compute d2q2 = d2*q2 which is of length 
      n1+n2-1
   */
   
   F_mpz_poly_init(d2q2);
   F_mpz_poly_mul(d2q2, d2, q2);  
   
   /*
      Compute dq2 = d1*q2*x^n2 + d2q2
      which is of length n1+n2-1
   */
   
   F_mpz_poly_init2(dq2, FLINT_MAX(d1q2->length + n2, d2q2->length));
   F_mpz_poly_left_shift(dq2, d1q2, n2);
   F_mpz_poly_clear(d1q2);
   _F_mpz_poly_add(dq2, dq2, d2q2);
   F_mpz_poly_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   
   F_mpz_poly_fit_length(Q, FLINT_MAX(q1->length + n2, q2->length));
   F_mpz_poly_left_shift(Q, q1, n2);
   F_mpz_poly_clear(q1);
   _F_mpz_poly_add(Q, Q, q2);
   F_mpz_poly_clear(q2);
   
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length n1+2*n2-1
      We truncate to length B->length - 1
   */
   
   F_mpz_poly_fit_length(BQ, FLINT_MAX(dq1->length + n2, dq2->length));
   F_mpz_poly_left_shift(BQ, dq1, n2);
   _F_mpz_poly_add(BQ, BQ, dq2);
   F_mpz_poly_truncate(BQ, B->length - 1);
   
   F_mpz_poly_clear(dq2);
   F_mpz_poly_clear(dq1);
}

void F_mpz_poly_div_divconquer(F_mpz_poly_t Q, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   if (B->length == 0)
   {
      printf("Exception : divide by zero in F_mpz_poly_div_divconquer\n");
   }
   
   if (A->length < B->length)
   {
      F_mpz_poly_zero(Q);
      
      return;
   }

   // A->length is now >= B->length
    
   ulong crossover = 16;
   
   if (A->length - B->length + 1  <= crossover) 
   {
      F_mpz_poly_div_basecase(Q, A, B);
      
      return;
   }
   
   F_mpz_poly_t d1, d2, d3, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
      
   ulong n1 = (B->length + 1)/2;
   ulong n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 
      d1 is of length n1 and
      d2 of length n2
   */

   _F_mpz_poly_attach_shift(d1, B, n2);
   _F_mpz_poly_attach_truncate(d2, B, n2);
   _F_mpz_poly_attach_shift(d3, B, n1);
   
   if (A->length < 2*B->length - 1)
   {
      /* We convert an unbalanced division into a 2*q by q division */
      
      ulong q = A->length - B->length + 1;
      
      F_mpz_poly_t t_A, t_B;
      _F_mpz_poly_attach_shift(t_B, B, B->length - q);
      _F_mpz_poly_attach_shift(t_A, A, A->length - 2*q + 1);

      F_mpz_poly_div_divconquer(Q, t_A, t_B);

      return; 
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length - 1
      // We call this polynomial p1
      
      ulong shift = A->length - 2*B->length + 1;
      _F_mpz_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = low(d1*q1) is length at most 2*B->length-1
         We discard the lower B->length-1 terms
      */
      
      F_mpz_poly_init(d1q1);
      F_mpz_poly_init(q1);
      
      F_mpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      F_mpz_poly_init(dq1);
      
      F_mpz_poly_fit_length(dq1, d1q1->length + shift);
      F_mpz_poly_left_shift(dq1, d1q1, shift);
      F_mpz_poly_clear(d1q1); 
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      F_mpz_poly_init(t);
      F_mpz_poly_sub(t, A, dq1);
      F_mpz_poly_clear(dq1);
      F_mpz_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      F_mpz_poly_init(q2);
      F_mpz_poly_div_divconquer(q2, t, B); 
      F_mpz_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      F_mpz_poly_fit_length(Q, FLINT_MAX(q1->length + shift, q2->length));
      
      F_mpz_poly_left_shift(Q, q1, shift);
      F_mpz_poly_clear(q1);
      _F_mpz_poly_add(Q, Q, q2);
      F_mpz_poly_clear(q2);
      
      return;
   }
   // We now have A->length == 2*B->length - 1
   
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length n1 and a2 is length n2 
      and a3 is length n1+n2-1 
   */
      
   // Set p1 to a1*x^(n1-1) + other terms
   // It has length 2*n1-1 and is normalised
      
   _F_mpz_poly_attach_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is a 2*n1-1 by n1 division so 
      q1 ends up being length n1
      d1q1 = low(d1*q1) is length n1-1
      Thus we have discarded the leading n1 terms 
   */
      
   F_mpz_poly_init(d1q1);
   F_mpz_poly_init(q1);
      
   F_mpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, d1); 
      
   /* 
      Compute d2q1 = d2*q1 with low n1 - 1 terms zeroed
      d2*q1 is length n1+n2-1 leaving
      n2 non-zero terms to the left
   */  
   
   F_mpz_poly_init(d2q1); 
   F_mpz_poly_mul_trunc_left(d2q1, d2, q1, n1 - 1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length 2*n1+n2-1
   */
   
   F_mpz_poly_init2(dq1, FLINT_MAX(d1q1->length + n2, d2q1->length));
   F_mpz_poly_left_shift(dq1, d1q1, n2);
   F_mpz_poly_clear(d1q1); 
   _F_mpz_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = a1*x^(2*n2-1) + a2*x^(n2-1) - dq1 
      after shifting dq1 to the right by (n1-n2)
      which has length 2*n1+n2-1, but we 
      n1 coefficients, so it has 
      effective length 2*n2-1 with the last n2-1
      coefficients ignored. Thus there are n2 
      significant coefficients
   */
   
   
   F_mpz_poly_init2(t, n1 + 2*n2 - 1);
   F_mpz_poly_right_shift(t, A, n1);
   _F_mpz_poly_attach_shift(temp, dq1, n1 - n2);
   _F_mpz_poly_sub(t, t, temp);
   F_mpz_poly_truncate(t, 2*n2-1);
     
   /*
      Compute q2 = t div d3
      It is a 2*n2-1 by n2 division, so
      the length of q2 will be n2 
   */
   
   F_mpz_poly_init(q2);
   F_mpz_poly_div_divconquer(q2, t, d3); 
   F_mpz_poly_clear(t);  
   F_mpz_poly_clear(dq1);
   F_mpz_poly_clear(d2q1);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   
   F_mpz_poly_fit_length(Q, q1->length + n2);
   F_mpz_poly_left_shift(Q, q1, n2);
   F_mpz_poly_clear(q1);
   _F_mpz_poly_add(Q, Q, q2);
   F_mpz_poly_clear(q2);   
}

/*===============================================================================

	Exact division

================================================================================*/

void F_mpz_poly_div_hensel(F_mpz_poly_t Q, const F_mpz_poly_t A, const ulong a_len, 
                                    const F_mpz_poly_t B, const ulong b_len)
{
   F_mpz_poly_t A_rev, B_rev;

   if (B->length == 0)
   {
      printf("Exception : divide by zero in F_mpz_poly_div_hensel\n");
      abort();
      return;
   }
   
   if (A->length == 0)
   {
      F_mpz_poly_zero(Q);
      return;
   }
   
   F_mpz_poly_init(A_rev);
   F_mpz_poly_init(B_rev);

   ulong q = a_len - b_len + 1;

   F_mpz_poly_reverse(B_rev, B, b_len);
   F_mpz_poly_reverse(A_rev, A, a_len);
   
   F_mpz_poly_div(Q, A_rev, B_rev);
   
   F_mpz_poly_clear(A_rev);
   F_mpz_poly_clear(B_rev);

   F_mpz_poly_reverse(Q, Q, q);
}

void F_mpz_poly_divexact(F_mpz_poly_t Q, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   if (B->length == 0)
   {
      printf("Exception : divide by zero in F_mpz_poly_div_hensel\n");
      abort();
      return;
   }
   
   if (A->length < B->length)
   {
      F_mpz_poly_zero(Q);
      return;
   }

   long q = A->length - B->length + 1;

   F_mpz_poly_fit_length(Q, q);
   Q->length = q;

   ulong a_len = A->length;
   ulong b_len = B->length;

   ulong i;
   
   for (i = 0; ; i++)
   {
      if (!F_mpz_is_zero(B->coeffs + i)) break;
   }

   ulong i2 = i;
   b_len -= i;

   for ( ; ; i++)
   {
      if (!F_mpz_is_zero(A->coeffs + i)) break;
      F_mpz_zero(Q->coeffs + i - i2);      
   }

   a_len -= i;

   q = a_len - b_len + 1;

   F_mpz_poly_t t_A, t_B, t_Q;

   if (q <= 1)
   {
      _F_mpz_poly_attach_shift(t_Q, Q, i - i2);
      _F_mpz_poly_attach_shift(t_A, A, i);
      _F_mpz_poly_attach_shift(t_B, B, i2);
      
      t_Q->alloc = q;

      F_mpz_poly_div(t_Q, t_A, t_B);
      
      return;
   }

   ulong q2 = (q + 1)/2;
   
   t_Q->alloc = q2;
   
   ulong q2b = FLINT_MIN(b_len, q2);
   _F_mpz_poly_attach_shift(t_Q, Q, q + i - i2 - q2);
   _F_mpz_poly_attach_shift(t_A, A, A->length - q2 - q2b + 1);
   _F_mpz_poly_attach_shift(t_B, B, B->length - q2b);

   F_mpz_poly_div(t_Q, t_A, t_B);
   
   q2 = q - q2;
   q2b = FLINT_MIN(b_len, q2);
   t_Q->coeffs = Q->coeffs + i - i2;
   t_Q->length = q2;
   t_A->coeffs = A->coeffs + i;
   t_A->length = q2 + q2b - 1;
   t_B->coeffs = B->coeffs + i2;
   t_B->length = q2b;
   
   F_mpz_poly_div_hensel(t_Q, t_A, q2 + q2b - 1, t_B, q2b);
}

/*===============================================================================

	Pseudo division

================================================================================*/

void F_mpz_poly_pseudo_divrem_basecase(F_mpz_poly_t Q, F_mpz_poly_t R, 
                            ulong * d, const F_mpz_poly_t A, const F_mpz_poly_t B)
{
   F_mpz_poly_t qB;
   
   if (B->length == 0)
   {
      printf("Exception : Divide by zero in F_mpz_poly_pseudo_divrem_basecase.\n");
      abort();
   }
  
   F_mpz * coeffs_A = A->coeffs;
   F_mpz * coeffs_B = B->coeffs;
   F_mpz * B_lead = coeffs_B + B->length - 1; 
   F_mpz * coeff_Q;
   F_mpz * coeff_R;
   F_mpz * coeffs_R;
   int scale;
   
   long m = A->length;
   long n = B->length;
   long q = m - n + 1;

   ulong size_B_lead = F_mpz_size(B_lead);
      
   F_mpz_poly_struct Rs;
   
   int want_rem = 1;
   if (R == NULL)
   {
      want_rem = 0;
      R = &Rs;
      F_mpz_poly_init(R);
      
   } else F_mpz_poly_set(R, A);

   coeffs_R = R->coeffs;
   
   *d = 0;
   
   if ((long) R->length >= (long) B->length)
   {
      F_mpz_poly_fit_length(Q, R->length - B->length + 1);
      
      ulong i;
      for (i = 0; i < R->length - B->length + 1; i++) 
         F_mpz_zero(Q->coeffs + i);

      Q->length = R->length - B->length+1;
   } else 
   {
      F_mpz_poly_zero(Q);
      return;
   }
   
   if (!want_rem) 
   {
      F_mpz_poly_fit_length(R, A->length);
      ulong i;
      for (i = A->length - q; i < A->length; i++)
          F_mpz_set(R->coeffs + i, A->coeffs + i);
   }

   F_mpz_poly_t Bm1;
   ulong Bsub_length = B->length;
   _F_mpz_poly_attach_truncate(Bm1, B, B->length - 1);

   coeff_R = coeffs_R + R->length - 1;

   F_mpz_t rem;
   F_mpz_init(rem);
   
   while ((long) R->length >= (long) B->length)
   {
      coeff_Q = Q->coeffs + R->length - Bsub_length;
          
      if (F_mpz_cmpabs(coeff_R, B_lead) >= 0)
      {
         F_mpz_fdiv_qr(coeff_Q, rem, coeff_R, B_lead);

      } else
      {
         F_mpz_zero(coeff_Q);
         if (F_mpz_is_zero(coeff_R)) F_mpz_zero(rem);
         else F_mpz_set_ui(rem, 1);
      }
      
      if (F_mpz_is_zero(rem))
      {
         scale = 0; 
      } else 
      {   
         F_mpz_poly_scalar_mul(Q, Q, B_lead);
         coeff_Q = Q->coeffs + R->length - B->length;
         F_mpz_set(coeff_Q, coeff_R);
         scale = 1;
         (*d)++;
      }
           
      if (B->length > 1)
      {
         F_mpz_poly_init2(qB, Bsub_length - 1);
         F_mpz_poly_scalar_mul(qB, Bm1, coeff_Q); 
      } else F_mpz_poly_init(qB); 
      
      if (scale)
      {
         coeffs_R = R->coeffs;
         F_mpz_poly_scalar_mul(R, R, B_lead);
      } else if (B->length > 1)
      {
         coeffs_R = R->coeffs;
      }
      
      F_mpz_poly_t R_sub;
      R_sub->coeffs = coeffs_R + R->length - Bsub_length;
      R_sub->length = Bsub_length - 1;
      
      if (B->length > 1)
      {
         _F_mpz_poly_sub(R_sub, R_sub, qB);
      }
      F_mpz_poly_clear(qB);
      
      ulong old_len = R->length;
      F_mpz_zero(R_sub->coeffs + Bsub_length - 1);
      
      _F_mpz_poly_normalise(R);
      coeff_R = coeffs_R + R->length - 1;

      if ((!want_rem) && (Bsub_length + B->length >= R->length + 1))
      {
         ulong diff = old_len - R->length;
         Bm1->coeffs += diff;
         Bm1->length -= diff;
         Bsub_length -= diff;
      }
   }
  
   F_mpz_clear(rem);
}

/*===============================================================================

   New Naive Standard Functions (without test code and written by Andy)

================================================================================*/

void F_mpz_poly_scalar_div_exact(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t d)
{

//check for d=+/-1?
   if (F_mpz_is_zero(d)){
      printf("FLINT Exception: Division by zero\n");
      abort();
   }

   if (F_mpz_is_one(d)){
      F_mpz_poly_set(res, f);
      return;
   }

   if (F_mpz_is_m1(d)){
      F_mpz_poly_neg(res, f);
      return;
   }

   F_mpz_poly_fit_length(res, f->length);

   res->length = f->length;

   long i;
   for (i = 0; i < f->length; i++){
      F_mpz_divexact(res->coeffs + i, f->coeffs + i, d);
   }
}

void F_mpz_poly_smod(F_mpz_poly_t res, F_mpz_poly_t f, F_mpz_t p)
{

   if (F_mpz_is_zero(p)){
      printf("FLINT Exception: Division by zero\n");
      abort();
   }

   if (F_mpz_is_one(p)){
      F_mpz_poly_zero(res);
      return;
   }

   F_mpz_t pdiv2;
   F_mpz_init(pdiv2);

   F_mpz_div_2exp(pdiv2, p, 1);

   F_mpz_poly_fit_length(res, f->length);

   res->length = f->length;

   long i;
   for (i = 0; i < f->length; i++){
      F_mpz_mod(f->coeffs + i, f->coeffs + i, p);

      if ( F_mpz_cmp( f->coeffs + i, pdiv2) > 0){
         F_mpz_sub(res->coeffs+i, f->coeffs + i, p);
      }
      else{
         F_mpz_set(res->coeffs + i, f->coeffs + i);
      }

   }

   _F_mpz_poly_normalise(res);

   F_mpz_clear(pdiv2);

}

void F_mpz_poly_derivative(F_mpz_poly_t der, F_mpz_poly_t poly)
{
	if (poly->length <= 1)
	{
		F_mpz_poly_zero(der);
		return;
	}
	

   F_mpz_poly_fit_length(der, poly->length - 1);

	der->length = poly->length - 1;

   ulong i;
   for (i = 0; i < poly->length - 1; i++)
	{
		F_mpz_mul_ui(der->coeffs + i, poly->coeffs + i + 1, i + 1);
	}

}

void F_mpz_poly_content(F_mpz_t c, const F_mpz_poly_t poly)
{
   unsigned long length = poly->length;

   if (length == 0) 
   {
      F_mpz_set_ui(c, 0L);
      return;
   }
   
   if (length == 1)
   {
      F_mpz_set(c, poly->coeffs);
//      if ((long) c[0] < 0L) c[0] = -c[0];
      return;
   }
   
   F_mpz_t coeff;
   F_mpz_init(coeff);

   F_mpz_set(coeff, poly->coeffs + length - 1);
   F_mpz_set(c, coeff);
   
   long i;
   for (i = length - 2; (i >= 0L) && !F_mpz_is_one(c); i--)
   {
      F_mpz_set(coeff, poly->coeffs + i);
      if (!F_mpz_is_zero(coeff))
         F_mpz_gcd(c, c, coeff);
   }

//   if (F_mpz_sgn(poly->coeffs + length -1) == -1)
//      F_mpz_neg(c, c);

   F_mpz_clear(coeff);

}

double F_mpz_poly_eval_horner_d(F_mpz_poly_t poly, double val){

   ulong n = poly->length;

   long exp;
   double temp;

   temp = F_mpz_get_d_2exp(&exp, poly->coeffs + n - 1);
   temp = temp*pow(2, exp);

   double ans = temp;

   long i;
   for (i = n - 2; i >= 0L; i--)
   {
      ans = ans * val;

      temp = F_mpz_get_d_2exp(&exp, poly->coeffs + i);
      temp = temp*pow(2, exp);

      ans = ans + temp;
   }
   return ans;
}

double F_mpz_poly_eval_horner_d_2exp(long * exp, F_mpz_poly_t poly, double val){

   ulong vbits = round( abs( log(val) / log(2.0) ) );
   long size_p = F_mpz_poly_max_bits(poly);
   ulong n = poly->length;
   ulong prec=(vbits*n) + FLINT_ABS(size_p) + 1;  
   mpf_set_default_prec(prec);

   mpz_t z_coeff_p;
   mpf_t f_coeff_p;
   mpf_t fval, output;

   mpf_init(output); 
   mpz_init(z_coeff_p);

   F_mpz_t coeff_p;
   F_mpz_init(coeff_p);
   F_mpz_set(coeff_p, poly->coeffs + n - 1);
   F_mpz_get_mpz(z_coeff_p,coeff_p);
   mpf_set_z(output,z_coeff_p); //Set output to top coeff

   mpz_clear(z_coeff_p);

   mpf_init(fval);

   mpf_set_d(fval,val);//set fval to mpf from the double val

   long i;
   for (i = n - 2; i >= 0L; i--)
   {

      mpf_mul(output,output,fval);

      F_mpz_set(coeff_p, poly->coeffs + i);

      mpz_init(z_coeff_p);
      mpf_init(f_coeff_p);

      F_mpz_get_mpz(z_coeff_p,coeff_p);//convert coeff from fmpz to mpz to mpf
      mpf_set_z(f_coeff_p,z_coeff_p);

      mpf_add(output,output,f_coeff_p);//add coeff to output then repeat

      mpf_clear(f_coeff_p);
      mpz_clear(z_coeff_p);
   }
   double res = mpf_get_d_2exp( exp, output);

   mpf_clear(output);
   mpf_clear(fval);

   return res;
}

void F_mpz_poly_scalar_abs(F_mpz_poly_t output, F_mpz_poly_t input){

   if (input == output){
      for( long i = 0; i < input->length; i++){
         F_mpz_abs(output->coeffs + i, input->coeffs + i);
      }
   }
   else{
      F_mpz_poly_fit_length(output, input->length);
      for( long i = 0; i < input->length; i++){
         F_mpz_abs(output->coeffs + i, input->coeffs + i);
      }
      output-> length = input-> length;      
   }
   return;
}

/*===========================================

   Some temporary global timing variables

==========================================*/

clock_t check_if_solve_start, check_if_solve_stop;
clock_t check_if_solve_total = 0;

clock_t local_factor_start, local_factor_stop;
clock_t local_factor_total = 0;

clock_t lll_start, lll_stop;
clock_t lll_total = 0;

clock_t hensel_start, hensel_stop;
clock_t hensel_total = 0;

clock_t linear_alg_start, linear_alg_stop;
clock_t linear_alg_total = 0;

clock_t cld_data_start, cld_data_stop;
clock_t cld_data_total = 0;


/*===========================================================================

   stupid F_mpz_poly functions which just wrap fmpz_poly functions

============================================================================*/

void F_mpz_poly_gcd(F_mpz_poly_t d, F_mpz_poly_t f, F_mpz_poly_t g){

   mpz_poly_t mpz_d, mpz_f, mpz_g;
   mpz_poly_init(mpz_d);
   mpz_poly_init(mpz_f);
   mpz_poly_init(mpz_g);

   fmpz_poly_t fmpz_d, fmpz_f, fmpz_g;
   fmpz_poly_init(fmpz_d);
   fmpz_poly_init(fmpz_f);
   fmpz_poly_init(fmpz_g);

   F_mpz_poly_to_mpz_poly(mpz_f, f);
   F_mpz_poly_to_mpz_poly(mpz_g, g);

   mpz_poly_to_fmpz_poly(fmpz_f, mpz_f);
   mpz_poly_to_fmpz_poly(fmpz_g, mpz_g);

   fmpz_poly_gcd(fmpz_d, fmpz_f, fmpz_g);

   fmpz_poly_to_mpz_poly(mpz_d, fmpz_d);

   mpz_poly_to_F_mpz_poly(d, mpz_d);
        
   fmpz_poly_clear(fmpz_d);
   fmpz_poly_clear(fmpz_f);
   fmpz_poly_clear(fmpz_g);

   mpz_poly_clear(mpz_d);
   mpz_poly_clear(mpz_f);
   mpz_poly_clear(mpz_g);
}

/*===========================================================================

   New Material for FLINT, computing fast/tight bounds for CLDs
      CLDs:= Coefficients of Logarithmic Derivatives.  f*g'/g

============================================================================*/

int _d_2exp_comp(double a, long ap, double b, long bp){
//assumes that if ap != 0 (or bp != 0) then a (or b) is in [1/2,1)
   if (ap == 0){
      if (bp == 0){
         if (a > 2*b)
            return 2;
         else if (b > 2*a)
            return -2;
         else if (a >= b)
            return 1;
         else
            return -1;
      }
// now we know that a is the number but b is in 1/2,1 with power bp
      long temp_ap = 1L + (long)log2(a);
      if (temp_ap >= bp + 2)
         return 2;
      else if (bp >= temp_ap + 2)
         return -2;
      else{
         double ta,tb;
         ta = a/4;
         tb = b*pow(2, (double)bp - 2);
         return _d_2exp_comp(ta, 0, tb, 0);
      }
   }
   else if (bp == 0){
      if (ap == 0){
         if (a > 2*b)
            return 2;
         else if (b > 2*a)
            return -2;
         else if (a >= b)
            return 1;
         else
            return -1;
      }
// now we know that b is the number but a is in 1/2,1 with power ap
      long temp_bp = 1L + (long)log2(b);
      if (ap >= temp_bp + 2)
         return 2;
      else if (temp_bp >= ap + 2)
         return -2;
      else{
         double ta,tb;
         ta = a*pow(2, (double)ap - 2);
         tb = b/4;
         return _d_2exp_comp(ta, 0, tb, 0);
      }
   }
   else{
// now we know that both are in 1/2, 1 with powers
      if (ap >= bp + 2)
         return 2;
      else if (bp >= ap + 2)
         return -2;
      else{
         double ta,tb;
         ta = a*pow(2, (double)ap - (double)bp);
         tb = b;
         return _d_2exp_comp(ta, 0, tb, 0);
      }
   }
}

void F_mpz_poly_CLD_bound(F_mpz_t res, F_mpz_poly_t f, ulong N){
   if ((N < 0) || (N >= f->length - 1)){
      printf("bad input\n");
      return;
   }
//Coded with n = c+1 and decided that the user would rather give 0,1,2,3 instead of 1,2,3,4
   ulong n = N + 1;
   F_mpz_poly_t low_f, up_f;
   F_mpz_poly_init(low_f);
   F_mpz_poly_init(up_f);
   F_mpz_poly_set(low_f, f);
   F_mpz_poly_truncate(low_f, n);
   F_mpz_poly_scalar_abs(low_f, low_f);

   F_mpz_poly_right_shift(up_f, f, n);
   F_mpz_poly_left_shift(up_f, up_f, n);
   F_mpz_poly_scalar_abs(up_f, up_f);
//Need to take scalar_abs of F_mpz_poly's here.
   double rpower = 0;
   double rshift = 1;
   double r = pow(2,rpower);
   double top_eval;
   double bottom_eval;
//OK right now we have a fast and loose double version, which suffices for the moment.  It is completely forseeable that a large polynomial input will require an F_mpz_t bound... I don't need a super precise output so I want to do the exponent trick to allow large numbers with only double bits of precision.  So need to make an F_mpz_poly_eval_horner_fast_d_2exp and some kind of an exponent size check on the flip side.  The eval program can use doubles and just track the exponents... in the meantime for small input polys this function works.
   long top_exp;
   long bot_exp;
   int dir = 1;
   double ans;
   int good_enough = 0;
   long size_p = F_mpz_poly_max_bits(f);
   ulong hn;// = poly->length;
   ulong vbits;// = round( abs( log(val) / log(2.0) ) );
   ulong prec;// =(vbits*n) + FLINT_ABS(size_p) + 1;
   int too_much = 0;
   while (good_enough == 0L){
      hn = up_f->length;
      vbits = round( abs( log2(r) ) );
      prec = (vbits*hn) + FLINT_ABS(size_p) + 1;
      //this is a rough bound for the number of bits of the answer...
      if ((prec > 950) || (too_much = 1)){
         top_eval = F_mpz_poly_eval_horner_d_2exp( &top_exp, up_f, r);
         // maybe I'll deal with this on it's own.  top_eval = top_eval*pow(2, top_exp);
      }      
      else{
         //Here we knew all along that doubles were good enough
         top_eval = F_mpz_poly_eval_horner_d( up_f, r);
         top_exp = 0;
      }
      hn = low_f->length;
      prec = (vbits*hn) + FLINT_ABS(size_p) + 1;
      if ((prec > 950) || (too_much = 1)){
         bottom_eval = F_mpz_poly_eval_horner_d_2exp( &bot_exp, low_f, r);
//         bottom_eval = bottom_eval*pow(2, bot_exp);
      }      
      else{
         bottom_eval = F_mpz_poly_eval_horner_d(low_f, r);
         bot_exp = 0;
      }
      if ((top_exp == 0) && (bot_exp == 0)){
         if ( 2*(bottom_eval) < (top_eval) ){
            if (dir == 1)
               rshift = rshift/2;
            dir = -1;
            rpower = rpower - rshift;
            r = pow(2, rpower);
         }
         else if (  (bottom_eval) > 2*(top_eval) ){
            if (dir == -1)
               rshift = rshift/2;
            dir = 1;
            rpower = rpower + rshift;
            r = pow(2, rpower);
         }
         else{
            good_enough = 1;

            if (isinf(top_eval) || isinf(bottom_eval))
            {
               too_much = 1; 
               good_enough = 0;
            }
            else
            {
               if (top_eval > bottom_eval)
                  ans = top_eval;
               else
                  ans = bottom_eval;
               ans = ans / pow(r, n);
               ans = ans*(f->length - 1);
               mpz_t temp;
               mpz_init(temp);
               mpz_set_d(temp, ans);
               F_mpz_set_mpz(res, temp);
               mpz_clear(temp);
            }
         }

      }
      else{
//here is trouble land, coeffs too big for doubles to handle.
// _d_2exp_comp should give 2 when 2*bottom < top and -2 when 2*top < bottom and 1 when bottom <= top and -1 when bottom > top
         int test_me = _d_2exp_comp(top_eval, top_exp, bottom_eval, bot_exp);

         if ( test_me == 2 ){
            if (dir == 1)
               rshift = rshift/2;
            dir = -1;
            rpower = rpower - rshift;
            r = pow(2, rpower);
         }
         else if ( test_me == -2 ){
            if (dir == -1)
               rshift = rshift/2;
            dir = 1;
            rpower = rpower + rshift;
            r = pow(2, rpower);
         }
         else{
            good_enough = 1;
            if ((test_me == 1L)){
// F_mpz_set_d_2exp and adjust and stuff using top_eval, top_exp 
               ans = top_eval;
               ans = ans / pow(r, n);
               ans = ans * f->length - 1;
               F_mpz_set_d_2exp(res, ans, top_exp);
               F_mpz_poly_clear(low_f);
               F_mpz_poly_clear(up_f);
               return;
            }
            else{
// F_mpz_set_d_2exp and adjust and junk using bottom_eval, bot_exp
               ans = bottom_eval;
               ans = ans / pow(r, n);
               ans = ans * f->length - 1;
               F_mpz_set_d_2exp(res, ans, bot_exp);
               F_mpz_poly_clear(low_f);
               F_mpz_poly_clear(up_f);
               return;
            }
         }
      }
   }
   F_mpz_poly_clear(low_f);
   F_mpz_poly_clear(up_f);
   return;
}

/*============================================================================

   Naive '_modp' ( := Large moduli ) F_mpz_poly functions

============================================================================*/

void F_mpz_poly_rem_modp_naive(F_mpz_poly_t R, F_mpz_poly_t A, F_mpz_poly_t B, F_mpz_t p){

   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      F_mpz_poly_set(R, A);      
      return;
   }

   F_mpz_poly_set(R,A);
   long coeff = A->length - 1;

   F_mpz_poly_t pre_inv_B, Bm1, qB;
   F_mpz_poly_init2(Bm1, B->length);

   F_mpz_poly_set(Bm1, B);

   F_mpz_t B_lead_inv;
   F_mpz_init(B_lead_inv);
   F_mpz_t coeff_Q;
   F_mpz_init( coeff_Q);

   F_mpz_invert(B_lead_inv, B->coeffs + (B->length - 1), p);

   long R_length;

   F_mpz_t temp;
   F_mpz_init( temp );

   while (coeff >= (long) Bm1->length - 1){
      while( (coeff >= (long) Bm1->length - 1) && F_mpz_is_zero( R->coeffs + coeff ) ){
         coeff--;
      }
      if (coeff >= (long) B->length - 1){
         F_mpz_mul2(temp, R->coeffs + coeff, B_lead_inv);
         F_mpz_mod(coeff_Q, temp, p);

         F_mpz_poly_init(qB);
         F_mpz_poly_fit_length(qB, Bm1->length);

         qB->length = Bm1->length;
         long i;
         for (i = 0; i < Bm1->length; i++){
            F_mpz_mul2(temp, Bm1->coeffs + i , coeff_Q );
            F_mpz_mod(qB->coeffs + i, temp, p);
         }
//for each coeff do mulmod by coeff_Q and write to temp_B
         F_mpz_poly_left_shift(qB, qB, coeff - B->length + 1);

         R_length = coeff;
         F_mpz_poly_sub(R, R, qB);
         R->length = R_length;

         F_mpz_poly_clear(qB);

         for (i = 0; i < R->length; i++){
            F_mpz_mod(R->coeffs + i, R->coeffs + i, p);
         }
         //subtract a shifted temp_B from R
         coeff--;
      }
   }
   R->length = B->length - 1;
   _F_mpz_poly_normalise(R);

   F_mpz_clear(B_lead_inv);

   F_mpz_poly_clear(Bm1);
   F_mpz_clear(coeff_Q);
   F_mpz_clear(temp);
}

void F_mpz_poly_mulmod_modp_naive(F_mpz_poly_t R, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t B, F_mpz_t p){
//Want to compute f*g mod B where all polynomials are considered modulo an F_mpz modulus p.  Need that the leading coeff of B is invertible mod p.
   if (B->length == 0){
      printf("FLINT Exception: Divide by zero\n");
      abort();
   }
   if ( (B->length ==1) || (f->length == 0) || (g->length == 0) ){
      F_mpz_poly_zero(R);
      return;
   }
   F_mpz_poly_t prod;
   F_mpz_poly_init(prod);
   F_mpz_poly_mul(prod, f, g);
   F_mpz_poly_rem_modp_naive(R, prod, B, p);
   F_mpz_poly_clear(prod);
}

void F_mpz_poly_div_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_t P, ulong n){
//assuming that g divides f mod P find the bottom n coeffs of f/g mod P.  (truncate now or no?) 
   if (g->length > f->length){
      if (n < f->length){
         for(ulong i = 0; i < n; i++)
            F_mpz_smod(res[i], f->coeffs + i, P);
         return;
      }
      else{
         for(ulong i = 0; i < f->length; i++)
            F_mpz_smod(res[i], f->coeffs + i, P);
         return;
      }
   }
   F_mpz_t temp, tc_inv;
   F_mpz_init(temp);
   F_mpz_init(tc_inv);

   F_mpz_poly_t t_f, t_g;
   F_mpz_poly_init(t_f);
   F_mpz_poly_init(t_g);

   F_mpz_poly_set(t_f, f);
   F_mpz_poly_set(t_g, g);
   F_mpz_poly_truncate(t_f, n);
   F_mpz_poly_truncate(t_g, n);
//now we have t_f, t_g truncated to the bottom n terms for speed reasons.
   F_mpz_set(temp, t_g->coeffs);
   F_mpz_invert(tc_inv, temp, P);

   if (F_mpz_is_zero(tc_inv)){
//Potential math problem here so we'll just get the answers the old fashioned way, hope this is rare, I should talk to somebody...
      F_mpz_poly_div(t_f, f, g);
      F_mpz_poly_smod(t_f, t_f, P);      

      if (n < t_f->length)
         for(ulong i = 0; i < n; i++){
            F_mpz_set(res[i], t_f->coeffs + i);
         }
      else
         for(ulong i = 0; i < t_f->length; i++){
            F_mpz_set(res[i], t_f->coeffs + i);
         }
      F_mpz_poly_clear(t_f);
      F_mpz_poly_clear(t_g);
   
      F_mpz_clear(tc_inv);
      F_mpz_clear(temp);
      return;
   }
   F_mpz_poly_t tempg;
   F_mpz_poly_init(tempg);

   ulong quo_length = n;

   if (n > f->length)
      quo_length = f->length - g->length + 1;

   for(ulong i = 0; (i < n) && (i < quo_length) ; i++){
//  This is a number which would cancel out the lowest term of f mod P
      F_mpz_mul2(temp, t_f->coeffs, tc_inv);
      F_mpz_smod(res[i], temp, P);

      F_mpz_poly_right_shift(tempg, t_g, 1UL);
      F_mpz_poly_right_shift(t_f, t_f, 1UL);

      F_mpz_poly_scalar_mul(tempg, tempg, temp);
      F_mpz_poly_sub(t_f, t_f, tempg);

      F_mpz_poly_smod(t_f, t_f, P);

      F_mpz_poly_truncate(t_f, n - i -1);
   }
   F_mpz_poly_clear(tempg);
   F_mpz_poly_clear(t_f);
   F_mpz_poly_clear(t_g);

   F_mpz_clear(tc_inv);
   F_mpz_clear(temp);
   return;
}

void F_mpz_poly_div_upper_trunc_modp( F_mpz_t *res, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_t P, ulong n){
//assuming that g divides f mod P find the top n coeffs of f/g mod P using only the top n coeffs of f and g.
   if (g->length > f->length){
      if (n < f->length){
         for(ulong i = 0; i < n; i++)
            F_mpz_smod(res[i], f->coeffs + i, P);
         return;
      }
      else{
         for(ulong i = 0; i < f->length; i++)
            F_mpz_smod(res[i], f->coeffs + i, P);
         return;
      }
   }
   F_mpz_t temp, lc_inv;
   F_mpz_init(temp);
   F_mpz_init(lc_inv);

   F_mpz_poly_t t_f, t_g;
   F_mpz_poly_init(t_f);
   F_mpz_poly_init(t_g);

   F_mpz_poly_set(t_f, f);
   F_mpz_poly_set(t_g, g);
   if (n < t_f->length){
      F_mpz_poly_right_shift(t_f, t_f, t_f->length - n);
   }
   if (n < t_g->length){
      F_mpz_poly_right_shift(t_g, t_g, t_g->length - n);
   }
//now we have t_f, t_g truncated to the bottom n terms.
   F_mpz_set(temp, t_g->coeffs + t_g->length - 1);
   F_mpz_invert(lc_inv, temp, P);

   if (F_mpz_is_zero(lc_inv)){
//Potential math problem here so we'll just get the answers the old fashioned way, hope this is rare, I should talk to somebody...
      F_mpz_poly_div(t_f, f, g);
      F_mpz_poly_smod(t_f, t_f, P);      

      if (n < t_f->length)
         for(ulong i = 0; i < n; i++){
            F_mpz_set(res[i], t_f->coeffs + t_f->length - 1 - i);
         }
      else
         for(ulong i = 0; i < t_f->length; i++){
            F_mpz_set(res[i], t_f->coeffs + t_f->length -1 - i);
         }
      F_mpz_poly_clear(t_f);
      F_mpz_poly_clear(t_g);   
      F_mpz_clear(lc_inv);
      F_mpz_clear(temp);
      return;
   }
   F_mpz_poly_t tempg;
   F_mpz_poly_init(tempg);

   ulong quo_length = n;
   if (n > f->length)
      quo_length = f->length - g->length + 1;
   for(ulong i = 0; (i < n) && (i < quo_length); i++){
      F_mpz_mul2(temp, t_f->coeffs + t_f->length - 1, lc_inv);
      F_mpz_smod(res[i], temp, P);

      long fg_diff;
      fg_diff = t_f->length - t_g->length;

      if (fg_diff>=0)
         F_mpz_poly_left_shift(tempg, t_g, fg_diff);
      else
         F_mpz_poly_right_shift(tempg, t_g, -fg_diff);

      F_mpz_poly_truncate(t_f, t_f->length - 1);
      F_mpz_poly_truncate(tempg, tempg->length - 1);
      F_mpz_poly_scalar_mul(tempg, tempg, res[i]);
      F_mpz_poly_sub(t_f, t_f, tempg);
      F_mpz_poly_smod(t_f, t_f, P);
   }
   F_mpz_poly_clear(tempg);
   F_mpz_poly_clear(t_f);
   F_mpz_poly_clear(t_g);

   F_mpz_clear(lc_inv);
   F_mpz_clear(temp);
   return;
}

/*============================================================================

   Square-Free Factorization

============================================================================*/

void F_mpz_poly_squarefree(F_mpz_poly_factor_t fac, F_mpz_t content, F_mpz_poly_t F)
{

   F_mpz_poly_content(content, F);

   F_mpz_poly_t f;
   F_mpz_poly_init(f);

   F_mpz_poly_scalar_div_exact(f, F, content);

   F_mpz_poly_factor_clear(fac);
   F_mpz_poly_factor_init(fac);

   if (f->length == 1)
      return;

   F_mpz_poly_t d, v, w, s, t1;
   F_mpz_poly_init(d);
   F_mpz_poly_init(v);
   F_mpz_poly_init(w);
   F_mpz_poly_init(s);
   F_mpz_poly_init(t1);

   F_mpz_poly_derivative(t1, f);
   F_mpz_poly_gcd(d, f, t1);

   if (d->length == 1){
      F_mpz_poly_factor_insert(fac, f, 1);

      F_mpz_poly_clear(d);
      F_mpz_poly_clear(v);
      F_mpz_poly_clear(w);
      F_mpz_poly_clear(s);
      F_mpz_poly_clear(t1);
      F_mpz_poly_clear(f);
      return;
   }

   F_mpz_poly_div(v, f, d);
   F_mpz_poly_div(w, t1, d);
   long i = 0;

   for ( ; ;){

      i = i + 1;
      F_mpz_poly_derivative(t1, v);
      F_mpz_poly_sub(s, w, t1);
      if (s->length == 0){
         if (v->length > 1)
            F_mpz_poly_factor_insert(fac, v, i);
         F_mpz_poly_clear(d);
         F_mpz_poly_clear(v);
         F_mpz_poly_clear(w);
         F_mpz_poly_clear(s);
         F_mpz_poly_clear(t1);
         F_mpz_poly_clear(f);
         return;
      }
      F_mpz_poly_gcd(d, v, s);
      F_mpz_poly_div(v, v, d);
      F_mpz_poly_div(w, s, d);

      if (d->length > 1)
         F_mpz_poly_factor_insert(fac, d, i);
   } 
   F_mpz_poly_clear(f);
   F_mpz_poly_clear(d);
   F_mpz_poly_clear(v);
   F_mpz_poly_clear(w);
   F_mpz_poly_clear(s);
   F_mpz_poly_clear(t1);
   return;
}


/*****************************************************************************

   NTL based Hensel Lifting procedures using F_mpz_mod_polys

*****************************************************************************/

void F_mpz_poly_build_hensel_tree(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, zmod_poly_factor_t fac)
{

   long r;
   unsigned long p;

   r = fac->num_factors;
   p = (fac->factors[0])->p;

   long i, j, s;
   long minp, mind;
   long tmp;

   zmod_poly_t V[2*r-2];
   zmod_poly_t W[2*r-2];

   for (i = 0; i < 2*r-2; i++)
      zmod_poly_init(V[i],p);

   for (i = 0; i < 2*r-2; i++)
      zmod_poly_init(W[i],p);

/* We will have five arrays: a v of fmpz_polys and a V of zmod_polys also a w and a W and link.  Here's the idea, we will sort each leaf and node of a factor tree by degree, in fact choosing to multiply the two smallest factors, then the next two smallest (factors or products) until a tree is made.  The tree will be stored in the v's.  The first two elements of v will be the smallest modular factors, the last two elements of v will multiply to form F itself.  Since v will be rearranging the original factors we will need to be able to recover the original order.  For this we use link which has nonnegative even numbers and negative numbers.  Link is an array of longs which aligns with V/v  if link has a negative number in spot j that means V[j] is an original modular factor which has been lifted, if link[j] is a nonnegative even number then V[j] stores a product of the two entries at V[ link[j] ] and V[ link[j]+1 ].  W/w plays the role of the extended GCD, at V[0], V[2], V[4], etc we have a new product, W[0], W[2], W[4], etc are the XGCD compliments of the V's.  So V[0]*W[0]+V[1]*W[1] = 1 mod p^(something)  these will be lifted along with the entries in V.  It's not enough to just lift each factor we have to lift the entire tree and the tree of XGCD inverses.*/

   for (i = 0; i < r; i++)
   {
      zmod_poly_set(V[i], fac->factors[i]);
      link[i] = -(i+1);
   }

   for (j = 0; j < 2*r - 4; j += 2)
   {
      minp = j;
      mind = zmod_poly_degree(V[j]);
      for (s = j+1; s < i; s++)
      {
         if (zmod_poly_degree(V[s]) < mind)
         {
            minp = s;
            mind = zmod_poly_degree(V[s]);
         }
      }
      zmod_poly_swap(V[j],V[minp]);

      tmp = link[j];
      link[j] = link[minp];
      link[minp] = tmp; 
      //swap link[j] and V[minp]

      minp = j+1;
      mind = zmod_poly_degree(V[j+1]);

      for ( s = j+2; s < i; s++)
      {
         if (zmod_poly_degree(V[s]) < mind )
         {
            minp = s;
            mind = zmod_poly_degree(V[s]);
         }
      }

      zmod_poly_swap(V[j+1],V[minp]);

      tmp = link[j+1];
      link[j+1] = link[minp];
      link[minp] = tmp; 
      //swap link[j+1] and V[minp]

      zmod_poly_mul(V[i], V[j], V[j+1]);
      link[i] = j;
      i++;
   }
   zmod_poly_t d;
   zmod_poly_init(d, p);
   for (j = 0; j < 2*r - 2; j += 2)
   {
      zmod_poly_xgcd(d, W[j], W[j+1], V[j], V[j+1]);
      //Make a check for d!=1
   }
   for (j = 0; j < 2*r-2; j++)
   {
      zmod_poly_to_F_mpz_poly(v[j], V[j]);
      zmod_poly_to_F_mpz_poly(w[j], W[j]);
   }
   for (i = 0; i < 2*r-2; i++)
      zmod_poly_clear(V[i]);
   for (i = 0; i < 2*r-2; i++)
      zmod_poly_clear(W[i]);
   zmod_poly_clear(d);
}

void F_mpz_poly_hensel_lift(F_mpz_poly_t Gout, F_mpz_poly_t Hout, F_mpz_poly_t Aout, F_mpz_poly_t Bout, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t h, F_mpz_poly_t a, F_mpz_poly_t b, F_mpz_t p, F_mpz_t p1, F_mpz_t big_P)
{

   F_mpz_poly_t c, g1, h1, G, H, A, B;
   F_mpz_poly_init(c);
   F_mpz_poly_init(g1);
   F_mpz_poly_init(h1);
   F_mpz_poly_init(G);
   F_mpz_poly_init(H);
   F_mpz_poly_init(A);
   F_mpz_poly_init(B);

   F_mpz_poly_mul(c, g, h);
   F_mpz_poly_sub(c, f, c);
   F_mpz_poly_scalar_div_exact(c, c, p);

   //Make a check that c is divisible by p
   F_mpz_mod_poly_t cc, gg, hh, aa, bb, tt, gg1, hh1;
   F_mpz_mod_poly_init(cc, p1);
   F_mpz_mod_poly_init(gg, p1);
   F_mpz_mod_poly_init(hh, p1);
   F_mpz_mod_poly_init(aa, p1);
   F_mpz_mod_poly_init(bb, p1);
   F_mpz_mod_poly_init(tt, p1);
   F_mpz_mod_poly_init(gg1, p1);
   F_mpz_mod_poly_init(hh1, p1);

   F_mpz_poly_to_F_mpz_mod_poly(cc, c);
   F_mpz_poly_to_F_mpz_mod_poly(gg, g);
   F_mpz_poly_to_F_mpz_mod_poly(hh, h);
   F_mpz_poly_to_F_mpz_mod_poly(aa, a);
   F_mpz_poly_to_F_mpz_mod_poly(bb, b);

//When I make a precomputing function use GG, HH instead

   F_mpz_mod_poly_rem(gg1, cc, gg);
   F_mpz_mod_poly_mulmod(gg1, gg1, bb, gg);

   F_mpz_mod_poly_rem(hh1, cc, hh);
   F_mpz_mod_poly_mulmod(hh1, hh1, aa, hh);

   F_mpz_mod_poly_to_F_mpz_poly(g1, gg1);
   F_mpz_poly_scalar_mul(g1, g1, p);

   F_mpz_mod_poly_to_F_mpz_poly(h1, hh1);
   F_mpz_poly_scalar_mul(h1, h1, p);

   F_mpz_poly_add(G, g, g1);
   F_mpz_poly_add(H, h, h1);
   
//Lifting the inverses now

   F_mpz_poly_t a1, b1, t1, t2, r, unity;
   F_mpz_poly_init(a1);
   F_mpz_poly_init(b1);
   F_mpz_poly_init(t1);
   F_mpz_poly_init(t2);
   F_mpz_poly_init(r);
   F_mpz_poly_init(unity);

   F_mpz_poly_set_coeff_si(unity, 0, -1);

   F_mpz_poly_mul(t1, a, G);
   F_mpz_poly_mul(t2, b, H);

   F_mpz_poly_add(t1, t1, t2);
   F_mpz_poly_add(t1, t1, unity);
   F_mpz_poly_neg(t1, t1);

   F_mpz_poly_scalar_div_exact(r, t1, p);
//Make a check that t1 is divisible by p
   F_mpz_mod_poly_t rr, aa1, bb1;
   F_mpz_mod_poly_init(rr, p1);
   F_mpz_mod_poly_init(aa1, p1);
   F_mpz_mod_poly_init(bb1, p1);

   F_mpz_poly_to_F_mpz_mod_poly(rr, r);

   F_mpz_mod_poly_rem(bb1, rr, gg);
   F_mpz_mod_poly_mulmod(bb1, bb1, bb, gg);

   F_mpz_mod_poly_rem(aa1, rr, hh);
   F_mpz_mod_poly_mulmod(aa1, aa1, aa, hh);   

   F_mpz_mod_poly_to_F_mpz_poly(a1, aa1);
   F_mpz_poly_scalar_mul(a1, a1, p);
   F_mpz_poly_add(A, a, a1);

   F_mpz_mod_poly_to_F_mpz_poly(b1, bb1);
   F_mpz_poly_scalar_mul(b1, b1, p);
   F_mpz_poly_add(B, b, b1);

   F_mpz_poly_set(Gout, G);
   F_mpz_poly_set(Hout, H);
   F_mpz_poly_set(Aout, A);
   F_mpz_poly_set(Bout, B);

   F_mpz_mod_poly_clear(rr);
   F_mpz_mod_poly_clear(aa1);
   F_mpz_mod_poly_clear(bb1);

   F_mpz_poly_clear(a1);
   F_mpz_poly_clear(b1);
   F_mpz_poly_clear(t1);
   F_mpz_poly_clear(t2);
   F_mpz_poly_clear(r);
   F_mpz_poly_clear(unity);

   F_mpz_mod_poly_clear(cc);
   F_mpz_mod_poly_clear(gg);
   F_mpz_mod_poly_clear(hh);
   F_mpz_mod_poly_clear(aa);
   F_mpz_mod_poly_clear(bb);
   F_mpz_mod_poly_clear(tt);
   F_mpz_mod_poly_clear(gg1);
   F_mpz_mod_poly_clear(hh1);

   F_mpz_poly_clear(c);
   F_mpz_poly_clear(g1);
   F_mpz_poly_clear(h1);
   F_mpz_poly_clear(G);
   F_mpz_poly_clear(H);
   F_mpz_poly_clear(A);
   F_mpz_poly_clear(B);
}

void F_mpz_poly_hensel_lift_without_inverse(F_mpz_poly_t Gout, F_mpz_poly_t Hout, F_mpz_poly_t f, F_mpz_poly_t g, F_mpz_poly_t h, F_mpz_poly_t a, F_mpz_poly_t b, F_mpz_t p, F_mpz_t p1, F_mpz_t big_P)
{

   F_mpz_poly_t c, g1, h1, G, H;
   F_mpz_poly_init(c);
   F_mpz_poly_init(g1);
   F_mpz_poly_init(h1);
   F_mpz_poly_init(G);
   F_mpz_poly_init(H);

   F_mpz_poly_mul(c, g, h);
   F_mpz_poly_sub(c, f, c);
   F_mpz_poly_scalar_div_exact(c, c, p);

   F_mpz_mod_poly_t cc, gg, hh, aa, bb, tt, gg1, hh1;
   F_mpz_mod_poly_init(cc, p1);
   F_mpz_mod_poly_init(gg, p1);
   F_mpz_mod_poly_init(hh, p1);
   F_mpz_mod_poly_init(aa, p1);
   F_mpz_mod_poly_init(bb, p1);
   F_mpz_mod_poly_init(tt, p1);
   F_mpz_mod_poly_init(gg1, p1);
   F_mpz_mod_poly_init(hh1, p1);

   F_mpz_poly_to_F_mpz_mod_poly(cc, c);
   F_mpz_poly_to_F_mpz_mod_poly(gg, g);
   F_mpz_poly_to_F_mpz_mod_poly(hh, h);
   F_mpz_poly_to_F_mpz_mod_poly(aa, a);
   F_mpz_poly_to_F_mpz_mod_poly(bb, b);

   F_mpz_mod_poly_rem(gg1, cc, gg);
   F_mpz_mod_poly_mulmod(gg1, gg1, bb, gg);

   F_mpz_mod_poly_rem(hh1, cc, hh);
   F_mpz_mod_poly_mulmod(hh1, hh1, aa, hh);

   F_mpz_mod_poly_to_F_mpz_poly(g1, gg1);
   F_mpz_poly_scalar_mul(g1, g1, p);

   F_mpz_mod_poly_to_F_mpz_poly(h1, hh1);
   F_mpz_poly_scalar_mul(h1, h1, p);

   F_mpz_poly_add(G, g, g1);
   F_mpz_poly_add(H, h, h1);
   
   F_mpz_poly_set(Gout, G);
   F_mpz_poly_set(Hout, H);

   F_mpz_mod_poly_clear(cc);
   F_mpz_mod_poly_clear(gg);
   F_mpz_mod_poly_clear(hh);
   F_mpz_mod_poly_clear(aa);
   F_mpz_mod_poly_clear(bb);
   F_mpz_mod_poly_clear(tt);
   F_mpz_mod_poly_clear(gg1);
   F_mpz_mod_poly_clear(hh1);

   F_mpz_poly_clear(c);
   F_mpz_poly_clear(g1);
   F_mpz_poly_clear(h1);
   F_mpz_poly_clear(G);
   F_mpz_poly_clear(H);
}

void F_mpz_poly_hensel_lift_only_inverse(F_mpz_poly_t Aout, F_mpz_poly_t Bout, F_mpz_poly_t f, F_mpz_poly_t G, F_mpz_poly_t H, F_mpz_poly_t a, F_mpz_poly_t b, F_mpz_t p, F_mpz_t p1, F_mpz_t big_P)
{

   F_mpz_poly_t A, B;
   F_mpz_poly_init(A);
   F_mpz_poly_init(B);

   F_mpz_poly_t a1, b1, t1, t2, r, unity;
   F_mpz_poly_init(a1);
   F_mpz_poly_init(b1);
   F_mpz_poly_init(t1);
   F_mpz_poly_init(t2);
   F_mpz_poly_init(r);
   F_mpz_poly_init(unity);

   //Make a check that c is divisible by p
   F_mpz_mod_poly_t gg, hh, aa, bb, g, h;
   F_mpz_mod_poly_init(gg, p1);
   F_mpz_mod_poly_init(hh, p1);
   F_mpz_mod_poly_init(g, p);
   F_mpz_mod_poly_init(h, p);
   F_mpz_mod_poly_init(aa, p1);
   F_mpz_mod_poly_init(bb, p1);

   F_mpz_mod_poly_t rr, aa1, bb1;
   F_mpz_mod_poly_init(rr, p1);
   F_mpz_mod_poly_init(aa1, p1);
   F_mpz_mod_poly_init(bb1, p1);

   F_mpz_poly_to_F_mpz_mod_poly(g, G);//reduce mod p
   F_mpz_poly_to_F_mpz_mod_poly(h, H);//reduce mod p
   F_mpz_mod_poly_set(gg, g);//embed mod p1
   F_mpz_mod_poly_set(hh, h);//embed mod p1
   F_mpz_poly_to_F_mpz_mod_poly(aa, a);
   F_mpz_poly_to_F_mpz_mod_poly(bb, b);
//Lifting the inverses now
   F_mpz_poly_set_coeff_si(unity, 0, -1);
   F_mpz_poly_mul(t1, a, G);
   F_mpz_poly_mul(t2, b, H);

   F_mpz_poly_add(t1, t1, t2);
   F_mpz_poly_add(t1, t1, unity);
   F_mpz_poly_neg(t1, t1);

   F_mpz_poly_scalar_div_exact(r, t1, p);
   F_mpz_poly_to_F_mpz_mod_poly(rr, r);

   F_mpz_mod_poly_rem(bb1, rr, gg);
   F_mpz_mod_poly_mulmod(bb1, bb1, bb, gg);

   F_mpz_mod_poly_rem(aa1, rr, hh);
   F_mpz_mod_poly_mulmod(aa1, aa1, aa, hh);   

   F_mpz_mod_poly_to_F_mpz_poly(a1, aa1);
   F_mpz_poly_scalar_mul(a1, a1, p);
   F_mpz_poly_add(A, a, a1);

   F_mpz_mod_poly_to_F_mpz_poly(b1, bb1);
   F_mpz_poly_scalar_mul(b1, b1, p);
   F_mpz_poly_add(B, b, b1);

   F_mpz_poly_set(Aout, A);
   F_mpz_poly_set(Bout, B);

   F_mpz_mod_poly_clear(rr);
   F_mpz_mod_poly_clear(aa1);
   F_mpz_mod_poly_clear(bb1);

   F_mpz_poly_clear(a1);
   F_mpz_poly_clear(b1);
   F_mpz_poly_clear(t1);
   F_mpz_poly_clear(t2);
   F_mpz_poly_clear(r);
   F_mpz_poly_clear(unity);

   F_mpz_mod_poly_clear(gg);
   F_mpz_mod_poly_clear(hh);
   F_mpz_mod_poly_clear(aa);
   F_mpz_mod_poly_clear(bb);
   F_mpz_mod_poly_clear(g);
   F_mpz_mod_poly_clear(h);

   F_mpz_poly_clear(A);
   F_mpz_poly_clear(B);
}

void F_mpz_poly_rec_tree_hensel_lift(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, F_mpz_t p, F_mpz_poly_t f, long j, long inv, F_mpz_t p1, F_mpz_t big_P)
{

   if (j < 0) return;

   if (inv == 1)
      F_mpz_poly_hensel_lift(v[j], v[j+1], w[j], w[j+1], f, v[j], v[j+1], w[j], w[j+1], p, p1, big_P);
   else if (inv == -1)
      F_mpz_poly_hensel_lift_only_inverse(w[j], w[j+1], f, v[j], v[j+1], w[j], w[j+1], p, p1, big_P);
   else
      F_mpz_poly_hensel_lift_without_inverse(v[j], v[j+1], f, v[j], v[j+1], w[j], w[j+1], p, p1, big_P);
   
//altered to check a bug, should be Hensel_Lift1
   F_mpz_poly_rec_tree_hensel_lift(link, v, w, p, v[j],   link[j],   inv, p1, big_P);
   F_mpz_poly_rec_tree_hensel_lift(link, v, w, p, v[j+1], link[j+1], inv, p1, big_P);
}

void F_mpz_poly_tree_hensel_lift(long *link, F_mpz_poly_t *v, F_mpz_poly_t *w, long e0, long e1, F_mpz_poly_t monic_f, long inv, long p, long r, F_mpz_t P)
{

   F_mpz_t temp, p0, p1;
   F_mpz_init(p0);
   F_mpz_init(p1);
   F_mpz_init(temp);

   F_mpz_set_ui(temp, p);

   F_mpz_pow_ui(p0, temp, e0);
   F_mpz_pow_ui(p1, temp, e1 - e0);

   F_mpz_mul2(P, p0, p1);

   F_mpz_poly_rec_tree_hensel_lift(link, v, w, p0, monic_f, 2*r-4, inv, p1, P);

   F_mpz_clear(temp);
   F_mpz_clear(p0);
   F_mpz_clear(p1);

}


/*
   This behemoth is going to take the local factors (in local_fac) and Hensel lift them until they are known mod p^target_exp.  These lifted factors will be stored (in the same ordering) in lifted_fac.  It is assumed that link, v, and w are initialized arrays of F_mpz_poly_t's with at least 2*r - 2 entries and that r >= 2. These are done outside of this function so that you can keep them for restarting Hensel lifting later.
*/

ulong _F_mpz_poly_start_hensel_lift(F_mpz_poly_factor_t lifted_fac, long * link, F_mpz_poly_t * v, F_mpz_poly_t * w, F_mpz_poly_t f, zmod_poly_factor_t local_fac, ulong target_exp)
{
   ulong r = local_fac->num_factors;
   ulong p = (local_fac->factors[0])->p;
   ulong i;
   F_mpz_t P, big_P, temp;
   F_mpz_init(temp);
   F_mpz_init(P);
   F_mpz_init(big_P);

   F_mpz_set_ui(P, p);

   F_mpz_poly_build_hensel_tree(link, v, w, local_fac);

   ulong current_exp = 1;

   F_mpz_poly_t monic_f;
   F_mpz_poly_init(monic_f);
   F_mpz_pow_ui(big_P, P, target_exp);

   if (F_mpz_is_one(f->coeffs + f->length - 1))
   {
      F_mpz_poly_set(monic_f, f);
   }
   else if (F_mpz_is_m1(f->coeffs + f->length - 1))
   {
      F_mpz_poly_neg(monic_f, f);
   }
   else
   {
      F_mpz_mod(temp, f->coeffs + f->length - 1, big_P);
      int OK = F_mpz_invert(temp, temp, big_P);
      if (OK == 0)
      {
         printf(" some problem\n");
         abort();
      }
      F_mpz_poly_scalar_mul(monic_f, f, temp);
      F_mpz_poly_smod(monic_f, monic_f, big_P);
   }
//later we're going to fine tune this powering process, in the meantime I want to have an array of exponents which we walk through
   ulong num_steps = 2 + (ulong) floor( log2( (double) (target_exp) ) );
   ulong exponents[num_steps];
   ulong pow = 1;
   for (i = 0; (i < num_steps) && (pow < target_exp); i++)
   { 
      exponents[i] = pow;
      pow = pow * 2;
   }
   num_steps = i;
   if (exponents[num_steps - 1] != target_exp ){
      exponents[num_steps] = target_exp;
      num_steps++;
   }
//here num_steps is actually the number of meaningful numbers in the array exponent so num_steps-2 means that the final time in the loop has 1 and one last
//time outside of the loop with 0

   for (i = 0; i < num_steps - 2; i++){
      F_mpz_poly_tree_hensel_lift(link, v, w, exponents[i], exponents[i+1], monic_f, 1, p, r, P);
   }
//Last run doesn't calculate the inverses
   F_mpz_poly_tree_hensel_lift(link, v, w, exponents[i], exponents[i+1], monic_f, 0, p, r, P);

   ulong prev_exp = exponents[i];

   F_mpz_poly_clear(monic_f);
//Now everything is lifted to p^target_exp just need to insert the factors into their correct places in lifted_fac.
   if (r > lifted_fac->alloc)
   {
      lifted_fac->factors = (F_mpz_poly_t *) flint_heap_realloc_bytes(lifted_fac->factors, sizeof(F_mpz_poly_t)*r);
      lifted_fac->exponents = (unsigned long *) flint_heap_realloc(lifted_fac->exponents, r);
      for (unsigned long i = lifted_fac->alloc; i < r; i++)
         F_mpz_poly_init(lifted_fac->factors[i]);
      lifted_fac->alloc = r;
   }

   for(i = 0; i < 2*r -2; i++){ 
      if (link[i] < 0){
//        Want to get F_mpz_mod_poly and want it to be smodded by P 
         F_mpz_poly_smod(lifted_fac->factors[-link[i]-1], v[i], P);
         lifted_fac->exponents[-link[i]-1] = 1L; 
      }
   }
   lifted_fac->num_factors = r;

   F_mpz_clear(temp);
   F_mpz_clear(big_P);
   F_mpz_clear(P);

   return prev_exp;
}

ulong _F_mpz_poly_continue_hensel_lift(F_mpz_poly_factor_t lifted_fac, long * link, F_mpz_poly_t * v, F_mpz_poly_t * w, F_mpz_poly_t f, ulong prev_exp, ulong current_exp, ulong target_exp, ulong p, ulong r)
{
   ulong i;
   F_mpz_t P, big_P, temp;
   F_mpz_init(temp);
   F_mpz_init(P);
   F_mpz_init(big_P);

   F_mpz_set_ui(P, p);

   F_mpz_poly_t monic_f;
   F_mpz_poly_init(monic_f);
   F_mpz_pow_ui(big_P, P, target_exp);

   if (F_mpz_is_one(f->coeffs + f->length - 1))
   {
      F_mpz_poly_set(monic_f, f);
   }
   else if (F_mpz_is_m1(f->coeffs + f->length - 1))
   {
      F_mpz_poly_neg(monic_f, f);
   }
   else
   {
      F_mpz_mod(temp, f->coeffs + f->length - 1, big_P);
      int OK = F_mpz_invert(temp, temp, big_P);
      if (OK == 0)
      {
         printf(" some problem\n");
         abort();
      }
      F_mpz_poly_scalar_mul(monic_f, f, temp);
//Note that this and one other smod could really be mod or scalar_mod once one of those exists
      F_mpz_poly_smod(monic_f, monic_f, big_P);
   }
//later we're going to fine tune this powering process, in the meantime I want to have an array of exponents which we walk through
//NEEDS SOME CHECKING!!!
   ulong num_steps = 2 + (ulong) floor( log2( (double) (target_exp) - (double) prev_exp ) );
   ulong exponents[num_steps];
   ulong pow = current_exp;
   exponents[0] = prev_exp;
   for (i = 1; (i < num_steps) && (pow < target_exp); i++)
   { 
      exponents[i] = pow;
      pow = pow * 2;
   }
   num_steps = i;
   if (exponents[num_steps - 1] != target_exp ){
      exponents[num_steps] = target_exp;
      num_steps++;
   }
//here num_steps is actually the number of meaningful numbers in the array exponent so num_steps-2 means that the final time in the loop has 1 and one last
//time outside of the loop with 0

   F_mpz_poly_tree_hensel_lift(link, v, w, exponents[0], exponents[1], monic_f, -1, p, r, P);
   for (i = 1; i < num_steps - 2; i++){
      F_mpz_poly_tree_hensel_lift(link, v, w, exponents[i], exponents[i+1], monic_f, 1, p, r, P);
   }
//Last run doesn't calculate the inverses
   F_mpz_poly_tree_hensel_lift(link, v, w, exponents[i], exponents[i+1], monic_f, 0, p, r, P);

   ulong new_prev_exp = exponents[i];

   F_mpz_poly_clear(monic_f);
//Now everything is lifted to p^target_exp just need to insert the factors into their correct places in lifted_fac.
   if (r > lifted_fac->alloc)
   {
      lifted_fac->factors = (F_mpz_poly_t *) flint_heap_realloc_bytes(lifted_fac->factors, sizeof(F_mpz_poly_t)*r);
      lifted_fac->exponents = (unsigned long *) flint_heap_realloc(lifted_fac->exponents, r);
      for (unsigned long i = lifted_fac->alloc; i < r; i++)
         F_mpz_poly_init(lifted_fac->factors[i]);
      lifted_fac->alloc = r;
   }

   for(i = 0; i < 2*r -2; i++){ 
      if (link[i] < 0){
//        Want to get F_mpz_mod_poly and want it to be smodded by P 
         F_mpz_poly_smod(lifted_fac->factors[-link[i]-1], v[i], P);
         lifted_fac->exponents[-link[i]-1] = 1L; 
      }
   }
   lifted_fac->num_factors = r;

   F_mpz_clear(temp);
   F_mpz_clear(big_P);
   F_mpz_clear(P);
   return new_prev_exp;
}

void F_mpz_poly_hensel_lift_once( F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, zmod_poly_factor_t local_fac, ulong target_exp)
{
   ulong r = local_fac->num_factors;

   F_mpz_poly_t v[2*r-2];
   F_mpz_poly_t w[2*r-2];   
   long link[2*r-2];

   for(long i = 0; i < 2*r - 2; i++)
   {
      F_mpz_poly_init(v[i]);
      F_mpz_poly_init(w[i]);
   }

   ulong trash = _F_mpz_poly_start_hensel_lift(lifted_fac, link, v, w, F, local_fac, target_exp);

   for (long i = 0; i < 2*r-2; i++){
      F_mpz_poly_clear(v[i]);
      F_mpz_poly_clear(w[i]);
   }
}

/***************************************************

Naive Zassenhaus

***************************/

void F_mpz_poly_zassenhaus_naive(F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc){

   ulong r = lifted_fac->num_factors;
   F_mpz_poly_t f;
   F_mpz_poly_init(f);

   F_mpz_poly_set(f, F);

   F_mpz_poly_t Q,R;
   F_mpz_poly_init(Q);
   F_mpz_poly_init(R);
   F_mpz_poly_t tryme;
   F_mpz_poly_init(tryme);
   F_mpz_t temp_lc;

   int k;
   int l;
   int indx;
   int used_arr[r];
   for(l = 0; l < r; l++){
      used_arr[l] = 0;
   }

   for (k = 1; k < r; k++){
      ulong count = 0;
      ulong sub_arr[k];
      for(l = 0; l < k; l++){
         sub_arr[l] = l;
      }
      indx = k-1;
      sub_arr[indx]--;
      while ((indx >= 0)){
         sub_arr[indx] = sub_arr[indx] + 1;
         for (l = indx + 1; l < k; l++){
            sub_arr[l] = sub_arr[l-1] + 1UL;
         }
         if (sub_arr[k-1] > r-1UL ){
            indx--;
         }
         else{
            for(l = 0; l < k; l++){
               if (used_arr[sub_arr[l]] == 1)
                  break;
            }
//Need to involve lc, perhaps set coeff 0 to lc and do lc * rest and check if under M_bits... here I'm using a trial division... hmm
            F_mpz_poly_fit_length(tryme, 1UL);
            tryme->length = 1UL;
            F_mpz_set(tryme->coeffs + 0, lc);
            for(l = 0; l < k; l++){
               F_mpz_poly_mul(tryme, tryme, lifted_fac->factors[sub_arr[l]]);
            }
            F_mpz_poly_smod(tryme, tryme, P);
            F_mpz_init(temp_lc);
            F_mpz_poly_content(temp_lc, tryme);
            F_mpz_poly_scalar_div_exact(tryme, tryme, temp_lc);
            F_mpz_poly_divrem(Q, R, f, tryme);
            if (R->length == 0){
//        FOUND ONE!!!!!
               F_mpz_poly_factor_insert(final_fac, tryme, exp);
               for(l = 0; l < k; l++){
                  used_arr[sub_arr[l]] = 1;
                  count++;
               }
               F_mpz_poly_set(f, Q);
               F_mpz_set(lc, Q->coeffs + Q->length - 1 );
//If r-count = k then the rest are irreducible.  But haven't added that
            }
            F_mpz_clear(temp_lc);
            indx = k-1;
         }
      }
//This is where we switch to the next loop for k.
//So we will have found all factors using <= k local factors
//We should/could update f to be the rest divided away (or multiply the remaining)
// could also adjust r.  It is the number of remaining factors  
//so if you update then test if r = k or k+1 in which case the remaining f is irreducible.
   }
   ulong test = 0;
   for (l = 0; l < r; l++){
      test = test + used_arr[l];
   }
   if (test == 0)
      F_mpz_poly_factor_insert(final_fac, f, exp);
   F_mpz_poly_clear(f);
   F_mpz_poly_clear(tryme);
   F_mpz_poly_clear(Q);
   F_mpz_poly_clear(R);
   return;
}

/*********

   Factoring wrapper after square free part

*********/

void F_mpz_poly_factor_sq_fr_prim( F_mpz_poly_factor_t final_fac, ulong exp, F_mpz_poly_t f ){

   if (f->length <= 1){
      return;
   }
   if (f->length == 2){
      F_mpz_poly_factor_insert( final_fac, f, exp);
      return;
   }
   ulong len = f->length;
   F_mpz_t lc;
   F_mpz_init(lc);   
   F_mpz_set(lc, f->coeffs + len - 1);
   ulong M_bits = 0;
   M_bits = M_bits + F_mpz_bits(lc) + FLINT_ABS(F_mpz_poly_max_bits(f)) + len + (long)ceil(log2((double) len));
   zmod_poly_t F, F_d, F_sbo;
   int tryme = 1;
   ulong p = 2UL;
   long i, num_primes;
   zmod_poly_factor_t fac;

   ulong min_p = p; 
   unsigned long lead_coeff;
   long r, min_r;
   min_r = len;
   i = 0;

   local_factor_start = clock();

for(num_primes = 1; num_primes < 3; num_primes++)
{
   tryme  = 1;
   for (; (i < 200) && (tryme == 1); i++){

      zmod_poly_init(F, p);
      zmod_poly_init(F_d, p);
      zmod_poly_init(F_sbo, p);

      F_mpz_poly_to_zmod_poly(F, f);
      if (F->length < f->length){
         p = z_nextprime( p, 0);
         zmod_poly_clear(F_d);
         zmod_poly_clear(F_sbo);
         zmod_poly_clear(F);
         continue;
      }

//Maybe faster some other way... checking if squarefee mod p
      zmod_poly_derivative(F_d, F);
      zmod_poly_gcd(F_sbo, F, F_d);


      if (zmod_poly_is_one(F_sbo)){
         tryme = 0;
      }
      else{
         p = z_nextprime( p, 0);
         zmod_poly_clear(F);
      }
      zmod_poly_clear(F_d);
      zmod_poly_clear(F_sbo);

   }
   if (i == 200){
      printf("wasn't square_free after 200 primes, maybe an error\n");
      zmod_poly_clear(F);
      F_mpz_clear(lc);
      return;
   }

   zmod_poly_factor_init(fac);
   lead_coeff = zmod_poly_factor(fac, F);

   r = fac->num_factors;
   zmod_poly_factor_clear(fac);
   zmod_poly_clear(F);

   printf("prime try r = %ld, p = %ld\n", r, p);

   if (r <= min_r)
   {
      min_r = r;
      min_p = p;
   }
   p = z_nextprime( p, 0);
}

   p = min_p;

   zmod_poly_init(F, p);
   F_mpz_poly_to_zmod_poly(F, f);
   zmod_poly_factor_init(fac);
   lead_coeff = zmod_poly_factor(fac, F);
   r = fac->num_factors;

   local_factor_stop = clock();
   local_factor_total = local_factor_total + local_factor_stop - local_factor_start;

   printf("primes chosen r = %ld, p = %ld total time = %f\n", r, p, (double) local_factor_total/(double) CLOCKS_PER_SEC);



   int use_Hoeij_Novocin = 0;
   int solved_yet = 0;
//doing away with mexpo for a while, might bring it back... don't forget intialization too
   int mexpo[4];
//   int mexpo[r + 2 * (f->length - 1)];
   F_mpz_mat_t M;
   long U_exp = r/4;

//In the near future we should go back and try some more primes might deduce irreducibility or find smaller r
   if (r > 6){
      use_Hoeij_Novocin = 1;
      F_mpz_mat_init_identity(M, r);

//Experimentally we are adjusting U = I*2^r this needs to match the precision in vHN approach
      F_mpz_mat_scalar_mul_2exp(M, M, U_exp);
      ulong i;
      for (i = 0; i < 4; i++)
         mexpo[i] = 0;
   }
   if (r == 0){
      printf("FLINT-exception: something broke\n");
      zmod_poly_clear(F);
      zmod_poly_factor_clear(fac);
      F_mpz_clear(lc);
      abort();
   }
   if (r == 1){
      F_mpz_poly_factor_insert(final_fac, f, exp);
      zmod_poly_clear(F);
      zmod_poly_factor_clear(fac);
      F_mpz_clear(lc);
      return;
   }
//Begin Hensel Lifting phase, make the tree in v, w, and link
// Later we'll do a check if use_Hoeij_Novocin (try for smaller a)
   ulong a;
   a = (long) ceil( (double) M_bits / log2( (double)p ) );
   a = (long) pow( (double) 2, ceil( log2( (double) a ) ) );

   if (use_Hoeij_Novocin == 1)
   {

   //Here we can add the cool new part...
      F_mpz_t lead_b, trail_b, mid_b;
      F_mpz_init(lead_b);
      F_mpz_init(trail_b);
      F_mpz_init(mid_b);
cld_data_start = clock();
      F_mpz_poly_CLD_bound(lead_b, f, len - 3);
//      F_mpz_poly_CLD_bound(mid_b, f, (len - 2)/2);
      F_mpz_poly_CLD_bound(trail_b, f, 0);
cld_data_stop = clock();
cld_data_total = cld_data_total + cld_data_stop - cld_data_start;

printf(" first two clds took %f seconds\n", (double) cld_data_total/ (double) CLOCKS_PER_SEC );

//reusing the lead_b to be the new average and trail_b to be 3 since there isn't an F_mpz_div_ui
/*      F_mpz_add(lead_b, lead_b, mid_b);
      F_mpz_add(lead_b, lead_b, trail_b);

      F_mpz_set_ui(trail_b, 3);
      F_mpz_cdiv_q(lead_b, lead_b, trail_b);
*/
      int cmp_b = F_mpz_cmpabs(lead_b, trail_b);

      ulong avg_b;
      if (cmp_b <= 0)
         avg_b = F_mpz_bits(lead_b);
      else
         avg_b = F_mpz_bits(trail_b);

//Trying to get (a-b)log(p)*(len-2) = .12*r^2 where b = avg_b/log2(p)...

//      long n_a = (long) (double)( 0.12 * r * r /(double)( (len - 2) ) + avg_b / log2( (double) p )  );

      long n_a = (long) (double)( (1.5)*r  + (double) avg_b / log2( (double) p )  );

//      n_a = (long) pow( (double) 2, ceil( log2( (double) n_a ) ) );

      a = FLINT_MIN(a, n_a);

      a = 50;

      printf(" new a = %ld\n", a);

      F_mpz_clear(lead_b);
      F_mpz_clear(trail_b);
      F_mpz_clear(mid_b);
   }

   F_mpz_poly_t v[2*r-2];
   F_mpz_poly_t w[2*r-2];   
   long link[2*r-2];
//P will be the F_mpz modulus
   F_mpz_poly_factor_t lifted_fac;
   F_mpz_poly_factor_init(lifted_fac);

   F_mpz_t P;
   F_mpz_init(P);


   for (i = 0; i < 2*r-2; i++)
      F_mpz_poly_init(v[i]);

   for (i = 0; i < 2*r-2; i++)
      F_mpz_poly_init(w[i]);
   hensel_start = clock();
   ulong prev_exp = _F_mpz_poly_start_hensel_lift(lifted_fac, link, v, w, f, fac, a);
   hensel_stop = clock();
   hensel_total = hensel_total + hensel_stop - hensel_start;

   printf("hensel lifted for %f seconds so far\n", (double) hensel_total / (double) CLOCKS_PER_SEC);
//clearing fac early, don't need them anymore
   zmod_poly_factor_clear(fac);

   while( solved_yet == 0 ){
//Have now Hensel lifted to p^a for the precalculated a, in the optimized version we will lift even less
//Here let's make a list of Hensel lifted factors for grabbing information and trial testing.

//Now we should undo the mystical link part to find the original local factors in original order, see the long explanation in the Hensel code
//Now we are ready to to the Zassenhaus testing... later the r > 20 (or even 10) test could go here
      F_mpz_set_ui(P, p);
      F_mpz_pow_ui(P, P, a);

      if (use_Hoeij_Novocin == 1){

         solved_yet = F_mpz_poly_factor_sq_fr_vHN(final_fac, lifted_fac, f, P, exp, M, mexpo, U_exp);
         if (solved_yet == 0){
//This is where we increase the Hensel Accuracy and go back
            hensel_start = clock();
            prev_exp = _F_mpz_poly_continue_hensel_lift(lifted_fac, link, v, w, f, prev_exp, a, 2*a, p, r);
            hensel_stop = clock();
            hensel_total = hensel_total + hensel_stop - hensel_start;

            printf("hensel lifted for %f seconds so far\n", (double) hensel_total / (double) CLOCKS_PER_SEC);

            a = 2*a;
         }
         else if (solved_yet == 5){
//This is where we increase the Hensel Accuracy and go back
            hensel_start = clock();
            prev_exp = _F_mpz_poly_continue_hensel_lift(lifted_fac, link, v, w, f, prev_exp, a, 4*a, p, r);
            hensel_stop = clock();
            hensel_total = hensel_total + hensel_stop - hensel_start;

            printf("hensel lifted for %f seconds so far\n", (double) hensel_total / (double) CLOCKS_PER_SEC);

            a = 4*a;
            solved_yet = 0;
         }
      }
      else{
         printf("called zassenhaus\n");
         F_mpz_poly_zassenhaus_naive(final_fac, lifted_fac, f, P, exp, lc);
         solved_yet = 1;
      }
   }
//Done factoring, just gotta clean house
   for (i = 0; i < 2*r-2; i++){
      F_mpz_poly_clear(v[i]);
      F_mpz_poly_clear(w[i]);
   }
   F_mpz_poly_factor_clear(lifted_fac);
   zmod_poly_clear(F);
   F_mpz_clear(P);
   if (use_Hoeij_Novocin == 1){
      F_mpz_mat_clear(M);
   }
   return;
}

void F_mpz_poly_factor(F_mpz_poly_factor_t final_fac, F_mpz_t cong, F_mpz_poly_t G){
   if (G->length == 0){
      F_mpz_set_ui(cong, 0UL);
      return;
   }
   if (G->length == 1){
      F_mpz_set(cong, G->coeffs);
      return;
   }
   F_mpz_poly_t g;
   F_mpz_poly_init(g);

   if (G->length == 2){
      F_mpz_poly_content(cong, G);
      F_mpz_poly_scalar_div_exact(g, G, cong);
      F_mpz_poly_factor_insert( final_fac, g, 1UL);
      return;
   }
//Does a presearch for a factor of form x^(x_pow)
   ulong x_pow = 0;
   while (F_mpz_is_zero(G->coeffs + x_pow)){
      x_pow++; 
   }
   if (x_pow != 0){
      F_mpz_poly_t temp_x;
      F_mpz_poly_init(temp_x);
      F_mpz_poly_set_coeff_ui(temp_x, 1, 1);
      F_mpz_poly_factor_insert(final_fac, temp_x, x_pow);
      F_mpz_poly_clear(temp_x);
   }
   F_mpz_poly_right_shift(g, G, x_pow);
   F_mpz_poly_factor_t sq_fr_fac;
//Could make other tests for x-1 or simple things 
// maybe take advantage of the composition algorithm
   F_mpz_poly_factor_init( sq_fr_fac );
   F_mpz_poly_squarefree(sq_fr_fac, cong, g);
//Now we can go through and factor each square free guy and add it to final factors.
   for(ulong j = 0; j < sq_fr_fac->num_factors; j++){
      F_mpz_poly_factor_sq_fr_prim(final_fac, sq_fr_fac->exponents[j], sq_fr_fac->factors[j]);
   }
   F_mpz_poly_factor_clear( sq_fr_fac );
   F_mpz_poly_clear(g);
}

/*============================================================================

   Hoeij/Novocin approach

============================================================================*/

void _F_mpz_poly_factor_CLD_mat(F_mpz_mat_t res, F_mpz_poly_t F, F_mpz_poly_factor_t lifted_fac, F_mpz_t P, ulong N)
{
//Here we want to take F such that the product of the monic lifted_facs is equiv to F mod P and find top N and 
//bottom N coeffs of each local Log Derivative.  Store the results in res with an extra row at the bottom for the bounds for that column
   ulong n;
   ulong d = lifted_fac->num_factors;
   F_mpz_poly_t gd,gcld;

   if (2*N >= F->length - 1){
      F_mpz_mat_clear(res);
      F_mpz_mat_init(res, d+1, F->length - 1);
      long i;
      for (i = 0; i < F->length - 1; i++)
         F_mpz_poly_CLD_bound(res->rows[d] + i, F, i);

      for (i = 0; i < d; i++){
         F_mpz_poly_init(gd);
         F_mpz_poly_init(gcld);
         F_mpz_poly_derivative(gd, lifted_fac->factors[i]);
         F_mpz_poly_mul(gcld, F, gd);
         F_mpz_poly_div(gcld, gcld, lifted_fac->factors[i]);
         F_mpz_poly_smod(gcld, gcld, P);
         long j;
         for (j = 0; j < F->length - 1; j++)
            F_mpz_set(res->rows[i] + j, gcld->coeffs + j);
         F_mpz_poly_clear(gd);
         F_mpz_poly_clear(gcld);
      }
      return;
   }
   n = N;
   F_mpz_mat_clear(res);
   F_mpz_mat_init(res, d+1, 2*n);

   long i;
   for (i = 0; i < n; i++){
      F_mpz_poly_CLD_bound(res->rows[d] + i, F, i);
      F_mpz_poly_CLD_bound(res->rows[d] + 2*n - 1 - i, F, F->length - 2 -i);
   }
   F_mpz_t temp[n];
   for (i = 0; i < n; i++)
      F_mpz_init(temp[i]);

   for (i = 0; i < d; i++){
      F_mpz_poly_init(gd);
      F_mpz_poly_init(gcld);

      F_mpz_poly_derivative(gd, lifted_fac->factors[i]);
//Should do upper and lower trunc multiplication soon, for speed sake
      F_mpz_poly_mul(gcld, F, gd);
      F_mpz_poly_div_trunc_modp(temp, gcld, lifted_fac->factors[i], P, n);

      long j;
      for (j = 0; j < n; j++)
         F_mpz_set(res->rows[i] + j, temp[j]);
      F_mpz_poly_div_upper_trunc_modp(temp, gcld, lifted_fac->factors[i], P, n);
      for (j = 0; j < n; j++)
         F_mpz_set(res->rows[i] + 2*n -1 - j, temp[j]);
      F_mpz_poly_clear(gd);
      F_mpz_poly_clear(gcld);
   }

   for (i = 0; i < n; i++)
      F_mpz_clear(temp[i]);
}

int _F_mpz_poly_try_to_solve(int num_facs, ulong * part, F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc)
{
   if (num_facs == 1){
      F_mpz_poly_factor_insert(final_fac, F, exp);
//      printf("Proved irreducibility\n");
      return 1;
   }
//we know that there is a 0-1 potential basis, let's make the potential factors and sort them by degree
   ulong r = lifted_fac->num_factors;
   if (r <= num_facs)
   {
      printf("highly likely not solved yet... \n");
      return 0;
   }


   F_mpz_poly_factor_t trial_factors;
   F_mpz_poly_factor_init(trial_factors);

   F_mpz_poly_t tryme;
   F_mpz_poly_init(tryme);
   F_mpz_t temp_lc;
   F_mpz_init(temp_lc);



   printf(" not yet r = %ld\n", r);
/*   ulong biggest, second;
   biggest = 0;
   second = 0;*/
   int i;
   for (i = 1; i <= num_facs; i++){
      F_mpz_poly_fit_length(tryme, 1UL);
      tryme->length = 1UL;
      F_mpz_set(tryme->coeffs + 0, lc);
      int j;
      for (j = 0; j < r; j++){
         if (part[j] == i)
         {
            F_mpz_poly_mul(tryme, tryme, lifted_fac->factors[j]);
            F_mpz_poly_smod(tryme, tryme, P);
         }
      }
      F_mpz_init(temp_lc);
      F_mpz_poly_content(temp_lc, tryme);
      F_mpz_poly_scalar_div_exact(tryme, tryme, temp_lc);
      F_mpz_poly_factor_insert(trial_factors, tryme, 1UL);
/*      if (tryme->length > biggest){
         second = biggest;
         biggest = tryme->length
      }
      else if (tryme-> length > second)
         second = tryme->length;*/
   }
//OK this was not optimal but we now have num_facs potential factors stored in trial_factors
//Would be more optimal if I did not multiply by lc yet, or if I trial divided along the way...
//BUT with the early termination we want to attempt our trial divisions from the bottom up so that
//if we find all but the largest factor and it matches our bound for the number of factors then 
//we've proven the factorization is legit
//maybe figure out what lengths I can prove this way and if I can prove biggest hooray go for it, if I can't then
//maybe I can prove second, in which case sort, if neither then try, maybe Hensel lift again
printf(" maybe here\n");
   for (i = 0; i < num_facs - 1; i++){
      int j;
      for (j = i+1; j <  num_facs; j++){
         if( (trial_factors->factors[i])->length < (trial_factors->factors[j])->length )
            F_mpz_poly_swap(trial_factors->factors[i], trial_factors->factors[j]);
      }
   }
printf("nope not there trial_factors->num_factors = %ld\n", trial_factors->num_factors);
   F_mpz_poly_t f,Q,R;
   F_mpz_poly_init(f);
   F_mpz_poly_init(Q);
   F_mpz_poly_init(R);
   F_mpz_poly_set(f, F);
   int j;
   for (i = 0; i < trial_factors->num_factors; i++){
      if (num_facs == 1){
         for (j = 0; j < i; j++)
            F_mpz_poly_factor_insert(final_fac, trial_factors->factors[j], exp);
         F_mpz_poly_factor_insert(final_fac, f, exp);
         return 1;
      }


      F_mpz_poly_divrem(Q, R, f, trial_factors->factors[i]);

      if (R->length == 0){
         //found one!!!! Don't insert just yet in case we find some but not all (which we handle suboptimally at the moment)
//         F_mpz_poly_factor_insert(final_fac, trial_factors->factors[i], exp);
         F_mpz_poly_set(f, Q);
         num_facs--;
      }
      else{
//Did not solve

/**** 
      Possibility of handling some things faster here by switching to Zassenhaus when we haven't solved but number of factors is small
      also we are not keeping the potentially valuable info in part (although indirectly we are...)
      Maybe it would be worth it to hand back the 0-1 version of U and start from scratch with new data (or recalculate old data)... 
      Sorry for all the comments but I took away the part/zassenhaus section for the moment.  My reasoning is that we might not have tons of Hensel 
      lifting in the final version in which case Zassenhaus could fail... which would suck to program that exception too or even figure out how to 
      use what it did find.  So for now the 'easy' route is just going back to van Hoeij regardless of the number of proven factors. 
 ****/

//         if (num_facs > 10){
// Still need van Hoeij and Hensel and stuff
printf("or thereere\n");
   F_mpz_poly_clear(f);
   F_mpz_poly_clear(Q);
   F_mpz_poly_clear(R);

   F_mpz_clear(temp_lc);
   F_mpz_poly_clear(tryme);
   F_mpz_poly_factor_clear(trial_factors);
         return 0;
//         }
//         else{
//Attempt part/zassenhaus.  We haven't found the right answers yet, but part could be used to help... maybe introduce a used array for keeping track of what's already been found
//Nah simplest would be to just pass the trial factors rather than lifted factors... 
//Isn't it concievable that we could reach this point with not enough Hensel Lifting?  In that case we'll need a check of landau bound...
//            printf("Will have solved once I write the partZassenhaus\n");
//            return 1;
//         }
      }
   }

   F_mpz_poly_clear(f);
   F_mpz_poly_clear(Q);
   F_mpz_poly_clear(R);

   F_mpz_clear(temp_lc);
   F_mpz_poly_clear(tryme);
   F_mpz_poly_factor_clear(trial_factors);
   return 0;
}

int _F_mpz_mat_check_if_solved(F_mpz_mat_t M, ulong r, F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_t lc){
   F_mpz_mat_t U;
   F_mpz_mat_init(U,0,0);
   F_mpz_mat_get_U(U, M, r);
   ulong part[U->c];
   ulong j;
   for (j = 0; j < U->c; j++)
      part[j] = 0;
   int ok = F_mpz_mat_check_0_1(part, U);

   if (ok == 0){
//not a 0-1 basis
      F_mpz_mat_clear(U);
      return 0;
   }

   if (ok > U->c){
      F_mpz_mat_clear(U);
      return 0;
   }

   check_if_solve_start = clock();
   int trym = _F_mpz_poly_try_to_solve(ok, part, final_fac, lifted_fac, F, P, exp, lc);
   check_if_solve_stop = clock();
   check_if_solve_total = check_if_solve_total + check_if_solve_stop - check_if_solve_start;
   printf("So far total checking time = %f\n", (double) check_if_solve_total /(double) CLOCKS_PER_SEC );
   F_mpz_mat_clear(U);
   return trym;
}

int F_mpz_poly_factor_sq_fr_vHN(F_mpz_poly_factor_t final_fac, F_mpz_poly_factor_t lifted_fac, F_mpz_poly_t F, F_mpz_t P, ulong exp, F_mpz_mat_t M, int * cexpo, long U_exp)
{
   int return_me = 0;
   ulong N = F->length - 1;
   F_mpz_t lc;
   F_mpz_init(lc);
   F_mpz_set(lc, F->coeffs + N);

   ulong r = lifted_fac->num_factors;
   ulong s = r; //M->r;
   F_mpz_mat_t col;
   F_mpz_mat_init(col, r, 1);
   ulong cur_col = 0;
   ulong worst_exp;
   ulong num_entries = M->c - r;

   F_mpz_t B;
   F_mpz_init(B);
// With U_exp B should be switched from r+1 to r+1 * 2^(2*U_exp) 
   F_mpz_set_ui(B, r + 1);
   if (U_exp >= 0)
      F_mpz_mul_2exp(B, B, (ulong) 2*U_exp);
   else
      F_mpz_div_2exp(B, B, (ulong) -2*U_exp);

//For the first run we'll only use 30 coeffs worth of data, should solve 99% of all 'random' polynomials
   ulong num_coeffs = 15UL;
   F_mpz_mat_t data;
   F_mpz_mat_init(data, 0, 0);


   cld_data_start = clock();
   _F_mpz_poly_factor_CLD_mat(data, F, lifted_fac, P, num_coeffs);
   cld_data_stop = clock();
   
   cld_data_total = cld_data_total + cld_data_stop - cld_data_start;
   printf(" spend a total of %f seconds on CLD stuff so far\n", (double) cld_data_total / (double)CLOCKS_PER_SEC);


   int all_coeffs = 0;
   if (data->c >= F->length - 1)
      all_coeffs = 1;
//assume that cexpo is correct for the first M->c entries, zero out the potential new entries
/* no cexpo for the moment   ulong i;
   for (i = M->c; i < M->c + 2*(F->length - 1); i++)
      cexpo[i] = 0;*/

   F_mpz_t temp;
   F_mpz_init(temp);
   ulong sqN;
//   printf("%ld sqN, %f sqrt(N)\n", sqN, sqrt( (double) (N) ) );
   int ok, col_cnt, solved,  since_last;
   ulong previously_checked;
   long newd;
   col_cnt = 0;
   solved = 0;
   since_last = 0;
   previously_checked = 0;
   while ((all_coeffs != 2) && (return_me == 0)){
      for (cur_col = previously_checked; cur_col < data->c - previously_checked; cur_col++){
         //aborting attempt to adjust sqN back to normal
         sqN = (ulong) sqrt( (double) N );         
         F_mpz_mul_ui(temp, data->rows[r] + cur_col, sqN);
         worst_exp = F_mpz_bits(temp);   
         for( ulong i = 0; i < r; i++)
            F_mpz_set(col->rows[i], data->rows[i] + cur_col);
   linear_alg_start = clock();
         ok = _F_mpz_mat_next_col(M, P, col, worst_exp, U_exp);
   linear_alg_stop = clock();
   
   linear_alg_total = linear_alg_total + linear_alg_stop - linear_alg_start;
         since_last++;
         if (ok != 0){
            since_last = 0;
            printf(" spend a total of %f seconds on Linear alg stuff since the last entry\n", (double) linear_alg_total / (double)CLOCKS_PER_SEC);
            num_entries++;
//            F_mpz_add_ui(B, B, r/2);
//            cexpo[r + col_cnt] = 0;
   lll_start = clock();
            newd = LLL_wrapper_with_removal(M, B);
   lll_stop = clock();
   
   lll_total = lll_total + lll_stop - lll_start;
   printf(" spend a total of %f seconds on LLL so far\n", (double) lll_total / (double)CLOCKS_PER_SEC);

            F_mpz_mat_resize(M, newd, M->c);
            col_cnt++;
//         This next line is what makes it 'gradual'... could try to prove that doing the same column twice won't add another P
//         But it's all the same
            cur_col--;
/*            if (M->r > 20)
            {
               newd = F_mpz_mat_check_rest(M, P, col, worst_exp);
               F_mpz_mat_resize(M, newd, M->c);               
            }
*/
            if (newd == 1){
               F_mpz_poly_factor_insert(final_fac, F, exp);
               return_me = 1;
               solved = 1;
               break;
            }
            solved =  _F_mpz_mat_check_if_solved(M, r, final_fac, lifted_fac, F, P, exp, lc);
            if (solved == 1){
               return_me = 1;
               break;
            }
         }
      }
      previously_checked = num_coeffs;
      if (solved == 0){
//This is reached when the data wasn't large enough to need LLL, this means that you had a super easy or super hard factorization
         solved =  _F_mpz_mat_check_if_solved(M, r, final_fac, lifted_fac, F, P, exp, lc);
         if (solved == 1){
//This is the easy case...
            return_me = 1;
         }
         else{
//This is the hard one...
            if (all_coeffs == 1){
               all_coeffs = 2;
//This is the worst case, all coeffs have been used and we still haven't solved the problem so more Hensel lifting needed
            }
            else{
//This condition should include a special case for when P is ridiculously large (for the sake of complexity proofs) although no example has ever needed it...
               if (since_last < N/2)
               {
                  num_coeffs = num_coeffs * 4;
   cld_data_start = clock();
                  _F_mpz_poly_factor_CLD_mat(data, F, lifted_fac, P, num_coeffs);
   cld_data_stop = clock();
   
   cld_data_total = cld_data_total + cld_data_stop - cld_data_start;
   printf(" spend a total of %f seconds on CLD stuff so far\n", (double) cld_data_total / (double)CLOCKS_PER_SEC);
                  if (data->c >= F->length - 1)
                     all_coeffs = 1;
               }
               else
               {
                  printf("maybe random would be nice but under estimated...\n");
                  num_coeffs = F->length - 1;
   cld_data_start = clock();
                  _F_mpz_poly_factor_CLD_mat(data, F, lifted_fac, P, num_coeffs);
   cld_data_stop = clock();
   
   cld_data_total = cld_data_total + cld_data_stop - cld_data_start;
   printf(" spend a total of %f seconds on CLD stuff so far\n", (double) cld_data_total / (double)CLOCKS_PER_SEC);
                  if (data->c >= F->length - 1)
                     all_coeffs = 1;               }
            }
         }
      }
   }
//   F_mpz_mat_print_pretty(M);
//   F_mpz_poly_factor_print(final_fac);
   F_mpz_clear(B);
   F_mpz_clear(lc);
   F_mpz_clear(temp);
   F_mpz_mat_clear(data);
   F_mpz_mat_clear(col);
   return return_me;
}

