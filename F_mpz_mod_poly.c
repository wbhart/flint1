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
/*****************************************************************************

   F_mpz_mod_poly.c: Polynomials over F_mpz mod p.
   
   Copyright (C) 2009, William Hart, Andy Novocin

*****************************************************************************/

#include "F_mpz_mod_poly.h"
#include "F_mpz_poly.h"
#include "long_extras.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "flint.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void F_mpz_mod_poly_init(F_mpz_mod_poly_t poly, const F_mpz_t P)
{
   poly->coeffs = NULL;
   F_mpz_init(poly->P);
   F_mpz_set(poly->P, P);
   poly->alloc = 0;
   poly->length = 0;
}

void F_mpz_mod_poly_init2(F_mpz_mod_poly_t poly, const F_mpz_t P, const ulong alloc)
{
   if (alloc) // allocate space for alloc small coeffs
   {
      poly->coeffs = (F_mpz *) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
   }
   else poly->coeffs = NULL;

   F_mpz_init(poly->P);
   F_mpz_set(poly->P, P);   
   poly->alloc = alloc;
   poly->length = 0;
}

void F_mpz_mod_poly_clear(F_mpz_mod_poly_t poly)
{
   for (ulong i = 0; i < poly->alloc; i++) // Clean up any mpz_t's
		_F_mpz_demote(poly->coeffs + i);
	if (poly->coeffs) flint_heap_free(poly->coeffs); // clean up ordinary coeffs
   F_mpz_clear(poly->P);
}

void F_mpz_mod_poly_realloc(F_mpz_mod_poly_t poly, const ulong alloc)
{
   if (!alloc) // alloc == 0, clear up
   {
         F_mpz_mod_poly_clear(poly);
			return;
   }  
   
	if (poly->alloc) // realloc
	{
		F_mpz_mod_poly_truncate(poly, alloc);

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

void F_mpz_mod_poly_fit_length(F_mpz_mod_poly_t poly, const ulong length)
{
   ulong alloc = length;
   
	if (alloc <= poly->alloc) return;

   // at least double number of allocated coeffs
	if (alloc < 2*poly->alloc) alloc = 2*poly->alloc; 
   
   F_mpz_mod_poly_realloc(poly, alloc);
}

/****************************************************************************

   Normalisation/truncation

****************************************************************************/

void _F_mpz_mod_poly_normalise(F_mpz_mod_poly_t poly)
{
   ulong length = poly->length;
	
	while ((length) && (!poly->coeffs[length - 1])) length--;

	poly->length = length;
}

/****************************************************************************

   Conversion

****************************************************************************/

void mpz_poly_to_F_mpz_mod_poly(F_mpz_mod_poly_t F_poly, const mpz_poly_t m_poly)
{
	F_mpz_mod_poly_fit_length(F_poly, m_poly->length);

	_F_mpz_mod_poly_set_length(F_poly, m_poly->length);
   
	for (ulong i = 0; i < m_poly->length; i++)
		F_mpz_set_mpz(F_poly->coeffs + i, m_poly->coeffs[i]);
}

void F_mpz_mod_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_mod_poly_t F_poly)
{
	mpz_poly_ensure_alloc(m_poly, F_poly->length);

   m_poly->length = F_poly->length;
   
   for (ulong i = 0; i < F_poly->length; i++)
	   F_mpz_get_mpz(m_poly->coeffs[i], F_poly->coeffs + i);
}

void zmod_poly_to_F_mpz_mod_poly(F_mpz_mod_poly_t fpol, const zmod_poly_t zpol)
{
   if (zpol->length == 0)
   {
      F_mpz_mod_poly_zero(fpol);
      return;
   }

   F_mpz_mod_poly_fit_length(fpol, zpol->length);
   
   for (ulong i = 0; i < zpol->length; i++)
      F_mpz_set_ui(fpol->coeffs + i, zpol->coeffs[i]);

   _F_mpz_mod_poly_set_length(fpol, zpol->length);
   _F_mpz_mod_poly_normalise(fpol);
}

/****************************************************************************

   Assignment/swap

****************************************************************************/

void F_mpz_mod_poly_swap(F_mpz_mod_poly_t poly1, F_mpz_mod_poly_t poly2)
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

/****************************************************************************

   Multiplication

****************************************************************************/

void _F_mpz_mod_poly_mul(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2)
{
   F_mpz_poly_t p1, p2, r;

   _F_mpz_poly_attach_F_mpz_mod_poly(p1, pol1);
   _F_mpz_poly_attach_F_mpz_mod_poly(p2, pol2);
   _F_mpz_poly_attach_F_mpz_mod_poly(r, res);

   _F_mpz_poly_mul(r, p1, p2);
   _F_mpz_poly_reduce_coeffs(r, res->P);

   _F_mpz_mod_poly_attach_F_mpz_poly(res, r);
   _F_mpz_mod_poly_normalise(res);
}

void F_mpz_mod_poly_mul(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2)
{
	if ((pol1->length == 0) || (pol2->length == 0)) // special case if either poly is zero
   {
      F_mpz_mod_poly_zero(res);
      return;
   }

	if ((pol1 == res) || (pol2 == res)) // aliased inputs
	{
		F_mpz_mod_poly_t output; // create temporary
		F_mpz_mod_poly_init2(output, pol1->P, pol1->length + pol2->length - 1);
		if (pol1->length >= pol2->length) _F_mpz_mod_poly_mul(output, pol1, pol2);
		else _F_mpz_mod_poly_mul(output, pol2, pol1);
		F_mpz_mod_poly_swap(output, res); // swap temporary with real output
		F_mpz_mod_poly_clear(output);
	} else // ordinary case
	{
		F_mpz_mod_poly_fit_length(res, pol1->length + pol2->length - 1);
      if (pol1->length >= pol2->length) _F_mpz_mod_poly_mul(res, pol1, pol2);
		else _F_mpz_mod_poly_mul(res, pol2, pol1);
	}		
}

