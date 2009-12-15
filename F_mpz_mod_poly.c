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

#include "F_mpz_poly.h"
#include "F_mpz_mod_poly.h"
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
   {
      F_mpz_set_mpz(F_poly->coeffs + i, m_poly->coeffs[i]);
      F_mpz_mod(F_poly->coeffs + i, F_poly->coeffs + i, F_poly->P);
   }
   _F_mpz_mod_poly_normalise(F_poly);
}

void F_mpz_mod_poly_to_mpz_poly(mpz_poly_t m_poly, const F_mpz_mod_poly_t F_poly)
{
	mpz_poly_ensure_alloc(m_poly, F_poly->length);

   m_poly->length = F_poly->length;
   
   for (ulong i = 0; i < F_poly->length; i++)
	   F_mpz_get_mpz(m_poly->coeffs[i], F_poly->coeffs + i);
}

void F_mpz_poly_to_F_mpz_mod_poly(F_mpz_mod_poly_t F_poly, const F_mpz_poly_t poly)
{
	F_mpz_mod_poly_fit_length(F_poly, poly->length);

	_F_mpz_mod_poly_set_length(F_poly, poly->length);
   
	for (ulong i = 0; i < poly->length; i++)
   {
      F_mpz_set(F_poly->coeffs + i, poly->coeffs + i);
      F_mpz_mod(F_poly->coeffs + i, F_poly->coeffs + i, F_poly->P);
   }
   _F_mpz_mod_poly_normalise(F_poly);
}

void F_mpz_mod_poly_to_F_mpz_poly(F_mpz_poly_t poly, const F_mpz_mod_poly_t F_poly)
{
	F_mpz_poly_fit_length(poly, F_poly->length);

   poly->length = F_poly->length;
   
   for (ulong i = 0; i < F_poly->length; i++)
	   F_mpz_set(poly->coeffs + i, F_poly->coeffs + i);
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

void F_mpz_mod_poly_set(F_mpz_mod_poly_t poly1, const F_mpz_mod_poly_t poly2)
{
	if (poly1 != poly2) // aliasing is trivial
	{
		ulong length = poly2->length;
	
	   F_mpz_mod_poly_fit_length(poly1, poly2->length);

		for (ulong i = 0; i < poly2->length; i++)
			F_mpz_set(poly1->coeffs + i, poly2->coeffs + i);
		
		_F_mpz_mod_poly_set_length(poly1, poly2->length);
	}
}

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

   Comparison

****************************************************************************/

int F_mpz_mod_poly_equal(const F_mpz_mod_poly_t poly1, const F_mpz_mod_poly_t poly2)
{
   if (poly1 == poly2) return 1; // same polynomial

	if (poly1->length != poly2->length) return 0; // check if lengths the same

	for (ulong i = 0; i < poly1->length; i++) // check if coefficients the same
		if (!F_mpz_equal(poly1->coeffs + i, poly2->coeffs + i)) 
		   return 0;

	return 1;
}

/*===============================================================================

	Shifting

================================================================================*/

void F_mpz_mod_poly_left_shift(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t poly, const ulong n)
{
   if (n == 0) // special case, no shift
	{
		if (res != poly) F_mpz_mod_poly_set(res, poly);
		return;
	}
	
	if (poly->length == 0) // nothing to shift
	{
		_F_mpz_mod_poly_set_length(res, 0);
		return;
	}
	
	F_mpz_mod_poly_fit_length(res, poly->length + n);
	
	// copy in reverse order to avoid writing over unshifted coeffs
	for (long i = poly->length - 1; i >= 0; i--) 
		F_mpz_set(res->coeffs + i + n, poly->coeffs + i);

   // insert n zeroes
	for (ulong i = 0; i < n; i++) F_mpz_zero(res->coeffs + i);
   
   _F_mpz_mod_poly_set_length(res, poly->length + n);
}

void F_mpz_mod_poly_right_shift(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t poly, const ulong n)
{
   if (poly->length <= n) 
   {
      F_mpz_mod_poly_zero(res);
      return;
   }

   F_mpz_mod_poly_fit_length(res, poly->length - n);
	
	// copy in forward order to avoid writing over unshifted coeffs
	for (ulong i = 0; i < poly->length - n; i++) 
		F_mpz_set(res->coeffs + i, poly->coeffs + i + n);
	
	_F_mpz_mod_poly_set_length(res, poly->length - n);
}

/****************************************************************************

   Add/sub

****************************************************************************/

void _F_mpz_mod_poly_add(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2)
{
   F_mpz_poly_t p1, p2, r;

   _F_mpz_poly_attach_F_mpz_mod_poly(p1, pol1);
   _F_mpz_poly_attach_F_mpz_mod_poly(p2, pol2);
   _F_mpz_poly_attach_F_mpz_mod_poly(r, res);

   _F_mpz_poly_add(r, p1, p2);
   for (ulong i = 0; i < r->length; i++)
   {
      if (F_mpz_cmpabs(r->coeffs + i, res->P) >= 0)
         F_mpz_sub(r->coeffs + i, r->coeffs + i, res->P);
   }
   
   _F_mpz_mod_poly_attach_F_mpz_poly(res, r);
   _F_mpz_mod_poly_normalise(res);
}

void _F_mpz_mod_poly_sub(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2)
{
   F_mpz_poly_t p1, p2, r;

   _F_mpz_poly_attach_F_mpz_mod_poly(p1, pol1);
   _F_mpz_poly_attach_F_mpz_mod_poly(p2, pol2);
   _F_mpz_poly_attach_F_mpz_mod_poly(r, res);

   _F_mpz_poly_sub(r, p1, p2);
   for (ulong i = 0; i < r->length; i++)
   {
      if (F_mpz_sgn(r->coeffs + i) < 0)
         F_mpz_add(r->coeffs + i, r->coeffs + i, res->P);
   }
   
   _F_mpz_mod_poly_attach_F_mpz_poly(res, r);
   _F_mpz_mod_poly_normalise(res);
}

void F_mpz_mod_poly_add(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t poly1, const F_mpz_mod_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);

	F_mpz_mod_poly_fit_length(res, longer);
	
	_F_mpz_mod_poly_add(res, poly1, poly2);
}

void F_mpz_mod_poly_sub(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t poly1, const F_mpz_mod_poly_t poly2)
{
   ulong longer = FLINT_MAX(poly1->length, poly2->length);

	F_mpz_mod_poly_fit_length(res, longer);
	
	_F_mpz_mod_poly_sub(res, poly1, poly2);
}

/****************************************************************************

   Scalar multiplication

****************************************************************************/

void F_mpz_mod_poly_scalar_mul(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_t x)
{
   F_mpz_poly_t p1, r;

   F_mpz_mod_poly_fit_length(res, pol1->length);
   
   _F_mpz_poly_attach_F_mpz_mod_poly(p1, pol1);
   _F_mpz_poly_attach_F_mpz_mod_poly(r, res);

   F_mpz_poly_scalar_mul(r, p1, x);
   _F_mpz_poly_reduce_coeffs(r, res->P);
   
   _F_mpz_mod_poly_attach_F_mpz_poly(res, r);
   _F_mpz_mod_poly_normalise(res);
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

void _F_mpz_mod_poly_mul_trunc_left(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t pol1, const F_mpz_mod_poly_t pol2, ulong trunc)
{
   F_mpz_poly_t p1, p2, r;

   if (trunc + 1 > pol1->length + pol2->length) trunc = pol1->length + pol2->length - 1;
   if (!pol1->length && !pol2->length) trunc = 0;

   _F_mpz_poly_attach_F_mpz_mod_poly(p1, pol1);
   _F_mpz_poly_attach_F_mpz_mod_poly(p2, pol2);
   _F_mpz_poly_attach_F_mpz_mod_poly(r, res);

   _F_mpz_poly_mul_trunc_left(r, p1, p2, trunc);
   _F_mpz_poly_reduce_coeffs(r, res->P);

   _F_mpz_mod_poly_attach_F_mpz_poly(res, r);
   _F_mpz_mod_poly_normalise(res);
}

void F_mpz_mod_poly_mul_trunc_left(F_mpz_mod_poly_t res, const F_mpz_mod_poly_t poly1, const F_mpz_mod_poly_t poly2, ulong trunc)
{
   if ((poly1->length == 0) || (poly2->length == 0) || (poly1->length + poly2->length <= trunc + 1)) // special case if either poly is zero
   {
      F_mpz_mod_poly_zero(res);
      return;
   }

	if ((poly1 == res) || (poly2 == res)) // aliased inputs
	{
		F_mpz_mod_poly_t output; // create temporary
		F_mpz_mod_poly_init2(output, res->P, poly1->length + poly2->length - 1);
		if (poly1->length >= poly2->length) _F_mpz_mod_poly_mul_trunc_left(output, poly1, poly2, trunc);
		else _F_mpz_mod_poly_mul_trunc_left(output, poly2, poly1, trunc);
		F_mpz_mod_poly_swap(output, res); // swap temporary with real output
		F_mpz_mod_poly_clear(output);
	} else // ordinary case
	{
		F_mpz_mod_poly_fit_length(res, poly1->length + poly2->length - 1);
      if (poly1->length >= poly2->length) _F_mpz_mod_poly_mul_trunc_left(res, poly1, poly2, trunc);
		else _F_mpz_mod_poly_mul_trunc_left(res, poly2, poly1, trunc);
	}		
}

/****************************************************************************

   Division

****************************************************************************/

void F_mpz_mod_poly_divrem_basecase(F_mpz_mod_poly_t Q, F_mpz_mod_poly_t R, const F_mpz_mod_poly_t A, const F_mpz_mod_poly_t B)
{
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      F_mpz_mod_poly_set(R, A);
      F_mpz_mod_poly_zero(Q);
      
      return;
   }

   F_mpz_t lead_inv;
   F_mpz_init(lead_inv);
   F_mpz_invert(lead_inv, B->coeffs + B->length - 1, B->P);
   
   F_mpz * coeff_Q;
   
   F_mpz_mod_poly_t qB;
   F_mpz_mod_poly_init2(qB, B->P, B->length);
   
   F_mpz_mod_poly_t Bm1;
   _F_mpz_mod_poly_attach_truncate(Bm1, B, B->length - 1);
   
   long coeff = A->length - 1;
   
   F_mpz_mod_poly_set(R, A);
   
   if (A->length >= B->length)
   {
      F_mpz_mod_poly_fit_length(Q, A->length - B->length + 1);
      _F_mpz_mod_poly_set_length(Q, A->length - B->length + 1);
   } else F_mpz_mod_poly_zero(Q); 

   coeff_Q = Q->coeffs - B->length + 1;
   
   while (coeff >= (long) B->length - 1)
   {
      while ((coeff >= (long) B->length - 1) && (F_mpz_is_zero(R->coeffs + coeff)))
      {
         F_mpz_zero(coeff_Q + coeff);
         coeff--;
      }
      
      if (coeff >= (long) B->length - 1)
      {
         F_mpz_mulmod2(coeff_Q + coeff, R->coeffs + coeff, lead_inv, B->P);
         
         F_mpz_mod_poly_scalar_mul(qB, Bm1, coeff_Q + coeff);
         
         F_mpz_mod_poly_t R_sub;
         F_mpz_init(R_sub->P);
         F_mpz_set(R_sub->P, B->P);
         R_sub->coeffs = R->coeffs + coeff - B->length + 1;
         R_sub->length = B->length - 1;
         _F_mpz_mod_poly_sub(R_sub, R_sub, qB);
         F_mpz_clear(R_sub->P);
         
         coeff--;
      }
   }
   
   _F_mpz_mod_poly_set_length(R, B->length - 1);
   _F_mpz_mod_poly_normalise(R);
   F_mpz_mod_poly_clear(qB);
   F_mpz_clear(lead_inv);
}

void F_mpz_mod_poly_div_divconquer_recursive(F_mpz_mod_poly_t Q, F_mpz_mod_poly_t BQ, const F_mpz_mod_poly_t A, const F_mpz_mod_poly_t B)
{
   if (A->length < B->length)
   {
      F_mpz_mod_poly_zero(Q);
      F_mpz_mod_poly_zero(BQ);

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
      
      F_mpz_mod_poly_t Rb;
      F_mpz_mod_poly_init(Rb, B->P);
      F_mpz_mod_poly_divrem_basecase(Q, Rb, A, B);
      F_mpz_mod_poly_fit_length(BQ, A->length);
      F_mpz_mod_poly_sub(BQ, A, Rb);
      F_mpz_mod_poly_clear(Rb);
      
      return;
   }
   
   F_mpz_mod_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   ulong n1 = (B->length + 1)/2;
   ulong n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _F_mpz_mod_poly_attach_shift(d1, B, n2);
   _F_mpz_mod_poly_attach_truncate(d2, B, n2);
   _F_mpz_mod_poly_attach_shift(d3, B, n1);
   _F_mpz_mod_poly_attach_truncate(d4, B, n1);
   
   if (A->length < 2*B->length - 1)
   {
      /* Convert unbalanced division into a 2*q - 1 by q division */
      F_mpz_mod_poly_t t_A, t_B, t_B2;
      
      ulong q = A->length - B->length + 1;
      ulong q2 = B->length - q;

      _F_mpz_mod_poly_attach_shift(t_A, A, A->length - 2*q + 1);
      _F_mpz_mod_poly_attach_shift(t_B, B, q2);
      _F_mpz_mod_poly_attach_truncate(t_B2, B, q2);
      
      F_mpz_mod_poly_init(d1q1, B->P);
      F_mpz_mod_poly_div_divconquer_recursive(Q, d1q1, t_A, t_B); 
      
      /*
         Compute d2q1 = Q*t_B2
         It is of length q2*q terms
      */
      
      F_mpz_mod_poly_init(d2q1, B->P);
      F_mpz_mod_poly_mul(d2q1, Q, t_B2);
      
      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length at most n1+n2-1
      */
      
      F_mpz_mod_poly_fit_length(BQ, FLINT_MAX(d1q1->length + B->length - q, d2q1->length));
      F_mpz_mod_poly_left_shift(BQ, d1q1, B->length - q);
      F_mpz_mod_poly_clear(d1q1);
      _F_mpz_mod_poly_add(BQ, BQ, d2q1);
      F_mpz_mod_poly_clear(d2q1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      ulong shift = A->length - 2*B->length + 1;
      _F_mpz_mod_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = d1*q1 is length at most 2*B->length-1
      */
      
      F_mpz_mod_poly_init(d1q1, B->P);
      F_mpz_mod_poly_init(q1, Q->P);
      
      F_mpz_mod_poly_div_divconquer_recursive(q1, d1q1, p1, B); 
      
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      F_mpz_mod_poly_init(dq1, B->P);
      
      F_mpz_mod_poly_fit_length(dq1, d1q1->length + shift);
      F_mpz_mod_poly_left_shift(dq1, d1q1, shift);
      F_mpz_mod_poly_clear(d1q1);
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      F_mpz_mod_poly_init(t, A->P);
      F_mpz_mod_poly_sub(t, A, dq1);
      F_mpz_mod_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      F_mpz_mod_poly_init(q2, Q->P);
      F_mpz_mod_poly_init(dq2, Q->P);
      F_mpz_mod_poly_div_divconquer_recursive(q2, dq2, t, B); 
      F_mpz_mod_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      F_mpz_mod_poly_fit_length(Q, FLINT_MAX(q1->length + shift, q2->length));
      
      F_mpz_mod_poly_left_shift(Q, q1, shift);
      F_mpz_mod_poly_clear(q1);
      F_mpz_mod_poly_add(Q, Q, q2);
      F_mpz_mod_poly_clear(q2);
      
      /*
         Write out BQ = dq1 + dq2
      */
      
      F_mpz_mod_poly_fit_length(BQ, FLINT_MAX(dq1->length, dq2->length));
      
      F_mpz_mod_poly_add(BQ, dq1, dq2);
      F_mpz_mod_poly_clear(dq1);
      F_mpz_mod_poly_clear(dq2);
      
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
      
   _F_mpz_mod_poly_attach_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length n1
      d1q1 = d1*q1 is length at most 2*n1-1
   */
      
   F_mpz_mod_poly_init(d1q1, B->P);
   F_mpz_mod_poly_init(q1, B->P);
   F_mpz_mod_poly_div_divconquer_recursive(q1, d1q1, p1, d1); 
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being at most length n1+n2-1
   */  
   
   F_mpz_mod_poly_init(d2q1, B->P);
   F_mpz_mod_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most 2*n1+n2-1
   */
   
   F_mpz_mod_poly_init2(dq1, B->P, FLINT_MAX(d1q1->length + n2, d2q1->length));
   F_mpz_mod_poly_left_shift(dq1, d1q1, n2);
   F_mpz_mod_poly_clear(d1q1);
   _F_mpz_mod_poly_add(dq1, dq1, d2q1);
   F_mpz_mod_poly_clear(d2q1);
   
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1
      which has length at most 2*n1+n2-1, but we are not interested 
      in up to the first n1 coefficients, so it has 
      effective length at most n1+n2-1
   */
   
   F_mpz_mod_poly_init2(t, A->P, FLINT_MAX(A->length-n2, dq1->length));
   F_mpz_mod_poly_right_shift(t, A, n2);
   _F_mpz_mod_poly_sub(t, t, dq1);
   F_mpz_mod_poly_truncate(t, B->length - 1);
   
   /*
      Compute q2 = t div d1
      It is at most an n1+n2-1 by n1 division, so
      the length of q2 will be at most n2
      Also compute d1q2 of length at most n1+n2-1
   */
   
   F_mpz_mod_poly_init(d1q2, B->P);
   F_mpz_mod_poly_init(q2, Q->P);
   F_mpz_mod_poly_div_divconquer_recursive(q2, d1q2, t, d1); 
   F_mpz_mod_poly_clear(t);
      
   /*
      Compute d2q2 = d2*q2 which is of length 
      at most n1+n2-1
   */
   
   F_mpz_mod_poly_init(d2q2, A->P);
   F_mpz_mod_poly_mul(d2q2, d2, q2);
   
   /*
      Compute dq2 = d1*q2*x^n2 + d2q2
      which is of length at most n1+2*n2-1
   */
   
   F_mpz_mod_poly_init2(dq2, A->P, FLINT_MAX(d1q2->length+n2, d2q2->length));
   F_mpz_mod_poly_left_shift(dq2, d1q2, n2);
   F_mpz_mod_poly_clear(d1q2);
   _F_mpz_mod_poly_add(dq2, dq2, d2q2);
   F_mpz_mod_poly_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length at most n1+n2
   */
   
   F_mpz_mod_poly_fit_length(Q, FLINT_MAX(q1->length+n2, q2->length));
   F_mpz_mod_poly_left_shift(Q, q1, n2);
   F_mpz_mod_poly_clear(q1);
   _F_mpz_mod_poly_add(Q, Q, q2);
   F_mpz_mod_poly_clear(q2);
   
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length at most 2*(n1+n2)-1
   */
   
   F_mpz_mod_poly_fit_length(BQ, FLINT_MAX(n2 + dq1->length, dq2->length));
   F_mpz_mod_poly_left_shift(BQ, dq1, n2);
   _F_mpz_mod_poly_add(BQ, BQ, dq2);
   
   F_mpz_mod_poly_clear(dq2);
   F_mpz_mod_poly_clear(dq1);
}

void F_mpz_mod_poly_divrem_divconquer(F_mpz_mod_poly_t Q, F_mpz_mod_poly_t R, const F_mpz_mod_poly_t A, const F_mpz_mod_poly_t B)
{
   F_mpz_mod_poly_t QB;
   
   F_mpz_mod_poly_init(QB, Q->P);
   
   F_mpz_mod_poly_div_divconquer_recursive(Q, QB, A, B);
   
   F_mpz_mod_poly_fit_length(R, A->length);
   _F_mpz_mod_poly_sub(R, A, QB);
   _F_mpz_mod_poly_normalise(R);
   
   F_mpz_mod_poly_clear(QB);
}

/*****************************************************************

   Wrappers to get this sheet working while Bill writes the full versions

******************************************************************/

void F_mpz_mod_poly_rem(F_mpz_mod_poly_t R, F_mpz_mod_poly_t A, F_mpz_mod_poly_t B)
{
   F_mpz_mod_poly_t Q;
   F_mpz_mod_poly_init(Q, A->P);

   F_mpz_mod_poly_divrem_divconquer(Q, R, A, B);

   F_mpz_mod_poly_clear(Q);
}

void F_mpz_mod_poly_mulmod( F_mpz_mod_poly_t res, F_mpz_mod_poly_t A, F_mpz_mod_poly_t B, F_mpz_mod_poly_t C)
{

   F_mpz_mod_poly_t Q;
   F_mpz_mod_poly_init(Q, A->P);

   F_mpz_mod_poly_mul(Q, A, B);
   F_mpz_mod_poly_rem(res, Q, C);

   F_mpz_mod_poly_clear(Q);
}


