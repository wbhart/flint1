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

/*
    F_mpz_mpoly.c: Multivariable polynomials over Z (FLINT 2.0 polynomials)

    Copyright (C) 2008, 2009 William Hart 
*/

#include <stdint.h>
#include <string.h>
#include <math.h>

#include "flint.h"
#include "F_mpz.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "packed_vec.h"
#include "F_mpz_mpoly.h"

/*===============================================================================

	Memory management

================================================================================*/

void F_mpz_mpoly_init(F_mpz_mpoly_t poly, ord_t ordering)
{
   poly->coeffs = NULL;
	poly->exps = NULL;
	poly->packed = NULL;
   
   poly->alloc = 0;
   poly->length = 0;
	poly->vars = 0;
	poly->ordering = ordering;
	poly->small = 1; // start with packed exponents
	poly->packed_bits = 8; // of 8 bits each
}

void F_mpz_mpoly_init2(F_mpz_mpoly_t poly, const ulong alloc, 
							                   const ulong vars, ord_t ordering)
{
   if (alloc) // allocate space for alloc small coeffs
   {
      poly->coeffs = (F_mpz *) flint_heap_alloc(alloc);
		poly->packed = (ulong *) flint_heap_alloc(alloc);
		F_mpn_clear(poly->coeffs, alloc);
		F_mpn_clear(poly->packed, alloc);
   }
   else poly->coeffs = NULL;

	/*if (vars)
	{
		poly->exps = (pv_s *) flint_heap_alloc_bytes(vars*sizeof(pv_s));
		ulong i;
		for (i = 0; i < vars; i++) pv_init(poly->exps + i, alloc, 0);
	} else poly->exps = NULL;*/
   
   poly->alloc = alloc;
	poly->vars = vars;
   poly->length = 0;
	poly->ordering = ordering;
	poly->small = 1; // start with packed exponents
   poly->packed_bits = 8; // of 8 bits each
}

void F_mpz_mpoly_realloc(F_mpz_mpoly_t poly, const ulong alloc, const ulong vars)
{
   if (!alloc || !vars) // alloc == 0, clear up
   {
         F_mpz_mpoly_clear(poly);
			return;
   }  
   
	if (poly->alloc != alloc)
	{
		if (poly->alloc) // realloc
	   {
		   _F_mpz_mpoly_truncate(poly, alloc);
         
		   poly->coeffs = (F_mpz *) flint_heap_realloc(poly->coeffs, alloc);
		   poly->packed = (ulong *) flint_heap_realloc(poly->packed, alloc);
		   if (alloc > poly->alloc)
			{
				F_mpn_clear(poly->coeffs + poly->alloc, alloc - poly->alloc);
            F_mpn_clear(poly->packed + poly->alloc, alloc - poly->alloc);
			}

	   } else // nothing allocated already so do it now
	   {
		   poly->coeffs = (F_mpz *) flint_heap_alloc(alloc);
		   poly->packed = (ulong *) flint_heap_alloc(alloc);
		   F_mpn_clear(poly->coeffs, alloc);
		   F_mpn_clear(poly->packed, alloc);
	   }
		
		if (!poly->small) 
			ulong i;
			for (i = 0; i < FLINT_MIN(vars, poly->vars); i++)
			   pv_realloc(poly->exps + i, alloc);
	}
   
   poly->alloc = alloc;

   if (vars != poly->vars)
	{
		if (poly->vars) // realloc
		{
			if (vars < poly->vars)
			   ulong i;
			   for (i = vars; i < poly->vars; i++)
					pv_clear(poly->exps + i);

			poly->exps = (pv_s *) flint_heap_realloc(poly->exps, vars*sizeof(pv_s));

			if (vars > poly->vars)
			   ulong i;
			   for (i = poly->vars; i < vars; i++)
				   pv_init(poly->exps + i, alloc, 0);
		} else // nothing allocated already so do it now
		{
			poly->exps = (pv_s *) flint_heap_alloc_bytes(vars*sizeof(pv_s));
			ulong i;
			for (i = 0; i < vars; i++)
				   pv_init(poly->exps, alloc, 0);
		}  
	}

	poly->vars = vars;
}

void F_mpz_mpoly_fit_length(F_mpz_mpoly_t poly, const ulong length)
{
   ulong alloc = length;
   
	if (alloc <= poly->alloc) return;

   // at least double number of allocated coeffs
	if (alloc < 2*poly->alloc) alloc = 2*poly->alloc; 
   
   F_mpz_mpoly_realloc(poly, alloc, poly->vars);
}

void F_mpz_mpoly_clear(F_mpz_mpoly_t poly)
{
   ulong i;
   for (i = 0; i < poly->alloc; i++) // Clean up any mpz_t's
		_F_mpz_demote(poly->coeffs + i);
	if (poly->coeffs) flint_heap_free(poly->coeffs); // clean up coeffs
   if (poly->coeffs) flint_heap_free(poly->packed); // and packed

	if (!poly->small)
	{
		ulong i;
		for (i = 0; i < poly->vars; i++) // and vars
		   pv_clear(poly->exps + i);
	   if (poly->vars) flint_heap_free(poly->exps);
	}
}

/*===============================================================================

	Set/get

================================================================================*/

void F_mpz_mpoly_set_coeff_ui(F_mpz_mpoly_t poly, const ulong n, const ulong x)
{
   F_mpz_mpoly_fit_length(poly, n+1);

   if (n + 1 > poly->length) // insert zeroes between end of poly and new coeff if needed
   {
      long i;
      for (i = poly->length; i + 1 < n; i++)
         F_mpz_zero(poly->coeffs + i); 
      poly->length = n+1;
   }

   F_mpz_set_ui(poly->coeffs + n, x);
}

ulong F_mpz_mpoly_get_coeff_ui(const F_mpz_mpoly_t poly, const ulong n)
{
   if (n + 1 > poly->length) // coefficient is beyond end of polynomial
      return 0;
   
	return F_mpz_get_ui(poly->coeffs + n);
}

void F_mpz_mpoly_set_var_exp_packed_grlex(F_mpz_mpoly_t poly, const ulong n, const ulong var, 
									                                        const ulong exp)
{
   int bits = poly->packed_bits;
	ulong shift = FLINT_BITS - bits*(var + 2); // first spot is for total degree, second spot is variable we are touching
	ulong mask = ~(((1L << bits) - 1L) << shift);
	poly->packed[n] &= mask; // clear the relevant bits
	poly->packed[n] += ((exp << (FLINT_BITS - bits)) + (exp << shift)); // update total degree and relevant exponent
}

void F_mpz_mpoly_set_var_exp(F_mpz_mpoly_t poly, const ulong n, const ulong var, 
									                                        const ulong exp)
{
   F_mpz_mpoly_fit_length(poly, n+1);
   F_mpz_mpoly_fit_vars(poly, var+1);

   if ((poly->small) && (FLINT_BIT_COUNT(exp) <= poly->packed_bits)) // monomials can be packed
	{
		switch (poly->ordering)
		{
		   case GRLEX:
			   F_mpz_mpoly_set_var_exp_packed_grlex(poly, n, var, exp);
			break;

			default:
			   abort(); // not implemented
		}
	} else if (poly->small)
	{
		abort(); // changing packed_bits or converting to big poly not yet implemented
	} else // poly is not small (exponents aren't packed)
	{
		int bits = pv_bit_fit(FLINT_BIT_COUNT(exp));
	   pv_s * expt = poly->exps + var;
	   if (bits > expt->bits) pv_set_bits(expt, bits);

	   PV_SET_ENTRY(*expt, n, exp);
	}
}

ulong F_mpz_mpoly_get_var_exp_packed_grlex(F_mpz_mpoly_t poly, const ulong n, const ulong var)
{
   int bits = poly->packed_bits;
	ulong shift = FLINT_BITS - bits*(var + 2); // first spot is for total degree, second spot is variable we are touching
	return ((poly->packed[n] >> shift) & ((1L << bits) - 1L));
}


ulong F_mpz_mpoly_get_var_exp(F_mpz_mpoly_t poly, const ulong n, const ulong var)
{
	if (poly->small)
	{
		switch (poly->ordering)
		{
		   case GRLEX: 
				return F_mpz_mpoly_get_var_exp_packed_grlex(poly, n, var);
			break;

			default:
				abort(); // not implemented yet
		}
	} else
	{
		if (var > poly->vars) return 0;
	   if (poly->exps[var].bits == 0) return 0;
	   if (n > poly->length) return 0;

      ulong val;
	   pv_s * expt = poly->exps + var;
	   PV_GET_ENTRY(val, *expt, n);

	   return val;
	}
}

/*===============================================================================

	Print/read

================================================================================*/

void F_mpz_mpoly_print_pretty(F_mpz_mpoly_t poly, char ** var_syms)
{
	if (poly->small)
	{
		if (poly->length == 0)
		{
			printf("0");
			return;
		}
		
		ulong bits = poly->packed_bits;
		ulong mask = ((1L << bits) - 1L);
		
		switch (poly->ordering)
		{
		   case GRLEX:
				long i;
				for (i = poly->length - 1; i >= 0; i--)
		      {
               if ((i != poly->length - 1) && (F_mpz_sgn(poly->coeffs + i) > 0)) 
						printf("+");

			      ulong exp = ((poly->packed[i] >> (FLINT_BITS - bits)) & mask);
               if ((poly->coeffs[i] == 1L) || (poly->coeffs[i] == -1L))
					{
						if (!exp) F_mpz_print(poly->coeffs + i);
					} else 
					{
						F_mpz_print(poly->coeffs + i);
						if (exp) printf("*");
					}

			      int first = 1;
					ulong j;
					for (j = 0; j < poly->vars; j++)
			      {
				      exp = ((poly->packed[i] >> (FLINT_BITS - (j+2)*bits)) & mask);
				      if (exp) 
						{
							if (!first) printf("*");
							first = 0;
						   if (exp == 1) printf("%s", var_syms[j]);
							else printf("%s^%ld", var_syms[j], exp);
						}
			      }
		      }
			break;

			default: abort(); // not implemented yet
		}
	} else
	{
		abort(); // printing not implemented for large polys
	}
}

/*===============================================================================

	Arithmetic

================================================================================*/

void _F_mpz_mpoly_add_inplace(F_mpz_mpoly_t res, ulong n, F_mpz_mpoly_t poly)
{	
	if ((poly->small) && (res->small))
	{
      ulong len = res->length - n;
		F_mpz * temp_c;
      ulong * temp_e;

		ulong z = len;
		if (len < 2100) z = 2100;
		if (len < 300) z = 300;
		if (len < 50) z = 50;
		if (z == 50)
		{
			temp_c = (F_mpz *) flint_stack_alloc_small(2*z);
         temp_e = (F_mpz *) temp_c + z;
		} else
      {
			temp_c = (F_mpz *) flint_stack_alloc(2*z);
         temp_e = (F_mpz *) temp_c + z;
		}

		F_mpz * res_c = res->coeffs;
		ulong * res_e = res->packed;

		F_mpz * poly_c = poly->coeffs;
		ulong * poly_e = poly->packed;

		ulong r_length = res->length;
      ulong p_length = poly->length;

		ulong extras = 0; // index of next temporary spot for a coefficient we have to pull out of res
      ulong start = 0; // first coefficient we haven't placed back in res
		ulong poly_i = 0; // index into poly
		ulong res_i = n; // index into res

      for ( ; poly_i < p_length; res_i++) 
		{
         //printf("res_i = %ld, poly_i = %ld, start = %ld, extras = %ld\n", res_i, poly_i, start, extras);
			if (extras != start) // we have coefficients to place
			{
            if (res_i < r_length) // check we actually have something to save
				{
					temp_c[extras] = res_c[res_i]; // we save the existing coefficient
			      res_c[res_i] = 0L;
					temp_e[extras] = res_e[res_i];
			      extras++;
				}

				if (temp_e[start] < poly_e[poly_i]) // temp coeff goes in next
				{
					res_c[res_i] = temp_c[start];
					res_e[res_i] = temp_e[start];
					start++;
				} else if (temp_e[start] == poly_e[poly_i]) // add coefficients
				{
               F_mpz_add(res_c + res_i, temp_c + start, poly_c + poly_i);
					_F_mpz_demote(poly_c + poly_i);
					_F_mpz_demote(temp_c + start);
					res_e[res_i] = temp_e[start];
					poly_i++;
               start++;

					if (res_c[res_i] == 0) // addition cancelled
					{
						if (start + 1 != extras) // we have more spare, so we can fix this
						{
					      extras--; // we will remove the last spare and place it back in the poly
							res_c[res_i] = temp_c[extras];
					      res_e[res_i] = temp_e[start];
						} else if (res_i + 1 < r_length) // we need to shift subsequent coeffs back one
						{
							F_mpn_copy_forward(res_c + res_i, res_c + res_i + 1, r_length - res_i - 1);
							res_c[r_length - 1] = 0L;
							F_mpn_copy_forward(res_e + res_i, res_e + res_i + 1, r_length - res_i - 1);
							res_e[r_length - 1] = 0L;
							r_length--;
						}
						// else we have no more coefficients from res to worry about
					   
						res_i--; // we need to revisit this coefficient
					}
				} else // poly coeff goes in next 
				{
					res_c[res_i] = poly_c[poly_i];
					poly_c[poly_i] = 0;
					//F_mpz_set(res_c + res_i, poly_c + poly_i);
					res_e[res_i] = poly_e[poly_i];
					poly_i++;
				}
			} else // it's either the poly coefficient or the existing one (if there is one)
			{
				//printf("%ld, %ld\n", poly_e[poly_i], res_e[res_i]);
				if (res_i >= r_length) // we've run out of coeffs in res
				{
               res_c[res_i] = poly_c[poly_i]; // copy poly coeff over
					poly_c[poly_i] = 0;
					//F_mpz_set(res_c + res_i, poly_c + poly_i);
					res_e[res_i] = poly_e[poly_i];
					poly_i++;
				} else if (poly_e[poly_i] < res_e[res_i]) // poly coeff goes in next
				{
					//printf("res_i = %ld, r_length = %ld\n", res_i, r_length);
					if (res_i < r_length) // check we have something to save
				   {
					   temp_c[extras] = res_c[res_i]; // save the existing coefficient
			         res_c[res_i] = 0L;
					   temp_e[extras] = res_e[res_i];
			         extras++;
					}

					res_c[res_i] = poly_c[poly_i];
					poly_c[poly_i] = 0;
					//F_mpz_set(res_c + res_i, poly_c + poly_i);
					res_e[res_i] = poly_e[poly_i];
					poly_i++;
				} else if (poly_e[poly_i] == res_e[res_i])  // add coefficients
				{
					F_mpz_add(res_c + res_i, res_c + res_i, poly_c + poly_i);
					_F_mpz_demote(poly_c + poly_i);
					poly_i++;

					if (res_c[res_i] == 0L) // addition cancelled
					{
						if (res_i + 1 < r_length) // we need to shift subsequent coeffs back one
						{
							F_mpn_copy_forward(res_c + res_i, res_c + res_i + 1, r_length - res_i - 1);
							res_c[r_length - 1] = 0L;
							F_mpn_copy_forward(res_e + res_i, res_e + res_i + 1, r_length - res_i - 1);
							res_e[r_length - 1] = 0L;
							r_length--;
						}
                  res_i--;  // we need to revisit this coefficient
					}
				}
			}
		}

		for ( ; extras != start; res_i++) // place remaining extras
		{
         if (res_i < r_length) // check we actually have something to save
			{
				temp_c[extras] = res_c[res_i]; // we save the existing coefficient
			   res_c[res_i] = 0;
				temp_e[extras] = res_e[res_i];
			   extras++;
			}

			res_c[res_i] = temp_c[start];
			res_e[res_i] = temp_e[start];
			start++;
		}

		if (z == 50)
		{
			flint_stack_release_small(); // temp_e
		} else
		{
			flint_stack_release(); // temp_e
		} 
		//flint_stack_release(); // temp_c

		res->length = FLINT_MAX(res_i, r_length); // All coeffs of poly may have added to res
	} else
	{
		abort(); // merge not implemented for large polynomials
	}
}

void F_mpz_mpoly_add_inplace(F_mpz_mpoly_t res, ulong n, F_mpz_mpoly_t poly)
{
   if (poly->length == 0) return; // nothing to merge
	
	F_mpz_mpoly_fit_length(res, res->length + poly->length);

	_F_mpz_mpoly_add_inplace(res, n, poly);
}
	
void _F_mpz_mpoly_mul_small_1(ulong * res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   ulong c1 = poly1->coeffs[n1];
   ulong c2 = poly2->coeffs[n2];
   
	umul_ppmm(res[1], res[0], c1, c2);
	res[2] = 0L;
}

void _F_mpz_mpoly_addmul_small_1(ulong * res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   ulong arr[2];
	ulong cry;
	
	ulong c1 = poly1->coeffs[n1];
   ulong c2 = poly2->coeffs[n2];
   
	umul_ppmm(arr[1], arr[0], c1, c2);
	add_ssaaaa(cry, res[0], 0, res[0], 0, arr[0]);
	add_ssaaaa(res[2], res[1], res[2], res[1], 0, cry);
	add_ssaaaa(res[2], res[1], res[2], res[1], 0, arr[1]);
}

void _F_mpz_mpoly_mul_small_1x1(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);
	
	res->length = 1;
}

void _F_mpz_mpoly_mul_small_2x1(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

   res_e[1] = poly1_e[1] + poly2_e[0];
	F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);

   res->length = 2;
}

void _F_mpz_mpoly_mul_small_2x2(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

	ulong p01 = poly1_e[0] + poly2_e[1];
	ulong p10 = poly1_e[1] + poly2_e[0];

	if (p01 < p10)
	{
		res_e[1] = p01;
		res_e[2] = p10;
		res_e[3] = poly1_e[1] + poly2_e[1];
		F_mpz_mul2(res_c + 1, poly1_c, poly2_c + 1);
		F_mpz_mul2(res_c + 2, poly1_c + 1, poly2_c);
		F_mpz_mul2(res_c + 3, poly1_c + 1, poly2_c + 1);
		res->length = 4;
	} else if (p01 > p10)
	{
		res_e[1] = p10;
		res_e[2] = p01;
		res_e[3] = poly1_e[1] + poly2_e[1];
		F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);
		F_mpz_mul2(res_c + 2, poly1_c, poly2_c + 1);
		F_mpz_mul2(res_c + 3, poly1_c + 1, poly2_c + 1);
		res->length = 4;
	} else
	{
		res_e[1] = p10;
		F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);
		F_mpz_addmul(res_c + 1, poly1_c, poly2_c + 1);
		if (res_c[1] == 0) // cancellation of middle terms occurred
		{
         res_e[1] = poly1_e[1] + poly2_e[1];
		   F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c + 1);
			res->length = 2;
		} else
		{
			res_e[2] = poly1_e[1] + poly2_e[1];
		   F_mpz_mul2(res_c + 2, poly1_c + 1, poly2_c + 1);
			res->length = 3;
		}
	}
}

void _F_mpz_mpoly_mul_small_3x1(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

	res_e[1] = poly1_e[1] + poly2_e[0];
	F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);

	res_e[2] = poly1_e[2] + poly2_e[0];
	F_mpz_mul2(res_c + 2, poly1_c + 2, poly2_c);

	res->length = 3;
}

void _F_mpz_mpoly_mul_small_3x2(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
	_F_mpz_mpoly_mul_small_3x1(res, poly1, n1, poly2, n2);

	F_mpz_mpoly_t temp;
	temp->packed = (ulong *) flint_stack_alloc(3);
	//F_mpn_clear(temp->packed, 3);
	temp->coeffs = (F_mpz *) flint_stack_alloc(3);
   F_mpn_clear(temp->coeffs, 3);
	temp->small = 1;
	
	_F_mpz_mpoly_mul_small_3x1(temp, poly1, n1, poly2, n2 + 1);

	_F_mpz_mpoly_add_inplace(res, 1, temp);

	//ulong i;
for (i = 0; i < temp->length; i++)
	//	_F_mpz_demote(temp->coeffs + i);

	flint_stack_release(); // coeffs
	flint_stack_release(); // packed
}

void _F_mpz_mpoly_mul_small_3x3(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

	res_e[1] = poly1_e[1] + poly2_e[0];
	F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);

	res_e[2] = poly1_e[2] + poly2_e[0];
	F_mpz_mul2(res_c + 2, poly1_c + 2, poly2_c);

	res_e[3] = poly1_e[2] + poly2_e[1];
	F_mpz_mul2(res_c + 3, poly1_c + 2, poly2_c + 1);

	res_e[4] = poly1_e[2] + poly2_e[2];
	F_mpz_mul2(res_c + 4, poly1_c + 2, poly2_c + 2);

	res->length = 5;

	ulong temp1[4];
	ulong temp2[4];
	F_mpz_mpoly_t temp;
	temp->packed = temp1;//(ulong *) flint_stack_alloc_small(4);
	//F_mpn_clear(temp->packed, 4);
	temp->coeffs = temp2;//(F_mpz *) flint_stack_alloc_small(4);
   F_mpn_clear(temp->coeffs, 4);
	temp->small = 1;
	
	_F_mpz_mpoly_mul_small_2x2(temp, poly1, n1, poly2, n2 + 1);

	_F_mpz_mpoly_add_inplace(res, 1, temp);

	//ulong i;
for (i = 0; i < temp->length; i++)
	//	_F_mpz_demote(temp->coeffs + i);

	//flint_stack_release_small(); // coeffs
	//flint_stack_release_small(); // packed
}

void _F_mpz_mpoly_mul_small_4x1(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

	res_e[1] = poly1_e[1] + poly2_e[0];
	F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);

	res_e[2] = poly1_e[2] + poly2_e[0];
	F_mpz_mul2(res_c + 2, poly1_c + 2, poly2_c);

	res_e[3] = poly1_e[3] + poly2_e[0];
	F_mpz_mul2(res_c + 3, poly1_c + 3, poly2_c);

	res->length = 4;
}

void _F_mpz_mpoly_mul_small_4x2(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
	_F_mpz_mpoly_mul_small_4x1(res, poly1, n1, poly2, n2);

	F_mpz_mpoly_t temp;
	temp->packed = (ulong *) flint_stack_alloc(4);
	//F_mpn_clear(temp->packed, 4);
	temp->coeffs = (F_mpz *) flint_stack_alloc(4);
   F_mpn_clear(temp->coeffs, 4);
	temp->small = 1;
	
	_F_mpz_mpoly_mul_small_4x1(temp, poly1, n1, poly2, n2 + 1);

	_F_mpz_mpoly_add_inplace(res, 1, temp);

	//ulong i;
for (i = 0; i < temp->length; i++)
	//	_F_mpz_demote(temp->coeffs + i);

	flint_stack_release(); // coeffs
	flint_stack_release(); // packed
}

void _F_mpz_mpoly_mul_small_4x3(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

	res_e[1] = poly1_e[1] + poly2_e[0];
	F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);

	res_e[2] = poly1_e[2] + poly2_e[0];
	F_mpz_mul2(res_c + 2, poly1_c + 2, poly2_c);

	res_e[3] = poly1_e[3] + poly2_e[0];
	F_mpz_mul2(res_c + 3, poly1_c + 3, poly2_c);

	res_e[4] = poly1_e[3] + poly2_e[1];
	F_mpz_mul2(res_c + 4, poly1_c + 3, poly2_c + 1);

	res_e[5] = poly1_e[3] + poly2_e[2];
	F_mpz_mul2(res_c + 5, poly1_c + 3, poly2_c + 2);

	res->length = 6;

	F_mpz_mpoly_t temp;
	temp->packed = (ulong *) flint_stack_alloc(6);
	//F_mpn_clear(temp->packed, 6);
	temp->coeffs = (F_mpz *) flint_stack_alloc(6);
   temp->small = 1;
	F_mpn_clear(temp->coeffs, 6);
	
	_F_mpz_mpoly_mul_small_3x2(temp, poly1, n1, poly2, n2 + 1);

	/*char * syms[4];
	syms[0] = "t";
	syms[1] = "x";
	syms[2] = "y";
   syms[3] = "z";
   temp->packed_bits = 8;
   temp->vars = 4;
   temp->ordering = GRLEX;
	temp->small = 1;
	printf("********-->length temp = %ld\n", temp->length);
	F_mpz_mpoly_print_pretty(temp, syms);printf("\n");
	ulong s;
	for (s = 0; s < res->length; s++)
	{
		F_mpz_print(res_c + s);printf("\n");
	}
	printf("\n");*/
      

	_F_mpz_mpoly_add_inplace(res, 1, temp);

	//ulong i;
for (i = 0; i < temp->length; i++)
	//	_F_mpz_demote(temp->coeffs + i);

	flint_stack_release(); // coeffs
	flint_stack_release(); // packed
}

void _F_mpz_mpoly_mul_small_4x4(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
   F_mpz * res_c = res->coeffs;
   ulong * res_e = res->packed;

   F_mpz * poly1_c = poly1->coeffs + n1;
   ulong * poly1_e = poly1->packed + n1;

   F_mpz * poly2_c = poly2->coeffs + n2;
   ulong * poly2_e = poly2->packed + n2;

	res_e[0] = poly1_e[0] + poly2_e[0];
	F_mpz_mul2(res_c, poly1_c, poly2_c);

	res_e[1] = poly1_e[1] + poly2_e[0];
	F_mpz_mul2(res_c + 1, poly1_c + 1, poly2_c);

	res_e[2] = poly1_e[2] + poly2_e[0];
	F_mpz_mul2(res_c + 2, poly1_c + 2, poly2_c);

	res_e[3] = poly1_e[3] + poly2_e[0];
	F_mpz_mul2(res_c + 3, poly1_c + 3, poly2_c);

	res_e[4] = poly1_e[3] + poly2_e[1];
	F_mpz_mul2(res_c + 4, poly1_c + 3, poly2_c + 1);

	res_e[5] = poly1_e[3] + poly2_e[2];
	F_mpz_mul2(res_c + 5, poly1_c + 3, poly2_c + 2);

	res_e[6] = poly1_e[3] + poly2_e[3];
	F_mpz_mul2(res_c + 6, poly1_c + 3, poly2_c + 3);

	res->length = 7;

	ulong temp1[9];
	ulong temp2[9];
	F_mpz_mpoly_t temp;
	temp->packed = temp1;//(ulong *) flint_stack_alloc_small(9);
	//F_mpn_clear(temp->packed, 9);
	temp->coeffs = temp2;//(F_mpz *) flint_stack_alloc_small(9);
   F_mpn_clear(temp->coeffs, 9);
	temp->small = 1;
	
	_F_mpz_mpoly_mul_small_3x3(temp, poly1, n1, poly2, n2 + 1);

	_F_mpz_mpoly_add_inplace(res, 1, temp);

	//ulong i;
for (i = 0; i < temp->length; i++)
	//	_F_mpz_demote(temp->coeffs + i);

	//flint_stack_release_small(); // coeffs
	//flint_stack_release_small(); // packed
}

void _F_mpz_mpoly_mul_small_mxn(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								         ulong n1, F_mpz_mpoly_t poly2, ulong n2)
{
	if (poly2->length - n2 == 1)
	{
		switch (poly1->length - n1)
		{
		   case 1: _F_mpz_mpoly_mul_small_1x1(res, poly1, n1, poly2, n2); break;
		   case 2: _F_mpz_mpoly_mul_small_2x1(res, poly1, n1, poly2, n2); break;
		   case 3: _F_mpz_mpoly_mul_small_3x1(res, poly1, n1, poly2, n2); break;
		   case 4: _F_mpz_mpoly_mul_small_4x1(res, poly1, n1, poly2, n2); break;
		}
	} else if (poly1->length - n1 == 1)
	{
		switch (poly2->length - n2)
		{
		   case 2: _F_mpz_mpoly_mul_small_2x1(res, poly2, n2, poly1, n1); break;
		   case 3: _F_mpz_mpoly_mul_small_3x1(res, poly2, n2, poly1, n1); break;
		   case 4: _F_mpz_mpoly_mul_small_4x1(res, poly2, n2, poly1, n1); break;
		}
	} else if (poly2->length - n2 == 2)
	{
		switch (poly1->length - n1)
		{
		   case 2: _F_mpz_mpoly_mul_small_2x2(res, poly1, n1, poly2, n2); break;
		   case 3: _F_mpz_mpoly_mul_small_3x2(res, poly1, n1, poly2, n2); break;
		   case 4: _F_mpz_mpoly_mul_small_4x2(res, poly1, n1, poly2, n2); break;
		}
	} else if (poly1->length - n1 == 2)
	{
		switch (poly2->length - n2)
		{
		   case 3: _F_mpz_mpoly_mul_small_3x2(res, poly2, n2, poly1, n1); break;
		   case 4: _F_mpz_mpoly_mul_small_4x2(res, poly2, n2, poly1, n1); break;
		}
	} else if (poly2->length - n2 == 3)
	{
		switch (poly1->length - n1)
		{
		   case 3: _F_mpz_mpoly_mul_small_3x3(res, poly1, n1, poly2, n2); break;
		   case 4: _F_mpz_mpoly_mul_small_4x3(res, poly1, n1, poly2, n2); break;
		}
	} else if (poly1->length - n1 == 3)
	{
		_F_mpz_mpoly_mul_small_4x3(res, poly2, n2, poly1, n1);
	} else _F_mpz_mpoly_mul_small_4x4(res, poly1, n1, poly2, n2);
}

void _F_mpz_mpoly_mul_small_mxn_2(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
	     ulong n1, F_mpz_mpoly_t poly2, ulong n2, ulong m, ulong n)
{
	if (n == 1)
	{
		switch (m)
		{
		   case 1: _F_mpz_mpoly_mul_small_1x1(res, poly1, n1, poly2, n2); break;
		   case 2: _F_mpz_mpoly_mul_small_2x1(res, poly1, n1, poly2, n2); break;
		   case 3: _F_mpz_mpoly_mul_small_3x1(res, poly1, n1, poly2, n2); break;
		   case 4: _F_mpz_mpoly_mul_small_4x1(res, poly1, n1, poly2, n2); break;
		}
	} else if (m == 1)
	{
		switch (n)
		{
		   case 2: _F_mpz_mpoly_mul_small_2x1(res, poly2, n2, poly1, n1); break;
		   case 3: _F_mpz_mpoly_mul_small_3x1(res, poly2, n2, poly1, n1); break;
		   case 4: _F_mpz_mpoly_mul_small_4x1(res, poly2, n2, poly1, n1); break;
		}
	} else if (n == 2)
	{
		switch (m)
		{
		   case 2: _F_mpz_mpoly_mul_small_2x2(res, poly1, n1, poly2, n2); break;
		   case 3: _F_mpz_mpoly_mul_small_3x2(res, poly1, n1, poly2, n2); break;
		   case 4: _F_mpz_mpoly_mul_small_4x2(res, poly1, n1, poly2, n2); break;
		}
	} else if (m == 2)
	{
		switch (n)
		{
		   case 3: _F_mpz_mpoly_mul_small_3x2(res, poly2, n2, poly1, n1); break;
		   case 4: _F_mpz_mpoly_mul_small_4x2(res, poly2, n2, poly1, n1); break;
		}
	} else if (n == 3)
	{
		switch (m)
		{
		   case 3: _F_mpz_mpoly_mul_small_3x3(res, poly1, n1, poly2, n2); break;
		   case 4: _F_mpz_mpoly_mul_small_4x3(res, poly1, n1, poly2, n2); break;
		}
	} else if (m == 3)
	{
		_F_mpz_mpoly_mul_small_4x3(res, poly2, n2, poly1, n1);
	} else _F_mpz_mpoly_mul_small_4x4(res, poly1, n1, poly2, n2);
}

void _F_mpz_mpoly_mul_small_Mxn(F_mpz_mpoly_t res, F_mpz_mpoly_t poly1, 
								           F_mpz_mpoly_t poly2, ulong n2, ulong n)
{
	ulong ptr = 0;
   ulong i, exp;

	F_mpz_mpoly_t temp;
	temp->coeffs = (F_mpz *) flint_stack_alloc(4*n);
   F_mpn_clear(temp->coeffs, 4*n);
	temp->packed = (ulong *) flint_stack_alloc(4*n);
	//F_mpn_clear(temp->packed, 4*n);
	temp->small = 1;
	ulong * poly1_e = poly1->packed;
	ulong * poly2_e = poly2->packed + n2;
	ulong * res_e = res->packed;

	if (poly1->length < 4) 
	{
		_F_mpz_mpoly_mul_small_mxn(res, poly1, 0, poly2, n2);
		return;
	}
	
	if (n == 4)
	{
		_F_mpz_mpoly_mul_small_4x4(res, poly1, 0, poly2, n2);
		for (i = 4; i + 4 <= poly1->length; i+=4)
	   {
			//printf("here-b %ld\n", i);
		   
			_F_mpz_mpoly_mul_small_4x4(temp, poly1, i, poly2, n2);
				
         //printf("here-c %ld\n", i);
		   exp = poly1_e[i] + poly2_e[0];
			while ((res_e[ptr] < exp) && (ptr < res->length)) ptr++;
         //printf("here-c2 %ld, %ld\n", i, temp->length);

			/*char * syms[4];
	      syms[0] = "t";
	      syms[1] = "x";
	      syms[2] = "y";
         syms[3] = "z";
         temp->packed_bits = 8;
         temp->vars = 4;
         temp->ordering = GRLEX;
      	temp->small = 1;
	      res->packed_bits = 8;
         res->vars = 4;
         res->ordering = GRLEX;
      	res->small = 1;
	      printf("********-->length temp = %ld\n", temp->length);
	      printf("********-->length res = %ld\n", res->length);
	      F_mpz_mpoly_print_pretty(temp, syms);printf("\n");
	      F_mpz_mpoly_print_pretty(res, syms);printf("\n");*/
	

		   _F_mpz_mpoly_add_inplace(res, ptr, temp);
			//F_mpz_mpoly_print_pretty(res, syms);printf("\n");
	
			//printf("here-d\n");
	   }
	} else if (n == 3)
	{
		_F_mpz_mpoly_mul_small_4x3(res, poly1, 0, poly2, n2);
		for (i = 4; i + 4 <= poly1->length; i+=4)
	   {
			_F_mpz_mpoly_mul_small_4x3(temp, poly1, i, poly2, n2);
         exp = poly1_e[i] + poly2_e[0];
			while ((res_e[ptr] < exp) && (ptr < res->length)) ptr++;
         _F_mpz_mpoly_add_inplace(res, ptr, temp);
	   }
	} else if (n == 2)
	{
		_F_mpz_mpoly_mul_small_4x2(res, poly1, 0, poly2, n2);
		for (i = 4; i + 4 <= poly1->length; i+=4)
	   {
			_F_mpz_mpoly_mul_small_4x2(temp, poly1, i, poly2, n2);
         exp = poly1_e[i] + poly2_e[0];
			while ((res_e[ptr] < exp) && (ptr < res->length)) ptr++;
         _F_mpz_mpoly_add_inplace(res, ptr, temp);
	   }
	} else // n == 1
	{
		_F_mpz_mpoly_mul_small_4x1(res, poly1, 0, poly2, n2);
		for (i = 4; i + 4 <= poly1->length; i+=4)
	   {
			ptr+=4;
			_F_mpz_mpoly_mul_small_4x1(temp, poly1, i, poly2, n2);
         _F_mpz_mpoly_add_inplace(res, ptr, temp);  	
	   }
	}

	if (poly1->length > i)
	{
		F_mpz_mpoly_t temp2;
		temp2->coeffs = (ulong *) flint_stack_alloc(n);
		F_mpn_clear(temp2->coeffs, n);
	   temp2->packed = (F_mpz *) flint_stack_alloc(n);
		//F_mpn_clear(temp2->packed, n);
	   temp2->small = 1;
		temp2->length = n;
      
		ulong j;
		for (j = 0; j < n; j++)
		{
			temp2->coeffs[j] = poly2->coeffs[n2 + j];
			temp2->packed[j] = poly2->packed[n2 + j];
		}

		_F_mpz_mpoly_mul_small_mxn(temp, poly1, i, temp2, 0);
		
		/*char * syms[4];
	   syms[0] = "t";
	   syms[1] = "x";
	   syms[2] = "y";
      syms[3] = "z";
		temp->packed_bits = 8;
		temp->vars = 4;
		temp->ordering = GRLEX;
		temp->small = 1;
	   temp2->packed_bits = 8;
		temp2->vars = 4;
		temp2->ordering = GRLEX;
		temp2->small = 1;
	   printf("-------->length temp = %ld\n", temp->length);
	   printf("-------->length temp2 = %ld\n", temp2->length);
	   printf("poly1->length = %ld, i = %ld\n", poly1->length, i);
		F_mpz_print(poly1->coeffs + poly1->length - 1); printf("\n");
		F_mpz_mpoly_print_pretty(temp, syms);printf("\n");
      F_mpz_mpoly_print_pretty(temp2, syms);printf("\n");*/
   
		_F_mpz_mpoly_add_inplace(res, ptr, temp);

		//ulong i;
for (i = 0; i < temp2->length; i++)
   		//_F_mpz_demote(temp2->coeffs + i);

		flint_stack_release(); // packed
		flint_stack_release(); // coeffs
	}

	//ulong i;
for (i = 0; i < temp->length; i++)
	//	_F_mpz_demote(temp->coeffs + i);

	flint_stack_release(); // packed
	flint_stack_release(); // coeffs
}

void F_mpz_mpoly_mul_small_recursive(F_mpz_mpoly_t res, 
				  F_mpz_mpoly_t poly1, F_mpz_mpoly_t poly2, ulong n2, ulong n)
{
	if (n <= 4) 
	{
		//printf("here\n");
		F_mpz_mpoly_fit_length(res, poly1->length * n);
		//printf("here2 %ld, %ld\n", n2, n);
		_F_mpz_mpoly_mul_small_Mxn(res, poly1, poly2, n2, n);
		//printf("here3\n");
		return;
	}

	ulong m1 = n/2;
	m1 = ((((m1 - 1)>>2)+1)<<2);
	ulong m2 = n - m1;

	F_mpz_mpoly_t temp;
	F_mpz_mpoly_init2(temp, 0, poly1->vars, poly1->ordering);

	//printf("1: %ld, %ld, %ld, %ld\n", poly1->length, poly2->length, n2, m1);
	F_mpz_mpoly_mul_small_recursive(res, poly1, poly2, n2, m1);
	//printf("2: %ld, %ld\n", n2 + m1, m2);
	F_mpz_mpoly_mul_small_recursive(temp, poly1, poly2, n2 + m1, m2);
   //printf("done1\n");

	// 70 -> 36 + 34 -> 20 + 16 + 20 + 14 -> 12 + 8 + 8 + 8 + 12 + 8 + 12 + 2
	/*if (m1 == 4)
	{
		printf("m1 = %ld, m2 = %ld\n", m1, m2);
		char * syms[4];
	   syms[0] = "t";
	   syms[1] = "x";
	   syms[2] = "y";
      syms[3] = "z";
	   F_mpz_mpoly_print_pretty(res, syms);printf("\n");
		F_mpz_mpoly_print_pretty(temp, syms);printf("\n");

	}*/

	F_mpz_mpoly_add_inplace(res, 0, temp);
      
	/*if (m1 == 4) 
	{
		char * syms[4];
	   syms[0] = "t";
	   syms[1] = "x";
	   syms[2] = "y";
      syms[3] = "z";
	   F_mpz_mpoly_print_pretty(res, syms); printf("\n");
	}*/

   //printf("done2\n");
	F_mpz_mpoly_clear(temp);
}

ulong F_mpz_mpoly_reheapify(ulong * heap, F_mpz_mpoly_heap_s * entries, ulong heap_bottom)
{
   ulong entry = 1;

	while (2*entry < heap_bottom)
	{
      if (entries[heap[2*entry]].packed <= entries[heap[2*entry + 1]].packed)
		{
         heap[entry] = heap[2*entry];
			entry = 2*entry;
		} else
		{
         heap[entry] = heap[2*entry + 1];
			entry = 2*entry + 1;
		}
	}

	if (2*entry == heap_bottom)
	{
		heap[entry] = heap[2*entry];
		entry = 2*entry;
	}

	return entry;
}

int insert_head(ulong * heap, F_mpz_mpoly_heap_s * entries, ulong i, ulong packed)
{
	ulong head = heap[1];
	if (packed == entries[head].packed)
	{
      entries[i].chain = head;
		heap[1] = i;
		return 1;
	}

	return 0;
}

ulong F_mpz_mpoly_heap_filter(ulong * heap, F_mpz_mpoly_heap_s * entries, ulong i, ulong packed, ulong heap_bottom)
{
	ulong entry = 1;

	if (heap_bottom == 0)
	{
		entries[i].chain = -1L;
		return i;
	}

	if (packed == entries[heap[1]].packed)
	{
		entries[i].chain = heap[1];
		heap[1] = i;
		return -1L;
	}

	if (packed > entries[heap[1]].packed)
	{

	while (2*entry < heap_bottom)
	{
      if (entries[heap[2*entry]].packed >= entries[heap[2*entry + 1]].packed) 
		{
         ulong packed2 = entries[heap[2*entry]].packed;
			if (packed <= packed2) 
			{
				if (packed < packed2)
				{
					entry = 2*entry;
				   goto place;
				} else
				{
               entries[i].chain = heap[entry*2];
			      heap[entry*2] = i;
			      return -1L;
			   }
			}
			
			entry = 2*entry;
		} else
		{
         ulong packed2 = entries[heap[2*entry + 1]].packed;
			if (packed <= packed2) 
			{
				if (packed < packed2)
				{
					entry = 2*entry + 1;
				   goto place;
				} else
				{
					entries[i].chain = heap[entry*2 + 1];
			      heap[entry*2 + 1] = i;
			      return -1L;
				}
			}
			entry = 2*entry + 1;
		}
	}

	if (2*entry == heap_bottom)
	{
		if (packed < entries[heap[2*entry]].packed) 
		{
			entry = 2*entry;
		} else if (packed == entries[heap[2*entry]].packed) 
		{
         entries[i].chain = heap[entry*2];
			heap[entry*2] = i;
			return -1L;
		}
	}

	}
   
place:

	entries[i].chain = -1L;

	ulong temp;
	ulong temp_i = i;

	while (entry <= heap_bottom)
	{
      temp = heap[entry];
		heap[entry] = temp_i;
		temp_i = temp;
      
		/*if (2*entry > heap_bottom) return temp_i;
		if (2*entry == heap_bottom) 
		{
			temp_i = heap[2*entry];
			heap[2*entry] = temp;
			return temp_i;
		}*/
		entry = 2*entry;
	}

	return temp_i;
}

ulong F_mpz_mpoly_heap_insert(ulong * heap, F_mpz_mpoly_heap_s * entries, ulong empty, ulong i, ulong packed, ulong heap_bottom)
{
	ulong entry = empty;
		
	while (entry > 1)
	{
		ulong heap2 = heap[entry/2];
		if (entries[heap2].packed > packed)
		{
			heap[entry] = heap2;
			entry /= 2;
		} else if (entries[heap2].packed == packed) 
		{
			entries[i].chain = heap2;
			heap[entry/2] = i;
			return entry;
		} else break;
	}
   
	while (2*entry < heap_bottom)
	{
      if (entries[heap[2*entry]].packed <= entries[heap[2*entry + 1]].packed) 
		{
         if (packed < entries[heap[2*entry]].packed) break;
			else if (packed == entries[heap[2*entry]].packed) 
			{
            entries[i].chain = heap[entry*2];
			   heap[entry*2] = i;
			   return entry;
			}
			heap[entry] = heap[2*entry];
			entry = 2*entry;
		} else
		{
         if (packed < entries[heap[2*entry + 1]].packed) break;
			else if (packed == entries[heap[2*entry + 1]].packed) 
			{
            entries[i].chain = heap[entry*2 + 1];
			   heap[entry*2 + 1] = i;
			   return entry;
			}
			heap[entry] = heap[2*entry + 1];
			entry = 2*entry + 1;
		}
	}
   
	if (2*entry == heap_bottom)
	{
		if (packed > entries[heap[2*entry]].packed) 
		{
			heap[entry] = heap[2*entry];
		   entry = 2*entry;
		} else if (packed == entries[heap[2*entry]].packed) 
		{
         entries[i].chain = heap[entry*2];
			heap[entry*2] = i;
			return entry;
		}
	}
   
	entries[i].chain = -1L;
	heap[entry] = i;
	return -1L;
}

ulong F_mpz_mpoly_heap_insert3(ulong * heap, F_mpz_mpoly_heap_s * entries, ulong empty, ulong i, ulong packed, ulong heap_bottom)
{
	ulong entry = 1;
   
	while (2*entry < heap_bottom)
	{
      if (entries[heap[2*entry]].packed <= entries[heap[2*entry + 1]].packed) 
		{
         if (packed < entries[heap[2*entry]].packed) break;
			else if (packed == entries[heap[2*entry]].packed) 
			{
            entries[i].chain = heap[entry*2];
			   heap[entry*2] = i;
			   return entry;
			}
			heap[entry] = heap[2*entry];
			entry = 2*entry;
		} else
		{
         if (packed < entries[heap[2*entry + 1]].packed) break;
			else if (packed == entries[heap[2*entry + 1]].packed) 
			{
            entries[i].chain = heap[entry*2 + 1];
			   heap[entry*2 + 1] = i;
			   return entry;
			}
			heap[entry] = heap[2*entry + 1];
			entry = 2*entry + 1;
		}
	}
   
	if (2*entry == heap_bottom)
	{
		if (packed > entries[heap[2*entry]].packed) 
		{
			heap[entry] = heap[2*entry];
		   entry = 2*entry;
		} else if (packed == entries[heap[2*entry]].packed) 
		{
         entries[i].chain = heap[entry*2];
			heap[entry*2] = i;
			return entry;
		}
	}
   
	entries[i].chain = -1L;
	heap[entry] = i;
	return -1L;
}

ulong F_mpz_mpoly_heap_insert2(ulong * heap, F_mpz_mpoly_heap_s * entries, ulong empty, ulong i, ulong packed, ulong heap_bottom)
{
	ulong entry = empty;
	ulong next;

	while (entry > 1)
	{
		ulong heap2 = heap[entry/2];
		if (entries[heap2].packed > packed)
		{
			heap[entry] = heap2;
			entry /= 2;
		} else if (entries[heap2].packed == packed) 
		{
			
			do {
				next = entries[i].chain;
				entries[i].chain = heap[entry/2];
			   heap[entry/2] = i;
				i = next;
			} while (i != -1L);

			return entry;
		} else break;
	}
   
	while (2*entry < heap_bottom)
	{
      if (entries[heap[2*entry]].packed <= entries[heap[2*entry + 1]].packed) 
		{
         if (packed <= entries[heap[2*entry]].packed) break;
			/*else if (packed == entries[heap[2*entry]].packed) 
			{
			   do {
				   next = entries[i].chain;
				   entries[i].chain = heap[entry*2];
			      heap[entry*2] = i;
				   i = next;
			   } while (i != -1L);

				return entry;
			}*/
			heap[entry] = heap[2*entry];
			entry = 2*entry;
		} else
		{
         if (packed <= entries[heap[2*entry + 1]].packed) break;
			/*else if (packed == entries[heap[2*entry + 1]].packed) 
			{
            do {
				   next = entries[i].chain;
				   entries[i].chain = heap[entry*2 + 1];
			      heap[entry*2 + 1] = i;
				   i = next;
			   } while (i != -1L);
				
				return entry;
			}*/
			heap[entry] = heap[2*entry + 1];
			entry = 2*entry + 1;
		}
	}
   
	if (2*entry == heap_bottom)
	{
		if (packed > entries[heap[2*entry]].packed) 
		{
			heap[entry] = heap[2*entry];
		   entry = 2*entry;
		} /*else if (packed == entries[heap[2*entry]].packed) 
		{
         do {
				next = entries[i].chain;
				entries[i].chain = heap[entry*2];
			   heap[entry*2] = i;
				i = next;
			} while (i != -1L);
			
			return entry;
		}*/
	}
   
	heap[entry] = i;
	return -1L;
}
	
void F_mpz_mpoly_mul_small_heap(F_mpz_mpoly_t res, 
				  F_mpz_mpoly_t poly1, F_mpz_mpoly_t poly2)
{
   ulong short_length = poly1->length;
	ulong short_length2 = poly2->length;
	ulong log_length = 1;
	while ((1L<<log_length) <= short_length) log_length++;
	ulong r = 0;
	
	ulong * saved = flint_heap_alloc(short_length);
	ulong saved_num;

	F_mpz_mpoly_heap_s * entries = (F_mpz_mpoly_heap_s *) flint_heap_alloc_bytes(short_length*sizeof(F_mpz_mpoly_heap_s));
		 
	ulong * heap = (ulong *) flint_heap_alloc(1L<<log_length);
   ulong heap_bottom = 1;//short_length;

   ulong i;
   for (i = 0; i < short_length; i++)
	{
		entries[i].j = 0;
		entries[i].chain = -1L;
   }

	entries[0].packed = poly1->packed[0] + poly2->packed[0];
	heap[1] = 0;	
	
	F_mpz_mpoly_t temp;
	F_mpz_mpoly_init2(temp, 16, poly1->vars, poly1->ordering);

	F_mpz_t temp_f;
	F_mpz_init(temp_f);

	ulong sum[3];

   ulong i, j, oldi, next;
	ulong packed, old_packed, empty;
	ulong arr[2];
	ulong cry;

	while (heap_bottom > 0)
	{
      saved_num = 0;
		
		i = heap[1];
		j = entries[i].j;
		packed = poly1->packed[i] + poly2->packed[j];
		
		do
		{
		   old_packed = packed;
			oldi = i;

			//_F_mpz_mpoly_mul_small_1(sum, poly1, i, poly2, j);
		   ulong c1 = poly1->coeffs[i];
         ulong c2 = poly2->coeffs[j];
   
	      umul_ppmm(sum[1], sum[0], c1, c2);
    	   i = entries[i].chain;
			
		   while (i != -1L)
		   {
		      j = entries[i].j;
   
	         ulong c1 = poly1->coeffs[i];
            ulong c2 = poly2->coeffs[j];
   
	         umul_ppmm(arr[1], arr[0], c1, c2);
	         add_ssaaaa(sum[1], sum[0], sum[1], sum[0], arr[1], arr[0]);
		      i = entries[i].chain;
		   } // end while (i != -1L)
      
		   if (!res->length)
		   {
			   F_mpz_mpoly_fit_length(res, r+1);
		      res->length++;
		   } else if (res->packed[r] < packed) 
		   {
			   r++;
			   res->length++;
			   F_mpz_mpoly_fit_length(res, r+1);
		      F_mpz_zero(res->coeffs + r);
		   }
      
		   if (sum[1])
			{
				F_mpz_set_ui(temp_f, sum[1]);
		      F_mpz_mul_2exp(temp_f, temp_f, FLINT_BITS);
		      F_mpz_add_ui(temp_f, temp_f, sum[0]);
			} else
            F_mpz_set_ui(temp_f, sum[0]);
			F_mpz_add(res->coeffs + r, res->coeffs + r, temp_f);
		   res->packed[r] = packed;

			ulong packed2 = (heap_bottom >= 2 ? entries[heap[2]].packed : -1L);
         ulong packed3 = (heap_bottom >= 3 ? entries[heap[3]].packed : -1L);

			if ((old_packed == packed2) || (old_packed == packed3))
			{
				empty = F_mpz_mpoly_reheapify(heap, entries, heap_bottom);

		      while (empty != -1L)
	         {	
		         if (empty != heap_bottom) 
			         empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			      else empty = -1L;
			      heap_bottom--;
	         } 

			   if (heap_bottom != 0)
			   {
				   i = heap[1];
		         j = entries[i].j;
		         packed = poly1->packed[i] + poly2->packed[j];
				}
			} else packed = -1L;
			
			saved[saved_num] = oldi;
			saved_num++;			

		} while (old_packed == packed);
 

		long k;
		for (k = saved_num - 1; k >= 0L; k--)
		{
		
		i = saved[k];

		if (k == saved_num - 1)
		{
		
		   next = entries[i].chain;
			
		   entries[i].j++;
		
		   if (entries[i].j == 1) // not exhausted this line yet
		   {
			   packed = poly1->packed[i] + poly2->packed[entries[i].j];
            entries[i].packed = packed;
 
			   empty = F_mpz_mpoly_heap_insert(heap, entries, 1, i, packed, heap_bottom);
			
			   while (empty != -1L)
	         {	
		         if (empty != heap_bottom) 
			         empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			      else empty = -1L;
			      heap_bottom--;
	         }

				if (i < short_length - 1)
				{
					packed = poly1->packed[i+1] + poly2->packed[0];
               entries[i+1].packed = packed;
 
				   ulong dislodged = F_mpz_mpoly_heap_filter(heap, entries, i+1, packed, heap_bottom);
			   
				   if (dislodged != -1L)
				   {
					   heap_bottom++;
				      empty = heap_bottom;

				      empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, dislodged, entries[dislodged].packed, heap_bottom);

		            while (empty != -1L)
	               {	
		               if (empty != heap_bottom) 
			               empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			            else empty = -1L;
			            heap_bottom--;
	               }
					}
				}
		   } else if (entries[i].j < short_length2) // not exhausted this line yet
		   {
			   packed = poly1->packed[i] + poly2->packed[entries[i].j];
            entries[i].packed = packed;
 
			   empty = F_mpz_mpoly_heap_insert3(heap, entries, 1, i, packed, heap_bottom);
			
			   while (empty != -1L)
	         {	
		         if (empty != heap_bottom) 
			         empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			      else empty = -1L;
			      heap_bottom--;
	         }
		   } else
		   {
            empty = F_mpz_mpoly_reheapify(heap, entries, heap_bottom);

		      while (empty != -1L)
	         {	
		         if (empty != heap_bottom) 
			         empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			      else empty = -1L;
			      heap_bottom--;
	         }
		   }

		   i = next;
		}

		while (i != -1L)
		{	
			next = entries[i].chain;
			
			entries[i].j++;
			if (entries[i].j == 1) // not exhausted this line yet
		   {
			   packed = poly1->packed[i] + poly2->packed[entries[i].j];
            entries[i].packed = packed;

				ulong dislodged = F_mpz_mpoly_heap_filter(heap, entries, i, packed, heap_bottom);
			   
				if (dislodged != -1L)
				{
					heap_bottom++;
				   empty = heap_bottom;

				   empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, dislodged, entries[dislodged].packed, heap_bottom);

		         while (empty != -1L)
	            {	
		            if (empty != heap_bottom) 
			            empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			         else empty = -1L;
			         heap_bottom--;
	            }
				}

				if (i < short_length - 1)
				{
					packed = poly1->packed[i+1] + poly2->packed[0];
               entries[i+1].packed = packed;
 
				   dislodged = F_mpz_mpoly_heap_filter(heap, entries, i+1, packed, heap_bottom);
			   
				   if (dislodged != -1L)
				   {
					   heap_bottom++;
				      empty = heap_bottom;

				      empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, dislodged, entries[dislodged].packed, heap_bottom);

		            while (empty != -1L)
	               {	
		               if (empty != heap_bottom) 
			               empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			            else empty = -1L;
			            heap_bottom--;
	               }
					}
				}
			} else if (entries[i].j < short_length2) // not exhausted this line yet
		   {
			   packed = poly1->packed[i] + poly2->packed[entries[i].j];
            entries[i].packed = packed;

				ulong dislodged = F_mpz_mpoly_heap_filter(heap, entries, i, packed, heap_bottom);
			   
				if (dislodged != -1L)
				{
					heap_bottom++;
				   empty = heap_bottom;

				   empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, dislodged, entries[dislodged].packed, heap_bottom);

		         while (empty != -1L)
	            {	
		            if (empty != heap_bottom) 
			            empty = F_mpz_mpoly_heap_insert2(heap, entries, empty, heap[heap_bottom], entries[heap[heap_bottom]].packed, heap_bottom - 1);
			         else empty = -1L;
			         heap_bottom--;
	            }
				}
			}


			i = next;
		}

		}

	}
	
   F_mpz_clear(temp_f);
	F_mpz_mpoly_clear(temp);
	flint_heap_free(saved);
	flint_heap_free(entries);
	flint_heap_free(heap);
}
