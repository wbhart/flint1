/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   fmpz_poly.c: Polynomials over "flat" multiprecision integers

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "fmpz_poly.h"
#include "fmpz.h"
#include "mpir.h"
#include "memory_manager.h"

/* ==============================================================================

   Memory management

===============================================================================*/

void fmpz_poly_init(fmpz_poly_t poly)
{
   poly->coeffs = NULL;
   
   poly->alloc = 0;
   poly->length = 0;
}

void fmpz_poly_init2(fmpz_poly_t poly, unsigned long alloc)
{
   if ((long) alloc > 0)
   {
      poly->coeffs = fmpz_init_array(alloc);
   }
   else poly->coeffs = NULL;
   
   poly->alloc = alloc;
   poly->length = 0;
}

/* 
   Shrink or expand a polynomial to "alloc" coefficients 
*/

void fmpz_poly_realloc(fmpz_poly_t poly, ulong alloc)
{
   if ((long) alloc > 0)
   {
      if (poly->alloc) poly->coeffs = fmpz_realloc_array(poly->coeffs, poly->alloc, alloc);
      else poly->coeffs = fmpz_init_array(alloc);
      poly->alloc = alloc;
   } else
   {
      if (poly->coeffs) fmpz_clear_array(poly->coeffs, poly->alloc);
      fmpz_poly_init(poly);
   }   
   
   // truncate actual data if necessary
   if (poly->length > alloc)
   {
      poly->length = alloc;
      _fmpz_poly_normalise(poly);
   }     
}

void fmpz_poly_fit_length(fmpz_poly_t poly, ulong alloc)
{
   if (alloc <= poly->alloc) return;

   if (alloc < 2*poly->alloc) alloc = 2*poly->alloc;
   
   fmpz_poly_realloc(poly, alloc);
}

void fmpz_poly_clear(fmpz_poly_t poly)
{
   if (poly->coeffs) fmpz_clear_array(poly->coeffs, poly->alloc);
}

/* ==============================================================================

   Conversion

===============================================================================*/

void mpz_poly_to_fmpz_poly(fmpz_poly_t res, mpz_poly_t poly)
{
   fmpz_poly_fit_length(res, poly->length);
   
   res->length = poly->length;
   if (poly->length == 0)
      return;

   for (ulong i = 0; i < poly->length; i++)
      mpz_to_fmpz(res->coeffs + i, poly->coeffs[i]);
}


void fmpz_poly_to_mpz_poly(mpz_poly_t res, fmpz_poly_t poly)
{
   mpz_poly_ensure_alloc(res, poly->length);

   res->length = poly->length;
   
   if (poly->length == 0) return;
   
   for (ulong i = 0; i < poly->length; i++)
   {
      fmpz_to_mpz(res->coeffs[i], poly->coeffs + i);
   }
}

/* ==============================================================================

   Addition/subtraction

===============================================================================*/

void fmpz_poly_add(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
{
   ulong shorter, longer;
   fmpz_t * coeffs1, * coeffs2, * coeffs_out;   
   
   if (input1->length > input2->length)
   {
      shorter = input2->length;
      longer = input1->length;
   } else
   {
      shorter = input1->length;
      longer = input2->length;
   }
   
   fmpz_poly_fit_length(output, longer);

   coeffs_out = output->coeffs;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   
   for (ulong i = 0; i < shorter; i++)
   {
      fmpz_add(coeffs_out+i, coeffs1+i, coeffs2+i);
   }    
   
   if (input1 != output)
   {
      for (ulong i = shorter; i < input1->length; i++)
      {
          fmpz_set(coeffs_out+i, coeffs1+i);
      }
   }

   if (input2 != output)
   {
      for (ulong i = shorter; i < input2->length; i++)
      {
         fmpz_set(coeffs_out+i, coeffs2+i);
      }
   }
   
   output->length = longer;
   if (input1->length == input2->length) _fmpz_poly_normalise(output);
}

void fmpz_poly_sub(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
{
   ulong shorter, longer;
   fmpz_t * coeffs1, * coeffs2, * coeffs_out;   
   
   if (input1->length > input2->length)
   {
      shorter = input2->length;
      longer = input1->length;
   } else
   {
      shorter = input1->length;
      longer = input2->length;
   }
   
   fmpz_poly_fit_length(output, longer);

   coeffs_out = output->coeffs;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   
   for (ulong i = 0; i < shorter; i++)
   {
      fmpz_sub(coeffs_out+i, coeffs1+i, coeffs2+i);
   }    
   
   if (input1 != output)
   {
      for (ulong i = shorter; i < input1->length; i++)
      {
          fmpz_set(coeffs_out+i, coeffs1+i);
      }
   }

   for (ulong i = shorter; i < input2->length; i++)
   {
      fmpz_neg(coeffs_out+i, coeffs2+i);
   }
   
   output->length = longer;
   if (input1->length == input2->length) _fmpz_poly_normalise(output);
}

// *************** end of file
