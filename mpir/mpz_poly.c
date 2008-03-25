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
/******************************************************************************

mpz_poly.c: Polynomials over Z, implemented as an array of mpz_t's

Copyright (C) 2007, William Hart and David Harvey

******************************************************************************/

#include <string.h>
#include <stdio.h>
#include "mpz_poly.h"
#include "memory_manager.h"
#include "mpir.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpz_poly_init(mpz_poly_t poly)
{
   poly->coeffs = (mpz_t*) mpir_alloc(sizeof(mpz_t));
   mpz_init(poly->coeffs[0]);
   poly->alloc = 1;
   poly->length = 0;
}


void mpz_poly_init2(mpz_poly_t poly, unsigned long alloc)
{
   if ((long) alloc <= 0) mpz_poly_init(poly);
   
   poly->coeffs = (mpz_t*) mpir_alloc(alloc * sizeof(mpz_t));
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init(poly->coeffs[i]);
   
   poly->alloc = alloc;
   poly->length = 0;
}


void mpz_poly_init3(mpz_poly_t poly, unsigned long alloc, unsigned long bits)
{
   if ((long) alloc <= 0) mpz_poly_init(poly);
   
   poly->coeffs = (mpz_t*) mpir_alloc(alloc * sizeof(mpz_t));
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init2(poly->coeffs[i], bits);
   
   poly->alloc = alloc;
   poly->length = 0;
}


void mpz_poly_clear(mpz_poly_t poly)
{
   for (long i = 0; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);

   mpir_free(poly->coeffs);
}


void mpz_poly_realloc(mpz_poly_t poly, unsigned long alloc)
{
   if ((long) alloc <= 0) alloc = 1;
   
   // clear any mpz_t's beyond the new array length
   for (unsigned long i = alloc; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);

   poly->coeffs = (mpz_t*) mpir_realloc(poly->coeffs,
                                              alloc * sizeof(mpz_t));
   
   // init any new mpz_t's required
   for (unsigned long i = poly->alloc; i < alloc; i++)
      mpz_init(poly->coeffs[i]);

   poly->alloc = alloc;
   
   // truncate poly if necessary
   if (poly->length > alloc)
   {
      poly->length = alloc;
      mpz_poly_normalise(poly);
   }
}


void mpz_poly_realloc2(mpz_poly_t poly, unsigned long alloc,
                       unsigned long bits)
{
   if ((long) alloc <= 0) alloc = 1;
   
   // clear any mpz_t's beyond the new array length
   for (unsigned long i = alloc; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);

   poly->coeffs = (mpz_t*) mpir_realloc(poly->coeffs,
                                              alloc * sizeof(mpz_t));
   
   // init any new mpz_t's required
   for (unsigned long i = poly->alloc; i < alloc; i++)
      mpz_init2(poly->coeffs[i], bits);

   poly->alloc = alloc;
   
   // truncate poly if necessary
   if (poly->length > alloc)
   {
      poly->length = alloc;
      mpz_poly_normalise(poly);
   }
}



void __mpz_poly_ensure_alloc(mpz_poly_t poly, unsigned long alloc)
{
   if (alloc < 2*poly->alloc)
      alloc = 2*poly->alloc;
   mpz_poly_realloc(poly, alloc);
}


/****************************************************************************

   Setting/retrieving coefficients

****************************************************************************/


void mpz_poly_set_coeff(mpz_poly_t poly, unsigned long n, mpz_t c)
{
   mpz_poly_ensure_alloc(poly, n+1);

   if (n+1 < poly->length)
      // set interior coefficient
      mpz_set(poly->coeffs[n], c);

   else if (n+1 == poly->length)
   {
      // set leading coefficient
      if (mpz_sgn(c))
         mpz_set(poly->coeffs[n], c);
      else
      {
         // set leading coefficient to zero
         poly->length--;
         mpz_poly_normalise(poly);
      }
   }
   
   else
   {
      // extend polynomial
      if (!mpz_sgn(c))
         return;

      for (ulong i = poly->length; i < n; i++)
         mpz_set_ui(poly->coeffs[i], 0);
         
      mpz_set(poly->coeffs[n], c);
      poly->length = n+1;
   }
}


void mpz_poly_set_coeff_ui(mpz_poly_t poly, ulong n, ulong c)
{
   mpz_poly_ensure_alloc(poly, n+1);
   
   if (n+1 < poly->length)
      // set interior coefficient
      mpz_set_ui(poly->coeffs[n], c);

   else if (n+1 == poly->length)
   {
      // set leading coefficient
      if (c)
         mpz_set_ui(poly->coeffs[n], c);
      else
      {
         // set leading coefficient to zero
         poly->length--;
         mpz_poly_normalise(poly);
      }
   }
   
   else
   {
      // extend polynomial
      if (!c)
         return;
      
      for (unsigned long i = poly->length; i < n; i++)
         mpz_set_ui(poly->coeffs[i], 0L);
         
      mpz_set_ui(poly->coeffs[n], c);
      poly->length = n+1;
   }
}


void mpz_poly_set_coeff_si(mpz_poly_t poly, unsigned long n, long c)
{
   mpz_poly_ensure_alloc(poly, n+1);
   
   if (n+1 < poly->length)
      // set interior coefficient
      mpz_set_si(poly->coeffs[n], c);

   else if (n+1 == poly->length)
   {
      // set leading coefficient
      if (c)
         mpz_set_si(poly->coeffs[n], c);
      else
      {
         // set leading coefficient to zero
         poly->length--;
         mpz_poly_normalise(poly);
      }
   }
   
   else
   {
      // extend polynomial
      if (!c)
         return;
      
      for (unsigned long i = poly->length; i < n; i++)
         mpz_set_ui(poly->coeffs[i], 0L);
         
      mpz_set_si(poly->coeffs[n], c);
      poly->length = n+1;
   }
}

/****************************************************************************

   Length and degree

****************************************************************************/


void mpz_poly_normalise(mpz_poly_t poly)
{
   while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]))
      poly->length--;
}

/****************************************************************************

   Comparison

****************************************************************************/


int mpz_poly_equal(mpz_poly_t poly1, mpz_poly_t poly2)
{
   if (poly1->length != poly2->length)
      return 0;

   for (long i = 0; i < poly1->length; i++)
      if (mpz_cmp(poly1->coeffs[i], poly2->coeffs[i]))
         return 0;

   return 1;
}

/****************************************************************************

   Addition/subtraction

****************************************************************************/


void mpz_poly_add(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   // rearrange parameters to make poly1 no longer than poly2
   if (poly1->length > poly2->length)
      SWAP_MPZ_POLY_PTRS(poly1, poly2);
      
   mpz_poly_ensure_alloc(res, poly2->length);

   unsigned long i;
   
   for (i = 0; i < poly1->length; i++)
      mpz_add(res->coeffs[i], poly1->coeffs[i], poly2->coeffs[i]);

   for (; i < poly2->length; i++)
      mpz_set(res->coeffs[i], poly2->coeffs[i]);

   res->length = poly2->length;
   mpz_poly_normalise(res);
}



void mpz_poly_sub(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   if (poly1 == poly2)
   {
      // equal operands
      res->length = 0;
      return;
   }

   // rearrange parameters to make poly1 no longer than poly2
   int swapped = 0;
   if (poly1->length > poly2->length)
   {
      swapped = 1;
      SWAP_MPZ_POLY_PTRS(poly1, poly2);
   }
      
   mpz_poly_ensure_alloc(res, poly2->length);

   unsigned long i;
   
   if (swapped)
   {
      for (i = 0; i < poly1->length; i++)
         mpz_sub(res->coeffs[i], poly2->coeffs[i], poly1->coeffs[i]);
      for (; i < poly2->length; i++)
         mpz_set(res->coeffs[i], poly2->coeffs[i]);
   }
   else
   {
      for (i = 0; i < poly1->length; i++)
         mpz_sub(res->coeffs[i], poly1->coeffs[i], poly2->coeffs[i]);
      for (; i < poly2->length; i++)
         mpz_neg(res->coeffs[i], poly2->coeffs[i]);
   }

   res->length = poly2->length;
   mpz_poly_normalise(res);
}


void mpz_poly_neg(mpz_poly_t res, mpz_poly_t poly)
{
   mpz_poly_ensure_alloc(res, poly->length);

   for (unsigned long i = 0; i < poly->length; i++)
      mpz_neg(res->coeffs[i], poly->coeffs[i]);
   
   res->length = poly->length;
}

// ********************* end of file 
