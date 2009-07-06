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
#include "flint.h"
#include "mpz_poly.h"
#include "mpz_poly-tuning.h"
#include "fmpz.h"
#include "fmpz_poly.h"


/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpz_poly_init(mpz_poly_t poly)
{
   poly->coeffs = (mpz_t*) flint_heap_alloc_bytes(sizeof(mpz_t));
   mpz_init(poly->coeffs[0]);
   poly->alloc = 1;
   poly->length = 0;
}


void mpz_poly_init2(mpz_poly_t poly, unsigned long alloc)
{
   if ((long) alloc <= 0) mpz_poly_init(poly);
   
   poly->coeffs = (mpz_t*) flint_heap_alloc_bytes(alloc * sizeof(mpz_t));
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init(poly->coeffs[i]);
   
   poly->alloc = alloc;
   poly->length = 0;
}


void mpz_poly_init3(mpz_poly_t poly, unsigned long alloc, unsigned long bits)
{
   if ((long) alloc <= 0) mpz_poly_init(poly);
   
   poly->coeffs = (mpz_t*) flint_heap_alloc_bytes(alloc * sizeof(mpz_t));
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init2(poly->coeffs[i], bits);
   
   poly->alloc = alloc;
   poly->length = 0;
}


void mpz_poly_clear(mpz_poly_t poly)
{
   for (long i = 0; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);

   flint_heap_free(poly->coeffs);
}


void mpz_poly_realloc(mpz_poly_t poly, unsigned long alloc)
{
   if ((long) alloc <= 0) alloc = 1;
   
   // clear any mpz_t's beyond the new array length
   for (unsigned long i = alloc; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);

   poly->coeffs = (mpz_t*) flint_heap_realloc_bytes(poly->coeffs,
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

   poly->coeffs = (mpz_t*) flint_heap_realloc_bytes(poly->coeffs,
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
   FLINT_ASSERT(alloc > poly->alloc);

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

      for (unsigned long i = poly->length; i < n; i++)
         mpz_set_ui(poly->coeffs[i], 0);
         
      mpz_set(poly->coeffs[n], c);
      poly->length = n+1;
   }
}


void mpz_poly_set_coeff_ui(mpz_poly_t poly, unsigned long n, unsigned long c)
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
         mpz_set_ui(poly->coeffs[i], 0);
         
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
         mpz_set_ui(poly->coeffs[i], 0);
         
      mpz_set_si(poly->coeffs[n], c);
      poly->length = n+1;
   }
}



/****************************************************************************

   String conversions and I/O

****************************************************************************/


int mpz_poly_from_string(mpz_poly_t poly, const char* s)
{
   const char* whitespace = " \t\n\r";
   
   // read poly length
   unsigned long length;
   if (!sscanf(s, "%ld", &length))
      return 0;
      
   // jump to next whitespace
   s += strcspn(s, whitespace);
   
   poly->length = 0;
   mpz_poly_ensure_alloc(poly, length);

   for (unsigned long i = 0; i < length; i++)
   {
      // skip whitespace
      s += strspn(s, whitespace);
      
      if (!gmp_sscanf(s, "%Zd", poly->coeffs[i]))
         return 0;
      poly->length++;

      // jump to next whitespace
      s += strcspn(s, whitespace);
   }
   
   mpz_poly_normalise(poly);
   
   return 1;
}


char* mpz_poly_to_string(mpz_poly_t poly)
{
   // estimate the size of the string
   // 20 = enough room for null terminator and length info
   unsigned long size = 20;
   for (unsigned long i = 0; i < poly->length; i++)
      // +2 is for the sign and a space
      size += mpz_sizeinbase(poly->coeffs[i], 10) + 2;

   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "%ld  ", poly->length);
   for (unsigned long i = 0; i < poly->length; i++)
   {
      mpz_get_str(ptr, 10, poly->coeffs[i]);
      ptr += strlen(ptr);
      *ptr = ' ';
      ptr++;
   }
   
   ptr--;
   *ptr = 0;
   
   return buf;
}

/***************************************************************************************************
**  LTOA.C
**
**  Converts a long integer to a string.
**
**  Copyright 1988-90 by Robert B. Stout dba MicroFirm
**
**  Released to public domain, 1991
**
**  Parameters: 1 - number to be converted
**              2 - buffer in which to build the converted string
**              3 - number base to use for conversion
**
**  Returns:  A character pointer to the converted string if
**            successful, a NULL pointer if the number base specified
**            is out of range.
***************************************************************************************************/

#define BUFSIZE (sizeof(long) * 8 + 1)

char *flint_ltoa(long N, char *str, int base)
{
      register int i = 2;
      long uarg;
      char *tail, *head = str, buf[BUFSIZE];

      if (36 < base || 2 > base)
            base = 10;                    /* can only use 0-9, A-Z        */
      tail = &buf[BUFSIZE - 1];           /* last character position      */
      *tail-- = '\0';

      if (10 == base && N < 0L)
      {
            *head++ = '-';
            uarg    = -N;
      }
      else  uarg = N;

      if (uarg)
      {
            for (i = 1; uarg; ++i)
            {
                  register ldiv_t r;

                  r       = ldiv(uarg, base);
                  *tail-- = (char)(r.rem + ((9L < r.rem) ?
                                  ('A' - 10L) : '0'));
                  uarg    = r.quot;
            }
      }
      else  *tail-- = '0';

      memcpy(head, ++tail, i);
      return str;
}

/*******************************************************************************************/

char* mpz_poly_to_string_pretty(mpz_poly_t poly, const char * x)
{
   if (poly->length == 0)
   {
      char* buf = (char*) malloc(2);
      *buf = '0';
      buf[1] = 0;
      return buf;
   }
   
   unsigned long x_len = strlen(x); // String length of the monomial
   unsigned long exp_len = FLINT_BIT_COUNT(poly->length)/3 + 1; // String length of largest degree
   // estimate the size of the string
   // 1 = enough room for null terminator
   unsigned long size = 1;
   long i;
   
   for (i = 0; i < poly->length; i++)
      // +3 is for the sign, a carot and a times
      size += mpz_sizeinbase(poly->coeffs[i], 10) + 3 + x_len + exp_len;

   // write the string
   char* buf = (char*) malloc(size);
   char* exp = (char*) malloc(exp_len+1);
   char* ptr = buf;
   
   for (i = poly->length - 1; i >=2; i--)
   {
      if ((mpz_sgn(poly->coeffs[i]) > 0L) && (i != poly->length - 1))
      {
         *ptr = '+';
         ptr++;
      }
      if (mpz_cmp_si(poly->coeffs[i], -1L) == 0)
      {
         *ptr = '-';
         ptr++;
      }
      if (mpz_sgn(poly->coeffs[i]) != 0L)
      {
         if ((mpz_cmp_si(poly->coeffs[i], -1L) != 0) && (mpz_cmp_ui(poly->coeffs[i], 1L) != 0))
         {
            mpz_get_str(ptr, 10, poly->coeffs[i]);
            ptr += strlen(ptr);
            *ptr = '*';
            ptr++;
         }
         strcpy(ptr, x);
         ptr += strlen(x);
         *ptr = '^';
         ptr++;
         flint_ltoa(i, exp, 10);
         strcpy(ptr, exp);
         ptr += strlen(exp);
      }
   }
   
   if (i == 1)
   {
      if ((mpz_sgn(poly->coeffs[i]) > 0L) && (i != poly->length - 1))
      {
         *ptr = '+';
         ptr++;
      }
      if (mpz_cmp_si(poly->coeffs[i], -1L) == 0)
      {
         *ptr = '-';
         ptr++;
      }
      if (mpz_sgn(poly->coeffs[i]) != 0L)
      {
         if ((mpz_cmp_si(poly->coeffs[i], -1L) != 0) && (mpz_cmp_ui(poly->coeffs[i], 1L) != 0))
         {
            mpz_get_str(ptr, 10, poly->coeffs[i]);
            ptr += strlen(ptr);
            *ptr = '*';
            ptr++;
         }
         strcpy(ptr, x);
         ptr += strlen(x);
      }
      i--;
   }
   
   if ((mpz_sgn(poly->coeffs[i]) > 0L) && (i != poly->length - 1))
   {
      *ptr = '+';
      ptr++;
   }
   if (mpz_sgn(poly->coeffs[i]) != 0L)
   {
      mpz_get_str(ptr, 10, poly->coeffs[i]);
      ptr += strlen(ptr);
   }

   *ptr = 0;
   
   return buf;
}


void mpz_poly_fprint(mpz_poly_t poly, FILE* f)
{
   char* s = mpz_poly_to_string(poly);
   fputs(s, f);
   free(s);
}

void mpz_poly_fprint_pretty(mpz_poly_t poly, FILE* f, const char * x)
{
   char* s = mpz_poly_to_string_pretty(poly, x);
   fputs(s, f);
   free(s);
}


void mpz_poly_print(mpz_poly_t poly)
{
   mpz_poly_fprint(poly, stdout);
}

void mpz_poly_print_pretty(mpz_poly_t poly, const char * x)
{
   mpz_poly_fprint_pretty(poly, stdout, x);
}


int mpz_poly_fread(mpz_poly_t poly, FILE* f)
{
   // read poly length
   unsigned long length;
   if (!fscanf(f, "%ld", &length))
      return 0;

   poly->length = 0;
   mpz_poly_ensure_alloc(poly, length);

   // read coefficients
   for (unsigned long i = 0; i < length; i++)
   {
      if (!mpz_inp_str(poly->coeffs[i], f, 10))
         return 0;
      poly->length++;
   }

   mpz_poly_normalise(poly);
   
   return 1;
}


int mpz_poly_read(mpz_poly_t poly)
{
   return mpz_poly_fread(poly, stdin);
}


/****************************************************************************

   Length and degree

****************************************************************************/


void mpz_poly_normalise(mpz_poly_t poly)
{
   while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]))
      poly->length--;
}


int mpz_poly_normalised(mpz_poly_t poly)
{
   return (poly->length == 0) || mpz_sgn(poly->coeffs[poly->length-1]);
}


void mpz_poly_pad(mpz_poly_t poly, unsigned long length)
{
   mpz_poly_ensure_alloc(poly, length);

   if (poly->length < length)
   {
      for (unsigned long i = poly->length; i < length; i++)
         mpz_set_ui(poly->coeffs[i], 0);
      poly->length = length;
   }
}


void mpz_poly_truncate(mpz_poly_t res, mpz_poly_t poly, unsigned long length)
{
   if (poly == res)
   {
      // inplace truncation

      if (length < poly->length)
         poly->length = length;
   }
   else
   {
      // copy and truncate

      if (length > poly->length)
      {
         mpz_poly_set(res, poly);
         return;
      }

      // todo: use mpz_init_set where appropriate
      
      mpz_poly_ensure_alloc(res, length);

      for (unsigned long i = 0; i < length; i++)
         mpz_set(res->coeffs[i], poly->coeffs[i]);
         
      res->length = length;
   }
   
   mpz_poly_normalise(res);
}



/****************************************************************************

   Assignment

****************************************************************************/


void mpz_poly_set(mpz_poly_t res, mpz_poly_t poly)
{
   if (res == poly)
      return;

   // todo: use mpz_init_set where appropriate

   mpz_poly_ensure_alloc(res, poly->length);
   
   for (unsigned long i = 0; i < poly->length; i++)
      mpz_set(res->coeffs[i], poly->coeffs[i]);
      
   res->length = poly->length;
}


/****************************************************************************

   Conversions

****************************************************************************/


// assumes coefficients are big enough, and alloc is big enough
void _mpz_poly_to_fmpz_poly(fmpz_poly_t res, mpz_poly_t poly)
{
   FLINT_ASSERT(res->alloc >= poly->length);

   res->length = poly->length;
   if (poly->length == 0)
      return;

   for (unsigned long i = 0; i < poly->length; i++)
   {
      FLINT_ASSERT(res->limbs >= mpz_size(poly->coeffs[i]));
      mpz_to_fmpz(res->coeffs + i*(res->limbs+1), poly->coeffs[i]);
   }
}


void mpz_poly_to_fmpz_poly(fmpz_poly_t res, mpz_poly_t poly)
{
   unsigned long limbs = mpz_poly_max_limbs(poly);

   // todo: there should be a single function that achieves both of the
   // following.... actually we don't even care in this case if the value
   // is preserved.
   fmpz_poly_fit_length(res, poly->length);
   fmpz_poly_fit_limbs(res, limbs);

   _mpz_poly_to_fmpz_poly(res, poly);
}


void fmpz_poly_to_mpz_poly(mpz_poly_t res, const fmpz_poly_t poly)
{
   mpz_poly_ensure_alloc(res, poly->length);

   res->length = poly->length;
   
   if (poly->length == 0) return;
   
   long i;
   mp_limb_t* ptr = poly->coeffs;

   for (i = 0; i < poly->length; i++, ptr += poly->limbs+1)
      fmpz_to_mpz(res->coeffs[i], ptr);
   
   mpz_poly_normalise(res);

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



/****************************************************************************

   Shifting

****************************************************************************/


void mpz_poly_lshift(mpz_poly_t res, mpz_poly_t poly, unsigned long k)
{
   mpz_poly_ensure_alloc(res, poly->length + k);

   if (poly == res)
   {
      // inplace; just shift the mpz_t's over
      for (long i = poly->length - 1; i >= 0; i--)
         mpz_swap(poly->coeffs[i], poly->coeffs[i+k]);
      
      for (unsigned long i = 0; i < k; i++)
         mpz_set_ui(poly->coeffs[i], 0);
   }
   else
   {
      // not inplace; need to copy data
      for (unsigned long i = 0; i < k; i++)
         mpz_set_ui(res->coeffs[i], 0);
      
      for (unsigned long i = 0; i < poly->length; i++)
         mpz_set(res->coeffs[i + k], poly->coeffs[i]);
   }
   
   res->length = poly->length + k;
}


void mpz_poly_rshift(mpz_poly_t res, mpz_poly_t poly, unsigned long k)
{
   if (k >= poly->length)
   {
      // shift all coefficients off the end
      res->length = 0;
      return;
   }

   if (poly == res)
   {
      // inplace; just shift the mpz_t's over

      for (unsigned long i = k; i < poly->length; i++)
         mpz_swap(poly->coeffs[i - k], poly->coeffs[i]);
   }
   else
   {
      // not inplace; need to copy data
      mpz_poly_ensure_alloc(res, poly->length - k);

      for (unsigned long i = k; i < poly->length; i++)
         mpz_set(res->coeffs[i - k], poly->coeffs[i]);
   }
   
   res->length = poly->length - k;
}


/****************************************************************************

   Norms

****************************************************************************/

void mpz_poly_2norm(mpz_t norm, mpz_poly_t poly)
{
   mpz_set_ui(norm, 0L);
   
   if (poly->length == 0)
   {
      return;
   }   

   mpz_t sqr;
   mpz_init(sqr);
   
   for (unsigned long i = 0; i < poly->length; i++)
   {
      mpz_mul(sqr, poly->coeffs[i], poly->coeffs[i]);
      mpz_add(norm, norm, sqr);
   }

   mpz_sqrtrem(norm, sqr, norm);
   if (mpz_sgn(sqr) != 0L) mpz_add_ui(norm, norm, 1L);

   mpz_clear(sqr);
}

/****************************************************************************

   Scalar multiplication and division

****************************************************************************/


void mpz_poly_scalar_mul(mpz_poly_t res, mpz_poly_t poly, mpz_t c)
{
   abort();
}

void mpz_poly_scalar_mul_ui(mpz_poly_t res, mpz_poly_t poly, unsigned long c)
{
   abort();
}


void mpz_poly_scalar_mul_si(mpz_poly_t res, mpz_poly_t poly, long c)
{
   abort();
}


void mpz_poly_scalar_div(mpz_poly_t res, mpz_poly_t poly, mpz_t c)
{
   abort();
}


void mpz_poly_scalar_div_ui(mpz_poly_t res, mpz_poly_t poly, unsigned long c)
{
   abort();
}


void mpz_poly_scalar_div_si(mpz_poly_t res, mpz_poly_t poly, long c)
{
   abort();
}


void mpz_poly_scalar_div_exact(mpz_poly_t res, mpz_poly_t poly, mpz_t c)
{
   abort();
}


void mpz_poly_scalar_div_exact_ui(mpz_poly_t res, mpz_poly_t poly,
                                  unsigned long c)
{
   abort();
}


void mpz_poly_scalar_div_exact_si(mpz_poly_t res, mpz_poly_t poly, long c)
{
   abort();
}


void mpz_poly_scalar_mod(mpz_poly_t res, mpz_poly_t poly, mpz_t c)
{
   abort();
}


void mpz_poly_scalar_mod_ui(mpz_poly_t res, mpz_poly_t poly, unsigned long c)
{
   abort();
}



/****************************************************************************

   Polynomial multiplication

****************************************************************************/


void mpz_poly_mul(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   // use naive KS for now
   mpz_poly_mul_naive_KS(res, poly1, poly2);
}

void mpz_poly_sqr(mpz_poly_t res, mpz_poly_t poly)
{
   // use naive KS for now
   mpz_poly_sqr_naive_KS(res, poly);
}


/*
 This is just like mpz_poly_mul_classical(), with the following restrictions:
 
  * assumes res does not alias poly1 and poly2
  * neither polynomial is zero
  * res->alloc >= poly1->length + poly2->length - 1
     (i.e. output has enough room for product)
*/
void _mpz_poly_mul_classical(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   FLINT_ASSERT(res != poly1);
   FLINT_ASSERT(res != poly2);
   FLINT_ASSERT(poly1->length && poly2->length);

   res->length = poly1->length + poly2->length - 1;
   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      mpz_set_ui(res->coeffs[i], 0);
   
   for (unsigned long i = 0; i < poly1->length; i++)
      for (unsigned long j = 0; j < poly2->length; j++)
         mpz_addmul(res->coeffs[i+j], poly1->coeffs[i], poly2->coeffs[j]);
}


void mpz_poly_mul_classical(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   if (!poly1->length || !poly2->length)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }
   
   if (poly1 == poly2)
   {
      // polys are identical, so call specialised squaring routine
      mpz_poly_sqr_classical(res, poly1);
      return;
   }
   
   unsigned long limbs = mpz_poly_product_max_limbs(poly1, poly2);
   unsigned long length = poly1->length + poly2->length - 1;
   
   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      mpz_poly_t temp;
      mpz_poly_init3(temp, length, FLINT_BITS * limbs);
      _mpz_poly_mul_classical(temp, poly1, poly2);
      mpz_poly_swap(temp, res);
      mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      mpz_poly_ensure_alloc(res, length);
      _mpz_poly_mul_classical(res, poly1, poly2);
   }
}


/*
 This is just like mpz_poly_sqr_classical(), with the following restrictions:
 
  * assumes res does not alias poly
  * poly is nonzero
  * res->alloc >= 2*poly->length - 1  (i.e. output has enough room for product)
*/
void _mpz_poly_sqr_classical(mpz_poly_t res, mpz_poly_t poly)
{
   FLINT_ASSERT(res != poly);
   FLINT_ASSERT(poly->length);

   res->length = 2*poly->length - 1;
   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      mpz_set_ui(res->coeffs[i], 0);
   
   // off-diagonal products
   for (unsigned long i = 1; i < poly->length; i++)
      for (unsigned long j = 0; j < i; j++)
         mpz_addmul(res->coeffs[i+j], poly->coeffs[i], poly->coeffs[j]);
         
   // double the off-diagonal products
   for (unsigned long i = 1; i < res->length - 1; i++)
      mpz_add(res->coeffs[i], res->coeffs[i], res->coeffs[i]);
      
   // add in diagonal products
   for (unsigned long i = 0; i < poly->length; i++)
      mpz_addmul(res->coeffs[2*i], poly->coeffs[i], poly->coeffs[i]);
}


void mpz_poly_sqr_classical(mpz_poly_t res, mpz_poly_t poly)
{
   if (!poly->length)
   {
      // input is zero
      res->length = 0;
      return;
   }
   
   unsigned long limbs = mpz_poly_product_max_limbs(poly, poly);
   unsigned long length = 2*poly->length - 1;
   
   if (res == poly)
   {
      // output is inplace, so need a temporary
      mpz_poly_t temp;
      mpz_poly_init3(temp, length, FLINT_BITS * limbs);
      _mpz_poly_sqr_classical(temp, poly);
      mpz_poly_swap(temp, res);
      mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      // allocate more coefficients if necessary
      mpz_poly_ensure_alloc(res, length);
      _mpz_poly_sqr_classical(res, poly);
   }
}



/*
Recursive portion of karatsuba multiplication.

Input arrays are "in1" of length "len1", staggered by "skip", ditto for in2.
Output array is "out" of length len1 + len2 - 1, also staggered by "skip".
scratch buffer should be length len1 + len2, also staggered by "skip".

All input/output/scratch space should be mpz_init'd, and shouldn't overlap.

Must have 1 <= len1 <= len2.

If len1*len2 <= crossover, it uses a classical multiplication algorithm. The
crossover parameter is passed down recursively to subproducts.
*/
void _mpz_poly_mul_kara_recursive(mpz_t* out,
                                  mpz_t* in1, unsigned long len1,
                                  mpz_t* in2, unsigned long len2,
                                  mpz_t* scratch, unsigned long skip,
                                  unsigned long crossover)
{
   FLINT_ASSERT(len1 >= 1);
   FLINT_ASSERT(len2 >= len1);

   // ==================== base cases

   if (len1 == 1)
   {
      // special case, just scalar multiplication
      for (unsigned long i = 0; i < len2; i++)
         mpz_mul(out[i*skip], in1[0], in2[i*skip]);
      return;
   }
   
   if (len1 * len2 < crossover)
   {
      // switch to classical multiplication
      for (unsigned long i = 0; i < len1 + len2 - 1; i++)
         mpz_set_ui(out[i*skip], 0);
   
      for (unsigned long i = 0; i < len1; i++)
         for (unsigned long j = 0; j < len2; j++)
            mpz_addmul(out[(i+j)*skip], in1[i*skip], in2[j*skip]);
            
      return;
   }

   // ==================== recursive case

   // Let in1 = A1(x^2) + x*B1(x^2) + x^(2*floor(len1/2))*C1,
   // where A1, B1 have length floor(len1/2),
   // and C1 is the leading term of in1 if len1 is odd

   // Similarly for in2 = A2(x^2) + x*B2(x^2) + x^(2*floor(len2/2))*C2.
   
   // Put A1 + B1 into even slots of scratch space
   // (uses len1/2 scratch slots)
   mpz_t* ptr = scratch;
   for (unsigned long i = 0; i < len1/2; i++, ptr += 2*skip)
      mpz_add(*ptr, in1[2*i*skip], in1[2*i*skip + skip]);

   // Put A2 + B2 into remaining even slots of scratch space
   // (uses len2/2 slots of scratch)
   mpz_t* scratch2 = ptr;
   for (unsigned long i = 0; i < len2/2; i++, ptr += 2*skip)
      mpz_add(*ptr, in2[2*i*skip], in2[2*i*skip + skip]);

   // The following three recursive calls all use the odd slots of the current
   // scratch array as the next layer's scratch space
   
   // Put product (A1+B1)*(A2+B2) into odd slots of output array
   _mpz_poly_mul_kara_recursive(out + skip, scratch, len1/2, scratch2, len2/2,
                                scratch + skip, 2*skip, crossover);

   // Put product x^2*(B1*B2) into even slots of output array
   // (except first slot, which is an implied zero)
   _mpz_poly_mul_kara_recursive(out + 2*skip, in1 + skip, len1/2, in2 + skip,
                                len2/2, scratch + skip, 2*skip, crossover);

   // Put product A1*A2 into even slots of scratch space
   _mpz_poly_mul_kara_recursive(scratch, in1, len1/2, in2, len2/2,
                                scratch + skip, 2*skip, crossover);
                            
   // Subtract A1*A2 and B1*B2 from (A1+B1)*(A2+B2) to get (A1*B2 + A2*B1)
   // in odd slots of output
   for (unsigned long i = 0; i < len1/2 + len2/2 - 1; i++)
   {
      mpz_sub(out[2*i*skip + skip], out[2*i*skip + skip], out[2*(i+1)*skip]);
      mpz_sub(out[2*i*skip + skip], out[2*i*skip + skip], scratch[2*i*skip]);
   }
      
   // Add A1*A2 to x^2*(B1*B2) into even slots of output
   mpz_set(out[0], scratch[0]);
   for (unsigned long i = 1; i < len1/2 + len2/2 - 1; i++)
      mpz_add(out[2*i*skip], out[2*i*skip], scratch[2*i*skip]);
   
   // Now we have the product (A1(x^2) + x*B1(x^2)) * (A2(x^2) + x*B2(x^2))
   // in the output array. Still need to handle C1 and C2 terms.
   
   if (len1 & 1)
   {
      if (len2 & 1)
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         mpz_t* term1 = in1 + skip*(len1-1);
         for (unsigned long i = 0; i < len2-2; i++)
            mpz_addmul(out[(i+len1-1)*skip], *term1, in2[i*skip]);
         mpz_mul(out[(len1+len2-3)*skip], *term1, in2[(len2-2)*skip]);

         // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
         mpz_t* term2 = in2 + skip*(len2-1);
         for (unsigned long i = 0; i < len1-1; i++)
            mpz_addmul(out[(i+len2-1)*skip], *term2, in1[i*skip]);
            
         // final C1*C2 term
         mpz_mul(out[(len1+len2-2)*skip], *term1, *term2);
      }
      else
      {
         // terms from x^(len1-1)*C1 * (A2(x^2) + x*B2(x^2))
         mpz_t* term = in1 + skip*(len1-1);
         for (unsigned long i = 0; i < len2-1; i++)
            mpz_addmul(out[(i+len1-1)*skip], *term, in2[i*skip]);
         mpz_mul(out[(len1+len2-2)*skip], *term, in2[(len2-1)*skip]);
      }
   }
   else if (len2 & 1)
   {
      // terms from x^(len2-1)*C2 * (A1(x^2) + x*B1(x^2))
      mpz_t* term = in2 + skip*(len2-1);
      for (unsigned long i = 0; i < len1-1; i++)
         mpz_addmul(out[(i+len2-1)*skip], *term, in1[i*skip]);
      mpz_mul(out[(len1+len2-2)*skip], *term, in1[(len1-1)*skip]);
   }
}


unsigned long _mpz_poly_mul_karatsuba_crossover(unsigned long limbs)
{
   unsigned long crossover;
   if (limbs >= mpz_poly_kara_crossover_table_size)
      crossover = 0;
   else
   {
      if (limbs == 0)
         limbs = 1;
      crossover = mpz_poly_kara_crossover_table[limbs - 1];
   }
   return crossover * crossover;
}


void mpz_poly_mul_karatsuba(mpz_poly_t res, mpz_poly_t poly1,
                            mpz_poly_t poly2)
{
   if (!poly1->length || !poly2->length)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }
   
   if (poly1 == poly2)
   {
      // polys are identical, so call specialised squaring routine
      mpz_poly_sqr_karatsuba(res, poly1);
      return;
   }

   // rearrange parameters to make poly1 no longer than poly2
   if (poly1->length > poly2->length)
      SWAP_MPZ_POLY_PTRS(poly1, poly2);

   // number of output coefficients, and a rough upper bound on the number
   // of limbs needed for each one
   unsigned long length = poly1->length + poly2->length - 1;
   unsigned long limbs = mpz_poly_product_max_limbs(poly1, poly2);
   
   // allocate scratch space for lower-level karatsuba routine
   mpz_t* scratch = (mpz_t*)
                     flint_stack_alloc_bytes((length+1) * sizeof(mpz_t));
   for (unsigned long i = 0; i <= length; i++)
      mpz_init2(scratch[i], limbs * FLINT_BITS);

   // look up crossover parameter (i.e. when to switch from classical to
   // karatsuba multiplication) based on coefficient size
   unsigned long crossover = _mpz_poly_mul_karatsuba_crossover(limbs/2);
   
   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      mpz_poly_t temp;
      mpz_poly_init3(temp, length, FLINT_BITS*limbs);

      _mpz_poly_mul_kara_recursive(
            temp->coeffs, poly1->coeffs, poly1->length,
            poly2->coeffs, poly2->length, scratch, 1, crossover);

      mpz_poly_swap(temp, res);
      mpz_poly_clear(temp);
   }
   else
   {
      // output not inplace
      
      // allocate more coefficients if necessary
      mpz_poly_ensure_alloc(res, length);

      _mpz_poly_mul_kara_recursive(
            res->coeffs, poly1->coeffs, poly1->length,
            poly2->coeffs, poly2->length, scratch, 1, crossover);
   }
   
   res->length = length;
   
   for (unsigned long i = 0; i <= length; i++)
      mpz_clear(scratch[i]);
   flint_stack_release();
}


void mpz_poly_mul_SS(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_sqr_SS(mpz_poly_t res, mpz_poly_t poly)
{
   abort();
}

void mpz_poly_sqr_karatsuba(mpz_poly_t res, mpz_poly_t poly)
{
   abort();
}


// ============================================================================
//
//   Naive KS multiplication and support routines

/*
   Sets y = \sum_{i=0}^{len-1} x[i] * 2^(ki)
*/

void mpz_poly_mul_naive_KS_pack(mpz_t y, mpz_t* x, unsigned long len,
                                unsigned long k)
{
   if (len == 1)
      mpz_set(y, x[0]);
   else
   {
      mpz_t temp;
      mpz_init(temp);
      unsigned long half = len/2;
      mpz_poly_mul_naive_KS_pack(temp, x, half, k);
      mpz_poly_mul_naive_KS_pack(y, x + half, len - half, k);
      mpz_mul_2exp(y, y, half*k);
      mpz_add(y, y, temp);
      mpz_clear(temp);
   }
}


/*
   Inverse operation of mpz_poly_mul_naive_KS_pack
   (note: y is destroyed)
*/

void mpz_poly_mul_naive_KS_unpack(mpz_t* x, unsigned long len, mpz_t y,
                                  unsigned long k)
{
   if (len == 1)
      mpz_set(x[0], y);
   else
   {
      mpz_t temp;
      mpz_init(temp);
      unsigned long half = len/2;
      if (mpz_tstbit(y, k*half - 1))
      {
         mpz_cdiv_q_2exp(temp, y, half*k);
         mpz_cdiv_r_2exp(y, y, half*k);
      }
      else
      {
         mpz_fdiv_q_2exp(temp, y, half*k);
         mpz_fdiv_r_2exp(y, y, half*k);
      }
      mpz_poly_mul_naive_KS_unpack(x, half, y, k);
      mpz_poly_mul_naive_KS_unpack(x + half, len - half, temp, k);
      mpz_clear(temp);
   }
}


/*
   Counts maximum number of bits in abs(x->coeffs[i])
   todo: isn't this subsumed into mpz_poly_max_bits()?
*/

unsigned long mpz_poly_mul_naive_KS_get_max_bits(mpz_poly_t x)
{
   unsigned long bits = 0, temp, i;
   for (i = 0; i < x->length; i++)
   {
      temp = mpz_sizeinbase(x->coeffs[i], 2);
      if (temp > bits)
         bits = temp;
   }
   return bits;
}


void mpz_poly_mul_naive_KS(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   if (poly1 == poly2)
   {
      mpz_poly_sqr_naive_KS(res, poly1);
      return;
   }

   if (!poly1->length || !poly2->length)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }
   
   mpz_t z1;
   mpz_t z2;
   mpz_init(z1);
   mpz_init(z2);

   unsigned long out_len = poly1->length + poly2->length - 1;
   unsigned long bits1 = mpz_poly_mul_naive_KS_get_max_bits(poly1);
   unsigned long bits2 = mpz_poly_mul_naive_KS_get_max_bits(poly2);
   unsigned long bits = bits1 + bits2 + 1 +
                        ceil_log2(FLINT_MIN(poly1->length, poly2->length));

   mpz_poly_mul_naive_KS_pack(z1, poly1->coeffs, poly1->length, bits);
   mpz_poly_mul_naive_KS_pack(z2, poly2->coeffs, poly2->length, bits);
   mpz_mul(z1, z1, z2);
   mpz_poly_ensure_alloc(res, out_len);
   mpz_poly_mul_naive_KS_unpack(res->coeffs, out_len, z1, bits);
   res->length = out_len;

   mpz_clear(z1);
   mpz_clear(z2);
}



void mpz_poly_sqr_naive_KS(mpz_poly_t res, mpz_poly_t poly)
{
   if (!poly->length)
   {
      // poly is zero
      res->length = 0;
      return;
   }
   
   mpz_t z;
   mpz_init(z);

   unsigned long out_len = 2*poly->length - 1;
   unsigned long bits = 2 * mpz_poly_mul_naive_KS_get_max_bits(poly)
                          + 1 + ceil_log2(poly->length);

   mpz_poly_mul_naive_KS_pack(z, poly->coeffs, poly->length, bits);
   mpz_mul(z, z, z);
   mpz_poly_ensure_alloc(res, out_len);
   mpz_poly_mul_naive_KS_unpack(res->coeffs, out_len, z, bits);
   res->length = out_len;
   
   mpz_clear(z);
}



/****************************************************************************

   Polynomial division

****************************************************************************/

/*
Input is a monic polynomial "poly" of degree n, and a nonzero polynomial Q1 of
degree k1 such that
    x^(k1+n) = poly*Q1 + R
where deg(R) < n.

Output is a nonzero polynomial Q2 of degree k2 such that
    x^(k2+n) = poly*Q2 + S
where deg(S) < n.

PRECONDITIONS:
   k2 >= k1
   poly and Q1 must be normalised
   Q1, Q2, poly must not alias each other

*/
void mpz_poly_monic_inverse_newton_extend(
             mpz_poly_t Q2, mpz_poly_t Q1, mpz_poly_t poly, unsigned long k2)
{
   FLINT_ASSERT(poly != Q1);
   FLINT_ASSERT(poly != Q2);
   FLINT_ASSERT(Q1 != Q2);
   FLINT_ASSERT(mpz_poly_normalised(poly));
   FLINT_ASSERT(mpz_poly_normalised(Q1));
   FLINT_ASSERT(Q1->length >= 1);
   
   unsigned long k1 = Q1->length - 1;
   FLINT_ASSERT(k2 >= k1);
   
   unsigned long n = poly->length - 1;

   if (k2 <= 2*k1)
   {
      // only one newton iteration is needed
      
      // temp := top k2+1 coefficients of Q1^2
      mpz_poly_t temp;
      mpz_poly_init(temp);
      mpz_poly_sqr(temp, Q1);
      mpz_poly_rshift(temp, temp, temp->length - (k2+1));
      
      // temp := top k2+1 coefficients of Q1^2*poly
      if (poly->length > k2+1)
      {
         // first get top k2+1 coefficients of poly
         mpz_poly_t top;
         mpz_poly_init(top);
         mpz_poly_rshift(top, poly, poly->length - (k2+1));

         // now get top k2+1 coefficients of Q1^2*poly
         mpz_poly_mul(temp, temp, top);
         mpz_poly_rshift(temp, temp, temp->length - (k2+1));
         
         mpz_poly_clear(top);
      }
      else
      {
         mpz_poly_mul(temp, temp, poly);
         mpz_poly_rshift(temp, temp, temp->length - (k2+1));
      }
      
      // Q2 = top k2+1 coefficients of 2*Q1*x^(k1+n) - Q1^2*poly
      mpz_poly_ensure_alloc(Q2, k2+1);
      mpz_t x;
      mpz_init(x);

      unsigned long i;
      for (i = 0; i <= k1; i++)
      {
         mpz_add(x, Q1->coeffs[k1-i], Q1->coeffs[k1-i]);
         mpz_sub(Q2->coeffs[k2-i], x, temp->coeffs[k2-i]);
      }
      for (; i <= k2; i++)
      {
         mpz_neg(Q2->coeffs[k2-i], temp->coeffs[k2-i]);
      }

      Q2->length = k2+1;

      mpz_clear(x);
      mpz_poly_clear(temp);
   }
   else
   {
      // more than one newton iteration is needed, so recurse
      mpz_poly_t temp;
      mpz_poly_init(temp);
      mpz_poly_monic_inverse_newton_extend(temp, Q1, poly, (k2+1)/2);
      mpz_poly_monic_inverse_newton_extend(Q2, temp, poly, k2);
      mpz_poly_clear(temp);
   }
}


void mpz_poly_monic_inverse(mpz_poly_t res, mpz_poly_t poly, unsigned long k)
{
   // todo: remove the following restrictions
   FLINT_ASSERT(k >= 2);
   FLINT_ASSERT(poly->length >= 2);
   FLINT_ASSERT(poly != res);

   // if poly is x^n + a*x^(n-1) + ..., then first approximation
   // to res is given by x - a
   mpz_poly_t temp;
   mpz_poly_init(temp);
   mpz_poly_pad(temp, 2);
   mpz_set_ui(temp->coeffs[1], 1);
   mpz_neg(temp->coeffs[0], poly->coeffs[poly->length-2]);
   temp->length = 2;

   // extend the approximation using newton's method
   mpz_poly_monic_inverse_newton_extend(res, temp, poly, k);
   mpz_poly_clear(temp);
}



void mpz_poly_pseudo_inverse(mpz_poly_t res, mpz_poly_t poly, unsigned long k)
{
   abort();
}

void mpz_poly_monic_div(mpz_poly_t quot, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_div(mpz_poly_t quot, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_rem(mpz_poly_t rem, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_rem(mpz_poly_t rem, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_div_rem(mpz_poly_t quot, mpz_poly_t rem,
                            mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_div_rem(mpz_poly_t quot, mpz_poly_t rem, 
                             mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_inverse_basecase(mpz_poly_t res, mpz_poly_t poly,
                                  unsigned long k)
{
   abort();
}

void mpz_poly_pseudo_inverse_basecase(mpz_poly_t res, mpz_poly_t poly,
                                   unsigned long k)
{
   abort();
}

void mpz_poly_monic_div_basecase(mpz_poly_t quot, mpz_poly_t poly1,
                              mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_div_basecase(mpz_poly_t quot, mpz_poly_t poly1,
                               mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_rem_basecase(mpz_poly_t rem, mpz_poly_t poly1,
                              mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_rem_basecase(mpz_poly_t rem, mpz_poly_t poly1,
                               mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_div_rem_basecase(mpz_poly_t quot, mpz_poly_t rem,
                                  mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_div_rem_basecase(mpz_poly_t quot, mpz_poly_t rem, 
                                   mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}



/****************************************************************************

   GCD and extended GCD

****************************************************************************/


void mpz_poly_content(mpz_t x, mpz_poly_t poly)
{
   abort();
}


unsigned long mpz_poly_content_ui(mpz_poly_t poly)
{
   abort();
}


void mpz_poly_gcd(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}


void mpz_poly_xgcd(mpz_poly_t res, mpz_poly_t a, mpz_poly_t b,
                   mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}



/****************************************************************************

   Miscellaneous

****************************************************************************/


unsigned long mpz_poly_max_limbs(mpz_poly_t poly)
{
   if (!poly->length)
      return 0;
   
   unsigned long temp, limbs = mpz_size(poly->coeffs[0]);
   
   for (unsigned long i = 1; i < poly->length; i++)
   {
      temp = mpz_size(poly->coeffs[i]);
      if (temp > limbs)
         limbs = temp;
   }

   return limbs;
}


unsigned long mpz_poly_max_bits(mpz_poly_t poly)
{
   abort();
}


unsigned long mpz_poly_product_max_limbs(mpz_poly_t poly1, mpz_poly_t poly2)
{
   unsigned long limbs1 = mpz_poly_max_limbs(poly1);
   unsigned long limbs2 = mpz_poly_max_limbs(poly2);

   // we're assuming poly lengths are at most 2^FLINT_BITS
   return limbs1 + limbs2 + 1;
}


unsigned long mpz_poly_product_max_bits(mpz_poly_t poly1, mpz_poly_t poly2)
{
   unsigned long bits1 = mpz_poly_max_bits(poly1);
   unsigned long bits2 = mpz_poly_max_bits(poly2);
   
   return bits1 + bits2 + ceil_log2(FLINT_MAX(poly1->length, poly2->length));
}


// *************** end of file
