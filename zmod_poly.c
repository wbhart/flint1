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

   zmod_poly->c: Polynomials over (unsigned) long mod p, for p prime.
   
   Copyright (C) 2007, David Howden.
   Copyright (C) 2007, 2008, William Hart
   Copyright (C) 2008, Richard Howell-Peak
   
*****************************************************************************/

#include "zmod_poly.h"
#include "zmod_mat.h"
#include "long_extras.h"
#include "flint.h"

#define PRINT_LIMB(a) print_limb(#a, a);
#define PRINT_VAR(a) print_var(#a, a);

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void zmod_poly_init(zmod_poly_t poly, unsigned long p)
{
   zmod_poly_init_precomp(poly, p, z_precompute_inverse(p));
}

void zmod_poly_init_precomp(zmod_poly_t poly, unsigned long p, double p_inv)
{
   poly->coeffs = (unsigned long*) flint_heap_alloc(1);
   
   poly->p = p;
   poly->p_inv = p_inv;
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) <= 32) poly->p32_inv = z_precompute_inverse32(p);
#endif
   poly->alloc = 1;
   poly->length = 0;
}

void zmod_poly_init2(zmod_poly_t poly, unsigned long p, unsigned long alloc)
{
   zmod_poly_init2_precomp(poly, p, z_precompute_inverse(p), alloc);
}

void zmod_poly_init2_precomp(zmod_poly_t poly, unsigned long p, double p_inv, unsigned long alloc)
{
   FLINT_ASSERT(alloc >= 1);

   poly->coeffs = (unsigned long*) flint_heap_alloc(alloc);
   
   poly->p = p;
   poly->p_inv = p_inv;
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) <= 32) poly->p32_inv = z_precompute_inverse32(p);
#endif
   
   poly->alloc = alloc;
   poly->length = 0;
}


void zmod_poly_clear(zmod_poly_t poly)
{
   flint_heap_free(poly->coeffs);
}


void zmod_poly_realloc(zmod_poly_t poly, unsigned long alloc)
{
   FLINT_ASSERT(alloc >= 1);

   // clear any mpz_t's beyond the new array length
   // for (unsigned long i = alloc; i < poly->alloc; i++)
   //    mpz_clear(poly->coeffs[i]);

   poly->coeffs = (unsigned long*) flint_heap_realloc(poly->coeffs,
                                              alloc);
   
   // init any new mpz_t's required
   // for (unsigned long i = poly->alloc; i < alloc; i++)
   //    mpz_init(poly->coeffs[i]);

   poly->alloc = alloc;
   
   // truncate poly if necessary
   if (poly->length > alloc)
   {
      poly->length = alloc;
      __zmod_poly_normalise(poly);
   }
}


void __zmod_poly_fit_length(zmod_poly_t poly, unsigned long alloc)
{
   FLINT_ASSERT(alloc > poly->alloc);

   if (alloc < 2*poly->alloc)
      alloc = 2*poly->alloc;
   zmod_poly_realloc(poly, alloc);
}

void zmod_poly_2x2_mat_init(zmod_poly_2x2_mat_t mat, ulong modulus)
{
   zmod_poly_init(mat->a, modulus);
   double p_inv = mat->a->p_inv;
	zmod_poly_init_precomp(mat->b, modulus, p_inv);
   zmod_poly_init_precomp(mat->c, modulus, p_inv);
   zmod_poly_init_precomp(mat->d, modulus, p_inv);
}

void zmod_poly_2x2_mat_clear(zmod_poly_2x2_mat_t mat)
{
   zmod_poly_clear(mat->a);
   zmod_poly_clear(mat->b);
   zmod_poly_clear(mat->c);
   zmod_poly_clear(mat->d);
}

/****************************************************************************

   Setting/retrieving coefficients

****************************************************************************/

void zmod_poly_set_coeff_ui(zmod_poly_t poly, unsigned long n, unsigned long c)
{
   c = z_mod_precomp(c, poly->p, poly->p_inv);
   
   zmod_poly_fit_length(poly, n+1);
   
   if (n+1 < poly->length)
      // set interior coefficient
      poly->coeffs[n] = c;

   else if (n+1 == poly->length)
   {
      // set leading coefficient
      if (c)
         poly->coeffs[n] = c;
      else
      {
         // set leading coefficient to zero
         poly->length--;
         __zmod_poly_normalise(poly);
      }
   }
   
   else
   {
      // extend polynomial
      if (!c)
         return;
      
      for (unsigned long i = poly->length; i < n; i++)
         poly->coeffs[i] = 0;
         
      poly->coeffs[n] = c;
      poly->length = n+1;
   }
}

/****************************************************************************

   String conversions and I/O

****************************************************************************/


/*
   Create a zmod_poly_t object from a string.
   
   Format: <Length> <Mod> <Coeffs>
*/

int zmod_poly_from_string(zmod_poly_t poly, char* s)
{
   const char* whitespace = " \t\n\r";

   unsigned long p, length;
   if (!sscanf(s, "%lx %lx", &length, &p))
      return 0;
      
   poly->p = p;
   poly->p_inv = z_precompute_inverse(p);
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) <= 32) poly->p32_inv = z_precompute_inverse32(p);
#endif

   // jump to next whitespace
   s += strcspn(s, whitespace);
   
   poly->length = 0;
   zmod_poly_fit_length(poly, length);
   
   for (unsigned long i = 0; i < length; i++)
   {
      // skip whitespace
      s += strspn(s, whitespace);
      
      if (!sscanf(s, "%ld", &poly->coeffs[i]))
         return 0;
      poly->length++;

      // jump to next whitespace
      s += strcspn(s, whitespace);
   }
   
   __zmod_poly_normalise(poly);
   
   return 1;
}


/*
   Convert a zmod_poly into a string.
   
   Format: <Length> <Mod> <Coeffs>
*/

char* zmod_poly_to_string(zmod_poly_t poly)
{
   // estimate the size of the string
   // 20 = enough room for null terminator and length info
   // and another 20 for p value...
   unsigned long size = 20*(2+poly->length);
   for (unsigned long i = 0; i < poly->length; i++)
   {
      // +2 is for the sign and a space
      if (poly->coeffs[i]) size += (unsigned long)ceil(log10(poly->coeffs[i])) + 2;
      else size += 3;
   }

   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "%ld  %ld  ", poly->length, poly->p);
   for (unsigned long i = 0; i < poly->length; i++)
   {
      ptr += sprintf(ptr, "%ld ", poly->coeffs[i]);
   }
   
   ptr--;
   *ptr = 0;
   
   return buf;
}


/*
   Convert a zmod_poly to a string and write it to the file f.
*/

void zmod_poly_fprint(zmod_poly_t poly, FILE* f)
{
   char* s = zmod_poly_to_string(poly);
   fputs(s, f);
   free(s);
}


/*
   Output the string representation of zmod_poly to stdout
*/

void zmod_poly_print(zmod_poly_t poly)
{
   zmod_poly_fprint(poly, stdout);
}


/*
   Create a zmod_poly from a string representation in file f
*/

int zmod_poly_fread(zmod_poly_t poly, FILE* f)
{
   // read poly length and mod
   unsigned long length, p;
   
   if (!fscanf(f, "%ld %ld", &length, &p))
      return 0;

   poly->length = 0;
   poly->p = p;
   poly->p_inv = z_precompute_inverse(p);
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) <= 32) poly->p32_inv = z_precompute_inverse32(p);
#endif
   
   zmod_poly_fit_length(poly, length);

   // read coefficients
   for (unsigned long i = 0; i < length; i++)
   {
      if (!fscanf(f, "%ld", &poly->coeffs[i]))
         return 0;
      poly->length++;
   }

   __zmod_poly_normalise(poly);
   
   return 1;
}


/*
   Create a zmod_poly from stdin
*/

int zmod_poly_read(zmod_poly_t poly)
{
   return zmod_poly_fread(poly, stdin);
}


/****************************************************************************

   Length and degree

****************************************************************************/


void __zmod_poly_normalise(zmod_poly_t poly)
{
   while (poly->length && (poly->coeffs[poly->length-1] == 0L))
      poly->length--;
}


int __zmod_poly_normalised(zmod_poly_t poly)
{
   return (poly->length == 0) || (poly->coeffs[poly->length-1] != 0L);
}


void zmod_poly_truncate(zmod_poly_t poly, unsigned long length)
{
   // inplace truncation

   if (length < poly->length)
         poly->length = length;
   
   __zmod_poly_normalise(poly);
}



/****************************************************************************

   Assignment

****************************************************************************/


void _zmod_poly_set(zmod_poly_t res, zmod_poly_t poly)
{
   if (res == poly)
      return;

   for (unsigned long i = 0; i < poly->length; i++)
      res->coeffs[i] = poly->coeffs[i];
      
   res->length = poly->length;
   
   res->p = poly->p;
   res->p_inv = poly->p_inv;
#if FLINT_BITS == 64
   res->p32_inv = poly->p32_inv;
#endif

}

void zmod_poly_set(zmod_poly_t res, zmod_poly_t poly)
{
   if (res == poly)
      return;

   zmod_poly_fit_length(res, poly->length);
   
   _zmod_poly_set(res, poly);
}


/****************************************************************************

   Comparison

****************************************************************************/


int zmod_poly_equal(zmod_poly_t poly1, zmod_poly_t poly2)
{
   if (poly1->p != poly2->p)
      return 0;
   
   if (poly1->length != poly2->length)
      return 0;

   for (unsigned long i = 0; i < poly1->length; i++)
      if (poly1->coeffs[i] != poly2->coeffs[i])
         return 0;

   return 1;
}

/****************************************************************************

   Reversal

****************************************************************************/


/* 
   Sets output to the reverse of input (i.e. reverse the order of the coefficients)
   assuming input to be a polynomial with _length_ coefficients (it may have a length
   that is less than _length_).
*/ 

void _zmod_poly_reverse(zmod_poly_t output, zmod_poly_t input, unsigned long length)
{
   long i;
   
   if (input != output)
   {
      for (i = 0; i < FLINT_MIN(length, input->length); i++)
      {
         output->coeffs[length - i - 1] = input->coeffs[i];
      }
      for ( ; i < length; i++)
      {
         output->coeffs[length - i - 1] = 0L;
      }
      output->length = length;
      __zmod_poly_normalise(output);
   } else
   {
      unsigned long temp;
      
      for (i = 0; i < length/2; i++)
      {
         if (i < input->length)
         {
            temp = input->coeffs[i];
         } else
         {
            temp = 0L;            
         }
         if (length - i - 1 < input->length)
         {
            input->coeffs[i] = input->coeffs[length - i - 1];
         } else
         {
            input->coeffs[i] = 0L;
         }
         input->coeffs[length - i - 1] =  temp;
      }
      if ((length & 1) && (i >= input->length)) input->coeffs[i] = 0L;

      output->length = length;
      __zmod_poly_normalise(output);
   }
}

void zmod_poly_reverse(zmod_poly_t output, zmod_poly_t input, unsigned long length)
{
   zmod_poly_fit_length(output, length);
   
   _zmod_poly_reverse(output, input, length);
}

/****************************************************************************

   Monic polys

****************************************************************************/

void zmod_poly_make_monic(zmod_poly_t output, zmod_poly_t pol)
{
   if (!pol->length) 
   {
      output->length = 0;
      return;
   }

   unsigned long lead_inv = pol->coeffs[pol->length-1];

   if (lead_inv == 1L) 
   {
      zmod_poly_set(output, pol);
      return;
   }

   lead_inv = z_invert(lead_inv, pol->p);

   zmod_poly_scalar_mul(output, pol, lead_inv);
}

/****************************************************************************

   Addition/subtraction

****************************************************************************/


void zmod_poly_add(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   // rearrange parameters to make poly1 no longer than poly2
   if (poly1->length > poly2->length)
      SWAP_ZMOD_POLY_PTRS(poly1, poly2);
      
   zmod_poly_fit_length(res, poly2->length);

   unsigned long i, neg1;
   /* The following standard technique was found in David Harvey's zn_poly */
   
   for (i = 0; i < poly1->length; i++)
   {
      res->coeffs[i] = z_addmod(poly1->coeffs[i], poly2->coeffs[i], poly1->p);
   }

   for (; i < poly2->length; i++)
      res->coeffs[i] = poly2->coeffs[i];

   res->length = poly2->length;
   __zmod_poly_normalise(res);
}

void _zmod_poly_add_without_mod(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   // rearrange parameters to make poly1 no longer than poly2
   if (poly1->length > poly2->length)
      SWAP_ZMOD_POLY_PTRS(poly1, poly2);
      
   unsigned long i, neg1;
   
   for (i = 0; i < poly1->length; i++)
   {
      res->coeffs[i] = poly1->coeffs[i] + poly2->coeffs[i];
   }

   for (; i < poly2->length; i++)
      res->coeffs[i] = poly2->coeffs[i];

   res->length = poly2->length;
   __zmod_poly_normalise(res);
}

void _zmod_poly_sub(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
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
      SWAP_ZMOD_POLY_PTRS(poly1, poly2);
   }
      
   unsigned long i;
   
   if (swapped)
   {
      for (i = 0; i < poly1->length; i++)
      {
         res->coeffs[i] = z_submod(poly2->coeffs[i], poly1->coeffs[i], poly2->p);
      }
         
      for (; i < poly2->length; i++)
         res->coeffs[i] = poly2->coeffs[i];
   }
   else
   {
      for (i = 0; i < poly1->length; i++)
      {
         res->coeffs[i] = z_submod(poly1->coeffs[i], poly2->coeffs[i], poly2->p);
      }
         
      for (; i < poly2->length; i++)
      {   
         res->coeffs[i] = poly2->p - poly2->coeffs[i];
         if (res->coeffs[i] == poly2->p) res->coeffs[i] = 0;
      }
   }

   res->length = poly2->length;
   __zmod_poly_normalise(res);
}

void zmod_poly_sub(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   if (poly1 == poly2)
   {
      // equal operands
      res->length = 0;
      return;
   }

   // rearrange parameters to make poly1 no longer than poly2
   if (poly1->length > poly2->length)
   {
      zmod_poly_fit_length(res, poly1->length);
   } else zmod_poly_fit_length(res, poly2->length);
   
   _zmod_poly_sub(res, poly1, poly2);
}
 
void zmod_poly_neg(zmod_poly_t res, zmod_poly_t poly)
{
   zmod_poly_fit_length(res, poly->length);

   for (unsigned long i = 0; i < poly->length; i++)
   {
      if (poly->coeffs[i]) res->coeffs[i] = poly->p - poly->coeffs[i];
      else res->coeffs[i] = 0L;
   }
   
   res->length = poly->length;
}



/****************************************************************************

   Shifting

****************************************************************************/


void zmod_poly_left_shift(zmod_poly_t res, zmod_poly_t poly, unsigned long k)
{
   zmod_poly_fit_length(res, poly->length + k);

   unsigned long temp;

   if (poly == res)
   {
      // inplace; just shift the coeffs over
      for (long i = poly->length - 1; i >= 0; i--)
      {
         poly->coeffs[i+k] = poly->coeffs[i];
      }
      
      for (unsigned long i = 0; i < k; i++)
         poly->coeffs[i] = 0L;
   }
   else
   {
      // not inplace; need to copy data
      for (unsigned long i = 0; i < k; i++)
         res->coeffs[i] = 0L;
      
      for (unsigned long i = 0; i < poly->length; i++)
         res->coeffs[i + k] = poly->coeffs[i];
         
      res->p = poly->p;
      res->p_inv = poly->p_inv;
#if FLINT_BITS == 64
      res->p32_inv = poly->p32_inv;
#endif

   }
   
   res->length = poly->length + k;
}


void zmod_poly_right_shift(zmod_poly_t res, zmod_poly_t poly, unsigned long k)
{
   if (k >= poly->length)
   {
      // shift all coefficients off the end
      res->length = 0;
      res->p = poly->p;
      res->p_inv = poly->p_inv;
#if FLINT_BITS == 64
      res->p32_inv = poly->p32_inv;
#endif

      return;
   }

   if (poly == res)
   {
      // inplace; just shift the mpz_t's over

      for (unsigned long i = k; i < poly->length; i++)
         poly->coeffs[i - k] = poly->coeffs[i];
   }
   else
   {
      // not inplace; need to copy data
      zmod_poly_fit_length(res, poly->length - k);

      for (unsigned long i = k; i < poly->length; i++)
         res->coeffs[i - k] = poly->coeffs[i];
         
      res->p = poly->p;
      res->p_inv = poly->p_inv;
#if FLINT_BITS == 64
      res->p32_inv = poly->p32_inv;
#endif

   }
   
   res->length = poly->length - k;
}



/*******************************************************************************

   Polynomial multiplication

********************************************************************************/


void zmod_poly_mul(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   if (poly1 == poly2)
   {
      zmod_poly_sqr(res, poly1);
      return;
   }
   
   if (poly1->length + poly2->length <= 6)
   {
      zmod_poly_mul_classical(res, poly1, poly2);
      return;
   }
   
   unsigned long bits = FLINT_BIT_COUNT(poly1->p);
   if ((bits <= 32) && (poly1->length + poly2->length <= 8))
   {
      zmod_poly_mul_classical(res, poly1, poly2);
      return;
   }

   zmod_poly_mul_KS(res, poly1, poly2, 0); 
}


void zmod_poly_sqr(zmod_poly_t res, zmod_poly_t poly)
{
   if (poly->length <= 4)
   {
      zmod_poly_sqr_classical(res, poly);
      return;
   }
   
   unsigned long bits = FLINT_BIT_COUNT(poly->p);
   if ((bits >= 32) && (bits <= 50) && (poly->length <= 10))
   {
      zmod_poly_sqr_classical(res, poly);
      return;
   }

   zmod_poly_mul_KS(res, poly, poly, 0); 
}


/*
 This is just like zmod_poly_mul_classical(), with the following restrictions:

  * assumes res does not alias poly1 and poly2
  * res->alloc >= poly1->length + poly2->length - 1
     (i.e. output has enough room for product)
*/

void _zmod_poly_mul_classical(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   FLINT_ASSERT(res != poly1);
   FLINT_ASSERT(res != poly2);
   
   if (!poly1->length || !poly2->length)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }

   res->length = poly1->length + poly2->length - 1;
   res->p = poly1->p;
   res->p_inv = poly1->p_inv;
#if FLINT_BITS == 64
   res->p32_inv = poly1->p32_inv;
#endif
   
   unsigned long length;
   
   if (poly1->length <= poly2->length)
   {
      length = poly1->length;
   }
   else
   {
      length = poly2->length;
   }
   
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   
   unsigned long bits = (FLINT_BIT_COUNT(poly1->p)<<1) + log_length;

   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = 0;

   if(bits < FLINT_BITS)
   {
      // the numbers of bits in the output of each coeff will be less than FLINT_BITS
      // so don't need to mod to stay in the single limb, hence can leave this for the
      // end...
      __zmod_poly_mul_classical_mod_last(res, poly1, poly2, bits);
   }
   else
   {
      bits = zmod_poly_bits(poly1) + zmod_poly_bits(poly2) + log_length;
      if (bits < FLINT_BITS)
      {
         __zmod_poly_mul_classical_mod_last(res, poly1, poly2, bits);
      }
      else
      {
         __zmod_poly_mul_classical_mod_throughout(res, poly1, poly2, bits);
      }
   }
      
   __zmod_poly_normalise(res);
}

/*
   Actually computes the classical multiplication, only applying mod at the end
   of the computations.
*/

void __zmod_poly_mul_classical_mod_last(zmod_poly_t res, zmod_poly_t poly1, 
                                             zmod_poly_t poly2, unsigned long bits)
{
   for (unsigned long i = 0; i < poly1->length; i++)
      for (unsigned long j = 0; j < poly2->length; j++)
         res->coeffs[i+j] = res->coeffs[i+j] + poly1->coeffs[i] * poly2->coeffs[j];
         
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   { 
      for (unsigned long i = 0; i < res->length; i++)
         res->coeffs[i] = z_mod_precomp(res->coeffs[i], res->p, res->p_inv);
   } else 
   { 
#endif
      for (unsigned long i = 0; i < res->length; i++)
         res->coeffs[i] = z_mod2_precomp(res->coeffs[i], res->p, res->p_inv);
#if FLINT_BITS == 64
   } 
#endif
}


/*
   Computes the classical multiplication, applying mods at each step.
*/

void __zmod_poly_mul_classical_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, 
                                            zmod_poly_t poly2, unsigned long bits)
{
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   {
      for (unsigned long i = 0; i < poly1->length; i++)
         for (unsigned long j = 0; j < poly2->length; j++)
            res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p);
   } else
   {
#endif
      for (unsigned long i = 0; i < poly1->length; i++)
         for (unsigned long j = 0; j < poly2->length; j++)
            res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod2_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p);
#if FLINT_BITS == 64
   }
#endif
}


void zmod_poly_mul_classical(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
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
      zmod_poly_sqr_classical(res, poly1);
      return;
   }

   unsigned long length = poly1->length + poly2->length - 1;

   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, poly1->p, length);
      _zmod_poly_mul_classical(temp, poly1, poly2);
      zmod_poly_swap(temp, res);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_fit_length(res, length);
      _zmod_poly_mul_classical(res, poly1, poly2);
   }
}


/*
 This is just like zmod_poly_sqr_classical(), with the following restrictions:

  * assumes res does not alias poly
  * res->alloc >= 2*poly->length - 1  (i.e. output has enough room for product)
*/
void _zmod_poly_sqr_classical(zmod_poly_t res, zmod_poly_t poly)
{
   FLINT_ASSERT(res != poly);
   FLINT_ASSERT(poly->length);

   if (!poly->length)
   {
      // input is zero
      res->length = 0;
      return;
   }

   res->length = 2*poly->length - 1;
   res->p = poly->p;
   res->p_inv = poly->p_inv;
#if FLINT_BITS == 64
   res->p32_inv = poly->p32_inv;
#endif
   FLINT_ASSERT(res->alloc >= res->length);
   
   unsigned long bits = FLINT_BIT_COUNT(poly->p);

   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = 0;

   // off-diagonal products
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   {
      for (unsigned long i = 1; i < poly->length; i++)
         for (unsigned long j = 0; j < i; j++)
            res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod_precomp(poly->coeffs[i], poly->coeffs[j], poly->p, poly->p_inv), poly->p);
   } else
   {
#endif
      for (unsigned long i = 1; i < poly->length; i++)
         for (unsigned long j = 0; j < i; j++)
            res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod2_precomp(poly->coeffs[i], poly->coeffs[j], poly->p, poly->p_inv), poly->p);
#if FLINT_BITS == 64
   }
#endif
   
   // double the off-diagonal products
   for (unsigned long i = 1; i < res->length - 1; i++)
      res->coeffs[i] = z_addmod(res->coeffs[i], res->coeffs[i], poly->p);

   // add in diagonal products
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   {
      for (unsigned long i = 0; i < poly->length; i++)
         res->coeffs[2*i] = z_addmod(res->coeffs[2*i], z_mulmod_precomp(poly->coeffs[i], poly->coeffs[i], poly->p, poly->p_inv), poly->p);
   } else
   {
#endif
      for (unsigned long i = 0; i < poly->length; i++)
         res->coeffs[2*i] = z_addmod(res->coeffs[2*i], z_mulmod2_precomp(poly->coeffs[i], poly->coeffs[i], poly->p, poly->p_inv), poly->p);
#if FLINT_BITS == 64
   }
#endif
      
   __zmod_poly_normalise(res);
}


void zmod_poly_sqr_classical(zmod_poly_t res, zmod_poly_t poly)
{
   if (!poly->length)
   {
      // input is zero
      res->length = 0;
      return;
   }

   unsigned long length = 2*poly->length - 1;

   if (res == poly)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, poly->p, length);
      _zmod_poly_sqr_classical(temp, poly);
      zmod_poly_swap(temp, res);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace

      // allocate more coefficients if necessary
      zmod_poly_fit_length(res, length);
      _zmod_poly_sqr_classical(res, poly);
   }
}

//=======================================================================

/*
 This is just like zmod_poly_mul_classical_trunc(), with the following restrictions:

  * assumes res does not alias poly1 and poly2
  * res->alloc >= MIN(trunc, poly1->length + poly2->length - 1)
     (i.e. output has enough room for truncated product)
*/

void _zmod_poly_mul_classical_trunc(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
{
   FLINT_ASSERT(res != poly1);
   FLINT_ASSERT(res != poly2);
   
   if (!poly1->length || !poly2->length || !trunc)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }
   
   if (trunc >= poly1->length + poly2->length - 1)
   {
      // there's no truncating to be done
      _zmod_poly_mul_classical(res, poly1, poly2);
      return;
   }

   res->length = trunc;
   res->p = poly1->p;
   res->p_inv = poly1->p_inv;
#if FLINT_BITS == 64
   res->p32_inv = poly1->p32_inv;
#endif
   
   unsigned long length;
   
   if (poly1->length <= poly2->length)
   {
      length = poly1->length;
   }
   else
   {
      length = poly2->length;
   }
   
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   
   unsigned long bits = (FLINT_BIT_COUNT(poly1->p)<<1) + log_length;

   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = 0;

   if(bits < FLINT_BITS)
   {
      // the numbers of bits in the output of each coeff will be less than FLINT_BITS
      // so don't need to mod to stay in the single limb, hence can leave this for the
      // end...
      __zmod_poly_mul_classical_trunc_mod_last(res, poly1, poly2, bits, trunc);
   }
   else
   {
      bits = zmod_poly_bits(poly1) + zmod_poly_bits(poly2) + log_length;
      if (bits < FLINT_BITS)
      {
         __zmod_poly_mul_classical_trunc_mod_last(res, poly1, poly2, bits, trunc);
      }
      else
      {
         __zmod_poly_mul_classical_trunc_mod_throughout(res, poly1, poly2, bits, trunc);
      }
   }
   
   __zmod_poly_normalise(res);
}

/*
   Actually computes the truncated classical multiplication, only applying mod at the end
   of the computations.
   
   Assumes neither poly length is zero and trunc is not zero
   Assumes res does not alias poly1 or poly2
   Assumes trunc < poly1->length + poly2->length - 1
*/

void __zmod_poly_mul_classical_trunc_mod_last(zmod_poly_t res, zmod_poly_t poly1, 
                                             zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
{
   for (unsigned long i = 0; i < poly1->length; i++)
      for (unsigned long j = 0; j < poly2->length; j++)
         if (i + j < trunc) 
            res->coeffs[i+j] = res->coeffs[i+j] + poly1->coeffs[i] * poly2->coeffs[j];
         
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   { 
      for (unsigned long i = 0; i < trunc; i++)
         res->coeffs[i] = z_mod_precomp(res->coeffs[i], res->p, res->p_inv);
   } else 
   { 
#endif
      for (unsigned long i = 0; i < trunc; i++)
         res->coeffs[i] = z_mod2_precomp(res->coeffs[i], res->p, res->p_inv);
#if FLINT_BITS == 64
   } 
#endif
}


/*
   Computes the classical multiplication, applying mods at each step.
   
   Assumes neither poly length is zero and trunc is not zero
   Assumes res does not alias poly1 or poly2
   Assumes trunc < poly1->length + poly2->length - 1

*/

void __zmod_poly_mul_classical_trunc_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, 
                                            zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
{
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   {
      for (unsigned long i = 0; i < poly1->length; i++)
         for (unsigned long j = 0; j < poly2->length; j++)
            if (i + j < trunc)
               res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p);
   } else
   {
#endif
      for (unsigned long i = 0; i < poly1->length; i++)
         for (unsigned long j = 0; j < poly2->length; j++)
            if (i + j < trunc)
               res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod2_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p);
#if FLINT_BITS == 64
   }
#endif
}


void zmod_poly_mul_classical_trunc(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
{   
   if (!poly1->length || !poly2->length || !trunc)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }

   /*if (poly1 == poly2)
   {
      // polys are identical, so call specialised truncated squaring routine
      zmod_poly_sqr_classical_trunc(res, poly1. trunc);
      return;
   }*/

   unsigned long length = poly1->length + poly2->length - 1;

   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, poly1->p, FLINT_MIN(length, trunc));
      _zmod_poly_mul_classical_trunc(temp, poly1, poly2, trunc);
      zmod_poly_swap(temp, res);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_fit_length(res, FLINT_MIN(length, trunc));
      _zmod_poly_mul_classical_trunc(res, poly1, poly2, trunc);
   }
}

//===================================================================================

void _zmod_poly_mul_classical_trunc_left(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
{
   FLINT_ASSERT(res != poly1);
   FLINT_ASSERT(res != poly2);
   
   if (!poly1->length || !poly2->length || (trunc >= poly1->length + poly2->length - 1))
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }
   
   if (trunc == 0)
   {
      _zmod_poly_mul_classical(res, poly1, poly2);
   }

   res->length = poly1->length + poly2->length - 1;
   res->p = poly1->p;
   res->p_inv = poly1->p_inv;
#if FLINT_BITS == 64
   res->p32_inv = poly1->p32_inv;
#endif
   
   unsigned long length;
   
   if (poly1->length <= poly2->length)
   {
      length = poly1->length;
   }
   else
   {
      length = poly2->length;
   }
   
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   
   unsigned long bits = (FLINT_BIT_COUNT(poly1->p)<<1) + log_length;

   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = 0;

   if(bits < FLINT_BITS)
   {
      // the numbers of bits in the output of each coeff will be less than FLINT_BITS
      // so don't need to mod to stay in the single limb, hence can leave this for the
      // end...
      __zmod_poly_mul_classical_trunc_left_mod_last(res, poly1, poly2, bits, trunc);
   }
   else
   {
      bits = zmod_poly_bits(poly1) + zmod_poly_bits(poly2) + log_length;
      if (bits < FLINT_BITS)
      {
         __zmod_poly_mul_classical_trunc_left_mod_last(res, poly1, poly2, bits, trunc);
      }
      else
      {
         __zmod_poly_mul_classical_trunc_left_mod_throughout(res, poly1, poly2, bits, trunc);
      }
   }
      
   __zmod_poly_normalise(res);
}

/*
   Actually computes the classical multiplication, only applying mod at the end
   of the computations.
*/

void __zmod_poly_mul_classical_trunc_left_mod_last(zmod_poly_t res, zmod_poly_t poly1, 
                                             zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
{
   for (unsigned long i = 0; i < poly1->length; i++)
      for (unsigned long j = 0; j < poly2->length; j++)
         if (i + j >= trunc)
            res->coeffs[i+j] = res->coeffs[i+j] + poly1->coeffs[i] * poly2->coeffs[j];
         
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   { 
      for (unsigned long i = trunc; i < res->length; i++)
         res->coeffs[i] = z_mod_precomp(res->coeffs[i], res->p, res->p_inv);
   } else 
   { 
#endif
      for (unsigned long i = trunc; i < res->length; i++)
         res->coeffs[i] = z_mod2_precomp(res->coeffs[i], res->p, res->p_inv);
#if FLINT_BITS == 64
   } 
#endif
}


/*
   Computes the classical multiplication, applying mods at each step.
*/

void __zmod_poly_mul_classical_trunc_left_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, 
                                            zmod_poly_t poly2, unsigned long bits, unsigned long trunc)
{
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   {
      for (unsigned long i = 0; i < poly1->length; i++)
         for (unsigned long j = 0; j < poly2->length; j++)
            if (i + j >= trunc)
               res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p);
   } else
   {
#endif
      for (unsigned long i = 0; i < poly1->length; i++)
         for (unsigned long j = 0; j < poly2->length; j++)
            if (i + j >= trunc)
               res->coeffs[i+j] = z_addmod(res->coeffs[i+j], z_mulmod2_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p);
#if FLINT_BITS == 64
   }
#endif
}


void zmod_poly_mul_classical_trunc_left(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
{   
   if (!poly1->length || !poly2->length)
   {
      // one of the polys is zero
      res->length = 0;
      return;
   }

   /*if (poly1 == poly2)
   {
      // polys are identical, so call specialised squaring routine
      zmod_poly_sqr_classical_trunc_left(res, poly1, trunc);
      return;
   }*/

   unsigned long length = poly1->length + poly2->length - 1;

   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, poly1->p, length);
      _zmod_poly_mul_classical_trunc_left(temp, poly1, poly2, trunc);
      zmod_poly_swap(temp, res);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_fit_length(res, length);
      _zmod_poly_mul_classical_trunc_left(res, poly1, poly2, trunc);
   }
}

//============================================================================================

/*
   Debugging function
*/

void print_var(char *name, unsigned long value)
{
   printf("%s = %d\n", name, value);
}

void zmod_poly_mul_KS(zmod_poly_t output, zmod_poly_p input1, zmod_poly_p input2, unsigned long bits_input)
{ 
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }
   
   unsigned long length = length1 + length2 - 1;
   
   zmod_poly_fit_length(output, length);
   
   if (output == input1 || output == input2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, input1->p, length);
      _zmod_poly_mul_KS(temp, input1, input2, bits_input);
      zmod_poly_swap(temp, output);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_fit_length(output, length);
      _zmod_poly_mul_KS(output, input1, input2, bits_input);
   }
} 
      
void _zmod_poly_mul_KS(zmod_poly_t output, zmod_poly_p input1, zmod_poly_p input2, unsigned long bits_input)
{   
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }

   unsigned long final_length = length1 + length2 - 1;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP_ZMOD_POLY_PTRS(input1, input2);
   }
      
   unsigned long bits1, bits2;
   
   bits1 = zmod_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : zmod_poly_bits(input2);
   
   unsigned long length = length2;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = bits1 + bits2 + log_length;
   
   if (bits_input) 
   {
      bits = bits_input;
   }
   
   mp_limb_t *mpn1, *mpn2, *res;

   unsigned long limbs1, limbs2;

   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = FLINT_MAX((long)((length2 * bits-1) / FLINT_BITS + 1), 0L);
      
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
   mpn2 = (input1 == input2) ? mpn1 : (mp_limb_t*) flint_stack_alloc(limbs2);

   _zmod_poly_bit_pack_mpn(mpn1, input1, bits, length1);
   
   if(input1 != input2)
      _zmod_poly_bit_pack_mpn(mpn2, input2, bits, length2);
   
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   res[limbs1+limbs2-1] = 0L;
   
   if (input1 != input2) F_mpn_mul(res, mpn1, limbs1, mpn2, limbs2);
   else F_mpn_mul(res, mpn1, limbs1, mpn1, limbs1);
   
   _zmod_poly_bit_unpack_mpn(output, res, length1 + length2 - 1, bits); 
   
   flint_stack_release();
   flint_stack_release();
   if(input1 != input2)
      flint_stack_release();
  
   output->length = final_length;

   /* The modulus may not be prime, so normalisation may be necessary */
   __zmod_poly_normalise(output);
}

//==========================================================================

void zmod_poly_mul_KS_trunc(zmod_poly_t output, zmod_poly_p input1, zmod_poly_p input2, unsigned long bits_input, unsigned long trunc)
{ 
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   if ((length1 == 0) || (length2 == 0) || (trunc == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }
   
   unsigned long length = length1 + length2 - 1;
   
   if (output == input1 || output == input2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, input1->p, FLINT_MIN(length, trunc));
      _zmod_poly_mul_KS_trunc(temp, input1, input2, bits_input, trunc);
      zmod_poly_swap(temp, output);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_fit_length(output, FLINT_MIN(length, trunc));
      _zmod_poly_mul_KS_trunc(output, input1, input2, bits_input, trunc);
   }
} 
      
void _zmod_poly_mul_KS_trunc(zmod_poly_t output, zmod_poly_p input1, zmod_poly_p input2, unsigned long bits_input, unsigned long trunc)
{   
   if ((input1->length == 0) || (input2->length == 0) || (trunc == 0))
   {
      output->length = 0;
      return;
   }

   unsigned long length1 = FLINT_MIN(input1->length, trunc);
   unsigned long length2 = FLINT_MIN(input2->length, trunc);
   
   while ((length1) && (input1->coeffs[length1-1] == 0)) length1--;
   while ((length2) && (input2->coeffs[length2-1] == 0)) length2--;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }
      
   unsigned long length = length1 + length2 - 1;
     
   if (trunc > length) trunc = length;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP_ZMOD_POLY_PTRS(input1, input2);
   }
         
   unsigned long bits1, bits2;
   
   bits1 = zmod_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : zmod_poly_bits(input2);
   
   unsigned long length_short = length2;
   unsigned log_length = 0;
   while ((1L<<log_length) < length_short) log_length++;
   unsigned long bits = bits1 + bits2 + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }
   
   mp_limb_t *mpn1, *mpn2, *res;

   unsigned long limbs1, limbs2;
   
   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = FLINT_MAX((long)((length2 * bits-1) / FLINT_BITS + 1), 0L);
      
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
   mpn2 = (input1 == input2) ? mpn1 : (mp_limb_t*) flint_stack_alloc(limbs2);
         
   _zmod_poly_bit_pack_mpn(mpn1, input1, bits, length1);
   
   if(input1 != input2)
      _zmod_poly_bit_pack_mpn(mpn2, input2, bits, length2);
         
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   res[limbs1+limbs2-1] = 0L;
   
   unsigned long output_length = FLINT_MIN(length1 + length2 - 1, trunc);
   
   if (input1 != input2) F_mpn_mul_trunc(res, mpn1, limbs1, mpn2, limbs2, (output_length*bits-1)/FLINT_BITS+1);
   else F_mpn_mul_trunc(res, mpn1, limbs1, mpn1, limbs1, (output_length*bits-1)/FLINT_BITS+1);
         
   _zmod_poly_bit_unpack_mpn(output, res, output_length, bits); 
   flint_stack_release(); //release res
   flint_stack_release(); //release mpn1 and mpn2
   if(input1 != input2)
      flint_stack_release();
  
   output->length = output_length;
   
   /* The modulus may not be prime, so normalisation may be necessary */
   __zmod_poly_normalise(output);
    
}

void zmod_poly_mul_trunc_n_precomp(zmod_poly_t output, zmod_poly_p input1, zmod_poly_precomp_t pre, unsigned long trunc)
{ 
   if ((input1->length == 0) || (pre->length2 == 0) || (trunc == 0))
   {
      output->length = 0;
      return;
   }
   
   unsigned long length = input1->length + pre->length2 - 1;

   zmod_poly_fit_length(output, FLINT_MIN(length, trunc));
   _zmod_poly_mul_KS_trunc_precomp(output, input1, pre, 0, trunc);
}

/*
   Prepare for caching of FFT for a truncated multiplicaton of input2, with the given number of 
   bits per output coefficient (0 if this is to be computed automatically) where the output will
   be truncated to the given length
*/
   
void zmod_poly_mul_trunc_n_precomp_init(zmod_poly_precomp_t pre, zmod_poly_p input2, unsigned long bits_input, unsigned long trunc)
{
   unsigned long length2 = FLINT_MIN(input2->length, trunc);
   while ((length2) && (input2->coeffs[length2-1] == 0)) length2--;

   pre->length2 = length2;
   if ((length2 == 0) || (trunc == 0))
   {
      F_mpn_mul_precomp_init(pre->precomp, NULL, 0, 0);
      return;
   }
  
   unsigned log_length = 0;
   while ((1L<<log_length) < length2) log_length++;
   unsigned long bits = 2*FLINT_BIT_COUNT(input2->p) + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }

   unsigned long limbs1 = FLINT_MAX((long)((trunc * bits-1) / FLINT_BITS + 1), 0L);
   unsigned long limbs2 = FLINT_MAX((long)((length2 * bits-1) / FLINT_BITS + 1), 0L);
   pre->limbs2 = limbs2;

   mp_limb_t* mpn2 = (mp_limb_t*) flint_stack_alloc(limbs2);
         
   _zmod_poly_bit_pack_mpn(mpn2, input2, bits, length2);

   F_mpn_mul_precomp_init(pre->precomp, mpn2, limbs2, limbs1);

   flint_stack_release(); // release mpn2 
}

void zmod_poly_mul_precomp_init(zmod_poly_precomp_t pre, zmod_poly_t input2, unsigned long bits_input, unsigned long length1)
{
   unsigned long length2 = input2->length;
   while ((length2) && (input2->coeffs[length2-1] == 0)) length2--;

   pre->length2 = length2;
   if (length2 == 0)
   {
      F_mpn_mul_precomp_init(pre->precomp, NULL, 0, 0);
      return;
   }
  
   unsigned log_length = 0;
   while ((1L<<log_length) < length2) log_length++;
   unsigned long bits = 2*FLINT_BIT_COUNT(input2->p) + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }

   unsigned long limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   unsigned long limbs2 = FLINT_MAX((long)((length2 * bits-1) / FLINT_BITS + 1), 0L);
   pre->limbs2 = limbs2;

   mp_limb_t* mpn2 = (mp_limb_t*) flint_stack_alloc(limbs2);
         
   _zmod_poly_bit_pack_mpn(mpn2, input2, bits, length2);

   F_mpn_mul_precomp_init(pre->precomp, mpn2, limbs2, limbs1);

   flint_stack_release(); // release mpn2 
}

void zmod_poly_precomp_clear(zmod_poly_precomp_t pre)
{
   F_mpn_mul_precomp_clear(pre->precomp);
}

void _zmod_poly_mul_KS_trunc_precomp(zmod_poly_t output, zmod_poly_t input1, zmod_poly_precomp_t pre, unsigned long bits_input, unsigned long trunc)
{   
   unsigned long length1 = FLINT_MIN(input1->length, trunc);
   
   while ((length1) && (input1->coeffs[length1-1] == 0)) length1--;
   
   unsigned long length2 = pre->length2;
   
   unsigned long length = length1 + length2 - 1;
     
   if (trunc > length) trunc = length;
         
   if ((length1 == 0) || (length2 == 0) || (trunc == 0))
   {
      output->length = 0;
      return;
   }
   unsigned log_length = 0;
   while ((1L<<log_length) < length2) log_length++;
   unsigned long bits = 2*FLINT_BIT_COUNT(input1->p) + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }
   
   mp_limb_t *mpn1, *res;

   unsigned long limbs1, limbs2;
   
   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = pre->limbs2;
   
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
         
   _zmod_poly_bit_pack_mpn(mpn1, input1, bits, length1);
         
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   F_mpn_clear(res, limbs1+limbs2);
   
   unsigned long output_length = FLINT_MIN(trunc, length);
   F_mpn_mul_precomp_trunc(res, mpn1, limbs1, pre->precomp, (output_length*bits-1)/FLINT_BITS+1);
        
   _zmod_poly_bit_unpack_mpn(output, res, output_length, bits); 
   flint_stack_release(); //release res
   flint_stack_release(); //release mpn1
   
   output->length = output_length;
   
   /* The modulus may not be prime, so normalisation may be necessary */
   __zmod_poly_normalise(output);
    
}

void _zmod_poly_mul_KS_precomp(zmod_poly_t output, zmod_poly_t input1, zmod_poly_precomp_t pre, unsigned long bits_input)
{   
   unsigned long length1 = input1->length;
   
   while ((length1) && (input1->coeffs[length1-1] == 0)) length1--;
   
   unsigned long length2 = pre->length2;
   
   unsigned long length = length1 + length2 - 1;
     
   if ((length1 == 0) || (length2 == 0))
   {
      output->length = 0;
      return;
   }
   unsigned log_length = 0;
   while ((1L<<log_length) < length2) log_length++;
   unsigned long bits = 2*FLINT_BIT_COUNT(input1->p) + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }
   
   mp_limb_t *mpn1, *res;

   unsigned long limbs1, limbs2;
   
   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = pre->limbs2;
   
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
         
   _zmod_poly_bit_pack_mpn(mpn1, input1, bits, length1);
         
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   res[limbs1+limbs2-1] = 0L;
   
   F_mpn_mul_precomp(res, mpn1, limbs1, pre->precomp);
        
   zmod_poly_fit_length(output, length);
   _zmod_poly_bit_unpack_mpn(output, res, length, bits); 
   flint_stack_release(); //release res
   flint_stack_release(); //release mpn1
   
   output->length = length;
   /* The modulus may not be prime, so normalisation may be necessary */
   __zmod_poly_normalise(output);
    
}

#if USE_MIDDLE_PRODUCT
/*
   This function can be used with _zmod_poly_mul_KS_precomp_init (which actually
	precomputes an FFT for a full product, not a middle product). The saving is that
	only the middle terms of the product are actually computed (and the precomputed
	FFT can be reused)
	The terms that are computed are [trunc/2, trunc), thus if the product is 2n x n
	and trunc is set to 2n, then terms [n, 2n) of the product will be computed
*/

void _zmod_poly_mul_KS_middle_precomp(zmod_poly_t output, zmod_poly_p input1, zmod_poly_precomp_t pre, unsigned long bits_input, unsigned long trunc)
{   
   unsigned long length1 = FLINT_MIN(input1->length, trunc);
   
   while ((length1) && (input1->coeffs[length1-1] == 0)) length1--;
   
   unsigned long length2 = pre->length2;
   
   unsigned long length = length1 + length2 - 1;
     
   if (trunc > length) trunc = length;
         
   unsigned log_length = 0;
   while ((1L<<log_length) < length2) log_length++;
   unsigned long bits = 2*FLINT_BIT_COUNT(input1->p) + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }
   
   mp_limb_t *mpn1, *res;

   unsigned long limbs1, limbs2;
   
   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = pre->limbs2;
   
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
         
   _zmod_poly_bit_pack_mpn(mpn1, input1, bits, length1);
         
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   res[limbs1+limbs2-1] = 0L;
   
   unsigned long output_length = FLINT_MIN(length1 + length2 - 1, trunc);
   unsigned long start = trunc/2;
   
   __F_mpn_mul_middle_precomp(res, mpn1, limbs1, pre->precomp, (start*bits)/FLINT_BITS, (output_length*bits-1)/FLINT_BITS+1);
        
   _zmod_poly_bit_unpack_mpn(output, res, output_length, bits); 
   flint_stack_release(); //release res
   flint_stack_release(); //release mpn1
   
   output->length = output_length;
   
   /* The modulus may not be prime, so normalisation may be necessary */
   __zmod_poly_normalise(output);
    
}

/*
   Computes a middle product of two polynomials
	The terms that are computed are [trunc/2, trunc), thus if the product is 2n x n
	and trunc is set to 2n, then terms [n, 2n) of the product will be computed
*/

void _zmod_poly_mul_KS_middle(zmod_poly_t output, zmod_poly_p input1, zmod_poly_p input2, unsigned long bits_input, unsigned long trunc)
{   
   unsigned long length1 = FLINT_MIN(input1->length, trunc);
   unsigned long length2 = FLINT_MIN(input2->length, trunc);
   
   while ((length1) && (input1->coeffs[length1-1] == 0)) length1--;
   while ((length2) && (input2->coeffs[length2-1] == 0)) length2--;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }
      
   unsigned long length = length1 + length2 - 1;
     
   if (trunc > length) trunc = length;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP_ZMOD_POLY_PTRS(input1, input2);
   }
         
   unsigned long bits1, bits2;
   
   bits1 = zmod_poly_bits(input1);
   bits2 = zmod_poly_bits(input2);
   
   unsigned long length_short = length2;
   unsigned log_length = 0;
   while ((1L<<log_length) < length_short) log_length++;
   unsigned long bits = bits1 + bits2 + log_length;
      
   if (bits_input) 
   {
      bits = bits_input;
   }
   
   mp_limb_t *mpn1, *mpn2, *res;

   unsigned long limbs1, limbs2;
   
   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = FLINT_MAX((long)((length2 * bits-1) / FLINT_BITS + 1), 0L);
      
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
   mpn2 = (input1 == input2) ? mpn1 : (mp_limb_t*) flint_stack_alloc(limbs2);
         
   _zmod_poly_bit_pack_mpn(mpn1, input1, bits, length1);
   
   _zmod_poly_bit_pack_mpn(mpn2, input2, bits, length2);
         
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   res[limbs1+limbs2-1] = 0L;
   
   unsigned long output_length = FLINT_MIN(length1 + length2 - 1, trunc);
   unsigned long start = trunc/2;
   
   __F_mpn_mul_middle(res, mpn1, limbs1, mpn2, limbs2, (start*bits)/FLINT_BITS, (output_length*bits-1)/FLINT_BITS+1);
         
   _zmod_poly_bit_unpack_mpn(output, res, output_length, bits); 
   flint_stack_release(); //release res
   flint_stack_release(); //release mpn1 and mpn2
   flint_stack_release();
  
   for (unsigned long i = 0; i < start; i++)
      output->coeffs[i] = 0L;
   output->length = output_length;
   
   /* The modulus may not be prime, so normalisation may be necessary */
   __zmod_poly_normalise(output);
    
}

void zmod_poly_mul_KS_middle(zmod_poly_t output, zmod_poly_p input1, zmod_poly_p input2, unsigned long bits_input, unsigned long trunc)
{ 
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   if ((length1 == 0) || (length2 == 0) || (trunc == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }
   
   unsigned long length = length1 + length2 - 1;
   
   if (output == input1 || output == input2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, input1->p, FLINT_MIN(length, trunc));
      _zmod_poly_mul_KS_middle(temp, input1, input2, bits_input, trunc);
      zmod_poly_swap(temp, output);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_fit_length(output, FLINT_MIN(length, trunc));
      _zmod_poly_mul_KS_middle(output, input1, input2, bits_input, trunc);
   }
} 
#endif

/*******************************************************************************

   Bitpacking functions

********************************************************************************/


/*
   Determines the maximum number of bits used in the coefficients of poly
*/

unsigned long zmod_poly_bits(zmod_poly_t poly)
{
   unsigned long bits = 0;
   unsigned long mask = -1L;
   for(unsigned long i = 0; i < poly->length; i++)
   {
      if(poly->coeffs[i])
      {
         if(poly->coeffs[i] & mask)
         {
            bits = FLINT_BIT_COUNT(poly->coeffs[i]);
            if(bits == FLINT_BITS) break;
            else mask = -1L - ((1L<<bits)-1);
         }
      }
   }
   return bits;
}


/*
   Debugging function
   
   Prints an unsigned long in binary of length len
*/

void print_binary(unsigned long n, unsigned long len)
{
   while(n || len)
   {
      if(n % 2)
      {
         printf("1");
      }
      else
      {
         printf("0");
      }
      n /=2;
      len--;
   }
}

void print_binary2(unsigned long n, unsigned long len, unsigned long space_bit)
{
   while(n || len)
   {
      if(len == space_bit) printf(" ");
      if(n % 2)
      {
         printf("1");
      }
      else
      {
         printf("0");
      }
      n /=2;
      len--;
   }
}


/*
   Debugging function
   
   Prints a limb, in the format "name = limb".
*/

void print_limb(char *name, unsigned long limb)
{
   printf("%s = ", name);
   print_binary(limb, FLINT_BITS);
   printf("\n");
}


/*
   Packs the zmod_poly into an mpn, using `bits` bits for each coefficient
*/

void _zmod_poly_bit_pack_mpn(mp_limb_t * res, zmod_poly_t poly, unsigned long bits, unsigned long length)
{  
   unsigned long current_limb = 0;
   unsigned int current_bit = 0;
   
   unsigned long temp_lower;
   unsigned long temp_upper;
   
   unsigned long total_limbs = FLINT_MAX((long)(((length * bits - 1)>>FLINT_LG_BITS_PER_LIMB) + 1), 0L);
   
   res[0] = 0L;
   
   if (bits < FLINT_BITS)
   {
      unsigned long boundary_limit_bit = FLINT_BITS - bits;

      //printf("Packing polynomial ****************************************\n");
      //print_limb("res[0]", res[0]);

      for(unsigned long i = 0; i < length; i++)
      {
         if (current_bit > boundary_limit_bit)
         {
            // the coefficient will be added accross a limb boundary,
            // so need the lower and upper parts (lower for limb with
            // lower index).
            //printf("coeff won't fit in limb...\n");
            //print_limb("poly->coeffs[i]                ", poly->coeffs[i]);

            // the part of the coeff that will be in the current limb
            temp_lower = (poly->coeffs[i] << current_bit);
            //print_limb("temp_lower                     ", temp_lower);
            // the part of the coeff that will be in the next limb
            temp_upper = (poly->coeffs[i] >> (FLINT_BITS - current_bit));
            //print_limb("temp_upper                     ", temp_upper);
            //print_limb("res[current_limb]              ", res[current_limb]);
            res[current_limb] |= temp_lower;
            //print_limb("res[current_limb] |= temp_lower", res[current_limb]);
            current_limb++;
            res[current_limb] = temp_upper;
            //print_limb("res[current_limb+1]            ", res[current_limb]);
            current_bit = bits + current_bit - FLINT_BITS;
         }
         else
         {
            // the coefficient will fit in the current limb
            //printf("coeff will fit in limb...\n");
            temp_lower = poly->coeffs[i] << current_bit;
            //print_limb("poly->coeffs[i]                ", poly->coeffs[i]);
            //print_limb("temp_lower                     ", temp_lower);
            //print_limb("res[current_limb]              ", res[current_limb]);
            res[current_limb] |= temp_lower;
            //print_limb("res[current_limb] |= temp_lower", res[current_limb]);
            current_bit += bits;
         }

         if (current_bit >= FLINT_BITS)
         {
            current_limb++;
            if (current_limb < total_limbs) res[current_limb] = 0L;
            current_bit -= FLINT_BITS;
         }
      }
   }
   else if (bits == FLINT_BITS)
   {
      for (unsigned long i = 0; i < length; i++)
      {
         res[i] = poly->coeffs[i];
      }
   }
   else if (bits == 2*FLINT_BITS)
   {
      for (unsigned long i = 0; i < length; i++)
      {
         res[current_limb] = poly->coeffs[i];
         current_limb++;
         res[current_limb] = 0L;
         current_limb++;
      }
   }
   else if (bits < 2*FLINT_BITS)
   {
      //printf("Packing Coeffs in Poly =============================================");
      
      for(unsigned long i = 0; i < length; i++)
      {
         //PRINT_VAR(current_bit);
         // the coefficient will be added accross a limb boundary,
         // so need the lower and upper parts (lower for limb with
         // lower index).
         //printf("coeff won't fit in limb... HERE\n");
         //print_limb("poly->coeffs[i]                ", poly->coeffs[i]);

         // the part of the coeff that will be in the current limb
         temp_lower = poly->coeffs[i] << current_bit;
         //print_limb("temp_lower                     ", temp_lower);
         // the part of the coeff that will be in the next limb
         if (current_bit)
         {
            //print_var("current_bit", current_bit);
            temp_upper = poly->coeffs[i] >> (FLINT_BITS - current_bit);
         }
         else
         {
            temp_upper = 0L;
         }
         //print_limb("temp_upper                     ", temp_upper);
         //print_limb("res[current_limb]              ", res[current_limb]);
         res[current_limb] |= temp_lower;
         //print_limb("res[current_limb] |= temp_lower", res[current_limb]);
         current_limb++;
         res[current_limb] = temp_upper;
         //print_limb("res[current_limb+1]            ", res[current_limb]);
         current_bit += bits - FLINT_BITS;
         //PRINT_VAR(current_bit);
         
         if (current_bit >= FLINT_BITS)
         {
            //printf("GOT HERE ****************\n");
            current_bit -= FLINT_BITS;
            current_limb++;
            if (current_limb < total_limbs) res[current_limb] = 0L;
         }
      }
   } else // 2*FLINT_BITS < bits < 3*FLINT_BITS
   {      
      for(unsigned long i = 0; i < length; i++)
      {
         // the part of the coeff that will be in the current limb
         temp_lower = poly->coeffs[i] << current_bit;
         // the part of the coeff that will be in the next limb
         if (current_bit)
         {
            temp_upper = poly->coeffs[i] >> (FLINT_BITS - current_bit);
         }
         else
         {
            temp_upper = 0L;
         }
         res[current_limb] |= temp_lower;
         current_limb++;
         res[current_limb] = temp_upper;
         current_limb++;
         if (current_limb < total_limbs) res[current_limb] = 0L;
         current_bit += bits - 2*FLINT_BITS;
         
         if (current_bit >= FLINT_BITS)
         {
            current_bit -= FLINT_BITS;
            current_limb++;
            if (current_limb < total_limbs) res[current_limb] = 0L;
         }
      }
   }
}


/*
   Unpacks a zmod_poly of length `length` from an mpn `mpn` with coeffs packed in `bits` bits.
*/

void _zmod_poly_bit_unpack_mpn(zmod_poly_t res, mp_limb_t * mpn, unsigned long length, unsigned long bits)
{
   unsigned long i;
   
#if FLINT_BITS == 64
   if ((bits <= 32) && (FLINT_BIT_COUNT(res->p) <= 32))
   {
      unsigned long current_limb = 0;
      unsigned long current_bit = 0;

      unsigned long boundary_limit_bit = FLINT_BITS - bits;

      unsigned long temp_lower;
      unsigned long temp_upper;
      
      unsigned long norm;
      count_lead_zeros(norm, res->p);
      norm -= 32;
      
      unsigned long mask;
      mask = 1L;
      i = bits - 1;
      while(i)
      {
         mask <<= 1;
         mask |= 1L;
         i--;
      }
      
      for (i = 0; i < length; i++)
      {
          if (current_bit > boundary_limit_bit)
          {
             // the coeff will be across a limb boundary...
             //printf("coeff won't only be in current limb...\n");

             // temp lower contains the part in the current limb
             temp_lower = mpn[current_limb];
             //print_limb("mpn[current_limb]       ", temp_lower);

             // need (bits - (FLINT_BITS - current_bit)) bits 
             // from the LSB side of this limb to complete the coeff...
             current_limb++;
             //print_limb("mpn[current_limb+1]     ", mpn[current_limb]);
             // so shift them up, OR with the lower part and apply the mask
             temp_upper = mpn[current_limb] << (FLINT_BITS - current_bit);
             //print_limb("temp_upper              ", temp_upper);
             temp_upper |= temp_lower;
             //print_limb("temp_upper |= temp_lower", temp_upper);
             temp_upper &= mask;
             //print_limb("temp_upper &= mask      ", temp_upper);
             _zmod_poly_set_coeff_ui(res, i, z_mod32_precomp(temp_upper<<norm, res->p<<norm, res->p32_inv)>>norm);
             
             current_bit = bits + current_bit - FLINT_BITS;
             mpn[current_limb] = mpn[current_limb] >> current_bit;
             //print_limb("mpn[current_limb+1]     ", mpn[current_limb]);
          }
          else
          {
             // the coeff will fit in the current limb...
             //printf("coeff will be in current limb...\n");
             //print_limb("mpn[current_limb]       ", mpn[current_limb]);
             temp_lower = mpn[current_limb] & mask;
             //print_limb("temp_lower              ", temp_lower);
             // less than a limb in size, so must be smaller than an unsigned long...

             //zmod_poly_set_coeff_ui(res, i, temp_lower);
             _zmod_poly_set_coeff_ui(res, i, z_mod32_precomp(temp_lower<<norm, res->p<<norm, res->p32_inv)>>norm);
                
             mpn[current_limb] = mpn[current_limb] >> bits;
             //print_limb("mpn[current_limb]       ", mpn[current_limb]);
             current_bit += bits;
          }

          if(current_bit == FLINT_BITS)
          {
             current_bit = 0;
             current_limb++;
          }
      }
   }
   else 
#endif
   if (bits < FLINT_BITS)
   {
      unsigned long current_limb = 0;
      unsigned long current_bit = 0;

      unsigned long boundary_limit_bit = FLINT_BITS - bits;

      unsigned long temp_lower;
      unsigned long temp_upper;
      
      
      unsigned long mask;
      mask = 1L;
      i = bits - 1;
      while(i)
      {
         mask <<= 1;
         mask |= 1L;
         i--;
      }
      
      for (i = 0; i < length; i++)
      {
          if (current_bit > boundary_limit_bit)
          {
             // the coeff will be across a limb boundary...
             //printf("coeff won't only be in current limb...\n");

             // temp lower contains the part in the current limb
             temp_lower = mpn[current_limb];
             //print_limb("mpn[current_limb]       ", temp_lower);

             // need (bits - (FLINT_BITS - current_bit)) bits 
             // from the LSB side of this limb to complete the coeff...
             current_limb++;
             //print_limb("mpn[current_limb+1]     ", mpn[current_limb]);
             // so shift them up, OR with the lower part and apply the mask
             temp_upper = mpn[current_limb] << (FLINT_BITS - current_bit);
             //print_limb("temp_upper              ", temp_upper);
             temp_upper |= temp_lower;
             //print_limb("temp_upper |= temp_lower", temp_upper);
             temp_upper &= mask;
             //print_limb("temp_upper &= mask      ", temp_upper);
#if FLINT_BITS == 64
             if (bits <= FLINT_D_BITS)
                _zmod_poly_set_coeff_ui(res, i, z_mod_precomp(temp_upper, res->p, res->p_inv));
             else 
#endif
                _zmod_poly_set_coeff_ui(res, i, z_mod2_precomp(temp_upper, res->p, res->p_inv));
             
             current_bit = bits + current_bit - FLINT_BITS;
             mpn[current_limb] = mpn[current_limb] >> current_bit;
             //print_limb("mpn[current_limb+1]     ", mpn[current_limb]);
          }
          else
          {
             // the coeff will fit in the current limb...
             //printf("coeff will be in current limb...\n");
             //print_limb("mpn[current_limb]       ", mpn[current_limb]);
             temp_lower = mpn[current_limb] & mask;
             //print_limb("temp_lower              ", temp_lower);
             // less than a limb in size, so must be smaller than an unsigned long...

             //zmod_poly_set_coeff_ui(res, i, temp_lower);
#if FLINT_BITS == 64
             if (bits <= FLINT_D_BITS)
                _zmod_poly_set_coeff_ui(res, i, z_mod_precomp(temp_lower, res->p, res->p_inv));
             else 
#endif
                _zmod_poly_set_coeff_ui(res, i, z_mod2_precomp(temp_lower, res->p, res->p_inv));
                
             mpn[current_limb] = mpn[current_limb] >> bits;
             //print_limb("mpn[current_limb]       ", mpn[current_limb]);
             current_bit += bits;
          }

          if(current_bit == FLINT_BITS)
          {
             current_bit = 0;
             current_limb++;
          }
      }
   }
   else if (bits == FLINT_BITS)
   {
      for (i = 0; i < length; i++)
      {
         _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(0L, mpn[i], res->p, res->p_inv));
      }
   }
   else if (bits == 2*FLINT_BITS)
   {
      unsigned long current_limb = 0;
      for (i = 0; i < length; i++)
      {
         _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(mpn[current_limb+1], mpn[current_limb], res->p, res->p_inv));
         current_limb+=2;
      }
   }
   else if (bits < 2*FLINT_BITS) // FLINT_BITS < bits < 2*FLINT_BITS
   {
      unsigned long current_limb = 0;
      unsigned long current_bit = 0;

      unsigned long double_boundary_limit_bit = bits - FLINT_BITS;

      unsigned long temp_lower;
      unsigned long temp_upper;
      
      for (i = 0; i < length; i++)
      {
         if(current_bit == 0)
         {
            // printf("Coeff across one boundary... current_bit == 0\n");
            temp_lower = mpn[current_limb];
            // PRINT_LIMB(temp_lower);
            current_limb++;
            temp_upper = (mpn[current_limb] << (2*FLINT_BITS - bits)) >> (2*FLINT_BITS - bits);
            // PRINT_LIMB(temp_upper);
            _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
            mpn[current_limb] >>= (bits - FLINT_BITS);
            // PRINT_LIMB(mpn[current_limb]);
            current_bit = 2*FLINT_BITS - bits;
            // PRINT_VAR(current_bit);
         }
         else if (current_bit < double_boundary_limit_bit)
         {
            // printf("Coeff across two boundaries...\n");
            // the coeff will be across two limb boundaries...
            temp_lower = mpn[current_limb];
            // PRINT_LIMB(temp_lower);
            // PRINT_VAR(current_bit);
            current_limb++;
            // PRINT_LIMB(mpn[current_limb] << current_bit);
            temp_lower |= (mpn[current_limb] << current_bit);
            // PRINT_LIMB(temp_lower);
            // FLINT_BITS - current_bit != FLINT_BITS as current_bit > double_boundary_limit_bit
            temp_upper = mpn[current_limb] >> (FLINT_BITS - current_bit);
            // PRINT_LIMB(temp_upper);
            current_limb++;
            // PRINT_LIMB(mpn[current_limb]);
            temp_upper |= (mpn[current_limb] << current_bit);
            temp_upper <<= 2*FLINT_BITS - bits;
            temp_upper >>= 2*FLINT_BITS - bits;
            // PRINT_LIMB(temp_upper)
            _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
            mpn[current_limb] >>= (bits - current_bit - FLINT_BITS);
            current_bit = 2*FLINT_BITS + current_bit - bits;
            // PRINT_VAR(current_bit);
         }
         else 
         {
            // the coeff will be across one limb boundary...
            // printf("Coeff across one boundary...\n");
            temp_lower = mpn[current_limb] | (mpn[current_limb+1] << current_bit);
            // PRINT_LIMB(mpn[current_limb]);
            //  PRINT_LIMB(mpn[current_limb+1]);
            //  PRINT_LIMB(temp_lower);

            current_limb++;
            
            //PRINT_LIMB(temp_lower);

            temp_upper = (mpn[current_limb] << (FLINT_BITS + current_bit - bits)) >> (2*FLINT_BITS - bits);

            // PRINT_LIMB(temp_upper);

            _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
            mpn[current_limb] >>= (bits - current_bit);
            // PRINT_LIMB(mpn[current_limb]);
            current_bit = FLINT_BITS + current_bit - bits;
            // PRINT_VAR(current_bit);
            if(!current_bit) current_limb++;
         }

         if(current_bit == FLINT_BITS)
         {
            current_bit = 0; 
            // PRINT_VAR(current_bit);
            current_limb++;  
         }     
      }
   } else // 2*FLINT_BITS < bits < 3*FLINT_BITS
   {
      unsigned long current_limb = 0;
      unsigned long current_bit = 0;

      unsigned long double_boundary_limit_bit = bits - 2*FLINT_BITS;

      unsigned long temp_lower;
      unsigned long temp_upper;
      unsigned long temp_upper2;
      
      for (i = 0; i < length; i++)
      {
         if(current_bit == 0)
         {
            // printf("Coeff across two boundaries... current_bit == 0\n");
            temp_lower = mpn[current_limb+1];
            temp_upper = (mpn[current_limb+2] << (3*FLINT_BITS - bits)) >> (3*FLINT_BITS - bits);
            temp_upper = z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv);
            temp_lower = mpn[current_limb];
            current_limb+=2;
            _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
            mpn[current_limb] >>= (bits - 2*FLINT_BITS);
            current_bit = 3*FLINT_BITS - bits;
         }
         else if (current_bit >= double_boundary_limit_bit)
         {
            // the coeff will be across two limb boundaries...
            temp_lower = mpn[current_limb];
            current_limb++;
            temp_lower |= (mpn[current_limb] << current_bit);
            
            temp_upper = mpn[current_limb] >> (FLINT_BITS - current_bit);
            current_limb++;
            temp_upper |= (mpn[current_limb] << current_bit);
            temp_upper2 = mpn[current_limb] >> (FLINT_BITS - current_bit);
            temp_upper2 <<= 3*FLINT_BITS - bits;
            temp_upper2 >>= 3*FLINT_BITS - bits;
            temp_upper = z_ll_mod_precomp(temp_upper2, temp_upper, res->p, res->p_inv);
            _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
            mpn[current_limb] >>= (bits - current_bit - FLINT_BITS);
            current_bit = 2*FLINT_BITS + current_bit - bits;
            if (!current_bit) current_limb++;
         }
         else 
         {
            // the coeff will be across three limb boundaries...
            temp_lower = mpn[current_limb];
            current_limb++;
            temp_lower |= (mpn[current_limb] << current_bit);
            
            temp_upper = mpn[current_limb] >> (FLINT_BITS - current_bit);
            current_limb++;
            temp_upper |= (mpn[current_limb] << current_bit);
            temp_upper2 = mpn[current_limb] >> (FLINT_BITS - current_bit);
            current_limb++;
            temp_upper2 |= (mpn[current_limb] << current_bit);
            temp_upper2 <<= 3*FLINT_BITS - bits;
            temp_upper2 >>= 3*FLINT_BITS - bits;
            temp_upper = z_ll_mod_precomp(temp_upper2, temp_upper, res->p, res->p_inv);
            _zmod_poly_set_coeff_ui(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
            mpn[current_limb] >>= (bits - current_bit - 2*FLINT_BITS);
            current_bit = 3*FLINT_BITS + current_bit - bits;
        }

         if(current_bit == FLINT_BITS)
         {
            current_bit = 0; 
            current_limb++;  
         }     
      }
   }
}

void zmod_poly_mul_trunc_n(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
{
   if (poly1->length + poly2->length <= 6)
   {
      zmod_poly_mul_classical_trunc(res, poly1, poly2, trunc);
      return;
   }
   
   if ((FLINT_BIT_COUNT(poly1->p) <= 30) && (poly1->length + poly2->length <= 16))
   {
      zmod_poly_mul_classical_trunc(res, poly1, poly2, trunc);
      return;
   }
   
   zmod_poly_mul_KS_trunc(res, poly1, poly2, 0, trunc);
}

void zmod_poly_mul_trunc_left_n(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, unsigned long trunc)
{
   if (poly1->length + poly2->length <= 10)
   {
      zmod_poly_mul_classical_trunc_left(res, poly1, poly2, trunc);
      return;
   }
   
   if ((FLINT_BIT_COUNT(poly1->p) <= 30) && (poly1->length + poly2->length < 30))
   {
      zmod_poly_mul_classical_trunc_left(res, poly1, poly2, trunc);
      return;
   }
   
   zmod_poly_mul_KS(res, poly1, poly2, 0);
}

/*******************************************************************************

   Scalar multiplication

********************************************************************************/

/* 
   Scalar multiplication
   
   Assumes the scalar is reduced modulo poly->p
*/

void __zmod_poly_scalar_mul_without_mod(zmod_poly_t res, zmod_poly_t poly, unsigned long scalar)
{
   if (poly != res)
      zmod_poly_fit_length(res, poly->length);

   if (scalar == 0) 
   {
      res->length = 0;
      return;
   }
   
   if (scalar == 1L) 
   {
      _zmod_poly_set(res, poly);
      return;
   }
   
   for (unsigned long i = 0; i < poly->length; i++)
   {
       res->coeffs[i] = poly->coeffs[i] * scalar;
   }
   
   res->length = poly->length;
   __zmod_poly_normalise(res);
}

void _zmod_poly_scalar_mul(zmod_poly_t res, zmod_poly_t poly, unsigned long scalar)
{
   if (scalar == 0) 
   {
      res->length = 0;
      return;
   }
   
   if (scalar == 1L) 
   {
      _zmod_poly_set(res, poly);
      return;
   }
   
   unsigned long bits = FLINT_BIT_COUNT(poly->p);
   
#if FLINT_BITS == 64
   if (bits <= FLINT_D_BITS)
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
          res->coeffs[i] = z_mulmod_precomp(poly->coeffs[i], scalar, poly->p, poly->p_inv);
      }
   } else
   {
#endif
      for (unsigned long i = 0; i < poly->length; i++)
      {
          res->coeffs[i] = z_mulmod2_precomp(poly->coeffs[i], scalar, poly->p, poly->p_inv);
      }
#if FLINT_BITS == 64
   }
#endif
   
   res->length = poly->length;
   __zmod_poly_normalise(res);
}

void zmod_poly_scalar_mul(zmod_poly_t res, zmod_poly_t poly, unsigned long scalar)
{
   if (poly != res)
      zmod_poly_fit_length(res, poly->length);
      
   _zmod_poly_scalar_mul(res, poly, scalar);
}

/* 
   Used to reduce a polynomial modulo its modulus if it has been left without reduction
   for a while. Assumes all the coefficients are positive and at most FLINT_D_BITS.
*/

void __zmod_poly_scalar_mod(zmod_poly_t poly)
{
   unsigned long p = poly->p;
   double p_inv = poly->p_inv;

   for (unsigned long i = 0; i < poly->length; i++)
   {
      poly->coeffs[i] = z_mod_precomp(poly->coeffs[i], p, p_inv);
   }

   __zmod_poly_normalise(poly);
}

/*
   Classical basecase division
   
   Requires that the leading coefficient be invertible modulo B->p
*/

void zmod_poly_divrem_classical(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
{
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      zmod_poly_set(R, A);

      return;
   }
   
   unsigned long p = B->p;
   if (2*FLINT_BIT_COUNT(p) + FLINT_BIT_COUNT(A->length - B->length + 1) <= FLINT_D_BITS)
   {
      __zmod_poly_divrem_classical_mod_last(Q, R, A, B);
      return;
   }

   double p_inv = B->p_inv;
   unsigned long lead_inv = z_invert(B->coeffs[B->length - 1], p);
   unsigned long * coeff_Q;
   
   zmod_poly_t qB;
   zmod_poly_init2(qB, p, B->length);
   
   zmod_poly_t Bm1;
   _zmod_poly_attach_truncate(Bm1, B, B->length - 1);
   
   long coeff = A->length - 1;
   
   zmod_poly_set(R, A);
   
   if (A->length >= B->length)
   {
      zmod_poly_fit_length(Q, A->length - B->length + 1);
      Q->length = A->length - B->length + 1;
   } else zmod_poly_zero(Q); 

   coeff_Q = Q->coeffs - B->length + 1;
   
#if FLINT_BITS == 64
   int small = (FLINT_BIT_COUNT(p) <= FLINT_D_BITS);
#endif   
      
   while (coeff >= (long) B->length - 1)
   {
      while ((coeff >= (long) B->length - 1) && (R->coeffs[coeff] == 0L))
      {
         coeff_Q[coeff] = 0L;
         coeff--;
      }
      
      if (coeff >= (long) B->length - 1)
      {
#if FLINT_BITS == 64
         if (small) coeff_Q[coeff] = z_mulmod_precomp(R->coeffs[coeff], lead_inv, p, p_inv); 
         else    
#endif
         coeff_Q[coeff] = z_mulmod2_precomp(R->coeffs[coeff], lead_inv, p, p_inv);
         
         zmod_poly_scalar_mul(qB, Bm1, coeff_Q[coeff]);
         
         zmod_poly_t R_sub;
         R_sub->p = p;
         R_sub->coeffs = R->coeffs + coeff - B->length + 1;
         R_sub->length = B->length - 1;
         _zmod_poly_sub(R_sub, R_sub, qB);
         
         coeff--;
      }
   }
   
   R->length = B->length - 1;
   __zmod_poly_normalise(R);
   zmod_poly_clear(qB);
}

void __zmod_poly_divrem_classical_mod_last(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
{
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      zmod_poly_set(R, A);

      return;
   }
   unsigned long p = B->p;
   double p_inv = B->p_inv;
   unsigned long lead_inv = z_invert(B->coeffs[B->length - 1], p);
   unsigned long * coeff_Q;
   
   zmod_poly_t qB;
   zmod_poly_init2(qB, p, B->length);
   
   zmod_poly_t Bm1;
   _zmod_poly_attach_truncate(Bm1, B, B->length - 1);
   
   long coeff = A->length - 1;
   
   zmod_poly_set(R, A);
   
   if (A->length >= B->length)
   {
      zmod_poly_fit_length(Q, A->length - B->length + 1);
      Q->length = A->length - B->length + 1;
   } else zmod_poly_zero(Q); 

   coeff_Q = Q->coeffs - B->length + 1;
   
   while (coeff >= (long) B->length - 1)
   {
      R->coeffs[coeff] = z_mod_precomp(R->coeffs[coeff], p, p_inv);
  
      while ((coeff >= (long) B->length - 1) && (R->coeffs[coeff] == 0L))
      {
         coeff_Q[coeff] = 0L;
         coeff--;
         if (coeff >= (long) B->length - 1) 
         {
               R->coeffs[coeff] = z_mod_precomp(R->coeffs[coeff], p, p_inv);
         }
      }
      
      if (coeff >= (long) B->length - 1)
      {
         coeff_Q[coeff] = z_mulmod_precomp(R->coeffs[coeff], lead_inv, p, p_inv); 
         
         __zmod_poly_scalar_mul_without_mod(qB, Bm1, z_negmod(coeff_Q[coeff], p));
         
         zmod_poly_t R_sub;
         R_sub->p = p;
         R_sub->coeffs = R->coeffs + coeff - B->length + 1;
         R_sub->length = B->length - 1;
         _zmod_poly_add_without_mod(R_sub, R_sub, qB);
         
         coeff--;
      }
   }
   
   R->length = B->length - 1;
   __zmod_poly_scalar_mod(R);
   __zmod_poly_normalise(R);
   zmod_poly_clear(qB);
}

/*
   Classical basecase division, without remainder
   
   Requires that the leading coefficient be invertible modulo B->p
*/

void zmod_poly_div_classical(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
{
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      
      return;
   }
   
   unsigned long p = B->p;
   if (2*FLINT_BIT_COUNT(p) + FLINT_BIT_COUNT(A->length - B->length + 1) <= FLINT_D_BITS)
   {
      __zmod_poly_div_classical_mod_last(Q, A, B);
      return;
   }

   double p_inv = B->p_inv;
   unsigned long lead_inv = z_invert(B->coeffs[B->length - 1], p);
   unsigned long * coeff_Q;
   
   zmod_poly_t qB;
   zmod_poly_init2(qB, p, B->length);
   zmod_poly_t R;
   zmod_poly_init(R, p);
   
   zmod_poly_t Bm1;
   _zmod_poly_attach_truncate(Bm1, B, B->length - 1);
   
   long coeff = A->length - 1;
   
   zmod_poly_set(R, A);
   
   if (A->length >= B->length)
   {
      zmod_poly_fit_length(Q, A->length - B->length + 1);
      Q->length = A->length - B->length + 1;
   } else zmod_poly_zero(Q); 

   coeff_Q = Q->coeffs - B->length + 1;
   
#if FLINT_BITS == 64
   int small = (FLINT_BIT_COUNT(p) <= FLINT_D_BITS);
#endif   
      
   while (coeff >= (long) B->length - 1)
   {
      while ((coeff >= (long) B->length - 1) && (R->coeffs[coeff] == 0L))
      {
         coeff_Q[coeff] = 0L;
         coeff--;
      }
      
      if (coeff >= (long) B->length - 1)
      {
#if FLINT_BITS == 64
         if (small) coeff_Q[coeff] = z_mulmod_precomp(R->coeffs[coeff], lead_inv, p, p_inv); 
         else    
#endif
         coeff_Q[coeff] = z_mulmod2_precomp(R->coeffs[coeff], lead_inv, p, p_inv);
         
         if (coeff >= (long) B->length)
         {
            zmod_poly_scalar_mul(qB, Bm1, coeff_Q[coeff]);
         
            zmod_poly_t R_sub;
            R_sub->p = p;
            R_sub->coeffs = R->coeffs + coeff - B->length + 1;
            R_sub->length = B->length - 1;
            _zmod_poly_sub(R_sub, R_sub, qB);
         }
         
         coeff--;
      }
   }
   
   zmod_poly_clear(R);
   zmod_poly_clear(qB);
}

void __zmod_poly_div_classical_mod_last(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
{
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();      
   }
   
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      
      return;
   }
   
   unsigned long p = B->p;
   double p_inv = B->p_inv;
   unsigned long lead_inv = z_invert(B->coeffs[B->length - 1], p);
   unsigned long * coeff_Q;
   
   zmod_poly_t qB;
   zmod_poly_init2(qB, p, B->length);
   zmod_poly_t R;
   zmod_poly_init(R, p);
   
   zmod_poly_t Bm1;
   _zmod_poly_attach_truncate(Bm1, B, B->length - 1);
   
   long coeff = A->length - 1;
   
   zmod_poly_set(R, A);
   
   if (A->length >= B->length)
   {
      zmod_poly_fit_length(Q, A->length - B->length + 1);
      Q->length = A->length - B->length + 1;
   } else zmod_poly_zero(Q); 

   coeff_Q = Q->coeffs - B->length + 1;
   
   while (coeff >= (long) B->length - 1)
   {
      R->coeffs[coeff] = z_mod_precomp(R->coeffs[coeff], p, p_inv);

      while ((coeff >= (long) B->length - 1) && (R->coeffs[coeff] == 0L))
      {
         coeff_Q[coeff] = 0L;
         coeff--;
         if (coeff >= (long) B->length - 1) R->coeffs[coeff] = z_mod_precomp(R->coeffs[coeff], p, p_inv);
      }
      
      if (coeff >= (long) B->length - 1)
      {
         coeff_Q[coeff] = z_mulmod_precomp(R->coeffs[coeff], lead_inv, p, p_inv); 
         
         if (coeff >= (long) B->length)
         {
            zmod_poly_scalar_mul(qB, Bm1, z_negmod(coeff_Q[coeff], p));
         
            zmod_poly_t R_sub;
            R_sub->p = p;
            R_sub->coeffs = R->coeffs + coeff - B->length + 1;
            R_sub->length = B->length - 1;
            _zmod_poly_add_without_mod(R_sub, R_sub, qB);
         }
         
         coeff--;
      }
   }
   
   zmod_poly_clear(R);
   zmod_poly_clear(qB);
}

/*
   Divide and conquer division
*/

void zmod_poly_div_divconquer_recursive(zmod_poly_t Q, zmod_poly_t BQ, zmod_poly_t A, zmod_poly_t B)
{
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      zmod_poly_zero(BQ);

      return;
   }
   
   // A->length is now >= B->length
   
   unsigned long p = A->p;
   unsigned long crossover = 16;
   unsigned long crossover2 = 128;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      /*
         Use the classical algorithm to compute the
         quotient and remainder, then use A-R to compute BQ
      */
      
      zmod_poly_t Rb;
      zmod_poly_init(Rb, p);
      zmod_poly_divrem_classical(Q, Rb, A, B);
      zmod_poly_sub(BQ, A, Rb);
      zmod_poly_clear(Rb);
      
      return;
   }
   
   zmod_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _zmod_poly_attach_shift(d1, B, n2);
   _zmod_poly_attach_truncate(d2, B, n2);
   _zmod_poly_attach_shift(d3, B, n1);
   _zmod_poly_attach_truncate(d4, B, n1);
   
   if (A->length <= n2 + B->length - 1)
   {
      /*
         If A->length <= B->length + n2 - 1
         then only a single quotient is needed
         We do a division of at most 2*n2 - 1
         terms by n2 terms yielding a quotient of
         at most n2 terms 
      */
      
      // Set p1 to be A without the last
      // n1 coefficients
      // 2*n2-1 >= p1->length > 0
      
      zmod_poly_init(p1, p);
      zmod_poly_right_shift(p1, A, n1);
      
      // Since A was normalised, then p1 will be
      // d3 is the leading terms of B and so must be normalised
      // d3 is length n2, so we get at most n2 terms in the quotient
      
      zmod_poly_init(d1q1, p);
      zmod_poly_div_divconquer_recursive(Q, d1q1, p1, d3); 
      zmod_poly_clear(p1);
      
      /*
         Compute d2q1 = Q*d4
         It is of length at most n1+n2-1 terms
      */
      
      zmod_poly_init(d2q1, p);
      zmod_poly_mul(d2q1, Q, d4);
      
      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length at most n1+2*n2-1
      */
      
      zmod_poly_left_shift(BQ, d1q1, n1);
      zmod_poly_clear(d1q1);
      zmod_poly_add(BQ, BQ, d2q1);
      zmod_poly_clear(d2q1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      unsigned long shift = A->length - 2*B->length + 1;
      _zmod_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = d1*q1 is length at most 2*B->length-1
      */
      
      zmod_poly_init(d1q1, p);
      zmod_poly_init(q1, p);
      
      zmod_poly_div_divconquer_recursive(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      zmod_poly_init(dq1, p);
      
      zmod_poly_left_shift(dq1, d1q1, shift);
      zmod_poly_clear(d1q1);
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      zmod_poly_init(t, p);
      zmod_poly_sub(t, A, dq1);
      zmod_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      zmod_poly_init(q2, p);
      zmod_poly_init(dq2, p);
      zmod_poly_div_divconquer_recursive(q2, dq2, t, B); 
      zmod_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      zmod_poly_left_shift(Q, q1, shift);
      zmod_poly_clear(q1);
      zmod_poly_add(Q, Q, q2);
      zmod_poly_clear(q2);
      
      /*
         Write out BQ = dq1 + dq2
      */
      
      zmod_poly_add(BQ, dq1, dq2);
      zmod_poly_clear(dq1);
      zmod_poly_clear(dq2);
      
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
      
   zmod_poly_init(p1, p);
   zmod_poly_right_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length n1
      d1q1 = d1*q1 is length at most 2*n1-1
   */
      
   zmod_poly_init(d1q1, p);
   zmod_poly_init(q1, p);
   zmod_poly_div_divconquer_recursive(q1, d1q1, p1, d1); 
   zmod_poly_clear(p1);   
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being at most length n1+n2-1
   */  
   
   zmod_poly_init(d2q1, p);
   zmod_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most 2*n1+n2-1
   */
   
   zmod_poly_init(dq1, p);
   zmod_poly_left_shift(dq1, d1q1, n2);
   zmod_poly_clear(d1q1);
   zmod_poly_add(dq1, dq1, d2q1);
   zmod_poly_clear(d2q1);
   
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1
      which has length at most 2*n1+n2-1, but we are not interested 
      in up to the first n1 coefficients, so it has 
      effective length at most n1+n2-1
   */
   
   zmod_poly_init(t, p);
   zmod_poly_right_shift(t, A, n2);
   zmod_poly_sub(t, t, dq1);
   zmod_poly_truncate(t, B->length - 1);
   
   /*
      Compute q2 = t div d1
      It is at most an n1+n2-1 by n1 division, so
      the length of q2 will be at most n2
      Also compute d1q2 of length at most n1+n2-1
   */
   
   zmod_poly_init(d1q2, p);
   zmod_poly_init(q2, p);
   zmod_poly_div_divconquer_recursive(q2, d1q2, t, d1); 
   zmod_poly_clear(t);
      
   /*
      Compute d2q2 = d2*q2 which is of length 
      at most n1+n2-1
   */
   
   zmod_poly_init(d2q2, p);
   zmod_poly_mul(d2q2, d2, q2);
   
   /*
      Compute dq2 = d1*q2*x^n2 + d2q2
      which is of length at most n1+2*n2-1
   */
   
   zmod_poly_init(dq2, p);
   zmod_poly_left_shift(dq2, d1q2, n2);
   zmod_poly_clear(d1q2);
   zmod_poly_add(dq2, dq2, d2q2);
   zmod_poly_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length at most n1+n2
   */
   
   zmod_poly_left_shift(Q, q1, n2);
   zmod_poly_clear(q1);
   zmod_poly_add(Q, Q, q2);
   zmod_poly_clear(q2);
   
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length at most 2*(n1+n2)-1
   */
   
   zmod_poly_left_shift(BQ, dq1, n2);
   zmod_poly_add(BQ, BQ, dq2);
   
   zmod_poly_clear(dq2);
   zmod_poly_clear(dq1);
}

void zmod_poly_div_divconquer(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
{
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      
      return;
   }

   // A->length is now >= B->length
    
   unsigned long crossover = 16;
   unsigned long crossover2 = 256;
   
   unsigned long p = B->p;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      zmod_poly_div_classical(Q, A, B);
      
      return;
   }
   
   // B->length is now >= crossover (16)
   
   zmod_poly_t d1, d2, d3, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
      
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   // n1 and n2 are at least 4
   
   /* We let B = d1*x^n2 + d2 
      d1 is of length n1 and
      d2 of length n2
   */

   _zmod_poly_attach_shift(d1, B, n2);
   _zmod_poly_attach_truncate(d2, B, n2);
   _zmod_poly_attach_shift(d3, B, n1);
   
   if (A->length <= n2 + B->length - 1)
   {
      /*
         If A->length <= B->length + n2 - 1
         then only a single quotient is needed
         We do a division of at most 2*n2 - 1
         terms by n2 terms yielding a quotient of
         at most n2 terms 
      */
      
      // Set p1 to be A without the last
      // n1 coefficients
      // 2*n2-1 >= p1->length > 0
      
      zmod_poly_init(p1, p);
      zmod_poly_right_shift(p1, A, n1);
      
      // Since A was normalised, then p1 will be
      // d3 is the leading terms of B and so must be normalised
      // d3 is length n2, so we get at most n2 terms in the quotient
      
      zmod_poly_div_divconquer(Q, p1, d3); 
      zmod_poly_clear(p1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      unsigned long shift = A->length - 2*B->length + 1;
      _zmod_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = low(d1*q1) is length at most 2*B->length-1
         We discard the lower B->length-1 terms
      */
      
      zmod_poly_init(d1q1, p);
      zmod_poly_init(q1, p);
      
      zmod_poly_div_divconquer_recursive(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      zmod_poly_init(dq1, p);
      
      zmod_poly_left_shift(dq1, d1q1, shift);
      zmod_poly_clear(d1q1); 
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      zmod_poly_init(t, p);
      zmod_poly_sub(t, A, dq1);
      zmod_poly_clear(dq1);
      zmod_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      zmod_poly_init(q2, p);
      zmod_poly_div_divconquer(q2, t, B); 
      zmod_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      zmod_poly_left_shift(Q, q1, shift);
      zmod_poly_clear(q1);
      zmod_poly_add(Q, Q, q2);
      zmod_poly_clear(q2);
      
      return;
   }
   // We now have n2 + B->length - 1 < A->length <= 2*B->length - 1
   
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length at most n1 and a2 is length n2 
      and a3 is length n1+n2-1 
   */
      
      // Set p1 to a1*x^(n1-1) + other terms
      // It has length at most 2*n1-1 and is normalised
      // A->length >= 2*n2
      
      zmod_poly_init(p1, p);
      zmod_poly_right_shift(p1, A, 2*n2);
      
      /* 
         Set q1 to p1 div d1 
         This is at most a 2*n1-1 by n1 division so 
         q1 ends up being at most length n1
         d1q1 = low(d1*q1) is length at most n1-1
         Thus we have discarded the leading n1 terms (at most)
      */
      
      zmod_poly_init(d1q1, p);
      zmod_poly_init(q1, p);
      
      zmod_poly_div_divconquer_recursive(q1, d1q1, p1, d1); 
      zmod_poly_clear(p1);
   
   /* 
      Compute d2q1 = d2*q1 with low n1 - 1 terms zeroed
      d2*q1 is length at most n1+n2-1 leaving at most
      n2 non-zero terms to the left
   */  
   
   zmod_poly_init(d2q1, p); 
   zmod_poly_mul_trunc_left_n(d2q1, d2, q1, n1 - 1);
       
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most 2*n1+n2-1
      but may have any length below this        
   */
   
   zmod_poly_init(dq1, p);
   zmod_poly_left_shift(dq1, d1q1, n2);
   zmod_poly_clear(d1q1); 
   zmod_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = a1*x^(2*n2-1) + a2*x^(n2-1) - dq1 
      after shifting dq1 to the right by (n1-n2)
      which has length at most 2*n1+n2-1, but we 
      discard up to n1 coefficients, so it has 
      effective length 2*n2-1 with the last n2-1
      coefficients ignored. Thus there are at most n2 
      significant coefficients
   */
   
   
   zmod_poly_init(t, p);
   zmod_poly_right_shift(t, A, n1);
   _zmod_poly_attach_shift(temp, dq1, n1-n2);
   zmod_poly_sub(t, t, temp);
   zmod_poly_truncate(t, 2*n2-1);
     
   /*
      Compute q2 = t div d3
      It is at most a 2*n2-1 by n2 division, so
      the length of q2 will be n2 at most
   */
   
   zmod_poly_init(q2, p);
   zmod_poly_div_divconquer(q2, t, d3); 
   zmod_poly_clear(t);  
   zmod_poly_clear(dq1);
   zmod_poly_clear(d2q1);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   
   zmod_poly_left_shift(Q, q1, n2);
   zmod_poly_clear(q1);
   zmod_poly_add(Q, Q, q2);
   zmod_poly_clear(q2);   
}

void zmod_poly_divrem_divconquer(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
{
   zmod_poly_t QB;
   zmod_poly_init(QB, B->p);
   
   zmod_poly_div_divconquer_recursive(Q, QB, A, B);   
   zmod_poly_sub(R, A, QB);
   
   zmod_poly_clear(QB);
}

/****************************************************************************

   Newton Inversion

****************************************************************************/

#define FLINT_ZMOD_NEWTON_INVERSE_BASECASE_CUTOFF 64 //32

/*
   Compute the polynomial X^{2n} / Q. 
   Used by Newton iteration to bootstrap power series inversion.
   Q must have length >= n and leading coefficient invertible with respect to the modulus.
*/

void zmod_poly_newton_invert_basecase(zmod_poly_t Q_inv, zmod_poly_t Q, unsigned long n)
{
   zmod_poly_t X2n, Qn;
   
   zmod_poly_init2(X2n, Q->p, 2*n-1);
   zmod_poly_set_coeff_ui(X2n, 2*n - 2, 1L);
   
   _zmod_poly_attach_shift(Qn, Q, Q->length - n);
   
   zmod_poly_div_divconquer(Q_inv, X2n, Qn);
      
   zmod_poly_clear(X2n);
}

/*
   Recursively compute 1 / Q mod x^n using Newton iteration
   Assumes Q is given as a power series to the full precision n required 
   with invertible constant term with respect to the modulus
*/

#define FLINT_ZMOD_NEWTON_INVERSE_CACHE_CUTOFF 2000

void zmod_poly_newton_invert(zmod_poly_t Q_inv, zmod_poly_t Q, unsigned long n)
{
   if (n < FLINT_ZMOD_NEWTON_INVERSE_BASECASE_CUTOFF)
   {
      zmod_poly_t Q_rev;
      zmod_poly_init2(Q_rev, Q->p, n);
		_zmod_poly_reverse(Q_rev, Q, n);
      zmod_poly_newton_invert_basecase(Q_inv, Q_rev, n);
      zmod_poly_reverse(Q_inv, Q_inv, n);
      zmod_poly_clear(Q_rev);
      
      return;
   }
   
   unsigned long m = (n+1)/2;
   unsigned long p = Q->p;
   
   zmod_poly_t g0, prod, prod2;
   zmod_poly_init(g0, p);
   zmod_poly_init(prod, p);
   zmod_poly_init(prod2, p);
   zmod_poly_newton_invert(g0, Q, m);
#if USE_MIDDLE_PRODUCT
   if (n < FLINT_ZMOD_NEWTON_INVERSE_CACHE_CUTOFF)
   {
#endif
      zmod_poly_mul_trunc_n(prod, Q, g0, n);
      prod->coeffs[0] = z_submod(prod->coeffs[0], 1L, p);
      //zmod_poly_mul_trunc_n(prod2, prod, g0, n);
      zmod_poly_t prod_s;
      zmod_poly_attach_shift(prod_s, prod, (n+1)/2);
      
      zmod_poly_t prod2_s;
      zmod_poly_fit_length(prod2, n);
      zmod_poly_attach_shift(prod2_s, prod2, (n+1)/2);
      _zmod_poly_mul_KS_trunc(prod2_s, prod_s, g0, 0, n - (n+1)/2);
      
      prod2->length = (n+1)/2 + prod2_s->length;
      for (unsigned long i = 0; i < (n+1)/2; i++)
         prod2->coeffs[i] = 0L;
#if USE_MIDDLE_PRODUCT
   } else
   {
      zmod_poly_precomp_t pre;
      zmod_poly_mul_trunc_n_precomp_init(pre, g0, 0, (n+1)/2);
      zmod_poly_fit_length(prod, n);
      _zmod_poly_mul_KS_middle_precomp(prod, Q, pre, 0, n);
      zmod_poly_t prod_s;
      zmod_poly_attach_shift(prod_s, prod, (n+1)/2);
      
      zmod_poly_t prod2_s;
      zmod_poly_fit_length(prod2, n);
      zmod_poly_attach_shift(prod2_s, prod2, (n+1)/2);
      _zmod_poly_mul_KS_trunc_precomp(prod2_s, prod_s, pre, 0, n - (n+1)/2);
      zmod_poly_precomp_clear(pre);
      prod2->length = (n+1)/2 + prod2_s->length;
      for (unsigned long i = 0; i < (n+1)/2; i++)
         prod2->coeffs[i] = 0L;
   }
#endif
   zmod_poly_sub(Q_inv, g0, prod2);
   
   zmod_poly_clear(prod2);
   zmod_poly_clear(prod);
   zmod_poly_clear(g0);
}

/****************************************************************************

   Newton Division

****************************************************************************/

/* 
   Yields a precision n power series quotient of A by B assuming A and B are both 
   given to precision n and B is normalised (i.e. constant coefficient is invertible).
*/

void zmod_poly_div_series(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B, unsigned long n)
{
   zmod_poly_t Ain, Bin;
   unsigned long p = B->p;
   
   if (A == Q)
   {
      zmod_poly_init(Ain, p);
      zmod_poly_set(Ain, A);
   } else _zmod_poly_attach(Ain, A);
   
   if (B == Q)
   {
      zmod_poly_init(Bin, p);
      zmod_poly_set(Bin, B);
   } else _zmod_poly_attach(Bin, B);

   zmod_poly_t B_inv;
   zmod_poly_init(B_inv, p);
   zmod_poly_newton_invert(B_inv, Bin, n);
   zmod_poly_mul_trunc_n(Q, B_inv, Ain, n);
   
   zmod_poly_clear(B_inv);

   if (A == Q) zmod_poly_clear(Ain);
   if (B == Q) zmod_poly_clear(Bin);
}

/*
   Polynomial division of A by B
   The remainder is not computed, to save time
*/

void zmod_poly_div_newton(zmod_poly_t Q, zmod_poly_t A, zmod_poly_t B)
{
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      return;
   }
   
   unsigned long p = B->p;
   
   zmod_poly_t A_rev, B_rev;
   zmod_poly_init2(A_rev, p, A->length);
   zmod_poly_init2(B_rev, p, B->length);
   
   zmod_poly_reverse(A_rev, A, A->length);
   zmod_poly_reverse(B_rev, B, B->length);
   
   zmod_poly_div_series(Q, A_rev, B_rev, A->length - B->length + 1);
   
   zmod_poly_reverse(Q, Q, A->length - B->length + 1);
   
   zmod_poly_clear(B_rev);
   zmod_poly_clear(A_rev);
}

void zmod_poly_divrem_newton(zmod_poly_t Q, zmod_poly_t R, zmod_poly_t A, zmod_poly_t B)
{
   if (A->length < B->length)
   {
      zmod_poly_zero(Q);
      zmod_poly_set(R, A);
      return;
   }
   
   zmod_poly_t QB, A_trunc;
   zmod_poly_init(QB, B->p);
   
   zmod_poly_div_newton(Q, A, B);
   zmod_poly_mul_trunc_n(QB, Q, B, B->length - 1);
   _zmod_poly_attach_truncate(A_trunc, A, B->length - 1);
   zmod_poly_sub(R, A_trunc, QB);
   
   zmod_poly_clear(QB); 
}

/****************************************************************************

   GCD/resultant

****************************************************************************/

void zmod_poly_gcd(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   zmod_poly_t Q, R, A, B;
   
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      zmod_poly_zero(res);
      return;
   }

   if ((poly1->length == 1) || (poly2->length == 1))
   {
      zmod_poly_set_coeff_ui(res, 0, 1L);
      res->length = 1;
      return;
   }
   
   unsigned long p = poly1->p;
   zmod_poly_init(Q, p);
   zmod_poly_init(R, p);    

   if (poly1->length > poly2->length)
   {
      _zmod_poly_attach(A, poly1);
      _zmod_poly_attach(B, poly2);
   } else
   {
      _zmod_poly_attach(A, poly2);
      _zmod_poly_attach(B, poly1);
   }

   int steps = 1;
      
   while (B->length > 1)
   {
      zmod_poly_divrem(Q, R, A, B);
      zmod_poly_swap(A, B);
      if (steps > 2) zmod_poly_clear(B);
      _zmod_poly_attach(B, R);
      zmod_poly_init(R, p); 
      steps++;   
   }
      
   if  (B->length == 1) 
   {
      zmod_poly_set_coeff_ui(res, 0, 1L);
      res->length = 1;
   }
   else zmod_poly_set(res, A);

   if (steps > 2) 
   {
      zmod_poly_clear(A);
   } 

   zmod_poly_clear(B);
   zmod_poly_clear(R);
   zmod_poly_clear(Q);
}

/* 
   Computes poly1^(-1) mod poly2
   Assumes poly1 is not zero and is already reduced mod poly2
*/

int zmod_poly_gcd_invert(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   zmod_poly_t Q, R, A, B, u1, u2, prod;
   long a;
   unsigned long coprime;

   if (poly1->length == 0)
   {
      printf("FLINT Exception: Divide by zero\n");
      abort();
   }
   
   if (poly1->length == 1)
   {
      z_gcd_invert(&a, poly1->coeffs[0], poly2->p);
      if (a < 0L) zmod_poly_set_coeff_ui(res, 0, poly2->p + a);
	  else zmod_poly_set_coeff_ui(res, 0, a);
      res->length = 1;
      return 1;
   }
   
   unsigned long p = poly1->p;
   zmod_poly_init(Q, p);
   zmod_poly_init(R, p);
   zmod_poly_init(u1, p);
   zmod_poly_init(u2, p);
   zmod_poly_init(prod, p);

   zmod_poly_set_coeff_ui(u2, 0, 1L);
   u2->length = 1;
   zmod_poly_zero(u1);    

   _zmod_poly_attach(A, poly2);
   _zmod_poly_attach(B, poly1);
   
   int steps = 1;
   
   while (B->length > 1)
   {
      zmod_poly_divrem(Q, R, A, B);
     
      zmod_poly_mul(prod, Q, u2);
      zmod_poly_swap(u1, u2);
      zmod_poly_sub(u2, u2, prod);

      zmod_poly_swap(A, B);
      if (steps > 2) zmod_poly_clear(B);
      _zmod_poly_attach(B, R);
      zmod_poly_init(R, p); 
      steps++;    
   }
   
   if (B->length == 1) 
   {
      zmod_poly_swap(u1, u2);     
      zmod_poly_set(res, u1);
      zmod_poly_scalar_mul(res, res, z_invert(B->coeffs[0], p));
      coprime = 1;
   } else
   {
      coprime = 0;
   }
   

   if (steps > 2) 
   {
      zmod_poly_clear(A);
   } 

   zmod_poly_clear(u1);
   zmod_poly_clear(u2);
   zmod_poly_clear(prod);
   zmod_poly_clear(B);
   zmod_poly_clear(R);
   zmod_poly_clear(Q);
   
   return coprime;
}

/*
   Compute res = gcd(poly1, poly2)
   Find s and t such that res = s*poly1 + t*poly2
*/

void zmod_poly_xgcd(zmod_poly_t res, zmod_poly_t s, zmod_poly_t t, zmod_poly_t poly1, zmod_poly_t poly2)
{
   zmod_poly_t Q, R, A, B, u1, u2, v1, v2, prod;
   unsigned long a;

   if ((poly1->length == 0) || (poly2->length == 0))
   {
      zmod_poly_zero(s);
      zmod_poly_zero(t);
      zmod_poly_zero(res);
      return;
   }
   
   if (poly1->length == 1)
   {
      a = z_invert(poly1->coeffs[0], poly2->p);
      zmod_poly_set_coeff_ui(s, 0, a);
      s->length = 1;
      zmod_poly_set_coeff_ui(res, 0, 1L);
      res->length = 1;
      zmod_poly_zero(t);
      return;
   }
   
   if (poly2->length == 1)
   {
      a = z_invert(poly2->coeffs[0], poly2->p);
      zmod_poly_set_coeff_ui(t, 0, a);
      t->length = 1;
      zmod_poly_set_coeff_ui(res, 0, 1L);
      res->length = 1;
      zmod_poly_zero(s);
      return;
   }
   
   unsigned long p = poly1->p;
   zmod_poly_init(Q, p);
   zmod_poly_init(R, p);
   zmod_poly_init(u1, p);
   zmod_poly_init(u2, p);
   zmod_poly_init(v1, p);
   zmod_poly_init(v2, p);
   zmod_poly_init(prod, p);

   zmod_poly_set_coeff_ui(u1, 0, 1L);
   u1->length = 1;
   zmod_poly_zero(u2);    
   zmod_poly_set_coeff_ui(v2, 0, 1L);
   v2->length = 1;
   zmod_poly_zero(v1);    

   if (poly1->length > poly2->length)
   {
      _zmod_poly_attach(A, poly1);
      _zmod_poly_attach(B, poly2);
   } else
   {
      _zmod_poly_attach(A, poly2);
      _zmod_poly_attach(B, poly1);
      zmod_poly_swap(u1, u2);
      zmod_poly_swap(v1, v2);
   }

   int steps = 1;
   
   while (B->length > 1)
   {
      zmod_poly_divrem(Q, R, A, B);
     
      zmod_poly_mul(prod, Q, u2);
      zmod_poly_swap(u1, u2);
      zmod_poly_sub(u2, u2, prod);

      zmod_poly_mul(prod, Q, v2);
      zmod_poly_swap(v1, v2);
      zmod_poly_sub(v2, v2, prod);

      zmod_poly_swap(A, B);
      if (steps > 2) zmod_poly_clear(B);
      _zmod_poly_attach(B, R);
      zmod_poly_init(R, p); 
      steps++;    
   }
   
   if (B->length == 1) 
   {      
      zmod_poly_swap(u1, u2);      
      zmod_poly_swap(v1, v2);
      
      zmod_poly_set(res, B);
   } else zmod_poly_set(res, A);  

   zmod_poly_set(s, u1);
   zmod_poly_set(t, v1);
   zmod_poly_scalar_mul(s, s, z_invert(res->coeffs[res->length-1], p));
   zmod_poly_scalar_mul(t, t, z_invert(res->coeffs[res->length-1], p));

   zmod_poly_make_monic(res, res);
   
   if (steps > 2) 
   {
      zmod_poly_clear(A);
   } 

   zmod_poly_clear(u1);
   zmod_poly_clear(u2);
   zmod_poly_clear(v1);
   zmod_poly_clear(v2);
   zmod_poly_clear(prod);
   zmod_poly_clear(B);
   zmod_poly_clear(R);
   zmod_poly_clear(Q);
}

unsigned long zmod_poly_resultant_euclidean(zmod_poly_t a, zmod_poly_t b)
{
   unsigned long res;
   
   if ((a->length == 0) || (b->length == 0))
   {
      return 0;
   }
   
   if ((a->length == 1) || (b->length == 1)) 
   {
      return 1;
   }

   unsigned long p = a->p;
   double p_inv = a->p_inv;
 
   unsigned long long l0, l1, l2;
   
   unsigned long lc;
   
   res = 1L;

   zmod_poly_t u, v, q;
   zmod_poly_init(u, p);
   zmod_poly_init(v, p);
   zmod_poly_init(q, p);
   
   zmod_poly_set(u, a);
   zmod_poly_set(v, b);

   for (;;) {
      l0 = u->length;
      l1 = v->length;
      lc = v->coeffs[v->length - 1];

      zmod_poly_divrem(q, u, u, v);
      
      zmod_poly_swap(u, v);

      l2 = v->length;
      if (l2 >= 1) 
      {
         lc = z_powmod2_precomp(lc, l0 - l2, p, p_inv);
         res = z_mulmod2_precomp(res, lc, p, p_inv);
         if (((l0 | l1) & 1) == 0) 
         {
            if (res) res = p - res;
         }  
      } else 
      {
         if (l1 == 1) {
            lc = z_powmod2_precomp(lc, l0 - 1, p, p_inv);
            res = z_mulmod2_precomp(res, lc, p, p_inv);
         } else
            res = 0L;
        
         break;
      }
   }

   zmod_poly_clear(q);
   zmod_poly_clear(u);
   zmod_poly_clear(v);

   return res;
}

long zmod_poly_half_gcd_iter(zmod_poly_2x2_mat_t res, zmod_poly_t a, zmod_poly_t b)
{
	ulong m = a->length/2;
	
	zmod_poly_fit_length(res->a, 1);
   zmod_poly_fit_length(res->d, 1);
   zmod_poly_zero(res->a);
   zmod_poly_zero(res->b);
   zmod_poly_zero(res->c);
   zmod_poly_zero(res->d);
   zmod_poly_set_coeff_ui(res->a, 0, 1);
   zmod_poly_set_coeff_ui(res->d, 0, 1);
   
	if (b->length < m+1) return 1;
	
	zmod_poly_t U, V, Q, temp;
	zmod_poly_init(U, a->p);
	zmod_poly_init(V, a->p);
   zmod_poly_init(Q, a->p);
   zmod_poly_init(temp, a->p);
  
	zmod_poly_set(U, a);
   zmod_poly_set(V, b);

	long sign = 1;

	while (V->length >= m+1)
	{
		zmod_poly_divrem(Q, U, U, V);
		zmod_poly_swap(U, V);
      
		zmod_poly_mul(temp, Q, res->c);
		zmod_poly_add(temp, res->d, temp);
      zmod_poly_swap(res->d, res->c);
		zmod_poly_swap(res->c, temp);

		zmod_poly_mul(temp, Q, res->a);
		zmod_poly_add(temp, res->b, temp);
      zmod_poly_swap(res->b, res->a);
		zmod_poly_swap(res->a, temp);
		sign = -sign;
	}

	zmod_poly_clear(temp);
	zmod_poly_clear(Q);
	zmod_poly_clear(U);
	zmod_poly_clear(V);

	return sign;
}

#define FLINT_ZMOD_POLY_HGCD_CUTOFF 75

long zmod_poly_half_gcd(zmod_poly_2x2_mat_t res, zmod_poly_t a, zmod_poly_t b)
{
   ulong m = a->length/2;

	if (b->length < m + 1)
	{
		zmod_poly_fit_length(res->a, 1);
      zmod_poly_fit_length(res->d, 1);
		zmod_poly_zero(res->a);
      zmod_poly_zero(res->b);
      zmod_poly_zero(res->c);
      zmod_poly_zero(res->d);
      zmod_poly_set_coeff_ui(res->a, 0, 1);
      zmod_poly_set_coeff_ui(res->d, 0, 1);
      return 1L;
	}
	//printf("hgcd: b->length = %ld\n", b->length);
	
   zmod_poly_t a0, b0;
	_zmod_poly_attach_shift(a0, a, m);
   _zmod_poly_attach_shift(b0, b, m);

	long R_sign;
	
	if (a0->length < FLINT_ZMOD_POLY_HGCD_CUTOFF) R_sign = zmod_poly_half_gcd_iter(res, a0, b0);
	else R_sign = zmod_poly_half_gcd(res, a0, b0);
   
	zmod_poly_t temp, a2, b2; 
	zmod_poly_init(temp, a->p);
   zmod_poly_init(b2, a->p);
   
	zmod_poly_mul(b2, res->c, a);
	zmod_poly_mul(temp, res->a, b);
	if (R_sign < 0L) zmod_poly_sub(b2, b2, temp);
	else zmod_poly_sub(b2, temp, b2);

	//printf("hgcd: b2->length = %ld\n", b2->length);
	
	if (b2->length < m + 1)
	{
	   zmod_poly_clear(temp);
	   zmod_poly_clear(b2);
		return R_sign;
	}
   
	zmod_poly_init(a2, a->p);
   
	zmod_poly_mul(a2, res->d, a);
   zmod_poly_mul(temp, res->b, b);
   if (R_sign < 0L) zmod_poly_sub(a2, temp, a2);
	else zmod_poly_sub(a2, a2, temp);
   
	zmod_poly_t q, d;
	zmod_poly_init(q, a->p);
	zmod_poly_init(d, a->p);

	zmod_poly_divrem(q, d, a2, b2);
	long k = 2*m - b2->length + 1;
	//printf("hgcd: d->length = %ld\n", d->length);
	
	zmod_poly_t c0, d0;
   _zmod_poly_attach_shift(c0, b2, k);
   _zmod_poly_attach_shift(d0, d, k);
   zmod_poly_2x2_mat_t S;
	zmod_poly_2x2_mat_init(S, a->p);
	
	long S_sign;
	if (c0->length < FLINT_ZMOD_POLY_HGCD_CUTOFF) S_sign = zmod_poly_half_gcd_iter(S, c0, d0);
	else R_sign = zmod_poly_half_gcd(S, c0, d0);
   
	zmod_poly_swap(S->a, S->c);
	zmod_poly_swap(S->b, S->d);
   zmod_poly_mul(temp, S->c, q);
	zmod_poly_add(S->a, S->a, temp);
   zmod_poly_mul(temp, S->d, q);
	zmod_poly_add(S->b, S->b, temp);

	zmod_poly_2x2_mat_mul(res, res, S);
   
	zmod_poly_2x2_mat_clear(S);
	zmod_poly_clear(temp);
	zmod_poly_clear(q);
	zmod_poly_clear(d);
   zmod_poly_clear(a2);
	zmod_poly_clear(b2);
	
	return -R_sign*S_sign;
}

#define FLINT_ZMOD_POLY_GCD_CUTOFF 100

void zmod_poly_gcd_hgcd(zmod_poly_t res, zmod_poly_t f, zmod_poly_t g)
{
   if (f->length == 0)
	{
		zmod_poly_set(res, g);
		return;
	}
	
	ulong p = f->p;
	
	zmod_poly_t h, j, q, r, temp;
	zmod_poly_init(q, p);
	zmod_poly_init(r, p);

	zmod_poly_divrem(q, r, f, g);
	if (r->length == 0)
	{
		zmod_poly_set(res, g);
	   zmod_poly_clear(q);
	   zmod_poly_clear(r);
      return;
	}

	if (g->length < FLINT_ZMOD_POLY_GCD_CUTOFF)
	{
		zmod_poly_gcd(res, g, r);
		return;
	}

	zmod_poly_2x2_mat_t R;
	zmod_poly_2x2_mat_init(R, p);

	zmod_poly_half_gcd(R, g, r);
	
	zmod_poly_init(j, p);
	zmod_poly_init(h, p);
	zmod_poly_init(temp, p);

	zmod_poly_mul(h, R->d, g);
	zmod_poly_mul(temp, R->b, r);
	zmod_poly_sub(h, h, temp);

	zmod_poly_mul(j, R->c, g);
	zmod_poly_mul(temp, R->a, r);
	zmod_poly_sub(j, temp, j);

	while (j->length != 0)
	{
      zmod_poly_divrem(q, r, h, j);
	
		if (r->length == 0)
	   {
		   zmod_poly_set(res, j);
	      zmod_poly_2x2_mat_clear(R);
	      zmod_poly_clear(temp);
	      zmod_poly_clear(j);
	      zmod_poly_clear(q);
	      zmod_poly_clear(r);
         return;
	   }

		if (j->length < FLINT_ZMOD_POLY_GCD_CUTOFF)
	   {
		   zmod_poly_gcd(res, j, r);
		   return;
	   }

      zmod_poly_half_gcd(R, j, r);
	   
      zmod_poly_mul(h, R->d, j);
	   zmod_poly_mul(temp, R->b, r);
	   zmod_poly_sub(h, h, temp);

	   zmod_poly_mul(j, R->c, j);
	   zmod_poly_mul(temp, R->a, r);
	   zmod_poly_sub(j, temp, j);
		printf("j->length = %ld\n", j->length);
	}

	zmod_poly_set(res, h);

	zmod_poly_2x2_mat_clear(R);
	zmod_poly_clear(temp);
	zmod_poly_clear(j);
	zmod_poly_clear(q);
	zmod_poly_clear(r);
}

/****************************************************************************

   Differentiation

****************************************************************************/

/**
 * Set's the polynomial x_primed to be the derivative of the polynomial x. It reduces all the coefficients
 * mod p after the differentiation so x_primed is in <b>Z</b>/p<b>Z</b>[X]
 * @param  x_primed  The derivative of the x param
 * @param  x The polynomial to be differentiated
 * @return void
 */
void zmod_poly_derivative(zmod_poly_t x_primed, zmod_poly_t x)
{
	//get the degree and modulus of the polynomial we are working with
	unsigned long length = zmod_poly_length(x);
	unsigned long p = zmod_poly_modulus(x);
	
	if (length <= 1) 
	{
	   zmod_poly_zero(x_primed);
	   return;
    }

    zmod_poly_fit_length(x_primed, length - 1);	
	pre_inv_t p_inv = z_precompute_inverse(p);
	
	unsigned long k = 1L;
	for (unsigned long j = 1; j < length; j++)
	{
		if (k == 0L) x_primed->coeffs[j-1] = 0L; 
		else if (k == 1L) x_primed->coeffs[j-1] = zmod_poly_get_coeff_ui(x, j);
		else x_primed->coeffs[j-1] = z_mulmod2_precomp(zmod_poly_get_coeff_ui(x, j), k, p, p_inv);
		if (k == p - 1) k = 0L;
		else k++;
	}

	x_primed->length = length - 1;
	__zmod_poly_normalise(x_primed); // final coefficients may be zero
}

/****************************************************************************

   Modular arithmetic (modulo a polynomial)

****************************************************************************/

/*
   Multiplies poly1 and poly2 modulo f
   Assumes poly1 and poly2 are reduced mod f
*/

void zmod_poly_mulmod(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2, zmod_poly_t f)
{
   unsigned long p = f->p;
   
   if (f->length == 0)
   {
      printf("FLINT Exception: Divide by zero\n");
      abort();
   }
   if ((f->length == 1) || (poly1->length == 0) || (poly2->length == 0))
   {
      zmod_poly_zero(res);
	  return;
   }

   zmod_poly_t prod, quot;
   zmod_poly_init(prod, p);
   zmod_poly_init(quot, p);
   zmod_poly_mul(prod, poly1, poly2);
   zmod_poly_divrem(quot, res, prod, f);
   zmod_poly_clear(quot);
   zmod_poly_clear(prod);
}

/*
   Returns poly^exp modulo f
   Assumes poly is reduced mod f
   There are no restrictions on exp, i.e. it can be zero or negative
   The leading coefficient of f must be invertible modulo the modulus
*/

void zmod_poly_powmod(zmod_poly_t res, zmod_poly_t pol, long exp, zmod_poly_t f)
{
   zmod_poly_t y;
   
   unsigned long e;
   unsigned long p = f->p;

   if (pol->length == 0) 
   {
      if (exp <= 0L) 
	  {
         printf("FLINT Exception: Divide by zero\n");
         abort();   
	  }
	  zmod_poly_zero(res);
	  return;
   }

   if (exp < 0L)
      e = (unsigned long) -exp;
   else
      e = exp;
   
   if (exp) 
   {
	  zmod_poly_init(y, p);
	  zmod_poly_set(y, pol);
   }

	zmod_poly_zero(res);
	zmod_poly_set_coeff_ui(res, 0, 1L);
   res->length = 1;
   
   while (e) {
      if (e & 1) zmod_poly_mulmod(res, res, y, f);
      e = e >> 1;
	   if (e) zmod_poly_mulmod(y, y, y, f);
   }
   
   if (exp < 0L) zmod_poly_gcd_invert(res, res, f);
   if (exp) zmod_poly_clear(y);
} 

void zmod_poly_powpowmod(zmod_poly_t res, zmod_poly_t pol, ulong exp, ulong exp2, zmod_poly_t f)
{
	zmod_poly_t pow;
	zmod_poly_init(pow, f->p);
	zmod_poly_powmod(pow, pol, exp, f);
	zmod_poly_set(res, pow);
	
   if (!zmod_poly_equal(pow, pol)) 
		for (ulong i = 0; i < exp2 - 1; i++)
         zmod_poly_powmod(res, res, exp, f);

	zmod_poly_clear(pow);
}

/**************************************************************************************************

   Factorisation/Irreducibility

**************************************************************************************************/

int zmod_poly_isirreducible(zmod_poly_t f)
{
    //printf("input = "); zmod_poly_print(f); printf("\n");
    //see if it is non linear
    if (zmod_poly_length(f) > 2)
    {
        unsigned long p = zmod_poly_modulus(f);
        unsigned long n = zmod_poly_degree(f);
        //Some polynomials we will need
        
        zmod_poly_t a, x, x_p, x_modf, Q;
        zmod_poly_init(a, p);
        zmod_poly_init(x, p);
        zmod_poly_init(x_p, p);
        //Set up the constant polynomials
        zmod_poly_set_coeff_ui(x, 1, 1);
		  //compute x^q mod f
        zmod_poly_powpowmod(x_p, x, p, n, f);
        zmod_poly_make_monic(x_p, x_p);
		  //now do the irreducibility test
        if (!zmod_poly_equal(x_p, x))
		     return 0;
        else
        {
            factor_t factors;
            z_factor(&factors, n);
            for (unsigned long i = 0; i < factors.num; i++)
            {
               zmod_poly_powpowmod(a, x, p, n/factors.p[i], f);
               zmod_poly_sub(a, a, x);
               zmod_poly_make_monic(a, a);
					zmod_poly_gcd(a, a, f);
               
					if (a->length != 1)
                  return 0;
				}
        }

        zmod_poly_clear(a);
        zmod_poly_clear(x);
        zmod_poly_clear(x_p); 	
    }

    return 1;
}

/**
 * Initialises an array of zmod_poly's
 */
void zmod_poly_factor_init(zmod_poly_factor_t fac)
{
   fac->alloc = 5;
   fac->num_factors = 0;
   fac->factors = (zmod_poly_t *) flint_heap_alloc_bytes(sizeof(zmod_poly_t)*5);
   fac->exponents = (unsigned long *) flint_heap_alloc(5);
   for (unsigned long i = 0; i < 5; i++)
   {
	  fac->factors[i]->coeffs = (unsigned long*) flint_heap_alloc(1);
	  fac->factors[i]->alloc = 1;
      fac->factors[i]->length = 0;
   }
}

/** 
 * Frees up memory being used
 */
void zmod_poly_factor_clear(zmod_poly_factor_t fac)
{
	for (unsigned long i = 0; i < fac->alloc; i++)
	   zmod_poly_clear(fac->factors[i]);
	free(fac->factors);
	free(fac->exponents);
}

/**
 * Adds an extra element to the array
 */
void zmod_poly_factor_add(zmod_poly_factor_t fac, zmod_poly_t poly)
{
   if (poly->length <= 1) return;
   
   // how much space left in the array?, 
   // if none make a new one twice as big (for efficiency) and copy contents across
   if(fac->alloc == fac->num_factors)
   {
      fac->factors = (zmod_poly_t *) flint_heap_realloc_bytes(fac->factors, sizeof(zmod_poly_t)*2*fac->alloc);
	  fac->exponents = (unsigned long *) flint_heap_realloc(fac->exponents, 2*fac->alloc);
	  for (unsigned long i = fac->alloc; i < 2*fac->alloc; i++)
      {
	     fac->factors[i]->coeffs = (unsigned long*) flint_heap_alloc(1);
	     fac->factors[i]->alloc = 1;
         fac->factors[i]->length = 0;
	  }
      fac->alloc = 2*fac->alloc;
   } 
	
   zmod_poly_set(fac->factors[fac->num_factors], poly);
   fac->exponents[fac->num_factors] = 1;
   fac->num_factors++;
}

/**
 * Concatenates array res to res + fac
 */
void zmod_poly_factor_concat(zmod_poly_factor_t res, zmod_poly_factor_t fac)
{
   for(unsigned long i = 0; i < fac->num_factors; i++)
	{  
		zmod_poly_factor_add(res, fac->factors[i]);
		res->exponents[res->num_factors - 1] = fac->exponents[i];
	}
}

/**
 * Dumps the array to stdout
 */
void zmod_poly_factor_print(zmod_poly_factor_t fac)
{
   for(unsigned long i = 0; i < fac->num_factors; i++)
   {
      zmod_poly_print(fac->factors[i]); 
      printf(" ^ %ld\n", fac->exponents[i]);
   }	
}

/**
 * Raise the exponents of a factor array to the given exponent
 */
void zmod_poly_factor_pow(zmod_poly_factor_t fac, unsigned long exp)
{
   for (unsigned long i = 0; i < fac->num_factors; i++)
      fac->exponents[i] *= exp;
}

/** 
 * Square-Free Algorithm, takes an arbitary polynomial in F_p[X] and returns an array of square free factors 
 * LOW MULTIPLICITIES 
 */ 
void zmod_poly_factor_square_free(zmod_poly_factor_t res, zmod_poly_t f)
{
   if (f->length <= 1) 
	{
	   res->num_factors = 0;
       return;
	}

	if (f->length == 2)
	{
	   zmod_poly_factor_add(res, f);
	   return;
	}

   unsigned long p = zmod_poly_modulus(f);    //order of the field
   unsigned long deg = zmod_poly_degree(f); //degree of the polynomial
    
	//Step 1, look at f', if it is zero then we are done since f = h(x)^p
	//for some particular h(x), clearly f(x) = sum a_k x^kp, k <= deg(f)
   zmod_poly_t f_d, g, g_1;

   zmod_poly_init(g_1, p);
   zmod_poly_init(f_d, p);
   zmod_poly_init(g, p);

   zmod_poly_derivative(f_d, f);

   //CASE 1:
   if(zmod_poly_is_zero(f_d))
   {
        zmod_poly_t h;
        zmod_poly_init(h, p);
        for(unsigned long i = 0; i <= deg/p; ++i)    //this will be an integer since f'=0
        {
            zmod_poly_set_coeff_ui(h, i, zmod_poly_get_coeff_ui(f, i*p));
        }
        
		  //now run square-free on h, and return it to the pth power
        zmod_poly_factor_t new_res;
        zmod_poly_factor_init(new_res);

        zmod_poly_factor_square_free(new_res, h);
        //now raise it to the power of p
        zmod_poly_factor_pow(new_res, p);
			    
        zmod_poly_factor_concat(res, new_res);    //note, concatenating is equivalent to multiplying
		  zmod_poly_clear(h);
		  zmod_poly_factor_clear(new_res);
   }
   else 
   { 
        zmod_poly_gcd(g, f, f_d);
        zmod_poly_div(g_1, f, g);
        unsigned long i = 1;

        zmod_poly_t h, z;
        zmod_poly_init(h, p);
        zmod_poly_init(z, p);
        //CASE 2:
        while (!zmod_poly_is_one(g_1)) 
        {
            
			   zmod_poly_gcd(h, g_1, g);
            zmod_poly_div(z, g_1, h);
            
				// out <- out.z
            if (z->length > 1)
				{
					zmod_poly_factor_add(res, z);
               zmod_poly_make_monic(res->factors[res->num_factors - 1], res->factors[res->num_factors - 1]);
					if (res->num_factors) res->exponents[res->num_factors-1] *= i;
				}
			   i++;
				zmod_poly_set(g_1, h);
            zmod_poly_div(g, g, h);
        }
		  zmod_poly_clear(h);
		  zmod_poly_clear(z);
        
		  zmod_poly_make_monic(g, g);
		  if(!zmod_poly_is_one(g))
        {
            //so now we multiply res with square-free(g^1/p) ^ p 
            zmod_poly_t g_p; //g^(1/p)
            zmod_poly_init(g_p, p);

            for(unsigned long i = 0; i <= zmod_poly_degree(g)/p; i++)    
                zmod_poly_set_coeff_ui(g_p, i, zmod_poly_get_coeff_ui(g, i*p));
            
            zmod_poly_factor_t new_res_2;
            zmod_poly_factor_init(new_res_2);
            //square-free(g^(1/p))
            zmod_poly_factor_square_free(new_res_2, g_p);
            //now raise it to the power of p
            zmod_poly_factor_pow(new_res_2, p);
			   zmod_poly_factor_concat(res, new_res_2);
				zmod_poly_clear(g_p);
				zmod_poly_factor_clear(new_res_2);
        }
   }
	zmod_poly_clear(g_1);
   zmod_poly_clear(f_d);
   zmod_poly_clear(g);
}

/*
   Berlekamp polynomial factoring algorithm
   It accepts a polynomial f and either factors it or returns f itself as the only factor
   (in which case it is irreducible)
*/

void zmod_poly_factor_berlekamp(zmod_poly_factor_t factors, zmod_poly_t f)
{
	if (f->length <= 2)
	{
		zmod_poly_factor_add(factors, f);
		return;
	}
	
	unsigned long p = zmod_poly_modulus(f); //order of the field
	unsigned long n = zmod_poly_degree(f);

	//Step 1, we compute x^p mod f in F_p[X]/<f>
	zmod_poly_t x, x_p;
	zmod_poly_init(x, p);
	zmod_poly_init(x_p, p);
	
	zmod_poly_set_coeff_ui(x, 1, 1);
	zmod_poly_powmod(x_p, x, p, f);
	zmod_poly_clear(x);
	
	//Step 2, compute the matrix for the Berlekamp Map
	zmod_mat_t matrix;
	zmod_mat_init(matrix, p, n, n); 
	zmod_poly_t x_pi, x_pi2;
	zmod_poly_init(x_pi, p);
	zmod_poly_init(x_pi2, p);
	zmod_poly_set_coeff_ui(x_pi, 0, 1);
    ulong coeff;

	int reducible = 0;
	for (long i = 0; i < n; i++)
	{
	    //Q - I
		zmod_poly_set(x_pi2, x_pi);
		coeff = zmod_poly_get_coeff_ui(x_pi2, i);
		if (coeff) zmod_poly_set_coeff_ui(x_pi2, i, coeff - 1);
		else zmod_poly_set_coeff_ui(x_pi2, i, p - 1);
		zmod_poly_to_zmod_mat_col(matrix, i, x_pi2);
        zmod_poly_mulmod(x_pi, x_pi, x_p, f); 
	}
    zmod_poly_clear(x_p);
    zmod_poly_clear(x_pi);
    zmod_poly_clear(x_pi2);

	//Now we do Gauss-Jordan on Q-I
	//todo: try implementing back substitution instead since it is potentially faster
	ulong nullity = n - zmod_mat_row_reduce_gauss_jordan(matrix);
	
	//Try and find a basis for the nullspace
	zmod_poly_t * basis = (zmod_poly_t *) flint_heap_alloc(nullity * sizeof(zmod_poly_t));
    ulong * shift = (ulong *) flint_heap_alloc(n);
    F_mpn_clear(shift, n);

	ulong col = 1; // first column is always zero
	ulong row = 0;
	shift[0] = 1;
	for (ulong i = 1; i < nullity; i++)
	{
	   zmod_poly_init(basis[i], p);
	   while (zmod_mat_get_coeff_ui(matrix, row, col)) 
	   {
	      row++;
		  col++;
	   }
	   zmod_mat_col_to_zmod_poly_shifted(basis[i], matrix, col, shift);
	   zmod_poly_set_coeff_ui(basis[i], col, p - 1);
	   shift[col] = 1;
	   col++;
	}
	flint_heap_free(shift);
	zmod_mat_clear(matrix);

	if (nullity == 1) //we are done
	{
		zmod_poly_factor_add(factors, f);
		flint_heap_free(basis);
	} else
	{		
		//generate random linear combinations
		zmod_poly_t factor, b, power, g;
		zmod_poly_init(factor, p); zmod_poly_init(b, p); 
		zmod_poly_init(power, p); zmod_poly_init(g, p);
		    
		while (1)
		{
		   do {
		      zmod_poly_zero(factor);
			  for(ulong i = 1; i < nullity; i++)
			  {
			     zmod_poly_scalar_mul(b, basis[i], z_randint(p));
			     zmod_poly_add(factor, factor, b);
			  }
			  zmod_poly_set_coeff_ui(factor, 0, z_randint(p));
			  zmod_poly_make_monic(factor, factor);
		   } while(zmod_poly_is_one(factor) || zmod_poly_is_zero(factor));
		
		   zmod_poly_gcd(g, f, factor);
		   if (zmod_poly_length(g) != 1) break;
		   if (p > 3) zmod_poly_powmod(power, factor, p>>1, f);
		   else zmod_poly_set(power, factor);
		   power->coeffs[0] = z_addmod(power->coeffs[0], p - 1, p);
		   zmod_poly_gcd(g, power, f);
		   if (zmod_poly_length(g) != 1) break;
		}

		for (ulong i = 1; i < nullity; i++)
		{
		   zmod_poly_clear(basis[i]);
		}
		flint_heap_free(basis);
        zmod_poly_clear(power);
		zmod_poly_clear(factor);
		zmod_poly_clear(b);

		zmod_poly_make_monic(g, g);
		zmod_poly_factor_t fac1, fac2;
	    zmod_poly_factor_init(fac1); zmod_poly_factor_init(fac2);
		zmod_poly_factor_berlekamp(fac1, g);
		zmod_poly_t Q;
		zmod_poly_init(Q, p);
		zmod_poly_div(Q, f, g);
		zmod_poly_make_monic(Q, Q);
		zmod_poly_factor_berlekamp(fac2, Q);
		zmod_poly_factor_concat(factors, fac1);
	    zmod_poly_factor_concat(factors, fac2);
		zmod_poly_factor_clear(fac1);
		zmod_poly_factor_clear(fac2);
		zmod_poly_clear(Q);
		zmod_poly_clear(g);
	}			
}

/**
 * This function takes an arbitary polynomial and factorises it. It first 
 * performs a square-free factorisation, then factorises all of the square 
 * free polynomails and returns the leading coefficient, all the factors will 
 * be monic. If the zero polynomial is passed, 0 is return. If a constant is 
 * passed, that constant is returned (and no factors).
 */
unsigned long zmod_poly_factor(zmod_poly_factor_t result, zmod_poly_t input)
{
   if (input->length == 0) return 0;
	
	zmod_poly_t monic_input;
   zmod_poly_init(monic_input, zmod_poly_modulus(input));
   
	//Now we must make sure the input polynomial is monic. Get the highest coeff and store it then call make monic
   unsigned long leading_coeff = zmod_poly_get_coeff_ui(input, zmod_poly_degree(input));
   if (input->length == 1) return leading_coeff;

	zmod_poly_make_monic(monic_input, input);
	
   //This will store all of the square-free factors
   zmod_poly_factor_t sqfree_factors;
   zmod_poly_factor_init(sqfree_factors);
   
	//Now we do the square-free factorisation of the input polynomial
   zmod_poly_factor_square_free(sqfree_factors, monic_input);
	zmod_poly_clear(monic_input);
	
   //space to store factors
   zmod_poly_factor_t factors;
   
   //Run berlekamp on each of the square-free factors
   for (unsigned long i = 0; i < sqfree_factors->num_factors; i++)
   {
      zmod_poly_factor_init(factors);
      
      zmod_poly_factor_berlekamp(factors, sqfree_factors->factors[i]);
      zmod_poly_factor_pow(factors, sqfree_factors->exponents[i]);

      //Now add all the factors to the final array
      zmod_poly_factor_concat(result, factors);

      //clean up the memory being used
      zmod_poly_factor_clear(factors);
   }

   //clean up memory
   zmod_poly_factor_clear(sqfree_factors);
   	
   return leading_coeff;   
}

/**************************************************************************************************

   zmod_poly matrix routines

**************************************************************************************************/

void zmod_poly_2x2_mat_mul_classical(zmod_poly_2x2_mat_t R, zmod_poly_2x2_mat_t A, zmod_poly_2x2_mat_t B)
{
	zmod_poly_t temp;
	zmod_poly_init_precomp(temp, A->a->p, A->a->p_inv);
   
	zmod_poly_mul(R->a, A->a, B->a);
	zmod_poly_mul(temp, A->b, B->c);
	zmod_poly_add(R->a, R->a, temp);

   zmod_poly_mul(R->b, A->a, B->b);
	zmod_poly_mul(temp, A->b, B->d);
	zmod_poly_add(R->b, R->b, temp);

   zmod_poly_mul(R->c, A->c, B->a);
	zmod_poly_mul(temp, A->d, B->c);
	zmod_poly_add(R->c, R->c, temp);

   zmod_poly_mul(R->d, A->c, B->b);
	zmod_poly_mul(temp, A->d, B->d);
	zmod_poly_add(R->d, R->d, temp);

	zmod_poly_clear(temp);
}

void zmod_poly_2x2_mat_mul_strassen(zmod_poly_2x2_mat_t R, zmod_poly_2x2_mat_t A, zmod_poly_2x2_mat_t B)
{
	zmod_poly_t x0;
	zmod_poly_t x1;
	zmod_poly_init_precomp(x0, A->a->p, A->a->p_inv);
	zmod_poly_init_precomp(x1, A->a->p, A->a->p_inv);

   zmod_poly_sub(x0, A->a, A->c);
	zmod_poly_sub(x1, B->d, B->b);
	zmod_poly_mul(R->c, x0, x1);

	zmod_poly_add(x0, A->c, A->d);
	zmod_poly_sub(x1, B->b, B->a);
	zmod_poly_mul(R->d, x0, x1);

   zmod_poly_sub(x0, x0, A->a);
	zmod_poly_sub(x1, B->d, x1);
	zmod_poly_mul(R->b, x0, x1);

	zmod_poly_sub(x0, A->b, x0);
	zmod_poly_mul(R->a, x0, B->d);

	zmod_poly_mul(x0, A->a, B->a);

	zmod_poly_add(R->b, x0, R->b);
	zmod_poly_add(R->c, R->b, R->c);
	zmod_poly_add(R->b, R->b, R->d);
	zmod_poly_add(R->d, R->c, R->d);
	zmod_poly_add(R->b, R->b, R->a);
	zmod_poly_sub(x1, x1, B->c);
	zmod_poly_mul(R->a, A->d, x1);

	zmod_poly_sub(R->c, R->c, R->a);
	zmod_poly_mul(R->a, A->b, B->c);

	zmod_poly_add(R->a, R->a, x0);

	zmod_poly_clear(x0);
	zmod_poly_clear(x1);
}

#define ZMOD_POLY_2X2_STRASSEN_CUTOFF 20

void zmod_poly_2x2_mat_mul(zmod_poly_2x2_mat_t R, zmod_poly_2x2_mat_t A, 
									                         zmod_poly_2x2_mat_t B)
{
	ulong min_A = FLINT_MIN(FLINT_MIN(A->a->length, A->b->length), FLINT_MIN(A->c->length, A->d->length));
   ulong min_B = FLINT_MIN(FLINT_MIN(B->a->length, B->b->length), FLINT_MIN(B->c->length, B->d->length));

	if ((R == A) || (R == B))
	{
		zmod_poly_2x2_mat_t T;
		zmod_poly_2x2_mat_init(T, A->a->p);

		if (FLINT_MIN(min_A, min_B) < ZMOD_POLY_2X2_STRASSEN_CUTOFF)
			zmod_poly_2x2_mat_mul_classical(T, A, B);
		else
         zmod_poly_2x2_mat_mul_strassen(T, A, B);
		
		zmod_poly_swap(T->a, R->a);
      zmod_poly_swap(T->b, R->b);
      zmod_poly_swap(T->c, R->c);
      zmod_poly_swap(T->d, R->d);
      
		zmod_poly_2x2_mat_clear(T);
	} else
	{
		if (FLINT_MIN(min_A, min_B) < ZMOD_POLY_2X2_STRASSEN_CUTOFF)
			zmod_poly_2x2_mat_mul_classical(R, A, B);
		else
         zmod_poly_2x2_mat_mul_strassen(R, A, B);
	}
}

