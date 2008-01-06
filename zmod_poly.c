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

   zmod_poly.c: Polynomials over (unsigned) long mod p, for p prime.
   
   Copyright (C) 2007, David Howden.
   
*****************************************************************************/

#include "zmod_poly.h"
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
      zmod_poly_normalise(poly);
   }
}


void __zmod_poly_ensure_alloc(zmod_poly_t poly, unsigned long alloc)
{
   FLINT_ASSERT(alloc > poly->alloc);

   if (alloc < 2*poly->alloc)
      alloc = 2*poly->alloc;
   zmod_poly_realloc(poly, alloc);
}


/****************************************************************************

   Setting/retrieving coefficients

****************************************************************************/

void zmod_poly_set_coeff(zmod_poly_t poly, unsigned long n, unsigned long c)
{
   c = z_mod_precomp(c, poly->p, poly->p_inv);
   
   zmod_poly_ensure_alloc(poly, n+1);
   
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
         zmod_poly_normalise(poly);
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

   // jump to next whitespace
   s += strcspn(s, whitespace);
   
   poly->length = 0;
   zmod_poly_ensure_alloc(poly, length);
   
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
   
   zmod_poly_normalise(poly);
   
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
      // +2 is for the sign and a space
      size += (unsigned long)log10(poly->coeffs[i]) + 2;

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
   
   zmod_poly_ensure_alloc(poly, length);

   // read coefficients
   for (unsigned long i = 0; i < length; i++)
   {
      if (!fscanf(f, "%ld", &poly->coeffs[i]))
         return 0;
      poly->length++;
   }

   zmod_poly_normalise(poly);
   
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


void zmod_poly_normalise(zmod_poly_t poly)
{
   while (poly->length && (poly->coeffs[poly->length-1] == 0L))
      poly->length--;
}


int zmod_poly_normalised(zmod_poly_t poly)
{
   return (poly->length == 0) || (poly->coeffs[poly->length-1] != 0L);
}


void zmod_poly_pad(zmod_poly_t poly, unsigned long length)
{
   zmod_poly_ensure_alloc(poly, length);

   if (poly->length < length)
   {
      for (unsigned long i = poly->length; i < length; i++)
         poly->coeffs[i] = 0L;
      poly->length = length;
   }
}


void zmod_poly_truncate(zmod_poly_t res, zmod_poly_t poly, unsigned long length)
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
         zmod_poly_set(res, poly);
         return;
      }

      // todo: use mpz_init_set where appropriate
      
      zmod_poly_ensure_alloc(res, length);

      for (unsigned long i = 0; i < length; i++)
         res->coeffs[i] = poly->coeffs[i];
         
      res->length = length;
      
      res->p = poly->p;
      res->p_inv = poly->p_inv;
   }
   
   zmod_poly_normalise(res);
}



/****************************************************************************

   Assignment

****************************************************************************/


void zmod_poly_set(zmod_poly_t res, zmod_poly_t poly)
{
   if (res == poly)
      return;

   // todo: use mpz_init_set where appropriate

   zmod_poly_ensure_alloc(res, poly->length);
   
   for (unsigned long i = 0; i < poly->length; i++)
      res->coeffs[i] = poly->coeffs[i];
      
   res->length = poly->length;
   
   res->p = poly->p;
   res->p_inv = poly->p_inv;
}


/****************************************************************************

   Conversions

****************************************************************************/


// assumes coefficients are big enough, and alloc is big enough
// void _zmod_poly_to_fzmod_poly(fzmod_poly_t res, zmod_poly_t poly)
// {
//    FLINT_ASSERT(res->alloc >= poly->length);
// 
//    res->length = poly->length;
//    if (poly->length == 0)
//       return;
// 
//    for (unsigned long i = 0; i < poly->length; i++)
//    {
//       FLINT_ASSERT(res->limbs >= mpz_size(poly->coeffs[i]));
//       mpz_to_fmpz(res->coeffs + i*(res->limbs+1), poly->coeffs[i]);
//    }
// }
// 
// 
// void zmod_poly_to_fzmod_poly(fzmod_poly_t res, zmod_poly_t poly)
// {
//    unsigned long limbs = zmod_poly_max_limbs(poly);
// 
//    // todo: there should be a single function that achieves both of the
//    // following.... actually we don't even care in this case if the value
//    // is preserved.
//    fzmod_poly_fit_length(res, poly->length);
//    fzmod_poly_fit_limbs(res, limbs);
// 
//    _zmod_poly_to_fzmod_poly(res, poly);
// }
// 
// 
// void fzmod_poly_to_zmod_poly(zmod_poly_t res, fzmod_poly_t poly)
// {
//    zmod_poly_ensure_alloc(res, poly->length);
// 
//    res->length = poly->length;
//    
//    // todo: is there a bug here if poly->coeffs is not actually allocated?
// 
//    unsigned long i;
//    mp_limb_t* ptr = poly->coeffs;
// 
//    for (i = 0; i < poly->length; i++, ptr += poly->limbs+1)
//       fmpz_to_mpz(res->coeffs[i], ptr);
//    
//    zmod_poly_normalise(res);
// 
// }


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

   Addition/subtraction

****************************************************************************/


void zmod_poly_add(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   // rearrange parameters to make poly1 no longer than poly2
   if (poly1->length > poly2->length)
      SWAP_ZMOD_POLY_PTRS(poly1, poly2);
      
   zmod_poly_ensure_alloc(res, poly2->length);

   unsigned long i;
   
   for (i = 0; i < poly1->length; i++)
      res->coeffs[i] = z_mod_precomp(poly1->coeffs[i] + poly2->coeffs[i], poly1->p, poly1->p_inv);

   for (; i < poly2->length; i++)
      res->coeffs[i] = poly2->coeffs[i];

   res->length = poly2->length;
   zmod_poly_normalise(res);
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
   int swapped = 0;
   if (poly1->length > poly2->length)
   {
      swapped = 1;
      SWAP_ZMOD_POLY_PTRS(poly1, poly2);
   }
      
   zmod_poly_ensure_alloc(res, poly2->length);

   unsigned long i;
   
   if (swapped)
   {
      for (i = 0; i < poly1->length; i++)
      {
         if (poly2->coeffs[i] < poly1->coeffs[i])
         {
            res->coeffs[i] = poly2->p + poly2->coeffs[i] - poly1->coeffs[i];
         }
         else
         {
            res->coeffs[i] = poly2->coeffs[i] - poly1->coeffs[i];
         }
      }
         
      for (; i < poly2->length; i++)
         res->coeffs[i] = poly2->coeffs[i];
   }
   else
   {
      for (i = 0; i < poly1->length; i++)
      {
         if (poly1->coeffs[i] < poly2->coeffs[i])
         {
            res->coeffs[i] = poly2->p + poly1->coeffs[i] - poly2->coeffs[i];
         }
         else
         {
            res->coeffs[i] = poly1->coeffs[i] - poly2->coeffs[i];
         }
      }
         
      for (; i < poly2->length; i++)
         res->coeffs[i] = poly2->p - poly2->coeffs[i];
   }

   res->length = poly2->length;
   zmod_poly_normalise(res);
}


void zmod_poly_neg(zmod_poly_t res, zmod_poly_t poly)
{
   zmod_poly_ensure_alloc(res, poly->length);

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


void zmod_poly_lshift(zmod_poly_t res, zmod_poly_t poly, unsigned long k)
{
   zmod_poly_ensure_alloc(res, poly->length + k);

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
   }
   
   res->length = poly->length + k;
}


void zmod_poly_rshift(zmod_poly_t res, zmod_poly_t poly, unsigned long k)
{
   if (k >= poly->length)
   {
      // shift all coefficients off the end
      res->length = 0;
      res->p = poly->p;
      res->p_inv = poly->p_inv;
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
      zmod_poly_ensure_alloc(res, poly->length - k);

      for (unsigned long i = k; i < poly->length; i++)
         res->coeffs[i - k] = poly->coeffs[i];
         
      res->p = poly->p;
      res->p_inv = poly->p_inv;
   }
   
   res->length = poly->length - k;
}



/*******************************************************************************

   Polynomial multiplication

********************************************************************************/


void zmod_poly_mul(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   // use naive for now
   zmod_poly_mul_naive(res, poly1, poly2);
}


void zmod_poly_sqr(zmod_poly_t res, zmod_poly_t poly)
{
   // use naive for now
   zmod_poly_sqr_naive(res, poly);
}


/*
 This is just like zmod_poly_mul_naive(), with the following restrictions:

  * assumes res does not alias poly1 and poly2
  * neither polynomial is zero
  * res->alloc >= poly1->length + poly2->length - 1
     (i.e. output has enough room for product)
*/

void _zmod_poly_mul_naive(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   FLINT_ASSERT(res != poly1);
   FLINT_ASSERT(res != poly2);
   FLINT_ASSERT(poly1->length && poly2->length);

   res->length = poly1->length + poly2->length - 1;
   res->p = poly1->p;
   res->p_inv = poly1->p_inv;
   
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
      __zmod_poly_mul_naive_mod_last(res, poly1, poly2);
   }
   else
   {
      bits = zmod_poly_bits(poly1) + zmod_poly_bits(poly2) + log_length;
      if (bits < FLINT_BITS)
      {
         __zmod_poly_mul_naive_mod_last(res, poly1, poly2);
      }
      else
      {
         __zmod_poly_mul_naive_mod_throughout(res, poly1, poly2);
      }
   }
}

void _zmod_poly_mul_naive_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   FLINT_ASSERT(res != poly1);
   FLINT_ASSERT(res != poly2);
   FLINT_ASSERT(poly1->length && poly2->length);

   res->length = poly1->length + poly2->length - 1;
   res->p = poly1->p;
   res->p_inv = poly1->p_inv;

   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = 0;

   __zmod_poly_mul_naive_mod_throughout(res, poly1, poly2);
}


/*
   Actually computes the naive multiplication, only applying mod at the end
   of the computations.
*/

void __zmod_poly_mul_naive_mod_last(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   for (unsigned long i = 0; i < poly1->length; i++)
      for (unsigned long j = 0; j < poly2->length; j++)
         res->coeffs[i+j] = res->coeffs[i+j] + poly1->coeffs[i] * poly2->coeffs[j];
         
   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = z_mod_precomp(res->coeffs[i], res->p, res->p_inv);
}


/*
   Computes the naive multiplication, applying mods at each step.
*/

void __zmod_poly_mul_naive_mod_throughout(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
{
   for (unsigned long i = 0; i < poly1->length; i++)
      for (unsigned long j = 0; j < poly2->length; j++)
         res->coeffs[i+j] = z_mod_precomp(res->coeffs[i+j] + z_mulmod_precomp(poly1->coeffs[i], poly2->coeffs[j], poly1->p, poly1->p_inv), poly1->p, poly1->p_inv);
}


void zmod_poly_mul_naive(zmod_poly_t res, zmod_poly_t poly1, zmod_poly_t poly2)
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
      zmod_poly_sqr_naive(res, poly1);
      return;
   }

   unsigned long length = poly1->length + poly2->length - 1;

   if (res == poly1 || res == poly2)
   {
      // output is inplace, so need a temporary
      zmod_poly_t temp;
      zmod_poly_init2(temp, poly1->p, length);
      _zmod_poly_mul_naive(temp, poly1, poly2);
      zmod_poly_swap(temp, res);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace
      zmod_poly_ensure_alloc(res, length);
      _zmod_poly_mul_naive(res, poly1, poly2);
   }
}


/*
 This is just like zmod_poly_sqr_naive(), with the following restrictions:

  * assumes res does not alias poly
  * poly is nonzero
  * res->alloc >= 2*poly->length - 1  (i.e. output has enough room for product)
*/
void _zmod_poly_sqr_naive(zmod_poly_t res, zmod_poly_t poly)
{
   FLINT_ASSERT(res != poly);
   FLINT_ASSERT(poly->length);

   res->length = 2*poly->length - 1;
   res->p = poly->p;
   res->p_inv = poly->p_inv;
   FLINT_ASSERT(res->alloc >= res->length);

   for (unsigned long i = 0; i < res->length; i++)
      res->coeffs[i] = 0;

   // off-diagonal products
   for (unsigned long i = 1; i < poly->length; i++)
      for (unsigned long j = 0; j < i; j++)
         res->coeffs[i+j] = z_mod_precomp(res->coeffs[i+j] + z_mulmod_precomp(poly->coeffs[i], poly->coeffs[j], poly->p, poly->p_inv), poly->p, poly->p_inv);

   // double the off-diagonal products
   for (unsigned long i = 1; i < res->length - 1; i++)
      res->coeffs[i] = z_mod_precomp(2*res->coeffs[i], poly->p, poly->p_inv);

   // add in diagonal products
   for (unsigned long i = 0; i < poly->length; i++)
      res->coeffs[2*i] = z_mod_precomp(res->coeffs[2*i] + z_mulmod_precomp(poly->coeffs[i], poly->coeffs[i], poly->p, poly->p_inv), poly->p, poly->p_inv);
}


void zmod_poly_sqr_naive(zmod_poly_t res, zmod_poly_t poly)
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
      _zmod_poly_sqr_naive(temp, poly);
      zmod_poly_swap(temp, res);
      zmod_poly_clear(temp);
   }
   else
   {
      // output not inplace

      // allocate more coefficients if necessary
      zmod_poly_ensure_alloc(res, length);
      _zmod_poly_sqr_naive(res, poly);
   }
}

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
   
   unsigned long final_length = length1 + length2 - 1;
   
   while ((length1) && (input1->coeffs[length1-1] == 0)) length1--;
   while ((length2) && (input2->coeffs[length2-1] == 0)) length2--;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      zmod_poly_zero(output);
      return;
   }

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
   
   // else
   // {
   //    printf("Bits: %d\n", bits);
   // }
   
   //print_var("bits", bits);
   // if (bits > FLINT_BITS)
   // {
   //    printf("Cannot multiply this\n");
   //    zmod_poly_zero(output);
   //    return;
   // }

   mp_limb_t *mpn1, *mpn2, *res;

   unsigned long limbs1, limbs2;

   limbs1 = (unsigned long) ((input1->length * bits-1) / FLINT_BITS + 1);
   limbs2 = (unsigned long) ((input2->length * bits-1) / FLINT_BITS + 1);
   
   mpn1 = (mp_limb_t*) flint_stack_alloc(limbs1);
   mpn2 = (input1 == input2) ? mpn1 : (mp_limb_t*) flint_stack_alloc(limbs2);

   zmod_poly_bit_pack_mpn(mpn1, input1, bits);
   
   // print_var("limbs1", limbs1);
   //  
   //  for (unsigned long i = 0; i < limbs1; i++)
   //  {
   //     print_binary(mpn1[i], FLINT_BITS);
   //     printf(" ");
   //  }
   //  printf("\n");
   
   if(input1 != input2)
      zmod_poly_bit_pack_mpn(mpn2, input2, bits);
   
   limbs1 = FLINT_MAX((long)((length1 * bits-1) / FLINT_BITS + 1), 0L);
   limbs2 = FLINT_MAX((long)((length2 * bits-1) / FLINT_BITS + 1), 0L);
   
   res = (mp_limb_t*) flint_stack_alloc(limbs1+limbs2);
   res[limbs1+limbs2-1] = 0L;
   
   F_mpn_mul(res, mpn1, limbs1, mpn2, limbs2);
   
   zmod_poly_ensure_alloc(output, final_length);
   
   zmod_poly_bit_unpack_mpn(output, res, length1 + length2 - 1, bits);
   
   for (unsigned long i = length1 + length2 - 1; i < final_length; i++)
   {
      _zmod_poly_set_coeff(output, i, 0L);
   } 
   
   flint_stack_release();
   flint_stack_release();
   if(input1 != input2)
      flint_stack_release();

   
   output->length = final_length;
}


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

void zmod_poly_bit_pack_mpn(mp_limb_t * res, zmod_poly_t poly, unsigned long bits)
{  
   unsigned long current_limb = 0;
   unsigned int current_bit = 0;
   
   unsigned long temp_lower;
   unsigned long temp_upper;
   
   unsigned long total_limbs = FLINT_MAX((long)(((poly->length * bits - 1)>>FLINT_LG_BITS_PER_LIMB) + 1), 0L);
   
   res[0] = 0L;
   
   if (bits < FLINT_BITS)
   {
      unsigned long boundary_limit_bit = FLINT_BITS - bits;

      //printf("Packing polynomial ****************************************\n");
      //print_limb("res[0]", res[0]);

      for(unsigned long i = 0; i < poly->length; i++)
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
      for (unsigned long i = 0; i < poly->length; i++)
      {
         res[i] = poly->coeffs[i];
      }
   }
   else if (bits == 2*FLINT_BITS)
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         res[current_limb] = poly->coeffs[i];
         current_limb++;
         res[current_limb] = 0L;
         current_limb++;
      }
   }
   else
   {
      //printf("Packing Coeffs in Poly =============================================");
      
      for(unsigned long i = 0; i < poly->length; i++)
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
   }
}


/*
   Unpacks a zmod_poly of length `length` from an mpn `mpn` with coeffs packed in `bits` bits.
*/

void zmod_poly_bit_unpack_mpn(zmod_poly_t res, mp_limb_t * mpn, unsigned long length, unsigned long bits)
{
   unsigned long i;
   
   //PRINT_VAR(bits);
   
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
             _zmod_poly_set_coeff(res, i, z_mod2_precomp(temp_upper, res->p, res->p_inv));
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

             //zmod_poly_set_coeff(res, i, temp_lower);
             _zmod_poly_set_coeff(res, i, z_mod2_precomp(temp_lower, res->p, res->p_inv));

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
         _zmod_poly_set_coeff(res, i, z_ll_mod_precomp(0L, mpn[i], res->p, res->p_inv));
      }
   }
   else if (bits == 2*FLINT_BITS)
   {
      unsigned long current_limb = 0;
      for (i = 0; i < length; i++)
      {
         _zmod_poly_set_coeff(res, i, z_ll_mod_precomp(mpn[current_limb+1], mpn[current_limb], res->p, res->p_inv));
         current_limb+=2;
      }
   }
   else  // FLINT_BITS < bits < 2*FLINT_BITS
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
            // printf("Coeff accross one boundary... current_bit == 0\n");
            temp_lower = mpn[current_limb];
            // PRINT_LIMB(temp_lower);
            current_limb++;
            temp_upper = (mpn[current_limb] << (2*FLINT_BITS - bits)) >> (2*FLINT_BITS - bits);
            // PRINT_LIMB(temp_upper);
            _zmod_poly_set_coeff(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
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
            _zmod_poly_set_coeff(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
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

            _zmod_poly_set_coeff(res, i, z_ll_mod_precomp(temp_upper, temp_lower, res->p, res->p_inv));
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
   }
   //printf("Unpacking poly **********************************************\n");
}
