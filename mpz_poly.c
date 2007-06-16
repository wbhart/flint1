/******************************************************************************

mpz_poly.c: Polynomials over Z, implemented as an array of mpz_t's

Copyright (C) 2007, William Hart and David Harvey

******************************************************************************/

#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "mpz_poly.h"


/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpz_poly_init(mpz_poly_t poly)
{
   poly->coeffs = (mpz_t*) flint_heap_alloc(sizeof(mpz_t));
   poly->alloc = 1;
   poly->init = 0;
   poly->length = 0;
}


void mpz_poly_clear(mpz_poly_t poly)
{
   for (unsigned long i = 0; i < poly->init; i++)
      mpz_clear(poly->coeffs[i]);

   flint_heap_free(poly->coeffs);
}


void mpz_poly_realloc(mpz_poly_t poly, unsigned long alloc)
{
   FLINT_ASSERT(alloc >= 1);

   // clear any mpz_t's beyond the new array length
   if (poly->init > alloc)
   {
      for (unsigned long i = alloc; i < poly->init; i++)
         mpz_clear(poly->coeffs[i]);
      poly->init = alloc;
   }

   poly->alloc = alloc;
   poly->coeffs = (mpz_t*) flint_heap_realloc(poly->coeffs,
                                              alloc * sizeof(mpz_t));
   
   // truncate poly if necessary
   if (poly->length > alloc)
   {
      poly->length = alloc;
      mpz_poly_normalise(poly);
   }
}


void mpz_poly_ensure_alloc(mpz_poly_t poly, unsigned long alloc)
{
   if (poly->alloc < alloc)
   {
      if (alloc < 2*poly->alloc)
         alloc = 2*poly->alloc;
      mpz_poly_realloc(poly, alloc);
   }
}


void mpz_poly_init_upto(mpz_poly_t poly, unsigned long init)
{
   mpz_poly_ensure_alloc(poly, init);

   if (poly->init < init)
   {
      for (unsigned long i = poly->init; i < init; i++)
         mpz_init(poly->coeffs[i]);
      poly->init = init;
   }
}


/****************************************************************************

   Setting/retrieving coefficients

****************************************************************************/


mpz_t* mpz_poly_get_coeff_ptr(mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return NULL;
   return &poly->coeffs[n];
}


void mpz_poly_get_coeff(mpz_t c, mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(c, 0);
   mpz_set(c, poly->coeffs[n]);
}


unsigned long mpz_poly_get_coeff_ui(mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return 0;
   return mpz_get_ui(poly->coeffs[n]);
}


long mpz_poly_get_coeff_si(mpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      return 0;
   return mpz_get_si(poly->coeffs[n]);
}


void mpz_poly_set_coeff(mpz_poly_t poly, unsigned long n, mpz_t c)
{
   if (n+1 <= poly->length)
   {
      if ((n+1 == poly->length) && !mpz_sgn(c))
      {
         // set last coefficient to zero (and normalise)
         do poly->length--;
         while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]));
      }
      else
         // set an existing coefficient
         mpz_set(poly->coeffs[n], c);
   }
   else
   {
      if (!mpz_sgn(c))
         // set zero coefficient beyond current length
         return;

      // zero extend, possibly init, and set new coefficient
      mpz_poly_ensure_alloc(poly, n+1);
      if (n < poly->init)
         mpz_set(poly->coeffs[n], c);
      else
      {
         for (; poly->init < n; poly->init++)
            mpz_init(poly->coeffs[poly->init]);
         mpz_init_set(poly->coeffs[n], c);
      }
      poly->length = n+1;
   }
}


void mpz_poly_set_coeff_ui(mpz_poly_t poly, unsigned long n, unsigned long c)
{
   if (n+1 <= poly->length)
   {
      if ((n+1 == poly->length) && !c)
      {
         // set last coefficient to zero (and normalise)
         do poly->length--;
         while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]));
      }
      else
         // set an existing coefficient
         mpz_set_ui(poly->coeffs[n], c);
   }
   else
   {
      if (!c)
         // set zero coefficient beyond current length
         return;

      // zero extend, possibly init, and set new coefficient
      mpz_poly_init_upto(poly, n+1);
      mpz_set_ui(poly->coeffs[n], c);
      poly->length = n+1;
   }
}


void mpz_poly_set_coeff_si(mpz_poly_t poly, unsigned long n, long c)
{
   if (n+1 <= poly->length)
   {
      if ((n+1 == poly->length) && !c)
      {
         // set last coefficient to zero (and normalise)
         do poly->length--;
         while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]));
      }
      else
         // set an existing coefficient
         mpz_set_si(poly->coeffs[n], c);
   }
   else
   {
      if (!c)
         // set zero coefficient beyond current length
         return;

      // zero extend, possibly init, and set new coefficient
      mpz_poly_init_upto(poly, n+1);
      mpz_set_si(poly->coeffs[n], c);
      poly->length = n+1;
   }
}



/****************************************************************************

   String conversions and I/O

****************************************************************************/


int mpz_poly_from_string(mpz_poly_t poly, char* s)
{
   const char* whitespace = " \t\n\r";
   
   // read poly length
   unsigned long length;
   if (!sscanf(s, "%d", &length))
      return 0;

   // jump to next whitespace
   s += strcspn(s, whitespace);
   
   poly->length = 0;
   mpz_poly_init_upto(poly, length);

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
   char* ptr = buf + sprintf(buf, "%d  ", poly->length);
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


void mpz_poly_fprint(mpz_poly_t poly, FILE* f)
{
   char* s = mpz_poly_to_string(poly);
   fputs(s, f);
   free(s);
}


void mpz_poly_print(mpz_poly_t poly)
{
   mpz_poly_fprint(poly, stdout);
}


int mpz_poly_fread(mpz_poly_t poly, FILE* f)
{
   // read poly length
   unsigned long length;
   if (!fscanf(f, "%d", &length))
      return 0;

   poly->length = 0;
   mpz_poly_init_upto(poly, length);

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
   mpz_poly_init_upto(poly, length);
   if (poly->length < length)
   {
      for (unsigned long i = poly->length; i < length; i++)
         mpz_set_ui(poly->coeffs[i], 0);
      poly->length = length;
   }
}



/****************************************************************************

   Assignment

****************************************************************************/


void mpz_poly_set(mpz_poly_t res, mpz_poly_t poly)
{
   mpz_poly_ensure_alloc(res, poly->length);
   
   // copy into coefficients that are already mpz_init'd
   unsigned long i, n = FLINT_MIN(poly->length, res->init);
   for (i = 0; i < n; i++)
      mpz_set(res->coeffs[i], poly->coeffs[i]);
      
   // copy into coefficients that need to be mpz_init'd
   if (i < poly->length)
   {
      for (; i < poly->length; i++)
         mpz_init_set(res->coeffs[i], poly->coeffs[i]);
      res->init = poly->length;
   }
   
   res->length = poly->length;
}



/****************************************************************************

   Comparison

****************************************************************************/


int mpz_poly_equal(mpz_poly_p poly1, mpz_poly_p poly2)
{
   if (poly1->length != poly2->length)
      return 0;

   for (unsigned long i = 0; i < poly1->length; i++)
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
   
   // first handle additions where target is already mpz_init'd
   unsigned long i, n = FLINT_MIN(poly1->length, res->init);
   for (i = 0; i < n; i++)
      mpz_add(res->coeffs[i], poly1->coeffs[i], poly2->coeffs[i]);

   // now handle additions where target is not yet mpz_init'd
   for (; i < poly1->length; i++)
   {
      mpz_t* x = poly1->coeffs + i;
      mpz_t* y = poly2->coeffs + i;
      // just take max of limbs; occasionally this will be too small
      unsigned long limbs = FLINT_MAX(mpz_size(*x), mpz_size(*y));
      mpz_init2(res->coeffs[i], FLINT_BITS * limbs);
      mpz_add(res->coeffs[i], *x, *y);
   }
   
   // now handle additions where target is already mpz_init'd, and one
   // input is zero
   n = FLINT_MIN(poly2->length, res->init);
   for (; i < n; i++)
      mpz_set(res->coeffs[i], poly2->coeffs[i]);

   // finally handle additions where target is not yet mpz_init'd, and one
   // input is zero
   for (; i < poly2->length; i++)
      mpz_init_set(res->coeffs[i], poly2->coeffs[i]);
   
   // update init and length
   res->length = poly2->length;
   if (res->init < poly2->length)
      res->init = poly2->length;

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

   unsigned long shorter, longer;

   if (poly1->length < poly2->length)
   {
      shorter = poly1->length;
      longer = poly2->length;
   }
   else
   {
      shorter = poly2->length;
      longer = poly1->length;
   }
   
   mpz_poly_ensure_alloc(res, longer);
   
   // first handle subtractions where target is already mpz_init'd
   unsigned long i, n = FLINT_MIN(shorter, res->init);
   for (i = 0; i < n; i++)
      mpz_sub(res->coeffs[i], poly1->coeffs[i], poly2->coeffs[i]);

   // now handle subtractions where target is not yet mpz_init'd
   for (; i < shorter; i++)
   {
      mpz_t* x = poly1->coeffs + i;
      mpz_t* y = poly2->coeffs + i;
      // just take max of limbs; occasionally this will be too small
      unsigned long limbs = FLINT_MAX(mpz_size(*x), mpz_size(*y));
      mpz_init2(res->coeffs[i], FLINT_BITS * limbs);
      mpz_sub(res->coeffs[i], *x, *y);
   }

   // now handle subtractions where one of the inputs is zero
   if (poly1->length <= poly2->length)
   {
      // target is already mpz_init'd
      n = FLINT_MIN(longer, res->init);
      for (; i < n; i++)
         mpz_neg(res->coeffs[i], poly2->coeffs[i]);

      // target is not yet mpz_init'd
      for (; i < longer; i++)
      {
         mpz_init_set(res->coeffs[i], poly2->coeffs[i]);
         mpz_neg(res->coeffs[i], res->coeffs[i]);
      }
   }
   else
   {
      // target is already mpz_init'd
      n = FLINT_MIN(longer, res->init);
      for (; i < n; i++)
         mpz_set(res->coeffs[i], poly1->coeffs[i]);

      // target is not yet mpz_init'd
      for (; i < longer; i++)
         mpz_init_set(res->coeffs[i], poly1->coeffs[i]);
   }
   
   // update init and length
   res->length = longer;
   if (res->init < longer)
      res->init = longer;

   mpz_poly_normalise(res);
}



void mpz_poly_neg(mpz_poly_t res, mpz_poly_t poly)
{
   if (poly == res)
   {
      // inplace case
      for (unsigned long i = 0; i < poly->length; i++)
         mpz_neg(poly->coeffs[i], poly->coeffs[i]);

      return;
   }
   
   // not inplace
   mpz_poly_ensure_alloc(res, poly->length);
   
   // first handle coefficients which are already mpz_init'd
   unsigned long i, n = FLINT_MIN(poly->length, res->init);
   for (i = 0; i < n; i++)
      mpz_neg(res->coeffs[i], poly->coeffs[i]);
      
   // copy into coefficients that need to be mpz_init'd
   if (i < poly->length)
   {
      for (; i < poly->length; i++)
      {
         mpz_init_set(res->coeffs[i], poly->coeffs[i]);
         mpz_neg(res->coeffs[i], res->coeffs[i]);
      }
      res->init = poly->length;
   }
   
   res->length = poly->length;
}



/****************************************************************************

   Shifting

****************************************************************************/


void mpz_poly_lshift(mpz_poly_t res, mpz_poly_t poly, unsigned long k)
{
   if (poly == res)
   {
      // inplace; just shift the mpz_t's over
      mpz_poly_init_upto(poly, poly->length + k);

      // todo: would probably be better to copy the data in large chunks,
      //       with memmove or something. Need to be careful though.....
      
      for (long i = poly->length - 1; i >= 0; i--)
         mpz_swap(poly->coeffs[i], poly->coeffs[i+k]);
      
      for (unsigned long i = 0; i < k; i++)
         mpz_set_ui(poly->coeffs[i], 0);
   }
   else
   {
      // not inplace; need to copy data
      mpz_poly_ensure_alloc(res, poly->length + k);

      // put zeroes at the bottom
      mpz_poly_init_upto(res, k);
      for (unsigned long i = 0; i < k; i++)
         mpz_set_ui(poly->coeffs[i], 0);
      
      // copy into coefficients that are already mpz_init'd
      unsigned long i, n = FLINT_MIN(poly->length, res->init - k);
      for (i = 0; i < n; i++)
         mpz_set(res->coeffs[i + k], poly->coeffs[i]);
         
      // copy into coefficients that need to be mpz_init'd
      if (i < poly->length)
      {
         for (; i < poly->length; i++)
            mpz_init_set(res->coeffs[i + k], poly->coeffs[i]);
         res->init = poly->length + k;
      }
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

      // todo: would probably be better to copy the data in large chunks,
      //       with memmove or something. Need to be careful though.....

      for (unsigned long i = k; i < poly->length; i++)
         mpz_swap(poly->coeffs[i-k], poly->coeffs[i]);
   }
   else
   {
      // not inplace; need to copy data
      mpz_poly_ensure_alloc(res, poly->length - k);

      // copy into coefficients that are already mpz_init'd
      unsigned long i, n = FLINT_MIN(poly->length, res->init + k);
      for (i = k; i < n; i++)
         mpz_set(res->coeffs[i - k], poly->coeffs[i]);
         
      // copy into coefficients that need to be mpz_init'd
      if (i < poly->length)
      {
         for (; i < poly->length; i++)
            mpz_init_set(res->coeffs[i - k], poly->coeffs[i]);
         res->init = poly->length - k;
      }
   }
   
   res->length = poly->length - k;
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

void mpz_poly_mul_naive(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_mul_karatsuba(mpz_poly_t res, mpz_poly_t poly1,
                            mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_mul_SS(mpz_poly_t res, mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_sqr_naive(mpz_poly_t res, mpz_poly_t poly)
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
   mpz_poly_init_upto(res, out_len);
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
   mpz_poly_init_upto(res, out_len);
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
      mpz_poly_init_upto(Q2, k2+1);
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

void mpz_poly_monic_inverse_naive(mpz_poly_t res, mpz_poly_t poly,
                                  unsigned long k)
{
   abort();
}

void mpz_poly_pseudo_inverse_naive(mpz_poly_t res, mpz_poly_t poly,
                                   unsigned long k)
{
   abort();
}

void mpz_poly_monic_div_naive(mpz_poly_t quot, mpz_poly_t poly1,
                              mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_div_naive(mpz_poly_t quot, mpz_poly_t poly1,
                               mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_rem_naive(mpz_poly_t rem, mpz_poly_t poly1,
                              mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_rem_naive(mpz_poly_t rem, mpz_poly_t poly1,
                               mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_monic_div_rem_naive(mpz_poly_t quot, mpz_poly_t rem,
                                  mpz_poly_t poly1, mpz_poly_t poly2)
{
   abort();
}

void mpz_poly_pseudo_div_rem_naive(mpz_poly_t quot, mpz_poly_t rem, 
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



// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//
// ======================== old code that hasn't been reviewed yet !!!
//
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

#if 0
                           
/* Naieve schoolboy polynomial multiplication routine */

void _mpz_poly_mul_naive(mpz_poly_t output, mpz_poly_t input1, mpz_poly_t input2)
{
   FLINT_ASSERT(output != input1);
   FLINT_ASSERT(output != input2);

   if (!input1->length || !input2->length)
   {
      // one of the inputs is zero
      output->length = 0;
      return;
   }
   
   output->length = input1->length + input2->length - 1;
   FLINT_ASSERT(output->alloc >= output->length);

   for (unsigned long i = 0; i < output->length; i++)
      mpz_set_ui(output->coeffs[i], 0);
   
   for (unsigned long i = 0; i < input1->length; i++)
      for (unsigned long j = 0; j < input2->length; j++)
         mpz_addmul(output->coeffs[i+j], input1->coeffs[i], input2->coeffs[j]);
}


#endif

// *************** end of file
