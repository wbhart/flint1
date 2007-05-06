/****************************************************************************

ZmodFpoly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "flint-manager.h"
#include "ZmodFpoly.h"


/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFpoly_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;

   // todo: change these to use the FLINT heap allocator
   
   poly->storage = (mp_limb_t*) malloc(bufs * (n+1) * sizeof(mp_limb_t));

   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) malloc(bufs * sizeof(ZmodF_t*));
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   for (unsigned long i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}


void ZmodFpoly_clear(ZmodFpoly_t poly)
{
   free(poly->coeffs);
   free(poly->storage);
}



/****************************************************************************

   Conversion Routines
   
****************************************************************************/

void ZmodFpoly_convert_in_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   abort();
}

void ZmodFpoly_convert_out_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   abort();
}

void ZmodFpoly_bit_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                            unsigned long bundle, unsigned long bits)
{
   abort();
}

void ZmodFpoly_bit_unpack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                              unsigned long bundle, unsigned long bits)
{
   abort();
}
     
void ZmodFpoly_byte_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                             unsigned long bundle, unsigned long bytes)
{
   abort();
}
     
void ZmodFpoly_byte_unpack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                               unsigned long bundle, unsigned long bytes)
{
   abort();
}

     
void ZmodFpoly_split_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                         unsigned long bundle, unsigned long limbs)
{
   abort();
}
     
void ZmodFpoly_unsplit_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                           unsigned long bundle, unsigned long limbs)
{
   abort();
}
     

/****************************************************************************

   Basic Arithmetic Routines
   
****************************************************************************/


void ZmodFpoly_set(ZmodFpoly_t x, ZmodFpoly_t y)
{
   abort();
}


void ZmodFpoly_mul(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   abort();
}


void ZmodFpoly_add(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   abort();
}


void ZmodFpoly_sub(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   abort();
}


void ZmodFpoly_normalise(ZmodFpoly_t poly)
{
   abort();
}


/****************************************************************************

   Fourier Transform Routines

****************************************************************************/


/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of buffers assumed to be nonzero
length = number of fourier coefficients requested
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_FFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                    unsigned long nonzero, unsigned long length,
                    unsigned long twist, unsigned long n,
                    ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> m;
   FLINT_ASSERT(twist < root);
   
   // ========================
   // base cases
   
   if (depth == 0)
      return;
      
   if (depth == 1)
   {
      if (length == 1)
      {
         if (nonzero == 2)
            ZmodF_add(x[0], x[0], x[skip], n);
      }
      else   // length == 2
      {
         if (nonzero == 1)
            ZmodF_mul_sqrt2exp(x+skip, x, scratch, twist, n);
         else  // nonzero == 2
            ZmodF_forward_butterfly_sqrt2exp(x, x+skip, scratch, twist, n);
      }
      
      return;
   }
   
   // ========================
   // factoring case

   unsigned long rows_depth = depth >> 1;
   unsigned long cols_depth = depth - rows_depth;
   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long length_whole_rows = length_cols ?
                                     (length_rows + 1) : length_rows;
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // column transforms
   for (i = 0, y = x, j = twist; i < nonzero_cols; i++, y += skip, j += root)
      _ZmodFpoly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                     length_whole_rows, j, n, scratch);

   if (nonzero_rows)
   {
      for (; i < cols; i++, y += skip, j += root)
         _ZmodFpoly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                        length_whole_rows, j, n, scratch);
      nonzero_cols = cols;
   }
   
   // row transforms
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, cols,
                     twist << cols_depth, n, scratch);

   if (length_cols)
      // The relevant portion of the last row:
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, length_cols,
                     twist << cols_depth, n, scratch);
}



void ZmodFpoly_FFT(ZmodFpoly_t poly, unsigned long length)
{
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         for (unsigned long i = 0; i < length; i++)
            ZmodF_clear(poly->coeffs[i], poly->n);
      }
      else
      {
         _ZmodFpoly_FFT(poly->coeffs, poly->depth, 1, poly->length,
                        length, 0, poly->n, poly->scratch);
      }
   }

   poly->length = length;
}


void ZmodFpoly_IFFT(ZmodFpoly_t poly)
{
   abort();
}


void ZmodFpoly_rescale(ZmodFpoly_t poly)
{
   abort();
}


void ZmodFpoly_convolution(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   abort();
}


/****************************************************************************

   Negacyclic Fourier Transform Routines
   
****************************************************************************/


void ZmodFpoly_negacyclic_FFT(ZmodFpoly_t poly, unsigned long length)
{
   abort();
}


void ZmodFpoly_negacyclic_IFFT(ZmodFpoly_t poly)
{
   abort();
}


void ZmodFpoly_negacyclic_convolution(ZmodFpoly_t res,
                                      ZmodFpoly_t x, ZmodFpoly_t y)
{
   abort();
}


// end of file ****************************************************************
