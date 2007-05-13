/****************************************************************************

ZmodFpoly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <sys/types.h>
#include "flint.h"
#include "flint-manager.h"
#include "ZmodFpoly.h"
#include "Zpoly_mpn.h"
#include "extras.h"


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
   unsigned long size_f = poly_f->n + 1;
   unsigned long size_m = poly_mpn->limbs+1;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;
   
   long size_j;
   
   for (unsigned long i = 0, j = 0; i < poly_mpn->length; i++, j += size_m)
   {
      if ((size_j = coeffs_m[j]) < 0)
      {
         negate_limbs(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         set_limbs(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      } else
      {
         copy_limbs(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         clear_limbs(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      }
   }
   poly_f->length = poly_mpn->length;   
}

void ZmodFpoly_convert_out_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f)
{
   unsigned long n = poly_f->n;
   unsigned long size_m = poly_mpn->limbs+1;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;

   for (unsigned long i = 0, j = 0; i < poly_f->length; i++, j += size_m)
   {
      ZmodF_normalise(coeffs_f[i], n);
      if (coeffs_f[i][n-1]>>(FLINT_BITS_PER_LIMB-1) || coeffs_f[i][n])
      {
         negate_limbs(coeffs_m + j + 1, coeffs_f[i], n);
         mpn_add_1(coeffs_m + j + 1, coeffs_m + j + 1, n, 1L);
         coeffs_m[j] = -n;
         NORM(coeffs_m + j);
      } else
      {
         copy_limbs(coeffs_m + j + 1, coeffs_f[i], n);
         coeffs_m[j] = n;
         NORM(coeffs_m + j);         
      }
   }
   poly_mpn->length = poly_f->length;   

}

static inline long __get_next_coeff(mp_limb_t * coeff_m, long * borrow, long * coeff, long mask)
{ 
   if ((long) coeff_m[0] == 0) *coeff = -*borrow;
   else if ((long) coeff_m[0] > 0) *coeff = coeff_m[1] - *borrow;
   else *coeff = (-coeff_m[1] - *borrow);
   *borrow = 0UL;
   if (*coeff < 0) 
   {
      *borrow = 1UL;
   }
   *coeff&=mask;
   
   return *coeff;
}

/*static inline long __get_next_coeff_unsigned(mp_limb_t * coeff_m, long * borrow, long * coeff, long mask)
{ 
   if ((long) coeff_m[0] == 0) *coeff = 0;
   else *coeff = coeff_m[1];
   *coeff&=mask;
   
   return *coeff;
}*/

void ZmodFpoly_bit_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                            unsigned long bundle, unsigned long bits)
{   
   unsigned long i, k, skip;
   unsigned long n = poly_f->n;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * array, * next_point;
    
   unsigned long temp;
   half_ulong lower;
   long coeff;
   long borrow;
   mp_limb_t extend;
   
   unsigned long coeffs_per_limb = FLINT_BITS_PER_LIMB/bits;

   const unsigned long mask = (1UL<<bits)-1;
      
   poly_f->length = 0;
   i=0;
      
   while (coeff_m < poly_mpn->coeffs + 2*poly_mpn->length)
   {
      k=0; skip=0;
      coeff = 0; borrow = 0L; temp = 0;
      array = poly_f->coeffs[i];
      i++;
   
      next_point = coeff_m + 2*bundle;
      if (next_point > poly_mpn->coeffs + 2*poly_mpn->length) 
         next_point = poly_mpn->coeffs + 2*poly_mpn->length;
         
      while (coeff_m < next_point)
      {
         // k is guaranteed to be less than FLINT_BITS_PER_LIMB at this point
         while ((k<HALF_FLINT_BITS_PER_LIMB)&&(coeff_m < next_point))
         {
            temp+=(__get_next_coeff(coeff_m, &borrow, &coeff, mask) << k);
            coeff_m+=2; k+=bits;
         }
         // k may exceed FLINT_BITS_PER_LIMB at this point but is less than 96

         if (k>FLINT_BITS_PER_LIMB)
         {
            // if k > FLINT_BITS_PER_LIMB write out a whole limb and read in remaining bits of coeff
            array[skip] = temp;
            skip++;
            temp=(coeff>>(bits+FLINT_BITS_PER_LIMB-k));
            k=(k-FLINT_BITS_PER_LIMB);
            // k < HALF_FLINT_BITS_PER_LIMB
         } else
         {
            // k <= FLINT_BITS_PER_LIMB
            if (k >= HALF_FLINT_BITS_PER_LIMB)
            {
               // if k >= HALF_FLINT_BITS_PER_LIMB store bottom HALF_FLINT_BITS_PER_LIMB bits
               lower = (half_ulong)temp;
               k-=HALF_FLINT_BITS_PER_LIMB;
               temp>>=HALF_FLINT_BITS_PER_LIMB;
               // k is now <= HALF_FLINT_BITS_PER_LIMB

               while ((k<HALF_FLINT_BITS_PER_LIMB)&&(coeff_m < next_point))
               {
                  temp+=(__get_next_coeff(coeff_m, &borrow, &coeff, mask) << k);
                  coeff_m+=2; k+=bits;
               }
               // k may again exceed FLINT_BITS_PER_LIMB bits but is less than 96
               if (k>FLINT_BITS_PER_LIMB)
               {
                  // if k > FLINT_BITS_PER_LIMB, write out bottom HALF_FLINT_BITS_PER_LIMB bits (along with HALF_FLINT_BITS_PER_LIMB bits from lower)
                  // read remaining bits from coeff and reduce k by HALF_FLINT_BITS_PER_LIMB
                  array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+(unsigned long)lower;
                  skip++;
                  temp>>=HALF_FLINT_BITS_PER_LIMB;
                  temp+=((coeff>>(bits+FLINT_BITS_PER_LIMB-k))<<HALF_FLINT_BITS_PER_LIMB);
                  k=(k-HALF_FLINT_BITS_PER_LIMB);
                  // k < FLINT_BITS_PER_LIMB and we are ready to read next coefficient if there is one
               } else if (k >= HALF_FLINT_BITS_PER_LIMB) 
               {
                  // k <= FLINT_BITS_PER_LIMB
                  // if k >= HALF_FLINT_BITS_PER_LIMB write out bottom HALF_FLINT_BITS_PER_LIMB bits (along with lower)
                  // and reduce k by HALF_FLINT_BITS_PER_LIMB
                  k-=HALF_FLINT_BITS_PER_LIMB;
                  array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
                  temp>>=HALF_FLINT_BITS_PER_LIMB;
                  skip++;
                  // k is now less than or equal to HALF_FLINT_BITS_PER_LIMB and we are now ready to read 
                  // the next coefficient if there is one
               } else
               {
                  // k < HALF_FLINT_BITS_PER_LIMB
                  // there isn't enough to write out a whole FLINT_BITS_PER_LIMB bits, so put it all 
                  // together in temp
                  temp = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
                  k+=HALF_FLINT_BITS_PER_LIMB;
                  // k is now guaranteed to be less than FLINT_BITS_PER_LIMB and we are ready for the
                  // next coefficient if there is one
               }
            } // if
         } // else
         poly_f->length++;
      } // while

      // sign extend the last FLINT_BITS_PER_LIMB bits we write out
      if (skip < n)
      {
        if (borrow) temp+=(-1UL << k);
        array[skip] = temp;
        skip++;
      } 
      // sign extend the remainder of the array, reducing modulo p 
      extend = 0;
      if (borrow) extend = -1L;
      while (skip < n+1) 
      {
         array[skip] = extend;
         skip++;
      }
   } // while
#if DEBUG
   for (unsigned long i = 0; i < n+1; i++)
   {
       printf("%lx ",poly_f->coeffs[0][i]);
   }
   printf("\n");
#endif

}


void ZmodFpoly_bit_unpack_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, 
                              unsigned long bundle, unsigned long bits)
{
   unsigned long k, skip;

   unsigned long temp2;
   unsigned long temp;
   unsigned long full_limb;
   unsigned long carry;
    
   mp_limb_t* array;
    
   const unsigned long mask = (1UL<<bits)-1;
   const unsigned long sign_mask = (1UL<<(bits-1));

   unsigned long s;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * next_point;
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long n = poly_f->n;
   
#if DEBUG
   for (unsigned long i = 0; i < n+1; i++)
   {
       printf("%lx ",poly_f->coeffs[0][i]);
   }
   printf("\n");
#endif
    
   for (unsigned long i = 0; coeff_m < poly_mpn->coeffs + poly_mpn->length*size_m; i++)
   {
      array = poly_f->coeffs[i];
      ZmodF_normalise(array, n);

      k=0; skip=0; carry = 0UL; temp2 = 0;
      next_point = coeff_m + size_m*bundle;
      if (next_point > poly_mpn->coeffs + poly_mpn->length*size_m) next_point = poly_mpn->coeffs + poly_mpn->length*size_m;
      
      while (coeff_m < next_point)
      {
         // read in a full limb
         full_limb = array[skip];
         temp2 += l_shift(full_limb,k);
         s=FLINT_BITS_PER_LIMB-k;
         k+=s;
         while ((k >= bits)&&(coeff_m < next_point))
         {
            if (!(temp2&sign_mask)) 
            {
               __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask)+carry);
               carry = 0UL;
            }  
            else
            {
               temp = ((-temp2)&mask)-carry;
               __Zpoly_mpn_sub_coeff_ui(coeff_m, temp);
               carry = 1UL;
            }
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         // read in remainder of full_limb
         temp2 += l_shift(r_shift(full_limb,s),k);
         k+=(FLINT_BITS_PER_LIMB-s);
       
         while ((k >= bits)&&(coeff_m < next_point))
         {
            if (!(temp2&sign_mask)) 
            {
               __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask)+carry);
               carry = 0UL;
            }
            else
            {
               temp = ((-temp2)&mask)-carry;
               __Zpoly_mpn_sub_coeff_ui(coeff_m, temp);
               carry = 1UL;
            }
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         skip++;
      }
   }
}
 
void ZmodFpoly_bit_unpack_unsigned_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, 
                              unsigned long bundle, unsigned long bits)
{
   unsigned long k, l, skip;

   unsigned long temp2;
   unsigned long temp;
   unsigned long full_limb;
    
   mp_limb_t* array;
    
   const unsigned long mask = (1UL<<bits)-1;

   unsigned long s;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * next_point;
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long n = poly_f->n;
       
   for (unsigned long i = 0; coeff_m < poly_mpn->coeffs + poly_mpn->length*size_m; i++)
   {
      array = poly_f->coeffs[i];
      ZmodF_normalise(array, n);

      k=0; skip=0; temp2 = 0;
      next_point = coeff_m + size_m*bundle;
      if (next_point > poly_mpn->coeffs + poly_mpn->length*size_m) next_point = poly_mpn->coeffs + poly_mpn->length*size_m;
      
      while (coeff_m < next_point)
      {
         // read in a full limb
         full_limb = array[skip];
         temp2 += l_shift(full_limb,k);
         s=FLINT_BITS_PER_LIMB-k;
         k+=s;
         while ((k >= bits)&&(coeff_m < next_point))
         {
            __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask));
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         // read in remainder of full_limb
         temp2 += l_shift(r_shift(full_limb,s),k);
         k+=(FLINT_BITS_PER_LIMB-s);
       
         while ((k >= bits)&&(coeff_m < next_point))
         {
            __Zpoly_mpn_add_coeff_ui(coeff_m, temp2&mask);
            coeff_m += size_m;
            temp2>>=bits;
            l++;
            k-=bits;
         }
         // k is now less than bits
         skip++;
      }
   }
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
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->n == y->n);

   for (unsigned long i = 0; i < y->length; i++)
      ZmodF_set(x->coeffs[i], y->coeffs[i], x->n);
   x->length = y->length;
}


void ZmodFpoly_pointwise_mul(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   // todo: use FLINT memory allocation
   
   mp_limb_t* scratch = (mp_limb_t*) malloc((2 * x->n) * sizeof(mp_limb_t));

   if (x != y)
      for (unsigned long i = 0; i < x->length; i++)
         ZmodF_mul(res->coeffs[i], x->coeffs[i], y->coeffs[i], scratch, x->n);
   else
      for (unsigned long i = 0; i < x->length; i++)
         ZmodF_sqr(res->coeffs[i], x->coeffs[i], scratch, x->n);

   res->length = x->length;
   
   free(scratch);
}


void ZmodFpoly_add(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   for (unsigned long i = 0; i < x->length; i++)
      ZmodF_add(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodFpoly_sub(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   for (unsigned long i = 0; i < x->length; i++)
      ZmodF_sub(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodFpoly_normalise(ZmodFpoly_t poly)
{
   for (unsigned long i = 0; i < poly->length; i++)
      ZmodF_normalise(poly->coeffs[i], poly->n);
}


void ZmodFpoly_rescale(ZmodFpoly_t poly)
{
   if (poly->depth == 0)
      return;

   for (unsigned long i = 0; i < poly->length; i++)
      ZmodF_short_div_2exp(poly->coeffs[i], poly->coeffs[i],
                           poly->depth, poly->n);
}


/****************************************************************************

   Forward fourier transforms (internal code)

****************************************************************************/


void _ZmodFpoly_FFT_iterative(
            ZmodF_t* x, unsigned long depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   FLINT_ASSERT(depth >= 1);

   unsigned long i, s, start;
   ZmodF_t* y, * z;

   // root is the (2^depth)-th root of unity for the current layer,
   // measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);

   // half = half the current block length
   unsigned long half = 1UL << (depth - 1);
   unsigned long half_skip = half * skip;

   unsigned long layer;

   // =========================================================================
   // Special case for first layer, if root and/or twist involve sqrt2
   // rotations.

   if ((root | twist) & 1)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length <= half)
      {
         // Only need first output for each butterfly,
         // i.e. (a, b) -> (a + b, ?)
         if (nonzero > half)
            for (i = 0, y = x; i < nonzero - half; i++, y += skip)
               ZmodF_add(y[0], y[0], y[half_skip], n);
      }
      else
      {
         // Need both outputs for each butterfly.
         if (nonzero <= half)
         {
            // The second half of each butterfly input are zeroes, so we just
            // computing (a, 0) -> (a, ra), where r is the appropriate root
            // of unity.
            for (i = 0, s = twist, y = x; i < nonzero;
                 i++, s += root, y += skip)
            {
               ZmodF_mul_sqrt2exp(y[half_skip], y[0], s, n);
            }
         }
         else
         {
            // If nonzero > half, then we need some full butterflies...
            for (i = 0, s = twist, y = x; i < nonzero - half;
                 i++, s += root, y += skip)
            {
               ZmodF_forward_butterfly_sqrt2exp(y, y + half_skip, 
                                                scratch, s, n);
            }
            // and also some partial butterflies (a, 0) -> (a, ra).
            for (; i < half; i++, s += root, y += skip)
               ZmodF_mul_sqrt2exp(y[half_skip], y[0], s, n);
         }
      }

      // Here we switch to measuring roots as powers of 2, but we also need
      // to update to the next layer's roots, and these two actions cancel
      // each other out :-)
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
      if (nonzero > 2*half)
         nonzero = 2*half;
   
      layer = 1;
   }
   else
   {
      // no special case for first layer
      layer = 0;

      // switch to measuring roots as powers of 2
      root >>= 1;
      twist >>= 1;
   }

   // =========================================================================
   // This section handles the layers where there are still zero coefficients
   // to take advantage of. In most cases, this will only happen for the
   // first layer or two, so we don't bother with specialised limbshift-only
   // code for these layers.

   // Note: from here on there are no sqrt2 rotations, and we measure all
   // roots as powers of 2.
   
   for (; (layer < depth) && (nonzero < 2*half); layer++)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length_remainder > half)
      {
         // If length overhangs by more than half the block, then we need to
         // perform full butterflies on the last block (i.e. the last block
         // doesn't get any special treatment).
         length_quantised += 2*half;
      }
      else if (length_remainder)
      {
         // If length overhangs the block by at most half the block size,
         // then we only need to compute the first output of each butterfly
         // for this block, i.e. (a, b) -> (a + b)
         if (nonzero > half)
         {
            y = x + skip * length_quantised;
            for (i = 0; i < nonzero - half; i++, y += skip)
               ZmodF_add(y[0], y[0], y[half_skip], n);
         }
      }

      if (nonzero <= half)
      {
         // If nonzero <= half, then the second half of each butterfly input
         // are zeroes, so we just computing (a, 0) -> (a, ra), where r is the
         // appropriate root of unity.
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, s = twist, z = y; i < nonzero;
                 i++, s += root, z += skip)
            {
               ZmodF_mul_2exp(z[half_skip], z[0], s, n);
            }
         }
      }
      else
      {
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            // If nonzero > half, then we need some full butterflies...
            for (i = 0, s = twist, z = y; i < nonzero - half;
                 i++, s += root, z += skip)
            {
               ZmodF_forward_butterfly_2exp(z, z + half_skip, scratch, s, n);
            }
            // and also some partial butterflies (a, 0) -> (a, ra).
            for (; i < half; i++, s += root, z += skip)
               ZmodF_mul_2exp(z[half_skip], z[0], s, n);
         }
      }

      // Update roots of unity
      twist <<= 1;
      root <<= 1;
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
      
      if (nonzero > 2*half)
         // no more zero coefficients to take advantage of:
         nonzero = 2*half;
   }

   // =========================================================================
   // Now we may assume there are no more zero coefficients.

   for (; layer < depth; layer++)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length_remainder > half)
      {
         // If length overhangs by more than half the block, then we need to
         // perform full butterflies on the last block (i.e. the last block
         // doesn't get any special treatment).
         length_quantised += 2*half;
      }
      else if (length_remainder)
      {
         // If length overhangs the block by at most half the block size,
         // then we only need to compute the first output of each butterfly
         // for this block, i.e. (a, b) -> (a + b)
         y = x + skip * length_quantised;
         for (i = 0; i < half; i++, y += skip)
            ZmodF_add(y[0], y[0], y[half_skip], n);
      }
      
      // To keep the inner loops long, we have two versions of the next loop.
      if (layer < depth/2)
      {
         // Version 1: only a few relatively long blocks.
         
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, s = twist, z = y; i < half; i++, s += root, z += skip)
               ZmodF_forward_butterfly_2exp(z, z + half_skip, scratch, s, n);
         }
      }
      else
      {
         // Version 2: lots of short blocks.
         
         // Two sub-versions, depending on whether the rotations are all by
         // a whole number of limbs.
         if ((root | twist) & (FLINT_BITS_PER_LIMB - 1))
         {
            // Version 2a: rotations still involve bitshifts.
            for (i = 0, s = twist, y = x; i < half; i++, s += root, y += skip)
               for (start = 0, z = y; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_forward_butterfly_2exp(z, z + half_skip,
                                               scratch, s, n);
               }
         }
         else
         {
            // Version 2b: rotations involve only limbshifts.
            unsigned long root_limbs = root / FLINT_BITS_PER_LIMB;

            if (twist == 0)
            {
               // special case, since ZmodF_forward_butterfly_Bexp doesn't
               // allow zero rotation count
               for (start = 0, z = x; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_simple_butterfly(z, z + half_skip, scratch, n);
               }
               i = 1;
               y = x + skip;
               s = root_limbs;
            }
            else
            {
               i = 0;
               y = x;
               s = twist / FLINT_BITS_PER_LIMB;
            }
            
            for (; i < half; i++, s += root_limbs, y += skip)
               for (start = 0, z = y; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_forward_butterfly_Bexp(z, z + half_skip,
                                               scratch, s, n);
               }
         }
      }

      // Update roots of unity
      twist <<= 1;
      root <<= 1;
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
   }
}


/*
Factors FFT of length 2^depth into length 2^rows_depth and length 2^cols_depth
transforms
*/
void _ZmodFpoly_FFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(rows_depth >= 1);
   FLINT_ASSERT(cols_depth >= 1);
   
   unsigned long depth = rows_depth + cols_depth;
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);

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
                     twist << rows_depth, n, scratch);

   if (length_cols)
      // The relevant portion of the last row:
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, length_cols,
                     twist << rows_depth, n, scratch);
}



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

   if (depth == 0)
      return;

   // If the data fits in L1 (2^depth coefficients of length n+1, plus a
   // scratch buffer), then use the iterative transform. Otherwise factor the
   // FFT into two chunks.
   if (depth == 1 ||
       ((1 << depth) + 1) * (n+1) <= ZMODFPOLY_FFT_FACTOR_THRESHOLD)
   {
      _ZmodFpoly_FFT_iterative(x, depth, skip, nonzero, length,
                               twist, n, scratch);
   }
   else
   {
      unsigned long rows_depth = depth >> 1;
      unsigned long cols_depth = depth - rows_depth;
      _ZmodFpoly_FFT_factor(x, rows_depth, cols_depth, skip, nonzero, length,
                            twist, n, scratch);
   }
}



/****************************************************************************

   Inverse fourier transforms (internal code)

****************************************************************************/


// temporary during development
void _ZmodFpoly_IFFT_basecase(ZmodF_t* x, unsigned long depth, unsigned long skip,
                     unsigned long nonzero, unsigned long length, int extra,
                     unsigned long twist, unsigned long n,
                     ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));


   FLINT_ASSERT(depth <= 1);
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);
   
   if (depth == 0)
      return;
      
   if (length == 0)
   {
      if (nonzero == 2)
         ZmodF_add(x[0], x[0], x[skip], n);
      ZmodF_short_div_2exp(x[0], x[0], 1, n);
   }
   else if (length == 1)
   {
      if (nonzero == 1)
      {
         if (extra)
            ZmodF_mul_sqrt2exp(x[skip], x[0], twist, n);
         ZmodF_add(x[0], x[0], x[0], n);
      }
      else  // nonzero == 2
      {
         if (extra)
         {
            ZmodF_sub(scratch[0], x[0], x[skip], n);
            ZmodF_add(x[0], x[0], scratch[0], n);
            ZmodF_mul_sqrt2exp(x[skip], scratch[0], twist, n);
         }
         else
         {
            ZmodF_add(x[0], x[0], x[0], n);
            ZmodF_sub(x[0], x[0], x[skip], n);
         }
      }
   }
   else  // length == 2
      ZmodF_inverse_butterfly_sqrt2exp(x, x + skip, scratch, twist, n);
}



/*
This one is for when there is no truncation.
*/
void _ZmodFpoly_IFFT_iterative(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   
   // for now just call the main routine back again, or the basecase :-)
   if (depth <= 1)
      _ZmodFpoly_IFFT_basecase(x, depth, skip, 1 << depth, 1 << depth, 0, twist, n, scratch);
   else
      _ZmodFpoly_IFFT(x, depth, skip, 1 << depth, 1 << depth, 0, twist, n, scratch);
}



/*
This one's for working in L1 when truncation is involved. It splits into
two halves.
*/
void _ZmodFpoly_IFFT_recursive(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long nonzero, unsigned long length, int extra,
               unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));

   if (depth == 0)
      return;

   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);

   long size = 1UL << depth;

   if (length == size)
   {
      // no truncation necessary
      _ZmodFpoly_IFFT_iterative(x, depth, skip, twist, n, scratch);
      return;
   }

   long cols = size >> 1;

   long i, j;
   ZmodF_t* y;

   // symbols in the following diagrams:
   // A = fully untransformed coefficient
   // a = fully untransformed coefficient (implied zero)
   // B = intermediate coefficient
   // b = intermediate coefficient (implied zero)
   // C = fully transformed coefficient
   // c = fully transformed coefficient (implied zero)
   // ? = garbage that we don't care about
   // * = the extra C coefficient, or "?" if no extra coefficient requested
   
   // the horizontal transforms convert between B and C
   // the vertical butterflies convert between A and B

   if ((length < cols) || (length == cols && !extra))
   {
      // The input could look like one of the following:
      // CCCCAAAA      CCCCAAAA      CCCCAAaa      CCCCaaaa
      // AAAAAAaa  or  AAaaaaaa  or  aaaaaaaa  or  aaaaaaaa

      long last_zero_forward_butterfly;
      long last_zero_cross_butterfly;
      long last_cross_butterfly;

      if (nonzero <= cols)
      {
         i = nonzero - 1;
         last_zero_forward_butterfly = length;
         last_zero_cross_butterfly = 0;
      }
      else
      {
         i = cols - 1;
         if (nonzero > length + cols)
         {
            last_zero_forward_butterfly = nonzero - cols;
            last_zero_cross_butterfly = length;
         }
         else
         {
            last_zero_forward_butterfly = length;
            last_zero_cross_butterfly = nonzero - cols;
         }
      }
      
      j = twist + root*i;
      y = x + skip*i;

      // First some forward butterflies ("Aa" => "B?") to make them look like:
      // CCCCAABB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAAAA??  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= last_zero_forward_butterfly; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 1, 0, 1, j, n, scratch);

      // Then some forward butterflies ("AA" => "B?") to make them look like:
      // CCCCBBBB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= (long)length; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 2, 0, 1, j, n, scratch);

      // Transform the first row to make them look like:
      // BBBB*???      BBBB*???      BBBB*???      BBBB*???
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      _ZmodFpoly_IFFT_recursive(x, depth - 1, skip,
                                (nonzero < cols) ? nonzero : cols,
                                length, extra, twist << 1, n, scratch);
      
      // Cross butterflies ("Ba" => "A?") to make them look like:
      // BBBB*???      BBAA*???      AAAA*???      AAAA*???
      // AAAA????  or  AA??????  or  ??????aa  or  ????aaaa
      for (; i >= last_zero_cross_butterfly; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 1, 1, 0, j, n, scratch);
         
      // Cross butterflies ("BA" => "A?") to make them look like:
      // AAAA*???      AAAA*???      AAAA*???      AAAA*???
      // ????????  or  ????????  or  ??????aa  or  ????aaaa
      for (; i >= 0; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 2, 1, 0, j, n, scratch);
   }
   else
   {
      // The input looks like one of these:
      // CCCCCCCC                   CCCCCCCC
      // AAAAaaaa (extra == 1)  or  CCCAAAaa
   
      // Transform first row (no truncation necessary) to make them look like:
      // BBBBBBBB                   BBBBBBBB
      // AAAAaaaa (extra == 1)  or  CCCAAAaa
      _ZmodFpoly_IFFT_iterative(x, depth - 1, skip, twist << 1, n, scratch);

      i = cols-1;
      j = twist + root*i;
      y = x + skip*i;
      
      long last_zero_cross_butterfly = nonzero - cols;
      long last_cross_butterfly = length - cols;
   
      // Cross butterflies ("Ba" => "AB") to make them look like:
      // BBBBAAAA                   BBBBBBAA
      // AAAABBBB (extra == 1)  or  CCCAAABB
      for (; i >= last_zero_cross_butterfly; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 1, 1, 1, j, n, scratch);
         
      // Cross butterflies ("BA" => "AB") to make them look like:
      // AAAAAAAA                   BBBAAAAA
      // BBBBBBBB (extra == 1)  or  CCCBBBBB
      for (; i >= last_cross_butterfly; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 2, 1, 1, j, n, scratch);
      
      // Transform second row to make them look like:
      // AAAAAAAA                   BBBAAAAA
      // *??????? (extra == 1)  or  BBB*????
      _ZmodFpoly_IFFT_recursive(x + skip*cols, depth - 1, skip, cols, length - cols, extra, twist << 1, n, scratch);

      // Inverse butterflies ("BB" => "AA") to make them look like:
      // AAAAAAAA                   AAAAAAAA
      // *??????? (extra == 1)  or  AAA*????
      for (; i >= 0; i--, j -= root, y -= skip)
         _ZmodFpoly_IFFT_basecase(y, 1, skip << (depth - 1), 2, 2, 0, j, n, scratch);
   }
}



/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_IFFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            int extra, unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(rows_depth >= 1);
   FLINT_ASSERT(cols_depth >= 1);

   unsigned long depth = rows_depth + cols_depth;
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);
   
   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // row transforms for the rows where we have all fourier coefficients
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodFpoly_IFFT(y, cols_depth, skip, cols, cols, 0,
                      twist << rows_depth, n, scratch);

   // column transforms where we have enough information
   for (i = length_cols, y = x + (skip * length_cols),
        j = twist + (root*length_cols);
        i < nonzero_cols; i++, y += skip, j += root)
   {
      _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                      length_rows, length_cols ? 1 : extra, j, n, scratch);
   }
   if (nonzero_rows)
      for (; i < cols; i++, y += skip, j += root)
         _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                         length_rows, length_cols ? 1 : extra, j, n, scratch);

   if (length_cols)
   {
      // a single switcheroo row transform
      _ZmodFpoly_IFFT(x + length_rows * (skip << cols_depth), cols_depth,
                      skip, (nonzero_rows ? cols : nonzero_cols),
                      length_cols, extra, twist << rows_depth, n, scratch);

      // remaining column transforms
      for (i = 0, y = x, j = twist; i < length_cols && i < nonzero_cols;
           i++, y += skip, j += root)
      {
         _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                         length_rows + 1, 0, j, n, scratch);
      }
      if (nonzero_rows)
      {
         for (; i < length_cols; i++, y += skip, j += root)
            _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                            length_rows + 1, 0, j, n, scratch);
      }
   }
   else if (extra)
   {
      // need one extra trivial fourier coefficient
      x += length_rows * (skip << cols_depth);
      for (i = 1, y = x + skip; i < (nonzero_rows ? cols : nonzero_cols);
           i++, y += skip)
      {
         ZmodF_add(x[0], x[0], y[0], n);
      }
      ZmodF_short_div_2exp(x[0], x[0], cols_depth, n);
   }
}


/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_IFFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                     unsigned long nonzero, unsigned long length, int extra,
                     unsigned long twist, unsigned long n,
                     ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));

   if (depth == 0)
      return;

/*
   // If the data fits in L1 (2^depth coefficients of length n+1, plus a
   // scratch buffer), then use the iterative transform. Otherwise factor the
   // FFT into two chunks.
   if (depth == 1 ||
       ((1 << depth) + 1) * (n+1) <= ZMODFPOLY_FFT_FACTOR_THRESHOLD)
*/
   if (depth <= 1)
   {
      _ZmodFpoly_IFFT_basecase(x, depth, skip, nonzero, length, extra,
                               twist, n, scratch);
   }
   else
   {
      unsigned long rows_depth = depth >> 1;
      unsigned long cols_depth = depth - rows_depth;
      _ZmodFpoly_IFFT_factor(x, rows_depth, cols_depth, skip, nonzero, length,
                             extra, twist, n, scratch);
   }
}


/****************************************************************************

   Fourier Transform Routines

****************************************************************************/


void ZmodFpoly_FFT(ZmodFpoly_t poly, unsigned long length)
{
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         for (unsigned long i = 0; i < length; i++)
            ZmodF_zero(poly->coeffs[i], poly->n);
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
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (poly->length != 0)
      _ZmodFpoly_IFFT(poly->coeffs, poly->depth, 1, poly->length,
                      poly->length, 0, 0, poly->n, poly->scratch);
}


void ZmodFpoly_convolution(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodFpoly_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodFpoly_FFT(y, length);
      
   ZmodFpoly_pointwise_mul(res, x, y);
   ZmodFpoly_IFFT(res);
   ZmodFpoly_rescale(res);
}


/****************************************************************************

   Negacyclic Fourier Transform Routines
   
****************************************************************************/


void ZmodFpoly_negacyclic_FFT(ZmodFpoly_t poly, unsigned long length)
{
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   // twist on the way in to make it negacyclic
   // todo: this needs to be cleaned up
   unsigned long twist = (2 * poly->n * FLINT_BITS_PER_LIMB) >> poly->depth;
   for (unsigned long i = 1; i < length; i++)
   {
      ZmodF_mul_sqrt2exp(*poly->scratch, poly->coeffs[i], i*twist, poly->n);
      ZmodF_swap(poly->scratch, poly->coeffs + i);
   }

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         for (unsigned long i = 0; i < length; i++)
            ZmodF_zero(poly->coeffs[i], poly->n);
      }
      else
      {
         _ZmodFpoly_FFT(poly->coeffs, poly->depth, 1, poly->length, length,
                        0, poly->n, poly->scratch);
      }
   }

   poly->length = length;
}


void ZmodFpoly_negacyclic_IFFT(ZmodFpoly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (poly->length != 0)
   {
      _ZmodFpoly_IFFT(poly->coeffs, poly->depth, 1, poly->length, poly->length,
                      0, 0, poly->n, poly->scratch);
   }

   // twist on the way out to make it negacyclic
   // todo: this needs to be cleaned up
   unsigned long twist = (2 * poly->n * FLINT_BITS_PER_LIMB) >> poly->depth;
   for (unsigned long i = 1; i < poly->length; i++)
   {
      ZmodF_mul_sqrt2exp(*poly->scratch, poly->coeffs[i],
                         2 * poly->n * FLINT_BITS_PER_LIMB - i*twist, poly->n);
      ZmodF_neg(poly->coeffs[i], *poly->scratch, poly->n);
   }
}


void ZmodFpoly_negacyclic_convolution(ZmodFpoly_t res,
                                      ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodFpoly_negacyclic_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodFpoly_negacyclic_FFT(y, length);
      
   ZmodFpoly_pointwise_mul(res, x, y);
   ZmodFpoly_negacyclic_IFFT(res);
   ZmodFpoly_rescale(res);
}


// end of file ****************************************************************
