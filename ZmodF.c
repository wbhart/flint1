/******************************************************************************

 ZmodF.c

 Copyright (C) 2007, David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
 These are currently used only in the ZmodFpoly module, which supplies the
 Schoenhage-Strassen FFT code.
 
******************************************************************************/

#include "ZmodF.h"


void ZmodF_normalise(ZmodF_t a, unsigned long n)
{
   mp_limb_t hi = a[n];

   if ((mp_limb_signed_t) hi < 0)
   {
      // If top limb (hi) is negative, we add -hi multiples of p
      a[n] = 0;
      mpn_add_1(a, a, n + 1, -hi);

      // If the result is >= p (very unlikely)...
      if (a[n] && a[0])
      {
         // ... need to subtract off p.
         a[n] = 0;
         a[0]--;
      }
   }
   else
   {
      // If top limb (hi) is non-negative, we subtract hi multiples of p
      a[n] = 0;
      mpn_sub_1(a, a, n + 1, hi);

      // If the result is negative (very unlikely)...
      if (a[n])
      {
         // ... need to add back p.
         a[n] = 0;
         mpn_add_1(a, a, n + 1, 1);
      }
   }
}


void ZmodF_mul(ZmodF_t res, ZmodF_t a, ZmodF_t b, mp_limb_t* scratch,
               unsigned long n)
{
   ZmodF_normalise(a, n);
   ZmodF_normalise(b, n);
      
   if (a[n])
   {
      // special case when a = -1 mod p
      ZmodF_neg(res, b, n);
   }
   else if (b[n])
   {
      // special case when b = -1 mod p
      ZmodF_neg(res, a, n);
   }
   else
   {
      // do the product
      mpn_mul_n(scratch, a, b, n);
      // reduce mod p
      res[n] = -mpn_sub_n(res, scratch, scratch + n, n);
   }
}


void ZmodF_sqr(ZmodF_t res, ZmodF_t a, mp_limb_t* scratch, unsigned long n)
{
   ZmodF_normalise(a, n);
      
   if (a[n])
   {
      // special case when a = -1 mod p
      if (a == res)
         res[n] = 0;
      else
         ZmodF_zero(res, n);
      res[0] = 1;
   }
   else
   {
      // do the product
      mpn_mul_n(scratch, a, a, n);
      // reduce mod p
      res[n] = -mpn_sub_n(res, scratch, scratch + n, n);
   }
}


void ZmodF_mul_2exp(ZmodF_t b, ZmodF_t a, unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(a != b);

   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
   if (bits)
   {
      // shift left by s+1 limbs...
      if (++s == n)
         ZmodF_neg(b, a, n);
      else
         ZmodF_mul_Bexp(b, a, s, n);

      // ... and then shift right by remaining bits
      ZmodF_short_div_2exp(b, b, FLINT_BITS_PER_LIMB - bits, n);
   }
   else
   {
      if (s)
         ZmodF_mul_Bexp(b, a, s, n);
      else
         ZmodF_set(b, a, n);
   }
}


void ZmodF_mul_sqrt2exp(ZmodF_t* b, ZmodF_t* a, ZmodF_t* z,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

   if (s & 1)
      // not implemented yet
      abort();
      
   ZmodF_mul_2exp(*b, *a, s >> 1, n);
}


void ZmodF_sub_mul_2exp(ZmodF_t c, ZmodF_t a, ZmodF_t b,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(c != a);
   FLINT_ASSERT(c != b);

   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
   if (bits)
   {
      // shift a-b left by s+1 limbs...
      if (++s == n)
         ZmodF_sub(c, b, a, n);
      else
         ZmodF_sub_mul_Bexp(c, a, b, s, n);
      
      // ... and then shift right by remaining bits
      ZmodF_short_div_2exp(c, c, FLINT_BITS_PER_LIMB - bits, n);
   }
   else
   {
      if (s)
         ZmodF_sub_mul_Bexp(c, a, b, s, n);
      else
         ZmodF_sub(c, a, b, n);
   }
}


void ZmodF_forward_butterfly_2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                  unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

   ZmodF_sub_mul_2exp(*z, *a, *b, s, n);
   ZmodF_add(*a, *a, *b, n);
   ZmodF_swap(b, z);
}


void ZmodF_forward_butterfly_sqrt2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                      unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

   if (s & 1)
      // not implemented yet
      abort();
      
   ZmodF_forward_butterfly_2exp(a, b, z, s >> 1, n);
}


void ZmodF_inverse_butterfly_2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                  unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   if (bits)
      // shift right by leftover bits
      ZmodF_short_div_2exp(*b, *b, bits, n);

   s /= FLINT_BITS_PER_LIMB;
   if (s)
   {
      ZmodF_div_Bexp_sub(*z, *a, *b, s, n);
      ZmodF_div_Bexp_add(*a, *a, *b, s, n);
   }
   else
   {
      ZmodF_sub(*z, *a, *b, n);
      ZmodF_add(*a, *a, *b, n);
   }

   ZmodF_swap(z, b);
}


void ZmodF_inverse_butterfly_sqrt2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                      unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

   if (s & 1)
      // not implemented yet
      abort();
      
   ZmodF_inverse_butterfly_2exp(a, b, z, s >> 1, n);
}

// end of file ****************************************************************
