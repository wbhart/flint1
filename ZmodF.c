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
         ZmodF_clear(res, n);
      res[0] = 1;
   }
   else
   {
      // do the product
      mpn_mul_n(scratch, x, x, n);
      // reduce mod p
      res[n] = -mpn_sub_n(res, scratch, scratch + n, n);
   }
}


// end of file ****************************************************************
