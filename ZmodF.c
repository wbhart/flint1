/******************************************************************************

 ZmodF.c

 Copyright (C) 2007, David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
 These are currently used only in the ZmodFpoly module, which supplies the
 Schoenhage-Strassen FFT code.
 
******************************************************************************/

#include "ZmodF.h"


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
