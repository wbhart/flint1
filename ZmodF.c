/******************************************************************************

 ZmodF.c

 Copyright (C) 2007, David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
 These are currently used only in the ZmodFpoly module, which supplies the
 Schoenhage-Strassen FFT code.
 
******************************************************************************/

#include "ZmodF.h"


/*
For odd s, finds "limbs" and "bits" such that 2^(s/2) is decomposed into
  2^(-bits) * B^limbs * (1 - B^(n/2)),
where 0 <= bits < FLINT_BITS_PER_LIMB, and 0 <= limbs < 2n.

i.e. we are decomposing a rotation involving a sqrt2 into a fractional
limbshift and a pseudosqrt2 call.

PRECONDIITONS:
   s must be odd
   0 <= s < 2n*FLINT_BITS_PER_LIMB
*/
void ZmodF_decompose_rotation(unsigned long* limbs, unsigned long* bits,
                              unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s & 1);
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);

   // first split into 2^r * (1 - B^(n/2))
   unsigned long r = (s >> 1) - 3*n*FLINT_BITS_PER_LIMB/4;
   if ((long)r < 0)
      r += 2*n*FLINT_BITS_PER_LIMB;

   // now split 2^r into 2^(-bits) and B^limbs
   unsigned long z = r & (FLINT_BITS_PER_LIMB - 1);
   r /= FLINT_BITS_PER_LIMB;
   if (z)
   {
      *bits = FLINT_BITS_PER_LIMB - z;
      if (++r == 2*n)
         r = 0;
   }
   else
      *bits = 0;

   *limbs = r;
}


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
      if (++s == n)
      {
         // special case if s == n-1
         ZmodF_neg(b, a, n);
         ZmodF_short_div_2exp(b, b, FLINT_BITS_PER_LIMB - bits, n);
         return;
      }

      // Need to shift left by s limbs and right by
      // (FLINT_BITS_PER_LIMB - bits) bits.
      bits = FLINT_BITS_PER_LIMB - bits;
      
      // Shift top part of input directly into bottom part of output
      ZmodF_fast_reduce(a, n);
      mp_limb_t carry1 = mpn_rshift(b, a+n-s, s+1, bits);
      mp_limb_t overlap = b[s];
      // complement the part we just shifted in
      long i = s-1;
      do b[i] = ~b[i]; while (--i >= 0);
      
      // shift bottom part of input directly into top part of output
      mp_limb_t carry2 = mpn_rshift(b+s, a, n-s, bits);
      b[n] = -1;  // compensate mod p for 1's complement
      
      // fiddle with carries
      mpn_add_1(b+n-1, b+n-1, 2, carry1);
      mpn_add_1(b+s-1, b+s-1, n-s+2, carry2);
      mpn_sub_1(b+s, b+s, n-s+1, overlap+1);
   }
   else
   {
      if (s)
         ZmodF_mul_Bexp(b, a, s, n);
      else
         ZmodF_set(b, a, n);
   }
}


void ZmodF_mul_sqrt2exp(ZmodF_t b, ZmodF_t a,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(a != b);

   if (s & 1)
   {
      unsigned long limbs, bits;
      ZmodF_decompose_rotation(&limbs, &bits, s, n);
      
      if (n & 1)
         ZmodF_mul_pseudosqrt2_n_odd(b, a, limbs, n);
      else
         ZmodF_mul_pseudosqrt2_n_even(b, a, limbs, n);
         
      if (bits)
         ZmodF_short_div_2exp(b, b, bits, n);
   }
   else
      ZmodF_mul_2exp(b, a, s >> 1, n);
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


void ZmodF_mul_pseudosqrt2_n_odd(ZmodF_t b, ZmodF_t a,
                                 unsigned long s, unsigned long n)
{
   // todo: this is not optimised yet

   FLINT_ASSERT(a != b);
   FLINT_ASSERT((n & 1) == 1);
   FLINT_ASSERT(s < 2*n);
   
   // first multiply by -B^(s+(n+1)/2) from a into b
   unsigned long limbs = s + (n+1)/2 + n;
   while (limbs >= 2*n)
      limbs -= 2*n;

   if (limbs == 0)
      ZmodF_set(b, a, n);
   else if (limbs < n)
      ZmodF_mul_Bexp(b, a, limbs, n);
   else if (limbs == n)
      ZmodF_neg(b, a, n);
   else
   {
      ZmodF_mul_Bexp(b, a, limbs - n, n);
      ZmodF_neg(b, b, n);
   }
   
   // Divide by B^(1/2)
   ZmodF_short_div_2exp(b, b, FLINT_BITS_PER_LIMB / 2, n);
   
   // Now add in B^s a
   if (s == 0)
      ZmodF_add(b, b, a, n);
   else if (s < n)
      ZmodF_div_Bexp_sub(b, b, a, n - s, n);
   else if (s == n)
      ZmodF_sub(b, b, a, n);
   else
      ZmodF_div_Bexp_add(b, b, a, 2*n - s, n);
}


void ZmodF_mul_pseudosqrt2_n_even(ZmodF_t b, ZmodF_t a,
                                  unsigned long s, unsigned long n)
{
   FLINT_ASSERT(a != b);
   FLINT_ASSERT((n & 1) == 0);
   FLINT_ASSERT(s < 2*n);
   
   mp_limb_t carry;

   if (s < n)
   {
      if (s < n/2)
      {
         if (s)
         {
            // We're computing B^s * (1 - B^(n/2)) * a.
            // If input is
            // 0                   n/2-s     n/2                   n-s       n
            // |         x0          |   y0   |         x1          |   y1   |
            // then output should be
            // 0        s                    n/2     n/2+s                   n
            // |  y0-y1 |        x0+x1        | y0+y1  |       -x0+x1        |

            // Store x1 - x0
            b[n] = -mpn_sub_n(b+n/2+s, a+n/2, a, n/2-s);
            // Store x0 + x1 and y0 + y1
            carry = mpn_add_n(b+s, a, a+n/2, n/2);
            signed_add_1(b+s+n/2, n/2-s+1, carry + a[n]);
            // Store y0 - y1
            carry = mpn_sub_n(b, a+n/2-s, a+n-s, s);
            signed_add_1(b+s, n-s+1, -a[n] - carry);
         }
         else
         {
            // Special case for s = 0.

            // Store x1 - x0
            b[n] = a[n] - mpn_sub_n(b+n/2, a+n/2, a, n/2);
            // Store x0 + x1
            carry = mpn_add_n(b, a, a+n/2, n/2);
            signed_add_1(b+n/2, n/2+1, carry + a[n]);
         }
      }
      else
      {
         s -= n/2;
      
         if (s)
         {
            // We're computing B^s * (1 + B^(n/2)) * a.
            // If input is
            // 0                   n/2-s     n/2                   n-s       n
            // |         x0          |   y0   |         x1          |   y1   |
            // then output should be
            // 0        s                    n/2     n/2+s                   n
            // | -y0-y1 |        x0-x1        | y0-y1  |       x0+x1         |
            
            // Store x0 + x1
            // (the -1 compensates mod p for the bottom bit of the complement)
            b[n] = mpn_add_n(b+s+n/2, a, a+n/2, n-s-n/2) - 1;
            // Store x0 - x1 and y0 - y1
            carry = mpn_sub_n(b+s, a, a+n/2, n/2);
            signed_add_1(b+s+n/2, n/2-s+1, -a[n] - carry);
            // Store -y0 - y1
            carry = mpn_add_n(b, a+n/2-s, a+n-s, s);
            // (the -1 compensates for the top bit of the complement)
            signed_add_1(b+s, n-s+1, -a[n] - carry - 1);
            long i = s-1;
            do b[i] = ~b[i]; while (--i >= 0);
         }
         else
         {
            // Special case for s = 0.

            // Store x0 + x1
            b[n] = a[n] + mpn_add_n(b+n/2, a, a+n/2, n/2);
            // Store x0 - x1
            carry = mpn_sub_n(b, a, a+n/2, n/2);
            signed_add_1(b+n/2, n/2+1, -a[n] - carry);
         }
      }
   }
   else
   {
      s -= n;

      if (s < n/2)
      {
         if (s)
         {
            // We're computing B^s * (-1 + B^(n/2)) * a.
            // If input is
            // 0                   n/2-s     n/2                   n-s       n
            // |         x0          |   y0   |         x1          |   y1   |
            // then output should be
            // 0        s                    n/2     n/2+s                   n
            // | -y0+y1 |       -x0-x1        | -y0-y1 |        x0-x1        |
            
            // Store x0 - x1
            b[n] = -mpn_sub_n(b+n/2+s, a, a+n/2, n/2-s);
            // Store -x0 - x1 and -y0 - y1
            carry = mpn_add_n(b+s, a, a+n/2, n/2);
            // (the -1 compensates for the top bit of the complement)
            signed_add_1(b+n/2+s, n/2-s+1, -a[n] - carry - 1);
            long i = n/2-1;
            do b[s+i] = ~b[s+i]; while (--i >= 0);
            // Store y1 - y0
            carry = mpn_sub_n(b, a+n-s, a+n/2-s, s);
            // (the +1 compensates for the bottom bit of the complement)
            signed_add_1(b+s, n-s+1, a[n] - carry + 1);
         }
         else
         {
            // Special case for s = 0.

            // Store x0 - x1
            // (the -1 compensates mod p for the bottom bit of the complement)
            b[n] = -mpn_sub_n(b+n/2, a, a+n/2, n/2) - a[n] - 1;
            // Store -x0 - x1
            carry = mpn_add_n(b, a, a+n/2, n/2);
            // (the -1 compensates for the top bit of the complement)
            signed_add_1(b+n/2, n/2+1, -a[n] - carry - 1);
            long i = n/2-1;
            do b[i] = ~b[i]; while (--i >= 0);
         }
      }
      else
      {
         s -= n/2;
      
         if (s)
         {
            // We're computing B^s * (-1 - B^(n/2)) * a.
            // If input is
            // 0                   n/2-s     n/2                   n-s       n
            // |         x0          |   y0   |         x1          |   y1   |
            // then output should be
            // 0        s                    n/2     n/2+s                   n
            // |  y0+y1 |       -x0+x1        | -y0+y1 |      -x0-x1         |

            // Store -x0 - x1
            // (the -1 compensates for the top bit of the complement)
            b[n] = -mpn_add_n(b+n/2+s, a, a+n/2, n/2-s) - 1;
            long i = n/2-s-1;
            do b[n/2+s+i] = ~b[n/2+s+i]; while (i-- >= 0);
            // Store x1 - x0 and y1 - y0
            carry = mpn_sub_n(b+s, a+n/2, a, n/2);
            // (the +1 compensates for the bottom bit of the complement)
            signed_add_1(b+n/2+s, n/2-s+1, a[n] - carry + 1);
            // Store y0 + y1
            carry = mpn_add_n(b, a+n/2-s, a+n-s, s);
            signed_add_1(b+s, n-s+1, a[n] + carry);
         }
         else
         {
            // Special case for s = 0.
            
            // Store -x0 - x1
            // (the -1 compensates for the top bit of the complement)
            b[n] = -mpn_add_n(b+n/2, a, a+n/2, n/2) - a[n] - 1;
            long i = n/2-1;
            do b[n/2+i] = ~b[n/2+i]; while (i-- >= 0);
            // Store x1 - x0
            carry = mpn_sub_n(b, a+n/2, a, n/2);
            // (the +1 compensates for the bottom bit of the complement)
            signed_add_1(b+n/2, n/2+1, a[n] - carry + 1);
         }
      }
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
   {
      unsigned long limbs, bits;
      ZmodF_decompose_rotation(&limbs, &bits, s, n);
      
      if (limbs == 0)
         ZmodF_sub(*z, *a, *b, n);
      else if (limbs < n)
         ZmodF_sub_mul_Bexp(*z, *a, *b, limbs, n);
      else if (limbs == n)
         ZmodF_sub(*z, *b, *a, n);
      else
         ZmodF_sub_mul_Bexp(*z, *b, *a, limbs - n, n);
      
      ZmodF_add(*a, *a, *b, n);
      
      if (n & 1)
         ZmodF_mul_pseudosqrt2_n_odd(*b, *z, 0, n);
      else
         ZmodF_mul_pseudosqrt2_n_even(*b, *z, 0, n);
         
      if (bits)
         ZmodF_short_div_2exp(*b, *b, bits, n);
   }
   else
      ZmodF_forward_butterfly_2exp(a, b, z, s >> 1, n);
}


void ZmodF_inverse_butterfly_2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                  unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

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
   {
      unsigned long limbs, bits;
      ZmodF_decompose_rotation(&limbs, &bits, 2*n*FLINT_BITS_PER_LIMB - s, n);
      
      if (n & 1)
         ZmodF_mul_pseudosqrt2_n_odd(*z, *b, 0, n);
      else
         ZmodF_mul_pseudosqrt2_n_even(*z, *b, 0, n);
      
      if (bits)
         ZmodF_short_div_2exp(*z, *z, bits, n);

      if (limbs == 0)
      {
         ZmodF_add(*b, *a, *z, n);
         ZmodF_sub(*a, *a, *z, n);
      }
      else if (limbs < n)
      {
         ZmodF_div_Bexp_sub(*b, *a, *z, n - limbs, n);
         ZmodF_div_Bexp_add(*a, *a, *z, n - limbs, n);
      }
      else if (limbs == n)
      {
         ZmodF_sub(*b, *a, *z, n);
         ZmodF_add(*a, *a, *z, n);
      }
      else
      {
         ZmodF_div_Bexp_add(*b, *a, *z, 2*n - limbs, n);
         ZmodF_div_Bexp_sub(*a, *a, *z, 2*n - limbs, n);
      }
   }
   else
      ZmodF_inverse_butterfly_2exp(a, b, z, s >> 1, n);
}


// end of file ****************************************************************
