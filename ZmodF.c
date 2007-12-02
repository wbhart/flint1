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

 ZmodF.c

 Copyright (C) 2007, David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS.
 
 These are currently used only in the ZmodF_poly module, which supplies the
 Schoenhage-Strassen FFT code.
 
******************************************************************************/

#include "ZmodF.h"
#include "longlong_wrapper.h"
#include "longlong.h"


/*
For odd s, finds "limbs" and "bits" such that 2^(s/2) is decomposed into
  2^(-bits) * B^limbs * (1 - B^(n/2)),
where 0 <= bits < FLINT_BITS, and 0 <= limbs < 2n.

i.e. we are decomposing a rotation involving a sqrt2 into a fractional
limbshift and a pseudosqrt2 call.

PRECONDIITONS:
   s must be odd
   0 <= s < 2n*FLINT_BITS
*/
void ZmodF_decompose_rotation(unsigned long* limbs, unsigned long* bits,
                              unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s & 1);
   FLINT_ASSERT(s < 2*n*FLINT_BITS);

   // first split into 2^r * (1 - B^(n/2))
   unsigned long r = (s >> 1) - 3*n*FLINT_BITS/4;
   if ((long)r < 0)
      r += 2*n*FLINT_BITS;

   // now split 2^r into 2^(-bits) and B^limbs
   unsigned long z = r & (FLINT_BITS - 1);
   r /= FLINT_BITS;
   if (z)
   {
      *bits = FLINT_BITS - z;
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


void ZmodF_mul_2exp(ZmodF_t b, ZmodF_t a, unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS);
   FLINT_ASSERT(a != b);

   unsigned long bits = s & (FLINT_BITS - 1);
   s /= FLINT_BITS;
   if (bits)
   {
      if (++s == n)
      {
         // special case if s == n-1
         ZmodF_neg(b, a, n);
         ZmodF_short_div_2exp(b, b, FLINT_BITS - bits, n);
         return;
      }

      // Need to shift left by s limbs and right by
      // (FLINT_BITS - bits) bits.
      bits = FLINT_BITS - bits;
      
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
   FLINT_ASSERT(s < 2*n*FLINT_BITS);
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
   FLINT_ASSERT(s < n*FLINT_BITS);
   FLINT_ASSERT(c != a);
   FLINT_ASSERT(c != b);

   unsigned long bits = s & (FLINT_BITS - 1);
   s /= FLINT_BITS;
   if (bits)
   {
      // shift a-b left by s+1 limbs...
      if (++s == n)
         ZmodF_sub(c, b, a, n);
      else
         ZmodF_sub_mul_Bexp(c, a, b, s, n);
      
      // ... and then shift right by remaining bits
      ZmodF_short_div_2exp(c, c, FLINT_BITS - bits, n);
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
   FLINT_ASSERT(a != b);
   FLINT_ASSERT((n & 1) == 1);
   FLINT_ASSERT(s < 2*n);

   // Let ss = s+(n+1)/2 mod n, in the range (0, n].
   unsigned long ss = s + (n+1)/2;
   if (ss > n)
      ss -= n;
   if (ss > n)
      ss -= n;

   // The next block of code has the following effect.
   // Pretend that the input is normalised to be divisible by B^(1/2)
   // (i.e. imagine that the bottom half-limb has been relocated mod p to the
   // overflow limb). Now write the input as
   //    a = (X + Y*B^(n-ss) + Z*B^n) * B^(1/2),
   // where X is exactly n-ss limbs long, Y is exactly ss limbs long,
   // and where Z is a signed quantity, just a few bits long.
   // This block computes Z, and sets b to Y + X*B^ss.
   // (It doesn't set the overflow limb of b to anything meaningful.)
   mp_limb_signed_t Z;
   ZmodF_fast_reduce(a, n);
   mpn_rshift(b, a+n-ss, ss+1, FLINT_BITS/2);
   mp_limb_t underflow = mpn_rshift(b+ss, a, n-ss+1, FLINT_BITS/2);
   sub_ddmmss(Z, b[ss-1], 0, b[ss-1], 0, underflow);
   
   mp_limb_t carry1, carry2;

   // Now we need to add in B^s*a, taking into account the fact that some of
   // b currently has the wrong sign. We split into various cases depending
   // on relative locations of s and ss, and depending on sign issues.
   if (s <= n)
   {
      if (s <= (n-1)/2)
      {
         carry1 = s ? -mpn_sub_n(b, b, a+n-s, s) : 0;
         carry2 = mpn_add_n(b+s, b+s, a, (n+1)/2);
         b[n] = (ss < n) ? -mpn_sub_n(b+ss, a+(n+1)/2, b+ss, n-ss) : 0;
         signed_add_1(b+s, n-s+1, carry1 - a[n]);
         signed_add_1(b+ss, n-ss+1, Z + carry2);
      }
      else
      {
         carry1 = mpn_add_n(b, b, a+n-s, ss);
         long i = ss-1;
         do b[i] = ~b[i]; while (--i >= 0);
         carry2 = (n > 1) ? mpn_sub_n(b+ss, b+ss, a+(n+1)/2, (n-1)/2) : 0;
         b[n] = (s < n) ? mpn_add_n(b+s, b+s, a, n-s) - 1 : -1;
         signed_add_1(b+ss, n-ss+1, -carry1 - Z - 1);
         signed_add_1(b+s, n-s+1, -carry2 - a[n]);
      }
   }
   else
   {
      s -= n;

      if (s <= (n-1)/2)
      {
         carry1 = s ? -mpn_sub_n(b, a+n-s, b, s) : 0;
         carry2 = mpn_add_n(b+s, b+s, a, (n+1)/2);
         long i = ss-1;
         do b[i] = ~b[i]; while (--i >= s);
         b[n] = (ss < n) ? -mpn_sub_n(b+ss, b+ss, a+(n+1)/2, n-ss) : 0;
         signed_add_1(b+s, n-s+1, carry1 + a[n] + 1);
         signed_add_1(b+ss, n-ss+1, -Z - carry2 - 1);
      }
      else
      {
         carry1 = mpn_add_n(b, b, a+n-s, ss);
         carry2 = (n > 1) ? mpn_sub_n(b+ss, a+(n+1)/2, b+ss, (n-1)/2) : 0;
         b[n] = mpn_add_n(b+s, b+s, a, n-s);
         long i = n;
         do b[i] = ~b[i]; while (--i >= s);
         signed_add_1(b+ss, n-ss+1, carry1 + Z);
         signed_add_1(b+s, n-s+1, -carry2 + a[n] + 1);
      }
   }
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
      if (s <= n/2)
      {
         // We're computing B^s * (1 - B^(n/2)) * a.
         // If input is
         // 0                   n/2-s     n/2                   n-s       n
         // |         x0          |   y0   |         x1          |   y1   |
         // then output should be
         // 0        s                    n/2     n/2+s                   n
         // |  y0-y1 |        x0+x1        | y0+y1  |       -x0+x1        |

         // Store x1 - x0
         b[n] = (s < n/2) ? -mpn_sub_n(b+n/2+s, a+n/2, a, n/2-s) : 0;
         // Store x0 + x1 and y0 + y1
         carry = mpn_add_n(b+s, a, a+n/2, n/2);
         signed_add_1(b+s+n/2, n/2-s+1, carry + a[n]);
         // Store y0 - y1
         carry = s ? mpn_sub_n(b, a+n/2-s, a+n-s, s) : 0;
         signed_add_1(b+s, n-s+1, -a[n] - carry);
      }
      else
      {
         s -= n/2;

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
   }
   else
   {
      s -= n;

      if (s < n/2)
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
         carry = s ? mpn_sub_n(b, a+n-s, a+n/2-s, s) : 0;
         // (the +1 compensates for the bottom bit of the complement)
         signed_add_1(b+s, n-s+1, a[n] - carry + 1);
      }
      else
      {
         s -= n/2;
      
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
         carry = s ? mpn_add_n(b, a+n/2-s, a+n-s, s) : 0;
         signed_add_1(b+s, n-s+1, a[n] + carry);
      }
   }
}


void ZmodF_forward_butterfly_2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                  unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS);
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
   FLINT_ASSERT(s < 2*n*FLINT_BITS);
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
   FLINT_ASSERT(s < n*FLINT_BITS);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

   unsigned long bits = s & (FLINT_BITS - 1);
   if (bits)
      // shift right by leftover bits
      ZmodF_short_div_2exp(*b, *b, bits, n);

   s /= FLINT_BITS;
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
   FLINT_ASSERT(s < 2*n*FLINT_BITS);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *z);
   FLINT_ASSERT(*a != *z);

   if (s & 1)
   {
      unsigned long limbs, bits;
      ZmodF_decompose_rotation(&limbs, &bits, 2*n*FLINT_BITS - s, n);
      
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


void ZmodF_divby3(ZmodF_t b, ZmodF_t a, unsigned long n)
{
   // make overflow limb nonnegative
   ZmodF_fast_reduce(a, n);

   // compute a "total" which is congruent to a mod 3
   unsigned long hi = 0, lo = 0;
   for (unsigned long i = 0; i <= n; i++)
      add_ssaaaa(hi, lo, hi, lo, 0, a[i]);

   unsigned long total = lo & ((1UL << (FLINT_BITS/2)) - 1);
   total += (lo >> (FLINT_BITS/2));
   total += hi;

   // add "total" times B^n + 1 (the latter is 2 mod 3),
   // so that a becomes exactly divisible by 3
   mpn_add_1(a, a, n+1, total);
   a[n] += total;
   
   unsigned long rem = mpn_divexact_by3(b, a, n+1);
   FLINT_ASSERT(!rem);
}


// end of file ****************************************************************
