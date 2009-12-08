/*
   mulmid_ks.c:  polynomial middle products by Kronecker substitution
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "zn_poly_internal.h"


/*
   In the routines below, we denote by f1(x) and f2(x) the input polynomials
   op1[0, n1) and op2[0, n2), and by h(x) their product in Z[x].
   
   We write h(x) = LO(x) + x^(n2-1) * g(x) + x^n1 * HI(x), where
   len(LO) = len(HI) = n2 - 1 and len(g) = n1 - n2 + 1. Our goal is to
   compute the middle segment g.
   
   The basic strategy is: if X is an evaluation point (i.e. X = 2^b, -2^b,
   2^(-b) or -2^(-b)) then g(X) corresponds roughly to the integer middle
   product of f1(X) and f2(X), and we will use mpn_mulmid() to compute the
   latter.
   
   Unfortunately there are some complications.
   
   First, mpn_mulmid() works in terms of whole limb counts, not bit counts,
   and moreover the first two and last two limbs of the output of mpn_mulmid()
   are always garbage. We handle this issue using zero-padding as follows.
   Suppose that we need s bits of g(X) starting at bit index r. We compute
   f2(X) as usual. Let k2 = number of limbs used to store f2(X). Instead of
   evaluating f1(X), we evaluate 2^p * f1(X), i.e. zero-pad by p bits, where

      p = (k2 + 1) * GMP_NUMB_BITS - r.

   (We will verify in each case that p >= 0.) This shifts g(X) left by p bits,
   and ensures that bit #r of g(X) starts exactly at the first bit of the
   third limb of the output of mpn_mulmid(). Let k1 = number of limbs used to
   store f1(X). To be guaranteed of obtaining s correct bits of g(X), we need
   to have
      
      (k1 - k2 - 1) * GMP_NUMB_BITS >= s,
      
   or equivalently
   
      k1 * GMP_NUMB_BITS >= p + r + s.           (*)
      
   We zero-pad 2^p * f1(X) on the right to ensure that (*) holds. In every
   case, it turns out that the total amount of zero-padding is O(1) bits.
   
   Second, in the "reciprocal" variants (KS3 and KS4) there is the problem of
   overlapping coefficients, e.g. when we compute the integer middle product,
   the low bits of g(X) are polluted by the high bits of LO(X). To deal with
   this we need to compute the low coefficient of g(X) separately, and remove
   its effect from the overlapping portion. Similarly at the other end. The
   diagonal() function accomplishes this.

*/



/*
   Middle product using Kronecker substitution at 2^b.
*/
void
zn_array_mulmid_KS1 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n1 <= ULONG_MAX);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   // length of g
   size_t n3 = n1 - n2 + 1;

   // bits in each output coefficient
   unsigned b = 2 * mod->bits + ceil_lg (n2);

   // number of ulongs required to store each output coefficient
   unsigned w = CEIL_DIV (b, ULONG_BITS);
   ZNP_ASSERT (w <= 3);

   // number of limbs needed to store f2(2^b)
   size_t k2 = CEIL_DIV (n2 * b, GMP_NUMB_BITS);

   // We need r = (n2 - 1) * b and s = (n1 - n2 + 1) * b. Note that p is
   // non-negative since k2 * GMP_NUMB_BITS >= n2 * b.
   unsigned p = GMP_NUMB_BITS * (k2 + 1) - (n2 - 1) * b;
   
   // For (*) to hold we need k1 * GMP_NUMB_BITS >= p + n1 * b.
   size_t k1 = CEIL_DIV (p + n1 * b, GMP_NUMB_BITS);
   
   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 2 * k1 + 3);
   mp_limb_t* v1 = limbs;       // k1 limbs
   mp_limb_t* v2 = v1 + k1;     // k2 limbs
   mp_limb_t* v3 = v2 + k2;     // k1 - k2 + 3 limbs

   // evaluate 2^p * f1(2^b) and f2(2^b)
   zn_array_pack (v1, op1, n1, 1, b, p, 0);
   zn_array_pack (v2, op2, n2, 1, b, 0, 0);

   // compute segment of f1(2^b) * f2(2^b) starting at bit index r
   ZNP_mpn_mulmid (v3, v1, k1, v2, k2);

   // unpack coefficients of g, and reduce mod m
   ZNP_FASTALLOC (z, ulong, 6624, n3 * w);
   zn_array_unpack_SAFE (z, v3 + 2, n3, b, 0, k1 - k2 - 1);
   array_reduce (res, 1, z, n3, w, redc, mod);

   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}



/*
   Middle product using Kronecker substitution at 2^b and -2^b.
*/
void
zn_array_mulmid_KS2 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n1 <= ULONG_MAX);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   if (n2 == 1)
   {
      // code below needs n2 > 1, so fall back on scalar multiplication
      _zn_array_scalar_mul (res, op1, n1, op2[0], redc, mod);
      return;
   }

   // bits in each output coefficient
   unsigned bits = 2 * mod->bits + ceil_lg (n2);
   
   // we're evaluating at x = B and -B, where B = 2^b, and b = ceil(bits / 2)
   unsigned b = (bits + 1) / 2;

   // number of ulongs required to store each output coefficient
   unsigned w = CEIL_DIV (2 * b, ULONG_BITS);
   ZNP_ASSERT (w <= 3);
   
   // Write f1(x) = f1e(x^2) + x * f1o(x^2)
   //       f2(x) = f2e(x^2) + x * f2o(x^2)
   //        h(x) =  he(x^2) + x *  ho(x^2)
   //        g(x) =  ge(x^2) + x *  go(x^2)
   // "e" = even, "o" = odd

   // When evaluating f2e(B^2) and B * f2o(B^2) the bit-packing routine needs
   // room for the last chunk of 2b bits, so we need to allow room for
   // (n2 + 1) * b bits.
   size_t k2 = CEIL_DIV ((n2 + 1) * b, GMP_NUMB_BITS);

   // We need r = (n2 - 2) * b + 1 and s = (n1 - n2 + 3) * b.
   // Note that p is non-negative (since k2 * GMP_NUMB_BITS >= (n2 + 1) * b
   // >= (n2 - 2) * b - 1).
   unsigned p = GMP_NUMB_BITS * (k2 + 1) - (n2 - 2) * b - 1;

   // For (*) to hold we need k1 * GMP_NUMB_BITS >= p + (n1 + 1) * b + 1.
   // Also, to ensure that there is enough room for bit-packing (as for k2
   // above), we need k1 * GMP_NUMB_BITS >= p + (n1 + 1) * b; this condition
   // is subsumed by the first one.
   size_t k1 = CEIL_DIV (p + (n1 + 1) * b + 1, GMP_NUMB_BITS);
   
   size_t k3 = k1 - k2 + 3;
   ZNP_ASSERT (k3 >= 5);

   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 3 * k3 + 5 * k2);
   mp_limb_t* v2_buf0 = limbs;             // k2 limbs
   mp_limb_t* v3_buf0 = v2_buf0 + k2;      // k3 limbs
   mp_limb_t* v2_buf1 = v3_buf0 + k3;      // k2 limbs
   mp_limb_t* v3_buf1 = v2_buf1 + k2;      // k3 limbs
   mp_limb_t* v2_buf2 = v3_buf1 + k3;      // k2 limbs
   mp_limb_t* v3_buf2 = v2_buf2 + k2;      // k3 limbs
   mp_limb_t* v2_buf3 = v3_buf2 + k3;      // k2 limbs
   mp_limb_t* v2_buf4 = v2_buf3 + k2;      // k2 limbs

   // arrange overlapping buffers to minimise memory use
   // "p" = plus, "m" = minus
   mp_limb_t* v1e = v2_buf0;
   mp_limb_t* v1o = v2_buf2;
   mp_limb_t* v1p = v2_buf1;
   mp_limb_t* v1m = v2_buf0;
   mp_limb_t* v2e = v2_buf2;
   mp_limb_t* v2o = v2_buf4;
   mp_limb_t* v2p = v2_buf3;
   mp_limb_t* v2m = v2_buf2;
   mp_limb_t* v3m = v3_buf2;
   mp_limb_t* v3p = v3_buf0;
   mp_limb_t* v3e = v3_buf1;
   mp_limb_t* v3o = v3_buf1;

   // length of g
   size_t n3 = n1 - n2 + 1;
   
   ZNP_FASTALLOC (z, ulong, 6624, w * ((n3 + 1) / 2));
   
   // evaluate 2^p * f1e(B^2) and 2^p * B * f1o(B^2)
   zn_array_pack (v1e, op1, (n1 + 1) / 2, 2, 2 * b, p, k1);
   zn_array_pack (v1o, op1 + 1, n1 / 2, 2, 2 * b, p + b, k1);

   // compute   2^p * f1(B)   =  2^p * (f1e(B^2) + B * f1o(B^2))
   //     and  |2^p * f1(-B)| = |2^p * (f1e(B^2) - B * f1o(B^2))|
   // v3m_neg is set if f1(-B) is negative
   ZNP_ASSERT_NOCARRY (mpn_add_n (v1p, v1e, v1o, k1));
   int v3m_neg = signed_mpn_sub_n (v1m, v1e, v1o, k1);
   
   // evaluate f2e(B^2) and B * f2o(B^2)
   zn_array_pack (v2e, op2, (n2 + 1) / 2, 2, 2 * b, 0, k2);
   zn_array_pack (v2o, op2 + 1, n2 / 2, 2, 2 * b, b, k2);

   // compute    f2(B)   =   f2e(B^2) + B * f2o(B^2)
   //     and   |f2(-B)| =  |f2e(B^2) - B * f2o(B^2)|
   // v3m_neg is set if f1(-B) and f2(-B) have opposite signs
   ZNP_ASSERT_NOCARRY (mpn_add_n (v2p, v2e, v2o, k2));
   v3m_neg ^= signed_mpn_sub_n (v2m, v2e, v2o, k2);

   // compute segment starting at bit index r of
   //           h(B)   =  f1(B)   *  f2(B)
   //    and   |h(-B)| = |f1(-B)| * |f2(-B)|
   // v3m_neg is set if h(-B) is negative
   ZNP_mpn_mulmid (v3m, v1m, k1, v2m, k2);
   ZNP_mpn_mulmid (v3p, v1p, k1, v2p, k2);

   // compute segment starting at bit index r of
   //         2     * he(B^2) = h(B) + h(-B)       (if n2 is odd)
   //    or   2 * B * ho(B^2) = h(B) - h(-B)       (if n2 is even)
   // i.e. the segment of he(B^2) or B * ho(B^2) starting at bit index
   // r - 1 = (n2 - 2) * b. This encodes the coefficients of ge(x).
   
   // Note that when we do the addition (resp. subtraction) below, we might
   // miss a carry (resp. borrow) from the unknown previous limbs. We arrange
   // so that the answers are either correct or one too big, by adding 1
   // appropriately.

   if (v3m_neg ^ (n2 & 1))
   {
      mpn_add_n (v3e, v3p + 2, v3m + 2, k3 - 4);    // miss carry?
      mpn_add_1 (v3e, v3e, k3 - 4, 1);
   }
   else
      mpn_sub_n (v3e, v3p + 2, v3m + 2, k3 - 4);    // miss borrow?

   // Now we extract ge(x). The first coefficient we want is the coefficient
   // of x^(n2 - 1) in h(x); this starts at bit b index in v3e. We want
   // ceil(n3 / 2) coefficients altogether, with 2b bits each. This accounts
   // for the definition of s.

   // Claim: if we committed a "one-too-big" error above, this does not affect
   // the coefficients we extract. Proof: the first b bits of v3e are the top
   // half of the coefficient of x^(n2 - 2) in h(x). The base-B digit in those
   // b bits has value at most B - 2. Therefore adding 1 to it will never
   // overflow those b bits.

   zn_array_unpack_SAFE (z, v3e, (n3 + 1) / 2, 2 * b, b, k3 - 4);
   array_reduce (res, 2, z, (n3 + 1) / 2, w, redc, mod);

   // Now repeat all the above for go(x).
   
   if (v3m_neg ^ (n2 & 1))
      mpn_sub_n (v3o, v3p + 2, v3m + 2, k3 - 4);
   else
   {
      mpn_add_n (v3o, v3p + 2, v3m + 2, k3 - 4);
      mpn_add_1 (v3o, v3o, k3 - 4, 1);
   }

   zn_array_unpack_SAFE (z, v3o, n3 / 2, 2 * b, 2 * b, k3 - 4);
   array_reduce (res + 1, 2, z, n3 / 2, w, redc, mod);

   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}



/*
   Computes the sum
   
      op1[0] * op2[n-1] + ... + op1[n-1] * op2[0]
      
   as an *integer*. The result is assumed to fit into w ulongs (where
   1 <= w <= 3), and is written to res[0, w). The return value is the result
   reduced modulo mod->m (using redc if requested).
*/
#define diagonal_sum \
    ZNP_diagonal_sum
ulong
diagonal_sum (ulong* res, const ulong* op1, const ulong* op2,
              size_t n, unsigned w, int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n >= 1);
   ZNP_ASSERT (w >= 1);
   ZNP_ASSERT (w <= 3);

   size_t i;
   
   if (w == 1)
   {
      ulong sum = op1[0] * op2[n - 1];

      for (i = 1; i < n; i++)
         sum += op1[i] * op2[n - 1 - i];

      res[0] = sum;
      return redc ? zn_mod_reduce_redc (sum, mod) : zn_mod_reduce (sum, mod);
   }
   else if (w == 2)
   {
      ulong lo, hi, sum0, sum1;
      
      ZNP_MUL_WIDE (sum1, sum0, op1[0], op2[n - 1]);
   
      for (i = 1; i < n; i++)
      {
         ZNP_MUL_WIDE (hi, lo, op1[i], op2[n - 1 - i]);
         ZNP_ADD_WIDE (sum1, sum0, sum1, sum0, hi, lo);
      }
      
      res[0] = sum0;
      res[1] = sum1;
      return redc ? zn_mod_reduce2_redc (sum1, sum0, mod)
                  : zn_mod_reduce2 (sum1, sum0, mod);
   }
   else    // w == 3
   {
      ulong lo, hi, sum0, sum1, sum2 = 0;
      
      ZNP_MUL_WIDE (sum1, sum0, op1[0], op2[n - 1]);

      for (i = 1; i < n; i++)
      {
         ZNP_MUL_WIDE (hi, lo, op1[i], op2[n - 1 - i]);
         ZNP_ADD_WIDE (sum1, sum0, sum1, sum0, hi, lo);
         // carry into third limb:
         if (sum1 <= hi)
            sum2 += (sum1 < hi || sum0 < lo);
      }
      
      res[0] = sum0;
      res[1] = sum1;
      res[2] = sum2;
      return redc ? zn_mod_reduce3_redc (sum2, sum1, sum0, mod)
                  : zn_mod_reduce3 (sum2, sum1, sum0, mod);
   }
}


/*
   Inplace subtract 2^i*x from res[0, n).
   x is an array of w ulongs, where 1 <= w <= 3.
   i may be any non-negative integer.
*/
#define subtract_ulongs \
    ZNP_subtract_ulongs
void
subtract_ulongs (mp_limb_t* res, size_t n, size_t i, ulong* x, unsigned w)
{
   ZNP_ASSERT (w >= 1);
   ZNP_ASSERT (w <= 3);

#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS
   size_t k = i / GMP_NUMB_BITS;
   
   if (k >= n)
      return;

   unsigned j = i % GMP_NUMB_BITS;
   
   if (j == 0)
      mpn_sub (res + k, res + k, n - k, x, ZNP_MIN (n - k, w));
   else
   {
      mp_limb_t y[4];
      y[w] = mpn_lshift (y, x, w, j);
      mpn_sub (res + k, res + k, n - k, y, ZNP_MIN (n - k, w + 1));
   }
#else
#error Not nails-safe yet
#endif
}


/*
   Middle product using Kronecker substitution at 2^b and 2^(-b).

   Note: this routine does not appear to be competitive in practice with the
   other KS routines. It's here just for fun.
*/
void
zn_array_mulmid_KS3 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n1 <= ULONG_MAX);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   // length of g
   size_t n3 = n1 - n2 + 1;

   // bits in each output coefficient
   unsigned bits = 2 * mod->bits + ceil_lg (n2);
   
   // we're evaluating at x = B and 1/B, where B = 2^b, and b = ceil(bits / 2)
   unsigned b = (bits + 1) / 2;

   // number of ulongs required to store each base-B digit
   unsigned w = CEIL_DIV (b, ULONG_BITS);
   ZNP_ASSERT (w <= 2);
   
   // number of ulongs needed to store each output coefficient
   unsigned ww = CEIL_DIV (2 * b, ULONG_BITS);
   ZNP_ASSERT (ww <= 3);

   // directly compute coefficient of x^0 in g(x)
   ulong dlo[3];
   res[0] = diagonal_sum (dlo, op1, op2, n2, ww, redc, mod);
   if (n3 == 1)
      return;      // only need one coefficient of output

   // directly compute coefficient of x^(n3-1) in g(x)
   ulong dhi[3];
   res[n3 - 1] = diagonal_sum (dhi, op1 + n3 - 1, op2, n2, ww, redc, mod);
   if (n3 == 2)
      return;      // only need two coefficients of output

   // limbs needed to store f2(B) and B^(n2-1) * f2(1/B)
   size_t k2 = CEIL_DIV (n2 * b, GMP_NUMB_BITS);
   
   // we need r = (n2 - 1) * b and s = (n1 - n2 + 1) * b, thus p is:
   unsigned p = GMP_NUMB_BITS * (k2 + 1) - (n2 - 1) * b;

   // for (*) we need k1 * GMP_NUMB_BITS >= p + n1 * b
   size_t k1 = CEIL_DIV (p + n1 * b, GMP_NUMB_BITS);

   size_t k3 = k1 - k2 + 3;
   ZNP_ASSERT (k3 >= 5);

   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 2 * k1 + 3);
   mp_limb_t* v1 = limbs;        // k1 limbs
   mp_limb_t* v2 = v1 + k1;      // k2 limbs
   mp_limb_t* v3 = v2 + k2;      // k1 - k2 + 3 limbs

   ZNP_FASTALLOC (z, ulong, 6624, 2 * w * n3);
   // "n" = normal order, "r" = reciprocal order
   ulong* zn = z;
   ulong* zr = z + w * n3;

   // -------------------------------------------------------------------------
   //     "normal" evaluation point

   // evaluate 2^p * f1(B) and f2(B)
   zn_array_pack (v1, op1, n1, 1, b, p, k1);
   zn_array_pack (v2, op2, n2, 1, b, 0, k2);
   
   // compute segment starting at bit index r of h(B) = f1(B) * f2(B)
   ZNP_mpn_mulmid (v3, v1, k1, v2, k2);
   
   // remove x^0 and x^(n3 - 1) coefficient of g(x)
   subtract_ulongs (v3 + 2, k3 - 4, 0, dlo, ww);
   subtract_ulongs (v3 + 2, k3 - 4, (n3 - 1) * b, dhi, ww);

   // decompose relevant portion of h(B) into base-B digits
   zn_array_unpack_SAFE (zn, v3 + 2, n3 - 1, b, b, k3 - 4);
   
   // At this stage zn contains (n3 - 1) base-B digits, representing the
   // integer g[1] + g[2]*B + ... + g[n3-2]*B^(n3-3)

   // -------------------------------------------------------------------------
   //     "reciprocal" evaluation point
   
   // evaluate 2^p * B^(n1-1) * f1(1/B) and B^(n2-1) * f2(B)
   zn_array_pack (v1, op1 + n1 - 1, n1, -1, b, p, k1);
   zn_array_pack (v2, op2 + n2 - 1, n2, -1, b, 0, k2);
   
   // compute segment starting at bit index r of B^(n1+n2-2) * h(1/B) =
   // (B^(n1-1) * f1(1/B)) * (B^(n2-1) * f2(1/B))
   ZNP_mpn_mulmid (v3, v1, k1, v2, k2);

   // remove x^0 and x^(n3 - 1) coefficient of g(x)
   subtract_ulongs (v3 + 2, k3 - 4, 0, dhi, ww);
   subtract_ulongs (v3 + 2, k3 - 4, (n3 - 1) * b, dlo, ww);
   
   // decompose relevant portion of B^(n1+n2-2) * h(1/B) into base-B digits
   zn_array_unpack_SAFE (zr, v3 + 2, n3 - 1, b, b, k3 - 4);

   // At this stage zr contains (n3 - 1) base-B digits, representing the
   // integer g[n3-2] + g[n3-3]*B + ... + g[1]*B^(n3-3)
   
   // -------------------------------------------------------------------------
   //     combine "normal" and "reciprocal" information

   zn_array_recover_reduce (res + 1, 1, zn, zr, n3 - 2, b, redc, mod);

   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}



/*
   Middle product using Kronecker substitution at 2^b, -2^b, 2^(-b)
   and -2^(-b).
*/
void
zn_array_mulmid_KS4 (ulong* res, 
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n1 <= ULONG_MAX);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   if (n2 == 1)
   {
      // code below needs n2 > 1, so fall back on scalar multiplication
      _zn_array_scalar_mul (res, op1, n1, op2[0], redc, mod);
      return;
   }

   // bits in each output coefficient
   unsigned bits = 2 * mod->bits + ceil_lg (n2);

   // we're evaluating at x = B, -B, 1/B, -1/B,
   // where B = 2^b, and b = ceil(bits / 4)
   unsigned b = (bits + 3) / 4;

   // number of ulongs required to store each base-B^2 digit
   unsigned w = CEIL_DIV (2 * b, ULONG_BITS);
   ZNP_ASSERT (w <= 2);

   // number of ulongs needed to store each output coefficient
   unsigned ww = CEIL_DIV (4 * b, ULONG_BITS);
   ZNP_ASSERT (ww <= 3);

   // mask = 2^c - 1, where c = number of bits used in high ulong of each
   // base-B^2 digit
   ulong mask;
   if (w == 1)
      mask = ((2 * b) < ULONG_BITS) ? ((1UL << (2 * b)) - 1) : (-1UL);
   else   // w == 2
      mask = (1UL << ((2 * b) - ULONG_BITS)) - 1;

   // Write f1(x) = f1e(x^2) + x * f1o(x^2)
   //       f2(x) = f2e(x^2) + x * f2o(x^2)
   //        h(x) =  he(x^2) + x *  ho(x^2)
   //        g(x) =  ge(x^2) + x *  go(x^2)
   // "e" = even, "o" = odd

   size_t n1o = n1 / 2;
   size_t n1e = n1 - n1o;

   size_t n2o = n2 / 2;
   size_t n2e = n2 - n2o;

   size_t n3 = n1 - n2 + 1;   // length of g
   size_t n3o = n3 / 2;
   size_t n3e = n3 - n3o;

   // directly compute coefficient of x^0 in ge(x)
   ulong delo[3];
   res[0] = diagonal_sum (delo, op1, op2, n2, ww, redc, mod);
   if (n3 == 1)
      return;      // only need one coefficient of output

   // directly compute coefficient of x^0 in go(x)
   ulong dolo[3];
   res[1] = diagonal_sum (dolo, op1 + 1, op2, n2, ww, redc, mod);
   if (n3 == 2)
      return;      // only need two coefficients of output

   // directly compute coefficient of x^(n3e - 1) in ge(x)
   ulong dehi[3];
   res[2*n3e - 2] = diagonal_sum (dehi, op1 + 2*n3e - 2, op2,
                                  n2, ww, redc, mod);
   if (n3 == 3)
      return;      // only need three coefficients of output

   // directly compute coefficient of x^(n3o - 1) in go(x)
   ulong dohi[3];
   res[2*n3o - 1] = diagonal_sum (dohi, op1 + 2*n3o - 1, op2,
                                  n2, ww, redc, mod);
   if (n3 == 4)
      return;      // only need four coefficients of output

   // In f2(B), the leading coefficient starts at bit position b * (n2 - 1)
   // and has length 2*b, and the coefficients overlap so we need an extra bit
   // for the carry: this gives (n2 + 1) * b + 1 bits.
   size_t k2 = CEIL_DIV ((n2 + 1) * b + 1, GMP_NUMB_BITS);

   // We need r = (n2 - 1) * b + 1 and s = (n3 + 1) * b.
   // Note that p is non-negative (since k2 * GMP_NUMB_BITS >= (n2 + 1) * b
   // >= (n2 - 1) * b - 1).
   unsigned p = GMP_NUMB_BITS * (k2 + 1) - (n2 - 1) * b - 1;

   // For (*) we need k1 * GMP_NUMB_BITS >= p + (n1 + 1) * b + 1.
   size_t k1 = CEIL_DIV (p + (n1 + 1) * b + 1, GMP_NUMB_BITS);

   size_t k3 = k1 - k2 + 3;
   ZNP_ASSERT (k3 >= 5);

   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 5 * (k2 + k3));
   mp_limb_t* v2_buf0 = limbs;             // k2 limbs
   mp_limb_t* v3_buf0 = v2_buf0 + k2;      // k3 limbs
   mp_limb_t* v2_buf1 = v3_buf0 + k3;      // k2 limbs
   mp_limb_t* v3_buf1 = v2_buf1 + k2;      // k3 limbs
   mp_limb_t* v2_buf2 = v3_buf1 + k3;      // k2 limbs
   mp_limb_t* v3_buf2 = v2_buf2 + k2;      // k3 limbs
   mp_limb_t* v2_buf3 = v3_buf2 + k3;      // k2 limbs
   mp_limb_t* v3_buf3 = v2_buf3 + k2;      // k3 limbs
   mp_limb_t* v2_buf4 = v3_buf3 + k3;      // k2 limbs
   mp_limb_t* v3_buf4 = v2_buf4 + k2;      // k3 limbs
   
   // arrange overlapping buffers to minimise memory use
   // "p" = plus, "m" = minus
   // "n" = normal order, "r" = reciprocal order
   mp_limb_t* v1en = v2_buf1;
   mp_limb_t* v1on = v2_buf2;
   mp_limb_t* v1pn = v2_buf0;
   mp_limb_t* v1mn = v2_buf1;
   mp_limb_t* v2en = v2_buf3;
   mp_limb_t* v2on = v2_buf4;
   mp_limb_t* v2pn = v2_buf2;
   mp_limb_t* v2mn = v2_buf3;
   mp_limb_t* v3mn = v3_buf2;
   mp_limb_t* v3pn = v3_buf3;
   mp_limb_t* v3en = v3_buf4;
   mp_limb_t* v3on = v3_buf3;

   mp_limb_t* v1er = v2_buf1;
   mp_limb_t* v1or = v2_buf2;
   mp_limb_t* v1pr = v2_buf0;
   mp_limb_t* v1mr = v2_buf1;
   mp_limb_t* v2er = v2_buf3;
   mp_limb_t* v2or = v2_buf4;
   mp_limb_t* v2pr = v2_buf2;
   mp_limb_t* v2mr = v2_buf3;
   mp_limb_t* v3mr = v3_buf2;
   mp_limb_t* v3pr = v3_buf1;
   mp_limb_t* v3er = v3_buf0;
   mp_limb_t* v3or = v3_buf1;

   ZNP_FASTALLOC (z, ulong, 6624, 2 * w * (n3e - 1));
   ulong* zn = z;
   ulong* zr = z + w * (n3e - 1);

   int v3m_neg;
   
   // -------------------------------------------------------------------------
   //     "normal" evaluation point

   // evaluate 2^p * f1e(B^2) and 2^p * B * f1o(B^2)
   zn_array_pack (v1en, op1, n1e, 2, 2 * b, p, k1);
   zn_array_pack (v1on, op1 + 1, n1o, 2, 2 * b, p + b, k1);

   // compute 2^p *   f1(B)  =  2^p * f1e(B^2) + 2^p * B * f1o(B^2)
   //    and  2^p * |f1(-B)| = |2^p * f1e(B^2) - 2^p * B * f1o(B^2)|
   ZNP_ASSERT_NOCARRY (mpn_add_n (v1pn, v1en, v1on, k1));
   v3m_neg = signed_mpn_sub_n (v1mn, v1en, v1on, k1);

   // evaluate f2e(B^2) and B * f2o(B^2)
   zn_array_pack (v2en, op2, n2e, 2, 2 * b, 0, k2);
   zn_array_pack (v2on, op2 + 1, n2o, 2, 2 * b, b, k2);

   // compute   f2(B)  =  f2e(B^2) + B * f2o(B^2)
   //    and  |f2(-B)| = |f2e(B^2) - B * f2o(B^2)|
   ZNP_ASSERT_NOCARRY (mpn_add_n (v2pn, v2en, v2on, k2));
   v3m_neg ^= signed_mpn_sub_n (v2mn, v2en, v2on, k2);

   // compute segment starting at bit index r of
   //            h(B)  =   f1(B)  *  f2(B)
   //    and   |h(-B)| = |f1(-B)| * |f2(-B)|
   // hn_neg is set if h(-B) is negative
   ZNP_mpn_mulmid (v3mn, v1mn, k1, v2mn, k2);
   ZNP_mpn_mulmid (v3pn, v1pn, k1, v2pn, k2);

   // compute segments starting at bit index r of
   //         2     * he(B^2) = h(B) + h(-B)
   //    and  2 * B * ho(B^2) = h(B) - h(-B)
   // ie. the segments of he(B^2) and B * ho(B^2) starting at bit index r - 1.

   // If n2 is odd, the former encodes ge(x) and the latter encodes go(x).
   // Otherwise the situation is reversed. We write the results to v3en/v3on
   // accordingly.
   
   // Note that when we do the addition (resp. subtraction) below, we might
   // miss a carry (resp. borrow) from the unknown previous limbs. We arrange
   // so that the answers are either correct or one too big, by adding 1
   // appropriately.

   if (v3m_neg ^ (n2 & 1))
   {
      mpn_add_n (v3en + 2, v3pn + 2, v3mn + 2, k3 - 4);    // miss carry?
      mpn_add_1 (v3en + 2, v3en + 2, k3 - 4, 1);
      mpn_sub_n (v3on + 2, v3pn + 2, v3mn + 2, k3 - 4);    // miss borrow?
   }
   else
   {
      mpn_sub_n (v3en + 2, v3pn + 2, v3mn + 2, k3 - 4);
      mpn_add_n (v3on + 2, v3pn + 2, v3mn + 2, k3 - 4);
      mpn_add_1 (v3on + 2, v3on + 2, k3 - 4, 1);
   }

   // remove x^0 and x^(n3e - 1) coefficients of ge(x),
   //    and x^0 and x^(n3o - 1) coefficients of go(x).
   subtract_ulongs (v3en + 2, k3 - 4, 0, delo, ww);
   subtract_ulongs (v3en + 2, k3 - 4, (2*n3e - 2) * b, dehi, ww);
   subtract_ulongs (v3on + 2, k3 - 4, b, dolo, ww);
   subtract_ulongs (v3on + 2, k3 - 4, (2*n3o - 1) * b, dohi, ww);

   // At this stage, the integer
   //   g[2] + g[4]*B^2 + ... + g[2*n3e - 4]*B^(2*n3e - 6)
   // appears in v3en + 2, starting at bit index 2*b, occupying
   // (2 * n3e - 2) * b bits. The integer
   //   g[3] + g[5]*B^2 + ... + g[2*n3o - 3]*B^(2*n3o - 6)
   // appears in v3on + 2, starting at bit index 3*b, occupying
   // (2 * n3o - 2) * b bits.

   // -------------------------------------------------------------------------
   //     "reciprocal" evaluation point

   // evaluate 2^p * B^(n1-1) * f1e(1/B^2) and 2^p * B^(n1-2) * f1o(B^2)
   zn_array_pack (v1er, op1 + 2*(n1e - 1), n1e, -2, 2 * b,
                  p + ((n1 & 1) ? 0 : b), k1);
   zn_array_pack (v1or, op1 + 1 + 2*(n1o - 1), n1o, -2, 2 * b,
                  p + ((n1 & 1) ? b : 0), k1);

   // compute 2^p * B^(n1-1) * f1(1/B) =
   //                2^p * B^(n1-1) * f1e(1/B^2) + 2^p * B^(n1-2) * f1o(B^2)
   //   and  |2^p * B^(n1-1) * f1(-1/B)| =
   //               |2^p * B^(n1-1) * f1e(1/B^2) - 2^p * B^(n1-2) * f1o(B^2)|
   ZNP_ASSERT_NOCARRY (mpn_add_n (v1pr, v1er, v1or, k1));
   v3m_neg = signed_mpn_sub_n (v1mr, v1er, v1or, k1);

   // evaluate B^(n2-1) * f2e(1/B^2) and B^(n2-2) * f2o(B^2)
   zn_array_pack (v2er, op2 + 2*(n2e - 1), n2e, -2, 2 * b,
                  (n2 & 1) ? 0 : b, k2);
   zn_array_pack (v2or, op2 + 1 + 2*(n2o - 1), n2o, -2, 2 * b,
                  (n2 & 1) ? b : 0, k2);

   // compute B^(n2-1) * f2(1/B) =
   //                B^(n2-1) * f2e(1/B^2) + B^(n2-2) * f2o(B^2)
   //   and  |B^(n2-1) * f2(-1/B)| =
   //               |B^(n2-1) * f2e(1/B^2) - B^(n2-2) * f2o(B^2)|
   ZNP_ASSERT_NOCARRY (mpn_add_n (v2pr, v2er, v2or, k2));
   v3m_neg ^= signed_mpn_sub_n (v2mr, v2er, v2or, k2);

   // compute segment starting at bit index r of
   //       B^(n3-1) * h(1/B)   = (B^(n1-1) * f1(1/B))  * (B^(n2-1) * f2(1/B))
   //  and |B^(n3-1) * h(-1/B)| = |B^(n1-1) * f1(-1/B)| * |B^(n2-1) * f2(-1/B)|
   // hr_neg is set if h(-1/B) is negative
   ZNP_mpn_mulmid (v3mr, v1mr, k1, v2mr, k2);
   ZNP_mpn_mulmid (v3pr, v1pr, k1, v2pr, k2);

   // compute segments starting at bit index r of
   //        2 * B^(n3-1) * he(1/B^2) = B^(n3-1) * h(1/B) + B^(n3-1) * h(-1/B)
   //   and  2 * B^(n3-2) * ho(1/B^2) = B^(n3-1) * h(1/B) - B^(n3-1) * h(-1/B)
   // ie. the segments of B^(n3-1) * he(1/B^2) and B^(n3-2) * ho(1/B^2)
   // starting at bit index r - 1.

   // If n2 is odd, the former encodes ge(x) and the latter encodes go(x).
   // Otherwise the situation is reversed. We write the results to v3er/v3or
   // accordingly.

   if (v3m_neg ^ (n2 & 1))
   {
      mpn_add_n (v3er + 2, v3pr + 2, v3mr + 2, k3 - 4);    // miss carry?
      mpn_add_1 (v3er + 2, v3er + 2, k3 - 4, 1);
      mpn_sub_n (v3or + 2, v3pr + 2, v3mr + 2, k3 - 4);    // miss borrow?
   }
   else
   {
      mpn_sub_n (v3er + 2, v3pr + 2, v3mr + 2, k3 - 4);
      mpn_add_n (v3or + 2, v3pr + 2, v3mr + 2, k3 - 4);
      mpn_add_1 (v3or + 2, v3or + 2, k3 - 4, 1);
   }
   
   unsigned s = (n3 & 1) ? 0 : b;

   // remove x^0 and x^(n3e - 1) coefficients of ge(x),
   //    and x^0 and x^(n3o - 1) coefficients of go(x).
   subtract_ulongs (v3er + 2, k3 - 4, s, dehi, ww);
   subtract_ulongs (v3er + 2, k3 - 4, (2*n3e - 2) * b + s, delo, ww);
   subtract_ulongs (v3or + 2, k3 - 4, b - s, dohi, ww);
   subtract_ulongs (v3or + 2, k3 - 4, (2*n3o - 2) * b + b - s, dolo, ww);

   // At this stage, the integer
   //   g[2*n3e - 4] + g[2*n3e - 6]*B^2 + ... + g[2]*B^(2*n3e - 6)
   // appears in v3er + 2, starting at bit index 2*b if n3 is odd, or 3*b
   // if n3 is even, and occupying (2 * n3e - 2) * b bits. The integer
   //   g[2*n3o - 3] + g[2*n3o - 5]*B^2 + ... + g[3]*B^(2*n3o - 6)
   // appears in v3or + 2, starting at bit index 3*b if n3 is odd, or 2*b
   // if n3 is even, and occupying (2 * n3o - 2) * b bits.

   // -------------------------------------------------------------------------
   //     combine "normal" and "reciprocal" information

   // decompose relevant portion of ge(B^2) and ge(1/B^2) into base-B^2 digits
   zn_array_unpack_SAFE (zn, v3en + 2, n3e - 1, 2 * b, 2 * b, k3 - 4);
   zn_array_unpack_SAFE (zr, v3er + 2, n3e - 1, 2 * b, 2 * b + s, k3 - 4);

   // combine ge(B^2) and ge(1/B^2) information to get even coefficients of g
   zn_array_recover_reduce (res + 2, 2, zn, zr, n3e - 2, 2 * b, redc, mod);

   // decompose relevant portion of go(B^2) and go(1/B^2) into base-B^2 digits
   zn_array_unpack_SAFE (zn, v3on + 2, n3o - 1, 2 * b, 3 * b, k3 - 4);
   zn_array_unpack_SAFE (zr, v3or + 2, n3o - 1, 2 * b, 3 * b - s, k3 - 4);

   // combine go(B^2) and go(1/B^2) information to get odd coefficients of g
   zn_array_recover_reduce (res + 3, 2, zn, zr, n3o - 2, 2 * b, redc, mod);

   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}


// end of file ****************************************************************
