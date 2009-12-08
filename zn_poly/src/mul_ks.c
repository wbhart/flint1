/*
   mul_ks.c:  polynomial multiplication by Kronecker substitution
   
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
*/


/*
   Multiplication/squaring using Kronecker substitution at 2^b.
*/
void
zn_array_mul_KS1 (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n1 <= ULONG_MAX);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   int sqr = (op1 == op2 && n1 == n2);

   // length of h
   size_t n3 = n1 + n2 - 1;
   
   // bits in each output coefficient
   unsigned b = 2 * mod->bits + ceil_lg (n2);
   
   // number of ulongs required to store each output coefficient
   unsigned w = CEIL_DIV (b, ULONG_BITS);
   ZNP_ASSERT (w <= 3);

   // number of limbs needed to store f1(2^b) and f2(2^b)
   size_t k1 = CEIL_DIV (n1 * b, GMP_NUMB_BITS);
   size_t k2 = CEIL_DIV (n2 * b, GMP_NUMB_BITS);

   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 2 * (k1 + k2));
   mp_limb_t* v1 = limbs;     // k1 limbs
   mp_limb_t* v2 = v1 + k1;   // k2 limbs
   mp_limb_t* v3 = v2 + k2;   // k1 + k2 limbs

   if (!sqr)
   {
      // multiplication version

      // evaluate f1(2^b) and f2(2^b)
      zn_array_pack (v1, op1, n1, 1, b, 0, 0);
      zn_array_pack (v2, op2, n2, 1, b, 0, 0);

      // compute h(2^b) = f1(2^b) * f2(2^b)
      ZNP_mpn_mul (v3, v1, k1, v2, k2);
   }
   else
   {
      // squaring version

      // evaluate f1(2^b)
      zn_array_pack (v1, op1, n1, 1, b, 0, 0);

      // compute h(2^b) = f1(2^b)^2
      ZNP_mpn_mul (v3, v1, k1, v1, k1);
   }

   // unpack coefficients of h, and reduce mod m
   ZNP_FASTALLOC (z, ulong, 6624, n3 * w);
   zn_array_unpack_SAFE (z, v3, n3, b, 0, k1 + k2);
   array_reduce (res, 1, z, n3, w, redc, mod);

   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}



/*
   Multiplication/squaring using Kronecker substitution at 2^b and -2^b.
*/
void
zn_array_mul_KS2 (ulong* res,
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

   int sqr = (op1 == op2 && n1 == n2);

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
   // "e" = even, "o" = odd

   size_t n1o = n1 / 2;
   size_t n1e = n1 - n1o;

   size_t n2o = n2 / 2;
   size_t n2e = n2 - n2o;

   size_t n3 = n1 + n2 - 1;    // length of h
   size_t n3o = n3 / 2;
   size_t n3e = n3 - n3o;

   // f1(B) and |f1(-B)| are at most ((n1 - 1) * b + mod->bits) bits long.
   // However, when evaluating f1e(B^2) and B * f1o(B^2) the bitpacking
   // routine needs room for the last chunk of 2b bits. Therefore we need to
   // allow room for (n1 + 1) * b bits. Ditto for f2.
   size_t k1 = CEIL_DIV ((n1 + 1) * b, GMP_NUMB_BITS);
   size_t k2 = CEIL_DIV ((n2 + 1) * b, GMP_NUMB_BITS);
   size_t k3 = k1 + k2;

   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 3 * k3);
   mp_limb_t* v1_buf0 = limbs;             // k1 limbs
   mp_limb_t* v2_buf0 = v1_buf0 + k1;      // k2 limbs
   mp_limb_t* v1_buf1 = v2_buf0 + k2;      // k1 limbs
   mp_limb_t* v2_buf1 = v1_buf1 + k1;      // k2 limbs
   mp_limb_t* v1_buf2 = v2_buf1 + k2;      // k1 limbs
   mp_limb_t* v2_buf2 = v1_buf2 + k1;      // k2 limbs

   // arrange overlapping buffers to minimise memory use
   // "p" = plus, "m" = minus
   mp_limb_t* v1e = v1_buf0;
   mp_limb_t* v2e = v2_buf0;
   mp_limb_t* v1o = v1_buf1;
   mp_limb_t* v2o = v2_buf1;
   mp_limb_t* v1p = v1_buf2;
   mp_limb_t* v2p = v2_buf2;
   mp_limb_t* v1m = v1_buf0;
   mp_limb_t* v2m = v2_buf0;
   mp_limb_t* v3m = v1_buf1;
   mp_limb_t* v3p = v1_buf0;
   mp_limb_t* v3e = v1_buf2;
   mp_limb_t* v3o = v1_buf0;
   
   ZNP_FASTALLOC (z, ulong, 6624, w * n3e);

   int v3m_neg;

   if (!sqr)
   {
      // multiplication version

      // evaluate f1e(B^2) and B * f1o(B^2)
      zn_array_pack (v1e, op1, n1e, 2, 2 * b, 0, k1);
      zn_array_pack (v1o, op1 + 1, n1o, 2, 2 * b, b, k1);

      // evaluate f2e(B^2) and B * f2o(B^2)
      zn_array_pack (v2e, op2, n2e, 2, 2 * b, 0, k2);
      zn_array_pack (v2o, op2 + 1, n2o, 2, 2 * b, b, k2);

      // compute f1(B) = f1e(B^2) + B * f1o(B^2)
      //     and f2(B) = f2e(B^2) + B * f2o(B^2)
      ZNP_ASSERT_NOCARRY (mpn_add_n (v1p, v1e, v1o, k1));
      ZNP_ASSERT_NOCARRY (mpn_add_n (v2p, v2e, v2o, k2));

      // compute |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
      //     and |f2(-B)| = |f2e(B^2) - B * f2o(B^2)|
      v3m_neg  = signed_mpn_sub_n (v1m, v1e, v1o, k1);
      v3m_neg ^= signed_mpn_sub_n (v2m, v2e, v2o, k2);

      // compute  h(B)   =  f1(B)   *  f2(B)
      // compute |h(-B)| = |f1(-B)| * |f2(-B)|
      // v3m_neg is set if h(-B) is negative
      ZNP_mpn_mul (v3m, v1m, k1, v2m, k2);
      ZNP_mpn_mul (v3p, v1p, k1, v2p, k2);
   }
   else
   {
      // squaring version

      // evaluate f1e(B^2) and B * f1o(B^2)
      zn_array_pack (v1e, op1, n1e, 2, 2 * b, 0, k1);
      zn_array_pack (v1o, op1 + 1, n1o, 2, 2 * b, b, k1);

      // compute f1(B) = f1e(B^2) + B * f1o(B^2)
      ZNP_ASSERT_NOCARRY (mpn_add_n (v1p, v1e, v1o, k1));

      // compute |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
      signed_mpn_sub_n (v1m, v1e, v1o, k1);

      // compute h(B)  = f1(B)^2
      // compute h(-B) = f1(-B)^2
      // v3m_neg is cleared (since f1(-B)^2 is never negative)
      ZNP_mpn_mul (v3m, v1m, k1, v1m, k1);
      ZNP_mpn_mul (v3p, v1p, k1, v1p, k1);
      v3m_neg = 0;
   }
   
   // he(B^2) and B * ho(B^2) are both at most b * (n3 + 1) bits long (since
   // the coefficients don't overlap). The buffers used below are at least
   // b * (n1 + n2 + 2) = b * (n3 + 3) bits long. So we definitely have
   // enough room for 2 * he(B^2) and 2 * B * ho(B^2).

   // compute 2 * he(B^2) = h(B) + h(-B)
   ZNP_ASSERT_NOCARRY (v3m_neg ? mpn_sub_n (v3e, v3p, v3m, k3)
                               : mpn_add_n (v3e, v3p, v3m, k3));

   // unpack coefficients of he, and reduce mod m
   zn_array_unpack_SAFE (z, v3e, n3e, 2 * b, 1, k3);
   array_reduce (res, 2, z, n3e, w, redc, mod);
   
   // compute 2 * b * ho(B^2) = h(B) - h(-B)
   ZNP_ASSERT_NOCARRY (v3m_neg ? mpn_add_n (v3o, v3p, v3m, k3)
                               : mpn_sub_n (v3o, v3p, v3m, k3));
   
   // unpack coefficients of ho, and reduce mod m
   zn_array_unpack_SAFE (z, v3o, n3o, 2 * b, b + 1, k3);
   array_reduce (res + 1, 2, z, n3o, w, redc, mod);

   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}



/*
   Multiplication/squaring using Kronecker substitution at 2^b and 2^(-b).
   
   Note: this routine does not appear to be competitive in practice with the
   other KS routines. It's here just for fun.
*/
void
zn_array_mul_KS3 (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n1 <= ULONG_MAX);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   int sqr = (op1 == op2 && n1 == n2);

   // length of h
   size_t n3 = n1 + n2 - 1;
   
   // bits in each output coefficient
   unsigned bits = 2 * mod->bits + ceil_lg (n2);
   
   // we're evaluating at x = B and 1/B, where B = 2^b, and b = ceil(bits / 2)
   unsigned b = (bits + 1) / 2;

   // number of ulongs required to store each base-B digit
   unsigned w = CEIL_DIV (b, ULONG_BITS);
   ZNP_ASSERT (w <= 2);
   
   // limbs needed to store f1(B) and B^(n1-1) * f1(1/B), ditto for f2
   size_t k1 = CEIL_DIV (n1 * b, GMP_NUMB_BITS);
   size_t k2 = CEIL_DIV (n2 * b, GMP_NUMB_BITS);
   
   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 2 * (k1 + k2));
   mp_limb_t* v1 = limbs;       // k1 limbs
   mp_limb_t* v2 = v1 + k1;     // k2 limbs
   mp_limb_t* v3 = v2 + k2;     // k1 + k2 limbs

   ZNP_FASTALLOC (z, ulong, 6624, 2 * w * (n3 + 1));
   // "n" = normal order, "r" = reciprocal order
   ulong* zn = z;
   ulong* zr = z + w * (n3 + 1);

   if (!sqr)
   {
      // multiplication version

      // evaluate f1(B) and f2(B)
      zn_array_pack (v1, op1, n1, 1, b, 0, k1);
      zn_array_pack (v2, op2, n2, 1, b, 0, k2);

      // compute h(B) = f1(B) * f2(B)
      ZNP_mpn_mul (v3, v1, k1, v2, k2);
   }
   else
   {
      // squaring version

      // evaluate f1(B)
      zn_array_pack (v1, op1, n1, 1, b, 0, k1);

      // compute h(B) = f1(B)^2
      ZNP_mpn_mul (v3, v1, k1, v1, k1);
   }

   // decompose h(B) into base-B digits
   zn_array_unpack_SAFE (zn, v3, n3 + 1, b, 0, k1 + k2);

   if (!sqr)
   {
      // multiplication version

      // evaluate B^(n1-1) * f1(1/B) and B^(n2-1) * f2(1/B)
      zn_array_pack (v1, op1 + n1 - 1, n1, -1, b, 0, k1);
      zn_array_pack (v2, op2 + n2 - 1, n2, -1, b, 0, k2);

      // compute B^(n1+n2-2) * h(1/B) =
      //                     (B^(n1-1) * f1(1/B)) * (B^(n2-1) * f2(1/B))
      ZNP_mpn_mul (v3, v1, k1, v2, k2);
   }
   else
   {
      // squaring version

      // evaluate B^(n1-1) * f1(1/B)
      zn_array_pack (v1, op1 + n1 - 1, n1, -1, b, 0, k1);

      // compute B^(2*n1-2) * h(1/B) = (B^(n1-1) * f1(1/B))^2
      ZNP_mpn_mul (v3, v1, k1, v1, k1);
   }

   // decompose h(1/B) into base-B digits
   zn_array_unpack_SAFE (zr, v3, n3 + 1, b, 0, k1 + k2);

   // recover h(x) from h(B) and h(1/B)
   // (note: need to check that the high digit of each output coefficient
   // is < B - 1; this follows from an estimate in section 3.2 of [Har07].)
   zn_array_recover_reduce (res, 1, zn, zr, n3, b, redc, mod);
   
   ZNP_FASTFREE(z);
   ZNP_FASTFREE(limbs);
}


/*
   Multiplication/squaring using Kronecker substitution at 2^b, -2^b,
   2^(-b) and -2^(-b).
*/
void
zn_array_mul_KS4 (ulong* res,
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

   int sqr = (op1 == op2 && n1 == n2);

   // bits in each output coefficient
   unsigned bits = 2 * mod->bits + ceil_lg (n2);
   
   // we're evaluating at x = B, -B, 1/B, -1/B,
   // where B = 2^b, and b = ceil(bits / 4)
   unsigned b = (bits + 3) / 4;

   // number of ulongs required to store each base-B^2 digit
   unsigned w = CEIL_DIV (2 * b, ULONG_BITS);
   ZNP_ASSERT (w <= 2);

   // Write f1(x) = f1e(x^2) + x * f1o(x^2)
   //       f2(x) = f2e(x^2) + x * f2o(x^2)
   //        h(x) =  he(x^2) + x *  ho(x^2)
   // "e" = even, "o" = odd

   size_t n1o = n1 / 2;
   size_t n1e = n1 - n1o;

   size_t n2o = n2 / 2;
   size_t n2e = n2 - n2o;

   size_t n3 = n1 + n2 - 1;   // length of h
   size_t n3o = n3 / 2;
   size_t n3e = n3 - n3o;

   // Put k1 = number of limbs needed to store f1(B) and |f1(-B)|.
   // In f1(B), the leading coefficient starts at bit position b * (n1 - 1)
   // and has length 2b, and the coefficients overlap so we need an extra bit
   // for the carry: this gives (n1 + 1) * b + 1 bits. Ditto for f2.
   size_t k1 = CEIL_DIV ((n1 + 1) * b + 1, GMP_NUMB_BITS);
   size_t k2 = CEIL_DIV ((n2 + 1) * b + 1, GMP_NUMB_BITS);
   size_t k3 = k1 + k2;

   // allocate space
   ZNP_FASTALLOC (limbs, mp_limb_t, 6624, 5 * k3);
   mp_limb_t* v1_buf0 = limbs;           // k1 limbs
   mp_limb_t* v2_buf0 = v1_buf0 + k1;    // k2 limbs
   mp_limb_t* v1_buf1 = v2_buf0 + k2;    // k1 limbs
   mp_limb_t* v2_buf1 = v1_buf1 + k1;    // k2 limbs
   mp_limb_t* v1_buf2 = v2_buf1 + k2;    // k1 limbs
   mp_limb_t* v2_buf2 = v1_buf2 + k1;    // k2 limbs
   mp_limb_t* v1_buf3 = v2_buf2 + k2;    // k1 limbs
   mp_limb_t* v2_buf3 = v1_buf3 + k1;    // k2 limbs
   mp_limb_t* v1_buf4 = v2_buf3 + k2;    // k1 limbs
   mp_limb_t* v2_buf4 = v1_buf4 + k1;    // k2 limbs

   // arrange overlapping buffers to minimise memory use
   // "p" = plus, "m" = minus
   // "n" = normal order, "r" = reciprocal order
   mp_limb_t* v1en = v1_buf0;
   mp_limb_t* v1on = v1_buf1;
   mp_limb_t* v1pn = v1_buf2;
   mp_limb_t* v1mn = v1_buf0;
   mp_limb_t* v2en = v2_buf0;
   mp_limb_t* v2on = v2_buf1;
   mp_limb_t* v2pn = v2_buf2;
   mp_limb_t* v2mn = v2_buf0;
   mp_limb_t* v3pn = v1_buf1;
   mp_limb_t* v3mn = v1_buf2;
   mp_limb_t* v3en = v1_buf0;
   mp_limb_t* v3on = v1_buf1;

   mp_limb_t* v1er = v1_buf2;
   mp_limb_t* v1or = v1_buf3;
   mp_limb_t* v1pr = v1_buf4;
   mp_limb_t* v1mr = v1_buf2;
   mp_limb_t* v2er = v2_buf2;
   mp_limb_t* v2or = v2_buf3;
   mp_limb_t* v2pr = v2_buf4;
   mp_limb_t* v2mr = v2_buf2;
   mp_limb_t* v3pr = v1_buf3;
   mp_limb_t* v3mr = v1_buf4;
   mp_limb_t* v3er = v1_buf2;
   mp_limb_t* v3or = v1_buf3;
   
   ZNP_FASTALLOC (z, ulong, 6624, 2 * w * (n3e + 1));
   ulong* zn = z;
   ulong* zr = z + w * (n3e + 1);

   int v3m_neg;

   // -------------------------------------------------------------------------
   //     "normal" evaluation points
   
   if (!sqr)
   {
      // multiplication version

      // evaluate f1e(B^2) and B * f1o(B^2)
      // We need max(2 * b*n1e, 2 * b*n1o + b) bits for this packing step,
      // which is safe since (n1 + 1) * b + 1 >= max(2 * b*n1e, 2 * b*n1o + b).
      // Ditto for f2 below.
      zn_array_pack (v1en, op1, n1e, 2, 2 * b, 0, k1);
      zn_array_pack (v1on, op1 + 1, n1o, 2, 2 * b, b, k1);

      // compute  f1(B)  =  f1e(B^2) + B * f1o(B^2)
      //    and |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
      ZNP_ASSERT_NOCARRY (mpn_add_n (v1pn, v1en, v1on, k1));
      v3m_neg = signed_mpn_sub_n (v1mn, v1en, v1on, k1);

      // evaluate f2e(B^2) and B * f2o(B^2)
      zn_array_pack (v2en, op2, n2e, 2, 2 * b, 0, k2);
      zn_array_pack (v2on, op2 + 1, n2o, 2, 2 * b, b, k2);
      
      // compute  f2(B)  =  f2e(B^2) + B * f2o(B^2)
      //    and |f2(-B)| = |f2e(B^2) - B * f2o(B^2)|
      ZNP_ASSERT_NOCARRY (mpn_add_n (v2pn, v2en, v2on, k2));
      v3m_neg ^= signed_mpn_sub_n (v2mn, v2en, v2on, k2);

      // compute  h(B)  =  f1(B)   *  f2(B)
      //    and |h(-B)| = |f1(-B)| * |f2(-B)|
      // hn_neg is set if h(-B) is negative
      ZNP_mpn_mul (v3pn, v1pn, k1, v2pn, k2);
      ZNP_mpn_mul (v3mn, v1mn, k1, v2mn, k2);
   }
   else
   {
      // squaring version

      // evaluate f1e(B^2) and B * f1o(B^2)
      zn_array_pack (v1en, op1, n1e, 2, 2 * b, 0, k1);
      zn_array_pack (v1on, op1 + 1, n1o, 2, 2 * b, b, k1);

      // compute  f1(B)  =  f1e(B^2) + B * f1o(B^2)
      //    and |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
      ZNP_ASSERT_NOCARRY (mpn_add_n (v1pn, v1en, v1on, k1));
      signed_mpn_sub_n (v1mn, v1en, v1on, k1);

      // compute h(B) =  f1(B)^2
      //    and h(-B) = |f1(-B)|^2
      // hn_neg is cleared since h(-B) is never negative
      ZNP_mpn_mul (v3pn, v1pn, k1, v1pn, k1);
      ZNP_mpn_mul (v3mn, v1mn, k1, v1mn, k1);
      v3m_neg = 0;
   }

   // Each coefficient of h(B) is up to 4b bits long, so h(B) needs at most
   // ((n1 + n2 + 2) * b + 1) bits. (The extra +1 is to accommodate carries
   // generated by overlapping coefficients.)  The buffer has at least
   // ((n1 + n2 + 2) * b + 2) bits. Therefore we can safely store 2*h(B) etc.

   // compute     2 * he(B^2) = h(B) + h(-B)
   // and     B * 2 * ho(B^2) = h(B) - h(-B)
   if (v3m_neg)
   {
      ZNP_ASSERT_NOCARRY (mpn_sub_n (v3en, v3pn, v3mn, k3));
      ZNP_ASSERT_NOCARRY (mpn_add_n (v3on, v3pn, v3mn, k3));
   }
   else
   {
      ZNP_ASSERT_NOCARRY (mpn_add_n (v3en, v3pn, v3mn, k3));
      ZNP_ASSERT_NOCARRY (mpn_sub_n (v3on, v3pn, v3mn, k3));
   }

   // -------------------------------------------------------------------------
   //     "reciprocal" evaluation points

   // correction factors to take into account that if a polynomial has even
   // length, its even and odd coefficients are swapped when the polynomial
   // is reversed
   unsigned a1 = (n1 & 1) ? 0 : b;
   unsigned a2 = (n2 & 1) ? 0 : b;
   unsigned a3 = (n3 & 1) ? 0 : b;

   if (!sqr)
   {
      // multiplication version
   
      // evaluate B^(n1-1) * f1e(1/B^2) and B^(n1-2) * f1o(1/B^2)
      zn_array_pack (v1er, op1 + 2*(n1e - 1), n1e, -2, 2 * b, a1, k1);
      zn_array_pack (v1or, op1 + 1 + 2*(n1o - 1), n1o, -2, 2 * b, b - a1, k1);

      // compute  B^(n1-1) * f1(1/B) =
      //              B^(n1-1) * f1e(1/B^2) + B^(n1-2) * f1o(1/B^2)
      //    and  |B^(n1-1) * f1(-1/B)| =
      //             |B^(n1-1) * f1e(1/B^2) - B^(n1-2) * f1o(1/B^2)|
      ZNP_ASSERT_NOCARRY (mpn_add_n (v1pr, v1er, v1or, k1));
      v3m_neg = signed_mpn_sub_n (v1mr, v1er, v1or, k1);

      // evaluate B^(n2-1) * f2e(1/B^2) and B^(n2-2) * f2o(1/B^2)
      zn_array_pack (v2er, op2 + 2*(n2e - 1), n2e, -2, 2 * b, a2, k2);
      zn_array_pack (v2or, op2 + 1 + 2*(n2o - 1), n2o, -2, 2 * b, b - a2, k2);

      // compute  B^(n2-1) * f2(1/B) =
      //              B^(n2-1) * f2e(1/B^2) + B^(n2-2) * f2o(1/B^2)
      //    and  |B^(n1-1) * f2(-1/B)| =
      //             |B^(n2-1) * f2e(1/B^2) - B^(n2-2) * f2o(1/B^2)|
      ZNP_ASSERT_NOCARRY (mpn_add_n (v2pr, v2er, v2or, k2));
      v3m_neg ^= signed_mpn_sub_n (v2mr, v2er, v2or, k2);

      // compute B^(n3-1) * h(1/B) =
      //                 (B^(n1-1) * f1(1/B)) * (B^(n2-1) * f2(1/B))
      //     and |B^(n3-1) * h(-1/B)| =
      //                 |B^(n1-1) * f1(-1/B)| * |B^(n2-1) * f2(-1/B)|
      // hr_neg is set if h(-1/B) is negative
      ZNP_mpn_mul (v3pr, v1pr, k1, v2pr, k2);
      ZNP_mpn_mul (v3mr, v1mr, k1, v2mr, k2);
   }
   else
   {
      // squaring version

      // evaluate B^(n1-1) * f1e(1/B^2) and B^(n1-2) * f1o(1/B^2)
      zn_array_pack (v1er, op1 + 2*(n1e - 1), n1e, -2, 2 * b, a1, k1);
      zn_array_pack (v1or, op1 + 1 + 2*(n1o - 1), n1o, -2, 2 * b, b - a1, k1);

      // compute  B^(n1-1) * f1(1/B) =
      //              B^(n1-1) * f1e(1/B^2) + B^(n1-2) * f1o(1/B^2)
      //    and  |B^(n1-1) * f1(-1/B)| =
      //             |B^(n1-1) * f1e(1/B^2) - B^(n1-2) * f1o(1/B^2)|
      ZNP_ASSERT_NOCARRY (mpn_add_n (v1pr, v1er, v1or, k1));
      signed_mpn_sub_n (v1mr, v1er, v1or, k1);

      // compute B^(n3-1) * h(1/B)  = (B^(n1-1) * f1(1/B))^2
      //     and B^(n3-1) * h(-1/B) = |B^(n1-1) * f1(-1/B)|^2
      // hr_neg is cleared since h(-1/B) is never negative
      ZNP_mpn_mul (v3pr, v1pr, k1, v1pr, k1);
      ZNP_mpn_mul (v3mr, v1mr, k1, v1mr, k1);
      v3m_neg = 0;
   }

   // compute 2 * B^(n3-1) * he(1/B^2)
   //                = B^(n3-1) * h(1/B) + B^(n3-1) * h(-1/B)
   //    and  2 * B^(n3-2) * ho(1/B^2)
   //                = B^(n3-1) * h(1/B) - B^(n3-1) * h(-1/B)
   if (v3m_neg)
   {
      ZNP_ASSERT_NOCARRY (mpn_sub_n (v3er, v3pr, v3mr, k3));
      ZNP_ASSERT_NOCARRY (mpn_add_n (v3or, v3pr, v3mr, k3));
   }
   else
   {
      ZNP_ASSERT_NOCARRY (mpn_add_n (v3er, v3pr, v3mr, k3));
      ZNP_ASSERT_NOCARRY (mpn_sub_n (v3or, v3pr, v3mr, k3));
   }

   // -------------------------------------------------------------------------
   //     combine "normal" and "reciprocal" information

   // decompose he(B^2) and B^(2*(n3e-1)) * he(1/B^2) into base-B^2 digits
   zn_array_unpack_SAFE (zn, v3en, n3e + 1, 2 * b, 1, k3);
   zn_array_unpack_SAFE (zr, v3er, n3e + 1, 2 * b, a3 + 1, k3);
   
   // combine he(B^2) and he(1/B^2) information to get even coefficients of h
   zn_array_recover_reduce (res, 2, zn, zr, n3e, 2 * b, redc, mod);

   // decompose ho(B^2) and B^(2*(n3o-1)) * ho(1/B^2) into base-B^2 digits
   zn_array_unpack_SAFE (zn, v3on, n3o + 1, 2 * b, b + 1, k3);
   zn_array_unpack_SAFE (zr, v3or, n3o + 1, 2 * b, b - a3 + 1, k3);

   // combine ho(B^2) and ho(1/B^2) information to get odd coefficients of h
   zn_array_recover_reduce (res + 1, 2, zn, zr, n3o, 2 * b, redc, mod);
   
   ZNP_FASTFREE (z);
   ZNP_FASTFREE (limbs);
}


// end of file ****************************************************************
