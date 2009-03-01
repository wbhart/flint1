/*
   mul_ks.c:  multiplication by Kronecker substitution
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
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


#ifdef ZNP_USE_FLINT

// use FLINT integer multiplication
#include <FLINT/mpn_extras.h>
#define ZNP_mpn_mul F_mpn_mul

#else

// use GMP integer multiplication
#define ZNP_mpn_mul mpn_mul

#endif



/*
   Sets res[i*skip] = reduction modulo _mod_ of the i-th entry of _op_,
   for 0 <= i < len. Each entry of _op_ is _words_ ulongs.
   
   If the _redc_ flag is set, the results are divided by -B mod n
   (only allowed if the modulus is odd).
   
   Must have 1 <= words <= 3.
*/
#define array_reduce \
    ZNP_array_reduce
void array_reduce(ulong* res, ptrdiff_t skip, const ulong* op, size_t len,
                  unsigned words, int redc, const zn_mod_t mod)
{
   ZNP_ASSERT(words >= 1 && words <= 3);
   ZNP_ASSERT((mod->n & 1) || !redc);

   if (words == 1)
   {
      if (redc)
      {
         for (; len; len--, res += skip, op++)
            *res = zn_mod_reduce_redc(*op, mod);
      }
      else
      {
         for (; len; len--, res += skip, op++)
            *res = zn_mod_reduce(*op, mod);
      }
   }
   else if (words == 2)
   {
      if (redc)
      {
         for (; len; len--, res += skip, op += 2)
            *res = zn_mod_reduce2_redc(op[1], op[0], mod);
      }
      else
      {
         for (; len; len--, res += skip, op += 2)
            *res = zn_mod_reduce2(op[1], op[0], mod);
      }
   }
   else    // words == 3
   {
      if (redc)
      {
         for (; len; len--, res += skip, op += 3)
            *res = zn_mod_reduce3_redc(op[2], op[1], op[0], mod);
      }
      else
      {
         for (; len; len--, res += skip, op += 3)
            *res = zn_mod_reduce3(op[2], op[1], op[0], mod);
      }
   }
}



/*
   void zn_array_recip_fix_reduce(
              ulong* res, ptrdiff_t skip, const ulong* op1, const ulong* op2,
              size_t len, unsigned bits, const zn_mod_t mod)
   
   This is a helper function for the variants of KS that evaluate at
   "reciprocal" evaluation points like 2^(-N); it implements essentially the
   algorithm of section 3.2 of [Har07], plus reductions mod n.

   It accepts two integers X and Y written in base M = 2^bits,
   where 2 <= 2*bits <= 3*ULONG_BITS. It assumes that

      X = a[0] + a[1]*M + ... + a[len-1]*M^(len-1),
      Y = a[len-1] + a[len-2]*M + ... + a[0]*M^(len-1),

   where each a[i] is two "digits" long, and where the high digit of a[i]
   is at most M-2 (i.e. may not equal M-1). It reconstructs the a[i],
   reduces them mod n, and stores the results in an array.
   
   The input is supplied as follows. X is in op1, Y is in op2. They are both
   arrays of values that are _bits_ bits wide, where bits <= 3*ULONG_BITS/2.
   Each value takes up one ulong if bits <= ULONG_BITS, otherwise two ulongs.
   There are len + 1 such values in each array. (i.e. each array consists of
   (len + 1) * ceil(bits / ULONG_BITS) ulongs.)
   
   The output (_len_ ulongs) is written to the array _res_, with consecutive
   outputs separated by _skip_ ulongs.
   
   _mod_ describes the modulus n.

   There are five versions:
   
   zn_array_recip_fix_reduce1(): requires 0 < 2*bits <= ULONG_BITS
   zn_array_recip_fix_reduce2(): requires ULONG_BITS < 2*bits < 2*ULONG_BITS
   zn_array_recip_fix_reduce2b(): requires bits == ULONG_BITS
   zn_array_recip_fix_reduce3(): requires 2*ULONG_BITS < 2*bits <= 3*ULONG_BITS
   zn_array_recip_fix_reduce(): dispatches to one of the above depending
                                on _bits_

   If the _redc_ flag is set, the modular reductions are performed using REDC,
   i.e. the result contain an extra factor of -1/B mod n
   (where B = 2^ULONG_BITS).
*/

#define zn_array_recip_fix_reduce1 \
    ZNP_zn_array_recip_fix_reduce1
void zn_array_recip_fix_reduce1(ulong* res, ptrdiff_t skip, const ulong* op1,
                                const ulong* op2, size_t len, unsigned bits,
                                int redc, const zn_mod_t mod)
{
   ZNP_ASSERT(bits >= 1 && 2*bits <= ULONG_BITS);

   ulong mask = (1UL << bits) - 1;

   // (x0, x1) and (y0, y1) are two-digit windows into X and Y.
   ulong x1, x0 = *op1++;

   op2 += len;
   ulong y0, y1 = *op2--;

   ulong borrow = 0;

   if (redc)
   {
      // REDC version
      for (; len; len--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT(y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce_redc(x0 + (y1 << bits), mod);
         res += skip;
         ZNP_ASSERT(y1 != mask);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = (y0 - x0) & mask;
         x0 = x1 & mask;
      }
   }
   else
   {
      // plain reduction version
      for (; len; len--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT(y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce(x0 + (y1 << bits), mod);
         res += skip;
         ZNP_ASSERT(y1 != mask);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = (y0 - x0) & mask;
         x0 = x1 & mask;
      }
   }
}


#define zn_array_recip_fix_reduce2 \
    ZNP_zn_array_recip_fix_reduce2
void zn_array_recip_fix_reduce2(ulong* res, ptrdiff_t skip, const ulong* op1,
                                const ulong* op2, size_t len, unsigned bits,
                                int redc, const zn_mod_t mod)
{
   ZNP_ASSERT(2*bits > ULONG_BITS  &&  bits < ULONG_BITS);

   // The main loop is the same as in zn_array_recip_fix_reduce1(), but the
   // modular reduction step needs to handle two input words.

   ulong mask = (1UL << bits) - 1;

   ulong x1, x0 = *op1++;

   op2 += len;
   ulong y0, y1 = *op2--;

   ulong borrow = 0;
   
   unsigned bits2 = ULONG_BITS - bits;

   if (redc)
   {
      // REDC version
      for (; len; len--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT(y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2_redc(y1 >> bits2, x0 + (y1 << bits), mod);
         res += skip;
         ZNP_ASSERT(y1 != mask);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = (y0 - x0) & mask;
         x0 = x1 & mask;
      }
   }
   else
   {
      // plain reduction version
      for (; len; len--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT(y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2(y1 >> bits2, x0 + (y1 << bits), mod);
         res += skip;
         ZNP_ASSERT(y1 != mask);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = (y0 - x0) & mask;
         x0 = x1 & mask;
      }
   }
}


#define zn_array_recip_fix_reduce2b \
    ZNP_zn_array_recip_fix_reduce2b
void zn_array_recip_fix_reduce2b(ulong* res, ptrdiff_t skip, const ulong* op1,
                                 const ulong* op2, size_t len, unsigned bits,
                                 int redc, const zn_mod_t mod)
{
   ZNP_ASSERT(bits == ULONG_BITS);

   // Basically the same code as zn_array_recip_fix_reduce2(), specialised
   // for bits == ULONG_BITS.

   ulong x1, x0 = *op1++;

   op2 += len;
   ulong y0, y1 = *op2--;

   ulong borrow = 0;

   if (redc)
   {
      // REDC version
      for (; len; len--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT(y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2_redc(y1, x0, mod);
         res += skip;
         ZNP_ASSERT(y1 != -1UL);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = y0 - x0;
         x0 = x1;
      }
   }
   else
   {
      // plain reduction version
      for (; len; len--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT(y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2(y1, x0, mod);
         res += skip;
         ZNP_ASSERT(y1 != -1UL);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = y0 - x0;
         x0 = x1;
      }
   }
}


#define zn_array_recip_fix_reduce3 \
    ZNP_zn_array_recip_fix_reduce3
void zn_array_recip_fix_reduce3(ulong* res, ptrdiff_t skip, const ulong* op1,
                                const ulong* op2, size_t len, unsigned bits,
                                int redc, const zn_mod_t mod)
{
   ZNP_ASSERT(bits > ULONG_BITS  &&  2*bits <= 3*ULONG_BITS);

   // The main loop is the same as in zn_array_recip_fix_reduce1(), but needs
   // to operate on double-word quantities everywhere, i.e. we basically
   // simulate double-word registers. The suffices L and H stand for low
   // and high words of each.

   ulong maskL = -1UL;
   ulong maskH = (1UL << (bits - ULONG_BITS)) - 1;

   ulong x1L, x0L = *op1++;
   ulong x1H, x0H = *op1++;

   op2 += 2 * len + 1;
   ulong y0H, y1H = *op2--;
   ulong y0L, y1L = *op2--;

   ulong borrow = 0;

   unsigned bits1 = bits - ULONG_BITS;
   unsigned bits2 = 2*ULONG_BITS - bits;

   if (redc)
   {
      // REDC version
      for (; len; len--)
      {
         y0H = *op2--;
         y0L = *op2--;
         x1L = *op1++;
         x1H = *op1++;
         if ((y0H < x0H) || (y0H == x0H  &&  y0L < x0L))
         {
            ZNP_ASSERT(y1H != 0 || y1L != 0);
            y1H -= (y1L-- == 0);
         }

         *res = zn_mod_reduce3_redc((y1H << bits1) + (y1L >> bits2),
                                    (y1L << bits1) + x0H, x0L, mod);
         res += skip;

         ZNP_ASSERT(y1L != maskL || y1H != maskH);
         if (borrow)
            y1H += (++y1L == 0);
         borrow = ((x1H < y1H) || (x1H == y1H  &&  x1L < y1L));
         ZNP_SUB_WIDE(x1H, x1L, x1H, x1L, y1H, y1L);
         ZNP_SUB_WIDE(y1H, y1L, y0H, y0L, x0H, x0L);
         y1H &= maskH;
         x0L = x1L;
         x0H = x1H & maskH;
      }
   }
   else
   {
      // plain reduction version
      for (; len; len--)
      {
         y0H = *op2--;
         y0L = *op2--;
         x1L = *op1++;
         x1H = *op1++;
         if ((y0H < x0H) || (y0H == x0H  &&  y0L < x0L))
         {
            ZNP_ASSERT(y1H != 0 || y1L != 0);
            y1H -= (y1L-- == 0);
         }

         *res = zn_mod_reduce3((y1H << bits1) + (y1L >> bits2),
                               (y1L << bits1) + x0H, x0L, mod);
         res += skip;

         ZNP_ASSERT(y1L != maskL || y1H != maskH);
         if (borrow)
            y1H += (++y1L == 0);
         borrow = ((x1H < y1H) || (x1H == y1H  &&  x1L < y1L));
         ZNP_SUB_WIDE(x1H, x1L, x1H, x1L, y1H, y1L);
         ZNP_SUB_WIDE(y1H, y1L, y0H, y0L, x0H, x0L);
         y1H &= maskH;
         x0L = x1L;
         x0H = x1H & maskH;
      }
   }
}


void zn_array_recip_fix_reduce(ulong* res, ptrdiff_t skip, const ulong* op1,
                               const ulong* op2, size_t len, unsigned bits,
                               int redc, const zn_mod_t mod)
{
   ZNP_ASSERT(bits > 0  &&  2*bits <= 3*ULONG_BITS);

   if (2*bits <= ULONG_BITS)
      zn_array_recip_fix_reduce1(res, skip, op1, op2, len, bits, redc, mod);
   else if (bits < ULONG_BITS)
      zn_array_recip_fix_reduce2(res, skip, op1, op2, len, bits, redc, mod);
   else if (bits == ULONG_BITS)
      zn_array_recip_fix_reduce2b(res, skip, op1, op2, len, bits, redc, mod);
   else
      zn_array_recip_fix_reduce3(res, skip, op1, op2, len, bits, redc, mod);
}



/*
   res := abs(op1 - op2).
   Returns 1 if op1 - op2 is negative, else zero.
*/
#define signed_mpn_sub_n \
    ZNP_signed_mpn_sub_n
ZNP_INLINE
int signed_mpn_sub_n(mp_limb_t* res, const mp_limb_t* op1,
                     const mp_limb_t* op2, size_t len)
{
   if (mpn_cmp(op1, op2, len) >= 0)
   {
      mpn_sub_n(res, op1, op2, len);
      return 0;
   }
   else
   {
      mpn_sub_n(res, op2, op1, len);
      return 1;
   }
}



/*
   Multiplication (or squaring) using kronecker substitution at 2^N.
*/
void zn_array_mul_KS1(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   ZNP_ASSERT(ULONG_MAX >= len1);
   ZNP_ASSERT((mod->n & 1) || !redc);

   int squaring = (op1 == op2 && len1 == len2);
   
   // Let f(x), g(x) be input polys (op1 and op2),
   // and h(x) their product in Z[x].

   size_t h_len = len1 + len2 - 1;
   
   // bits in each output coefficient (i.e. if we were really in Z[x])
   unsigned bits = 2 * mod->bits + ceil_lg(len2);
   
   // number of ulongs required to store each output coefficient
   unsigned words = CEIL_DIV(bits, ULONG_BITS);
   ZNP_ASSERT(words <= 3);

   // number of limbs needed to store f(2^bits) and g(2^bits)
   size_t f_limbs = CEIL_DIV(len1 * bits, GMP_NUMB_BITS);
   size_t g_limbs = CEIL_DIV(len2 * bits, GMP_NUMB_BITS);

   // allocate space
   ZNP_FASTALLOC(limbs, mp_limb_t, 6624, 2 * (f_limbs + g_limbs));
   mp_limb_t* t = limbs;
   mp_limb_t* f_eval = t; t += f_limbs;
   mp_limb_t* g_eval = t; t += g_limbs;
   mp_limb_t* h_eval = t; // t += f_limbs + g_limbs;

   if (!squaring)
   {
      // multiplication version

      // evaluate f(2^bits) and g(2^bits)
      zn_array_pack(f_eval, op1, len1, 1, bits, 0, 0);
      zn_array_pack(g_eval, op2, len2, 1, bits, 0, 0);

      // compute h(2^bits) = f(2^bits) * g(2^bits)
      ZNP_mpn_mul(h_eval, f_eval, f_limbs, g_eval, g_limbs);
   }
   else
   {
      // squaring version

      // evaluate f(2^bits)
      zn_array_pack(f_eval, op1, len1, 1, bits, 0, 0);

      // compute h(2^bits) = f(2^bits)^2
      ZNP_mpn_mul(h_eval, f_eval, f_limbs, f_eval, f_limbs);
   }

   // unpack coefficients of h, and reduce mod n
   ZNP_FASTALLOC(unpack, ulong, 6624, h_len * words);
   zn_array_unpack(unpack, h_eval, h_len, bits, 0);
   array_reduce(res, 1, unpack, h_len, words, redc, mod);

   ZNP_FASTFREE(unpack);
   ZNP_FASTFREE(limbs);
}



/*
   Multiplication (or squaring) using kronecker substitution at 2^N and -2^N.
*/
void zn_array_mul_KS2(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   ZNP_ASSERT(ULONG_MAX >= len1);
   ZNP_ASSERT((mod->n & 1) || !redc);
   
   // Let f(x), g(x) be input polys (op1 and op2),
   // and h(x) their product in Z[x].

   if (len2 == 1)
   {
      // code requires len2 > 1, so fall back on scalar multiplication routine
      _zn_array_scalar_mul(res, op1, len1, op2[0], redc, mod);
      return;
   }

   int squaring = (op1 == op2 && len1 == len2);

   // bits in each output coefficient (i.e. if we were really in Z[x])
   unsigned bits = 2 * mod->bits + ceil_lg(len2);
   
   // we're evaluating at x = B and -B,
   // where B = 2^N, and N = ceil(bits / 2)
   unsigned N = (bits + 1) / 2;

   // number of ulongs required to store each output coefficient
   unsigned words = CEIL_DIV(bits, ULONG_BITS);
   
   // Write f(x) = f0(x^2) + x * f1(x^2)
   //   and g(x) = g0(x^2) + x * g1(x^2)

   size_t f1_len = len1 / 2;
   size_t f0_len = len1 - f1_len;

   size_t g1_len = len2 / 2;
   size_t g0_len = len2 - g1_len;

   size_t h_len = len1 + len2 - 1;
   size_t h1_len = h_len / 2;
   size_t h0_len = h_len - h1_len;

   // Put f_limbs = number of limbs needed to store f(B) and |f(-B)|.
   // Actually it's slightly more, since when computing f0(B^2) and
   // B*f1(B^2) the bitpacking routine needs room for the last chunk of 2N
   // bits. So we need max(2*N*f0_len, 2*N*f1_len + N) = N*(len1 + 1) bits.
   size_t f_limbs = CEIL_DIV((len1 + 1)*N, GMP_NUMB_BITS);
   // ditto for g
   size_t g_limbs = CEIL_DIV((len2 + 1)*N, GMP_NUMB_BITS);
   size_t h_limbs = f_limbs + g_limbs;

   // allocate space
   // "p" = plus, "m" = minus
   ZNP_FASTALLOC(limbs, mp_limb_t, 6624, 3 * h_limbs);
   mp_limb_t* t = limbs;
   mp_limb_t* f_buf0 = t;  t += f_limbs;
   mp_limb_t* g_buf0 = t;  t += g_limbs;
   mp_limb_t* f_buf1 = t;  t += f_limbs;
   mp_limb_t* g_buf1 = t;  t += g_limbs;
   mp_limb_t* f_buf2 = t;  t += f_limbs;
   mp_limb_t* g_buf2 = t;  // t += g_limbs;

   // arrange overlapping buffers to minimise memory use
   mp_limb_t* f0_eval = f_buf0;
   mp_limb_t* g0_eval = g_buf0;
   mp_limb_t* f1_eval = f_buf1;
   mp_limb_t* g1_eval = g_buf1;
   mp_limb_t* fp_eval = f_buf2;
   mp_limb_t* gp_eval = g_buf2;
   mp_limb_t* fm_eval = f_buf0;
   mp_limb_t* gm_eval = g_buf0;
   mp_limb_t* hm_eval = f_buf1;
   mp_limb_t* hp_eval = f_buf0;
   mp_limb_t* h0_eval = f_buf2;
   mp_limb_t* h1_eval = f_buf2;
   
   ZNP_FASTALLOC(unpack, ulong, 6624, words * h0_len);

   int h_neg;
   mp_limb_t carry;

   if (!squaring)
   {
      // multiplication version

      // evaluate f0(B^2) and B * f1(B^2)
      zn_array_pack(f0_eval, op1, f0_len, 2, 2*N, 0, f_limbs);
      zn_array_pack(f1_eval, op1 + 1, f1_len, 2, 2*N, N, f_limbs);

      // evaluate g0(B^2) and B * g1(B^2)
      zn_array_pack(g0_eval, op2, g0_len, 2, 2*N, 0, g_limbs);
      zn_array_pack(g1_eval, op2 + 1, g1_len, 2, 2*N, N, g_limbs);

      // compute f(B) = f0(B^2) + B * f1(B^2)
      //     and g(B) = g0(B^2) + B * g1(B^2)
      carry = mpn_add_n(fp_eval, f0_eval, f1_eval, f_limbs);
      ZNP_ASSERT(carry == 0);
      carry = mpn_add_n(gp_eval, g0_eval, g1_eval, g_limbs);
      ZNP_ASSERT(carry == 0);

      // compute f(-B) = f0(B^2) - B * f1(B^2)
      //     and g(-B) = g0(B^2) - B * g1(B^2)
      h_neg  = signed_mpn_sub_n(fm_eval, f0_eval, f1_eval, f_limbs);
      h_neg ^= signed_mpn_sub_n(gm_eval, g0_eval, g1_eval, g_limbs);

      // compute h(B) = f(B) * g(B)
      // compute h(-B) = f(-B) * g(-B)
      // h_neg is set if h(-B) is negative
      ZNP_mpn_mul(hm_eval, fm_eval, f_limbs, gm_eval, g_limbs);
      ZNP_mpn_mul(hp_eval, fp_eval, f_limbs, gp_eval, g_limbs);
   }
   else
   {
      // squaring version

      // evaluate f0(B^2) and B * f1(B^2)
      zn_array_pack(f0_eval, op1, f0_len, 2, 2*N, 0, f_limbs);
      zn_array_pack(f1_eval, op1 + 1, f1_len, 2, 2*N, N, f_limbs);

      // compute f(B) = f0(B^2) + B * f1(B^2)
      carry = mpn_add_n(fp_eval, f0_eval, f1_eval, f_limbs);
      ZNP_ASSERT(carry == 0);

      // compute f(-B) = f0(B^2) - B * f1(B^2)
      signed_mpn_sub_n(fm_eval, f0_eval, f1_eval, f_limbs);

      // compute h(B) = f(B)^2
      // compute h(-B) = f(-B)^2
      // h_neg is cleared (since f(-B)^2 is never negative)
      ZNP_mpn_mul(hm_eval, fm_eval, f_limbs, fm_eval, f_limbs);
      ZNP_mpn_mul(hp_eval, fp_eval, f_limbs, fp_eval, f_limbs);
      h_neg = 0;
   }

   // compute 2 * h0(B^2) = h(B) + h(-B)
   // Note that the buffer is at least (len1 + len2 + 2)*N bits long, while
   // 2*h(B) needs at most ((len1 + len2 - 2)*N + bits + 2) bits. Thus we
   // have at least (4N - bits - 2) >= (bits - 2) >= 0 spare bits.
   if (h_neg)
      carry = mpn_sub_n(h0_eval, hp_eval, hm_eval, h_limbs);
   else
      carry = mpn_add_n(h0_eval, hp_eval, hm_eval, h_limbs);
   ZNP_ASSERT(carry == 0);

   // unpack coefficients of h0, and reduce mod n
   zn_array_unpack(unpack, h0_eval, h0_len, 2*N, 1);
   array_reduce(res, 2, unpack, h0_len, words, redc, mod);
   
   // compute 2*B * h1(B^2) = h(B) - h(-B)
   if (h_neg)
      carry = mpn_add_n(h0_eval, hp_eval, hm_eval, h_limbs);
   else
      carry = mpn_sub_n(h0_eval, hp_eval, hm_eval, h_limbs);
   ZNP_ASSERT(carry == 0);
   
   // unpack coefficients of h1, and reduce mod n
   zn_array_unpack(unpack, h1_eval, h1_len, 2*N, N + 1);
   array_reduce(res + 1, 2, unpack, h1_len, words, redc, mod);

   ZNP_FASTFREE(unpack);
   ZNP_FASTFREE(limbs);
}



/*
   Multiplication (or squaring) using kronecker substitution at 2^N and 2^(-N).
   
   Note: this routine does not appear to be competitive in practice with the
   other KS routines. It's here just for fun.
*/
void zn_array_mul_KS3(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   ZNP_ASSERT(ULONG_MAX >= len1);
   ZNP_ASSERT((mod->n & 1) || !redc);

   int squaring = (op1 == op2 && len1 == len2);

   // Let f(x), g(x) be input polys (op1 and op2),
   // and h(x) their product in Z[x].

   size_t h_len = len1 + len2 - 1;
   
   // bits in each output coefficient (i.e. if we were really in Z[x])
   unsigned bits = 2 * mod->bits + ceil_lg(len2);
   
   // we're evaluating at x = B and 1/B,
   // where B = 2^N, and N = ceil(bits / 2)
   unsigned N = (bits + 1) / 2;

   // number of ulongs required to store each base-B digit
   unsigned words = CEIL_DIV(N, ULONG_BITS);
   ZNP_ASSERT(words <= 2);
   
   // number of limbs needed to store f(B) and B^(len1-1) * f(1/B)
   size_t f_limbs = CEIL_DIV(len1 * N, GMP_NUMB_BITS);
   // ditto for g
   size_t g_limbs = CEIL_DIV(len2 * N, GMP_NUMB_BITS);
   
   // allocate space
   ZNP_FASTALLOC(limbs, mp_limb_t, 6624, 2 * (f_limbs + g_limbs));
   mp_limb_t* t = limbs;
   mp_limb_t* f_eval = t;  t += f_limbs;
   mp_limb_t* g_eval = t;  t += g_limbs;
   mp_limb_t* h_eval = t;  //  t1 += f_limbs + g_limbs;

   ZNP_FASTALLOC(digits, ulong, 6624, 2 * words * (h_len + 1));
   // "n" = normal order, "r" = reciprocal order
   ulong* hn_digits = digits;
   ulong* hr_digits = digits + words * (h_len + 1);

   if (!squaring)
   {
      // evaluate f(B) and g(B)
      zn_array_pack(f_eval, op1, len1, 1, N, 0, f_limbs);
      zn_array_pack(g_eval, op2, len2, 1, N, 0, g_limbs);

      // compute h(B) = f(B) * g(B)
      ZNP_mpn_mul(h_eval, f_eval, f_limbs, g_eval, g_limbs);
   }
   else
   {
      // evaluate f(B)
      zn_array_pack(f_eval, op1, len1, 1, N, 0, f_limbs);

      // compute h(B) = f(B)^2
      ZNP_mpn_mul(h_eval, f_eval, f_limbs, f_eval, f_limbs);
   }

   // decompose h(B) into base-B digits
   zn_array_unpack(hn_digits, h_eval, h_len + 1, N, 0);

   if (!squaring)
   {
      // multiplication version

      // evaluate B^(len1-1) * f(1/B) and B^(len2-1) * g(1/B)
      zn_array_pack(f_eval, op1 + len1 - 1, len1, -1, N, 0, f_limbs);
      zn_array_pack(g_eval, op2 + len2 - 1, len2, -1, N, 0, g_limbs);

      // compute B^(len1+len2-1) * h(1/B) =
      //                     (B^(len1-1) * f(1/B)) * (B^(len1-2) * g(1/B))
      ZNP_mpn_mul(h_eval, f_eval, f_limbs, g_eval, g_limbs);
   }
   else
   {
      // squaring version

      // evaluate B^(len1-1) * f(1/B)
      zn_array_pack(f_eval, op1 + len1 - 1, len1, -1, N, 0, f_limbs);

      // compute B^(len1+len2-1) * h(1/B) = (B^(len1-1) * f(1/B))^2
      ZNP_mpn_mul(h_eval, f_eval, f_limbs, f_eval, f_limbs);
   }

   // decompose h(1/B) into base-B digits
   zn_array_unpack(hr_digits, h_eval, h_len + 1, N, 0);

   // combine h(B) and h(1/B) information
   // (note: need to check that the high digit of each output coefficient
   // is < B - 1; this follows from an estimate in section 3.2 of [Har07].
   zn_array_recip_fix_reduce(res, 1, hn_digits, hr_digits, h_len, N,
                             redc, mod);
   
   ZNP_FASTFREE(digits);
   ZNP_FASTFREE(limbs);
}


/*
   Multiplication (or squaring) using kronecker substitution at 2^N, -2^N,
   2^(-N) and -2^(-N).
*/
void zn_array_mul_KS4(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   ZNP_ASSERT(ULONG_MAX >= len1);
   ZNP_ASSERT((mod->n & 1) || !redc);

   // Let f(x), g(x) be input polys (op1 and op2),
   // and h(x) their product in Z[x].

   if (len2 == 1)
   {
      // code requires len2 > 1, so fall back on scalar multiplication routine
      _zn_array_scalar_mul(res, op1, len1, op2[0], redc, mod);
      return;
   }

   int squaring = (op1 == op2 && len1 == len2);

   // bits in each output coefficient (i.e. if we were really in Z[x])
   unsigned bits = 2 * mod->bits + ceil_lg(len2);
   
   // we're evaluating at x = B, -B, 1/B, -1/B,
   // where B = 2^N, and N = ceil(bits / 4)
   unsigned N = (bits + 3) / 4;

   // Write f(x) = f0(x^2) + x * f1(x^2)
   //   and g(x) = g0(x^2) + x * g1(x^2)

   size_t f_len = len1;
   size_t f1_len = f_len / 2;
   size_t f0_len = f_len - f1_len;

   size_t g_len = len2;
   size_t g1_len = g_len / 2;
   size_t g0_len = g_len - g1_len;

   size_t h_len = f_len + g_len - 1;
   size_t h1_len = h_len / 2;
   size_t h0_len = h_len - h1_len;

   // number of ulongs required to store each base-B^2 digit
   unsigned words = CEIL_DIV(2*N, ULONG_BITS);
   ZNP_ASSERT(words <= 2);

   // Put f_limbs = number of limbs needed to store f(B) and |f(-B)|.
   // In f(B), the leading coefficient starts at bit position N*(f_len-1) and
   // has length 2N, and the coefficients overlap so we need on extra bit for
   // the carry: this gives (f_len + 1) * N + 1.
   size_t f_limbs = CEIL_DIV((f_len + 1) * N + 1, GMP_NUMB_BITS);
   // ditto for g
   size_t g_limbs = CEIL_DIV((g_len + 1) * N + 1, GMP_NUMB_BITS);
   size_t h_limbs = f_limbs + g_limbs;

   // allocate space
   // "p" = plus, "m" = minus
   // "n" = normal order, "r" = reciprocal order
   ZNP_FASTALLOC(limbs, mp_limb_t, 6624, 5 * h_limbs);
   mp_limb_t* t = limbs;
   mp_limb_t* f_buf0 = t;  t += f_limbs;
   mp_limb_t* g_buf0 = t;  t += g_limbs;
   mp_limb_t* f_buf1 = t;  t += f_limbs;
   mp_limb_t* g_buf1 = t;  t += g_limbs;
   mp_limb_t* f_buf2 = t;  t += f_limbs;
   mp_limb_t* g_buf2 = t;  t += g_limbs;
   mp_limb_t* f_buf3 = t;  t += f_limbs;
   mp_limb_t* g_buf3 = t;  t += g_limbs;
   mp_limb_t* f_buf4 = t;  t += f_limbs;
   mp_limb_t* g_buf4 = t;  // t += g_limbs;

   // arrange overlapping buffers to minimise memory use
   mp_limb_t* f0n = f_buf0;
   mp_limb_t* f1n = f_buf1;
   mp_limb_t* fpn = f_buf2;
   mp_limb_t* fmn = f_buf0;
   mp_limb_t* g0n = g_buf0;
   mp_limb_t* g1n = g_buf1;
   mp_limb_t* gpn = g_buf2;
   mp_limb_t* gmn = g_buf0;
   mp_limb_t* hpn = f_buf1;
   mp_limb_t* hmn = f_buf2;
   mp_limb_t* h0n = f_buf0;
   mp_limb_t* h1n = f_buf1;

   mp_limb_t* f0r = f_buf2;
   mp_limb_t* f1r = f_buf3;
   mp_limb_t* fpr = f_buf4;
   mp_limb_t* fmr = f_buf2;
   mp_limb_t* g0r = g_buf2;
   mp_limb_t* g1r = g_buf3;
   mp_limb_t* gpr = g_buf4;
   mp_limb_t* gmr = g_buf2;
   mp_limb_t* hpr = f_buf3;
   mp_limb_t* hmr = f_buf4;
   mp_limb_t* h0r = f_buf2;
   mp_limb_t* h1r = f_buf3;
   
   ZNP_FASTALLOC(digits, ulong, 6624, words * 2 * (h0_len + 1));
   ulong* hn_digits = digits;
   ulong* hr_digits = digits + words * (h0_len + 1);

   mp_limb_t carry;
   int hmn_neg, hmr_neg;

   // -------------------------------------------------------------------------
   //     "normal" evaluation points
   
   if (!squaring)
   {
      // multiplication version

      // evaluate f0(B^2) and B * f1(B^2)
      // Note: we need max(2*N*f0_len, 2*N*f1_len + N) bits for this packing
      // step, which is safely covered.
      zn_array_pack(f0n, op1, f0_len, 2, 2*N, 0, f_limbs);
      zn_array_pack(f1n, op1 + 1, f1_len, 2, 2*N, N, f_limbs);

      // compute f(B) = f0(B^2) + B * f1(B^2)
      //    and f(-B) = f0(B^2) - B * f1(B^2)
      carry = mpn_add_n(fpn, f0n, f1n, f_limbs);
      ZNP_ASSERT(carry == 0);
      hmn_neg = signed_mpn_sub_n(fmn, f0n, f1n, f_limbs);

      // evaluate g0(B^2) and B * g1(B^2)
      zn_array_pack(g0n, op2, g0_len, 2, 2*N, 0, g_limbs);
      zn_array_pack(g1n, op2 + 1, g1_len, 2, 2*N, N, g_limbs);
      
      // compute g(B) = g0(B^2) + B * g1(B^2)
      //    and g(-B) = g0(B^2) - B * g1(B^2)
      carry = mpn_add_n(gpn, g0n, g1n, g_limbs);
      ZNP_ASSERT(carry == 0);
      hmn_neg ^= signed_mpn_sub_n(gmn, g0n, g1n, g_limbs);

      // compute h(B) = f(B)  * g(B)
      //    and h(-B) = f(-B) * g(-B)
      // hmn_neg is set if h(-B) is negative
      ZNP_mpn_mul(hpn, fpn, f_limbs, gpn, g_limbs);
      ZNP_mpn_mul(hmn, fmn, f_limbs, gmn, g_limbs);
   }
   else
   {
      // squaring version

      // evaluate f0(B^2) and B * f1(B^2)
      zn_array_pack(f0n, op1, f0_len, 2, 2*N, 0, f_limbs);
      zn_array_pack(f1n, op1 + 1, f1_len, 2, 2*N, N, f_limbs);

      // compute f(B) = f0(B^2) + B * f1(B^2)
      //    and f(-B) = f0(B^2) - B * f1(B^2)
      carry = mpn_add_n(fpn, f0n, f1n, f_limbs);
      ZNP_ASSERT(carry == 0);
      signed_mpn_sub_n(fmn, f0n, f1n, f_limbs);

      // compute h(B) = f(B)^2
      //    and h(-B) = f(-B)^2
      // hmn_neg is cleared since h(-B) is never negative
      ZNP_mpn_mul(hpn, fpn, f_limbs, fpn, f_limbs);
      ZNP_mpn_mul(hmn, fmn, f_limbs, fmn, f_limbs);
      hmn_neg = 0;
   }

   // compute     2 * h0(B^2) = h(B) + h(-B)
   // and     B * 2 * h1(B^2) = h(B) - h(-B)
   // Note that the buffer is at least (f_len + g_len + 2) * N + 2 bits long,
   // while 2*h(B) needs at most (f_len + g_len - 2) * N + bits + 2 bits.
   // Thus we have at least 4N - bits >= 0 spare bits.
   if (hmn_neg)
   {
      carry = mpn_sub_n(h0n, hpn, hmn, h_limbs);
      ZNP_ASSERT(carry == 0);
      carry = mpn_add_n(h1n, hpn, hmn, h_limbs);
      ZNP_ASSERT(carry == 0);
   }
   else
   {
      carry = mpn_add_n(h0n, hpn, hmn, h_limbs);
      ZNP_ASSERT(carry == 0);
      carry = mpn_sub_n(h1n, hpn, hmn, h_limbs);
      ZNP_ASSERT(carry == 0);
   }

   // -------------------------------------------------------------------------
   //     "reciprocal" evaluation points

   if (!squaring)
   {
      // multiplication version
   
      // evaluate B^(len1-1) * f0(1/B^2)
      //      and B^(len1-2) * f1(1/B^2)
      zn_array_pack(f0r, op1 + 2*(f0_len - 1), f0_len, -2, 2*N,
                    (len1 & 1) ? 0 : N, f_limbs);
      zn_array_pack(f1r, op1 + 1 + 2*(f1_len - 1), f1_len, -2, 2*N,
                    (len1 & 1) ? N : 0, f_limbs);

      // compute B^(len1-1) * f(1/B) =
      //             B^(len1-1) * f0(1/B^2) + B^(len1-2) * f1(1/B^2)
      //    and  B^(len1-1) * f(-1/B) =
      //             B^(len1-1) * f0(1/B^2) - B^(len1-2) * f1(1/B^2)
      carry = mpn_add_n(fpr, f0r, f1r, f_limbs);
      ZNP_ASSERT(carry == 0);
      hmr_neg = signed_mpn_sub_n(fmr, f0r, f1r, f_limbs);

      // evaluate B^(len2-1) * g0(1/B^2)
      //      and B^(len2-2) * g1(1/B^2)
      zn_array_pack(g0r, op2 + 2*(g0_len - 1), g0_len, -2, 2*N,
                    (len2 & 1) ? 0 : N, g_limbs);
      zn_array_pack(g1r, op2 + 1 + 2*(g1_len - 1), g1_len, -2, 2*N,
                    (len2 & 1) ? N : 0, g_limbs);

      // compute B^(len2-1) * g(1/B) =
      //             B^(len2-1) * g0(1/B^2) + B^(len2-2) * g1(1/B^2)
      //    and  B^(len1-1) * g(-1/B) =
      //             B^(len2-1) * g0(1/B^2) - B^(len2-2) * g1(1/B^2)
      carry = mpn_add_n(gpr, g0r, g1r, g_limbs);
      ZNP_ASSERT(carry == 0);
      hmr_neg ^= signed_mpn_sub_n(gmr, g0r, g1r, g_limbs);

      // compute B^(h_len-1) * h(1/B) =
      //                 (B^(len1-1) * f(1/B)) * (B^(len2-1) * g(1/B))
      //     and B^(h_len-1) * h(-1/B) =
      //                 (B^(len1-1) * f(-1/B)) * (B^(len2-1) * g(-1/B))
      // hmr_neg is set if h(-1/B) is negative
      ZNP_mpn_mul(hpr, fpr, f_limbs, gpr, g_limbs);
      ZNP_mpn_mul(hmr, fmr, f_limbs, gmr, g_limbs);
   }
   else
   {
      // squaring version

      // evaluate B^(len1-1) * f0(1/B^2)
      //      and B^(len1-2) * f1(1/B^2)
      zn_array_pack(f0r, op1 + 2*(f0_len - 1), f0_len, -2, 2*N,
                    (len1 & 1) ? 0 : N, f_limbs);
      zn_array_pack(f1r, op1 + 1 + 2*(f1_len - 1), f1_len, -2, 2*N,
                    (len1 & 1) ? N : 0, f_limbs);

      // compute B^(len1-1) * f(1/B) =
      //             B^(len1-1) * f0(1/B^2) + B^(len1-2) * f1(1/B^2)
      //    and  B^(len1-1) * f(-1/B) =
      //             B^(len1-1) * f0(1/B^2) - B^(len1-2) * f1(1/B^2)
      carry = mpn_add_n(fpr, f0r, f1r, f_limbs);
      ZNP_ASSERT(carry == 0);
      signed_mpn_sub_n(fmr, f0r, f1r, f_limbs);

      // compute B^(h_len-1) * h(1/B) = (B^(len1-1) * f(1/B))^2
      //     and B^(h_len-1) * h(-1/B) = (B^(len1-1) * f(-1/B))^2
      // hmr_neg is cleared since h(-1/B) is never negative
      ZNP_mpn_mul(hpr, fpr, f_limbs, fpr, f_limbs);
      ZNP_mpn_mul(hmr, fmr, f_limbs, fmr, f_limbs);
      hmr_neg = 0;
   }

   // compute 2 * B^(h_len-1) * h0(1/B^2)
   //                = B^(h_len-1) * h(1/B) + B^(h_len-1) * h(-1/B)
   //    and  2 * B^(h_len-2) * h1(1/B^2)
   //                = B^(h_len-1) * h(1/B) - B^(h_len-1) * h(-1/B)
   if (hmr_neg)
   {
      carry = mpn_sub_n(h0r, hpr, hmr, h_limbs);
      ZNP_ASSERT(carry == 0);
      carry = mpn_add_n(h1r, hpr, hmr, h_limbs);
      ZNP_ASSERT(carry == 0);
   }
   else
   {
      carry = mpn_add_n(h0r, hpr, hmr, h_limbs);
      ZNP_ASSERT(carry == 0);
      carry = mpn_sub_n(h1r, hpr, hmr, h_limbs);
      ZNP_ASSERT(carry == 0);
   }

   // -------------------------------------------------------------------------
   //     combine "normal" and "reciprocal" information

   // decompose h0(B^2) and B^(2*(h0_len-1)) * h0(1/B^2) into base-B^2 digits
   zn_array_unpack(hn_digits, h0n, h0_len + 1, 2*N, 1);
   zn_array_unpack(hr_digits, h0r, h0_len + 1, 2*N, (h_len & 1) ? 1 : (N + 1));
   
   // combine h0(B^2) and h0(1/B^2) information to get even coefficients of h
   zn_array_recip_fix_reduce(res, 2, hn_digits, hr_digits, h0_len, 2*N,
                             redc, mod);

   // decompose h1(B^2) and B^(2*(h1_len-1)) * h1(1/B^2) into base-B^2 digits
   zn_array_unpack(hn_digits, h1n, h1_len + 1, 2*N, N + 1);
   zn_array_unpack(hr_digits, h1r, h1_len + 1, 2*N, (h_len & 1) ? (N + 1) : 1);

   // combine h1(B^2) and h1(1/B^2) information to get odd coefficients of h
   zn_array_recip_fix_reduce(res + 1, 2, hn_digits, hr_digits,
                             h1_len, 2*N, redc, mod);
   
   ZNP_FASTFREE(digits);
   ZNP_FASTFREE(limbs);
}



// end of file ****************************************************************
