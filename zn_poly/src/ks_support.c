/*
   ks_support.c:  support routines for algorithms based on Kronecker
                  substitution
   
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


void
array_reduce (ulong* res, ptrdiff_t s, const ulong* op, size_t n, unsigned w,
              int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (w >= 1 && w <= 3);
   ZNP_ASSERT ((mod->m & 1) || !redc);

   if (w == 1)
   {
      if (redc)
      {
         for (; n; n--, res += s, op++)
            *res = zn_mod_reduce_redc (*op, mod);
      }
      else
      {
         for (; n; n--, res += s, op++)
            *res = zn_mod_reduce (*op, mod);
      }
   }
   else if (w == 2)
   {
      if (redc)
      {
         for (; n; n--, res += s, op += 2)
            *res = zn_mod_reduce2_redc (op[1], op[0], mod);
      }
      else
      {
         for (; n; n--, res += s, op += 2)
            *res = zn_mod_reduce2 (op[1], op[0], mod);
      }
   }
   else    // w == 3
   {
      if (redc)
      {
         for (; n; n--, res += s, op += 3)
            *res = zn_mod_reduce3_redc (op[2], op[1], op[0], mod);
      }
      else
      {
         for (; n; n--, res += s, op += 3)
            *res = zn_mod_reduce3 (op[2], op[1], op[0], mod);
      }
   }
}



/*
   Same as zn_array_recover_reduce(), but requires 0 < 2 * b <= ULONG_BITS
*/
#define zn_array_recover_reduce1 \
    ZNP_zn_array_recover_reduce1
void
zn_array_recover_reduce1 (ulong* res, ptrdiff_t s, const ulong* op1,
                          const ulong* op2, size_t n, unsigned b,
                          int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (b >= 1 && 2 * b <= ULONG_BITS);

   ulong mask = (1UL << b) - 1;

   // (x0, x1) and (y0, y1) are two-digit windows into X and Y.
   ulong x1, x0 = *op1++;

   op2 += n;
   ulong y0, y1 = *op2--;

   ulong borrow = 0;

   if (redc)
   {
      // REDC version
      for (; n; n--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT (y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce_redc (x0 + (y1 << b), mod);
         res += s;
         ZNP_ASSERT (y1 != mask);
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
      for (; n; n--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT (y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce (x0 + (y1 << b), mod);
         res += s;
         ZNP_ASSERT (y1 != mask);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = (y0 - x0) & mask;
         x0 = x1 & mask;
      }
   }
}


/*
   Same as zn_array_recover_reduce(), but requires
   ULONG_BITS < 2 * b < 2*ULONG_BITS
*/
#define zn_array_recover_reduce2 \
    ZNP_zn_array_recover_reduce2
void
zn_array_recover_reduce2 (ulong* res, ptrdiff_t s, const ulong* op1,
                          const ulong* op2, size_t n, unsigned b,
                          int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (2 * b > ULONG_BITS  &&  b < ULONG_BITS);

   // The main loop is the same as in zn_array_recover_reduce1(), but the
   // modular reduction step needs to handle two input words.

   ulong mask = (1UL << b) - 1;

   ulong x1, x0 = *op1++;

   op2 += n;
   ulong y0, y1 = *op2--;

   ulong borrow = 0;
   
   unsigned b2 = ULONG_BITS - b;

   if (redc)
   {
      // REDC version
      for (; n; n--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT (y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2_redc (y1 >> b2, x0 + (y1 << b), mod);
         res += s;
         ZNP_ASSERT (y1 != mask);
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
      for (; n; n--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT (y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2 (y1 >> b2, x0 + (y1 << b), mod);
         res += s;
         ZNP_ASSERT (y1 != mask);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = (y0 - x0) & mask;
         x0 = x1 & mask;
      }
   }
}


/*
   Same as zn_array_recover_reduce(), but requires b == ULONG_BITS
*/
#define zn_array_recover_reduce2b \
    ZNP_zn_array_recover_reduce2b
void
zn_array_recover_reduce2b (ulong* res, ptrdiff_t s, const ulong* op1,
                           const ulong* op2, size_t n, unsigned b,
                           int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (b == ULONG_BITS);

   // Basically the same code as zn_array_recover_reduce2(), specialised
   // for b == ULONG_BITS.

   ulong x1, x0 = *op1++;

   op2 += n;
   ulong y0, y1 = *op2--;

   ulong borrow = 0;

   if (redc)
   {
      // REDC version
      for (; n; n--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT (y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2_redc (y1, x0, mod);
         res += s;
         ZNP_ASSERT (y1 != -1UL);
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
      for (; n; n--)
      {
         y0 = *op2--;
         x1 = *op1++;
         if (y0 < x0)
         {
            ZNP_ASSERT (y1 != 0);
            y1--;
         }
         *res = zn_mod_reduce2 (y1, x0, mod);
         res += s;
         ZNP_ASSERT (y1 != -1UL);
         y1 += borrow;
         borrow = (x1 < y1);
         x1 -= y1;
         y1 = y0 - x0;
         x0 = x1;
      }
   }
}


/*
   Same as zn_array_recover_reduce(), but requires
   2 * ULONG_BITS < 2 * b <= 3 * ULONG_BITS.
*/
#define zn_array_recover_reduce3 \
    ZNP_zn_array_recover_reduce3
void
zn_array_recover_reduce3 (ulong* res, ptrdiff_t s, const ulong* op1,
                          const ulong* op2, size_t n, unsigned b,
                          int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (b > ULONG_BITS  &&  2 * b <= 3 * ULONG_BITS);

   // The main loop is the same as in zn_array_recover_reduce1(), but needs
   // to operate on double-word quantities everywhere, i.e. we simulate
   // double-word registers. The suffixes L and H stand for low and high words
   // of each.

   ulong maskL = -1UL;
   ulong maskH = (1UL << (b - ULONG_BITS)) - 1;

   ulong x1L, x0L = *op1++;
   ulong x1H, x0H = *op1++;

   op2 += 2 * n + 1;
   ulong y0H, y1H = *op2--;
   ulong y0L, y1L = *op2--;

   ulong borrow = 0;

   unsigned b1 = b - ULONG_BITS;
   unsigned b2 = 2 * ULONG_BITS - b;

   if (redc)
   {
      // REDC version
      for (; n; n--)
      {
         y0H = *op2--;
         y0L = *op2--;
         x1L = *op1++;
         x1H = *op1++;
         if ((y0H < x0H) || (y0H == x0H  &&  y0L < x0L))
         {
            ZNP_ASSERT (y1H != 0 || y1L != 0);
            y1H -= (y1L-- == 0);
         }

         *res = zn_mod_reduce3_redc ((y1H << b1) + (y1L >> b2),
                                     (y1L << b1) + x0H, x0L, mod);
         res += s;

         ZNP_ASSERT (y1L != maskL || y1H != maskH);
         if (borrow)
            y1H += (++y1L == 0);
         borrow = ((x1H < y1H) || (x1H == y1H  &&  x1L < y1L));
         ZNP_SUB_WIDE (x1H, x1L, x1H, x1L, y1H, y1L);
         ZNP_SUB_WIDE (y1H, y1L, y0H, y0L, x0H, x0L);
         y1H &= maskH;
         x0L = x1L;
         x0H = x1H & maskH;
      }
   }
   else
   {
      // plain reduction version
      for (; n; n--)
      {
         y0H = *op2--;
         y0L = *op2--;
         x1L = *op1++;
         x1H = *op1++;
         if ((y0H < x0H) || (y0H == x0H  &&  y0L < x0L))
         {
            ZNP_ASSERT (y1H != 0 || y1L != 0);
            y1H -= (y1L-- == 0);
         }

         *res = zn_mod_reduce3 ((y1H << b1) + (y1L >> b2),
                                (y1L << b1) + x0H, x0L, mod);
         res += s;

         ZNP_ASSERT (y1L != maskL || y1H != maskH);
         if (borrow)
            y1H += (++y1L == 0);
         borrow = ((x1H < y1H) || (x1H == y1H  &&  x1L < y1L));
         ZNP_SUB_WIDE (x1H, x1L, x1H, x1L, y1H, y1L);
         ZNP_SUB_WIDE (y1H, y1L, y0H, y0L, x0H, x0L);
         y1H &= maskH;
         x0L = x1L;
         x0H = x1H & maskH;
      }
   }
}


/*
   Dispatches to one of the above routines depending on b.
*/
void
zn_array_recover_reduce (ulong* res, ptrdiff_t s, const ulong* op1,
                         const ulong* op2, size_t n, unsigned b,
                         int redc, const zn_mod_t mod)
{
   ZNP_ASSERT (b > 0  &&  2 * b <= 3 * ULONG_BITS);

   if (2 * b <= ULONG_BITS)
      zn_array_recover_reduce1 (res, s, op1, op2, n, b, redc, mod);
   else if (b < ULONG_BITS)
      zn_array_recover_reduce2 (res, s, op1, op2, n, b, redc, mod);
   else if (b == ULONG_BITS)
      zn_array_recover_reduce2b (res, s, op1, op2, n, b, redc, mod);
   else
      zn_array_recover_reduce3 (res, s, op1, op2, n, b, redc, mod);
}


// end of file ****************************************************************
