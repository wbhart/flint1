/*
   array.c:  simple operations on arrays mod m
   
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


int
zn_array_cmp (const ulong* op1, const ulong* op2, size_t n)
{
   for (; n > 0; n--)
      if (*op1++ != *op2++)
         return 1;

   return 0;
}


void
zn_array_copy (ulong* res, const ulong* op, size_t n)
{
   for (; n > 0; n--)
      *res++ = *op++;
}


void
zn_array_neg (ulong* res, const ulong* op, size_t n, const zn_mod_t mod)
{
   for (; n > 0; n--)
      *res++ = zn_mod_neg (*op++, mod);
}



void
zn_array_scalar_mul_or_copy (ulong* res, const ulong* op, size_t n,
                             ulong x, const zn_mod_t mod)
{
   if (x != 1)
      zn_array_scalar_mul (res, op, n, x, mod);
   else if (res != op)
      zn_array_copy (res, op, n);
}



/*
   Same as zn_array_scalar_mul, but:
      * always uses REDC reduction (requires modulus is odd);
      * requires that residues fit into half a word.
*/
#define _zn_array_scalar_mul_redc_v1 \
    ZNP__zn_array_scalar_mul_redc_v1
void
_zn_array_scalar_mul_redc_v1 (ulong* res, const ulong* op, size_t n, ulong x,
                              const zn_mod_t mod)
{
   ZNP_ASSERT (mod->bits <= ULONG_BITS/2);
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (x < mod->m);

   for (; n; n--, op++, res++)
      *res = zn_mod_reduce_redc ((*op) * x, mod);
}


/*
   Same as zn_array_scalar_mul, but:
      * always uses REDC reduction (requires modulus is odd);
      * requires that modulus is slim.
*/
#define _zn_array_scalar_mul_redc_v2 \
    ZNP__zn_array_scalar_mul_redc_v2
void
_zn_array_scalar_mul_redc_v2 (ulong* res, const ulong* op, size_t n, ulong x,
                              const zn_mod_t mod)
{
   ZNP_ASSERT (zn_mod_is_slim (mod));
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (x < mod->m);

   for (; n; n--, op++, res++)
   {
      ulong hi, lo;
      ZNP_MUL_WIDE (hi, lo, *op, x);
      *res = zn_mod_reduce_wide_redc_slim (hi, lo, mod);
   }
}


/*
   Same as zn_array_scalar_mul, but:
      * always uses REDC reduction (requires modulus is odd).
*/
#define _zn_array_scalar_mul_redc_v3 \
    ZNP__zn_array_scalar_mul_redc_v3
void
_zn_array_scalar_mul_redc_v3 (ulong* res, const ulong* op, size_t n, ulong x,
                              const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (x < mod->m);

   for (; n; n--, op++, res++)
   {
      ulong hi, lo;
      ZNP_MUL_WIDE (hi, lo, *op, x);
      *res = zn_mod_reduce_wide_redc (hi, lo, mod);
   }
}


/*
   Same as zn_array_scalar_mul, but always uses REDC reduction (requires that
   modulus is odd).
   
   Dispatches to one of the three versions above, depending on modulus size.
*/
#define _zn_array_scalar_mul_redc \
    ZNP__zn_array_scalar_mul_redc
void
_zn_array_scalar_mul_redc (ulong* res, const ulong* op, size_t n, ulong x,
                           const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (x < mod->m);
   
   if (mod->bits <= ULONG_BITS/2)
      _zn_array_scalar_mul_redc_v1 (res, op, n, x, mod);
   else if (zn_mod_is_slim (mod))
      _zn_array_scalar_mul_redc_v2 (res, op, n, x, mod);
   else
      _zn_array_scalar_mul_redc_v3 (res, op, n, x, mod);
}


/*
   Same as zn_array_scalar_mul, but:
      * always uses plain reduction;
      * requires that residues fit into half a word.
*/
#define _zn_array_scalar_mul_plain_v1 \
    ZNP__zn_array_scalar_mul_plain_v1
void
_zn_array_scalar_mul_plain_v1 (ulong* res, const ulong* op, size_t n, ulong x,
                               const zn_mod_t mod)
{
   ZNP_ASSERT (mod->bits <= ULONG_BITS/2);
   ZNP_ASSERT (x < mod->m);

   for (; n; n--, op++, res++)
      *res = zn_mod_reduce ((*op) * x, mod);
}


/*
   Same as zn_array_scalar_mul, but:
      * always uses plain reduction.
*/
#define _zn_array_scalar_mul_plain_v2 \
    ZNP__zn_array_scalar_mul_plain_v2
void
_zn_array_scalar_mul_plain_v2 (ulong* res, const ulong* op, size_t n, ulong x,
                               const zn_mod_t mod)
{
   ZNP_ASSERT (x < mod->m);

   for (; n; n--, op++, res++)
   {
      ulong hi, lo;
      ZNP_MUL_WIDE (hi, lo, *op, x);
      *res = zn_mod_reduce_wide (hi, lo, mod);
   }
}


/*
   Same as zn_array_scalar_mul, but always uses plain reduction.
   
   Dispatches to one of the versions above, depending on modulus size.
*/
#define _zn_array_scalar_mul_plain \
    ZNP__zn_array_scalar_mul_plain
void
_zn_array_scalar_mul_plain (ulong* res, const ulong* op, size_t n, ulong x,
                            const zn_mod_t mod)
{
   ZNP_ASSERT (x < mod->m);
   
   if (mod->bits <= ULONG_BITS/2)
      _zn_array_scalar_mul_plain_v1 (res, op, n, x, mod);
   else
      _zn_array_scalar_mul_plain_v2 (res, op, n, x, mod);
}


void
_zn_array_scalar_mul (ulong* res, const ulong* op, size_t n, ulong x,
                      int redc, const zn_mod_t mod)
{
   if (redc)
      _zn_array_scalar_mul_redc (res, op, n, x, mod);
   else
      _zn_array_scalar_mul_plain (res, op, n, x, mod);
}


void zn_array_scalar_mul (ulong* res, const ulong* op, size_t n, ulong x,
                          const zn_mod_t mod)
{
   ZNP_ASSERT (x < mod->m);

   // Do plain reduction if the vector is really short, or if the modulus
   // is even (in which case REDC reduction is not available).
   if (n < 5  ||  !(mod->m & 1))
   {
      _zn_array_scalar_mul_plain (res, op, n, x, mod);
   }
   else
   {
      // modulus is odd, and vector is not too short, so we can go faster
      // by adjusting the multiplier and using REDC reduction
      _zn_array_scalar_mul_redc (res, op, n,
                                 zn_mod_mul_redc (x, mod->B2, mod), mod);
   }
}



void
zn_array_sub (ulong* res, const ulong* op1, const ulong* op2, size_t n,
              const zn_mod_t mod)
{
   if (zn_mod_is_slim (mod))
      for (; n; n--)
         *res++ = zn_mod_sub_slim (*op1++, *op2++, mod);
   else
      for (; n; n--)
         *res++ = zn_mod_sub (*op1++, *op2++, mod);
}


ulong*
zn_skip_array_signed_add (ulong* res, ptrdiff_t s, size_t n,
                          const ulong* op1, int neg1,
                          const ulong* op2, int neg2,
                          const zn_mod_t mod)
{
   if (zn_mod_is_slim (mod))
   {
      // slim version
      if (neg1)
      {
         if (neg2)
            // res = -(op1 + op2)
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_neg (zn_mod_add_slim (*op1, *op2, mod), mod);
         else
            // res = op2 - op1
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_sub_slim (*op2, *op1, mod);
      }
      else
      {
         if (neg2)
            // res = op1 - op2
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_sub_slim (*op1, *op2, mod);
         else
            // res = op1 + op2
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_add_slim (*op1, *op2, mod);
      }
   }
   else
   {
      // non-slim version
      if (neg1)
      {
         if (neg2)
            // res = -(op1 + op2)
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_neg (zn_mod_add (*op1, *op2, mod), mod);
         else
            // res = op2 - op1
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_sub (*op2, *op1, mod);
      }
      else
      {
         if (neg2)
            // res = op1 - op2
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_sub (*op1, *op2, mod);
         else
            // res = op1 + op2
            for (; n > 0; n--, res += s, op1++, op2++)
               *res = zn_mod_add (*op1, *op2, mod);
      }
   }
   
   return res;
}


// end of file ****************************************************************
