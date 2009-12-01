/*
   pack.c:  bit-packing/unpacking for Kronecker substitution routines
   
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
   Same as zn_array_pack(), but requires b <= ULONG_BITS.
*/
#define zn_array_pack1 \
    ZNP_zn_array_pack1
void
zn_array_pack1 (mp_limb_t* res, const ulong* op, size_t n, ptrdiff_t s,
                unsigned b, unsigned k, size_t r)
{
   ZNP_ASSERT (b > 0 && b <= ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // where to write the next limb
   mp_limb_t* dest = res;
   
   // write leading zero-padding
   while (k >= ULONG_BITS)
   {
      *dest++ = 0;
      k -= ULONG_BITS;
   }

   // limb currently being filled
   mp_limb_t buf = 0;
   // number of bits used in buf; always in [0, ULONG_BITS)
   unsigned buf_b = k;
   unsigned buf_b_old;
   
   for (; n > 0; n--, op += s)
   {
      ZNP_ASSERT (b >= ULONG_BITS  ||  *op < (1UL << b));
      
      // put low bits of current input into buffer
      buf += *op << buf_b;
      buf_b_old = buf_b;
      buf_b += b;
      if (buf_b >= ULONG_BITS)
      {
         // buffer is full; flush it
         *dest++ = buf;
         buf_b -= ULONG_BITS;
         // put remaining bits of current input into buffer
         buf = buf_b_old ? (*op >> (ULONG_BITS - buf_b_old)) : 0;
      }
   }
   
   // write last limb if it's non-empty
   if (buf_b)
      *dest++ = buf;

   // zero-pad up to requested length
   if (r)
   {
      size_t written = dest - res;
      ZNP_ASSERT (written <= r);
      for (; written < r; written++)
         *dest++ = 0;
   }
#else
#error Not nails-safe yet
#endif
}



void
zn_array_pack (mp_limb_t* res, const ulong* op, size_t n, ptrdiff_t s,
               unsigned b, unsigned k, size_t r)
{
   ZNP_ASSERT (b > 0 && b < 3 * ULONG_BITS);
   
   if (b <= ULONG_BITS)
   {
      // use specialised version if b is small enough
      zn_array_pack1 (res, op, n, s, b, k, r);
      return;
   }
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // where to write the next limb
   mp_limb_t* dest = res;
   
   // write leading zero-padding
   while (k >= ULONG_BITS)
   {
      *dest++ = 0;
      k -= ULONG_BITS;
   }

   // limb currently being filled
   mp_limb_t buf = 0;
   // number of bits used in buf; always in [0, ULONG_BITS)
   unsigned buf_b = k;
   unsigned buf_b_old;
   
   for (; n > 0; n--, op += s)
   {
      ZNP_ASSERT (b >= ULONG_BITS  ||  *op < (1UL << b));
      
      // put low bits of current input into buffer
      buf += *op << buf_b;
      buf_b_old = buf_b;
      buf_b += b;
      if (buf_b >= ULONG_BITS)
      {
         // buffer is full; flush it
         *dest++ = buf;
         buf_b -= ULONG_BITS;
         // put remaining bits of current input into buffer
         buf = buf_b_old ? (*op >> (ULONG_BITS - buf_b_old)) : 0;

         // write as many extra zeroes as necessary
         if (buf_b >= ULONG_BITS)
         {
            *dest++ = buf;
            buf = 0;
            buf_b -= ULONG_BITS;
            if (buf_b >= ULONG_BITS)
            {
               *dest++ = 0;
               buf_b -= ULONG_BITS;
               ZNP_ASSERT (buf_b < ULONG_BITS);
            }
         }
      }
   }
   
   // write last limb if it's non-empty
   if (buf_b)
      *dest++ = buf;

   // zero-pad up to requested length
   if (r)
   {
      size_t written = dest - res;
      ZNP_ASSERT (written <= r);
      for (; written < r; written++)
         *dest++ = 0;
   }
#else
#error Not nails-safe yet
#endif
}



/*
   Same as zn_array_unpack(), but requires b <= ULONG_BITS
   (i.e. writes one word per coefficient)
*/
#define zn_array_unpack1 \
    ZNP_zn_array_unpack1
void
zn_array_unpack1 (ulong* res, const mp_limb_t* op, size_t n, unsigned b,
                  unsigned k)
{
   ZNP_ASSERT (b <= ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // limb we're currently extracting bits from
   mp_limb_t buf = 0;
   // number of bits currently in buf; always in [0, ULONG_BITS)
   unsigned buf_b = 0;
   
   // skip over k leading bits
   while (k >= GMP_NUMB_BITS)
   {
      k -= GMP_NUMB_BITS;
      op++;
   }

   if (k)
   {
      buf = *op++;
      buf >>= k;
      buf_b = ULONG_BITS - k;
   }

   if (b == ULONG_BITS)
   {
      // various special cases
      if (buf_b)
      {
         for (; n > 0; n--)
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_b);
            buf >>= (ULONG_BITS - buf_b);
         }
      }
      else
      {
         for (; n > 0; n--)
            *res++ = *op++;
      }
   }
   else
   {
      ulong mask = (1UL << b) - 1;

      for (; n > 0; n--)
      {
         if (b <= buf_b)
         {
            // buf contains all the bits we need
            *res++ = buf & mask;
            buf >>= b;
            buf_b -= b;
         }
         else
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + ((buf << buf_b) & mask);
            buf >>= (b - buf_b);
            buf_b = ULONG_BITS - (b - buf_b);
         }
      }
   }
   
#else
#error Not nails-safe yet
#endif
}



/*
   Same as zn_array_unpack(), but requires ULONG_BITS < b <= 2 * ULONG_BITS
   (i.e. writes two words per coefficient)
*/
#define zn_array_unpack2 \
    ZNP_zn_array_unpack2
void
zn_array_unpack2 (ulong* res, const mp_limb_t* op, size_t n, unsigned b,
                  unsigned k)
{
   ZNP_ASSERT (b > ULONG_BITS && b <= 2 * ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // limb we're currently extracting bits from
   mp_limb_t buf = 0;
   // number of bits currently in buf; always in [0, ULONG_BITS)
   unsigned buf_b = 0;
   
   // skip over k leading bits
   while (k >= GMP_NUMB_BITS)
   {
      k -= GMP_NUMB_BITS;
      op++;
   }

   if (k)
   {
      buf = *op++;
      buf >>= k;
      buf_b = ULONG_BITS - k;
   }

   if (b == 2 * ULONG_BITS)
   {
      n *= 2;
      
      // various special cases
      if (buf_b)
      {
         for (; n > 0; n--)
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_b);
            buf >>= (ULONG_BITS - buf_b);
         }
      }
      else
      {
         for (; n > 0; n--)
            *res++ = *op++;
      }
   }
   else
   {
      b -= ULONG_BITS;
      ulong mask = (1UL << b) - 1;

      for (; n > 0; n--)
      {
         // shunt one whole limb through first
         if (buf_b)
         {
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_b);
            buf >>= (ULONG_BITS - buf_b);
         }
         else
            *res++ = *op++;
      
         // now handle the fractional limb
         if (b <= buf_b)
         {
            // buf contains all the bits we need
            *res++ = buf & mask;
            buf >>= b;
            buf_b -= b;
         }
         else
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + ((buf << buf_b) & mask);
            buf >>= (b - buf_b);
            buf_b = ULONG_BITS - (b - buf_b);
         }
      }
   }
   
#else
#error Not nails-safe yet
#endif
}



/*
   Same as zn_array_unpack(), but requires 2 * ULONG_BITS < b < 3 * ULONG_BITS
   (i.e. writes three words per coefficient)
*/
#define zn_array_unpack3 \
    ZNP_zn_array_unpack3
void
zn_array_unpack3 (ulong* res, const mp_limb_t* op, size_t n, unsigned b,
                  unsigned k)
{
   ZNP_ASSERT (b > 2 * ULONG_BITS && b < 3 * ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // limb we're currently extracting bits from
   mp_limb_t buf = 0;
   // number of bits currently in buf; always in [0, ULONG_BITS)
   unsigned buf_b = 0;

   // skip over k leading bits
   while (k >= GMP_NUMB_BITS)
   {
      k -= GMP_NUMB_BITS;
      op++;
   }

   if (k)
   {
      buf = *op++;
      buf >>= k;
      buf_b = ULONG_BITS - k;
   }

   b -= 2 * ULONG_BITS;
   ulong mask = (1UL << b) - 1;

   for (; n > 0; n--)
   {
      // shunt two whole limbs through first
      if (buf_b)
      {
         ulong temp = buf;
         buf = *op++;
         *res++ = temp + (buf << buf_b);
         buf >>= (ULONG_BITS - buf_b);

         temp = buf;
         buf = *op++;
         *res++ = temp + (buf << buf_b);
         buf >>= (ULONG_BITS - buf_b);
      }
      else
      {
         *res++ = *op++;
         *res++ = *op++;
      }
   
      // now handle the fractional limb
      if (b <= buf_b)
      {
         // buf contains all the bits we need
         *res++ = buf & mask;
         buf >>= b;
         buf_b -= b;
      }
      else
      {
         // we need bits from both sides of a limb boundary
         ulong temp = buf;
         buf = *op++;
         *res++ = temp + ((buf << buf_b) & mask);
         buf >>= (b - buf_b);
         buf_b = ULONG_BITS - (b - buf_b);
      }
   }
   
#else
#error Not nails-safe yet
#endif
}


void
zn_array_unpack (ulong* res, const mp_limb_t* op, size_t n, unsigned b,
                 unsigned k)
{
   ZNP_ASSERT (b >= 1 && b <= 3 * ULONG_BITS);

   if (b <= ULONG_BITS)
      zn_array_unpack1 (res, op, n, b, k);
   else if (b <= 2 * ULONG_BITS)
      zn_array_unpack2 (res, op, n, b, k);
   else    // b < 3 * ULONG_BITS
      zn_array_unpack3 (res, op, n, b, k);
}



// end of file ****************************************************************
