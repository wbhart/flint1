/*
   pack.c:  bit-packing/unpacking for Kronecker substitution routines
   
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


/*
   Same as zn_array_pack(), but requires bits <= ULONG_BITS.
*/
#define zn_array_pack1 \
    ZNP_zn_array_pack1
void zn_array_pack1(mp_limb_t* res, const ulong* op, size_t len,
                    ptrdiff_t skip, unsigned bits, unsigned lead,
                    size_t requested)
{
   ZNP_ASSERT(bits > 0 && bits <= ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // where to write the next limb
   mp_limb_t* dest = res;
   
   // write leading zero-padding
   while (lead >= ULONG_BITS)
   {
      *dest++ = 0;
      lead -= ULONG_BITS;
   }

   // limb currently being filled
   mp_limb_t buf = 0;
   // number of bits used in buf; always in [0, ULONG_BITS)
   unsigned buf_bits = lead;
   unsigned old_buf_bits;
   
   for (; len > 0; len--, op += skip)
   {
      ZNP_ASSERT(bits >= ULONG_BITS  ||  *op < (1UL << bits));
      
      // put low bits of current input into buffer
      buf += *op << buf_bits;
      old_buf_bits = buf_bits;
      buf_bits += bits;
      if (buf_bits >= ULONG_BITS)
      {
         // buffer is full; flush it
         *dest++ = buf;
         buf_bits -= ULONG_BITS;
         // put remaining bits of current input into buffer
         buf = old_buf_bits ? (*op >> (ULONG_BITS - old_buf_bits)) : 0;
      }
   }
   
   // write last limb if it's non-empty
   if (buf_bits)
      *dest++ = buf;

   // zero-pad up to requested length
   if (requested)
   {
      size_t written = dest - res;
      ZNP_ASSERT(written <= requested);
      for (; written < requested; written++)
         *dest++ = 0;
   }
#else
#error Not nails-safe yet
#endif
}



void zn_array_pack(mp_limb_t* res, const ulong* op, size_t len,
                   ptrdiff_t skip, unsigned bits, unsigned lead,
                   size_t requested)
{
   ZNP_ASSERT(bits > 0 && bits < 3*ULONG_BITS);
   
   if (bits <= ULONG_BITS)
   {
      // use specialised version if bits is small enough
      zn_array_pack1(res, op, len, skip, bits, lead, requested);
      return;
   }
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // where to write the next limb
   mp_limb_t* dest = res;
   
   // write leading zero-padding
   while (lead >= ULONG_BITS)
   {
      *dest++ = 0;
      lead -= ULONG_BITS;
   }

   // limb currently being filled
   mp_limb_t buf = 0;
   // number of bits used in buf; always in [0, ULONG_BITS)
   unsigned buf_bits = lead;
   unsigned old_buf_bits;
   
   for (; len > 0; len--, op += skip)
   {
      ZNP_ASSERT(bits >= ULONG_BITS  ||  *op < (1UL << bits));
      
      // put low bits of current input into buffer
      buf += *op << buf_bits;
      old_buf_bits = buf_bits;
      buf_bits += bits;
      if (buf_bits >= ULONG_BITS)
      {
         // buffer is full; flush it
         *dest++ = buf;
         buf_bits -= ULONG_BITS;
         // put remaining bits of current input into buffer
         buf = old_buf_bits ? (*op >> (ULONG_BITS - old_buf_bits)) : 0;

         // write as many extra zeroes as necessary
         if (buf_bits >= ULONG_BITS)
         {
            *dest++ = buf;
            buf = 0;
            buf_bits -= ULONG_BITS;
            if (buf_bits >= ULONG_BITS)
            {
               *dest++ = 0;
               buf_bits -= ULONG_BITS;
               ZNP_ASSERT(buf_bits < ULONG_BITS);
            }
         }
      }
   }
   
   // write last limb if it's non-empty
   if (buf_bits)
      *dest++ = buf;

   // zero-pad up to requested length
   if (requested)
   {
      size_t written = dest - res;
      ZNP_ASSERT(written <= requested);
      for (; written < requested; written++)
         *dest++ = 0;
   }
#else
#error Not nails-safe yet
#endif
}



/*
   Same as zn_array_unpack(), but requires bits <= ULONG_BITS
   (i.e. writes one word per coefficient)
*/
#define zn_array_unpack1 \
    ZNP_zn_array_unpack1
void zn_array_unpack1(ulong* res, const mp_limb_t* op, size_t len,
                      unsigned bits, unsigned lead)
{
   ZNP_ASSERT(bits <= ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // limb we're currently extracting bits from
   mp_limb_t buf = 0;
   // number of bits currently in buf; always in [0, ULONG_BITS)
   unsigned buf_bits = 0;
   
   // skip over _lead_ leading bits
   while (lead >= GMP_NUMB_BITS)
   {
      lead -= GMP_NUMB_BITS;
      op++;
   }

   if (lead)
   {
      buf = *op++;
      buf >>= lead;
      buf_bits = ULONG_BITS - lead;
   }

   if (bits == ULONG_BITS)
   {
      // various special cases
      if (buf_bits)
      {
         for (; len > 0; len--)
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_bits);
            buf >>= (ULONG_BITS - buf_bits);
         }
      }
      else
      {
         for (; len > 0; len--)
            *res++ = *op++;
      }
   }
   else
   {
      ulong mask = (1UL << bits) - 1;

      for (; len > 0; len--)
      {
         if (bits <= buf_bits)
         {
            // buf contains all the bits we need
            *res++ = buf & mask;
            buf >>= bits;
            buf_bits -= bits;
         }
         else
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + ((buf << buf_bits) & mask);
            buf >>= (bits - buf_bits);
            buf_bits = ULONG_BITS - (bits - buf_bits);
         }
      }
   }
   
#else
#error Not nails-safe yet
#endif
}



/*
   Same as zn_array_unpack(), but requires ULONG_BITS < bits <= 2*ULONG_BITS
   (i.e. writes two words per coefficient)
*/
#define zn_array_unpack2 \
    ZNP_zn_array_unpack2
void zn_array_unpack2(ulong* res, const mp_limb_t* op, size_t len,
                      unsigned bits, unsigned lead)
{
   ZNP_ASSERT(bits > ULONG_BITS && bits <= 2*ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // limb we're currently extracting bits from
   mp_limb_t buf = 0;
   // number of bits currently in buf; always in [0, ULONG_BITS)
   unsigned buf_bits = 0;
   
   // skip over _lead_ leading bits
   while (lead >= GMP_NUMB_BITS)
   {
      lead -= GMP_NUMB_BITS;
      op++;
   }

   if (lead)
   {
      buf = *op++;
      buf >>= lead;
      buf_bits = ULONG_BITS - lead;
   }

   if (bits == 2*ULONG_BITS)
   {
      len = 2*len;
      
      // various special cases
      if (buf_bits)
      {
         for (; len > 0; len--)
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_bits);
            buf >>= (ULONG_BITS - buf_bits);
         }
      }
      else
      {
         for (; len > 0; len--)
            *res++ = *op++;
      }
   }
   else
   {
      bits -= ULONG_BITS;
      ulong mask = (1UL << bits) - 1;

      for (; len > 0; len--)
      {
         // shunt one whole limb through first
         if (buf_bits)
         {
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + (buf << buf_bits);
            buf >>= (ULONG_BITS - buf_bits);
         }
         else
            *res++ = *op++;
      
         // now handle the fractional limb
         if (bits <= buf_bits)
         {
            // buf contains all the bits we need
            *res++ = buf & mask;
            buf >>= bits;
            buf_bits -= bits;
         }
         else
         {
            // we need bits from both sides of a limb boundary
            ulong temp = buf;
            buf = *op++;
            *res++ = temp + ((buf << buf_bits) & mask);
            buf >>= (bits - buf_bits);
            buf_bits = ULONG_BITS - (bits - buf_bits);
         }
      }
   }
   
#else
#error Not nails-safe yet
#endif
}



/*
   Same as zn_array_unpack(), but requires 2*ULONG_BITS < bits < 3*ULONG_BITS
   (i.e. writes three words per coefficient)
*/
#define zn_array_unpack3 \
    ZNP_zn_array_unpack3
void zn_array_unpack3(ulong* res, const mp_limb_t* op, size_t len,
                      unsigned bits, unsigned lead)
{
   ZNP_ASSERT(bits > 2*ULONG_BITS && bits < 3*ULONG_BITS);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // limb we're currently extracting bits from
   mp_limb_t buf = 0;
   // number of bits currently in buf; always in [0, ULONG_BITS)
   unsigned buf_bits = 0;

   // skip over _lead_ leading bits
   while (lead >= GMP_NUMB_BITS)
   {
      lead -= GMP_NUMB_BITS;
      op++;
   }

   if (lead)
   {
      buf = *op++;
      buf >>= lead;
      buf_bits = ULONG_BITS - lead;
   }

   bits -= 2*ULONG_BITS;
   ulong mask = (1UL << bits) - 1;

   for (; len > 0; len--)
   {
      // shunt two whole limbs through first
      if (buf_bits)
      {
         ulong temp = buf;
         buf = *op++;
         *res++ = temp + (buf << buf_bits);
         buf >>= (ULONG_BITS - buf_bits);

         temp = buf;
         buf = *op++;
         *res++ = temp + (buf << buf_bits);
         buf >>= (ULONG_BITS - buf_bits);
      }
      else
      {
         *res++ = *op++;
         *res++ = *op++;
      }
   
      // now handle the fractional limb
      if (bits <= buf_bits)
      {
         // buf contains all the bits we need
         *res++ = buf & mask;
         buf >>= bits;
         buf_bits -= bits;
      }
      else
      {
         // we need bits from both sides of a limb boundary
         ulong temp = buf;
         buf = *op++;
         *res++ = temp + ((buf << buf_bits) & mask);
         buf >>= (bits - buf_bits);
         buf_bits = ULONG_BITS - (bits - buf_bits);
      }
   }
   
#else
#error Not nails-safe yet
#endif
}


void zn_array_unpack(ulong* res, const mp_limb_t* op, size_t len,
                     unsigned bits, unsigned lead)
{
   ZNP_ASSERT(bits >= 1 && bits <= 3 * ULONG_BITS);

   if (bits <= ULONG_BITS)
      zn_array_unpack1(res, op, len, bits, lead);
   else if (bits <= 2 * ULONG_BITS)
      zn_array_unpack2(res, op, len, bits, lead);
   else    // bits < 3 * ULONG_BITS
      zn_array_unpack3(res, op, len, bits, lead);
}



// end of file ****************************************************************
