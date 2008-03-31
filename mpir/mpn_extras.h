/*============================================================================

    Copyright (C) 2007, William Hart, David Harvey

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

#include "mpir.h"

/*============================================================================

"mpn-wannabe" code.

These are functions that we wish were in GMP's mpn layer.

=============================================================================*/

/*
Computes the negation of a multiple-precision integer in 2's complement.

Input is count limbs stored at src. Output is stored at dest.

src and dest can be the same buffer. If they're not, they should be disjoint.
*/

static inline 
void F_mpn_com(mp_limb_t* dest, mp_limb_t* src, ulong count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = ~src[i];
}

static inline 
void F_mpn_negate(mp_limb_t* dest, mp_limb_t* src, ulong count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = ~src[i];
   mpn_add_1(dest, dest, count, 1);
}

/*
Copies a bunch of limbs from one buffer to another.

Input is count limbs stored at src. Output is stored at dest.

src and dest can be the same buffer. If they're not, they should be disjoint.

 */
static inline 
void F_mpn_copy(mp_limb_t* dest, const mp_limb_t* src, ulong count)
{
   for (long i = count - 1; i >= 0; i--)
   {
      dest[i] = src[i];
   }
}

static inline 
void F_mpn_copy_forward(mp_limb_t* dest, const mp_limb_t* src, ulong count)
{
   for (long i = 0; i < count; i++)
   {
      dest[i] = src[i];
   }
}

/*
   Sets a bunch of limbs to zero.
*/

static inline 
void F_mpn_clear(mp_limb_t* dest, ulong count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = 0;
}

/*
Sets a bunch of limbs to 0xfff....
 */
static inline 
void F_mpn_set(mp_limb_t* dest, ulong count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = (mp_limb_t)(-1L);
}

#endif

