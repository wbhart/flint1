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
/****************************************************************************

NTL-interface.cpp: Functions for conversion between NTL and FLINT format

Copyright (C) 2007, William Hart

*****************************************************************************/

#include <cstdio>

#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/lip.h>
#include <NTL/ctools.h>
#include <NTL/g_lip.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "NTL-interface.h"
#include "NTL-interface_impl.h"

#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

NTL_CLIENT

unsigned long ZZ_limbs(ZZ z)
{
   if (z.rep) return SIZE(z.rep);
   else return 0;
}

void ZZ_to_fmpz(fmpz_t output, ZZ z)
{
   _ntl_gbigint x = z.rep;
   
   if (!x) 
   {
      output[0] = 0L;
      return;
   }
   
   long lw = SIZE(x);
   mp_limb_t *xp = DATA(x);

   F_mpn_copy(output + 1, xp, lw);
   
   output[0] = lw;
}

void fmpz_to_ZZ(ZZ& output, fmpz_t z)
{
   mp_limb_t *xp;
   _ntl_gbigint *x = &output.rep;
   long lw = FLINT_ABS(z[0]);;
   
   if (lw == 0) {
      if (*x) SIZE(*x) = 0;
      return;
   }

   _ntl_gsetlength(x, lw); 
   xp = DATA(*x);

   F_mpn_copy(xp, z + 1, lw);
   
   SIZE(*x) = lw; 
}

void ZZX_to_fmpz_poly(ZZX& output, fmpz_poly_t poly)
{
}

void fmpz_poly_to_ZZX(fmpz_poly_t output, ZZX poly)
{
}

 
     
