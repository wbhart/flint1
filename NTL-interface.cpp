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
#include <NTL/mat_ZZ.h>
#include <NTL/lip.h>
#include <NTL/ctools.h>
#include <NTL/g_lip.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"
#include "fmpz_poly.h"
#include "NTL-interface.h"

#define SIZE(p) (((long *) (p))[1])
#define DATA(p) ((mp_limb_t *) (((long *) (p)) + 2))

NTL_CLIENT

unsigned long ZZ_limbs(const ZZ& z)
{
   if (z.rep) return FLINT_ABS(SIZE(z.rep));
   else return 0;
}

unsigned long ZZX_maxlimbs(const ZZX& z)
{
   unsigned long length = deg(z)+1;
   unsigned long maxlimbs = 0;
   unsigned long newlimbs, i;
   const ZZ *ap; 

   if (length == 0) return 0;
   
   for (i = 0, ap = z.rep.elts(); i < length; ap++, i++)
   {
      newlimbs = ZZ_limbs(*ap);
      if (newlimbs > maxlimbs) maxlimbs = newlimbs;
   }
   return maxlimbs;
}

void ZZ_to_fmpz(fmpz_t output, const ZZ& z)
{
   _ntl_gbigint x = z.rep;
   
   if (!x) 
   {
      output[0] = 0L;
      return;
   }
   
   unsigned long lw = ZZ_limbs(z);
   mp_limb_t *xp = DATA(x);

   F_mpn_copy(output + 1, xp, lw);
   
   if (z < 0L) output[0] = -lw;
   else output[0] = lw;
}

void ZZ_to_F_mpz(F_mpz_t output, const ZZ& z)
{
   _ntl_gbigint x = z.rep;
   
   if (!x) 
   {
      F_mpz_zero(output);
      return;
   }
   
   unsigned long lw = ZZ_limbs(z);
   mp_limb_t * xp = DATA(x);

   F_mpz_set_limbs(output, xp, lw);
   
   if (z < 0L) F_mpz_neg(output, output);
}

void fmpz_to_ZZ(ZZ& output, const fmpz_t z)
{
   mp_limb_t *xp;
   _ntl_gbigint *x = &output.rep;
   long lw = FLINT_ABS(z[0]);;
   
   if (lw == 0) 
	{
      if (*x) SIZE(*x) = 0;
      return;
   }

   _ntl_gsetlength(x, lw); 
   xp = DATA(*x);

   F_mpn_copy(xp, z + 1, lw);
   
   if ((long) z[0] < 0L) SIZE(*x) = -lw;
   else SIZE(*x) = lw;
}

void F_mpz_to_ZZ(ZZ& output, const F_mpz_t z)
{
   mp_limb_t *xp;
   _ntl_gbigint *x = &output.rep;
   long lw = F_mpz_size(z);
   
   if (lw == 0) 
	{
      if (*x) SIZE(*x) = 0;
      return;
   }

   _ntl_gsetlength(x, lw); 
   xp = DATA(*x);

	F_mpz_get_limbs(xp, z);
	
   if (F_mpz_sgn(z) < 0) SIZE(*x) = -lw;
   else SIZE(*x) = lw;
}

void fmpz_poly_to_ZZX(ZZX& output, const fmpz_poly_t poly)
{
   unsigned long length = poly->length;
   unsigned long i;
   fmpz_t coeff;
   ZZ *ap; 
   
   if (length == 0)
   {
      output = 0;
      return;
   }
   
   output.rep.SetLength(length);
   
   for (i = 0, ap = output.rep.elts(); i < length; ap++, i++)
   {
      coeff = fmpz_poly_get_coeff_ptr(poly, i);
      fmpz_to_ZZ(*ap, coeff);
   }
}

void ZZX_to_fmpz_poly(fmpz_poly_t output, const ZZX& poly)
{
  
   unsigned long length = deg(poly) + 1;
   unsigned long limbs = ZZX_maxlimbs(poly);
   unsigned long i;
   
   fmpz_t coeff_f;
   const ZZ *ap; 
   
   if (length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   
   fmpz_poly_fit_length(output, length);
   fmpz_poly_fit_limbs(output, limbs);
   
   output->length = length;
   for (i = 0, ap = poly.rep.elts(); i < length; ap++, i++)
   {
      coeff_f = fmpz_poly_get_coeff_ptr(output, i);
      ZZ_to_fmpz(coeff_f, *ap);
   }
}

void mat_ZZ_to_F_mpz_mat(F_mpz_mat_t output, const mat_ZZ& mat)
{
	ulong r = mat.NumRows();
	ulong c = mat.NumCols();

	for (ulong i = 0; i < r; i++)
		for (ulong j = 0; j < c; j++)
		{
			ZZ_to_F_mpz(output->rows[i] + j, mat[i][j]);
		}
}

void F_mpz_mat_to_mat_ZZ(mat_ZZ& output, const F_mpz_mat_t mat)
{
	ulong r = mat->r;
	ulong c = mat->c;

	for (ulong i = 0; i < r; i++)
		for (ulong j = 0; j < c; j++)
		{
			F_mpz_to_ZZ(output[i][j], mat->rows[i] + j);
		}
}

 
     
