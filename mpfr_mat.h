/*
  Copyright 2005, 2006 Damien Stehlé.
  Copyright 2009, 2010 William Hart.
  Copyright 2009, 2010 Andy Novocin.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#ifndef MPFR_MAT_H
#define MPFR_MAT_H

#ifdef __cplusplus
 extern "C" {
#endif

#include <gmp.h>
#include <mpfr.h>
#include "flint.h"

__mpfr_struct ** mpfr_mat_init(int d, int n);

__mpfr_struct ** mpfr_mat_init2(int d, int n, mp_prec_t prec);

static inline
int _mpfr_vec_equal(__mpfr_struct * vec1, __mpfr_struct * vec2, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
      if (mpfr_cmp(vec1 + i, vec2 + i) != 0)
	     return 0;

   return 1;
}

void mpfr_mat_clear(__mpfr_struct ** B, int d, int n);

void mpfr_mat_print(__mpfr_struct ** B, int d, int n);

int _mpfr_vec_scalar_product(mpfr_t sp, __mpfr_struct * vec1, __mpfr_struct * vec2, int n);

void _mpfr_vec_norm(mpfr_t norm, __mpfr_struct * vec, int n);

int _mpfr_vec_scalar_product2(mpfr_t sp, __mpfr_struct * vec1, __mpfr_struct * vec2, int n, mp_prec_t prec);

void _mpfr_vec_norm2(mpfr_t norm, __mpfr_struct * vec, int n, mp_prec_t prec);

#ifdef __cplusplus
 }
#endif
 
#endif


