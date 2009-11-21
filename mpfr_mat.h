/*
  Copyright 2005, 2006 Damien Stehlé.
  Copyright 2009, William Hart.

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

#include <gmp.h>
#include <mpfr.h>

mpfr_t ** mpfr_mat_init(int d, int n);

void mpfr_mat_clear(mpfr_t ** B, int d, int n);

void mpfr_mat_print(mpfr_t ** B, int d, int n);

int mpfr_vec_scalar_product(mpfr_t sp, mpfr_t * vec1, mpfr_t * vec2, int n);

void mpfr_vec_norm(mpfr_t norm, mpfr_t * vec, int n);



