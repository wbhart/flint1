/*
  Copyright 2005, 2006 Damien Stehlé.

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

#include "F_mpz_mat.h"

double ** d_mat_init(int d, int n);

void d_mat_clear(double ** B);

void d_mat_print(double ** B, int * expo, int d, int n);

double d_vec_scalar_product(double * vec1, double * vec2, int n);

double d_vec_scalar_product_heuristic(double * vec1, double * vec2, int n, F_mpz_mat_t B, ulong kappa, ulong j, long exp_adj);

double d_vec_norm(double * vec, int n);

double d_2exp_vec_scalar_product(double * vec1, double * vec2, int n, int *cexpo, F_mpz_mat_t B, ulong kappa, ulong j);

double d_2exp_vec_norm(double * vec, int n, int *cexpo);



