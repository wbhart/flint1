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

/* 
   Determines if the entries in r1 are equal to r2 to within about 2^-exp
*/
static inline 
int d_mat_row_equal(double * r1, double * r2, ulong start, ulong n, ulong exp)
{
   ulong i;
   double eps = ldexp(1.0, -exp);

   for (i = start; i < start + n; i++)
      if (fabs(r1[i] - r2[i]) > eps) return 0;

   return 1;
}

void d_mat_row_add(double * r1, double * r2, double * r3, ulong start, ulong n);

void d_mat_row_sub(double * r1, double * r2, double * r3, ulong start, ulong n);

void d_mat_print(double ** B, int * expo, int d, int n);

double d_vec_scalar_product(double * vec1, double * vec2, int n);

double d_vec_norm(double * vec, int n);

double d_vec_scalar_product_heuristic(double * vec1, double * vec2, int n, F_mpz_mat_t B, ulong kappa, ulong j, long exp_adj);



