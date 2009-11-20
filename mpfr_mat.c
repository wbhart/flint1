/*
   Copyright 2005, 2006 Damien Stehlé.
   Copyright 2009, William Hart

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctype.h>
#include <limits.h>
#include <gmp.h>
#include <mpfr.h>

#include "mpfr_mat.h"

mpfr_t ** mpfr_mat_init(int d, int n)
{
   mpfr_t ** B;

   B = (mpfr_t **) malloc (d*sizeof(mpfr_t *) + n*d*sizeof(mpfr_t));
   B[0] = (mpfr_t *) (B + d);
	for (long i = 1; i < d; i++) B[i] = B[i-1] + n;

   for (long i = 0; i < d; i++)
      for (long j = 0; j < n; j++)
         mpfr_init(B[i][j]);

	return B;
}

void mpfr_mat_clear(mpfr_t ** B, int d, int n)
{
   for (long i = 0; i < d; i++)
      for (long j = 0; j < n; j++)
         mpfr_clear(B[i][j]);

   free(B);
}

void mpfr_mat_print(mpfr_t ** B, int d, int n) 
{
   long i, j; 

   printf("[");
   for (i = 0; i < d; i++) 
   {
      printf("[");
      for (j = 0; j < n; j++) 
      { 
         mpfr_printf("%.128Rf", B[i][j]); 
         if (j < n-1) printf(" "); 
      }
      printf("]"); 
   }  
   printf("]\n"); 
}

int mpfr_vec_scalar_product(mpfr_t sp, mpfr_t * vec1, mpfr_t * vec2, int n)
{
  mpfr_t tmp, tmp2;
  mpfr_init(tmp);
  mpfr_init(tmp2);
  
  int res;

  mpfr_mul(sp, vec1[0], vec2[0], GMP_RNDN);
  
  for (long i = 1; i < n; i++)
  {
     mpfr_mul(tmp, vec1[i], vec2[i], GMP_RNDN);
     mpfr_add(sp, sp, tmp, GMP_RNDN);
  }

  mpfr_vec_norm(tmp, vec1, n);
  mpfr_vec_norm(tmp2, vec2, n);
  mpfr_mul(tmp, tmp, tmp2, GMP_RNDN);
  mpfr_div_2ui(tmp, tmp, 70, GMP_RNDN);
  mpfr_mul(tmp2, sp, sp, GMP_RNDN);

  if (mpfr_cmp(tmp2, tmp) <= 0) res = 0;
  else res = 1;

  mpfr_clear(tmp);
  mpfr_clear(tmp2);

  return res;
} 

void mpfr_vec_norm(mpfr_t norm, mpfr_t * vec, int n)
{
  mpfr_t tmp;
  mpfr_init(tmp);
  
  mpfr_mul(norm, vec[0], vec[0], GMP_RNDN);
  
  for (long i = 1 ; i < n ; i++)
  {
     mpfr_mul(tmp, vec[i], vec[i], GMP_RNDN);
     mpfr_add(norm, norm, tmp, GMP_RNDN);
  }

  mpfr_clear(tmp);

  return;
} 
