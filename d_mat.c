/*
   Copyright 2005, 2006 Damien Stehlé.
   Copyright 2010 William Hart.

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

/****************************************************************************

   d_mat.c: Vectors and matrices with double precision floating point entries.

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "d_mat.h"

double ** d_mat_init(ulong r, ulong c)
{
   double ** B;

   B = (double **) malloc (r*sizeof(double*) + r*c*sizeof(double));
   B[0] = (double *) (B + r);
	long i;
	for (i = 1; i < r; i++) B[i] = B[i-1] + c;

	return B;
}

void d_mat_clear(double ** B)
{
   free(B);
}

void d_mat_print(double ** B, int * expo, ulong r, ulong c) 
{
   long i, j; 

   printf("[");
   for (i = 0; i < r; i++) 
   {
      printf("[");
      for (j = 0; j < c; j++) 
      { 
         printf("%E", B[i][j]); 
         if (j < c-1) printf(" "); 
      }
      printf("] * 2^%d\n", expo[i]); 
   }  
   printf("]\n"); 
}

void _d_vec_add(double * r1, double * r2, double * r3, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
	  r1[i] = r2[i] + r3[i];
}

void _d_vec_sub(double * r1, double * r2, double * r3, ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
	  r1[i] = r2[i] - r3[i];
}

double _d_vec_scalar_product(double * vec1, double * vec2, ulong n)
{
  double sum;

  sum = vec1[0] * vec2[0];
  long i;
  for (i = 1; i < n; i++)
     sum += vec1[i] * vec2[i];

  return sum;
} 

double _d_vec_norm(double * vec, ulong n)
{
  double sum;

  sum = vec[0] * vec[0];
  long i;
  for (i = 1 ; i < n ; i++)
     sum += vec[i] * vec[i];

  return sum;
} 


