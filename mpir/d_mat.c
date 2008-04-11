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

#include "d_mat.h"
#include "fmpz.h"

double* init_vecf (int n)
{
   double *a;
   a = (double*) malloc (n * sizeof (double));
   return a;
}

void init_matrixf (double **B, int d, int n)
{
   int i;

   B = (double **) malloc (d * sizeof (double*));
   for (i=0; i<d; i++) B[i] = (double*) malloc (n * sizeof (double));
}

void clear_matrixf (double **B, int d)
{
   int i;
   for (i=0; i<d; i++) free (B[i]);
   free (B);
}

void Print_matf (double **B, int *expo, int d, int n) 
{
   int i, j; 

   printf("[");
   for (i=0;i<d;i++) 
   {
      printf("[");
      for (j=0;j<n;j++) 
      { 
         printf("%E", B[i][j]); 
         if (j < n-1) printf(" "); 
      }
      printf("] * 2^%d\n", expo[i]); 
   }  
   printf("]\n"); 
}


int set_line(double *appv, fmpz_t * v, int n)
{
   long *exp, i, maxexp = 0L;
   exp = (long *) malloc(n * sizeof(long)); 
  
   for (i=0; i<n; i++)
   {
      appv[i] = fmpz_get_d_2exp(&exp[i], v + i);
      if (exp[i]>maxexp) maxexp=exp[i];
   }

   for (i=0; i<n; i++) appv[i] = ldexp( appv[i], exp[i]-maxexp);

   free(exp);
   return maxexp;
}


