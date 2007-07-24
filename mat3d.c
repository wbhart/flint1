/****************************************************************************

   mat3d.c: functions for dealing with 3D matrices over double floating 
                 point numbers

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mat3d.h"

/* 
   Compute the Gram-Schmidt Orthogonalisation of a set of 3 vectors b_i (given as the
   columns of the input matrix B_in). The orthogonal vectors b*_i are returned as the 
   columns of the output matrix B_out and a matrix of the Gram coefficients :
                   Q_ij = <b_i, b*_j>/<b*_j, b*_j>; 1 <= j < i <= 3
   are returned as the entries below the diagonal in the matrix Q. 
   
   Aliasing is permitted.
*/

void mat3dc_gram_schmidt(mat3dc_t Q, mat3dc_t B_out, mat3dc_t B_in)
{
   vec_3d_t temp;
   
   vec3d_set(B_out->c1, B_in->c1);
   Q->c1->x2 = vec3d_scalar_proj(B_in->c2, B_in->c1);
   Q->c1->x3 = vec3d_scalar_proj(B_in->c3, B_in->c1);
   
   vec3d_sub_scalar_mul(B_out->c2, B_in->c2, B_in->c1, Q->c1->x2);
   Q->c2->x3 = vec3d_scalar_proj(B_in->c3, B_out->c2);
   
   vec3d_sub_scalar_mul(B_out->c3, B_in->c3, B_out->c2, Q->c2->x3);
   vec3d_sub_scalar_mul(B_out->c3, B_out->c3, B_out->c1, Q->c1->x3);
}  

/* 
   Performs LLL reduction on the matrix B_in with the given delta value.
   We assume 1/4 < delta <= 1.
   
   Note B_out must be a mat3dc_p for efficiency.
   
   Not finished yet....
*/

void mat3dc_LLL(mat3dc_p B_out, mat3dc_p B_in, double delta)
{
   mat3dc_t Q;
   
   mat3dc_gram_schmidt(Q, B_out, B_in);
   
   // k = 2
   if (abs(Q->c1->x2) <= 0.5) 
   {
      unsigned long r = round(Q->c1->x2);
      vec3d_sub_scalar_mul(B_out->c2, B_out->c2, B_out->c1, (double) r);
      Q->c1->x2 -= (double) r;
   }
   
   
   // k = 3
   
   // k = 4 : terminate
}  

