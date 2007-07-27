/****************************************************************************

   mat3d.c: functions for dealing with 3D matrices over double floating 
                 point numbers

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mat3d.h"
#include "vec3d.h"
#include "memory-manager.h"

#define ROW(R, x) ((R)[x-1])
#define COL(R, x) ((R)[x-1])
#define MAT_R(R, x, y) ((R)[x-1][y-1])
#define MAT_C(R, x, y) ((R)[y-1][x-1])

void z_mat3dr_stack_init(z_mat3dr_t * C)
{
   (*C) = (unsigned long **) flint_stack_alloc(12);
   (*C)[0] = (unsigned long *) (C+3);
   (*C)[1] = (unsigned long *) (C+6);
   (*C)[2] = (unsigned long *) (C+9); 
}

void z_mat3dr_stack_clear(void)
{
   flint_stack_release();
}

void z_mat3dc_stack_init(z_mat3dc_t * C)
{
   (*C) = (unsigned long **) flint_stack_alloc(12);
   (*C)[0] = (unsigned long *) (C+3);
   (*C)[1] = (unsigned long *) (C+6);
   (*C)[2] = (unsigned long *) (C+9); 
}

void z_mat3dc_stack_clear(void)
{
   flint_stack_release();
}

void z_mat3dr_swap12(z_mat3dr_t C)
{
   z_vec3d temp;
   temp = ROW(C, 1);
   ROW(C, 1) = ROW(C, 2);
   ROW(C, 2) = temp;
} 

void z_mat3dr_swap13(z_mat3dr_t C)
{
   z_vec3d temp;
   temp = ROW(C, 1);
   ROW(C, 1) = ROW(C, 3);
   ROW(C, 3) = temp;
} 

void z_mat3dr_swap23(z_mat3dr_t C)
{
   z_vec3d temp;
   temp = ROW(C, 2);
   ROW(C, 2) = ROW(C, 3);
   ROW(C, 3) = temp;
} 

void z_mat3dr_set_identity(z_mat3dr_t C)
{
   MAT_R(C, 1, 1) = 1;
   MAT_R(C, 1, 2) = 0;
   MAT_R(C, 1, 3) = 0;
   MAT_R(C, 2, 1) = 0;
   MAT_R(C, 2, 2) = 1;
   MAT_R(C, 2, 3) = 0;
   MAT_R(C, 3, 1) = 0;
   MAT_R(C, 3, 2) = 0;
   MAT_R(C, 3, 3) = 1;
}

void d_mat3dr_stack_init(d_mat3dr_t * C)
{
   (*C) = (double **) flint_stack_alloc(12);
   (*C)[0] = (double *) (C+3);
   (*C)[1] = (double *) (C+6);
   (*C)[2] = (double *) (C+9); 
}

void d_mat3dr_stack_clear()
{
   flint_stack_release();
}

void d_mat3dc_stack_init(d_mat3dc_t * C)
{
   (*C) = (double **) flint_stack_alloc(12);
   (*C)[0] = (double *) (C+3);
   (*C)[1] = (double *) (C+6);
   (*C)[2] = (double *) (C+9); 
}

void d_mat3dc_stack_clear()
{
   flint_stack_release();
}

/* 
   Compute the Gram-Schmidt Orthogonalisation of a set of 3 vectors b_i (given as the
   columns of the input matrix B_in). The orthogonal vectors b*_i are returned as the 
   columns of the output matrix B_out and a matrix of the Gram coefficients :
                   Q_ij = <b_i, b*_j>/<b*_j, b*_j>; 1 <= j < i <= 3
   are returned as the entries below the diagonal in the matrix Q. 
   
   Aliasing is permitted.
*/

void d_mat3dc_gram_schmidt(d_mat3dc_t Q, d_mat3dc_t B_out, d_mat3dc_t B_in)
{
   d_vec3d_set(COL(B_out, 1), COL(B_in, 1));
   MAT_C(Q, 2, 1) = d_vec3d_scalar_proj(COL(B_in, 2), COL(B_in, 1));
   MAT_C(Q, 3, 1) = d_vec3d_scalar_proj(COL(B_in, 3), COL(B_in, 1));
   
   d_vec3d_sub_scalar_mul(COL(B_out, 2), COL(B_in, 2), COL(B_in, 1), MAT_C(Q, 2, 1));
   MAT_C(Q, 3, 2) = d_vec3d_scalar_proj(COL(B_in, 3), COL(B_out, 2));
   
   d_vec3d_sub_scalar_mul(COL(B_out, 3), COL(B_in, 3), COL(B_out, 2), MAT_C(Q, 3, 2));
   d_vec3d_sub_scalar_mul(COL(B_out, 3), COL(B_out, 3), COL(B_out, 1), MAT_C(Q, 3, 1));
}  

/* 
   Performs LLL reduction on the matrix B_in with a given delta value.
   We assume 1/4 < delta <= 1 (for now delta doesn't do anything and is
   set to 1/4).
   
   B_in is a matrix whose columns are a basis specifying a subspace of R^3.
   
   The matrix C will become a set of row vectors, each of which specifies a 
   vector of the reduced basis as a linear combination of the original 
   basis vectors.
*/

void d_mat3dc_LLL(z_mat3dr_t C, d_mat3dc_t B_in, double delta)
{
   d_mat3dc_t Q;
   d_mat3dc_t B;
   double B1, B2, B3, B_temp, mu;
   
   d_mat3dc_stack_init(&Q);
   d_mat3dc_stack_init(&B);
   
   z_mat3dr_set_identity(C);
   
   d_mat3dc_gram_schmidt(Q, B, B_in);
   
   B1 = d_vec3d_norm(COL(B, 1));
   B2 = d_vec3d_norm(COL(B, 2));
   B3 = d_vec3d_norm(COL(B, 3));
   
   // k = 2
k2:
   mu = MAT_C(Q, 2, 1);
   if (abs(mu) > 0.5) 
   {
      unsigned long r = round(mu);
      z_vec3d_sub_scalar_mul(ROW(C, 2), ROW(C, 2), ROW(C, 1), r);
      MAT_C(Q, 2, 1) -= (double) r;
   }
   
   mu = MAT_C(Q, 2, 1);
   if (B2 < (0.75 - mu*mu)*B1)
   {
      B_temp = B2 + mu*B1;
      MAT_C(Q, 2, 1) = mu*B1 / B_temp; 
      B2 = B1*B2 / B_temp;
      B1 = B_temp;
      z_mat3dr_swap12(C);
      mu = MAT_C(Q, 3, 1) - mu*MAT_C(Q, 3, 2);
      MAT_C(Q, 3, 1) = MAT_C(Q, 3, 2) + mu*MAT_C(Q, 2, 1);
      MAT_C(Q, 3, 2) = mu;
      goto k2;
   }
   
   // k = 3
k3:
   mu = MAT_C(Q, 3, 2);
   if (abs(mu) > 0.5) 
   {
      unsigned long r = round(mu);
      z_vec3d_sub_scalar_mul(ROW(C, 3), ROW(C, 3), ROW(C, 2), r);
      MAT_C(Q, 3, 1) -= r*MAT_C(Q, 2, 1);
      MAT_C(Q, 3, 2) -= (double) r;
   }
   
   mu = MAT_C(Q, 3, 2);
   if (B3 >= (0.75 - mu*mu)*B2)
   {
      mu = MAT_C(Q, 3, 1);
      if (abs(mu) > 0.5)
      {
         unsigned long r = round(mu);
         z_vec3d_sub_scalar_mul(ROW(C, 3), ROW(C, 3), ROW(C, 1), r);
         MAT_C(Q, 3, 1) -= (double) r;
      }
   } else
   {
      B_temp = B3 + mu*mu*B2;
      MAT_C(Q, 3, 2) = mu*B2 / B_temp;
      B3 = B2*B3 / B_temp;
      B2 = B_temp;
      z_mat3dr_swap23(C);
      mu = MAT_C(Q, 2, 1);
      MAT_C(Q, 2, 1) = MAT_C(Q, 3, 1);
      MAT_C(Q, 3, 1) = mu;
      goto k2;
   } 
    
   // k = 4 : terminate
   d_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
} 
