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

#include "qd/c_dd.h"

#define ROW(R, x) ((R)[x-1])
#define COL(R, x) ((R)[x-1])
#define MAT_R(R, x, y) ((R)[x-1][y-1])
#define MAT_C(R, x, y) ((R)[y-1][x-1])

/************************************************************************************

   long int matrices
   
************************************************************************************/

void z_mat3dr_stack_init(z_mat3dr_t * C)
{
   (*C) = (long **) flint_stack_alloc(12);
   (*C)[0] = (long *) (*C+3);
   (*C)[1] = (long *) (*C+6);
   (*C)[2] = (long *) (*C+9); 
}

void z_mat3dr_stack_clear(void)
{
   flint_stack_release();
}

void z_mat3dc_stack_init(z_mat3dc_t * C)
{
   (*C) = (long **) flint_stack_alloc(12);
   (*C)[0] = (long *) (*C+3);
   (*C)[1] = (long *) (*C+6);
   (*C)[2] = (long *) (*C+9); 
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

void z_mat3dc_swap12(z_mat3dc_t C)
{
   z_vec3d temp;
   temp = COL(C, 1);
   COL(C, 1) = COL(C, 2);
   COL(C, 2) = temp;
} 

void z_mat3dc_swap13(z_mat3dc_t C)
{
   z_vec3d temp;
   temp = COL(C, 1);
   COL(C, 1) = ROW(C, 3);
   COL(C, 3) = temp;
} 

void z_mat3dc_swap23(z_mat3dc_t C)
{
   z_vec3d temp;
   temp = COL(C, 2);
   COL(C, 2) = COL(C, 3);
   COL(C, 3) = temp;
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

void z_mat3dc_set_identity(z_mat3dc_t C)
{
   MAT_C(C, 1, 1) = 1;
   MAT_C(C, 1, 2) = 0;
   MAT_C(C, 1, 3) = 0;
   MAT_C(C, 2, 1) = 0;
   MAT_C(C, 2, 2) = 1;
   MAT_C(C, 2, 3) = 0;
   MAT_C(C, 3, 1) = 0;
   MAT_C(C, 3, 2) = 0;
   MAT_C(C, 3, 3) = 1;
}

void z_mat3dc_printf(z_mat3dc_t mat)
{
   printf("\n[%ld %ld %ld]\n", MAT_C(mat, 1, 1), MAT_C(mat, 1, 2), MAT_C(mat, 1, 3));
   printf("[%ld %ld %ld]\n", MAT_C(mat, 2, 1), MAT_C(mat, 2, 2), MAT_C(mat, 2, 3));
   printf("[%ld %ld %ld]\n", MAT_C(mat, 3, 1), MAT_C(mat, 3, 2), MAT_C(mat, 3, 3)); 
}

void z_mat3dc_scanf(z_mat3dc_t mat)
{
   scanf("%ld %ld %ld", &MAT_C(mat, 1, 1), &MAT_C(mat, 1, 2), &MAT_C(mat, 1, 3));
   scanf("%ld %ld %ld", &MAT_C(mat, 2, 1), &MAT_C(mat, 2, 2), &MAT_C(mat, 2, 3));
   scanf("%ld %ld %ld", &MAT_C(mat, 3, 1), &MAT_C(mat, 3, 2), &MAT_C(mat, 3, 3)); getchar();
}

void z_mat3dr_printf(z_mat3dr_t mat)
{
   printf("\n[%ld %ld %ld]\n", MAT_R(mat, 1, 1), MAT_R(mat, 1, 2), MAT_R(mat, 1, 3));
   printf("[%ld %ld %ld]\n", MAT_R(mat, 2, 1), MAT_R(mat, 2, 2), MAT_R(mat, 2, 3));
   printf("[%ld %ld %ld]\n", MAT_R(mat, 3, 1), MAT_R(mat, 3, 2), MAT_R(mat, 3, 3)); 
}

void z_mat3dr_scanf(z_mat3dr_t mat)
{
   scanf("%ld %ld %ld", &MAT_R(mat, 1, 1), &MAT_R(mat, 1, 2), &MAT_R(mat, 1, 3));
   scanf("%ld %ld %ld", &MAT_R(mat, 2, 1), &MAT_R(mat, 2, 2), &MAT_R(mat, 2, 3));
   scanf("%ld %ld %ld", &MAT_R(mat, 3, 1), &MAT_R(mat, 3, 2), &MAT_R(mat, 3, 3)); getchar();
}

/************************************************************************************

   Double precision matrices
   
************************************************************************************/

void d_mat3dc_set_identity(d_mat3dc_t C)
{
   MAT_C(C, 1, 1) = 1.0;
   MAT_C(C, 1, 2) = 0.0;
   MAT_C(C, 1, 3) = 0.0;
   MAT_C(C, 2, 1) = 0.0;
   MAT_C(C, 2, 2) = 1.0;
   MAT_C(C, 2, 3) = 0.0;
   MAT_C(C, 3, 1) = 0.0;
   MAT_C(C, 3, 2) = 0.0;
   MAT_C(C, 3, 3) = 1.0;
}

void d_mat3dc_set_zero(d_mat3dc_t C)
{
   MAT_C(C, 1, 1) = 0.0;
   MAT_C(C, 1, 2) = 0.0;
   MAT_C(C, 1, 3) = 0.0;
   MAT_C(C, 2, 1) = 0.0;
   MAT_C(C, 2, 2) = 0.0;
   MAT_C(C, 2, 3) = 0.0;
   MAT_C(C, 3, 1) = 0.0;
   MAT_C(C, 3, 2) = 0.0;
   MAT_C(C, 3, 3) = 0.0;
}

void d_mat3dr_stack_init(d_mat3dr_t * C)
{
   (*C) = (double **) flint_stack_alloc(12);
   (*C)[0] = (double *) (*C+3);
   (*C)[1] = (double *) (*C+6);
   (*C)[2] = (double *) (*C+6); 
}

void d_mat3dr_stack_clear()
{
   flint_stack_release();
}

void d_mat3dc_stack_init(d_mat3dc_t * C)
{
   (*C) = (double **) flint_stack_alloc(12);
   (*C)[0] = (double *) (*C+3);
   (*C)[1] = (double *) (*C+6);
   (*C)[2] = (double *) (*C+9); 
}

void d_mat3dc_stack_clear()
{
   flint_stack_release();
}

void d_mat3dc_printf(d_mat3dc_t mat)
{
   printf("\n[%.18f %.18f %.18f]\n", MAT_C(mat, 1, 1), MAT_C(mat, 1, 2), MAT_C(mat, 1, 3));
   printf("[%.18f %.18f %.18f]\n", MAT_C(mat, 2, 1), MAT_C(mat, 2, 2), MAT_C(mat, 2, 3));
   printf("[%.18f %.18f %.18f]\n", MAT_C(mat, 3, 1), MAT_C(mat, 3, 2), MAT_C(mat, 3, 3)); 
}

void d_mat3dc_scanf(d_mat3dc_t mat)
{
   scanf("%lf %lf %lf", &MAT_C(mat, 1, 1), &MAT_C(mat, 1, 2), &MAT_C(mat, 1, 3));
   scanf("%lf %lf %lf", &MAT_C(mat, 2, 1), &MAT_C(mat, 2, 2), &MAT_C(mat, 2, 3));
   scanf("%lf %lf %lf", &MAT_C(mat, 3, 1), &MAT_C(mat, 3, 2), &MAT_C(mat, 3, 3)); getchar();
}

double d_mat3dc_det(d_mat3dc_t B)
{
   return MAT_C(B, 1, 1)*(MAT_C(B, 3, 3)*MAT_C(B, 2, 2) - MAT_C(B, 3, 2)*MAT_C(B, 2, 3))
        - MAT_C(B, 2, 1)*(MAT_C(B, 3, 3)*MAT_C(B, 1, 2) - MAT_C(B, 3, 2)*MAT_C(B, 1, 3))
        + MAT_C(B, 3, 1)*(MAT_C(B, 2, 3)*MAT_C(B, 1, 2) - MAT_C(B, 2, 2)*MAT_C(B, 1, 3));
   
}

void d_mat3dc_invert(d_mat3dc_t B_inv, d_mat3dc_t B)
{
   double det = d_mat3dc_det(B);
   
   MAT_C(B_inv, 1, 1) = (MAT_C(B, 3, 3)*MAT_C(B, 2, 2) - MAT_C(B, 3, 2)*MAT_C(B, 2, 3))/det;
   MAT_C(B_inv, 1, 2) = -(MAT_C(B, 3, 3)*MAT_C(B, 1, 2) - MAT_C(B, 3, 2)*MAT_C(B, 1, 3))/det;
   MAT_C(B_inv, 1, 3) = (MAT_C(B, 2, 3)*MAT_C(B, 1, 2) - MAT_C(B, 2, 2)*MAT_C(B, 1, 3))/det;
   MAT_C(B_inv, 2, 1) = -(MAT_C(B, 3, 3)*MAT_C(B, 2, 1) - MAT_C(B, 3, 1)*MAT_C(B, 2, 3))/det;
   MAT_C(B_inv, 2, 2) = (MAT_C(B, 3, 3)*MAT_C(B, 1, 1) - MAT_C(B, 3, 1)*MAT_C(B, 1, 3))/det;
   MAT_C(B_inv, 2, 3) = -(MAT_C(B, 2, 3)*MAT_C(B, 1, 1) - MAT_C(B, 2, 1)*MAT_C(B, 1, 3))/det;
   MAT_C(B_inv, 3, 1) = (MAT_C(B, 3, 2)*MAT_C(B, 2, 1) - MAT_C(B, 3, 1)*MAT_C(B, 2, 2))/det;
   MAT_C(B_inv, 3, 2) = -(MAT_C(B, 3, 2)*MAT_C(B, 1, 1) - MAT_C(B, 3, 1)*MAT_C(B, 1, 2))/det;
   MAT_C(B_inv, 3, 3) = (MAT_C(B, 2, 2)*MAT_C(B, 1, 1) - MAT_C(B, 2, 1)*MAT_C(B, 1, 2))/det;
}

void d_mat3dc_mul_z_mat3dc(d_mat3dc_t B_out, d_mat3dc_t B, z_mat3dc_t Z)
{
   MAT_C(B_out, 1, 1) = MAT_C(B, 1, 1)*MAT_C(Z, 1, 1) + MAT_C(B, 1, 2)*MAT_C(Z, 2, 1) + MAT_C(B, 1, 3)*MAT_C(Z, 3, 1);
   MAT_C(B_out, 1, 2) = MAT_C(B, 1, 1)*MAT_C(Z, 1, 2) + MAT_C(B, 1, 2)*MAT_C(Z, 2, 2) + MAT_C(B, 1, 3)*MAT_C(Z, 3, 2);
   MAT_C(B_out, 1, 3) = MAT_C(B, 1, 1)*MAT_C(Z, 1, 3) + MAT_C(B, 1, 2)*MAT_C(Z, 2, 3) + MAT_C(B, 1, 3)*MAT_C(Z, 3, 3);
   MAT_C(B_out, 2, 1) = MAT_C(B, 2, 1)*MAT_C(Z, 1, 1) + MAT_C(B, 2, 2)*MAT_C(Z, 2, 1) + MAT_C(B, 2, 3)*MAT_C(Z, 3, 1);
   MAT_C(B_out, 2, 2) = MAT_C(B, 2, 1)*MAT_C(Z, 1, 2) + MAT_C(B, 2, 2)*MAT_C(Z, 2, 2) + MAT_C(B, 2, 3)*MAT_C(Z, 3, 2);
   MAT_C(B_out, 2, 3) = MAT_C(B, 2, 1)*MAT_C(Z, 1, 3) + MAT_C(B, 2, 2)*MAT_C(Z, 2, 3) + MAT_C(B, 2, 3)*MAT_C(Z, 3, 3);
   MAT_C(B_out, 3, 1) = MAT_C(B, 3, 1)*MAT_C(Z, 1, 1) + MAT_C(B, 3, 2)*MAT_C(Z, 2, 1) + MAT_C(B, 3, 3)*MAT_C(Z, 3, 1);
   MAT_C(B_out, 3, 2) = MAT_C(B, 3, 1)*MAT_C(Z, 1, 2) + MAT_C(B, 3, 2)*MAT_C(Z, 2, 2) + MAT_C(B, 3, 3)*MAT_C(Z, 3, 2);
   MAT_C(B_out, 3, 3) = MAT_C(B, 3, 1)*MAT_C(Z, 1, 3) + MAT_C(B, 3, 2)*MAT_C(Z, 2, 3) + MAT_C(B, 3, 3)*MAT_C(Z, 3, 3);
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

int d_mat3dc_LLL(z_mat3dc_t C, d_mat3dc_t B_out, d_mat3dc_t B_in, double delta)
{
   d_mat3dc_t Q;
   d_mat3dc_t B;
   double B1, B2, B3, B_temp, mu;
   unsigned long i = 0;
   double alpha = 0.500001;
   int oldexp, newexp;
   int loss = 0;
   
   d_mat3dc_stack_init(&Q);
   d_mat3dc_stack_init(&B);
   
   z_mat3dc_set_identity(C);
   
   d_mat3dc_gram_schmidt(Q, B, B_in);
   
   B1 = d_vec3d_norm(COL(B, 1));
   B2 = d_vec3d_norm(COL(B, 2));
   B3 = d_vec3d_norm(COL(B, 3));
   
   // k = 2
k2:
   mu = MAT_C(Q, 2, 1);
   if (fabsl(mu) > alpha) 
   {
      long r = (long) round(mu);
      z_vec3d_sub_scalar_mul(COL(C, 2), COL(C, 2), COL(C, 1), r);
      //frexpl(MAT_C(Q, 2, 1), &oldexp);
      MAT_C(Q, 2, 1) -= r;
      //frexpl(MAT_C(Q, 2, 1), &newexp);
      //if (oldexp > newexp) loss += (oldexp - newexp);
      /*if (loss > 15) 
      {
         d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         d_mat3dc_gram_schmidt(Q, B, B_out);
   
         B1 = d_vec3d_norm(COL(B, 1));
         B2 = d_vec3d_norm(COL(B, 2));
         B3 = d_vec3d_norm(COL(B, 3));
         loss = 0;
      }*/
   }
   
   mu = MAT_C(Q, 2, 1);
   if (B2 < (delta - mu*mu)*B1)
   {
      z_mat3dc_swap12(C);
      //frexp(B2, &oldexp);
      B_temp = B2 + mu*mu*B1;
      MAT_C(Q, 2, 1) = mu*B1 / B_temp; 
      B2 = B1*B2 / B_temp;
      B1 = B_temp;
      //frexp(B1, &newexp);
      //if (oldexp > newexp) loss += (oldexp - newexp);
      mu = MAT_C(Q, 3, 1) - mu*MAT_C(Q, 3, 2);
      //frexp(MAT_C(Q, 3, 1), &oldexp);
      MAT_C(Q, 3, 1) = MAT_C(Q, 3, 2) + mu*MAT_C(Q, 2, 1); //6-10 lost occasionally
      //frexp(MAT_C(Q, 3, 1), &newexp);
      //if (oldexp > newexp) loss += (oldexp - newexp);
      /*if (loss > 15) 
      {
         d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         d_mat3dc_gram_schmidt(Q, B, B_out);
   
         B1 = d_vec3d_norm(COL(B, 1));
         B2 = d_vec3d_norm(COL(B, 2));
         B3 = d_vec3d_norm(COL(B, 3));
         loss = 0;
      } else */
      MAT_C(Q, 3, 2) = mu;
      i++;
      if (i > 5)
      {
         d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         d_mat3dc_gram_schmidt(Q, B, B_out);
   
         B1 = d_vec3d_norm(COL(B, 1));
         B2 = d_vec3d_norm(COL(B, 2));
         B3 = d_vec3d_norm(COL(B, 3));
         alpha += 0.001;
         if (alpha > 0.6) 
         {
            d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
            d_mat3dc_stack_clear();
            d_mat3dc_stack_clear();
            return 0;
         }
         i=0;
      }
      goto k2;
   }
   
   // k = 3
k3:
   mu = MAT_C(Q, 3, 2);
   if (fabsl(mu) > alpha) 
   {
      long r = (long) round(mu);
      z_vec3d_sub_scalar_mul(COL(C, 3), COL(C, 3), COL(C, 2), r);
      //frexp(MAT_C(Q, 3, 1), &oldexp);
      MAT_C(Q, 3, 1) -= r*MAT_C(Q, 2, 1); 
      //frexp(MAT_C(Q, 3, 1), &newexp);
      //if (oldexp > newexp) loss += (oldexp - newexp);
      //frexp(MAT_C(Q, 3, 2), &oldexp);
      MAT_C(Q, 3, 2) -= r; // 20 - 21 bits lost often
      //frexp(MAT_C(Q, 3, 2), &newexp);
      //if (oldexp > newexp) loss += (oldexp - newexp);
      /*if (loss > 15) 
      {
         d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         d_mat3dc_gram_schmidt(Q, B, B_out);
   
         B1 = d_vec3d_norm(COL(B, 1));
         B2 = d_vec3d_norm(COL(B, 2));
         B3 = d_vec3d_norm(COL(B, 3));
         loss = 0;
      }*/
   }
   
   mu = MAT_C(Q, 3, 2);
   if (B3 >= (delta - mu*mu)*B2)
   {
      mu = MAT_C(Q, 3, 1);
      if (fabs(mu) > alpha)
      {
         long r = (long) round(mu);
         z_vec3d_sub_scalar_mul(COL(C, 3), COL(C, 3), COL(C, 1), r);
         //frexp(MAT_C(Q, 3, 1), &oldexp);
         MAT_C(Q, 3, 1) -= r; // 6-7 bits lost occasionally
         //frexp(MAT_C(Q, 3, 1), &newexp);
         //if (oldexp > newexp) loss += (oldexp - newexp);
         /*if (loss > 15) 
         {
            d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
            d_mat3dc_gram_schmidt(Q, B, B_out);
   
            B1 = d_vec3d_norm(COL(B, 1));
            B2 = d_vec3d_norm(COL(B, 2));
            B3 = d_vec3d_norm(COL(B, 3));
            loss = 0;
         }*/
      }
   } else
   {
      z_mat3dc_swap23(C);  
      B_temp = B3 + mu*mu*B2;
      MAT_C(Q, 3, 2) = mu*B2 / B_temp;
      B3 = B2*B3 / B_temp;
      B2 = B_temp;
      mu = MAT_C(Q, 2, 1);
      MAT_C(Q, 2, 1) = MAT_C(Q, 3, 1);
      MAT_C(Q, 3, 1) = mu;
      i++;
      if (i > 5)
      {
         d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         d_mat3dc_gram_schmidt(Q, B, B_out);
   
         B1 = d_vec3d_norm(COL(B, 1));
         B2 = d_vec3d_norm(COL(B, 2));
         B3 = d_vec3d_norm(COL(B, 3));
         alpha += 0.001;
         if (alpha > 0.6) 
         {
            d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
            d_mat3dc_stack_clear();
            d_mat3dc_stack_clear();
            return 0;
         }
         i=0;
      }
      
      goto k2;
   } 
    
   // k = 4 : terminate
   d_mat3dc_mul_z_mat3dc(B_out, B_in, C);
   d_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
   return 1;
} 

/************************************************************************************

   Double-double precision matrices
   
************************************************************************************/

#define D_COORD(A, x) ((A) + 2*(x-1))
#define D_COL(A, y) ((A)[y-1])
#define D_ROW(A, x) ((A)[x-1])
#define D_MAT_R(R, x, y) D_COORD(((R)[x-1]),y)
#define D_MAT_C(R, x, y) D_COORD(((R)[y-1]),x)

void dd_mat3dr_stack_init(dd_mat3dr_t * C)
{
   (*C) = (double **) flint_stack_alloc(21);
   (*C)[0] = (double *) (*C+3);
   (*C)[1] = (double *) (*C+9);
   (*C)[2] = (double *) (*C+15); 
}

void dd_mat3dr_stack_clear()
{
   flint_stack_release();
}

void dd_mat3dc_stack_init(dd_mat3dc_t * C)
{
   (*C) = (double **) flint_stack_alloc(21);
   (*C)[0] = (double *) (*C+3);
   (*C)[1] = (double *) (*C+9);
   (*C)[2] = (double *) (*C+15); 
}

void dd_mat3dc_stack_clear()
{
   flint_stack_release();
}

void dd_mat3dc_printf(dd_mat3dc_t mat)
{
   char d_str1[64];
   char d_str2[64];
   char d_str3[64];

   c_dd_swrite(D_MAT_C(mat, 1, 1), d_str1);
   c_dd_swrite(D_MAT_C(mat, 1, 2), d_str2);
   c_dd_swrite(D_MAT_C(mat, 1, 3), d_str3);
   printf("\n[%s, %s, %s]\n", d_str1, d_str2, d_str3); 
   c_dd_swrite(D_MAT_C(mat, 2, 1), d_str1);
   c_dd_swrite(D_MAT_C(mat, 2, 2), d_str2);
   c_dd_swrite(D_MAT_C(mat, 2, 3), d_str3);
   printf("\n[%s, %s, %s]\n", d_str1, d_str2, d_str3); 
   c_dd_swrite(D_MAT_C(mat, 3, 1), d_str1);
   c_dd_swrite(D_MAT_C(mat, 3, 2), d_str2);
   c_dd_swrite(D_MAT_C(mat, 3, 3), d_str3);
   printf("\n[%s, %s, %s]\n", d_str1, d_str2, d_str3); 
}

void dd_mat3dc_scanf(dd_mat3dc_t mat)
{
   char str1[64];
   char str2[64];
   char str3[64];
   
   scanf("%s %s %s", str1, str2, str3);
   c_dd_read(str1, D_MAT_C(mat, 1, 1));
   c_dd_read(str2, D_MAT_C(mat, 1, 2));
   c_dd_read(str3, D_MAT_C(mat, 1, 3));
   scanf("%s %s %s", str1, str2, str3);
   c_dd_read(str1, D_MAT_C(mat, 2, 1));
   c_dd_read(str2, D_MAT_C(mat, 2, 2));
   c_dd_read(str3, D_MAT_C(mat, 2, 3));
   scanf("%s %s %s", str1, str2, str3); getchar();
   c_dd_read(str1, D_MAT_C(mat, 3, 1));
   c_dd_read(str2, D_MAT_C(mat, 3, 2));
   c_dd_read(str3, D_MAT_C(mat, 3, 3));
}


void dd_mat3dc_set_identity(dd_mat3dc_t C)
{
   c_dd_copy_d(1.0, D_MAT_C(C, 1, 1));
   c_dd_copy_d(0.0, D_MAT_C(C, 1, 2));
   c_dd_copy_d(0.0, D_MAT_C(C, 1, 3));
   c_dd_copy_d(0.0, D_MAT_C(C, 2, 1));
   c_dd_copy_d(1.0, D_MAT_C(C, 2, 2));
   c_dd_copy_d(0.0, D_MAT_C(C, 2, 3));
   c_dd_copy_d(0.0, D_MAT_C(C, 3, 1));
   c_dd_copy_d(0.0, D_MAT_C(C, 3, 2));
   c_dd_copy_d(1.0, D_MAT_C(C, 3, 3));
}

void dd_mat3dc_set_zero(dd_mat3dc_t C)
{
   c_dd_copy_d(0.0, D_MAT_C(C, 1, 1));
   c_dd_copy_d(0.0, D_MAT_C(C, 1, 2));
   c_dd_copy_d(0.0, D_MAT_C(C, 1, 3));
   c_dd_copy_d(0.0, D_MAT_C(C, 2, 1));
   c_dd_copy_d(0.0, D_MAT_C(C, 2, 2));
   c_dd_copy_d(0.0, D_MAT_C(C, 2, 3));
   c_dd_copy_d(0.0, D_MAT_C(C, 3, 1));
   c_dd_copy_d(0.0, D_MAT_C(C, 3, 2));
   c_dd_copy_d(0.0, D_MAT_C(C, 3, 3));
}

/* 
   Compute the Gram-Schmidt Orthogonalisation of a set of 3 vectors b_i (given as the
   columns of the input matrix B_in). The orthogonal vectors b*_i are returned as the 
   columns of the output matrix B_out and a matrix of the Gram coefficients :
                   Q_ij = <b_i, b*_j>/<b*_j, b*_j>; 1 <= j < i <= 3
   are returned as the entries below the diagonal in the matrix Q. 
   
   Aliasing is permitted.
*/

void dd_mat3dc_gram_schmidt(dd_mat3dc_t Q, dd_mat3dc_t B_out, dd_mat3dc_t B_in)
{
   dd_vec3d_set(D_COL(B_out, 1), D_COL(B_in, 1));
   dd_vec3d_scalar_proj(D_MAT_C(Q, 2, 1), D_COL(B_in, 2), D_COL(B_in, 1));
   dd_vec3d_scalar_proj(D_MAT_C(Q, 3, 1), D_COL(B_in, 3), D_COL(B_in, 1));
   
   dd_vec3d_sub_scalar_mul(D_COL(B_out, 2), D_COL(B_in, 2), D_COL(B_in, 1), D_MAT_C(Q, 2, 1));
   dd_vec3d_scalar_proj(D_MAT_C(Q, 3, 2), D_COL(B_in, 3), D_COL(B_out, 2));
   
   dd_vec3d_sub_scalar_mul(D_COL(B_out, 3), D_COL(B_in, 3), D_COL(B_out, 2), D_MAT_C(Q, 3, 2));
   dd_vec3d_sub_scalar_mul(D_COL(B_out, 3), D_COL(B_out, 3), D_COL(B_out, 1), D_MAT_C(Q, 3, 1));
}  

void dd_mat3dc_det(double * out, dd_mat3dc_t B)
{
   double t1[2];
   double t2[2];
   double t3[2];
   
   c_dd_mul(D_MAT_C(B, 3, 3), D_MAT_C(B, 2, 2), t1);
   c_dd_mul(D_MAT_C(B, 3, 2), D_MAT_C(B, 2, 3), t2);
   c_dd_sub(t1, t2, out);
   c_dd_mul(out, D_MAT_C(B, 1, 1), out);
   
   c_dd_mul(D_MAT_C(B, 3, 3), D_MAT_C(B, 1, 2), t1);
   c_dd_mul(D_MAT_C(B, 3, 2), D_MAT_C(B, 1, 3), t2);
   c_dd_sub(t1, t2, t3);
   c_dd_mul(t3, D_MAT_C(B, 2, 1), t3);
   
   c_dd_sub(out, t3, out);
   
   c_dd_mul(D_MAT_C(B, 2, 3), D_MAT_C(B, 1, 2), t1);
   c_dd_mul(D_MAT_C(B, 2, 2), D_MAT_C(B, 1, 3), t2);
   c_dd_sub(t1, t2, t3);
   c_dd_mul(t3, D_MAT_C(B, 3, 1), t3);
   
   c_dd_add(out, t3, out);
}

void dd_mat3dc_invert(dd_mat3dc_t B_inv, dd_mat3dc_t B)
{
   double t1[2];
   double t2[2];
   double det[2];
   dd_mat3dc_det(det, B);
   
   c_dd_mul(D_MAT_C(B, 3, 3), D_MAT_C(B, 3, 2), t1);
   c_dd_mul(D_MAT_C(B, 3, 2), D_MAT_C(B, 2, 3), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 1, 1));
   
   c_dd_mul(D_MAT_C(B, 3, 2), D_MAT_C(B, 1, 3), t1);
   c_dd_mul(D_MAT_C(B, 3, 3), D_MAT_C(B, 1, 2), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 1, 2));
   
   c_dd_mul(D_MAT_C(B, 2, 3), D_MAT_C(B, 1, 2), t1);
   c_dd_mul(D_MAT_C(B, 2, 2), D_MAT_C(B, 1, 3), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 1, 3));
   
   c_dd_mul(D_MAT_C(B, 3, 1), D_MAT_C(B, 2, 3), t1);
   c_dd_mul(D_MAT_C(B, 3, 3), D_MAT_C(B, 2, 1), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 2, 1));
   
   c_dd_mul(D_MAT_C(B, 3, 3), D_MAT_C(B, 1, 1), t1);
   c_dd_mul(D_MAT_C(B, 3, 1), D_MAT_C(B, 1, 3), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 2, 2));
   
   c_dd_mul(D_MAT_C(B, 2, 1), D_MAT_C(B, 1, 3), t1);
   c_dd_mul(D_MAT_C(B, 2, 3), D_MAT_C(B, 1, 1), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 2, 3));
   
   c_dd_mul(D_MAT_C(B, 3, 2), D_MAT_C(B, 2, 1), t1);
   c_dd_mul(D_MAT_C(B, 3, 1), D_MAT_C(B, 2, 2), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 3, 1));
   
   c_dd_mul(D_MAT_C(B, 3, 1), D_MAT_C(B, 1, 2), t1);
   c_dd_mul(D_MAT_C(B, 3, 2), D_MAT_C(B, 1, 1), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 3, 2));
   
   c_dd_mul(D_MAT_C(B, 2, 2), D_MAT_C(B, 1, 1), t1);
   c_dd_mul(D_MAT_C(B, 2, 1), D_MAT_C(B, 1, 2), t2);
   c_dd_sub(t1, t2, t1);
   c_dd_div(t1, det, D_MAT_C(B_inv, 3, 3));
}

void dd_mat3dc_mul_z_mat3dc(dd_mat3dc_t B_out, dd_mat3dc_t B, z_mat3dc_t Z)
{
   double t1[2];
   double t2[2];
   
   c_dd_mul_dd_d(D_MAT_C(B, 1, 1), MAT_C(Z, 1, 1), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 1, 2), MAT_C(Z, 2, 1), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 1, 3), MAT_C(Z, 3, 1), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 1, 1));
   
   c_dd_mul_dd_d(D_MAT_C(B, 1, 1), MAT_C(Z, 1, 2), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 1, 2), MAT_C(Z, 2, 2), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 1, 3), MAT_C(Z, 3, 2), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 1, 2));
   
   c_dd_mul_dd_d(D_MAT_C(B, 1, 1), MAT_C(Z, 1, 3), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 1, 2), MAT_C(Z, 2, 3), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 1, 3), MAT_C(Z, 3, 3), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 1, 3));
   
   c_dd_mul_dd_d(D_MAT_C(B, 2, 1), MAT_C(Z, 1, 1), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 2, 2), MAT_C(Z, 2, 1), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 2, 3), MAT_C(Z, 3, 1), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 2, 1));
   
   c_dd_mul_dd_d(D_MAT_C(B, 2, 1), MAT_C(Z, 1, 2), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 2, 2), MAT_C(Z, 2, 2), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 2, 3), MAT_C(Z, 3, 2), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 2, 2));
   
   c_dd_mul_dd_d(D_MAT_C(B, 2, 1), MAT_C(Z, 1, 3), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 2, 2), MAT_C(Z, 2, 3), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 2, 3), MAT_C(Z, 3, 3), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 2, 3));
   
   c_dd_mul_dd_d(D_MAT_C(B, 3, 1), MAT_C(Z, 1, 1), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 3, 2), MAT_C(Z, 2, 1), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 3, 3), MAT_C(Z, 3, 1), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 3, 1));
   
   c_dd_mul_dd_d(D_MAT_C(B, 3, 1), MAT_C(Z, 1, 2), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 3, 2), MAT_C(Z, 2, 2), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 3, 3), MAT_C(Z, 3, 2), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 3, 2));
   
   c_dd_mul_dd_d(D_MAT_C(B, 3, 1), MAT_C(Z, 1, 3), t1);
   c_dd_mul_dd_d(D_MAT_C(B, 3, 2), MAT_C(Z, 2, 3), t2);
   c_dd_add(t1, t2, t1);
   c_dd_mul_dd_d(D_MAT_C(B, 3, 3), MAT_C(Z, 3, 3), t2);
   c_dd_add(t1, t2, D_MAT_C(B_out, 3, 3));
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

int dd_mat3dc_LLL(z_mat3dc_t C, dd_mat3dc_t B_out, dd_mat3dc_t B_in, double del)
{
   dd_mat3dc_t Q;
   dd_mat3dc_t B;
   double B1[2];
   double B2[2];
   double B3[2];
   double B_temp[2];
   double mu[2];
   double temp[2];
   double temp2[2];
   double delta[2];
   
   int test;
   
   unsigned long i = 0;
   double alpha = 0.500001;
   c_dd_copy_d(del, delta);
   int oldexp, newexp;
   int loss = 0;
   
   dd_mat3dc_stack_init(&Q);
   dd_mat3dc_stack_init(&B);
   
   z_mat3dc_set_identity(C);
   
   dd_mat3dc_gram_schmidt(Q, B, B_in);
   
   dd_vec3d_norm(B1, D_COL(B, 1));
   dd_vec3d_norm(B2, D_COL(B, 2));
   dd_vec3d_norm(B3, D_COL(B, 3));
   
   // k = 2
k2:
   c_dd_copy(D_MAT_C(Q, 2, 1), mu);
   c_dd_abs(mu, temp);
   c_dd_comp_dd_d(temp, alpha, &test);
   if (test > 0) 
   {
      c_dd_nint(mu, temp);
      long r = (long) temp[0];
      z_vec3d_sub_scalar_mul(COL(C, 2), COL(C, 2), COL(C, 1), r);
      frexp(D_MAT_C(Q, 2, 1)[0], &oldexp);
      c_dd_sub_dd_d(D_MAT_C(Q, 2, 1), r, D_MAT_C(Q, 2, 1));
      frexp(D_MAT_C(Q, 2, 1)[0], &newexp);
      if (oldexp > newexp) loss += (oldexp - newexp);
      if (loss > 15) 
      {
         dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         dd_mat3dc_gram_schmidt(Q, B, B_out);
   
         dd_vec3d_norm(B1, D_COL(B, 1));
         dd_vec3d_norm(B2, D_COL(B, 2));
         dd_vec3d_norm(B3, D_COL(B, 3));
         loss = 0;
      }
   }
   
   c_dd_copy(D_MAT_C(Q, 2, 1), mu);
   c_dd_sqr(mu, temp);
   c_dd_sub(delta, temp, temp2);
   c_dd_mul(temp2, B1, temp2);
   c_dd_comp(B2, temp2, &test);
   if (test < 0)
   {
      z_mat3dc_swap12(C);
      frexp(B2[0], &oldexp);
      c_dd_mul(temp, B1, temp2);
      c_dd_add(B2, temp2, B_temp);
      c_dd_mul(mu, B1, temp);
      c_dd_div(temp, B_temp, D_MAT_C(Q, 2, 1));
      c_dd_mul(B1, B2, temp);
      c_dd_div(temp, B_temp, B2);
      c_dd_copy(B_temp, B1);
      frexp(B1[0], &newexp);
      if (oldexp > newexp) loss += (oldexp - newexp);
      c_dd_mul(mu, D_MAT_C(Q, 3, 2), temp);
      c_dd_sub(D_MAT_C(Q, 3, 1), temp, mu);
      frexp(D_MAT_C(Q, 3, 1)[0], &oldexp);
      c_dd_mul(mu, D_MAT_C(Q, 2, 1), temp);
      c_dd_add(D_MAT_C(Q, 3, 2), temp, D_MAT_C(Q, 3, 1));
      frexp(D_MAT_C(Q, 3, 1)[0], &newexp);
      if (oldexp > newexp) loss += (oldexp - newexp);
      if (loss > 15) 
      {
         dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         dd_mat3dc_gram_schmidt(Q, B, B_out);
   
         dd_vec3d_norm(B1, D_COL(B, 1));
         dd_vec3d_norm(B2, D_COL(B, 2));
         dd_vec3d_norm(B3, D_COL(B, 3));
         loss = 0;
      } else
         c_dd_copy(mu, D_MAT_C(Q, 3, 2));
      i++;
      if (i > 30)
      {
         alpha += 0.001;
         if (alpha > 0.6) 
         {
            dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
            dd_mat3dc_stack_clear();
            dd_mat3dc_stack_clear();
            return 0;
         }
         i=0;
      }
      goto k2;
   }
   
   // k = 3
k3:
   c_dd_copy(D_MAT_C(Q, 3, 2), mu);
   c_dd_abs(mu, temp);
   c_dd_comp_dd_d(temp, alpha, &test);
   if (test > 0) 
   {
      c_dd_nint(mu, temp);
      long r = (long) temp[0];
      z_vec3d_sub_scalar_mul(COL(C, 3), COL(C, 3), COL(C, 2), r);
      frexp(D_MAT_C(Q, 3, 1)[0], &oldexp);
      c_dd_mul_dd_d(D_MAT_C(Q, 2, 1), r, temp);
      c_dd_sub(D_MAT_C(Q, 3, 1), temp, D_MAT_C(Q, 3, 1));
      frexp(D_MAT_C(Q, 3, 1)[0], &newexp);
      if (oldexp > newexp) loss += (oldexp - newexp);
      frexp(D_MAT_C(Q, 3, 2)[0], &oldexp);
      c_dd_sub_dd_d(D_MAT_C(Q, 3, 2), r, D_MAT_C(Q, 3, 2));
      frexp(D_MAT_C(Q, 3, 2)[0], &newexp);
      if (oldexp > newexp) loss += (oldexp - newexp);
      if (loss > 15) 
      {
         dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
         dd_mat3dc_gram_schmidt(Q, B, B_out);
   
         dd_vec3d_norm(B1, D_COL(B, 1));
         dd_vec3d_norm(B2, D_COL(B, 2));
         dd_vec3d_norm(B3, D_COL(B, 3));
         loss = 0;
      }
   }
   
   c_dd_copy(D_MAT_C(Q, 3, 2), mu);
   c_dd_sqr(mu, temp);
   c_dd_sub(delta, temp, temp2);
   c_dd_mul(temp2, B2, temp2);
   c_dd_comp(B3, temp2, &test);
   if (test >= 0)
   {
      c_dd_copy(D_MAT_C(Q, 3, 1), mu);
      c_dd_abs(mu, temp);
      c_dd_comp_dd_d(temp, alpha, &test);
      if (test > 0) 
      {
         c_dd_nint(mu, temp);
         long r = (long) temp[0];
         z_vec3d_sub_scalar_mul(COL(C, 3), COL(C, 3), COL(C, 1), r);
         frexp(D_MAT_C(Q, 3, 1)[0], &oldexp);
         c_dd_sub_dd_d(D_MAT_C(Q, 3, 1), r, D_MAT_C(Q, 3, 1));
         frexp(D_MAT_C(Q, 3, 1)[0], &newexp);
         if (oldexp > newexp) loss += (oldexp - newexp);
         if (loss > 15) 
         {
            dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
            dd_mat3dc_gram_schmidt(Q, B, B_out);
   
            dd_vec3d_norm(B1, D_COL(B, 1));
            dd_vec3d_norm(B2, D_COL(B, 2));
            dd_vec3d_norm(B3, D_COL(B, 3));
            loss = 0;
         }
      }
   } else
   {
      z_mat3dc_swap23(C);  
      c_dd_sqr(mu, temp);
      c_dd_mul(temp, B2, temp2);
      c_dd_add(B3, temp2, B_temp);
      c_dd_mul(mu, B2, temp);
      c_dd_div(temp, B_temp, D_MAT_C(Q, 3, 2));
      c_dd_mul(B2, B3, temp);
      c_dd_div(temp, B_temp, B3);
      c_dd_copy(B_temp, B2);
      c_dd_copy(D_MAT_C(Q, 2, 1), mu);
      c_dd_copy(D_MAT_C(Q, 3, 1), D_MAT_C(Q, 2, 1));
      c_dd_copy(mu, D_MAT_C(Q, 3, 1));
      
      i++;
      if (i > 30)
      {
         alpha += 0.001;
         if (alpha > 0.6) 
         {
            dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
            dd_mat3dc_stack_clear();
            dd_mat3dc_stack_clear();
            return 0;
         }
         i=0;
      }
      
      goto k2;
   } 
    
   // k = 4 : terminate
   dd_mat3dc_mul_z_mat3dc(B_out, B_in, C);
   dd_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
   return 1;
} 
