/****************************************************************************

   vecmat.c: simple driver file to play with fucntions in vec3d and mat3d

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"
#include "matrix.h"
#include "vec3d.h"
#include "mat3d.h"

#define COORD(A,x) ((A)[x-1])
#define SQ(A,x) (((A)[x-1])*((A)[x-1]))
#define ROW(R, x) ((R)[x-1])
#define COL(R, x) ((R)[x-1])
#define MAT_R(R, x, y) ((R)[x-1][y-1])
#define MAT_C(R, x, y) ((R)[y-1][x-1])

int main(void)
{
   double a1, a2, a3, scalar;
   d_vec3d avec;
   d_vec3d avec2;
   d_vec3d avec3;
   d_mat3dc_t mat, mat2, mat4;
   z_mat3dc_t mat3;
   
   d_vec3d_stack_init(&avec);
   d_vec3d_stack_init(&avec2);
   d_vec3d_stack_init(&avec3);
   d_mat3dc_stack_init(&mat);
   d_mat3dc_stack_init(&mat2);
   d_mat3dc_stack_init(&mat4);
   z_mat3dc_stack_init(&mat3);
   
   printf("Input a 3d vector A, e.g. 1.20 3.45 -4.55 : ");
   d_vec3d_scanf(avec);
   printf("You entered A = "); d_vec3d_printf(avec); printf("\n");
   printf("The length of A is %f\n", d_vec3d_length(avec));
   printf("The norm of A is %f\n", d_vec3d_norm(avec));
   d_vec3d_normalise(avec2, avec);
   printf("A normalised is "); d_vec3d_printf(avec2); printf("\n");
   printf("Enter a scalar : ");
   scanf("%lf", &scalar); getchar();
   d_vec3d_mul_scalar(avec3, avec, scalar);
   printf("A multiplied by %f is ", scalar); d_vec3d_printf(avec3); printf("\n");
   printf("Enter another 3d vector B : ");
   d_vec3d_scanf(avec3);
   printf("The scalar product of A and B is %f\n", d_vec3d_scalar_prod(avec, avec2));
   d_vec3d_vector_proj(avec2, avec, avec3);
   printf("The vector projection of A onto B is "); d_vec3d_printf(avec2); printf("\n");
   printf("The scalar projection of A onto B is %f\n", d_vec3d_scalar_proj(avec, avec3));
   printf("Setting C equal to A....\n");
   d_vec3d_set(avec2, avec);
   printf("Checking that "); d_vec3d_printf(avec); printf(" equals "); d_vec3d_printf(avec2); printf("\n");
   d_vec3d_add_scalar_mul(avec2, avec, avec3, scalar);
   printf("A+%f*B = ", scalar);d_vec3d_printf(avec2); printf("\n"); 
   d_vec3d_sub_scalar_mul(avec2, avec, avec3, scalar);
   printf("A-%f*B = ", scalar);d_vec3d_printf(avec2); printf("\n\n"); 
   
   printf("Setting D to the identity matrix....\n");
   d_mat3dc_set_identity(mat);
   printf("D is:"); d_mat3dc_printf(mat);
   printf("Enter nine floating point numbers for the entries of E, a 3x3 matrix: ");
   d_mat3dc_scanf(mat2);
   printf("E is:"); d_mat3dc_printf(mat2);
   printf("Enter nine integers for the entries of F, a 3x3 matrix: ");
   z_mat3dc_scanf(mat3);
   printf("F is:"); z_mat3dc_printf(mat3);
   printf("Swapping rows 1 and 2 of F....\n");
   z_mat3dc_swap12(mat3);
   printf("F is now:"); z_mat3dc_printf(mat3);
   printf("The Gram-Schmidt orthogonalisation of E is: ");
   d_mat3dc_set_zero(mat4);
   d_mat3dc_gram_schmidt(mat4, mat, mat2);
   d_mat3dc_printf(mat);
   printf("The Gram coefficients of E are: ");
   d_mat3dc_printf(mat4);
   printf("The determinant of E is: %f\n", d_mat3dc_det(mat2));
   printf("The inverse of E is:"); 
   d_mat3dc_invert(mat, mat2);
   d_mat3dc_printf(mat);
   
   printf("E is:"); d_mat3dc_printf(mat2);
   printf("The LLL coefficients of E are: ");
   d_mat3dc_LLL(mat3, mat4, mat2, 0.75);
   z_mat3dc_printf(mat3);
   d_mat3dc_mul_z_mat3dc(mat, mat2, mat3);
   printf("The reduced matrix is: ");
   d_mat3dc_printf(mat);
   printf("The lengths are %f, %f, %f\n", d_vec3d_length(COL(mat, 1)), d_vec3d_length(COL(mat, 2)), d_vec3d_length(COL(mat, 3)));
   printf("The scalar product of C1 and C2 is %f\n", d_vec3d_scalar_prod(COL(mat,1), COL(mat,2)));
   printf("The scalar product of C1 and C3 is %f\n", d_vec3d_scalar_prod(COL(mat,1), COL(mat,3)));
   printf("The scalar product of C2 and C3 is %f\n", d_vec3d_scalar_prod(COL(mat,2), COL(mat,3)));
   
   z_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
   d_vec3d_stack_clear();
   d_vec3d_stack_clear();
   d_vec3d_stack_clear();
}
