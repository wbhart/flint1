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

#include "qd/c_dd.h"

#define D_COORD(A,x) ((A) + 2*(x-1))
#define D_COL(A, y) ((A)[y-1])

int main(void)
{
   double a1[2];
   double scalar[2];
   dd_vec3d avec;
   dd_vec3d avec2;
   dd_vec3d avec3;
   dd_mat3dc_t mat, mat2, mat4;
   z_mat3dc_t mat3;
   
   char str[64];
   
   dd_vec3d_stack_init(&avec);
   dd_vec3d_stack_init(&avec2);
   dd_vec3d_stack_init(&avec3);
   dd_mat3dc_stack_init(&mat);
   dd_mat3dc_stack_init(&mat2);
   dd_mat3dc_stack_init(&mat4);
   z_mat3dc_stack_init(&mat3);
   
   printf("Input a 3d vector A, e.g. 1.20 3.45 -4.55 : ");
   dd_vec3d_scanf(avec);
   printf("You entered A = "); dd_vec3d_printf(avec); printf("\n");
   dd_vec3d_norm(a1, avec);
   c_dd_swrite(a1, str);
   printf("The norm of A is %s\n", str);
   dd_vec3d_length(a1, avec);
   c_dd_swrite(a1, str);
   printf("The length of A is %s\n", str);
   dd_vec3d_normalise(avec2, avec);
   printf("A normalised is "); dd_vec3d_printf(avec2); printf("\n");
   printf("Enter a scalar : ");
   scanf("%s", &str); getchar(); c_dd_read(str, scalar);
   dd_vec3d_mul_scalar(avec3, avec, scalar);
   printf("A multiplied by %s is ", str); dd_vec3d_printf(avec3); printf("\n");
   printf("Enter another 3d vector B : ");
   dd_vec3d_scanf(avec3);
   dd_vec3d_scalar_prod(a1, avec, avec3);
   c_dd_swrite(a1, str);
   printf("The scalar product of A and B is %s\n", str);
   dd_vec3d_vector_proj(avec2, avec, avec3);
   printf("The vector projection of A onto B is "); dd_vec3d_printf(avec2); printf("\n");
   dd_vec3d_scalar_proj(a1, avec, avec3);
   c_dd_swrite(a1, str);
   printf("The scalar projection of A onto B is %s\n", str);
   printf("Setting C equal to A....\n");
   dd_vec3d_set(avec2, avec);
   printf("Checking that "); dd_vec3d_printf(avec); printf(" equals "); dd_vec3d_printf(avec2); printf("\n");
   dd_vec3d_add_scalar_mul(avec2, avec, avec3, scalar);
   c_dd_swrite(scalar, str);
   printf("A+%s*B = ", str); dd_vec3d_printf(avec2); printf("\n"); 
   dd_vec3d_sub_scalar_mul(avec2, avec, avec3, scalar);
   printf("A-%s*B = ", str); dd_vec3d_printf(avec2); printf("\n\n"); 
   
   printf("Setting D to the identity matrix....\n");
   dd_mat3dc_set_identity(mat);
   printf("D is:"); dd_mat3dc_printf(mat);
   printf("Enter nine floating point numbers for the entries of E, a 3x3 matrix: ");
   dd_mat3dc_scanf(mat2);
   printf("E is:"); dd_mat3dc_printf(mat2);
   printf("Enter nine integers for the entries of F, a 3x3 matrix: ");
   z_mat3dc_scanf(mat3);
   printf("F is:"); z_mat3dc_printf(mat3);
   printf("Swapping rows 1 and 2 of F....\n");
   z_mat3dc_swap12(mat3);
   printf("F is now:"); z_mat3dc_printf(mat3);
   printf("The Gram-Schmidt orthogonalisation of E is: ");
   dd_mat3dc_set_zero(mat4);
   dd_mat3dc_gram_schmidt(mat4, mat, mat2);
   dd_mat3dc_printf(mat);
   printf("The Gram coefficients of E are: ");
   dd_mat3dc_printf(mat4);
   dd_mat3dc_det(a1, mat2);
   c_dd_swrite(a1, str);
   printf("The determinant of E is: %s\n", str);
   printf("The inverse of E is:"); 
   dd_mat3dc_invert(mat, mat2);
   dd_mat3dc_printf(mat);
   
   printf("E is:"); dd_mat3dc_printf(mat2);
   printf("The LLL coefficients of E are: ");
   d_mat3dc_LLL(mat3, mat, mat2, 0.99);
   z_mat3dc_printf(mat3);
   printf("The reduced matrix is: ");
   dd_mat3dc_printf(mat);
   dd_vec3d_length(a1, D_COL(mat, 1));
   c_dd_swrite(a1, str);
   printf("The lengths are %s, ", str);
   dd_vec3d_length(a1, D_COL(mat, 2));
   c_dd_swrite(a1, str);
   printf("%s, ", str);
   dd_vec3d_length(a1, D_COL(mat, 3));
   c_dd_swrite(a1, str);
   printf("%s\n", str);
   dd_vec3d_scalar_prod(a1, D_COL(mat,1), D_COL(mat,2));
   c_dd_swrite(a1, str);
   printf("The scalar product of C1 and C2 is %s\n", str);
   dd_vec3d_scalar_prod(a1, D_COL(mat,1), D_COL(mat,3));
   c_dd_swrite(a1, str);
   printf("The scalar product of C1 and C3 is %s\n", str);
   dd_vec3d_scalar_prod(a1, D_COL(mat,2), D_COL(mat,3));
   c_dd_swrite(a1, str);
   printf("The scalar product of C2 and C3 is %s\n", str);

   z_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
   dd_vec3d_stack_clear();
   dd_vec3d_stack_clear();
   dd_vec3d_stack_clear();
}
