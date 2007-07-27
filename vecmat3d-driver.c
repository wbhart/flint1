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
   
   d_vec3d_stack_init(&avec);
   d_vec3d_stack_init(&avec2);
   d_vec3d_stack_init(&avec3);
   
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
    
   
   d_vec3d_stack_clear();
   d_vec3d_stack_clear();
   d_vec3d_stack_clear();
}
