/****************************************************************************

   vec3d.c: functions for dealing with 3D vectors of double floating 
                 point numbers

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vec3d.h"
#include "memory-manager.h"

#define COORD(A,x) ((A)[x-1])
#define SQ(A,x) (((A)[x-1])*((A)[x-1]))

void z_vec3d_stack_init(z_vec3d * vec)
{
   *vec = (long *) flint_stack_alloc_bytes(3*sizeof(long));
}

void d_vec3d_stack_init(d_vec3d * vec)
{
   *vec = (double *) flint_stack_alloc_bytes(3*sizeof(double));
}

void z_vec3d_stack_clear()
{
   flint_stack_release();
}

void d_vec3d_stack_clear()
{
   flint_stack_release();
}

/*
   Returns the Euclidean length of a 3d vector
*/

double d_vec3d_length(d_vec3d v)
{
   return sqrt(SQ(v,1) + SQ(v,2) + SQ(v,3));
}

/*
   Returns the square of the Euclidean length of a 3d vector
*/

double d_vec3d_norm(d_vec3d v)
{
   return SQ(v,1) + SQ(v,2) + SQ(v,3);
}

/*
   Normalises the given input vector (i.e. divides it by its length)
*/

void d_vec3d_normalise(d_vec3d v_out, d_vec3d v_in)
{
   double length = d_vec3d_length(v_in);
   
   COORD(v_out, 1) = COORD(v_in, 1) / length;
   COORD(v_out, 2) = COORD(v_in, 2) / length;
   COORD(v_out, 3) = COORD(v_in, 3) / length;
}

/* 
   Multiplies the given vector by the given scalar
*/

void d_vec3d_mul_scalar(d_vec3d v_out, d_vec3d v_in, double scalar)
{
   COORD(v_out, 1) = COORD(v_in, 1)*scalar;
   COORD(v_out, 2) = COORD(v_in, 2)*scalar;
   COORD(v_out, 3) = COORD(v_in, 3)*scalar;
}

/*
   Computes the scalar product of the two input vectors
*/

double d_vec3d_scalar_prod(d_vec3d v1, d_vec3d v2)
{
   return COORD(v1, 1)*COORD(v2, 1) + COORD(v1, 2)*COORD(v2, 2) + COORD(v1, 3)*COORD(v2, 3);
}

/*
   Computes the vector projection of _v_ on _u_
*/

void d_vec3d_vector_proj(d_vec3d v_out, d_vec3d v, d_vec3d u)
{
   double proj_scalar = d_vec3d_scalar_prod(v, u) / d_vec3d_scalar_prod(u, u);
   printf("Scalar is %lf\n", proj_scalar);
   d_vec3d_mul_scalar(v_out, u, proj_scalar);   
}

/*
   Returns the scalar projection of _v_ on _u_
*/

double d_vec3d_scalar_proj(d_vec3d v, d_vec3d u)
{
   return d_vec3d_scalar_prod(v, u) / d_vec3d_scalar_prod(u, u);
}

/* 
   Sets the vector v_out to be equal to v_in
*/

void d_vec3d_set(d_vec3d v_out, d_vec3d v_in)
{
   COORD(v_out, 1) = COORD(v_in, 1);
   COORD(v_out, 2) = COORD(v_in, 2);
   COORD(v_out, 3) = COORD(v_in, 3);
}

/*
   Sets v_out to v_in1 - scalar * v_in2
*/

void d_vec3d_sub_scalar_mul(d_vec3d v_out, d_vec3d v_in1, d_vec3d v_in2, double scalar)
{
   COORD(v_out, 1) = COORD(v_in1, 1) - scalar * COORD(v_in2, 1);
   COORD(v_out, 2) = COORD(v_in1, 2) - scalar * COORD(v_in2, 2);
   COORD(v_out, 3) = COORD(v_in1, 3) - scalar * COORD(v_in2, 3);
}

/*
   Sets v_out to v_in1 + scalar * v_in2
*/

void d_vec3d_add_scalar_mul(d_vec3d v_out, d_vec3d v_in1, d_vec3d v_in2, double scalar)
{
   COORD(v_out, 1) = COORD(v_in1, 1) + scalar * COORD(v_in2, 1);
   COORD(v_out, 2) = COORD(v_in1, 2) + scalar * COORD(v_in2, 2);
   COORD(v_out, 3) = COORD(v_in1, 3) + scalar * COORD(v_in2, 3);
}

/*
   Sets v_out to v_in1 - scalar * v_in2
*/

void z_vec3d_sub_scalar_mul(z_vec3d v_out, z_vec3d v_in1, z_vec3d v_in2, long scalar)
{
   COORD(v_out, 1) = COORD(v_in1, 1) - scalar * COORD(v_in2, 1);
   COORD(v_out, 2) = COORD(v_in1, 2) - scalar * COORD(v_in2, 2);
   COORD(v_out, 3) = COORD(v_in1, 3) - scalar * COORD(v_in2, 3);
}

/*
   Sets v_out to v_in1 + scalar * v_in2
*/

void z_vec3d_add_scalar_mul(z_vec3d v_out, z_vec3d v_in1, z_vec3d v_in2, long scalar)
{
   COORD(v_out, 1) = COORD(v_in1, 1) + scalar * COORD(v_in2, 1);
   COORD(v_out, 2) = COORD(v_in1, 2) + scalar * COORD(v_in2, 2);
   COORD(v_out, 3) = COORD(v_in1, 3) + scalar * COORD(v_in2, 3);
}

void d_vec3d_scanf(d_vec3d vec)
{
   scanf("%lf %lf %lf", &COORD(vec, 1), &COORD(vec, 2), &COORD(vec, 3)); getchar(); 
}

void d_vec3d_printf(d_vec3d vec)
{
   printf("[%f %f %f]", COORD(vec, 1), COORD(vec, 2), COORD(vec, 3)); 
}
