/****************************************************************************

   vec3d.c: functions for dealing with 3D vectors of double floating 
                 point numbers

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "vec3d.h"
#include "memory-manager.h"

#include "qd/c_dd.h"

#define COORD(A,x) ((A)[x-1])
#define SQ(A,x) (((A)[x-1])*((A)[x-1]))

/************************************************************************************

   long int vectors
   
************************************************************************************/

void z_vec3d_stack_init(z_vec3d * vec)
{
   *vec = (long *) flint_stack_alloc_bytes(3*sizeof(long));
}

void z_vec3d_stack_clear()
{
   flint_stack_release();
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

/************************************************************************************

   Double precision vectors
   
************************************************************************************/

void d_vec3d_stack_init(d_vec3d * vec)
{
   *vec = (double *) flint_stack_alloc_bytes(3*sizeof(double));
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

void d_vec3d_scanf(d_vec3d vec)
{
   scanf("%lf %lf %lf", &COORD(vec, 1), &COORD(vec, 2), &COORD(vec, 3)); getchar(); 
}

void d_vec3d_printf(d_vec3d vec)
{
   printf("[%f, %f, %f]", COORD(vec, 1), COORD(vec, 2), COORD(vec, 3)); 
}

/************************************************************************************

   Double-double precision vectors
   
************************************************************************************/

#define D_COORD(A,x) ((A) + 2*(x-1))

void dd_vec3d_stack_init(dd_vec3d * vec)
{
   *vec = (double *) flint_stack_alloc_bytes(3*2*sizeof(double));
}

void dd_vec3d_stack_clear()
{
   flint_stack_release();
}

void dd_vec3d_scanf(dd_vec3d vec)
{
   char d_str1[64];
   char d_str2[64];
   char d_str3[64];
   
   scanf("%s %s %s", &d_str1, &d_str2, &d_str3); getchar(); 
   c_dd_read(d_str1, D_COORD(vec, 1));  
   c_dd_read(d_str2, D_COORD(vec, 2));
   c_dd_read(d_str3, D_COORD(vec, 3));  
}

void dd_vec3d_printf(dd_vec3d vec)
{
   char d_str1[64];
   char d_str2[64];
   char d_str3[64];
   
   c_dd_swrite(D_COORD(vec, 1), d_str1);
   c_dd_swrite(D_COORD(vec, 2), d_str2);
   c_dd_swrite(D_COORD(vec, 3), d_str3);
   
   printf("[%s, %s, %s]", d_str1, d_str2, d_str3); 
}

/*
   Returns the Euclidean norm of a 3d vector
*/

void dd_vec3d_norm(double * out, d_vec3d v)
{
   double temp[2];
   
   c_dd_sqr(D_COORD(v, 1), out);
   c_dd_sqr(D_COORD(v, 2), temp);
   c_dd_add(temp, out, out); 
   c_dd_sqr(D_COORD(v, 3), temp);
   c_dd_add(temp, out, out); 
}

/*
   Returns the Euclidean norm of a 3d vector
*/

void dd_vec3d_length(double * out, d_vec3d v)
{
   double temp[2];
   
   c_dd_sqr(D_COORD(v, 1), out);
   c_dd_sqr(D_COORD(v, 2), temp);
   c_dd_add(temp, out, out); 
   c_dd_sqr(D_COORD(v, 3), temp);
   c_dd_add(temp, out, out);
   c_dd_sqrt(out, out); 
}

/*
   Normalises the given input vector (i.e. divides it by its length)
*/

void dd_vec3d_normalise(dd_vec3d v_out, dd_vec3d v_in)
{
   double length[2];
   
   dd_vec3d_length(length, v_in);
   
   c_dd_div(D_COORD(v_in, 1), length, D_COORD(v_out, 1));
   c_dd_div(D_COORD(v_in, 2), length, D_COORD(v_out, 2));
   c_dd_div(D_COORD(v_in, 3), length, D_COORD(v_out, 3));
}

/* 
   Multiplies the given vector by the given scalar
*/

void dd_vec3d_mul_scalar(dd_vec3d v_out, dd_vec3d v_in, double * scalar)
{
   c_dd_mul(D_COORD(v_in, 1), scalar, D_COORD(v_out, 1));
   c_dd_mul(D_COORD(v_in, 2), scalar, D_COORD(v_out, 2));
   c_dd_mul(D_COORD(v_in, 3), scalar, D_COORD(v_out, 3));
}

/*
   Computes the scalar product of the two input vectors
*/

double dd_vec3d_scalar_prod(double * out, dd_vec3d v1, dd_vec3d v2)
{
   double temp[2];
   
   c_dd_mul(D_COORD(v1, 1), D_COORD(v2, 1), out);
   c_dd_mul(D_COORD(v1, 2), D_COORD(v2, 2), temp);
   c_dd_add(temp, out, out); 
   c_dd_mul(D_COORD(v1, 3), D_COORD(v2, 3), temp);
   c_dd_add(temp, out, out); 
}

/*
   Computes the vector projection of _v_ on _u_
*/

void dd_vec3d_vector_proj(dd_vec3d v_out, dd_vec3d v, dd_vec3d u)
{
   double proj1[2];
   double proj2[2];
   double proj_scalar[2];
   
   dd_vec3d_scalar_prod(proj1, v, u);
   dd_vec3d_norm(proj2, u);
   c_dd_div(proj1, proj2, proj_scalar);
   dd_vec3d_mul_scalar(v_out, u, proj_scalar);   
}

/*
   Returns the scalar projection of _v_ on _u_
*/

double dd_vec3d_scalar_proj(double * proj_scalar, dd_vec3d v, dd_vec3d u)
{
   double proj1[2];
   double proj2[2];
   
   dd_vec3d_scalar_prod(proj1, v, u);
   dd_vec3d_norm(proj2, u);
   c_dd_div(proj1, proj2, proj_scalar);
}

/* 
   Sets the vector v_out to be equal to v_in
*/

void dd_vec3d_set(dd_vec3d v_out, dd_vec3d v_in)
{
   c_dd_copy(D_COORD(v_in, 1), D_COORD(v_out, 1));
   c_dd_copy(D_COORD(v_in, 2), D_COORD(v_out, 2));
   c_dd_copy(D_COORD(v_in, 3), D_COORD(v_out, 3));
}

/*
   Sets v_out to v_in1 - scalar * v_in2
*/

void dd_vec3d_sub_scalar_mul(dd_vec3d v_out, dd_vec3d v_in1, dd_vec3d v_in2, double * scalar)
{
   double temp[2];
   
   c_dd_mul(D_COORD(v_in2, 1), scalar, temp);
   c_dd_sub(D_COORD(v_in1, 1), temp, D_COORD(v_out, 1));
   c_dd_mul(D_COORD(v_in2, 2), scalar, temp);
   c_dd_sub(D_COORD(v_in1, 2), temp, D_COORD(v_out, 2));
   c_dd_mul(D_COORD(v_in2, 3), scalar, temp);
   c_dd_sub(D_COORD(v_in1, 3), temp, D_COORD(v_out, 3));
}

/*
   Sets v_out to v_in1 + scalar * v_in2
*/

void dd_vec3d_add_scalar_mul(dd_vec3d v_out, dd_vec3d v_in1, dd_vec3d v_in2, double * scalar)
{
   double temp[2];
   
   c_dd_mul(D_COORD(v_in2, 1), scalar, temp);
   c_dd_add(D_COORD(v_in1, 1), temp, D_COORD(v_out, 1));
   c_dd_mul(D_COORD(v_in2, 2), scalar, temp);
   c_dd_add(D_COORD(v_in1, 2), temp, D_COORD(v_out, 2));
   c_dd_mul(D_COORD(v_in2, 3), scalar, temp);
   c_dd_add(D_COORD(v_in1, 3), temp, D_COORD(v_out, 3));
}

/*
   Sets v_out to v_in1 - scalar * v_in2
*/

void dd_vec3d_sub_scalar_mul_d(dd_vec3d v_out, dd_vec3d v_in1, dd_vec3d v_in2, double scalar)
{
   double temp[2];
   
   c_dd_mul_dd_d(D_COORD(v_in2, 1), scalar, temp);
   c_dd_sub(D_COORD(v_in1, 1), temp, D_COORD(v_out, 1));
   c_dd_mul_dd_d(D_COORD(v_in2, 2), scalar, temp);
   c_dd_sub(D_COORD(v_in1, 2), temp, D_COORD(v_out, 2));
   c_dd_mul_dd_d(D_COORD(v_in2, 3), scalar, temp);
   c_dd_sub(D_COORD(v_in1, 3), temp, D_COORD(v_out, 3));
}

/*
   Sets v_out to v_in1 + scalar * v_in2
*/

void dd_vec3d_add_scalar_mul_d(dd_vec3d v_out, dd_vec3d v_in1, dd_vec3d v_in2, double scalar)
{
   double temp[2];
   
   c_dd_mul_dd_d(D_COORD(v_in2, 1), scalar, temp);
   c_dd_add(D_COORD(v_in1, 1), temp, D_COORD(v_out, 1));
   c_dd_mul_dd_d(D_COORD(v_in2, 2), scalar, temp);
   c_dd_add(D_COORD(v_in1, 2), temp, D_COORD(v_out, 2));
   c_dd_mul_dd_d(D_COORD(v_in2, 3), scalar, temp);
   c_dd_add(D_COORD(v_in1, 3), temp, D_COORD(v_out, 3));
}

/* 
   Multiplies the given vector by the given scalar
*/

void dd_vec3d_mul_scalar_d(dd_vec3d v_out, dd_vec3d v_in, double scalar)
{
   c_dd_mul_dd_d(D_COORD(v_in, 1), scalar, D_COORD(v_out, 1));
   c_dd_mul_dd_d(D_COORD(v_in, 2), scalar, D_COORD(v_out, 2));
   c_dd_mul_dd_d(D_COORD(v_in, 3), scalar, D_COORD(v_out, 3));
}
