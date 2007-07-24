/******************************************************************************

 vec3d.h
 Header file for vec3d.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef VEC3D_H
#define VEC3D_H

#include <math.h>

/*
   Defines a 3D vector whose coordinates are doubles
*/

typedef struct vec3d
{
   double x1;
   double x2;
   double x3;
} vec3d;

typedef vec3d vec3d_t[1];

/*
   Returns the Euclidean length of a 3d vector
*/

static inline double vec3d_length(vec3d_t v)
{
   return sqrt(v->x1*v->x2 + v->x2*v->x2 + v->x3*v->x3);
}

/*
   Returns the square of the Euclidean length of a 3d vector
*/

static inline double vec3d_norm(vec3d_t v)
{
   return v->x1*v->x2 + v->x2*v->x2 + v->x3*v->x3;
}

/*
   Normalises the given input vector (i.e. divides it by its length)
*/

static inline void vec3d_normalise(vec3d_t v_out, vec3d_t v_in)
{
   double length = vec3d_length(v_in);
   
   v_out->x1 = v_in->x1 / length;
   v_out->x2 = v_in->x2 / length;
   v_out->x3 = v_in->x3 / length;
}

/* 
   Multiplies the given vector by the given scalar
*/

static inline double vec3d_mul_scalar(vec3d_t v_out, vec3d_t v_in, double scalar)
{
   v_out->x1 = v_in->x1*scalar;
   v_out->x2 = v_in->x2*scalar;
   v_out->x3 = v_in->x3*scalar;  
}

/*
   Computes the scalar product of the two input vectors
*/

static inline double vec3d_scalar_prod(vec3d_t v1, vec3d_t v2)
{
   return v1->x1*v2->x1 + v1->x2*v2->x2 + v1->x3*v2->x3;
}

/*
   Computes the vector projection of _v_ on _u_
*/

static inline void vec3d_vector_proj(vec3d_t v_out, vec3d_t v, vec3d_t u)
{
   double proj_scalar = vec3d_scalar_prod(v, u) / vec3d_scalar_prod(u, u);
   vec3d_mul_scalar(v_out, u, proj_scalar);   
}

/*
   Returns the scalar projection of _v_ on _u_
*/

static inline void vec3d_scalar_proj(vec3d_t v, vec3d_t u)
{
   return vec3d_scalar_prod(v, u) / vec3d_scalar_prod(u, u);
}

/* 
   Sets the vector v_out to be equal to v_in
*/

static inline void vec3d_set(vec3d_t v_out, vec3d_t v_in)
{
   v_out->x1 = v_in->x1;
   v_out->x2 = v_in->x2;
   v_out->x3 = v_in->x3;   
}

/*
   Sets v_out to v_in1 - scalar * v_in2
*/

static inline void vec3d_sub_scalar_mul(vec3d_t v_out, vec3d_t v_in1, vec3d_t v_in2, double scalar)
{
   v_out->x1 = v_in1->x1 - scalar * v_in2->x1;
   v_out->x2 = v_in1->x2 - scalar * v_in2->x2;
   v_out->x3 = v_in1->x3 - scalar * v_in2->x3;
}

/*
   Sets v_out to v_in1 + scalar * v_in2
*/

static inline void vec3d_add_scalar_mul(vec3d_t v_out, vec3d_t v_in1, vec3d_t v_in2, double scalar)
{
   v_out->x1 = v_in1->x1 + scalar * v_in2->x1;
   v_out->x2 = v_in1->x2 + scalar * v_in2->x2;
   v_out->x3 = v_in1->x3 + scalar * v_in2->x3;
}


#endif
