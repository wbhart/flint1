/******************************************************************************

 vec3d.h
 Header file for vec3d.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef VEC3D_H
#define VEC3D_H

#include <math.h>

#include "vector.h"

/*
   Defines a 3D vector whose coordinates are doubles (the length is implicitly 3)
*/

typedef d_vec d_vec3d;

/*
   Defines a 3D vector whose coordinates are longs (the length is implicitly 3)
*/

typedef z_vec z_vec3d;

void z_vec3d_stack_init(z_vec3d * vec);

void d_vec3d_stack_init(d_vec3d * vec);

void z_vec3d_stack_clear();

void d_vec3d_stack_clear();

void d_vec3d_scanf(d_vec3d vec);

void d_vec3d_printf(d_vec3d vec);

double d_vec3d_length(d_vec3d v);

double d_vec3d_norm(d_vec3d v);

void d_vec3d_normalise(d_vec3d v_out, d_vec3d v_in);

void d_vec3d_mul_scalar(d_vec3d v_out, d_vec3d v_in, double scalar);

double d_vec3d_scalar_prod(d_vec3d v1, d_vec3d v2);

void d_vec3d_vector_proj(d_vec3d v_out, d_vec3d v, d_vec3d u);

double d_vec3d_scalar_proj(d_vec3d v, d_vec3d u);

void d_vec3d_set(d_vec3d v_out, d_vec3d v_in);

void d_vec3d_sub_scalar_mul(d_vec3d v_out, d_vec3d v_in1, d_vec3d v_in2, double scalar);

void d_vec3d_add_scalar_mul(d_vec3d v_out, d_vec3d v_in1, d_vec3d v_in2, double scalar);

void z_vec3d_sub_scalar_mul(z_vec3d v_out, z_vec3d v_in1, z_vec3d v_in2, long scalar);

void z_vec3d_add_scalar_mul(z_vec3d v_out, z_vec3d v_in1, z_vec3d v_in2, long scalar);

#endif
