/******************************************************************************

 mat3d.h
 Header file for mat3d.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MAT3D_H
#define MAT3D_H

#include <math.h>

#include "vec3d.h"

/*
   mat3d: There are numerous different types here:
          
   1) d_mat3dc_t : a 3x3 matrix of doubles whose columns are 3D vectors (vec3d_t's). 
   
   2) d_mat3dr_t : a 3x3 matrix of doubles whose rows are 3D vectors (vec3d_t's). 
   
   3) z_mat3dc_t : a 3x3 matrix of longs whose columns are 3D vectors (vec3u_t's). 
   
   4) z_mat3dr_t : a 3x3 matrix of longs whose rows are 3D vectors (vec3u_t's). 
      
   These types need to be initialised before being used. 
*/

// Defines (implicitly) a matrix of dimension 3 whose columns are 3d vectors of doubles

typedef d_vec3d * d_mat3dc_t;

// Defines (implicitly) a matrix of dimension 3 whose rows are 3d vectors of doubles

typedef d_vec3d * d_mat3dr_t;

// Defines (implicitly) a matrix of dimension 3 whose columns are 3d vectors of longs

typedef z_vec3d * z_mat3dc_t;

// Defines (implicitly) a matrix of dimension 3 whose rows are 3d vectors of longs

typedef z_vec3d * z_mat3dr_t;

void z_mat3dr_stack_init(z_mat3dr_t * C);

void z_mat3dr_stack_clear(void);

void z_mat3dc_stack_init(z_mat3dc_t * C);

void z_mat3dc_stack_clear(void);

void d_mat3dr_stack_init(d_mat3dr_t * C);

void d_mat3dr_stack_clear();

void d_mat3dc_stack_init(d_mat3dc_t * C);

void d_mat3dc_stack_clear();

void z_mat3dr_swap12(z_mat3dr_t C);

void z_mat3dr_swap13(z_mat3dr_t C);

void z_mat3dr_swap23(z_mat3dr_t C);

void z_mat3dc_swap12(z_mat3dc_t C);

void z_mat3dc_swap13(z_mat3dc_t C);

void z_mat3dc_swap23(z_mat3dc_t C);

void z_mat3dr_set_identity(z_mat3dr_t C);

void z_mat3dc_set_identity(z_mat3dc_t C);

void d_mat3dc_set_identity(d_mat3dc_t C);

void d_mat3dc_set_zero(d_mat3dc_t C);

void d_mat3dc_mul_z_mat3dc(d_mat3dc_t B_out, d_mat3dc_t B, z_mat3dc_t Z);

double d_mat3dc_det(d_mat3dc_t B);

void d_mat3dc_invert(d_mat3dc_t B_inv, d_mat3dc_t B);

void d_mat3dc_printf(d_mat3dc_t mat);

void d_mat3dc_scanf(d_mat3dc_t mat);

void z_mat3dc_printf(z_mat3dc_t mat);

void z_mat3dc_scanf(z_mat3dc_t mat);

void z_mat3dr_printf(z_mat3dr_t mat);

void z_mat3dr_scanf(z_mat3dr_t mat);

void d_mat3dc_gram_schmidt(d_mat3dc_t Q, d_mat3dc_t B_out, d_mat3dc_t B_in);

void d_mat3dc_LLL(z_mat3dr_t C, d_mat3dc_t B_in, double delta);

#endif
