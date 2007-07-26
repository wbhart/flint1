/******************************************************************************

 mat3d.h
 Header file for mat3d.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MAT3D_H
#define MAT3D_H

#include <math.h>

/*
   mat3d: There are numerous different types here:
          
   1) d_mat3dc_t : a 3x3 matrix of doubles whose columns are 3D vectors (vec3d_t's). 
   
   2) d_mat3dr_t : a 3x3 matrix of doubles whose rows are 3D vectors (vec3d_t's). 
   
   3) z_mat3dc_t : a 3x3 matrix of unsigned longs whose columns are 3D vectors (vec3u_t's). 
   
   4) z_mat3dr_t : a 3x3 matrix of unsigned longs whose rows are 3D vectors (vec3u_t's). 
      
   These types need to be initialised before being used. 
*/

// Defines (implicitly) a matrix of dimension 3 whose columns are 3d vectors of doubles

typedef d_mat d_mat3dc_t;

// Defines (implicitly) a matrix of dimension 3 whose rows are 3d vectors of doubles

typedef d_mat d_mat3dr_t;

// Defines (implicitly) a matrix of dimension 3 whose columns are 3d vectors of unsigned longs

typedef z_mat z_mat3dc_t;

// Defines (implicitly) a matrix of dimension 3 whose rows are 3d vectors of unsigned longs

typedef z_mat z_mat3dr_t;


void d_mat3dc_gram_schmidt(d_mat3dc_t Q, d_mat3dc_t B_out, d_mat3dc_t B_in);

void d_mat3dc_LLL(z_mat3dr_t B_out, d_mat3dc_t B_in, double delta);

#endif
