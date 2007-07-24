/******************************************************************************

 mat3d.h
 Header file for mat3d.c.

 (C) 2006 William Hart

******************************************************************************/

#ifndef MAT3D_H
#define MAT3D_H

#include <math.h>

/*
   mat3d: There are five different types here:
          
   1) mat3d_t  : a 3x3 matrix implemented as a double subscripted array
                 to access row 1 column 2 of such an object B, use B[1][2]
   2) mat3dc_t : a 3x3 matrix whose columns are 3D vectors (vec3d_t's). To
                 access row 1 column 2, use B->c2->x1
   3) mat3dr_t : a 3x3 matrix whose rows are 3D vectors (vec3d_t's). To
                 access row 1 column 2, use B->r1->x2
   4) mat3dc_p : same as 2, except these need initialising before use. The
                 advantage is that since there is a layer of indirection,
                 swapping of columns can be done efficiently
   5) mat3dr_p : same as 3, except these need initialising before use. The
                 advantage is that since there is a layer of indirection,
                 swapping of rows can be done efficiently
                 
   Note that functions which accept objects of types 2 also accept those of type
   4 and those which accept objects of type 3 also accept those of type 5 respectively. 
   However, more efficient versions of some function exist solely for types 4 and 5. 
   
   Some functions accept matrices of one type and output a different type. This
   enables carefully thought out algorithms to involve a mix of different types
   where this is more efficient than using a single one of the types.
   
   The disadvantage of types 4 and 5 is that they need to be initialised and then
   cleared after use. There is some for the memory allocation associated with this.
   
*/

/*
   Define a square matrix of dimension 3 implemented as a double subscripted array
*/

typedef double mat3d[3][3];

/*
   Defines a square matrix of dimension 3 whose entries are doubles. The columns are 
   implemented as 3d vectors
*/

typedef struct mat3dc_s
{
   vec_3d_t c1;
   vec_3d_t c2;
   vec_3d_t c3;
} mat3dc;

typedef mat3dc mat3dc_t[1];

// Defines a matrix of dimension 3 whose columns are 3d vectors

typedef struct mat3dc_p
{
   vec_3d * c1;
   vec_3d * c2;
   vec_3d * c3;
} mat3dc_p;

/*
   Defines a square matrix of dimension 3 whose entries are doubles. The rows are 
   implemented as 3d vectors
*/

typedef struct mat3dr_s
{
   vec_3d_t r1;
   vec_3d_t r2;
   vec_3d_t r3;
} mat3dr;

typedef mat3dr mat3dr_t[1];

// Defines a matrix of dimension 3 whose rows are 3d vectors

typedef struct mat3dr_p
{
   vec_3d * r1;
   vec_3d * r2;
   vec_3d * r3;
} mat3dr_p;

void mat3dc_gram_schmidt(mat3dc_t Q, mat3dc_t B_out, mat3dc_t B_out);

void mat3dc_LLL(mat3dc_p B_out, mat3dc_p B_in, double delta);

#endif
