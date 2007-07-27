/******************************************************************************

 matrix.h
 
 Header file for vector types

 (C) 2006 William Hart

******************************************************************************/

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>

/*
   Vector of unsigned longs
*/

typedef unsigned long * z_vec;

typedef struct z_vec_s
{
   z_vec coords;
   unsigned long length;
} z_vec_s;

typedef z_vec_s z_vec_t[1];

/*
   Vector of doubles
*/

typedef double * d_vec;

typedef struct d_vec_s
{
   d_vec coords;
   unsigned long length;
} d_vec_s;

typedef d_vec_s d_vec_t[1];

#endif
