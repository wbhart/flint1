/******************************************************************************

 matrix.h
 
 Header file for matrix types

 (C) 2006 William Hart

******************************************************************************/

#ifndef MATRIX_H
#define MATRIX_H

#include <math.h>

/*
   Matrix of unsigned longs
*/

typedef z_vec * z_mat;

typedef struct z_mat_s
{
   z_mat entries;
   unsigned long * data; // Internal use, will point to a single block of memory containing all the coefficients 
   unsigned long rows; // Number of rows
   unsigned long cols; // Number of columns
   int format;  // If 0 the rows of this matrix are vectors, if 1 the columns are vectors
} z_mat_s;

typedef z_mat_s z_mat_t[1];

/*
   Matrix of doubles
*/

typedef d_vec * d_mat;

typedef struct d_mat_s
{
   d_mat entries; 
   double * data; // Internal use, will point to a single block of memory containing all the coefficients 
   unsigned long rows; // Number of rows
   unsigned long cols; // Number of columns
   int format;  // If 0 the rows of this matrix are vectors, if 1 the columns are vectors
} d_mat_s;

typedef d_mat_s d_mat_t[1];

#endif
