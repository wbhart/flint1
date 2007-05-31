/******************************************************************************

 ZmodF_mul.h

 Copyright (C) 2007, David Harvey
 
******************************************************************************/

#ifndef FLINT_ZMODF_MUL_H
#define FLINT_ZMODF_MUL_H

#include <stdlib.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "ZmodFpoly.h"


/*
This struct stores info used to speed up multiplications mod p
for a specific n.
*/
typedef struct
{
   unsigned long n;

   // 0 means use plain old mpn_mul_n; 1 means use the negacyclic FFT.
   int use_fft;
   
   // ------------------ fields used only for mpn_mul_n
   
   // scratch buffer of length 2*n
   mp_limb_t* scratch;
   
   // ------------------ fields used only for negacyclic FFT

   ZmodFpoly_t polys[2];

} ZmodF_mul_precomp_struct;


// ZmodF_mul_precomp_t allows reference-like semantics for
// ZmodF_mul_precomp_struct:
typedef ZmodF_mul_precomp_struct ZmodF_mul_precomp_t[1];


// Input is a requested coefficient length n. Output is n' >= n which is
// optimal with respect to running a negacyclic FFT.
// (i.e. If n is sufficiently small that mpn_mul_n should be used directly,
// then it returns n. Otherwise it might round n up slightly to a multiple
// of a power of two, to make a negacyclic convolution possible.)
// If depth != NULL, it will also write the negacyclic transform depth to
// that location.
unsigned long ZmodF_mul_precomp_get_feasible_n(unsigned long* depth,
                                               unsigned long n);


// initialises ZmodF_mul_precomp_t for a given n
// (The squaring flag is 1 if you are going to use this object for squaring.
// This doesn't affect whether the struct can be used for squaring or
// multiplying, the only effect is to possibly modify the tuning parameters.)
void ZmodF_mul_precomp_init(ZmodF_mul_precomp_t info, unsigned long n,
                            int squaring);
                            
// releases resources
void ZmodF_mul_precomp_clear(ZmodF_mul_precomp_t info);

// sets res := a * b, assuming they all have the same n as the given
// ZmodF_mul_precomp_t.
void ZmodF_mul_precomp(ZmodF_mul_precomp_t, ZmodF_t res, ZmodF_t a, ZmodF_t b);
// res := a * b
void ZmodF_sqr_precomp(ZmodF_mul_precomp_t, ZmodF_t res, ZmodF_t a);


/*
   res := a * b
   
   PRECONDITIONS:
      Any combination of aliasing among res, a, b is allowed.
      scratch must be a buffer of length 2*n, and must NOT alias a, b, res.

*/
void ZmodF_mul(ZmodF_t res, ZmodF_t a, ZmodF_t b, mp_limb_t* scratch,
               unsigned long n);


/*
   res := a * a
   
   PRECONDITIONS:
      a may alias res.
      scratch must be a buffer of length 2*n, and must NOT overlap a or res.
*/
void ZmodF_sqr(ZmodF_t res, ZmodF_t a, mp_limb_t* scratch, unsigned long n);


#endif

// end of file ****************************************************************
