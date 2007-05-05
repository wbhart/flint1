/******************************************************************************

 ZmodFermatVec.h

 Copyright (C) 2007
 David Harvey
 
 Vectors over ZmodFermat, routines for performing truncated Schoenhage-Strassen
 FFTs.
  
******************************************************************************/

#include "ZmodFermat.h"


/*
ZmodFermatVec_t represents a vector with coefficients in Z/pZ, where
p = B^n + 1, B = 2^FLINT_BITS_PER_LIMB. Coefficients are represented in the
format described in ZmodFermat.h.

Each vector has a fixed transform length M = 2^m, specified at creation time,
where m >= 0.

M must divide 4*n*FLINT_BITS_PER_LIMB. This ensures that Z/pZ has Mth roots
of unity (since we have available a sqrt2).

A vector may be in either "natural" or "transformed" state (i.e. if the vector
represents a polynomial, then "natural" means just the list of coefficients,
and "transformed" means the fourier transform). The distinction is only
conceptual, i.e. the vector doesn't keep track of which state it is in, but
it does have some bearing on the semantics of ZmodFermatVec attributes.

If a vector is "natural", then the length attribute indicates how many
coefficients are actually used. Coefficients beyond "length" are assumed to
be *zero*.

If a vector is "transformed", then the length attribute indicates how many
fourier coefficients are *known*. Coefficients beyond "length" might not
be zero, their values are simply unknown.

*/

typedef struct
{
   unsigned long m, M;
   unsigned long n;
   unsigned long length;

   // Single chunk of memory where all coefficients live.
   mp_limb_t* storage;
   
   // Pointers to M coefficients. These may get permuted during FFTs.
   ZmodFermat* coeffs;
   
   // A scratch buffer. It may get permuted along with the other coefficients.
   // (NOTE: in future, this may be an array of scratch buffers.)
   ZmodFermat* scratch;

} ZmodFermatVec_struct;

// ZmodFermatVec_t allows reference-like semantics for ZpolyFermatVec_struct:
typedef ZmodFermatVec_struct ZmodFermatVec_t[1];


/*
Initialises a ZmodFermatVec_t with supplied m and n, and length = 0.
Coefficients are not zeroed out.
*/
void ZmodFermatVec_init(ZmodFermatVec_t x, unsigned long m, unsigned long n);


/*
Frees resources for supplied ZmodFermatVec_t.
*/
void ZmodFermatVec_clear(ZmodFermatVec_t x);


/*
Copies y into x. Only y.length coefficients are copied.

PRECONDITIONS:
   x and y must have the same m and n.
*/
void ZmodFermatVec_set(ZmodFermatVec_t x, ZmodFermatVec_t y);


/*
Converts from "natural" to "transformed" state. The number of fourier
coefficients desired is in output_length.

Output is inplace. (Note that in general *all* M coefficients will get
overwritten in intermediate steps.)

x.length is set to output_length at the end.

PRECONDITIONS:
   0 <= output_length <= x.M
*/
void ZmodFermatVec_FFT(ZmodFermatVec_t x, unsigned long output_length);


/*
Converts from "transformed" to "natural" state.

That is, it runs an inverse fourier transform, *assuming* that the supplied
vector is actually the fourier transform of a vector which has all zeroes
beyond index x.length.

Result is inplace, x.length is not modified. (Note: after it's finished, the
coefficients beyond x.length will contain garbage.)
*/
void ZmodFermatVec_IFFT(ZmodFermatVec_t x);



/*
todo: need to think about if/how to expose the more sophisticated truncated
IFFT functionality, where you can tell it the values of some of the output
coefficients, which are not necessarily zero. Not sure how useful this would
be.
*/


/*
todo: need some functions for converting in and out of ZmodFermat_t format.
But probably this belongs in the ZmodFermat module.
*/


/*
Pointwise multiplication of two ZmodFermatVecs. Only coefficients up to
x.length are multiplied.

PRECONDITIONS:
   Any combination of aliasing among res, x, y is allowed.
   x and y must have the same m, n, and length.
   res must have the same m, n as x and y.
*/
void ZmodFermatVec_mul(ZmodFermatVec_t res,
                       ZmodFermatVec_t x, ZmodFermatVec_t y);


/*
Convolution of two ZmodFermatVecs.

The resulting length will be x.length + y.length - 1. If this is more than M,
then the resulting length is M, and the convolution is actually cyclic of
length M.

PRECONDITIONS:
   Any combination of aliasing among res, x, y is allowed.
   x and y must have the same m, n, and length.
   res must have the same m, n as x and y.

NOTE:
   x and y will both be converted to "transformed" state. If you don't like
   it, make a copy.
*/
void ZmodFermatVec_convolve(ZmodFermatVec_t res,
                            ZmodFermatVec_t x, ZmodFermatVec_t y);
{
}


/*
Pointwise addition of two ZmodFermatVecs. Only coefficients up to x.length are
multiplied.

PRECONDITIONS:
   Any combination of aliasing among res, x, y is allowed.
   x and y must have the same m, n, and length.
   res must have the same m, n as x and y.
*/
void ZmodFermatVec_add(ZmodFermatVec_t res,
                       ZmodFermatVec_t x, ZmodFermatVec_t y);


/* subtraction too I guess */


/*
Divides all coefficients by M mod p. This should be used after running an
inverse transform.
*/
void ZmodFermatVec_scale(ZmodFermatVec_t x);


/*
Normalises all coefficients (up to x.length) to be in the range [0, p).
*/
void ZmodFermatVec_normalise(ZmodFermatVec_t x);


// end of file ****************************************************************
