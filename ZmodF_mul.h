/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/******************************************************************************

 ZmodF_mul.h

 Copyright (C) 2007, David Harvey and William Hart
 
 Routines for multiplication of elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS.
 
******************************************************************************/

#ifndef FLINT_ZMODF_MUL_H
#define FLINT_ZMODF_MUL_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "ZmodF_poly.h"


/*
   Several algorithms for multiplication mod p:

   ZMODF_MUL_PLAIN: use mpn_mul_n and then reduce mod p

   ZMODF_MUL_THREEWAY: if n = 3m, do multiplication mod B^m + 1 and mod
     B^2m - B^m + 1 separately, combine using CRT to get answer mod B^n + 1.

   ZMODF_MUL_FFT: split into polynomial of length 2^depth with coefficients
     mod (B^m + 1)*B^k, where k = 0, 1, or 2. Do negacyclic convolution mod
     B^m + 1 using FFT; do convolution mod B^k using naive algorithm; combine
     results by CRT.

*/
#define ZMODF_MUL_ALGO_PLAIN 0
#define ZMODF_MUL_ALGO_THREEWAY 1
#define ZMODF_MUL_ALGO_FFT 2


/*
   This struct manages temporary allocation and tuning parameters for a
   specific n.
*/
typedef struct
{
   unsigned long n;

   // possible values for algo are the ZMOD_MUL_ALGO_xyz constants above
   int algo;
   
   // this flag indicates that this struct is being used for squaring, in
   // which case less memory will be allocated in the init routines
   int squaring;
   
   // scratch buffer:
   // of length 2n           (for ZMODF_MUL_ALGO_PLAIN)
   // or length 3n+1         (for ZMODF_MUL_ALGO_THREEWAY)
   // or length 3*k*2^depth  (for ZMODF_MUL_ALGO_FFT)
   // or NULL if not allocated
   mp_limb_t* scratch;

   // for ZMODF_MUL_ALGO_THREEWAY, m = n/3.
   // for ZMODF_MUL_ALGO_FFT, the FFT operates mod B^m + 1, and the naive
   // algorithm works mod B^k
   unsigned long m, k;
   
   // used only for ZMODF_MUL_ALGO_FFT
   // (if squaring == 1, only the first one is init'd)
   ZmodF_poly_t polys[2];

} ZmodF_mul_info_struct;


// ZmodF_mul_info_t allows reference-like semantics for
// ZmodF_mul_info_struct:
typedef ZmodF_mul_info_struct ZmodF_mul_info_t[1];


/*
   Initialises ZmodF_mul_info_t for a given n. This function automatically
   selects the best underlying multiplication algorithm for the given n.

   The squaring flag is 1 if you intend to use this object for squaring. The
   object will then use up less memory. If squaring == 0, it's still possible
   to use this object for squaring, but if squaring == 1, you can't use it for
   arbitrary multiplication.

   WARNING: the multiplication time does NOT increase monotonically with n.
   * If n is divisible by 3, the "threeway" algorithm is available, which is
     faster than the "plain" algorithm for n >= 24 (on our opteron test
     machine).
   * For large n, the "fft" algorithm is available, but there are conditions on
     the 2-divisibility of n (not very onerous though).

   NOTE: this function uses stack based memory management.
*/
void ZmodF_mul_info_init(ZmodF_mul_info_t info, unsigned long n, int squaring);


// the following functions initialise with a specific algorithm:
void ZmodF_mul_info_init_plain(ZmodF_mul_info_t info, unsigned long n,
                               int squaring);

void ZmodF_mul_info_init_threeway(ZmodF_mul_info_t info, unsigned long n,
                                  int squaring);

void ZmodF_mul_info_init_fft(ZmodF_mul_info_t info, unsigned long n,
                             unsigned long depth, unsigned long m,
                             unsigned long k, int squaring);

                            
// releases resources
void ZmodF_mul_info_clear(ZmodF_mul_info_t info);

// sets res := a * b using the given ZmodF_mul_info_t object
void ZmodF_mul_info_mul(ZmodF_mul_info_t, ZmodF_t res, ZmodF_t a, ZmodF_t b);


// the following functions are for standalone multiplying/squaring, without
// the use of the precomputation stuff above:

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



// ============================================================================
// the following functions are exported for testing purposes:

void _ZmodF_mul_fft_split(ZmodF_poly_t poly, ZmodF_t x, unsigned long n);

void _ZmodF_mul_fft_combine(ZmodF_t x, ZmodF_poly_t poly,
                            unsigned long m, unsigned long k, unsigned long n);

void _ZmodF_mul_fft_convolve_modB(
            mp_limb_t* out, mp_limb_t* in1, mp_limb_t* in2, unsigned long len);

void _ZmodF_mul_fft_convolve_modB2(
            mp_limb_t* out, mp_limb_t* in1, mp_limb_t* in2, unsigned long len);



void _ZmodF_mul_threeway_reduce1(ZmodF_t res, ZmodF_t a, unsigned long m);

void _ZmodF_mul_threeway_reduce2(mp_limb_t* res, ZmodF_t a, unsigned long m);

#ifdef __cplusplus
 }
#endif

#endif

// end of file ****************************************************************
