/*
   zn_poly_internal.h:  main header file #included internally by zn_poly
                        modules
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


/*

   IMPORTANT NOTE!!!!!!
   
   Everything in this file is internal, and may change incompatibly in
   future versions of zn_poly. You have been warned!

*/


#ifndef ZN_POLY_INTERNAL_H
#define ZN_POLY_INTERNAL_H

#include <gmp.h>
#include <stddef.h>
#include <stdio.h>
#include "zn_poly.h"


#ifdef __cplusplus
extern "C" {
#endif
 

/*
   For integers a >= 1 and b >= 1, returns ceil(a / b).
*/
#define CEIL_DIV(a, b) ((((a) - 1) / (b)) + 1)


/*
   For integers a >= 1 and 0 <= r < ULONG BITS, returns ceil(a / 2^r).
*/
#define CEIL_DIV_2EXP(a, r) ((((a) - 1) >> (r)) + 1)


#define ZNP_MIN(aaa, bbb) (((aaa) < (bbb)) ? (aaa) : (bbb))
#define ZNP_MAX(aaa, bbb) (((aaa) > (bbb)) ? (aaa) : (bbb))


/*
   Estimate of the L1 cache size, in bytes.
   If this is a bit on the small side, it's probably not a big deal.
   If it's on the big side, that might start to seriously degrade performance.
*/
#define ZNP_CACHE_SIZE 32768


/*
   Returns ceil(log2(x)).
   x must be >= 1.
*/
#define ceil_lg \
    ZNP_ceil_lg
int ceil_lg(ulong x);


/*
   Returns floor(log2(x)).
   Returns -1 for x == 0.
*/
#define floor_lg \
    ZNP_floor_lg
int floor_lg(ulong x);


/*
   The ZNP_FASTALLOC and ZNP_FASTFREE macros are used for allocating memory
   which is taken off the stack if the request is small enough, or off the
   heap if not.
   
   Example usage:
   
   ZNP_FASTALLOC(stuff, int, 100, n);

   This does two things. It allocates an array of 100 ints on the stack.
   It also declares a pointer "int* stuff", which points to a block of ints
   of length n. If n <= 100, the block will be the one just allocated on the
   stack. If n > 100, the block will be found using malloc.
   
   Then afterwards, you need to do:
   
   ZNP_FASTFREE(stuff);
   
   This will call free() if the block was originally taken off the heap.
*/

#define ZNP_FASTALLOC(ptr, type, reserve, request)        \
   size_t __FASTALLOC_request_##ptr = (request);          \
   type* ptr;                                             \
   type __FASTALLOC_##ptr [reserve];                      \
   if (__FASTALLOC_request_##ptr <= (reserve))            \
      ptr = __FASTALLOC_##ptr;                            \
   else                                                   \
      ptr = (type*) malloc(sizeof(type) * __FASTALLOC_request_##ptr);
      

#define ZNP_FASTFREE(ptr)                                 \
   if (ptr != __FASTALLOC_##ptr)                          \
      free(ptr);



/*
   Stores tuning data for moduli of a specific bitsize.
*/
#define tuning_info_t \
    ZNP_tuning_info_t
typedef struct
{
   // Crossovers from KS1 -> KS2 and KS2 -> KS4 for Kronecker substitution
   // multiplication algorithms (see mul_ks.c)
   size_t mul_KS2_crossover;
   size_t mul_KS4_crossover;
   // Crossover from KS4 -> fft multiplication algorithm
   size_t mul_fft_crossover;

   // As above, but for squaring.
   size_t sqr_KS2_crossover;
   size_t sqr_KS4_crossover;
   size_t sqr_fft_crossover;

   // Crossover from plain -> fft middle product
   size_t midmul_fft_crossover;
   
   // for negacyclic multiplications, switch from KS to Nussbaumer FFT
   // when length reaches 2^nuss_mul_crossover
   unsigned nuss_mul_crossover;
   // ditto for nussbaumer squaring
   unsigned nuss_sqr_crossover;
   
} tuning_info_t;


/*
   Global array of tuning_info_t's, one for each bitsize.
*/
#define tuning_info \
    ZNP_tuning_info
extern tuning_info_t tuning_info[];



/* ============================================================================

     stuff from pack.c

============================================================================ */

/*
   Computes
   res := 2^lead * (op[0] + op[skip]*2^bits + op[2*skip]*2^(2*bits) + ...
                              + op[(len-1)*skip]*2^((len-1)*bits) ).
                              
   Assumes each op[i] satisfies 0 <= op[i] < 2^bits.

   Must have 0 < bits < 3*ULONG_BITS.
   
   If _requested_ is zero, then exactly ceil(len * bits / GMP_NUMB_BITS)
   limbs are written. Otherwise, the output will be zero-padded up to exactly
   _requested_ limbs, which must be at least ceil(len * bits / GMP_NUMB_BITS).
   
*/
#define zn_array_pack \
    ZNP_zn_array_pack
void zn_array_pack(mp_limb_t* res, const ulong* op, size_t len,
                   ptrdiff_t skip, unsigned bits, unsigned lead,
                   size_t requested);


/*
   Let op be an integer of the form

      2^lead * (a[0] + a[1]*2^bits + ... + a[len-1]*2^((len-1)*bits)),

   where 0 <= a[i] < 2^bits for each i.
   
   This function reads off the a[i]'s and stores them at res. Each output
   coefficient occupies exactly ceil(bits / ULONG_BITS) words. The input
   should be exactly ceil((lead + bits*len) / GMP_NUMB_BITS) limbs long.
   
   Must have 0 < bits < 3*ULONG_BITS.
*/
#define zn_array_unpack \
    ZNP_zn_array_unpack
void zn_array_unpack(ulong* res, const mp_limb_t* op,
                     size_t len, unsigned bits, unsigned lead);



/* ============================================================================

     stuff from mul.c

============================================================================ */

/*
   Identical to zn_array_mul(), except for the _fastred_ flag.
   
   If fastred is cleared, the output is the same as for zn_array_mul().
   
   If fastred is set, the routine uses the fastest modular reduction strategy
   available for the given parameters. The result will come out divided by a
   fudge factor, which can be recovered via _zn_array_mul_get_fudge().
*/
#define _zn_array_mul \
    ZNP__zn_array_mul
void _zn_array_mul(ulong* res, const ulong* op1, size_t len1,
                   const ulong* op2, size_t len2,
                   int fastred, const zn_mod_t mod);

#define _zn_array_mul_get_fudge \
    ZNP__zn_array_mul_get_fudge
ulong _zn_array_mul_get_fudge(size_t len1, size_t len2, int squaring,
                              const zn_mod_t mod);



/* ============================================================================

     stuff from mul_ks.c

============================================================================ */

/*
   These are the same as zn_array_mul(). They use four different types of
   Kronecker substitution. They automatically use a faster algorithm for
   squaring (if the inputs are identical buffers). See mul_ks.c.

   Aliasing of all operands allowed.

   Must have len1 >= len2 >= 1.
   
   If the _redc_ flag is set, the outputs will be divided by -B mod n.
   (Only allowed if the modulus is odd.)
*/
#define zn_array_mul_KS1 \
    ZNP_zn_array_mul_KS1
void zn_array_mul_KS1(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod);

#define zn_array_mul_KS2 \
    ZNP_zn_array_mul_KS2
void zn_array_mul_KS2(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod);

#define zn_array_mul_KS3 \
    ZNP_zn_array_mul_KS3
void zn_array_mul_KS3(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod);

#define zn_array_mul_KS4 \
    ZNP_zn_array_mul_KS4
void zn_array_mul_KS4(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int redc,
                      const zn_mod_t mod);




#define zn_array_recip_fix_reduce \
    ZNP_zn_array_recip_fix_reduce
void zn_array_recip_fix_reduce(ulong* res, ptrdiff_t skip, const ulong* op1,
                               const ulong* op2, size_t len, unsigned bits,
                               int redc, const zn_mod_t mod);



/* ============================================================================

     zn_pmf_t stuff

============================================================================ */


/*
   Let R = Z/nZ. A zn_pmf_t ("pmf" = "polynomial modulo fermat") represents an
   element of S = R[Y]/(Y^M + 1). This is used as the coefficient ring in
   the Schonhage and Nussbaumer FFT routines.

   It's an array of ulongs of length M + 1, where M = 2^lgM is a power of two.
   (The value M is not stored, the caller needs to keep track of it.)

   The first value in the array is an integer b called the "bias" of the
   representation.

   The remaining M values are coefficients a_0, ..., a_{M-1}.

   These together represent the polynomial

     Y^b (a_0 + a_1 Y + ... + a_{M-1} Y^{M-1}).
     
   Note that elements of S do not have a unique representation in this form;
   in fact they have one possible representation for each value of b in
   [0, 2M). (By allowing nonzero bias, we get more efficient in-place FFT
   butterflies.) The stored bias value need not be in [0, 2M), but it is
   interpreted mod 2M.

   Currently the values a_i are always normalised into [0, n). Later we might
   drop that restriction to obtain faster butterflies...

*/
#define zn_pmf_t \
    ZNP_zn_pmf_t
typedef ulong* zn_pmf_t;

#define zn_pmf_const_t \
    ZNP_zn_pmf_const_t
typedef const ulong* zn_pmf_const_t;


#define zn_pmf_scalar_mul \
    ZNP_zn_pmf_scalar_mul
ZNP_INLINE
void zn_pmf_scalar_mul(zn_pmf_t op, ulong M, ulong x, const zn_mod_t mod)
{
   zn_array_scalar_mul(op + 1, op + 1, M, x, mod);
}


/*
   op := 0, with bias reset to zero too
*/
#define zn_pmf_zero \
    ZNP_zn_pmf_zero
ZNP_INLINE
void zn_pmf_zero(zn_pmf_t op, ulong M)
{
   for (M++; M > 0; M--)
      *op++ = 0;
}


/*
   op := op / 2

   Modulus must be odd.
*/
#define zn_pmf_divby2 \
    ZNP_zn_pmf_divby2
ZNP_INLINE
void zn_pmf_divby2(zn_pmf_t op, ulong M, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);

   for (op++; M > 0; M--, op++)
      *op = zn_mod_divby2(*op, mod);
}


/*
   res := op
*/
#define zn_pmf_set \
    ZNP_zn_pmf_set
ZNP_INLINE
void zn_pmf_set(zn_pmf_t res, zn_pmf_t op, ulong M)
{
   for (M++; M > 0; M--)
      *res++ = *op++;
}


/*
   op := Y^r * op
*/
#define zn_pmf_rotate \
    ZNP_zn_pmf_rotate
ZNP_INLINE
void zn_pmf_rotate(zn_pmf_t op, ulong r)
{
   op[0] += r;
}


/*
   op1 := op2 + op1
   op2 := op2 - op1

   Inputs must be [0, n); outputs will be in [0, n).
*/
#define zn_pmf_bfly \
    ZNP_zn_pmf_bfly
void zn_pmf_bfly(zn_pmf_t op1, zn_pmf_t op2, ulong M, const zn_mod_t mod);

/*
   op1 := op1 + op2

   Inputs must be [0, n); outputs will be in [0, n).
*/
#define zn_pmf_add \
    ZNP_zn_pmf_add
void zn_pmf_add(zn_pmf_t op1, const zn_pmf_t op2, ulong M, const zn_mod_t mod);

/*
   op1 := op1 - op2

   Inputs must be [0, n); outputs will be in [0, n).
*/
#define zn_pmf_sub \
    ZNP_zn_pmf_sub
void zn_pmf_sub(zn_pmf_t op1, const zn_pmf_t op2, ulong M, const zn_mod_t mod);




/*
   These functions are exported just for profiling purposes:
*/

#define zn_array_bfly_inplace \
    ZNP_zn_array_bfly_inplace
void zn_array_bfly_inplace(ulong* op1, ulong* op2, ulong len,
                           const zn_mod_t mod);

#define zn_array_add_inplace \
    ZNP_zn_array_add_inplace
void zn_array_add_inplace(ulong* op1, const ulong* op2, ulong len,
                          const zn_mod_t mod);

#define zn_array_sub_inplace \
    ZNP_zn_array_sub_inplace
void zn_array_sub_inplace(ulong* op1, const ulong* op2, ulong len,
                          const zn_mod_t mod);



/* ============================================================================

     zn_pmf_vec_t stuff

============================================================================ */


/*
   A zn_pmf_vec_t stores a vector of length K = 2^lgK of elements of S.
   
   Used to represent an element of S[Z]/(Z^K + 1) or S[Z]/(Z^K - 1), or some
   other quotient like that.
   
   The functions zn_pmf_vec_init/clear should be used to allocate storage for
   this type. Also sometimes fake ones get created temporarily to point at
   sub-vectors of existing vectors.
*/
#define zn_pmf_vec_struct \
    ZNP_zn_pmf_vec_struct
typedef struct
{
   // points to the first coefficient
   zn_pmf_t data;
   
   // number of coefficients
   ulong K;
   unsigned lgK;      // lg2(K)
   
   // length of coefficients (see definition of zn_pmf_t)
   ulong M;
   unsigned lgM;      // lg2(M)

   // distance between adjacent coefficients, measured in ulongs
   // (this is at least M + 1, might be more)
   ptrdiff_t skip;
   
   // associated modulus
   const zn_mod_struct* mod;
}
zn_pmf_vec_struct;

#define zn_pmf_vec_t \
    ZNP_zn_pmf_vec_t
typedef zn_pmf_vec_struct zn_pmf_vec_t[1];


/*
   Checks that vec1 and vec2 have compatible data,
   i.e. have the same K, M, mod.
*/
#define zn_pmf_vec_compatible \
    ZNP_zn_pmf_vec_compatible
ZNP_INLINE
int zn_pmf_vec_compatible(const zn_pmf_vec_t vec1, const zn_pmf_vec_t vec2)
{
   return (vec1->K == vec2->K) && (vec1->M == vec2->M) &&
          (vec1->mod == vec2->mod);
}


/*
   Initialises _res_ with given parameters, allocates memory.
*/
#define zn_pmf_vec_init \
    ZNP_zn_pmf_vec_init
void zn_pmf_vec_init(zn_pmf_vec_t res, unsigned lgK, ptrdiff_t skip,
                     unsigned lgM, const zn_mod_t mod);


/*
   Initialises _res_ in preparation for a Nussbaumer multiplication of
   length 2^lgL.
*/
#define zn_pmf_vec_init_nussbaumer \
    ZNP_zn_pmf_vec_init_nussbaumer
void zn_pmf_vec_init_nussbaumer(zn_pmf_vec_t res, unsigned lgL,
                                const zn_mod_t mod);


/*
   Destroys op, frees all associated memory.
*/
#define zn_pmf_vec_clear \
    ZNP_zn_pmf_vec_clear
void zn_pmf_vec_clear(zn_pmf_vec_t op);


/*
   res := op
*/
#define zn_pmf_vec_set \
    ZNP_zn_pmf_vec_set
void zn_pmf_vec_set(zn_pmf_vec_t res, const zn_pmf_vec_t op);


/*
   Multiplies first _len_ coefficients of _op_ by _x_.
*/
#define zn_pmf_vec_scalar_mul \
    ZNP_zn_pmf_vec_scalar_mul
void zn_pmf_vec_scalar_mul(zn_pmf_vec_t op, ulong len, ulong x);


/*
   Multiplies the first _length_ coefficients of op1 and op2, puts result
   in res.
   
   It's okay for res to alias op1 or op2. The modulus must be odd.
   
   If the special_first_two flag is set, the routine assumes that the first
   two coefficients are of length only M/2 (this is the typical situation
   after performing the FFT), and multiplies them more quickly accordingly.
   
   The routine automatically selects KS or Nussbaumer multiplication
   depending on the bitsize and on M.
   
   The output will be divided by a fudge factor, which can be retrieved
   via zn_pmf_vec_mul_get_fudge().
   
   Automatically uses specialised squaring algorithm if the inputs are the
   same zn_pmf_vec_t object.
*/
#define zn_pmf_vec_mul \
    ZNP_zn_pmf_vec_mul
void zn_pmf_vec_mul(zn_pmf_vec_t res, const zn_pmf_vec_t op1,
                    const zn_pmf_vec_t op2, ulong length,
                    int special_first_two);

#define zn_pmf_vec_mul_get_fudge \
    ZNP_zn_pmf_vec_mul_get_fudge
ulong zn_pmf_vec_mul_get_fudge(unsigned lgM, int squaring, const zn_mod_t mod);


/*
   Modifies the _data_ and _skip_ members of op to make it look as if the
   first coefficient is the one at index (length-1), and the last coefficient
   is the one at index 0. Calling this function again undoes the reversal.
   Note that this function *must* be called a second time before calling
   zn_pmf_vec_clear(), so that free() is not called on the wrong pointer!
*/
#define zn_pmf_vec_reverse \
    ZNP_zn_pmf_vec_reverse
void zn_pmf_vec_reverse(zn_pmf_vec_t op, ulong length);



/* ============================================================================

     stuff in array.c

============================================================================ */


#define zn_skip_array_signed_add \
    ZNP_zn_skip_array_signed_add
ulong* zn_skip_array_signed_add(ulong* res, ptrdiff_t skip, size_t len,
                                const ulong* op1, int neg1,
                                const ulong* op2, int neg2,
                                const zn_mod_t mod);


/*
   Same as zn_array_scalar_mul, but has a _redc_ flag. If the flag is set,
   then REDC reduction is used (in which case the modulus must be odd),
   otherwise ordinary reduction is used.
*/
#define _zn_array_scalar_mul \
    ZNP__zn_array_scalar_mul
void _zn_array_scalar_mul(ulong* res, const ulong* op, size_t len,
                          ulong x, int redc, const zn_mod_t mod);



/* ============================================================================

     stuff in nussbaumer.c

============================================================================ */


#define nussbaumer_params \
    ZNP_nussbaumer_params
void nussbaumer_params(unsigned* lgK, unsigned* lgM, unsigned lgL);

#define nussbaumer_mul \
    ZNP_nussbaumer_mul
void nussbaumer_mul(ulong* res, const ulong* op1, const ulong* op2,
                    zn_pmf_vec_t vec1, zn_pmf_vec_t vec2);

#define nussbaumer_mul_get_fudge \
    ZNP_nussbaumer_mul_get_fudge
ulong nussbaumer_mul_get_fudge(unsigned lgL, int squaring, const zn_mod_t mod);



/* ============================================================================

     stuff from mul_fft.c

============================================================================ */


/*
   See mul_fft.c for documentation of the following FFT and IFFT functions.
*/

#define zn_pmf_vec_fft \
    ZNP_zn_pmf_vec_fft
void zn_pmf_vec_fft(zn_pmf_vec_t op, ulong length, ulong nonzero, ulong twist);

#define zn_pmf_vec_fft_factor \
    ZNP_zn_pmf_vec_fft_factor
void zn_pmf_vec_fft_factor(zn_pmf_vec_t op, unsigned lgT,
                           ulong length, ulong nonzero, ulong twist);

#define zn_pmf_vec_fft_small \
    ZNP_zn_pmf_vec_fft_small
void zn_pmf_vec_fft_small(zn_pmf_vec_t op, ulong length, ulong nonzero,
                          ulong twist);

#define zn_pmf_vec_fft_notrunc_iterative \
    ZNP_zn_pmf_vec_fft_notrunc_iterative
void zn_pmf_vec_fft_notrunc_iterative(zn_pmf_vec_t op, ulong twist);


#define zn_pmf_vec_ifft \
    ZNP_zn_pmf_vec_ifft
void zn_pmf_vec_ifft(zn_pmf_vec_t op, ulong length, int forward,
                     ulong nonzero, ulong twist);

#define zn_pmf_vec_ifft_factor \
    ZNP_zn_pmf_vec_ifft_factor
void zn_pmf_vec_ifft_factor(zn_pmf_vec_t op, unsigned lgT, ulong length,
                            int forward, ulong nonzero, ulong twist);

#define zn_pmf_vec_ifft_small \
    ZNP_zn_pmf_vec_ifft_small
void zn_pmf_vec_ifft_small(zn_pmf_vec_t op, ulong length, int forward,
                           ulong nonzero, ulong twist);

#define zn_pmf_vec_ifft_notrunc_iterative \
    ZNP_zn_pmf_vec_ifft_notrunc_iterative
void zn_pmf_vec_ifft_notrunc_iterative(zn_pmf_vec_t op, ulong twist);



/*
   Same as zn_array_mul(), but uses the Schonhage/Nussbaumer FFT algorithm.
   
   Uses faster algorithm for squaring if inputs are identical buffers.

   The modulus must be odd.

   Output may overlap the inputs.

   The output will come out divided by a fudge factor, which can be recovered
   via zn_array_mul_fft_get_fudge().
   
   If scale != 1, the output is further multiplied by _scale_.
*/
#define zn_array_mul_fft \
    ZNP_zn_array_mul_fft
void zn_array_mul_fft(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, ulong scale,
                      const zn_mod_t mod);

#define zn_array_mul_fft_get_fudge \
    ZNP_zn_array_mul_fft_get_fudge
ulong zn_array_mul_fft_get_fudge(size_t len1, size_t len2, int squaring,
                                 const zn_mod_t mod);


/*
   Computes the best lgK, lgM, coeffs1, coeffs2 such that polynomials of
   length len1 and len2 may be multiplied with fourier transform parameters
   lgK and lgM, and where the polynomials get split into coeffs1 (resp.
   coeffs2) chunks of length M/2.
   
   More precisely, the outputs satisfy:
   
   * lgM + 1 >= lgK (i.e. there are enough roots of unity)
   * coeffs1 + coeffs2 - 1 <= K,
     where coeffs1 = ceil(len1 / (M/2)) and coeffs2 = ceil(len2 / (M/2))
     (i.e. the transform has enough room for the answer)
   * lgM >= 1
   * lgM is minimal subject to the above conditions.
*/
#define mul_fft_params \
    ZNP_mul_fft_params
void mul_fft_params(unsigned* lgK, unsigned* lgM,
                    ulong* coeffs1, ulong* coeffs2,
                    size_t len1, size_t len2);



/*
   Splits op[-lead, len) into pieces of length M/2, where M = res->M, and
   where the first _lead_ coefficients are assumed to be zero. The pieces
   are written to the first ceil((len + lead) / (M/2)) coefficients of res.
   The last fractional piece is treated as if zero-padded up to length M/2.
   The second half of each target zn_pmf_t is zeroed out, and the bias fields
   are all set to _bias_.
   
   If scale != 1, then all entries are multiplied by _scale_ mod n.
*/
#define fft_split \
    ZNP_fft_split
void fft_split(zn_pmf_vec_t res, const ulong* op, size_t len, size_t lead,
               ulong scale, ulong bias);


/*
   Performs the substitution back from S[Z]/(Z^K - 1) to a polynomial in X,
   i.e. mapping Y -> X, Z -> X^(M/2). It only looks at the first _nonzero_
   coefficients of op; it assumes the rest are zero. It writes exactly
   _len_ coefficients of output. If skip_first is set, it ignores the first
   M/2 coefficients of output, and then writes the *next* M/2 coefficients
   of output.
   
   NOTE: this routine is not threadsafe: it temporarily modifies the bias
   field of the first coefficient of op.
*/
#define fft_combine \
    ZNP_fft_combine
void fft_combine(ulong* res, size_t len, const zn_pmf_vec_t op, ulong nonzero,
                 int skip_first);


/* ============================================================================

     stuff from midmul_fft.c

============================================================================ */


/*
   See midmul_fft.c for documentation of the following transposed FFT and
   IFFT functions.
*/

#define zn_pmf_vec_fft_transposed \
    ZNP_zn_pmf_vec_fft_transposed
void zn_pmf_vec_fft_transposed(zn_pmf_vec_t op, ulong length, ulong nonzero,
                               ulong twist);

#define zn_pmf_vec_fft_transposed_factor \
    ZNP_zn_pmf_vec_fft_transposed_factor
void zn_pmf_vec_fft_transposed_factor(zn_pmf_vec_t op, unsigned lgT,
                                      ulong length, ulong nonzero,
                                      ulong twist);

#define zn_pmf_vec_fft_transposed_small \
    ZNP_zn_pmf_vec_fft_transposed_small
void zn_pmf_vec_fft_transposed_small(zn_pmf_vec_t op, ulong length,
                                     ulong nonzero, ulong twist);

#define zn_pmf_vec_fft_transposed_notrunc_iterative \
    ZNP_zn_pmf_vec_fft_transposed_notrunc_iterative
void zn_pmf_vec_fft_transposed_notrunc_iterative(zn_pmf_vec_t op, ulong twist);


#define zn_pmf_vec_ifft_transposed \
    ZNP_zn_pmf_vec_ifft_transposed
void zn_pmf_vec_ifft_transposed(zn_pmf_vec_t op, ulong length, int forward,
                                ulong nonzero, ulong twist);

#define zn_pmf_vec_ifft_transposed_factor \
    ZNP_zn_pmf_vec_ifft_transposed_factor
void zn_pmf_vec_ifft_transposed_factor(zn_pmf_vec_t op, unsigned lgT,
                                       ulong length, int forward,
                                       ulong nonzero, ulong twist);

#define zn_pmf_vec_ifft_transposed_small \
    ZNP_zn_pmf_vec_ifft_transposed_small
void zn_pmf_vec_ifft_transposed_small(zn_pmf_vec_t op, ulong length,
                                      int forward, ulong nonzero, ulong twist);

#define zn_pmf_vec_ifft_transposed_notrunc_iterative \
    ZNP_zn_pmf_vec_ifft_transposed_notrunc_iterative
void zn_pmf_vec_ifft_transposed_notrunc_iterative(zn_pmf_vec_t op,
                                                  ulong twist);


/*
   Stores precomputed information for performing an FFT-based middle product
   where the first input array is invariant, and the length of the second
   input array is invariant.
*/
#define zn_array_midmul_fft_precomp1_struct \
    ZNP_zn_array_midmul_fft_precomp1_struct
struct zn_array_midmul_fft_precomp1_struct
{
   // these parameters are as described at the top of midmul_fft.c.
   size_t len1, len2;
   ulong coeffs1, coeffs2, pad;
   
   // holds the transposed IFFT of the input array
   zn_pmf_vec_t vec1;
};

#define zn_array_midmul_fft_precomp1_t \
    ZNP_zn_array_midmul_fft_precomp1_t
typedef struct zn_array_midmul_fft_precomp1_struct
               zn_array_midmul_fft_precomp1_t[1];


/*
   Initialises res to perform middle product of op1[0, len1) by operands of
   size len2.

   If scale != 1, the data is multiplied by _scale_. Since middle products
   are linear, this has the effect of multiplying the output of subsequent
   calls to zn_array_midmul_fft_precomp1_execute() by _scale_.
*/
#define zn_array_midmul_fft_precomp1_init \
    ZNP_zn_array_midmul_fft_precomp1_init
void zn_array_midmul_fft_precomp1_init(
           zn_array_midmul_fft_precomp1_t res, const ulong* op1,
           size_t len1, size_t len2, ulong scale, const zn_mod_t mod);

/*
   Performs middle product of op1[0, len1) by op2[0, len2), stores result
   at res[0, len1 - len2 + 1).
   
   The output will come out divided by a fudge factor, which can be recovered
   via zn_array_midmul_fft_precomp1_get_fudge().
   
   If scale != 1, the output is further multiplied by _scale_.
*/
#define zn_array_midmul_fft_precomp1_execute \
    ZNP_zn_array_midmul_fft_precomp1_execute
void zn_array_midmul_fft_precomp1_execute(
            ulong* res, const ulong* op2, ulong scale,
            const zn_array_midmul_fft_precomp1_t precomp);

#define zn_array_midmul_fft_precomp1_get_fudge \
    ZNP_zn_array_midmul_fft_precomp1_get_fudge
ulong zn_array_midmul_fft_precomp1_get_fudge(size_t len1, size_t len2,
                                             const zn_mod_t mod);


/*
   Deallocates op.
*/
#define zn_array_midmul_fft_precomp1_clear \
    ZNP_zn_array_midmul_fft_precomp1_clear
void zn_array_midmul_fft_precomp1_clear(zn_array_midmul_fft_precomp1_t op);


/*
   Same as zn_array_midmul(), but uses the Schonhage/Nussbaumer FFT algorithm.

   The modulus must be odd.

   Output may overlap the inputs.
   
   The output will come out divided by a fudge factor, which can be recovered
   via zn_array_midmul_fft_get_fudge().
   
   If scale != 1, the output is further multiplied by _scale_.
*/
#define zn_array_midmul_fft \
    ZNP_zn_array_midmul_fft
void zn_array_midmul_fft(ulong* res, const ulong* op1, size_t len1,
                         const ulong* op2, size_t len2,
                         ulong scale, const zn_mod_t mod);

#define zn_array_midmul_fft_get_fudge \
    ZNP_zn_array_midmul_fft_get_fudge
ulong zn_array_midmul_fft_get_fudge(size_t len1, size_t len2,
                                    const zn_mod_t mod);


/*
   Computes the best lgK, lgM, coeffs1, coeffs2, pad such that the middle
   product of polynomials of length len1 and len2 may be computed with
   fourier transform parameters lgK and lgM, where the first polynomial is
   padded on the left by _pad_ zeroes, and where the polynomials get split
   into coeffs1 (resp. coeffs2) chunks of length M/2.
   
   More precisely, the outputs satisfy (see midmul_fft.c for further
   discussion):
   
   * lgM >= 1
   * lgM + 1 >= lgK (i.e. there are enough roots of unity)
   * 1 <= pad <= M/2  and  len2 + pad - 1 is divisible by M/2
   * coeffs1 = ceil((len1 + pad) / (M/2))
   * coeffs2 = ceil(len2 / (M/2))
   * coeffs1 <= K
   * lgM is minimal subject to the above conditions.
*/
#define midmul_fft_params \
    ZNP_midmul_fft_params
void midmul_fft_params(unsigned* lgK, unsigned* lgM,
                       ulong* coeffs1, ulong* coeffs2, ulong* pad,
                       size_t len1, size_t len2);



/* ============================================================================

     stuff from midmul.c

============================================================================ */


/*
   Same as zn_array_midmul(), but always falls back on doing the full product
   via _zn_array_mul().
   
   If fastred is cleared, the output is the same as for zn_array_midmul().
   
   If fastred is set, the routine uses the fastest modular reduction strategy
   available for the given parameters. The result will come out divided by a
   fudge factor, which can be recovered via
   _zn_array_midmul_fallback_get_fudge().
*/
#define zn_array_midmul_fallback \
    ZNP_zn_array_midmul_fallback
void zn_array_midmul_fallback(ulong* res, const ulong* op1, size_t len1,
                              const ulong* op2, size_t len2, int fastred,
                              const zn_mod_t mod);

#define zn_array_midmul_fallback_get_fudge \
    ZNP_zn_array_midmul_fallback_get_fudge
ulong zn_array_midmul_fallback_get_fudge(size_t len1, size_t len2,
                                         const zn_mod_t mod);


/*
   Identical to zn_array_midmul(), except for the _fastred_ flag.
   
   If fastred is cleared, the output is the same as for zn_array_midmul().
   
   If fastred is set, the routine uses the fastest modular reduction strategy
   available for the given parameters. The result will come out divided by a
   fudge factor, which can be recovered via _zn_array_midmul_get_fudge().
*/
#define _zn_array_midmul \
    ZNP__zn_array_midmul
void _zn_array_midmul(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int fastred,
                      const zn_mod_t mod);
                      
#define _zn_array_midmul_get_fudge \
    ZNP__zn_array_midmul_get_fudge
ulong _zn_array_midmul_get_fudge(size_t len1, size_t len2, const zn_mod_t mod);



/* ============================================================================

     other random stuff

============================================================================ */


/*
   Compute 2^k mod n.
   
   Modulus must be odd.
   
   Must have -ULONG_BITS < k < ULONG_BITS.
*/
#define zn_mod_pow2 \
    ZNP_zn_mod_pow2
ulong zn_mod_pow2(int k, const zn_mod_t mod);



/*
   Constants describing algorithms for precomputed middle products
*/
#define ZNP_MIDMUL_ALGO_FALLBACK 0
#define ZNP_MIDMUL_ALGO_FFT 1



#ifdef __cplusplus
}
#endif

#endif

// end of file ****************************************************************
