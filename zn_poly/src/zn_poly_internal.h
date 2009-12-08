/*
   zn_poly_internal.h:  main header file #included internally by zn_poly
                        modules
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
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
 


#include "../../mpn_extras.h"
#define ZNP_mpn_mul F_mpn_mul



/*
   Executes "stuff", and checks that it returns 0. Usage is e.g.:
   
      ZNP_ASSERT_NOCARRY (mpn_add_n (z, x, y, 42));

   This is used where we want the test suite to check that the return value
   is zero, but the production version should skip the check, and we don't
   want to stuff around with temporaries everywhere.
*/
#define ZNP_ASSERT_NOCARRY(stuff)   \
    do {                            \
       mp_limb_t __xyz_cy;          \
       __xyz_cy = (stuff);          \
       ZNP_ASSERT (__xyz_cy == 0);  \
    } while (0)


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
int ceil_lg (ulong x);


/*
   Returns floor(log2(x)).
   Returns -1 for x == 0.
*/
#define floor_lg \
    ZNP_floor_lg
int
floor_lg (ulong x);


/*
   res := abs(op1 - op2).

   Returns 1 if op1 - op2 is negative, else zero.
*/
#define signed_mpn_sub_n \
    ZNP_signed_mpn_sub_n
ZNP_INLINE int
signed_mpn_sub_n (mp_limb_t* res, const mp_limb_t* op1, const mp_limb_t* op2,
                  size_t n)
{
   if (mpn_cmp (op1, op2, n) >= 0)
   {
      mpn_sub_n (res, op1, op2, n);
      return 0;
   }
   else
   {
      mpn_sub_n (res, op2, op1, n);
      return 1;
   }
}


/*
   The ZNP_FASTALLOC and ZNP_FASTFREE macros are used for allocating memory
   which is taken off the stack if the request is small enough, or off the
   heap if not.
   
   Example usage:
   
   ZNP_FASTALLOC (stuff, int, 100, n);

   This does two things. It allocates an array of 100 ints on the stack.
   It also declares a pointer "int* stuff", which points to a block of ints
   of length n. If n <= 100, the block will be the one just allocated on the
   stack. If n > 100, the block will be found using malloc.
   
   Then afterwards, you need to do:
   
   ZNP_FASTFREE (stuff);
   
   This will call free() if the block was originally taken off the heap.
*/

#define ZNP_FASTALLOC(ptr, type, reserve, request)        \
   size_t __FASTALLOC_request_##ptr = (request);          \
   type* ptr;                                             \
   type __FASTALLOC_##ptr [reserve];                      \
   if (__FASTALLOC_request_##ptr <= (reserve))            \
      ptr = __FASTALLOC_##ptr;                            \
   else                                                   \
      ptr = (type*) malloc (sizeof (type) * __FASTALLOC_request_##ptr);
      

#define ZNP_FASTFREE(ptr)                                 \
   if (ptr != __FASTALLOC_##ptr)                          \
      free (ptr);



extern size_t
ZNP_mpn_smp_kara_thresh;

extern size_t
ZNP_mpn_mulmid_fallback_thresh;


/*
   Stores tuning data for moduli of a specific bitsize.
*/
#define tuning_info_t \
    ZNP_tuning_info_t
typedef struct
{
   // thresholds for array multiplication
   size_t mul_KS2_thresh;    // KS1 -> KS2 threshold
   size_t mul_KS4_thresh;    // KS2 -> KS4 threshold
   size_t mul_fft_thresh;    // KS4 -> fft threshold

   // as above, but for squaring
   size_t sqr_KS2_thresh;
   size_t sqr_KS4_thresh;
   size_t sqr_fft_thresh;

   // as above, but for middle products
   size_t mulmid_KS2_thresh;
   size_t mulmid_KS4_thresh;
   size_t mulmid_fft_thresh;
   
   // for negacyclic multiplications, switch from KS to Nussbaumer FFT
   // when length reaches 2^nuss_mul_thresh
   unsigned nuss_mul_thresh;
   // ditto for nussbaumer squaring
   unsigned nuss_sqr_thresh;
   
}
tuning_info_t;


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

      res := 2^k * ( op[0] + op[s]*2^b + ... + op[(n-1)*s]*2^((n-1)*b) ).
                              
   Assumes each op[i] satisfies 0 <= op[i] < 2^b.

   Must have 0 < b < 3 * ULONG_BITS.
   
   If r == 0, then exactly
       ceil((k + n * b) / GMP_NUMB_BITS)
   limbs are written. Otherwise, the output will be zero-padded up to exactly
   r limbs, which must be at least the above number of limbs.
   
*/
#define zn_array_pack \
    ZNP_zn_array_pack
void
zn_array_pack (mp_limb_t* res, const ulong* op, size_t n, ptrdiff_t s,
               unsigned b, unsigned k, size_t r);


/*
   Let op be an integer of the form

      2^k * (a[0] + a[1]*2^b + ... + a[n-1]*2^((n-1)*b)) + junk,

   where 0 <= a[i] < 2^b for each i, and where 0 <= junk < 2^k.
   
   This function reads off the a[i]'s and stores them at res. Each output
   coefficient occupies exactly ceil(b / ULONG_BITS) words. The input
   should be exactly ceil((k + n * b) / GMP_NUMB_BITS) limbs long.
   
   Must have 0 < b < 3 * ULONG_BITS.
*/
#define zn_array_unpack \
    ZNP_zn_array_unpack
void
zn_array_unpack (ulong* res, const mp_limb_t* op, size_t n, unsigned b,
                 unsigned k);


/*
   Same as zn_array_unpack, but adds an assertion to check that the unpacking
   routine will not read beyond the first r limbs of op.
*/
#define zn_array_unpack_SAFE(res, op, n, b, k, r)                     \
do                                                                    \
{                                                                     \
   ZNP_ASSERT((n) * (b) + (k) <= (r) * GMP_NUMB_BITS);                \
   zn_array_unpack(res, op, n, b, k);                                 \
} while (0)



/* ============================================================================

     stuff from mul.c

============================================================================ */

/*
   Identical to zn_array_mul(), except for the fastred flag.
   
   If fastred is cleared, the output is the same as for zn_array_mul().
   
   If fastred is set, the routine uses the fastest modular reduction strategy
   available for the given parameters. The result will come out divided by a
   fudge factor, which can be recovered via _zn_array_mul_fudge().
*/
#define _zn_array_mul \
    ZNP__zn_array_mul
void
_zn_array_mul (ulong* res,
               const ulong* op1, size_t n1,
               const ulong* op2, size_t n2,
               int fastred, const zn_mod_t mod);

#define _zn_array_mul_fudge \
    ZNP__zn_array_mul_fudge
ulong
_zn_array_mul_fudge (size_t n1, size_t n2, int sqr, const zn_mod_t mod);



/* ============================================================================

     stuff from ks_support.c

============================================================================ */


/*
   Sets res[i * s] = reduction modulo mod of the i-th entry of op,
   for 0 <= i < n. Each entry of op is w ulongs.
   
   If the redc flag is set, the results are divided by -B mod m
   (only allowed if the modulus is odd).
   
   Must have 1 <= w <= 3.
*/
#define array_reduce \
    ZNP_array_reduce
void
array_reduce (ulong* res, ptrdiff_t s, const ulong* op, size_t n, unsigned w,
              int redc, const zn_mod_t mod);



/*
   This is a helper function for the variants of KS that evaluate at
   "reciprocal" evaluation points like 2^(-N); it implements essentially the
   algorithm of section 3.2 of [Har07], plus reductions mod n.

   It accepts two integers X and Y written in base M = 2^b, where
   1 <= b <= 3 * ULONG_BITS / 2. It assumes that

      X = a[0] + a[1]*M + ... + a[n-1]*M^(n-1),
      Y = a[n-1] + a[n-2]*M + ... + a[0]*M^(n-1),

   where each a[i] is two "digits" long, and where the high digit of a[i]
   is at most M-2 (i.e. may not equal M-1). It reconstructs the a[i],
   reduces them mod m, and stores the results in an array.
   
   The input is supplied as follows. X is in op1, Y is in op2. They are both
   arrays of values that are b bits wide, where b <= 3 * ULONG_BITS / 2.
   Each value takes up one ulong if b <= ULONG_BITS, otherwise two ulongs.
   There are n + 1 such values in each array (i.e. each array consists of
   (n + 1) * ceil(b / ULONG_BITS) ulongs).
   
   The output (n ulongs) is written to the array res, with consecutive
   outputs separated by s ulongs.
   
   mod describes the modulus m.

   If the redc flag is set, the modular reductions are performed using REDC,
   i.e. the result contain an extra factor of -1/B mod m (where
   B = 2^ULONG_BITS).
*/

#define zn_array_recover_reduce \
    ZNP_zn_array_recover_reduce
void
zn_array_recover_reduce (ulong* res, ptrdiff_t s, const ulong* op1,
                         const ulong* op2, size_t n, unsigned b, int redc,
                         const zn_mod_t mod);



/* ============================================================================

     stuff from mul_ks.c

============================================================================ */

/*
   These are the same as zn_array_mul(). They use four different types of
   Kronecker substitution. They automatically use a faster algorithm for
   squaring (if the inputs are identical buffers).

   Aliasing of all operands allowed.

   Must have n1 >= n2 >= 1.
   
   If the redc flag is set, the outputs will be divided by -B mod m.
   (Only allowed if the modulus is odd.)
*/
#define zn_array_mul_KS1 \
    ZNP_zn_array_mul_KS1
void
zn_array_mul_KS1 (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int redc, const zn_mod_t mod);

#define zn_array_mul_KS2 \
    ZNP_zn_array_mul_KS2
void
zn_array_mul_KS2 (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int redc, const zn_mod_t mod);

#define zn_array_mul_KS3 \
    ZNP_zn_array_mul_KS3
void
zn_array_mul_KS3 (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int redc, const zn_mod_t mod);

#define zn_array_mul_KS4 \
    ZNP_zn_array_mul_KS4
void
zn_array_mul_KS4 (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int redc, const zn_mod_t mod);


/* ============================================================================

     stuff from mulmid_ks.c

============================================================================ */


/*
   These are the same as zn_array_mulmid(). They use four different types of
   Kronecker substitution.

   Aliasing of all operands allowed.

   Must have n1 >= n2 >= 1.
   
   If the redc flag is set, the outputs will be divided by -B mod m.
   (Only allowed if the modulus is odd.)
*/
#define zn_array_mulmid_KS1 \
    ZNP_zn_array_mulmid_KS1
void
zn_array_mulmid_KS1 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod);

#define zn_array_mulmid_KS2 \
    ZNP_zn_array_mulmid_KS2
void
zn_array_mulmid_KS2 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod);

#define zn_array_mulmid_KS3 \
    ZNP_zn_array_mulmid_KS3
void
zn_array_mulmid_KS3 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod);

#define zn_array_mulmid_KS4 \
    ZNP_zn_array_mulmid_KS4
void
zn_array_mulmid_KS4 (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     int redc, const zn_mod_t mod);



/* ============================================================================

     pmf_t stuff

============================================================================ */


/*
   Let R = Z/mZ. A pmf_t ("pmf" = "polynomial modulo fermat") represents an
   element of S = R[Y]/(Y^M + 1). This is used as the coefficient ring in
   the Schonhage/Nussbaumer FFT routines.

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

   Currently the values a_i are always normalised into [0, m). Later we might
   drop that restriction to obtain faster butterflies...

*/
#define pmf_t \
    ZNP_pmf_t
typedef ulong*  pmf_t;

#define pmf_const_t \
    ZNP_pmf_const_t
typedef const ulong*  pmf_const_t;


/*
   op := op * x
*/
#define pmf_scalar_mul \
    ZNP_pmf_scalar_mul
ZNP_INLINE void
pmf_scalar_mul (pmf_t op, ulong M, ulong x, const zn_mod_t mod)
{
   zn_array_scalar_mul (op + 1, op + 1, M, x, mod);
}


/*
   op := 0, with bias reset to zero too
*/
#define pmf_zero \
    ZNP_pmf_zero
ZNP_INLINE void
pmf_zero (pmf_t op, ulong M)
{
   for (M++; M > 0; M--)
      *op++ = 0;
}


/*
   op := op / 2

   Modulus must be odd.
*/
#define pmf_divby2 \
    ZNP_pmf_divby2
ZNP_INLINE void
pmf_divby2 (pmf_t op, ulong M, const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);

   for (op++; M > 0; M--, op++)
      *op = zn_mod_divby2 (*op, mod);
}


/*
   res := op
*/
#define pmf_set \
    ZNP_pmf_set
ZNP_INLINE void
pmf_set (pmf_t res, pmf_t op, ulong M)
{
   for (M++; M > 0; M--)
      *res++ = *op++;
}


/*
   op := Y^r * op
*/
#define pmf_rotate \
    ZNP_pmf_rotate
ZNP_INLINE void
pmf_rotate (pmf_t op, ulong r)
{
   op[0] += r;
}


/*
   op1 := op2 + op1
   op2 := op2 - op1

   Inputs must be [0, m); outputs will be in [0, m).
*/
#define pmf_bfly \
    ZNP_pmf_bfly
void
pmf_bfly (pmf_t op1, pmf_t op2, ulong M, const zn_mod_t mod);

/*
   op1 := op1 + op2

   Inputs must be [0, m); outputs will be in [0, m).
*/
#define pmf_add \
    ZNP_pmf_add
void
pmf_add (pmf_t op1, const pmf_t op2, ulong M, const zn_mod_t mod);

/*
   op1 := op1 - op2

   Inputs must be [0, m); outputs will be in [0, m).
*/
#define pmf_sub \
    ZNP_pmf_sub
void
pmf_sub (pmf_t op1, const pmf_t op2, ulong M, const zn_mod_t mod);




/*
   These functions are exported just for profiling purposes:
*/

#define zn_array_bfly_inplace \
    ZNP_zn_array_bfly_inplace
void
zn_array_bfly_inplace (ulong* op1, ulong* op2, ulong n, const zn_mod_t mod);

#define zn_array_add_inplace \
    ZNP_zn_array_add_inplace
void
zn_array_add_inplace (ulong* op1, const ulong* op2, ulong n,
                      const zn_mod_t mod);

#define zn_array_sub_inplace \
    ZNP_zn_array_sub_inplace
void
zn_array_sub_inplace (ulong* op1, const ulong* op2, ulong n,
                      const zn_mod_t mod);



/* ============================================================================

     pmfvec_t stuff

============================================================================ */


/*
   A pmfvec_t stores a vector of length K = 2^lgK of elements of S.
   
   Used to represent an element of S[Z]/(Z^K + 1) or S[Z]/(Z^K - 1), or some
   other quotient like that.
   
   The functions pmfvec_init/clear should be used to allocate storage for
   this type. Also sometimes fake ones get created temporarily to point at
   sub-vectors of existing vectors.
*/
#define pmfvec_struct \
    ZNP_pmfvec_struct
typedef struct
{
   // points to the first coefficient
   pmf_t data;
   
   // number of coefficients
   ulong K;
   unsigned lgK;      // lg2(K)
   
   // length of coefficients (see definition of pmf_t)
   ulong M;
   unsigned lgM;      // lg2(M)

   // distance between adjacent coefficients, measured in ulongs
   // (this is at least M + 1, might be more)
   ptrdiff_t skip;
   
   // associated modulus
   const zn_mod_struct* mod;
}
pmfvec_struct;

#define pmfvec_t \
    ZNP_pmfvec_t
typedef pmfvec_struct  pmfvec_t[1];


/*
   Checks that vec1 and vec2 have compatible data,
   i.e. have the same K, M, mod.
*/
#define pmfvec_compatible \
    ZNP_pmfvec_compatible
ZNP_INLINE int
pmfvec_compatible (const pmfvec_t vec1, const pmfvec_t vec2)
{
   return (vec1->K == vec2->K) && (vec1->M == vec2->M) &&
          (vec1->mod == vec2->mod);
}


/*
   Initialises res with given parameters, allocates memory.
*/
#define pmfvec_init \
    ZNP_pmfvec_init
void
pmfvec_init (pmfvec_t res, unsigned lgK, ptrdiff_t skip, unsigned lgM,
             const zn_mod_t mod);


/*
   Initialises res in preparation for a Nussbaumer multiplication of
   length 2^lgL.
*/
#define pmfvec_init_nuss \
    ZNP_pmfvec_init_nuss
void
pmfvec_init_nuss (pmfvec_t res, unsigned lgL, const zn_mod_t mod);


/*
   Destroys op, frees all associated memory.
*/
#define pmfvec_clear \
    ZNP_pmfvec_clear
void
pmfvec_clear (pmfvec_t op);


/*
   res := op
*/
#define pmfvec_set \
    ZNP_pmfvec_set
void
pmfvec_set (pmfvec_t res, const pmfvec_t op);


/*
   Multiplies first n coefficients of op by x.
*/
#define pmfvec_scalar_mul \
    ZNP_pmfvec_scalar_mul
void
pmfvec_scalar_mul (pmfvec_t op, ulong n, ulong x);


/*
   Multiplies pointwise the first n coefficients of op1 and op2, puts result
   in res.
   
   It's okay for res to alias op1 or op2. The modulus must be odd.
   
   If the special_first_two flag is set, the routine assumes that the first
   two coefficients are of length only M/2 (this is the typical situation
   after performing the FFT), and multiplies them more quickly accordingly.
   
   The routine automatically selects KS or Nussbaumer multiplication
   depending on the modulus bitsize and on M.
   
   The output will be divided by a fudge factor, which can be retrieved
   via pmfvec_mul_fudge().
   
   Automatically uses specialised squaring algorithm if the inputs are the
   same pmfvec_t object.
*/
#define pmfvec_mul \
    ZNP_pmfvec_mul
void
pmfvec_mul (pmfvec_t res, const pmfvec_t op1, const pmfvec_t op2, ulong n,
            int special_first_two);

#define pmfvec_mul_fudge \
    ZNP_pmfvec_mul_fudge
ulong
pmfvec_mul_fudge (unsigned lgM, int sqr, const zn_mod_t mod);


/*
   Modifies the op->data and op->skip to make it look as if the first
   coefficient is the one at index n - 1, and the last coefficient is
   the one at index 0. Calling this function again undoes the reversal.
   Note that this function *must* be called a second time before calling
   pmfvec_clear(), so that free() is not called on the wrong pointer!
*/
#define pmfvec_reverse \
    ZNP_pmfvec_reverse
void
pmfvec_reverse (pmfvec_t op, ulong n);



/* ============================================================================

     stuff in pmfvec_fft.c

============================================================================ */

/*
   ---------- forward FFTs ----------

   The functions
   
      pmfvec_fft()
      pmfvec_fft_basecase()
      pmfvec_fft_dc()
      pmfvec_fft_huge()
      
   operate on a pmfvec_t, and compute inplace the FFT:

      b_k = Y^{t k'} \sum_{i=0}^{K-1} Y^{2 M i k' / K} a_i,

   where 0 <= t < 2M / K is a twist parameter. The notation k' indicates the
   length lgK bit-reversal of k.
   
   All except the "basecase" version have an n parameter; they only compute the
   first n outputs. The remaining buffers are used in intermediate steps, and
   contain junk at the end. For "basecase", all K outputs are computed.

   All except the "basecase" version have a z parameter; they assume that the
   input coefficients are zero from index z and beyond. They never read from
   those coefficients. For "basecase", all the inputs are used.
   
   The four versions use different algorithms as follows:
   
      * pmfvec_fft(): main entry point, delegates to one of the other routines
        based on the size of the transform.
        
      * pmfvec_fft_basecase(): low-overhead iterative FFT, no truncation logic.
        
      * pmfvec_fft_dc(): divide-and-conquer. It handles the top layer of
        butterflies, and then recurses into the two halves. This is intended
        for fairly small transforms where locality is not a big issue. The
        algorithm implemented here is essentially van der Hoeven's "truncated
        Fourier transform" [vdH04], [vdH05].
        
      * pmfvec_fft_huge(): intended for large transforms, where locality
        is an issue. It factors the FFT into U = 2^lgU transforms of length
        T = 2^lgT followed by T transforms of length U, where K = T * U. This
        is done recursively until we're in L1 cache (or as small as possible),
        at which point we switch to pmfvec_fft_dc().
        
        The algorithm is straightforward, but I believe it to be new. It is
        simultaneously a generalisation of van der Hoeven's truncated Fourier
        transform and Bailey's FFT algorithm [Bai89]. (I used something
        similar in the ZmodF_poly module in FLINT.)


   ---------- inverse FFTs ----------

   The functions

      pmfvec_ifft()
      pmfvec_ifft_basecase()
      pmfvec_ifft_dc()
      pmfvec_ifft_huge()

   compute the corresponding inverse FFT. They are a little more complicated
   than the FFTs owing to the truncation.
   
   Let a_i and b_k be as described above for the FFTs. The IFFT functions
   take as input the array

      b_{0'}, b_{1'}, ..., b_{(n-1)'}, K*a_n, K*a_{n+1}, ..., K*a_{K-1},
      
   for some 0 <= n <= K. If the fwd flag is zero, then the output of the IFFT
   routine is:
   
      K*a_0, K*a_1, ..., K*a_{n-1}
      
   followed by K - n junk coefficients. If the fwd flag is set, then we
   require that 0 <= n <= K - 1, and the output is

      K*a_0, K*a_1, ..., K*a_{n-1}, b_n
      
   followed by K - n - 1 junk coefficients; i.e. it also computes one
   coefficient of the *forward* FFT.

   The "basecase" version has no n parameter, and assumes that n = K. Here
   the routine becomes the (non-truncated) IFFT as usually understood, with
   inputs in bit-reversed order and outputs in normal order.

   All except the "basecase" version have a z parameter (with z >= n); they
   assume that the input coefficients are zero from index z and beyond. They
   never read from those coefficients. For the "basecase" version, all of the
   inputs are used.
   
   The four versions use different algorithms as follows:
   
      * pmfvec_ifft(): main entry point, delegates to one of the other routines
        based on the size of the transform.
        
      * pmfvec_ifft_basecase(): low-overhead iterative IFFT, no truncation
        logic.
        
      * pmfvec_ifft_dc(): divide-and-conquer. It recurses into the two halves,
        and handles the top layer of butterflies. This is intended for fairly
        small transforms where locality is not a big issue. The algorithm
        implemented here is essentially van der Hoeven's beautiful "truncated
        inverse Fourier transform".
        
      * pmfvec_ifft_huge(): intended for large transforms, where locality
        is an issue. It factors the FFT into U = 2^lgU transforms of length
        T = 2^lgT and T transforms of length U, where K = T * U. This is done
        recursively until we're in L1 cache (if possible), at which point we
        switch to pmfvec_ifft_dc().
        
        The algorithm is not as simple as the FFT version; it is necessary
        to alternate between "row" and "column" transforms in a slightly
        complicated way. I believe the algorithm to be new.


   ---------- transposed forward and inverse FFTs ----------

   The functions

      pmfvec_tpfft_basecase()
      pmfvec_tpfft_dc()
      pmfvec_tpfft_huge()
      pmfvec_tpfft()
      
   are *transposed* versions of the corresponding pmfvec_fft() routines. This
   means that if the FFT computes an S-linear map from S^z to S^n, the
   transposed version computes the transpose map from S^n to S^z.
   
   Similarly, the functions
       
      pmfvec_tpifft_basecase()
      pmfvec_tpifft_dc()
      pmfvec_tpifft_huge()
      pmfvec_tpifft()
      
   are transposed versions of the IFFT routines. If the IFFT computes an
   S-linear map from S^z to S^(n + fwd), the transposed version computes the
   transpose map from S^(n + fwd) to S^z.
   
   The algorithms are transposed essentially by reversing them, and
   transposing every step of the algorithm; see for example [BLS03] for how
   this is done. We don't have comments in these routines; see the comments
   on the corresponding FFT/IFFT routines.


   ---------- notes for all the above functions ----------

   These functions all perform O(n * lgK + K) operations in S (an "operation"
   being an addition/subtraction/copy with possibly an implied rotation by a
   power of Y). In particular the running time varies fairly smoothly with
   n instead of jumping with K.

   For all four algorithms, the "dc" version is essentially equivalent to the
   "huge" version with lgT = 1.

   These functions are not thread-safe. Apart from modifying the input inplace,
   they also temporarily modify the pmfvec_t structs themselves.

   Our approach to improving cache performance is certainly not ideal. The
   main problem is that we can get address conflicts, especially since
   everything gets spread out by powers of two. Mitigating factors:
   associativity in the cache; the extra bias word scrambles the addresses
   somewhat; when the transforms gets large, so do the coefficients, so we
   don't expect to fit that many in cache anyway.

*/


#define pmfvec_fft \
    ZNP_pmfvec_fft
void
pmfvec_fft (pmfvec_t op, ulong n, ulong z, ulong t);

#define pmfvec_fft_huge \
    ZNP_pmfvec_fft_huge
void
pmfvec_fft_huge (pmfvec_t op, unsigned lgT, ulong n, ulong z, ulong t);

#define pmfvec_fft_dc \
    ZNP_pmfvec_fft_dc
void
pmfvec_fft_dc (pmfvec_t op, ulong n, ulong z, ulong t);

#define pmfvec_fft_basecase \
    ZNP_pmfvec_fft_basecase
void
pmfvec_fft_basecase (pmfvec_t op, ulong t);

#define pmfvec_ifft \
    ZNP_pmfvec_ifft
void
pmfvec_ifft (pmfvec_t op, ulong n, int fwd, ulong z, ulong t);

#define pmfvec_ifft_huge \
    ZNP_pmfvec_ifft_huge
void
pmfvec_ifft_huge (pmfvec_t op, unsigned lgT, ulong n, int fwd, ulong z,
                  ulong t);

#define pmfvec_ifft_dc \
    ZNP_pmfvec_ifft_dc
void
pmfvec_ifft_dc (pmfvec_t op, ulong n, int fwd, ulong z, ulong t);

#define pmfvec_ifft_basecase \
    ZNP_pmfvec_ifft_basecase
void
pmfvec_ifft_basecase (pmfvec_t op, ulong t);

#define pmfvec_tpfft \
    ZNP_pmfvec_tpfft
void
pmfvec_tpfft (pmfvec_t op, ulong n, ulong z, ulong t);

#define pmfvec_tpfft_huge \
    ZNP_pmfvec_tpfft_huge
void
pmfvec_tpfft_huge (pmfvec_t op, unsigned lgT, ulong n, ulong z, ulong t);

#define pmfvec_tpfft_dc \
    ZNP_pmfvec_tpfft_dc
void
pmfvec_tpfft_dc (pmfvec_t op, ulong n, ulong z, ulong t);

#define pmfvec_tpfft_basecase \
    ZNP_pmfvec_tpfft_basecase
void
pmfvec_tpfft_basecase (pmfvec_t op, ulong t);

#define pmfvec_tpifft \
    ZNP_pmfvec_tpifft
void
pmfvec_tpifft (pmfvec_t op, ulong n, int fwd, ulong z, ulong t);

#define pmfvec_tpifft_huge \
    ZNP_pmfvec_tpifft_huge
void
pmfvec_tpifft_huge (pmfvec_t op, unsigned lgT, ulong n, int fwd, ulong z,
                    ulong t);

#define pmfvec_tpifft_dc \
    ZNP_pmfvec_tpifft_dc
void
pmfvec_tpifft_dc (pmfvec_t op, ulong n, int fwd, ulong z, ulong t);

#define pmfvec_tpifft_basecase \
    ZNP_pmfvec_tpifft_basecase
void
pmfvec_tpifft_basecase (pmfvec_t op, ulong t);



/* ============================================================================

     stuff in array.c

============================================================================ */


/*
   Computes

      res = sign1*op1 + sign2*op2,

   where sign1 = -1 if neg1 is set, otherwise +1; ditto for sign2.
   
   op1 and op2 are arrays of length n.
   
   res is a staggered array, entries separated by s.
   
   Return value is res + s*n, i.e. points beyond the written array.
*/
#define zn_skip_array_signed_add \
    ZNP_zn_skip_array_signed_add
ulong*
zn_skip_array_signed_add (ulong* res, ptrdiff_t skip, size_t n,
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
void
_zn_array_scalar_mul (ulong* res, const ulong* op, size_t n, ulong x,
                      int redc, const zn_mod_t mod);


/*
   Behaves just like zn_array_scalar_mul, except it uses the obvious
   optimisation if x == 1.
*/
#define zn_array_scalar_mul_or_copy \
    ZNP_zn_array_scalar_mul_or_copy
void
zn_array_scalar_mul_or_copy (ulong* res, const ulong* op, size_t n,
                             ulong x, const zn_mod_t mod);


/* ============================================================================

     stuff in nuss.c

============================================================================ */


/*
   Performs negacyclic multiplication using Nussbaumer's algorithm.
   
   vec1 and vec2 must be pre-initialised pmfvec_t's with the same
   modulus and the same lgM and lgK, satisfying lgM + 1 >= lgK (i.e. there
   are enough roots of unity). These are used for scratch space.
   
   The convolution length is L = 2^lgL, where lgL = lgM + lgK - 1.

   Inputs are op1[0, L) and op2[0, L), output is res[0, L). It's okay for res
   to alias op1 or op2.
   
   If op1 == op2, then a faster squaring version is used. In this case
   vec2 is ignored.

   The result comes out divided by a fudge factor, which can be recovered
   via nuss_mul_fudge().
*/
#define nuss_mul \
    ZNP_nuss_mul
void
nuss_mul (ulong* res, const ulong* op1, const ulong* op2,
                pmfvec_t vec1, pmfvec_t vec2);

#define nuss_mul_fudge \
    ZNP_nuss_mul_fudge
ulong
nuss_mul_fudge (unsigned lgL, int sqr, const zn_mod_t mod);


/*
   Computes optimal lgK and lgM for given lgL, as described above for
   nuss_mul().
*/
#define nuss_params \
    ZNP_nuss_params
void
nuss_params (unsigned* lgK, unsigned* lgM, unsigned lgL);



/* ============================================================================

     stuff from mul_fft.c

============================================================================ */


/*
   Splits op[-k, n) into pieces of length M/2, where M = res->M, and where
   the first k coefficients are assumed to be zero. The pieces are written to
   the first ceil((n + k) / (M/2)) coefficients of res. The last fractional
   piece is treated as if zero-padded up to length M/2. The second half of
   each target pmf_t is zeroed out, and the bias fields are all set to b.
   
   If x != 1, then all entries are multiplied by x mod m.
*/
#define fft_split \
    ZNP_fft_split
void
fft_split (pmfvec_t res, const ulong* op, size_t n, size_t k, ulong x,
           ulong b);


/*
   Performs the substitution back from S[Z]/(Z^K - 1) to a polynomial in X,
   i.e. mapping Y -> X, Z -> X^(M/2). It only looks at the first z coefficients
   of op; it assumes the rest are zero. It writes exactly n coefficients of
   output. If skip_first is set, it ignores the first M/2 coefficients of
   output, and then writes the *next* M/2 coefficients of output.
   
   NOTE: this routine is not threadsafe: it temporarily modifies the bias
   field of the first coefficient of op.
*/
#define fft_combine \
    ZNP_fft_combine
void
fft_combine (ulong* res, size_t n, const pmfvec_t op, ulong z,
             int skip_first);



/*
   Same as zn_array_mul(), but uses the Schonhage/Nussbaumer FFT algorithm.
   
   Uses faster algorithm for squaring if inputs are identical buffers.

   The modulus must be odd.

   Output may overlap the inputs.

   The output will come out divided by a fudge factor, which can be recovered
   via zn_array_mul_fft_fudge().
   
   If x != 1, the output is further multiplied by x.
*/
#define zn_array_mul_fft \
    ZNP_zn_array_mul_fft
void
zn_array_mul_fft (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  ulong x, const zn_mod_t mod);

#define zn_array_mul_fft_fudge \
    ZNP_zn_array_mul_fft_fudge
ulong
zn_array_mul_fft_fudge (size_t n1, size_t n2, int sqr, const zn_mod_t mod);


/*
   Computes the best lgK, lgM, m1, m2 such that polynomials of length n1 and
   n2 may be multiplied with fourier transform parameters lgK and lgM, and
   where the polynomials get split into m1 (resp. m2) chunks of length M/2.
   
   More precisely, the outputs satisfy:
   
   * lgM + 1 >= lgK (i.e. there are enough roots of unity)
   * m1 + m2 - 1 <= K, where m1 = ceil(n1 / (M/2)) and m2 = ceil(n2 / (M/2))
     (i.e. the transform has enough room for the answer)
   * lgM >= 1
   * lgM is minimal subject to the above conditions.
*/
#define mul_fft_params \
    ZNP_mul_fft_params
void
mul_fft_params (unsigned* lgK, unsigned* lgM, ulong* m1, ulong* m2,
                size_t n1, size_t n2);


/*
   Same as zn_array_mulmid(), but uses the Schonhage/Nussbaumer FFT algorithm.

   The modulus must be odd.

   Output may overlap the inputs.
   
   The output will come out divided by a fudge factor, which can be recovered
   via zn_array_mulmid_fft_fudge().
   
   If x != 1, the output is further multiplied by x.
*/
#define zn_array_mulmid_fft \
    ZNP_zn_array_mulmid_fft
void
zn_array_mulmid_fft (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     ulong x, const zn_mod_t mod);

#define zn_array_mulmid_fft_fudge \
    ZNP_zn_array_mulmid_fft_fudge
ulong
zn_array_mulmid_fft_fudge (size_t n1, size_t n2, const zn_mod_t mod);


/*
   Computes the best lgK, lgM, m1, m2, p such that the middle product of
   polynomials of length n1 and n2 may be computed with fourier transform
   parameters lgK and lgM, where the first polynomial is padded on the left by
   p zeroes, and where the polynomials get split into m1 (resp. m2) chunks of
   length M/2.
   
   More precisely, the outputs satisfy (see mul_fft.c for further discussion):
   
   * lgM >= 1
   * lgM + 1 >= lgK (i.e. there are enough roots of unity)
   * 1 <= p <= M/2  and  n2 + p - 1 is divisible by M/2
   * m1 = ceil((n1 + p) / (M/2))
   * m2 = ceil(n2 / (M/2))
   * m1 <= K
   * lgM is minimal subject to the above conditions.
*/
#define mulmid_fft_params \
    ZNP_mulmid_fft_params
void
mulmid_fft_params (unsigned* lgK, unsigned* lgM, ulong* m1, ulong* m2,
                   ulong* p, size_t n1, size_t n2);


/*
   Stores precomputed information for performing an FFT-based middle product
   where the first input array is invariant, and the length of the second
   input array is invariant.
*/
#define zn_array_mulmid_fft_precomp1_struct \
    ZNP_zn_array_mulmid_fft_precomp1_struct
struct zn_array_mulmid_fft_precomp1_struct
{
   // these parameters are as described at the top of mul_fft.c.
   size_t n1, n2;
   ulong m1, m2, p;
   
   // holds the transposed IFFT of the input array
   pmfvec_t vec1;
};

#define zn_array_mulmid_fft_precomp1_t \
    ZNP_zn_array_mulmid_fft_precomp1_t
typedef struct zn_array_mulmid_fft_precomp1_struct
               zn_array_mulmid_fft_precomp1_t[1];


/*
   Initialises res to perform middle product of op1[0, n1) by operands of
   size n2.

   If x != 1, the data is multiplied by x. Since middle products are linear,
   this has the effect of multiplying the output of subsequent calls to
   zn_array_mulmid_fft_precomp1_execute() by x.
*/
#define zn_array_mulmid_fft_precomp1_init \
    ZNP_zn_array_mulmid_fft_precomp1_init
void
zn_array_mulmid_fft_precomp1_init (zn_array_mulmid_fft_precomp1_t res,
                                   const ulong* op1, size_t n1, size_t n2,
                                   ulong x, const zn_mod_t mod);

/*
   Performs middle product of op1[0, n1) by op2[0, n2), stores result at
   res[0, n1 - n2 + 1).
   
   The output will come out divided by a fudge factor, which can be recovered
   via zn_array_mulmid_fft_precomp1_fudge().
   
   If x != 1, the output is further multiplied by x.
*/
#define zn_array_mulmid_fft_precomp1_execute \
    ZNP_zn_array_mulmid_fft_precomp1_execute
void
zn_array_mulmid_fft_precomp1_execute
                     (ulong* res, const ulong* op2, ulong x,
                      const zn_array_mulmid_fft_precomp1_t precomp);

#define zn_array_mulmid_fft_precomp1_fudge \
    ZNP_zn_array_mulmid_fft_precomp1_fudge
ulong
zn_array_mulmid_fft_precomp1_fudge (size_t n1, size_t n2, const zn_mod_t mod);


/*
   Deallocates op.
*/
#define zn_array_mulmid_fft_precomp1_clear \
    ZNP_zn_array_mulmid_fft_precomp1_clear
void
zn_array_mulmid_fft_precomp1_clear (zn_array_mulmid_fft_precomp1_t op);



/* ============================================================================

     stuff from mpn_mulmid.c

============================================================================ */


/*
   Let n1 >= n2 >= 1, and let

      a = \sum_{i=0}^{n1-1} a_i B^i
      b = \sum_{j=0}^{n2-1} b_j B^j

   be integers with n1 and n2 limbs respectively. We define SMP(a, b), the
   *simple* middle product of a and b, to be the integer

      \sum_{0 <= i < n1, 0 <= j < n2, n2-1 <= i+j < n1} a_i b_j B^(i+j-(n2-1)).

   In other words, it's as if we treat a and b as polynomials in Z[B] of
   length n1 and n2 respectively, compute the polynomial middle product over
   Z, and then propagate the high words and subsequent carries.
   
   Note that SMP(a, b) is at most n1 - n2 + 3 limbs long (we assume throughout
   that n1 is less than the maximum value stored in a limb).
*/


/*
   Computes SMP(op1[0, n1), op2[0, n2)).
   
   Stores result at res[0, n1 - n2 + 3).
*/
void
ZNP_mpn_smp (mp_limb_t* res,
             const mp_limb_t* op1, size_t n1,
             const mp_limb_t* op2, size_t n2);


/*
   Same as mpn_smp(), but always uses basecase (quadratic-time) algorithm.
   
   res[0, n1 - n2 + 3) must not overlap op1[0, n1) or op2[0, n2).
*/
void
ZNP_mpn_smp_basecase (mp_limb_t* res,
                      const mp_limb_t* op1, size_t n1,
                      const mp_limb_t* op2, size_t n2);


/*
   Computes SMP(op1[0, 2*n - 1), op2[0, n)). Algorithm is selected depending
   on size of n.
   
   Stores result at res[0, n + 2). Output must not overlap inputs.
*/
void
ZNP_mpn_smp_n (mp_limb_t* res, const mp_limb_t* op1,
               const mp_limb_t* op2, size_t n);


/*
   Computes SMP(op1[0, 2*n - 1), op2[0, n)), using Karatsuba algorithm.
   
   Must have n >= 2.

   Stores result at res[0, n + 2). Output must not overlap inputs.
*/
void
ZNP_mpn_smp_kara (mp_limb_t* res, const mp_limb_t* op1, const mp_limb_t* op2,
                  size_t n);



/*
   Computes the *true* middle product of op1[0, n1) and op2[0, n2).
   
   More precisely, let P be the product op1 * op2. This function computes
   P[n2 + 1, n1), and stores this at res[2, n1 - n2 + 1).

   Must have n1 >= n2 >= 1.

   The output buffer res *must* have room for n1 - n2 + 3 limbs, but the first
   two limbs and the last two limbs of the output will be *garbage*. In
   particular, only n1 - n2 - 1 limbs of useful output are produced. If
   n1 <= n2 + 1, then no useful output is produced.
*/
void
ZNP_mpn_mulmid (mp_limb_t* res,
                const mp_limb_t* op1, size_t n1,
                const mp_limb_t* op2, size_t n2);


/*
   Same as mpn_mulmid, but always just falls back on using mpn_mul.
*/
void
ZNP_mpn_mulmid_fallback (mp_limb_t* res,
                         const mp_limb_t* op1, size_t n1,
                         const mp_limb_t* op2, size_t n2);




/* ============================================================================

     stuff from mulmid.c

============================================================================ */


/*
   Same as zn_array_mulmid(), but always falls back on doing the full product
   via _zn_array_mul().
   
   If fastred is cleared, the output is the same as for zn_array_mulmid().
   
   If fastred is set, the routine uses the fastest modular reduction strategy
   available for the given parameters. The result will come out divided by a
   fudge factor, which can be recovered via
   _zn_array_mulmid_fallback_fudge().
*/
#define zn_array_mulmid_fallback \
    ZNP_zn_array_mulmid_fallback
void
zn_array_mulmid_fallback (ulong* res,
                          const ulong* op1, size_t n1,
                          const ulong* op2, size_t n2,
                          int fastred, const zn_mod_t mod);

#define zn_array_mulmid_fallback_fudge \
    ZNP_zn_array_mulmid_fallback_fudge
ulong
zn_array_mulmid_fallback_fudge (size_t n1, size_t n2, const zn_mod_t mod);


/*
   Identical to zn_array_mulmid(), except for the fastred flag.
   
   If fastred is cleared, the output is the same as for zn_array_mulmid().
   
   If fastred is set, the routine uses the fastest modular reduction strategy
   available for the given parameters. The result will come out divided by a
   fudge factor, which can be recovered via _zn_array_mulmid_fudge().
*/
#define _zn_array_mulmid \
    ZNP__zn_array_mulmid
void
_zn_array_mulmid (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  int fastred, const zn_mod_t mod);

#define _zn_array_mulmid_fudge \
    ZNP__zn_array_mulmid_fudge
ulong
_zn_array_mulmid_fudge (size_t n1, size_t n2, const zn_mod_t mod);


/* ============================================================================

     other random stuff

============================================================================ */


/*
   Compute 2^k mod m.
   
   Modulus must be odd.
   
   Must have -ULONG_BITS < k < ULONG_BITS.
*/
#define zn_mod_pow2 \
    ZNP_zn_mod_pow2
ulong
zn_mod_pow2 (int k, const zn_mod_t mod);



/*
   Constants describing algorithms for precomputed middle products
*/
#define ZNP_MULMID_ALGO_FALLBACK  0
#define ZNP_MULMID_ALGO_KS        1
#define ZNP_MULMID_ALGO_FFT       2



#ifdef __cplusplus
}
#endif

#endif

// end of file ****************************************************************
