/*
   zn_poly.h:  main header file to be #included by zn_poly users
   
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

#ifndef ZN_POLY_H
#define ZN_POLY_H

#ifdef __cplusplus
extern "C" {
#endif
 

#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <assert.h>


/*
   why not, let's have our own assertion macro just to be difficult
*/
#define ZNP_ASSERT assert


/*
   hmmm.... I'm not totally sure if static inline is portable enough.
*/
#define ZNP_INLINE static inline


/*
   Returns a string like "3.1"
*/
extern const char* zn_poly_version_string();


/*
   Three components of "version x.y.z"
*/
#define ZNP_VERSION_MAJOR 0
#define ZNP_VERSION_MINOR 8
#define ZNP_VERSION_REVISION 0


/*
   ULONG_BITS = number of bits per unsigned long
*/
#if ULONG_MAX == 4294967295U
#define ULONG_BITS 32
#elif ULONG_MAX == 18446744073709551615U
#define ULONG_BITS 64
#else
#error zn_poly requires that unsigned long is either 32 bits or 64 bits
#endif

#include "wide_arith.h"

#include "../flint.h"


/* ============================================================================

     zn_mod_t stuff

============================================================================ */

/*
   zn_mod_t stores precomputed information about a modulus.
   
   The modulus can be any integer in the range 2 <= n < 2^ULONG_BITS.
   
   A modulus n is called "slim" if n <= 2^(ULONG_BITS - 1), i.e. the residues
   never occupy the top bit of the word. Many routines are much faster for
   slim moduli.
*/
typedef struct
{
   // the modulus, must be >= 2
   ulong n;
   
   // ceil(log2(n)) = number of bits in a non-negative residue
   int bits;
   
   // reduction of B and B^2 mod n (where B = 2^ULONG_BITS)
   ulong B, B2;
   
   // sh1 and inv1 are respectively ell-1 and m'
   // from Figure 4.1 of [GM94]
   unsigned sh1;
   ulong inv1;
   
   // sh2, sh3, inv2 and n_norm are respectively N-ell, ell-1, m', d_norm
   // from Figure 8.1 of [GM94]
   unsigned sh2, sh3;
   ulong inv2, n_norm;
   
   // inv3 = n^(-1) mod B (only valid if n is odd)
   ulong inv3;
}
zn_mod_struct;

typedef zn_mod_struct zn_mod_t[1];


/*
   Initialises zn_mod_t with given modulus, does some (fairly cheap)
   precomputations.
*/
void zn_mod_init(zn_mod_t mod, ulong n);


/*
   Must be called when the modulus object goes out of scope.
*/
void zn_mod_clear(zn_mod_t mod);


/*
   Return nonzero if mod is a slim modulus.
*/
ZNP_INLINE
int zn_mod_is_slim(const zn_mod_t mod)
{
   return (long) mod->n >= 0;
}


/*
   Returns op1 + op2 mod n.
   
   Both op1 and op2 must be in [0, n).
*/
ZNP_INLINE
ulong zn_mod_add(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   ulong temp = mod->n - op2;
   if (op1 < temp)
      return op1 + op2;
   else
      return op1 - temp;
}


/*
   Same as zn_mod_add, but only for slim moduli.
   
   This is usually several times faster than zn_mod_add, depending on the
   context; see the bfly-profile target for examples.
*/
ZNP_INLINE
ulong zn_mod_add_slim(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(zn_mod_is_slim(mod));
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   ulong temp = op1 + op2;
   if (temp >= mod->n)
      temp -= mod->n;
   return temp;
}


/*
   Returns op1 - op2 mod n.
   
   Both op1 and op2 must be in [0, n).
*/
ZNP_INLINE
ulong zn_mod_sub(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   ulong x = op1 - op2;
   if (op1 < op2)
      x += mod->n;
   return x;
}


/*
   Same as zn_mod_sub, but only for slim moduli.
*/
ZNP_INLINE
ulong zn_mod_sub_slim(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(zn_mod_is_slim(mod));
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   long temp = op1 - op2;
   temp += (temp < 0) ? mod->n : 0;
   return temp;
}


/*
   Returns -op mod n.

   op must be in [0, n).
*/
ZNP_INLINE
ulong zn_mod_neg(ulong op, const zn_mod_t mod)
{
   ZNP_ASSERT(op < mod->n);
   return op ? (mod->n - op) : op;
}


/*
   Return op/2 mod n.

   op must be in [0, n).
   
   If the modulus is even, op must be even too.
*/
ZNP_INLINE
ulong zn_mod_divby2(ulong op, const zn_mod_t mod)
{
   ZNP_ASSERT(op < mod->n);
   ZNP_ASSERT((mod->n & 1) || !(op & 1));

   return (op >> 1) + ((-(op & 1)) & ((mod->n >> 1) + 1));
}


/*
   Returns floor(x / n).

   No restrictions on x.
   
   Algorithm is essentially Figure 4.1 of [GM94].
*/
ZNP_INLINE
ulong zn_mod_quotient(ulong x, const zn_mod_t mod)
{
   ulong t;
   ZNP_MUL_HI(t, x, mod->inv1);
   return (t + ((x - t) >> 1)) >> mod->sh1;
}


/*
   Returns x mod n.
   
   No restrictions on x.
*/
ZNP_INLINE
ulong zn_mod_reduce(ulong x, const zn_mod_t mod)
{
   return x - zn_mod_quotient(x, mod) * mod->n;
}


/*
   Returns -x/B mod n.
   
   n must be odd.
   
   No restrictions on x.
*/
ZNP_INLINE
ulong zn_mod_reduce_redc(ulong x, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);

   ulong y = x * mod->inv3;
   ulong z;
   ZNP_MUL_HI(z, y, mod->n);
   return z;
}



/*
   Returns x1*B + x0 mod n.
   
   Assumes x1 is already in [0, n).

   Algorithm is essentially Figure 8.1 of [GM94].
*/
ZNP_INLINE
ulong zn_mod_reduce_wide(ulong x1, ulong x0, const zn_mod_t mod)
{
   ZNP_ASSERT(x1 < mod->n);

   ulong y1 = (x1 << mod->sh2) + ((x0 >> 1) >> mod->sh3);
   ulong y0 = (x0 << mod->sh2);
   
   ulong sign = y0 >> (ULONG_BITS - 1);
   ulong z0 = y0 + (mod->n_norm & -sign);
   
   ulong a1, a0;
   ZNP_MUL_WIDE(a1, a0, mod->inv2, y1 + sign);
   ZNP_ADD_WIDE(a1, a0, a1, a0, y1, z0);
   
   ulong b1, b0;
   ZNP_MUL_WIDE(b1, b0, (-a1 - 1), mod->n);
   ZNP_ADD_WIDE(b1, b0, b1, b0, x1, x0);
   b1 -= mod->n;
   
   return b0 + (b1 & mod->n);
}


/*
   Returns -(x1*B + x0)/B mod n.
   
   Assumes x1 is already in [0, n), and that n is odd.
   
   Uses essentially Montgomery's REDC algorithm [Mon85].
*/
ZNP_INLINE
ulong zn_mod_reduce_wide_redc(ulong x1, ulong x0, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(x1 < mod->n);
   
   ulong y = x0 * mod->inv3;
   ulong z;
   ZNP_MUL_HI(z, y, mod->n);
   return zn_mod_sub(z, x1, mod);
}


/*
   Returns -(x1*B + x0)/B mod n.
   
   Assumes x1 is already in [0, n), and that n is odd, and that n is slim.
   
   Uses essentially Montgomery's REDC algorithm [Mon85].
*/
ZNP_INLINE
ulong zn_mod_reduce_wide_redc_slim(ulong x1, ulong x0, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(x1 < mod->n);
   
   ulong y = x0 * mod->inv3;
   ulong z;
   ZNP_MUL_HI(z, y, mod->n);
   return zn_mod_sub_slim(z, x1, mod);
}


/*
   Returns x1*B + x0 mod n.

   No restrictions on x0 and x1.
*/
ZNP_INLINE
ulong zn_mod_reduce2(ulong x1, ulong x0, const zn_mod_t mod)
{
   // first reduce into [0, Bn)
   ulong c0, c1;
   ZNP_MUL_WIDE(c1, c0, x1, mod->B);
   ZNP_ADD_WIDE(c1, c0, c1, c0, 0, x0);
   // (must still have c1 < n)

   return zn_mod_reduce_wide(c1, c0, mod);
}



/*
   Returns -(x1*B + x0)/B mod n.
   
   n must be odd.

   No restrictions on x0 and x1.
*/
ZNP_INLINE
ulong zn_mod_reduce2_redc(ulong x1, ulong x0, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);

   // first reduce into [0, Bn)
   ulong c0, c1;
   ZNP_MUL_WIDE(c1, c0, x1, mod->B);
   ZNP_ADD_WIDE(c1, c0, c1, c0, 0, x0);
   // (must still have c1 < n)

   return zn_mod_reduce_wide_redc(c1, c0, mod);
}


/*
   Returns x2*B^2 + x1*B + x0 mod n.

   No restrictions on x0, x1 or x2.
*/
ZNP_INLINE
ulong zn_mod_reduce3(ulong x2, ulong x1, ulong x0, const zn_mod_t mod)
{
   // reduce B^2*x2 and B*x1 into [0, Bn)
   ulong c0, c1, d0, d1;
   ZNP_MUL_WIDE(c1, c0, x2, mod->B2);
   ZNP_MUL_WIDE(d1, d0, x1, mod->B);
   
   // add B^2*x2 and B*x1 and x0 mod Bn
   ZNP_ADD_WIDE(c1, c0, c1, c0, 0, d0);
   // (must still have c1 < n)
   ZNP_ADD_WIDE(c1, c0, c1, c0, 0, x0);
   if (c1 >= mod->n)
      c1 -= mod->n;

   c1 = zn_mod_add(c1, d1, mod);

   // finally reduce it mod n
   return zn_mod_reduce_wide(c1, c0, mod);
}



/*
   Returns -(x2*B^2 + x1*B + x0)/B mod n.

   n must be odd.

   No restrictions on x0, x1 or x2.
*/
ZNP_INLINE
ulong zn_mod_reduce3_redc(ulong x2, ulong x1, ulong x0, const zn_mod_t mod)
{
   // reduce B^2*x2 and B*x1 into [0, Bn)
   ulong c0, c1, d0, d1;
   ZNP_MUL_WIDE(c1, c0, x2, mod->B2);
   ZNP_MUL_WIDE(d1, d0, x1, mod->B);
   
   // add B^2*x2 and B*x1 and x0 mod Bn
   ZNP_ADD_WIDE(c1, c0, c1, c0, 0, d0);
   // (must still have c1 < n)
   ZNP_ADD_WIDE(c1, c0, c1, c0, 0, x0);
   if (c1 >= mod->n)
      c1 -= mod->n;

   c1 = zn_mod_add(c1, d1, mod);

   // finally reduce it mod n
   return zn_mod_reduce_wide_redc(c1, c0, mod);
}



/*
   Returns op1 * op2 mod n.
   
   op1 and op2 must be in [0, n).
*/
ZNP_INLINE
ulong zn_mod_mul(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   ulong hi, lo;
   ZNP_MUL_WIDE(hi, lo, op1, op2);
   return zn_mod_reduce_wide(hi, lo, mod);
}


/*
   Returns -(op1 * op2)/B mod n.
   
   op1 and op2 must be in [0, n), and n must be odd.
*/
ZNP_INLINE
ulong zn_mod_mul_redc(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   ulong hi, lo;
   ZNP_MUL_WIDE(hi, lo, op1, op2);
   return zn_mod_reduce_wide_redc(hi, lo, mod);
}


/*
   Returns -(op1 * op2)/B mod n.
   
   op1 and op2 must be in [0, n), and n must be odd and slim.
*/
ZNP_INLINE
ulong zn_mod_mul_redc_slim(ulong op1, ulong op2, const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(op1 < mod->n && op2 < mod->n);

   ulong hi, lo;
   ZNP_MUL_WIDE(hi, lo, op1, op2);
   return zn_mod_reduce_wide_redc_slim(hi, lo, mod);
}


/*
   Returns x^k mod n.
   
   x must be in [0, n).
   
   Negative indices are not supported (yet).
*/
ulong zn_mod_pow(ulong x, long k, const zn_mod_t mod);


/*
   Returns 1/x mod n, or 0 if x is not invertible mod n.
   
   x must be in [0, n).
*/
ulong zn_mod_invert(ulong x, const zn_mod_t mod);



/* ============================================================================

     scalar multiplication on raw arrays
     
============================================================================ */

/*
   Multiplies each element of op[0, len) by x, stores result at res[0, len).
   
   res and op must be either identical or disjoint buffers.
*/
void zn_array_scalar_mul(ulong* res, const ulong* op, size_t len,
                         ulong x, const zn_mod_t mod);



/* ============================================================================

     polynomial multiplication on raw arrays
     
============================================================================ */


/*
   Multiplies op1[0, len1) by op2[0, len2), stores result in
   res[0, len1 + len2 - 1).
   
   op1 and op2 may alias each other, but neither may overlap res.
   
   Must have len1 >= len2 >= 1.
   
   Automatically selects best multiplication algorithm based on modulus
   and input lengths. Automatically uses specialised squaring code if inputs
   buffers are identical.
*/
void zn_array_mul(ulong* res, const ulong* op1, size_t len1,
                  const ulong* op2, size_t len2, const zn_mod_t mod);


/*
   Middle product of op1[0, len1) and op2[0, len2), stores result in
   res[0, len1 - len2 + 1).
   
   (i.e. this is the subarray of the ordinary product op1 * op2 consisting of
   those coefficients with indices in the range [len2 - 1, len1).)

   Must have len1 >= len2 >= 1.

   op1 and op2 may alias each other, but neither may overlap res.

   Performance note: for large inputs (in the FFT range), we use a "correct"
   implementation of the middle product (i.e. a 2n*n middle product takes the
   same time as an n*n full product). For small inputs (in the KS range), we
   currently use an ordinary product and extract the relevant coefficients;
   this will be improved in the future.
*/
void zn_array_midmul(ulong* res, const ulong* op1, size_t len1,
                     const ulong* op2, size_t len2, const zn_mod_t mod);


// forward declaration (see zn_poly_internal.h)
struct ZNP_zn_array_midmul_fft_precomp1_struct;


/*
   Stores precomputed information for performing a middle product where the
   first input array op1[0, len1) is invariant, and the *length* of the second
   input array op2[0, len2) is invariant.
*/
typedef struct
{
   // Determines which middle product algorithm we're using.
   // One of the constants:
   //    ZNP_MIDMUL_ALGO_FALLBACK: fall back on zn_array_midmul
   //    ZNP_MIDMUL_ALGO_FFT: use zn_array_midmul_fft_precomp1
   int algo;

   size_t len1, len2;
   const zn_mod_struct* mod;

   // stores a copy of op1[0, len1) if we're using ZNP_MIDMUL_ALGO_FALLBACK
   ulong* op1;
   
   // precomputed data if we're using ZNP_MIDMUL_ALGO_FFT
   struct ZNP_zn_array_midmul_fft_precomp1_struct* precomp_fft;
}
zn_array_midmul_precomp1_struct;

typedef zn_array_midmul_precomp1_struct zn_array_midmul_precomp1_t[1];


/*
   Initialises res to perform middle product of op1[0, len1) by operands of
   size len2.
*/
void zn_array_midmul_precomp1_init(zn_array_midmul_precomp1_t res,
                                   const ulong* op1, size_t len1,
                                   size_t len2, const zn_mod_t mod);

/*
   Performs middle product of op1[0, len1) by op2[0, len2), stores result
   at res[0, len1 - len2 + 1).
*/
void zn_array_midmul_precomp1_execute(
            ulong* res, const ulong* op2,
            const zn_array_midmul_precomp1_t precomp);


/*
   Deallocates op.
*/
void zn_array_midmul_precomp1_clear(zn_array_midmul_precomp1_t op);



/*
   Same as zn_array_mul(), but uses the Schonhage/Nussbaumer FFT algorithm,
   with a few layers of naive DFT to save memory.
   
   lgT is the number of layers of DFT. Larger values of lgT save more memory,
   as long as lgT doesn't get too close to lg2(sqrt(len1 + len2)). Larger
   values also make the function slower. Probably you never want to make lgT
   bigger than 4; after that the savings are marginal.

   The modulus must be odd.

   Output may *not* overlap inputs.
   
   NOTE: this interface is preliminary and may change in future versions.
*/
void zn_array_mul_fft_dft(ulong* res, const ulong* op1, size_t len1,
                          const ulong* op2, size_t len2, unsigned lgT,
                          const zn_mod_t mod);



/* ============================================================================

     polynomial division on raw arrays
     
============================================================================ */


/*
   Computes len terms of power series inverse of op[0, len).
   
   Must have len >= 1.
   
   Must have op[0] == 1 (todo: this will be fixed later)
   
   Output may not overlap input.
*/
void zn_array_invert(ulong* res, const ulong* op, size_t len,
                     const zn_mod_t mod);



/* ============================================================================

     other miscellaneous zn_array stuff

============================================================================ */


/*
   res := -op
   
   Inputs and outputs in [0, n).
*/
void zn_array_neg(ulong* res, const ulong* op, size_t len, const zn_mod_t mod);


/*
   res := op1 - op2.
   
   Inputs and outputs in [0, n).
*/
void zn_array_sub(ulong* res, const ulong* op1, const ulong* op2, size_t len,
                  const zn_mod_t mod);


/*
   Returns zero if op1[0, len) and op2[0, len) are equal, otherwise nonzero.
*/
int zn_array_cmp(const ulong* op1, const ulong* op2, size_t len);


/*
   Copies op[0, len) to res[0, len).
   
   Buffers must not overlap.
*/
void zn_array_copy(ulong* res, const ulong* op, size_t len);


/*
   Sets res[0, len) to zero.
*/
ZNP_INLINE
void zn_array_zero(ulong* res, size_t len)
{
   for (; len; len--)
      *res++ = 0;
}



#ifdef __cplusplus
}
#endif

#endif

// end of file ****************************************************************
