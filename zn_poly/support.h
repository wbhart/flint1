/*
   support.h:  various support routines for test, profiling and tuning code
   
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

#ifndef ZNP_SUPPORT_H
#define ZNP_SUPPORT_H


#include <stdio.h>
#include <gmp.h>
#include "zn_poly_internal.h"


#ifdef __cplusplus
extern "C" {
#endif



/*
   single global random state for test/profile modules
*/
extern gmp_randstate_t randstate;



/*
   An array of modulus bitsizes, used by several test functions.
*/
extern unsigned test_bitsizes[];
extern unsigned num_test_bitsizes;    // how big the array is


/*
   Exports abs(op) to res, storing exactly len limbs
   (zero-padded if necessary).

   Sign of op is ignored.

   abs(op) must fit into len limbs.
*/
void mpz_to_mpn(mp_limb_t* res, size_t len, const mpz_t op);


/*
   Converts mpn buffer (exactly len limbs) to mpz.

   Output is always non-negative.
*/
void mpn_to_mpz(mpz_t res, const mp_limb_t* op, size_t len);


/*
   Returns random unsigned long in [0, max).
*/
ulong random_ulong(ulong max);


/*
   Returns random unsigned long in [0, 2^bits).
*/
ulong random_ulong_bits(unsigned bits);


/*
   Returns random modulus with exactly _bits_ bits, i.e. in the range
   [2^(bits-1), 2^bits).
   
   If require_odd is set, the returned modulus will be odd.
*/
ulong random_modulus(unsigned bits, int require_odd);


/*
   Prints array to stdout, in format e.g. "[2 3 7]".
*/
void zn_array_print(const ulong* x, size_t len);



void ref_zn_array_mul(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, const zn_mod_t mod);

void ref_zn_array_scalar_mul(ulong* res, const ulong* op, size_t len,
                             ulong x, const zn_mod_t mod);

void ref_zn_array_midmul(ulong* res, const ulong* op1, size_t len1,
                         const ulong* op2, size_t len2, const zn_mod_t mod);

void ref_zn_array_negamul(ulong* res, const ulong* op1, const ulong* op2,
                          size_t len, const zn_mod_t mod);


#if DEBUG

/*
   Prints op to standard output (in normalised form).
*/
void zn_pmf_print(const zn_pmf_t op, ulong M, const zn_mod_t mod);

/*
   Prints op to standard output.
*/
void zn_pmf_vec_print(const zn_pmf_vec_t op);

/*
   Prints first _length_ coefficients of op to standard output.
*/
void zn_pmf_vec_print_trunc(const zn_pmf_vec_t op, ulong length);


#endif



/* ============================================================================

     tuning routines

============================================================================ */


#define tune_mul_KS \
    ZNP_tune_mul_KS
void tune_mul_KS(FILE* flog, int squaring, int verbose);

#define tune_mul_nussbaumer \
    ZNP_tune_mul_nussbaumer
void tune_nussbaumer(FILE* flog, int squaring, int verbose);

#define tune_mul \
    ZNP_tune_mul
void tune_mul(FILE* flog, int squaring, int verbose);




/* ============================================================================

     structs used in profiling routines

============================================================================ */


enum
{
   ALGO_MUL_BEST,
   ALGO_MUL_KS1,
   ALGO_MUL_KS1_REDC,
   ALGO_MUL_KS2,
   ALGO_MUL_KS2_REDC,
   ALGO_MUL_KS3,
   ALGO_MUL_KS3_REDC,
   ALGO_MUL_KS4,
   ALGO_MUL_KS4_REDC,
   ALGO_MUL_FFT,
   ALGO_MUL_NTL,
};

/*
   struct for passing information to profile_mul
*/
typedef struct
{
   size_t len;         // length of polynomials to multiply
   ulong n;            // the modulus
   int algo;           // one of the ALGO_MUL_* values
   int squaring;       // whether to profile squaring or multiplication
}
profile_mul_info_struct;

typedef profile_mul_info_struct profile_mul_info_t[1];

/*
   Profiles one of the multiplication routines.
   
   _arg_ should point to a profile_mul_info_t describing what to profile.

   Returns total cycle count for _count_ calls.
*/
double profile_mul(void* arg, unsigned long count);

/*
   As above, but assumes that the algorithm is ALGO_MUL_NTL.
*/
double profile_mul_ntl(void* arg, unsigned long count);



enum
{
   // fall back on calling zn_array_mul and reducing negacyclically
   ALGO_NEGAMUL_FALLBACK,

   // use Nussbaumer convolution
   ALGO_NEGAMUL_NUSSBAUMER,
};

/*
   struct for passing information to profile_negamul
*/
typedef struct
{
   unsigned lg_len;    // lg2 of polynomial length
   ulong n;            // the modulus
   int algo;           // one of the ALGO_NEGAMUL_* values
   int squaring;       // whether to profile squaring or multiplication
}
profile_negamul_info_struct;

typedef profile_negamul_info_struct profile_negamul_info_t[1];

/*
   Profiles one of the negacyclic multiplication routines.
   
   _arg_ should point to a profile_negamul_info_t describing what to profile.

   Returns total cycle count for _count_ calls.
*/
double profile_negamul(void* arg, unsigned long count);



enum
{
   ALGO_MIDMUL_BEST,
   ALGO_MIDMUL_FALLBACK,
   ALGO_MIDMUL_KS1,
   ALGO_MIDMUL_KS2,
   ALGO_MIDMUL_KS3,
   ALGO_MIDMUL_KS4,
   ALGO_MIDMUL_FFT,
};

/*
   struct for passing information to profile_midmul
*/
typedef struct
{
   size_t len;         // we're doing a (2*len) x len middle product
   ulong n;            // the modulus
   int algo;           // one of the ALGO_MIDMUL_* values
}
profile_midmul_info_struct;

typedef profile_midmul_info_struct profile_midmul_info_t[1];

/*
   Profiles one of the middle product routines.
   
   _arg_ should point to a profile_midmul_info_t describing what to profile.

   Returns total cycle count for _count_ calls.
*/
double profile_midmul(void* arg, unsigned long count);



enum
{
   ALGO_INVERT_BEST,
   ALGO_INVERT_NTL,
};

/*
   struct for passing information to profile_invert
*/
typedef struct
{
   size_t len;         // we're doing a length len inversion
   ulong n;            // the modulus
   int algo;           // one of the ALGO_INVERT_* values
}
profile_invert_info_struct;

typedef profile_invert_info_struct profile_invert_info_t[1];

/*
   Profiles one of the series inversion routines.
   
   _arg_ should point to a profile_invert_info_t describing what to profile.

   Returns total cycle count for _count_ calls.
*/
double profile_invert(void* arg, unsigned long count);

/*
   As above, but assumes that the algorithm is ALGO_INVERT_NTL.
*/
double profile_invert_ntl(void* arg, unsigned long count);




double profile_bfly(void* arg, unsigned long count);
double profile_mpn_aors(void* arg, unsigned long count);
double profile_scalar_mul(void* arg, unsigned long count);

void prof_main(int argc, char* argv[]);


#ifdef __cplusplus
}
#endif

#endif

// end of file ****************************************************************
