/*
   support.h:  various support routines for test, profiling and tuning code
   
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
   Exports abs(op) to res, storing exactly n limbs (zero-padded if necessary).

   Sign of op is ignored.

   abs(op) must fit into n limbs.
*/
void
mpz_to_mpn (mp_limb_t* res, size_t n, const mpz_t op);


/*
   Converts mpn buffer (exactly n limbs) to mpz.

   Output is always non-negative.
*/
void
mpn_to_mpz (mpz_t res, const mp_limb_t* op, size_t n);


/*
   Returns random unsigned long in [0, max).
*/
ulong
random_ulong (ulong max);


/*
   Returns random unsigned long in [0, 2^b).
*/
ulong
random_ulong_bits (unsigned b);


/*
   Returns random modulus with exactly b bits, i.e. in [2^(b-1), 2^b).
   
   If require_odd is set, the returned modulus will be odd.
*/
ulong
random_modulus (unsigned b, int require_odd);


/*
   Prints array to stdout, in format e.g. "[2 3 7]".
*/
void
zn_array_print (const ulong* x, size_t n);


/*
   Similar to mpn_random2, but flips all the bits with probability 1/2.
   This deals with the annoying way that mpn_random2 always writes 1's to
   the high bits of the buffer.
*/
void
ZNP_mpn_random2 (mp_limb_t* res, size_t n);




void
ref_zn_array_mul (ulong* res,
                  const ulong* op1, size_t n1,
                  const ulong* op2, size_t n2,
                  const zn_mod_t mod);

void
ref_zn_array_scalar_mul (ulong* res, const ulong* op, size_t n, ulong x,
                         const zn_mod_t mod);

void
ref_zn_array_mulmid (ulong* res,
                     const ulong* op1, size_t n1,
                     const ulong* op2, size_t n2,
                     const zn_mod_t mod);

void
ref_zn_array_negamul (ulong* res, const ulong* op1, const ulong* op2,
                      size_t n, const zn_mod_t mod);


void
ref_mpn_mulmid (mp_limb_t* res, const mp_limb_t* op1, size_t n1,
                const mp_limb_t* op2, size_t n2);
                

void
ref_mpn_smp (mp_limb_t* res,
             const mp_limb_t* op1, size_t n1,
             const mp_limb_t* op2, size_t n2);


#if DEBUG


/*
   Sets res to a uniformly random pmf. Bias is uniformly random in [0, 2M).
*/
void
pmf_rand (pmf_t res, ulong M, const zn_mod_t mod);

/*
   Compares op1 and op2, returns 0 if equal.
*/
int
pmf_cmp (const pmf_t op1, const pmf_t op2, ulong M, const zn_mod_t mod);

/*
   Prints op to standard output (in normalised form).
*/
void
pmf_print (const pmf_t op, ulong M, const zn_mod_t mod);

/*
   Prints op to standard output.
*/
void
pmfvec_print (const pmfvec_t op);

/*
   Prints first n coefficients of op to standard output.
*/
void
pmfvec_print_trunc (const pmfvec_t op, ulong n);


#endif



/* ============================================================================

     tuning routines

============================================================================ */


#define tune_mul_KS \
    ZNP_tune_mul_KS
void
tune_mul_KS (FILE* flog, int sqr, int verbose);

#define tune_mulmid_KS \
    ZNP_tune_mulmid_KS
void
tune_mulmid_KS (FILE* flog, int verbose);

#define tune_mul_nuss \
    ZNP_tune_mul_nuss
void
tune_nuss (FILE* flog, int sqr, int verbose);

#define tune_mul \
    ZNP_tune_mul
void
tune_mul (FILE* flog, int sqr, int verbose);

#define tune_mulmid \
    ZNP_tune_mulmid
void
tune_mulmid (FILE* flog, int verbose);

#define tune_mpn_smp_kara \
    ZNP_tune_mpn_smp_kara
void
tune_mpn_smp_kara (FILE* flog, int verbose);

#define tune_mpn_mulmid_fallback \
    ZNP_tune_mpn_mulmid_fallback
void
tune_mpn_mulmid_fallback (FILE* flog, int verbose);



/* ============================================================================

     structs used in profiling routines

============================================================================ */

/*
   Struct for passing parameters to various profiling routines. Not all
   members are used by all routines, and they may have different meanings
   for different routines.
*/
typedef struct
{
   // modulus
   ulong m;
   
   // length of input polynomials
   size_t n;
   // lengths of input polynomials for routines taking two input lengths
   size_t n1, n2;
   // for negacyclic multiplication, log2 of the convolution length
   unsigned lgL;

   // which algorithm to use. Meaning depends on routine selected.
   int algo;
   
   // for routines profiling multiplication, indicates whether we want to
   // profile squaring
   int sqr;
}
profile_info_struct;

typedef profile_info_struct  profile_info_t[1];



/*
   legal algo values for profile_mul
*/
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
   Profiles one of the multiplication routines.
   
   arg points to a profile_info_t with parameters m, n1, n2, sqr, algo.

   Returns total cycle count for count calls.
*/
double
profile_mul (void* arg, unsigned long count);

/*
   As above, but assumes that the algorithm is ALGO_MUL_NTL.
*/
double
profile_mul_ntl (void* arg, unsigned long count);



/*
   legal algo values for profile_negamul
*/
enum
{
   // fall back on calling zn_array_mul and reducing negacyclically
   ALGO_NEGAMUL_FALLBACK,

   // use Nussbaumer convolution
   ALGO_NEGAMUL_NUSS,
};

/*
   Profiles one of the negacyclic multiplication routines.
   
   arg points to a profile_info_t with parameters m, lgL, sqr, algo.

   Returns total cycle count for count calls.
*/
double
profile_negamul (void* arg, unsigned long count);



/*
   legal algo values for profile_mulmid
*/
enum
{
   ALGO_MULMID_BEST,
   ALGO_MULMID_FALLBACK,
   ALGO_MULMID_KS1,
   ALGO_MULMID_KS1_REDC,
   ALGO_MULMID_KS2,
   ALGO_MULMID_KS2_REDC,
   ALGO_MULMID_KS3,
   ALGO_MULMID_KS3_REDC,
   ALGO_MULMID_KS4,
   ALGO_MULMID_KS4_REDC,
   ALGO_MULMID_FFT,
};

/*
   Profiles one of the middle product routines.
   
   arg points to a profile_info_t with parameters m, n1, n2, algo.

   Returns total cycle count for count calls.
*/
double
profile_mulmid (void* arg, unsigned long count);



/*
   legal algo values for profile_invert
*/
enum
{
   ALGO_INVERT_BEST,
   ALGO_INVERT_NTL,
};

/*
   Profiles one of the series inversion routines.
   
   arg points to a profile_info_t with parameters m, n, algo.

   Returns total cycle count for count calls.
*/
double
profile_invert (void* arg, unsigned long count);

/*
   As above, but assumes that the algorithm is ALGO_INVERT_NTL.
*/
double
profile_invert_ntl (void* arg, unsigned long count);




/*
   Profiles mpn_mul.
   
   arg points to a profile_info_t with parameters n1, n2.

   Returns total cycle count for count calls.
*/
double
profile_mpn_mul (void* arg, unsigned long count);

/*
   Profiles mpn_smp.
   
   arg points to a profile_info_t with parameters n1, n2.

   Returns total cycle count for count calls.
*/
double
profile_mpn_smp (void* arg, unsigned long count);

/*
   As above, for mpn_mulmid_fallback.
*/
double
profile_mpn_mulmid_fallback (void* arg, unsigned long count);

/*
   As above, for mpn_smp_basecase.
*/
double
profile_mpn_smp_basecase (void* arg, unsigned long count);

/*
   As above, for mpn_smp_kara, except that the n parameter is used
   instead of n1, n2.
*/
double
profile_mpn_smp_kara (void* arg, unsigned long count);



double
profile_bfly (void* arg, unsigned long count);

double
profile_mpn_aors (void* arg, unsigned long count);

double
profile_scalar_mul (void* arg, unsigned long count);


void 
prof_main (int argc, char* argv[]);


#ifdef __cplusplus
}
#endif

#endif

// end of file ****************************************************************
