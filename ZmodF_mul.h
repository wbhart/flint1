/******************************************************************************

 ZmodF_mul.h

 Copyright (C) 2007, David Harvey
 
 Routines for multiplication of elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS.
 
******************************************************************************/

#ifndef FLINT_ZMODF_MUL_H
#define FLINT_ZMODF_MUL_H

#include <stdlib.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "ZmodFpoly.h"


// whether to use negacyclic/threeway algorithms in the automatic
// algorithm selection
// todo: these flags are just for testing, they should eventually
// be removed
#define ZMODF_MUL_ENABLE_NEGACYCLIC 1
#define ZMODF_MUL_ENABLE_THREEWAY 1


/*
   Three possible algorithms for multiplication mod p.
   ZMOD_MUL_PLAIN: use mpn_mul_n and then reduce mod p
   ZMOD_MUL_THREEWAY: use B^(3n) + 1 = (B^n + 1)(B^2n - B^n + 1)
   ZMOD_MUL_NEGACYCLIC: use negacyclic FFT
*/
#define ZMODF_MUL_ALGO_PLAIN 0
#define ZMODF_MUL_ALGO_THREEWAY 1
#define ZMODF_MUL_ALGO_NEGACYCLIC 2


/*
This struct stores info used to speed up multiplications mod p
for a specific n.
*/
typedef struct
{
   unsigned long n;

   // possible values for algo are the ZMOD_MUL_ALGO_xyz constants above
   int algo;
   
   // scratch buffer of length 2n (for ZMODF_MUL_ALGO_PLAIN)
   // or length 3n+1 (for ZMODF_MUL_ALGO_THREEWAY)
   mp_limb_t* scratch;

   // m = n/3; used only in ZMODF_MUL_ALGO_THREEWAY
   unsigned long m;
   
   // used only for ZMODF_MUL_ALGO_NEGACYCLIC
   ZmodFpoly_t polys[2];

} ZmodF_mul_info_struct;


// ZmodF_mul_info_t allows reference-like semantics for
// ZmodF_mul_info_struct:
typedef ZmodF_mul_info_struct ZmodF_mul_info_t[1];


/*
Initialises ZmodF_mul_info_t for a given n. This function automatically selects
the best underlying multiplication algorithm for the given n.

The squaring flag is 1 if you intend to use this object for squaring. (This
doesn't affect whether the struct can be used for squaring or multiplying,
the only effect is to select the best algorithm using different tuning
parameters.)

WARNING: the multiplication time does NOT increase monotonically with n.
* If n is divisible by 3, the "threeway" algorithm is available, which is
  faster than the "plain" algorithm for n >= 36 (on our opteron test machine).
* For large n, the "negacyclic" algorithm is available, but there are
  conditions on the 2-divisibility of n (not very onerous though).

NOTE: this function uses stack based memory management.
*/
void ZmodF_mul_info_init(ZmodF_mul_info_t info, unsigned long n, int squaring);


// the following functions initialise with a specific algorithm:
void ZmodF_mul_info_init_plain(ZmodF_mul_info_t info, unsigned long n);
void ZmodF_mul_info_init_threeway(ZmodF_mul_info_t info, unsigned long n);
void ZmodF_mul_info_init_negacyclic(ZmodF_mul_info_t info, unsigned long n,
                                    unsigned long depth);

                            
// releases resources
void ZmodF_mul_info_clear(ZmodF_mul_info_t info);

// sets res := a * b using the given ZmodF_mul_info_t object
void ZmodF_mul_info_mul(ZmodF_mul_info_t, ZmodF_t res, ZmodF_t a, ZmodF_t b);
// sets res := a * a using the given ZmodF_mul_info_t object
void ZmodF_mul_info_sqr(ZmodF_mul_info_t, ZmodF_t res, ZmodF_t a);


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

void _ZmodF_mul_negacyclic_split(ZmodFpoly_t poly, ZmodF_t x, unsigned long n);

void _ZmodF_mul_negacyclic_combine(ZmodF_t x, ZmodFpoly_t poly,
                                   unsigned long n);


FLINT_INLINE
void _ZmodF_mul_threeway_reduce1(ZmodF_t res, ZmodF_t a, unsigned long m);

FLINT_INLINE
void _ZmodF_mul_threeway_reduce2(mp_limb_t* res, ZmodF_t a, unsigned long m);

FLINT_INLINE
void _ZmodF_mul_threeway_crt(mp_limb_t* res, ZmodF_t a, mp_limb_t* b,
                             unsigned long m);


#endif

// end of file ****************************************************************
