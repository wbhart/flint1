/******************************************************************************

 ZmodF.h

 Copyright (C) 2007, David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
 These are currently used only in the ZmodFpoly module, which supplies the
 Schoenhage-Strassen FFT code.
 
******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "mpn_extras.h"


/*
A ZmodF_t is stored as a *signed* value in two's complement format, using
n+1 limbs. The value is not normalised into any particular range, so the top
limb may pick up overflow bits. Of course the arithmetic functions in this
module may implicitly reduce mod p whenever they like.

More precisely, suppose that the first n limbs are x[0], ..., x[n-1] (unsigned)
and the last limb is x[n] (signed). Then the value being represented is
 x[0] + x[1]*B + ... + x[n-1]*B^(n-1) - x[n]   (mod p).

*/
typedef mp_limb_t* ZmodF_t;


static inline
void ZmodF_swap(ZmodF_t* a, ZmodF_t* b)
{
   ZmodF_t temp = *a;
   *a = *b;
   *b = temp;
}


/*
   Normalises a into the range [0, p).
   (Note that the top limb will be set if and only if a = -1 mod p.)
*/
void ZmodF_normalise(ZmodF_t a, unsigned long n);


/*
   Adjusts a mod p so that the top limb is in the interval [0, 2].

   This in general will be faster then ZmodF_normalise(); in particular
   the branching is much more predictable.
*/
static inline
void ZmodF_fast_reduce(ZmodF_t a, unsigned long n)
{
   mp_limb_t hi = a[n];
   a[n] = 1;
   signed_add_1(a, n+1, 1-hi);
}


/*
   a := 0
*/
static inline
void ZmodF_zero(ZmodF_t a, unsigned long n)
{
   long i = n;
   do a[i] = 0; while (--i >= 0);
}


/*
   b := a
*/
static inline
void ZmodF_set(ZmodF_t b, ZmodF_t a, unsigned long n)
{
   long i = n;
   do b[i] = a[i]; while (--i >= 0);
}


/*
   b := -a
   
   PRECONDITIONS:
      a and b may alias each other
*/
static inline
void ZmodF_neg(ZmodF_t b, ZmodF_t a, unsigned long n)
{
   b[n] = ~a[n] - 1;     // -1 is to make up mod p for 2's complement negation
   long i = n-1;
   do b[i] = ~a[i]; while (--i >= 0);
}


/*
   res := a + b
   
   PRECONDITIONS:
      Any combination of aliasing among res, a, b is allowed.
*/
static inline
void ZmodF_add(ZmodF_t res, ZmodF_t a, ZmodF_t b, unsigned long n)
{
   mpn_add_n(res, a, b, n+1);
}


/*
   res := a - b
   
   PRECONDITIONS:
      Any combination of aliasing among res, a, b is allowed.
*/
static inline
void ZmodF_sub(ZmodF_t res, ZmodF_t a, ZmodF_t b, unsigned long n)
{
   mpn_sub_n(res, a, b, n+1);
}


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


/*
   b := 2^(-s) a

   PRECONDITIONS:
      0 < s < FLINT_BITS_PER_LIMB
      b may alias a
*/
static inline
void ZmodF_short_div_2exp(ZmodF_t b, ZmodF_t a,
                          unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0 && s < FLINT_BITS_PER_LIMB);
   
   // quick adjustment mod p to ensure a is non-negative
   ZmodF_fast_reduce(a, n);

   // do the rotation, and push the overflow back to the top limb
   mp_limb_t overflow = mpn_rshift(b, a, n+1, s);
   mpn_sub_1(b+n-1, b+n-1, 2, overflow);
}


/*
   b := B^s a

   PRECONDITIONS:
      0 < s < n
      b must not alias a
*/
static inline
void ZmodF_mul_Bexp(ZmodF_t b, ZmodF_t a, unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(b != a);

   // let a = ex*B^n + hi*B^(n-s) + lo,
   // where  0 <= lo < B^(n-s)  and  0 <= hi < B^s  and  abs(ex) < B/2.
   // Then the output should be  -ex*B^s + lo*B^s - hi  (mod p).

   long i;
   
   // Put B^s - hi - 1 into b
   i = s-1;
   do b[i] = ~a[n-s+i]; while (--i >= 0);

   // Put lo*B^s into b
   i = n-s-1;
   do b[i+s] = a[i]; while (--i >= 0);

   // Put -B^n into b (to compensate mod p for -1 added in first loop)
   b[n] = (mp_limb_t)(-1L);
   
   // Add (-ex-1)*B^s to b
   signed_add_1(b+s, n-s+1, -a[n]-1);
}


/*
   c := a - 2^(-Bs) b

PRECONDITIONS:
   0 < s < n
   b must not alias c
   a may alias b or c
*/
static inline
void ZmodF_div_Bexp_sub(ZmodF_t c, ZmodF_t a, ZmodF_t b,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(b != c);

   // add low limbs of b to high limbs of a
   c[n] = a[n] + mpn_add_n(c+n-s, b, a+n-s, s);
   // subtract high limbs of b from low limbs of a
   mp_limb_t overflow = b[n] + mpn_sub_n(c, a, b+s, n-s);
   // propagate overflow
   signed_add_1(c+n-s, s+1, -overflow);
}


/*
   c := a + 2^(-Bs) b

   PRECONDITIONS:
      0 < s < n
      b must not alias c
      a may alias b or c
*/
static inline
void ZmodF_div_Bexp_add(ZmodF_t c, ZmodF_t a, ZmodF_t b,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(b != c);

   // subtract low limbs of b from high limbs of a
   c[n] = a[n] - mpn_sub_n(c+n-s, a+n-s, b, s);
   // add high limbs of b to low limbs of a
   mp_limb_t overflow = b[n] + mpn_add_n(c, a, b+s, n-s);
   // propagate overflow
   signed_add_1(c+n-s, s+1, overflow);
}


/*
   c := B^s (a - b)

PRECONDITIONS:   
   c must not alias a or b
   0 < s < n
*/
static inline
void ZmodF_sub_mul_Bexp(ZmodF_t c, ZmodF_t a, ZmodF_t b,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(c != a);
   FLINT_ASSERT(c != b);

   // get low limbs of a - b into high limbs of c
   c[n] = -mpn_sub_n(c+s, a, b, n-s);
   // get high limbs of b - a into low limbs of c
   mp_limb_t overflow = b[n] - a[n] - mpn_sub_n(c, b+n-s, a+n-s, s);
   // propagate overflow
   signed_add_1(c+s, n+1-s, overflow);
}


/*
   b := 2^s a

   PRECONDITIONS:
      0 <= s < n*FLINT_BITS_PER_LIMB
      b may not alias a
*/
void ZmodF_mul_2exp(ZmodF_t b, ZmodF_t a, unsigned long s, unsigned long n);


/*
   b := 2^(s/2) a
   z := destroyed

   PRECONDITIONS:
      0 <= s < 2*n*FLINT_BITS_PER_LIMB
      a may alias z (in which case a gets destroyed)
      b may not alias a or z

   NOTE: a, b, z may get permuted
*/
void ZmodF_mul_sqrt2exp(ZmodF_t* b, ZmodF_t* a, ZmodF_t* z,
                        unsigned long s, unsigned long n);


/*
   c := 2^s (a - b)

   PRECONDITIONS:
      c must not alias a or b
      0 <= s < n*FLINT_BITS_PER_LIMB
*/
void ZmodF_sub_mul_2exp(ZmodF_t c, ZmodF_t a, ZmodF_t b,
                        unsigned long s, unsigned long n);

/*
   a := a + b
   b := 2^s (a - b)
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
void ZmodF_forward_butterfly_2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                  unsigned long s, unsigned long n);


/*
   a := a + b
   b := 2^(s/2) (a - b)
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
void ZmodF_forward_butterfly_sqrt2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                      unsigned long s, unsigned long n);


/*
   a := a + 2^(-s/2) b
   b := a - 2^(-s/2) b
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
void ZmodF_inverse_butterfly_sqrt2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                      unsigned long s, unsigned long n);


/*
   a := a + b
   b := a - b
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
void ZmodF_simple_butterfly(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                            unsigned long n);


// end of file ****************************************************************
