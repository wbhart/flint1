/******************************************************************************

 Schonhage-Strassen FFT/IFFT.
 
 Uses van der Hoeven's truncated transforms, modified to use a Bailey-style
 factoring algorithm in order to improve cache efficiency.

 Copyright (C) 2006
 David Harvey

todo: need to acknowledge source of sqrt2 idea
todo: I saw Bill added some constant like FLINT_LG_BITS_PER_LIMB...
this could be very handy throughout this file?

 Throughout the documentation for this file, B denotes FLINT_BITS_PER_LIMB.
 
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include "flint.h"


#define ENABLE_THREAD_CODE 1


#if ENABLE_THREAD_CODE
#include <pthread.h>
#define NUM_THREADS 2
#endif


// todo: get these assertion macros happening uniformly throughout FLINT
#if 0
#define FLINT_ASSERT assert
#else
#define FLINT_ASSERT(condition)
#endif


// ****************************************************
// This block is temporarily here until we move it somewhere else
// (probably into Zvec.c?)


/*
Sets y = \sum_{i=0}^{len-1} x[i] * 2^(ki)
Running time should be O(k*len*log(len))
*/
void naive_KS_pack_coeffs(mpz_t y, mpz_t* x, unsigned long len,
                          unsigned long k)
{
   if (len == 1)
      mpz_set(y, x[0]);
   else
   {
      mpz_t temp;
      mpz_init(temp);
      unsigned long half = len/2;
      naive_KS_pack_coeffs(temp, x, half, k);
      naive_KS_pack_coeffs(y, x + half, len - half, k);
      mpz_mul_2exp(y, y, half*k);
      mpz_add(y, y, temp);
      mpz_clear(temp);
   }
}


/*
Inverse operation of naive_KS_pack_coeffs
(note: y is destroyed)
*/
void naive_KS_unpack_coeffs(mpz_t* x, unsigned long len, mpz_t y,
                            unsigned long k)
{
   if (len == 1)
      mpz_set(x[0], y);
   else
   {
      mpz_t temp;
      mpz_init(temp);
      unsigned long half = len/2;
      if (mpz_tstbit(y, k*half - 1))
      {
         mpz_cdiv_q_2exp(temp, y, half*k);
         mpz_cdiv_r_2exp(y, y, half*k);
      }
      else
      {
         mpz_fdiv_q_2exp(temp, y, half*k);
         mpz_fdiv_r_2exp(y, y, half*k);
      }
      naive_KS_unpack_coeffs(x, half, y, k);
      naive_KS_unpack_coeffs(x + half, len - half, temp, k);
      mpz_clear(temp);
   }
}


/*
Counts maximum number of bits in abs(x[i]) for 0 <= i < len
*/
unsigned long naive_KS_get_max_bits(mpz_t* x, unsigned long len)
{
   unsigned long bits = 0, temp, i;
   for (i = 0; i < len; i++)
   {
      temp = mpz_sizeinbase(x[i], 2);
      if (temp > bits)
         bits = temp;
   }
   return bits;
}


/*
A very silly implementation of KS multiplication.
x1 and x2 should be polys of length len1 and len2
y should be exactly len1+len2-1 long, and already mpz_init'd
*/
void naive_KS_mul(mpz_t* y,
                  mpz_t* x1, unsigned long len1,
                  mpz_t* x2, unsigned long len2)
{
   mpz_t z1;
   mpz_t z2;
   mpz_init(z1);
   mpz_init(z2);

   unsigned long output_len = len1 + len2 - 1;
   unsigned long bits1 = naive_KS_get_max_bits(x1, len1);
   unsigned long bits2 = naive_KS_get_max_bits(x2, len2);
   unsigned long bits = bits1 + bits2 + 2 + ceil_log2(output_len);

   naive_KS_pack_coeffs(z1, x1, len1, bits);
   naive_KS_pack_coeffs(z2, x2, len2, bits);
   mpz_mul(z1, z1, z2);
   naive_KS_unpack_coeffs(y, output_len, z1, bits);
   
   mpz_clear(z1);
   mpz_clear(z2);
}


// ************************ end temporary block




// tests whether a rotation by 2^(zzz/2) is a rotation by a whole
// number of limbs
#define IS_WHOLE_LIMB_ROTATION(zzz) (!((zzz) & (2*FLINT_BITS_PER_LIMB - 1)))

// tests whether a rotation by 2^(zzz/2) is a rotation by a whole
// number of bits
#define IS_WHOLE_BIT_ROTATION(zzz) (!((zzz) & 1))
  

/*
Add the given *signed* limb to the buffer [x, x+count), much like
mpn_add_1 and mpn_sub_1 (except it's always inplace).

PRECONDITIONS:
   count >= 1
*/
inline void ssfft_signed_add_1(mp_limb_t* x, unsigned long count,
                               mp_limb_signed_t limb)
{
   FLINT_ASSERT(count >= 1);
   
   // If the high bit of *x doesn't change when we add "limb" to it,
   // then there's no possibility of overflow. If we assume that abs(limb)
   // is relatively small (which will always be the case here, because we only
   // ever use this function to add an overflow limb to something), then
   // this will be a very predictable branch.
   mp_limb_t temp = *x + limb;
   if ((mp_limb_signed_t)(temp ^ *x) >= 0)
      // the likely case
      *x = temp;
   else
   {
      // the unlikely case; here we need to branch based on the sign of
      // the limb being added
      if (limb >= 0)
         mpn_add_1(x, x, count, limb);
      else
         mpn_sub_1(x, x, count, -limb);
   }
}


/* ============================================================================

   Basic coefficient operations.

These are the building blocks of the butterflies and hard-wired short
transforms. They operate on coefficient buffers n+1 limbs long. Each
coefficient is interpreted as an integer in 2's complement form (i.e. the high
bit of x[n] is the sign bit). It represents an integer mod p = 2^(Bn) + 1, but
need not be normalised; we generally let bits build up happily in the high
limb as long as we can get away with it. All the functions below describe
exactly what is the worst number of bits of the high limb of the output that
could end up occupied.

Note: some of the cross relations appear to use up TWO bits of the overflow
limb instead of just one. I haven't analysed this carefully yet, but for now,
just to be safe, we do a fast approximation reduction after any cross relation.

=============================================================================*/


/*
Approximate reduction of x mod p.

POSTCONDITIONS:
   0 <= x[n] <= 2.
*/
inline void basic_fast_reduce(mp_limb_t* x, unsigned long n)
{
   mp_limb_t high = x[n];
   x[n] = 1;
   ssfft_signed_add_1(x, n+1, 1-high);
}


/*
   b := 2^(-s) a   mod p

PRECONDITIONS:
   0 < s < B
   b may alias a

POSTCONDITIONS:
   -1 <= b[n] <= 1
    0 <= a[n] <= 2
    
   (Note: a may get modified, but it won't change mod p.)
*/
inline void basic_unrotate_bits(mp_limb_t* b, mp_limb_t* a,
                                unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < FLINT_BITS_PER_LIMB);
   
   // quick adjustment mod p to ensure a is positive
   basic_fast_reduce(a, n);

   // do the rotation, and push the overflow back to the top limb
   mp_limb_t overflow = mpn_rshift(b, a, n+1, s);
   mpn_sub_1(b+n-1, b+n-1, 2, overflow);
}


/*
   b := 2^s a   mod p

PRECONDITIONS:
   0 < s < B
   b may alias a

POSTCONDITIONS:
   Number of bits of b[n] is s more than number of bits in a[n].
   (i.e. no normalisation is performed at all.)
*/
inline void basic_rotate_bits(mp_limb_t* b, mp_limb_t* a,
                              unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < FLINT_BITS_PER_LIMB);
   
   mpn_lshift(b, a, n+1, s);
}


/*
   b := 2^(Bs) a   mod p

PRECONDITIONS:
   b must not alias a
   0 < s < n

POSTCONDITIONS:
   -2 <= b[n] <= 0.
*/
inline void basic_rotate_limbs(mp_limb_t* b, mp_limb_t* a,
                               unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(b != a);

   // let a = ex*2^(Bn) + hi*2^((n-s)B) + lo,
   // where 0 <= lo < 2^((n-s)B) and 0 <= hi < 2^(Bs) and abs(ex) < 2^(B-1).
   // Then the output should be -ex*2^(Bs) + lo*2^(Bs) - hi (mod p).

   unsigned long i;
   
   // Put 2^(Bs) - hi - 1 into b
   i = s;
   do b[i-1] = ~a[n-s+i-1]; while (--i);

   // Put lo*2^(Bs) into b
   i = n-s;
   do b[i+s-1] = a[i-1]; while (--i);

   // Put -2^(Bn) into b (to compensate mod p for -1 added in first loop)
   b[n] = (mp_limb_t)(-1L);
   
   // Add (-ex-1)*2^(Bs) to b
   ssfft_signed_add_1(b+s, n-s+1, -a[n]-1);
}


/*
   b := 2^(-Bs) a   mod p

PRECONDITIONS:
   b must not alias a
   0 < s < n

POSTCONDITIONS:
   -2 <= b[n] <= 0
*/
inline void basic_unrotate_limbs(mp_limb_t* b, mp_limb_t* a,
                                 unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(b != a);

   // let a = ex*2^(Bn) + hi*2^(sB) + lo,
   // where 0 <= lo < 2^(sB) and 0 <= hi < 2^(B(n-s)) and abs(ex) < 2^(B-1).
   // Then the output should be ex*2^(B(n-s)) + hi - lo*2^(B(n-s)).
   
   unsigned long i;
   
   // Put hi into b
   i = n-s;
   do b[i-1] = a[s+i-1]; while (--i);
      
   // Put 2^(B(n-s)) (2^(Bs) - lo - 1) into b
   i = s;
   do b[n-s+i-1] = ~a[i-1]; while (--i);
      
   // Put -2^(Bn) into b (to compensate for extra 2^(Bn) in previous loop)
   b[n] = (mp_limb_t)(-1L);

   // Add (ex+1)*2^(B(n-s)) to b
   ssfft_signed_add_1(b+n-s, s+1, a[n]+1);
}


/*
   c := 2^(Bs) (a - b)    mod p

PRECONDITIONS:   
   c must not alias a or b
   0 < s < n

POSTCONDITIONS:
   -2 <= c[n] <= 1
*/
inline void basic_sub_rotate_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
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
   ssfft_signed_add_1(c+s, n+1-s, overflow);
}


/*
   c := a - 2^(Bs) b   mod p

PRECONDITIONS:
   0 < s < n
   b must not alias c
   a may alias b or c

POSTCONDITIONS:
   -2 <= (c[n] - a[n]) <= 1
   i.e. c[n] is very close to a[n]
*/
inline void basic_rotate_sub_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                                   unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(b != c);

   // subtract low limbs of b from high limbs of a
   c[n] = a[n] - mpn_sub_n(c+s, a+s, b, n-s);
   // add high limbs of b to low limbs of a
   mp_limb_t overflow = b[n] + mpn_add_n(c, a, b+n-s, s);
   // propagate overflow
   ssfft_signed_add_1(c+s, n+1-s, overflow);
}


/*
   c := a - 2^(-Bs) b   mod p

PRECONDITIONS:
   0 < s < n
   b must not alias c
   a may alias b or c

POSTCONDITIONS:
   -1 <= (c[n] - a[n]) <= 2
   i.e. c[n] is very close to a[n]
*/
inline void basic_unrotate_sub_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
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
   ssfft_signed_add_1(c+n-s, s+1, -overflow);
}


/*
   c := a + 2^(-Bs) b   mod p

PRECONDITIONS:
   0 < s < n
   b must not alias c
   a may alias b or c

POSTCONDITIONS:
   -2 <= (c[n] - a[n]) <= 1
   i.e. c[n] is very close to a[n]
*/
inline void basic_unrotate_add_limbs(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
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
   ssfft_signed_add_1(c+n-s, s+1, overflow);
}


/*
   c := a + b

PRECONDITIONS:   
   All combinations of aliasing are permitted.

POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k, then abs(c) < 2^(k+1).
*/
inline void basic_add(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                      unsigned long n)
{
   mpn_add_n(c, a, b, n+1);
}


/*
   c := a - b

PRECONDITIONS:   
   All combinations of aliasing are permitted.

POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k, then abs(c) < 2^(k+1).
*/
inline void basic_sub(mp_limb_t* c, mp_limb_t* a, mp_limb_t* b,
                      unsigned long n)
{
   mpn_sub_n(c, a, b, n+1);
}


/*
   b := a
*/
inline void basic_copy(mp_limb_t* b, mp_limb_t* a, unsigned long n)
{
   unsigned long i = n+1;
   do b[i-1] = a[i-1]; while (--i);
}


/*
   b := -a  mod p

PRECONDITIONS:
   a may alias b

POSTCONDITIONS:
   b[n] = -a[n] - 2
   i.e. b[n] is very close to a[n]
*/
inline void basic_negate(mp_limb_t* b, mp_limb_t* a, unsigned long n)
{
   b[n] = ~a[n] - 1;     // -1 is to make up mod p for 2's complement negation
   unsigned long i = n;
   do b[i-1] = ~a[i-1]; while (--i);
}



/* ============================================================================

   Compound coefficient operations.

These are some more complicated coefficient operations which are built up
of the basic operations of the previous section.

They operate on *pointers* to coefficient buffers. Most of them take a scratch
buffer. They may permute the buffers passed to them.

=============================================================================*/


/* 
Swaps arrays of limbs efficiently
*/
inline void swap_limb_ptrs(mp_limb_t** x, mp_limb_t** y)
{
   mp_limb_t* temp;
   temp = *x;
   *x = *y;
   *y = temp;
}


/*
   a := a + b
   b := a - b
   scratch := destroyed
   
PRECONDITIONS:
   a, b, scratch may not alias each other
   
POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k, then both outputs have abs < 2^(k+1).
*/
void coeff_forward_simple_butterfly(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch, unsigned long n)
{
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   basic_sub(*scratch, *a, *b, n);
   basic_add(*a, *a, *b, n);
   swap_limb_ptrs(scratch, b);
}


/*
Inverse simple butterfly is same as forward simple butterfly
*/
#define coeff_inverse_simple_butterfly \
        coeff_forward_simple_butterfly

/*
   a := a + b
   b := 2^(Bs) (a - b)
   scratch := destroyed
   
PRECONDITIONS:
   0 < s < n
   a, b, scratch may not alias each other

POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k, then abs(a_out) < 2^(k+1).
   -2 <= b_out[n] <= 1
*/
void coeff_forward_butterfly_limbs(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   basic_sub_rotate_limbs(*scratch, *a, *b, s, n);
   basic_add(*a, *a, *b, n);
   swap_limb_ptrs(scratch, b);
}


/*
   a := a + 2^(-Bs) b
   b := a - 2^(-Bs) b
   scratch := destroyed
   
PRECONDITIONS:
   0 < s < n
   a, b, scratch may not alias each other
   
POSTCONDITIONS:
   -2 <= a_out[n] - a_in[n] <= 1
   -1 <= b_out[n] - a_in[n] <= 2
   i.e. both outputs are very close to a_input
*/
void coeff_inverse_butterfly_limbs(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s > 0);
   FLINT_ASSERT(s < n);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   basic_unrotate_sub_limbs(*scratch, *a, *b, s, n);
   basic_unrotate_add_limbs(*a, *a, *b, s, n);
   swap_limb_ptrs(scratch, b);
}


/*
   c := 2^s (a - b)    mod p

PRECONDITIONS:   
   c must not alias a or b
   0 <= s < nB

POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k and k >= Bn, then abs(c) < 2^(k+1).
*/
void coeff_sub_rotate_bits(mp_limb_t** c, mp_limb_t** a, mp_limb_t** b,
                           unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*c != *a);
   FLINT_ASSERT(*c != *b);

   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
   if (bits)
   {
      // we will shift a-b left by s+1 limbs...
      if (++s == n)
         basic_sub(*c, *b, *a, n);
      else
         basic_sub_rotate_limbs(*c, *a, *b, s, n);
      
      // ... and then shift right by remaining bits
      basic_unrotate_bits(*c, *c, FLINT_BITS_PER_LIMB - bits, n);
   }
   else
   {
      if (s)
         basic_sub_rotate_limbs(*c, *a, *b, s, n);
      else
         basic_sub(*c, *a, *b, n);
   }
}


/*
   a := a + b
   b := 2^s (a - b)
   scratch := destroyed
   
PRECONDITIONS:
   0 <= s < Bn
   a, b, scratch may not alias each other

POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k and k >= Bn,
   then abs(a_out) < 2^(k+1) and abs(b_out) < 2^(k+1).
*/
void coeff_forward_butterfly_bits(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   coeff_sub_rotate_bits(scratch, a, b, s, n);
   basic_add(*a, *a, *b, n);
   swap_limb_ptrs(scratch, b);
}


/*
   a := a + 2^(-s) b
   b := a - 2^(-s) b
   scratch := destroyed
   
PRECONDITIONS:
   0 <= s < nB
   a, b, scratch may not alias each other
   
POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k and k >= Bn,
   then abs(a_out) < 2^(k+1) and abs(b_out) < 2^(k+1).
*/
void coeff_inverse_butterfly_bits(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   if (bits)
      // shift right by leftover bits
      basic_unrotate_bits(*b, *b, bits, n);
   
   s /= FLINT_BITS_PER_LIMB;
   if (s)
   {
      basic_unrotate_sub_limbs(*scratch, *a, *b, s, n);
      basic_unrotate_add_limbs(*a, *a, *b, s, n);
   }
   else
   {
      basic_sub(*scratch, *a, *b, n);
      basic_add(*a, *a, *b, n);
   }
   swap_limb_ptrs(scratch, b);
}


/*
   b := 2^s a
   
PRECONDITIONS:
   0 <= s < Bn
   a, b may not alias each other
   
POSTCONDITIONS:
   If abs(a) < 2^k and k >= Bn+1, then abs(b) < 2^k.

todo: we can do better than this. Should bit-rotate the whole array into the
target, and then negate only part of the array (with appropriate bit-fiddling
at the end). This will reduce the amount of copying by 50% on average.
*/
void coeff_rotate_bits(mp_limb_t** b, mp_limb_t** a,
                       unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);

   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
   if (bits)
   {
      // we will shift left by s+1 limbs...
      if (++s == n)
         basic_negate(*b, *a, n);
      else
         basic_rotate_limbs(*b, *a, s, n);

      // ... and then shift right by remaining bits
      basic_unrotate_bits(*b, *b, FLINT_BITS_PER_LIMB - bits, n);
   }
   else
   {
      if (s)
         basic_rotate_limbs(*b, *a, s, n);
      else
         basic_copy(*b, *a, n);
   }
}


/*
   a := 2*a - b
   b := 2^s (a - b)
   scratch := destroyed
   
PRECONDITIONS:
   0 <= s < nB
   a, b, scratch may not alias each other
   
POSTCONDITIONS:
   If abs(a) < 2^k and abs(b) < 2^k and k >= Bn, then abs(b_out) < 2^(k+1)
   0 <= a_out[n] <= 2

NOTE:
   This is one of the basic "cross relations" in the truncated transform.
   It maps (x+y, 2y) to (2x, 2^s (x-y)).
*/
void coeff_cross_butterfly_bits(
         mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
         unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*a != *b);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*a != *scratch);

   coeff_sub_rotate_bits(scratch, a, b, s, n);
   basic_add(*a, *a, *a, n);
   basic_sub(*a, *a, *b, n);
   basic_fast_reduce(*a, n);    // just to be safe
   swap_limb_ptrs(b, scratch);
}


/*
   a := 2*a - b

POSTCONDITIONS:
   0 <= a_out[n] <= 2
   
NOTE:
   This is not quite a "cross relation", but it's used frequently enough to
   deserve its own function.
*/
void coeff_cross_butterfly2(mp_limb_t** a, mp_limb_t** b, unsigned long n)
{
   basic_add(*a, *a, *a, n);
   basic_sub(*a, *a, *b, n);
   basic_fast_reduce(*a, n);    // just to be safe
}


/*
   b := (1 - 2^(nB/2)) * a
   a := destroyed

PRECONDITIONS:
   a must not alias b
   
POSTCONDITIONS:
   If abs(a) < 2^k and k >= nB+1, then abs(b) < 2^(k+1).
   
NOTE:
   This is a helper function for rotations involving sqrt2.
   See coeff_rotate_arbitrary().
*/
void coeff_sqrt2_helper(mp_limb_t** b, mp_limb_t** a, unsigned long n)
{
   FLINT_ASSERT(*a != *b);

   if (n & 1)
   {
      basic_unrotate_bits(*b, *a, FLINT_BITS_PER_LIMB/2, n);
      if (n == 1)
         // special case since basic_rotate_sub_limbs requires s < n
         basic_add(*b, *b, *a, n);
      else
      {
         basic_rotate_sub_limbs(*a, *a, *b, (n >> 1) + 1, n);
         swap_limb_ptrs(a, b);
      }
   }
   else
      basic_rotate_sub_limbs(*b, *a, *a, n/2, n);
}


/*
For odd s, finds j and k such that 2^(s/2) is decomposed into
  2^(-j) * 2^(k*B) * (1 - 2^(nB/2)),
where 0 <= j < B, and 0 <= k < 2n.

i.e. we are decomposition a rotation involving a sqrt2 into a full limb shift,
a fractional limb shift, and a call to coeff_sqrt2_helper.

PRECONDIITONS:
   s must be odd
   0 <= s < 2nB
   
todo: there are some cases here (when n is odd) where we might sometimes be
able to avoid one of the bit rotations... although it's pretty rare, maybe
not worth worrying about
*/
void decompose_rotation(unsigned long* k, unsigned long* j,
                        unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s & 1);
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);

   // first split into 2^r * (1 - 2^(nB/2))
   unsigned long r = (s >> 1) - 3*n*FLINT_BITS_PER_LIMB/4;
   if ((long)r < 0)
      r += 2*n*FLINT_BITS_PER_LIMB;

   unsigned long bits = r & (FLINT_BITS_PER_LIMB - 1);
   r /= FLINT_BITS_PER_LIMB;
   if (bits)
   {
      *j = FLINT_BITS_PER_LIMB - bits;
      if (++r == 2*n)
         r = 0;
   }
   else
      *j = 0;

   *k = r;
}


/*
   b := 2^(s/2) a
   scratch := destroyed
   
PRECONDITIONS:
   b must not alias a or scratch
   a MAY alias scratch; in that case obviously a is destroyed
   0 <= s < 2nB
   
POSTCONDITIONS:
   If abs(a) < 2^k and k >= Bn+1, then abs(b) < 2^(k+1).
*/
void coeff_rotate_arbitrary(mp_limb_t** b, mp_limb_t** a, mp_limb_t** scratch,
                            unsigned long s, unsigned long n)
{
   FLINT_ASSERT(s < 2*n*FLINT_BITS_PER_LIMB);
   FLINT_ASSERT(*b != *scratch);
   FLINT_ASSERT(*b != *a);

   if (s & 1)
   {
      // involves sqrt2
      unsigned long limbs;
      unsigned long bits;
      decompose_rotation(&limbs, &bits, s, n);
      
      // try to do full limb rotation first
      if (limbs)
      {
         if (limbs < n)
            basic_rotate_limbs(*b, *a, limbs, n);
         else if (limbs > n)
            basic_unrotate_limbs(*b, *a, 2*n-limbs, n);
         else   // limbs == n
         {
            FLINT_ASSERT(limbs == n);
            basic_negate(*b, *a, n);
         }
         
         // do sqrt2 portion
         coeff_sqrt2_helper(scratch, b, n);
         
         // fractional limb component
         if (bits)
            basic_unrotate_bits(*scratch, *scratch, bits, n);
         
         swap_limb_ptrs(b, scratch);
      }
      else
      {
         // no full limb rotation available, do fraction limb rotation
         // directly into scratch
         if (bits)
            basic_unrotate_bits(*scratch, *a, bits, n);
         else
         {
            // damn, no full limb or fractional limb rotations available,
            // forced to make a copy
            basic_copy(*scratch, *a, n);
         }

         // do sqrt2 portion
         coeff_sqrt2_helper(b, scratch, n);
      }
   }
   else
   {
      // no sqrt2 involved
      coeff_rotate_bits(b, a, s >> 1, n);
   }
}



/* ============================================================================

   Forward transform.

It's a fast fourier transform with Schoenhage-Strassen-style coefficients,
using a square root of 2 mod p to extend the transform range, "twisted" roots
of unity to enable the recursive Bailey-style FFT factoring, and some fairly
involved truncation logic to keep performance smooth.

The transform takes several parameters, described below. There are several
versions of the transform, optimised for particular choices of parameters.

Let M = 2^m be the transform length. We are working modulo the Fermat "prime"
p = 2^(Bn) + 1.

Let U = some power of two. The transform requires a primitive (UM)-th
root of unity mod p, denoted by w = 2^(r/2). We allow r to be odd, in which
case we take 2^(1/2) to be
   2^(3nB/4) - 2^(nB/4),
which is a square root of 2 mod p.

Denote the input coefficients by a_0, ..., a_{M-1} and the output coefficients
by b_0, ..., b_{M-1}. Then the transform being computed is:

    b_i = \sum_{k=0}^{M-1} w^{i' (u + kU)} a_k,
    
where 0 <= u < U, 0 <= i < M, and i' is the (length m) bit-reversal of i.

(If U = 1 and u = 0, this is the usual "untwisted" fourier transform of
length M with bit-reversed output.)

The coefficients being operated on are stored at x[0], x[t], x[2*t], ...,
x[(M-1)*t]. The buffers are each n+1 limbs long. The transform may permute
these buffers, and and the permutation may include a scratch buffer (also n+1
limbs) which is passed in separately.

The transforms are truncated in the following sense: given 1 <= z <= M and
1 <= g <= M, they only compute the first g coefficients of the output, and
they assume that all inputs from z onwards are zero.

More precisely, when calling the transform, the input a_i should be stored at
x + i*t for 0 <= i < z (the remaining input coefficients are ignored, and a_i
is assumed to be zero for z <= i < M), and when the transform returns, the
output b_i is stored at x + i*t for 0 <= i < g (the remaining output
coefficients will be garbage). Note that in general ALL of the M coefficient
buffers may be used in intermediate steps, regardless of the values of z and g.

The transforms take parameters ru and rU, which are r*u and r*U in the above
notation. Some of them take instead ru_bits, which means ru/2, or ru_limbs,
which means ru/(2*FLINT_BITS_PER_LIMB). Similarly for rU.

=============================================================================*/


// ----------------------------------------------------------------------------
// length 4 transforms

/*
Assumes M = z = g = 4, all rotations by full limbs
*/
void ssfft_fft_size4_z4g4_limbs(
      mp_limb_t** x, unsigned long t,
      unsigned long ru_limbs, unsigned long rU_limbs,
      unsigned long n, mp_limb_t** scratch)
{
   FLINT_ASSERT(ru_limbs < rU_limbs);

   if (ru_limbs)
   {
      coeff_forward_butterfly_limbs(x, x+2*t, scratch, ru_limbs, n);
      coeff_forward_butterfly_limbs(x+t, x+3*t, scratch,
                                    ru_limbs + rU_limbs, n);
      coeff_forward_butterfly_limbs(x, x+t, scratch, 2*ru_limbs, n);
      coeff_forward_butterfly_limbs(x+2*t, x+3*t, scratch, 2*ru_limbs, n);
   }
   else
   {
      coeff_forward_simple_butterfly(x, x+2*t, scratch, n);
      coeff_forward_butterfly_limbs(x+t, x+3*t, scratch, rU_limbs, n);
      coeff_forward_simple_butterfly(x, x+t, scratch, n);
      coeff_forward_simple_butterfly(x+2*t, x+3*t, scratch, n);
   }
}


/*
The following versions all assume M = 4, and z and g are given in the name
of the function. Rotations are by an integral number of bits.
*/

void ssfft_fft_size4_z1g1_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
}

void ssfft_fft_size4_z2g1_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   basic_add(x[0], x[0], x[t], n);
}

void ssfft_fft_size4_z3g1_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   basic_add(x[0], x[0], x[t], n);
   basic_add(x[0], x[0], x[2*t], n);
}

void ssfft_fft_size4_z4g1_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   basic_add(x[0], x[0], x[t], n);
   basic_add(x[0], x[0], x[2*t], n);
   basic_add(x[0], x[0], x[3*t], n);
}

void ssfft_fft_size4_z1g2_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_rotate_bits(x+t, x, 2*ru_bits, n);
}

void ssfft_fft_size4_z2g2_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_forward_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
}

void ssfft_fft_size4_z3g2_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   basic_add(x[0], x[0], x[2*t], n);
   coeff_forward_butterfly_bits(x, x+t, x+2*t, 2*ru_bits, n);
}

void ssfft_fft_size4_z4g2_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   basic_add(x[0], x[0], x[2*t], n);
   basic_add(x[t], x[t], x[3*t], n);
   coeff_forward_butterfly_bits(x, x+t, x+3*t, 2*ru_bits, n);
}

void ssfft_fft_size4_z1g3_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_rotate_bits(x+t, x, 2*ru_bits, n);
   coeff_rotate_bits(x+2*t, x, ru_bits, n);
}

void ssfft_fft_size4_z2g3_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_rotate_bits(x+2*t, x, ru_bits, n);
   coeff_rotate_bits(x+3*t, x+t, ru_bits + rU_bits, n);
   basic_add(x[2*t], x[2*t], x[3*t], n);
   coeff_forward_butterfly_bits(x, x+t, x+3*t, 2*ru_bits, n);
}

void ssfft_fft_size4_z3g3_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_forward_butterfly_bits(x, x+2*t, x+3*t, ru_bits, n);
   coeff_rotate_bits(x+3*t, x+t, ru_bits + rU_bits, n);
   basic_add(x[2*t], x[2*t], x[3*t], n);
   coeff_forward_butterfly_bits(x, x+t, x+3*t, 2*ru_bits, n);
}

void ssfft_fft_size4_z4g3_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_forward_butterfly_bits(x, x+2*t, scratch, ru_bits, n);
   coeff_forward_butterfly_bits(x+t, x+3*t, scratch, ru_bits + rU_bits, n);
   coeff_forward_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   basic_add(x[2*t], x[2*t], x[3*t], n);
}

void ssfft_fft_size4_z1g4_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_rotate_bits(x+t, x, 2*ru_bits, n);
   coeff_rotate_bits(x+2*t, x, ru_bits, n);
   coeff_rotate_bits(x+3*t, x+2*t, 2*ru_bits, n);
}

void ssfft_fft_size4_z2g4_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_rotate_bits(x+2*t, x, ru_bits, n);
   coeff_rotate_bits(x+3*t, x+t, ru_bits + rU_bits, n);
   coeff_forward_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   coeff_forward_butterfly_bits(x+2*t, x+3*t, scratch, 2*ru_bits, n);
}

void ssfft_fft_size4_z3g4_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_forward_butterfly_bits(x, x+2*t, scratch, ru_bits, n);
   coeff_rotate_bits(x+3*t, x+t, ru_bits + rU_bits, n);
   coeff_forward_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   coeff_forward_butterfly_bits(x+2*t, x+3*t, scratch, 2*ru_bits, n);
}

void ssfft_fft_size4_z4g4_bits(mp_limb_t** x, unsigned long t,
                          unsigned long ru_bits, unsigned long rU_bits,
                          unsigned long n, mp_limb_t** scratch)
{
   coeff_forward_butterfly_bits(x, x+2*t, scratch, ru_bits, n);
   coeff_forward_butterfly_bits(x+t, x+3*t, scratch, ru_bits + rU_bits, n);
   coeff_forward_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   coeff_forward_butterfly_bits(x+2*t, x+3*t, scratch, 2*ru_bits, n);
}

// ssfft_fft_size4_bits_func_t is the type of the 16 functions above
typedef void (*ssfft_fft_size4_bits_func_t)
                        (mp_limb_t**, unsigned long, unsigned long,
                         unsigned long, unsigned long, mp_limb_t**);


// A lookup table for the 16 functions above, indexed by z-1 and g-1
ssfft_fft_size4_bits_func_t
ssfft_fft_size4_bits_table[4][4] = {
   {ssfft_fft_size4_z1g1_bits, ssfft_fft_size4_z1g2_bits,
    ssfft_fft_size4_z1g3_bits, ssfft_fft_size4_z1g4_bits},
   {ssfft_fft_size4_z2g1_bits, ssfft_fft_size4_z2g2_bits,
    ssfft_fft_size4_z2g3_bits, ssfft_fft_size4_z2g4_bits},
   {ssfft_fft_size4_z3g1_bits, ssfft_fft_size4_z3g2_bits,
    ssfft_fft_size4_z3g3_bits, ssfft_fft_size4_z3g4_bits},
   {ssfft_fft_size4_z4g1_bits, ssfft_fft_size4_z4g2_bits,
    ssfft_fft_size4_z4g3_bits, ssfft_fft_size4_z4g4_bits}};


/*
This version assumes M = 4, and rotation by an integral number of bits.
Any 1 <= z <= 4 and 1 <= g <= 4 is allowed.
*/
inline void ssfft_fft_size4_bits(
       mp_limb_t** x, unsigned long t, unsigned long z, unsigned long g,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch)
{
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= 4);
   FLINT_ASSERT(g >= 1);
   FLINT_ASSERT(g <= 4);
   FLINT_ASSERT(ru_bits < rU_bits);

   ssfft_fft_size4_bits_table[z-1][g-1](x, t, ru_bits, rU_bits, n, scratch);
}


// ----------------------------------------------------------------------------
// length 8 transforms

/*
Assumes M = z = g = 8, all rotations by full limbs
*/
void ssfft_fft_size8_z8g8_limbs(
      mp_limb_t** x, unsigned long t,
      unsigned long ru_limbs, unsigned long rU_limbs,
      unsigned long n, mp_limb_t** scratch)
{
   FLINT_ASSERT(ru_limbs < rU_limbs);

   // First row of butterflies.
   if (ru_limbs)
      coeff_forward_butterfly_limbs(x, x+4*t, scratch, ru_limbs, n);
   else
      coeff_forward_simple_butterfly(x, x+4*t, scratch, n);
   coeff_forward_butterfly_limbs(x+t, x+5*t, scratch, ru_limbs + rU_limbs, n);
   coeff_forward_butterfly_limbs(x+2*t, x+6*t, scratch,
                                 ru_limbs + 2*rU_limbs, n);
   coeff_forward_butterfly_limbs(x+3*t, x+7*t, scratch,
                                 ru_limbs + 3*rU_limbs, n);

   // recurse into two halves
   ssfft_fft_size4_z4g4_limbs(x, t, 2*ru_limbs, 2*rU_limbs, n, scratch);
   ssfft_fft_size4_z4g4_limbs(x+4*t, t, 2*ru_limbs, 2*rU_limbs, n, scratch);
}


/*
This version assumes M = 8, and rotation by an integral number of bits.
Any 1 <= z <= 8 and 1 <= g <= 8 is allowed.
*/
void ssfft_fft_size8_bits(
       mp_limb_t** x, unsigned long t, unsigned long z, unsigned long g,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch)
{
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= 8);
   FLINT_ASSERT(g >= 1);
   FLINT_ASSERT(g <= 8);
   FLINT_ASSERT(ru_bits < rU_bits);

   if (g <= 4)
   {
      // Only left four coefficients required.
      
      switch (z)
      {
         // If more than half of the coefficients are nonzero, need to
         // butterfly them over to the left half....
      
         case 8:
            basic_add(x[3*t], x[3*t], x[7*t], n);
         case 7:
            basic_add(x[2*t], x[2*t], x[6*t], n);
         case 6:
            basic_add(x[t], x[t], x[5*t], n);
         case 5:
            basic_add(x[0], x[0], x[4*t], n);
            // ... and then run length 4 on the left half.
            ssfft_fft_size4_bits(x, t, 4, g, 2*ru_bits, 2*rU_bits, n, scratch);
            break;
            
         // Otherwise, just run left half immediately:
         default:
            ssfft_fft_size4_bits(x, t, z, g, 2*ru_bits, 2*rU_bits, n, scratch);
            break;
      }
   }
   else
   {
      // g >= 5, so need the whole left half, and at least one coefficient
      // in the right half.

      switch (z)
      {
         // If the right half of the coefficients are zero, butterfly the
         // left coefficients over to the right,
      
         case 4:
            coeff_rotate_bits(x+7*t, x+3*t, ru_bits + 3*rU_bits, n);
         case 3:
            coeff_rotate_bits(x+6*t, x+2*t, ru_bits + 2*rU_bits, n);
         case 2:
            coeff_rotate_bits(x+5*t, x+t, ru_bits + rU_bits, n);
         case 1:
            coeff_rotate_bits(x+4*t, x, ru_bits, n);
            // and run length 4 transforms on both halves:
            ssfft_fft_size4_bits(x+4*t, t, z, g-4,
                                 2*ru_bits, 2*rU_bits, n, scratch);
            ssfft_fft_size4_bits(x, t, z, 4, 2*ru_bits, 2*rU_bits, n, scratch);
            return;

         // If there are nonzero coefficients on the right, we need to
         // run a bunch of butterflies/rotations to prepare for length 4
         // transforms on both halves.
            
         case 5:
            coeff_forward_butterfly_bits(x, x+4*t, scratch, ru_bits, n);
            coeff_rotate_bits(x+5*t, x+t, ru_bits + rU_bits, n);
            coeff_rotate_bits(x+6*t, x+2*t, ru_bits + 2*rU_bits, n);
            coeff_rotate_bits(x+7*t, x+3*t, ru_bits + 3*rU_bits, n);
            break;

         case 6:
            coeff_forward_butterfly_bits(x, x+4*t, scratch, ru_bits, n);
            coeff_forward_butterfly_bits(x+t, x+5*t, scratch,
                                         ru_bits + rU_bits, n);
            coeff_rotate_bits(x+6*t, x+2*t, ru_bits + 2*rU_bits, n);
            coeff_rotate_bits(x+7*t, x+3*t, ru_bits + 3*rU_bits, n);
            break;

         case 7:
            coeff_forward_butterfly_bits(x, x+4*t, scratch, ru_bits, n);
            coeff_forward_butterfly_bits(x+t, x+5*t, scratch,
                                         ru_bits + rU_bits, n);
            coeff_forward_butterfly_bits(x+2*t, x+6*t, scratch,
                                         ru_bits + 2*rU_bits, n);
            coeff_rotate_bits(x+7*t, x+3*t, ru_bits + 3*rU_bits, n);
            break;

         default:  // case 8:
            coeff_forward_butterfly_bits(x, x+4*t, scratch, ru_bits, n);
            coeff_forward_butterfly_bits(x+t, x+5*t, scratch,
                                         ru_bits + rU_bits, n);
            coeff_forward_butterfly_bits(x+2*t, x+6*t, scratch,
                                         ru_bits + 2*rU_bits, n);
            coeff_forward_butterfly_bits(x+3*t, x+7*t, scratch,
                                         ru_bits + 3*rU_bits, n);
            break;
      }

      // two halves:
      ssfft_fft_size4_bits(x+4*t, t, 4, g-4, 2*ru_bits, 2*rU_bits, n, scratch);
      ssfft_fft_size4_bits(x, t, 4, 4, 2*ru_bits, 2*rU_bits, n, scratch);
   }
}



// ----------------------------------------------------------------------------
// longer transforms


/*
*/
void ssfft_fft_size2(mp_limb_t** x, unsigned long t,
                     unsigned long z, unsigned long g,
                     unsigned long ru, unsigned long n,
                     mp_limb_t** scratch)
{
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= 2);
   FLINT_ASSERT(g >= 1);
   FLINT_ASSERT(g <= 2);
   
   if (g == 1)
   {
      if (z == 2)
         basic_add(x[0], x[0], x[t], n);
      return;
   }
   
   // g == 2 case
   
   if (z == 1)
      coeff_rotate_arbitrary(x+t, x, scratch, ru, n);
   else
   {
      if (ru & 1)
      {
         // sqrt2 involved
         basic_sub(*scratch, x[0], x[t], n);
         basic_add(x[0], x[0], x[t], n);
         coeff_rotate_arbitrary(x+t, scratch, scratch, ru, n);
      }
      else
         // no sqrt2 involved
         coeff_forward_butterfly_bits(x, x+t, scratch, ru >> 1, n);
   }
}


/*
This is the most general forward transform.
It can handle any parameters you throw at it.
*/
void ssfft_fft(mp_limb_t** x, unsigned long t,
               unsigned long m, unsigned long z, unsigned long g,
               unsigned long ru, unsigned long rU, unsigned long n,
               mp_limb_t** scratch)
{
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= (1UL << m));
   FLINT_ASSERT(g >= 1);
   FLINT_ASSERT(g <= (1UL << m));
   FLINT_ASSERT(t >= 1);
   FLINT_ASSERT(ru < rU);
   
   // base cases

   if (m == 2)
   {
      // length 4 cases
      if (z == 4 && g == 4 && IS_WHOLE_LIMB_ROTATION(ru | rU))
      {
         ssfft_fft_size4_z4g4_limbs(x, t, ru / (2*FLINT_BITS_PER_LIMB),
                               rU / (2*FLINT_BITS_PER_LIMB), n, scratch);
         return;
      }
      // todo: would be cleaner to replace this condition with a macro
      else if (IS_WHOLE_BIT_ROTATION(ru | rU))
      {
         ssfft_fft_size4_bits(x, t, z, g, ru >> 1, rU >> 1, n, scratch);
         return;
      }
   }
   
   if (m == 3)
   {
      // length 8 cases
      if (z == 8 && g == 8 && IS_WHOLE_LIMB_ROTATION(ru | rU))
      {
         ssfft_fft_size8_z8g8_limbs(x, t, ru / (2*FLINT_BITS_PER_LIMB),
                                    rU / (2*FLINT_BITS_PER_LIMB), n, scratch);
         return;
      }
      else if (IS_WHOLE_BIT_ROTATION(ru | rU))
      {
         ssfft_fft_size8_bits(x, t, z, g, ru >> 1, rU >> 1, n, scratch);
         return;
      }
   }

   if (m <= 1)
   {
      if (m == 1)
         ssfft_fft_size2(x, t, z, g, ru, n, scratch);
      return;
   }

   // factoring case

   unsigned long k, s;
   mp_limb_t** y;

   unsigned long m2 = m >> 1;
   unsigned long m1 = m - m2;
   unsigned long M2 = 1UL << m2;

   unsigned long g_rows = g >> m2;
   unsigned long z_rows = z >> m2;

   unsigned long m2_mask = M2 - 1;
   unsigned long z_cols = z & m2_mask;
   unsigned long g_cols = g & m2_mask;
   
   unsigned long t_M2 = t << m2;
   unsigned long rU_M2 = rU << m2;
   unsigned long ru_M1 = ru << m1;
   unsigned long rU_M1 = rU << m1;

   if (g_cols)
   {
      for (k = 0, y = x, s = ru; k < z_cols; k++, y += t, s += rU)
         ssfft_fft(y, t_M2, m1, z_rows + 1, g_rows + 1, s, rU_M2, n, scratch);
      
      if (z_rows)
      {
         for (k = z_cols; k < M2; k++, y += t, s += rU)
            ssfft_fft(y, t_M2, m1, z_rows, g_rows + 1, s, rU_M2, n, scratch);

         z_cols = M2;
      }
      
      for (k = 0, y = x; k < g_rows; k++, y += t_M2)
         ssfft_fft(y, t, m2, z_cols, M2, ru_M1, rU_M1, n, scratch);
      ssfft_fft(y, t, m2, z_cols, g_cols, ru_M1, rU_M1, n, scratch);
   }
   else
   {
      for (k = 0, y = x, s = ru; k < z_cols; k++, y += t, s += rU)
         ssfft_fft(y, t_M2, m1, z_rows + 1, g_rows, s, rU_M2, n, scratch);

      if (z_rows)
      {
         for (k = z_cols; k < M2; k++, y += t, s += rU)
            ssfft_fft(y, t_M2, m1, z_rows, g_rows, s, rU_M2, n, scratch);

         z_cols = M2;
      }
      
      for (k = 0, y = x; k < g_rows; k++, y += t_M2)
         ssfft_fft(y, t, m2, z_cols, M2, ru_M1, rU_M1, n, scratch);
   }
}


/* ============================================================================

   Experimental multithreaded version of ssfft_fft.
   
=============================================================================*/

#if ENABLE_THREAD_CODE
typedef struct job_list_t
{
   pthread_mutex_t lock;
   void** jobs;
   void (*func)(void*);
   unsigned long num_jobs;
   unsigned long next_job;
} job_list_t;


void* execute_parallel_worker(void* arg)
{
   job_list_t* job_list = (job_list_t*) arg;
   
   while (1)
   {
      // get a lock on the job queue
      pthread_mutex_lock(&job_list->lock);
      // if a job is available, grab it, update the queue and run the job
      if (job_list->next_job < job_list->num_jobs)
      {
         void* my_job = job_list->jobs[job_list->next_job];
         job_list->next_job++;
         pthread_mutex_unlock(&job_list->lock);
         job_list->func(my_job);
      }
      else
      {
         // no jobs left, so quit this thread
         pthread_mutex_unlock(&job_list->lock);
         return NULL;
      }
   }
}
#endif

/*
Calls func(arg[0]), func(arg[1]), ..., func(arg[n-1]) in parallel.
These are all assumed to be independent tasks, and can be called in any order.

Currently this is just a stub, which calls them in sequence.
*/
void execute_parallel(void (*func)(void*), void** arg, unsigned long n)
{
   unsigned long i;

#if !ENABLE_THREAD_CODE
   // plain boring sequential version
   for (i = 0; i < n; i++)
      func(arg[i]);
#else
   // multithreaded version

   job_list_t job_list;
   pthread_mutex_init(&job_list.lock, NULL);
   job_list.jobs = arg;
   job_list.func = func;
   job_list.num_jobs = n;
   job_list.next_job = 0;

   // start up a bunch of worker threads
   pthread_t threads[NUM_THREADS-1];
   for (i = 0; i < NUM_THREADS-1; i++)
      pthread_create(&threads[i], NULL, execute_parallel_worker,
                     (void*)(&job_list));

   // make this thread a worker too
   execute_parallel_worker((void*)(&job_list));

   // wait for the other threads to catch up
   for (i = 0; i < NUM_THREADS-1; i++)
      pthread_join(threads[i], NULL);

   // clean up
   pthread_mutex_destroy(&job_list.lock);
#endif
}



typedef struct fft_job_t
{
   mp_limb_t** x;
   mp_limb_t** scratch;
   unsigned long m, t, z, g, ru, rU, n;
} fft_job_t;


void ssfft_fft_run_job(void* arg)
{
   fft_job_t* job = (fft_job_t*) arg;
   ssfft_fft(job->x, job->t, job->m, job->z, job->g, job->ru, job->rU, job->n,
             job->scratch);
}


/*
Length must be at least 4 (i.e. m >= 2).
*/
void ssfft_fft_threaded(
               mp_limb_t** x, unsigned long t,
               unsigned long m, unsigned long z, unsigned long g,
               unsigned long ru, unsigned long rU, unsigned long n)
{
   FLINT_ASSERT(m >= 2);

   unsigned long i, j, k, s;
   mp_limb_t** y;

   unsigned long M = 1 << m;
   unsigned long m2 = m >> 1;
   unsigned long m1 = m - m2;
   unsigned long M2 = 1UL << m2;

   unsigned long g_rows = g >> m2;
   unsigned long z_rows = z >> m2;

   unsigned long m2_mask = M2 - 1;
   unsigned long z_cols = z & m2_mask;
   unsigned long g_cols = g & m2_mask;
   
   unsigned long t_M2 = t << m2;
   unsigned long rU_M2 = rU << m2;
   unsigned long ru_M1 = ru << m1;
   unsigned long rU_M1 = rU << m1;
   
   // allocate space for scratch buffers
   unsigned long dim = (m1 > m2) ? (1UL << m1) : (1UL << m2);
   mp_limb_t* scratch_bufs =
              (mp_limb_t*) malloc((n+1) * dim * sizeof(mp_limb_t));
   mp_limb_t* scratch_bufs_end = scratch_bufs + (n+1) * dim;

   mp_limb_t** scratch_buf_ptrs =
              (mp_limb_t**) malloc(dim * sizeof(mp_limb_t*));
   scratch_buf_ptrs[0] = scratch_bufs;
   for (i = 1; i < dim; i++)
      scratch_buf_ptrs[i] = scratch_buf_ptrs[i-1] + (n+1);

   // allocate space for job descriptions and pointers to them
   fft_job_t* jobs = (fft_job_t*) malloc(dim * sizeof(fft_job_t));
   fft_job_t** job_ptrs = (fft_job_t**) malloc(dim * sizeof(fft_job_t*));
   for (i = 0; i < dim; i++)
      job_ptrs[i] = jobs + i;
   unsigned long job_count = 0;

   // generate jobs for column transforms
   unsigned long k_count = z_rows ? M2 : z_cols;
   for (k = 0, y = x, s = ru; k < k_count; k++, y += t, s += rU)
   {
      FLINT_ASSERT(job_count < dim);
      fft_job_t* job = jobs + job_count;
      job->x = y;
      job->t = t_M2;
      job->m = m1;
      job->z = (k < z_cols) ? (z_rows + 1) : z_rows;
      job->g = g_cols ? (g_rows + 1) : g_rows;
      job->ru = s;
      job->rU = rU_M2;
      job->n = n;
      job->scratch = scratch_buf_ptrs + job_count;
      job_count++;
   }
   if (z_rows)
      z_cols = M2;

   // run column transforms
   execute_parallel(ssfft_fft_run_job, (void**) job_ptrs, job_count);

   // generate jobs for row transforms
   job_count = 0;
   k_count = g_cols ? (g_rows + 1) : g_rows;
   for (k = 0, y = x; k < k_count; k++, y += t_M2)
   {
      FLINT_ASSERT(job_count < dim);
      fft_job_t* job = jobs + job_count;
      job->x = y;
      job->t = t;
      job->m = m2;
      job->z = z_cols;
      job->g = (k < g_rows) ? M2 : g_cols;
      job->ru = ru_M1;
      job->rU = rU_M1;
      job->n = n;
      job->scratch = scratch_buf_ptrs + job_count;
      job_count++;
   }

   // run row transforms
   execute_parallel(ssfft_fft_run_job, (void**) job_ptrs, job_count);

   // rearrange data so that all outputs are back in the input buffers, and
   // all scratch buffers are unused again
   i = j = 0;
   while (1)
   {
      // find next input buffer that's in a scratch slot
      for (; i < M && !(x[i] >= scratch_bufs && x[i] < scratch_bufs_end); i++);
      // bail out if all buffers have been examined
      if (i == M)
         break;
      // find next scratch buffer that's in an input slot
      for (; j < dim && (scratch_buf_ptrs[j] >= scratch_bufs &&
                         scratch_buf_ptrs[j] < scratch_bufs_end); j++);

      // copy the data from the scratch slot to the input slot
      swap_limb_ptrs(x + i, scratch_buf_ptrs + j);
      basic_copy(x[i], scratch_buf_ptrs[j], n);
      
      i++; j++;
   }
   
   // clean up
   free(job_ptrs);
   free(jobs);
   free(scratch_buf_ptrs);
   free(scratch_bufs);
}



/* ============================================================================

   Inverse transform.

This is pretty much the inverse of the above forward transform, with certain
complications arising from the truncation business.

Let a_i and b_i be as above. The inverse transform requires the input array
to contain

   b_0, b_1, ..., b_{g-1}, 2^m a_g, 2^m a_{g+1}, ... , 2^m a_{z-1}.

Coefficients from z onwards are ignored, and are assumed to be zero.
It is required that 0 < z <= 2^m and also that z >= g.

The output consists of either

   2^m a_0, 2^m a_1, ..., 2^m a_{g-1}        (if the flag e is not set),
or
   2^m a_0, 2^m a_1, ..., 2^m a_{g-1}, b_g   (if the flag e is set);

the remainder of the array becomes garbage. Note that g is allowed to be zero
(this still makes sense if e is set). If g == 2^m then e is ignored.

=============================================================================*/

// ----------------------------------------------------------------------------
// length 4 transforms

/*
Assumes M = z = g = 4, all rotations by full limbs
*/
void ssfft_ifft_size4_z4g4_limbs(
      mp_limb_t** x, unsigned long t,
      unsigned long ru_limbs, unsigned long rU_limbs,
      unsigned long n, mp_limb_t** scratch)
{
   FLINT_ASSERT(ru_limbs < rU_limbs);

   if (ru_limbs)
   {
      coeff_inverse_butterfly_limbs(x+2*t, x+3*t, scratch, 2*ru_limbs, n);
      coeff_inverse_butterfly_limbs(x, x+t, scratch, 2*ru_limbs, n);
      coeff_inverse_butterfly_limbs(x+t, x+3*t, scratch,
                                    ru_limbs + rU_limbs, n);
      coeff_inverse_butterfly_limbs(x, x+2*t, scratch, ru_limbs, n);
   }
   else
   {
      coeff_inverse_simple_butterfly(x+2*t, x+3*t, scratch, n);
      coeff_inverse_simple_butterfly(x, x+t, scratch, n);
      coeff_inverse_butterfly_limbs(x+t, x+3*t, scratch, rU_limbs, n);
      coeff_inverse_simple_butterfly(x, x+2*t, scratch, n);
   }
}


/*
Assumes M = z = g = 8, all rotations by full limbs
*/
void ssfft_ifft_size8_z8g8_limbs(
      mp_limb_t** x, unsigned long t,
      unsigned long ru_limbs, unsigned long rU_limbs,
      unsigned long n, mp_limb_t** scratch)
{
   FLINT_ASSERT(ru_limbs < rU_limbs);

   // recurse into two halves
   ssfft_ifft_size4_z4g4_limbs(x+4*t, t, 2*ru_limbs, 2*rU_limbs, n, scratch);
   ssfft_ifft_size4_z4g4_limbs(x, t, 2*ru_limbs, 2*rU_limbs, n, scratch);

   // Last row of butterflies.
   if (ru_limbs)
      coeff_inverse_butterfly_limbs(x, x+4*t, scratch, ru_limbs, n);
   else
      coeff_inverse_simple_butterfly(x, x+4*t, scratch, n);
   coeff_inverse_butterfly_limbs(x+t, x+5*t, scratch, ru_limbs + rU_limbs, n);
   coeff_inverse_butterfly_limbs(x+2*t, x+6*t, scratch,
                                 ru_limbs + 2*rU_limbs, n);
   coeff_inverse_butterfly_limbs(x+3*t, x+7*t, scratch,
                                 ru_limbs + 3*rU_limbs, n);
}


/*
The following versions all assume M = 4, and z and g and e are given in the
name of the function. Rotations are by an integral number of bits.
*/

// input  [4*a0]
// output [b0]
void ssfft_ifft_size4_z1g0e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := x0 / 4
   basic_unrotate_bits(x[0], x[0], 2, n);
}

// input  [4*a0, 4*a1]
// output [b0]
void ssfft_ifft_size4_z2g0e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := x0 + x1  = 4*(a0 + a1)
   basic_add(x[0], x[0], x[t], n);
   // x0 := x0 / 4  = a0 + a1   = b0
   basic_unrotate_bits(x[0], x[0], 2, n);
}

// input  [4*a0, 4*a1, 4*a2]
// output [b0]
void ssfft_ifft_size4_z3g0e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := x0 + x1 + x2  = 4*(a0 + a1 + a2)
   basic_add(x[0], x[0], x[t], n);
   basic_add(x[0], x[0], x[2*t], n);
   // x0 := x0 / 4   = a0 + a1 + a2   = b0
   basic_unrotate_bits(x[0], x[0], 2, n);
}

// input  [4*a0, 4*a1, 4*a2, 4*a3]
// output [b0]
void ssfft_ifft_size4_z4g0e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := x0 + x1 + x2 + x3  = 4*(a0 + a1 + a2 + a3)
   basic_add(x[0], x[0], x[t], n);
   basic_add(x[0], x[0], x[2*t], n);
   basic_add(x[0], x[0], x[3*t], n);
   // x0 := x0 / 4  = a0 + a1 + a2 + a3  = b0
   basic_unrotate_bits(x[0], x[0], 2, n);
}

// input  [b0]
// output [4*a0]
void ssfft_ifft_size4_z1g1e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*x0  = 4*a0
   basic_rotate_bits(x[0], x[0], 2, n);
}

// input  [b0]
// output [4*a0, b1]
void ssfft_ifft_size4_z1g1e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x1 := w^(2u) x0  = w^(2u) a0  = b1
   coeff_rotate_bits(x+t, x, 2*ru_bits, n);
   // x0 := 4*x0  = 4*a0
   basic_rotate_bits(x[0], x[0], 2, n);
}

// input  [b0, 4*a1]
// output [4*a0]
void ssfft_ifft_size4_z2g1e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*x0  = 4*b0  = 4*(a0 + a1)
   basic_rotate_bits(x[0], x[0], 2, n);
   // x0 := x0 - x1  = 4*a0
   basic_sub(x[0], x[0], x[t], n);
   basic_fast_reduce(x[0], n);    // just to be safe
}

// input  [b0, 4*a1]
// output [4*a0, b1]
void ssfft_ifft_size4_z2g1e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // scr := x1 / 4  = a1
   basic_unrotate_bits(*scratch, x[t], 2, n);
   // x0 := x0 - scr  = (a0 + a1) - a1  = a0
   basic_sub(x[0], x[0], *scratch, n);
   // x1 := w^(2u) (x0 - scr)  = b1
   coeff_sub_rotate_bits(x+t, x, scratch, 2*ru_bits, n);
   // x0 := 4*x0  = 4*a0
   basic_rotate_bits(x[0], x[0], 2, n);
   basic_fast_reduce(x[0], n);    // just to be safe
}

// input  [b0, 4*a1, 4*a2]
// output [4*a0]
void ssfft_ifft_size4_z3g1e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*x0  = 4*b0  = 4*(a0 + a1 + a2)
   basic_rotate_bits(x[0], x[0], 2, n);
   // x0 := x0 - x1 - x2  = 4*a0
   basic_sub(x[0], x[0], x[t], n);
   basic_sub(x[0], x[0], x[2*t], n);
   basic_fast_reduce(x[0], n);    // just to be safe
}

// input  [b0, 4*a1, 4*a2]
// output [4*a0, b1]
void ssfft_ifft_size4_z3g1e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // scr := x1 / 4  = a1
   basic_unrotate_bits(*scratch, x[t], 2, n);
   // x0 := x0 - scr  = a0 + a2
   basic_sub(x[0], x[0], *scratch, n);
   // x1 := w^(2u) (x0 - scr)  = w^(2u) (a0 + a2 - a1)  = b1
   coeff_sub_rotate_bits(x+t, x, scratch, 2*ru_bits, n);
   // x0 := 4*x0  = 4*a0 + 4*a2
   basic_rotate_bits(x[0], x[0], 2, n);
   // x0 := x0 - x2  = 4*a0
   basic_sub(x[0], x[0], x[2*t], n);
   basic_fast_reduce(x[0], n);    // just to be safe
}

// input  [b0, 4*a1, 4*a2, 4*a3]
// output [4*a0]
void ssfft_ifft_size4_z4g1e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*x0  = 4*b0  = 4*(a0 + a1 + a2 + a3)
   basic_rotate_bits(x[0], x[0], 2, n);
   // x0 := x0 - x1 - x2 - x3  = 4*a0
   basic_sub(x[0], x[0], x[t], n);
   basic_sub(x[0], x[0], x[2*t], n);
   basic_sub(x[0], x[0], x[3*t], n);
   basic_fast_reduce(x[0], n);    // just to be safe
}

// input  [b0, 4*a1, 4*a2, 4*a3]
// output [4*a0, b1]
void ssfft_ifft_size4_z4g1e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x3 := x1 + x3  = 4*a1 + 4*a3
   basic_add(x[3*t], x[t], x[3*t], n);
   // x3 := x3 / 4  = a1 + a3
   basic_unrotate_bits(x[3*t], x[3*t], 2, n);
   // x0 := x0 - x3  = (a0 + a1 + a2 + a3) - (a1 + a3)  = (a0 + a2)
   basic_sub(x[0], x[0], x[3*t], n);
   // x1 := w^(2u) (x0 - x3)  = w^(2u) (a0 + a2 - a1 - a3)
   coeff_sub_rotate_bits(x+t, x, x+3*t, 2*ru_bits, n);
   // x0 := 4*x0  = 4*a0 + 4*a2
   basic_rotate_bits(x[0], x[0], 2, n);
   // x0 := x0 - x2  = 4*a0
   basic_sub(x[0], x[0], x[2*t], n);
   basic_fast_reduce(x[0], n);    // just to be safe
}

// input  [b0, b1]
// output [4*a0, 4*a1]
void ssfft_ifft_size4_z2g2e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0
   // x1 := 2*a1
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   // x0 := 2*x0 = 4*a0
   basic_add(x[0], x[0], x[0], n);
   // x1 := 2*x1 = 4*a1
   basic_add(x[t], x[t], x[t], n);
}

// input  [b0, b1]
// output [4*a0, 4*a1, b2]
void ssfft_ifft_size4_z2g2e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*a0
   // x1 := 4*a1
   ssfft_ifft_size4_z2g2e0(x, t, ru_bits, rU_bits, n, scratch);

   // x2 := 1/4 w^u x0  = w^u a0
   if (ru_bits >= 2)
      coeff_rotate_bits(x+2*t, x, ru_bits - 2, n);
   else
      basic_unrotate_bits(x[2*t], x[0], 2 - ru_bits, n);
   
   // x3 := 1/4 w^(u+U) x1  = w^(u+U) a1
   unsigned long s = ru_bits + rU_bits;
   if (s >= 2)
      coeff_rotate_bits(x+3*t, x+t, s-2, n);
   else
      basic_unrotate_bits(x[3*t], x[t], 2-s, n);
      
   // x2 := x2 + x3  = w^u a0 + w^(u+U) a1  = b2
   basic_add(x[2*t], x[2*t], x[3*t], n);
   basic_fast_reduce(x[2*t], n);    // just to be safe
}

// input  [b0, b1, 4*a2]
// output [4*a0, 4*a1]
void ssfft_ifft_size4_z3g2e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0 + 2*a2
   // x1 := 2*a1
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   // x0 := 2*x0 - x2  = 4*a0
   coeff_cross_butterfly2(x, x+2*t, n);
   // x1 := 2*x1  = 4*a1
   basic_add(x[t], x[t], x[t], n);
}

// input  [b0, b1, 4*a2]
// output [4*a0, 4*a1, b2]
void ssfft_ifft_size4_z3g2e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*a0
   // x1 := 4*a1
   ssfft_ifft_size4_z3g2e0(x, t, ru_bits, rU_bits, n, scratch);

   // scr := 1/4 w^u (x0 - x2)  = w^u (a0 - a2)
   if (ru_bits >= 2)
      coeff_sub_rotate_bits(scratch, x, x+2*t, ru_bits - 2, n);
   else
   {
      basic_sub(*scratch, x[0], x[2*t], n);
      if (ru_bits <= 1)
         basic_unrotate_bits(*scratch, *scratch, 2 - ru_bits, n);
   }
   
   // x2 := 1/4 w^(u+U) x1  = w^(u+U) a1
   FLINT_ASSERT(rU_bits >= 2);
   coeff_rotate_bits(x+2*t, x+t, ru_bits + rU_bits - 2, n);

   // x2 := x2 + scr  = w^u (a0 - a2) + w^(u+U) a1   = b2
   basic_add(x[2*t], x[2*t], *scratch, n);
   basic_fast_reduce(x[2*t], n);    // just to be safe
}

// input  [b0, b1, 4*a2, 4*a3]
// output [4*a0, 4*a1]
void ssfft_ifft_size4_z4g2e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0 + 2*a2
   // x1 := 2*a1 + 2*a3
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   // x0 := 2*x0 - x2 = 4*a0
   coeff_cross_butterfly2(x, x+2*t, n);
   // x1 := 2*x1 - x3  = 4*a1
   coeff_cross_butterfly2(x+t, x+3*t, n);
}

// input  [b0, b1, 4*a2, 4*a3]
// output [4*a0, 4*a1, b2]
void ssfft_ifft_size4_z4g2e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 4*a0
   // x1 := 4*a1
   ssfft_ifft_size4_z4g2e0(x, t, ru_bits, rU_bits, n, scratch);
   
   // scr := 1/4 w^u (x0 - x2)  = w^u (a0 - a2)
   if (ru_bits >= 2)
      coeff_sub_rotate_bits(scratch, x, x+2*t, ru_bits - 2, n);
   else
   {
      basic_sub(*scratch, x[0], x[2*t], n);
      if (ru_bits <= 1)
         basic_unrotate_bits(*scratch, *scratch, 2 - ru_bits, n);
   }
   
   // x2 := 1/4 w^(u+U) (x1 - x3)  = w^(u+U) (a1 - a3)
   FLINT_ASSERT(rU_bits >= 2);
   coeff_sub_rotate_bits(x+2*t, x+t, x+3*t, ru_bits + rU_bits - 2, n);
   
   // x2 := x2 + scr  = b2
   basic_add(x[2*t], x[2*t], *scratch, n);
   basic_fast_reduce(x[2*t], n);    // just to be safe
}

// input  [b0, b1, b2]
// output [4*a0, 4*a1, 4*a2]
void ssfft_ifft_size4_z3g3e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0 + 2*a2
   // x1 := 2*a1
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);

   // scr := 1/2 w^(u+U) x1  = w^(u+U) a1
   FLINT_ASSERT(rU_bits >= 1);
   coeff_rotate_bits(scratch, x+t, ru_bits + rU_bits - 1, n);
   // x2 := x2 - scr  = w^u (a0 - a2)
   basic_sub(x[2*t], x[2*t], *scratch, n);
   
   // x0 := 4*a0
   // x2 := 4*a2
   if (ru_bits)
      coeff_inverse_butterfly_bits(x, x+2*t, scratch, ru_bits - 1, n);
   else
   {
      basic_add(x[2*t], x[2*t], x[2*t], n);
      coeff_inverse_simple_butterfly(x, x+2*t, scratch, n);
   }

   basic_fast_reduce(x[0], n);    // just to be safe
   basic_fast_reduce(x[2*t], n);    // just to be safe
   
   // x1 := 2*x1  = 4*a1
   basic_add(x[t], x[t], x[t], n);
}

// input  [b0, b1, b2]
// output [4*a0, 4*a1, 4*a2, b3]
void ssfft_ifft_size4_z3g3e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0 + 2*a2
   // x1 := 2*a1
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);

   // scr := 1/2 w^(u+U) x1  = w^(u+U) a1
   FLINT_ASSERT(rU_bits >= 1);
   coeff_rotate_bits(scratch, x+t, ru_bits + rU_bits - 1, n);
   // x2 := x2 - scr  = w^u (a0 - a2)
   basic_sub(x[2*t], x[2*t], *scratch, n);

   // x3 := w^(2u) (x2 - scr)  = w^(3u) (a0 - a2) - w^(3u+U) a1  = b3
   coeff_sub_rotate_bits(x+3*t, x+2*t, scratch, 2*ru_bits, n);
   basic_fast_reduce(x[3*t], n);    // just to be safe

   // x0 := 4*a0
   // x2 := 4*a2
   if (ru_bits)
      coeff_inverse_butterfly_bits(x, x+2*t, scratch, ru_bits - 1, n);
   else
   {
      basic_add(x[2*t], x[2*t], x[2*t], n);
      coeff_inverse_simple_butterfly(x, x+2*t, scratch, n);
   }

   basic_fast_reduce(x[0], n);    // just to be safe
   basic_fast_reduce(x[2*t], n);    // just to be safe
   
   // x1 := 2*x1  = 4*a1
   basic_add(x[t], x[t], x[t], n);
}

// input  [b0, b1, b2, 4*a3]
// output [4*a0, 4*a1, 4*a2]
void ssfft_ifft_size4_z4g3e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0 + 2*a2
   // x1 := 2*a1 + 2*a3
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   
   // scr := 1/2 w^(u+U) (x1 - x3)  = w^(u+U) (a1 - a3)
   FLINT_ASSERT(rU_bits >= 1);
   coeff_sub_rotate_bits(scratch, x+t, x+3*t, ru_bits + rU_bits - 1, n);

   // x2 := x2 - scr  = w^u (a0 - a2)
   basic_sub(x[2*t], x[2*t], *scratch, n);
   basic_fast_reduce(x[2*t], n);    // just to be safe

   // x0 := 4*a0
   // x2 := 4*a2
   if (ru_bits)
      coeff_inverse_butterfly_bits(x, x+2*t, scratch, ru_bits - 1, n);
   else
   {
      basic_add(x[2*t], x[2*t], x[2*t], n);
      coeff_inverse_simple_butterfly(x, x+2*t, scratch, n);
   }

   // x1 := 2*x1 - x3  = 4*a1
   coeff_cross_butterfly2(x+t, x+3*t, n);
}

// input  [b0, b1, b2, 4*a3]
// output [4*a0, 4*a1, 4*a2, b3]
void ssfft_ifft_size4_z4g3e1(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   // x0 := 2*a0 + 2*a2
   // x1 := 2*a1 + 2*a3
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   
   // scr := 1/2 w^(u+U) (x1 - x3)  = w^(u+U) (a1 - a3)
   FLINT_ASSERT(rU_bits >= 1);
   coeff_sub_rotate_bits(scratch, x+t, x+3*t, ru_bits + rU_bits - 1, n);

   // x2 := x2 - scr  = w^u (a0 - a2)
   basic_sub(x[2*t], x[2*t], *scratch, n);
   basic_fast_reduce(x[2*t], n);    // just to be safe

   // x1 := 2*x1 - x3  = 4*a1
   coeff_cross_butterfly2(x+t, x+3*t, n);

   // x3 := w^(2u) (x2 - scr)  = w^(3u) (a0 - a2) - w^(3u+U) (a1 - a3)  = b3
   coeff_sub_rotate_bits(x+3*t, x+2*t, scratch, 2*ru_bits, n);
   basic_fast_reduce(x[3*t], n);    // just to be safe

   // x0 := 4*a0
   // x2 := 4*a2
   if (ru_bits)
      coeff_inverse_butterfly_bits(x, x+2*t, scratch, ru_bits - 1, n);
   else
   {
      basic_add(x[2*t], x[2*t], x[2*t], n);
      coeff_inverse_simple_butterfly(x, x+2*t, scratch, n);
   }
}

// input  [b0, b1, b2, b3]
// output [4*a0, 4*a1, 4*a2, 4*a3]
void ssfft_ifft_size4_z4g4e0(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   coeff_inverse_butterfly_bits(x+2*t, x+3*t, scratch, 2*ru_bits, n);
   coeff_inverse_butterfly_bits(x, x+t, scratch, 2*ru_bits, n);
   coeff_inverse_butterfly_bits(x+t, x+3*t, scratch, ru_bits + rU_bits, n);
   coeff_inverse_butterfly_bits(x, x+2*t, scratch, ru_bits, n);
}

// ssfft_ifft_size4_func_t is the type of the functions above
typedef void (*ssfft_ifft_size4_func_t)
                        (mp_limb_t**, unsigned long, unsigned long,
                         unsigned long, unsigned long, mp_limb_t**);


// a placeholder to pick up invalid combinations of z, g, e during debugging
void ssfft_ifft_size4_broken(mp_limb_t** x, unsigned long t,
                           unsigned long ru_bits, unsigned long rU_bits,
                           unsigned long n, mp_limb_t** scratch)
{
   FLINT_ASSERT(0);
}


// A lookup table for the functions above, indexed by e and z-1 and g
ssfft_ifft_size4_func_t
ssfft_ifft_size4_table[2][4][8] = {
   {
      {ssfft_ifft_size4_broken, ssfft_ifft_size4_z1g1e0,
       ssfft_ifft_size4_broken, ssfft_ifft_size4_broken,
       ssfft_ifft_size4_broken},

      {ssfft_ifft_size4_broken, ssfft_ifft_size4_z2g1e0,
       ssfft_ifft_size4_z2g2e0, ssfft_ifft_size4_broken,
       ssfft_ifft_size4_broken},

      {ssfft_ifft_size4_broken, ssfft_ifft_size4_z3g1e0,
       ssfft_ifft_size4_z3g2e0, ssfft_ifft_size4_z3g3e0,
       ssfft_ifft_size4_broken},
       
      {ssfft_ifft_size4_broken, ssfft_ifft_size4_z4g1e0,
       ssfft_ifft_size4_z4g2e0, ssfft_ifft_size4_z4g3e0,
       ssfft_ifft_size4_z4g4e0}
   },
   {
      {ssfft_ifft_size4_z1g0e1, ssfft_ifft_size4_z1g1e1,
       ssfft_ifft_size4_broken, ssfft_ifft_size4_broken,
       ssfft_ifft_size4_broken},
      
      {ssfft_ifft_size4_z2g0e1, ssfft_ifft_size4_z2g1e1,
       ssfft_ifft_size4_z2g2e1, ssfft_ifft_size4_broken,
       ssfft_ifft_size4_broken},
       
      {ssfft_ifft_size4_z3g0e1, ssfft_ifft_size4_z3g1e1,
       ssfft_ifft_size4_z3g2e1, ssfft_ifft_size4_z3g3e1,
       ssfft_ifft_size4_broken},
       
      {ssfft_ifft_size4_z4g0e1, ssfft_ifft_size4_z4g1e1,
       ssfft_ifft_size4_z4g2e1, ssfft_ifft_size4_z4g3e1,
       ssfft_ifft_size4_broken}
   }
};


/*
This version assumes M = 4, and rotation by an integral number of bits.
Any 1 <= z <= 4, 0 <= g <= 4, z >= g, and 0 <= e <= 1 is allowed.

todo: those conditions are wrong, fix them
*/
inline void ssfft_ifft_size4_bits(
       mp_limb_t** x, unsigned long t,
       unsigned long z, unsigned long g, int e,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch)
{
   FLINT_ASSERT(ru_bits < rU_bits);
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= 4);
   FLINT_ASSERT(g <= 4);
   FLINT_ASSERT(z >= g);
   FLINT_ASSERT(e <= 1);
   FLINT_ASSERT((g == 0 && e) || (g == 4 && !e) || (g >= 1 && g <= 3));
   
   ssfft_ifft_size4_table[e][z-1][g](x, t, ru_bits, rU_bits, n, scratch);
}


// ----------------------------------------------------------------------------
// length 8 transforms

/*
This version assumes M = 8, and rotation by an integral number of bits.
Any 1 <= z <= 8, z >= g, 0 <= g <= 8 is allowed. (todo: what about e?)

The code here is basically a fully unrolled version of the general
ssfft_ifft, for m1 = 1 and m2 = 2.
*/
void ssfft_ifft_size8_bits(
       mp_limb_t** x, unsigned long t,
       unsigned long z, unsigned long g, int e,
       unsigned long ru_bits, unsigned long rU_bits, unsigned long n,
       mp_limb_t** scratch)
{
   FLINT_ASSERT(ru_bits < rU_bits);
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= 8);
   FLINT_ASSERT(g <= 8);
   FLINT_ASSERT(z >= g);
   FLINT_ASSERT(e <= 1);
   FLINT_ASSERT((g == 0 && e) || (g == 8 && !e) || (g >= 1 && g <= 7));

   unsigned long k, s;
   mp_limb_t** X;
   mp_limb_t** Y;

   if (g >= 4)
   {
      ssfft_ifft_size4_z4g4e0(x, t, 2*ru_bits, 2*rU_bits, n, scratch);
      if (g == 8)
         ssfft_ifft_size4_z4g4e0(x+4*t, t, 2*ru_bits, 2*rU_bits, n, scratch);
   }

   if (g <= 3)
   {
      if (g == 0)
      {
         FLINT_ASSERT(e);
         for (k = 1; k < z; k++)
            basic_add(x[0], x[0], x[k*t], n);
         basic_unrotate_bits(x[0], x[0], 3, n);
      }
      else if (z <= 3)
      {
         // g_rows == 0, z_rows == 0
         for (k = g, X = x+g*t; k < z; k++, X += t)
            basic_unrotate_bits(*X, *X, 1, n);

         ssfft_ifft_size4_bits(x, t, z, g, e,
                               2*ru_bits, 2*rU_bits, n, scratch);

         for (k = 0, X = x; k < g; k++, X += t)
            basic_add(*X, *X, *X, n);
      }
      else if (z <= 7)
      {
         // g_rows == 0, z_rows == 1
         z -= 4;
      
         for (k = g, X = x+g*t, Y = X+4*t; k < z; k++, X += t, Y += t)
         {
            basic_add(*X, *X, *Y, n);
            basic_unrotate_bits(*X, *X, 1, n);
         }
         for (; k < 4; k++, X += t)
            basic_unrotate_bits(*X, *X, 1, n);

         ssfft_ifft_size4_bits(x, t, 4, g, e,
                               2*ru_bits, 2*rU_bits, n, scratch);
                    
         for (k = 0, X = x, Y = X+4*t; k < g && k < z; k++, X += t, Y += t)
            coeff_cross_butterfly2(X, Y, n);
         for (; k < g; k++, X += t)
            basic_add(*X, *X, *X, n);
      }
      else
      {
         // g_rows == 0, z_rows == 2
         for (k = g, X = x+g*t, Y = X+4*t; k < 4; k++, X += t, Y += t)
         {
            basic_add(*X, *X, *Y, n);
            basic_unrotate_bits(*X, *X, 1, n);
         }

         ssfft_ifft_size4_bits(x, t, 4, g, e,
                               2*ru_bits, 2*rU_bits, n, scratch);
                    
         for (k = 0, X = x, Y = X+4*t; k < g; k++, X += t, Y += t)
            coeff_cross_butterfly2(X, Y, n);
      }
   }
   else if (g <= 7)
   {
      if (z <= 7)
      {
         // g_rows == 1, z_rows == 1
         g -= 4;
         z -= 4;

         if (g)
         {
            for (k = g, X = x+g*t, Y = X+4*t, s = ru_bits + g*rU_bits;
                 k < z; k++, X += t, Y += t, s += rU_bits)
            {
               coeff_cross_butterfly_bits(X, Y, scratch, s, n);
            }
            for (; k < 4; k++, X += t, Y += t, s += rU_bits)
            {
               coeff_rotate_bits(Y, X, s, n);
               basic_add(*X, *X, *X, n);
            }

            ssfft_ifft_size4_bits(x+4*t, t, 4, g, e,
                                  2*ru_bits, 2*rU_bits, n, scratch);
                       
            for (k = 0, X = x, Y = X+4*t, s = ru_bits; k < g && k < z;
                 k++, X += t, Y += t, s += rU_bits)
            {
               coeff_inverse_butterfly_bits(X, Y, scratch, s, n);
            }

            for (; k < g; k++, X += t, Y += t)
               coeff_cross_butterfly2(X, Y, n);
         }
         else if (e)
         {
            for (k = 0, X = x, Y = X+4*t, s = ru_bits; k < z;
                 k++, X += t, Y += t, s += rU_bits)
            {
               coeff_cross_butterfly_bits(X, Y, scratch, s, n);
            }
            for (; k < 4; k++, X += t, Y += t, s += rU_bits)
            {
               coeff_rotate_bits(Y, X, s, n);
               basic_add(*X, *X, *X, n);
            }
            
            x += 4*t;
            for (k = 1, X = x+t; k < 4; k++, X += t)
               basic_add(*x, *x, *X, n);
            basic_unrotate_bits(*x, *x, 2, n);
         }
         else
         {
            for (k = 0, X = x, Y = X+4*t; k < z; k++, X += t, Y += t)
               coeff_cross_butterfly2(X, Y, n);
            for (; k < 4; k++, X += t)
               basic_add(*X, *X, *X, n);
         }
      }
      else
      {
         // g_rows == 1, z_rows == 2
         g -= 4;

         if (g)
         {
            for (k = g, X = x+g*t, Y = X+4*t, s = ru_bits + g*rU_bits;
                 k < 4; k++, X += t, Y += t, s += rU_bits)
            {
               coeff_cross_butterfly_bits(X, Y, scratch, s, n);
            }

            ssfft_ifft_size4_bits(x+4*t, t, 4, g, e,
                                  2*ru_bits, 2*rU_bits, n, scratch);
                       
            for (k = 0, X = x, Y = X+4*t, s = ru_bits; k < g;
                 k++, X += t, Y += t, s += rU_bits)
            {
               coeff_inverse_butterfly_bits(X, Y, scratch, s, n);
            }
         }
         else if (e)
         {
            for (k = 0, X = x, Y = X+4*t, s = ru_bits; k < 4;
                 k++, X += t, Y += t, s += rU_bits)
            {
               coeff_cross_butterfly_bits(X, Y, scratch, s, n);
            }
            
            x += 4*t;
            for (k = 1, X = x+t; k < 4; k++, X += t)
               basic_add(*x, *x, *X, n);
            basic_unrotate_bits(*x, *x, 2, n);
         }
         else
         {
            for (k = 0, X = x, Y = X+4*t; k < 4; k++, X += t, Y += t)
               coeff_cross_butterfly2(X, Y, n);
         }
      }
   }
   else
   {
      // g_rows == 2 and z_rows == 2
      for (k = 0, X = x, Y = X+4*t, s = ru_bits; k < 4;
           k++, X += t, Y += t, s += rU_bits)
      {
         coeff_inverse_butterfly_bits(X, Y, scratch, s, n);
      }
   }
}


void ssfft_ifft_size2(mp_limb_t** x, unsigned long t,
                      unsigned long z, unsigned long g, int e,
                      unsigned long ru, unsigned long n,
                      mp_limb_t** scratch)
{
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= 2);
   FLINT_ASSERT(g <= 2);
   FLINT_ASSERT(z >= g);
   FLINT_ASSERT(e <= 1);
   FLINT_ASSERT((g == 0 && e) || (g == 2 && !e) || (g == 1));

   if (g == 0)
   {
      FLINT_ASSERT(e);
      if (z == 2)
         basic_add(x[0], x[0], x[t], n);
      basic_unrotate_bits(x[0], x[0], 1, n);
   }
   else if (g == 1)
   {
      if (z == 1)
      {
         if (e)
            coeff_rotate_arbitrary(x+t, x, scratch, ru, n);
         basic_add(x[0], x[0], x[0], n);
      }
      else  // z == 2
      {
         if (e)
         {
            basic_sub(*scratch, x[0], x[t], n);
            basic_add(x[0], x[0], *scratch, n);
            coeff_rotate_arbitrary(x+t, scratch, scratch, ru, n);
            basic_fast_reduce(x[0], n);    // just to be safe
            basic_fast_reduce(x[t], n);    // just to be safe
         }
         else
            coeff_cross_butterfly2(x, x+t, n);
      }
   }
   else  // g == 2
   {
      if (ru & 1)
      {
         // sqrt2 involved
         coeff_rotate_arbitrary(scratch, x+t, x+t,
                                2*n*FLINT_BITS_PER_LIMB - ru, n);
         basic_fast_reduce(*scratch, n);     // just to be safe
         basic_add(x[t], x[0], *scratch, n);
         basic_sub(x[0], x[0], *scratch, n);
      }
      else
         // no sqrt2 involved
         coeff_inverse_butterfly_bits(x, x+t, scratch, ru >> 1, n);
   }
}


/*
This is the most general inverse transform.
It can handle any parameters you throw at it.
*/
void ssfft_ifft(mp_limb_t** x, unsigned long t,
                unsigned long m, unsigned long z, unsigned long g, int e,
                unsigned long ru, unsigned long rU, unsigned long n,
                mp_limb_t** scratch)
{
   FLINT_ASSERT(ru < rU);
   FLINT_ASSERT(z >= 1);
   FLINT_ASSERT(z <= (1UL << m));
   FLINT_ASSERT(g <= (1UL << m));
   FLINT_ASSERT(z >= g);
   FLINT_ASSERT(e <= 1);
   FLINT_ASSERT((g == 0 && e) || (g == (1UL << m) && !e) ||
                (g >= 1 && g < (1UL << m)));
   
   
   // base cases
   if (m == 2)
   {
      if (z == 4 && g == 4 && IS_WHOLE_LIMB_ROTATION(ru | rU))
      {
         ssfft_ifft_size4_z4g4_limbs(x, t, ru / (2*FLINT_BITS_PER_LIMB),
                                     rU / (2*FLINT_BITS_PER_LIMB), n, scratch);
         return;
      }
      else if (IS_WHOLE_BIT_ROTATION(ru | rU))
      {
         ssfft_ifft_size4_bits(x, t, z, g, e, ru >> 1, rU >> 1, n, scratch);
         return;
      }
   }
   
   if (m == 3)
   {
      if (z == 8 && g == 8 && IS_WHOLE_LIMB_ROTATION(ru | rU))
      {
         ssfft_ifft_size8_z8g8_limbs(x, t, ru / (2*FLINT_BITS_PER_LIMB),
                                     rU / (2*FLINT_BITS_PER_LIMB), n, scratch);
         return;
      }
      else if (IS_WHOLE_BIT_ROTATION(ru | rU))
      {
         ssfft_ifft_size8_bits(x, t, z, g, e, ru >> 1, rU >> 1, n, scratch);
         return;
      }
   }
      
   if (m <= 1)
   {
      if (m == 1)
         ssfft_ifft_size2(x, t, z, g, e, ru, n, scratch);
      return;
   }
   
   // factoring case
      
   unsigned long k, s;
   mp_limb_t** y;

   unsigned long m1 = m >> 1;
   unsigned long m2 = m - m1;
   unsigned long M2 = 1UL << m2;
   
   unsigned long g_rows = g >> m2;
   unsigned long z_rows = z >> m2;

   unsigned long m2_mask = M2 - 1;
   unsigned long z_cols = z & m2_mask;
   unsigned long g_cols = g & m2_mask;
   
   unsigned long t_M2 = t << m2;
   unsigned long rU_M2 = rU << m2;
   unsigned long ru_M1 = ru << m1;
   unsigned long rU_M1 = rU << m1;

   for (k = 0, y = x; k < g_rows; k++, y += t_M2)
      ssfft_ifft(y, t, m2, M2, M2, 0, ru_M1, rU_M1, n, scratch);

   if (g_cols)
   {
      for (k = g_cols, y = x + g_cols*t, s = ru + g_cols*rU;
           k < z_cols; k++, y += t, s += rU)
      {
         ssfft_ifft(y, t_M2, m1, z_rows + 1, g_rows, 1, s, rU_M2, n, scratch);
      }

      if (z_rows)
         for (; k < M2; k++, y += t, s += rU)
            ssfft_ifft(y, t_M2, m1, z_rows, g_rows, 1, s, rU_M2, n, scratch);

      ssfft_ifft(x + (g_rows << m2)*t, t, m2, (z_rows ? M2 : z_cols),
                 g_cols, e, ru_M1, rU_M1, n, scratch);
                 
      for (k = 0, y = x, s = ru; k < g_cols && k < z_cols;
           k++, y += t, s += rU)
      {
         ssfft_ifft(y, t_M2, m1, z_rows + 1, g_rows + 1, 0,
                    s, rU_M2, n, scratch);
      }
      if (z_rows)
      {
         for (; k < g_cols; k++, y += t, s += rU)
            ssfft_ifft(y, t_M2, m1, z_rows, g_rows + 1, 0,
                       s, rU_M2, n, scratch);
      }
   }
   else
   {
      for (k = 0, y = x, s = ru; k < z_cols; k++, y += t, s += rU)
         ssfft_ifft(y, t_M2, m1, z_rows + 1, g_rows, e, s, rU_M2, n, scratch);

      if (z_rows)
         for (; k < M2; k++, y += t, s += rU)
            ssfft_ifft(y, t_M2, m1, z_rows, g_rows, e, s, rU_M2, n, scratch);
      
      if (e)
      {
         x += (g_rows << m2) * t;
         for (k = 1, y = x + t; k < (z_rows ? M2 : z_cols); k++, y += t)
            basic_add(*x, *x, *y, n);
         basic_unrotate_bits(*x, *x, m2, n);
      }
   }
}


// end of file ****************************************************************
