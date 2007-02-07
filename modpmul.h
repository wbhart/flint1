/******************************************************************************

 DEVELOPMENT CODE

 modpmul.h

 Copyright (C) 2006, David Harvey

 Throughout the documentation in this file, R denotes 2^FLINT_BITS_PER_LIMB.

******************************************************************************/

#ifndef FLINT_MODPMUL_H
#define FLINT_MODPMUL_H

#include "longlong_wrapper.h"


// for debugging:
void print_array(unsigned long* data, unsigned long length);


/*
Returns x * y mod p.

This uses a SLOW division instruction, doesn't use any helpful precomputation.

todo: perhaps this belongs in Z.h, with a better name?
*/
static inline
unsigned long mod_mul(unsigned long x, unsigned long y, unsigned long p)
{
   unsigned long hi, lo, result, dummy;
   
   umul_ppmm(hi, lo, x, y);
   udiv_qrnnd(dummy, result, hi, lo, p);

   return result;
}



/*
In the matrix transpose operations, we switch from recursive subdivision to
a basecase loop when we get down to MODP_FFT_TRANSPOSE_THRESHOLD. The value
is a tradeoff between a few things:

* Shouldn't be too small or we incur too much loop overhead.
* Need that THRESHOLD * THRESHOLD coefficients fit safely in L1.
* Preferably want THRESHOLD >= cache line length (in words).

todo: select this threshold at tuning time
 */
#define FLINT_MODP_FFT_TRANSPOSE_THRESHOLD 32


/*
We switch from "cache-friendly" FFT algorithm to basecase iterative algorithm
when we get down to length 2^MODP_FFT_BASECASE_THRESHOLD. The best value
is around log2(L1 cache size in words).

 (One would expect a value one less to be optimal, because the precomputed
 roots of unity need to go in cache too. But in practice it seems better to
 choose exactly the value above.)

todo: get best value determined at tuning time
*/
#define FLINT_MODP_FFT_BASECASE_THRESHOLD 13



/******************************************************************************

 REDC stuff

******************************************************************************/

// number of bits in the modulus in REDC arithmetic
#define FLINT_REDC_BITS (FLINT_BITS_PER_LIMB - 1)
 
/*
Precomputed information for doing REDC operations, for an odd (not necessarily
prime) modulus. Modulus must fit into FLINT_REDC_BITS bits.
*/
typedef struct redc_precomp_t
{
   unsigned long p;       // the modulus
   unsigned long pinv;    // -p^(-1) mod R
   unsigned long R;       // R mod p
   unsigned long RR;      // R^2 mod p
} redc_precomp_t;

void redc_precomp_init(redc_precomp_t* redc_info, unsigned long p);


/*
Applies Montgomery's REDC reduction to x = R*hi + lo.
It returns x/R mod p, in the range [0, p).
The input x must be in the range [0, Rp).

p and pinv should be as in redc_precomp_t.
 */
static inline
unsigned long reduce_redc(unsigned long hi, unsigned long lo,
                          unsigned long p, unsigned long pinv)
{
   unsigned long d, e;
   umul_ppmm(d, e, lo * pinv, p);
   add_ssaaaa(d, e, d, e, hi, lo);
   if (d >= p)
      d -= p;
   return d;
}


/*
Convenience wrapper for reduce_redc.
 */
static inline
unsigned long reduce_redc2(unsigned long hi, unsigned long lo,
                           redc_precomp_t* info)
{
   return reduce_redc(hi, lo, info->p, info->pinv);
}

/*
Modular multiplication via REDC.

Computes (x*y)/R mod p. The product x*y must be in [0, Rp).

p and pinv should be as in redc_precomp_t.
*/
static inline
unsigned long mul_redc(unsigned long x, unsigned long y,
                       unsigned long p, unsigned long pinv)
{
   unsigned long hi, lo;
   umul_ppmm(hi, lo, x, y);
   return reduce_redc(hi, lo, p, pinv);
}

/*
Convenience wrapper for mul_redc.
 */
static inline
unsigned long mul_redc2(unsigned long x, unsigned long y, redc_precomp_t* info)
{
   return mul_redc(x, y, info->p, info->pinv);
}

/*
Converts input x into REDC format, i.e. returns Rx mod p.

This is achieved by multiplying by R^2 and then using REDC reduction
(So effectively this is an alias for mul_redc :-))

RR, p, pinv should be as in redc_precomp_t.
 */
static inline
unsigned long convert_to_redc(unsigned long x, unsigned long RR,
                              unsigned long p, unsigned long pinv)
{
   return mul_redc(x, RR, p, pinv);
}

/*
Convenience wrapper for convert_to_redc.
*/
static inline
unsigned long convert_to_redc2(unsigned long x, redc_precomp_t* info)
{
   return convert_to_redc(x, info->RR, info->p, info->pinv);
}

/*
Convert input x from REDC format, i.e. returns x/R mod p.

This is achieved by simply calling REDC reduction on x.

p, pinv should be as in redc_precomp_t.
*/
static inline
unsigned long convert_from_redc(unsigned long x, unsigned long p,
                                unsigned long pinv)
{
   return reduce_redc(0, x, p, pinv);
}

/*
Convenience wrapper for convert_from_redc.
 */
static inline
unsigned long convert_from_redc2(unsigned long x, redc_precomp_t* info)
{
   return convert_from_redc(x, info->p, info->pinv);
}



/******************************************************************************

 FFT stuff

******************************************************************************/

/*
Struct for holding precomputed information about an FFT, including precomputed
roots of unity, and information about how to split up the transform to
achieve cache friendliness.
 */
typedef struct modp_fft_precomp_t
{
   // Associated precomputed REDC info
   redc_precomp_t* redc_info;
   
   // N = 2^logN = length of transform
   unsigned long logN;
   unsigned long N;
   
   // primitive N-th root of unity (in REDC format)
   unsigned long w;
   
   // flag indicating whether to use basecase for this transform
   int basecase;

   // If basecase == 1, roots has length N/2, and holds the roots of unity
   // 1, w, w^2, ..., w^(N/2 - 1), in REDC format.

   // If basecase == 0, roots has length 2*logJ, where J is the length of
   // transforms corresponding to child1 (i.e. J = 2^(child1->logN)).
   // It contains w, w^2, w^4, ..., w^(J/2),
   // followed by w^-1, w^-2, w^-4, ..., w^(-J/2), all in REDC format.
   unsigned long* roots;

   // These are precomputed information for subtransforms, only valid
   // if basecase == 0. Must satisfy child1->logN + child2->logN == this->logN.
   struct modp_fft_precomp_t* child1;
   struct modp_fft_precomp_t* child2;

} modp_fft_precomp_t;


void modp_fft_precomp_init(modp_fft_precomp_t* info, redc_precomp_t* redc_info,
                           unsigned long logN, unsigned long w,
                           unsigned long basecase_threshold);

void modp_fft_precomp_clear(modp_fft_precomp_t* info);


void modp_fft(modp_fft_precomp_t* info, unsigned long* data, int first);
void modp_ifft(modp_fft_precomp_t* info, unsigned long* data);


void modp_convolution(redc_precomp_t* redc_info, unsigned long n,
                      unsigned long w,
                      unsigned long* data1, unsigned long* data2);

#endif // #ifndef FLINT_MODPMUL_H
