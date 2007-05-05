/******************************************************************************

 ZmodFermat.h

 Copyright (C) 2007
 David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
 These are currently used only in the ZmodFermatVec module, which supplies
 the Schoenhage-Strassen FFT code.
 
******************************************************************************/


/*
A ZmodFermat_t is stored as a *signed* value in two's complement format, using
n+1 limbs. The value is not normalised into any particular range, so the top
limb may pick up overflow bits. Of course the arithmetic functions in this
module may implicitly reduce mod p whenever they like.

More precisely, suppose that the first n limbs are x[0], ..., x[n-1] (unsigned)
and the last limb is x[n] (signed). Then the value being represented is
 x[0] + x[1]*B + ... + x[n-1]*B^(n-1) - x[n]   (mod p).

*/
typedef mp_limb_t* ZmodFermat_t;


/*
Add the given *signed* limb to the buffer [x, x+count), much like
mpn_add_1 and mpn_sub_1 (except it's always inplace).

PRECONDITIONS:
   count >= 1
   
NOTE:
   The branch predictability of this function is optimised for the case that
   abs(limb) is relatively small and that the first limb of x is randomly
   distributed, which should be the normal usage in the FFT routines.
   
todo: perhaps this should go in mpn_extras instead of here. And perhaps it
      shouldn't start with "mpn", that's not very nice.
*/
static inline
void mpn_signed_add_1(mp_limb_t* x, unsigned long count, mp_limb_signed_t limb)
{
}


/*
Normalises x into the range [0, p).
(Note that the top limb will be set if and only if x = -1 mod p.)
*/
static inline
void ZmodFermat_normalise(ZmodFermat_t x, unsigned long n)
{
}


/*
Adjusts x mod p so that the top limb is in the interval [-1, 1].

This in general will be faster then ZmodFermat_normalise(); in particular
the branching is much more predictable.
*/
static inline
void ZmodFermat_fast_reduce(ZmodFermat_t x, unsigned long n)
{
}


/*
   b := 2^(-s) a   mod p

PRECONDITIONS:
   0 < s < logB
   b may alias a

*/
static inline
void ZmodFermat_short_div_2exp(ZmodFermat_t b, ZmodFermat_t a,
                               unsigned long s, unsigned long n)
{
}



/*
   b := B^s a   mod p

PRECONDITIONS:
   0 < s < n
   b must not alias a

*/
static inline
void ZmodFermat_mul_Bexp(ZmodFermat_t b, ZmodFermat_t a,
                         unsigned long s, unsigned long n)
{
}


/*
   a := a + b
   b := a - b
   z := destroyed

PRECONDITIONS:
   a, b, z may not alias each other
   
NOTE: a, b, z may get permuted

*/
void ZmodFermat_simple_butterfly(ZmodFermat_t* a, ZmodFermat_t* b,
                                 ZmodFermat_t* z, unsigned long n)
{
}


// end of file ****************************************************************
