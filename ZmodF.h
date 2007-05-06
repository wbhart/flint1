/******************************************************************************

 ZmodF.h

 Copyright (C) 2007, David Harvey
 
 Routines for arithmetic on elements of Z/pZ where p = B^n + 1,
 B = 2^FLINT_BITS_PER_LIMB.
 
 These are currently used only in the ZmodFpoly module, which supplies the
 Schoenhage-Strassen FFT code.
 
******************************************************************************/


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


/*
   x := y
*/
static inline
void ZmodF_set(ZmodF_t x, ZmodF_t y, unsigned long n)
{
   abort();
}


/*
   x := -y
   
   PRECONDITIONS:
      x and y may alias each other
*/
static inline
void ZmodF_neg(ZmodF_t x, ZmodF_t y, unsigned long n)
{
   abort();
}


/*
   res := x*y
   
   PRECONDITIONS:
      Any combination of aliasing among res, x, y is allowed.
      scratch must be a buffer of length 2*n, and must NOT alias x, y, res.
*/
static inline
void ZmodF_mul(ZmodF_t res, ZmodF_t x, ZmodF_t y, mp_limb_t* scratch,
               unsigned long n)
{
   abort();
}


/*
   Normalises x into the range [0, p).
   (Note that the top limb will be set if and only if x = -1 mod p.)
*/
static inline
void ZmodF_normalise(ZmodF_t x, unsigned long n)
{
   abort();
}


/*
   Adjusts x mod p so that the top limb is in the interval [-1, 1].

   This in general will be faster then ZmodF_normalise(); in particular
   the branching is much more predictable.
*/
static inline
void ZmodF_fast_reduce(ZmodF_t x, unsigned long n)
{
   abort();
}


/*
   b := 2^(-s) a

   PRECONDITIONS:
      0 < s < logB
      b may alias a

*/
static inline
void ZmodF_short_div_2exp(ZmodF_t b, ZmodF_t a,
                          unsigned long s, unsigned long n)
{
   abort();
}


/*
   b := 2^(s/2) a
   z := destroyed

   PRECONDITIONS:
      0 < s < logB
      a may alias z
      b may not alias a or z

   NOTE: a, b, z may get permuted
*/
static inline
void ZmodF_mul_sqrt2exp(ZmodF_t* b, ZmodF_t* a, ZmodF_t* z,
                        unsigned long s, unsigned long n)
{
   abort();
}


/*
   a := a + b
   b := 2^(s/2) (a - b)
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
static inline
void ZmodF_forward_butterfly_sqrt2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                      unsigned long s, unsigned long n)
{
   abort();
}


/*
   a := a + 2^(-s/2) b
   b := a - 2^(-s/2) b
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
static inline
void ZmodF_inverse_butterfly_sqrt2exp(ZmodF_t* a, ZmodF_t* b, ZmodF_t* z,
                                      unsigned long s, unsigned long n)
{
   abort();
}


/*
   b := B^s a

   PRECONDITIONS:
      0 < s < n
      b must not alias a
*/
static inline
void ZmodF_mul_Bexp(ZmodF_t b, ZmodF_t a,
                    unsigned long s, unsigned long n)
{
   abort();
}


/*
   a := a + b
   b := a - b
   z := destroyed

   PRECONDITIONS:
      a, b, z may not alias each other
      
   NOTE: a, b, z may get permuted
*/
void ZmodF_simple_butterfly(ZmodF_t* a, ZmodF_t* b,
                            ZmodF_t* z, unsigned long n)
{
   abort();
}


// end of file ****************************************************************
