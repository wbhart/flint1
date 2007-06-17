/****************************************************************************

   fmpz.h: "flat" integer format

   Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_FMPZ_H
#define FLINT_FMPZ_H

#include <gmp.h>


typedef mp_limb_t* fmpz_t;


#define ABS(x) (((long) x < 0) ? -x : x)


static inline
unsigned long fmpz_size(fmpz_t x)
{
   long limb = (long) x[0];
   return (unsigned long)  ((limb < 0) ? -limb : limb);
}


// returns positive, negative or zero according to sign of x
static inline
long fmpz_sgn(fmpz_t x)
{
   return (long) x[0];
}


// res must have enough space for x
void mpz_to_fmpz(fmpz_t res, mpz_t x);


void fmpz_to_mpz(mpz_t res, fmpz_t x);



#endif

// *************** end of file
