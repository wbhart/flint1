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


// res := x
// if x == 0, then res needs room only for the control limb
// if x != 0, res needs room for one limb beyond control limb
static inline
void fmpz_set_ui(fmpz_t res, unsigned long x)
{
   if (x) 
   {
      res[0] = 1UL;
      res[1] = x;
   }
   else
      res[0] = 0UL;
}


// same as fmpz_set_ui
static inline
void fmpz_set_si(fmpz_t res, long x)
{
   if (x > 0)
   {
      res[0] = 1L;
      res[1] = x;
   }
   else if (x < 0)
   {
      res[0] = -1L;
      res[1] = -x;
   }
   else
      res[0] = 0UL;
}


// res must have enough space for x
void mpz_to_fmpz(fmpz_t res, mpz_t x);

void fmpz_to_mpz(mpz_t res, fmpz_t x);



#endif

// *************** end of file
