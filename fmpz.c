/****************************************************************************

   fmpz.c: "flat" integer format

   Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "fmpz.h"


void mpz_to_fmpz(fmpz_t res, mpz_t x)
{
   if (mpz_sgn(x))
   {
      size_t countp;
      mpz_export(res + 1, &countp, -1, sizeof(mp_limb_t), 0, 0, x);
      res[0] = (mpz_sgn(x) > 0) ? countp : -countp;
   }
   else
      res[0] = 0;
}


void fmpz_to_mpz(mpz_t res, fmpz_t x)
{
   long size = x[0];
   
   if (size == 0)
      mpz_set_ui(res, 0);
   else
   {
      mpz_import(res, ABS(size), -1, sizeof(mp_limb_t), 0, 0, x + 1);

      if (size < 0)
         mpz_neg(res, res);
   }
}


// *************** end of file
