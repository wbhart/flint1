/****************************************************************************

   fmpz.c: "flat" integer format

   Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>

#include "fmpz.h"
#include "flint.h"
#include "memory-manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "extras.h"

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

void _fmpz_mul_ui(fmpz_t output, fmpz_t input, unsigned long x)
{
   if (x == 0)
   {
      output[0] = 0;
      return;
   }
   
   mp_limb_t mslimb;
   
   if (output[0] = input[0]) // This isn't a typo
   {
      mslimb = mpn_mul_1(output+1, input+1, FLINT_ABS(input[0]), x);
      if (mslimb) 
      {
         output[FLINT_ABS(input[0])+1] = mslimb; 
         if ((long) output[0] > 0) output[0]++;
         else output[0]--;  
      }
   }
}

void _fmpz_div_ui(fmpz_t output, fmpz_t input, unsigned long x)
{
   output[0] = input[0];
   unsigned long size = FLINT_ABS(input[0]);
      
   if (size > 2) 
   {
      unsigned long norm;
      mp_limb_t xinv;
      
      count_lead_zeros(norm, x);
      x <<= norm;
      invert_limb(xinv, x);
      x >>= norm;
      
      mpn_divmod_1_preinv(output+1, input+1, size, x, xinv, norm);
   } else
   {
      mpn_divmod_1(output+1, input+1, size, x);
   }
   
   NORM(output);
}

/*
   Raise input to the power exp
   Very simplistic at this point. It just converts to an mpz_t and uses
   GMP's mpz_pow_ui function
   This function allocates space for output
*/

void fmpz_pow_ui(fmpz_t * output, fmpz_t input, unsigned long exp)
{
   mpz_t power;
   mpz_init(power);
   fmpz_to_mpz(power, input);
   mpz_pow_ui(power, power, exp);
   *output = (fmpz_t) flint_heap_alloc(mpz_size(power) + 1);
   mpz_to_fmpz(*output, power);
   mpz_clear(power);
}

/*
   Raise input to the power exp
   Very simplistic at this point. It just converts to an mpz_t and uses
   GMP's mpz_pow_ui function
*/

void _fmpz_pow_ui(fmpz_t output, fmpz_t input, unsigned long exp)
{
   mpz_t power;
   mpz_init(power);
   fmpz_to_mpz(power, input);
   mpz_pow_ui(power, power, exp);
   mpz_to_fmpz(output, power);
   mpz_clear(power);
}

// *************** end of file
