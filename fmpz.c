/****************************************************************************

   fmpz.c: "flat" integer format

   Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "fmpz.h"
#include "flint.h"
#include "memory-manager.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "extras.h"
#include "Z_mpn.h"
#include "Z_mpn_mul-tuning.h"

void mpz_to_fmpz(fmpz_t res, const mpz_t x)
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

void fmpz_to_mpz(mpz_t res, const fmpz_t x)
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

/*
    Adds two fmpz's together
*/

void fmpz_add(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2)
{
   fmpz_t coeffs1 = in1;
   fmpz_t coeffs2 = in2;
   
   long carry;
   unsigned long size1 = ABS(coeffs1[0]);
   unsigned long size2 = ABS(coeffs2[0]);
   
   if (size1 < size2) 
   {
      SWAP_PTRS(coeffs1, coeffs2);
      size1 = ABS(coeffs1[0]);
      size2 = ABS(coeffs2[0]);
   } 
   
   if (!size1)
   {
      if (!size2) coeffs_out[0] = 0L;
      else
      {
         if (coeffs_out != coeffs2) copy_limbs(coeffs_out, coeffs2, size2+1);
      }
   } else if (!size2)
   {
      if (coeffs_out != coeffs1) copy_limbs(coeffs_out, coeffs1, size1+1);
   } else if ((long) (coeffs1[0] ^ coeffs2[0]) >= 0L)
   {
      coeffs_out[0] = coeffs1[0];
      carry = mpn_add(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
      if (carry) 
      {
         coeffs_out[size1+1] = carry;
         if ((long) coeffs_out[0] < 0L) coeffs_out[0]--;
         else coeffs_out[0]++;
      }
   } else
   {
      carry = 0;
      if (size1 != size2) carry = 1;
      else carry = mpn_cmp(coeffs1+1, coeffs2+1, size1); 
          
      if (carry == 0) coeffs_out[0] = 0L;
      else if (carry > 0) 
      {
         mpn_sub(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
         coeffs_out[0] = coeffs1[0];
         NORM(coeffs_out);
      }
      else
      {
         mpn_sub_n(coeffs_out+1, coeffs2+1, coeffs1+1, size1);
         coeffs_out[0] = -coeffs1[0];
         NORM(coeffs_out);
      }
   }
}

/* 
     Add an unsigned long to an fmpz, inplace
*/

void fmpz_add_ui_inplace(fmpz_t output, const unsigned long x)
{
   unsigned long carry;
   
   if (x)
   {
      if (!output[0])
      {
         output[1] = x;
         output[0] = 1;
      } else if ((long) output[0] > 0)
      {
         carry = mpn_add_1(output + 1, output + 1, output[0], x); 
         if (carry)
         {
            output[output[0]+1] = carry;
            output[0]++;
         }
      } else if ((long) output[0] < -1L)
      {
         mpn_sub_1(output + 1, output + 1, ABS(output[0]), x); 
         NORM(output);
      } else
      {
         if (x <= output[1]) 
         {
            output[1] -= x;
            if (!output[1]) output[0] = 0;
         } else
         {
            output[1] = x - output[1];
            output[0] = 1;
         } 
      }
   }
}

/*
   Add an unsigned long to a coefficient. 
   Assumes the output coefficient is non-negative.   
*/

void __fmpz_add_ui_inplace(fmpz_t output, const unsigned long x)
{
   unsigned long carry;
   
   if (x)
   {
      if (!output[0])
      {
         output[1] = x;
         output[0] = 1;
      } else 
      {
         carry = mpn_add_1(output + 1, output + 1, output[0], x); 
         if (carry)
         {
            output[output[0]+1] = carry;
            output[0]++;
         }
      } 
   }
}

void fmpz_sub(fmpz_t coeffs_out, const fmpz_t in1, const fmpz_t in2)
{
   fmpz_t coeffs1 = in1;
   fmpz_t coeffs2 = in2;
   
   long carry;
   unsigned long size1 = ABS(coeffs1[0]);
   unsigned long size2 = ABS(coeffs2[0]);
   int in_order = 1;
   
   if (size1 < size2) 
   {
      SWAP_PTRS(coeffs1, coeffs2);
      size1 = ABS(coeffs1[0]);
      size2 = ABS(coeffs2[0]);
      in_order = 0;
   } 
   
   if (!size1)
   {
      if (!size2) coeffs_out[0] = 0L;
      else
      {
         if (coeffs2 != coeffs_out) copy_limbs(coeffs_out, coeffs2, size2+1);
         if (in_order) coeffs_out[0] = -coeffs_out[0];
      }
   } else if (!size2)
   {
      if (coeffs1 != coeffs_out) copy_limbs(coeffs_out, coeffs1, size1+1);
      if (!in_order) coeffs_out[0] = -coeffs_out[0];
   } else if ((long) (coeffs1[0] ^ coeffs2[0]) < 0)
   {
      if (in_order) coeffs_out[0] = coeffs1[0];
      else coeffs_out[0] = -coeffs1[0];
      carry = mpn_add(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
      if (carry) 
      {
         coeffs_out[size1+1] = carry;
         if ((long) coeffs_out[0] < 0) coeffs_out[0]--;
         else coeffs_out[0]++;
      }
   } else
   {
      carry = 0;
      if (size1 != size2) carry = 1;
      else carry = mpn_cmp(coeffs1+1, coeffs2+1, size1); 
          
      if (carry == 0) coeffs_out[0] = 0L;
      else if (carry > 0) 
      {
         mpn_sub(coeffs_out+1, coeffs1+1, size1, coeffs2+1, size2);
         if (in_order) coeffs_out[0] = coeffs1[0];
         else coeffs_out[0] = -coeffs1[0];
         NORM(coeffs_out);
      }
      else
      {
         mpn_sub_n(coeffs_out+1, coeffs2+1, coeffs1+1, size1);
         if (in_order) coeffs_out[0] = -coeffs1[0];
         else coeffs_out[0] = coeffs1[0];
         NORM(coeffs_out);
      }
   }
}

void fmpz_sub_ui_inplace(fmpz_t output, const unsigned long x)
{
   unsigned long carry;
   
   if (x)
   {
      if (!output[0])
      {
         output[1] = x;
         output[0] = -1L;
      } else if ((long) output[0] < 0)
      {
         carry = mpn_add_1(output + 1, output + 1, ABS(output[0]), x); 
         if (carry)
         {
            output[ABS(output[0])+1] = carry;
            output[0]--;
         }
      } else if ((long) output[0] > 1L)
      {
         mpn_sub_1(output + 1, output + 1, output[0], x); 
         NORM(output);
      } else
      {
         if (x <= output[1]) 
         {
            output[1] -= x;
            if (!output[1]) output[0] = 0;
         } else
         {
            output[1] = x - output[1];
            output[0] = -1L;
         } 
      }
   }
}

/* 
   Multiplies two fmpz's
   Assumes no overlap
*/

void fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b) 
{
      long a0 = a[0];
      long b0 = b[0];
      unsigned long sizea = FLINT_ABS(a0);
      unsigned long sizeb = FLINT_ABS(b0);
      while ((!a[sizea]) && (sizea)) sizea--;
      while ((!b[sizeb]) && (sizeb)) sizeb--;
      
      mp_limb_t mslimb;
      fmpz_t temp;
      
      if ((sizea == 0) || (sizeb == 0))
      {
        res[0] = 0;
      } else if (sizea + sizeb < 100)
      {
         temp = (fmpz_t) flint_stack_alloc_small(sizea + sizeb + 1);
         if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
         else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
         temp[0] = sizea + sizeb - (mslimb == 0);
         copy_limbs(res, temp, temp[0]+1);
         if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
         flint_stack_release_small();     
      } else if (sizea + sizeb < 2*FLINT_FFT_LIMBS_CROSSOVER)
      {
         temp = (fmpz_t) flint_stack_alloc(sizea + sizeb + 1);
         if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
         else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
         temp[0] = sizea + sizeb - (mslimb == 0);
         copy_limbs(res, temp, temp[0]+1);
         if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
         flint_stack_release();   
      } else
      {
         if (sizea >= sizeb) mslimb = Z_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = Z_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea+sizeb - (mslimb == 0);
         if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
      }
}      

/* 
   Used internally by fmpz_poly.
   Multiplies two fmpz's but assumes res has enough space to contain the 
   number of limbs in _a_ plus the number of limbs of _b_, whenever 
   this sum is less than 2*FLINT_FFT_LIMBS_CROSSOVER  
   Assumes no overlap
*/

void __fmpz_mul(fmpz_t res, const fmpz_t a, const fmpz_t b) 
{
      long a0 = a[0];
      long b0 = b[0];
      unsigned long sizea = FLINT_ABS(a0);
      unsigned long sizeb = FLINT_ABS(b0);
      while ((!a[sizea]) && (sizea)) sizea--;
      while ((!b[sizeb]) && (sizeb)) sizeb--;
      
      mp_limb_t mslimb;
      fmpz_t temp;
      
      if ((sizea == 0) || (sizeb == 0))
      {
        res[0] = 0;
      } 
      else if (sizea + sizeb < 100)
      {
         if (sizea >= sizeb) mslimb = mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea + sizeb - (mslimb == 0);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
      } else
      {
         if (sizea >= sizeb) mslimb = Z_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = Z_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea+sizeb - (mslimb == 0);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
      }
}      


void fmpz_mul_ui(fmpz_t output, const fmpz_t input, const unsigned long x)
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

/* 
   Sets res to res+a*b
   Assumes no overlap
*/

void fmpz_addmul(fmpz_t res, const fmpz_t a, const fmpz_t b) 
{
      long a0 = a[0];
      long b0 = b[0];
      unsigned long sizea = FLINT_ABS(a0);
      unsigned long sizeb = FLINT_ABS(b0);
      while ((!a[sizea]) && (sizea)) sizea--;
      while ((!b[sizeb]) && (sizeb)) sizeb--;
      
      fmpz_t temp;
      mp_limb_t mslimb;
      
      if (sizea && sizeb)
      {
         if (sizea + sizeb < 100)
         {
            temp = (fmpz_t) flint_stack_alloc_small(sizea + sizeb + 1);
            if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
            else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
            temp[0] = sizea + sizeb - (mslimb == 0);
            if ((long) (a[0] ^ b[0]) < 0) temp[0] = -temp[0];
            fmpz_add(res, res, temp);
            flint_stack_release_small();
         } else
         {
            temp = (fmpz_t) flint_stack_alloc(sizea + sizeb + 1);
            if (sizea >= sizeb) mslimb = Z_mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
            else mslimb = Z_mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
            temp[0] = sizea + sizeb - (mslimb == 0);
            if ((long) (a[0] ^ b[0]) < 0) temp[0] = -temp[0];
            fmpz_add(res, res, temp);
            flint_stack_release();
         }        
      }     
} 

/* 
   Sets res to a / b
   Assumes no overlap
*/

void fmpz_div(fmpz_t res, const fmpz_t a, const fmpz_t b) 
{
   long a0 = a[0];
   long b0 = b[0];
   unsigned long sizea = FLINT_ABS(a0);
   unsigned long sizeb = FLINT_ABS(b0);
   while ((!a[sizea]) && (sizea)) sizea--;
   while ((!b[sizeb]) && (sizeb)) sizeb--;
      

   mp_limb_t mslimb;
   fmpz_t temp;
      
   if (sizeb == 0)
   {
      printf("Error: division by zero!\n");
      abort();          
   } else if (sizea < sizeb) // Todo: make this deal with sizea == sizeb but a < b
   {
      res[0] = 0;
   } else 
   {
      temp = (fmpz_t) flint_stack_alloc(sizeb);
      mpn_tdiv_qr(res+1, temp, 0, a+1, sizea, b+1, sizeb);
      res[0] = sizea - sizeb + 1;
      if ((long) (a0 ^ b0) < 0) res[0] = -res[0];
      flint_stack_release();  
   }
}     

void fmpz_div_ui(fmpz_t output, const fmpz_t input, const unsigned long x)
{
   output[0] = input[0];
   unsigned long size = FLINT_ABS(input[0]);
   unsigned long xnorm;
     
   if (size > 2) 
   {
      unsigned long norm;
      mp_limb_t xinv;
      
      count_lead_zeros(norm, x);
      xnorm = (x<<norm);
      invert_limb(xinv, xnorm);
      
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
*/

void fmpz_pow_ui(fmpz_t output, const fmpz_t input, const unsigned long exp)
{
   mpz_t power;
   mpz_init(power);
   fmpz_to_mpz(power, input);
   mpz_pow_ui(power, power, exp);
   mpz_to_fmpz(output, power);
   mpz_clear(power);
}

// *************** end of file
