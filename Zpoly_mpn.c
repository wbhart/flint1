/****************************************************************************

Zpoly_mpn.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"
#include "mpn_extras.h"
#include "extras.h"
#include "longlong_wrapper.h"
#include "memory-manager.h"
#include "ZmodFpoly.h"
#include "Z_mpn.h"
#include "ZmodF_mul.h"

/****************************************************************************

   _Zpoly_mpn_* layer

****************************************************************************/

/* 
   Create a polynomial of length zero with "alloc" allocated coefficients
   each with enough space for limbs limbs
*/

void _Zpoly_mpn_stack_init(Zpoly_mpn_t poly, unsigned long alloc, unsigned long limbs)
{
   FLINT_ASSERT(alloc >= 1);
   FLINT_ASSERT(limbs >= 1);

   poly->coeffs = (mp_limb_t *) flint_stack_alloc(alloc*(limbs+1));
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

void _Zpoly_mpn_stack_clear(Zpoly_mpn_t poly)
{
   flint_stack_release();
}

static inline
void _convert_coeff_to_mpz(mpz_t x, mp_limb_t* coeff)
{
   long size = coeff[0];
   
   if (size == 0) 
      mpz_set_ui(x, 0);
   else
   {
      mpz_import(x, ABS(size), -1, sizeof(mp_limb_t), 0, 0, coeff + 1);

      if (size < 0)
         mpz_neg(x, x);
   }
}


// retrieves coefficient #n as an mpz, no bounds checking
void _Zpoly_mpn_get_coeff_mpz(mpz_t x, Zpoly_mpn_t poly, unsigned long n)
{
   FLINT_ASSERT(n < poly->length);
   _convert_coeff_to_mpz(x, poly->coeffs + n*(poly->limbs+1));
}


// converts Zpoly_mpn to Zpoly, assumes enough space allocated for output
void _Zpoly_mpn_convert_out(Zpoly_t poly_mpz, Zpoly_mpn_t poly_mpn)
{
   FLINT_ASSERT(poly_mpz->alloc >= poly_mpn->length);

   poly_mpz->length = poly_mpn->length;
   
   for (unsigned long i = 0; i < poly_mpn->length; i++)
      _convert_coeff_to_mpz(poly_mpz->coeffs[i],
                            poly_mpn->coeffs + i*(poly_mpn->limbs+1));
}


void _Zpoly_mpn_convert_in(Zpoly_mpn_t poly_mpn, Zpoly_t poly_mpz)
{
   FLINT_ASSERT(poly_mpn->alloc >= poly_mpz->length);

   poly_mpn->length = poly_mpz->length;
   if (poly_mpz->length == 0) return;
   for (unsigned long i = 0; i < poly_mpz->length; i++)
   {
      if (mpz_sgn(poly_mpz->coeffs[i]) == 0) 
          poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = 0L;
      else
      {
          size_t countp;
          
          FLINT_ASSERT(poly_mpn->limbs >= mpz_size(poly_mpz->coeffs[i]));
          
          mpz_export(poly_mpn->coeffs + i*(poly_mpn->limbs+1) + 1, &countp, 
                     -1, sizeof(mp_limb_t), 0, 0, poly_mpz->coeffs[i]);
                     
          if (mpz_sgn(poly_mpz->coeffs[i]) < 0) 
          {
             poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = -countp;
          } else
             poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = countp;
      }
   }
}

/* 
   Set a coefficient to the given unsigned value.
   Clears dirty limbs unless the coefficient is set to zero. 
   Sets the sign to 1 if x is positive, else to zero.
   Assumes the polynomial length is greater than n.
*/
void _Zpoly_mpn_set_coeff_ui(Zpoly_mpn_t poly, unsigned long n, unsigned long x)
{
   FLINT_ASSERT(poly->alloc > n);
   if (x == 0) 
   {
      poly->coeffs[n*(poly->limbs+1)] = 0UL;
      return;
   }
   poly->coeffs[n*(poly->limbs+1)] = 1UL;
   poly->coeffs[n*(poly->limbs+1)+1] = x;
}

/* 
   Set a coefficient to the given signed value.
   Clears dirty limbs unless the coefficient is set to zero. 
   Sets the sign to 1 if x is positive, -1 if negative, else to zero.
   Assumes the polynomial length is greater than n.
*/

void _Zpoly_mpn_set_coeff_si(Zpoly_mpn_t poly, unsigned long n, long x)
{
   FLINT_ASSERT(poly->length > n);
   if (x == 0)
   {
      poly->coeffs[n*(poly->limbs+1)] = 0UL;
      return;
   }

   if (x > 0)
   {
      poly->coeffs[n*(poly->limbs+1)] = 1L;
      poly->coeffs[n*(poly->limbs+1)+1] = x;
   }
   else
   {
      poly->coeffs[n*(poly->limbs+1)] = -1L;
      poly->coeffs[n*(poly->limbs+1)+1] = -x;
   } 
}

void _Zpoly_mpn_normalise(Zpoly_mpn_t poly)
{
   while (poly->length && poly->coeffs[(poly->length-1)*(poly->limbs+1)] == 0)
      poly->length--;
}

/* 
   Sets the output poly to equal the input poly 
   Assumes the output poly is big enough to hold the nonzero limbs of the input poly
*/

void _Zpoly_mpn_set(Zpoly_mpn_t output, Zpoly_mpn_t input)
{
   if (output->coeffs != input->coeffs) 
   {
      unsigned long input_size = input->limbs + 1;
      unsigned long output_size = output->limbs + 1;
      if ((output->coeffs < input->coeffs) || (output->coeffs >= input->coeffs + input->length*(input->limbs+1)))
      {
         for (unsigned long i = 0; i < input->length; i++)
         {
            if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
            else copy_limbs(output->coeffs+i*output_size, input->coeffs+i*input_size, ABS(input->coeffs[i*input_size])+1);
         }
      } else
      {
         for (long i = input->length - 1; i >= 0; i--)
         {
            if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
            else copy_limbs(output->coeffs+i*output_size, input->coeffs+i*input_size, ABS(input->coeffs[i*input_size])+1);
         }
      }
   }
   output->length = input->length;
}

void _Zpoly_mpn_swap(Zpoly_mpn_t x, Zpoly_mpn_t y)
{
   mp_limb_t * temp_p;
   mp_limb_t temp_l;
   
   temp_p = x->coeffs;
   x->coeffs = y->coeffs;
   y->coeffs = temp_p;
   
   temp_l = x->alloc;
   x->alloc = y->alloc;
   y->alloc = temp_l;
   
   temp_l = x->length;
   x->length = y->length;
   y->length = temp_l;
   
   temp_l = x->limbs;
   x->limbs = y->limbs;
   y->limbs = temp_l;
}

/*
   Determines the maximum number of bits in a coefficient of poly_mpn. This
   function assumes every coefficient fits in a limb. The returned value is
   negative if any of the coefficients was negative.
*/
long _Zpoly_mpn_bits1(Zpoly_mpn_t poly_mpn)
{
   unsigned long mask = -1L;
   long bits = 0;
   long sign = 1;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   unsigned long i, j;
   
   for (i = 0, j = 0; i < poly_mpn->length; i++, j += 2)
   {
      if (i&3 == 0) FLINT_PREFETCH(coeffs_m+j,64);
      if ((long) coeffs_m[j] < 0) sign = -1L;
      if (coeffs_m[j])
      {
         if (coeffs_m[j+1] & mask)
         {
            bits = FLINT_BIT_COUNT(coeffs_m[j+1]);   
            if (bits == FLINT_BITS_PER_LIMB) break;
            else mask = -1L - ((1L<<bits)-1);
         }
      }
   }
   
   if (sign == 1)
   {
      for ( ; i < poly_mpn->length; i++, j += 2)
      { 
         if ((long) coeffs_m[j] < 0) 
         {
            sign = -1L;
            break;
         }
      }
   }
   
   return sign*bits;
}

/*
   Determines the maximum number of bits in a coefficient of poly_mpn. 
   The returned value is negative if any of the coefficients was negative.
*/
long _Zpoly_mpn_bits(Zpoly_mpn_t poly_mpn)
{
   unsigned long mask = -1L;
   long bits = 0;
   long sign = 1;
   long limbs = 0;
   long size_j;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long i, j;
   
   for (i = 0, j = 0; i < poly_mpn->length; i++, j += size_m)
   {
      size_j = (long) coeffs_m[j];
      if (size_j < 0) sign = -1L;
      if (ABS(size_j) > limbs + 1)
      {
         limbs = ABS(size_j) - 1;
         bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]);   
         if (bits == FLINT_BITS_PER_LIMB) mask = 0L;
         else mask = -1L - ((1L<<bits)-1);
      } else if (ABS(size_j) == limbs+1)
      {
         if (coeffs_m[j+ABS(size_j)] & mask)
         {
            bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]);   
            if (bits == FLINT_BITS_PER_LIMB) mask = 0L;
            else mask = -1L - ((1L<<bits)-1);
         }
      }       
   }
   
   if (sign == 1)
   {
      for ( ; i < poly_mpn->length; i++, j += size_m)
      { 
         if ((long) coeffs_m[j] < 0) 
         {
            sign = -1L;
            break;
         }
      }
   }
   
   return sign*(FLINT_BITS_PER_LIMB*limbs+bits);
}

int _Zpoly_mpn_equal(Zpoly_mpn_p input1, Zpoly_mpn_p input2)
{
   int shorter_poly;
   unsigned long shorter_length;
   unsigned long i,j;
   
   if (input1->length > input2->length)
   {
      SWAP(input1, input2);
   }  
   
   for (unsigned long i = input1->length; i < input2->length; i++)
      if (input2->coeffs[i*(input2->limbs+1)] != 0L) return 0;
      
   for (i = 0; i < input1->length; i++)
   {
      for (j = 0; j < ABS(input1->coeffs[i*(input1->limbs+1)])+1; j++)
      {
         if (input1->coeffs[i*(input1->limbs+1)+j] != input2->coeffs[i*(input2->limbs+1)+j]) 
            return 0; 
      }
   }
      
   return 1;

}

void _Zpoly_mpn_negate(Zpoly_mpn_t output, Zpoly_mpn_t input)
{
   if (input->coeffs == output->coeffs)
   {
      for (unsigned long i = 0; i < input->length; i++)
         output->coeffs[i*(output->limbs+1)] = -output->coeffs[i*(output->limbs+1)];
   } else
   {
      unsigned long input_size = input->limbs + 1;
      unsigned long output_size = output->limbs + 1;
      for (unsigned long i = 0; i < input->length; i++)
      {
         if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
         else 
         {
            output->coeffs[i*output_size] = -input->coeffs[i*input_size];
            copy_limbs(output->coeffs+i*output_size+1, input->coeffs+i*input_size+1, ABS(input->coeffs[i*input_size]));
         }
      }
   }
   output->length = input->length;
}

/* 
   Multiplies input by x^n and sets output to the result
   Assumes output is large enough to contain the result
*/

void _Zpoly_mpn_left_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, 
                                                 unsigned long n)
{
   Zpoly_mpn_t part;   
   
   part->length = input->length;
   part->limbs = output->limbs;
   part->coeffs = output->coeffs + n*(output->limbs+1);
      
   _Zpoly_mpn_set(part, input);
   for (unsigned long i = 0; i < n; i++) output->coeffs[i*(output->limbs+1)] = 0;
   
   output->length = input->length + n;
}

/* 
   Divides input by x^n losing the remainder and sets output to the result
   Assumes output is large enough to contain the result
*/

void _Zpoly_mpn_right_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, unsigned long n)
{
   if (input->length <= n) 
   {
      _Zpoly_mpn_zero(output);
      return;
   }
   Zpoly_mpn_t part;
   part->length = input->length - n;
   part->limbs = input->limbs;
   part->coeffs = input->coeffs + n*(input->limbs + 1);
   _Zpoly_mpn_set(output, part);
}

/*
    Adds two coefficients together
*/

void __Zpoly_mpn_add_coeffs(mp_limb_t * coeffs_out, mp_limb_t * coeffs1, mp_limb_t * coeffs2)
{
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
         copy_limbs(coeffs_out, coeffs2, size2+1);
      }
   } else if (!size2)
   {
      copy_limbs(coeffs_out, coeffs1, size1+1);
   } else if ((long) (coeffs1[0] ^ coeffs2[0]) >= 0)
   {
      coeffs_out[0] = coeffs1[0];
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

void __Zpoly_mpn_add_coeff_ui(mp_limb_t * output, unsigned long x)
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

void __Zpoly_mpn_add_coeff2_ui(mp_limb_t * output, unsigned long x)
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

void __Zpoly_mpn_sub_coeff_ui(mp_limb_t * output, unsigned long x)
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
    Add two polynomials together 
*/

void _Zpoly_mpn_add(Zpoly_mpn_t output, Zpoly_mpn_p input1, Zpoly_mpn_p input2)
{
   unsigned long size1, size2, shorter, size_out;
   mp_limb_t * coeffs1, * coeffs2, * coeffs_out;
   
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   
   size_out = output->limbs+1;
   coeffs_out = output->coeffs;
   
   shorter = (input1->length > input2->length) ? input2->length : input1->length;
   
   for (unsigned long i = 0; i < shorter; i++)
   {
      __Zpoly_mpn_add_coeffs(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2+i*size2);
   }
   
   for (unsigned long i = shorter; i < input1->length; i++)
   {
       copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, ABS(coeffs1[i*size1])+1);
   }
   for (unsigned long i = shorter; i < input2->length; i++)
   {
       copy_limbs(coeffs_out+i*size_out, coeffs2+i*size2, ABS(coeffs2[i*size2])+1);
   }
   
   output->length = (input1->length > input2->length) ? input1->length : input2->length;
}

void __Zpoly_mpn_sub_coeffs(mp_limb_t * coeffs_out, mp_limb_t * coeffs1, mp_limb_t * coeffs2)
{
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
         copy_limbs(coeffs_out, coeffs2, size2+1);
         if (in_order) coeffs_out[0] = -coeffs_out[0];
      }
   } else if (!size2)
   {
      copy_limbs(coeffs_out, coeffs1, size1+1);
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

/* 
    Add two polynomials together 
*/

void _Zpoly_mpn_sub(Zpoly_mpn_t output, Zpoly_mpn_p input1, Zpoly_mpn_p input2)
{
   unsigned long size1, size2, shorter, size_out;
   mp_limb_t * coeffs1, * coeffs2, * coeffs_out;
   
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   
   size_out = output->limbs+1;
   coeffs_out = output->coeffs;
   
   shorter = (input1->length > input2->length) ? input2->length : input1->length;
   
   for (unsigned long i = 0; i < shorter; i++)
   {
      __Zpoly_mpn_sub_coeffs(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2+i*size2);
   }
   
   for (unsigned long i = shorter; i < input1->length; i++)
   {
       copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, ABS(coeffs1[i*size1])+1);
   }
   for (unsigned long i = shorter; i < input2->length; i++)
   {
       copy_limbs(coeffs_out+i*size_out+1, coeffs2+i*size2+1, ABS(coeffs2[i*size2]));
       coeffs_out[i*size_out] = -coeffs2[i*size2];
   }
   
   output->length = (input1->length > input2->length) ? input1->length : input2->length;
}

/* 
   Multiplies two coefficient
   Assumes no overlap
*/

void __Zpoly_mpn_mul_coeffs(mp_limb_t * res, mp_limb_t * a, mp_limb_t * b) 
{
      NORM(a);
      NORM(b);
      unsigned long sizea = ABS(a[0]);
      unsigned long sizeb = ABS(b[0]);
      mp_limb_t mslimb;
      if ((sizea == 0) || (sizeb == 0))
      {
        res[0] = 0;
      } else
      {
         if (sizea >= sizeb) mslimb = Z_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = Z_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea+sizeb - (mslimb == 0);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
      }
}      

void _Zpoly_mpn_scalar_mul_ui(Zpoly_mpn_t output, Zpoly_mpn_t poly, unsigned long x)
{
     if (x == 0) 
     {
        output->length = 0;
        return;
     }
     
     mp_limb_t * coeffs1 = poly->coeffs;
     mp_limb_t * coeffs_out = output->coeffs;
     unsigned long size1 = poly->limbs+1;
     unsigned long size_out = output->limbs+1;
     mp_limb_t mslimb;
     
     for (unsigned long i = 0; i < poly->length; i++)
     {
        if ((coeffs_out[i*size_out] = coeffs1[i*size1]))
        {
           mslimb = mpn_mul_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
           if (mslimb) 
           {
              coeffs_out[i*size_out+ABS(coeffs1[i*size1])+1] = mslimb; 
              if ((long) coeffs_out[i*size_out] > 0) coeffs_out[i*size_out]++;
              else coeffs_out[i*size_out]--;  
           }
        }
     }
     output->length = poly->length;
}

void _Zpoly_mpn_scalar_mul_si(Zpoly_mpn_t output, Zpoly_mpn_t poly, long x)
{
     if (x == 0) 
     {
        output->length = 0;
        return;
     }
     
     mp_limb_t * coeffs1 = poly->coeffs;
     mp_limb_t * coeffs_out = output->coeffs;
     unsigned long size1 = poly->limbs+1;
     unsigned long size_out = output->limbs+1;
     mp_limb_t mslimb;
     
     for (unsigned long i = 0; i < poly->length; i++)
     {
        if (x < 0)
        {
           if ((coeffs_out[i*size_out] = -coeffs1[i*size1]))
           {
              mslimb = mpn_mul_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), -x);
              if (mslimb) 
              {
                 coeffs_out[i*size_out+ABS(coeffs1[i*size1])+1] = mslimb; 
                 if ((long) coeffs_out[i*size_out] > 0) coeffs_out[i*size_out]++;
                 else coeffs_out[i*size_out]--;
              }    
           }
        } else 
        {
           if ((coeffs_out[i*size_out] = coeffs1[i*size1]))
           {
              mslimb = mpn_mul_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
              if (mslimb) 
              {
                 coeffs_out[i*size_out+ABS(coeffs1[i*size1])+1] = mslimb; 
                 if ((long) coeffs_out[i*size_out] > 0) coeffs_out[i*size_out]++;
                 else coeffs_out[i*size_out]--;
              }    
           }
        }
     }
     output->length = poly->length;
}

void _Zpoly_mpn_scalar_div_exact_ui(Zpoly_mpn_t output, Zpoly_mpn_t poly, unsigned long x)
{
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs1 = poly->coeffs;
      
   if (size_out != size1)
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
         coeffs_out[i*size_out] = coeffs1[i*size1];
         NORM(coeffs_out+i*size_out);
      }
   } else
   {
      if (coeffs_out != coeffs1)
      {
         coeffs_out[0] = 0;
         for (unsigned long i = 0; i < poly->length-1; i++)
         {
            copy_limbs(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]));
            clear_limbs(coeffs_out+i*size_out+ABS(coeffs1[i*size1])+1, size_out-ABS(coeffs1[i*size1]));
         } 
         copy_limbs(coeffs_out+(poly->length-1)*size_out+1, coeffs1+(poly->length-1)*size1+1, ABS(coeffs1[(poly->length-1)*size1]));
         if (size_out > ABS(coeffs1[(poly->length-1)*size1])+1) clear_limbs(coeffs_out+(poly->length-1)*size_out+ABS(coeffs1[(poly->length-1)*size1])+1, size_out-ABS(coeffs1[(poly->length-1)*size1])-1);
         
         mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            coeffs_out[i*size_out] = coeffs1[i*size1];
            NORM(coeffs_out+i*size_out);
         }
      } else
      {
         mp_limb_t * signs = (mp_limb_t *) flint_stack_alloc(poly->length);
         signs[0] = coeffs1[0];
         coeffs_out[0] = 0;
         for (unsigned long i = 0; i < poly->length-1; i++)
         {
             signs[i+1] = coeffs1[(i+1)*size1];
             clear_limbs(coeffs_out+i*size_out+ABS(signs[i])+1, size_out-ABS(signs[i]));
         } 
         if (size_out > ABS(signs[poly->length-1])+1) clear_limbs(coeffs_out+(poly->length-1)*size_out+ABS(signs[poly->length-1])+1, size_out-ABS(signs[poly->length-1])-1);
         mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            coeffs_out[i*size_out] = signs[i];
            NORM(coeffs_out+i*size_out);
         }
         flint_stack_release();
      }
   }
   output->length = poly->length;
}

void _Zpoly_mpn_scalar_div_exact_si(Zpoly_mpn_t output, Zpoly_mpn_t poly, long x)
{
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs1 = poly->coeffs;
      
   if (size_out != size1)
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (x < 0)
         {
            mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), -x);
            coeffs_out[i*size_out] = -coeffs1[i*size1];
         } else
         {
            mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
            coeffs_out[i*size_out] = coeffs1[i*size1];
         }
         NORM(coeffs_out+i*size_out);
      }
   } else
   {
      if (coeffs_out != coeffs1)
      {
         coeffs_out[0] = 0;
         for (unsigned long i = 0; i < poly->length-1; i++)
         {
            copy_limbs(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]));
            clear_limbs(coeffs_out+i*size_out+ABS(coeffs1[i*size1])+1, size_out-ABS(coeffs1[i*size1]));
         } 
         copy_limbs(coeffs_out+(poly->length-1)*size_out+1, coeffs1+(poly->length-1)*size1+1, ABS(coeffs1[(poly->length-1)*size1]));
         if (size_out > ABS(coeffs1[(poly->length-1)*size1])+1) clear_limbs(coeffs_out+(poly->length-1)*size_out+ABS(coeffs1[(poly->length-1)*size1])+1, size_out-ABS(coeffs1[(poly->length-1)*size1])-1);
         
         if (x < 0) mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, -x);
         else mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            if (x < 0) coeffs_out[i*size_out] = -coeffs1[i*size1];
            else coeffs_out[i*size_out] = coeffs1[i*size1];
            NORM(coeffs_out+i*size_out);
         }
      } else
      {
         mp_limb_t * signs = (mp_limb_t *) flint_stack_alloc(poly->length);
         signs[0] = coeffs1[0];
         coeffs_out[0] = 0;
         for (unsigned long i = 0; i < poly->length-1; i++)
         {
             signs[i+1] = coeffs1[(i+1)*size1];
             clear_limbs(coeffs_out+i*size_out+ABS(signs[i])+1, size_out-ABS(signs[i]));
         } 
         if (size_out > ABS(signs[poly->length-1])+1) clear_limbs(coeffs_out+(poly->length-1)*size_out+ABS(signs[poly->length-1])+1, size_out-ABS(signs[poly->length-1])-1);
         if (x < 0) mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, -x);
         else mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            if (x < 0) coeffs_out[i*size_out] = -signs[i];
            else coeffs_out[i*size_out] = signs[i];
            NORM(coeffs_out+i*size_out);
         }
         flint_stack_release();
      }
   }
   output->length = poly->length;
}

/* 
    Does scalar division of a polynomial by a limb x. Currently rounding is done towards
    zero.
*/

void _Zpoly_mpn_scalar_div_ui(Zpoly_mpn_t output, Zpoly_mpn_t poly, unsigned long x)
{
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs1 = poly->coeffs;
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      
      count_leading_zeros (norm, x);
      x <<= norm;
      invert_limb(xinv,x);
      x >>= norm;
      
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1_preinv(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x, xinv, norm);
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
      }
   }
   
   output->length = poly->length;
}

void _Zpoly_mpn_scalar_div_si(Zpoly_mpn_t output, Zpoly_mpn_t poly, long x)
{
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs1 = poly->coeffs;
   int sign = (x < 0);
   if (sign) x = -x; 
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      
      count_leading_zeros (norm, (unsigned long) x);
      x = ((unsigned long) x) << norm;
      invert_limb(xinv,(unsigned long) x);
      x = ((unsigned long) x) >> norm;
      
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (sign) coeffs_out[i*size_out] = -coeffs1[i*size1];
         else coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1_preinv(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x, xinv, norm);
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (sign) coeffs_out[i*size_out] = -coeffs1[i*size1];
         else coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
      }
   }
   
   output->length = poly->length;
}

void _Zpoly_mpn_mul_naive(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   mp_limb_t * coeffs_out = output->coeffs;
   unsigned long size_out = output->limbs+1;
   mp_limb_t * coeffs1, * coeffs2;
   unsigned long size1, size2;
   unsigned long len1, len2; 
   unsigned long lenm1;
      
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   lenm1 = input1->length-1;
   len1 = input1->length;
   len2 = input2->length;
      
   mp_limb_t * temp;
   
   if (len1 + len2 > 2) temp = (mp_limb_t *) flint_stack_alloc(size1+size2-1);
         
   for (unsigned long i = 0; i < len1; i++)
   {
      /* Set out[i] = in1[i]*in2[0] */
      if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
      {
         coeffs_out[i*size_out]=0;
      } else
      {
         __Zpoly_mpn_mul_coeffs(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
      }
   }
   for (unsigned long i = 1; i < len2; i++)
   {
      /* Set out[i+in1->length-1] = in1[in1->length-1]*in2[i] */
      if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
      {
         coeffs_out[(i+lenm1)*size_out]=0;
      } else
      {
         __Zpoly_mpn_mul_coeffs(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
      }      
   }
   for (unsigned long i = 0; i < lenm1; i++)
   {      
      for (unsigned long j = 1; j < len2; j++)
      {
         /* out[i+j] += in1[i]*in2[j] */
         if ((coeffs1[i*size1] != 0) && (coeffs2[j*size2] != 0))
         {
            if (!coeffs_out[(i+j)*size_out])
            {
               __Zpoly_mpn_mul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
            } else 
            {
               __Zpoly_mpn_mul_coeffs(temp, coeffs1+i*size1, coeffs2+j*size2);
               __Zpoly_mpn_add_coeffs(coeffs_out+(i+j)*size_out, temp, coeffs_out+(i+j)*size_out);
            } 
         }
      }
   } 
   output->length = len1 + len2 - 1;
   if (len1 + len2 > 2) flint_stack_release();
}

void __Zpoly_mpn_karamul_recursive(Zpoly_mpn_t res, Zpoly_mpn_t a, Zpoly_mpn_t b, Zpoly_mpn_t scratch, Zpoly_mpn_t scratchb)
{
   Zpoly_mpn_t temp;
   
   if ((a->length <= 1) || (b->length <= 1)) 
   {
      _Zpoly_mpn_mul_naive(res, a, b);
      
      return;
   }
   
   if (a->length ==2 && b->length == 2) {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      __Zpoly_mpn_mul_coeffs(res->coeffs, a->coeffs, b->coeffs); 
      __Zpoly_mpn_add_coeffs(scratchb->coeffs, a->coeffs, a->coeffs+asize);
      __Zpoly_mpn_mul_coeffs(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
      __Zpoly_mpn_add_coeffs(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
      __Zpoly_mpn_mul_coeffs(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
      __Zpoly_mpn_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
      __Zpoly_mpn_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
      
      res->length = a->length + b->length - 1;
      
      return;
   }
   
   /*if ((a->length <= 4) || (b->length <= 4)) 
   {
      _Zpoly_mpn_mul_naive(res, a, b);
      
      return;
   }*/


   /* 
      As we may have dirty limbs in our res polynomial (it might be part of the
      original scratch polynomial) we need to clean it before using it 
   */
   
   for (unsigned long i = 0, j = 0; i < a->length + b->length - 1; i++, j+= (scratch->limbs+1))
   {
      res->coeffs[j] = 0;
   }
      
   Zpoly_mpn_t a1,a2,b1,b2;
      
   unsigned long l2 = 0;
      
   a1->length = (a->length+1)/2;
   a2->length = a->length-a1->length;
   a1->coeffs = a->coeffs;
   a2->coeffs = a->coeffs+a1->length*(a->limbs+1);
   a1->limbs = a->limbs;
   a2->limbs = a->limbs;
   
   if (a1->length < b->length) //ordinary case
   {
      /*
         (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
      */
      
      b1->length = a1->length;
      b2->length = b->length - b1->length;
      b1->coeffs = b->coeffs;
      b2->coeffs = b->coeffs + b1->length*(b->limbs+1);
      b1->limbs = b->limbs;
      b2->limbs = b->limbs;
      
      Zpoly_mpn_t asum, bsum, prodsum, scratch2, scratch3;
     
      asum->length = a1->length;
      asum->coeffs = scratchb->coeffs;
      asum->limbs = scratchb->limbs;
      bsum->length = a1->length;
      bsum->coeffs = scratchb->coeffs + a1->length*(scratchb->limbs+1);
      bsum->limbs = scratchb->limbs;
      prodsum->length = (a1->length<<1)-1;
      prodsum->coeffs = scratch->coeffs + (a1->length<<1)*(scratch->limbs+1);
      prodsum->limbs = scratch->limbs;
      
      // res_lo = a1*b1
      scratch2->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs+((a1->length<<2)-1)*(scratch->limbs+1);
      __Zpoly_mpn_karamul_recursive(res, a1, b1, scratch2, scratchb);
     
      // res_hi = a2*b2
      temp->coeffs = res->coeffs+(a1->length<<1)*(res->limbs+1);
      temp->limbs = res->limbs;
      __Zpoly_mpn_karamul_recursive(temp, a2, b2, scratch2, scratchb);
      
      // asum = a1+a2
      _Zpoly_mpn_add(asum, a1, a2);
      // bsum = b1+b2
      _Zpoly_mpn_add(bsum, b1, b2);
      // prodsum = asum*bsum
      scratch3->coeffs = scratchb->coeffs+(a1->length<<1)*(scratchb->limbs+1);
      scratch3->limbs = scratchb->limbs;
      
      __Zpoly_mpn_karamul_recursive(prodsum, asum, bsum, scratch2, scratch3);
      
      // prodsum = prodsum - res_lo
      temp->coeffs = res->coeffs;
      temp->length = (a1->length<<1)-1;
      _Zpoly_mpn_sub(prodsum, prodsum, temp);
       
      // prodsum = prodsum - res_hi
      temp->coeffs = res->coeffs + (a1->length<<1)*(res->limbs+1);
      temp->length = a2->length+b2->length-1;
      _Zpoly_mpn_sub(prodsum, prodsum, temp);
      
      // res_mid += prodsum
      temp->coeffs = res->coeffs + a1->length*(res->limbs+1);
      temp->length = prodsum->length;
      _Zpoly_mpn_add(temp, temp, prodsum);
      
      res->length = a->length + b->length - 1;
     
   } else 
   {
      Zpoly_mpn_t scratch2, temp1; 

      while ((1<<l2)<a1->length) l2++;
      if ((1<<l2) < a->length) a1->length = (1<<l2);
      a2->length = a->length-a1->length;
      a1->coeffs = a->coeffs;
      a2->coeffs = a->coeffs+a1->length*(a->limbs+1);

      // res_lo = a1*b
      __Zpoly_mpn_karamul_recursive(res,a1,b,scratch,scratchb);
      
      //temp = a2*b
      temp->coeffs = scratch->coeffs;
      temp->length = a2->length + b->length - 1;
      temp->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs+temp->length*(scratch->limbs+1);
      scratch2->limbs = scratch->limbs;
      if (b->length <= a2->length) __Zpoly_mpn_karamul_recursive(temp,a2,b,scratch2,scratchb);
      else __Zpoly_mpn_karamul_recursive(temp,b,a2,scratch2,scratchb);
      
      // res_mid += temp
      temp1->coeffs = res->coeffs+a1->length*(res->limbs+1);
      temp1->length = temp->length;
      temp1->limbs = res->limbs;
      _Zpoly_mpn_add(temp1,temp1,temp);
  
      res->length = a->length + b->length - 1;
   } 
}

void _Zpoly_mpn_mul_karatsuba(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   unsigned long limbs = input1->limbs+input2->limbs+2;
   unsigned long log_length = 0;
   Zpoly_mpn_t scratch, scratchb;
   scratch->coeffs = (mp_limb_t *) flint_stack_alloc(5*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs;
   scratchb->limbs = FLINT_MAX(input1->limbs,input2->limbs)+1;
   scratchb->coeffs = (mp_limb_t *) flint_stack_alloc(5*FLINT_MAX(input1->length,input2->length)*(scratchb->limbs+1));
   
   if (input1->length >= input2->length)
       __Zpoly_mpn_karamul_recursive(output, input1, input2, scratch, scratchb);
   else
       __Zpoly_mpn_karamul_recursive(output, input2, input1, scratch, scratchb);
   
   flint_stack_release(); flint_stack_release();
}

void _Zpoly_mpn_mul_KS(Zpoly_mpn_t output, Zpoly_mpn_p input1, Zpoly_mpn_p input2)
{
   long sign1 = 1L;
   long sign2 = 1L;

   _Zpoly_mpn_normalise(input1);
   _Zpoly_mpn_normalise(input2);
   
   if (input2->length > input1->length) SWAP(input1, input2);
   
   if ((input1->length == 0) || (input2->length == 0)) 
   {
      _Zpoly_mpn_zero(output);
      return;
   }
   
              if ((long) input1->coeffs[(input1->length-1)*(input1->limbs+1)] < 0)
   {
      _Zpoly_mpn_negate(input1, input1);
      sign1 = -1L;
   }
   
   if (input1 != input2)
   {
      if ((long) input2->coeffs[(input2->length-1)*(input2->limbs+1)] < 0)
      {
         _Zpoly_mpn_negate(input2, input2);
         sign2 = -1L;
      }
   } else sign2 = sign1;
   
   long bits1, bits2;
   int bitpack = 0;
   
   if ((input1->limbs == 1) && (input2->limbs == 1))
   {
      bits1 = _Zpoly_mpn_bits1(input1);
      bits2 = (input1 == input2) ? bits1 : _Zpoly_mpn_bits1(input2);
   } else
   {
      bits1 = _Zpoly_mpn_bits(input1);
      bits2 = (input1 == input2) ? bits1 : _Zpoly_mpn_bits(input2);
   }
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = input2->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   unsigned long limbs = (bits-1)/FLINT_BITS_PER_LIMB + 1;
   
   if ((bits < 64) && (input1->limbs == 1) && (input2->limbs == 1) && (output->limbs == 1)) bitpack = 1;
   
   unsigned long bytes = ((bits-1)>>3)+1;
   
   ZmodFpoly_t poly1, poly2, poly3;
   if (bitpack)
   {
      ZmodFpoly_stack_init(poly1, 0, (bits*input1->length-1)/FLINT_BITS_PER_LIMB+1, 0);
      if (input1 != input2)
         ZmodFpoly_stack_init(poly2, 0, (bits*input2->length-1)/FLINT_BITS_PER_LIMB+1, 0);

      if (sign) bits = -1L*bits;
      ZmodFpoly_bit_pack_mpn(poly2, input2, input2->length, bits);
      if (input1 != input2)
         ZmodFpoly_bit_pack_mpn(poly1, input1, input1->length, bits);

      bits=ABS(bits);
   } else
   {
      ZmodFpoly_stack_init(poly1, 0, ((bytes*input1->length-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);
      if (input1 != input2)
         ZmodFpoly_stack_init(poly2, 0, ((bytes*input2->length-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);

      ZmodFpoly_byte_pack_mpn(poly1, input1, input1->length, bytes);
      if (input1 != input2)
         ZmodFpoly_byte_pack_mpn(poly2, input2, input2->length, bytes);
   }
   
   if (input1 == input2)
   {
      poly2->coeffs = poly1->coeffs;
      poly2->n = poly1->n;
   }
   
   ZmodFpoly_stack_init(poly3, 0, poly1->n + poly2->n, 0);
           
   Z_mpn_mul(poly3->coeffs[0], poly1->coeffs[0], poly1->n, poly2->coeffs[0], poly2->n);
   
   poly3->coeffs[0][poly1->n+poly2->n] = 0;
   poly3->length = 1;
   
   output->length = input1->length+input2->length-1;
  
   for (unsigned long i = 0; i < output->length; i++)
      output->coeffs[i*(output->limbs+1)] = 0;
      
   if (bitpack)
   {
      if (sign) ZmodFpoly_bit_unpack_mpn(output, poly3, input1->length+input2->length-1, bits);  
      else ZmodFpoly_bit_unpack_unsigned_mpn(output, poly3, input1->length+input2->length-1, bits);  
   } else
   {
      if (sign) ZmodFpoly_byte_unpack_mpn(output, poly3->coeffs[0], input1->length+input2->length-1, bytes);        
      else ZmodFpoly_byte_unpack_unsigned_mpn(output, poly3->coeffs[0], input1->length+input2->length-1, bytes);  
   }
   
   ZmodFpoly_stack_clear(poly3);
   if (input1 != input2)
      ZmodFpoly_stack_clear(poly2);
   ZmodFpoly_stack_clear(poly1);
     
   if ((long) (sign1 ^ sign2) < 0) _Zpoly_mpn_negate(output, output);
   
   if (sign1 < 0) _Zpoly_mpn_negate(input1, input1);
   if ((sign2 < 0) && (input1 != input2)) _Zpoly_mpn_negate(input2, input2);
}


void _Zpoly_mpn_mul_SS(Zpoly_mpn_t output, Zpoly_mpn_p input1, Zpoly_mpn_p input2)
{
   _Zpoly_mpn_normalise(input1);
   _Zpoly_mpn_normalise(input2);
   
   if (input1->length < input2->length) SWAP(input1, input2);
   
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      _Zpoly_mpn_zero(output);
      return;
   }
   
   unsigned long size1 = input1->limbs;
   unsigned long size2 = input2->limbs;
   
   unsigned long log_length = 0;
   while ((1<<log_length) < length1) log_length++;
   unsigned long log_length2 = 0;
   while ((1<<log_length2) < length2) log_length2++;
   
   /* Start with an upper bound on the number of bits needed */
   
   unsigned long output_bits = FLINT_BITS_PER_LIMB * (size1 + size2) + log_length2 + 2;

   if (output_bits <= length1)
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   else 
      output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
      
   unsigned long n = (output_bits - 1) / FLINT_BITS_PER_LIMB + 1;
   
   n = ZmodF_mul_precomp_get_feasible_n(NULL, n);

   ZmodFpoly_t poly1, poly2, res;
   long bits1, bits2;
   unsigned long sign = 0;
   
   ZmodFpoly_stack_init(poly1, log_length + 1, n, 1);
   ZmodFpoly_stack_init(poly2, log_length + 1, n, 1);
   ZmodFpoly_stack_init(res, log_length + 1, n, 1);
   
   bits1 = ZmodFpoly_convert_in_mpn(poly1, input1);
   bits2 = ZmodFpoly_convert_in_mpn(poly2, input2);
   
   if ((bits1 < 0) || (bits2 < 0)) 
   {
      sign = 1;  
      bits1 = ABS(bits1);
      bits2 = ABS(bits2);
   }
   
   /* Recompute the length of n now that we know how large everything really is */
   
   output_bits = bits1 + bits2 + log_length2 + sign;
   
   if (output_bits <= length1)
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   else 
      output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
      
   n = (output_bits - 1) / FLINT_BITS_PER_LIMB + 1;
   
   n = ZmodF_mul_precomp_get_feasible_n(NULL, n);
   
   ZmodFpoly_decrease_n(poly1, n);
   ZmodFpoly_decrease_n(poly2, n);
   ZmodFpoly_decrease_n(res, n);
                    
   ZmodFpoly_convolution(res, poly1, poly2);
   ZmodFpoly_normalise(res);
          
   output->length = length1 + length2 - 1;
   
   ZmodFpoly_convert_out_mpn(output, res, sign);
   
   ZmodFpoly_stack_clear(res);
   ZmodFpoly_stack_clear(poly2);
   ZmodFpoly_stack_clear(poly1);
}

void _Zpoly_mpn_mul(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   if ((input1->length <= 2) && (input2->length <= 2)) 
   {
      _Zpoly_mpn_mul_naive(output, input1, input2);
      return;
   }
   
   if (input1->limbs + input2->limbs <= 512/FLINT_BITS_PER_LIMB)
   {
      _Zpoly_mpn_mul_KS(output, input1, input2);
      return;
   }
   
   if (input1->length + input2->length <= 32) 
   {
      _Zpoly_mpn_mul_karatsuba(output, input1, input2);
      return;
   }
   
   unsigned long bits1 = _Zpoly_mpn_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _Zpoly_mpn_bits(input2);
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _Zpoly_mpn_mul_SS(output, input1, input2);
      return;
   } 
   
   _Zpoly_mpn_mul_KS(output, input1, input2);     
}

/****************************************************************************

   Zpoly_mpn_* layer

****************************************************************************/

/* 
   Create a polynomial of length zero with "alloc" allocated coefficients
   each with enough space for limbs limbs
*/

void Zpoly_mpn_init(Zpoly_mpn_t poly, unsigned long alloc, unsigned long limbs)
{
   FLINT_ASSERT(alloc >= 1);
   FLINT_ASSERT(limbs >= 1);

   poly->coeffs = (mp_limb_t *) flint_heap_alloc(alloc*(limbs+1));
   
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

/* Shrink or expand a polynomial to "alloc" coefficients */

void Zpoly_mpn_realloc(Zpoly_mpn_t poly, unsigned long alloc)
{
   if (alloc <= 0) alloc = 1;
   poly->coeffs = (mp_limb_t*) flint_heap_realloc(poly->coeffs, alloc*(poly->limbs+1));
   
   poly->alloc = alloc;
   
   // truncate actual data if necessary
   if (poly->length > alloc)
      poly->length = alloc;
}

void Zpoly_mpn_clear(Zpoly_mpn_t poly)
{
   flint_heap_free(poly->coeffs);
}

long Zpoly_mpn_degree(Zpoly_mpn_t poly)
{
   _Zpoly_mpn_normalise(poly);
   return poly->length - 1;
}

unsigned long Zpoly_mpn_length(Zpoly_mpn_t poly)
{
   _Zpoly_mpn_normalise(poly);
   return poly->length;
}

void Zpoly_mpn_ensure_space2(Zpoly_mpn_t poly, unsigned long alloc)
{
   FLINT_ASSERT(alloc >= poly->alloc);

   if (alloc < 2*poly->alloc)
      alloc = 2*poly->alloc;
   Zpoly_mpn_realloc(poly, alloc);
}

void Zpoly_mpn_ensure_length2(Zpoly_mpn_t poly, unsigned long length)
{
   FLINT_ASSERT(length >= poly->length);
   
   Zpoly_mpn_ensure_space(poly, length);

   for (unsigned long i = poly->length; i < length; i++)
      poly->coeffs[i*(poly->limbs+1)] = 0;

   poly->length = length;
}

void Zpoly_mpn_get_coeff_mpz(mpz_t x, Zpoly_mpn_t poly, unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(x, 0);
   else
      _Zpoly_mpn_get_coeff_mpz(x, poly, n);
}
