/****************************************************************************

fmpz_poly.c: Polynomials over Z, implemented as contiguous block of fmpz_t's

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include "mpz_poly.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"
#include "extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "Z_mpn.h"
#include "ZmodF_mul.h"
#include "Z_mpn_mul-tuning.h"

#include "mpz_poly.h"

/****************************************************************************

   _fmpz_poly_* layer

****************************************************************************/

/* 
   Create a polynomial of length zero with "alloc" allocated coefficients
   each with enough space for limbs limbs
*/

void _fmpz_poly_stack_init(fmpz_poly_t poly, unsigned long alloc, unsigned long limbs)
{
   FLINT_ASSERT(alloc >= 1);
   FLINT_ASSERT(limbs >= 1);

   poly->coeffs = (mp_limb_t *) flint_stack_alloc(alloc*(limbs+1));
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

void _fmpz_poly_stack_clear(fmpz_poly_t poly)
{
   flint_stack_release();
}

/*
   Used for debugging polynomial code
   Checks that length <= alloc and that both are positive or zero
   Checks that limbs >= 0 otherwise
   Checks that each coefficient has at most _limbs_ limbs 
*/

void _fmpz_poly_check(fmpz_poly_t poly)
{
   if ((long) poly->alloc < 0)
   {
      printf("Error: Poly alloc < 0\n");
      abort();
   }
   if ((long) poly->length < 0)
   {
      printf("Error: Poly length < 0\n");
      abort();
   }
   if (poly->length > poly->alloc) 
   {
      printf("Error: Poly length = %ld > alloc = %ld\n", poly->length, poly->alloc);
      abort();
   }
   if ((long) poly->limbs < 0) 
   {
      printf("Error: Poly limbs < 0\n");
      abort();
   }
   for (unsigned long i = 0; i < poly->length; i++)
   {
      if (FLINT_ABS(poly->coeffs[i*(poly->limbs+1)]) > poly->limbs)
      {
         printf("Error: coefficient %ld is too large (%ld limbs vs %ld limbs)\n", 
                        i, FLINT_ABS(poly->coeffs[i*(poly->limbs+1)]), poly->limbs);
         abort();
      }
   }
}


// retrieves coefficient #n as an mpz, no bounds checking
void _fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, unsigned long n)
{
   FLINT_ASSERT(n < poly->length);
   fmpz_to_mpz(x, poly->coeffs + n*(poly->limbs + 1));
}


/* 
   Set a coefficient to the given unsigned value.
   If x is nonzero, poly->limbs must be positive.
   Assumes the polynomial length is greater than n.
*/
void _fmpz_poly_set_coeff_ui(fmpz_poly_t poly, unsigned long n, unsigned long x)
{
   FLINT_ASSERT(poly->alloc > n);
   fmpz_set_ui(poly->coeffs + n*(poly->limbs + 1), x);
}

/* 
   Set a coefficient to the given signed value.
   If x is nonzero, poly->limbs must be positive.
   Assumes the polynomial length is greater than n.
*/

void _fmpz_poly_set_coeff_si(fmpz_poly_t poly, unsigned long n, long x)
{
   FLINT_ASSERT(poly->length > n);
   fmpz_set_si(poly->coeffs + n*(poly->limbs + 1), x);
}

void _fmpz_poly_normalise(fmpz_poly_t poly)
{
   while (poly->length && poly->coeffs[(poly->length-1)*(poly->limbs+1)] == 0)
      poly->length--;
}

/* 
   Sets the output poly to equal the input poly 
   Assumes the output poly is big enough to hold the nonzero limbs of the input poly
*/

void _fmpz_poly_set(fmpz_poly_t output, fmpz_poly_t input)
{
   if (input->length == 0) 
   {
      output->length = 0;
      return;
   }
   
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

void _fmpz_poly_swap(fmpz_poly_t x, fmpz_poly_t y)
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
long _fmpz_poly_bits1(fmpz_poly_t poly_mpn)
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
            if (bits == FLINT_BITS) break;
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
long _fmpz_poly_bits(fmpz_poly_t poly_mpn)
{
   if (poly_mpn->limbs == 0) return 0;
   if (poly_mpn->limbs == 1) return _fmpz_poly_bits1(poly_mpn);
   
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
         if (bits == FLINT_BITS) mask = 0L;
         else mask = -1L - ((1L<<bits)-1);
      } else if (ABS(size_j) == limbs+1)
      {
         if (coeffs_m[j+ABS(size_j)] & mask)
         {
            bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]);   
            if (bits == FLINT_BITS) mask = 0L;
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
   
   return sign*(FLINT_BITS*limbs+bits);
}

int _fmpz_poly_equal(fmpz_poly_p input1, fmpz_poly_p input2)
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

void _fmpz_poly_neg(fmpz_poly_t output, fmpz_poly_t input)
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

void _fmpz_poly_left_shift(fmpz_poly_t output, fmpz_poly_t input, 
                                                 unsigned long n)
{
   fmpz_poly_t part;   
   
   part->length = input->length;
   part->limbs = output->limbs;
   part->coeffs = output->coeffs + n*(output->limbs+1);
      
   _fmpz_poly_set(part, input);
   for (unsigned long i = 0; i < n; i++) output->coeffs[i*(output->limbs+1)] = 0;
   
   output->length = input->length + n;
}

/* 
   Divides input by x^n losing the remainder and sets output to the result
   Assumes output is large enough to contain the result
*/

void _fmpz_poly_right_shift(fmpz_poly_t output, fmpz_poly_t input, unsigned long n)
{
   if (input->length <= n) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   fmpz_poly_t part;
   part->length = input->length - n;
   part->limbs = input->limbs;
   part->coeffs = input->coeffs + n*(input->limbs + 1);
   _fmpz_poly_set(output, part);
}

/* 
   Sets output to the reverse of input (i.e. reverse the order of the coefficients)
   assuming input to be a polynomial with _length_ coefficients (it may have a length
   that is less than _length_).
*/ 

void _fmpz_poly_reverse(fmpz_poly_t output, fmpz_poly_t input, unsigned long length)
{
   unsigned long coeff_limbs;
   unsigned long size_in = input->limbs + 1;
   unsigned long size_out = output->limbs + 1;
   long i;
   
   if (input != output)
   {
      for (i = 0; i < FLINT_MIN(length, input->length); i++)
      {
         coeff_limbs = ABS(input->coeffs[i*size_in]) + 1;
         copy_limbs(output->coeffs + (length - i - 1)*size_out, input->coeffs + i*size_in, coeff_limbs);
      }
      for ( ; i < length; i++)
      {
         output->coeffs[(length - i - 1)*size_out] = 0;
      }
      output->length = length;
      _fmpz_poly_normalise(output);
   } else
   {
      mp_limb_t * temp = (mp_limb_t *) flint_stack_alloc(size_in);
      unsigned long coeff_limbs2;
      
      for (i = 0; i < length/2; i++)
      {
         if (i < input->length)
         {
            coeff_limbs = ABS(input->coeffs[i*size_in]) + 1;
            copy_limbs(temp, input->coeffs + i*size_in, coeff_limbs);
         } else
         {
            coeff_limbs = 1;
            temp[0] = 0;            
         }
         if (length - i - 1 < input->length)
         {
            coeff_limbs2 = ABS(input->coeffs[(length - i - 1)*size_in]) + 1;
            copy_limbs(input->coeffs + i*size_in, input->coeffs + (length - i - 1)*size_in, coeff_limbs2);
         } else
         {
            input->coeffs[i*size_in] = 0;
         }
         copy_limbs(input->coeffs + (length - i - 1)*size_in, temp, coeff_limbs);
      }
      if ((length & 1) && (i >= input->length)) input->coeffs[i*size_in] = 0;

      input->length = length;
      _fmpz_poly_normalise(input);
   }
}

/*
    Adds two coefficients together
*/

void __fmpz_poly_add_coeffs(mp_limb_t * coeffs_out, mp_limb_t * coeffs1, mp_limb_t * coeffs2)
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
     Add an unsigned long to a polynomial coefficient
*/

void __fmpz_poly_add_coeff_ui(mp_limb_t * output, unsigned long x)
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

void __fmpz_poly_add_coeff2_ui(mp_limb_t * output, unsigned long x)
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

void __fmpz_poly_sub_coeff_ui(mp_limb_t * output, unsigned long x)
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

void _fmpz_poly_add(fmpz_poly_t output, fmpz_poly_p input1, fmpz_poly_p input2)
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
      __fmpz_poly_add_coeffs(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2+i*size2);
   }    
   
   if (input1 != output)
   {
      for (unsigned long i = shorter; i < input1->length; i++)
      {
          copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, ABS(coeffs1[i*size1])+1);
      }
   }
   if (input2 != output)
   {
      for (unsigned long i = shorter; i < input2->length; i++)
      {
         copy_limbs(coeffs_out+i*size_out, coeffs2+i*size2, ABS(coeffs2[i*size2])+1);
      }
   }
   
   output->length = (input1->length > input2->length) ? input1->length : input2->length;
}

void __fmpz_poly_sub_coeffs(mp_limb_t * coeffs_out, mp_limb_t * coeffs1, mp_limb_t * coeffs2)
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

/* 
    Subtract two polynomials
*/

void _fmpz_poly_sub(fmpz_poly_t output, fmpz_poly_p input1, fmpz_poly_p input2)
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
      __fmpz_poly_sub_coeffs(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2+i*size2);
   }
   
   if (input1 != output)
   {
      for (unsigned long i = shorter; i < input1->length; i++)
      {
         copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, ABS(coeffs1[i*size1])+1);
      }
   }
   if (input2 != output)
   {
      for (unsigned long i = shorter; i < input2->length; i++)
      {
         copy_limbs(coeffs_out+i*size_out+1, coeffs2+i*size2+1, ABS(coeffs2[i*size2]));
         coeffs_out[i*size_out] = -coeffs2[i*size2];
      }
   } else
   {
      for (unsigned long i = shorter; i < input2->length; i++)
      {
         coeffs_out[i*size_out] = -coeffs2[i*size2];
      }
   }

   output->length = (input1->length > input2->length) ? input1->length : input2->length;
}

/* 
   Multiplies two coefficients
   Assumes no overlap
*/

void __fmpz_poly_mul_coeffs(mp_limb_t * res, mp_limb_t * a, mp_limb_t * b) 
{
      NORM(a);
      NORM(b);
      unsigned long sizea = ABS(a[0]);
      unsigned long sizeb = ABS(b[0]);
      mp_limb_t mslimb;
      mp_limb_t * temp;
      
      if ((sizea == 0) || (sizeb == 0))
      {
        res[0] = 0;
      } else if (sizea + sizeb < 100)
      {
         temp = (mp_limb_t *) flint_stack_alloc_small(sizea + sizeb + 1);
         if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
         else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
         temp[0] = sizea + sizeb - (mslimb == 0);
         copy_limbs(res, temp, temp[0]+1);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
         flint_stack_release_small();     
      } else if (sizea + sizeb < 2*FLINT_FFT_LIMBS_CROSSOVER)
      {
         temp = (mp_limb_t *) flint_stack_alloc(sizea + sizeb + 1);
         if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
         else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
         temp[0] = sizea + sizeb - (mslimb == 0);
         copy_limbs(res, temp, temp[0]+1);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
         flint_stack_release();   
      } else
      {
         if (sizea >= sizeb) mslimb = Z_mpn_mul(res+1, a+1, sizea, b+1, sizeb);
         else mslimb = Z_mpn_mul(res+1, b+1, sizeb, a+1, sizea);
         res[0] = sizea+sizeb - (mslimb == 0);
         if ((long) (a[0] ^ b[0]) < 0) res[0] = -res[0];
      }
}      

/* 
   Multiplies two coefficients but assumes that there is space in each coeff
   of res for the number of limbs of _a_ plus the number of limbs of _b_, whenever 
   this sum is less than 2*FLINT_FFT_LIMBS_CROSSOVER  
.
   Assumes no overlap
*/

void __fmpz_poly_mul_coeffs2(mp_limb_t * res, mp_limb_t * a, mp_limb_t * b) 
{
      NORM(a);
      NORM(b);
      unsigned long sizea = ABS(a[0]);
      unsigned long sizeb = ABS(b[0]);
      mp_limb_t mslimb;
      mp_limb_t * temp;
      
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

/* 
   Sets out to out+in1*in2
   Assumes no overlap
   Assumes all coefficients are only a small number of limbs (< 100 say)
*/

void __fmpz_poly_addmul_coeffs(mp_limb_t * res, mp_limb_t * a, mp_limb_t * b) 
{
      NORM(a);
      NORM(b);
      unsigned long sizea = ABS(a[0]);
      unsigned long sizeb = ABS(b[0]);
      mp_limb_t * temp;
      mp_limb_t mslimb;
      
      if (sizea && sizeb)
      {
         if (sizea + sizeb < 100)
         {
            temp = (mp_limb_t *) flint_stack_alloc_small(sizea + sizeb + 1);
            if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
            else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
            temp[0] = sizea + sizeb - (mslimb == 0);
            if ((long) (a[0] ^ b[0]) < 0) temp[0] = -temp[0];
            __fmpz_poly_add_coeffs(res, res, temp);
            flint_stack_release_small();
         } else
         {
            temp = (mp_limb_t *) flint_stack_alloc(sizea + sizeb + 1);
            if (sizea >= sizeb) mslimb = mpn_mul(temp+1, a+1, sizea, b+1, sizeb);
            else mslimb = mpn_mul(temp+1, b+1, sizeb, a+1, sizea);
            temp[0] = sizea + sizeb - (mslimb == 0);
            if ((long) (a[0] ^ b[0]) < 0) temp[0] = -temp[0];
            __fmpz_poly_add_coeffs(res, res, temp);
            flint_stack_release();
         }        
      }     
}      

void _fmpz_poly_scalar_mul_ui(fmpz_poly_t output, fmpz_poly_t poly, unsigned long x)
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

void _fmpz_poly_scalar_mul_si(fmpz_poly_t output, fmpz_poly_t poly, long x)
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

void _fmpz_poly_scalar_div_exact_ui(fmpz_poly_t output, fmpz_poly_t poly, unsigned long x)
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

void _fmpz_poly_scalar_div_exact_si(fmpz_poly_t output, fmpz_poly_t poly, long x)
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

void _fmpz_poly_scalar_div_ui(fmpz_poly_t output, fmpz_poly_t poly, unsigned long x)
{
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs1 = poly->coeffs;
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      
      count_lead_zeros(norm, x);
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

void _fmpz_poly_scalar_div_si(fmpz_poly_t output, fmpz_poly_t poly, long x)
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
      
      count_lead_zeros(norm, (unsigned long) x);
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

/*
   Multiply two polynomials using the naive technique.
   Currently doesn't allow aliasing
*/

void _fmpz_poly_mul_naive(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
{
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs1, * coeffs2;
   unsigned long len1, len2; 
   unsigned long lenm1;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   len1 = input1->length;
   len2 = input2->length;
      
   // Special case if the length of both inputs is 1
   if ((len1 == 1) && (len2 == 1))
   {
      if ((coeffs1[0] == 0) || (coeffs2[0] == 0))
      {
         coeffs_out[0] = 0;
      } else
      {
         __fmpz_poly_mul_coeffs(coeffs_out, coeffs1, coeffs2);
      }      
   }         
   // Ordinary case
   else
   {
      unsigned long size_out = output->limbs+1;
      unsigned long size1, size2;
      size1 = input1->limbs+1;
      size2 = input2->limbs+1;
      lenm1 = input1->length-1;
      
      mp_limb_t * temp;
      long i, j;
      
      for (i = 0; i < len1; i++)
      {
         /* Set out[i] = in1[i]*in2[0] */
         if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
         {
            coeffs_out[i*size_out] = 0;
         } else
         {
            __fmpz_poly_mul_coeffs2(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
         }
      }
      for (i = 1; i < len2 - 1; i++)
      {
         /* Set out[i+in1->length-1] = in1[in1->length-1]*in2[i] */
         if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
         {
            coeffs_out[(i+lenm1)*size_out]=0;
         } else
         {
            __fmpz_poly_mul_coeffs2(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
         }      
      }
      /* 
         The above coefficient multiplications overwrite the first limb of the next coefficient
         in each case, using the function __fmpz_poly_mul_coeffs2. The final multiplication 
         cannot do this however.
      */
      if ((coeffs1[lenm1*size1] == 0) || (coeffs2[(len2-1)*size2] == 0))
      {
         coeffs_out[(len2+lenm1-1)*size_out]=0;
      } else
      {
         __fmpz_poly_mul_coeffs(coeffs_out+(len2+lenm1-1)*size_out, coeffs1+lenm1*size1, coeffs2+(len2-1)*size2);
      }      
      
      for (i = 0; i < lenm1; i++)
      {      
         for (j = 1; j < len2; j++)
         {
            /* out[i+j] += in1[i]*in2[j] */
            if ((coeffs1[i*size1] != 0) && (coeffs2[j*size2] != 0))
            {
               if (!coeffs_out[(i+j)*size_out])
               {
                  __fmpz_poly_mul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } else 
               {
                  __fmpz_poly_addmul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } 
            }
         }
      }
   } 
   
   output->length = len1 + len2 - 1;
}

void _fmpz_poly_truncate(fmpz_poly_t poly, unsigned long trunc)
{
   if (poly->length > trunc) poly->length = trunc;
}

/*
   Multiply two polynomials using the naive technique truncating the result to trunc terms.
   Currently doesn't allow aliasing
*/

void _fmpz_poly_mul_naive_trunc(fmpz_poly_t output, fmpz_poly_t input1, 
                                          fmpz_poly_t input2, unsigned long trunc)
{
   mp_limb_t * coeffs_out = output->coeffs;
   unsigned long size_out = output->limbs+1;
   mp_limb_t * coeffs1, * coeffs2;
   unsigned long size1, size2;
   unsigned long len1, len2; 
   unsigned long lenm1;
   
   if (trunc == 0) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   if ((input1->length == 0) || (input2->length == 0)) 
   {
      for (unsigned long i = 0; i < trunc; i++)
      {
         coeffs_out[i*size_out] = 0;
      }
      _fmpz_poly_zero(output);
      return;      
   }
      
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   lenm1 = input1->length-1;
   len1 = input1->length;
   len2 = input2->length;
   
   long i, j;
      
   mp_limb_t * temp;
            
   // Special case if the length of both inputs is 1
   if ((len1 == 1) && (len2 == 1))
   {
      if ((coeffs1[0] == 0) || (coeffs2[0] == 0))
      {
         coeffs_out[0] = 0;
      } else
      {
         __fmpz_poly_mul_coeffs(coeffs_out, coeffs1, coeffs2);
      }      
   }
   // Ordinay case
   else
   {
      for (i = 0; (i < len1) && (i < trunc - 1); i++)
      {
         /* Set out[i] = in1[i]*in2[0] */
         if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
         {
            coeffs_out[i*size_out] = 0;
         } else
         {
            __fmpz_poly_mul_coeffs2(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
         }
      }
      if (i != len1)
      {
         if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
         {
            coeffs_out[i*size_out] = 0;
         } else
         {
            __fmpz_poly_mul_coeffs(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
         }
      } else
      {
         for (i = 1; (i < len2 - 1) && (i + lenm1 < trunc - 1); i++)
         {
            /* Set out[i+in1->length-1] = in1[in1->length-1]*in2[i] */
            if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
            {
               coeffs_out[(i+lenm1)*size_out] = 0;
            } else
            {
               __fmpz_poly_mul_coeffs2(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
            }      
         }
         
         /* 
            The above coefficient multiplications overwrite the first limb of the next coefficient
            in each case, using the function __fmpz_poly_mul_coeffs2. The final multiplication 
            cannot do this however.
         */
         if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
         {
            coeffs_out[(i+lenm1)*size_out] = 0;
         } else
         {
            __fmpz_poly_mul_coeffs(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
         }     
      }
         
      for (i = 0; i < lenm1; i++)
      {      
         for (j = 1; (j < len2) && (i + j < trunc); j++)
         {
            /* out[i+j] += in1[i]*in2[j] */
            if ((coeffs1[i*size1] != 0) && (coeffs2[j*size2] != 0))
            {
               if (!coeffs_out[(i+j)*size_out])
               {
                  __fmpz_poly_mul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } else 
               {
                  __fmpz_poly_addmul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } 
            }
         }
      }
   } 
   
   output->length = FLINT_MIN(len1 + len2 - 1, trunc);
}

/*
   Multiply two polynomials using the naive technique truncating the result 
   so that the first trunc terms are zero.
   Currently doesn't allow aliasing
*/

void _fmpz_poly_mul_naive_trunc_left(fmpz_poly_t output, fmpz_poly_t input1, 
                                          fmpz_poly_t input2, unsigned long trunc)
{
   mp_limb_t * coeffs_out = output->coeffs;
   unsigned long size_out = output->limbs+1;
   mp_limb_t * coeffs1, * coeffs2;
   unsigned long size1, size2;
   unsigned long len1, len2; 
   unsigned long lenm1;
   
   if ((input1->length == 0) || (input2->length == 0) || (trunc >= input1->length + input2->length - 1)) 
   {
      for (unsigned long i = 0; i < input1->length + input2->length - 1; i++)
      {
         coeffs_out[i*size_out] = 0;
      }
      _fmpz_poly_zero(output);
      return;      
   }
      
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   lenm1 = input1->length-1;
   len1 = input1->length;
   len2 = input2->length;
   
   long i, j;
      
   mp_limb_t * temp;
            
   // Special case if the length of both inputs is 1
   if ((len1 == 1) && (len2 == 1))
   {
      if ((coeffs1[0] == 0) || (coeffs2[0] == 0))
      {
         coeffs_out[0] = 0;
      } else
      {
         __fmpz_poly_mul_coeffs(coeffs_out, coeffs1, coeffs2);
      }      
   }
   // Ordinay case
   else
   {
      for (i = trunc; (i < len1); i++)
      {
         /* Set out[i] = in1[i]*in2[0] */
         if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
         {
            coeffs_out[i*size_out] = 0;
         } else
         {
            __fmpz_poly_mul_coeffs2(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
         }
      }
      for (i = 1; i < len2 - 1; i++)
      {
         if (i + lenm1 >= trunc)
         {
            /* Set out[i+in1->length-1] = in1[in1->length-1]*in2[i] */
            if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
            {
               coeffs_out[(i+lenm1)*size_out] = 0;
            } else
            {
               __fmpz_poly_mul_coeffs2(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
            }
         }      
      }
      if (len2 == 1) i = 0;   
      /* 
         The above coefficient multiplications overwrite the first limb of the next coefficient
         in each case, using the function __fmpz_poly_mul_coeffs2. The final multiplication 
         cannot do this however.
      */
      if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
      {
         coeffs_out[(i+lenm1)*size_out] = 0;
      } else
      {
         __fmpz_poly_mul_coeffs(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
      }     
         
      for (i = 0; i < lenm1; i++)
      {      
         for (j = 1; j < len2; j++)
         {
            /* out[i+j] += in1[i]*in2[j] */
            if ((coeffs1[i*size1] != 0) && (coeffs2[j*size2] != 0) && (i + j >= trunc))
            {
               if (!coeffs_out[(i+j)*size_out])
               {
                  __fmpz_poly_mul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } else 
               {
                  __fmpz_poly_addmul_coeffs(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } 
            }
         }
      }
   } 
   
   for (i = 0; (i < trunc) && (i < len1 + len2 - 1); i++)
   {
      coeffs_out[i*size_out] = 0;
   }
   
   output->length = len1 + len2 - 1;
}

unsigned long _fmpz_poly_max_limbs(fmpz_poly_t poly)
{
   unsigned long limbs = poly->limbs;
   unsigned long max_limbs = 0;
   unsigned long next_limbs;
   
   for (unsigned long i = 0; (i < poly->length) && (max_limbs != limbs); i++)
   {
       next_limbs = ABS(poly->coeffs[i*(limbs+1)]);
       if (next_limbs > max_limbs) max_limbs = next_limbs;
   } 
   return max_limbs;
}

void __fmpz_poly_karamul_recursive(fmpz_poly_t res, fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, unsigned long crossover)
{
   fmpz_poly_t temp;
   
   if ((crossover < 4) && (a->length == 2 && b->length == 2)) {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      __fmpz_poly_mul_coeffs2(res->coeffs, a->coeffs, b->coeffs); 
      __fmpz_poly_add_coeffs(scratchb->coeffs, a->coeffs, a->coeffs+asize);
      __fmpz_poly_mul_coeffs(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
      __fmpz_poly_add_coeffs(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
      __fmpz_poly_mul_coeffs(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
      __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
      __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
      res->length = a->length + b->length - 1;
      
      return;
   }
   
   if ((a->length+b->length <= crossover) ||  (a->length <= 1) || (b->length <= 1) ||  ((a->length == 2) || (b->length == 2)))
   {
      _fmpz_poly_mul_naive(res, a, b);
      
      return;
   }  
        
   fmpz_poly_t a1,a2,b1,b2;
      
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
      
      /* 
         from 0 for 2 * a1->length - 1, from 2 * a1->length for a2->length + b2->length - 1
         will be written directly to, so we need to clean the coefficient in between
      */
      res->coeffs[((a1->length<<1)-1)*(res->limbs+1)] = 0;
      
      fmpz_poly_t asum, bsum, prodsum, scratch2, scratch3;
     
      asum->length = a1->length;
      asum->coeffs = scratchb->coeffs;
      asum->limbs = scratchb->limbs;
      bsum->length = a1->length;
      bsum->coeffs = scratchb->coeffs + a1->length*(scratchb->limbs+1);
      bsum->limbs = scratchb->limbs;
      prodsum->length = (a1->length<<1)-1;
      prodsum->coeffs = scratch->coeffs;// + (a1->length<<1)*(scratch->limbs+1);
      prodsum->limbs = scratch->limbs;
      
      // res_lo = a1*b1
      __fmpz_poly_karamul_recursive(res, a1, b1, scratch, scratchb, crossover);
      
      // res_hi = a2*b2
      temp->coeffs = res->coeffs+(a1->length<<1)*(res->limbs+1);
      temp->limbs = res->limbs;
      __fmpz_poly_karamul_recursive(temp, a2, b2, scratch, scratchb, crossover);
      
      // asum = a1+a2
      _fmpz_poly_add(asum, a1, a2);
      // bsum = b1+b2
      _fmpz_poly_add(bsum, b1, b2);
      // prodsum = asum*bsum
      scratch3->coeffs = scratchb->coeffs+(a1->length<<1)*(scratchb->limbs+1);
      scratch3->limbs = scratchb->limbs;
      
      scratch2->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs+((a1->length<<1)-1)*(scratch->limbs+1);
      __fmpz_poly_karamul_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover);
      
      // prodsum = prodsum - res_lo
      temp->coeffs = res->coeffs;
      temp->length = (a1->length<<1)-1;
      _fmpz_poly_sub(prodsum, prodsum, temp);
       
      // prodsum = prodsum - res_hi
      temp->coeffs = res->coeffs + (a1->length<<1)*(res->limbs+1);
      temp->length = a2->length+b2->length-1;
      _fmpz_poly_sub(prodsum, prodsum, temp);
      
      // res_mid += prodsum
      temp->coeffs = res->coeffs + a1->length*(res->limbs+1);
      temp->length = prodsum->length;
      _fmpz_poly_add(temp, temp, prodsum);
      
      res->length = a->length + b->length - 1;
     
   } else 
   {
      fmpz_poly_t scratch2, temp1; 

      while ((1<<l2)<a1->length) l2++;
      if ((1<<l2) < a->length) a1->length = (1<<l2);
      a2->length = a->length-a1->length;
      a1->coeffs = a->coeffs;
      a2->coeffs = a->coeffs+a1->length*(a->limbs+1);

      /* 
         The first a1->length + b->length - 1 coefficients will be written to directly, 
         so we need to clean the remaining coefficients
      */
      for (unsigned long i = a1->length + b->length - 1; i < a->length + b->length - 1; i++)
         res->coeffs[i*(res->limbs+1)] = 0;
      
      // res_lo = a1*b
      __fmpz_poly_karamul_recursive(res, a1, b, scratch, scratchb, crossover);
      
      //temp = a2*b
      temp->coeffs = scratch->coeffs;
      temp->length = a2->length + b->length - 1;
      temp->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs+temp->length*(scratch->limbs+1);
      scratch2->limbs = scratch->limbs;
      if (b->length <= a2->length) __fmpz_poly_karamul_recursive(temp, a2, b, scratch2, scratchb, crossover);
      else __fmpz_poly_karamul_recursive(temp, b, a2, scratch2, scratchb, crossover);
      
      // res_mid += temp
      temp1->coeffs = res->coeffs+a1->length*(res->limbs+1);
      temp1->length = temp->length;
      temp1->limbs = res->limbs;
      _fmpz_poly_add(temp1,temp1,temp);
  
      res->length = a->length + b->length - 1;
   } 
}

void _fmpz_poly_mul_karatsuba(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
{
   unsigned long limbs = output->limbs;
   unsigned long log_length = 0;
   unsigned long crossover;
   
   fmpz_poly_t scratch, scratchb, temp;
   scratch->coeffs = (mp_limb_t *) flint_stack_alloc(5*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs + 1;
   scratchb->limbs = FLINT_MAX(input1->limbs, input2->limbs) + 1;
   scratchb->coeffs = (mp_limb_t *) flint_stack_alloc(5*FLINT_MAX(input1->length, input2->length)*(scratchb->limbs+1));
   
   if (_fmpz_poly_max_limbs(input1) + _fmpz_poly_max_limbs(input2) >= 19) crossover = 0;
   else crossover = 19 - _fmpz_poly_max_limbs(input1) - _fmpz_poly_max_limbs(input2);
   
   if (input1->length >= input2->length)
       __fmpz_poly_karamul_recursive(output, input1, input2, scratch, scratchb, crossover);
   else
       __fmpz_poly_karamul_recursive(output, input2, input1, scratch, scratchb, crossover);
   
   flint_stack_release(); 
   flint_stack_release();
}

void __fmpz_poly_karatrunc_recursive(fmpz_poly_t res, fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, unsigned long crossover, unsigned long trunc)
{
   fmpz_poly_t temp, temp2;
   
   if ((a->length <= 1) || (b->length <= 1)) 
   {
      trunc = FLINT_MIN(trunc, a->length + b->length - 1);
      _fmpz_poly_mul_naive_trunc(res, a, b, trunc);
      
      return;
   }
   
   if (((a->length == 2 && b->length == 2) && (crossover < 4)) || (trunc == 1)) {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      if (trunc > 1)
      {
         __fmpz_poly_mul_coeffs2(res->coeffs, a->coeffs, b->coeffs); 
         
         __fmpz_poly_add_coeffs(scratchb->coeffs, a->coeffs, a->coeffs+asize);
         __fmpz_poly_add_coeffs(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
         
         if (trunc > 2) __fmpz_poly_mul_coeffs(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
         else __fmpz_poly_mul_coeffs(scratch->coeffs, a->coeffs+asize, b->coeffs+bsize); 
         
         __fmpz_poly_mul_coeffs(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
         
         __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
         if (trunc > 2) __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
         else __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, scratch->coeffs);
      } else
      {
         __fmpz_poly_mul_coeffs(res->coeffs, a->coeffs, b->coeffs); 
      }
            
      res->length = FLINT_MIN(a->length + b->length - 1, trunc);
      
      return;
   }
   if ((a->length+b->length <= crossover) || ((a->length == 2) && (b->length == 2)))
   {
      trunc = FLINT_MIN(trunc, a->length + b->length - 1);
      _fmpz_poly_mul_naive_trunc(res, a, b, trunc);
      
      return;
   }   
        
   fmpz_poly_t a1, a2, b1, b2;
      
   unsigned long l2 = 0;
   
   unsigned long sa = FLINT_MIN(a->length, trunc);
   unsigned long sb = FLINT_MIN(b->length, trunc);
   unsigned long old_length;
     
   a1->length = (sa+1)/2;
   a2->length = sa-a1->length;
   a1->coeffs = a->coeffs;
   a2->coeffs = a->coeffs+a1->length*(a->limbs+1);
   a1->limbs = a->limbs;
   a2->limbs = a->limbs;
   
   if (a1->length < sb) //ordinary case
   {
      /*
         (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
      */
      
      b1->length = a1->length;
      b2->length = sb - b1->length;
      b1->coeffs = b->coeffs;
      b2->coeffs = b->coeffs + b1->length*(b->limbs+1);
      b1->limbs = b->limbs;
      b2->limbs = b->limbs;
      
      /* 
         from 0 for 2 * a1->length - 1, from 2 * a1->length for a2->length + b2->length - 1
         will be written directly to, so we need to clean the coefficient in between
      */
      if ((a1->length<<1)-1 < trunc) res->coeffs[((a1->length<<1)-1)*(res->limbs+1)] = 0;
  
      fmpz_poly_t asum, bsum, prodsum, scratch2, scratch3;
     
      asum->length = a1->length;
      asum->coeffs = scratchb->coeffs;
      asum->limbs = scratchb->limbs;
      bsum->length = a1->length;
      bsum->coeffs = scratchb->coeffs + a1->length*(scratchb->limbs+1);
      bsum->limbs = scratchb->limbs;
      prodsum->length = (a1->length<<1)-1;
      prodsum->coeffs = scratch->coeffs+(a2->length+b2->length-1)*(scratch->limbs+1);
      prodsum->limbs = scratch->limbs;
      
      // res_lo = a1*b1
      __fmpz_poly_karatrunc_recursive(res, a1, b1, scratch, scratchb, crossover, trunc);
      
      // res_hi = a2*b2
      temp->coeffs = scratch->coeffs;
      temp->limbs = scratch->limbs;
      scratch2->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs + (a2->length+b2->length-1)*(scratch->limbs+1);
      __fmpz_poly_karatrunc_recursive(temp, a2, b2, scratch2, scratchb, crossover, trunc - a1->length);
      
      temp2->limbs = res->limbs;
      if (trunc > (a1->length<<1))
      {
         old_length = temp->length;
         temp->length = FLINT_MIN(old_length, trunc-(a1->length<<1));
         temp2->coeffs = res->coeffs+(a1->length<<1)*(res->limbs+1);
         _fmpz_poly_set(temp2, temp);
         temp->length = old_length;
      }
      
      // asum = a1+a2
      _fmpz_poly_add(asum, a1, a2);
      // bsum = b1+b2
      _fmpz_poly_add(bsum, b1, b2);
      // prodsum = asum*bsum
      scratch3->coeffs = scratchb->coeffs+(a1->length<<1)*(scratchb->limbs+1);
      scratch3->limbs = scratchb->limbs;
      
      scratch2->coeffs = scratch->coeffs + (sa+sb-2)*(scratch->limbs+1);
      __fmpz_poly_karatrunc_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover, trunc - a1->length);
      
      // prodsum = prodsum - res_lo
      temp2->coeffs = res->coeffs;
      temp2->length = (a1->length<<1)-1;
      _fmpz_poly_sub(prodsum, prodsum, temp2);
       
      // prodsum = prodsum - res_hi
      _fmpz_poly_sub(prodsum, prodsum, temp);
      
      // res_mid += prodsum
      prodsum->length = FLINT_MIN(prodsum->length, trunc - a1->length);
      temp2->coeffs = res->coeffs + a1->length*(res->limbs+1);
      temp2->length = prodsum->length;
      _fmpz_poly_add(temp2, temp2, prodsum);
      
      res->length = FLINT_MIN(a->length + b->length - 1, trunc);
      
   } else 
   {
      fmpz_poly_t scratch2, temp1; 

      while ((1<<l2)<a1->length) l2++;
      if ((1<<l2) < sa) a1->length = (1<<l2);
      a2->length = sa - a1->length;
      a1->coeffs = a->coeffs;
      a2->coeffs = a->coeffs+a1->length*(a->limbs+1);
      
      /* 
         The first a1->length + b->length - 1 coefficients will be written to directly, 
         so we need to clean the remaining coefficients
      */
      if (trunc > a1->length + sb - 1)
         for (unsigned long i = a1->length + sb - 1; i < FLINT_MIN(sa + sb - 1, trunc); i++)
            res->coeffs[i*(res->limbs+1)] = 0;
            
      temp->coeffs = b->coeffs;
      temp->limbs = b->limbs;
      temp->length = sb;
      
      // res_lo = a1*b
      if (sb <= a1->length) __fmpz_poly_karatrunc_recursive(res, a1, temp, scratch, scratchb, crossover, trunc);
      else __fmpz_poly_karatrunc_recursive(res, temp, a1, scratch, scratchb, crossover, trunc);
      
      //temp2 = a2*b
      temp2->coeffs = scratch->coeffs;
      temp2->length = FLINT_MIN(a2->length + sb - 1, trunc - a1->length);
      temp2->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs+temp2->length*(scratch->limbs+1);
      scratch2->limbs = scratch->limbs;
      if (sb <= a2->length) __fmpz_poly_karatrunc_recursive(temp2, a2, temp, scratch2, scratchb, crossover, trunc - a1->length);
      else __fmpz_poly_karatrunc_recursive(temp2, temp, a2, scratch2, scratchb, crossover, trunc - a1->length);
      
      // res_mid += 2temp
      temp1->coeffs = res->coeffs+a1->length*(res->limbs+1);
      temp1->length = FLINT_MIN(temp2->length, trunc - a1->length);
      temp1->limbs = res->limbs;
      _fmpz_poly_add(temp1, temp1, temp2);
      
      res->length = FLINT_MIN(sa + sb - 1, trunc);
   } 
}

void _fmpz_poly_mul_karatsuba_trunc(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2, unsigned long trunc)
{
   if (trunc == 0) 
   {
      output->length = 0;
      return;
   }
   
   unsigned long limbs = output->limbs;
   unsigned long log_length = 0;
   unsigned long crossover;
   
   fmpz_poly_t scratch, scratchb, temp;
   scratch->coeffs = (mp_limb_t *) flint_stack_alloc(6*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs+1;
   scratchb->limbs = FLINT_MAX(input1->limbs,input2->limbs)+1;
   scratchb->coeffs = (mp_limb_t *) flint_stack_alloc(6*FLINT_MAX(input1->length,input2->length)*(scratchb->limbs+1));
   
   if (_fmpz_poly_max_limbs(input1) + _fmpz_poly_max_limbs(input2) >= 19) crossover = 0;
   else crossover = 19 - _fmpz_poly_max_limbs(input1) - _fmpz_poly_max_limbs(input2);
   
   if (input1->length >= input2->length)
       __fmpz_poly_karatrunc_recursive(output, input1, input2, scratch, scratchb, crossover, trunc);
   else
       __fmpz_poly_karatrunc_recursive(output, input2, input1, scratch, scratchb, crossover, trunc);
   
   flint_stack_release(); flint_stack_release();
}

void __fmpz_poly_karatrunc_left_recursive(fmpz_poly_t res, fmpz_poly_t a, fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, unsigned long crossover, unsigned long trunc)
{
   fmpz_poly_t temp, temp2;
   
   long non_zero = a->length + b->length - trunc - 1;
   if (non_zero <= 0)
   {
      for (unsigned long i = 0; i < a->length + b->length - 1; i++)
      {
         res->coeffs[i*(res->limbs+1)] = 0;
      }
      return;
   }
   
   if ((a->length <= 1) || (b->length <= 1) || (non_zero == 1)) 
   {
      _fmpz_poly_mul_naive_trunc_left(res, a, b, trunc);
      
      return;
   }
   
   if ((a->length == 2 && b->length == 2) && (crossover < 4) && (!trunc)) 
   {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      __fmpz_poly_mul_coeffs2(res->coeffs, a->coeffs, b->coeffs); 
         
      __fmpz_poly_add_coeffs(scratchb->coeffs, a->coeffs, a->coeffs+asize);
      __fmpz_poly_add_coeffs(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
         
      __fmpz_poly_mul_coeffs(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
         
      __fmpz_poly_mul_coeffs(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
         
      __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
      __fmpz_poly_sub_coeffs(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
      
            
      res->length = a->length + b->length - 1;
      
      return;
   }
   
   if ((a->length+b->length <= crossover) || ((a->length == 2) && (b->length == 2)))
   {
      _fmpz_poly_mul_naive_trunc_left(res, a, b, trunc);
      
      return;
   }   
        
   fmpz_poly_t a1, a2, b1, b2;
      
   unsigned long l2 = 0;
   
   unsigned long old_length;
     
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
      
      /* 
         from 0 for 2 * a1->length - 1, from 2 * a1->length for a2->length + b2->length - 1
         will be written directly to, so we need to clean the coefficient in between
      */
      res->coeffs[((a1->length<<1)-1)*(res->limbs+1)] = 0;
  
      fmpz_poly_t asum, bsum, prodsum, scratch2, scratch3;
     
      asum->length = a1->length;
      asum->coeffs = scratchb->coeffs;
      asum->limbs = scratchb->limbs;
      bsum->length = a1->length;
      bsum->coeffs = scratchb->coeffs + a1->length*(scratchb->limbs+1);
      bsum->limbs = scratchb->limbs;
      prodsum->length = (a1->length<<1)-1;
      prodsum->coeffs = scratch->coeffs;
      prodsum->limbs = scratch->limbs;
      
      /*
         (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
      */
      
      // res_lo = a1*b1
      if (trunc > a1->length) __fmpz_poly_karatrunc_left_recursive(res, a1, b1, scratch, scratchb, crossover, trunc - a1->length);
      else __fmpz_poly_karatrunc_left_recursive(res, a1, b1, scratch, scratchb, crossover, 0);
      
      // res_hi = a2*b2
      temp->coeffs = res->coeffs+(a1->length<<1)*(res->limbs+1);
      temp->limbs = res->limbs;
      if (trunc > a1->length*2) __fmpz_poly_karatrunc_left_recursive(temp, a2, b2, scratch, scratchb, crossover, trunc - a1->length*2);
      else __fmpz_poly_karatrunc_left_recursive(temp, a2, b2, scratch, scratchb, crossover, 0);
      
      if (trunc < 3*a1->length - 1)
      {
         // asum = a1+a2
         _fmpz_poly_add(asum, a1, a2);
         // bsum = b1+b2
         _fmpz_poly_add(bsum, b1, b2);
         // prodsum = asum*bsum
         scratch3->coeffs = scratchb->coeffs+(a1->length<<1)*(scratchb->limbs+1);
         scratch3->limbs = scratchb->limbs;
      
         scratch2->limbs = scratch->limbs;
         scratch2->coeffs = scratch->coeffs+((a1->length<<1)-1)*(scratch->limbs+1);
         if (trunc > a1->length) __fmpz_poly_karatrunc_left_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover, trunc - a1->length);
         else __fmpz_poly_karatrunc_left_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover, 0);
         
         // prodsum = prodsum - res_lo
         temp->coeffs = res->coeffs;
         temp->length = (a1->length<<1)-1;
         _fmpz_poly_sub(prodsum, prodsum, temp);
       
         // prodsum = prodsum - res_hi
         temp->coeffs = res->coeffs + (a1->length<<1)*(res->limbs+1);
         temp->length = a2->length+b2->length-1;
         _fmpz_poly_sub(prodsum, prodsum, temp);
      
         // res_mid += prodsum
         temp->coeffs = res->coeffs + a1->length*(res->limbs+1);
         temp->length = prodsum->length;
         _fmpz_poly_add(temp, temp, prodsum);
      
      }
      
      res->length = a->length + b->length - 1;
      
   } else 
   {
      fmpz_poly_t scratch2, temp1; 

      while ((1<<l2)<a1->length) l2++;
      if ((1<<l2) < a->length) a1->length = (1<<l2);
      a2->length = a->length-a1->length;
      a1->coeffs = a->coeffs;
      a2->coeffs = a->coeffs+a1->length*(a->limbs+1);

      /* 
         The first a1->length + b->length - 1 coefficients will be written to directly, 
         so we need to clean the remaining coefficients
      */
      for (unsigned long i = a1->length + b->length - 1; i < a->length + b->length - 1; i++)
         res->coeffs[i*(res->limbs+1)] = 0;
      
      // res_lo = a1*b
      if (trunc < a1->length + b->length - 1) __fmpz_poly_karatrunc_left_recursive(res, a1, b, scratch, scratchb, crossover, trunc);
      
      //temp = a2*b
      temp->coeffs = scratch->coeffs;
      temp->length = a2->length + b->length - 1;
      temp->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs+temp->length*(scratch->limbs+1);
      scratch2->limbs = scratch->limbs;
      if (trunc > a1->length)
      {
         if (b->length <= a2->length) __fmpz_poly_karatrunc_left_recursive(temp, a2, b, scratch2, scratchb, crossover, trunc - a1->length);
         else __fmpz_poly_karatrunc_left_recursive(temp, b, a2, scratch2, scratchb, crossover, trunc - a1->length);
      } else
      {
         if (b->length <= a2->length) __fmpz_poly_karatrunc_left_recursive(temp, a2, b, scratch2, scratchb, crossover, 0);
         else __fmpz_poly_karatrunc_left_recursive(temp, b, a2, scratch2, scratchb, crossover, 0);
      }
      
      // res_mid += temp
      temp1->coeffs = res->coeffs+a1->length*(res->limbs+1);
      temp1->length = temp->length;
      temp1->limbs = res->limbs;
      _fmpz_poly_add(temp1,temp1,temp);
  
      res->length = a->length + b->length - 1;
   } 
   
   for (unsigned long i = 0; i < trunc; i++)
   {
      res->coeffs[i*(res->limbs+1)] = 0;
   }     
}

void _fmpz_poly_mul_karatsuba_trunc_left(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2, unsigned long trunc)
{
   unsigned long limbs = output->limbs;
   unsigned long log_length = 0;
   unsigned long crossover;
   
   fmpz_poly_t scratch, scratchb, temp;
   scratch->coeffs = (mp_limb_t *) flint_stack_alloc(5*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs + 1;
   scratchb->limbs = FLINT_MAX(input1->limbs, input2->limbs) + 1;
   scratchb->coeffs = (mp_limb_t *) flint_stack_alloc(5*FLINT_MAX(input1->length, input2->length)*(scratchb->limbs+1));
   
   if (_fmpz_poly_max_limbs(input1) + _fmpz_poly_max_limbs(input2) >= 19) crossover = 0;
   else crossover = 19 - _fmpz_poly_max_limbs(input1) - _fmpz_poly_max_limbs(input2);
   
   if (input1->length >= input2->length)
       __fmpz_poly_karatrunc_left_recursive(output, input1, input2, scratch, scratchb, crossover, trunc);
   else
       __fmpz_poly_karatrunc_left_recursive(output, input2, input1, scratch, scratchb, crossover, trunc);
   
   flint_stack_release(); 
   flint_stack_release();
}

void _fmpz_poly_mul_KS(fmpz_poly_t output, fmpz_poly_p input1, fmpz_poly_p input2)
{
   long sign1 = 1L;
   long sign2 = 1L;
   
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   unsigned long final_length = length1 + length2 - 1;
   
   while ((input1->coeffs[(length1-1)*(input1->limbs+1)] == 0) && (length1)) length1--;
   while ((input2->coeffs[(length2-1)*(input2->limbs+1)] == 0) && (length2)) length2--;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }

   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP(input1, input2);
   }
   
   if ((long) input1->coeffs[(length1-1)*(input1->limbs+1)] < 0)
   {
      _fmpz_poly_neg(input1, input1);
      sign1 = -1L;
   }
   
   if (input1 != input2)
   {
      if ((long) input2->coeffs[(length2-1)*(input2->limbs+1)] < 0)
      {
         _fmpz_poly_neg(input2, input2);
         sign2 = -1L;
      }
   } else sign2 = sign1;
   
   long bits1, bits2;
   int bitpack = 0;
   
   bits1 = _fmpz_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = length2;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   unsigned long limbs = (bits-1)/FLINT_BITS + 1;
   
   if ((bits < 64) && (input1->limbs == 1) && (input2->limbs == 1) && (output->limbs == 1)) bitpack = 1;
   
   unsigned long bytes = ((bits-1)>>3)+1;
   
   ZmodF_poly_t poly1, poly2, poly3;
   if (bitpack)
   {
      ZmodF_poly_stack_init(poly1, 0, (bits*length1-1)/FLINT_BITS+1, 0);
      if (input1 != input2)
         ZmodF_poly_stack_init(poly2, 0, (bits*length2-1)/FLINT_BITS+1, 0);

      if (sign) bits = -1L*bits;
      if (input1 != input2)
         ZmodF_poly_bit_pack_mpn(poly2, input2, length2, bits, length2);
      ZmodF_poly_bit_pack_mpn(poly1, input1, length1, bits, length1);

      bits=ABS(bits);
   } else
   {
      ZmodF_poly_stack_init(poly1, 0, ((bytes*length1-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);
      if (input1 != input2)
         ZmodF_poly_stack_init(poly2, 0, ((bytes*length2-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);

      ZmodF_poly_byte_pack_mpn(poly1, input1, length1, bytes, length1);
      if (input1 != input2)
         ZmodF_poly_byte_pack_mpn(poly2, input2, length2, bytes, length2);
   }
   
   if (input1 == input2)
   {
      poly2->coeffs = poly1->coeffs;
      poly2->n = poly1->n;
   }
   
   ZmodF_poly_stack_init(poly3, 0, poly1->n + poly2->n, 0);
           
   mp_limb_t msl = Z_mpn_mul(poly3->coeffs[0], poly1->coeffs[0], poly1->n, poly2->coeffs[0], poly2->n);
   
   poly3->coeffs[0][poly1->n+poly2->n-1] = msl;
   poly3->coeffs[0][poly1->n+poly2->n] = 0;
   poly3->length = 1;
   
   output->length = length1+length2-1;
  
   for (unsigned long i = 0; i < output->length; i++)
      output->coeffs[i*(output->limbs+1)] = 0;
      
   if (bitpack)
   {
      if (sign) ZmodF_poly_bit_unpack_mpn(output, poly3, length1+length2-1, bits);  
      else ZmodF_poly_bit_unpack_unsigned_mpn(output, poly3, length1+length2-1, bits);  
   } else
   {
      if (sign) ZmodF_poly_byte_unpack_mpn(output, poly3->coeffs[0], length1+length2-1, bytes);        
      else ZmodF_poly_byte_unpack_unsigned_mpn(output, poly3->coeffs[0], length1+length2-1, bytes);  
   }
   
   ZmodF_poly_stack_clear(poly3);
   if (input1 != input2)
      ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
     
   if ((long) (sign1 ^ sign2) < 0) _fmpz_poly_neg(output, output);
   
   if (sign1 < 0) _fmpz_poly_neg(input1, input1);
   if ((sign2 < 0) && (input1 != input2)) _fmpz_poly_neg(input2, input2);
   for (unsigned long i = output->length; i < input1->length + input2->length - 1; i++)
   {
      output->coeffs[i*(output->limbs+1)] = 0;
   }
   output->length = final_length;
}

void _fmpz_poly_mul_KS_trunc(fmpz_poly_t output, fmpz_poly_p input1, 
                                        fmpz_poly_p input2, unsigned long trunc)
{
   long sign1 = 1L;
   long sign2 = 1L;
   
   unsigned long length1 = FLINT_MIN(input1->length, trunc);
   unsigned long length2 = FLINT_MIN(input2->length, trunc);
   
   while ((input1->coeffs[(length1-1)*(input1->limbs+1)] == 0) && (length1)) length1--;
   while ((input2->coeffs[(length2-1)*(input2->limbs+1)] == 0) && (length2)) length2--;
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP(input1, input2);
   }
   
   if ((long) input1->coeffs[(length1-1)*(input1->limbs+1)] < 0)
   {
      _fmpz_poly_neg(input1, input1);
      sign1 = -1L;
   }
   
   if (input1 != input2)
   {
      if ((long) input2->coeffs[(length2-1)*(input2->limbs+1)] < 0)
      {
         _fmpz_poly_neg(input2, input2);
         sign2 = -1L;
      }
   } else sign2 = sign1;
   
   long bits1, bits2;
   int bitpack = 0;
   
   bits1 = _fmpz_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = length2;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   unsigned long limbs = (bits-1)/FLINT_BITS + 1;
   
   if ((bits < FLINT_BITS) && (input1->limbs == 1) && (input2->limbs == 1) && (output->limbs == 1)) bitpack = 1;
   
   unsigned long bytes = ((bits-1)>>3)+1;
   
   ZmodF_poly_t poly1, poly2, poly3;
   if (bitpack)
   {
      ZmodF_poly_stack_init(poly1, 0, (bits*length1-1)/FLINT_BITS+1, 0);
      if (input1 != input2)
         ZmodF_poly_stack_init(poly2, 0, (bits*length2-1)/FLINT_BITS+1, 0);

      if (sign) bits = -bits;
      if (input1 != input2)
         ZmodF_poly_bit_pack_mpn(poly2, input2, length2, bits, length2);
      ZmodF_poly_bit_pack_mpn(poly1, input1, length1, bits, length1);

      bits = ABS(bits);
   } else
   {
      ZmodF_poly_stack_init(poly1, 0, ((bytes*length1-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);
      if (input1 != input2)
         ZmodF_poly_stack_init(poly2, 0, ((bytes*length2-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);

      ZmodF_poly_byte_pack_mpn(poly1, input1, length1, bytes, length1);
      if (input1 != input2)
         ZmodF_poly_byte_pack_mpn(poly2, input2, length2, bytes, length2);
   }
   
   if (input1 == input2)
   {
      poly2->coeffs = poly1->coeffs;
      poly2->n = poly1->n;
   }
   
   ZmodF_poly_stack_init(poly3, 0, poly1->n + poly2->n, 0);
           
   output->length = FLINT_MIN(length1+length2-1, trunc);

   mp_limb_t msl;
   
   if (bitpack)
   {
      msl = Z_mpn_mul_trunc(poly3->coeffs[0], poly1->coeffs[0], poly1->n, poly2->coeffs[0], poly2->n, (output->length*bits-1)/FLINT_BITS+1);
   } else
   {
      msl = Z_mpn_mul_trunc(poly3->coeffs[0], poly1->coeffs[0], poly1->n, poly2->coeffs[0], poly2->n, ((output->length*bytes-1)>>FLINT_LG_BYTES_PER_LIMB) + 1);
   }
   
   poly3->coeffs[0][poly1->n+poly2->n-1] = msl;
   poly3->coeffs[0][poly1->n+poly2->n] = 0;
   poly3->length = 1;
      
   for (unsigned long i = 0; i < trunc; i++)
      output->coeffs[i*(output->limbs+1)] = 0;
      
   if (bitpack)
   {
      if (sign) ZmodF_poly_bit_unpack_mpn(output, poly3, output->length, bits);  
      else ZmodF_poly_bit_unpack_unsigned_mpn(output, poly3, output->length, bits);  
   } else
   {
      if (sign) ZmodF_poly_byte_unpack_mpn(output, poly3->coeffs[0], output->length, bytes);        
      else ZmodF_poly_byte_unpack_unsigned_mpn(output, poly3->coeffs[0], output->length, bytes);  
   }
   
   ZmodF_poly_stack_clear(poly3);
   if (input1 != input2)
      ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
     
   if ((long) (sign1 ^ sign2) < 0) _fmpz_poly_neg(output, output);
   
   if (sign1 < 0) _fmpz_poly_neg(input1, input1);
   if ((sign2 < 0) && (input1 != input2)) _fmpz_poly_neg(input2, input2);
}


void _fmpz_poly_mul_SS(fmpz_poly_t output, fmpz_poly_p input1, fmpz_poly_p input2)
{
   unsigned long length1 = input1->length;
   unsigned long length2 = input2->length;
   
   while ((input1->coeffs[(length1-1)*(input1->limbs+1)] == 0) && (length1)) length1--;
   while ((input2->coeffs[(length2-1)*(input2->limbs+1)] == 0) && (length2)) length2--;
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP(input1, input2);
   }
   
   unsigned long size1 = input1->limbs;
   unsigned long size2 = input2->limbs;
   
   unsigned long log_length = 0;
   while ((1<<log_length) < length1) log_length++;
   unsigned long log_length2 = 0;
   while ((1<<log_length2) < length2) log_length2++;
   
   /* Start with an upper bound on the number of bits needed */
   
   unsigned long output_bits = FLINT_BITS * (size1 + size2) + log_length2 + 2;

   if (output_bits <= length1)
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   else 
      output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
      
   unsigned long n = (output_bits - 1) / FLINT_BITS + 1;
   
   ZmodF_poly_t poly1, poly2, res;
   long bits1, bits2;
   unsigned long sign = 0;
   
   ZmodF_poly_stack_init(poly1, log_length + 1, n, 1);
   ZmodF_poly_stack_init(poly2, log_length + 1, n, 1);
   ZmodF_poly_stack_init(res, log_length + 1, n, 1);
   
   bits1 = ZmodF_poly_convert_in_mpn(poly1, input1, length1);
   bits2 = ZmodF_poly_convert_in_mpn(poly2, input2, length2);
   
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
      
   n = (output_bits - 1) / FLINT_BITS + 1;
   
   ZmodF_poly_decrease_n(poly1, n);
   ZmodF_poly_decrease_n(poly2, n);
   ZmodF_poly_decrease_n(res, n);
                    
   ZmodF_poly_convolution(res, poly1, poly2);
   ZmodF_poly_normalise(res);
          
   output->length = length1 + length2 - 1;
   
   ZmodF_poly_convert_out_mpn(output, res, sign);
   
   ZmodF_poly_stack_clear(res);
   ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
}

void _fmpz_poly_mul_SS_trunc(fmpz_poly_t output, fmpz_poly_p input1, 
                                        fmpz_poly_p input2, unsigned long trunc)
{
   unsigned long length1 = FLINT_MIN(input1->length, trunc);
   unsigned long length2 = FLINT_MIN(input2->length, trunc);
   
   while ((input1->coeffs[(length1-1)*(input1->limbs+1)] == 0) && (length1)) length1--;
   while ((input2->coeffs[(length2-1)*(input2->limbs+1)] == 0) && (length2)) length2--;
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      SWAP(input1, input2);
   }
   
   unsigned long size1 = input1->limbs;
   unsigned long size2 = input2->limbs;
   
   unsigned long log_length = 0;
   while ((1<<log_length) < length1) log_length++;
   unsigned long log_length2 = 0;
   while ((1<<log_length2) < length2) log_length2++;
   
   /* Start with an upper bound on the number of bits needed */
   
   unsigned long output_bits = FLINT_BITS * (size1 + size2) + log_length2 + 2;

   if (output_bits <= length1)
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
   else 
      output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
      
   unsigned long n = (output_bits - 1) / FLINT_BITS + 1;
   
   ZmodF_poly_t poly1, poly2, res;
   long bits1, bits2;
   unsigned long sign = 0;
   
   ZmodF_poly_stack_init(poly1, log_length + 1, n, 1);
   ZmodF_poly_stack_init(poly2, log_length + 1, n, 1);
   ZmodF_poly_stack_init(res, log_length + 1, n, 1);
   
   bits1 = ZmodF_poly_convert_in_mpn(poly1, input1, length1);
   bits2 = ZmodF_poly_convert_in_mpn(poly2, input2, length2);
   
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
      
   n = (output_bits - 1) / FLINT_BITS + 1;
   
   ZmodF_poly_decrease_n(poly1, n);
   ZmodF_poly_decrease_n(poly2, n);
   ZmodF_poly_decrease_n(res, n);
                    
   ZmodF_poly_convolution_trunc(res, poly1, poly2, trunc);
   res->length = FLINT_MIN(res->length, trunc);
   ZmodF_poly_normalise(res);
          
   output->length = FLINT_MIN(length1 + length2 - 1, trunc);
   
   ZmodF_poly_convert_out_mpn(output, res, sign);
   
   ZmodF_poly_stack_clear(res);
   ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
}

/*
   A truncating polynomial multiplication.
   The number of terms require, _trunc_ can be any value, but the function is
   tuned for truncation to length n where both inputs have length approximately n.
*/

void _fmpz_poly_mul_trunc_n(fmpz_poly_t output, fmpz_poly_t input1, 
                                fmpz_poly_t input2, unsigned long trunc)
{
   if ((input1->length == 0) || (input2->length == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }

   if ((input1->length <= 3) && (input2->length <= 3)) 
   {
      _fmpz_poly_mul_karatsuba_trunc(output, input1, input2, trunc);
      return;
   }
   
   unsigned long bits1 = _fmpz_poly_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
   bits1 = ABS(bits1);
   bits2 = ABS(bits2);
   
   if ((bits1 + bits2 >= 64) && (input1->length + input2->length <= 10)) 
   {
      _fmpz_poly_mul_karatsuba_trunc(output, input1, input2, trunc);
      return;
   }
   
   if ((bits1 + bits2 >= 370) && (input1->length + input2->length <= 32)) 
   {
      _fmpz_poly_mul_karatsuba_trunc(output, input1, input2, trunc);
      return;
   }   
   
   if (bits1 + bits2 < 512)
   {
      _fmpz_poly_mul_KS_trunc(output, input1, input2, trunc);
      return;
   } 
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _fmpz_poly_mul_SS_trunc(output, input1, input2, trunc);
      return;
   } 
   
   _fmpz_poly_mul_KS_trunc(output, input1, input2, trunc);     
}

void _fmpz_poly_mul(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
{
   if ((input1->length == 0) || (input2->length == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }

   if ((input1->length <= 2) && (input2->length <= 2)) 
   {
      _fmpz_poly_mul_karatsuba(output, input1, input2);
      return;
   }
   
   if ((input1->limbs <= 256/FLINT_BITS) && (input1->limbs >= 200/FLINT_BITS) && (input1->length == 256)) 
   {
      _fmpz_poly_mul_SS(output, input1, input2);
      return;
   } 
   
   if (input1->limbs + input2->limbs <= 512/FLINT_BITS)
   {
      _fmpz_poly_mul_KS(output, input1, input2);
      return;
   }
   
   if (input1->length + input2->length <= 32) 
   {
      _fmpz_poly_mul_karatsuba(output, input1, input2);
      return;
   }
   
   unsigned long bits1 = _fmpz_poly_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
   bits1 = ABS(bits1);
   bits2 = ABS(bits2);
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _fmpz_poly_mul_SS(output, input1, input2);
      return;
   } 
   
   _fmpz_poly_mul_KS(output, input1, input2);     
}

/*
   A truncating polynomial multiplication which ignores the first trunc coeffs of
   the output (which can end up being anything - often zero).
   The number of zero terms, _trunc_ can be any value, but the function is
   tuned for truncation of length n-1 terms, where both inputs have length approximately n.
*/

void _fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, fmpz_poly_t input1, 
                                fmpz_poly_t input2, unsigned long trunc)
{
   if ((input1->length == 0) || (input2->length == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }

   if ((input1->length <= 3) && (input2->length <= 3)) 
   {
      _fmpz_poly_mul_karatsuba_trunc_left(output, input1, input2, trunc);
      return;
   }
   
   unsigned long bits1 = _fmpz_poly_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
   bits1 = ABS(bits1);
   bits2 = ABS(bits2);
   
   if ((bits1 + bits2 >= 64) && (input1->length + input2->length <= 10)) 
   {
      _fmpz_poly_mul_karatsuba_trunc_left(output, input1, input2, trunc);
      return;
   }
   
   if ((bits1 + bits2 >= 370) && (input1->length + input2->length <= 32)) 
   {
      _fmpz_poly_mul_karatsuba_trunc_left(output, input1, input2, trunc);
      return;
   }   
   
   if (bits1 + bits2 < 512)
   {
      _fmpz_poly_mul_KS(output, input1, input2);
      return;
   } 
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _fmpz_poly_mul_SS(output, input1, input2);
      return;
   } 
   
   _fmpz_poly_mul_KS(output, input1, input2);     
}


/*
   Scalar multiplication of a polynomial by a scalar
*/

void _fmpz_poly_scalar_mul(fmpz_poly_t output, fmpz_poly_t poly, mp_limb_t * x)
{
   NORM(x);
   
   unsigned long limbs1 = ABS(x[0]);
   unsigned long limbs2 = poly->limbs;
   unsigned long total_limbs;
   unsigned long msl;
   unsigned long limbs_out = output->limbs+1;
   mp_limb_t * coeffs_out = output->coeffs;
   mp_limb_t * coeffs2 = poly->coeffs;
   long sign1 = x[0];
   
   if (limbs1 == 1)
   {
      for (long i = 0; i < poly->length; i++)
      {
          total_limbs = 1 + ABS(coeffs2[i*(limbs2+1)]);
          if (total_limbs != 1)
          {
             msl = mpn_mul_1(coeffs_out + i*limbs_out + 1, coeffs2 + i*(limbs2+1) + 1, ABS(coeffs2[i*(limbs2+1)]), x[1]);
             if (msl) coeffs_out[i*limbs_out+ABS(coeffs2[i*(limbs2+1)])+1] = msl;
             if (((long) coeffs2[i*(limbs2+1)] ^ sign1) < 0) coeffs_out[i*limbs_out] = -total_limbs + (msl == 0L);
             else coeffs_out[i*limbs_out] = total_limbs - (msl == 0L);
          } else coeffs_out[i*limbs_out] = 0;
      }
   } else if (limbs1 + limbs2 > 1000)//FLINT_FFT_LIMBS_CROSSOVER*2)
   {
      Z_mpn_precomp_t precomp;
   
      Z_mpn_mul_precomp_init(precomp, x+1, limbs1, limbs2);   
      
      for (long i = 0; i < poly->length; i++)
      {
          total_limbs = limbs1 + ABS(coeffs2[i*(limbs2+1)]);
          if (total_limbs != limbs1)
          {
             msl = Z_mpn_mul_precomp(coeffs_out + i*limbs_out + 1, coeffs2 + i*(limbs2+1) + 1, ABS(coeffs2[i*(limbs2+1)]), precomp);
             if (((long) coeffs2[i*(limbs2+1)] ^ sign1) < 0) coeffs_out[i*limbs_out] = -total_limbs + (msl == 0L);
             else coeffs_out[i*limbs_out] = total_limbs - (msl == 0L);
          } else coeffs_out[i*limbs_out] = 0;
      }
      Z_mpn_mul_precomp_clear(precomp);
   } else
   {
      if (poly != output)
      {
         for (long i = 0; i < poly->length - 1; i++)
         {
            __fmpz_poly_mul_coeffs2(coeffs_out + i*limbs_out, coeffs2 + i*(limbs2+1), x);
         }
         __fmpz_poly_mul_coeffs(coeffs_out + (poly->length - 1)*limbs_out, coeffs2 + (poly->length - 1)*(limbs2+1), x);
      } else
      {
         for (long i = 0; i < poly->length; i++)
         {
            __fmpz_poly_mul_coeffs(coeffs_out + i*limbs_out, coeffs2 + i*(limbs2+1), x);
         }
      }
   } 
   output->length = poly->length;
}

/* 
   Set n of the coefficients of poly to zero.
*/

void _fmpz_poly_zero_coeffs(fmpz_poly_t poly, unsigned long n)
{
   unsigned long size = poly->limbs+1;
   mp_limb_t * coeff = poly->coeffs;
   for (unsigned long i = 0; i < n; i++)
   {
      coeff[0] = 0;
      coeff+=size;
   }
}

/****************************************************************************

   fmpz_poly_* layer

****************************************************************************/

/* 
   Create a polynomial of length zero with zero allocated coefficients
*/

void fmpz_poly_init(fmpz_poly_t poly)
{
   poly->coeffs = NULL;
   
   poly->alloc = 0;
   poly->length = 0;
   poly->limbs = 0;
}

/* 
   Create a polynomial of length zero with "alloc" allocated coefficients
   each with enough space for "limbs" limbs
*/

void fmpz_poly_init2(fmpz_poly_t poly, unsigned long alloc, unsigned long limbs)
{
   if ((alloc > 0) && (limbs > 0))
      poly->coeffs = (mp_limb_t *) flint_heap_alloc(alloc*(limbs+1));
   else poly->coeffs = NULL;
   
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

/* Shrink or expand a polynomial to "alloc" coefficients */

void fmpz_poly_realloc(fmpz_poly_t poly, unsigned long alloc)
{
   if (poly->limbs > 0)
   {
      if (alloc > 0)
      {
         poly->coeffs = (mp_limb_t*) flint_heap_realloc(poly->coeffs, alloc*(poly->limbs+1));
      } else
      {
         if (poly->coeffs) flint_heap_free(poly->coeffs);
      }   
      poly->alloc = alloc;
   
      // truncate actual data if necessary
      if (poly->length > alloc)
         poly->length = alloc;
   } else
   {
      poly->alloc = alloc;
   }
}

void fmpz_poly_fit_length(fmpz_poly_t poly, unsigned long alloc)
{
   if (alloc <= poly->alloc) return;

   if (alloc < 2*poly->alloc) alloc = 2*poly->alloc;
   
   fmpz_poly_realloc(poly, alloc);
}

void fmpz_poly_resize_limbs(fmpz_poly_t poly, unsigned long limbs)
{
   if (limbs > 0)
   {
      if (limbs == poly->limbs) return;
      
      unsigned long i;
      mp_limb_t * coeff_i;
      mp_limb_t * coeff_i_old = poly->coeffs;
      
      if (limbs < poly->limbs)
      {
         coeff_i = poly->coeffs + limbs+1;
         coeff_i_old += (poly->limbs+1);
         for (i = 1; i < poly->length; i++)
         {
            forward_copy_limbs(coeff_i, coeff_i_old, limbs+1);
            FLINT_ASSERT(ABS(coeff_i[0]) > limbs); 
            coeff_i += (limbs+1);
            coeff_i_old += (poly->limbs+1);
         } 
      } else
      {
         mp_limb_t * temp_coeffs = (mp_limb_t*) flint_heap_alloc(poly->alloc*(limbs+1));
         coeff_i = temp_coeffs;
         for (i = 0; i < poly->length; i++)
         {
            copy_limbs(coeff_i, coeff_i_old, poly->limbs+1);
            coeff_i += (limbs+1);
            coeff_i_old += (poly->limbs+1);
         } 
         flint_heap_free(poly->coeffs);
         poly->coeffs = temp_coeffs;
      }
      for ( ; i < poly->alloc; i++)
      {
         coeff_i[0] = 0;
         coeff_i += (limbs+1);
      } 
      poly->limbs = limbs;
   } else
   {
      if (poly->coeffs) flint_heap_free(poly->coeffs);
      poly->length = 0;
      poly->limbs = 0;
   }
}

void fmpz_poly_clear(fmpz_poly_t poly)
{
   if (poly->coeffs) flint_heap_free(poly->coeffs);
}

long fmpz_poly_degree(fmpz_poly_t poly)
{
   _fmpz_poly_normalise(poly);
   return poly->length - 1;
}

unsigned long fmpz_poly_length(fmpz_poly_t poly)
{
   _fmpz_poly_normalise(poly);
   return poly->length;
}

void fmpz_poly_set_length(fmpz_poly_t poly, unsigned long length)
{
   fmpz_poly_fit_length(poly, length);

   for (unsigned long i = poly->length; i < length; i++)
      poly->coeffs[i*(poly->limbs+1)] = 0;

   poly->length = length;
}

void fmpz_poly_get_coeff_mpz(mpz_t x, fmpz_poly_t poly, unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(x, 0);
   else
      _fmpz_poly_get_coeff_mpz(x, poly, n);
}

void fmpz_poly_mul(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2)
{
   if ((input1->length == 0) || (input2->length == 0))
   {
      fmpz_poly_fit_length(output, 1);
      fmpz_poly_fit_limbs(output, 1);
      _fmpz_poly_zero(output);
      return;
   }
   
   unsigned long limbs = input1->limbs + input2->limbs;
   unsigned long total_length = input1->length + input2->length - 1;
   
   long bits1, bits2;
      
   bits1 = _fmpz_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = (input1->length > input2->length) ? input2->length : input1->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   
   fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+1);
   fmpz_poly_fit_length(output, input1->length + input2->length - 1);
   
   _fmpz_poly_mul(output, input1, input2);
   output->length = total_length;
   
}

void fmpz_poly_mul_trunc_n(fmpz_poly_t output, fmpz_poly_t input1, 
                                          fmpz_poly_t input2, unsigned long trunc)
{
   long bits1, bits2;
      
   bits1 = _fmpz_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
     
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = (input1->length > input2->length) ? input2->length : input1->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   
   fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+1);
   fmpz_poly_fit_length(output, FLINT_MIN(input1->length + input2->length - 1, trunc));
   
   _fmpz_poly_mul_trunc_n(output, input1, input2, FLINT_MIN(input1->length + input2->length - 1, trunc));
   fmpz_poly_set_length(output, FLINT_MIN(input1->length + input2->length - 1, trunc));
}

void fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, fmpz_poly_t input1, fmpz_poly_t input2, unsigned long trunc)
{
   unsigned long limbs = input1->limbs + input2->limbs;
   
   long bits1, bits2;
      
   bits1 = _fmpz_poly_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_bits(input2);
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = (input1->length > input2->length) ? input2->length : input1->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   
   fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+1);
   fmpz_poly_fit_length(output, input1->length + input2->length - 1);
   
   _fmpz_poly_mul_trunc_left_n(output, input1, input2, trunc);
   fmpz_poly_set_length(output, input1->length + input2->length - 1);
}


void fmpz_poly_divrem_naive(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
   _fmpz_poly_normalise(A);
   _fmpz_poly_normalise(B);
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   long coeff = A->length-1;
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   mp_limb_t * coeffs_A = A->coeffs;
   mp_limb_t * coeffs_B = B->coeffs;
   mp_limb_t * coeff_i = coeffs_A + coeff*size_A;
   mp_limb_t * B_lead = coeffs_B + (B->length-1)*size_B; 
   mp_limb_t * coeff_Q;
   mp_limb_t * coeffs_R;

   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   mp_limb_t sign_quot;
   
   while (1)
   {
      if (coeff < (long) B->length - 1) break;
      NORM(coeff_i);
      if (ABS(coeff_i[0]) < size_B_lead)
      {
         coeff--;
         coeff_i -= size_A;
      } else if (ABS(coeff_i[0]) > size_B_lead) break;
      else if (mpn_cmp(coeff_i+1, B_lead+1, size_B_lead) >= 0) break;
      else 
      {
         coeff--;
         coeff_i -= size_A;         
      }    
   }
   
   mp_limb_t * rem = (mp_limb_t *) flint_heap_alloc(size_B_lead);
   
   fmpz_poly_fit_length(R, A->length);
   fmpz_poly_fit_limbs(R, A->limbs);
   R->length = A->length;
   _fmpz_poly_set(R, A);
   coeffs_R = R->coeffs;
   size_R = R->limbs+1;
      
   if (coeff >= (long) B->length - 1)
   {
      fmpz_poly_fit_length(Q, coeff-B->length+2);
      fmpz_poly_fit_limbs(Q, 1);
      Q->length = coeff-B->length+2;
      size_Q = Q->limbs+1;
   } else _fmpz_poly_zero(Q);
   
   while (coeff >= (long) B->length - 1)
   {
      coeff_Q = Q->coeffs+(coeff-B->length+1)*size_Q;
      while (1)
      {
         if (coeff < (long) B->length - 1) break;
         NORM(coeffs_R+coeff*size_R);
         if (ABS(coeffs_R[coeff*size_R]) < size_B_lead)
         {
            coeff_Q[0] = 0;
            coeff_Q -= size_Q;
            coeff--;
         } else if (ABS(coeffs_R[coeff*size_R]) > size_B_lead) break;
         else if (mpn_cmp(coeffs_R+coeff*size_R+1, B_lead+1, size_B_lead) >= 0) break;
         else 
         {
            coeff_Q[0] = 0;
            coeff_Q -= size_Q;
            coeff--;
         }    
      }
      
      if (coeff >= (long) B->length - 1)
      {
         limbs_Q = ABS(coeffs_R[coeff*size_R]) - size_B_lead + 1;
         fmpz_poly_fit_limbs(Q, limbs_Q);
         size_Q = Q->limbs+1;
         coeff_Q = Q->coeffs+(coeff - B->length+1)*size_Q;
         sign_quot = ABS(coeffs_R[coeff*size_R]) - size_B_lead + 1;
         if (((long) (sign_B_lead ^ coeffs_R[coeff*size_R])) < 0) 
         {
            mpn_tdiv_qr(coeff_Q+1, rem, 0, coeffs_R+coeff*size_R+1, ABS(coeffs_R[coeff*size_R]), B_lead+1, size_B_lead);
            coeff_Q[0] = -sign_quot;
            for (unsigned long i = 0; i < size_B_lead; i++)
            {
               if (rem[i])
               {
                  __fmpz_poly_sub_coeff_ui(coeff_Q,1);
                  break;
               }
            }
         }
         else 
         {
            mpn_tdiv_qr(coeff_Q+1, rem, 0, coeffs_R+coeff*size_R+1, ABS(coeffs_R[coeff*size_R]), B_lead+1, size_B_lead);
            coeff_Q[0] = sign_quot;
         }
         NORM(coeff_Q);
         
         fmpz_poly_init2(qB, B->length, B->limbs+ABS(coeff_Q[0]));
         _fmpz_poly_scalar_mul(qB, B, coeff_Q); 
         
         fmpz_poly_fit_limbs(R, qB->limbs+1);
         coeffs_R = R->coeffs;
         size_R = R->limbs+1;
      
         fmpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R+(coeff - B->length + 1)*size_R;
         R_sub->limbs = R->limbs;
         R_sub->length = B->length;
         _fmpz_poly_sub(R_sub, R_sub, qB);
         
         coeff--;
         fmpz_poly_clear(qB);
      }
   }
   
   _fmpz_poly_normalise(R);
   flint_heap_free(rem);
}

/* 
   Divides A by B and returns the quotient Q, but only the low half of the remainder R
*/

void fmpz_poly_divrem_naive_low(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
   _fmpz_poly_normalise(A);
   _fmpz_poly_normalise(B);
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   long coeff = A->length-1;
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   mp_limb_t * coeffs_A = A->coeffs;
   mp_limb_t * coeffs_B = B->coeffs;
   mp_limb_t * coeff_i = coeffs_A + coeff*size_A;
   mp_limb_t * B_lead = coeffs_B + (B->length-1)*size_B; 
   mp_limb_t * coeff_Q;
   mp_limb_t * coeffs_R;

   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   mp_limb_t sign_quot;
   
   while (1)
   {
      if (coeff < (long) B->length - 1) break;
      NORM(coeff_i);
      if (ABS(coeff_i[0]) < size_B_lead)
      {
         coeff--;
         coeff_i -= size_A;
      } else if (ABS(coeff_i[0]) > size_B_lead) break;
      else if (mpn_cmp(coeff_i+1, B_lead+1, size_B_lead) >= 0) break;
      else 
      {
         coeff--;
         coeff_i -= size_A;         
      }    
   }
   
   mp_limb_t * rem = (mp_limb_t *) flint_heap_alloc(size_B_lead);
   
   fmpz_poly_fit_length(R, A->length);
   fmpz_poly_fit_limbs(R, A->limbs);
   R->length = A->length;
   _fmpz_poly_set(R, A);
   coeffs_R = R->coeffs;
   size_R = R->limbs+1;
      
   if (coeff >= (long) B->length - 1)
   {
      fmpz_poly_fit_length(Q, coeff-B->length+2);
      fmpz_poly_fit_limbs(Q, 1);
      Q->length = coeff-B->length+2;
      size_Q = Q->limbs+1;
   } else _fmpz_poly_zero(Q);
   
   while (coeff >= (long) B->length - 1)
   {
      coeff_Q = Q->coeffs+(coeff-B->length+1)*size_Q;
      while (1)
      {
         if (coeff < (long) B->length - 1) break;
         NORM(coeffs_R+coeff*size_R);
         if (ABS(coeffs_R[coeff*size_R]) < size_B_lead)
         {
            coeff_Q[0] = 0;
            coeff_Q -= size_Q;
            coeff--;
         } else if (ABS(coeffs_R[coeff*size_R]) > size_B_lead) break;
         else if (mpn_cmp(coeffs_R+coeff*size_R+1, B_lead+1, size_B_lead) >= 0) break;
         else 
         {
            coeff_Q[0] = 0;
            coeff_Q -= size_Q;
            coeff--;
         }    
      }
      
      if (coeff >= (long) B->length - 1)
      {
         limbs_Q = ABS(coeffs_R[coeff*size_R]) - size_B_lead + 1;
         fmpz_poly_fit_limbs(Q, limbs_Q);
         size_Q = Q->limbs+1;
         coeff_Q = Q->coeffs+(coeff - B->length+1)*size_Q;
         sign_quot = ABS(coeffs_R[coeff*size_R]) - size_B_lead + 1;
         if (((long) (sign_B_lead ^ coeffs_R[coeff*size_R])) < 0) 
         {
            mpn_tdiv_qr(coeff_Q+1, rem, 0, coeffs_R+coeff*size_R+1, ABS(coeffs_R[coeff*size_R]), B_lead+1, size_B_lead);
            coeff_Q[0] = -sign_quot;
            for (unsigned long i = 0; i < size_B_lead; i++)
            {
               if (rem[i])
               {
                  __fmpz_poly_sub_coeff_ui(coeff_Q,1);
                  break;
               }
            }
         }
         else 
         {
            mpn_tdiv_qr(coeff_Q+1, rem, 0, coeffs_R+coeff*size_R+1, ABS(coeffs_R[coeff*size_R]), B_lead+1, size_B_lead);
            coeff_Q[0] = sign_quot;
         }
         NORM(coeff_Q);
         
         fmpz_poly_t temp;
         fmpz_poly_init2(qB, B->length, B->limbs+ABS(coeff_Q[0]));
         temp->coeffs = B->coeffs;
         temp->length = B->length - 1;
         temp->limbs = B->limbs;
         _fmpz_poly_scalar_mul(qB, temp, coeff_Q); 
         
         fmpz_poly_fit_limbs(R, qB->limbs+1);
         coeffs_R = R->coeffs;
         size_R = R->limbs+1;
      
         fmpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R+(coeff - B->length + 1)*size_R;
         R_sub->limbs = R->limbs;
         R_sub->length = B->length - 1;
         _fmpz_poly_sub(R_sub, R_sub, qB);
         coeffs_R[coeff*size_R] = 0;
         
         coeff--;
         fmpz_poly_clear(qB);
      }
   }
   
   _fmpz_poly_normalise(R);
   flint_heap_free(rem);
}

/*
   Divide the polynomial A by the polynomial B but do not compute the remainder
*/

void fmpz_poly_div_naive(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)
{
   fmpz_poly_t qB, R;
   
   fmpz_poly_init(R);
   
   _fmpz_poly_normalise(A);
   _fmpz_poly_normalise(B);
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   long coeff = A->length-1;
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   mp_limb_t * coeffs_A = A->coeffs;
   mp_limb_t * coeffs_B = B->coeffs;
   mp_limb_t * coeff_i = coeffs_A + coeff*size_A;
   mp_limb_t * B_lead = coeffs_B + (B->length-1)*size_B; 
   mp_limb_t * coeff_Q;
   mp_limb_t * coeffs_R;

   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   mp_limb_t sign_quot;
   
   // Find the first coefficient greater than B_lead
   while (1)
   {
      if (coeff < (long) B->length - 1) break;
      NORM(coeff_i);
      if (ABS(coeff_i[0]) < size_B_lead)
      {
         coeff--;
         coeff_i -= size_A;
      } else if (ABS(coeff_i[0]) > size_B_lead) break;
      else if (mpn_cmp(coeff_i+1, B_lead+1, size_B_lead) >= 0) break;
      else 
      {
         coeff--;
         coeff_i -= size_A;         
      }    
   }
   
   mp_limb_t * rem = (mp_limb_t *) flint_heap_alloc(size_B_lead);
   
   // Set R to A
   fmpz_poly_fit_length(R, A->length);
   fmpz_poly_fit_limbs(R, A->limbs);
   R->length = A->length;
   _fmpz_poly_set(R, A);
   coeffs_R = R->coeffs;
   size_R = R->limbs+1;
      
   // Set the quotient to zero if R is shorter than B
   if (coeff >= (long) B->length - 1)
   {
      fmpz_poly_fit_length(Q, coeff-B->length+2);
      fmpz_poly_fit_limbs(Q, 1);
      Q->length = coeff-B->length+2;
      size_Q = Q->limbs+1;
   } else _fmpz_poly_zero(Q);
   
   while (coeff >= (long) B->length - 1)
   {
      coeff_Q = Q->coeffs+(coeff-B->length+1)*size_Q;
      // Set quotient coefficients to 0 if the R coefficients are already smaller than B
      while (1)
      {
         if (coeff < (long) B->length - 1) break;
         NORM(coeffs_R+coeff*size_R);
         if (ABS(coeffs_R[coeff*size_R]) < size_B_lead)
         {
            coeff_Q[0] = 0;
            coeff_Q -= size_Q;
            coeff--;
         } else if (ABS(coeffs_R[coeff*size_R]) > size_B_lead) break;
         else if (mpn_cmp(coeffs_R+coeff*size_R+1, B_lead+1, size_B_lead) >= 0) break;
         else 
         {
            coeff_Q[0] = 0;
            coeff_Q -= size_Q;
            coeff--;
         }    
      }
      
      if (coeff >= (long) B->length - 1)
      {
         // else compute the quotient of the coefficient by B_lead
         limbs_Q = ABS(coeffs_R[coeff*size_R]) - size_B_lead + 1;
         fmpz_poly_fit_limbs(Q, limbs_Q);
         size_Q = Q->limbs+1;
         coeff_Q = Q->coeffs+(coeff - B->length+1)*size_Q;
         sign_quot = ABS(coeffs_R[coeff*size_R]) - size_B_lead + 1;
         if (((long) (sign_B_lead ^ coeffs_R[coeff*size_R])) < 0) 
         {
            mpn_tdiv_qr(coeff_Q+1, rem, 0, coeffs_R+coeff*size_R+1, ABS(coeffs_R[coeff*size_R]), B_lead+1, size_B_lead);
            coeff_Q[0] = -sign_quot;
            for (unsigned long i = 0; i < size_B_lead; i++)
            {
               if (rem[i])
               {
                  __fmpz_poly_sub_coeff_ui(coeff_Q,1);
                  break;
               }
            }
         }
         else 
         {
            mpn_tdiv_qr(coeff_Q+1, rem, 0, coeffs_R+coeff*size_R+1, ABS(coeffs_R[coeff*size_R]), B_lead+1, size_B_lead);
            coeff_Q[0] = sign_quot;
         }
         NORM(coeff_Q);
         
         if (coeff >= (long) B->length)
         {
            // Now multiply B by this new quotient coefficient and subtract from R
            fmpz_poly_t R_sub;
            unsigned long length = FLINT_MIN(coeff - B->length + 2, B->length);
            
            fmpz_poly_init2(qB, length, B->limbs+ABS(coeff_Q[0]));
            R_sub->coeffs = B->coeffs + (B->length - length)*(B->limbs + 1);
            R_sub->limbs = B->limbs;
            R_sub->length = length;
            _fmpz_poly_scalar_mul(qB, R_sub, coeff_Q); 
         
            fmpz_poly_fit_limbs(R, qB->limbs+1);
            coeffs_R = R->coeffs;
            size_R = R->limbs+1;
      
            R_sub->coeffs = coeffs_R+(coeff - length + 1)*size_R;
            R_sub->limbs = R->limbs;
            _fmpz_poly_sub(R_sub, R_sub, qB);

            fmpz_poly_clear(qB);
         }
         coeff--;
      }
   }
   
   fmpz_poly_clear(R);
   flint_heap_free(rem);
}

/* 
   Integer polynomial division using a Karatsuba-like algorithm.
   We assume for now that the length of B is 2*n = 2^l, and that the length of A is 4*n-1 
   Note BQ is not the remainder but B*Q, so the remainder R = A-BQ
*/
 
void fmpz_poly_div_karatsuba_recursive(fmpz_poly_t Q, fmpz_poly_t BQ, fmpz_poly_t A, fmpz_poly_t B)
{
   _fmpz_poly_normalise(A);
   _fmpz_poly_normalise(B);
         
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      _fmpz_poly_zero(BQ);

      return;
   }
   
   if ((B->length <= 16) || (A->length > 2*B->length - 1))
   {
      fmpz_poly_t Rb;
      fmpz_poly_init(Rb);
      fmpz_poly_divrem_naive(Q, Rb, A, B);
      fmpz_poly_fit_length(BQ, A->length);
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(A->limbs, Rb->limbs)+1);
      _fmpz_poly_sub(BQ, A, Rb);
      fmpz_poly_clear(Rb);
      
      return;
   }
   
   fmpz_poly_t d1, d2, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   d1->length = n1;
   d2->length = n2;
   d1->limbs = B->limbs;
   d2->limbs = B->limbs;
   d1->coeffs = B->coeffs + n2*(B->limbs+1);
   d2->coeffs = B->coeffs;
   
   if (A->length <= n1+2*n2-1)
   {
      temp->length = A->length - (n1+n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
      
      fmpz_poly_init(d1q1);
      fmpz_poly_div_karatsuba_recursive(Q, d1q1, p1, d1); //******************************
      _fmpz_poly_stack_clear(p1);

      _fmpz_poly_stack_init(d2q1, d2->length+Q->length-1, d2->limbs+Q->limbs+1); 
      _fmpz_poly_mul(d2q1, d2, Q);
      
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(d1q1->limbs, d2q1->limbs));
      fmpz_poly_fit_length(BQ, n2+d1q1->length);
      BQ->length = n2+d1q1->length;
      _fmpz_poly_zero_coeffs(BQ, n2);
   
      temp->length = d1q1->length;
      temp->limbs = BQ->limbs;
      temp->coeffs = BQ->coeffs + n2*(BQ->limbs+1);
      _fmpz_poly_set(temp, d1q1);  
      _fmpz_poly_add(BQ, BQ, d2q1);
   
      _fmpz_poly_stack_clear(d2q1);
      fmpz_poly_clear(d1q1);
            
      return;   
   } else
   {
   /* 
      We let A = a1*x^(3n-1) + a2*x^(2n-1) + a3 
      where a1 and a2 are length n and a3 is length 2n-1 
      We set p1 = a1*x^(n-1), so it has length 2n-1
   */
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length n1 and a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n2-1), so it has length n1+n2-1
   */
      temp->length = A->length - (n1+2*n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+2*n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
   
      /* 
         Set q1 to p1 div d1 
         This is an 2n-1 by n division so 
         q1 ends up being length n
         d1q1 = d1*q1 is length 2n-1
      */
      /* 
         Set q1 to p1 div d1 
         This is a 2*n1-1 by n1 division so 
         q1 ends up being length n1
         d1q1 = d1*q1 is length 2*n1-1
      */
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
   
      fmpz_poly_div_karatsuba_recursive(q1, d1q1, p1, d1); //******************************
      _fmpz_poly_stack_clear(p1);   
   }
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length 2n-1
   */
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length n1+n2-1
   */  
   
   _fmpz_poly_stack_init(d2q1, d2->length+q1->length-1, d2->limbs+q1->limbs+1); 
   _fmpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n + d2*q1
      dq1 is then of length 3n-1
   */
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length 2*n1+n2-1
   */
   
   
   _fmpz_poly_stack_init(dq1, d1q1->length + n2, B->limbs+q1->limbs+1);
   dq1->length = d1q1->length + n2;
   
   _fmpz_poly_zero_coeffs(dq1, n2);
   temp->length = d1q1->length;
   temp->limbs = dq1->limbs;
   temp->coeffs = dq1->coeffs + n2*(dq1->limbs+1);
   _fmpz_poly_set(temp, d1q1);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = p1*x^(2n-1) + p2*x^(n-1) - dq1
      which has length 3*n-1, but the first
      n coefficients will be 0, so it has 
      effective length 2n-1
   */
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1 (shifted left by 1 if n1 > n2)
      which has length 2*n1+n2-1, but we are not interested 
      in the first n1 coefficients, so it has 
      effective length n1+n2-1
   */
   
   temp->length = A->length - (n1+n2-1);
   temp->limbs = A->limbs;
   temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
   _fmpz_poly_stack_init(t, 2*n1+n2-1, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_left_shift(t, temp, n1-1);
   t->length = temp->length+n1-1;
   _fmpz_poly_sub(t, t, dq1);
   _fmpz_poly_normalise(t); 
     
   /*
      Compute q2 = t div d1
      It is a 2*n-1 by n division, so
      the length of q2 will be n
      Also compute d1q2 of length 2n-1
   */
   /*
      Compute q2 = t div d1
      It is a n1+n2-1 by n1 division, so
      the length of q2 will be n2
      Also compute d1q2 of length n1+n2-1
   */
   fmpz_poly_init(d1q2);
   fmpz_poly_init(q2);
   fmpz_poly_div_karatsuba_recursive(q2, d1q2, t, d1); //******************************
   _fmpz_poly_stack_clear(t);
      
   /*
      Compute d1*q2*x^n2
   */
   _fmpz_poly_stack_init(dq2, d1q2->length+n2, B->limbs+q2->limbs+2);
   dq2->length = d1q2->length+n2;
   _fmpz_poly_zero_coeffs(dq2, n2);
   temp->length = d1q2->length;
   temp->limbs = dq2->limbs;
   temp->coeffs = dq2->coeffs + n2*(dq2->limbs+1);
   _fmpz_poly_set(temp, d1q2);
   fmpz_poly_clear(d1q2);
   
   /* 
      Compute dq2 = d1*q2*x^n2 + d2*q2
      dq2 has length 2*n1+n2-1
   */
   _fmpz_poly_stack_init(d2q2, d2->length+q2->length-1, d2->limbs+q2->limbs+1);
   _fmpz_poly_mul(d2q2, d2, q2);
   _fmpz_poly_add(dq2, dq2, d2q2);
   _fmpz_poly_stack_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   fmpz_poly_fit_length(Q, q1->length+n2);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(q1->limbs, q2->limbs));
   _fmpz_poly_set(Q, q2);
   fmpz_poly_clear(q2);
   Q->length = q1->length + n2;
   temp->length = q1->length;
   temp->limbs = Q->limbs;
   temp->coeffs = Q->coeffs + n2*(Q->limbs+1);
   _fmpz_poly_set(temp, q1);
   fmpz_poly_clear(q1);
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length 4n-1
   */
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length 2*(n1+n2)-1
   */
   fmpz_poly_fit_limbs(BQ, FLINT_MAX(dq1->limbs, dq2->limbs)+1);
   fmpz_poly_fit_length(BQ, n2+dq1->length);
   BQ->length = n2+dq1->length;
   _fmpz_poly_zero_coeffs(BQ, n2);
   temp->length = dq1->length;
   temp->limbs = BQ->limbs;
   temp->coeffs = BQ->coeffs + n2*(BQ->limbs+1);
   _fmpz_poly_set(temp, dq1);  
   _fmpz_poly_add(BQ, BQ, dq2);
   
   _fmpz_poly_stack_clear(dq2);
   _fmpz_poly_stack_clear(dq1);
   _fmpz_poly_stack_clear(d2q1);
}

/****************************************************************************

   String conversions and I/O

****************************************************************************/


int fmpz_poly_from_string(fmpz_poly_t poly, char* s)
{
   // const char* whitespace = " \t\n\r";
   // 
   // // read poly length
   // unsigned long length;
   // if (!sscanf(s, "%ld", &length))
   //    return 0;
   //    
   // printf("length = %ld\n", length);
   // 
   // // jump to next whitespace
   // s += strcspn(s, whitespace);
   // 
   // unsigned long limbs = 1;
   // unsigned long size;
   // 
   // //mpz_t* coeffs = (mpz_t*) malloc(sizeof(mpz_t) * length);
   // mpz_t coeff;
   // mpz_init(coeff);
   // 
   // char* ptr;
   // ptr = s;
   // 
   // // count how many limbs needed to store all coefficients!
   // for (unsigned long i = 0; i < length; i++)
   // {
   //    // skip whitespace
   //    s += strspn(s, whitespace);
   //    
   //    if (!gmp_sscanf(s, "%Zd", coeff))
   //       return 0;
   //       
   //    size = mpz_sizeinbase(coeff, 2) / FLINT_BITS;
   //    
   //    printf("coeff = %s\n", mpz_get_str(NULL, 10, coeff));
   //    
   //    if(size > limbs)
   //       limbs = size;
   //       
   //    // skip to next whitespace
   //    s += strcspn(s, whitespace);
   // }
   // 
   // printf("\n\nDone!  limbs = %i\n\n", limbs);
   // 
   // s = ptr;
   // 
   // fmpz_poly_init2(poly, length, limbs);
   // fmpz_poly_set_length(poly, length);
   // 
   // for (unsigned long i = 0; i < length; i++)
   // {
   //    s += strspn(s, whitespace);
   //    if (!gmp_sscanf(s, "%Zd", coeff))
   //       return 0;
   //    mpz_to_fmpz((fmpz_t)(poly->coeffs + i*(limbs+1)), coeff);
   //    s += strcspn(s, whitespace);
   // }
   // 
   // //free(coeffs);
   // 
   // _fmpz_poly_normalise(poly);
   
   mpz_poly_t p;
   mpz_poly_init(p);
   mpz_poly_from_string(p, s);
   mpz_poly_to_fmpz_poly(poly, p);
   mpz_poly_clear(p);
   
   return 1;
}

// taken from mpz_poly_to_string, and altered using fmpz_to_mpz
char* fmpz_poly_to_string(fmpz_poly_t poly)
{
   // estimate the size of the string
   // 20 = enough room for null terminator and length info
   unsigned long size = 20;
   mpz_t temp;
   mpz_init(temp);
   for (unsigned long i = 0; i < poly->length; i++) {
      // +2 is for the sign and a space
		fmpz_to_mpz(temp, poly->coeffs + i*(poly->limbs+1));
      size += mpz_sizeinbase(temp, 10) + 2;      
   }   
   
   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "%ld  ", poly->length);
   for (unsigned long i = 0; i < poly->length; i++)
   {
		fmpz_to_mpz(temp, poly->coeffs + i*(poly->limbs+1));
		mpz_get_str(ptr, 10, temp);
      ptr += strlen(ptr);
      *ptr = ' ';
      ptr++;
   }
   
   mpz_clear(temp);
   
   ptr--;
   *ptr = 0;
   
   return buf;
}


void fmpz_poly_fprint(fmpz_poly_t poly, FILE* f)
{
   char* s = fmpz_poly_to_string(poly);
   fputs(s, f);
   free(s);
}


void fmpz_poly_print(fmpz_poly_t poly)
{
   fmpz_poly_fprint(poly, stdout);
}


int fmpz_poly_fread(fmpz_poly_t poly, FILE* f)
{
   // // read poly length
   // unsigned long length;
   // if (!fscanf(f, "%ld", &length))
   //    return 0;
   // 
   // poly->length = 0;
   // fmpz_poly_init_upto(poly, length);
   // 
   // // read coefficients
   // for (unsigned long i = 0; i < length; i++)
   // {
   //    if (!fmpz_inp_str(poly->coeffs[i], f, 10))
   //       return 0;
   //    poly->length++;
   // }
   // 
   // fmpz_poly_normalise(poly);
   
   return 1;
}


int fmpz_poly_read(fmpz_poly_t poly)
{
   return fmpz_poly_fread(poly, stdin);
}

/*
   Divide and conquer division of A by B but only computing the low half of Q*B
*/

void fmpz_poly_div_karatsuba_recursive_low(fmpz_poly_t Q, fmpz_poly_t BQ, fmpz_poly_t A, fmpz_poly_t B)
{
   _fmpz_poly_normalise(A);
   _fmpz_poly_normalise(B);
         
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      _fmpz_poly_zero(BQ);

      return;
   }
   
   if (A->length > 2*B->length - 1)
   {
      fmpz_poly_t Rb;
      fmpz_poly_init(Rb);
      fmpz_poly_divrem_naive(Q, Rb, A, B);
      fmpz_poly_fit_length(BQ, A->length);
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(A->limbs, Rb->limbs)+1);
      _fmpz_poly_sub(BQ, A, Rb);
      fmpz_poly_clear(Rb);
      
      for (unsigned long i = B->length - 1; i < BQ->length; i++)
      {
         BQ->coeffs[i*(BQ->limbs+1)] = 0;
      }
      _fmpz_poly_normalise(BQ);
      
      return;
   }
   
   fmpz_poly_t temp;
   
   if (B->length <= 16)
   {
      fmpz_poly_t Rb;
      fmpz_poly_init(Rb);
      fmpz_poly_divrem_naive_low(Q, Rb, A, B);
      fmpz_poly_fit_length(BQ, B->length - 1);
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(A->limbs, Rb->limbs)+1);
      temp->limbs = A->limbs;
      temp->length = B->length - 1;
      temp->coeffs = A->coeffs;
      _fmpz_poly_sub(BQ, temp, Rb);
      fmpz_poly_clear(Rb);
      _fmpz_poly_normalise(BQ);
      
      return;
   }
   
   fmpz_poly_t d1, d2, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   d1->length = n1;
   d2->length = n2;
   d1->limbs = B->limbs;
   d2->limbs = B->limbs;
   d1->coeffs = B->coeffs + n2*(B->limbs+1);
   d2->coeffs = B->coeffs;
   
   if (A->length <= n1+2*n2-1)
   {
      temp->length = A->length - (n1+n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
      
      fmpz_poly_init(d1q1);
      fmpz_poly_div_karatsuba_recursive_low(Q, d1q1, p1, d1); //******************************
      _fmpz_poly_stack_clear(p1);

      _fmpz_poly_stack_init(d2q1, d2->length+Q->length-1, d2->limbs+Q->limbs+1); 
      _fmpz_poly_mul(d2q1, d2, Q);
      
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(d1q1->limbs, d2q1->limbs));
      fmpz_poly_fit_length(BQ, n2+d1q1->length);
      BQ->length = n2+d1q1->length;
      _fmpz_poly_zero_coeffs(BQ, n2);
   
      temp->length = d1q1->length;
      temp->limbs = BQ->limbs;
      temp->coeffs = BQ->coeffs + n2*(BQ->limbs+1);
      _fmpz_poly_set(temp, d1q1);  
      _fmpz_poly_add(BQ, BQ, d2q1);
   
      _fmpz_poly_stack_clear(d2q1);
      fmpz_poly_clear(d1q1);
            
      return;   
   } else
   {
   /* 
      We let A = a1*x^(3n-1) + a2*x^(2n-1) + a3 
      where a1 and a2 are length n and a3 is length 2n-1 
      We set p1 = a1*x^(n-1), so it has length 2n-1
   */
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length n1 and a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n2-1), so it has length n1+n2-1
   */
      temp->length = A->length - (n1+2*n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+2*n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
   
      /* 
         Set q1 to p1 div d1 
         This is an 2n-1 by n division so 
         q1 ends up being length n
         d1q1 = d1*q1 is length 2n-1
      */
      /* 
         Set q1 to p1 div d1 
         This is a 2*n1-1 by n1 division so 
         q1 ends up being length n1
         d1q1 = d1*q1 is length 2*n1-1
      */
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
   
      fmpz_poly_div_karatsuba_recursive_low(q1, d1q1, p1, d1); //******************************
      _fmpz_poly_stack_clear(p1);   
   }
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length 2n-1
   */
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length n1+n2-1
   */  
   
   _fmpz_poly_stack_init(d2q1, d2->length+q1->length-1, d2->limbs+q1->limbs+1); 
   _fmpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n + d2*q1
      dq1 is then of length 3n-1
   */
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length 2*n1+n2-1
   */
   
   
   _fmpz_poly_stack_init(dq1, d1q1->length + n2, B->limbs+q1->limbs+1);
   dq1->length = d1q1->length + n2;
   
   _fmpz_poly_zero_coeffs(dq1, n2);
   temp->length = d1q1->length;
   temp->limbs = dq1->limbs;
   temp->coeffs = dq1->coeffs + n2*(dq1->limbs+1);
   _fmpz_poly_set(temp, d1q1);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = p1*x^(2n-1) + p2*x^(n-1) - dq1
      which has length 3*n-1, but the first
      n coefficients will be 0, so it has 
      effective length 2n-1
   */
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1 (shifted left by 1 if n1 > n2)
      which has length 2*n1+n2-1, but we are not interested 
      in the first n1 coefficients, so it has 
      effective length n1+n2-1
   */
   
   temp->length = A->length - (n1+n2-1);
   temp->limbs = A->limbs;
   temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
   _fmpz_poly_stack_init(t, 2*n1+n2-1, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_left_shift(t, temp, n1-1);
   t->length = temp->length+n1-1;
   _fmpz_poly_sub(t, t, dq1);
   t->length = n1+n2-1;
   _fmpz_poly_normalise(t); 
     
   /*
      Compute q2 = t div d1
      It is a 2*n-1 by n division, so
      the length of q2 will be n
      Also compute d1q2 of length 2n-1
   */
   /*
      Compute q2 = t div d1
      It is a n1+n2-1 by n1 division, so
      the length of q2 will be n2
      Also compute d1q2 of length n1+n2-1
   */
   fmpz_poly_init(d1q2);
   fmpz_poly_init(q2);
   fmpz_poly_div_karatsuba_recursive_low(q2, d1q2, t, d1); //******************************
   _fmpz_poly_stack_clear(t);
      
   /*
      Compute d1*q2*x^n2
   */
   _fmpz_poly_stack_init(dq2, d1q2->length+n2, B->limbs+q2->limbs+2);
   dq2->length = d1q2->length+n2;
   _fmpz_poly_zero_coeffs(dq2, n2);
   temp->length = d1q2->length;
   temp->limbs = dq2->limbs;
   temp->coeffs = dq2->coeffs + n2*(dq2->limbs+1);
   _fmpz_poly_set(temp, d1q2);
   fmpz_poly_clear(d1q2);
   
   /* 
      Compute dq2 = d1*q2*x^n2 + d2*q2
      dq2 has length 2*n1+n2-1
   */
   _fmpz_poly_stack_init(d2q2, d2->length+q2->length-1, d2->limbs+q2->limbs+1);
   _fmpz_poly_mul(d2q2, d2, q2);
   _fmpz_poly_add(dq2, dq2, d2q2);
   _fmpz_poly_stack_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   fmpz_poly_fit_length(Q, q1->length+n2);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(q1->limbs, q2->limbs));
   _fmpz_poly_set(Q, q2);
   fmpz_poly_clear(q2);
   Q->length = q1->length + n2;
   temp->length = q1->length;
   temp->limbs = Q->limbs;
   temp->coeffs = Q->coeffs + n2*(Q->limbs+1);
   _fmpz_poly_set(temp, q1);
   fmpz_poly_clear(q1);
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length 4n-1
   */
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length 2*(n1+n2)-1
   */
   fmpz_poly_fit_limbs(BQ, FLINT_MAX(dq1->limbs, dq2->limbs)+1);
   fmpz_poly_fit_length(BQ, FLINT_MAX(n2+dq1->length, dq2->length));
   BQ->length = FLINT_MAX(n2+dq1->length, dq2->length);
   _fmpz_poly_zero_coeffs(BQ, n2);
   temp->length = dq1->length;
   temp->limbs = BQ->limbs;
   temp->coeffs = BQ->coeffs + n2*(BQ->limbs+1);
   _fmpz_poly_set(temp, dq1);  
   _fmpz_poly_add(BQ, BQ, dq2);
   
   for (unsigned long i = B->length - 1; i < BQ->length; i++)
   {
      BQ->coeffs[i*(BQ->limbs+1)] = 0;
   }
   _fmpz_poly_normalise(BQ);
   
   _fmpz_poly_stack_clear(dq2);
   _fmpz_poly_stack_clear(dq1);
   _fmpz_poly_stack_clear(d2q1);
}


void fmpz_poly_div_karatsuba(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)
{
   _fmpz_poly_normalise(A);
   _fmpz_poly_normalise(B);
         
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      
      return;
   }
   
   if ((B->length <= 16) || (A->length > 2*B->length - 1))
   {
      fmpz_poly_div_naive(Q, A, B);
      
      return;
   }
   
   fmpz_poly_t d1, d2, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
      
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   d1->length = n1;
   d2->length = n2;
   d1->limbs = B->limbs;
   d2->limbs = B->limbs;
   d1->coeffs = B->coeffs + n2*(B->limbs+1);
   d2->coeffs = B->coeffs;
   
   if (A->length <= n1+2*n2-1)
   {
      temp->length = A->length - (n1+n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
      
      fmpz_poly_init(d1q1);
      fmpz_poly_div_karatsuba_recursive_low(Q, d1q1, p1, d1); //******************************
      fmpz_poly_clear(d1q1);
      _fmpz_poly_stack_clear(p1);
            
      return;   
   } else
   {
   /* 
      We let A = a1*x^(3n-1) + a2*x^(2n-1) + a3 
      where a1 and a2 are length n and a3 is length 2n-1 
      We set p1 = a1*x^(n-1), so it has length 2n-1
   */
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length n1 and a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n2-1), so it has length n1+n2-1
   */
      temp->length = A->length - (n1+2*n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+2*n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
   
      /* 
         Set q1 to p1 div d1 
         This is an 2n-1 by n division so 
         q1 ends up being length n
         d1q1 = d1*q1 is length 2n-1
      */
      /* 
         Set q1 to p1 div d1 
         This is a 2*n1-1 by n1 division so 
         q1 ends up being length n1
         d1q1 = d1*q1 is length 2*n1-1
      */
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
   
      fmpz_poly_div_karatsuba_recursive_low(q1, d1q1, p1, d1); //******************************
      _fmpz_poly_stack_clear(p1);
   }
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length 2n-1
   */
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length n1+n2-1
   */  
   
   fmpz_poly_t temp2;
   _fmpz_poly_stack_init(d2q1, d2->length+q1->length-1, d2->limbs+q1->limbs+1); 
   _fmpz_poly_mul_trunc_left_n(d2q1, d2, q1, n1 - 1);
     
   /* 
      Compute dq1 = d1*q1*x^n + d2*q1
      dq1 is then of length 3n-1
   */
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length 2*n1+n2-1
   */
   
   
   _fmpz_poly_stack_init(dq1, d1q1->length + n2, B->limbs+q1->limbs+1);
   dq1->length = d1q1->length + n2;
   
   _fmpz_poly_zero_coeffs(dq1, n2);
   temp->length = d1q1->length;
   temp->limbs = dq1->limbs;
   temp->coeffs = dq1->coeffs + n2*(dq1->limbs+1);
   _fmpz_poly_set(temp, d1q1);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);

   
   /*
      Compute t = p1*x^(2n-1) + p2*x^(n-1) - dq1
      which has length 3*n-1, but the first
      n coefficients will be 0, so it has 
      effective length 2n-1
   */
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1 (shifted left by 1 if n1 > n2)
      which has length 2*n1+n2-1, but we are not interested 
      in the first n1 coefficients, so it has 
      effective length n1+n2-1
   */
   
   temp->length = A->length - (n1+n2-1);
   temp->limbs = A->limbs;
   temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
   _fmpz_poly_stack_init(t, 2*n1+n2-1, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_left_shift(t, temp, n1-1);
   t->length = temp->length+n1-1;
   _fmpz_poly_sub(t, t, dq1);
   t->length = n1+n2-1; 
   _fmpz_poly_normalise(t); 
     
   /*
      Compute q2 = t div d1
      It is a 2*n-1 by n division, so
      the length of q2 will be n
      Also compute d1q2 of length 2n-1
   */
   /*
      Compute q2 = t div d1
      It is a n1+n2-1 by n1 division, so
      the length of q2 will be n2
      Also compute d1q2 of length n1+n2-1
   */
   fmpz_poly_init(q2);
   fmpz_poly_div_karatsuba(q2, t, d1); 
   _fmpz_poly_stack_clear(t);  
   _fmpz_poly_stack_clear(dq1);
   _fmpz_poly_stack_clear(d2q1);
      
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   fmpz_poly_fit_length(Q, q1->length+n2);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(q1->limbs, q2->limbs));
   _fmpz_poly_set(Q, q2);
   fmpz_poly_clear(q2);
   Q->length = q1->length + n2;
   temp->length = q1->length;
   temp->limbs = Q->limbs;
   temp->coeffs = Q->coeffs + n2*(Q->limbs+1);
   _fmpz_poly_set(temp, q1);
   fmpz_poly_clear(q1);
}

void fmpz_poly_divrem_karatsuba(fmpz_poly_t Q, fmpz_poly_t R, fmpz_poly_t A, fmpz_poly_t B)
{
   fmpz_poly_t QB;
   
   fmpz_poly_init(QB);
   
   fmpz_poly_div_karatsuba_recursive(Q, QB, A, B);
   
   fmpz_poly_fit_limbs(R, FLINT_MAX(QB->limbs, A->limbs)+1);
   fmpz_poly_fit_length(R, A->length);
   _fmpz_poly_sub(R, A, QB);
   _fmpz_poly_normalise(R);
   
   fmpz_poly_clear(QB);
}

/*
   Compute the polynomial X^{2n} / Q. 
   Used by Newton iteration to bootstrap power series inversion.
   Q must be monic and have length >= n.
*/

void fmpz_poly_newton_invert_basecase(fmpz_poly_t Q_inv, fmpz_poly_t Q, unsigned long n)
{
   fmpz_poly_t X2n, Qn;
   
   fmpz_poly_init2(X2n, 2*n-1, 1);
   _fmpz_poly_zero_coeffs(X2n, 2*n - 2);
   _fmpz_poly_set_coeff_ui(X2n, 2*n - 2, 1);
   X2n->length = 2*n-1;
   
   Qn->coeffs = Q->coeffs + (Q->length - n)*(Q->limbs + 1);
   Qn->limbs = Q->limbs;
   Qn->length = n;
   
   fmpz_poly_div_karatsuba(Q_inv, X2n, Qn); 
   
   fmpz_poly_clear(X2n);
}

#define FLINT_NEWTON_INVERSE_BASECASE_CUTOFF 32

/*
   Recursively compute 1 / Q mod x^n using Newton iteration
   Assumes Q is given to the full precision n required and has constant term 1
*/

void fmpz_poly_newton_invert(fmpz_poly_t Q_inv, fmpz_poly_t Q, unsigned long n)
{
   if (n < FLINT_NEWTON_INVERSE_BASECASE_CUTOFF)
   {
      fmpz_poly_t Q_rev;
      fmpz_poly_init(Q_rev);
      fmpz_poly_fit_length(Q_rev, n);
      fmpz_poly_fit_limbs(Q_rev, Q->limbs);
      _fmpz_poly_reverse(Q_rev, Q, n);
      fmpz_poly_newton_invert_basecase(Q_inv, Q_rev, n);
      fmpz_poly_fit_length(Q_inv, n);
      _fmpz_poly_reverse(Q_inv, Q_inv, n);
      fmpz_poly_clear(Q_rev);
      
      return;
   }
   
   unsigned long m = (n+1)/2;
   
   fmpz_poly_t g0, prod, prod2;
   fmpz_poly_init(g0);
   fmpz_poly_init(prod);
   fmpz_poly_init(prod2);
   fmpz_poly_newton_invert(g0, Q, m);
   fmpz_poly_mul_trunc_n(prod, Q, g0, n);
   __fmpz_poly_sub_coeff_ui(prod->coeffs, 1);
   fmpz_poly_mul_trunc_n(prod2, prod, g0, n);
   fmpz_poly_fit_length(Q_inv, n);
   fmpz_poly_fit_limbs(Q_inv, FLINT_MAX(prod2->limbs, g0->limbs)+1);
   _fmpz_poly_sub(Q_inv, g0, prod2);
   
   fmpz_poly_clear(prod2);
   fmpz_poly_clear(prod);
   fmpz_poly_clear(g0);
}

/* 
   Yields a precision n power series quotient of A by B assuming A and B are both 
   given to precision n and B is normalised (i.e. constant coefficient is 1).
*/

void fmpz_poly_div_series(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B, unsigned long n)
{
   fmpz_poly_t B_inv;
   fmpz_poly_init(B_inv);
   fmpz_poly_newton_invert(B_inv, B, n);
   fmpz_poly_mul_trunc_n(Q, B_inv, A, n);
   
   fmpz_poly_clear(B_inv);
}

/*
   Polynomial division of A by B
   The remainder is not computed, to save time
   B is assumed to be monic
*/

void fmpz_poly_div_newton(fmpz_poly_t Q, fmpz_poly_t A, fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      fmpz_poly_set_coeff_si(Q, 0, 0);
      _fmpz_poly_normalise(Q);
      return;
   }
   
   fmpz_poly_t A_rev, B_rev;
   fmpz_poly_init2(A_rev, A->length, A->limbs);
   fmpz_poly_init2(B_rev, B->length, B->limbs);
   
   _fmpz_poly_reverse(A_rev, A, A->length);
   _fmpz_poly_reverse(B_rev, B, B->length);
   
   fmpz_poly_div_series(Q, A_rev, B_rev, A->length - B->length + 1);
   
   fmpz_poly_fit_length(Q, A->length - B->length + 1);
   _fmpz_poly_reverse(Q, Q, A->length - B->length + 1);
   _fmpz_poly_normalise(Q);
   
   fmpz_poly_clear(B_rev);
   fmpz_poly_clear(A_rev);
}

void fmpz_poly_power(fmpz_poly_t output, fmpz_poly_t poly, unsigned long exp)
{
   fmpz_poly_t power, temp;
   
   fmpz_poly_init(power);
   fmpz_poly_init(temp);
   
   if (exp == 0) 
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      fmpz_poly_set_coeff_ui(output, 0, 1);
      output->length = 1;
      return;
   }
   if ((poly->length == 1) && (poly->coeffs[0] == 1) && (poly->coeffs[1] == 1))
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      fmpz_poly_set_coeff_ui(output, 0, 1);
      output->length = 1;
      return;
   } 
   if (poly->length == 0)
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      output->length = 0;
      return;      
   }
   
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, poly->limbs);
   _fmpz_poly_set(output, poly);
    
   while (!(exp & 1L))
   {
      fmpz_poly_mul(temp, output, output);
      fmpz_poly_fit_length(output, temp->length);
      fmpz_poly_fit_limbs(output, temp->limbs);
      _fmpz_poly_set(output, temp);
      exp >>= 1;
   }
   
   exp >>= 1;
   if (exp)
   {
      fmpz_poly_fit_length(power, output->length);
      fmpz_poly_fit_limbs(power, output->limbs);
      _fmpz_poly_set(power, output);
      
      while (exp)
      {
         fmpz_poly_mul(temp, power, power);
         fmpz_poly_fit_length(power, temp->length);
         fmpz_poly_fit_limbs(power, temp->limbs);
         _fmpz_poly_set(power, temp);
         if (exp & 1) 
         {
            fmpz_poly_mul(temp, output, power);
            fmpz_poly_fit_length(output, temp->length);
            fmpz_poly_fit_limbs(output, temp->limbs);
            _fmpz_poly_set(output, temp);
         }
         exp >>= 1;
      }
   }
}
