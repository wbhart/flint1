/****************************************************************************

Zpoly_mpn.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "flint-manager.h"

/****************************************************************************

   _Zpoly_mpn_* layer

****************************************************************************/

/* 
   Set a coefficient to the given unsigned value.
   Clears dirty limbs. Sets the sign to 1 if x is positive, else to zero.
*/

void _Zpoly_mpn_set_coeff_ui(Zpoly_mpn_t poly, unsigned long n, unsigned long x)
{
   if (x == 0) poly->coeffs[n*(poly->coeff_size+1)] = 0UL;
   else poly->coeffs[n*(poly->coeff_size+1)] = 1UL;
   poly->coeffs[n*(poly->coeff_size+1)+1] = x;
   if (poly->coeff_size > 1) 
          clear_limbs(poly->coeffs+n*(poly->coeff_size+1)+2, poly->coeff_size - 1);
}

void _Zpoly_mpn_set_coeff_si(Zpoly_mpn_t poly, unsigned long n, long x)
{
   if (x == 0)
   {
      clear_limbs(poly->coeffs+n*(poly->coeff_size+1), poly->coeff_size+1);
      return;
   }
   if (x > 0)
   {
      poly->coeffs[n*(poly->coeff_size+1)] = 1L;
      poly->coeffs[n*(poly->coeff_size+1)+1] = x;
   } else if (x < 0)
   {
      poly->coeffs[n*(poly->coeff_size+1)] = -1L;
      poly->coeffs[n*(poly->coeff_size+1)+1] = -x;
   } 
   if (poly->coeff_size > 1) 
          clear_limbs(poly->coeffs+n*(poly->coeff_size+1)+2, poly->coeff_size - 1);
}

/* 
   Sets the output poly to equal the input poly 
   Assumes the output poly is at least as big as the input poly
   Assumes polynomials don't overlap in a bad way
*/

void _Zpoly_mpn_set(Zpoly_mpn_t output, Zpoly_mpn_t input)
{
   if (output->coeffs != input->coeffs) 
   {
      if (input->coeff_size == output->coeff_size)
      {
         copy_limbs(output->coeffs, input_coeffs, input->length*(input->coeff_size+1));
      } else
      {
         unsigned long diff_size = output->coeff_size - input->coeff->size;
         unsigned long input_size = input->coeff_size + 1;
         unsigned long output_size = output->coeff_size + 1;
         for (unsigned long i = 0; i < input->length; i++)
         {
             copy_limbs(output->coeffs+i*output_size, input->coeffs+i*input_size, input_size);
             clear_limbs(output->coeffs+i*output_size+input_size, diff_size);
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
   y->coeffs = temp->coeffs;
   
   temp_l = x->alloc;
   x->alloc = y->alloc;
   y->alloc = temp_l;
   
   temp_l = x->length;
   x->length = y->length;
   y->length = temp_l;
   
   temp_l = x->coeff_size;
   x->coeff_size = y->coeff_size;
   y->coeff_size = temp_l;
}

int _Zpoly_mpn_equal(Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   if (input1->length != input2-> length) return 0;
   if (input1->coeff_size == input2->coeff_size)
   {
      for (unsigned long i = 0; i < input1->length*(input1->coeff_size+1); i++)
         if (input1->coeffs[i] != input2->coeffs[i]) return 0;
      return 1;
   }
    
   unsigned long i,j;
   if (input1->coeff_size > input2->coeff_size)
   {
      for (i = 0; i < input1->length; i++)
      {
         for (j = 0; j < input2->coeff_size+1; j++)
         {
            if (input1[i*(input1->coeff_size+1)+j] != input2[i*(input2->coeff_size+1)+j])
               return 0;
         }
         for (j = input2->coeff_size + 1; j < input1->coeff_size + 1; j++)
         {
            if (input1[i*(input1->coeff_size+1)+j] != 0) return 0;
         }
      }
      return 1;
   } else
   {
      for (i = 0; i < input1->length; i++)
      {
         for (j = 0; j < input1->coeff_size+1; j++)
         {
            if (input1[i*(input1->coeff_size+1)+j] != input2[i*(input2->coeff_size+1)+j])
               return 0;
         }
         for (j = input1->coeff_size + 1; j < input2->coeff_size + 1; j++)
         {
            if (input2[i*(input1->coeff_size+1)+j] != 0) return 0;
         }
      }
      return 1;
   }
}

void _Zpoly_mpn_negate(Zpoly_mpn_t output, Zpoly_mpn_t input)
{
   if (input->coeffs == output->coeffs)
   {
      for (unsigned long i = 0; i < input->length; i++)
         output->coeffs[i*(output->coeff_size+1)] = -output->coeffs[i*(output->coeff_size+1)];
   } else if (input->coeff_size == output->coeff_size)
   {
      copy_limbs(output->coeffs, input->coeffs, input->length*(input->coeff_size+1));
      for (unsigned long i = 0; i < input->length; i++)
         output->coeffs[i*(input->coeff_size+1)] = -input->coeffs[i*(input->coeff_size+1)];
   } else
   {
      unsigned long diff_size = output->coeff_size - input->coeff->size;
      unsigned long input_size = input->coeff_size + 1;
      unsigned long output_size = output->coeff_size + 1;
      for (unsigned long i = 0; i < input->length; i++)
      {
         output->coeffs[i*output_size] = -input->coeffs[i*input_size];
         copy_limbs(output->coeffs+i*output_size+1, input->coeffs+i*input_size+1, input_size-1);
         clear_limbs(output->coeffs+i*output_size+input_size, diff_size);
      }
   }
   output->length = input->length;
}

/* 
   Multiplies input by x^n and sets output to the result
   Assumes output is large enough to contain the result
   If input and output are part of the same polynomial, but 
   not starting at the same point, we assume they do not 
   overlap in a bad way
*/
void _Zpoly_mpn_left_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, 
                                                 unsigned long n)
{
   if (n == 0) return;
   
   Zpoly_mpn_t part;   
   if ((input->coeffs == output->coeffs) && (n < output_length))
   {
      Zpoly part2;
      part->length = input->length - n;
      part->coeff_size = input->coeff_size;
      part->coeffs = input->coeffs + n*(input->coeff_size+1);
      part2->length = input->length - n;
      part2->coeff_size = input->coeff_size;
      part2->coeffs = input->coeffs + 2*n*(input->coeff_size+1);
      _Zpoly_mpn_set(part2, part);
      part->length = n;
      part2->length = n;
      part2->coeffs = input->coeffs;
      _Zpoly_mpn_set(part2, part);
      clear_limbs(output->coeffs, n*(output->coeff_size+1);
   } else
   {
      part->length = input->length;
      part->coeff_size = output->coeff_size;
      part->coeffs = output->coeffs + n*(output->coeff_size+1);
      _Zpoly_mpn_set(part, input);
      clear_limbs(output->coeffs, n*(output->coeff_size+1);
   }
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
   part->coeff_size = output->coeff_size;
   part->coeffs = input->coeffs + n*(output->coeff_size + 1);
   _Zpoly_mpn_set(part, input);
   clear_limbs(output->coeffs, n*(output->coeff_size+1);  
   output->length = part->length; 
}

