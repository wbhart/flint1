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
