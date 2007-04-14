/****************************************************************************

Zpoly_mpn.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"
#include "flint-manager.h"

/****************************************************************************

   _Zpoly_mpn_* layer

****************************************************************************/

void _Zpoly_mpn_convert_out(Zpoly_t poly_mpz, Zpoly_mpn_t poly_mpn)
{
   poly_mpz->length = poly_mpn->length;
   for (unsigned long i = 0; i < poly_mpn->length; i++)
   {
       if (poly_mpn->coeffs[i*(poly_mpn->limbs+1)] == 0) 
          mpz_set_ui(poly_mpz->coeffs[i],0);
       else 
       {
          mpz_import(poly_mpz->coeffs[i], poly_mpn->limbs, -1, 
          sizeof(mp_limb_t), 0, 0, poly_mpn->coeffs+i*(poly_mpn->limbs+1)+1);
          if (poly_mpn->coeffs[i*(poly_mpn->limbs+1)] < 0)
             mpz_neg(poly_mpz->coeffs[i],poly_mpz->coeffs[i]);
       }
   }
}

void _Zpoly_mpn_convert_in(Zpoly_mpn_t poly_mpn, Zpoly_t poly_mpz)
{
   poly_mpn->length = poly_mpz->length;
   for (unsigned long i = 0; i < poly_mpz->length; i++)
   {
      if (mpz_sgn(poly_mpz->coeffs[i]) == 0) 
          _Zpoly_mpn_set_coeff_ui(poly_mpn,i,0);
      else
      {
          size_t countp;
          mpz_export(poly_mpn->coeffs+i*(poly_mpn->limbs+1)+1, &countp, 
             -1, sizeof(mp_limb_t), 0, 0, poly_mpz->coeffs[i]);
          if (mpz_sgn(poly_mpz->coeffs[i]) < 0) 
             poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = -1L;
          else poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = 1L;
      }
   }
}

/* 
   Set a coefficient to the given unsigned value.
   Clears dirty limbs unless the coefficient is set to zero. 
   Sets the sign to 1 if x is positive, else to zero.
*/

void _Zpoly_mpn_set_coeff_ui(Zpoly_mpn_t poly, unsigned long n, unsigned long x)
{
   if (x == 0) 
   {
      poly->coeffs[n*(poly->limbs+1)] = 0UL;
      return;
   }
   poly->coeffs[n*(poly->limbs+1)] = 1UL;
   poly->coeffs[n*(poly->limbs+1)+1] = x;
   if (poly->limbs > 1) 
          clear_limbs(poly->coeffs+n*(poly->limbs+1)+2, poly->limbs - 1);
}

void _Zpoly_mpn_set_coeff_si(Zpoly_mpn_t poly, unsigned long n, long x)
{
   if (x == 0)
   {
      poly->coeffs[n*(poly->limbs+1)] = 0UL;
      return;
   }
   if (x > 0)
   {
      poly->coeffs[n*(poly->limbs+1)] = 1L;
      poly->coeffs[n*(poly->limbs+1)+1] = x;
   } else if (x < 0)
   {
      poly->coeffs[n*(poly->limbs+1)] = -1L;
      poly->coeffs[n*(poly->limbs+1)+1] = -x;
   } 
   if (poly->limbs > 1) 
          clear_limbs(poly->coeffs+n*(poly->limbs+1)+2, poly->limbs - 1);
}

void _Zpoly_mpn_normalise(Zpoly_mpn_t poly)
{
   while (poly->coeffs[(poly->length-1)*(poly->limbs+1)] == 0)
   {
      poly->length--;
   }
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
      if (input->limbs == output->limbs)
      {
         copy_limbs(output->coeffs, input->coeffs, input->length*(input->limbs+1));
      } else
      {
         unsigned long diff_size = output->limbs - input->limbs;
         unsigned long input_size = input->limbs + 1;
         unsigned long output_size = output->limbs + 1;
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

int _Zpoly_mpn_equal(Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   if (input1->length != input2-> length) return 0;
   if (input1->limbs == input2->limbs)
   {
      for (unsigned long i = 0; i < input1->length*(input1->limbs+1); i++)
         if (input1->coeffs[i] != input2->coeffs[i]) return 0;
      return 1;
   }
    
   unsigned long i,j;
   if (input1->limbs > input2->limbs)
   {
      for (i = 0; i < input1->length; i++)
      {
         for (j = 0; j < input2->limbs+1; j++)
         {
            if (input1->coeffs[i*(input1->limbs+1)+j] != input2->coeffs[i*(input2->limbs+1)+j])
               return 0;
         }
         for (j = input2->limbs + 1; j < input1->limbs + 1; j++)
         {
            if (input1->coeffs[i*(input1->limbs+1)+j] != 0) return 0;
         }
      }
      return 1;
   } else
   {
      for (i = 0; i < input1->length; i++)
      {
         for (j = 0; j < input1->limbs+1; j++)
         {
            if (input1->coeffs[i*(input1->limbs+1)+j] != input2->coeffs[i*(input2->limbs+1)+j])
               return 0;
         }
         for (j = input1->limbs + 1; j < input2->limbs + 1; j++)
         {
            if (input2->coeffs[i*(input1->limbs+1)+j] != 0) return 0;
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
         output->coeffs[i*(output->limbs+1)] = -output->coeffs[i*(output->limbs+1)];
   } else if (input->limbs == output->limbs)
   {
      copy_limbs(output->coeffs, input->coeffs, input->length*(input->limbs+1));
      for (unsigned long i = 0; i < input->length; i++)
         output->coeffs[i*(input->limbs+1)] = -input->coeffs[i*(input->limbs+1)];
   } else
   {
      unsigned long diff_size = output->limbs - input->limbs;
      unsigned long input_size = input->limbs + 1;
      unsigned long output_size = output->limbs + 1;
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
   if ((input->coeffs == output->coeffs) && (n < output->length))
   {
      Zpoly_mpn_t part2;
      part->length = input->length - n;
      part->limbs = input->limbs;
      part->coeffs = input->coeffs + n*(input->limbs+1);
      part2->length = input->length - n;
      part2->limbs = input->limbs;
      part2->coeffs = input->coeffs + 2*n*(input->limbs+1);
      _Zpoly_mpn_set(part2, part);
      part->length = n;
      part2->length = n;
      part2->coeffs = input->coeffs;
      _Zpoly_mpn_set(part2, part);
      clear_limbs(output->coeffs, n*(output->limbs+1));
   } else
   {
      part->length = input->length;
      part->limbs = output->limbs;
      part->coeffs = output->coeffs + n*(output->limbs+1);
      _Zpoly_mpn_set(part, input);
      clear_limbs(output->coeffs, n*(output->limbs+1));
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
   part->limbs = output->limbs;
   part->coeffs = input->coeffs + n*(output->limbs + 1);
   _Zpoly_mpn_set(part, input);
   clear_limbs(output->coeffs, n*(output->limbs+1));  
   output->length = part->length; 
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

   poly->coeffs = (mp_limb_t *) flint_malloc(sizeof(mp_limb_t)*alloc*(limbs+1));
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

/* Shrink or expand a polynomial to "alloc" coefficients */

void Zpoly_mpn_realloc(Zpoly_mpn_t poly, unsigned long alloc)
{
   if (alloc <= 0) alloc = 1;
   poly->coeffs = (mp_limb_t*) flint_realloc(poly->coeffs, sizeof(mp_limb_t)*alloc*(poly->limbs+1));
   
   poly->alloc = alloc;
   
   // truncate actual data if necessary
   if (poly->length > alloc)
      poly->length = alloc;
}

void Zpoly_mpn_clear(Zpoly_mpn_t poly)
{
   flint_free(poly->coeffs);
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
