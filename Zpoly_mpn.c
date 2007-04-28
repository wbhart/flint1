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
#include "flint-manager.h"

/****************************************************************************

   _Zpoly_mpn_* layer

****************************************************************************/

void _Zpoly_mpn_convert_out(Zpoly_t poly_mpz, Zpoly_mpn_t poly_mpn)
{
   FLINT_ASSERT(poly_mpz->alloc >= poly_mpn->length);

   poly_mpz->length = poly_mpn->length;
   for (unsigned long i = 0; i < poly_mpn->length; i++)
   {
       if (poly_mpn->coeffs[i*(poly_mpn->limbs+1)] == 0) 
          mpz_set_ui(poly_mpz->coeffs[i], 0);
       else 
       {
          mpz_import(poly_mpz->coeffs[i], poly_mpn->limbs, -1, 
                     sizeof(mp_limb_t), 0, 0,
                     poly_mpn->coeffs + i*(poly_mpn->limbs+1) + 1);

          if (poly_mpn->coeffs[i*(poly_mpn->limbs+1)] == -1L)
          {
             mpz_neg(poly_mpz->coeffs[i], poly_mpz->coeffs[i]);
          }
       }
   }
}

void _Zpoly_mpn_convert_in(Zpoly_mpn_t poly_mpn, Zpoly_t poly_mpz)
{
   FLINT_ASSERT(poly_mpn->alloc >= poly_mpz->length);

   poly_mpn->length = poly_mpz->length;
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
                     
          for (unsigned long pad = countp; pad < poly_mpn->limbs; pad++)
          {
              poly_mpn->coeffs[i*(poly_mpn->limbs+1)+pad+1] = 0L;
          }

          if (mpz_sgn(poly_mpz->coeffs[i]) < 0) 
          {
             poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = -1L;
          } else
             poly_mpn->coeffs[i*(poly_mpn->limbs+1)] = 1L;
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
   if (poly->limbs > 1) 
      clear_limbs(poly->coeffs + n*(poly->limbs+1) + 2, poly->limbs - 1);
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
   if (poly->limbs > 1) 
      clear_limbs(poly->coeffs + n*(poly->limbs+1) + 2, poly->limbs - 1);
}

void _Zpoly_mpn_normalise(Zpoly_mpn_t poly)
{
   while (poly->length && poly->coeffs[(poly->length-1)*(poly->limbs+1)] == 0)
      poly->length--;
}

/* 
   Sets the output poly to equal the input poly 
   Assumes the output poly is at least as big as the input poly
*/

void _Zpoly_mpn_set(Zpoly_mpn_t output, Zpoly_mpn_t input)
{
   if (output->coeffs != input->coeffs) 
   {
      if ((output->coeffs < input->coeffs) || (output->coeffs >= input->coeffs + input->length*(input->limbs+1)))
      {
         if (input->limbs == output->limbs)
         {
            forward_copy_limbs(output->coeffs, input->coeffs, input->length*(input->limbs+1));
         } else
         {
            unsigned long diff_size = output->limbs - input->limbs;
            unsigned long input_size = input->limbs + 1;
            unsigned long output_size = output->limbs + 1;
            for (unsigned long i = 0; i < input->length; i++)
            {
               forward_copy_limbs(output->coeffs+i*output_size, input->coeffs+i*input_size, input_size);
               clear_limbs(output->coeffs+i*output_size+input_size, diff_size);
            }
         }
      } else
      {
         unsigned long shift_limbs = output->coeffs - input->coeffs;
         copy_limbs(output->coeffs+shift_limbs, output->coeffs, input->length*(input->limbs+1)-shift_limbs);
         copy_limbs(output->coeffs, input->coeffs, shift_limbs);
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
   int shorter_poly;
   unsigned long shorter_length;
   unsigned long i,j;
   
   if (input1->length > input2->length)
   {
      shorter_length = input2->length;
      shorter_poly = 2;
   } else
   {
      shorter_length = input1->length;
      shorter_poly = 1;
   }
   
   if (input1->limbs == input2->limbs)
   {
      for (i = 0; i < shorter_length; i++)
      {  
         if (input1->coeffs[i*(input1->limbs+1)] == 0)
         {
             if (input2->coeffs[i*(input2->limbs+1)] != 0) 
                return 0;
         } else                                    
         {  
            for (j = 0; j < input1->limbs+1; j++)
               if (input1->coeffs[i*(input1->limbs+1)+j] != input2->coeffs[i*(input1->limbs+1)+j]) 
                  return 0;
         }
      }
      if (shorter_poly == 1) 
      {   
         for (i = input1->length; i < input2->length; i++)
            if (input2->coeffs[i*(input2->limbs+1)] != 0L) return 0;
      } else
      {   
         for (i = input2->length; i < input1->length; i++)
            if (input1->coeffs[i*(input1->limbs+1)] != 0L) return 0;
      }
      
      return 1;
   }
    
   if (input1->limbs > input2->limbs)
   {
      for (i = 0; i < shorter_length; i++)
      {
         if (input1->coeffs[i*(input1->limbs+1)] == 0)
         {
             if (input2->coeffs[i*(input2->limbs+1)] != 0) 
                return 0; 
         } else
         {
            for (j = 0; j < input2->limbs+1; j++)
            {
               if (input1->coeffs[i*(input1->limbs+1)+j] != input2->coeffs[i*(input2->limbs+1)+j]) 
                  return 0; 
            }
            for (j = input2->limbs + 1; j < input1->limbs + 1; j++)
            {
               if (input1->coeffs[i*(input1->limbs+1)+j] != 0)
                  return 0; 
            }
         }
      }
      if (shorter_poly == 1) 
      {   
         for (unsigned long i = input1->length; i < input2->length; i++)
            if (input2->coeffs[i*(input2->limbs+1)] != 0L) return 0;
      } else
      {   
         for (unsigned long i = input2->length; i < input1->length; i++)
            if (input1->coeffs[i*(input1->limbs+1)] != 0L) return 0;
      }
      
      return 1;
   } else
   {
      for (i = 0; i < shorter_length; i++)
      {
         if (input1->coeffs[i*(input1->limbs+1)] == 0)
         {
             if (input2->coeffs[i*(input2->limbs+1)] != 0) return 0;
         } else
         {
            for (j = 0; j < input1->limbs+1; j++)
            {
               if (input1->coeffs[i*(input1->limbs+1)+j] != input2->coeffs[i*(input2->limbs+1)+j])
                  return 0;
            }
            for (j = input1->limbs + 1; j < input2->limbs + 1; j++)
            {
               if (input2->coeffs[i*(input2->limbs+1)+j] != 0) return 0;
            }
         }
      }
      if (shorter_poly == 1) 
      {   
         for (unsigned long i = input1->length; i < input2->length; i++)
            if (input2->coeffs[i*(input2->limbs+1)] != 0L) return 0;
      } else
      {   
         for (unsigned long i = input2->length; i < input1->length; i++)
            if (input1->coeffs[i*(input1->limbs+1)] != 0L) return 0;
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
*/

void _Zpoly_mpn_left_shift(Zpoly_mpn_t output, Zpoly_mpn_t input, 
                                                 unsigned long n)
{
   Zpoly_mpn_t part;   
   
   part->length = input->length;
   part->limbs = output->limbs;
   part->coeffs = output->coeffs + n*(output->limbs+1);
      
   _Zpoly_mpn_set(part, input);
   clear_limbs(output->coeffs, n*(output->limbs+1));
   
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
    Add two polynomials together 
*/

void _Zpoly_mpn_add(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   unsigned long size1, size2, shorter, size_out;
   mp_limb_t * coeffs1, * coeffs2, * coeffs_out;
   
   if (input1->limbs > input2->limbs)
   {
      size1 = input1->limbs+1;
      size2 = input2->limbs+1;
      coeffs1 = input1->coeffs;
      coeffs2 = input2->coeffs;
   } else
   {
      size1 = input2->limbs+1;
      size2 = input1->limbs+1;
      coeffs1 = input2->coeffs;
      coeffs2 = input1->coeffs;
   }
   
   size_out = output->limbs+1;
   coeffs_out = output->coeffs;
   
   long carry;
   
   shorter = (input1->length > input2->length) ? input2->length : input1->length;
   
   for (unsigned long i = 0; i < shorter; i++)
   {
       if (!coeffs1[i*size1])
       {
          if (!coeffs2[i*size2]) coeffs_out[i*size_out] = 0L;
          else
          {
              copy_limbs(coeffs_out+i*size_out, coeffs2+i*size2, size2);
              if (size_out > size2) clear_limbs(coeffs_out+i*size_out+size2, size_out - size2);
          }
       } else if (!coeffs2[i*size2])
       {
          copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, size1);
          if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
       } else if (coeffs1[i*size1] == coeffs2[i*size2])
       {
          coeffs_out[i*size_out] = coeffs1[i*size1];
          carry = mpn_add(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, coeffs2+i*size2+1, size2-1);
          if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
          if (carry) coeffs_out[i*size_out+size1] = carry;
       } else
       {
          carry = 0;
          for (unsigned long j = size2; (!carry) && (j < size1); j++) carry |= coeffs1[i*size1+j];
          if (carry) carry = 1;
          else carry = mpn_cmp(coeffs1+i*size1+1, coeffs2+i*size2+1, size2-1); 
          
          if (carry == 0) coeffs_out[i*size_out] = 0L;
          else if (carry > 0) 
          {
             mpn_sub(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, coeffs2+i*size2+1, size2-1);
             coeffs_out[i*size_out] = coeffs1[i*size1];
             if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
          }
          else
          {
             mpn_sub(coeffs_out+i*size_out+1, coeffs2+i*size2+1, size2-1, coeffs1+i*size1+1, size2-1);
             coeffs_out[i*size_out] = -coeffs1[i*size1];
             if (size_out > size2) clear_limbs(coeffs_out+i*size_out+size2, size_out - size2);
          }
       }
   }
   
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   
   for (unsigned long i = shorter; i < input1->length; i++)
   {
       if (!coeffs1[i*size1]) coeffs_out[i*size_out] = 0L;
       else
       {
           copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, size1);
           if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);   
       }
   }
   for (unsigned long i = shorter; i < input2->length; i++)
   {
       if (!coeffs2[i*size2]) coeffs_out[i*size_out] = 0L;
       else
       {
           copy_limbs(coeffs_out+i*size_out, coeffs2+i*size2, size2);
           if (size_out > size2) clear_limbs(coeffs_out+i*size_out+size2, size_out - size2);   
       }
   }
   
   output->length = (input1->length > input2->length) ? input1->length : input2->length;
}

/* 
    Set output poly to input1 - input2 
*/

void _Zpoly_mpn_sub(Zpoly_mpn_t output, Zpoly_mpn_t input1, Zpoly_mpn_t input2)
{
   unsigned long size1, size2, shorter, size_out;
   int in_order = 1;
   mp_limb_t * coeffs1, * coeffs2, * coeffs_out;
   
   if (input1->limbs > input2->limbs)
   {
      size1 = input1->limbs+1;
      size2 = input2->limbs+1;
      coeffs1 = input1->coeffs;
      coeffs2 = input2->coeffs;
   } else
   {
      size1 = input2->limbs+1;
      size2 = input1->limbs+1;
      coeffs1 = input2->coeffs;
      coeffs2 = input1->coeffs;
      in_order = 0;
   }
   
   size_out = output->limbs+1;
   coeffs_out = output->coeffs;
   
   long carry;
   
   shorter = (input1->length > input2->length) ? input2->length : input1->length;
   
   for (unsigned long i = 0; i < shorter; i++)
   {
       if (!coeffs1[i*size1])
       {
          if (!coeffs2[i*size2]) coeffs_out[i*size_out] = 0L;
          else
          {
              copy_limbs(coeffs_out+i*size_out, coeffs2+i*size2, size2);
              if (size_out > size2) clear_limbs(coeffs_out+i*size_out+size2, size_out - size2);
              if (in_order) coeffs_out[i*size_out] = -coeffs_out[i*size_out];
          }
       } else if (!coeffs2[i*size2])
       {
          copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, size1);
          if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
          if (!in_order) coeffs_out[i*size_out] = -coeffs_out[i*size_out];
       } else if (coeffs1[i*size1] != coeffs2[i*size2])
       {
          if (in_order) coeffs_out[i*size_out] = coeffs1[i*size1];
          else coeffs_out[i*size_out] = -coeffs1[i*size1];
          carry = mpn_add(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, coeffs2+i*size2+1, size2-1);
          if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
          if (carry) coeffs_out[i*size_out+size1] = carry;
       } else
       {
          carry = 0;
          for (unsigned long j = size2; (!carry) && (j < size1); j++) carry |= coeffs1[i*size1+j];
          if (carry) carry = 1;
          else carry = mpn_cmp(coeffs1+i*size1+1, coeffs2+i*size2+1, size2-1); 
          
          if (carry == 0) coeffs_out[i*size_out] = 0L;
          else if (carry > 0) 
          {
             mpn_sub(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, coeffs2+i*size2+1, size2-1);
             if (in_order) coeffs_out[i*size_out] = coeffs1[i*size1];
             else coeffs_out[i*size_out] = -coeffs1[i*size1];
             if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
          }
          else
          {
             mpn_sub(coeffs_out+i*size_out+1, coeffs2+i*size2+1, size2-1, coeffs1+i*size1+1, size2-1);
             if (in_order) coeffs_out[i*size_out] = -coeffs1[i*size1];
             else coeffs_out[i*size_out] = coeffs1[i*size1];
             if (size_out > size2) clear_limbs(coeffs_out+i*size_out+size2, size_out - size2);
          }
       }
   }
   
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   
   for (unsigned long i = shorter; i < input1->length; i++)
   {
       if (!coeffs1[i*size1]) coeffs_out[i*size_out] = 0L;
       else
       {
           copy_limbs(coeffs_out+i*size_out, coeffs1+i*size1, size1);
           if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);   
       }
   }
   for (unsigned long i = shorter; i < input2->length; i++)
   {
       if (!coeffs2[i*size2]) coeffs_out[i*size_out] = 0L;
       else
       {
           copy_limbs(coeffs_out+i*size_out, coeffs2+i*size2, size2);
           if (size_out > size2) clear_limbs(coeffs_out+i*size_out+size2, size_out - size2);   
           coeffs_out[i*size_out] = -coeffs_out[i*size_out];
       }
   }
   
   output->length = (input1->length > input2->length) ? input1->length : input2->length;
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
           mslimb = mpn_mul_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, x);
           if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
           if (mslimb) coeffs_out[i*size_out+size1] = mslimb;     
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
              mslimb = mpn_mul_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, -x);
              if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
              if (mslimb) coeffs_out[i*size_out+size1] = mslimb;     
           }
        } else 
        {
           if ((coeffs_out[i*size_out] = coeffs1[i*size1]))
           {
              mslimb = mpn_mul_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, x);
              if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out - size1);
              if (mslimb) coeffs_out[i*size_out+size1] = mslimb;     
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
         coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, x);
         if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out-size1);
      }
   } else
   {
      if (coeffs_out != coeffs1)
      {
         copy_limbs(coeffs_out, coeffs1, size1*poly->length); 
         for (unsigned long i = 0; i < poly->length; i++)
         {
            if (!coeffs_out[i*size_out]) clear_limbs(coeffs_out+i*size_out+1, size_out-1);
            else coeffs_out[i*size_out] = 0L;
         }
         mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            coeffs_out[i*size_out] = coeffs1[i*size1];
         }
      } else
      {
         mp_limb_t * signs = (mp_limb_t *) flint_malloc_limbs(poly->length);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            signs[i] = coeffs_out[i*size_out];
            if (!coeffs_out[i*size_out]) clear_limbs(coeffs_out+i*size_out+1, size_out-1);
            else coeffs_out[i*size_out] = 0L;
         }
         mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            coeffs_out[i*size_out] = signs[i];
         }
         flint_free(signs);
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
            coeffs_out[i*size_out] = -coeffs1[i*size1];
            mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, -x);
         }
         else 
         {
            coeffs_out[i*size_out] = coeffs1[i*size1];
            mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, x);
         }
         if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out-size1);
      }
   } else
   {
      if (coeffs_out != coeffs1)
      {
         copy_limbs(coeffs_out, coeffs1, size1*poly->length); 
         for (unsigned long i = 0; i < poly->length; i++)
         {
            if (!coeffs_out[i*size_out]) clear_limbs(coeffs_out+i*size_out+1, size_out-1);
            else coeffs_out[i*size_out] = 0L;
         }
         if (x < 0) 
         {
            mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, -x);
            for (unsigned long i = 0; i < poly->length; i++)
            {
               coeffs_out[i*size_out] = -coeffs1[i*size1];
            }
         }
         else 
         {
            mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
            for (unsigned long i = 0; i < poly->length; i++)
            {
               coeffs_out[i*size_out] = coeffs1[i*size1];
            }
         }
      } else
      {
         mp_limb_t * signs = (mp_limb_t *) flint_malloc_limbs(poly->length);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            signs[i] = coeffs_out[i*size_out];
            if (!coeffs_out[i*size_out]) clear_limbs(coeffs_out+i*size_out+1, size_out-1);
            else coeffs_out[i*size_out] = 0L;
         }
         if (x < 0) 
         {
            mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, -x);
            for (unsigned long i = 0; i < poly->length; i++)
            {
               coeffs_out[i*size_out] = -signs[i];
            }
         }
         else 
         {
            mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
            for (unsigned long i = 0; i < poly->length; i++)
            {
               coeffs_out[i*size_out] = signs[i];
            }
         }
         flint_free(signs);
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
         mpn_divmod_1_preinv(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, x, xinv);
         if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out-size1);
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, size1-1, x);
         if (size_out > size1) clear_limbs(coeffs_out+i*size_out+size1, size_out-size1);
      }
   }
   
   output->length = poly->length;
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
