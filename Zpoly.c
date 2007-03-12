/****************************************************************************

Zpoly.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include <string.h>
#include "flint.h"
#include "Zpoly.h"


//////////////////////////////////////////////////////////
/*
Temporary baby implementations of memory management functions
*/


/*
This function gets called if malloc fails (i.e. returns NULL).
For now we just die via abort(); if we change this to do nothing, then
anyone who calls flint_malloc() etc will get a NULL pointer returned, and
perhaps get a chance to clean up before completely dying...
*/
void flint_memory_failure()
{
   abort();
}

void* flint_malloc_limbs(unsigned long limbs)
{
   void* buf = malloc(limbs * sizeof(mp_limb_t));
   if (!buf)
      flint_memory_failure();
   return buf;
}

void* flint_malloc(unsigned long bytes)
{
   void* buf = malloc(bytes);
   if (!buf)
      flint_memory_failure();
   return buf;
}

void* flint_realloc_limbs(void* block, unsigned long limbs)
{
   void* buf = realloc(block, limbs * sizeof(mp_limb_t));
   if (!buf)
      flint_memory_failure();
   return buf;
}

void* flint_realloc(void* block, unsigned long bytes)
{
   void* buf = realloc(block, bytes);
   if (!buf)
      flint_memory_failure();
   return buf;
}

void flint_free(void* block)
{
   free(block);
}

//////////////////////////////////////////////////////////


/****************************************************************************

   Zpoly_mpz_raw_* layer

****************************************************************************/


void Zpoly_mpz_raw_normalise(Zpoly_mpz_t poly)
{
   while (poly->length && !mpz_sgn(poly->coeffs[poly->length-1]))
      poly->length--;
}


void Zpoly_mpz_raw_set(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   FLINT_ASSERT(output->alloc >= input->length);
    
   for (unsigned long i = 0; i < input->length; i++)
      mpz_set(output->coeffs[i], input->coeffs[i]);
   
   output->length = input->length;
}


int Zpoly_mpz_raw_equal(Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   int is_input1_longer = (input1->length > input2->length);

   unsigned long short_length =
           is_input1_longer ? input2->length : input1->length;

   unsigned long i;
   for (i = 0; i < short_length; i++)
      if (mpz_cmp(input1->coeffs[i], input2->coeffs[i]))
         return 0;

   if (is_input1_longer)
      for (; i < input1->length; i++)
         if (mpz_sgn(input1->coeffs[i]))
            return 0;
   else
      for (; i < input2->length; i++)
         if (mpz_sgn(input2->coeffs[i]))
            return 0;

   return 1;
}


void Zpoly_mpz_raw_add(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   FLINT_ASSERT(output->alloc >= input1->length);
   FLINT_ASSERT(output->alloc >= input2->length);

   int is_input1_longer = (input1->length > input2->length);

   unsigned long short_length =
           is_input1_longer ? input2->length : input1->length;
   
   unsigned long i;
   for (i = 0; i < short_length; i++)
      mpz_add(output->coeffs[i], input1->coeffs[i], input2->coeffs[i]);
      
   if (is_input1_longer)
      for (; i < input1->length; i++)
         mpz_set(output->coeffs[i], input1->coeffs[i]);
   else
      for (; i < input2->length; i++)
         mpz_set(output->coeffs[i], input2->coeffs[i]);

   output->length = i;
}


void Zpoly_mpz_raw_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   FLINT_ASSERT(output->alloc >= input1->length);
   FLINT_ASSERT(output->alloc >= input2->length);

   int is_input1_longer = (input1->length > input2->length);

   unsigned long short_length =
           is_input1_longer ? input2->length : input1->length;
   
   unsigned long i;
   for (i = 0; i < short_length; i++)
      mpz_sub(output->coeffs[i], input1->coeffs[i], input2->coeffs[i]);
      
   if (is_input1_longer)
      for (; i < input1->length; i++)
         mpz_set(output->coeffs[i], input1->coeffs[i]);
   else
      for (; i < input2->length; i++)
         mpz_neg(output->coeffs[i], input2->coeffs[i]);

   output->length = i;
}


void Zpoly_mpz_raw_negate(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   FLINT_ASSERT(output->alloc >= input->length);
   
   for (unsigned long i = 0; i < input->length; i++)
      mpz_neg(output->coeffs[i], input->coeffs[i]);
   
   output->length = input->length;
}


void Zpoly_mpz_raw_scalar_mul(Zpoly_mpz_t poly, mpz_t x)
{
   abort();
}

void Zpoly_mpz_raw_scalar_mul_ui(Zpoly_mpz_t poly, unsigned long x)
{
   abort();
}

void Zpoly_mpz_raw_scalar_mul_si(Zpoly_mpz_t poly, long x)
{
   abort();
}

void Zpoly_mpz_raw_scalar_div(Zpoly_mpz_t poly, mpz_t x)
{
   abort();
}

void Zpoly_mpz_raw_scalar_div_ui(Zpoly_mpz_t poly, unsigned long x)
{
   abort();
}

void Zpoly_mpz_raw_mul(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   abort();
}
                           
void Zpoly_mpz_raw_mul_naive(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                             Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_raw_mul_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                                 Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_raw_sqr(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_raw_sqr_naive(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_raw_sqr_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_raw_left_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                              unsigned long n)
{
   abort();
}

void Zpoly_mpz_raw_right_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                               unsigned long n)
{
   abort();
}


void Zpoly_mpz_raw_div(Zpoly_mpz_t quotient, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_raw_rem(Zpoly_mpz_t remainder, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_raw_div_rem(Zpoly_mpz_t quotient, Zpoly_mpz_t remainder,
                           Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_raw_gcd(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_raw_xgcd(Zpoly_mpz_t a, Zpoly_mpz_t b, Zpoly_mpz_t output,
                        Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   abort();
}



/****************************************************************************

   Zpoly_mpz_* layer  ("non-raw" functions)

****************************************************************************/


void Zpoly_mpz_init(Zpoly_mpz_t poly)
{
   poly->coeffs = (mpz_t*) flint_malloc(sizeof(mpz_t));
   mpz_init(poly->coeffs[0]);
   poly->alloc = 1;
   poly->length = 0;
}


void Zpoly_mpz_init2(Zpoly_mpz_t poly, unsigned long alloc)
{
   poly->coeffs = (mpz_t*) flint_malloc(sizeof(mpz_t) * alloc);
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init(poly->coeffs[i]);
   poly->alloc = alloc;
   poly->length = 0;
}


void Zpoly_mpz_init3(Zpoly_mpz_t poly, unsigned long alloc,
                     unsigned long coeff_bits)
{
   poly->coeffs = (mpz_t*) flint_malloc(sizeof(mpz_t) * alloc);
   for (unsigned long i = 0; i < alloc; i++)
      mpz_init2(poly->coeffs[i], coeff_bits);
   poly->alloc = alloc;
   poly->length = 0;
}


void Zpoly_mpz_realloc(Zpoly_mpz_t poly, unsigned long alloc)
{
   if (alloc < poly->alloc)
   {
      // shrink available space; need to mpz_clear mpz's in the tail
      for (unsigned long i = alloc; i < poly->alloc; i++)
         mpz_clear(poly->coeffs[i]);
   }
   
   poly->coeffs = (mpz_t*) flint_realloc(poly->coeffs, sizeof(mpz_t) * alloc);

   // create new mpz's if necessary
   for (unsigned long i = poly->alloc; i < alloc; i++)
      mpz_init(poly->coeffs[i]);

   poly->alloc = alloc;
   
   // truncate actual data if necessary
   if (poly->length > alloc)
      poly->length = alloc;
}


void Zpoly_mpz_ensure_space_IMPL(Zpoly_mpz_t poly, unsigned long alloc)
{
   if (alloc < 2*poly->alloc)
      alloc = 2*poly->alloc;
   Zpoly_mpz_realloc(poly, alloc);
}


void Zpoly_mpz_clear(Zpoly_mpz_t poly)
{
   for (unsigned long i = 0; i < poly->alloc; i++)
      mpz_clear(poly->coeffs[i]);
   flint_free(poly->coeffs);
}


mpz_t* Zpoly_mpz_get_coeff_ptr(Zpoly_mpz_t poly, unsigned long n)
{
   if (n >= poly->length)
      return NULL;
   return &poly->coeffs[n];
}


void Zpoly_mpz_get_coeff(mpz_t output, Zpoly_mpz_t poly,
                         unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(output, 0);
   else
      mpz_set(output, poly->coeffs[n]);
}


int Zpoly_mpz_set_from_string(Zpoly_mpz_t output, char* s)
{
   const char* whitespace = " \t\n\r";
   
   output->length = 0;
   
   while (1)
   {
      // skip whitespace
      s += strspn(s, whitespace);
      if (*s == 0)
         return 1;
         
      // ensure sufficient space, and grab coefficient
      Zpoly_mpz_ensure_space(output, output->length + 1);
      if (!gmp_sscanf(s, "%Zd", output->coeffs[output->length]))
         return 0;
      output->length++;
      
      // jump to next whitespace
      s += strcspn(s, whitespace);
   }
}


unsigned long Zpoly_mpz_get_string_size(Zpoly_mpz_t poly)
{
   unsigned long size = 0;
   for (unsigned long i = 0; i < poly->length; i++)
      // (+2 is for the sign and a space)
      size += mpz_sizeinbase(poly->coeffs[i], 10) + 2;
   // (+1 is for the null terminator)
   return size + 1;
}


void Zpoly_mpz_get_as_string(char* output, Zpoly_mpz_t poly)
{
   if (!poly->length)
   {
      output[0] = 0;
      return;
   }
   
   for (unsigned long i = 0; i < poly->length; i++)
   {
      mpz_get_str(output, 10, poly->coeffs[i]);
      output += strlen(output);
      *output = ' ';
      output++;
   }
   
   output--;
   *output = ' ';
}


unsigned long Zpoly_mpz_get_coeff_ui(Zpoly_mpz_t poly, unsigned long n)
{
   abort();
}

long Zpoly_mpz_get_coeff_si(Zpoly_mpz_t poly, unsigned long n)
{
   abort();
}

void Zpoly_mpz_set_coeff(Zpoly_mpz_t poly, unsigned long n, mpz_t x)
{
   abort();
}

void Zpoly_mpz_set_coeff_ui(Zpoly_mpz_t poly, unsigned long n, unsigned long x)
{
   abort();
}

void Zpoly_mpz_set_coeff_si(Zpoly_mpz_t poly, unsigned long n, long x)
{
   abort();
}

void Zpoly_mpz_normalise(Zpoly_mpz_t poly)
{
   abort();
}

void Zpoly_mpz_set(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_swap(Zpoly_mpz_t x, Zpoly_mpz_t y)
{
   abort();
}

void Zpoly_mpz_add(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   Zpoly_mpz_ensure_space(output, FLINT_MAX(input1->length, input2->length));
   Zpoly_mpz_raw_add(output, input1, input2);
}

void Zpoly_mpz_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   Zpoly_mpz_ensure_space(output, FLINT_MAX(input1->length, input2->length));
   Zpoly_mpz_raw_sub(output, input1, input2);
}

void Zpoly_mpz_negate(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_scalar_mul(Zpoly_mpz_t poly, mpz_t x)
{
   abort();
}

void Zpoly_mpz_scalar_mul_ui(Zpoly_mpz_t poly, unsigned long x)
{
   abort();
}

void Zpoly_mpz_scalar_div(Zpoly_mpz_t poly, mpz_t x)
{
   abort();
}

void Zpoly_mpz_scalar_div_ui(Zpoly_mpz_t poly, unsigned long x)
{
   abort();
}

void Zpoly_mpz_mul(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_mul_naive(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                         Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_mul_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                             Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_sqr(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_sqr_naive(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_sqr_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
   abort();
}

void Zpoly_mpz_left_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                          unsigned long n)
{
   abort();
}

void Zpoly_mpz_right_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                           unsigned long n)
{
   abort();
}

void Zpoly_mpz_div(Zpoly_mpz_t quotient, Zpoly_mpz_t input1,
                   Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_rem(Zpoly_mpz_t remainder, Zpoly_mpz_t input1,
                   Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_div_rem(Zpoly_mpz_t quotient, Zpoly_mpz_t remainder,
                       Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_gcd(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   abort();
}

void Zpoly_mpz_xgcd(Zpoly_mpz_t a, Zpoly_mpz_t b, Zpoly_mpz_t output,
                    Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
   abort();
}

// *************** end of file
