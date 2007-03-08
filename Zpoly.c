/****************************************************************************

Zpoly.c: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

!!!!! this code is not even at the "compile" stage.... it's strictly for
amusement for the moment

*****************************************************************************/


//////////////////////////////////////////////////////////
/*
Temporary baby implementations of memory management functions
*/

void* flint_malloc_limbs(unsigned long limbs)
{
    return malloc(limbs * sizeof(mp_limb_t));
}

void* flint_malloc(unsigned long bytes)
{
    return malloc(bytes);
}

void* flint_realloc_limbs(void* block, unsigned long limbs)
{
    return realloc(block, limbs * sizeof(mp_limb_t));
}

void* flint_realloc(void* block, unsigned long bytes)
{
    return realloc(block, bytes);
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
}


void Zpoly_mpz_raw_add(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   FLINT_ASSERT(output->alloc >= input1->length);
   FLINT_ASSERT(output->alloc >= input2->length);

   Zpoly_mpz_struct* short_poly, long_poly;
   if (input1->length < input2->length)
   {
      short_poly = input1;
      long_poly = input2;
   }
   else
   {
      short_poly = input2;
      long_poly = input1;
   }
    
   unsigned long i;
   for (i = 0; i < short_poly->length; i++)
      mpz_add(output->coeffs[i], input1->coeffs[i], input2->coeffs[i]);
   for (; i < long_poly->length; i++)
      mpz_set(output->coeffs[i], long_poly->coeffs[i]);

   output->length = long_poly->length;
}


void Zpoly_mpz_raw_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
   FLINT_ASSERT(output->alloc >= input1->length);
   FLINT_ASSERT(output->alloc >= input2->length);
   
   Zpoly_mpz_struct* short_poly, long_poly;
   if (input1->length < input2->length)
   {
      short_poly = input1;
      long_poly = input2;
   }
   else
   {
      short_poly = input2;
      long_poly = input1;
   }
   
   unsigned long i;
   for (i = 0; i < short_poly->length; i++)
      mpz_sub(output->coeffs[i], input1->coeffs[i], input2->coeffs[i]);
   if (long_poly == input1)
   {
      for (; i < long_poly->length; i++)
         mpz_set(output->coeffs[i], long_poly->coeffs[i]);
   }
   else
   {
      for (; i < long_poly->length; i++)
         mpz_neg(output->coeffs[i], long_poly->coeffs[i]);
   }
   
   output->length = long_poly->length;
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
}

void Zpoly_mpz_raw_scalar_mul_ui(Zpoly_mpz_t poly, unsigned long x)
{
}

void Zpoly_mpz_raw_scalar_div(Zpoly_mpz_t poly, mpz_t x)
{
}

void Zpoly_mpz_raw_scalar_div_ui(Zpoly_mpz_t poly, unsigned long x)
{
}

void Zpoly_mpz_raw_mul(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}
                           
void Zpoly_mpz_raw_mul_naive(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                             Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_mul_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                                 Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_sqr(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_raw_sqr_naive(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_raw_sqr_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_raw_left_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                              unsigned long n)
{
}

void Zpoly_mpz_raw_right_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                               unsigned long n)
{
}


void Zpoly_mpz_raw_div(Zpoly_mpz_t quotient, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_rem(Zpoly_mpz_t remainder, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_div_rem(Zpoly_mpz_t quotient, Zpoly_mpz_t remainder,
                           Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_gcd(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_xgcd(Zpoly_mpz_t a, Zpoly_mpz_t b, Zpoly_mpz_t output,
                        Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
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
      // shrink available space; need to clear mpz's in the tail
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
   return poly->coeffs[n];
}

void Zpoly_mpz_get_coeff(mpz_t output, Zpoly_mpz_t poly,
                         unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(output, 0);
   else
      mpz_set(output, poly->coeffs[n]);
}

void Zpoly_mpz_set_from_string(Zpoly_mpz_t output, char* s)
{
}

unsigned long Zpoly_mpz_get_coeff_ui(Zpoly_mpz_t poly, unsigned long n)
{
}

void Zpoly_mpz_set_coeff(Zpoly_mpz_t poly, unsigned long n, mpz_t x)
{
}

void Zpoly_mpz_set_coeff_ui(Zpoly_mpz_t poly, unsigned long n, unsigned long x)
{
}

void Zpoly_mpz_set_coeff_si(Zpoly_mpz_t poly, unsigned long n, long x)
{
}

void Zpoly_mpz_normalise(Zpoly_mpz_t poly)
{
}

void Zpoly_mpz_set(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_swap(Zpoly_mpz_t x, Zpoly_mpz_t y)
{
}

void Zpoly_mpz_add(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_negate(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_scalar_mul(Zpoly_mpz_t poly, mpz_t x)
{
}

void Zpoly_mpz_scalar_mul_ui(Zpoly_mpz_t poly, unsigned long x)
{
}

void Zpoly_mpz_scalar_div(Zpoly_mpz_t poly, mpz_t x)
{
}

void Zpoly_mpz_scalar_div_ui(Zpoly_mpz_t poly, unsigned long x)
{
}

void Zpoly_mpz_mul(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_mul_naive(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                         Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_mul_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                             Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_sqr(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_sqr_naive(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_sqr_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_left_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                          unsigned long n)
{
}

void Zpoly_mpz_right_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                           unsigned long n)
{
}

void Zpoly_mpz_div(Zpoly_mpz_t quotient, Zpoly_mpz_t input1,
                   Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_rem(Zpoly_mpz_t remainder, Zpoly_mpz_t input1,
                   Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_div_rem(Zpoly_mpz_t quotient, Zpoly_mpz_t remainder,
                       Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_gcd(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_xgcd(Zpoly_mpz_t a, Zpoly_mpz_t b, Zpoly_mpz_t output,
                    Zpoly_mpz_t input1, Zpoly_mpz_t input2)
{
}

// end of file
