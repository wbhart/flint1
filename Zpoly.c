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


mpz_t* Zpoly_mpz_raw_get_coeff_ptr(Zpoly_mpz_t poly, unsigned long n)
{
}

void Zpoly_mpz_raw_get_coeff(mpz_t output, Zpoly_mpz_t poly,
                             unsigned long n)
{
}

unsigned long Zpoly_mpz_raw_get_coeff_ui(Zpoly_mpz_t poly, unsigned long n)
{
}

void Zpoly_mpz_raw_set_coeff(Zpoly_mpz_t poly, unsigned long n, mpz_t x)
{
}

void Zpoly_mpz_raw_set_coeff_ui(Zpoly_mpz_t poly, unsigned long n,
                                unsigned long x)
{
}

void Zpoly_mpz_raw_set_coeff_si(Zpoly_mpz_t poly, unsigned long n, long x)
{
}

void Zpoly_mpz_raw_normalise(Zpoly_mpz_t poly)
{
}

void Zpoly_mpz_raw_set(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
}

void Zpoly_mpz_raw_swap(Zpoly_mpz_t x, Zpoly_mpz_t y)
{
}

void Zpoly_mpz_raw_add(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2)
{
}

void Zpoly_mpz_raw_negate(Zpoly_mpz_t output, Zpoly_mpz_t input)
{
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
}

void Zpoly_mpz_init2(Zpoly_mpz_t poly, unsigned long alloc)
{
}

void Zpoly_mpz_init3(Zpoly_mpz_t poly, unsigned long alloc,
                     unsigned long coeff_size)
{
}

void Zpoly_mpz_clear(Zpoly_mpz_t poly)
{
}

mpz_t* Zpoly_mpz_get_coeff_ptr(Zpoly_mpz_t poly, unsigned long n)
{
}

void Zpoly_mpz_get_coeff(mpz_t output, Zpoly_mpz_t poly,
                         unsigned long n)
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
