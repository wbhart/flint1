/****************************************************************************

fmpz_poly.h: Polynomials over Z, implemented as contiguous block of fmpz_t's

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_FMPZ_POLY_H
#define FLINT_FMPZ_POLY_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "mpn_extras.h"
#include "fmpz.h"

/****************************************************************************

   fmpz_poly_t
   -----------

fmpz_poly_t represents a dense polynomial in Z[x] using a single block of
memory to hold all the coefficients.

This type is better suited to handling very dense polynomials with relatively
small coefficients, where the memory management overhead of Zpoly_t would
be too expensive.

"coeffs" is an array of limbs of length (alloc * (limbs+1)). Each
coefficient uses limbs+1 limbs. For each coefficient, the first limb is
a sign/size limb: the number of limbs of the absolute value of the 
coefficient is given by the absolute value of this limb, and the sign 
of this limb is the sign of the coefficient. (Zero is stored as a 
sign/size of zero followed by arbitrary data.) The remaining "limbs" 
limbs represent the qbsolute value of the coefficient, stored in 
GMP's mpn format.

Only the first "length" coefficients actually represent coefficients 
of the polynomial; i.e. it's a polynomial of degree at most length-1. 
If length == 0, this is the zero polynomial. All functions normalise
so that the (length-1)-th coefficient is non-zero.
Obviously always alloc >= length.

There are two classes of functions operating on fmpz_poly_t:

-- The _fmpz_poly_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough space for
   the result. They also NEVER modify the limbs attribute (since this
   would screw up the block size).

-- The fmpz_poly_* functions ASSUME that "coeffs" was allocated via
   flint_heap_alloc, and they MAY free or reallocate "coeffs" using 
   flint_heap_realloc, flint_heap_free etc, whenever they feel the need.

*/
 
typedef struct
{
   mp_limb_t* coeffs;
   unsigned long alloc;
   unsigned long length;
   unsigned long limbs;
} fmpz_poly_struct;

// fmpz_poly_t allows reference-like semantics for fmpz_poly_struct:
typedef fmpz_poly_struct fmpz_poly_t[1];
typedef fmpz_poly_struct * fmpz_poly_p;

#define SWAP(x_dummy, y_dummy) \
do { \
   fmpz_poly_p swap_temp = x_dummy; \
   x_dummy = y_dummy; \
   y_dummy = swap_temp; \
} while(0);

#define SWAP_PTRS(x_dummy_p, y_dummy_p) \
do { \
   fmpz_t swap_temp_p = x_dummy_p; \
   x_dummy_p = y_dummy_p; \
   y_dummy_p = swap_temp_p; \
} while(0);

/****************************************************************************

   Conversion Routines
   
****************************************************************************/

/* 
   Converts fmpz_poly_t "poly_mpn" to a ZmodF_poly.
   
   Each coefficient of poly_mpn is assumed to fit into a coefficient 
   of poly_f.
   
   The maximum number of *bits* that any coefficient takes is returned
   and made negative if any of the coefficients was negative.
   
   Only _length_ coefficients are converted.
   
   Assumes 0 < length <= poly_mpn->length 
*/

long fmpz_poly_to_ZmodF_poly(ZmodF_poly_t poly_f, const fmpz_poly_t poly_mpn, 
                                                         const unsigned long length);


/* 
   Normalise and converts ZmodF_poly "poly_f" to a fmpz_poly_t. 
   
   Each coefficient of poly_f is assumed to fit into a coefficient 
   of poly_mpn. 
   
   The normalisation ensures that this function is the inverse of 
   ZmodF_poly_convert_in_mpn.
   
   Assumes 0 < poly_f->length 
*/

void ZmodF_poly_to_fmpz_poly(fmpz_poly_t poly_mpn, const ZmodF_poly_t poly_f, const long sign);

/* 
   Packs bundle coefficients, each padded out to the given number of limbs, into 
   the first coefficient of poly_f.
   
   Assumes 0 < bundle <= poly_mpn->length
*/
void fmpz_poly_limb_pack(ZmodF_poly_t poly_f, const fmpz_poly_t poly_mpn,
                                           const unsigned long bundle, const long limbs);

/* 
   Unpacks bundle coefficients from the first coefficient of poly_f, each 
   assumed to be stored in a field of the given number of limbs.
   
   Assumes 0 < bundle <= poly_f->length
*/
void fmpz_poly_limb_unpack(fmpz_poly_t poly_mpn, const ZmodF_poly_t poly_f, 
                                  const unsigned long bundle, const unsigned long limbs);

/* 
   Unpacks bundle coefficients from the first coefficient of poly_f, each 
   assumed to be stored in a field of the given number of limbs. Assumes the
   coefficients are unsigned.

   Assumes 0 < bundle <= poly_f->length
*/
void fmpz_poly_limb_unpack_unsigned(fmpz_poly_t poly_mpn, const ZmodF_poly_t poly_f, 
                                  const unsigned long bundle, const unsigned long limbs);


/*
   Packs poly_mpn down to the bit into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it. Each of the original
   coefficients is packed into a bitfield "bits" bits wide including one 
   bit for a sign bit. 
   
   "bits" is assumed to be less than FLINT_BITS. If bits is 
   negative, the input poly is assumed to have signed coefficients.

   Assumes 0 < bundle and 0 < poly_mpn->length
*/ 
   
void fmpz_poly_bit_pack(ZmodF_poly_t poly_f, const fmpz_poly_t poly_mpn,
     const unsigned long bundle, const long bits, const unsigned long length, const long negate);


/*
   Unpacks poly_f into poly_mpn. This is the inverse of ZmodF_poly_bitpack_mpn, so
   long as the final coefficient in the polynomial is positive.
   Each coeff of poly_f is assumed to contain "bundle" coefficients, each stored 
   in a bitfield "bits" bits wide with the most significant bit being reserved for
   the sign. 
   
   The total number of coefficients to be unpacked is given by the length of 
   poly_mpn. One must ensure each of the coefficients of poly_mpn are set to zero
   before calling this function for the first time since it adds to existing 
   coefficients of poly_mpn, rather than overwriting them.
   
   "bits" is assumed to be less than FLINT_BITS. 

   Assumes 0 < bundle and 0 < poly_f->length
*/ 
   
void fmpz_poly_bit_unpack(fmpz_poly_t poly_mpn, const ZmodF_poly_t poly_f, 
                              const unsigned long bundle, const unsigned long bits);
void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly_mpn, const ZmodF_poly_t poly_f, 
                              const unsigned long bundle, const unsigned long bits);


     
/*
   Packs poly_mpn down to the byte into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it, each packed into a field
   "bytes" bytes wide.
   
   "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
   coefficients are assumed to be at least a limb wide.

   Assumes 0 < bundle and 0 < poly_mpn->length
*/ 
   
void fmpz_poly_byte_pack(ZmodF_poly_t poly_f, const fmpz_poly_t poly_mpn,
                   const unsigned long bundle, const unsigned long coeff_bytes,
                                const unsigned long length, const long negate);

     
/*
   Unpacks array into poly_mpn. Each coefficient stored in array will have "bundle" 
   coefficients, each packed into a field "bytes" bytes wide.
   
   The total number of coefficients to be unpacked is given by the length of 
   poly_mpn.
   
   "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
   coefficients are assumed to be at least a limb wide.

   Assumes 0 < bundle and poly->mpn->length > 0
*/ 
   
void fmpz_poly_byte_unpack_unsigned(fmpz_poly_t poly_m, const mp_limb_t* array,
                               const unsigned long bundle, const unsigned long coeff_bytes);

void fmpz_poly_byte_unpack(fmpz_poly_t poly_m, const mp_limb_t* array,
                               const unsigned long bundle, const unsigned long coeff_bytes);


     
/*
   Splits each coefficient of poly_mpn into pieces "limbs" limbs long and 
   stores each piece into bundle coefficients of poly_f.
*/
 
void fmpz_poly_split(ZmodF_poly_t poly_f, fmpz_poly_t poly_mpn,
     unsigned long bundle, unsigned long limbs);

     
/*
   Combines each "bundle" coefficients of poly_f, each taken to be "limbs" 
   limbs long, into a coefficient of poly_mpn. 
   
   This function is used for testing purposed only, and is the exact inverse
   of ZmodF_poly_split_mpn.
   
   The number of coefficients extracted is given by the length of poly_mpn.
*/
 
void fmpz_poly_unsplit(ZmodF_poly_t poly_f, fmpz_poly_t poly_mpn,
     unsigned long bundle, unsigned long limbs);


/*============================================================================
  
    Functions in _fmpz_poly_* layer
    
===============================================================================*/

void _fmpz_poly_stack_init(fmpz_poly_t poly, const unsigned long alloc, const unsigned long limbs);

void _fmpz_poly_stack_clear(const fmpz_poly_t poly);

void _fmpz_poly_check(const fmpz_poly_t poly);

static inline
fmpz_t _fmpz_poly_get_coeff_ptr(fmpz_poly_t poly, const unsigned long n)
{
   return poly->coeffs+n*(poly->limbs+1);
}

/* 
   Set "output" to the given coefficient and return the sign
   Assumes length of output is poly->limbs limbs long.  
*/
   
static inline
long _fmpz_poly_get_coeff(fmpz_t output, const fmpz_poly_t poly,
                          const unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) F_mpn_clear(output, poly->limbs);
   F_mpn_copy(output, poly->coeffs+n*(poly->limbs+1)+1, poly->limbs);
   return poly->coeffs[n*(poly->limbs+1)];
}

static inline
unsigned long _fmpz_poly_get_coeff_ui(fmpz_poly_t poly, const unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   else return poly->coeffs[n*(poly->limbs+1)+1];
}

static inline
long _fmpz_poly_get_coeff_si(fmpz_poly_t poly, const unsigned long n)
{
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   if (poly->coeffs[n*(poly->limbs+1)] == 1L) 
                                 return poly->coeffs[n*(poly->limbs+1)+1];
   else return -poly->coeffs[n*(poly->limbs+1)+1];
}

void _fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, const unsigned long n);

void _fmpz_poly_normalise(fmpz_poly_t poly);

/* 
   Set a coefficient to the given value having "size" limbs.
   Assumes that the poly->limbs is at least "size" and
   that n < poly->length
*/

static inline void _fmpz_poly_set_coeff(fmpz_poly_t poly, const unsigned long n, 
                                  const fmpz_t x, const long sign, const unsigned long size)
{
   FLINT_ASSERT(poly->limbs >= size);
   F_mpn_copy(poly->coeffs+n*(poly->limbs+1)+1, x, size);
   poly->coeffs[n*(poly->limbs+1)] = sign;
   if (poly->limbs > size) 
     F_mpn_clear(poly->coeffs+n*(poly->limbs+1)+size+1, poly->limbs-size);
   if ((x[0] == 0) && (n == poly->length - 1)) _fmpz_poly_normalise(poly);
}

void _fmpz_poly_set_coeff_ui(fmpz_poly_t poly, const unsigned long n, const unsigned long x);

void _fmpz_poly_set_coeff_si(fmpz_poly_t poly, const unsigned long n, const long x);

static inline long _fmpz_poly_degree(const fmpz_poly_t poly)
{
   return poly->length - 1;
}

static inline unsigned long _fmpz_poly_length(const fmpz_poly_t poly)
{
   return poly->length;
}

static inline unsigned long _fmpz_poly_limbs(const fmpz_poly_t poly)
{
   return poly->limbs;
}

void _fmpz_poly_set(fmpz_poly_t output, const fmpz_poly_t input);

/* 
   Zero the polynomial by setting the length to zero.
   Does not set the actual limbs to zero.
*/

static inline void _fmpz_poly_zero(fmpz_poly_t output)
{
   output->length = 0;
}

/* Zero first n coefficients of poly, regardless of what length is */

void _fmpz_poly_zero_coeffs(fmpz_poly_t poly, const unsigned long n);

static inline void _fmpz_poly_transfer(fmpz_poly_t output, const fmpz_poly_t input)
{
   output->length = input->length;
   output->limbs = input->limbs;
   output->coeffs = input->coeffs;
}

long _fmpz_poly_bits1(const fmpz_poly_t poly_mpn);

long _fmpz_poly_bits(const fmpz_poly_t poly_mpn);

unsigned long _fmpz_poly_max_limbs(const fmpz_poly_t poly);

int _fmpz_poly_equal(const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_neg(fmpz_poly_t output, const fmpz_poly_t input);

void _fmpz_poly_truncate(fmpz_poly_t poly, const unsigned long trunc);

void _fmpz_poly_reverse(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long length);

void _fmpz_poly_left_shift(fmpz_poly_t output, const fmpz_poly_t input, 
                                                 const unsigned long n);

void _fmpz_poly_right_shift(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long n);

void _fmpz_poly_add(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_sub(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_scalar_mul(fmpz_poly_t output, const fmpz_poly_t poly, const fmpz_t x);

void _fmpz_poly_scalar_mul_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x);

void _fmpz_poly_scalar_mul_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x);

void _fmpz_poly_scalar_div(fmpz_poly_t output, const fmpz_poly_t poly, const fmpz_t x);

void _fmpz_poly_scalar_div_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x);

void _fmpz_poly_scalar_div_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x);

void _fmpz_poly_scalar_div_exact_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x);

void _fmpz_poly_scalar_div_exact_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x);

void _fmpz_poly_mul_classical(fmpz_poly_t output, const fmpz_poly_t input1, 
                                                 const fmpz_poly_t input2);
                                                 
void _fmpz_poly_mul_classical_trunc(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc);
                                          
void _fmpz_poly_mul_classical_trunc_left(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc);

void __fmpz_poly_karamul_recursive(fmpz_poly_t res, const fmpz_poly_t a, const fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, const unsigned long crossover);

void _fmpz_poly_mul_karatsuba(fmpz_poly_t output, const fmpz_poly_t input1, 
                                                 const fmpz_poly_t input2);
                                                 
void _fmpz_poly_mul_karatsuba_trunc(fmpz_poly_t output, const fmpz_poly_t input1, 
                                           const fmpz_poly_t input2, const unsigned long trunc);
                                                 
void _fmpz_poly_mul_karatsuba_trunc_left(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2, const unsigned long trunc);

/*
   Multiply two polynomials together using the Kronecker segmentation method.
   Currently assumes that the number of output bits per coefficient is <= 64 and
   is supplied by the parameter "bits"
*/

void _fmpz_poly_mul_KS(fmpz_poly_t output, const fmpz_poly_t input1, 
                                       const fmpz_poly_t input2);
                                       
void _fmpz_poly_mul_KS_trunc(fmpz_poly_t output, const fmpz_poly_t input1, 
                                        const fmpz_poly_t input2, const unsigned long trunc);

void _fmpz_poly_mul_SS(fmpz_poly_t output, const fmpz_poly_t input1, 
                                                      const fmpz_poly_t input2);
                                                            
void _fmpz_poly_mul_SS_trunc(fmpz_poly_t output, const fmpz_poly_t input1, 
                           const fmpz_poly_t input2, const unsigned long trunc);
                                        
void _fmpz_poly_mul_trunc_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                const fmpz_poly_t input2, const unsigned long trunc);
                                
void _fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                const fmpz_poly_t input2, const unsigned long trunc);
                                
void _fmpz_poly_mul(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_sqr(fmpz_poly_t output, const fmpz_poly_t input);

void _fmpz_poly_sqr_naive(fmpz_poly_t output, const fmpz_poly_t input);

void _fmpz_poly_sqr_karatsuba(fmpz_poly_t output, const fmpz_poly_t input);

void _fmpz_poly_content(fmpz_t content, const fmpz_poly_t a);

/*============================================================================
  
    Functions in fmpz_poly_* layer
    
===============================================================================*/

void fmpz_poly_init(fmpz_poly_t poly);

void fmpz_poly_init2(fmpz_poly_t poly, const unsigned long alloc, const unsigned long limbs);
                                              
void fmpz_poly_realloc(fmpz_poly_t poly, const unsigned long alloc);

void fmpz_poly_fit_length(fmpz_poly_t poly, const unsigned long alloc);

void fmpz_poly_resize_limbs(fmpz_poly_t poly, const unsigned long limbs);

static inline void fmpz_poly_fit_limbs(fmpz_poly_t poly, const unsigned long limbs)
{
   if ((long) limbs > (long) poly->limbs) fmpz_poly_resize_limbs(poly, limbs);
}

void fmpz_poly_clear(fmpz_poly_t poly);

void fmpz_poly_check(const fmpz_poly_t poly);

// ------------------------------------------------------
// String conversions and I/O

int fmpz_poly_from_string(fmpz_poly_t poly, const char* s);
char* fmpz_poly_to_string(const fmpz_poly_t poly);
void fmpz_poly_print(const fmpz_poly_t poly);
void fmpz_poly_fprint(const fmpz_poly_t poly, FILE* f);
int fmpz_poly_fread(fmpz_poly_t poly, FILE* f);

static inline unsigned long fmpz_poly_limbs(const fmpz_poly_t poly)
{
   return poly->limbs;
}

static inline 
int fmpz_poly_read(fmpz_poly_t poly)
{
   return fmpz_poly_fread(poly, stdin);
}

static inline 
long fmpz_poly_degree(const fmpz_poly_t poly)
{
   return poly->length - 1;
}

static inline 
unsigned long fmpz_poly_length(const fmpz_poly_t poly)
{
   return poly->length;
}

static inline 
long fmpz_poly_bits(const fmpz_poly_t poly)
{
   return _fmpz_poly_bits(poly);
}

static inline
long fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
   return _fmpz_poly_max_limbs(poly);
}

static inline void fmpz_poly_truncate(fmpz_poly_t poly, const unsigned long length)
{
   FLINT_ASSERT(poly->length >= length);
   poly->length = length;
   _fmpz_poly_normalise(poly);
}

void fmpz_poly_swap(fmpz_poly_t x, fmpz_poly_t y);

static inline
fmpz_t fmpz_poly_get_coeff_ptr(const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length)
   {
      return NULL;
   }
   return poly->coeffs+n*(poly->limbs+1);
}

static inline
unsigned long fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length)
   {
      return 0;
   }
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   else return poly->coeffs[n*(poly->limbs+1)+1];
}

static inline
long fmpz_poly_get_coeff_si(const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length)
   {
      return 0;
   }
   if (poly->coeffs[n*(poly->limbs+1)] == 0) return 0;
   if (poly->coeffs[n*(poly->limbs+1)] == 1L) 
                                 return poly->coeffs[n*(poly->limbs+1)+1];
   else return -poly->coeffs[n*(poly->limbs+1)+1];
}

void fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, const unsigned long n);

static inline void fmpz_poly_set_coeff_si(fmpz_poly_t poly, const unsigned long n, const long x)
{
   fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, 1);
   if (n+1 > poly->length) 
   {
      for (long i = poly->length; i + 1 < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0L;
      } 
      poly->length = n+1;
   }
   _fmpz_poly_set_coeff_si(poly, n, x);
   _fmpz_poly_normalise(poly);
}

static inline void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, const unsigned long n, const unsigned long x)
{
   fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, 1);
   if (n+1 > poly->length) 
   {
      for (long i = poly->length; i + 1 < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0L;
      } 
      poly->length = n+1;
   }
   _fmpz_poly_set_coeff_ui(poly, n, x);
   _fmpz_poly_normalise(poly);
}

static inline
void fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, const unsigned long n, const mpz_t x)
{
   fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, mpz_size(x));
   
   if (poly->limbs == 0) return;
   
   if (n+1 > poly->length) 
   {
      for (long i = poly->length; i + 1 < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0;
      } 
      poly->length = n+1;
   }
   mpz_to_fmpz(poly->coeffs + n*(poly->limbs+1), x);
   _fmpz_poly_normalise(poly);
}

static inline
void fmpz_poly_set(fmpz_poly_t output, const fmpz_poly_t input)
{
   fmpz_poly_fit_length(output, input->length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_set(output, input);
}

static inline void fmpz_poly_zero(fmpz_poly_t output)
{
   output->length = 0;
}

static inline
void fmpz_poly_neg(fmpz_poly_t output, const fmpz_poly_t input)
{
   fmpz_poly_fit_length(output, input->length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_neg(output, input);
}

static inline
void fmpz_poly_reverse(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long length)
{
   fmpz_poly_fit_length(output, length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_reverse(output, input, length);
}

void fmpz_poly_scalar_mul(fmpz_poly_t output, 
                          const fmpz_poly_t input, const fmpz_t x);
                          
void fmpz_poly_scalar_mul_ui(fmpz_poly_t output, 
                          const fmpz_poly_t input, unsigned long x);
                          
void fmpz_poly_scalar_mul_si(fmpz_poly_t output, 
                          const fmpz_poly_t input, long x);

static inline
void fmpz_poly_add(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   unsigned long bits1 = FLINT_ABS(_fmpz_poly_bits(input1));
   unsigned long bits2 = FLINT_ABS(_fmpz_poly_bits(input2));
   
   fmpz_poly_fit_length(output, FLINT_MAX(input1->length, input2->length));
   fmpz_poly_fit_limbs(output, FLINT_MAX(bits1, bits2)/FLINT_BITS + 1);
   
   _fmpz_poly_add(output, input1, input2);
}

static inline
void fmpz_poly_sub(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   unsigned long bits1 = FLINT_ABS(_fmpz_poly_bits(input1));
   unsigned long bits2 = FLINT_ABS(_fmpz_poly_bits(input2));
   
   fmpz_poly_fit_length(output, FLINT_MAX(input1->length, input2->length));
   fmpz_poly_fit_limbs(output, FLINT_MAX(bits1, bits2)/FLINT_BITS + 1);
   
   _fmpz_poly_sub(output, input1, input2);
}

void fmpz_poly_mul(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void fmpz_poly_mul_trunc_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc);
                                          
void fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc);

void fmpz_poly_div_classical(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_divrem_classical(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_div_bisection_recursive(fmpz_poly_t Q, fmpz_poly_t DQ, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_divrem_bisection(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_div_bisection(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_div_mulders(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_newton_invert_basecase(fmpz_poly_t Q_inv, const fmpz_poly_t Q, const unsigned long n);

void fmpz_poly_newton_invert(fmpz_poly_t Q_inv, const fmpz_poly_t Q, const unsigned long n);

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, const unsigned long n);

void fmpz_poly_div_newton(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

static inline
void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_div_mulders(Q, A, B);
}

static inline
void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_divrem_bisection(Q, R, A, B);
}

void fmpz_poly_power(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long exp);

void fmpz_poly_power_trunc_n(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long exp, const unsigned long n);

void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_pseudo_divrem_schoup(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                               unsigned long * d, const fmpz_poly_t A, const fmpz_poly_t B);
                               
void fmpz_poly_pseudo_div_basecase(fmpz_poly_t Q, unsigned long * d, 
                                               const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_pseudo_divrem_recursive(fmpz_poly_t Q, fmpz_poly_t R, 
                               unsigned long * d, const fmpz_poly_t A, const fmpz_poly_t B);

static inline
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, unsigned long * d, 
                                               const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_pseudo_divrem_recursive(Q, R, d, A, B);
}

void fmpz_poly_pseudo_div_recursive(fmpz_poly_t Q, unsigned long * d, 
                                      const fmpz_poly_t A, const fmpz_poly_t B);
                                      
static inline
void fmpz_poly_pseudo_div(fmpz_poly_t Q, unsigned long * d, 
                                      const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_pseudo_div_recursive(Q, d, A, B);
}
                                      


// *************** end of file

#ifdef __cplusplus
 }
#endif
 
#endif
