/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
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
#include "zmod_poly.h"

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
limbs represent the absolute value of the coefficient, stored in 
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
   Converts fmpz_poly_t "poly_fmpz" to a ZmodF_poly.
   
   Each coefficient of poly_fmpz is assumed to fit into a coefficient 
   of poly_f.
   
   The maximum number of *bits* that any coefficient takes is returned
   and made negative if any of the coefficients was negative.
   
   Only _length_ coefficients are converted.
   
   Assumes 0 < length <= poly_fmpz->length 
*/

long fmpz_poly_to_ZmodF_poly(ZmodF_poly_t poly_f, const fmpz_poly_t poly_fmpz, 
                                                         const unsigned long length);

/* 
   Normalise and converts ZmodF_poly "poly_f" to a fmpz_poly_t. 
   
   Each coefficient of poly_f is assumed to fit into a coefficient 
   of poly_fmpz. 
   
   The normalisation ensures that this function is the inverse of 
   ZmodF_poly_convert_in_mpn.
   
   Assumes 0 < poly_f->length 
*/

void ZmodF_poly_to_fmpz_poly(fmpz_poly_t poly_fmpz, const ZmodF_poly_t poly_f, const long sign);

/* 
   Packs poly_fmpz->length coefficients, each assumed to fit into a limb, into
   the array of limbs coeffs_f. Coefficients are stored in B-adic symmetric 
	format. The neg variant changes the sign of the coefficients before packing 
	(the input polynomial is not modified). 
   
   Assumes 0 < poly_fmpz->length
*/
void fmpz_poly_limb_pack_1(mp_limb_t * array, const fmpz_poly_t poly_fmpz);
void fmpz_poly_limb_pack_neg_1(mp_limb_t * array, const fmpz_poly_t poly_fmpz);

/* 
   Unpacks length signed coefficients from the array, each assumed to be 
	stored in precisely one limb, into the polynomial poly_fmpz. 
	Coefficients are assumed to be stored in B-adic symmetric format. 
   
   Assumes 0 < length <= poly_f->length
*/
void fmpz_poly_limb_unpack_1(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                                  const unsigned long length);

/* 
   Packs length coefficients of poly_fmpz into the array into fields with the 
	given number of limbs. The neg variant negates the coefficients.
   
   Assumes 0 < length.
*/
void fmpz_poly_limb_pack(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
                                           const unsigned long length, const long limbs);
void fmpz_poly_limb_pack_neg(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
                                           const unsigned long length, const long limbs);

/* 
   Unpacks length signed coefficients from the array, each assumed to be 
	stored in a field of the given number of limbs, into the polynomial 
	poly_fmpz. Coefficients are assumed to be stored in B-adic symmetric 
	format. 
   
   Assumes 0 < length <= poly_f->length
*/
void fmpz_poly_limb_unpack(fmpz_poly_t poly_fmpz, mp_limb_t * array, 
                                  const unsigned long length, const unsigned long limbs);

/* 
   Unpacks bundle coefficients from the given array, each assumed to be 
	stored in a field of the given number of limbs. Assumes the coefficients 
	are unsigned.

   Assumes 0 < bundle <= poly_f->length
*/
void fmpz_poly_limb_unpack_unsigned(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                                  const unsigned long length, const unsigned long limbs);


/*
   Packs length coefficients of poly_fmpz down to the bit into the given 
	array. Each coefficient is packed into a bitfield "bits" bits wide 
	including one bit for a sign bit. 
   
   "bits" is assumed to be less than FLINT_BITS. If bits is 
   negative, the input poly is assumed to have signed coefficients.

   Assumes 0 < bundle and 0 < poly_fmpz->length
*/ 
   
void fmpz_poly_bit_pack(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
     const unsigned long length, const long bits, const long negate);

/*
   Unpacks the given array into poly_fmpz. This is the inverse of fmpz_poly_bit_pack, 
	however array must have one more (zero) limb than is necessary to store all the
	packed coefficients. 
   
   The total number of coefficients to be unpacked is given by length. One must 
	ensure each of the coefficients of poly_fmpz are set to zero before calling 
	this function for the first time since it adds to existing coefficients of 
	poly_fmpz, rather than overwriting them.
   
	The unsigned variant assumes the coefficients are unsigned, otherwise the bitfield
	is assumed to be one bit wider to allow for a sign bit.

   "bits" is assumed to be less than FLINT_BITS. 

   Assumes 0 < bundle and 0 < poly_f->length
*/ 
   
void fmpz_poly_bit_unpack(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                              const unsigned long length, const unsigned long bits);
void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                              const unsigned long length, const unsigned long bits);
     
/*
   Packs coefficients of poly_fmpz down to the byte into array, each packed 
	into a field "bytes" bytes wide.
   
   "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
   coefficients are assumed to be at least a limb wide.

   Assumes 0 < length and 0 < poly_fmpz->length
*/ 
   
void fmpz_poly_byte_pack(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
                   const unsigned long length, const unsigned long coeff_bytes,
                                const long negate);

     
/*
   Unpacks coefficients from array into poly_fmpz. Each coefficient stored is
	assumed to be packed into a field "bytes" bytes wide, with one bit reserved
	for a sign bit except for the unsigned variant of the function where the 
	coefficients are assumed to be unsigned. 
	
	It is also assumed that array has one extra (zero) limb beyond what is 
	required to store the packed coefficients.
   
   The total number of coefficients to be unpacked is given by length.
   
   "coeff_bytes" is assumed to be at least FLINT_BITS/8, i.e. the
   coefficients are assumed to be at least a limb wide.

   Assumes 0 < length.
*/ 
   
void fmpz_poly_byte_unpack_unsigned(fmpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long bytes);

void fmpz_poly_byte_unpack(fmpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long bytes);

void fmpz_poly_pack_bytes(fmpz_poly_t res, fmpz_poly_t poly, ulong n, ulong bytes);

void fmpz_poly_unpack_bytes(fmpz_poly_t res, fmpz_poly_t poly, ulong n, ulong bytes);

/*
   Reduce coefficients of the given fmpz_poly fpol, modulo the modulus of the given 
   zmod_poly zpol, and store the result in zpol.
*/

void fmpz_poly_to_zmod_poly(zmod_poly_t zpol, const fmpz_poly_t fpol);

void fmpz_poly_to_zmod_poly_no_red(zmod_poly_t zpol, const fmpz_poly_t fpol);

/*
   Store the unsigned long coefficients of the zmod_poly zpol in the given fmpz_poly fpol.
   The unsigned version normalised to [0, p) the other version to [-(p-1)/2, p/2]
*/

void zmod_poly_to_fmpz_poly(fmpz_poly_t fpol, zmod_poly_t zpol);

void zmod_poly_to_fmpz_poly_unsigned(fmpz_poly_t fpol, zmod_poly_t zpol);

/*
   Given an fmpz_poly_t fpol representing the reduction modulo oldmod of a polynomial 
   and a zmod_poly zpol with modulus p, use Chinese remaindering to reconstruct the 
   polynomial modulo newmod, with newmod = p*oldmod, where each new coefficient fpol[i] 
   is set to the unique non-negative integer in [0, newmod) which is fpol[i] mod oldmod 
   and zpol[i] mod p. 

   Assumes p and oldmod are coprime.

   Returns 1 if the CRT has stabilised, i.e. if the new output equals the old input
   else it returns 0.

   The unsigned version normalises the output to [0, newmod) the main version
   normalises to [-nemod/2, newmod/2].

   Allows aliasing of fpol and res.
*/

int fmpz_poly_CRT(fmpz_poly_t res, fmpz_poly_t fpol, zmod_poly_t zpol, fmpz_t newmod, fmpz_t oldmod);

int fmpz_poly_CRT_unsigned(fmpz_poly_t res, fmpz_poly_t fpol, zmod_poly_t zpol, fmpz_t newmod, fmpz_t oldmod);

/*============================================================================
  
    Functions in _fmpz_poly_* layer
    
===============================================================================*/

void _fmpz_poly_stack_init(fmpz_poly_t poly, const unsigned long alloc, const unsigned long limbs);

void _fmpz_poly_stack_clear(fmpz_poly_t poly);

void _fmpz_poly_check(const fmpz_poly_t poly);

void _fmpz_poly_normalise(fmpz_poly_t poly);

void _fmpz_poly_check_normalisation(const fmpz_poly_t poly);

static inline
fmpz_t _fmpz_poly_get_coeff_ptr(const fmpz_poly_t poly, const unsigned long n)
{
   return poly->coeffs+n*(poly->limbs+1);
}

static inline
fmpz_t _fmpz_poly_lead(const fmpz_poly_t poly)
{
   if (poly->length == 0) return NULL;
   return poly->coeffs+(poly->length-1)*(poly->limbs+1);
}

static inline
fmpz_t fmpz_poly_lead(const fmpz_poly_t poly)
{
   if (poly->length == 0) return NULL;
   return poly->coeffs+(poly->length-1)*(poly->limbs+1);
}

/* 
   Set "output" to the given coefficient and return the sign
   Assumes length of output is poly->limbs limbs long.  
*/
   
static inline
long _fmpz_poly_get_coeff(mp_limb_t * output, const fmpz_poly_t poly,
                          const unsigned long n)
{
   F_mpn_clear(output, poly->limbs);
   if (poly->coeffs[n*(poly->limbs+1)] != 0L)  
   {
      F_mpn_copy(output, poly->coeffs+n*(poly->limbs+1)+1, ABS(poly->coeffs[n*(poly->limbs+1)]));
   }
      
   return poly->coeffs[n*(poly->limbs+1)];
}

static inline
unsigned long _fmpz_poly_get_coeff_ui(fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length) return 0L;
   if (poly->coeffs[n*(poly->limbs+1)] == 0L) return 0L;
   else return poly->coeffs[n*(poly->limbs+1)+1];
}

static inline
long _fmpz_poly_get_coeff_si(fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length) return 0L;
   if (poly->coeffs[n*(poly->limbs+1)] == 0L) return 0L;
   if ((long) poly->coeffs[n*(poly->limbs+1)] > 0L) 
                                 return poly->coeffs[n*(poly->limbs+1)+1];
   else return -poly->coeffs[n*(poly->limbs+1)+1];
}

void _fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, const unsigned long n);

void _fmpz_poly_get_coeff_mpz_read_only(mpz_t x, const fmpz_poly_t poly, const unsigned long n);

static inline
void _fmpz_poly_set_coeff_mpz(fmpz_poly_t poly, const unsigned long n, const mpz_t x)
{
   if (poly->limbs == 0) return;
   
   if (n+1 > poly->length) 
   {
      long i;
		for (i = poly->length; i + 1 < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0;
      } 
      poly->length = n+1;
   }
   mpz_to_fmpz(poly->coeffs + n*(poly->limbs+1), x);
   _fmpz_poly_normalise(poly);
}

void _fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, const unsigned long n, fmpz_t x);

void _fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, const unsigned long n);

/* 
   Set a coefficient to the given value having "size" limbs.
   Assumes that the poly->limbs is at least "size" and
   that n < poly->length
*/

static inline 
void _fmpz_poly_set_coeff(fmpz_poly_t poly, const unsigned long n, 
                                  const mp_limb_t * x, const long sign, const unsigned long size)
{
   FLINT_ASSERT(poly->limbs >= size);
   F_mpn_copy(poly->coeffs+n*(poly->limbs+1)+1, x, size);
   poly->coeffs[n*(poly->limbs+1)] = sign;
   if (poly->limbs > size) 
     F_mpn_clear(poly->coeffs+n*(poly->limbs+1)+size+1, poly->limbs-size);
   _fmpz_poly_normalise(poly);
}

void _fmpz_poly_set_coeff_ui(fmpz_poly_t poly, const unsigned long n, const unsigned long x);

void _fmpz_poly_set_coeff_si(fmpz_poly_t poly, const unsigned long n, const long x);

static inline 
long _fmpz_poly_degree(const fmpz_poly_t poly)
{
   return poly->length - 1;
}

static inline 
unsigned long _fmpz_poly_length(const fmpz_poly_t poly)
{
   return poly->length;
}

static inline 
unsigned long _fmpz_poly_limbs(const fmpz_poly_t poly)
{
   return poly->limbs;
}

void _fmpz_poly_set(fmpz_poly_t output, const fmpz_poly_t input);

/* 
   Zero the polynomial by setting the length to zero.
   Does not set the actual limbs to zero.
*/

static inline 
void _fmpz_poly_zero(fmpz_poly_t output)
{
   output->length = 0;
}

/* Zero first n coefficients of poly, regardless of what length is */

void _fmpz_poly_zero_coeffs(fmpz_poly_t poly, const unsigned long n);

static inline 
void _fmpz_poly_attach(fmpz_poly_t output, const fmpz_poly_t input)
{
   output->length = input->length;
   output->limbs = input->limbs;
   output->coeffs = input->coeffs;
}

/*
   Attach input shifted right by n to output
*/

static inline 
void _fmpz_poly_attach_shift(fmpz_poly_t output, 
             const fmpz_poly_t input, unsigned long n)
{
   if (input->length >= n) output->length = input->length - n;
   else output->length = 0;
   output->limbs = input->limbs;
   output->coeffs = input->coeffs + n*(input->limbs+1);
}

/*
   Attach input to first n coefficients of input
*/

static inline 
void _fmpz_poly_attach_truncate(fmpz_poly_t output, 
             const fmpz_poly_t input, unsigned long n)
{
   if (input->length < n) output->length = input->length;
   else output->length = n;
   output->limbs = input->limbs;
   output->coeffs = input->coeffs;
   _fmpz_poly_normalise(output);
}

long _fmpz_poly_max_bits1(const fmpz_poly_t poly);

long _fmpz_poly_max_bits(const fmpz_poly_t poly);

unsigned long _fmpz_poly_max_limbs(const fmpz_poly_t poly);

int _fmpz_poly_equal(const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_neg(fmpz_poly_t output, const fmpz_poly_t input);

void _fmpz_poly_scalar_abs(fmpz_poly_t output, const fmpz_poly_t input);

void _fmpz_poly_truncate(fmpz_poly_t poly, const unsigned long trunc);

void _fmpz_poly_reverse(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long length);

void _fmpz_poly_left_shift(fmpz_poly_t output, const fmpz_poly_t input, 
                                                 const unsigned long n);

void _fmpz_poly_right_shift(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long n);

void _fmpz_poly_add(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_sub(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void _fmpz_poly_scalar_mul_fmpz(fmpz_poly_t output, const fmpz_poly_t poly, const fmpz_t x);

void _fmpz_poly_scalar_mul_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x);

void _fmpz_poly_scalar_mul_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x);

void _fmpz_poly_scalar_div_fmpz(fmpz_poly_t output, const fmpz_poly_t poly, const fmpz_t x);

void _fmpz_poly_scalar_tdiv_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x);

void _fmpz_poly_scalar_div_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x);

void _fmpz_poly_scalar_div_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x);

void _fmpz_poly_scalar_tdiv_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x);

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

void _fmpz_poly_mul_modular(fmpz_poly_t output, const fmpz_poly_t input1,
                                                const fmpz_poly_t input2, const ulong bits);

void _fmpz_poly_mul_KS(fmpz_poly_t output, const fmpz_poly_t input1, 
                                       const fmpz_poly_t input2, const long bits);
                                       
void _fmpz_poly_mul_KS_trunc(fmpz_poly_t output, const fmpz_poly_t in1, 
                                        const fmpz_poly_t in2, const unsigned long trunc, long bits);

void _fmpz_poly_mul_SS(fmpz_poly_t output, const fmpz_poly_t input1, 
                                                      const fmpz_poly_t input2);
                                                            
void _fmpz_poly_mul_SS_trunc(fmpz_poly_t output, const fmpz_poly_t input1, 
                           const fmpz_poly_t input2, const unsigned long trunc);
                                        
void _fmpz_poly_mul_trunc_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                const fmpz_poly_t input2, const unsigned long trunc);
                                
void _fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                const fmpz_poly_t input2, const unsigned long trunc);
                                
void _fmpz_poly_mul(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

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

void fmpz_poly_check_normalisation(const fmpz_poly_t poly);

// ------------------------------------------------------
// String conversions and I/O

int fmpz_poly_from_string(fmpz_poly_t poly, const char* s);

char* fmpz_poly_to_string(const fmpz_poly_t poly);

void fmpz_poly_print(const fmpz_poly_t poly);

void fmpz_poly_fprint(const fmpz_poly_t poly, FILE* f);

int fmpz_poly_fread(fmpz_poly_t poly, FILE* f);

char* fmpz_poly_to_string_pretty(const fmpz_poly_t poly, const char * x);

void fmpz_poly_fprint_pretty(const fmpz_poly_t poly, FILE* f, const char * x);

void fmpz_poly_print_pretty(const fmpz_poly_t poly, const char * x);

static inline 
int fmpz_poly_read(fmpz_poly_t poly)
{
   return fmpz_poly_fread(poly, stdin);
}

static inline 
unsigned long fmpz_poly_limbs(const fmpz_poly_t poly)
{
   return poly->limbs;
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
long fmpz_poly_max_bits1(const fmpz_poly_t poly)
{
   return _fmpz_poly_max_bits1(poly);
}

static inline 
long fmpz_poly_max_bits(const fmpz_poly_t poly)
{
   return _fmpz_poly_max_bits(poly);
}

static inline
long fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
   return _fmpz_poly_max_limbs(poly);
}

static inline 
void fmpz_poly_truncate(fmpz_poly_t poly, const unsigned long length)
{
   poly->length = FLINT_MIN(length, poly->length);
   _fmpz_poly_normalise(poly);
}

static inline
void fmpz_poly_swap(fmpz_poly_t x, fmpz_poly_t y)
{
   if (x == y) return;
   fmpz_t temp_p;
   
   mp_limb_t temp_l;
   
   temp_l = x->alloc;
   x->alloc = y->alloc;
   y->alloc = temp_l;
   
   temp_p = x->coeffs;
   x->coeffs = y->coeffs;
   y->coeffs = temp_p;
   
   temp_l = x->length;
   x->length = y->length;
   y->length = temp_l;
   
   temp_l = x->limbs;
   x->limbs = y->limbs;
   y->limbs = temp_l;
}


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
long fmpz_poly_get_coeff(mp_limb_t * output, const fmpz_poly_t poly,
                                                  const unsigned long n)
{
   if (n >= poly->length)
   {
       F_mpn_clear(output, poly->limbs);
       return 0;  
   }
   return _fmpz_poly_get_coeff(output, poly, n);
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
   if ((long) poly->coeffs[n*(poly->limbs+1)] >= 1L) 
                                 return poly->coeffs[n*(poly->limbs+1)+1];
   else return -poly->coeffs[n*(poly->limbs+1)+1];
}

void fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, const unsigned long n);

void fmpz_poly_get_coeff_mpz_read_only(mpz_t x, const fmpz_poly_t poly, const unsigned long n);

static inline
void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, const unsigned long n, fmpz_t x) 
{
   if ((poly->length == 0) && (x[0] == 0)) return;
	fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, fmpz_size(x));
   if (n+1 > poly->length) 
   {
      long i;
		for (i = poly->length; i < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0L;
      } 
      poly->length = n+1;
   }
   _fmpz_poly_set_coeff_fmpz(poly, n, x);
   _fmpz_poly_normalise(poly);
}

static inline
void fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, const unsigned long n) 
{
   if (n >= poly->length)
   {
      x[0] = 0;
      return;
   }
   _fmpz_poly_get_coeff_fmpz(x, poly, n);
} 

static inline
void fmpz_poly_set_coeff(fmpz_poly_t poly, const unsigned long n, 
                        const mp_limb_t * x, const long sign, const unsigned long size)
{
   if ((poly->length == 0) && (size == 0)) return;
	fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, size);
   if (n+1 > poly->length) 
   {
      long i;
		for (i = poly->length; i + 1 < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0L;
      } 
      poly->length = n+1;
   }
   _fmpz_poly_set_coeff(poly, n, x, sign, size);
   _fmpz_poly_normalise(poly);
}    

static inline 
void fmpz_poly_set_coeff_si(fmpz_poly_t poly, const unsigned long n, const long x)
{
   if ((poly->length == 0) && (x == 0)) return;
	fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, 1);
   if (n+1 > poly->length) 
   {
      long i;
		for (i = poly->length; i < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0L;
      } 
      poly->length = n+1;
   }
   _fmpz_poly_set_coeff_si(poly, n, x);
   _fmpz_poly_normalise(poly);
}

static inline
void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, const unsigned long n, const unsigned long x)
{
   if ((poly->length == 0) && (x == 0)) return;
	fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, 1);
   if (n+1 > poly->length) 
   {
      long i;
		for (i = poly->length; i < n; i++)
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
   if ((poly->length == 0) && (mpz_size(x) == 0)) return;
	fmpz_poly_fit_length(poly, n+1);
   fmpz_poly_fit_limbs(poly, mpz_size(x));
   if (n+1 > poly->length) 
   {
      long i;
		for (i = poly->length; i < n; i++)
      {
         poly->coeffs[i*(poly->limbs+1)] = 0L;
      } 
      poly->length = n+1;
   }
   
   _fmpz_poly_set_coeff_mpz(poly, n, x);
}

static inline
void fmpz_poly_set(fmpz_poly_t output, const fmpz_poly_t input)
{
   fmpz_poly_fit_length(output, input->length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_set(output, input);
}

static inline
int fmpz_poly_equal(const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   return _fmpz_poly_equal(input1, input2);
}

static inline
void fmpz_poly_zero(fmpz_poly_t output)
{
   output->length = 0;
}

static inline 
void fmpz_poly_zero_coeffs(fmpz_poly_t poly, const unsigned long n)
{
   fmpz_poly_fit_length(poly, n);
   if (n >= poly->length) 
   {
      fmpz_poly_zero(poly);
      return;
   }
   _fmpz_poly_zero_coeffs(poly, n);
}

static inline
void fmpz_poly_neg(fmpz_poly_t output, const fmpz_poly_t input)
{
   fmpz_poly_fit_length(output, input->length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_neg(output, input);
}

static inline
void fmpz_poly_scalar_abs(fmpz_poly_t output, const fmpz_poly_t input)
{
   fmpz_poly_fit_length(output, input->length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_scalar_abs(output, input);
}

static inline
void fmpz_poly_reverse(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long length)
{
   fmpz_poly_fit_length(output, length);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_reverse(output, input, length);
}

static inline
void fmpz_poly_left_shift(fmpz_poly_t output, const fmpz_poly_t input, 
                                                 const unsigned long n)
{
   if (input->length + n == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   fmpz_poly_fit_length(output, input->length + n);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_left_shift(output, input, n);  
}

static inline
void fmpz_poly_right_shift(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long n)
{
   if ((long)(input->length - n) <= 0L)
   {
      fmpz_poly_zero(output);
      return;
   }
   fmpz_poly_fit_length(output, input->length - n);
   fmpz_poly_fit_limbs(output, input->limbs);
   _fmpz_poly_right_shift(output, input, n);
}

void fmpz_poly_2norm(fmpz_t norm, const fmpz_poly_t pol);

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const fmpz_t x);
                          
void fmpz_poly_scalar_mul_ui(fmpz_poly_t output, 
                          const fmpz_poly_t input, unsigned long x);
                          
void fmpz_poly_scalar_mul_si(fmpz_poly_t output, 
                          const fmpz_poly_t input, long x);
                          
void fmpz_poly_scalar_mul_mpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const mpz_t x);
                          
static inline
void fmpz_poly_scalar_div_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs = fmpz_poly_max_limbs(poly);
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, limbs);
   _fmpz_poly_scalar_div_ui(output, poly, x);
}
                          
static inline
void fmpz_poly_scalar_div_si(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs = fmpz_poly_max_limbs(poly);
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, limbs);
   _fmpz_poly_scalar_div_si(output, poly, x);
}
                          
static inline
void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs = fmpz_poly_max_limbs(poly);
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, limbs);
   _fmpz_poly_scalar_tdiv_ui(output, poly, x);
}
                          
static inline
void fmpz_poly_scalar_tdiv_si(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs = fmpz_poly_max_limbs(poly);
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, limbs);
   _fmpz_poly_scalar_tdiv_si(output, poly, x);
}

static inline
void fmpz_poly_scalar_div_exact_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs = fmpz_poly_max_limbs(poly);
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, limbs);
   _fmpz_poly_scalar_div_exact_ui(output, poly, x);
}
                          
static inline
void fmpz_poly_scalar_div_exact_si(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs = fmpz_poly_max_limbs(poly);
   fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, limbs);
   _fmpz_poly_scalar_div_exact_si(output, poly, x);
}
                          
void fmpz_poly_scalar_div_fmpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const fmpz_t x);
                          
void fmpz_poly_scalar_div_mpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const mpz_t x);

static inline
void fmpz_poly_add(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   if (input1 == input2) 
   {
      fmpz_poly_scalar_mul_ui(output, input1, 2UL);
      
      return;
   }

   unsigned long bits1 = FLINT_ABS(_fmpz_poly_max_bits(input1));
   unsigned long bits2 = FLINT_ABS(_fmpz_poly_max_bits(input2));
   
   fmpz_poly_fit_length(output, FLINT_MAX(input1->length, input2->length));
   fmpz_poly_fit_limbs(output, FLINT_MAX(bits1, bits2)/FLINT_BITS + 1);
   
   _fmpz_poly_add(output, input1, input2);
}

static inline
void fmpz_poly_sub(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   if (input1 == input2) 
   {
      fmpz_poly_zero(output);
      
      return;
   }
   
   unsigned long bits1 = FLINT_ABS(_fmpz_poly_max_bits(input1));
   unsigned long bits2 = FLINT_ABS(_fmpz_poly_max_bits(input2));
   
   fmpz_poly_fit_length(output, FLINT_MAX(input1->length, input2->length));
   fmpz_poly_fit_limbs(output, FLINT_MAX(bits1, bits2)/FLINT_BITS + 1);
   
   _fmpz_poly_sub(output, input1, input2);
}

void fmpz_poly_mul(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2);

void fmpz_poly_mul_modular(fmpz_poly_t output, const fmpz_poly_t input1,
                                                const fmpz_poly_t input2, const ulong bits);

void fmpz_poly_mul_trunc_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc);
                                          
void fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc);

void fmpz_poly_div_classical(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_divrem_classical(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_div_divconquer_recursive(fmpz_poly_t Q, fmpz_poly_t DQ, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_div_divconquer(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_div_mulders(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_newton_invert_basecase(fmpz_poly_t Q_inv, const fmpz_poly_t Q, const unsigned long n);

void fmpz_poly_newton_invert(fmpz_poly_t Q_inv, const fmpz_poly_t Q, const unsigned long n);

static inline
void fmpz_poly_invert_series(fmpz_poly_t Q_inv, const fmpz_poly_t Q, const unsigned long n)
{
	fmpz_poly_newton_invert(Q_inv, Q, n);
}

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, const unsigned long n);

void fmpz_poly_div_newton(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B);

static inline
void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A == B)
   {
      fmpz_poly_fit_length(Q, 1);
      fmpz_poly_fit_limbs(Q, 1);
      fmpz_poly_zero(Q);
      fmpz_poly_set_coeff_ui(Q, 0, 1UL);
      
      return;
   }
   
   fmpz_poly_t Ain, Bin;
   
   if (A == Q)
   {
      _fmpz_poly_stack_init(Ain, A->length, A->limbs);
      _fmpz_poly_set(Ain, A);
   } else _fmpz_poly_attach(Ain, A);
   
   if (B == Q)
   {
      _fmpz_poly_stack_init(Bin, B->length, B->limbs);
      _fmpz_poly_set(Bin, B);
   } else _fmpz_poly_attach(Bin, B);

   fmpz_poly_div_mulders(Q, Ain, Bin);

   if (A == Q) _fmpz_poly_stack_clear(Ain);
   if (B == Q) _fmpz_poly_stack_clear(Bin);
}

static inline
void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A == B)
   {
      fmpz_poly_fit_length(Q, 1);
      fmpz_poly_fit_limbs(Q, 1);
      fmpz_poly_zero(Q);
      fmpz_poly_zero(R);
      fmpz_poly_set_coeff_ui(Q, 0, 1UL);
      
      return;
   }
   
   fmpz_poly_t Ain, Bin;
   
   if ((A == R) || (A == Q))
   {
      _fmpz_poly_stack_init(Ain, A->length, A->limbs);
      _fmpz_poly_set(Ain, A);
   } else _fmpz_poly_attach(Ain, A);
   
   if ((B == R) || (B == Q))
   {
      _fmpz_poly_stack_init(Bin, B->length, B->limbs);
      _fmpz_poly_set(Bin, B);
   } else _fmpz_poly_attach(Bin, B);

   fmpz_poly_divrem_divconquer(Q, R, Ain, Bin);

   if ((A == R) || (A == Q)) _fmpz_poly_stack_clear(Ain);
   if ((B == R) || (B == Q)) _fmpz_poly_stack_clear(Bin);
}

/*
   Returns 1 and the quotient if B divides A
   else returns 0
*/

static inline
int fmpz_poly_divides_divconquer(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t R;
   int divides = 0;

   fmpz_poly_init(R);

   fmpz_poly_divrem(Q, R, A, B);
   if (R->length == 0) divides = 1;
   fmpz_poly_clear(R);
   
   return divides;
}

int fmpz_poly_divides_modular(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, const ulong bound);

static inline
int fmpz_poly_divides(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A == B)
	{
		fmpz_poly_set_coeff_ui(Q, 0, 1L);
		Q->length = 1;
		return 1;
	}
	
	if ((Q == A) || (Q == B))
	{
		fmpz_poly_t Qa;
		fmpz_poly_init(Qa);
      int divides = fmpz_poly_divides_modular(Qa, A, B, 0);
		fmpz_poly_swap(Q, Qa);
		fmpz_poly_clear(Qa);
		return divides;
	} else
	   return fmpz_poly_divides_modular(Q, A, B, 0);
}


void fmpz_poly_power(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long exp);

void fmpz_poly_power_trunc_n(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long exp, const unsigned long n);

void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

void fmpz_poly_pseudo_divrem_shoup(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B);

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
   if (A == B)
   {
      fmpz_poly_fit_length(Q, 1);
      fmpz_poly_fit_limbs(Q, 1);
      fmpz_poly_zero(Q);
      fmpz_poly_zero(R);
      d = 0;
      fmpz_poly_set_coeff_ui(Q, 0, 1UL);
      
      return;
   }

   fmpz_poly_t Ain, Bin;
   
   if ((A == R) || (A == Q))
   {
      _fmpz_poly_stack_init(Ain, A->length, A->limbs);
      _fmpz_poly_set(Ain, A);
   } else _fmpz_poly_attach(Ain, A);
   
   if ((B == R) || (B == Q))
   {
      _fmpz_poly_stack_init(Bin, B->length, B->limbs);
      _fmpz_poly_set(Bin, B);
   } else _fmpz_poly_attach(Bin, B);

   fmpz_poly_pseudo_divrem_recursive(Q, R, d, Ain, Bin);
   
   if ((A == R) || (A == Q)) _fmpz_poly_stack_clear(Ain);
   if ((B == R) || (B == Q)) _fmpz_poly_stack_clear(Bin);
}

void fmpz_poly_pseudo_div_recursive(fmpz_poly_t Q, unsigned long * d, 
                                      const fmpz_poly_t A, const fmpz_poly_t B);
                                      
static inline
void fmpz_poly_pseudo_div(fmpz_poly_t Q, unsigned long * d, 
                                      const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A == B)
   {
      fmpz_poly_fit_length(Q, 1);
      fmpz_poly_fit_limbs(Q, 1);
      fmpz_poly_zero(Q);
      d = 0;
      fmpz_poly_set_coeff_ui(Q, 0, 1UL);
      
      return;
   }
   
   fmpz_poly_t Ain, Bin;
   
   if (A == Q)
   {
      _fmpz_poly_stack_init(Ain, A->length, A->limbs);
      _fmpz_poly_set(Ain, A);
   } else _fmpz_poly_attach(Ain, A);
   
   if (B == Q)
   {
      _fmpz_poly_stack_init(Bin, B->length, B->limbs);
      _fmpz_poly_set(Bin, B);
   } else _fmpz_poly_attach(Bin, B);

   fmpz_poly_pseudo_div_recursive(Q, d, Ain, Bin);
  
   if (A == Q) _fmpz_poly_stack_clear(Ain);
   if (B == Q) _fmpz_poly_stack_clear(Bin);
}
                                      
void fmpz_poly_content(fmpz_t c, const fmpz_poly_t poly);

static inline
void _fmpz_poly_primitive_part(fmpz_poly_t prim, fmpz_poly_t poly)
{
   if (poly->length == 0)
   {
      _fmpz_poly_zero(prim);
      return;
   }

   fmpz_t c = fmpz_init(poly->limbs);
   fmpz_poly_content(c, poly);
   _fmpz_poly_scalar_div_fmpz(prim, poly, c);
   fmpz_clear(c);
}

static inline
void fmpz_poly_primitive_part(fmpz_poly_t prim, fmpz_poly_t poly)
{
   if (poly->length == 0)
   {
      fmpz_poly_zero(prim);
      return;
   }

   fmpz_poly_fit_length(prim, poly->length);
	fmpz_t c = fmpz_init(poly->limbs);
   fmpz_poly_content(c, poly);
   fmpz_poly_fit_limbs(prim, poly->limbs - FLINT_ABS(c[0]) + 1);
	fmpz_poly_scalar_div_fmpz(prim, poly, c);
   fmpz_clear(c);
}

void fmpz_poly_gcd_subresultant(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_gcd_modular(fmpz_poly_t H, const fmpz_poly_t poly1, 
					    const fmpz_poly_t poly2, const ulong bits1, const ulong bits2);

int fmpz_poly_gcd_heuristic(fmpz_poly_t H, const fmpz_poly_t poly1, 
						 const fmpz_poly_t poly2, const ulong bits1, const ulong bits2);

void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2);

void fmpz_poly_invmod_modular(fmpz_t d, fmpz_poly_t H, fmpz_poly_t poly1, fmpz_poly_t poly2);

static inline
void fmpz_poly_invmod(fmpz_t d, fmpz_poly_t H, fmpz_poly_t poly1, fmpz_poly_t poly2)
{
   if ((H == poly1) || (H == poly2))
	{
		fmpz_poly_t res;
		fmpz_poly_init(res);
      fmpz_poly_invmod_modular(d, res, poly1, poly2);
		fmpz_poly_swap(H, res);
		fmpz_poly_clear(res);
	} else
		fmpz_poly_invmod_modular(d, H, poly1, poly2);
}

unsigned long fmpz_poly_resultant_bound(fmpz_poly_t a, fmpz_poly_t b);

void fmpz_poly_resultant(fmpz_t res, fmpz_poly_t a, fmpz_poly_t b);

void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, fmpz_poly_t a, fmpz_poly_t b);

static inline
void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, fmpz_poly_t a, fmpz_poly_t b)
{
   fmpz_poly_t T, S;
	
	if ((t == a) || (t == b))
	{
		fmpz_poly_init(T);
		if ((s == a) || (s == b))
		{
			fmpz_poly_init(S);
         fmpz_poly_xgcd_modular(r, S, T, a, b);
		   fmpz_poly_swap(S, s);
		   fmpz_poly_clear(S);
		} else fmpz_poly_xgcd_modular(r, s, T, a, b);
      fmpz_poly_swap(T, t);
		fmpz_poly_clear(T);
	} else if ((s == a) || (s == b))
	{
		fmpz_poly_init(S);
      fmpz_poly_xgcd_modular(r, S, t, a, b);
		fmpz_poly_swap(S, s);
		fmpz_poly_clear(S);
	} else
	   fmpz_poly_xgcd_modular(r, s, t, a, b);
}

void fmpz_poly_derivative(fmpz_poly_t der, fmpz_poly_t poly);

void fmpz_poly_evaluate_horner_range(fmpz_t output, fmpz_poly_t poly, fmpz_t val, ulong start, ulong n);

static inline
void fmpz_poly_evaluate_horner(fmpz_t output, fmpz_poly_t poly, fmpz_t val)
{
	fmpz_poly_evaluate_horner_range(output, poly, val, 0, poly->length);
}

void fmpz_poly_evaluate_divconquer(fmpz_t output, fmpz_poly_t poly, fmpz_t val);

void fmpz_poly_evaluate(fmpz_t output, fmpz_poly_t poly, fmpz_t val);

void fmpz_poly_compose(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val);

void fmpz_poly_array_compose_divconquer(fmpz_poly_t output, fmpz_poly_t * poly, ulong length, fmpz_poly_t val);

void fmpz_poly_compose_divconquer(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val);

void fmpz_poly_compose_horner_range(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val, ulong start, ulong n);

static inline
void fmpz_poly_compose_horner(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val)
{
	if (poly->length == 0) 
	{
		fmpz_poly_zero(output);
		return;
	}
	
	if (val->length == 0)
	{
		fmpz_poly_zero(output);
		fmpz_poly_set_coeff_fmpz(output, 0, poly->coeffs); 
		return;
	}

	if (poly->length == 1)
	{
		fmpz_poly_set(output, poly);
		return;
	}

	fmpz_poly_compose_horner_range(output, poly, val, 0, poly->length);
}

// *************** end of file

#ifdef __cplusplus
 }
#endif
 
#endif
