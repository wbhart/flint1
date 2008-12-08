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

fmpz_poly.c: Polynomials over Z, implemented as contiguous block of fmpz_t's

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mpz_poly.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "long_extras.h"
#include "zmod_poly.h"

#ifdef HAVE_ZNPOLY
#include "zn_poly.h"
#endif

/****************************************************************************

   Conversion Routines
   
****************************************************************************/

/* 
   Convert length coefficients of an fmpz_poly_t to an
   already initialised ZmodF_poly_t. Each coefficient will
   be represented mod p = 2^Bn+1 where n is given by the field
   n of the ZmodF_poly_t. Coefficients will be assumed to 
   be in the range [-p/2, p/2].
   Assumes 0 < length <= poly_fmpz->length 
*/
   
long fmpz_poly_to_ZmodF_poly(ZmodF_poly_t poly_f, const fmpz_poly_t poly_fmpz, 
                                                        const unsigned long length)
{
   unsigned long size_f = poly_f->n + 1;
   unsigned long size_m = poly_fmpz->limbs+1;
   mp_limb_t * coeffs_m = poly_fmpz->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;
   
   unsigned long mask = -1L;
   long bits = 0;
   long limbs = 0;
   long sign = 1;
   
   long size_j;
   
   for (unsigned long i = 0, j = 0; i < length; i++, j += size_m)
   {
      size_j = coeffs_m[j];
      if ((long) size_j < 0) sign = -1L;
      if (ABS(size_j) > limbs + 1)
      {
         limbs = ABS(size_j) - 1;
         bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]); 
         if (bits == FLINT_BITS) mask = 0L;
         else mask = -1L - ((1L<<bits)-1);  
      } else if (ABS(size_j) == limbs + 1)
      {
         if (coeffs_m[j+ABS(size_j)] & mask)
         {
            bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]);   
            if (bits == FLINT_BITS) mask = 0L;
            else mask = -1L - ((1L<<bits)-1);
         }
      }
      if (size_j < 0)
      {
         F_mpn_negate(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         F_mpn_set(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      } else
      {
         F_mpn_copy(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         F_mpn_clear(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      }
   }
   poly_f->length = length; 
   
   return sign*(FLINT_BITS*limbs+bits);  
}

/* 
   Convert a ZmodF_poly_t to an fmpz_poly_t. Coefficients will
   be taken to be in the range [-p/2, p/2] where p = 2^nB+1.
   Assumes 0 < poly_f->length 
*/

void ZmodF_poly_to_fmpz_poly(fmpz_poly_t poly_fmpz, const ZmodF_poly_t poly_f, const long sign)
{
   unsigned long n = poly_f->n;
   unsigned long size_m = poly_fmpz->limbs+1;
   unsigned long limbs = FLINT_MIN(n, size_m-1);
   
   mp_limb_t * coeffs_m = poly_fmpz->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;

   if (sign)
   {
      for (unsigned long i = 0, j = 0; i < poly_f->length; i++, j += size_m)
      {
         ZmodF_normalise(coeffs_f[i], n);
         if (coeffs_f[i][n-1]>>(FLINT_BITS-1) || coeffs_f[i][n])
         {
            F_mpn_negate(coeffs_m + j + 1, coeffs_f[i], limbs);
            mpn_add_1(coeffs_m + j + 1, coeffs_m + j + 1, limbs, 1L);
            coeffs_m[j] = -limbs;
            NORM(coeffs_m + j);
         } else
         {
            F_mpn_copy(coeffs_m + j + 1, coeffs_f[i], limbs);
            coeffs_m[j] = limbs;
            NORM(coeffs_m + j);         
         }
      }
   } else 
   {
      for (unsigned long i = 0, j = 0; i < poly_f->length; i++, j += size_m)
      {
         ZmodF_normalise(coeffs_f[i], n);
         F_mpn_copy(coeffs_m + j + 1, coeffs_f[i], limbs);
         coeffs_m[j] = limbs;
         NORM(coeffs_m + j);         
      }
   }
   
   poly_fmpz->length = poly_f->length;
   _fmpz_poly_normalise(poly_fmpz);   
}

static inline long __get_next_coeff(const mp_limb_t * coeff_m, long * borrow, long * coeff, const long mask, const long negate)
{ 
   if ((long) coeff_m[0] == 0) *coeff = -*borrow;
   else if ((((long) coeff_m[0]) ^ negate) >= 0L) *coeff = coeff_m[1] - *borrow;
   else *coeff = (-coeff_m[1] - *borrow);
   *borrow = 0UL;
   if (*coeff < 0) 
   {
      *borrow = 1UL;
   }
   *coeff&=mask;
   
   return *coeff;
}

static inline long __get_next_coeff_unsigned(const mp_limb_t * coeff_m, long * coeff)
{ 
   if ((long) coeff_m[0] == 0) *coeff = 0;
   else *coeff = coeff_m[1];
   
   return *coeff;
}

void fmpz_poly_bit_pack(fmpz_t array, const fmpz_poly_t poly_fmpz,
                            const unsigned long length, const long bitwidth, 
                            const long negate)
{   
   ulong i, k, skip;
   mp_limb_t * coeff_m = poly_fmpz->coeffs;
   mp_limb_t * last_point;
    
   ulong temp;
   half_ulong lower;
   long coeff;
   long borrow;
   mp_limb_t extend;
   
   long bits = bitwidth;
   int sign = (bits < 0);
   if (sign) bits = ABS(bits);
   
   ulong coeffs_per_limb = FLINT_BITS/bits;

   const ulong mask = (1UL<<bits)-1;
      
	k = 0; skip = 0;
   coeff = 0; borrow = 0L; temp = 0;
	ulong size = poly_fmpz->limbs + 1;
      
   last_point = coeff_m + size*length;
         
   while (coeff_m < last_point)
   {
      if ((ulong) coeff_m & 7 == 0) FLINT_PREFETCH(coeff_m, 64);
      // k is guaranteed to be less than FLINT_BITS at this point
      while ((k < HALF_FLINT_BITS) && (coeff_m < last_point))
      {
         if (sign) temp+=(__get_next_coeff(coeff_m, &borrow, &coeff, mask, negate) << k);
         else temp+=(__get_next_coeff_unsigned(coeff_m, &coeff) << k);
         coeff_m+=size; k+=bits;
      }
      // k may exceed FLINT_BITS at this point but is less than 96

      if (k > FLINT_BITS)
      {
         // if k > FLINT_BITS write out a whole limb and read in remaining bits of coeff
			array[skip] = temp;
         skip++;
         temp = (coeff>>(bits+FLINT_BITS-k));
         k = (k-FLINT_BITS);
         // k < HALF_FLINT_BITS
      } else
      {
         // k <= FLINT_BITS
         if (k >= HALF_FLINT_BITS)
         {
            // if k >= HALF_FLINT_BITS store bottom HALF_FLINT_BITS bits
            lower = (half_ulong)temp;
            k-=HALF_FLINT_BITS;
            temp>>=HALF_FLINT_BITS;
            // k is now <= HALF_FLINT_BITS

            while ((k<HALF_FLINT_BITS) && (coeff_m < last_point))
            {
               if (sign) temp+=(__get_next_coeff(coeff_m, &borrow, &coeff, mask, negate) << k);
               else temp+=(__get_next_coeff_unsigned(coeff_m, &coeff) << k);
               coeff_m+=size; k+=bits;
            }
            // k may again exceed FLINT_BITS bits but is less than 96
            if (k>FLINT_BITS)
            {
               // if k > FLINT_BITS, write out bottom HALF_FLINT_BITS bits (along with HALF_FLINT_BITS bits from lower)
               // read remaining bits from coeff and reduce k by HALF_FLINT_BITS
               array[skip] = (temp<<HALF_FLINT_BITS) + (ulong) lower;
               skip++;
               temp>>=HALF_FLINT_BITS;
               temp+=((coeff>>(bits+FLINT_BITS-k))<<HALF_FLINT_BITS);
               k = (k-HALF_FLINT_BITS);
               // k < FLINT_BITS and we are ready to read next coefficient if there is one
            } else if (k >= HALF_FLINT_BITS) 
            {
               // k <= FLINT_BITS
               // if k >= HALF_FLINT_BITS write out bottom HALF_FLINT_BITS bits (along with lower)
               // and reduce k by HALF_FLINT_BITS
               k-=HALF_FLINT_BITS;
               array[skip] = (temp<<HALF_FLINT_BITS)+lower;
               temp>>=HALF_FLINT_BITS;
               skip++;
               // k is now less than or equal to HALF_FLINT_BITS and we are now ready to read 
               // the next coefficient if there is one
            } else
            {
               // k < HALF_FLINT_BITS
               // there isn't enough to write out a whole FLINT_BITS bits, so put it all 
               // together in temp
               temp = (temp<<HALF_FLINT_BITS)+lower;
               k += HALF_FLINT_BITS;
               // k is now guaranteed to be less than FLINT_BITS and we are ready for the
               // next coefficient if there is one
            } // else
         } // if
      } // else
	} // while
       
   // write out final coefficient/limb
   array[skip] = temp;
}

void fmpz_poly_bit_unpack(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                              const ulong length, const ulong bits)
{
   ulong k, skip;

   ulong temp2;
   ulong temp;
   ulong full_limb;
   ulong carry;
     
   const ulong mask = (1UL<<bits)-1;
   const ulong sign_mask = (1UL<<(bits-1));

   ulong s;

	fmpz_poly_fit_length(poly_fmpz, length + 1);
	fmpz_poly_fit_limbs(poly_fmpz, 1);
   mp_limb_t * coeff_m = poly_fmpz->coeffs;
   mp_limb_t * last_point;
   ulong size_m = poly_fmpz->limbs+1;
    
   k=0; skip=0; carry = 0UL; temp2 = 0;
   last_point = coeff_m + size_m*length;
      
   while (coeff_m < last_point)
   {
      // read in a full limb
      full_limb = array[skip];
      temp2 += l_shift(full_limb,k);
      s = FLINT_BITS-k;
      k+=s;
      while ((k >= bits)&&(coeff_m < last_point))
      {
         if (!(temp2 & sign_mask)) 
         {
            fmpz_add_ui_inplace(coeff_m, (temp2 & mask) + carry);
            carry = 0UL;
         }  
         else
         {
            temp = ((-temp2) & mask) - carry;
            fmpz_sub_ui_inplace(coeff_m, temp);
            carry = 1UL;
         }
         coeff_m += size_m;
         temp2>>=bits;
         k-=bits;
      }
      // k is now less than bits
      // read in remainder of full_limb
      temp2 += l_shift(r_shift(full_limb, s), k);
      k+=(FLINT_BITS - s);
       
      while ((k >= bits)&&(coeff_m < last_point))
      {
         if (!(temp2 & sign_mask)) 
         {
            fmpz_add_ui_inplace(coeff_m, (temp2 & mask) + carry);
            carry = 0UL;
         }
         else
         {
            temp = ((-temp2) & mask) - carry;
            fmpz_sub_ui_inplace(coeff_m, temp);
            carry = 1UL;
         }
         coeff_m += size_m;
         temp2>>=bits;
         k-=bits;
      }
      // k is now less than bits
      skip++;
   }

	poly_fmpz->length = length;
   _fmpz_poly_normalise(poly_fmpz);
}
 
void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                              const unsigned long length, const unsigned long bits)
{
   unsigned long k, l, skip;

   unsigned long temp2;
   unsigned long temp;
   unsigned long full_limb;
     
   const unsigned long mask = (1UL<<bits)-1;

   unsigned long s;
   fmpz_t coeff_m = poly_fmpz->coeffs;
   fmpz_t next_point;
   unsigned long size_m = poly_fmpz->limbs+1;
         
   k=0; skip=0; temp2 = 0;
   next_point = coeff_m + size_m*length;
      
   while (coeff_m < next_point)
   {
      if (skip & 7 == 0) FLINT_PREFETCH(array + skip, 64);
      // read in a full limb
      full_limb = array[skip];
      temp2 += l_shift(full_limb,k);
      s = FLINT_BITS - k;
      k += s;
      while ((k >= bits) && (coeff_m < next_point))
      {
         __fmpz_add_ui_inplace(coeff_m, (temp2 & mask));
         coeff_m += size_m;
         temp2>>=bits;
         k-=bits;
      }
      // k is now less than bits
      // read in remainder of full_limb
      temp2 += l_shift(r_shift(full_limb,s),k);
      k += (FLINT_BITS-s);
       
      while ((k >= bits)&&(coeff_m < next_point))
      {
         __fmpz_add_ui_inplace(coeff_m, temp2 & mask);
         coeff_m += size_m;
         temp2 >>= bits;
         l++;
         k -= bits;
      }
      // k is now less than bits
      skip++;
   }

   poly_fmpz->length = length;
	_fmpz_poly_normalise(poly_fmpz);
}

void fmpz_poly_limb_pack(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
                                           const unsigned long length, const long limbs)
{
   unsigned long size_m = poly_fmpz->limbs + 1;
   long size_j;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
   long carry = 0;
   
   for (unsigned long i = 0, j = 0, k = 0; i < length; i++, j += size_m, k += limbs)
   {
      size_j = (long) coeffs_m[j];
      if (size_j < 0)
      {
         F_mpn_negate(array + k, coeffs_m + j + 1, ABS(size_j)); 
         F_mpn_set(array + k + ABS(size_j), limbs - ABS(size_j));
         if (carry) mpn_sub_1(array + k, array + k, limbs, 1L);
         carry = 1L;
      } else if (size_j > 0)
      {
         F_mpn_copy(array + k, coeffs_m + j + 1, ABS(size_j)); 
         F_mpn_clear(array + k + ABS(size_j), limbs - ABS(size_j)); 
         if (carry) mpn_sub_1(array + k, array + k, limbs, 1L);
         carry = 0L;
      } else
      {
         if (carry) 
         {
            F_mpn_set(array + k, limbs);
            carry = 1L;
         } else 
         {
            F_mpn_clear(array + k, limbs);
            carry = 0L;
         }
      }
   }
}

void fmpz_poly_limb_pack_neg(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
                                           const unsigned long length, const long limbs)
{
   unsigned long size_m = poly_fmpz->limbs + 1;
   long size_j;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
   long carry = 0;
   
   for (unsigned long i = 0, j = 0, k = 0; i < length; i++, j += size_m, k += limbs)
   {
      size_j = (long) coeffs_m[j];
      if (size_j > 0)
      {
         F_mpn_negate(array + k, coeffs_m + j + 1, ABS(size_j)); 
         F_mpn_set(array + k + ABS(size_j), limbs - ABS(size_j));
         if (carry) mpn_sub_1(array + k, array + k, limbs, 1L);
         carry = 1L;
      } else if (size_j < 0)
      {
         F_mpn_copy(array + k, coeffs_m + j + 1, ABS(size_j)); 
         F_mpn_clear(array + k + ABS(size_j), limbs - ABS(size_j)); 
         if (carry) mpn_sub_1(array + k, array + k, limbs, 1L);
         carry = 0L;
      } else
      {
         if (carry) 
         {
            F_mpn_set(array + k, limbs);
            carry = 1L;
         } else 
         {
            F_mpn_clear(array + k, limbs);
            carry = 0L;
         }
      }
   }
}

void fmpz_poly_limb_unpack(fmpz_poly_t poly_fmpz, mp_limb_t * array, 
                                  const unsigned long length, const unsigned long limbs)
{
   unsigned long carry = 0L;
   fmpz_poly_fit_length(poly_fmpz, length + 1);
	fmpz_poly_fit_limbs(poly_fmpz, limbs);
	fmpz_t coeffs_m = poly_fmpz->coeffs;
   unsigned long size_m = poly_fmpz->limbs + 1;
   ulong i, j, k;
   
   for (i = 0, j = 0, k = 0; i < length; i++, j += size_m, k += limbs)
   {
      if (array[k+limbs-1]>>(FLINT_BITS-1))
      {
         F_mpn_negate(coeffs_m + j + 1, array + k, limbs);
         coeffs_m[j] = -limbs;
         NORM(coeffs_m + j);
			if (carry) mpn_sub_1(coeffs_m + j + 1, coeffs_m + j + 1, -coeffs_m[j], 1L);
			if (!coeffs_m[j - coeffs_m[j]]) coeffs_m[j]++;
			carry = 1L;
      } else
      {
         F_mpn_copy(coeffs_m + j + 1, array + k, limbs);
         coeffs_m[j] = limbs;
         NORM(coeffs_m + j); 
			mp_limb_t cry = 0L;
			if (carry) 
			{
				if (coeffs_m[j]) cry = mpn_add_1(coeffs_m + j + 1, coeffs_m + j + 1, coeffs_m[j], 1L);
				else cry = 1L;
			}
			if (cry) 
			{
				coeffs_m[j + coeffs_m[j] + 1] = cry;
				coeffs_m[j]++;
			}
			carry = 0L;
      }
   }

	if (carry) // check if there was a carried 1 from the final mpn_add_1
	{
		coeffs_m[j+1] = 1;
      coeffs_m[j] = 1;
		poly_fmpz->length = length + 1;
	} else
	{
		poly_fmpz->length = length;
	   _fmpz_poly_normalise(poly_fmpz);
	}
}

void fmpz_poly_limb_pack_1(mp_limb_t * array, const fmpz_poly_t poly_fmpz)
{
   unsigned long size_m = poly_fmpz->limbs + 1;
   long size_j;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
	unsigned long length = poly_fmpz->length;
   long carry = 0;
   
   for (unsigned long i = 0, j = 0, k = 0; i < length; i++, j += size_m, k++)
   {
      size_j = (long) coeffs_m[j];
      if (size_j < 0L)
      {
         if (carry) array[k] = -coeffs_m[j + 1] - 1; 
         else array[k] = -coeffs_m[j + 1];
			carry = 1L;
      } else if (size_j > 0L)
      {
         if (carry) array[k] = coeffs_m[j + 1] - 1; 
         else array[k] = coeffs_m[j + 1];
			carry = 0L;
      } else
      {
         if (carry) 
         {
            array[k] = -1L;
				carry = 1L;
         } else 
         {
            array[k] = 0L;
				carry = 0L;
         }
      }
   }
}

void fmpz_poly_limb_pack_neg_1(mp_limb_t * array, const fmpz_poly_t poly_fmpz)
{
   unsigned long size_m = poly_fmpz->limbs + 1;
   long size_j;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
	unsigned long length = poly_fmpz->length;
   long carry = 0;
   
   for (unsigned long i = 0, j = 0, k = 0; i < length; i++, j += size_m, k++)
   {
      size_j = (long) coeffs_m[j];
      if (size_j > 0L)
      {
         if (carry) array[k] = -coeffs_m[j + 1] - 1; 
         else array[k] = -coeffs_m[j + 1];
			carry = 1L;
      } else if (size_j < 0L)
      {
         if (carry) array[k] = coeffs_m[j + 1] - 1; 
         else array[k] = coeffs_m[j + 1];
			carry = 0L;
      } else
      {
         if (carry) 
         {
            array[k] = -1L;
				carry = 1L;
         } else 
         {
            array[k] = 0L;
				carry = 0L;
         }
      }
   }
}

void fmpz_poly_limb_unpack_1(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                                  const unsigned long length)
{
   fmpz_poly_fit_length(poly_fmpz, length + 1);
	fmpz_poly_fit_limbs(poly_fmpz, 1);
	
	unsigned long size_m = poly_fmpz->limbs + 1;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
   unsigned long carry = 0L;
   ulong i, j, k;

   for (i = 0, j = 0, k = 0; i < length; i++, j += size_m, k++)
   {
      if ((long) array[k] < 0L)
      {
         coeffs_m[j+1] = -array[k] - carry;
			if (coeffs_m[j+1] == 0L)
			   coeffs_m[j] = 0L;
		   else
			   coeffs_m[j] = -1;
			carry = 1L;
      } else if (array[k] == 0L)
		{
			if (carry)
			{
				coeffs_m[j+1] = 1L;
				coeffs_m[j] = 1L;
				carry = 0L;
			} else
			{
				coeffs_m[j] = 0L;
				carry = 0L;
			}
		} else
      {
         coeffs_m[j+1] = array[k] + carry;
			coeffs_m[j] = 1;
         NORM(coeffs_m + j); 
         carry = 0L;        
      }
   }

	if (carry) 
	{
		coeffs_m[j+1] = 1;
      coeffs_m[j] = 1;
		poly_fmpz->length = length + 1;
	} else
	{
		poly_fmpz->length = length;
	   _fmpz_poly_normalise(poly_fmpz);
	}
}

void fmpz_poly_limb_unpack_unsigned(fmpz_poly_t poly_fmpz, const mp_limb_t * array, 
                                  const unsigned long length, const unsigned long limbs)
{
   unsigned long size_m = poly_fmpz->limbs + 1;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
   
   for (unsigned long i = 0, j = 0, k = 0; i < length; i++, j += size_m, k += limbs)
   {
         F_mpn_copy(coeffs_m + j + 1, array + k, limbs);
         coeffs_m[j] = limbs;
         NORM(coeffs_m + j);          
   }
   _fmpz_poly_normalise(poly_fmpz);
}

void __fmpz_poly_write_next_limb(fmpz_t array, unsigned long * temp, unsigned long * offset_limb, 
             const unsigned long next_limb, const unsigned long shift_1, const unsigned long shift_2)
{
   *temp += l_shift(next_limb, shift_1);
   array[*offset_limb] = *temp + ((l_shift(1UL,shift_1)-1)&array[*offset_limb]);
   (*offset_limb)++;
   *temp = r_shift(next_limb, shift_2);
}

void __fmpz_poly_write_whole_limb(fmpz_t array, unsigned long * temp, unsigned long * offset_limb, 
             const unsigned long next_limb, const unsigned long shift_1, const unsigned long shift_2)
{
   *temp += l_shift(next_limb,shift_1);
   array[*offset_limb] = *temp;
   (*offset_limb)++;
   *temp = r_shift(next_limb,shift_2);
}

void fmpz_poly_byte_pack(mp_limb_t * array, const fmpz_poly_t poly_fmpz,
                   const unsigned long length, const unsigned long coeff_bytes, 
                                                                const long negate)
{
   unsigned long size_m = poly_fmpz->limbs+1;
   fmpz_t coeff_m = poly_fmpz->coeffs;
    
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   // Start limb of the current coefficient within array
   unsigned long coeff_limb;
   // Additional offset in bytes of current coefficient within array
   unsigned long coeff_byte;
    
   // Where we are up to in the current coefficient: limbs + bytes
   unsigned long offset_limb;
    
   // Next limb to be written to bytes
   unsigned long next_limb;
   unsigned long temp;
   unsigned long extend;
   
   unsigned long shift_1, shift_2;
   
   fmpz_t scratch = (fmpz_t) flint_stack_alloc(size_m);
   
   fmpz_t co;
    
   // when a coefficient is negative, we need to borrow from the next coefficient
   int borrow;
    
   unsigned long j;
    
   coeff_limb = 0;
   coeff_byte = 0;
   offset_limb = 0;
   temp = 0;
   borrow = 0;
            
   fmpz_t next_point = coeff_m + size_m*length;
         
   while (coeff_m < next_point)
   {
       // compute shifts to be used
       shift_1 = coeff_byte << 3;
       shift_2 = FLINT_BITS - shift_1;
                     
       /* Coefficient is negative after borrow */
       if (((negate > 0L) && ((long) coeff_m[0] - (long) borrow < 0L)) || ((negate < 0L) && ((long) -coeff_m[0] - (long) borrow < 0L)))
       {
          // mpz_t's store the absolute value only, so add 1 then complement
          if (borrow) 
          {
             if (coeff_m[0] == 0) next_limb = ~0L;
             else next_limb = ~coeff_m[1];
             co = coeff_m;
          } else 
          {
             if (negate > 0L) fmpz_add_ui(scratch, coeff_m, 1L);
             else fmpz_sub_ui(scratch, coeff_m, 1L);
             if (scratch[0] == 0) next_limb = ~0L;
             else next_limb = ~scratch[1];
             co = scratch;
          }
          // deal with first limb of coefficient
          if (limbs_per_coeff == 0) 
          {
             if (coeff_m == next_point - size_m) 
             {
                __fmpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                temp += l_shift(-1UL, shift_1);
                array[offset_limb] = temp;
                offset_limb++;
                extend = -1L;
             } else 
             {
                next_limb &= ((1UL<<(extra_bytes_per_coeff<<3))-1);
                __fmpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                array[offset_limb] = temp;
             }
          } else
          {
             __fmpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             // deal with remaining limbs
             for (j = 1; j < ABS(co[0]); j++)
             {
                next_limb = ~co[j+1];
                __fmpz_poly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             }
             // write remaining part of coefficient and fill 
             // remainder of coeff_bytes with binary 1's
             if ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte) 
             {
                temp += l_shift(-1UL,shift_1);
                array[offset_limb] = temp;
                offset_limb++;
             }
             for ( ; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = -1UL;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = -1UL;
                offset_limb++;
             }
             extend = -1L;
          }
          temp = 0;
             
          borrow = 1;
       }        
       else 
       {
          if (borrow) 
          {
             if (negate > 0L) fmpz_sub_ui(scratch, coeff_m, 1L);
             else fmpz_add_ui(scratch, coeff_m, 1L); 
             co = scratch;
          } else 
          {
             co = coeff_m;
          }     
          
          /* Coefficient is positive after borrow */
          if (co[0] != 0L)
          {
             // deal with first limb of coefficient
             next_limb = co[1];
             __fmpz_poly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             if (shift_2 == FLINT_BITS) temp = 0;
             // deal with remaining limbs
             for (j = 1; j < ABS(co[0]); j++)
             {
                next_limb = co[j+1];
                __fmpz_poly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             }
             // write remaining part of coefficient
             array[offset_limb] = temp;
             offset_limb++;
             for (; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = 0UL;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = 0UL;
                offset_limb++;
             }
             extend = 0L;
             temp = 0;
             borrow = 0;
          }
          /* Coefficient is zero after borrow */
          else 
          {
             array[offset_limb] = ((l_shift(1UL,shift_1)-1)&array[offset_limb]);
             offset_limb++;
             for ( ; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = 0UL;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                 limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = 0UL;
                offset_limb++;
             }

             extend = 0L;
             temp = 0;
             borrow = 0;
          }
       }
       // Extend sign of final input coefficient to end of output coefficient 
       /*if (coeff_m == next_point-size_m)
       {
          while (offset_limb < total_limbs)
          {
             array[offset_limb] = extend;
             offset_limb++;
          } 
       }*/
       // update information for next coefficient
       coeff_limb += limbs_per_coeff;
       coeff_byte += extra_bytes_per_coeff;
       if (coeff_byte > FLINT_BYTES_PER_LIMB) 
       {
          coeff_byte -= FLINT_BYTES_PER_LIMB;
          coeff_limb++;
       }
       offset_limb = coeff_limb;
          
       coeff_m += size_m;
   }

	flint_stack_release();
}
     
static inline void __fmpz_poly_unpack_bytes(mp_limb_t * output, const mp_limb_t * array, 
            const unsigned long limb_start, const unsigned long byte_start, 
            const unsigned long num_bytes)
{
    const unsigned long limbs_to_extract = (num_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const unsigned long extra_bytes_to_extract = num_bytes 
                                    - (limbs_to_extract<<FLINT_LG_BYTES_PER_LIMB);
                                    
    unsigned long next_limb;
    unsigned long temp = 0;
    
    // the limb we are up to in the array and output respectively
    unsigned long coeff_limb = limb_start;
    unsigned long output_limb = 0;

    unsigned long shift_1, shift_2;
    
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS - shift_1;
    
    temp = array[coeff_limb];
    coeff_limb++;
    while (output_limb < limbs_to_extract)
    {
       next_limb = r_shift(temp,shift_1);
       temp = array[coeff_limb];
       coeff_limb++;
       next_limb += l_shift(temp,shift_2);
       output[output_limb] = next_limb;
       output_limb++;
    }
    if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
    {
       next_limb = r_shift(temp,shift_1);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    } else
    {
       next_limb = r_shift(temp,shift_1);
       temp = array[coeff_limb];
       next_limb += l_shift(temp,shift_2);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    }
}

static inline unsigned long __fmpz_poly_unpack_signed_bytes(mp_limb_t * output, 
			               const mp_limb_t * array, const unsigned long limb_start, 
				         const unsigned long byte_start, const unsigned long num_bytes)
{
    const unsigned long limbs_to_extract = (num_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const unsigned long extra_bytes_to_extract = num_bytes 
                                    - (limbs_to_extract<<FLINT_LG_BYTES_PER_LIMB);
                                    
    unsigned long next_limb;
    unsigned long temp = 0;
    
    // the limb we are up to in the array and output respectively
    unsigned long coeff_limb = limb_start;
    unsigned long output_limb = 0;

    unsigned long shift_1, shift_2;
    
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS - shift_1;
    
    unsigned long sign;

    if (byte_start + extra_bytes_to_extract > FLINT_BYTES_PER_LIMB)
    {
       sign = array[limb_start+limbs_to_extract+1]&(1UL<<(((byte_start 
            + extra_bytes_to_extract - FLINT_BYTES_PER_LIMB)<<3)-1));
    } else if (byte_start + extra_bytes_to_extract == FLINT_BYTES_PER_LIMB)
    {
       sign = array[limb_start+limbs_to_extract]&(1UL<<(FLINT_BITS-1));
    } else if (byte_start + extra_bytes_to_extract == 0)
    {
       sign = array[limb_start+limbs_to_extract-1]&(1UL<<(FLINT_BITS-1));
    } else
    {
       sign = array[limb_start+limbs_to_extract]&(1UL<<(((byte_start 
            + extra_bytes_to_extract)<<3)-1));
    }
    
    if (sign)
    {
        temp = ~array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract)
        {
           next_limb = r_shift(temp,shift_1);
           temp = ~array[coeff_limb];
           coeff_limb++;
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
        {
           next_limb = r_shift(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else
        {
           next_limb = r_shift(temp,shift_1);
           temp = ~array[coeff_limb];
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    else
    {
        temp = array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract)
        {
           next_limb = r_shift(temp,shift_1);
           temp = array[coeff_limb];
           coeff_limb++;
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
        {
           next_limb = r_shift(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else
        {
           next_limb = r_shift(temp,shift_1);
           temp = array[coeff_limb];
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    return sign;
}   

void fmpz_poly_byte_unpack_unsigned(fmpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long coeff_bytes)
{
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   const unsigned long limbs = ((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
   
   mp_limb_t* temp = (mp_limb_t*) flint_stack_alloc(limbs+2);
   
   unsigned long limb_upto = 0;
   unsigned long byte_offset = 0;
   
   fmpz_t coeff_m = poly_m->coeffs;
   unsigned long size_m = poly_m->limbs+1;
   poly_m->length = length;
   
   for (unsigned long i = 0; i < length; i++)
   {
       F_mpn_clear(temp, limbs+2);
       __fmpz_poly_unpack_bytes(temp + 1, array, limb_upto, 
                                             byte_offset, coeff_bytes);
       temp[0] = limbs;
       NORM(temp);
       
       fmpz_add(coeff_m, coeff_m, temp);
       
       limb_upto += limbs_per_coeff;
       
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
       coeff_m += size_m;
   }                        
   flint_stack_release();
   _fmpz_poly_normalise(poly_m);
}

void fmpz_poly_byte_unpack(fmpz_poly_t poly_m, const mp_limb_t * array,
                               const unsigned long length, const unsigned long coeff_bytes)
{
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   const unsigned long limbs = ((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
   
   mp_limb_t* temp = (mp_limb_t*) flint_stack_alloc(limbs+2);
   
   unsigned long limb_upto = 0;
   unsigned long byte_offset = 0;
   
   unsigned long sign;
   unsigned long borrow = 0;
   
   fmpz_t coeff_m = poly_m->coeffs;
   unsigned long size_m = poly_m->limbs+1;
   poly_m->length = length;
   
   for (unsigned long i = 0; i < length; i++)
   {
       F_mpn_clear(temp,limbs+2);
       sign = __fmpz_poly_unpack_signed_bytes(temp + 1, array, limb_upto, 
                                             byte_offset, coeff_bytes);
       if (sign) temp[0] = -limbs;
       else temp[0] = limbs;
       NORM(temp);
       
       if (sign)
       {
          fmpz_sub_ui_inplace(temp, 1);
       }
       if (borrow) fmpz_add_ui_inplace(temp, 1);
       
       fmpz_add(coeff_m, coeff_m, temp);
       
       borrow = 0;
       if (sign) borrow = 1;
       
       limb_upto += limbs_per_coeff;
       
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
       coeff_m += size_m;
   }                        
   flint_stack_release();
   _fmpz_poly_normalise(poly_m);
}
  
void fmpz_poly_to_zmod_poly(zmod_poly_t zpol, const fmpz_poly_t fpol)
{
   unsigned long p = zpol->p;

   if (fpol->length == 0) 
   {
      zmod_poly_zero(zpol);
      return;
   } 
   
   zmod_poly_fit_length(zpol, fpol->length);
   
   unsigned long sizef = fpol->limbs+1;
   fmpz_t fcoeff = fpol->coeffs;
   unsigned long * zcoeff = zpol->coeffs;

   for (unsigned long i = 0; i < fpol->length; i++)
   {
      zcoeff[i] = fmpz_mod_ui(fcoeff, p);
      fcoeff += sizef;
   }

   zpol->length = fpol->length;
   __zmod_poly_normalise(zpol);
}

void fmpz_poly_to_zmod_poly_no_red(zmod_poly_t zpol, const fmpz_poly_t fpol)
{
   unsigned long p = zpol->p;

   if (fpol->length == 0) 
   {
      zmod_poly_zero(zpol);
      return;
   } 
   
   zmod_poly_fit_length(zpol, fpol->length);
   
   unsigned long sizef = fpol->limbs+1;
   fmpz_t fcoeff = fpol->coeffs;
   unsigned long * zcoeff = zpol->coeffs;
  
   for (unsigned long i = 0; i < fpol->length; i++)
   {
      if (fcoeff[0] == 0) zcoeff[i] = 0;
      else if ((long)fcoeff[0] < 0L) zcoeff[i] = p - fcoeff[1];
      else zcoeff[i] = fcoeff[1];
      fcoeff += sizef;
   }

   zpol->length = fpol->length;
   __zmod_poly_normalise(zpol);
}

void zmod_poly_to_fmpz_poly_unsigned(fmpz_poly_t fpol, zmod_poly_t zpol)
{
   unsigned long p = zpol->p;

   if (zpol->length == 0) 
   {
      fmpz_poly_zero(fpol);
      return;
   } 

   fmpz_poly_fit_length(fpol, zpol->length);
   fmpz_poly_fit_limbs(fpol, 1);

   unsigned long sizef = fpol->limbs+1;
   fmpz_t fcoeff = fpol->coeffs;
   unsigned long * zcoeff = zpol->coeffs;

   for (unsigned long i = 0; i < zpol->length; i++)
   {
      if (zcoeff[i])
      {
         fcoeff[0] = 1L;
         fcoeff[1] = zcoeff[i];
      } else fcoeff[0] = 0L;
      fcoeff += sizef;
   }
   
   fpol->length = zpol->length;
}

void zmod_poly_to_fmpz_poly(fmpz_poly_t fpol, zmod_poly_t zpol)
{
   unsigned long p = zpol->p;
   unsigned long pdiv2 = p/2;

   if (zpol->length == 0) 
   {
      fmpz_poly_zero(fpol);
      return;
   } 

   fmpz_poly_fit_length(fpol, zpol->length);
   fmpz_poly_fit_limbs(fpol, 1);

   unsigned long sizef = fpol->limbs+1;
   fmpz_t fcoeff = fpol->coeffs;
   unsigned long * zcoeff = zpol->coeffs;

   for (unsigned long i = 0; i < zpol->length; i++)
   {
      if (zcoeff[i] == 0)
      {
         fcoeff[0] = 0L;
      } else if (zcoeff[i] > pdiv2)
      {
         fcoeff[0] = -1L;
         fcoeff[1] = p - zcoeff[i];
      } else
      {
         fcoeff[0] = 1L;
         fcoeff[1] = zcoeff[i];
      }
      fcoeff += sizef;
   }
   
   fpol->length = zpol->length;
}

int __fmpz_poly_CRT_unsigned(fmpz_poly_t out, fmpz_poly_t fpol, zmod_poly_t zpol, fmpz_t newmod, fmpz_t oldmod)
{
   unsigned long p = zpol->p;
   double pre = zpol->p_inv;
   unsigned long c, r1;

   c = fmpz_mod_ui(oldmod, p);
   c = z_invert(c, p);

   unsigned long shortest = (fpol->length < zpol->length) ? fpol->length : zpol->length;

   unsigned limbs = FLINT_ABS(oldmod[0]) + 1;
   
   fmpz_poly_fit_length(out, FLINT_MAX(fpol->length, zpol->length));
   fmpz_poly_fit_limbs(out, limbs);

   unsigned long sizef = fpol->limbs+1;
   unsigned long sizeo = out->limbs+1;
   fmpz_t fcoeff = fpol->coeffs;
   fmpz_t ocoeff = out->coeffs;
   unsigned long * zcoeff = zpol->coeffs;
  
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) > FLINT_D_BITS-1)
   {
      for (unsigned long i = 0; i < shortest; i++)
      {
         fmpz_CRT_ui2_precomp(ocoeff, fcoeff, oldmod, zcoeff[i], p, c, pre);
         fcoeff += sizef;
         ocoeff += sizeo;
      }
   } else
#endif 
   for (unsigned long i = 0; i < shortest; i++)
   {
      fmpz_CRT_ui_precomp(ocoeff, fcoeff, oldmod, zcoeff[i], p, c, pre);
      fcoeff += sizef;
      ocoeff += sizeo;
   }

   /* fpol is longer */
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) > FLINT_D_BITS-1)
   {
      for (unsigned long i = shortest; i < fpol->length; i++)
      {
         fmpz_CRT_ui2_precomp(ocoeff, fcoeff, oldmod, 0L, p, c, pre);
         fcoeff += sizef;
         ocoeff += sizeo;
      } 
   } else
#endif 
   for (unsigned long i = shortest; i < fpol->length; i++)
   {
      fmpz_CRT_ui_precomp(ocoeff, fcoeff, oldmod, 0L, p, c, pre);
      fcoeff += sizef;
      ocoeff += sizeo;
   }

   /* zpol is longer */
   unsigned long s;
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) > FLINT_D_BITS-1)
   {
      for (unsigned long i = shortest; i < zpol->length; i++)
      {
         s = z_mulmod2_precomp(zcoeff[i], c, p, pre);
         fmpz_mul_ui(ocoeff, oldmod, s);
         ocoeff += sizeo;
      }
   } else
#endif 
   for (unsigned long i = shortest; i < zpol->length; i++)
   {
      s = z_mulmod_precomp(zcoeff[i], c, p, pre);
      fmpz_mul_ui(ocoeff, oldmod, s);
      ocoeff += sizeo;
   }

   out->length = FLINT_MAX(fpol->length, zpol->length);
   _fmpz_poly_normalise(out);
	fmpz_mul_ui(newmod, oldmod, p);
   return fmpz_poly_equal(fpol, out);
}

int fmpz_poly_CRT_unsigned(fmpz_poly_t res, fmpz_poly_t fpol, zmod_poly_t zpol, fmpz_t newmod, fmpz_t oldmod)
{
	fmpz_poly_t out;

	if (res == fpol)
	{
		fmpz_poly_init(out);
		int same = __fmpz_poly_CRT_unsigned(out, fpol, zpol, newmod, oldmod);
		fmpz_poly_swap(out, res);
		fmpz_poly_clear(out);
		return same;
	} else
		return __fmpz_poly_CRT_unsigned(res, fpol, zpol, newmod, oldmod);
}

int __fmpz_poly_CRT(fmpz_poly_t out, fmpz_poly_t fpol, zmod_poly_t zpol, fmpz_t newmod, fmpz_t oldmod)
{
   unsigned long p = zpol->p;
   double pre = zpol->p_inv;
   unsigned long c, r1;

   c = fmpz_mod_ui(oldmod, p);
   c = z_invert(c, p);

   unsigned long shortest = (fpol->length < zpol->length) ? fpol->length : zpol->length;
   
   fmpz_t newm;
	
	if (newmod == oldmod)
	   newm = fmpz_init(FLINT_ABS(newm[0]) + 1);
	else 
	   newm = newmod;
		
	fmpz_mul_ui(newm, oldmod, p);

   unsigned limbs = FLINT_ABS(newm[0]);
	
	fmpz_poly_fit_length(out, FLINT_MAX(fpol->length, zpol->length));
   fmpz_poly_fit_limbs(out, limbs);

   unsigned long sizef = fpol->limbs+1;
   unsigned long sizeo = out->limbs+1;
   fmpz_t fcoeff = fpol->coeffs;
   fmpz_t ocoeff = out->coeffs;
   unsigned long * zcoeff = zpol->coeffs;
   fmpz_t moddiv2 = fmpz_init(newm[0]);

   fmpz_div_2exp(moddiv2, newm, 1);
  
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) > FLINT_D_BITS-1)
   {
      for (unsigned long i = 0; i < shortest; i++)
      {
         fmpz_CRT_ui2_precomp(ocoeff, fcoeff, oldmod, zcoeff[i], p, c, pre);
         if (fmpz_cmpabs(ocoeff, moddiv2) > 0)
            fmpz_sub(ocoeff, ocoeff, newm);
         fcoeff += sizef;
         ocoeff += sizeo;
      }
   } else
#endif 
   for (unsigned long i = 0; i < shortest; i++)
   {
      fmpz_CRT_ui_precomp(ocoeff, fcoeff, oldmod, zcoeff[i], p, c, pre);
      if (fmpz_cmpabs(ocoeff, moddiv2) > 0)
         fmpz_sub(ocoeff, ocoeff, newm);
      fcoeff += sizef;
      ocoeff += sizeo;
   }
   
   /* fpol is longer */
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) > FLINT_D_BITS-1)
   {
      for (unsigned long i = shortest; i < fpol->length; i++)
      {
         fmpz_CRT_ui2_precomp(ocoeff, fcoeff, oldmod, 0L, p, c, pre);
         if (fmpz_cmpabs(ocoeff, moddiv2) > 0)
            fmpz_sub(ocoeff, ocoeff, newm);
         fcoeff += sizef;
         ocoeff += sizeo;
      } 
   } else
#endif 
   for (unsigned long i = shortest; i < fpol->length; i++)
   {
      fmpz_CRT_ui_precomp(ocoeff, fcoeff, oldmod, 0L, p, c, pre);
      if (fmpz_cmpabs(ocoeff, moddiv2) > 0)
         fmpz_sub(ocoeff, ocoeff, newm);
      fcoeff += sizef;
      ocoeff += sizeo;
   }

   /* zpol is longer */
   unsigned long s;
#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) > FLINT_D_BITS-1)
   {
      for (unsigned long i = shortest; i < zpol->length; i++)
      {
         s = z_mulmod2_precomp(zcoeff[i], c, p, pre);
         fmpz_mul_ui(ocoeff, oldmod, s);
         if (fmpz_cmpabs(ocoeff, moddiv2) > 0)
            fmpz_sub(ocoeff, ocoeff, newm);
         ocoeff += sizeo;
      }
   } else
#endif 
   for (unsigned long i = shortest; i < zpol->length; i++)
   {
      s = z_mulmod_precomp(zcoeff[i], c, p, pre);
      fmpz_mul_ui(ocoeff, oldmod, s);
      if (fmpz_cmpabs(ocoeff, moddiv2) > 0)
         fmpz_sub(ocoeff, ocoeff, newm);
      ocoeff += sizeo;
   }

   out->length = FLINT_MAX(fpol->length, zpol->length);
   _fmpz_poly_normalise(out);
   
   fmpz_clear(moddiv2); 

   if (newmod == oldmod)
	{
		fmpz_set(newmod, newm);
	   fmpz_clear(newm);
	} 

	return fmpz_poly_equal(fpol, out);
}

int fmpz_poly_CRT(fmpz_poly_t res, fmpz_poly_t fpol, zmod_poly_t zpol, fmpz_t newmod, fmpz_t oldmod)
{
	fmpz_poly_t out;

	if (res == fpol)
	{
		fmpz_poly_init(out);
		int same = __fmpz_poly_CRT(out, fpol, zpol, newmod, oldmod);
		fmpz_poly_swap(out, res);
		fmpz_poly_clear(out);
		return same;
	} else
		return __fmpz_poly_CRT(res, fpol, zpol, newmod, oldmod);
}

/****************************************************************************

   _fmpz_poly_* layer
   
   All inputs are assumed to be normalised. All outputs are
   normalised given this assumption.

****************************************************************************/

/* 
   Create a polynomial of length zero with "alloc" allocated coefficients
   each with enough space for limbs limbs
*/

void _fmpz_poly_stack_init(fmpz_poly_t poly, const unsigned long alloc, const unsigned long limbs)
{
   if ((alloc) && (limbs)) poly->coeffs = (fmpz_t) flint_stack_alloc(alloc*(limbs+1));
   else poly->coeffs = NULL;
   
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

void _fmpz_poly_stack_clear(fmpz_poly_t poly)
{
   if (poly->coeffs) flint_stack_release();
   poly->coeffs = NULL;
}

void _fmpz_poly_check(const fmpz_poly_t poly)
{
  if ((long) poly->length < 0)
   {
      printf("Error: Poly length < 0\n");
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

void _fmpz_poly_check_normalisation(const fmpz_poly_t poly)
{
   if (poly->length)
   {
      if (!poly->coeffs[(poly->length-1)*(poly->limbs+1)])
      {
         printf("Error: Poly not normalised\n");
         abort();
      }
   }
   if ((long) poly->length < 0)
   {
      printf("Error: Poly length < 0\n");
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


// Retrieves the n-th coefficient as an mpz
void _fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length) 
   {
      mpz_set_ui(x, 0);
      return;
   }
   fmpz_to_mpz(x, poly->coeffs + n*(poly->limbs + 1));
}


/* 
   Set a coefficient to the given unsigned value.
   If x is nonzero, poly->limbs must be positive.
   Assumes the polynomial length is greater than n.
*/
void _fmpz_poly_set_coeff_ui(fmpz_poly_t poly, const unsigned long n, const unsigned long x)
{
   FLINT_ASSERT(poly->length > n);
   fmpz_set_ui(poly->coeffs + n*(poly->limbs + 1), x);
   if ((x==0L) && (poly->length == n+1)) _fmpz_poly_normalise(poly);
}

/* 
   Set a coefficient to the given signed value.
   If x is nonzero, poly->limbs must be positive.
   Assumes the polynomial length is greater than n.
   Normalises only if the leading coefficient is 
   set to zero.
*/

void _fmpz_poly_set_coeff_si(fmpz_poly_t poly, const unsigned long n, const long x)
{
   FLINT_ASSERT(poly->length > n);
   fmpz_set_si(poly->coeffs + n*(poly->limbs + 1), x);
   if ((x==0L) && (poly->length == n+1)) _fmpz_poly_normalise(poly);
}

/* 
   Set a coefficient to the given fmpz_t.
   Assumes the polynomial length is greater than n.
   Normalises only if the leading coefficient is 
   set to zero.
*/

void _fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, const unsigned long n, fmpz_t x)
{
   FLINT_ASSERT(poly->length > n);
   fmpz_set(poly->coeffs + n*(poly->limbs + 1), x);
   if (fmpz_is_zero(x) && (poly->length == n+1)) _fmpz_poly_normalise(poly);
}

// Retrieves the n-th coefficient as an fmpz
void _fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length) 
   {
      fmpz_set_ui(x, 0);
      return;
   }
   fmpz_set(x, poly->coeffs + n*(poly->limbs + 1));
}

void _fmpz_poly_get_coeff_mpz_read_only(mpz_t x, const fmpz_poly_t poly, const unsigned long n)
{
   mp_limb_t * coeff = _fmpz_poly_get_coeff_ptr(poly, n);
   if (poly->length) 
   {
      x->_mp_d = coeff + 1;
      x->_mp_size = coeff[0];
      x->_mp_alloc = poly->limbs;
   } else 
   {
      x->_mp_d = (mp_limb_t *) &poly; // We need to point to something, and at least this exists
      x->_mp_size = 0;
      x->_mp_alloc = FLINT_MAX(1, poly->limbs);
   }
}

void _fmpz_poly_normalise(fmpz_poly_t poly)
{
   while (poly->length && poly->coeffs[(poly->length-1)*(poly->limbs+1)] == 0)
      poly->length--;
}

/* 
   Sets the output poly to equal the input poly 
   Assumes the output poly is big enough to hold 
   the nonzero limbs of the input poly
*/

void _fmpz_poly_set(fmpz_poly_t output, const fmpz_poly_t input)
{
   if (input->length == 0) 
   {
      output->length = 0;
      return;
   }
   
   if (output != input) 
   {
      unsigned long input_size = input->limbs + 1;
      unsigned long output_size = output->limbs + 1;
      if ((output->coeffs < input->coeffs) || (output->coeffs >= input->coeffs + input->length*(input->limbs+1)))
      {
         for (long i = 0; i < input->length; i++)
         {
            if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
            else F_mpn_copy(output->coeffs+i*output_size, input->coeffs+i*input_size, ABS(input->coeffs[i*input_size])+1);
         }
      } else
      {
         for (long i = input->length - 1; i >= 0; i--)
         {
            if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
            else F_mpn_copy(output->coeffs+i*output_size, input->coeffs+i*input_size, ABS(input->coeffs[i*input_size])+1);
         }
      }
   }
   output->length = input->length;
}

/*
   Determines the maximum number of bits in any coefficient of poly_fmpz. This
   function assumes every coefficient fits in a limb. The returned value is
   negative if any of the coefficients was negative.
*/
long _fmpz_poly_max_bits1(const fmpz_poly_t poly_fmpz)
{
   unsigned long mask = -1L;
   long bits = 0;
   long sign = 1;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
   long i, j;
   
   for (i = 0, j = 0; i < poly_fmpz->length; i++, j += 2)
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
      for ( ; i < poly_fmpz->length; i++, j += 2)
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
   Determines the maximum number of bits in a coefficient of poly_fmpz. 
   The returned value is negative if any of the coefficients was negative.
*/

long _fmpz_poly_max_bits(const fmpz_poly_t poly_fmpz)
{
   if (poly_fmpz->limbs == 0) return 0;
   if (poly_fmpz->limbs == 1) return _fmpz_poly_max_bits1(poly_fmpz);
   
   unsigned long mask = -1L;
   long bits = 0;
   long sign = 1;
   long limbs = 0;
   long size_j;
   fmpz_t coeffs_m = poly_fmpz->coeffs;
   unsigned long size_m = poly_fmpz->limbs+1;
   long i, j;
   
   for (i = 0, j = 0; i < poly_fmpz->length; i++, j += size_m)
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
      for ( ; i < poly_fmpz->length; i++, j += size_m)
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

/*
   Returns the maximum number of limbs of any coefficient of poly.
   Does not count the sign/size limb.
*/

unsigned long _fmpz_poly_max_limbs(const fmpz_poly_t poly)
{
   unsigned long limbs = poly->limbs;
   unsigned long max_limbs = 0;
   unsigned long next_limbs;
   
   for (long i = 0; (i < poly->length) && (max_limbs != limbs); i++)
   {
       next_limbs = ABS(poly->coeffs[i*(limbs+1)]);
       if (next_limbs > max_limbs) max_limbs = next_limbs;
   } 
   return max_limbs;
}

/* 
   Checks if two polynomials are arithmetically equal and
   returns 1 if they are, 0 otherwise.
*/

int _fmpz_poly_equal(const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   if (input1 == input2) return 1;
   if (input1->length != input2->length) return 0;
   
   long i,j;
      
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

/*
   Truncate the polynomial to trunc terms and normalise
*/

void _fmpz_poly_truncate(fmpz_poly_t poly, const unsigned long trunc)
{
   if (poly->length > trunc) poly->length = trunc;
   _fmpz_poly_normalise(poly);
}

/*
   Set output to -input
*/

void _fmpz_poly_neg(fmpz_poly_t output, const fmpz_poly_t input)
{
   if (input == output)
   {
      for (long i = 0; i < input->length; i++)
         output->coeffs[i*(output->limbs+1)] = -output->coeffs[i*(output->limbs+1)];
   } else
   {
      unsigned long input_size = input->limbs + 1;
      unsigned long output_size = output->limbs + 1;
      for (long i = 0; i < input->length; i++)
      {
         if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
         else 
         {
            output->coeffs[i*output_size] = -input->coeffs[i*input_size];
            F_mpn_copy(output->coeffs+i*output_size+1, input->coeffs+i*input_size+1, ABS(input->coeffs[i*input_size]));
         }
      }
   }
   output->length = input->length;
}

void _fmpz_poly_scalar_abs(fmpz_poly_t output, const fmpz_poly_t input)
{
   if (input == output)
   {
      for (long i = 0; i < input->length; i++)
         output->coeffs[i*(output->limbs+1)] = FLINT_ABS(output->coeffs[i*(output->limbs+1)]);
   } else
   {
      unsigned long input_size = input->limbs + 1;
      unsigned long output_size = output->limbs + 1;
      for (long i = 0; i < input->length; i++)
      {
         if (!input->coeffs[i*input_size]) output->coeffs[i*output_size] = 0;
         else 
         {
            output->coeffs[i*output_size] = FLINT_ABS(input->coeffs[i*input_size]);
            F_mpn_copy(output->coeffs+i*output_size+1, input->coeffs+i*input_size+1, ABS(input->coeffs[i*input_size]));
         }
      }
   }
   output->length = input->length;
}

/* 
   Set n of the coefficients of poly to zero starting with
   the constant term *regardless of the original length*. 
   Normalises if n >= poly->length-1
*/

void _fmpz_poly_zero_coeffs(fmpz_poly_t poly, const unsigned long n)
{
   unsigned long size = poly->limbs+1;
   fmpz_t coeff = poly->coeffs;
   for (long i = 0; i < n; i++)
   {
      coeff[0] = 0;
      coeff+=size;
   }
   if (n >= poly->length-1) _fmpz_poly_normalise(poly);
}

/* 
   Multiplies input by x^n and sets output to the result
   Assumes output is large enough to contain the result
*/

void _fmpz_poly_left_shift(fmpz_poly_t output, const fmpz_poly_t input, 
                                                 const unsigned long n)
{
   fmpz_poly_t part;   
   
   part->length = input->length;
   part->limbs = output->limbs;
   part->coeffs = output->coeffs + n*(output->limbs+1);
      
   _fmpz_poly_set(part, input);
	for (long i = 0; i < n; i++) 
	   output->coeffs[i*(output->limbs+1)] = 0;
   
   if (input->length > 0) output->length = input->length + n;
   else (output->length = 0);
}

/* 
   Divides input by x^n losing the remainder and sets output to the result
   Assumes output is large enough to contain the result
*/

void _fmpz_poly_right_shift(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long n)
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

void _fmpz_poly_reverse(fmpz_poly_t output, const fmpz_poly_t input, const unsigned long length)
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
         F_mpn_copy(output->coeffs + (length - i - 1)*size_out, input->coeffs + i*size_in, coeff_limbs);
      }
      for ( ; i < length; i++)
      {
         output->coeffs[(length - i - 1)*size_out] = 0L;
      }
      output->length = length;
      _fmpz_poly_normalise(output);
   } else
   {
      fmpz_t temp = (fmpz_t) flint_stack_alloc(size_in);
      unsigned long coeff_limbs2;
      
      for (i = 0; i < length/2; i++)
      {
         if (i < input->length)
         {
            coeff_limbs = ABS(input->coeffs[i*size_in]) + 1;
            F_mpn_copy(temp, input->coeffs + i*size_in, coeff_limbs);
         } else
         {
            coeff_limbs = 1;
            temp[0] = 0;            
         }
         if (length - i - 1 < input->length)
         {
            coeff_limbs2 = ABS(input->coeffs[(length - i - 1)*size_in]) + 1;
            F_mpn_copy(input->coeffs + i*size_in, input->coeffs + (length - i - 1)*size_in, coeff_limbs2);
         } else
         {
            input->coeffs[i*size_in] = 0;
         }
         F_mpn_copy(input->coeffs + (length - i - 1)*size_in, temp, coeff_limbs);
      }
      if ((length & 1) && (i >= input->length)) input->coeffs[i*size_in] = 0;

      output->length = length;
      _fmpz_poly_normalise(output);
      flint_stack_release();
   }
}


/* 
    Add two polynomials together 
*/

void _fmpz_poly_add(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   if (input1 == input2) 
   {
      _fmpz_poly_scalar_mul_ui(output, input1, 2UL);
      
      return;
   }
   
   unsigned long size1, size2, shorter, size_out;
   fmpz_t coeffs1, coeffs2, coeffs_out;
   
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   
   size_out = output->limbs+1;
   coeffs_out = output->coeffs;
   
   shorter = (input1->length > input2->length) ? input2->length : input1->length;
   
   for (long i = 0; i < shorter; i++)
   {
      fmpz_add(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2+i*size2);
   }    
   
   if (input1 != output)
   {
      for (long i = shorter; i < input1->length; i++)
      {
          F_mpn_copy(coeffs_out+i*size_out, coeffs1+i*size1, ABS(coeffs1[i*size1])+1);
      }
   }
   if (input2 != output)
   {
      for (unsigned long i = shorter; i < input2->length; i++)
      {
         F_mpn_copy(coeffs_out+i*size_out, coeffs2+i*size2, ABS(coeffs2[i*size2])+1);
      }
   }
   
   if (input1->length == input2->length)
   {
      output->length = input1->length;
      _fmpz_poly_normalise(output);
   } else
   {
      output->length = (input1->length > input2->length) ? input1->length : input2->length;
   }
}

/* 
    Subtract two polynomials
*/

void _fmpz_poly_sub(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
{
   if (input1 == input2) 
   {
      _fmpz_poly_zero_coeffs(output, input1->length);
      _fmpz_poly_zero(output);
      
      return;
   }
   
   unsigned long size1, size2, shorter, size_out;
   fmpz_t coeffs1, coeffs2, coeffs_out;
   
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   
   size_out = output->limbs+1;
   coeffs_out = output->coeffs;
   
   shorter = (input1->length > input2->length) ? input2->length : input1->length;
   
   for (long i = 0; i < shorter; i++)
   {
      fmpz_sub(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2+i*size2);
   }
   
   if (input1 != output)
   {
      for (long i = shorter; i < input1->length; i++)
      {
         F_mpn_copy(coeffs_out+i*size_out, coeffs1+i*size1, ABS(coeffs1[i*size1])+1);
      }
   }
   if (input2 != output)
   {
      for (long i = shorter; i < input2->length; i++)
      {
         F_mpn_copy(coeffs_out+i*size_out+1, coeffs2+i*size2+1, ABS(coeffs2[i*size2]));
         coeffs_out[i*size_out] = -coeffs2[i*size2];
      }
   } else
   {
      for (long i = shorter; i < input2->length; i++)
      {
         coeffs_out[i*size_out] = -coeffs2[i*size2];
      }
   }

   if (input1->length == input2->length)
   {
      output->length = input1->length;
      _fmpz_poly_normalise(output);
   } else
   {
      output->length = (input1->length > input2->length) ? input1->length : input2->length;
   }
}

void _fmpz_poly_scalar_mul_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
     if (x == 0)  
     {
        unsigned long size = output->limbs + 1;
        for (long i = 0; i < poly->length; i++)
        {
           output->coeffs[i*size] = 0L;
        }
        output->length = 0;
        return;
     }
     
     fmpz_t coeffs1 = poly->coeffs;
     fmpz_t coeffs_out = output->coeffs;
     unsigned long size1 = poly->limbs+1;
     unsigned long size_out = output->limbs+1;
     mp_limb_t mslimb;
     
     for (long i = 0; i < poly->length; i++)
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

void _fmpz_poly_scalar_mul_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x)
{
     if (x == 0)  
     {
        unsigned long size = output->limbs + 1;
        for (long i = 0; i < poly->length; i++)
        {
           output->coeffs[i*size] = 0L;
        }
        output->length = 0;
        return;
     }
     
     fmpz_t coeffs1 = poly->coeffs;
     fmpz_t coeffs_out = output->coeffs;
     unsigned long size1 = poly->limbs+1;
     unsigned long size_out = output->limbs+1;
     mp_limb_t mslimb;
     
     for (long i = 0; i < poly->length; i++)
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

/*
   Scalar multiplication of a polynomial by a scalar
*/

void _fmpz_poly_scalar_mul_fmpz(fmpz_poly_t output, const fmpz_poly_t poly, const fmpz_t x)
{
   if (poly->length == 0)
   {
      output->length = 0;
      return;
   }   
   if (x[0] == 0)
   {
      unsigned long size = output->limbs + 1;
      for (long i = 0; i < poly->length; i++)
      {
         output->coeffs[i*size] = 0L;
      }
      output->length = 0;
      return;
   }   
   
   unsigned long x0 = ABS(x[0]);
   while ((!x[x0]) && (x0)) x0--;
   unsigned long limbs1 = x0;
   unsigned long limbs2 = poly->limbs;
   unsigned long total_limbs;
   unsigned long msl;
   unsigned long limbs_out = output->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs2 = poly->coeffs;
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
          } else coeffs_out[i*limbs_out] = 0L;
      }
   } else if (limbs1 + limbs2 > 1000)
   {
      F_mpn_precomp_t precomp;
   
      F_mpn_mul_precomp_init(precomp, x+1, limbs1, limbs2);   
             
      for (long i = 0; i < poly->length; i++)
      {
          total_limbs = limbs1 + ABS(coeffs2[i*(limbs2+1)]);
          if (total_limbs != limbs1)
          {
             msl = F_mpn_mul_precomp(coeffs_out + i*limbs_out + 1, coeffs2 + i*(limbs2+1) + 1, ABS(coeffs2[i*(limbs2+1)]), precomp);
             if (((long) coeffs2[i*(limbs2+1)] ^ sign1) < 0) coeffs_out[i*limbs_out] = -total_limbs + (msl == 0L);
             else coeffs_out[i*limbs_out] = total_limbs - (msl == 0L);
          } else coeffs_out[i*limbs_out] = 0L;
      }
      F_mpn_mul_precomp_clear(precomp);
   } else
   {
      if (poly != output)
      {
         for (long i = 0; i < poly->length - 1; i++)
         {
            __fmpz_mul(coeffs_out + i*limbs_out, coeffs2 + i*(limbs2+1), x);
         }
         fmpz_mul(coeffs_out + (poly->length - 1)*limbs_out, coeffs2 + (poly->length - 1)*(limbs2+1), x);
      } else
      {
         for (long i = 0; i < poly->length; i++)
         {


            fmpz_mul(coeffs_out + i*limbs_out, coeffs2 + i*(limbs2+1), x);
         }
      }
   } 
   output->length = poly->length;
}

void _fmpz_poly_scalar_div_exact_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0) 
   {
      output->length = 0;
      return;
   }
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1 = poly->coeffs;
      
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
            F_mpn_copy(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]));
            F_mpn_clear(coeffs_out+i*size_out+ABS(coeffs1[i*size1])+1, size_out-ABS(coeffs1[i*size1]));
         } 
         F_mpn_copy(coeffs_out+(poly->length-1)*size_out+1, coeffs1+(poly->length-1)*size1+1, ABS(coeffs1[(poly->length-1)*size1]));
         if (size_out > ABS(coeffs1[(poly->length-1)*size1])+1) F_mpn_clear(coeffs_out+(poly->length-1)*size_out+ABS(coeffs1[(poly->length-1)*size1])+1, size_out-ABS(coeffs1[(poly->length-1)*size1])-1);
         
         mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (unsigned long i = 0; i < poly->length; i++)
         {
            coeffs_out[i*size_out] = coeffs1[i*size1];
            NORM(coeffs_out+i*size_out);
         }
      } else
      {
         fmpz_t signs = (fmpz_t) flint_stack_alloc(poly->length);
         signs[0] = coeffs1[0];
         coeffs_out[0] = 0;
         for (unsigned long i = 0; i < poly->length-1; i++)
         {
             signs[i+1] = coeffs1[(i+1)*size1];
             F_mpn_clear(coeffs_out+i*size_out+ABS(signs[i])+1, size_out-ABS(signs[i]));
         } 
         if (size_out > ABS(signs[poly->length-1])+1) F_mpn_clear(coeffs_out+(poly->length-1)*size_out+ABS(signs[poly->length-1])+1, size_out-ABS(signs[poly->length-1])-1);
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

void _fmpz_poly_scalar_div_exact_si(fmpz_poly_t output, const fmpz_poly_t poly, const long x)
{
   if (poly->length == 0) 
   {
      output->length = 0;
      return;
   }
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1 = poly->coeffs;
      
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
            F_mpn_copy(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]));
            F_mpn_clear(coeffs_out+i*size_out+ABS(coeffs1[i*size1])+1, size_out-ABS(coeffs1[i*size1]));
         } 
         F_mpn_copy(coeffs_out+(poly->length-1)*size_out+1, coeffs1+(poly->length-1)*size1+1, ABS(coeffs1[(poly->length-1)*size1]));
         if (size_out > ABS(coeffs1[(poly->length-1)*size1])+1) F_mpn_clear(coeffs_out+(poly->length-1)*size_out+ABS(coeffs1[(poly->length-1)*size1])+1, size_out-ABS(coeffs1[(poly->length-1)*size1])-1);
         
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
         fmpz_t signs = (fmpz_t) flint_stack_alloc(poly->length);
         signs[0] = coeffs1[0];
         coeffs_out[0] = 0;
         for (long i = 0; i < poly->length-1; i++)
         {
             signs[i+1] = coeffs1[(i+1)*size1];
             F_mpn_clear(coeffs_out+i*size_out+ABS(signs[i])+1, size_out-ABS(signs[i]));
         } 
         if (size_out > ABS(signs[poly->length-1])+1) F_mpn_clear(coeffs_out+(poly->length-1)*size_out+ABS(signs[poly->length-1])+1, size_out-ABS(signs[poly->length-1])-1);
         if (x < 0) mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, -x);
         else mpn_divmod_1(coeffs_out, coeffs_out, size_out*poly->length, x);
         for (long i = 0; i < poly->length; i++)
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
    Does scalar division of a polynomial by a limb x. Rounding is done towards zero.
*/

void _fmpz_poly_scalar_tdiv_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0) 
   {
      output->length = 0;
      return;
   }
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1 = poly->coeffs;
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      unsigned long xnorm;
      
      count_lead_zeros(norm, x);
      xnorm = (x<<norm);
      xinv = F_mpn_precompute_inverse(xnorm);
      
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         F_mpn_divrem_ui_precomp(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x, xinv);
         NORM(coeffs_out+i*size_out);
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
         NORM(coeffs_out+i*size_out);
      }
   }
   
   output->length = poly->length;
   _fmpz_poly_normalise(output);
}

/* 
    Does scalar division of a polynomial by a limb x. Rounding is done towards 
    minus infinity so that the remainder is positive.
*/

void _fmpz_poly_scalar_div_ui(fmpz_poly_t output, const fmpz_poly_t poly, const unsigned long x)
{
   if (poly->length == 0) 
   {
      output->length = 0;
      return;
   }
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1 = poly->coeffs;
   
   mp_limb_t rem;
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      unsigned long xnorm;
      
      count_lead_zeros(norm, x);
      xnorm = (x<<norm);
      xinv = F_mpn_precompute_inverse(xnorm);
      
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         rem = F_mpn_divrem_ui_precomp(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x, xinv);
         if (((long) coeffs_out[i*size_out] < 0L) && (rem))
         {
            NORM(coeffs_out+i*size_out);
            fmpz_sub_ui_inplace(coeffs_out+i*size_out, 1UL);
         } else
         {
            NORM(coeffs_out+i*size_out);
         }     
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         coeffs_out[i*size_out] = coeffs1[i*size1];
         rem = mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
         if (((long) coeffs_out[i*size_out] < 0L) && (rem))
         {
            NORM(coeffs_out+i*size_out);
            fmpz_sub_ui_inplace(coeffs_out+i*size_out, 1UL);
         } else
         {
            NORM(coeffs_out+i*size_out);
         }     
      }
   }
   
   output->length = poly->length;
   _fmpz_poly_normalise(output);
}

/*
   Divide each coefficient by the signed scalar, rounding the quotient towards zero
*/

void _fmpz_poly_scalar_tdiv_si(fmpz_poly_t output, const fmpz_poly_t poly, const long scalar)
{
   long x = scalar;
   
   if (poly->length == 0) 
   {
      output->length = 0;
      return;
   }
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1 = poly->coeffs;
   int sign = (x < 0);
   if (sign) x = -x; 
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      unsigned long xnorm;
      
      count_lead_zeros(norm, (unsigned long) x);
      xnorm = ((unsigned long) x<<norm);
      xinv = F_mpn_precompute_inverse(xnorm);
      
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (sign) coeffs_out[i*size_out] = -coeffs1[i*size1];
         else coeffs_out[i*size_out] = coeffs1[i*size1];
         F_mpn_divrem_ui_precomp(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x, xinv);
         NORM(coeffs_out+i*size_out);
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (sign) coeffs_out[i*size_out] = -coeffs1[i*size1];
         else coeffs_out[i*size_out] = coeffs1[i*size1];
         mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
         NORM(coeffs_out+i*size_out);
      }
   }
   
   output->length = poly->length;
   _fmpz_poly_normalise(output);
}

/*
   Divide each coefficient by the signed scalar, rounding the quotient towards minus infinity
*/

void _fmpz_poly_scalar_div_si(fmpz_poly_t output, const fmpz_poly_t poly, const long scalar)
{
   long x = scalar;
   
   if (poly->length == 0) 
   {
      output->length = 0;
      return;
   }
   unsigned long size_out = output->limbs+1;
   unsigned long size1 = poly->limbs+1;
   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1 = poly->coeffs;
   int sign = (x < 0L);
   if (sign) x = -x; 
   mp_limb_t rem;
      
   if (poly->length > FLINT_POL_DIV_1_LENGTH)
   {
      unsigned long norm;
      mp_limb_t xinv;
      unsigned long xnorm;
      
      count_lead_zeros(norm, (unsigned long) x);
      xnorm = ((unsigned long) x<<norm);
      xinv = F_mpn_precompute_inverse(xnorm);
      
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (sign) coeffs_out[i*size_out] = -coeffs1[i*size1];
         else coeffs_out[i*size_out] = coeffs1[i*size1];
         rem = F_mpn_divrem_ui_precomp(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x, xinv);
         if (((long) coeffs_out[i*size_out] < 0L) && (rem))
         {
            NORM(coeffs_out+i*size_out);
            fmpz_sub_ui_inplace(coeffs_out+i*size_out, 1UL);
         } else
         {
            NORM(coeffs_out+i*size_out);
         }     
      }
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         if (sign) coeffs_out[i*size_out] = -coeffs1[i*size1];
         else coeffs_out[i*size_out] = coeffs1[i*size1];
         rem = mpn_divmod_1(coeffs_out+i*size_out+1, coeffs1+i*size1+1, ABS(coeffs1[i*size1]), x);
         if (((long) coeffs_out[i*size_out] < 0L) && (rem))
         {
            NORM(coeffs_out+i*size_out);
            fmpz_sub_ui_inplace(coeffs_out+i*size_out, 1UL);
         } else
         {
            NORM(coeffs_out+i*size_out);
         }
      }
   }
   
   output->length = poly->length;
   _fmpz_poly_normalise(output);
}

/*
   Divide each coefficient of poly by scalar. Rounds towards minus infinity.
*/

void _fmpz_poly_scalar_div_fmpz(fmpz_poly_t output, const fmpz_poly_t poly, const fmpz_t scalar)
{
   if (scalar[0] == 1L) 
   {
      _fmpz_poly_scalar_div_ui(output, poly, scalar[1]);
      return;
   }
   
   if ((scalar[0] == -1L) && (fmpz_bits(scalar) < FLINT_BITS))
   {
      _fmpz_poly_scalar_div_si(output, poly, -scalar[1]);
      return;
   }
   
   if (poly == output)
   {
      fmpz_poly_t temp;
      fmpz_poly_init(temp);
      fmpz_poly_set(temp, poly);
      for (unsigned long i = 0; i < temp->length; i++)
      {
         fmpz_fdiv(output->coeffs+i*(output->limbs+1), temp->coeffs+i*(temp->limbs+1), scalar);
      }
      fmpz_poly_clear(temp);
   } else
   {
      for (unsigned long i = 0; i < poly->length; i++)
      {
         fmpz_fdiv(output->coeffs+i*(output->limbs+1), poly->coeffs+i*(poly->limbs+1), scalar);
      }
   }
   
   output->length = poly->length;
   _fmpz_poly_normalise(output);
}

/*
   Multiply two polynomials using the classical technique.
   Currently doesn't allow aliasing
*/

void _fmpz_poly_mul_classical(fmpz_poly_t output, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   fmpz_poly_t input1, input2;
   
   if (output == poly1)
   {
      _fmpz_poly_stack_init(input1, poly1->length, poly1->limbs);
      _fmpz_poly_set(input1, poly1);
      if (output == poly2)
      {
         _fmpz_poly_attach(input2, input1);
      } else _fmpz_poly_attach(input2, poly2);
   } else if (output == poly2)
   {
      _fmpz_poly_stack_init(input2, poly2->length, poly2->limbs);
      _fmpz_poly_set(input2, poly2);
      _fmpz_poly_attach(input1, poly1);
   } else
   {
      _fmpz_poly_attach(input1, poly1);
      _fmpz_poly_attach(input2, poly2);
   }

   fmpz_t coeffs_out = output->coeffs;
   fmpz_t coeffs1, coeffs2;
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
         fmpz_mul(coeffs_out, coeffs1, coeffs2);
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
      
      fmpz_t temp;
      long i, j;
      
      for (i = 0; i < len1; i++)
      {
         /* Set out[i] = in1[i]*in2[0] */
         if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
         {
            coeffs_out[i*size_out] = 0;
         } else
         {
            __fmpz_mul(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
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
            __fmpz_mul(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
         }      
      }
      /* 
         The above coefficient multiplications overwrite the first limb of the next coefficient
         in each case, using the function __fmpz_mul. The final multiplication 
         cannot do this however.
      */
      if ((coeffs1[lenm1*size1] == 0) || (coeffs2[(len2-1)*size2] == 0))
      {
         coeffs_out[(len2+lenm1-1)*size_out]=0;
      } else
      {
         fmpz_mul(coeffs_out+(len2+lenm1-1)*size_out, coeffs1+lenm1*size1, coeffs2+(len2-1)*size2);
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
                  fmpz_mul(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } else 
               {
                  fmpz_addmul(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } 
            }
         }
      }
   } 
   
   output->length = len1 + len2 - 1;
   
   if (poly1 == output) _fmpz_poly_stack_clear(input1);
   else if (poly2 == output) _fmpz_poly_stack_clear(input2);
}

#ifdef HAVE_ZNPOLY
/*
   Multiply two polynomials using the multimodular technique.
   This version requires an appropriate fmpz_comb_t struct, already initialized with primes
   This function allows aliasing
*/
void __fmpz_poly_mul_modular_comb(fmpz_poly_t output, const fmpz_poly_t poly1, const fmpz_poly_t poly2,
        fmpz_comb_t comb)
{
    if ((poly1->length == 0) || (poly2->length == 0))
    {
        _fmpz_poly_zero(output);
        return;
    }

    // Multi-reduce poly1, place result into block1
    unsigned long len1 = poly1->length;
    unsigned long * block1 = flint_heap_alloc(len1 << comb->n);
    for(int i=0; i<len1; i++) {
        fmpz_multi_mod_ui(block1 + (i << comb->n),
                _fmpz_poly_get_coeff_ptr(poly1, i), comb);
    }
    
    // Multi-reduce poly2, place result into block2
    unsigned long len2 = poly2->length;
    unsigned long * block2 = flint_heap_alloc(len2 << comb->n);
    for(int i=0; i<len2; i++) {
        fmpz_multi_mod_ui(block2 + (i << comb->n),
                _fmpz_poly_get_coeff_ptr(poly2, i), comb);
    }
    
    unsigned long len_out = len1 + len2 - 1;
    unsigned long * block_out = flint_heap_alloc(len_out << comb->n);

    // Inputs for zn_array_mul
    unsigned long * in1 = flint_heap_alloc(len1);
    unsigned long * in2 = flint_heap_alloc(len2);

    // Output for zn_array_mul
    unsigned long * out = flint_heap_alloc(len_out);

    unsigned long numprimes = comb->num_primes;

    // FIXME: reorganize this loop to optimize cache line usage
    
    for(int i=0; i<numprimes; i++) {

        // in1 := poly1 % comb->primes[i]
        for(int j=0; j<len1; j++) {
            in1[j] = block1[i + (j << comb->n)];
        }

        // in2 := poly2 % comb->primes[i]
        for(int j=0; j<len2; j++) {
            in2[j] = block2[i + (j << comb->n)];
        }

        // multiply using zn_poly (requires len1>=len2>=1)
        if(len1>=len2)
            zn_array_mul(out, in1, len1, in2, len2, comb->mod[i]);
        else
            zn_array_mul(out, in2, len2, in1, len1, comb->mod[i]);
        // place result in block_out with proper spacing
        for(int j=0; j<len_out; j++) {
            block_out[i + (j << comb->n)] = out[j];
        }
        
    }

    // Reconstruct output from data in block_out
    for(int i=0; i<len_out; i++) {
        fmpz_t coeff = _fmpz_poly_get_coeff_ptr(output, i);
        fmpz_multi_crt_ui(coeff, block_out + (i << comb->n), comb);
		fmpz_multi_crt_sign(coeff, coeff, comb);
    }
    
    output->length = len_out;
    _fmpz_poly_normalise(output);
    // Free all stuff
    flint_heap_free(block1);
    flint_heap_free(block2);
    flint_heap_free(block_out);
    flint_heap_free(in1);
    flint_heap_free(in2);
    flint_heap_free(out);
}

/*
   Multiply two polynomials using the multimodular technique.
   This function allows aliasing
*/
void _fmpz_poly_mul_modular(fmpz_poly_t output, const fmpz_poly_t poly1, 
									 const fmpz_poly_t poly2, const ulong bits_in)
{

#if FLINT_BITS == 32
    // start from nextprime(2^31)
    // WARNING: there are only about 2^26 primes between 2^31 and 2^32
    // this limits the size of output coefficients to about 65M limbs
    unsigned long p0 = z_nextprime(1UL << 31);
    // primes_per_limb = 32/log2(p0)
    // = 1.0322580642700560413378333946892625831
    double primes_per_limb = 1.0322580642701;
#elif FLINT_BITS == 64
    // start from nextprime(2^64 - 2^56)
    // There should be about 2^50 primes starting there, which should
    // be enough to coefficients of almost 2^50 limbs.
    unsigned long p0 = z_nextprime((unsigned long) -(1L << 56));
    // primes_per_limb = 64/log2(p0)
    // = 1.0000882353338675941356105657797766168
    double primes_per_limb = 1.0000882353339;
#else
#error FLINT_BITS must be either 32 or 64
#endif

    // estimated bound for the size of output coefficients
    unsigned long length = FLINT_MIN(poly1->length, poly2->length);
	 unsigned long output_bits;
	 long bits1, bits2;
	 if (bits_in) output_bits = bits_in;
	 else
	 {
		 bits1= fmpz_poly_max_bits(poly1);
	    bits2 = fmpz_poly_max_bits(poly2);
	    unsigned long log_length = 0;
       while (length > (1L<<log_length)) log_length++;
       output_bits = FLINT_ABS(bits1) + FLINT_ABS(bits2) + log_length + 1;
	 }

    // round up number of primes to a power of two;
    unsigned long numprimes = (output_bits * primes_per_limb)/FLINT_BITS + 1;
	
	 unsigned long* primes = flint_heap_alloc(numprimes);

    unsigned long p = p0;
    for(unsigned long i = 0; i < numprimes; i++) {
        primes[i] = p;
        p = z_nextprime(p);
    }

    // precomputation space
    fmpz_comb_t comb;

    fmpz_comb_init(comb, primes, numprimes);
    __fmpz_poly_mul_modular_comb(output, poly1, poly2, comb);
    fmpz_comb_clear(comb);

    // Free allocated stuff
    flint_heap_free(primes);
}
#endif

/*
   Multiply two polynomials using the classical technique truncating the result to trunc terms.
   Currently doesn't allow aliasing
*/

void _fmpz_poly_mul_classical_trunc(fmpz_poly_t output, const fmpz_poly_t poly1, 
                                          const fmpz_poly_t poly2, const unsigned long trunc)
{
   fmpz_t coeffs_out = output->coeffs;
   unsigned long size_out = output->limbs+1;
   fmpz_t coeffs1, coeffs2;
   unsigned long size1, size2;
   unsigned long len1, len2; 
   unsigned long lenm1;
   
   if (trunc == 0) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      for (unsigned long i = 0; i < trunc; i++)
      {
         coeffs_out[i*size_out] = 0;
      }
      _fmpz_poly_zero(output);
      return;      
   }
      
   fmpz_poly_t input1, input2;
   
   if (output == poly1)
   {
      _fmpz_poly_stack_init(input1, poly1->length, poly1->limbs);
      _fmpz_poly_set(input1, poly1);
      if (output == poly2)
      {
         _fmpz_poly_attach(input2, input1);
      } else _fmpz_poly_attach(input2, poly2);
   } else if (output == poly2)
   {
      _fmpz_poly_stack_init(input2, poly2->length, poly2->limbs);
      _fmpz_poly_set(input2, poly2);
      _fmpz_poly_attach(input1, poly1);
   } else
   {
      _fmpz_poly_attach(input1, poly1);
      _fmpz_poly_attach(input2, poly2);
   }

   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   lenm1 = input1->length-1;
   len1 = input1->length;
   len2 = input2->length;
   
   long i, j;
      
   fmpz_t temp;
            
   // Special case if the length of both inputs is 1
   if ((len1 == 1) && (len2 == 1))
   {
      if ((coeffs1[0] == 0) || (coeffs2[0] == 0))
      {
         coeffs_out[0] = 0;
      } else
      {
         fmpz_mul(coeffs_out, coeffs1, coeffs2);
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
            __fmpz_mul(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
         }
      }
      if (i != len1)
      {
         if ((coeffs1[i*size1] == 0) || (coeffs2[0] == 0))
         {
            coeffs_out[i*size_out] = 0;
         } else
         {
            fmpz_mul(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
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
               __fmpz_mul(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
            }      
         }
         
         /* 
            The above coefficient multiplications overwrite the first limb of the next coefficient
            in each case, using the function __fmpz_mul. The final multiplication 
            cannot do this however.
         */
         if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
         {
            coeffs_out[(i+lenm1)*size_out] = 0;
         } else
         {
            fmpz_mul(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
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
                  fmpz_mul(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } else 
               {
                  fmpz_addmul(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } 
            }
         }
      }
   } 
   
   output->length = FLINT_MIN(len1 + len2 - 1, trunc);
   _fmpz_poly_normalise(output);
   
   if (poly1 == output) _fmpz_poly_stack_clear(input1);
   else if (poly2 == output) _fmpz_poly_stack_clear(input2);
}

/*
   Multiply two polynomials using the classical technique truncating the result 
   so that the first trunc terms are zero.
   Currently doesn't allow aliasing
*/

void _fmpz_poly_mul_classical_trunc_left(fmpz_poly_t output, const fmpz_poly_t poly1, 
                                          const fmpz_poly_t poly2, const unsigned long trunc)
{
   fmpz_t coeffs_out = output->coeffs;
   unsigned long size_out = output->limbs+1;
   fmpz_t coeffs1, coeffs2;
   unsigned long size1, size2;
   unsigned long len1, len2; 
   unsigned long lenm1;
   
   if ((poly1->length == 0) || (poly2->length == 0) || (trunc >= poly1->length + poly2->length - 1)) 
   {
      for (long i = 0; i < (long) (poly1->length + poly2->length - 1); i++)
      {
         coeffs_out[i*size_out] = 0L;
      }
      _fmpz_poly_zero(output);
      return;      
   }
      
   fmpz_poly_t input1, input2;
   
   if (output == poly1)
   {
      _fmpz_poly_stack_init(input1, poly1->length, poly1->limbs);
      _fmpz_poly_set(input1, poly1);
      if (output == poly2)
      {
         _fmpz_poly_attach(input2, input1);
      } else _fmpz_poly_attach(input2, poly2);
   } else if (output == poly2)
   {
      _fmpz_poly_stack_init(input2, poly2->length, poly2->limbs);
      _fmpz_poly_set(input2, poly2);
      _fmpz_poly_attach(input1, poly1);
   } else
   {
      _fmpz_poly_attach(input1, poly1);
      _fmpz_poly_attach(input2, poly2);
   }

   coeffs1 = input1->coeffs;
   coeffs2 = input2->coeffs;
   size1 = input1->limbs+1;
   size2 = input2->limbs+1;
   lenm1 = input1->length-1;
   len1 = input1->length;
   len2 = input2->length;
   
   long i, j;
      
   fmpz_t temp;
            
   // Special case if the length of both inputs is 1
   if ((len1 == 1) && (len2 == 1))
   {
      if ((coeffs1[0] == 0) || (coeffs2[0] == 0))
      {
         coeffs_out[0] = 0;
      } else
      {
         fmpz_mul(coeffs_out, coeffs1, coeffs2);
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
            __fmpz_mul(coeffs_out+i*size_out, coeffs1+i*size1, coeffs2);
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
               __fmpz_mul(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
            }
         }      
      }
      if (len2 == 1) i = 0;   
      /* 
         The above coefficient multiplications overwrite the first limb of the next coefficient
         in each case, using the function __fmpz_mul. The final multiplication 
         cannot do this however.
      */
      if ((coeffs1[lenm1*size1] == 0) || (coeffs2[i*size2] == 0))
      {
         coeffs_out[(i+lenm1)*size_out] = 0;
      } else
      {
         fmpz_mul(coeffs_out+(i+lenm1)*size_out, coeffs1+lenm1*size1, coeffs2+i*size2);
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
                  fmpz_mul(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
               } else 
               {
                  fmpz_addmul(coeffs_out+(i+j)*size_out, coeffs1+i*size1, coeffs2+j*size2);
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
   if (trunc >= output->length) _fmpz_poly_normalise(output);
   
   if (poly1 == output) _fmpz_poly_stack_clear(input1);
   else if (poly2 == output) _fmpz_poly_stack_clear(input2);
}

void __fmpz_poly_karamul_recursive(fmpz_poly_t res, const fmpz_poly_t a, const fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, const unsigned long crossover)
{
   fmpz_poly_t temp;
   
   if ((crossover < 4) && (a->length == 2 && b->length == 2)) {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      __fmpz_mul(res->coeffs, a->coeffs, b->coeffs); 
      fmpz_add(scratchb->coeffs, a->coeffs, a->coeffs+asize);
      fmpz_mul(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
      fmpz_add(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
      fmpz_mul(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
      fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
      fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
      res->length = a->length + b->length - 1;
      
      return;
   }
   
   if ((a->length+b->length <= crossover) ||  (a->length <= 1) || (b->length <= 1) ||  ((a->length == 2) || (b->length == 2)))
   {
      _fmpz_poly_mul_classical(res, a, b);
      
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
      if (asum->length > bsum->length) __fmpz_poly_karamul_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover);
      else __fmpz_poly_karamul_recursive(prodsum, bsum, asum, scratch2, scratch3, crossover);
      for (long i = prodsum->length; i < (a1->length<<1)-1; i++)
          prodsum->coeffs[i*(prodsum->limbs+1)] = 0L;
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
         res->coeffs[i*(res->limbs+1)] = 0L;
      
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
      _fmpz_poly_add(temp1, temp1, temp);
  
      res->length = a->length + b->length - 1;
   } 
}

void _fmpz_poly_mul_karatsuba(fmpz_poly_t output, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   unsigned long limbs = output->limbs;
   unsigned long log_length = 0;
   unsigned long crossover;
   
   fmpz_poly_t input1, input2;
   
   if (output == poly1)
   {
      _fmpz_poly_stack_init(input1, poly1->length, poly1->limbs);
      _fmpz_poly_set(input1, poly1);
      if (output == poly2)
      {
         _fmpz_poly_attach(input2, input1);
      } else _fmpz_poly_attach(input2, poly2);
   } else if (output == poly2)
   {
      _fmpz_poly_stack_init(input2, poly2->length, poly2->limbs);
      _fmpz_poly_set(input2, poly2);
      _fmpz_poly_attach(input1, poly1);
   } else
   {
      _fmpz_poly_attach(input1, poly1);
      _fmpz_poly_attach(input2, poly2);
   }

   fmpz_poly_t scratch, scratchb, temp;
   scratch->coeffs = (fmpz_t) flint_stack_alloc(5*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs + 1;
   scratchb->limbs = FLINT_MAX(input1->limbs, input2->limbs) + 1;
   scratchb->coeffs = (fmpz_t) flint_stack_alloc(5*FLINT_MAX(input1->length, input2->length)*(scratchb->limbs+1));
   
   if (_fmpz_poly_max_limbs(input1) + _fmpz_poly_max_limbs(input2) >= 19) crossover = 0;
   else crossover = 19 - _fmpz_poly_max_limbs(input1) - _fmpz_poly_max_limbs(input2);
   
   if (input1->length >= input2->length)
       __fmpz_poly_karamul_recursive(output, input1, input2, scratch, scratchb, crossover);
   else
       __fmpz_poly_karamul_recursive(output, input2, input1, scratch, scratchb, crossover);
   
   flint_stack_release(); 
   flint_stack_release();
   
   if (poly1 == output) _fmpz_poly_stack_clear(input1);
   else if (poly2 == output) _fmpz_poly_stack_clear(input2);
}

/*
   The recursive part of fmpz_poly_mul_karatsuba_trunc. It multiplies polynomials a and b
   storing the first trunc coefficients of the result in res.
   To prevent valgrind complaining, leading coefficients, up to length trunc are
   zeroed and res is normalised.
   If a->length + b->length <= crossover, the algorithm drops into the truncated classical 
   technique. 
   Optimised versions of the algorithm are implemented for less than or equal to
   2x2 multiplications, should crossover be low enough to allow these.
   There are two lots of scratch space; scratch has space for coefficients with size
   equal to those of the final output coefficients of the entire multiplication 
   whilst scratchb has space for coefficients one more limb in size than either 
   of the input coefficients (to allow for addition of such coefficients).
   trunc may be bigger than a full product of polynomials a and b.
   Assumes trunc is at least 1 and a->length >= b->length
   a and b need not be normalised
   warning: the code to zero up to trunc should remain before any branches using 
   scratch space
   todo: have this function only deal with normalised inputs for speed
   todo: currently if trunc_res is much bigger than 2*a1->length-1 the recursive 
   call to karatrunc_recursive in the main branch zeroes all the higher coefficients, 
   some of which get overwritten, a similar thing happens in the degenerate branch
   instead the call to karatrunc_recursive should be truncated earlier and only those
   coefficients which won't get overwritten should be zeroed
*/

void __fmpz_poly_karatrunc_recursive(fmpz_poly_t res, const fmpz_poly_t a, const fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, const unsigned long crossover, const unsigned long trunc)
{
   // Allow for trunc to be bigger than the length of a full product of a and b
   unsigned long trunc_res = FLINT_MIN(trunc, a->length + b->length - 1);
   for (unsigned long i = trunc_res; i < trunc; i++)
   {
	  res->coeffs[i*(res->limbs+1)] = 0L;
   }
         
   // Drop to truncated classical algorithm if either of the polynomials has length <= 1
   if ((a->length <= 1) || (b->length <= 1)) 
   {
      // todo: a and b are not necessarily normalised; they should be normalised
      // for speed
	  _fmpz_poly_mul_classical_trunc(res, a, b, trunc_res);
      
	  // fixme: this zeroing is unnecessary if classical_trunc zeroes up to trunc_res
	  for (unsigned long i = res->length; i < trunc_res; i++)
	  {
	     res->coeffs[i*(res->limbs+1)] = 0L;
	  }
      
      return;
   }
   
   // Special code for trunc_res == 1 (which the main branch cannot handle)
   // recall trunc_res == 0 is not allowed and a and b have length at least 2 by now
   if (trunc_res == 1)
   {
      fmpz_mul(res->coeffs, a->coeffs, b->coeffs); 
	  res->length = 1;

	  // the bottom coefficient of either a or b may be zero, so normalise
	  _fmpz_poly_normalise(res);
	  
	  return;
   }

   // todo: check to see if (crossover < 3) would be better here as the next 
   // case would seem to imply
   // Special code for 2 x 2 in the case where the crossover is <= 3
   if ((a->length == 2 && b->length == 2) && (crossover < 4)) 
   {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      // trunc_res is at least 2 so multiply lowest coefficients
	  // we can use fast multiplication since output has at least one
	  // more coefficient
	  __fmpz_mul(res->coeffs, a->coeffs, b->coeffs); 
         
      // scratchb[0] = a[0] + a[1]
	  // scratchb[1] = b[0] + b[1]
	  fmpz_add(scratchb->coeffs, a->coeffs, a->coeffs+asize);
      fmpz_add(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
         
      // if trunc_res > 2 multiply out top coefficient and store in result
	  // otherwise store it in scratch[0]
	  if (trunc_res > 2) fmpz_mul(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
      else fmpz_mul(scratch->coeffs, a->coeffs+asize, b->coeffs+bsize); 
         
      // trunc_res is at least 2 so set res[1] = scratchb[0]*scratchb[1]
	  fmpz_mul(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
      
	  // subtract a[0]*b[0] from res[1]
      fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
      
	  // subtract a[1]*b[1] from res[1]
	  if (trunc_res > 2) fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
      else fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, scratch->coeffs);
            
      res->length = trunc_res;
      _fmpz_poly_normalise(res);
	  
      return;
   }

   // If we have reached the crossover or if the 2 x 2 case wasn't dealt with
   if ((a->length + b->length <= crossover) || ((a->length == 2) && (b->length == 2)))
   {
      // todo: a and b are not necessarily normalised; they should be normalised
      // for speed
	  _fmpz_poly_mul_classical_trunc(res, a, b, trunc_res);
      
	  // fixme: this zeroing is unnecessary if classical_trunc zeroes up to trunc_res
	  for (unsigned long i = res->length; i < trunc_res; i++)
	  {
	     res->coeffs[i*(res->limbs+1)] = 0L;
	  }
      
      return;
   }   
        
   fmpz_poly_t a1, a2, b1, b2, temp, temp2;
      
   unsigned long log_length = 0; 
   unsigned long old_length;
   
   // we only consider the first trunc_res coefficients of a and b if this is 
   // less than their current length
   // note sa and sb are both at least 2
   unsigned long sa = FLINT_MIN(a->length, trunc_res);
   unsigned long sb = FLINT_MIN(b->length, trunc_res);
     
   // split a into a1 and a2 with a1->length >= a2->length >= 1
   a1->length = (sa+1)/2;
   a2->length = sa-a1->length;
   a1->coeffs = a->coeffs;
   a2->coeffs = a->coeffs+a1->length*(a->limbs+1);
   a1->limbs = a->limbs;
   a2->limbs = a->limbs;
   
   if (a1->length < sb) //ordinary branch
   {
      /*
         (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
      */
      
      // sb is at least 2 and sa >= sb so we split b into b1 and b2 with b1 the 
      // same length as a1 and b2 equal to the rest
      // as sb > a1->length, b2 has length at least 1
      // note, none of a1, a2, b1, b2 is necessarily normalised
	  b1->length = a1->length;
      b2->length = sb - b1->length;
      b1->coeffs = b->coeffs;
      b2->coeffs = b->coeffs + b1->length*(b->limbs+1);
      b1->limbs = b->limbs;
      b2->limbs = b->limbs;
      
      fmpz_poly_t asum, bsum, prodsum, scratch2, scratch3;
     
      asum->length = a1->length;
      asum->coeffs = scratchb->coeffs;
      asum->limbs = scratchb->limbs;
      bsum->length = a1->length;
      bsum->coeffs = scratchb->coeffs + a1->length*(scratchb->limbs+1);
      bsum->limbs = scratchb->limbs;
      prodsum->length = 2*a1->length-1;
      prodsum->coeffs = scratch->coeffs+(a2->length+b2->length-1)*(scratch->limbs+1);
      prodsum->limbs = scratch->limbs;
      
      // res_lo = a1*b1
	  // note everything gets zeroed or written to up to trunc_res at this point
      __fmpz_poly_karatrunc_recursive(res, a1, b1, scratch, scratchb, crossover, trunc_res);
      
      // res_hi = a2*b2
	  // write to scratch[0] and following, with length no greater than
	  // a->length+b->length-1
	  // extras zeroes will only be written at the start of the algorithm and so 
	  // will not overwrite the scratch space for the main branches of the algorithm
	  // which use the scratch space
	  // note trunc_res > a1->length
      temp->coeffs = scratch->coeffs;
      temp->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs + (a2->length+b2->length-1)*(scratch->limbs+1);
      scratch2->limbs = scratch->limbs;
      __fmpz_poly_karatrunc_recursive(temp, a2, b2, scratch2, scratchb, crossover, trunc_res - a1->length);
      
      // temp2 will be used as an alias for parts of res
	  temp2->limbs = res->limbs;
      
	  // write any of the coefficients just computed into res according to whether 
	  // they occur before trunc_res and zero any unwritten coefficients
	  if (trunc_res > 2*a1->length)
      {
         old_length = temp->length; // save temp->length
         temp->length = FLINT_MIN(old_length, trunc_res-2*a1->length);
         temp2->coeffs = res->coeffs+2*a1->length*(res->limbs+1);
         _fmpz_poly_set(temp2, temp);

		 temp->length = old_length; // restore temp->length
      }
      
      // asum = a1+a2
      _fmpz_poly_add(asum, a1, a2);
      // bsum = b1+b2
      _fmpz_poly_add(bsum, b1, b2);

      // asum and bsum are currently in scratchb, both have length a1->length
	  // prodsum of length 2*a1->length - 1 and temp of length 
	  // a2->length + b2->length - 1 are in scratch
	  scratch3->coeffs = scratchb->coeffs+2*a1->length*(scratchb->limbs+1);
      scratch3->limbs = scratchb->limbs;   
      scratch2->coeffs = scratch->coeffs + (sa+sb-2)*(scratch->limbs+1);
      
	  // prodsum = asum*bsum
      if (asum->length > bsum->length) __fmpz_poly_karatrunc_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover, trunc_res - a1->length);
      else __fmpz_poly_karatrunc_recursive(prodsum, bsum, asum, scratch2, scratch3, crossover, trunc_res - a1->length);
      
	  // prodsum is normalised, so zero remaining coefficients as these will have been
	  // used as scratch space
	  for (long i = prodsum->length; i < trunc_res - a1->length; i++)
          prodsum->coeffs[i*(prodsum->limbs+1)] = 0L;
      
      // prodsum = prodsum - res_lo
      temp2->coeffs = res->coeffs;
      temp2->length = FLINT_MIN(trunc_res - a1->length, 2*a1->length - 1);
      _fmpz_poly_sub(prodsum, prodsum, temp2);
       
      // prodsum = prodsum - res_hi
      temp->length = FLINT_MIN(trunc_res - a1->length, temp->length); 
      _fmpz_poly_sub(prodsum, prodsum, temp);
      
	  prodsum->length = FLINT_MIN(prodsum->length, trunc_res - a1->length);
      
      // res_mid += prodsum
      temp2->coeffs = res->coeffs + a1->length*(res->limbs+1);
      temp2->length = prodsum->length;
      _fmpz_poly_add(temp2, temp2, prodsum);
      
      res->length = FLINT_MIN(a->length + b->length - 1, trunc_res);
      _fmpz_poly_normalise(res);

   } else // a1->length <= sb
   {
      fmpz_poly_t scratch2, temp1; 

      while ((1<<log_length)<a1->length) log_length++;
      if ((1<<log_length) < sa) a1->length = (1<<log_length);
      a2->length = sa - a1->length;
      a1->coeffs = a->coeffs;
      a2->coeffs = a->coeffs+a1->length*(a->limbs+1);
      
      temp->coeffs = b->coeffs;
      temp->limbs = b->limbs;
      temp->length = sb;
      
      // res_lo = a1*b
	  // note everything is zeroed or written to up to trunc_res at this point
      if (sb <= a1->length) __fmpz_poly_karatrunc_recursive(res, a1, temp, scratch, scratchb, crossover, trunc_res);
      else __fmpz_poly_karatrunc_recursive(res, temp, a1, scratch, scratchb, crossover, trunc_res);
      
      //temp2 = a2*b
      temp2->coeffs = scratch->coeffs;
      temp2->length = FLINT_MIN(a2->length + sb - 1, trunc_res - a1->length);
      temp2->limbs = scratch->limbs;
      scratch2->coeffs = scratch->coeffs + temp2->length*(scratch->limbs+1);
      scratch2->limbs = scratch->limbs;
      if (sb <= a2->length) __fmpz_poly_karatrunc_recursive(temp2, a2, temp, scratch2, scratchb, crossover, trunc_res - a1->length);
      else __fmpz_poly_karatrunc_recursive(temp2, temp, a2, scratch2, scratchb, crossover, trunc_res - a1->length);
      
      // res_mid += temp2
      temp1->coeffs = res->coeffs+a1->length*(res->limbs+1);
      temp1->length = temp2->length;
      temp1->limbs = res->limbs;
      _fmpz_poly_add(temp1, temp1, temp2);
      
      res->length = FLINT_MIN(sa + sb - 1, trunc_res);
	  _fmpz_poly_normalise(res);
   } 
}

void _fmpz_poly_mul_karatsuba_trunc(fmpz_poly_t output, const fmpz_poly_t poly1, const fmpz_poly_t poly2, const unsigned long trunc)
{
   if ((poly1->length == 0) || (poly2->length == 0) || (trunc == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   unsigned long limbs = output->limbs;
   unsigned long log_length = 0;
   unsigned long crossover;
   
    fmpz_poly_t input1, input2;
   
   if (output == poly1)
   {
      _fmpz_poly_stack_init(input1, poly1->length, poly1->limbs);
      _fmpz_poly_set(input1, poly1);
      if (output == poly2)
      {
         _fmpz_poly_attach(input2, input1);
      } else _fmpz_poly_attach(input2, poly2);
   } else if (output == poly2)
   {
      _fmpz_poly_stack_init(input2, poly2->length, poly2->limbs);
      _fmpz_poly_set(input2, poly2);
      _fmpz_poly_attach(input1, poly1);
   } else
   {
      _fmpz_poly_attach(input1, poly1);
      _fmpz_poly_attach(input2, poly2);
   }

   fmpz_poly_t scratch, scratchb, temp;
   scratch->coeffs = (fmpz_t) flint_stack_alloc(6*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs;
   scratchb->limbs = FLINT_MAX(input1->limbs, input2->limbs)+1;
   scratchb->coeffs = (fmpz_t) flint_stack_alloc(6*FLINT_MAX(input1->length,input2->length)*(scratchb->limbs+1));
   
   if (_fmpz_poly_max_limbs(input1) + _fmpz_poly_max_limbs(input2) >= 19) crossover = 0;
   else crossover = 19 - _fmpz_poly_max_limbs(input1) - _fmpz_poly_max_limbs(input2);
   
   unsigned long trunc_res = FLINT_MIN(trunc, input1->length + input2->length - 1);

   if (input1->length >= input2->length)
       __fmpz_poly_karatrunc_recursive(output, input1, input2, scratch, scratchb, crossover, trunc_res);
   else
       __fmpz_poly_karatrunc_recursive(output, input2, input1, scratch, scratchb, crossover, trunc_res);
   
   flint_stack_release(); flint_stack_release();
   
   if (poly1 == output) _fmpz_poly_stack_clear(input1);
   else if (poly2 == output) _fmpz_poly_stack_clear(input2);
}

void __fmpz_poly_karatrunc_left_recursive(fmpz_poly_t res, const fmpz_poly_t a, const fmpz_poly_t b, fmpz_poly_t scratch, fmpz_poly_t scratchb, const unsigned long crossover, const unsigned long trunc)
{
   fmpz_poly_t temp, temp2;
   
   long non_zero = a->length + b->length - trunc - 1;
   if (non_zero <= 0)
   {
      for (long i = 0; i < (long) (a->length + b->length - 1); i++)
      {
         res->coeffs[i*(res->limbs+1)] = 0;
      }
      res->length = 0;
      return;
   }
   
   if ((a->length <= 1) || (b->length <= 1) || (non_zero == 1)) 
   {
      _fmpz_poly_mul_classical_trunc_left(res, a, b, trunc);
         
      return;
   }
   
   if ((a->length == 2 && b->length == 2) && (crossover < 4) && (!trunc)) 
   {
      const unsigned long asize = a->limbs+1;
      const unsigned long bsize = b->limbs+1;
      const unsigned long rsize = res->limbs+1;
      const unsigned long ssize = scratchb->limbs+1;
      
      __fmpz_mul(res->coeffs, a->coeffs, b->coeffs); 
         
      fmpz_add(scratchb->coeffs, a->coeffs, a->coeffs+asize);
      fmpz_add(scratchb->coeffs+ssize, b->coeffs, b->coeffs+bsize);
         
      fmpz_mul(res->coeffs+2*rsize, a->coeffs+asize, b->coeffs+bsize); 
         
      fmpz_mul(res->coeffs+rsize, scratchb->coeffs, scratchb->coeffs+ssize); 
         
      fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, res->coeffs);
      fmpz_sub(res->coeffs+rsize, res->coeffs+rsize, res->coeffs+2*rsize);
      
            
      res->length = a->length + b->length - 1;
      
      return;
   }
   
   if ((a->length+b->length <= crossover) || ((a->length == 2) && (b->length == 2)))
   {
      _fmpz_poly_mul_classical_trunc_left(res, a, b, trunc);
      
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
         if (trunc > a1->length) 
         {
            if (asum->length > bsum->length) __fmpz_poly_karatrunc_left_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover, trunc - a1->length);
            else __fmpz_poly_karatrunc_left_recursive(prodsum, bsum, asum, scratch2, scratch3, crossover, trunc - a1->length);
         } else 
         {
            if (asum->length > bsum->length) __fmpz_poly_karatrunc_left_recursive(prodsum, asum, bsum, scratch2, scratch3, crossover, 0);
            else __fmpz_poly_karatrunc_left_recursive(prodsum, bsum, asum, scratch2, scratch3, crossover, 0);
         }
         for (long i = prodsum->length; i < (a1->length<<1)-1; i++)
            prodsum->coeffs[i*(prodsum->limbs+1)] = 0L;
               
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
      for (unsigned long i = 0; i < a->length + b->length - 1; i++) // <- Bug (should only need a1->length + b->length to a->length+b->length-1
         res->coeffs[i*(res->limbs+1)] = 0L;
      
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
      _fmpz_poly_add(temp1, temp1, temp); 
  
      res->length = a->length + b->length - 1;
   } 
   
   for (unsigned long i = 0; i < trunc; i++)
   {
      res->coeffs[i*(res->limbs+1)] = 0L;
   }     
}

void _fmpz_poly_mul_karatsuba_trunc_left(fmpz_poly_t output, const fmpz_poly_t poly1, const fmpz_poly_t poly2, const unsigned long trunc)
{
   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }

   unsigned long limbs = output->limbs;
   unsigned long log_length = 0;
   unsigned long crossover;
   
   fmpz_poly_t input1, input2;
   
   if (output == poly1)
   {
      _fmpz_poly_stack_init(input1, poly1->length, poly1->limbs);
      _fmpz_poly_set(input1, poly1);
      if (output == poly2)
      {
         _fmpz_poly_attach(input2, input1);
      } else _fmpz_poly_attach(input2, poly2);
   } else if (output == poly2)
   {
      _fmpz_poly_stack_init(input2, poly2->length, poly2->limbs);
      _fmpz_poly_set(input2, poly2);
      _fmpz_poly_attach(input1, poly1);
   } else
   {
      _fmpz_poly_attach(input1, poly1);
      _fmpz_poly_attach(input2, poly2);
   }

   fmpz_poly_t scratch, scratchb, temp;
   scratch->coeffs = (fmpz_t) flint_stack_alloc(5*FLINT_MAX(input1->length,input2->length)*(limbs+1));
   scratch->limbs = limbs + 1;
   scratchb->limbs = FLINT_MAX(input1->limbs, input2->limbs) + 1;
   scratchb->coeffs = (fmpz_t) flint_stack_alloc(5*FLINT_MAX(input1->length, input2->length)*(scratchb->limbs+1));
   
   if (_fmpz_poly_max_limbs(input1) + _fmpz_poly_max_limbs(input2) >= 19) crossover = 0;
   else crossover = 19 - _fmpz_poly_max_limbs(input1) - _fmpz_poly_max_limbs(input2);
   
   if (input1->length >= input2->length)
       __fmpz_poly_karatrunc_left_recursive(output, input1, input2, scratch, scratchb, crossover, trunc);
   else
       __fmpz_poly_karatrunc_left_recursive(output, input2, input1, scratch, scratchb, crossover, trunc);
   
   flint_stack_release(); 
   flint_stack_release();
   if (trunc >= input1->length+input2->length-1) _fmpz_poly_normalise(output);
   
   if (poly1 == output) _fmpz_poly_stack_clear(input1);
   else if (poly2 == output) _fmpz_poly_stack_clear(input2);
}

void _fmpz_poly_mul_KS(fmpz_poly_t output, const fmpz_poly_t in1, const fmpz_poly_t in2, const long bits_in)
{
   long sign1 = 1L;
   long sign2 = 1L;
   
   unsigned long length1 = in1->length;
   unsigned long length2 = in2->length;
   
   unsigned long final_length = length1 + length2 - 1;
   
   while ((length1) && (in1->coeffs[(length1-1)*(in1->limbs+1)] == 0L)) length1--;
   while ((length2) && (in2->coeffs[(length2-1)*(in2->limbs+1)] == 0L)) length2--;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }

   fmpz_poly_t input1, input2;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      _fmpz_poly_attach(input1, in2);
      _fmpz_poly_attach(input2, in1);
   } else
   {
      _fmpz_poly_attach(input1, in1);
      _fmpz_poly_attach(input2, in2);
   }
   
   long bits1, bits2;
   int bitpack = 0;
   
   unsigned long bits;
	unsigned long sign;
	if (bits_in) 
	{
		bits = FLINT_ABS(bits_in);
		sign = (bits_in < 0L);
	}
	else
	{
		bits1 = _fmpz_poly_max_bits(input1);
      bits2 = (in1 == in2) ? bits1 : _fmpz_poly_max_bits(input2);
      
      sign = ((bits1 < 0L) || (bits2 < 0L));
      unsigned long length = length2;
      unsigned long log_length = 0L;
      while ((1<<log_length) < length) log_length++;
      bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
	}
	
	unsigned long limbs = (bits-1)/FLINT_BITS + 1;
   
   if ((bits < FLINT_BITS) && (input1->limbs == 1) && (input2->limbs == 1) && (output->limbs == 1)) bitpack = 1;
   
   unsigned long bytes = ((bits-1)>>3)+1;
   
   if ((long) input1->coeffs[(length1-1)*(input1->limbs+1)] < 0L) sign1 = -1L;
   
   if (in1 != in2)
   {
      if ((long) input2->coeffs[(length2-1)*(input2->limbs+1)] < 0L) sign2 = -1L;
   } else sign2 = sign1;
   
   mp_limb_t * array1, * array2, * array3;
	ulong p1n, p2n;
   if (bitpack)
   {
      p1n = (bits*length1 - 1)/FLINT_BITS + 1;
		array1 = flint_stack_alloc(p1n + 1); //extra limb required
      if (in1 != in2)
		{
			p2n = (bits*length2 - 1)/FLINT_BITS + 1;
		   array2 = flint_stack_alloc(p2n + 1); //extra limb required
		}
      
      if (sign) bits = -1L*bits;
      if (in1 != in2)
         fmpz_poly_bit_pack(array2, input2, length2, bits, sign2);
      fmpz_poly_bit_pack(array1, input1, length1, bits, sign1);

      bits=ABS(bits);
   } else
   {
      p1n = ((bytes*length1-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
		array1 = flint_stack_alloc(p1n + 1); // extra limb required
      if (in1 != in2)
		{
         p2n = ((bytes*length2-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
			array2 = flint_stack_alloc(p2n + 1); // extra limb required
		}

      fmpz_poly_byte_pack(array1, input1, length1, bytes, sign1);
      if (in1 != in2)
         fmpz_poly_byte_pack(array2, input2, length2, bytes, sign2);
   }
   
   if (in1 == in2)
   {
      array2 = array1;
      p2n = p1n;
   }
   
   array3 = flint_stack_alloc(p1n + p2n + 1);
           
   mp_limb_t msl = F_mpn_mul(array3, array1, p1n, array2, p2n);
   
   array3[p1n+p2n-1] = msl;
   array3[p1n + p2n] = 0L;

   output->length = length1 + length2 - 1;
  
   for (unsigned long i = 0; i < output->length; i++)
   {
      output->coeffs[i*(output->limbs+1)] = 0L;
   }
      
   if (bitpack)
   {
      if (sign) fmpz_poly_bit_unpack(output, array3, length1 + length2 - 1, bits);  
      else fmpz_poly_bit_unpack_unsigned(output, array3, length1 + length2 - 1, bits);  
   } else
   {
      if (sign) fmpz_poly_byte_unpack(output, array3, length1 + length2 - 1, bytes);        
      else fmpz_poly_byte_unpack_unsigned(output, array3, length1 + length2 - 1, bytes);  
   }
   
   flint_stack_release(); // array3
   if (in1 != in2)
      flint_stack_release(); // array2
   flint_stack_release(); // array1
     
   if ((long) (sign1 ^ sign2) < 0L) _fmpz_poly_neg(output, output);
   
   output->length = length1 + length2 - 1;
}

void _fmpz_poly_mul_KS_trunc(fmpz_poly_t output, const fmpz_poly_t in1, 
                                        const fmpz_poly_t in2, const unsigned long trunc, long bits_in)
{
   long sign1 = 1L;
   long sign2 = 1L;
   
   unsigned long length1 = FLINT_MIN(in1->length, trunc);
   unsigned long length2 = FLINT_MIN(in2->length, trunc);
   
   while ((length1) && (in1->coeffs[(length1-1)*(in1->limbs+1)] == 0)) length1--;
   while ((length2) && (in2->coeffs[(length2-1)*(in2->limbs+1)] == 0)) length2--;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   fmpz_poly_t input1, input2;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      _fmpz_poly_attach(input1, in2);
      _fmpz_poly_attach(input2, in1);
   } else
   {
      _fmpz_poly_attach(input1, in1);
      _fmpz_poly_attach(input2, in2);
   }
   
   long bits1, bits2;
   int bitpack = 0;
   
   unsigned long sign;
   unsigned long length;
   unsigned log_length;

   if (!bits_in)
   {
	  bits1 = _fmpz_poly_max_bits(input1);
      bits2 = (in1 == in2) ? bits1 : _fmpz_poly_max_bits(input2);
      sign = ((bits1 < 0) || (bits2 < 0));
   } else sign = ((long) bits_in < 0L);

   length = length2;
   log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits;
   
   if (bits_in) bits = FLINT_ABS(bits_in);
   else bits = ABS(bits1) + ABS(bits2) + log_length + sign; 

   unsigned long limbs = (bits-1)/FLINT_BITS + 1;
   
   if ((bits < FLINT_BITS) && (input1->limbs == 1) && (input2->limbs == 1) && (output->limbs == 1)) bitpack = 1;
   
   unsigned long bytes = ((bits-1)>>3)+1;
   
   if ((long) input1->coeffs[(length1-1)*(input1->limbs+1)] < 0) sign1 = -1L;
   
   if (in1 != in2)
   {
      if ((long) input2->coeffs[(length2-1)*(input2->limbs+1)] < 0) sign2 = -1L;
   } else sign2 = sign1;
   
   mp_limb_t * array1, * array2, * array3;
	ulong p1n, p2n;
   if (bitpack)
   {
      p1n = (bits*length1-1)/FLINT_BITS+1;
		array1 = flint_stack_alloc(p1n + 1); // extra limb required
      if (in1 != in2)
		{
			p2n = (bits*length2-1)/FLINT_BITS+1;
		   array2 = flint_stack_alloc(p2n + 1); // extra limb required
		}

      if (sign) bits = -bits;
      if (in1 != in2)
         fmpz_poly_bit_pack(array2, input2, length2, bits, sign2);
      fmpz_poly_bit_pack(array1, input1, length1, bits, sign1);

      bits = ABS(bits);
   } else
   {
      p1n = ((bytes*length1-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
		array1 = flint_stack_alloc(p1n + 1); // extra limb required
      if (in1 != in2)
		{
			p2n = ((bytes*length2-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
		   array2 = flint_stack_alloc(p2n + 1); // extra limb required
		}
      fmpz_poly_byte_pack(array1, input1, length1, bytes, sign1);
      if (in1 != in2)
         fmpz_poly_byte_pack(array2, input2, length2, bytes, sign2);
   }
   
   if (in1 == in2)
   {
      array2 = array1;
      p2n = p1n;
   }
   
   array3 = flint_stack_alloc(p1n + p2n + 1);
           
   output->length = FLINT_MIN(length1 + length2 - 1, trunc);

   mp_limb_t msl;
   
   if (bitpack)
   {
      msl = F_mpn_mul_trunc(array3, array1, p1n, array2, p2n, (output->length*bits-1)/FLINT_BITS+1);
   } else
   {
      msl = F_mpn_mul_trunc(array3, array1, p1n, array2, p2n, ((output->length*bytes-1)>>FLINT_LG_BYTES_PER_LIMB) + 1);
   }
   
   array3[p1n + p2n - 1] = msl;
   array3[p1n + p2n] = 0;
      
   for (unsigned long i = 0; i < trunc; i++)
      output->coeffs[i*(output->limbs+1)] = 0;
      
   if (bitpack)
   {
      if (sign) fmpz_poly_bit_unpack(output, array3, output->length, bits);  
      else fmpz_poly_bit_unpack_unsigned(output, array3, output->length, bits);  
   } else
   {
      if (sign) fmpz_poly_byte_unpack(output, array3, output->length, bytes);        
      else fmpz_poly_byte_unpack_unsigned(output, array3, output->length, bytes);  
   }
   
   flint_stack_release(); // array3
   if (in1 != in2)
      flint_stack_release(); // array2
   flint_stack_release(); // array1
     
   if ((long) (sign1 ^ sign2) < 0) _fmpz_poly_neg(output, output);
   _fmpz_poly_normalise(output);
}

void _fmpz_poly_mul_SS(fmpz_poly_t output, const fmpz_poly_t in1, const fmpz_poly_t in2)
{
   unsigned long length1 = in1->length;
   while ((length1) && (in1->coeffs[(length1-1)*(in1->limbs+1)] == 0)) length1--;
   
   unsigned long length2;
   if (in1 != in2)
   {
       length2= in2->length;
      while ((length2) && (in2->coeffs[(length2-1)*(in2->limbs+1)] == 0)) length2--;
   } else length2 = length1;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   fmpz_poly_t input1, input2;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      _fmpz_poly_attach(input1, in2);
      _fmpz_poly_attach(input2, in1);
   } else
   {
      _fmpz_poly_attach(input1, in1);
      _fmpz_poly_attach(input2, in2);
   }
   
   unsigned long size1 = input1->limbs;
   unsigned long size2 = input2->limbs;
   
   unsigned long log_length = 0;
   while ((1<<log_length) < length1) log_length++;
   unsigned long log_length2 = 0;
   
   if (in1 != in2) while ((1<<log_length2) < length2) log_length2++;
   else (log_length2 = log_length);
   
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
   if (in1 != in2) ZmodF_poly_stack_init(poly2, log_length + 1, n, 1);
   ZmodF_poly_stack_init(res, log_length + 1, n, 1);
   
   bits1 = fmpz_poly_to_ZmodF_poly(poly1, input1, length1);
   if (in1 != in2) bits2 = fmpz_poly_to_ZmodF_poly(poly2, input2, length2);
   else bits2 = bits1;
   
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
   if (in1 != in2) ZmodF_poly_decrease_n(poly2, n);
   ZmodF_poly_decrease_n(res, n);
                    
   if (in1 != in2) ZmodF_poly_convolution(res, poly1, poly2);
   else ZmodF_poly_convolution(res, poly1, poly1);
   ZmodF_poly_normalise(res);
          
   output->length = length1 + length2 - 1;
   
   ZmodF_poly_to_fmpz_poly(output, res, sign);
   
   ZmodF_poly_stack_clear(res);
   if (in1 != in2) ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
}

void _fmpz_poly_mul_SS_trunc(fmpz_poly_t output, const fmpz_poly_t in1, 
                                        const fmpz_poly_t in2, const unsigned long trunc)
{
   unsigned long length1 = FLINT_MIN(in1->length, trunc);
   while ((length1) && (in1->coeffs[(length1-1)*(in1->limbs+1)] == 0)) length1--;
   
   unsigned long length2;
   if (in1 != in2)
   {
      length2 = FLINT_MIN(in2->length, trunc);
      while ((length2) && (in2->coeffs[(length2-1)*(in2->limbs+1)] == 0)) length2--;
   } else length2 = length1;
   
   if ((length1 == 0) || (length2 == 0)) 
   {
      _fmpz_poly_zero(output);
      return;
   }
   
   fmpz_poly_t input1, input2;
   
   if (length2 > length1) 
   {
      unsigned long temp = length1;
      length1 = length2;
      length2 = temp;
      _fmpz_poly_attach(input1, in2);
      _fmpz_poly_attach(input2, in1);
   } else
   {
      _fmpz_poly_attach(input1, in1);
      _fmpz_poly_attach(input2, in2);
   }
   
   unsigned long size1 = input1->limbs;
   unsigned long size2 = input2->limbs;
   
   unsigned long log_length = 0;
   while ((1<<log_length) < length1) log_length++;
   unsigned long log_length2 = 0;
   if (in1 != in2) while ((1<<log_length2) < length2) log_length2++;
   else log_length2 = log_length;
   
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
   if (in1 != in2) ZmodF_poly_stack_init(poly2, log_length + 1, n, 1);
   ZmodF_poly_stack_init(res, log_length + 1, n, 1);
   
   bits1 = fmpz_poly_to_ZmodF_poly(poly1, input1, length1);
   if (in1 != in2) bits2 = fmpz_poly_to_ZmodF_poly(poly2, input2, length2);
   else (bits2 = bits1);
   
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
   if (in1 != in2) ZmodF_poly_decrease_n(poly2, n);
   ZmodF_poly_decrease_n(res, n);
                    
   if (in1 != in2) ZmodF_poly_convolution_range(res, poly1, poly2, 0, trunc);
   else ZmodF_poly_convolution_range(res, poly1, poly1, 0, trunc);

   res->length = FLINT_MIN(res->length, trunc);
   ZmodF_poly_normalise(res);
          
   output->length = FLINT_MIN(length1 + length2 - 1, trunc);
   
   ZmodF_poly_to_fmpz_poly(output, res, sign);
   
   ZmodF_poly_stack_clear(res);
   if (in1 != in2) ZmodF_poly_stack_clear(poly2);
   ZmodF_poly_stack_clear(poly1);
   _fmpz_poly_normalise(output);
}

void _fmpz_poly_mul(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
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
      _fmpz_poly_mul_KS(output, input1, input2, 0);
      return;
   }
   
   if (input1->length + input2->length <= 32) 
   {
      _fmpz_poly_mul_karatsuba(output, input1, input2);
      return;
   }
   
   unsigned long bits1 = _fmpz_poly_max_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _fmpz_poly_max_bits(input2);
   bits1 = ABS(bits1);
   bits2 = ABS(bits2);
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _fmpz_poly_mul_SS(output, input1, input2);
      return;
   } 
   
   _fmpz_poly_mul_KS(output, input1, input2, 0);     
}

/*
   A truncating polynomial multiplication.
   The number of terms require, _trunc_ can be any value, but the function is
   tuned for truncation to length n where both inputs have length approximately n.
*/

void _fmpz_poly_mul_trunc_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                const fmpz_poly_t input2, const unsigned long trunc)
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
   
   unsigned long bits1 = _fmpz_poly_max_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _fmpz_poly_max_bits(input2);
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
      _fmpz_poly_mul_KS_trunc(output, input1, input2, trunc, 0);
      return;
   } 
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _fmpz_poly_mul_SS_trunc(output, input1, input2, trunc);
      return;
   } 
   
   _fmpz_poly_mul_KS_trunc(output, input1, input2, trunc, 0);     
}

/*
   A truncating polynomial multiplication which ignores the first trunc coeffs of
   the output (which can end up being anything - often zero).
   The number of zero terms, _trunc_ can be any value, but the function is
   tuned for truncation of length n-1 terms, where both inputs have length approximately n.
*/

void _fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                const fmpz_poly_t input2, const unsigned long trunc)
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
   
   unsigned long bits1 = _fmpz_poly_max_bits(input1);
   unsigned long bits2 = (input1 == input2) ? bits1 : _fmpz_poly_max_bits(input2);
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
      _fmpz_poly_mul_KS(output, input1, input2, 0);
      return;
   } 
   
   if (3*(bits1 + bits2) >= input1->length + input2->length)
   {
      _fmpz_poly_mul_SS(output, input1, input2);
      return;
   } 
   
   _fmpz_poly_mul_KS(output, input1, input2, 0); 
}

/*===========================================================================

   fmpz_poly_* layer

===========================================================================*/

/****************************************************************************

   Memory management

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

void fmpz_poly_init2(fmpz_poly_t poly, const unsigned long alloc, const unsigned long limbs)
{
   if (((long)alloc > 0) && ((long)limbs > 0))
   {
      poly->coeffs = (fmpz_t) flint_heap_alloc(alloc*(limbs+1));
   }
   else poly->coeffs = NULL;
   
   poly->alloc = alloc;
   poly->length = 0;
   poly->limbs = limbs;
}

/* 
   Shrink or expand a polynomial to "alloc" coefficients 
*/

void fmpz_poly_realloc(fmpz_poly_t poly, const unsigned long alloc)
{
   if (poly->limbs > 0)
   {
      if ((long)alloc > 0)
      {
         if (poly->alloc) poly->coeffs = (mp_limb_t*) flint_heap_realloc(poly->coeffs, alloc*(poly->limbs+1));
         else poly->coeffs = (mp_limb_t*) flint_heap_alloc(alloc*(poly->limbs+1));
      } else
      {
         if (poly->coeffs) flint_heap_free(poly->coeffs);
         poly->coeffs = NULL;
         poly->limbs = 0;
      }   
      poly->alloc = alloc;
   
      // truncate actual data if necessary
      if (poly->length > alloc)
      {
         poly->length = alloc;
         _fmpz_poly_normalise(poly);
      }     
   } else
   {
      poly->alloc = alloc;
   }
}

void fmpz_poly_fit_length(fmpz_poly_t poly, const unsigned long alc)
{
   unsigned long alloc = alc;
   if (alloc <= poly->alloc) return;

   if (alloc < 2*poly->alloc) alloc = 2*poly->alloc;
   
   fmpz_poly_realloc(poly, alloc);
}

void fmpz_poly_resize_limbs(fmpz_poly_t poly, const unsigned long limbs)
{
   if ((long)limbs > 0)
   {
      if (limbs == poly->limbs) return;
      
      unsigned long i = 0;
      fmpz_t coeff_i;
      fmpz_t coeff_i_old = poly->coeffs;
      
      if (limbs < poly->limbs)
      {
         coeff_i = poly->coeffs + limbs+1;
         coeff_i_old += (poly->limbs+1);
         for (i = 1; i < poly->length; i++)
         {
            F_mpn_copy_forward(coeff_i, coeff_i_old, limbs+1);
            FLINT_ASSERT(ABS(coeff_i[0]) > limbs); 
            coeff_i += (limbs+1);
            coeff_i_old += (poly->limbs+1);
         } 
      } else
      {
         if (poly->alloc)
         {
            fmpz_t temp_coeffs = (mp_limb_t*) flint_heap_alloc(poly->alloc*(limbs+1));
            coeff_i = temp_coeffs;
            for (i = 0; i < poly->length; i++)
            {
               F_mpn_copy(coeff_i, coeff_i_old, poly->limbs+1);
               coeff_i += (limbs+1);
               coeff_i_old += (poly->limbs+1);
            } 
            if (poly->coeffs) flint_heap_free(poly->coeffs);
            poly->coeffs = temp_coeffs;
         }
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

/****************************************************************************

   Polynomial Checking

****************************************************************************/

/*
   Used for debugging polynomial code
   Checks that length <= alloc and that both are positive or zero
   Checks that limbs >= 0 otherwise
   Checks that each coefficient has at most _limbs_ limbs 
*/

void fmpz_poly_check(const fmpz_poly_t poly)
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

void fmpz_poly_check_normalisation(const fmpz_poly_t poly)
{
   if (poly->length)
   {
      if (!poly->coeffs[(poly->length-1)*(poly->limbs+1)])
      {
         printf("Error: Poly not normalised\n");
         abort();
      }
   }
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

/****************************************************************************

   Coefficient setting and retrieval

****************************************************************************/

void fmpz_poly_get_coeff_mpz(mpz_t x, const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length)
      mpz_set_ui(x, 0);
   else
      _fmpz_poly_get_coeff_mpz(x, poly, n);
}

void fmpz_poly_get_coeff_mpz_read_only(mpz_t x, const fmpz_poly_t poly, const unsigned long n)
{
   if (n >= poly->length)
   {
      x->_mp_alloc = 1;
      x->_mp_d = (mp_limb_t *) &poly; // We need to point to something, and at least this exists
      x->_mp_size = 0;  
   } else 
      _fmpz_poly_get_coeff_mpz_read_only(x, poly, n);
}

/****************************************************************************

   String conversions and I/O

****************************************************************************/

int fmpz_poly_from_string(fmpz_poly_t poly, const char* s)
{
   int ok;
   
   mpz_poly_t p;
   mpz_poly_init(p);
   ok = mpz_poly_from_string(p, s);
   if (ok)
   {
      mpz_poly_to_fmpz_poly(poly, p);
   }
   mpz_poly_clear(p);
   
   return ok;
}

char* fmpz_poly_to_string(const fmpz_poly_t poly)
{
   char* buf;
   mpz_poly_t m_poly;
   mpz_poly_init(m_poly);
   fmpz_poly_to_mpz_poly(m_poly, poly);
   buf = mpz_poly_to_string(m_poly);
   mpz_poly_clear(m_poly);
   return buf;
}

char* fmpz_poly_to_string_pretty(const fmpz_poly_t poly, const char * x)
{
   char* buf;
   mpz_poly_t m_poly;
   mpz_poly_init(m_poly);
   fmpz_poly_to_mpz_poly(m_poly, poly);
   buf = mpz_poly_to_string_pretty(m_poly, x);
   mpz_poly_clear(m_poly);
   return buf;
}


void fmpz_poly_fprint(const fmpz_poly_t poly, FILE* f)
{
   char* s = fmpz_poly_to_string(poly);
   fputs(s, f);
   free(s);
}

void fmpz_poly_fprint_pretty(const fmpz_poly_t poly, FILE* f, const char * x)
{
   char* s = fmpz_poly_to_string_pretty(poly, x);
   fputs(s, f);
   free(s);
}


void fmpz_poly_print(const fmpz_poly_t poly)
{
   fmpz_poly_fprint(poly, stdout);
}

void fmpz_poly_print_pretty(const fmpz_poly_t poly, const char * x)
{
   fmpz_poly_fprint_pretty(poly, stdout, x);
}


int fmpz_poly_fread(fmpz_poly_t poly, FILE* f)
{
   int ok;
   
   mpz_poly_t p;
   mpz_poly_init(p);
   ok = mpz_poly_fread(p, f);
   if (ok)
   {
      mpz_poly_to_fmpz_poly(poly, p);
   }
   mpz_poly_clear(p);
   
   return ok;
}

/****************************************************************************

   Scalar multiplications and divisions

****************************************************************************/

void fmpz_poly_scalar_mul_ui(fmpz_poly_t output, 
                          const fmpz_poly_t input, unsigned long x)
{
   if ((input->length == 0) || (x == 0))
   {
      _fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs, bits, max_limbs, max_bits, x_bits, top_bits;
   
   max_bits = (input->limbs << FLINT_LG_BITS_PER_LIMB);
   
   max_limbs = 0;
   bits = 0;
   top_bits = 0;
   x_bits = FLINT_BIT_COUNT(x);
   fmpz_t next_coeff = input->coeffs;
   unsigned long size = input->limbs+1; 
   unsigned long i;
   
   for (i = 0; (i < input->length) && (top_bits + x_bits <= max_bits); i++)
   {
      limbs = ABS(next_coeff[0]);
      if ((limbs >= max_limbs) && (limbs))
      {
         max_limbs = limbs;
         bits = ((limbs - 1) << FLINT_LG_BITS_PER_LIMB) + FLINT_BIT_COUNT(next_coeff[limbs]);
         if (bits > top_bits) top_bits = bits;
      }
      next_coeff += size;
   }
   
   fmpz_poly_fit_length(output, input->length);
   if (i < input->length)
   {
      fmpz_poly_fit_limbs(output, input->limbs + 1);
   } else
   {
      fmpz_poly_fit_limbs(output, ((top_bits + x_bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1);
   }
   
   _fmpz_poly_scalar_mul_ui(output, input, x);
}

void fmpz_poly_scalar_mul_si(fmpz_poly_t output, 
                          const fmpz_poly_t input, long x)
{
   if ((input->length == 0) || (x == 0))
   {
      _fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs, bits, max_limbs, max_bits, x_bits, top_bits;
   
   max_bits = (input->limbs << FLINT_LG_BITS_PER_LIMB);
   
   max_limbs = 0;
   bits = 0;
   top_bits = 0;
   x_bits = FLINT_BIT_COUNT(FLINT_ABS(x));
   fmpz_t next_coeff = input->coeffs;
   unsigned long size = input->limbs+1; 
   unsigned long i;
   
   for (i = 0; (i < input->length) && (top_bits + x_bits <= max_bits); i++)
   {
      limbs = ABS(next_coeff[0]);
      if ((limbs >= max_limbs) && (limbs))
      {
         max_limbs = limbs;
         bits = ((limbs - 1) << FLINT_LG_BITS_PER_LIMB) + FLINT_BIT_COUNT(next_coeff[limbs]);
         if (bits > top_bits) top_bits = bits;
      }
      next_coeff += size;
   }
   
   fmpz_poly_fit_length(output, input->length);
   if (i < input->length)
   {
      fmpz_poly_fit_limbs(output, input->limbs + 1);
   } else
   {
      fmpz_poly_fit_limbs(output, ((top_bits + x_bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1);
   }
   
   _fmpz_poly_scalar_mul_si(output, input, x);
}

void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const fmpz_t x)
{
   if ((input->length == 0) || (x[0] == 0))
   {
      _fmpz_poly_zero(output);
      return;
   }
   unsigned long limbs, bits, max_limbs, max_bits, x_bits, top_bits;
   
   max_bits = ((input->limbs + ABS(x[0]) - 1) << FLINT_LG_BITS_PER_LIMB);
   
   max_limbs = 0;
   bits = 0;
   top_bits = 0;
   x_bits = ABS(fmpz_bits(x));
   fmpz_t next_coeff = input->coeffs;
   unsigned long size = input->limbs+1; 
   unsigned long i;
   
   for (i = 0; (i < input->length) && (top_bits + x_bits <= max_bits); i++)
   {
      limbs = ABS(next_coeff[0]);
      if ((limbs >= max_limbs) && (limbs))
      {
         max_limbs = limbs;
         bits = ((limbs - 1) << FLINT_LG_BITS_PER_LIMB) + FLINT_BIT_COUNT(next_coeff[limbs]);
         if (bits > top_bits) top_bits = bits;
      }
      next_coeff += size;
   }
   
   fmpz_poly_fit_length(output, input->length);
   if (i < input->length)
   {
      fmpz_poly_fit_limbs(output, input->limbs + ABS(x[0]));
   } else
   {
      fmpz_poly_fit_limbs(output, ((top_bits + x_bits - 1) >> FLINT_LG_BITS_PER_LIMB) + 1);
   }
   
   _fmpz_poly_scalar_mul_fmpz(output, input, x);
}

void fmpz_poly_scalar_mul_mpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const mpz_t x)
{
   if ((input->length == 0) || (mpz_sgn(x) == 0))
   {
      _fmpz_poly_zero(output);
      return;
   }
   fmpz_t x_fmpz = fmpz_init(mpz_size(x));
   mpz_to_fmpz(x_fmpz, x);   
   fmpz_poly_scalar_mul_fmpz(output, input, x_fmpz);
   fmpz_clear(x_fmpz);
}

void fmpz_poly_scalar_div_fmpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const fmpz_t x)
{
   if (input->length == 0)
   {
      _fmpz_poly_zero(output);
      return;
   }
    
   fmpz_poly_fit_length(output, input->length);
   
   unsigned long inlimbs = _fmpz_poly_max_limbs(input);
   unsigned long xlimbs = ABS(x[0]);
   
   if (inlimbs >= xlimbs) fmpz_poly_fit_limbs(output, inlimbs - xlimbs + 1);
   else fmpz_poly_fit_limbs(output, 1);
   
   _fmpz_poly_scalar_div_fmpz(output, input, x);
}

void fmpz_poly_scalar_div_mpz(fmpz_poly_t output, 
                          const fmpz_poly_t input, const mpz_t x)
{
   if (input->length == 0)
   {
      _fmpz_poly_zero(output);
      return;
   }
   fmpz_t x_fmpz = fmpz_init(mpz_size(x));
   mpz_to_fmpz(x_fmpz, x);   
   fmpz_poly_scalar_div_fmpz(output, input, x_fmpz);
   fmpz_clear(x_fmpz);
}


/****************************************************************************

   Multiplication

****************************************************************************/

void fmpz_poly_mul(fmpz_poly_t output, const fmpz_poly_t input1, const fmpz_poly_t input2)
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
      
   bits1 = _fmpz_poly_max_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_max_bits(input2);
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = (input1->length > input2->length) ? input2->length : input1->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   
   fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+2);
   fmpz_poly_fit_length(output, input1->length + input2->length - 1);
   
   _fmpz_poly_mul(output, input1, input2);
}

void fmpz_poly_mul_trunc_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                                          const fmpz_poly_t input2, const unsigned long trunc)
{
   long bits1, bits2;
      
   bits1 = _fmpz_poly_max_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_max_bits(input2);
     
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = (input1->length > input2->length) ? input2->length : input1->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   
   if (bits == 0)
   {
      _fmpz_poly_zero(output);
      return;
   }
   fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+1);
   fmpz_poly_fit_length(output, FLINT_MIN(input1->length + input2->length - 1, trunc));
   
   _fmpz_poly_mul_trunc_n(output, input1, input2, FLINT_MIN(input1->length + input2->length - 1, trunc));
}

void fmpz_poly_mul_trunc_left_n(fmpz_poly_t output, const fmpz_poly_t input1, 
                             const fmpz_poly_t input2, const unsigned long trunc)
{
   unsigned long limbs = input1->limbs + input2->limbs;
   
   long bits1, bits2;
      
   bits1 = _fmpz_poly_max_bits(input1);
   bits2 = (input1 == input2) ? bits1 : _fmpz_poly_max_bits(input2);
      
   unsigned long sign = ((bits1 < 0) || (bits2 < 0));
   unsigned long length = (input1->length > input2->length) ? input2->length : input1->length;
   unsigned log_length = 0;
   while ((1<<log_length) < length) log_length++;
   unsigned long bits = ABS(bits1) + ABS(bits2) + log_length + sign; 
   
   fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+1);
   if (input1->length + input2->length) fmpz_poly_fit_length(output, input1->length + input2->length - 1);
   
   _fmpz_poly_mul_trunc_left_n(output, input1, input2, trunc);
}

/****************************************************************************

   Division

****************************************************************************/

void fmpz_poly_divrem_classical(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
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
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t coeff_i = coeffs_A + coeff*size_A;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeffs_R;

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
   
   fmpz_t rem = (fmpz_t) flint_heap_alloc(size_B_lead);
   
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
                  fmpz_sub_ui_inplace(coeff_Q,1);
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
         _fmpz_poly_scalar_mul_fmpz(qB, B, coeff_Q); 
         
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

void fmpz_poly_divrem_classical_low(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
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
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t coeff_i = coeffs_A + coeff*size_A;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeffs_R;

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
   
   fmpz_t rem = (fmpz_t) flint_heap_alloc(size_B_lead);
   
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
                  fmpz_sub_ui_inplace(coeff_Q,1);
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
         fmpz_poly_init2(qB, B->length-1, B->limbs+ABS(coeff_Q[0]));
         temp->coeffs = B->coeffs;
         temp->length = B->length - 1;
         temp->limbs = B->limbs;
         _fmpz_poly_scalar_mul_fmpz(qB, temp, coeff_Q); 
         
         fmpz_poly_fit_limbs(R, qB->limbs+1);
         coeffs_R = R->coeffs;
         size_R = R->limbs+1;
      
         fmpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R+(coeff - B->length + 1)*size_R;
         R_sub->limbs = R->limbs;
         R_sub->length = B->length - 1;
         _fmpz_poly_sub(R_sub, R_sub, qB);
         coeffs_R[coeff*size_R] = 0L;
         
         coeff--;
         fmpz_poly_clear(qB);
      }
   }
   
   R->length = B->length - 1;
   _fmpz_poly_normalise(R);
   flint_heap_free(rem);
}

/*
   Divide the polynomial A by the polynomial B but do not compute the remainder
*/

void fmpz_poly_div_classical(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB, R;
   
   fmpz_poly_init(R);
   
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
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t coeff_i = coeffs_A + coeff*size_A;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeffs_R;

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
   
   fmpz_t rem = (fmpz_t) flint_heap_alloc(size_B_lead);
   
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
                  fmpz_sub_ui_inplace(coeff_Q,1);
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
            
            fmpz_poly_init2(qB, length, B->limbs+ABS(coeff_Q[0])+1);
            R_sub->coeffs = B->coeffs + (B->length - length)*(B->limbs + 1);
            R_sub->limbs = B->limbs;
            R_sub->length = length;
            _fmpz_poly_scalar_mul_fmpz(qB, R_sub, coeff_Q); 
            
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
   Integer polynomial division using a divide and conquer algorithm.
   Note BQ is not the remainder but it is B*Q, so the remainder R = A-BQ
*/
 
void fmpz_poly_div_divconquer_recursive(fmpz_poly_t Q, fmpz_poly_t BQ, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      _fmpz_poly_zero(BQ);

      return;
   }
   
   // A->length is now >= B->length
   
   unsigned long crossover = 16;
   unsigned long crossover2 = 128;
   
   if (B->limbs > 16) crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      /*
         Use the classical algorithm to compute the
         quotient and remainder, then use A-R to compute BQ
      */
      
      fmpz_poly_t Rb;
      fmpz_poly_init(Rb);
      fmpz_poly_divrem_classical(Q, Rb, A, B);
      fmpz_poly_fit_length(BQ, A->length);
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(A->limbs, Rb->limbs)+1);
      _fmpz_poly_sub(BQ, A, Rb);
      fmpz_poly_clear(Rb);
      
      return;
   }
   
   fmpz_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _fmpz_poly_attach_shift(d1, B, n2);
   _fmpz_poly_attach_truncate(d2, B, n2);
   _fmpz_poly_attach_shift(d3, B, n1);
   _fmpz_poly_attach_truncate(d4, B, n1);
   
   if (A->length <= n2 + B->length - 1)
   {
      /*
         If A->length <= B->length + n2 - 1
         then only a single quotient is needed
         We do a division of at most 2*n2 - 1
         terms by n2 terms yielding a quotient of
         at most n2 terms 
      */
      
      // Set p1 to be A without the last
      // n1 coefficients
      // 2*n2-1 >= p1->length > 0
      
      fmpz_poly_init(p1);
      fmpz_poly_fit_length(p1, A->length-n1);
      fmpz_poly_fit_limbs(p1, A->limbs);
      _fmpz_poly_right_shift(p1, A, n1);
      
      // Since A was normalised, then p1 will be
      // d3 is the leading terms of B and so must be normalised
      // d3 is length n2, so we get at most n2 terms in the quotient
      
      fmpz_poly_init(d1q1);
      fmpz_poly_div_divconquer_recursive(Q, d1q1, p1, d3); 
      fmpz_poly_clear(p1);
      
      /*
         Compute d2q1 = Q*d4
         It is of length at most n1+n2-1 terms
      */
      
      fmpz_poly_init(d2q1);
      fmpz_poly_mul(d2q1, Q, d4);
      
      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length at most n1+2*n2-1
      */
      
      fmpz_poly_fit_length(BQ, FLINT_MAX(d1q1->length+n1, d2q1->length));
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(d1q1->limbs, d2q1->limbs)+1);
      _fmpz_poly_left_shift(BQ, d1q1, n1);
      fmpz_poly_clear(d1q1);
      _fmpz_poly_add(BQ, BQ, d2q1);
      fmpz_poly_clear(d2q1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      unsigned long shift = A->length - 2*B->length + 1;
      _fmpz_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = d1*q1 is length at most 2*B->length-1
      */
      
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
      
      fmpz_poly_div_divconquer_recursive(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      fmpz_poly_init(dq1);
      
      fmpz_poly_fit_length(dq1, d1q1->length + shift);
      fmpz_poly_fit_limbs(dq1, d1q1->limbs);
      _fmpz_poly_left_shift(dq1, d1q1, shift);
      fmpz_poly_clear(d1q1);
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      fmpz_poly_init(t);
      fmpz_poly_sub(t, A, dq1);
      _fmpz_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      fmpz_poly_init(q2);
      fmpz_poly_init(dq2);
      fmpz_poly_div_divconquer_recursive(q2, dq2, t, B); 
      fmpz_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      fmpz_poly_fit_length(Q, FLINT_MAX(q1->length+shift, q2->length));
      fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   
      _fmpz_poly_left_shift(Q, q1, shift);
      fmpz_poly_clear(q1);
      _fmpz_poly_add(Q, Q, q2);
      fmpz_poly_clear(q2);
      
      /*
         Write out BQ = dq1 + dq2
      */
      
      fmpz_poly_fit_length(BQ, FLINT_MAX(dq1->length, dq2->length));
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(dq1->limbs, dq2->limbs)+1);
      
      _fmpz_poly_add(BQ, dq1, dq2);
      fmpz_poly_clear(dq1);
      fmpz_poly_clear(dq2);
      
      return;
   } 
   
   // n2 + B->length - 1 < A->length <= n1 + n2 + B->length - 1
    
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length at most n1 (and at least 1), 
      a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n1-1)+ other terms, so it has 
      length at most 2*n1-1 
   */
      
   _fmpz_poly_stack_init(p1, A->length-2*n2, A->limbs);
   _fmpz_poly_right_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length n1
      d1q1 = d1*q1 is length at most 2*n1-1
   */
      
   fmpz_poly_init(d1q1);
   fmpz_poly_init(q1);
   fmpz_poly_div_divconquer_recursive(q1, d1q1, p1, d1); 
   _fmpz_poly_stack_clear(p1);   
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being at most length n1+n2-1
   */  
   
   fmpz_poly_init(d2q1);
   fmpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most 2*n1+n2-1
   */
   
   _fmpz_poly_stack_init(dq1, FLINT_MAX(d1q1->length + n2, d2q1->length), FLINT_MAX(d1q1->limbs, d2q1->limbs)+1);
   _fmpz_poly_left_shift(dq1, d1q1, n2);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);
   fmpz_poly_clear(d2q1);
   
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n1-1) - dq1
      which has length at most 2*n1+n2-1, but we are not interested 
      in up to the first n1 coefficients, so it has 
      effective length at most n1+n2-1
   */
   
   _fmpz_poly_stack_init(t, FLINT_MAX(A->length-n2, dq1->length), FLINT_MAX(A->limbs, dq1->limbs)+1);
   _fmpz_poly_right_shift(t, A, n2);
   _fmpz_poly_sub(t, t, dq1);
   _fmpz_poly_truncate(t, B->length - 1);
   
   /*
      Compute q2 = t div d1
      It is at most an n1+n2-1 by n1 division, so
      the length of q2 will be at most n2
      Also compute d1q2 of length at most n1+n2-1
   */
   
   fmpz_poly_init(d1q2);
   fmpz_poly_init(q2);
   fmpz_poly_div_divconquer_recursive(q2, d1q2, t, d1); 
   _fmpz_poly_stack_clear(t);
      
   /*
      Compute d2q2 = d2*q2 which is of length 
      at most n1+n2-1
   */
   
   fmpz_poly_init(d2q2);
   fmpz_poly_mul(d2q2, d2, q2);
   
   /*
      Compute dq2 = d1*q2*x^n2 + d2q2
      which is of length at most n1+2*n2-1
   */
   
   _fmpz_poly_stack_init(dq2, FLINT_MAX(d1q2->length+n2, d2q2->length), FLINT_MAX(d1q2->limbs, d2q2->limbs)+1);
   _fmpz_poly_left_shift(dq2, d1q2, n2);
   fmpz_poly_clear(d1q2);
   _fmpz_poly_add(dq2, dq2, d2q2);
   fmpz_poly_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length at most n1+n2
   */
   
   fmpz_poly_fit_length(Q, FLINT_MAX(q1->length+n2, q2->length));
   fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   _fmpz_poly_left_shift(Q, q1, n2);
   fmpz_poly_clear(q1);
   _fmpz_poly_add(Q, Q, q2);
   fmpz_poly_clear(q2);
   
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length at most 2*(n1+n2)-1
   */
   
   fmpz_poly_fit_length(BQ, FLINT_MAX(n2+dq1->length, dq2->length));
   fmpz_poly_fit_limbs(BQ, FLINT_MAX(dq1->limbs, dq2->limbs)+1);
   _fmpz_poly_left_shift(BQ, dq1, n2);
   _fmpz_poly_add(BQ, BQ, dq2);
   
   _fmpz_poly_stack_clear(dq2);
   _fmpz_poly_stack_clear(dq1);
}

/*
   Divide and conquer division of A by B but only computing the low half of Q*B
*/

void fmpz_poly_div_divconquer_recursive_low(fmpz_poly_t Q, fmpz_poly_t BQ, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      _fmpz_poly_zero(BQ);
      
      return;
   }
   
   // A->length is now >= B->length
   
   unsigned long crossover = 16;
   unsigned long crossover2 = 128;
   
   if (B->limbs > 16) crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      /*
         Use the classical algorithm to compute the
         quotient and low half of the remainder, then 
         truncate A-R to compute BQ
      */
      
      fmpz_poly_t Rb;
      fmpz_poly_init(Rb);
      fmpz_poly_divrem_classical_low(Q, Rb, A, B);
      fmpz_poly_fit_length(BQ, A->length);
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(A->limbs, Rb->limbs)+1);
      _fmpz_poly_sub(BQ, A, Rb);
      fmpz_poly_clear(Rb);
      _fmpz_poly_truncate(BQ, B->length - 1);
      
      return;
   }
   
   fmpz_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _fmpz_poly_attach_shift(d1, B, n2);
   _fmpz_poly_attach_truncate(d2, B, n2);
   _fmpz_poly_attach_shift(d3, B, n1);
   _fmpz_poly_attach_truncate(d4, B, n1);
   
   if (A->length <= n2 + B->length - 1)
   {
      /*
         If A->length <= B->length + n2 - 1
         then only a single quotient is needed
         We do a division of at most 2*n2 - 1
         terms by n2 terms yielding a quotient of
         at most n2 terms 
      */
      
      // Set p1 to be A without the last
      // n1 coefficients
      // 2*n2-1 >= p1->length > 0
      
      fmpz_poly_init(p1);
      fmpz_poly_fit_length(p1, A->length-n1);
      fmpz_poly_fit_limbs(p1, A->limbs);
      _fmpz_poly_right_shift(p1, A, n1);
      
      // Since A was normalised, then p1 will be
      // d3 is the leading terms of B and so must be normalised
      // d3 is length n2, so we get at most n2 terms in the quotient
      // We compute only the low n2-1 terms of the product d1q1
      
      fmpz_poly_init(d1q1);
      fmpz_poly_div_divconquer_recursive_low(Q, d1q1, p1, d3); 
      fmpz_poly_clear(p1);
      
      /*
         Compute d2q1 = Q*d4
         It is of length at most n1+n2-1 terms
      */
      
      fmpz_poly_init(d2q1);
      fmpz_poly_mul(d2q1, Q, d4);
      
      /*
         Compute BQ = d1q1*x^n1 + d2q1
         It has length at most n1+n2-1
      */
      
      fmpz_poly_fit_length(BQ, FLINT_MAX(d1q1->length+n1, d2q1->length));
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(d1q1->limbs, d2q1->limbs)+1);
      _fmpz_poly_left_shift(BQ, d1q1, n1);
      fmpz_poly_clear(d1q1);
      _fmpz_poly_add(BQ, BQ, d2q1);
      fmpz_poly_clear(d2q1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length - 1
      // We call this polynomial p1
      
      unsigned long shift = A->length - 2*B->length + 1;
      _fmpz_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = d1*q1 is truncated to length at most B->length-1
      */
      
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
      
      fmpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length - B->length
         dq1 is normalised since d1q1 was
      */
   
      fmpz_poly_init(dq1);
      
      fmpz_poly_fit_length(dq1, d1q1->length + shift);
      fmpz_poly_fit_limbs(dq1, d1q1->limbs);
      _fmpz_poly_left_shift(dq1, d1q1, shift);
      fmpz_poly_clear(d1q1);
      
      /*
         Compute t = A - dq1 
         We truncate, leaving at most A->length - B->length 
         significant terms
      */
   
      fmpz_poly_init(t);
      fmpz_poly_sub(t, A, dq1);
      _fmpz_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
         dq2 has length at most B->length - 1
      */
   
      fmpz_poly_init(q2);
      fmpz_poly_init(dq2);
      fmpz_poly_div_divconquer_recursive_low(q2, dq2, t, B); 
      fmpz_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      fmpz_poly_fit_length(Q, FLINT_MAX(q1->length+shift, q2->length));
      fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   
      _fmpz_poly_left_shift(Q, q1, shift);
      fmpz_poly_clear(q1);
      _fmpz_poly_add(Q, Q, q2);
      fmpz_poly_clear(q2);
      
      /*
         Write out BQ = dq1 + dq2
      */
      
      fmpz_poly_fit_length(BQ, FLINT_MAX(dq1->length, dq2->length));
      fmpz_poly_fit_limbs(BQ, FLINT_MAX(dq1->limbs, dq2->limbs)+1);
      
      _fmpz_poly_add(BQ, dq1, dq2);
      _fmpz_poly_truncate(BQ, B->length - 1);
      fmpz_poly_clear(dq1);
      fmpz_poly_clear(dq2);
      
      return;
   } 
   
   // n2 + B->length - 1 < A->length <= n1 + n2 + B->length - 1
    
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length at most n1 (and at least 1), 
      a2 is length n2 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n1-1)+ other terms, so it has 
      length at most 2*n1-1 
   */
      
   _fmpz_poly_stack_init(p1, A->length-2*n2, A->limbs);
   _fmpz_poly_right_shift(p1, A, 2*n2);
      
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length n1
      d1q1 = d1*q1 is truncated to length at most n1-1
   */
      
   fmpz_poly_init(d1q1);
   fmpz_poly_init(q1);
   fmpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, d1); 
   _fmpz_poly_stack_clear(p1);   
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being at most length n1+n2-1
   */  
   
   fmpz_poly_init(d2q1);
   fmpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most n1+n2-1
   */
   
   _fmpz_poly_stack_init(dq1, FLINT_MAX(d1q1->length + n2, d2q1->length), FLINT_MAX(d1q1->limbs, d2q1->limbs)+1);
   _fmpz_poly_left_shift(dq1, d1q1, n2);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);
   fmpz_poly_clear(d2q1);
   
   /*
      Compute t = a1*x^(n1+n2-1) + a2*x^(n1-1) - dq1
      which has length at most 2*n1+n2-1, but we are not interested 
      in up to the first n1 coefficients, so it has 
      effective length at most n1+n2-1
   */
   
   _fmpz_poly_stack_init(t, FLINT_MAX(A->length-n2, dq1->length), FLINT_MAX(A->limbs, dq1->limbs)+1);
   _fmpz_poly_right_shift(t, A, n2);
   _fmpz_poly_sub(t, t, dq1);
   _fmpz_poly_truncate(t, B->length - 1);
   
   /*
      Compute q2 = t div d1
      It is at most an n1+n2-1 by n1 division, so
      the length of q2 will be at most n2
      Also compute d1q2 truncated to length at most n1-1
   */
   
   fmpz_poly_init(d1q2);
   fmpz_poly_init(q2);
   fmpz_poly_div_divconquer_recursive_low(q2, d1q2, t, d1); 
   _fmpz_poly_stack_clear(t);
      
   /*
      Compute d2q2 = d2*q2 which is of length 
      at most n1+n2-1
   */
   
   fmpz_poly_init(d2q2);
   fmpz_poly_mul(d2q2, d2, q2);  
   
   /*
      Compute dq2 = d1*q2*x^n2 + d2q2
      which is of length at most n1+n2-1
   */
   
   _fmpz_poly_stack_init(dq2, FLINT_MAX(d1q2->length+n2, d2q2->length), FLINT_MAX(d1q2->limbs, d2q2->limbs)+1);
   _fmpz_poly_left_shift(dq2, d1q2, n2);
   fmpz_poly_clear(d1q2);
   _fmpz_poly_add(dq2, dq2, d2q2);
   fmpz_poly_clear(d2q2);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length at most n1+n2
   */
   
   fmpz_poly_fit_length(Q, FLINT_MAX(q1->length+n2, q2->length));
   fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   _fmpz_poly_left_shift(Q, q1, n2);
   fmpz_poly_clear(q1);
   _fmpz_poly_add(Q, Q, q2);
   fmpz_poly_clear(q2);
   
   /*
      Write out BQ = dq1*x^n2 + dq2
      BQ has length at most n1+2*n2-1
      We truncate to at most length B->length - 1
   */
   
   fmpz_poly_fit_length(BQ, FLINT_MAX(n2+dq1->length, dq2->length));
   fmpz_poly_fit_limbs(BQ, FLINT_MAX(dq1->limbs, dq2->limbs)+1);
   _fmpz_poly_left_shift(BQ, dq1, n2);
   _fmpz_poly_add(BQ, BQ, dq2);
   _fmpz_poly_truncate(BQ, B->length - 1);
   
   _fmpz_poly_stack_clear(dq2);
   _fmpz_poly_stack_clear(dq1);
}

void fmpz_poly_div_divconquer(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      
      return;
   }

   // A->length is now >= B->length
    
   unsigned long crossover = 16;
   unsigned long crossover2 = 256;
   
   if (B->limbs > 16) crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      fmpz_poly_div_classical(Q, A, B);
      
      return;
   }
   
   // B->length is now >= crossover (8 or 16)
   
   fmpz_poly_t d1, d2, d3, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
      
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   // n1 and n2 are at least 4
   
   /* We let B = d1*x^n2 + d2 
      d1 is of length n1 and
      d2 of length n2
   */

   _fmpz_poly_attach_shift(d1, B, n2);
   _fmpz_poly_attach_truncate(d2, B, n2);
   _fmpz_poly_attach_shift(d3, B, n1);
   
   if (A->length <= n2 + B->length - 1)
   {
      /*
         If A->length <= B->length + n2 - 1
         then only a single quotient is needed
         We do a division of at most 2*n2 - 1
         terms by n2 terms yielding a quotient of
         at most n2 terms 
      */
      
      // Set p1 to be A without the last
      // n1 coefficients
      // 2*n2-1 >= p1->length > 0
      
      fmpz_poly_init(p1);
      fmpz_poly_fit_length(p1, A->length-n1);
      fmpz_poly_fit_limbs(p1, A->limbs);
      _fmpz_poly_right_shift(p1, A, n1);
      
      // Since A was normalised, then p1 will be
      // d3 is the leading terms of B and so must be normalised
      // d3 is length n2, so we get at most n2 terms in the quotient
      
      fmpz_poly_div_divconquer(Q, p1, d3); 
      fmpz_poly_clear(p1);
            
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length -1
      // We call this polynomial p1
      
      unsigned long shift = A->length - 2*B->length + 1;
      _fmpz_poly_attach_shift(p1, A, shift);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         d1q1 = low(d1*q1) is length at most 2*B->length-1
         We discard the lower B->length-1 terms
      */
      
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
      
      fmpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, B); 
       
      /* 
         Compute dq1 = d1*q1*x^shift
         dq1 is then of length at most A->length
         dq1 is normalised since d1q1 was
      */
   
      fmpz_poly_init(dq1);
      
      fmpz_poly_fit_length(dq1, d1q1->length + shift);
      fmpz_poly_fit_limbs(dq1, d1q1->limbs);
      _fmpz_poly_left_shift(dq1, d1q1, shift);
      fmpz_poly_clear(d1q1); 
      
      /*
         Compute t = A - dq1 
         The first B->length coefficients cancel
         if the division is exact, leaving
          A->length - B->length significant terms
         otherwise we truncate at this length 
      */
   
      fmpz_poly_init(t);
      fmpz_poly_sub(t, A, dq1);
      fmpz_poly_clear(dq1);
      _fmpz_poly_truncate(t, A->length - B->length);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length-B->length
      */
   
      fmpz_poly_init(q2);
      fmpz_poly_div_divconquer(q2, t, B); 
      fmpz_poly_clear(t);  
      
      /*
         Write out Q = q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
      */
   
      fmpz_poly_fit_length(Q, FLINT_MAX(q1->length+shift, q2->length));
      fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   
      _fmpz_poly_left_shift(Q, q1, shift);
      fmpz_poly_clear(q1);
      _fmpz_poly_add(Q, Q, q2);
      fmpz_poly_clear(q2);
      
      return;
   }
   // We now have n2 + B->length - 1 < A->length <= 2*B->length - 1
   
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length at most n1 and a2 is length n2 
      and a3 is length n1+n2-1 
   */
      
      // Set p1 to a1*x^(n1-1) + other terms
      // It has length at most 2*n1-1 and is normalised
      // A->length >= 2*n2
      
      fmpz_poly_init(p1);
      fmpz_poly_fit_length(p1, A->length - 2*n2);
      fmpz_poly_fit_limbs(p1, A->limbs);
      _fmpz_poly_right_shift(p1, A, 2*n2);
      
      /* 
         Set q1 to p1 div d1 
         This is at most a 2*n1-1 by n1 division so 
         q1 ends up being at most length n1
         d1q1 = low(d1*q1) is length at most n1-1
         Thus we have discarded the leading n1 terms (at most)
      */
      
      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
      
      fmpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, d1); 
      fmpz_poly_clear(p1);
   
   /* 
      Compute d2q1 = d2*q1 with low n1 - 1 terms zeroed
      d2*q1 is length at most n1+n2-1 leaving at most
      n2 non-zero terms to the left
   */  
   
   _fmpz_poly_stack_init(d2q1, d2->length+q1->length-1, d2->limbs+q1->limbs+1); 
   _fmpz_poly_mul_trunc_left_n(d2q1, d2, q1, n1 - 1);
       
   /* 
      Compute dq1 = d1*q1*x^n2 + d2*q1
      dq1 is then of length at most 2*n1+n2-1
      but may have any length below this        
   */
   
   _fmpz_poly_stack_init(dq1, FLINT_MAX(d1q1->length + n2, d2q1->length), B->limbs+q1->limbs+1);
   _fmpz_poly_left_shift(dq1, d1q1, n2);
   fmpz_poly_clear(d1q1); 
   _fmpz_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = a1*x^(2*n2-1) + a2*x^(n2-1) - dq1 
      after shifting dq1 to the right by (n1-n2)
      which has length at most 2*n1+n2-1, but we 
      discard up to n1 coefficients, so it has 
      effective length 2*n2-1 with the last n2-1
      coefficients ignored. Thus there are at most n2 
      significant coefficients
   */
   
   
   _fmpz_poly_stack_init(t, n1+2*n2-1, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_right_shift(t, A, n1);
   _fmpz_poly_attach_shift(temp, dq1, n1-n2);
   _fmpz_poly_sub(t, t, temp);
   _fmpz_poly_truncate(t, 2*n2-1);
     
   /*
      Compute q2 = t div d3
      It is at most a 2*n2-1 by n2 division, so
      the length of q2 will be n2 at most
   */
   
   fmpz_poly_init(q2);
   fmpz_poly_div_divconquer(q2, t, d3); 
   _fmpz_poly_stack_clear(t);  
   _fmpz_poly_stack_clear(dq1);
   _fmpz_poly_stack_clear(d2q1);
   
   /*
      Write out Q = q1*x^n2 + q2
      Q has length n1+n2
   */
   
   fmpz_poly_fit_length(Q, q1->length+n2);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   _fmpz_poly_left_shift(Q, q1, n2);
   fmpz_poly_clear(q1);
   _fmpz_poly_add(Q, Q, q2);
   fmpz_poly_clear(q2);   
}

void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t QB;
   
   fmpz_poly_init(QB);
   
   fmpz_poly_div_divconquer_recursive(Q, QB, A, B);
   
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

void fmpz_poly_newton_invert_basecase(fmpz_poly_t Q_inv, const fmpz_poly_t Q, unsigned long n)
{
   fmpz_poly_t X2n, Qn;
   
   fmpz_poly_init2(X2n, 2*n-1, 1);
   _fmpz_poly_zero_coeffs(X2n, 2*n - 2);
   _fmpz_poly_set_coeff_ui(X2n, 2*n - 2, 1);
   X2n->length = 2*n-1;
   
   Qn->coeffs = Q->coeffs + (Q->length - n)*(Q->limbs + 1);
   Qn->limbs = Q->limbs;
   Qn->length = n;
   
   fmpz_poly_div_mulders(Q_inv, X2n, Qn);
      
   fmpz_poly_clear(X2n);
}

#define FLINT_NEWTON_INVERSE_BASECASE_CUTOFF 32

/*
   Recursively compute 1 / Q mod x^n using Newton iteration
   Assumes Q is given to the full precision n required and has constant term 1
*/

void fmpz_poly_newton_invert(fmpz_poly_t Q_inv, const fmpz_poly_t Q, const unsigned long n)
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
   fmpz_sub_ui_inplace(prod->coeffs, 1);
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

void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, const unsigned long n)
{
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

   fmpz_poly_t B_inv;
   fmpz_poly_init(B_inv);
   fmpz_poly_newton_invert(B_inv, Bin, n);
   fmpz_poly_mul_trunc_n(Q, B_inv, Ain, n);
   
   fmpz_poly_clear(B_inv);

   if (A == Q) _fmpz_poly_stack_clear(Ain);
   if (B == Q) _fmpz_poly_stack_clear(Bin);
}

/*
   Polynomial division of A by B
   The remainder is not computed, to save time
   B is assumed to be monic
*/

void fmpz_poly_div_newton(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
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
   
   fmpz_poly_clear(B_rev);
   fmpz_poly_clear(A_rev);
}

/*===================================================================================

   Mulder's short division algorithm
   
====================================================================================*/

// Mulders algorithm without improvements

/*void fmpz_poly_div_mulders(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      
      return;
   }
   
   unsigned long crossover = 16;
   
   if (B->limbs > 16)  crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) || (A->length > 2*B->length - 1))
   {
      fmpz_poly_div_classical(Q, A, B);
      return;
   }
   
   unsigned long k;
   
   k = 0;
   /*if (B->length <= 100) k = B->length/5;
   if (B->length <= 20) k = B->length/4;
   if (B->length == 10) k = B->length/3;*/
   
/*   fmpz_poly_t d1, d2, g1, g2, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
      
#if MULDERS_NEGATIVE
   unsigned long n1 = (B->length+1)/2 - k; 
#else
   unsigned long n1 = (B->length+1)/2 + k; 
#endif
   unsigned long n2 = B->length - n1; 
   
   /* We let B = d1*x^n2 + d2 */
/*   d1->length = n1;
   d2->length = n2;
   g1->length = n2;
   g2->length = n1;
   d1->limbs = B->limbs;
   d2->limbs = B->limbs;
   g1->limbs = B->limbs;
   g2->limbs = B->limbs;
   d1->coeffs = B->coeffs + n2*(B->limbs+1);
   d2->coeffs = B->coeffs;
   g1->coeffs = B->coeffs + n1*(B->limbs+1);
   g2->coeffs = B->coeffs;
      
   if (A->length <= 2*n1+n2-1) 
   {
      temp->length = A->length - (n1+n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n1-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n1-1);
      p1->length = temp->length+n1-1;
      
      fmpz_poly_init(d1q1);
      _fmpz_poly_normalise(p1);
      fmpz_poly_div_divconquer_recursive(Q, d1q1, p1, d1); //******************************
      fmpz_poly_clear(d1q1);
      _fmpz_poly_stack_clear(p1);
      
      return;   
   } else
   {
   /* 
      We let A = a1*x^(2*n1+n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is length n2 and a2 is length n1 and a3 is length n1+n2-1 
      We set p1 = a1*x^(n2-1), so it has length 2*n2-1
   */
      
/*      temp->length = A->length - (2*n1+n2-1);
      temp->limbs = A->limbs;
      temp->coeffs = A->coeffs + (2*n1+n2-1)*(A->limbs+1);
      _fmpz_poly_stack_init(p1, temp->length+n2-1, A->limbs);
      _fmpz_poly_left_shift(p1, temp, n2-1);
      p1->length = temp->length+n2-1;
      
      /* 
         Set q1 to p1 div g1 
         This is a 2*n2-1 by n2 division so 
         q1 ends up being length n2
         g1q1 = g1*q1 is length 2*n2-1 but we retrieve only the low n2-1 terms
      */
/*      fmpz_poly_init(d1q1);
      fmpz_poly_init(q1);
   
      _fmpz_poly_normalise(p1);
      fmpz_poly_div_divconquer_recursive(q1, d1q1, p1, g1); //******************************
      _fmpz_poly_stack_clear(p1);
   }
   
   /* 
      Compute g2q1 = g2*q1 
      which ends up being length n1+n2-1 but we set the right most n2-1 terms to zero
   */  
   
/*   _fmpz_poly_stack_init(d2q1, g2->length+q1->length-1, g2->limbs+q1->limbs+1); 
   _fmpz_poly_mul_trunc_left_n(d2q1, g2, q1, n2 - 1);
     
   /* 
      Compute dq1 = g1*q1*x^n1 + g2*q1
      dq1 is then of length n1+2*n2-1 but we have only the rightmost n1+n2-1 terms
   */
   
/*   _fmpz_poly_stack_init(dq1, FLINT_MAX(d1q1->length + n1, d2q1->length), B->limbs+q1->limbs+1);
   
   _fmpz_poly_zero_coeffs(dq1, n1);
   dq1->length = d1q1->length + n1;
   temp->length = d1q1->length;
   temp->limbs = dq1->limbs;
   temp->coeffs = dq1->coeffs + n1*(dq1->limbs+1);
   _fmpz_poly_set(temp, d1q1);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n2-1) - dq1 
      which has length 2*n1+n2-1, but we are not interested 
      in the first n1 coefficients, so it has 
      effective length n1+n2-1
   */
   
/*   temp->length = A->length - (n1+n2-1);
   temp->limbs = A->limbs;
   temp->coeffs = A->coeffs + (n1+n2-1)*(A->limbs+1);
#if MULDERS_NEGATIVE
   _fmpz_poly_stack_init(t, n1+2*n2-1, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_left_shift(t, temp, n2-1);
   t->length = temp->length+n2-1;
   _fmpz_poly_sub(t, t, dq1);
   _fmpz_poly_right_shift(t, t, n2-n1);
#else
   _fmpz_poly_stack_init(t, 2*n1+n2-1, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_left_shift(t, temp, n1-1);
   t->length = temp->length+n1-1;
   temp->length = dq1->length;
   temp->limbs = t->limbs;
   temp->coeffs = t->coeffs + (n1-n2)*(t->limbs+1);
   _fmpz_poly_sub(temp, temp, dq1);
#endif
   t->length = 2*n1-1; 
   _fmpz_poly_normalise(t); 
     
   /*
      Compute q2 = t div d1
      It is a 2*n1-1 by n1 division, so
      the length of q2 will be n1
   */
/*   fmpz_poly_init(q2);
   _fmpz_poly_normalise(t);
   fmpz_poly_div_mulders(q2, t, d1); 
   _fmpz_poly_stack_clear(t);  
   _fmpz_poly_stack_clear(dq1);
   _fmpz_poly_stack_clear(d2q1);
      
   /*
      Write out Q = q1*x^n1 + q2
      Q has length n1+n2
   */
/*   fmpz_poly_fit_length(Q, q1->length+n1);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(q1->limbs, q2->limbs));
   _fmpz_poly_set(Q, q2);
   fmpz_poly_clear(q2);
   Q->length = q1->length + n1;
   temp->length = q1->length;
   temp->limbs = Q->limbs;
   temp->coeffs = Q->coeffs + n1*(Q->limbs+1);
   _fmpz_poly_set(temp, q1);
   fmpz_poly_clear(q1);
}*/

void fmpz_poly_div_mulders(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      
      return;
   }
   
   // Crossover must be at least 8 so that n2 is not zero
   
   unsigned long crossover = 16;
   unsigned long crossover2 = 256;
   
   if (B->limbs > 16) crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      fmpz_poly_div_classical(Q, A, B);
      return;
   }
   
   unsigned long k;
   
   k = 0;
   if (B->length <= 100) k = B->length/5;
   if (B->length <= 20) k = B->length/4;
   if (B->length == 10) k = B->length/3;
   
   fmpz_poly_t d1, d2, g1, g2, p1, q1, q2, dq1, dq2, d1q1, d2q1, d2q2, d1q2, t, temp;
      
   // We demand n2 not be zero, this holds since 
   // crossover is at least 8
   
   unsigned long n1 = (B->length+1)/2 + k; 
   unsigned long n2 = B->length - n1; 
   
   /* We let B = d1*x^n2 + d2 */
   
   _fmpz_poly_attach_shift(d1, B, n2);
   _fmpz_poly_attach_truncate(d2, B, n2);
   _fmpz_poly_attach_shift(g1, B, n1);
   _fmpz_poly_attach_truncate(g2, B, n1);
      
   if (A->length <= n1 + B->length - 1) 
   {
      /*
         We only need a single division so we
         shift and make a recursive call
         Since n2 is not zero the size has been reduced
      */
      
      _fmpz_poly_stack_init(p1, A->length - n2, A->limbs);
      _fmpz_poly_right_shift(p1, A, n2);
      
      fmpz_poly_div_mulders(Q, p1, d1); 
      _fmpz_poly_stack_clear(p1);
      
      return;   
   } 
   
   if (A->length > 2*B->length - 1)
   {
      fmpz_poly_div_divconquer(Q, A, B);
      
      return;              
   }              
   
   /* 
      We let A = a1*x^(2*n1+n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is at most length n2 and a2 is length n1 
      and a3 is length n1+n2-1 
      We set p1 = a1*x^(n2-1) + other terms, so it has 
      length at most 2*n2-1
      A->length is at least 2*n1 + n2 - 1 which is at 
      least the requisite 2*n1
   */
      
   _fmpz_poly_stack_init(p1, A->length-2*n2, A->limbs);
   _fmpz_poly_right_shift(p1, A, 2*n1);
      
   /* 
      Set q1 to p1 div g1 
      This is at most a 2*n2-1 by n2 division so 
      q1 ends up being at most length n2
      d1q1 = g1*q1 is length 2*n2-1 but we retrieve 
      only the low n2-1 terms
   */
   
   fmpz_poly_init(d1q1);
   fmpz_poly_init(q1);
      
   fmpz_poly_div_divconquer_recursive_low(q1, d1q1, p1, g1); 
   _fmpz_poly_stack_clear(p1);
   
   /* 
      Compute d2q1 = g2*q1 
      which ends up being at most length n1+n2-1 but we set 
      the right most n2-1 terms to zero
      g2->length cannot be zero since it is d1
   */  
   
   _fmpz_poly_stack_init(d2q1, g2->length+q1->length-1, g2->limbs+q1->limbs+1); 
   _fmpz_poly_mul_trunc_left_n(d2q1, g2, q1, n2 - 1);
     
   /* 
      Compute dq1 = g1*q1*x^n1 + g2*q1
      dq1 is then of length n1+2*n2-1 but we have only 
      the rightmost n1+n2-1 terms, the last n2-1 of 
      which are irrelevant
   */
   
   _fmpz_poly_stack_init(dq1, FLINT_MAX(d1q1->length + n1, d2q1->length), FLINT_MAX(d1q1->limbs, d2q1->limbs)+1);
   
   _fmpz_poly_left_shift(dq1, d1q1, n1);
   fmpz_poly_clear(d1q1);
   _fmpz_poly_add(dq1, dq1, d2q1);
   
   /*
      Compute t = p1*x^(n1+n2-1) + p2*x^(n2-1) - dq1 
      where dq1 has been shifted left by (n1-n2),
      which has length 2*n1+n2-1, but we are not interested 
      in the first n2 coefficients, so it has 
      effective length at most 2*n1-1
   */
   
   _fmpz_poly_stack_init(t, n1+B->length, FLINT_MAX(A->limbs,dq1->limbs)+1);
   _fmpz_poly_right_shift(t, A, n2);
   _fmpz_poly_attach_shift(temp, t, n1-n2);
   _fmpz_poly_sub(temp, temp, dq1);
   _fmpz_poly_truncate(t, 2*n1-1); 
     
   /*
      Compute q2 = t div d1
      It is at most a 2*n1-1 by n1 division, so
      the length of q2 will be at most n1
   */
   
   fmpz_poly_init(q2);
   fmpz_poly_div_mulders(q2, t, d1); 
   _fmpz_poly_stack_clear(t);  
   _fmpz_poly_stack_clear(dq1);
   _fmpz_poly_stack_clear(d2q1);
      
   /*
      Write out Q = q1*x^n1 + q2
      Q has length at most n1+n2
   */
   
   fmpz_poly_fit_length(Q, FLINT_MAX(q1->length+n1, q2->length));
   fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs, q2->limbs), 1));
   _fmpz_poly_left_shift(Q, q1, n1);
   fmpz_poly_clear(q1);
   _fmpz_poly_add(Q, Q, q2);
   fmpz_poly_clear(q2);
   
}

/* 
   Computes the 2norm of H and divides by the leading coefficient and returns
   an upper bound on the number of bits of the result
	Used for computing the Landau-Mignotte bound
*/
unsigned long fmpz_poly_2norm_bits_normalised(const fmpz_poly_t H)
{
   fmpz_t norm = fmpz_init(H->limbs + 1);
	fmpz_poly_2norm(norm, H);
	unsigned long bits = fmpz_bits(norm);
	fmpz_clear(norm);
   unsigned long bits_lc = fmpz_bits(_fmpz_poly_lead(H));
   return bits - bits_lc + 1;
}

int fmpz_poly_divides_modular(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, const ulong bound_in)
{
	if (B->length == 0)
   {
		printf("FLINT Exception: divide by zero!\n");
		abort();
   }

	if (A->length == 0)
	{
      fmpz_poly_zero(Q);
		return 1;
	}

	if (A->length < B->length)
	   return 0;

	ulong size_B = B->limbs + 1;
	ulong size_A = A->limbs + 1;
	fmpz_t coeff_B = B->coeffs;
	fmpz_t coeff_A = A->coeffs;

	// special case, check if B has trailing zeroes that A doesn't 
	// or B has all terms zero except leading term
	long i;
	for (i = 0; i < B->length - 1; i++)
	{
		if (!fmpz_is_zero(coeff_B)) break;
		if (!fmpz_is_zero(coeff_A)) return 0;
		coeff_A += size_A;
		coeff_B += size_B;
	}

	if (i == B->length - 1)
	   return fmpz_poly_divides_divconquer(Q, A, B);
   
   int divides = 0;
	
	fmpz_t lead_A = _fmpz_poly_lead(A);
   fmpz_t lead_B = _fmpz_poly_lead(B);

	fmpz_t q_lead = fmpz_init(FLINT_MAX((long)(FLINT_ABS(lead_A[0]) - FLINT_ABS(lead_B[0]) + 1), 0L));
   if (!fmpz_divides(q_lead, lead_A, lead_B)) 
	{
		fmpz_clear(q_lead);
		return 0;
	}
   
	int q_one = fmpz_is_one(q_lead);
	
   unsigned long p;
   unsigned long pbits;

	unsigned long bits1 = FLINT_ABS(fmpz_poly_max_bits(A));
   unsigned long bits2 = FLINT_ABS(fmpz_poly_max_bits(B));
   
   long bits = bits1 - bits2;
	if (bits < 2L) bits = 2;
	
#if FLINT_BITS == 64
	if (bits == FLINT_D_BITS - 1) 
   {
      pbits = FLINT_D_BITS;
      p = (1L<<FLINT_D_BITS)-200;
   } else if (bits < FLINT_D_BITS) 
   {
      pbits = bits;
      p = (1L<<bits);
   } else 
   {
#endif
		pbits = FLINT_BITS - 1;
      p = (1L<<(FLINT_BITS-2));
#if FLINT_BITS == 64
   }
#endif
   
   zmod_poly_t a, b, q, r;
   
   int first = 1;

   unsigned long n = B->length;

   unsigned long modsize = FLINT_MAX(A->limbs, B->limbs);
   fmpz_t modulus = fmpz_init(modsize);
   modulus[0] = 0L;
   
	fmpz_poly_t P;
	fmpz_poly_init(P);

	ulong bound;
	
	if (bound_in) bound = bound_in;
	else
	{
		/*
		   Bound is easily derived from Thm 5.3 of 
			http://compalg.inf.elte.hu/~tony/Informatikai-Konyvtar/03-Algorithms%20of%20Informatics%201,%202,%203/CompAlg29May.pdf
			The + 1 is to allow for signed coefficients after Chinese Remaindering
		*/
	   ulong na = fmpz_poly_2norm_bits_normalised(A);
      bound = A->length - B->length + na + fmpz_bits(lead_A) + 1;
	}

   for (;;)
   {
      if (!first)
      {
         zmod_poly_clear(a);
         zmod_poly_clear(b);
         zmod_poly_clear(q); 
         zmod_poly_clear(r); 
      } else first = 0;
      do { p = z_nextprime(p); }
      while ((!fmpz_mod_ui(lead_A, p)) || (!fmpz_mod_ui(lead_B, p))); 

		zmod_poly_init(a, p);
      zmod_poly_init(b, p);
      zmod_poly_init(q, p);
      zmod_poly_init(r, p);
      
      if (bits1 + 1 < pbits) fmpz_poly_to_zmod_poly_no_red(a, A);
      else fmpz_poly_to_zmod_poly(a, A);
      if (bits2 + 1 < pbits) fmpz_poly_to_zmod_poly_no_red(b, B);
      else fmpz_poly_to_zmod_poly(b, B);
      
      zmod_poly_divrem(q, r, a, b);

		if (r->length != 0) // doesn't divide
		   break;

      if (q_one) zmod_poly_make_monic(q, q);
      else
      {
         unsigned long q_inv = z_invert(q->coeffs[q->length-1], q->p);
         unsigned long q_mod = fmpz_mod_ui(q_lead, q->p);
         q_inv = z_mulmod2_precomp(q_inv, q_mod, q->p, q->p_inv);
         zmod_poly_scalar_mul(q, q, q_inv);
      }
         
		if (modulus[0] == 0) 
		{
			zmod_poly_to_fmpz_poly(Q, q);
		
		   if (bits <= pbits)
			{
				fmpz_poly_mul(P, Q, B);
				if (fmpz_poly_equal(P, A))
				{
					divides = 1;
					break;
				} else if (pbits > bound)
				{
				   divides = 0;
					break;
				}
			}
				
			fmpz_set_ui(modulus, p);
         continue;
      }
      
      fmpz_t newmod = fmpz_init(modulus[0] + 1);
      
      ulong n_bits;
		if ((fmpz_poly_CRT(Q, Q, q, newmod, modulus)) || (bits <= (n_bits = fmpz_bits(newmod))) || (bound < n_bits))
		{
			fmpz_poly_mul(P, Q, B);
			if (fmpz_poly_equal(P, A))
			{
				divides = 1;
				fmpz_clear(newmod);
				break;
			} else if (bound < n_bits)
			{  
				divides = 0;
			   break;
			}
      } 
      
      if (newmod[0] >= modsize) 
      {
         modulus = fmpz_realloc(modulus, (modsize+8));
         modsize += 8;
      }
      
      fmpz_set(modulus, newmod);
      fmpz_clear(newmod); // release newmod
   }

	
   fmpz_clear(q_lead); // release g
   
	zmod_poly_clear(a);
   zmod_poly_clear(b);
   zmod_poly_clear(q); 
   zmod_poly_clear(r); 
        
   fmpz_poly_clear(P);
   fmpz_clear(modulus);
   
	return divides;
}

/*===================================================================================

   Pseudo-division algorithm
   
====================================================================================*/

//Pseudo-division a la Cohen. 

void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeffs_R;
   
   long m = A->length;
   long n = B->length;
   long e = m - n + 1;
   
   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   unsigned long bits_B_lead = fmpz_bits(B_lead);
   
   mp_limb_t sign_quot;
      
   fmpz_poly_fit_length(R, A->length);
   fmpz_poly_fit_limbs(R, A->limbs);  
   R->length = A->length;
   _fmpz_poly_set(R, A);
   coeffs_R = R->coeffs;
   size_R = R->limbs+1;
   
   if ((long) R->length >= (long) B->length)
   {
      fmpz_poly_fit_length(Q, R->length-B->length+1);
      fmpz_poly_fit_limbs(Q, ABS(A->coeffs[(A->length-1)*(A->limbs+1)]));
      for (unsigned long i = 0; i < R->length-B->length+1; i++) Q->coeffs[i*(Q->limbs+1)] = 0;
      Q->length = R->length-B->length+1;
      size_Q = Q->limbs+1;
   } else 
   {
      _fmpz_poly_zero(Q);
      return;
   }
   fmpz_poly_t Bm1;
   Bm1->length = B->length - 1;
   Bm1->limbs = B->limbs;
   Bm1->coeffs = B->coeffs;
   
   while ((long) R->length >= (long) B->length)
   {
      _fmpz_poly_scalar_mul_fmpz(Q, Q, B_lead);
      coeff_Q = coeffs_R + (R->length-1)*size_R;
      fmpz_add(Q->coeffs + (R->length-B->length)*size_Q, Q->coeffs + (R->length-B->length)*size_Q, coeff_Q);
      
      if (B->length > 1)
      {
         fmpz_poly_init2(qB, B->length-1, B->limbs+ABS(coeff_Q[0]));
         _fmpz_poly_scalar_mul_fmpz(qB, Bm1, coeff_Q); 
         fmpz_poly_fit_limbs(R, FLINT_MAX(R->limbs + size_B_lead, qB->limbs) + 1);
      } else
      {
         fmpz_poly_fit_limbs(R, R->limbs + size_B_lead);
      }  
      
      coeffs_R = R->coeffs;
      size_R = R->limbs+1;
      _fmpz_poly_scalar_mul_fmpz(R, R, B_lead);
      
      fmpz_poly_t R_sub;
      R_sub->coeffs = coeffs_R+(R->length-B->length)*size_R;
      R_sub->limbs = R->limbs;
      R_sub->length = B->length-1;
      if (B->length > 1)
      {
         _fmpz_poly_sub(R_sub, R_sub, qB);
         
         fmpz_poly_clear(qB);
      }
      R_sub->coeffs[(B->length-1)*(R_sub->limbs+1)] = 0;
      _fmpz_poly_normalise(R);

      if (R->length) fmpz_poly_fit_limbs(Q, FLINT_MAX(Q->limbs + size_B_lead, ABS(R->coeffs[(R->length-1)*(R->limbs+1)])) + 1);
      size_Q = Q->limbs+1;
      e--;
   }
   fmpz_t pow = (fmpz_t) flint_stack_alloc((bits_B_lead*e)/FLINT_BITS+2);
   fmpz_pow_ui(pow, B_lead, e);
   fmpz_poly_fit_limbs(Q, Q->limbs+ABS(pow[0]));
   fmpz_poly_fit_limbs(R, R->limbs+ABS(pow[0]));
   _fmpz_poly_scalar_mul_fmpz(Q, Q, pow);
   _fmpz_poly_scalar_mul_fmpz(R, R, pow);   
   flint_stack_release();
}

void fmpz_poly_pseudo_divrem_shoup(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeffs_R;
   
   long m = A->length;
   long n = B->length;
   long e = m - n + 1;
   
   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   mp_limb_t sign_quot;

   int lead_is_one = fmpz_is_one(B_lead);
   
   fmpz_t coeff_R, coeff_A;
   if ((long) A->length >= (long) B->length)
   {
      fmpz_poly_fit_length(R, A->length);
      fmpz_poly_fit_limbs(R, A->limbs + (m - n)*size_B_lead);  
      R->length = A->length;
      coeffs_R = R->coeffs;
      size_R = R->limbs + 1;
      fmpz_poly_fit_length(Q, A->length-B->length+1);
      fmpz_poly_fit_limbs(Q, ABS(A->coeffs[(A->length-1)*(A->limbs+1)]));
      for (unsigned long i = 0; i < A->length-B->length+1; i++) Q->coeffs[i*(Q->limbs+1)] = 0;
      Q->length = A->length-B->length+1;
      size_Q = Q->limbs+1;
   } else 
   {
      fmpz_poly_set(R, A);
      _fmpz_poly_zero(Q);
      return;
   }
   
   if (!lead_is_one)
   {
      fmpz_t pow = (fmpz_t) flint_stack_alloc(size_B_lead*(m-n+1)+1);
      fmpz_set(pow, B_lead);
      coeff_R = coeffs_R + (m-n-1)*size_R; 
      coeff_A = A->coeffs + (m-n-1)*(A->limbs+1);
      
      for (long i = m - n - 1; i >= 0; i--)
      {
         fmpz_mul(coeff_R, coeff_A, pow);
         if (i > 0) fmpz_mul(pow, pow, B_lead);
         coeff_R -= size_R;
         coeff_A -= (A->limbs+1);
      }
      coeff_R = coeffs_R + (m-n)*size_R;
      coeff_A = A->coeffs + (m-n)*(A->limbs+1);
      for (long i = m - n; i < m; i++)
      {
         fmpz_set(coeff_R, coeff_A);
         coeff_R += size_R;
         coeff_A += (A->limbs+1);
      }
      
      flint_stack_release();
      
   } else
   {
      coeff_R = coeffs_R;
      coeff_A = A->coeffs;
      for (unsigned long i = 0; i < m; i++)
      {
         fmpz_set(coeff_R, coeff_A);
         coeff_R += size_R;
         coeff_A += (A->limbs+1);
      }
   } 
   unsigned long coeff = R->length;
   
   fmpz_poly_t Bm1;
   Bm1->length = B->length - 1;
   Bm1->limbs = B->limbs;
   Bm1->coeffs = B->coeffs;
   
   while ((long) coeff >= (long) B->length)
   {
      coeff_Q = coeffs_R + (coeff-1)*size_R;
      fmpz_set(Q->coeffs + (coeff-B->length)*size_Q, coeff_Q);
      coeff_Q = Q->coeffs + (coeff-B->length)*size_Q;
      
      if (B->length > 1) 
      {
         fmpz_poly_init2(qB, Bm1->length, Bm1->limbs+ABS(coeff_Q[0]));
         _fmpz_poly_scalar_mul_fmpz(qB, Bm1, coeff_Q); 
         
         fmpz_poly_fit_limbs(R, FLINT_MAX(R->limbs + size_B_lead, qB->limbs) + 1);
         coeffs_R = R->coeffs;
         size_R = R->limbs+1;
      
         fmpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R+(coeff-B->length)*size_R;
         R_sub->limbs = R->limbs;
         R_sub->length = B->length-1;
         if (!lead_is_one) _fmpz_poly_scalar_mul_fmpz(R_sub, R_sub, B_lead);
         _fmpz_poly_sub(R_sub, R_sub, qB);
      
         fmpz_poly_clear(qB);
      }   
      
      coeff--;
      if (coeff >= B->length) 
      {
         fmpz_poly_fit_limbs(Q, ABS(R->coeffs[(coeff-1)*size_R]));
         size_Q = Q->limbs+1;
      }
   }
   
   R->length = B->length - 1;
   _fmpz_poly_normalise(R);
   if (!lead_is_one)
   {
      size_Q = Q->limbs + 1;
      coeff_Q = Q->coeffs; 
      fmpz_t pow = (fmpz_t) flint_stack_alloc(size_B_lead*(m-n+1)+1);
      fmpz_set(pow, B_lead);
      
      for (long i = 1; i <= m-n; i++)
      {
         coeff_Q = Q->coeffs + i*size_Q;
         fmpz_poly_fit_limbs(Q, ABS(coeff_Q[0]) + ABS(pow[0]));
         size_Q = Q->limbs + 1;
         coeff_Q = Q->coeffs + i*size_Q;
         fmpz_mul(coeff_Q, coeff_Q, pow);
         if (i < m-n) fmpz_mul(pow, pow, B_lead);
      }  
       
      flint_stack_release();
   }
}

/*
   Pseudo division of A by B. Returns Q, R and d such that l^d A = QB + R where l 
   is the leading coefficient of B. This is faster than the pseudo divisions above
   when there is no coefficient explosion, but is slower otherwise. It is usually
   desirable to use this version unless you know specifically that coefficient
   explosion will occur.
*/

void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, 
                               unsigned long * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB;
   
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeff_R;
   fmpz_t coeffs_R;
   int scale;
   
   long m = A->length;
   long n = B->length;
   
   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   mp_limb_t sign_quot;
      
   fmpz_poly_fit_length(R, A->length);
   fmpz_poly_fit_limbs(R, A->limbs);  
   R->length = A->length;
   _fmpz_poly_set(R, A);
   coeffs_R = R->coeffs;
   size_R = R->limbs+1;
   
   *d = 0;
   
   if ((long) R->length >= (long) B->length)
   {
      fmpz_poly_fit_length(Q, R->length-B->length+1);
      fmpz_poly_fit_limbs(Q, ABS(A->coeffs[(A->length-1)*(A->limbs+1)]));
      for (unsigned long i = 0; i < R->length-B->length+1; i++) Q->coeffs[i*(Q->limbs+1)] = 0;
      Q->length = R->length-B->length+1;
      size_Q = Q->limbs+1;
   } else 
   {
      _fmpz_poly_zero(Q);
      return;
   }
   
   fmpz_poly_t Bm1;
   Bm1->length = B->length - 1;
   Bm1->limbs = B->limbs;
   Bm1->coeffs = B->coeffs;
   
   coeff_R = coeffs_R + (R->length-1)*size_R;
   
   fmpz_t rem = (fmpz_t) flint_heap_alloc(size_B_lead+1);
      
   while ((long) R->length >= (long) B->length)
   {
      coeff_Q = Q->coeffs+(R->length - B->length)*size_Q;
          
      __fmpz_normalise(coeff_R);
      sign_quot = ABS(coeff_R[0]) - size_B_lead + 1;
         
      if (((long) sign_quot > 1) || ((sign_quot == 1) && (mpn_cmp(coeff_R+1, B_lead+1, size_B_lead) >= 0)))
      {
         mpn_tdiv_qr(coeff_Q+1, rem+1, 0, coeff_R+1, ABS(coeff_R[0]), B_lead+1, size_B_lead);
      
         rem[0] = size_B_lead;
         
         __fmpz_normalise(rem);
      } else
      {
         coeff_Q[0] = 0;
         if (coeff_R[0] == 0) rem[0] = 0;
         else 
         {
            rem[0] = 1;
            rem[1] = 1;
         }
      }
      if (fmpz_is_zero(rem))
      {
         if (((long) (sign_B_lead ^ coeff_R[0])) < 0) 
         {
            coeff_Q[0] = -sign_quot;
            for (unsigned long i = 0; i < size_B_lead; i++)
            {
               if (rem[i])
               {
                  fmpz_sub_ui_inplace(coeff_Q,1);
                  break;
               }
            }
         }
         else 
         {
            coeff_Q[0] = sign_quot;
         }
         NORM(coeff_Q);  
         scale = 0; 
      } else 
      {   
         _fmpz_poly_scalar_mul_fmpz(Q, Q, B_lead);
         fmpz_set(coeff_Q, coeff_R);
         scale = 1;
         (*d)++;
      }
           
      if (B->length > 1)
      {
         fmpz_poly_init2(qB, B->length-1, B->limbs+ABS(coeff_Q[0]));
         _fmpz_poly_scalar_mul_fmpz(qB, Bm1, coeff_Q); 
      }   
      
      if (scale)
      {
         fmpz_poly_fit_limbs(R, FLINT_MAX(R->limbs + size_B_lead, qB->limbs) + 1);
         coeffs_R = R->coeffs;
         size_R = R->limbs+1;
         _fmpz_poly_scalar_mul_fmpz(R, R, B_lead);
      } else if (B->length > 1)
      {
         fmpz_poly_fit_limbs(R, FLINT_MAX(R->limbs, qB->limbs) + 1);
         coeffs_R = R->coeffs;
         size_R = R->limbs+1;
      }
      
      fmpz_poly_t R_sub;
      R_sub->coeffs = coeffs_R+(R->length-B->length)*size_R;
      R_sub->limbs = R->limbs;
      R_sub->length = B->length-1;
      
      if (B->length > 1)
      {
         _fmpz_poly_sub(R_sub, R_sub, qB);
         
         fmpz_poly_clear(qB);
      }
      
      R_sub->coeffs[(B->length-1)*(R_sub->limbs+1)] = 0;
      
      _fmpz_poly_normalise(R);
      coeff_R = coeffs_R + (R->length-1)*size_R;
      
      if (R->length) fmpz_poly_fit_limbs(Q, R->limbs);
      size_Q = Q->limbs+1;
   }
  
   flint_heap_free(rem);
}

void fmpz_poly_pseudo_div_basecase(fmpz_poly_t Q, unsigned long * d, 
                                     const fmpz_poly_t A, const fmpz_poly_t B)
{
   fmpz_poly_t qB, R;
   
   fmpz_poly_init(R);
   
   if (B->length == 0)
   {
      printf("Error: Divide by zero\n");
      abort();
   }
   
   unsigned long size_A = A->limbs + 1;
   unsigned long size_B = B->limbs + 1;
   unsigned long size_R;
   unsigned long limbs_R;
   unsigned long size_Q;
   unsigned long limbs_Q;
   fmpz_t coeffs_A = A->coeffs;
   fmpz_t coeffs_B = B->coeffs;
   fmpz_t B_lead = coeffs_B + (B->length-1)*size_B; 
   fmpz_t coeff_Q;
   fmpz_t coeff_R;
   fmpz_t coeffs_R;
   int scale;
   
   long m = A->length;
   long n = B->length;
   
   NORM(B_lead);
   
   unsigned long size_B_lead = ABS(B_lead[0]);
   mp_limb_t sign_B_lead = B_lead[0];
   mp_limb_t sign_quot;
      
   fmpz_poly_fit_length(R, A->length);
   fmpz_poly_fit_limbs(R, A->limbs);  
   R->length = A->length;
   _fmpz_poly_set(R, A);
   coeffs_R = R->coeffs;
   size_R = R->limbs+1;
   
   *d = 0;
   
   if ((long) R->length >= (long) B->length)
   {
      fmpz_poly_fit_length(Q, R->length-B->length+1);
      fmpz_poly_fit_limbs(Q, ABS(A->coeffs[(A->length-1)*(A->limbs+1)]));
      for (unsigned long i = 0; i < R->length-B->length+1; i++) Q->coeffs[i*(Q->limbs+1)] = 0;
      Q->length = R->length-B->length+1;
      size_Q = Q->limbs+1;
   } else 
   {
      _fmpz_poly_zero(Q);
      return;
   }
   
   fmpz_poly_t Bm1;
   Bm1->length = B->length - 1;
   Bm1->limbs = B->limbs;
   Bm1->coeffs = B->coeffs;
   
   coeff_R = coeffs_R + (R->length-1)*size_R;
   
   fmpz_t rem = (fmpz_t) flint_heap_alloc(size_B_lead+1);
      
   while ((long) R->length >= (long) B->length)
   {
      coeff_Q = Q->coeffs+(R->length - B->length)*size_Q;
          
      __fmpz_normalise(coeff_R);
      sign_quot = ABS(coeff_R[0]) - size_B_lead + 1;
         
      if (((long) sign_quot > 1) || ((sign_quot == 1) && (mpn_cmp(coeff_R+1, B_lead+1, size_B_lead) >= 0)))
      {
         mpn_tdiv_qr(coeff_Q+1, rem+1, 0, coeff_R+1, ABS(coeff_R[0]), B_lead+1, size_B_lead);
      
         rem[0] = size_B_lead;
         
         __fmpz_normalise(rem);
      } else
      {
         coeff_Q[0] = 0;
         if (coeff_R[0] == 0) rem[0] = 0;
         else 
         {
            rem[0] = 1;
            rem[1] = 1;
         }
      }
      if (fmpz_is_zero(rem))
      {
         if (((long) (sign_B_lead ^ coeff_R[0])) < 0) 
         {
            coeff_Q[0] = -sign_quot;
            for (unsigned long i = 0; i < size_B_lead; i++)
            {
               if (rem[i])
               {
                  fmpz_sub_ui_inplace(coeff_Q,1);
                  break;
               }
            }
         }
         else 
         {
            coeff_Q[0] = sign_quot;
         }
         NORM(coeff_Q);  
         scale = 0; 
      } else 
      {   
         _fmpz_poly_scalar_mul_fmpz(Q, Q, B_lead);
         fmpz_set(coeff_Q, coeff_R);
         scale = 1;
         (*d)++;
      }
               
      if (R->length != B->length)
      {
         if (B->length > 1)
         {
            fmpz_poly_init2(qB, B->length-1, B->limbs+ABS(coeff_Q[0]));
            _fmpz_poly_scalar_mul_fmpz(qB, Bm1, coeff_Q); 
         }   
      
         if (scale)
         {
            fmpz_poly_fit_limbs(R, FLINT_MAX(R->limbs + size_B_lead, qB->limbs) + 1);
            coeffs_R = R->coeffs;
            size_R = R->limbs+1;
            _fmpz_poly_scalar_mul_fmpz(R, R, B_lead);
         }
      
         fmpz_poly_t R_sub;
         R_sub->coeffs = coeffs_R+(R->length-B->length)*size_R;
         R_sub->limbs = R->limbs;
         R_sub->length = B->length-1;
         if (B->length > 1)
         {
            _fmpz_poly_sub(R_sub, R_sub, qB);
         
            fmpz_poly_clear(qB);
         }
         R_sub->coeffs[(B->length-1)*(R_sub->limbs+1)] = 0;
      
         _fmpz_poly_normalise(R);
         coeff_R = coeffs_R + (R->length-1)*size_R;
      
         if (R->length) fmpz_poly_fit_limbs(Q, R->limbs);
         size_Q = Q->limbs+1;
      } else R->length = 0;
   }
   
   fmpz_poly_clear(R);
   flint_heap_free(rem);
}

/* 
   Pseudo division using a divide and conquer algorithm.
*/
 
void fmpz_poly_pseudo_divrem_recursive(fmpz_poly_t Q, fmpz_poly_t R, unsigned long * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      fmpz_poly_fit_length(R, A->length);
      fmpz_poly_fit_limbs(R, A->limbs);  
      _fmpz_poly_set(R, A);
      _fmpz_poly_zero(Q);
      *d = 0;
      
      return;
   }
   
   unsigned long crossover = 16;
   unsigned long crossover2 = 128;
   
   if (B->limbs > 16) crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      fmpz_poly_pseudo_divrem_basecase(Q, R, d, A, B);
      
      return;
   }
   fmpz_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, r1, d2q1, d2q2, r2, t, u, temp;
   
   fmpz_t B_lead;
   unsigned long size_B_lead;
   unsigned long bits_B_lead;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _fmpz_poly_attach_shift(d1, B, n2);
   _fmpz_poly_attach_truncate(d2, B, n2);
   _fmpz_poly_attach_shift(d3, B, n1);
   _fmpz_poly_attach_truncate(d4, B, n1);
   
   /* We need the leading coefficient of B */
   
   B_lead = B->coeffs + (B->length-1)*(B->limbs+1);
   size_B_lead = ABS(B_lead[0]);
   bits_B_lead = fmpz_bits(B_lead);
      
   if (A->length <= n2 + B->length - 1)
   {
      /*
         A is greater than length n1+n2-1 and at most 
         length n1+2*n2-1
         We shift right by n1 and zero the last n2-1
         coefficients, leaving at at most n2 significant
         terms
      */
      
      _fmpz_poly_stack_init(p1, A->length-n1, A->limbs);
      _fmpz_poly_right_shift(p1, A, n1);
      _fmpz_poly_zero_coeffs(p1, n2-1);
      
      /* 
         We compute p1 div d3 which is at most 
         a 2*n2-1 by n2 division, leaving n2 terms 
         in the quotient. Since we are doing pseudo
         division, the remainder will have at most
         n2-1 terms
      */
      
      fmpz_poly_init(r1);
      fmpz_poly_pseudo_divrem_recursive(Q, r1, d, p1, d3); 
      _fmpz_poly_stack_clear(p1);

      /*
         We compute d2q1 = Q*d4
         It will have at most n1+n2-1 terms
      */
      
      _fmpz_poly_stack_init(d2q1, d4->length+Q->length-1, d4->limbs+Q->limbs+1); 
      _fmpz_poly_mul(d2q1, d4, Q);
      
      /*
         Compute R = lead(B)^n * R' where R' is 
         the terms of A we haven't dealt with, 
         of which there are at most n1+n2-1
      */
      
      fmpz_poly_fit_length(R, n1+n2-1);
      fmpz_poly_fit_limbs(R, FLINT_MAX(FLINT_MAX(A->limbs+((*d)*bits_B_lead)/FLINT_BITS+1, r1->limbs), d2q1->limbs)+1);
      fmpz_t pow = (fmpz_t) flint_stack_alloc((bits_B_lead*(*d))/FLINT_BITS+2);
      fmpz_pow_ui(pow, B_lead, *d);
      _fmpz_poly_attach_truncate(temp, A, n1+n2-1);
      _fmpz_poly_scalar_mul_fmpz(R, temp, pow);
      flint_stack_release();
      
      /*
         Compute the original remainder from the
         first pseudo division r' = r1^n1 - d2q1.
         This should be thought of as r'/lead(B)^n
         We add this to the remainder R', first
         multiplying everything through by 
         lead(B)^n. This gives the remainder 
         R + r'. We note r' will have at most 
         n1+n2-1 terms.
      */
      
      fmpz_poly_fit_length(r1, FLINT_MAX(r1->length+n1, d2q1->length));
      _fmpz_poly_left_shift(r1, r1, n1); 
      _fmpz_poly_sub(r1, r1, d2q1);
      _fmpz_poly_stack_clear(d2q1);
      _fmpz_poly_add(R, R, r1);
      fmpz_poly_clear(r1);
      
      return;   
   } 
   
   unsigned long s1, s2;
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length - 1
      // We call this polynomial p1. Zero the final B->length-1
      // coefficients. Note A->length > 2*B->length - 1
      
      unsigned long shift = A->length - 2*B->length + 1;
      _fmpz_poly_stack_init(p1, 2*B->length - 1, A->limbs);
      _fmpz_poly_right_shift(p1, A, shift);
      _fmpz_poly_zero_coeffs(p1, B->length - 1);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         r1 is length at most B->length-1
      */
      
      fmpz_poly_init(r1);
      fmpz_poly_init(q1);
      
      fmpz_poly_pseudo_divrem_recursive(q1, r1, &s1, p1, B); 
      _fmpz_poly_stack_clear(p1);
       
      /* 
         Compute t = (lead(B)^s1) * a2 + r1*x^shift
         which ends up being at most length A->length - B->length 
         since r1 is at most length B->length-1 
         Here a2 is what remains of A after the first R->length
         coefficients are removed.
      */  
   
      _fmpz_poly_stack_init(t, A->length - B->length, FLINT_MAX(A->limbs+(bits_B_lead*s1)/FLINT_BITS+1, r1->limbs)+1);
      _fmpz_poly_attach_truncate(temp, A, A->length - B->length);
      
      fmpz_t pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s1)/FLINT_BITS+2);
      fmpz_pow_ui(pow, B_lead, s1);
      _fmpz_poly_scalar_mul_fmpz(t, temp, pow);
      flint_stack_release();
   
      fmpz_poly_fit_length(r1, r1->length+shift);
      _fmpz_poly_left_shift(r1, r1, shift);
      _fmpz_poly_add(t, t, r1);
      fmpz_poly_clear(r1);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length - B->length
         r2 has length at most B->length - 1
      */
   
      fmpz_poly_init(q2);
      fmpz_poly_pseudo_divrem_recursive(q2, R, &s2, t, B); 
      _fmpz_poly_stack_clear(t);  
      
      /*
         Write out Q = lead(B)^s2*q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
         q1 cannot have length zero since we
         are doing pseudo division
      */
   
      fmpz_poly_fit_length(Q, q1->length+shift);
      fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs + (s2*bits_B_lead)/FLINT_BITS+1, q2->limbs), 1));
   
      pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s2)/FLINT_BITS+2);
      fmpz_pow_ui(pow, B_lead, s2);
      _fmpz_poly_scalar_mul_fmpz(Q, q1, pow);
      fmpz_poly_clear(q1);
      flint_stack_release();
      _fmpz_poly_left_shift(Q, Q, shift);
      _fmpz_poly_add(Q, Q, q2);
      fmpz_poly_clear(q2);
   
      /* 
         Set d to the power of lead(B) that everything
         must be multiplied by
      */
      
      *d = s1 + s2;
      
      return;
   } 
   
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is at most length n1 and a2 is length n2 
      and a3 is length n1+n2-1 
      We set p1 = a1*x^(n1-1), so it has length at most
      2*n1-1. We note A is at least length n1+2*n2-1
   */
      
   _fmpz_poly_stack_init(p1, A->length-2*n2, A->limbs);
   _fmpz_poly_right_shift(p1, A, 2*n2);
   _fmpz_poly_zero_coeffs(p1, n1-1);
   
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length length n1
      r1 is length n1-1
   */
   
   fmpz_poly_init(r1);
   fmpz_poly_init(q1);
   fmpz_poly_pseudo_divrem_recursive(q1, r1, &s1, p1, d1); 
   _fmpz_poly_stack_clear(p1);   
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length n1+n2-1
      Note q1->length is at least 1 since we are doing 
      pseudo division
   */  
   
   _fmpz_poly_stack_init(d2q1, d2->length+q1->length-1, d2->limbs+q1->limbs+1); 
   _fmpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute t = (lead(B)^s1) * (a2*x^(n1+n2-1)+a3) 
                               + r1*x^(2*n2) - d2q1*x^n2
      which ends up being at most length n2+B->length-1 
      since r1 is at most length n1-1 and d2q1 is at 
      most length n1+n2-1
   */  
   
   _fmpz_poly_stack_init(t, n2+B->length-1, FLINT_MAX(FLINT_MAX(A->limbs+(bits_B_lead*s1)/FLINT_BITS+1, r1->limbs), d2q1->limbs)+1);
   _fmpz_poly_attach_truncate(temp, A, n2+B->length-1);
   fmpz_t pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s1)/FLINT_BITS+2);
   fmpz_pow_ui(pow, B_lead, s1);
   _fmpz_poly_scalar_mul_fmpz(t, temp, pow);
   flint_stack_release();
   
   fmpz_poly_fit_length(r1, FLINT_MAX(r1->length+2*n2, d2q1->length+n2));
   _fmpz_poly_left_shift(r1, r1, n2);
   _fmpz_poly_sub(r1, r1, d2q1);
   _fmpz_poly_left_shift(r1, r1, n2);
   _fmpz_poly_add(t, t, r1);
   fmpz_poly_clear(r1);
   
   /*
      Compute q2 = t div B and set R to the remainder
      It is at most a n2+B->length-1 by n1+n2 division, 
      so the length of q2 will be at most n2 .
      R will have length at most n1+n2-1 since we are
      doing pseudo division
   */
   
   fmpz_poly_init(q2);
   fmpz_poly_pseudo_divrem_recursive(q2, R, &s2, t, B); 
   _fmpz_poly_stack_clear(t);
   _fmpz_poly_stack_clear(d2q1);
      
   /*
      Write out Q = lead(B)^s2 * q1*x^n2 + q2
      Q has length n1+n2
      Note q1->length is not zero since we are doing
      pseudo division
   */
   
   fmpz_poly_fit_length(Q, q1->length+n2);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs + (s2*bits_B_lead)/FLINT_BITS+1, q2->limbs), 1));
   pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s2)/FLINT_BITS+2);
   fmpz_pow_ui(pow, B_lead, s2);
   _fmpz_poly_scalar_mul_fmpz(Q, q1, pow);
   fmpz_poly_clear(q1);
   flint_stack_release();
   _fmpz_poly_left_shift(Q, Q, n2);
   _fmpz_poly_add(Q, Q, q2);
   fmpz_poly_clear(q2);
   
   /* 
      Set d to the power of lead(B) which everything 
      has been raised to
   */
   
   *d = s1+s2;  
}

void fmpz_poly_pseudo_div_recursive(fmpz_poly_t Q, unsigned long * d, const fmpz_poly_t A, const fmpz_poly_t B)
{
   if (A->length < B->length)
   {
      _fmpz_poly_zero(Q);
      *d = 0;
      
      return;
   }
   
   unsigned long crossover = 16;
   unsigned long crossover2 = 256;
   
   if (B->limbs > 16) crossover = 8;
   if ((B->length <= 12) && (B->limbs > 8)) crossover = 8;
   
   if ((B->length <= crossover) 
   || ((A->length > 2*B->length - 1) && (A->length < crossover2)))
   {
      fmpz_poly_pseudo_div_basecase(Q, d, A, B);
      
      return;
   }
   
   fmpz_poly_t d1, d2, d3, d4, p1, q1, q2, dq1, dq2, r1, d2q1, d2q2, r2, t, u, temp;
   
   fmpz_t B_lead;
   unsigned long size_B_lead;
   unsigned long bits_B_lead;
   
   unsigned long n1 = (B->length+1)/2;
   unsigned long n2 = B->length - n1;
   
   /* We let B = d1*x^n2 + d2 */
   
   _fmpz_poly_attach_shift(d1, B, n2);
   _fmpz_poly_attach_truncate(d2, B, n2);
   _fmpz_poly_attach_shift(d3, B, n1);
   _fmpz_poly_attach_truncate(d4, B, n1);
   
   /* We need the leading coefficient of B */
   
   B_lead = B->coeffs + (B->length-1)*(B->limbs+1);
   size_B_lead = ABS(B_lead[0]);
   bits_B_lead = fmpz_bits(B_lead);
      
   if (A->length <= n2 + B->length - 1)
   {
      /*
         A is greater than length n1+n2-1 and at most 
         length n1+2*n2-1
         We shift right by n1 and zero the last n2-1
         coefficients, leaving at at most n2 significant
         terms
      */
      
      _fmpz_poly_stack_init(p1, A->length-n1, A->limbs);
      _fmpz_poly_right_shift(p1, A, n1);
      _fmpz_poly_zero_coeffs(p1, n2-1);
      
      /* 
         We compute p1 div d3 which is at most 
         a 2*n2-1 by n2 division, leaving n2 terms 
         in the quotient. 
      */
      
      fmpz_poly_pseudo_div_recursive(Q, d, p1, d3); 
      _fmpz_poly_stack_clear(p1);
      
      return;   
   } 
   
   unsigned long s1, s2;
   
   if (A->length > 2*B->length - 1)
   {
      // We shift A right until it is length 2*B->length - 1
      // We call this polynomial p1. Zero the final B->length-1
      // coefficients. Note A->length > 2*B->length - 1
      unsigned long shift = A->length - 2*B->length + 1;
      _fmpz_poly_stack_init(p1, 2*B->length - 1, A->limbs);
      _fmpz_poly_right_shift(p1, A, shift);
      _fmpz_poly_zero_coeffs(p1, B->length - 1);
      
      /* 
         Set q1 to p1 div B 
         This is a 2*B->length-1 by B->length division so 
         q1 ends up being at most length B->length
         r1 is length at most B->length-1
      */
      
      fmpz_poly_init(r1);
      fmpz_poly_init(q1);
      
      fmpz_poly_pseudo_divrem_recursive(q1, r1, &s1, p1, B); 
      _fmpz_poly_stack_clear(p1);
       
      /* 
         Compute t = (lead(B)^s1) * a2 + r1*x^shift
         which ends up being at most length A->length - B->length 
         since r1 is at most length B->length-1 
         Here a2 is what remains of A after the first R->length
         coefficients are removed.
      */  
   
      _fmpz_poly_stack_init(t, A->length - B->length, FLINT_MAX(A->limbs+(bits_B_lead*s1)/FLINT_BITS+1, r1->limbs)+1);
      _fmpz_poly_attach_truncate(temp, A, A->length - B->length);
      
      fmpz_t pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s1)/FLINT_BITS+2);
      fmpz_pow_ui(pow, B_lead, s1);
      _fmpz_poly_scalar_mul_fmpz(t, temp, pow);
      flint_stack_release();
   
      fmpz_poly_fit_length(r1, r1->length+shift);
      _fmpz_poly_left_shift(r1, r1, shift);
      _fmpz_poly_add(t, t, r1);
      fmpz_poly_clear(r1);
      
      /*
         Compute q2 = t div B
         It is a smaller division than the original 
         since t->length <= A->length - B->length
      */
   
      fmpz_poly_init(q2);
      fmpz_poly_pseudo_div_recursive(q2, &s2, t, B); 
      _fmpz_poly_stack_clear(t);  
      
      /*
         Write out Q = lead(B)^s2*q1*x^shift + q2
         Q has length at most B->length+shift
         Note q2 has length at most shift since 
         at most it is an A->length-B->length 
         by B->length division
         q1 cannot have length zero since we
         are doing pseudo division
      */
   
      fmpz_poly_fit_length(Q, q1->length+shift);
      fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs + (s2*bits_B_lead)/FLINT_BITS+1, q2->limbs), 1));
   
      pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s2)/FLINT_BITS+2);
      fmpz_pow_ui(pow, B_lead, s2);
      _fmpz_poly_scalar_mul_fmpz(Q, q1, pow);
      flint_stack_release();
      fmpz_poly_clear(q1);
      _fmpz_poly_left_shift(Q, Q, shift);
      _fmpz_poly_add(Q, Q, q2);
      fmpz_poly_clear(q2);
   
      /* 
         Set d to the power of lead(B) that everything
         must be multiplied by
      */
      
      *d = s1 + s2;
      
      return;
   } 
   
   /* 
      We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
      where a1 is at most length n1 and a2 is length n2 
      and a3 is length n1+n2-1 
      We set p1 = a1*x^(n1-1), so it has length at most
      2*n1-1. We note A is at least length n1+2*n2-1
   */
      
   _fmpz_poly_stack_init(p1, A->length-2*n2, A->limbs);
   _fmpz_poly_right_shift(p1, A, 2*n2);
   _fmpz_poly_zero_coeffs(p1, n1-1);
   
   /* 
      Set q1 to p1 div d1 
      This is at most a 2*n1-1 by n1 division so 
      q1 ends up being at most length length n1
      r1 is length n1-1
   */
   
   fmpz_poly_init(r1);
   fmpz_poly_init(q1);
   fmpz_poly_pseudo_divrem_recursive(q1, r1, &s1, p1, d1); 
   _fmpz_poly_stack_clear(p1);   
   
   /* 
      Compute d2q1 = d2*q1 
      which ends up being length n1+n2-1
      Note q1->length is at least 1 since we are doing 
      pseudo division
   */  
   
   _fmpz_poly_stack_init(d2q1, d2->length+q1->length-1, d2->limbs+q1->limbs+1); 
   _fmpz_poly_mul(d2q1, d2, q1);
   
   /* 
      Compute t = (lead(B)^s1) * (a2*x^(n1+n2-1)+a3) 
                               + r1*x^(2*n2) - d2q1*x^n2
      which ends up being at most length n2+B->length-1 
      since r1 is at most length n1-1 and d2q1 is at 
      most length n1+n2-1
   */  
   
   _fmpz_poly_stack_init(t, n2+B->length-1, FLINT_MAX(FLINT_MAX(A->limbs+(bits_B_lead*s1)/FLINT_BITS+1, r1->limbs), d2q1->limbs)+1);
   _fmpz_poly_attach_truncate(temp, A, n2+B->length-1);
   fmpz_t pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s1)/FLINT_BITS+2);
   fmpz_pow_ui(pow, B_lead, s1);
   _fmpz_poly_scalar_mul_fmpz(t, temp, pow);
   flint_stack_release();
   
   fmpz_poly_fit_length(r1, FLINT_MAX(r1->length+2*n2, d2q1->length+n2));
   _fmpz_poly_left_shift(r1, r1, n2);
   _fmpz_poly_sub(r1, r1, d2q1);
   _fmpz_poly_left_shift(r1, r1, n2);
   _fmpz_poly_add(t, t, r1);
   fmpz_poly_clear(r1);
   
   /*
      Compute q2 = t div B and set R to the remainder
      It is at most a n2+B->length-1 by n1+n2 division, 
      so the length of q2 will be at most n2 .
      R will have length at most n1+n2-1 since we are
      doing pseudo division
   */
   
   fmpz_poly_init(q2);
   fmpz_poly_pseudo_div_recursive(q2, &s2, t, B); 
   _fmpz_poly_stack_clear(t);
   _fmpz_poly_stack_clear(d2q1);
      
   /*
      Write out Q = lead(B)^s2 * q1*x^n2 + q2
      Q has length n1+n2
      Note q1->length is not zero since we are doing
      pseudo division
   */
   
   fmpz_poly_fit_length(Q, q1->length+n2);
   fmpz_poly_fit_limbs(Q, FLINT_MAX(FLINT_MAX(q1->limbs + (s2*bits_B_lead)/FLINT_BITS+1, q2->limbs), 1));
   pow = (fmpz_t) flint_stack_alloc((bits_B_lead*s2)/FLINT_BITS+2);
   fmpz_pow_ui(pow, B_lead, s2);
   _fmpz_poly_scalar_mul_fmpz(Q, q1, pow);
   fmpz_poly_clear(q1);
   flint_stack_release();
   _fmpz_poly_left_shift(Q, Q, n2);
   _fmpz_poly_add(Q, Q, q2);
   fmpz_poly_clear(q2);
   
   /* 
      Set d to the power of lead(B) which everything 
      has been raised to
   */
   
   *d = s1+s2;
}

/****************************************************************************

   Powering

****************************************************************************/

void fmpz_poly_power(fmpz_poly_t output, const fmpz_poly_t poly_in, const unsigned long exp)
{
   if (exp == 0) 
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      fmpz_poly_set_coeff_ui(output, 0, 1);
      output->length = 1;
      return;
   }
   if (poly_in->length == 0)
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      output->length = 0;
      return;      
   }
   
   fmpz_poly_t poly;
	ulong trailing = 0L;
	while (!poly_in->coeffs[trailing*(poly_in->limbs + 1)]) trailing++;
	
   _fmpz_poly_attach_shift(poly, poly_in, trailing);
	
	if ((poly->length == 1) && ((poly->coeffs[0] == 1L) || (poly->coeffs[0] == -1L)) && (poly->coeffs[1] == 1L))
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1 + trailing*exp);
      output->length = 0;
		if ((exp & 1) && (poly->coeffs[0] == -1L))
			   fmpz_poly_set_coeff_si(output, trailing*exp, -1L);
		else fmpz_poly_set_coeff_ui(output, trailing*exp, 1L);
      output->length = 1 + trailing*exp;
      return;
   } 
   if (poly->length == 1)
   {
      fmpz_poly_fit_length(output, 1 + trailing*exp);
      fmpz_poly_fit_limbs(output, fmpz_size(poly->coeffs)*exp);
      _fmpz_poly_attach_shift(poly, poly_in, trailing);
		fmpz_pow_ui(output->coeffs + trailing*exp*(output->limbs + 1), poly->coeffs, exp);
      for (ulong j = 0; j < trailing*exp; j++)
			output->coeffs[j*(output->limbs + 1)] = 0L;
		output->length = 1 + trailing*exp;
      return;
   }
   
   //==================================================================================
   
   if (poly->length == 2) // Compute using binomial expansion
   {
      fmpz_t coeff1, coeff2;
      
      if (poly_in == output)
      {
         coeff1 = fmpz_init(poly->limbs);
         fmpz_set(coeff1, poly->coeffs);
         coeff2 = fmpz_init(poly->limbs);
         fmpz_set(coeff2, poly->coeffs + poly->limbs + 1);  
      } else 
      {
         coeff1 = poly->coeffs;
         coeff2 = poly->coeffs + poly->limbs + 1;
      }
             
      fmpz_poly_fit_length(output, exp + 1);     
      
      unsigned long bits2 = fmpz_bits(coeff2);
      
      // A rough estimate of the max number of limbs needed for a coefficient 
      unsigned long bits1 = fmpz_bits(coeff1);
      unsigned long bits = FLINT_MAX(bits1, bits2);
      
      fmpz_t pow;
      if (!(fmpz_is_one(coeff1) && fmpz_is_one(coeff2)))
      {
         bits = exp*(bits+1);
      } else
      {
         bits = exp;
      }
      
      fmpz_poly_fit_limbs(output, (bits-1)/FLINT_BITS+2);
      
      long i;
      unsigned long cbits;
      fmpz_t coeff_out = output->coeffs;
      fmpz_t last_coeff;
      
      if (fmpz_is_one(coeff2))
      {
         if (fmpz_is_one(coeff1))
         {
            coeff_out = output->coeffs;
            coeff_out[0] = 1;
            coeff_out[1] = 1;
               
            for (i = 1; i <= exp; i++)
            {
               last_coeff = coeff_out;
               coeff_out += (output->limbs+1);
               __fmpz_binomial_next(coeff_out, last_coeff, exp, i);
            }
         } else
         {
            fmpz_poly_set_coeff_ui(output, exp, 1);
            coeff_out = output->coeffs + exp*(output->limbs+1);   
            for (i = exp-1; i >= 0; i--)
            {
               coeff_out = output->coeffs + i*(output->limbs+1);   
               last_coeff = coeff_out + output->limbs+1;
               __fmpz_binomial_next(coeff_out, last_coeff, exp, exp - i);
               fmpz_mul(coeff_out, coeff_out, coeff1);
            }
         }
      } else
      {
         if (fmpz_is_one(coeff1))
         {
            coeff_out = output->coeffs;
            coeff_out[0] = 1;
            coeff_out[1] = 1;
               
            for (i = 1; i <= exp; i++)
            {
               output->length++;
               coeff_out = output->coeffs + i*(output->limbs+1);   
               last_coeff = coeff_out - output->limbs - 1;
               __fmpz_binomial_next(coeff_out, last_coeff, exp, i);
               fmpz_mul(coeff_out, coeff_out, coeff2);
            }
         } else
         {
            coeff_out = output->coeffs;
            fmpz_pow_ui(coeff_out, coeff1, exp);
               
            for (i = 1; i <= exp; i++)
            {
               output->length++;
               coeff_out = output->coeffs + i*(output->limbs+1);   
               last_coeff = coeff_out - output->limbs - 1;
               fmpz_tdiv(coeff_out, last_coeff, coeff1);
               __fmpz_binomial_next(coeff_out, coeff_out, exp, i);
               fmpz_mul(coeff_out, coeff_out, coeff2);
            }
         }
      }
      
      output->length = exp + 1;
      
      if (poly_in == output)
      {
         fmpz_clear(coeff1);
         fmpz_clear(coeff2);
      }      
        
		fmpz_poly_left_shift(output, output, trailing*exp);
      return;
   }
   //===================================================================================
   fmpz_poly_t temp;
   fmpz_poly_init(temp);
   
   unsigned long bits = FLINT_BIT_COUNT(exp);
   
   fmpz_poly_t polycopy;
   
   if (poly_in == output)
   {
      fmpz_poly_init(polycopy);
      fmpz_poly_set(polycopy, poly);
   } else _fmpz_poly_attach(polycopy, poly);
   
	fmpz_poly_fit_length(output, poly->length);
   fmpz_poly_fit_limbs(output, poly->limbs);
   
	_fmpz_poly_set(output, polycopy);
   
   while (bits > 1)
   {
      fmpz_poly_mul(output, output, output);
      if ((1L<<(bits-2)) & exp)
      {
         fmpz_poly_mul(output, output, polycopy);
      }
      bits--;
   } 
   
   if (poly_in == output) fmpz_poly_clear(polycopy);

	fmpz_poly_left_shift(output, output, trailing*exp);
}

void fmpz_poly_power_trunc_n(fmpz_poly_t output, const fmpz_poly_t poly_i, const unsigned long exponent, const unsigned long n)
{
   unsigned long exp = exponent;
   fmpz_poly_t power, poly_in;

	if (exp == 0) 
   {
		fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      if ((n) && (poly_i->length))
		{
			fmpz_poly_set_coeff_ui(output, 0, 1);
         output->length = 1;
		} else output->length = 0;
		return;
   }

	_fmpz_poly_attach_truncate(poly_in, poly_i, n);
   
   if ((poly_in->length == 0) || (n == 0))
   {
      fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      output->length = 0;
      return;      
   }
   
	fmpz_poly_t poly;
	ulong trailing = 0L;
	while (!poly_in->coeffs[trailing*(poly_in->limbs + 1)]) trailing++;
	
   if (trailing*exp >= n) // the first n coeffs are zero
	{
	   fmpz_poly_fit_limbs(output, 1);
      fmpz_poly_fit_length(output, 1);
      output->length = 0;
      return;  
   }

	ulong nn = n - exp*trailing;
	_fmpz_poly_attach_shift(poly, poly_in, trailing);
	_fmpz_poly_truncate(poly, nn);
   
	if (poly->length == 1)
	{    
	   if (((poly->coeffs[0] == 1L) || (poly->coeffs[0] == -1L)) && (poly->coeffs[1] == 1L))
      {		
         fmpz_poly_fit_limbs(output, 1);
         fmpz_poly_fit_length(output, 1 + trailing*exp);
         output->length = 0;
		   if ((exp & 1) && (poly->coeffs[0] == -1L))
			   fmpz_poly_set_coeff_si(output, trailing*exp, -1L);
			else fmpz_poly_set_coeff_ui(output, trailing*exp, 1L);
         output->length = 1 + trailing*exp;
         return;
      } 

	   fmpz_poly_fit_length(output, 1 + trailing*exp);
      fmpz_poly_fit_limbs(output, fmpz_size(poly->coeffs)*exp);
      _fmpz_poly_attach_truncate(poly_in, poly_i, n);
      _fmpz_poly_attach_shift(poly, poly_in, trailing);
		fmpz_pow_ui(output->coeffs + trailing*exp*(output->limbs + 1), poly->coeffs, exp);
      for (ulong j = 0; j < trailing*exp; j++)
			output->coeffs[j*(output->limbs + 1)] = 0L;
		output->length = 1 + trailing*exp;
      return;
   }
   
   fmpz_poly_fit_length(output, n);  // Set output to poly
   fmpz_poly_fit_limbs(output, poly->limbs);
   _fmpz_poly_attach_truncate(poly_in, poly_i, n);
   _fmpz_poly_attach_shift(poly, poly_in, trailing);
	_fmpz_poly_truncate(poly, nn);

	fmpz_poly_t polycopy;
	if (poly_in == output)
	{
		fmpz_poly_init(polycopy);
		fmpz_poly_set(polycopy, poly);
	} else _fmpz_poly_attach(polycopy, poly);
	
   fmpz_poly_set(output, polycopy);

	ulong expcpy = exp;
   
   while (!(exp & 1L))  // Square until we get to the first binary 1 in the exponent
   {
      fmpz_poly_mul_trunc_n(output, output, output, nn);  
      exp >>= 1;
   }
   
   exp >>= 1;
   if (exp) // Exponent is not just a power of 2, so keep multiplying by higher powers
   {
      fmpz_poly_init(power);
      fmpz_poly_fit_length(power, nn);
      fmpz_poly_fit_limbs(power, output->limbs);
      _fmpz_poly_set(power, output);
      
      while (exp)
      {
         fmpz_poly_mul_trunc_n(power, power, power, nn);
         if (exp & 1) 
         {
            fmpz_poly_mul_trunc_n(output, output, power, nn);
         }
         exp >>= 1;
      }
	  fmpz_poly_clear(power);
   }
   
	if (poly_in == output) 
		fmpz_poly_clear(polycopy);
   
	fmpz_poly_left_shift(output, output, expcpy*trailing);
}

/****************************************************************************

   Content

****************************************************************************/

void fmpz_poly_content(fmpz_t c, const fmpz_poly_t poly)
{
   unsigned long length = poly->length;
   
   if (length == 0) 
   {
      fmpz_set_ui(c, 0L);
      return;
   }
   
   if (length == 1)
   {
      fmpz_set(c, poly->coeffs);
      if ((long) c[0] < 0L) c[0] = -c[0];
      return;
   }
   
   fmpz_t coeff = fmpz_poly_get_coeff_ptr(poly, length - 1);
   fmpz_set(c, coeff);
   
   for (long i = length - 2; (i >= 0L) && !fmpz_is_one(c); i--)
   {
      coeff = fmpz_poly_get_coeff_ptr(poly, i);
      fmpz_gcd(c, c, coeff);
   }
}

/****************************************************************************

   GCD

****************************************************************************/

void fmpz_poly_gcd_subresultant(fmpz_poly_t D, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   fmpz_poly_t Ain, Bin;
     
   if (poly2->length > poly1->length)
   {
      _fmpz_poly_attach(Ain, poly2);
      _fmpz_poly_attach(Bin, poly1);
   } else
   {
      _fmpz_poly_attach(Ain, poly1);
      _fmpz_poly_attach(Bin, poly2);      
   }
   if (Bin->length == 0)
   {
      fmpz_poly_set(D, Ain);
      return;
   }
   
   fmpz_t a, b, d;
   a = fmpz_init(Ain->limbs);
   b = fmpz_init(Bin->limbs);
   fmpz_poly_content(a, Ain);
   fmpz_poly_content(b, Bin);
   
	d = fmpz_init(FLINT_MIN(fmpz_size(a), fmpz_size(b)));
   fmpz_gcd(d, a, b);
   
   fmpz_poly_t A, B, Q, R;
   fmpz_poly_init(A);
   fmpz_poly_init(B);
   fmpz_poly_init(Q);
   fmpz_poly_init(R);
   unsigned long s;
   
   fmpz_poly_scalar_div_fmpz(A, Ain, a);
   fmpz_poly_scalar_div_fmpz(B, Bin, b);
   fmpz_clear(b); //release b
   fmpz_clear(a); //release a
   
   int done = 0;
   
   fmpz_t g;
   fmpz_t h = fmpz_init(1);
   fmpz_t one = fmpz_init(1);
   fmpz_set_ui(h, 1UL);
   fmpz_set_ui(one, 1UL);
   g = one;
   unsigned long olddelta = 1;
        
   while (!done)
   {
      
		unsigned long delta = A->length - B->length;
      fmpz_poly_pseudo_divrem(Q, R, &s, A, B);
      
		if (R->length > 1)
      {
         fmpz_poly_swap(A, B);
         fmpz_t r;
         if (olddelta == 1) 
         {
            r = fmpz_init((delta+1)*fmpz_size(g)+1);
            fmpz_pow_ui(r, g, delta+1);
         } else
         {
            r = fmpz_init(fmpz_size(g) + delta*fmpz_size(h)+1);
            fmpz_pow_ui(r, h, delta);
            fmpz_mul(r, r, g);
         }
         
         g = fmpz_poly_get_coeff_ptr(A, A->length - 1);
         fmpz_t temp = fmpz_init((delta-s+1)*fmpz_size(g)+1);
         fmpz_pow_ui(temp, g, delta-s+1);
         fmpz_poly_scalar_mul_fmpz(R, R, temp);
         fmpz_clear(temp); // release temp
         
         fmpz_poly_scalar_div_fmpz(B, R, r);
         fmpz_clear(r); // release r
         
         olddelta = delta;
         if (delta == 0)
         {
            fmpz_clear(h);
            h = fmpz_init(delta*fmpz_size(g)+1);
            fmpz_pow_ui(h, g, delta);
         } else if (delta == 1)
         {
            olddelta = 1;
            fmpz_clear(h);
            h = fmpz_init(fmpz_size(g));
            fmpz_set(h, g);
         } else
         {
            temp = fmpz_init((delta-1)*fmpz_size(h)+1);
            fmpz_pow_ui(temp, h, delta - 1);
            fmpz_clear(h);
            h = fmpz_init(delta*fmpz_size(g)+1);
            fmpz_t temp2 = fmpz_init(delta*fmpz_size(g)+1);
            fmpz_pow_ui(temp2, g, delta);
            fmpz_fdiv(h, temp2, temp);
            fmpz_clear(temp2); // release temp2
            fmpz_clear(temp); // release temp
         }
      } else
      {
         if (R->length == 1)
         {
            fmpz_poly_zero(B);
            fmpz_poly_set_coeff_ui(B, 0, 1UL);
         }
         done = 1;
      }
   }
   
   b = fmpz_init(B->limbs+1);
   fmpz_poly_content(b, B);
   fmpz_poly_scalar_div_fmpz(D, B, b);
   fmpz_poly_scalar_mul_fmpz(D, D, d);
   fmpz_clear(b); // release b
         
   if ((long) (_fmpz_poly_lead(D)[0]) < 0L) fmpz_poly_neg(D, D);
   
   fmpz_clear(h);
   fmpz_clear(one);
   fmpz_poly_clear(A);
   fmpz_poly_clear(B);
   fmpz_poly_clear(Q);
   fmpz_poly_clear(R);
   fmpz_clear(d); //release d
}

void fmpz_poly_gcd_modular(fmpz_poly_t H, const fmpz_poly_t poly1, 
					    const fmpz_poly_t poly2, const ulong bits1_in1, const ulong bits2_in1)
{
   fmpz_poly_t Ain, Bin;
	ulong bits1_in = bits1_in1;
	ulong bits2_in = bits2_in1;
      
   if (poly2->length > poly1->length)
   {
      _fmpz_poly_attach(Ain, poly2);
      _fmpz_poly_attach(Bin, poly1);
      bits1_in = bits2_in1;
	   bits2_in = bits1_in1;
   } else
   {
      _fmpz_poly_attach(Ain, poly1);
      _fmpz_poly_attach(Bin, poly2);      
   }
   if (Bin->length == 0)
   {
      fmpz_poly_set(H, Ain);
      return;
   }
       
   fmpz_t ac, bc, d;
   ac = fmpz_init(Ain->limbs);
   bc = fmpz_init(Bin->limbs);
   fmpz_poly_content(ac, Ain);
   fmpz_poly_content(bc, Bin);
   
	d = fmpz_init(FLINT_MIN(fmpz_size(ac), fmpz_size(bc)));
   fmpz_gcd(d, ac, bc);

   if (Bin->length == 1)
   {
      fmpz_poly_set_coeff_fmpz(H, 0, d);
      H->length = 1;
      fmpz_clear(d);
      fmpz_clear(ac); // release ac
      fmpz_clear(bc); // release bc
      return;
   }
   
   fmpz_poly_t A, B;
   fmpz_poly_init(A);
   fmpz_poly_init(B);
   
   fmpz_poly_scalar_div_fmpz(A, Ain, ac);
   fmpz_poly_scalar_div_fmpz(B, Bin, bc);
   fmpz_clear(bc); //release bc
   fmpz_clear(ac); //release ac
	
	ulong nb1, nb2, bits1, bits2, bound;

	if (bits1_in) bits1 = bits1_in;
	else bits1 = FLINT_ABS(fmpz_poly_max_bits(A));
	if (bits2_in) bits2 = bits2_in;
	else bits2 = FLINT_ABS(fmpz_poly_max_bits(B));

	fmpz_t lead_A = _fmpz_poly_lead(A);
   fmpz_t lead_B = _fmpz_poly_lead(B);

	if ((A->length < 64) && (B->length < 64))
	{
		nb1 = fmpz_poly_2norm_bits_normalised(A);
		nb2 = fmpz_poly_2norm_bits_normalised(B);
	} else // approximate to save time 
	{
		nb1 = (2*bits1 + FLINT_BIT_COUNT(A->length) + 1)/2 - fmpz_bits(lead_A) + 1;
		nb2 = (2*bits2 + FLINT_BIT_COUNT(B->length) + 1)/2 - fmpz_bits(lead_B) + 1;
	}
   
   fmpz_t g = fmpz_init(FLINT_MIN(FLINT_ABS(lead_A[0]), FLINT_ABS(lead_B[0])));
   fmpz_gcd(g, lead_A, lead_B);
   ulong gbits = fmpz_bits(g);
   int g_pm1 = 0;
   if ((FLINT_ABS(g[0]) == 1L) && (g[1] == 1L)) g_pm1 = 1;

   fmpz_t eval_A = fmpz_init(A->limbs + 1);
   fmpz_t eval_B = fmpz_init(B->limbs + 1);
	fmpz_t coeff_A = A->coeffs;
   fmpz_t coeff_B = B->coeffs;
   ulong size_A = A->limbs + 1;
	ulong size_B = B->limbs + 1;
   eval_A[0] = 0;
	eval_B[0] = 0;

	for (ulong i = 0; i < A->length; i++)
	{
      if (i&1) fmpz_add(eval_A, eval_A, coeff_A);
		else fmpz_sub(eval_A, eval_A, coeff_A);
		coeff_A += size_A;
	}

   for (ulong i = 0; i < B->length; i++)
	{
      if (i&1) fmpz_add(eval_B, eval_B, coeff_B);
		else fmpz_sub(eval_B, eval_B, coeff_B);
		coeff_B += size_B;
	}

   fmpz_t eval_GCD = fmpz_init(FLINT_MAX(FLINT_ABS(eval_A[0]), FLINT_ABS(eval_B[0])));
	fmpz_gcd(eval_GCD, eval_A, eval_B);
   long bits_small = FLINT_MAX(fmpz_bits(eval_GCD), fmpz_bits(g));
	if (bits_small < 2L) bits_small = 2;
	fmpz_clear(eval_GCD);
	fmpz_clear(eval_A);
	fmpz_clear(eval_B);

	ulong p;
   ulong pbits;
   
#if FLINT_BIT == 64
   if (bits_small < FLINT_D_BITS) 
   {
      pbits = FLINT_D_BITS;
      p = (1L<<FLINT_D_BITS) - 200;
   } else 
   {
#endif
		pbits = FLINT_BITS - 1;
      p = (1L<<(FLINT_BITS-2));
#if FLINT_BIT == 64
   }
#endif

   zmod_poly_t a, b, h;
   
   int first = 1;

   unsigned long n = B->length;

   unsigned long modsize = FLINT_MAX(A->limbs, B->limbs);
   fmpz_t modulus = fmpz_init(modsize);
   modulus[0] = 0L;

   fmpz_poly_t Q;
   fmpz_poly_init(Q);
      
   for (;;)
   {
      if (!first)
      {
         zmod_poly_clear(a);
         zmod_poly_clear(b);
         zmod_poly_clear(h); 
      } else first = 0;
      do { p = z_nextprime(p); }
      while (!fmpz_mod_ui(g, p)); 
      zmod_poly_init(a, p);
      zmod_poly_init(b, p);
      zmod_poly_init(h, p);
		
      if (bits1 + 1 < pbits) fmpz_poly_to_zmod_poly_no_red(a, A);
      else fmpz_poly_to_zmod_poly(a, A);
      if (bits2 + 1 < pbits) fmpz_poly_to_zmod_poly_no_red(b, B);
      else fmpz_poly_to_zmod_poly(b, B);
      
      zmod_poly_gcd(h, a, b);
      
      if (h->length == 1) // gcd is 1
      {
         fmpz_poly_set_coeff_ui(H, 0, 1L);
         H->length = 1;
			break; 
      }
      
      if (h->length - 1 > n) // discard
         continue;      
      
      if (g_pm1) zmod_poly_make_monic(h, h);
      else
      {
         unsigned long h_inv = z_invert(h->coeffs[h->length-1], h->p);
         unsigned long g_mod = fmpz_mod_ui(g, h->p);
         h_inv = z_mulmod2_precomp(h_inv, g_mod, h->p, h->p_inv);
         zmod_poly_scalar_mul(h, h, h_inv);
      }
      
      if (h->length - 1 < n)
      {
         zmod_poly_to_fmpz_poly(H, h);
         
         if (g_pm1)
         {
            if (fmpz_poly_divides_modular(Q, A, H, 0) && fmpz_poly_divides_modular(Q, B, H, 0)) 
				{
				   break;
				}
         } else
         {
            /*
				   Bound is easily derived from Thm 5.3 and Cor 5.4 of 
				   http://compalg.inf.elte.hu/~tony/Informatikai-Konyvtar/03-Algorithms%20of%20Informatics%201,%202,%203/CompAlg29May.pdf
					The + 1 is to allow for signed coefficients after Chinese Remaindering
				*/
				bound = h->length + FLINT_MIN(nb1, nb2) + gbits + 1;
            if (pbits > bound)
            { 
               fmpz_t hc = fmpz_init(H->limbs);
               fmpz_poly_content(hc, H);
               fmpz_poly_scalar_div_fmpz(H, H, hc);
               fmpz_clear(hc); // release hc
               break;
            }

				if (pbits >= bits_small)
				{
               fmpz_t hc = fmpz_init(H->limbs);
               fmpz_poly_content(hc, H);
               fmpz_poly_scalar_div_fmpz(H, H, hc);
               if (fmpz_poly_divides_modular(Q, A, H, 0) && fmpz_poly_divides_modular(Q, B, H, 0)) 
					{
						fmpz_clear(hc); // release hc
						break;
					}
					fmpz_poly_scalar_mul_fmpz(H, H, hc);
					fmpz_clear(hc); // release hc
				}
         }

         fmpz_set_ui(modulus, p);
         n = h->length - 1;
         continue;
      }
      
      fmpz_t newmod = fmpz_init(modulus[0] + 1);
      
      if (g_pm1)
      {
         if ((fmpz_poly_CRT(H, H, h, newmod, modulus)) || (fmpz_bits(newmod) > bits_small))
			{
				if (fmpz_poly_divides_modular(Q, A, H, 0) && fmpz_poly_divides_modular(Q, B, H, 0)) 
            {
               fmpz_clear(newmod); // release newmod
               break;
            }
			}
      } else
      {
         if (fmpz_poly_CRT(H, H, h, newmod, modulus) || (fmpz_bits(newmod) > bound) || (fmpz_bits(newmod) >= bits_small))
         { 
            fmpz_t hc = fmpz_init(H->limbs);
            fmpz_poly_content(hc, H);
            fmpz_poly_scalar_div_fmpz(H, H, hc);
            if ((fmpz_bits(newmod) > bound) || (fmpz_poly_divides_modular(Q, A, H, 0) && fmpz_poly_divides_modular(Q, B, H, 0))) 
            {
               fmpz_clear(newmod); // release newmod
               fmpz_clear(hc); // release hc
               break;
            }
				fmpz_poly_scalar_mul_fmpz(H, H, hc);
				fmpz_clear(hc);
         }
      }
      
      if (newmod[0] >= modsize) 
      {
         modulus = fmpz_realloc(modulus, (modsize+8));
         modsize += 8;
      }
      
      fmpz_set(modulus, newmod);
      fmpz_clear(newmod); // release newmod
   }
   fmpz_clear(g); // release g
   zmod_poly_clear(a);
   zmod_poly_clear(b);
   zmod_poly_clear(h); 
        
   fmpz_poly_scalar_mul_fmpz(H, H, d);

   fmpz_poly_clear(A);
   fmpz_poly_clear(B);
   fmpz_poly_clear(Q);
   fmpz_clear(modulus);
   fmpz_clear(d); //release d
}

void fmpz_poly_limb_unpack_wrap(fmpz_poly_t H, fmpz_t arrayg, ulong pack_bits)
{
#if FLINT_BITS == 64
	if (pack_bits <= 32)
	{
		ulong length = (arrayg[0]*FLINT_BITS)/pack_bits + 1; // may have one extra coeff
		                                        // due to 1 0 -x being packed as 0 -1 -x
		//ulong rem = length*pack_bits % FLINT_BITS;
		//if ((rem) && (arrayg[arrayg[0]]>>rem)) length ++;
		fmpz_poly_fit_length(H, length);
		fmpz_poly_fit_limbs(H, 1);
		for (ulong i = 0; i < length; i++) H->coeffs[i*(H->limbs+1)] = 0;
		fmpz_poly_bit_unpack(H, arrayg + 1, length, pack_bits);
	} else
#endif
 	if (pack_bits == FLINT_BITS) 
	{
	   fmpz_poly_limb_unpack_1(H, arrayg + 1, arrayg[0]);
	}
   else
	{
	   ulong pack_limbs = (pack_bits >> FLINT_LG_BITS_PER_LIMB);
		fmpz_poly_limb_unpack(H, arrayg + 1, (arrayg[0] - 1)/pack_limbs + 1, pack_limbs);
	}
}

int fmpz_poly_gcd_heuristic(fmpz_poly_t H, const fmpz_poly_t poly1, 
					const fmpz_poly_t poly2, const ulong bits1_in, const ulong bits2_in)
{
	if (poly2->length == 0)
   {
      fmpz_poly_set(H, poly1);
      return 1;
   }

	if (poly1->length == 0)
	{
		fmpz_poly_set(H, poly2);
      return 1;
	}

	fmpz_t ac, bc, d;
   ac = fmpz_stack_init(poly1->limbs);
   bc = fmpz_stack_init(poly2->limbs);
   
	fmpz_poly_content(ac, poly1);
   fmpz_poly_content(bc, poly2);

	d = fmpz_init(FLINT_MIN(fmpz_size(ac), fmpz_size(bc)));
   fmpz_gcd(d, ac, bc);

   if ((poly1->length == 1) || (poly2->length == 1))
   {
      fmpz_poly_set_coeff_fmpz(H, 0, d);
      H->length = 1;
      fmpz_stack_release(); //release bc
      fmpz_stack_release(); //release ac
	   fmpz_clear(d);
      return 1;
   }
   
   fmpz_poly_t A, B;
   fmpz_poly_init(A);
   fmpz_poly_init(B);
   
   fmpz_poly_scalar_div_fmpz(A, poly1, ac);
   fmpz_poly_scalar_div_fmpz(B, poly2, bc);
   fmpz_stack_release(); //release bc
   fmpz_stack_release(); //release ac

	if (A->length == 2) 
	{
		fmpz_poly_t Q;
		fmpz_poly_init(Q);
		if (fmpz_poly_divides_modular(Q, B, A, 0))
		{
			fmpz_poly_scalar_mul_fmpz(H, A, d);
		}
		else  
		{
			fmpz_poly_set_coeff_fmpz(H, 0, d);
         H->length = 1;
		}
		fmpz_clear(d);
		fmpz_poly_clear(A);
      fmpz_poly_clear(B);
      fmpz_poly_clear(Q);
      return 1;
	}
	
	if (B->length == 2) 
	{
		fmpz_poly_t Q;
		fmpz_poly_init(Q);
		if (fmpz_poly_divides_modular(Q, A, B, 0))
		   fmpz_poly_scalar_mul_fmpz(H, B, d);
		else  
		{
			fmpz_poly_set_coeff_fmpz(H, 0, d);
         H->length = 1;
		}
		fmpz_clear(d);
		fmpz_poly_clear(A);
      fmpz_poly_clear(B);
      fmpz_poly_clear(Q);
      return 1;
	}
	
	ulong bits1, bits2;
	if (bits1_in) bits1 = bits1_in;
	else bits1 = FLINT_ABS(fmpz_poly_max_bits(A));
	if (bits2_in) bits2 = bits2_in;
	else bits2 = FLINT_ABS(fmpz_poly_max_bits(B));
	
	ulong max_bits = FLINT_MAX(bits1, bits2);
	ulong pack_limbs = (max_bits - 1)/FLINT_BITS + 1;
   
	long sign1 = ((long) fmpz_poly_lead(A)[0] < 0L) ? -1L : 1L;
	long sign2 = ((long) fmpz_poly_lead(B)[0] < 0L) ? -1L : 1L;
   ulong pack_bits;

	/*
	   This bound ensures that if H | A and H | B with H primitive then H is the 
		gcd of A and B. The bound is taken from 
		http://arxiv.org/abs/cs/0206032v1
	*/
			
	ulong bound_bits = FLINT_MIN(bits1, bits2) + 6; // can use 
	
	long bound_limbs = (bound_bits - 1)/FLINT_BITS + 1; 
	if (bound_limbs > pack_limbs) pack_limbs = bound_limbs;

	fmpz_t array1 = fmpz_stack_init(poly1->length*pack_limbs + 1);
   fmpz_t array2 = fmpz_stack_init(poly2->length*pack_limbs + 1);
   fmpz_t arrayg = fmpz_stack_init(pack_limbs*FLINT_MAX(poly1->length, poly2->length) + 1);

	F_mpn_clear(array1, poly1->length*pack_limbs + 2);
   F_mpn_clear(array2, poly2->length*pack_limbs + 2);
   F_mpn_clear(arrayg, pack_limbs*FLINT_MAX(poly1->length, poly2->length) + 2);
	
#if FLINT_BITS == 64
	if ((bits1 < 32) && (bits2 < 32) && (bound_bits < 32))
	{
	   pack_bits = FLINT_MAX(bits1, bits2) + 6;
		if (pack_bits > 32) pack_bits = 32;
		if (pack_bits < bound_bits) pack_bits = bound_bits;
		
		fmpz_poly_bit_pack(array1 + 1, A, A->length, -pack_bits, sign1);
	   array1[0] = (pack_bits*A->length - 1)/FLINT_BITS + 1;

      fmpz_poly_bit_pack(array2 + 1, B, B->length, -pack_bits, sign2);
	   array2[0] = (pack_bits*B->length - 1)/FLINT_BITS + 1;
	} else 
#endif
	if ((bits1 < FLINT_BITS) && (bits2 < FLINT_BITS) && (bound_bits < FLINT_BITS))
	{
		pack_bits = FLINT_BITS;
		if (sign1 < 0L) fmpz_poly_limb_pack_neg_1(array1 + 1, A);
	   else fmpz_poly_limb_pack_1(array1 + 1, A);
	   array1[0] = poly1->length;

      if (sign2 < 0L) fmpz_poly_limb_pack_neg_1(array2 + 1, B);
	   else fmpz_poly_limb_pack_1(array2 + 1, B);
      array2[0] = poly2->length;
   } else
	{
		pack_bits = FLINT_BITS*pack_limbs;
		if (sign1 < 0L) fmpz_poly_limb_pack_neg(array1 + 1, A, A->length, pack_limbs);
	   else fmpz_poly_limb_pack(array1 + 1, A, A->length, pack_limbs);
	   array1[0] = poly1->length*pack_limbs;
      
      if (sign2 < 0L) fmpz_poly_limb_pack_neg(array2 + 1, B, B->length, pack_limbs);
	   else fmpz_poly_limb_pack(array2 + 1, B, B->length, pack_limbs);
	   array2[0] = poly2->length*pack_limbs;
	}

   NORM(array1);
	NORM(array2);
		
	fmpz_gcd(arrayg, array1, array2);
	
   fmpz_poly_limb_unpack_wrap(H, arrayg, pack_bits);
   
	fmpz_poly_t Q, Q1;
	fmpz_poly_init(Q);
   fmpz_poly_init(Q1);

	fmpz_t hc = fmpz_stack_init(H->limbs);
	fmpz_t temp = fmpz_stack_init(FLINT_ABS(arrayg[0]) + 1);
	fmpz_poly_content(hc, H);

	if (!fmpz_is_one(hc)) 
	  if (hc[0] == 1) fmpz_tdiv_ui(temp, arrayg, hc[1]);
	  else fmpz_tdiv(temp, arrayg, hc);
	else fmpz_set(temp, arrayg);

	ulong qlimbs = FLINT_MAX(FLINT_ABS(array1[0]), FLINT_ABS(array2[0])) + pack_limbs + 1;
	fmpz_t q = fmpz_stack_init(qlimbs);
   
	int divides = 0;

	if (fmpz_divides(q, array1, temp))
	{
      F_mpn_clear(q + FLINT_ABS(q[0]) + 1, FLINT_MIN(qlimbs - FLINT_ABS(q[0]) - 1, pack_limbs + 1)); // clear additional limbs for limb_unpack
		fmpz_poly_limb_unpack_wrap(Q1, q, pack_bits);
      ulong bits_H = FLINT_ABS(fmpz_poly_max_bits(H));
		ulong bits_Q1 = FLINT_ABS(fmpz_poly_max_bits(Q1));
		ulong log1 = ceil_log2(H->length);
		ulong log_length = FLINT_MIN(log1, ceil_log2(Q1->length));
		int ok = 0;
		if (bits_H + bits_Q1 + log_length < pack_bits)
		{
			ok = 1;
		} else
		{
			fmpz_poly_scalar_div_fmpz(H, H, hc);
		   fmpz_poly_mul(Q1, Q1, H);
		   if (sign1 < 0L) fmpz_poly_neg(Q1, Q1);
		}
		if (ok || fmpz_poly_equal(Q1, A))
		{
         if (fmpz_divides(q, array2, temp))
	      {
            F_mpn_clear(q + FLINT_ABS(q[0]) + 1, FLINT_MIN(qlimbs - FLINT_ABS(q[0]) - 1, pack_limbs)); // clear additional limbs for limb_unpack
		      if (ok) fmpz_poly_scalar_div_fmpz(H, H, hc);
		      fmpz_poly_limb_unpack_wrap(Q1, q, pack_bits);
            bits_Q1 = FLINT_ABS(fmpz_poly_max_bits(Q1));
				log_length = FLINT_MIN(log1, ceil_log2(Q1->length));
				ok = 0;
				if (bits_H + bits_Q1 + log_length < pack_bits)
		      {
			      ok = 1;
		      } else
				{
					fmpz_poly_mul(Q1, Q1, H);
		         if (sign2 < 0L) fmpz_poly_neg(Q1, Q1);
				}
				if (ok || fmpz_poly_equal(Q1, B)) divides = 1;
			} 
		} 
	}

	fmpz_stack_release(); // release q
	fmpz_stack_release(); // release temp
	fmpz_stack_release(); // release hc
	fmpz_stack_release(); // release arrayg
	fmpz_stack_release(); // release array2
	fmpz_stack_release(); // release array1
	fmpz_poly_clear(Q1);
   fmpz_poly_clear(Q);
	fmpz_poly_clear(A);
	fmpz_poly_clear(B);
			
	if (divides)
	{
	   fmpz_poly_scalar_mul_fmpz(H, H, d);
		fmpz_clear(d);
		return 1;
	} else
	{
	   fmpz_clear(d);
	   return 0;
	}
}

void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
   if (poly1 == poly2)
   {
      if (res != poly1)
         fmpz_poly_set(res, poly1);
      return;
   }

	ulong max_length = FLINT_MAX(poly1->length, poly2->length);
	ulong limbs1 = fmpz_poly_max_limbs(poly1);
	ulong limbs2 = fmpz_poly_max_limbs(poly2);
   ulong max_limbs = FLINT_MAX(limbs1, limbs2);
   
	if ((max_length <= 80) && (max_limbs <= 16))
	{
		if (fmpz_poly_gcd_heuristic(res, poly1, poly2, 0, 0))
		   return;

		if (max_length <= 1) 
	   {
		   fmpz_poly_gcd_modular(res, poly1, poly2, 0, 0);
            return;
      }
		 
		if ((max_length <= 4) || ((max_limbs <= 1) && (max_length <= 12)))
      {
         fmpz_poly_gcd_subresultant(res, poly1, poly2);
         return;
      }

      fmpz_poly_gcd_modular(res, poly1, poly2, 0, 0);
         return;  
	}

   if (max_length <= 6) 
   {
      fmpz_poly_gcd_subresultant(res, poly1, poly2);
      return;
   }
   
   if (max_limbs > 16)
   {
      fmpz_poly_gcd_modular(res, poly1, poly2, 0, 0);
      return;
   }

	ulong bits1 = FLINT_ABS(fmpz_poly_max_bits(poly1));
	ulong bits2 = FLINT_ABS(fmpz_poly_max_bits(poly2));
   ulong max_bits = FLINT_MAX(bits1, bits2);

	if (max_limbs <= 2)
	{
		if (max_bits <= 32)
		{
			if (max_length*max_bits < 1650000)
			{
				if (fmpz_poly_gcd_heuristic(res, poly1, poly2, bits1, bits2))
		         return;
			}
		} else if (max_limbs <= 1)
		{
			if (max_length < 9500)
		   {
			   if (fmpz_poly_gcd_heuristic(res, poly1, poly2, bits1, bits2))
		         return;
		   }
	   } else if (max_length < 3000)
		{
			if (fmpz_poly_gcd_heuristic(res, poly1, poly2, bits1, bits2))
		      return;
		}
	} else if (max_length * max_limbs < 4000)
	{
		if (fmpz_poly_gcd_heuristic(res, poly1, poly2, bits1, bits2))
		   return;
	}

   fmpz_poly_gcd_modular(res, poly1, poly2, bits1, bits2);
}

/*
   Invert poly1 modulo poly2 with denominator (not guaranteed to be the resultant)
   i.e. H*poly1 = d modulo poly2
   Assumes poly1 is reduced modulo poly2, which is monic and irreducible
   Assumes d has enough space to store fmpz_poly_resultant_bound(poly1, poly2)/FLINT_BITS + 2 limbs
*/

void fmpz_poly_invmod_modular(fmpz_t d, fmpz_poly_t H, fmpz_poly_t poly1, fmpz_poly_t poly2)
{
   FLINT_ASSERT(poly2->length > poly1->length);
   
   if ((poly1->length == 1) && (poly2->length > 1))
   {
      fmpz_set(d, poly1->coeffs);
      fmpz_poly_set_coeff_ui(H, 0, 1L);
      H->length = 1;
      return;
   }

   if ((poly1->length == 0) || (poly2->length == 0)) 
   {
      printf("Error: divide by zero!\n");
      abort();
   }
 
   fmpz_poly_t A, B;
   fmpz_poly_init(A);
   fmpz_poly_init(B);
   
   fmpz_poly_set(A, poly1);
   fmpz_poly_set(B, poly2);

   fmpz_poly_t prod, quot, rem;
   fmpz_poly_init(prod);
   fmpz_poly_init(quot);
   fmpz_poly_init(rem);
   
   unsigned long p = (1L<<(FLINT_BITS-2));
   
   zmod_poly_t a, b, h;
   
   int first = 1;

   unsigned long n = B->length;

   unsigned long modsize = FLINT_MAX(A->limbs, B->limbs);
   fmpz_t modulus = fmpz_init(modsize);
   modulus[0] = 0L;

   fmpz_poly_t Q;
   fmpz_poly_init(Q);
   for (;;)
   {
      if (!first)
      {
         zmod_poly_clear(a);
         zmod_poly_clear(b);
         zmod_poly_clear(h); 
      } 
      p = z_nextprime(p);  
      zmod_poly_init(a, p);
      zmod_poly_init(b, p);
      zmod_poly_init(h, p);

      fmpz_poly_to_zmod_poly(a, A);
      fmpz_poly_to_zmod_poly(b, B);
      
      unsigned long r = zmod_poly_resultant(a, b);

      if ((fmpz_mod_ui(_fmpz_poly_lead(A), p) == 0L) || (r == 0L))
      {
         continue;
      }

      unsigned long coprime = zmod_poly_gcd_invert(h, a, b);
      if (!coprime) 
      {
         continue;
      }

      zmod_poly_scalar_mul(h, h, r);

      if (first)
      {
         zmod_poly_to_fmpz_poly(H, h);
         fmpz_set_ui(modulus, p);
         first = 0;
         continue;
      } 
      fmpz_t newmod = fmpz_init(modulus[0] + 1);
      
      if (fmpz_poly_CRT(H, H, h, newmod, modulus))
      {
         fmpz_t hc = fmpz_init(H->limbs);
         fmpz_poly_content(hc, H);
         fmpz_poly_scalar_div_fmpz(H, H, hc);
         fmpz_clear(hc); // release hc
         fmpz_clear(newmod); // release newmod
         fmpz_poly_mul(prod, H, poly1);
         fmpz_poly_divrem(quot, rem, prod, poly2);
         if (rem->length == 1)
         {
            fmpz_set(d, rem->coeffs);
            break;
         }
      }
      
      if (newmod[0] >= modsize) 
      {
         modulus = fmpz_realloc(modulus, (modsize+8));
         modsize += 8;
      }
      
      fmpz_set(modulus, newmod);
      fmpz_clear(newmod); // release newmod
   }
   zmod_poly_clear(a);
   zmod_poly_clear(b);
   zmod_poly_clear(h); 
       
   fmpz_poly_clear(quot);
   fmpz_poly_clear(rem);
   fmpz_poly_clear(prod);
   
   fmpz_poly_clear(A);
   fmpz_poly_clear(B);
   fmpz_poly_clear(Q);
   fmpz_clear(modulus);
}

void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, fmpz_poly_t a, fmpz_poly_t b)
{
   fmpz_poly_resultant(r, a, b);

   if (r[0] == 0) 
   {
      return;
   }

   int stabilised = 0;

   fmpz_t prod = fmpz_init(a->limbs + 1);
   unsigned long modsize = a->limbs + 1;
   
   fmpz_set_ui(prod, 1L);
   
   fmpz_poly_zero(s);
   fmpz_poly_zero(t);

   unsigned long p = (1L<<(FLINT_BITS-2));
   
   int first = 1;
   
   for (;;) {
      p = z_nextprime(p);

      unsigned long R = fmpz_mod_ui(r, p);

      if ((fmpz_mod_ui(_fmpz_poly_lead(a), p) == 0L) || (fmpz_mod_ui(_fmpz_poly_lead(b), p) == 0L)
         || (R == 0))
         continue;

      zmod_poly_t D, S, T, A, B;
      zmod_poly_init(D, p);
      zmod_poly_init(S, p);
      zmod_poly_init(T, p);
      zmod_poly_init(A, p);
      zmod_poly_init(B, p);
      
      fmpz_poly_to_zmod_poly(A, a);
      fmpz_poly_to_zmod_poly(B, b);

      if (stabilised) {
         fmpz_poly_to_zmod_poly(S, s);
         fmpz_poly_to_zmod_poly(T, t);
         zmod_poly_t t1, t2;
         zmod_poly_init(t1, p);
         zmod_poly_init(t2, p);
         zmod_poly_mul(t1, A, S); 
         zmod_poly_mul(t2, B, T);
         zmod_poly_add(t1, t1, t2);
         
         if ((t1->length == 1) && (t1->coeffs[0] == R))
            fmpz_mul_ui(prod, prod, p);
         else
            stabilised = 0;

         zmod_poly_clear(t1);
		 zmod_poly_clear(t2);
		 if (prod[0] >= modsize - 1) 
         {
            prod = fmpz_realloc(prod, modsize+8);
            modsize += 8;
         }   
      }

      if (!stabilised) {
         zmod_poly_xgcd(D, S, T, A, B);
         zmod_poly_scalar_mul(S, S, R);
         zmod_poly_scalar_mul(T, T, R);
   
         if (first)
         {
            zmod_poly_to_fmpz_poly(s, S);
            zmod_poly_to_fmpz_poly(t, T);
            fmpz_set_ui(prod, p);
            stabilised = 1;
            first = 0;
         } else
         {
            if (prod[0] >= modsize - 2) 
            {
               modsize += 8;
            }

            fmpz_t tmp = fmpz_init(modsize);
         
            int S_stabilised = fmpz_poly_CRT(s, s, S, tmp, prod);
            int T_stabilised = fmpz_poly_CRT(t, t, T, tmp, prod);

            fmpz_clear(prod);
            prod = tmp;

            stabilised = S_stabilised && T_stabilised;
         }
      }

      if (stabilised) {
         unsigned long bound1 = FLINT_BIT_COUNT(FLINT_MIN(a->length, s->length)) 
                      + FLINT_ABS(_fmpz_poly_max_bits(a)) + FLINT_ABS(_fmpz_poly_max_bits(s));
         unsigned long bound2 = FLINT_BIT_COUNT(FLINT_MIN(b->length, t->length)) 
                      + FLINT_ABS(_fmpz_poly_max_bits(b)) + FLINT_ABS(_fmpz_poly_max_bits(t));

         unsigned long bound = 4 + FLINT_MAX(fmpz_bits(r), FLINT_MAX(bound1, bound2));

         if (modsize < bound/FLINT_BITS + 2) 
         {
            prod = fmpz_realloc(prod, bound/FLINT_BITS + 2);
            modsize = bound/FLINT_BITS + 2;
         }

         if (fmpz_bits(prod) > bound)
		 {
            zmod_poly_clear(D);
			zmod_poly_clear(S);
			zmod_poly_clear(T);
			zmod_poly_clear(A);
			zmod_poly_clear(B);
			break;
		 }
      }
	  zmod_poly_clear(D);
	  zmod_poly_clear(S);
	  zmod_poly_clear(T);
	  zmod_poly_clear(A);
	  zmod_poly_clear(B);		
   }
   fmpz_clear(prod);
}

/****************************************************************************

   Resultant

****************************************************************************/

void fmpz_poly_2norm(fmpz_t norm, const fmpz_poly_t pol)
{
   if (pol->length == 0)
   {
      norm[0] = 0L;
      return;
   }

   fmpz_t sqr = fmpz_init(2*pol->limbs);
   fmpz_t sum = fmpz_init(2*pol->limbs + 1);
   fmpz_t temp = fmpz_init(2*pol->limbs + 1);
   
   unsigned long size = pol->limbs+1;
   fmpz_t coeff = pol->coeffs;

   fmpz_set_ui(sum, 0L);
   for (unsigned long i = 0; i < pol->length; i++)
   {
      fmpz_mul(sqr, coeff, coeff);
      fmpz_add(sum, sum, sqr);
      coeff += size;
   }
   
   fmpz_sqrtrem(norm, temp, sum);
   if (temp[0]) fmpz_add_ui(norm, norm, 1L);

   fmpz_clear(temp); // release temp
   fmpz_clear(sum); // release sum
   fmpz_clear(sqr); // release sqr
} 

unsigned long fmpz_poly_resultant_bound(fmpz_poly_t a, fmpz_poly_t b)
{
   if (b->length == 0) return 0;
   if (a->length == 0) return 0;
   fmpz_t t1, t2, tt;
   t1 = fmpz_init(b->length*(a->limbs+1));
   t2 = fmpz_init(a->length*(b->limbs+1));
   fmpz_poly_2norm(t1, a);
   fmpz_poly_2norm(t2, b);
   fmpz_pow_ui(t1, t1, b->length - 1);
   fmpz_pow_ui(t2, t2, a->length - 1);
   tt = fmpz_init(fmpz_size(t1)+fmpz_size(t2));
   fmpz_mul(tt, t1, t2);
   fmpz_clear(t1);
   fmpz_clear(t2);

   unsigned long bound = fmpz_bits(tt);
   fmpz_clear(tt);

   return bound;
}

void fmpz_poly_resultant(fmpz_t res, fmpz_poly_t a, fmpz_poly_t b)
{
   if ((a->length == 0) || (b->length == 0)) 
   {
      res[0] = 0L;
      return;
   }

   unsigned long bound = fmpz_poly_resultant_bound(a, b)+2;
   
   fmpz_t prod = fmpz_init(bound/FLINT_BITS + 2);
   
   fmpz_set_ui(res, 0L);
   fmpz_set_ui(prod, 1L);

   unsigned long p = (1L<<(FLINT_BITS-2));
   unsigned long t;

   int first = 1;

   for (;;) {
      if (fmpz_bits(prod) > bound)
         break;

      p = z_nextprime(p);

      if ((fmpz_mod_ui(_fmpz_poly_lead(a), p) == 0L) || (fmpz_mod_ui(_fmpz_poly_lead(b), p) == 0L))
         continue;

      zmod_poly_t A, B;
      zmod_poly_init(A, p);
      zmod_poly_init(B, p);

      fmpz_poly_to_zmod_poly(A, a);
      fmpz_poly_to_zmod_poly(B, b);

      t = zmod_poly_resultant(A, B);
      
      if (first)
      {
         fmpz_set_ui(prod, p);
         fmpz_set_ui(res, t);
         first = 0;
      } else
      {
         unsigned long c = fmpz_mod_ui(prod, p);
         c = z_invert(c, p);
         fmpz_CRT_ui2_precomp(res, res, prod, t, p, c, A->p_inv);
      
         fmpz_mul_ui(prod, prod, p);
      }

      zmod_poly_clear(A);
      zmod_poly_clear(B);
   }

   fmpz_t proddiv2 = fmpz_init(prod[0]);
   fmpz_div_2exp(proddiv2, prod, 1);
   if (fmpz_cmpabs(res, proddiv2) > 0L) 
      fmpz_sub(res, res, prod);
   fmpz_clear(proddiv2);
   fmpz_clear(prod);
}

/****************************************************************************

   Derivative

****************************************************************************/

void fmpz_poly_derivative(fmpz_poly_t der, fmpz_poly_t poly)
{
	if (poly->length <= 1)
	{
		fmpz_poly_zero(der);
		return;
	}
	
	ulong bits = FLINT_ABS(_fmpz_poly_max_bits(poly));
	bits += ceil_log2(poly->length);
	fmpz_poly_fit_length(der, poly->length - 1);
   fmpz_poly_fit_limbs(der, (bits - 1)/FLINT_BITS + 1);

	ulong size_d = der->limbs + 1;
	ulong size_p = poly->limbs + 1;
   fmpz_t coeff_d = der->coeffs;
	fmpz_t coeff_p = poly->coeffs + size_p;
	
	for (ulong i = 0; i < poly->length - 1; i++)
	{
		fmpz_mul_ui(coeff_d, coeff_p, i + 1);
		coeff_d += size_d;
		coeff_p += size_p;
	}

	der->length = poly->length - 1;
}

/****************************************************************************

   Evaluation

****************************************************************************/

/*
    Treat the n coefficients starting at the given coefficient of poly as a polynomial
	 and evaluate this polynomial at the given value
*/

void fmpz_poly_evaluate_horner_range(fmpz_t output, fmpz_poly_t poly, fmpz_t val, ulong start, ulong n)
{
    if ((n == 0) || (val[0] == 0L))
	 {
		 output[0] = 0L;
		 return;
	 }

	 if (n == 1)
	 {
		 fmpz_set(output, poly->coeffs + start*(poly->limbs+1));
		 return;
	 }

	 ulong size_p = poly->limbs + 1;
	 fmpz_t coeff_p = poly->coeffs + (start + n - 1)* size_p;
	 
	 fmpz_set(output, coeff_p); // set output to top coefficient

	 if (val[0] == 1L)
	 {
		 if (val[1] == 1L) // special case, value == 1
	    {
		    for (long i = n - 2; i >= 0L; i--)
		    {
             coeff_p -= size_p;
			    fmpz_add(output, output, coeff_p);
		    }

		    return;
	    }

		 // value is a positive single limb
       ulong value = val[1];
		 
		 for (long i = n - 2; i >= 0L; i--)
		 {
          coeff_p -= size_p;
			 fmpz_mul_ui(output, output, value);
			 fmpz_add(output, output, coeff_p);
		 }

		 return;
	 }

	 if (val[0] == -1L)
	 {
		 if (val[1] == 1L) // special case, value == -1
	    {
		    long i;

		    for (i = n - 2; i >= 1L; i-=2)
		    {
             coeff_p -= size_p;
			    fmpz_sub(output, output, coeff_p);
             coeff_p -= size_p;
			    fmpz_add(output, output, coeff_p);
		    }

          if (i == 0) 
		    {
			    coeff_p -= size_p;
			    fmpz_sub(output, output, coeff_p);
			    output[0] = -output[0];
		    }
		    return;
	    }

		 // value is a negative single limb
       ulong value = val[1];
		 long i;

		 for (i = n - 2; i >= 1L; i-=2)
		 {
          coeff_p -= size_p;
			 fmpz_mul_ui(output, output, value);
			 fmpz_sub(output, output, coeff_p);
          coeff_p -= size_p;
			 fmpz_mul_ui(output, output, value);
			 fmpz_add(output, output, coeff_p);
		 }
		 
		 if (i == 0)
		 {
			 coeff_p -= size_p;
			 fmpz_mul_ui(output, output, value);
			 fmpz_sub(output, output, coeff_p);
          output[0] = -output[0];
		 }

		 return;
	 }
	 
	 // value is more than one limb
	 for (long i = n - 2; i >= 0L; i--)
    {
       coeff_p -= size_p;
		 fmpz_mul(output, output, val);
		 fmpz_add(output, output, coeff_p);
	 }
}

void fmpz_poly_evaluate_divconquer(fmpz_t output, fmpz_poly_t poly, fmpz_t val)
{
	if ((poly->length == 0) || (val[0] == 0))
	{
		fmpz_set_ui(output, 0L);
		return;
	}

	if (((FLINT_ABS(val[0]) == 1) && (val[1] == 1)) || (poly->length == 2))
	{
		fmpz_poly_evaluate_horner(output, poly, val);
		return;
	}

	if (poly->length == 1)
	{
		fmpz_set(output, poly->coeffs);
		return;
	}

	fmpz_poly_t half, temp;

	ulong val_bits = fmpz_bits(val);
	ulong bits = FLINT_ABS(fmpz_poly_max_bits(poly)) + val_bits;
	fmpz_poly_init2(temp, (poly->length + 1)/2, bits/FLINT_BITS + 1);

	ulong size_t = poly->limbs + 1;
   ulong size_h = temp->limbs + 1;
	fmpz_t coeff_t = poly->coeffs;
	fmpz_t coeff_h = temp->coeffs;

	for (ulong i = 0; i < poly->length/2; i++)
	{
		fmpz_mul(coeff_h, coeff_t + size_t, val);
		fmpz_add(coeff_h, coeff_h, coeff_t);
		coeff_t += (size_t*2);
		coeff_h += size_h;
	}

	if (poly->length & 1) fmpz_set(coeff_h, coeff_t);
	temp->length = (poly->length + 1)/2;

	ulong log_iter = 1;
	fmpz_t val_pow = fmpz_init((val_bits*(1L<<ceil_log2(poly->length)))/FLINT_BITS + 1);
	fmpz_mul(val_pow, val, val);

	while (temp->length > 2)
	{
		bits = FLINT_ABS(fmpz_poly_max_bits(temp)) + (val_bits<<log_iter);
		fmpz_poly_init2(half, (temp->length + 1)/2, bits/FLINT_BITS + 1);

		size_t = temp->limbs + 1;
      size_h = half->limbs + 1;
	   coeff_t = temp->coeffs;
	   coeff_h = half->coeffs;

		for (ulong i = 0; i < temp->length/2; i++)
	   {
		   fmpz_mul(coeff_h, coeff_t + size_t, val_pow);
		   fmpz_add(coeff_h, coeff_h, coeff_t);
		   coeff_t += (size_t*2);
		   coeff_h += size_h;
	   }

		half->length = (temp->length + 1)/2;

		fmpz_mul(val_pow, val_pow, val_pow);

		if (temp->length & 1) fmpz_set(coeff_h, coeff_t);
	
		fmpz_poly_swap(half, temp);
		fmpz_poly_clear(half);
		log_iter++;
	}
   
	fmpz_mul(output, temp->coeffs + temp->limbs + 1, val_pow);
   fmpz_add(output, output, temp->coeffs);
		
	fmpz_poly_clear(temp);
	fmpz_clear(val_pow);
}

void fmpz_poly_evaluate(fmpz_t output, fmpz_poly_t poly, fmpz_t value)
{
   fmpz_t val;
	
	if ((poly->length == 0) || (value[0] == 0)) 
	{
      output[0] = 0L;
		return;
	}

	if (output == value)
	{
		val = fmpz_init(FLINT_ABS(value[0]));
		fmpz_set(val, value);
	} else val = value;
	
	if (((FLINT_ABS(val[0]) == 1) && (val[1] == 1)) || (poly->length == 2))
	{
		fmpz_poly_evaluate_horner(output, poly, val);
		if (output == value) fmpz_clear(val);
	   return;
	}

	if (poly->length == 1)
	{
		fmpz_set(output, poly->coeffs);
		if (output == value) fmpz_clear(val);
	   return;
	}

	ulong eval_length;
	ulong val_bits = fmpz_bits(val);
	if (val_bits <= 6) eval_length = 256;
	else if (val_bits <= 12) eval_length = 128;
	else if (val_bits <= 128) eval_length = 64;
	else 
	{
		fmpz_poly_evaluate_divconquer(output, poly, val);
		if (output == value) fmpz_clear(val);
	   return;
	}

	ulong bits = FLINT_ABS(fmpz_poly_max_bits(poly)) + val_bits*eval_length;
	ulong short_length = (poly->length - 1)/eval_length + 1;

	fmpz_poly_t temp;
	fmpz_poly_init2(temp, short_length, bits/FLINT_BITS + 1);

	fmpz_t coeff_t = temp->coeffs;
	ulong size_t = temp->limbs + 1;

	long i;
	for (i = 0; i < short_length - 1; i++)
	{
		fmpz_poly_evaluate_horner_range(coeff_t, poly, val, i*eval_length, eval_length);
		coeff_t += size_t;
	}
   if (short_length == 1)
		fmpz_poly_evaluate_horner(output, poly, val);
	else
	{
		fmpz_poly_evaluate_horner_range(coeff_t, poly, val, i*eval_length, poly->length - i*eval_length);

	   temp->length = short_length;

	   fmpz_t val_pow = fmpz_init((val_bits*eval_length - 1)/FLINT_BITS + 1);
	   fmpz_pow_ui(val_pow, val, eval_length);

	   fmpz_poly_evaluate_divconquer(output, temp, val_pow);

	   fmpz_clear(val_pow);
	}
	
	fmpz_poly_clear(temp);
	if (output == value) fmpz_clear(val);
}

/****************************************************************************

   Composition

****************************************************************************/

void fmpz_poly_compose_horner_range(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val, ulong start, ulong n)
{
	if (n == 0)
	{
		fmpz_poly_zero(output);
		return;
	}
	
	ulong size_p = poly->limbs + 1;
   fmpz_t coeff_p = poly->coeffs + (start + n - 1)* size_p;
   fmpz_poly_zero(output);
	fmpz_poly_set_coeff_fmpz(output, 0, coeff_p);
   fmpz_poly_fit_length(output, 1);
	fmpz_poly_fit_limbs(output, size_p);
	for (long i = n - 2; i >= 0L; i--)
   {
      coeff_p -= size_p;
		fmpz_poly_mul(output, output, val);
		if (output->length) fmpz_poly_fit_limbs(output, FLINT_ABS(output->coeffs[0]) + 1);
	   else 
		{
			output->coeffs[0] = 0;
			output->length = 1;
		}
		fmpz_add(output->coeffs, output->coeffs, coeff_p);
		_fmpz_poly_normalise(output);
	}
}

void fmpz_poly_compose_divconquer(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val)
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

	if (poly->length == 2)
	{
		fmpz_poly_compose_horner(output, poly, val);
		return;
	}

	fmpz_poly_t * half = (fmpz_poly_t *) flint_heap_alloc(((poly->length + 1)/2)*sizeof(fmpz_poly_t));
	fmpz_poly_t * temp  = (fmpz_poly_t *) flint_heap_alloc(((poly->length + 1)/2)*sizeof(fmpz_poly_t));

	for (ulong i = 0; i < (poly->length + 1)/2; i++) fmpz_poly_init(temp[i]);
   for (ulong i = 0; i < (poly->length + 1)/2; i++) fmpz_poly_init(half[i]);

	ulong size_t = poly->limbs + 1;
   fmpz_t coeff_t = poly->coeffs;
	
	ulong i;
	for (i = 0; i < poly->length/2; i++)
	{
		fmpz_poly_scalar_mul_fmpz(temp[i], val, coeff_t + size_t);
		fmpz_poly_fit_limbs(temp[i], size_t);
		fmpz_poly_fit_length(temp[i], 1);
		if (temp[i]->length) fmpz_poly_fit_limbs(temp[i], FLINT_ABS(temp[i]->coeffs[0])+1);
	   else 
		{
			temp[i]->coeffs[0] = 0;
			temp[i]->length = 1;
		}
		fmpz_add(temp[i]->coeffs, temp[i]->coeffs, coeff_t);
		_fmpz_poly_normalise(temp[i]);
		coeff_t += (size_t*2);
	}

	if (poly->length & 1) 
	{
		fmpz_poly_zero(temp[i]);
		fmpz_poly_set_coeff_fmpz(temp[i], 0, coeff_t);
	}

   ulong length = (poly->length + 1)/2;
   
	fmpz_poly_t val_pow;
	fmpz_poly_init(val_pow);
	fmpz_poly_mul(val_pow, val, val);

	while (length > 2)
	{
		ulong i;
		for (i = 0; i < length/2; i++)
	   {
		   fmpz_poly_mul(half[i], temp[2*i+1], val_pow);
		   fmpz_poly_add(half[i], half[i], temp[2*i]);
	   }

		if (length & 1) fmpz_poly_set(half[i], temp[2*i]);
	
		fmpz_poly_mul(val_pow, val_pow, val_pow);

		for (ulong i = 0; i < length; i++) fmpz_poly_swap(half[i], temp[i]);

		length = (length + 1)/2;
	}
   
	fmpz_poly_mul(output, temp[1], val_pow);
   fmpz_poly_add(output, output, temp[0]);
		
	for (ulong i = 0; i < (poly->length + 1)/2; i++) fmpz_poly_clear(temp[i]);
   for (ulong i = 0; i < (poly->length + 1)/2; i++) fmpz_poly_clear(half[i]);

	fmpz_poly_clear(val_pow);

	flint_heap_free(temp);
	flint_heap_free(half);
}

void fmpz_poly_array_compose_divconquer(fmpz_poly_t output, fmpz_poly_t * poly, ulong length_in, fmpz_poly_t val)
{
	if (length_in == 0)
	{
		fmpz_poly_zero(output);
		return;
	}
	
	if (length_in == 1)
	{
		fmpz_poly_set(output, poly[0]);
		return;
	}

	if (length_in == 2)
	{
		fmpz_poly_mul(output, poly[1], val);
		fmpz_poly_add(output, output, poly[0]);
		return;
	}

	ulong length = length_in;
	
	fmpz_poly_t * half = (fmpz_poly_t *) flint_heap_alloc(((length + 1)/2)*sizeof(fmpz_poly_t));
	fmpz_poly_t * temp  = (fmpz_poly_t *) flint_heap_alloc(((length + 1)/2)*sizeof(fmpz_poly_t));

	for (ulong i = 0; i < (length + 1)/2; i++) fmpz_poly_init(temp[i]);
   for (ulong i = 0; i < (length + 1)/2; i++) fmpz_poly_init(half[i]);

	ulong i;
	for (i = 0; i < length/2; i++)
	{
		fmpz_poly_mul(temp[i], val, poly[2*i+1]);
		fmpz_poly_add(temp[i], temp[i], poly[2*i]);
	}

	if (length & 1) 
	   fmpz_poly_set(temp[i], poly[2*i]);

   length = (length + 1)/2;

	fmpz_poly_t val_pow;
	fmpz_poly_init(val_pow);
	fmpz_poly_mul(val_pow, val, val);

	while (length > 2)
	{
		ulong i;
		for (i = 0; i < length/2; i++)
	   {
		   fmpz_poly_mul(half[i], temp[2*i+1], val_pow);
		   fmpz_poly_add(half[i], half[i], temp[2*i]);
	   }

		if (length & 1) fmpz_poly_set(half[i], temp[2*i]);
	
		fmpz_poly_mul(val_pow, val_pow, val_pow);

		for (ulong i = 0; i < (length+1)/2; i++) fmpz_poly_swap(half[i], temp[i]);

		length = (length + 1)/2;
	}
   
	fmpz_poly_mul(output, temp[1], val_pow);
   fmpz_poly_add(output, output, temp[0]);
		
	for (ulong i = 0; i < (length_in + 1)/2; i++) fmpz_poly_clear(temp[i]);
   for (ulong i = 0; i < (length_in + 1)/2; i++) fmpz_poly_clear(half[i]);

	fmpz_poly_clear(val_pow);

	flint_heap_free(temp);
	flint_heap_free(half);
}

void fmpz_poly_compose_combined(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val)
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

	if (poly->length <= 4) 
	{
		fmpz_poly_compose_divconquer(output, poly, val);
		return;
	}
   ulong eval_length = 4;

	ulong short_length = (poly->length - 1)/eval_length + 1;

	fmpz_poly_t * temp  = (fmpz_poly_t *) flint_heap_alloc(short_length*sizeof(fmpz_poly_t));

	for (ulong i = 0; i < short_length; i++) fmpz_poly_init(temp[i]);
   
	long i;
	for (i = 0; i < short_length - 1; i++)
	{
		fmpz_poly_compose_horner_range(temp[i], poly, val, i*eval_length, eval_length);
	}

   if (short_length == 1)
		fmpz_poly_compose_horner_range(output, poly, val, 0, poly->length);
	else
	{
		fmpz_poly_compose_horner_range(temp[i], poly, val, i*eval_length, poly->length - i*eval_length);

		fmpz_poly_t val_pow;
		fmpz_poly_init(val_pow);
		
		fmpz_poly_power(val_pow, val, eval_length);

	   fmpz_poly_array_compose_divconquer(output, temp, short_length, val_pow);

	   fmpz_poly_clear(val_pow);
	}
	
	for (ulong i = 0; i < short_length; i++) fmpz_poly_clear(temp[i]);
}

void fmpz_poly_compose(fmpz_poly_t output, fmpz_poly_t poly, fmpz_poly_t val)
{
	fmpz_poly_t out;

	if ((output == poly) || (output == val))
	{
		fmpz_poly_init(out);
		fmpz_poly_compose_combined(out, poly, val);
		fmpz_poly_swap(out, output);
		fmpz_poly_clear(out);
	} else 
      fmpz_poly_compose_combined(output, poly, val);
}
