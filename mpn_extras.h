/*============================================================================

    Copyright (C) 2007, William Hart, David Harvey

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
#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

#include "flint.h"
#include "ZmodF_poly.h"

/*============================================================================

"mpn-wannabe" code.

These are functions that I wish were in GMP's mpn layer.

=============================================================================*/

/*
Computes the negation of a multiple-precision integer in 2's complement.

Input is count limbs stored at src. Output is stored at dest.

src and dest can be the same buffer. If they're not, they should be disjoint.

todo: currently this code will make only 1 pass over the data, EXCEPT in the
      case where all limbs are zero, in which case it will make two passes.
      FIX THIS!

todo: try writing another version that makes a block of zeroes and then
      uses mpn_sub_n repeatedly. This could be faster, if GMP's assembler
      is better than what gcc can come up with.

todo: consider writing this in assembly

todo: write test function for this

todo: consider using GMP's mpn_com_n (undocumented)

*/

static inline void F_mpn_negate(mp_limb_t* dest, mp_limb_t* src, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = ~src[i];
   mpn_add_1(dest, dest, count, 1);
}

/*
Copies a bunch of limbs from one buffer to another.

Input is count limbs stored at src. Output is stored at dest.

src and dest can be the same buffer. If they're not, they should be disjoint.

todo: it's completely insane that we're not using memcpy. But memcpy seems
      to have crazy overhead and is slow!! Why is this?

todo: GMP has code to do limb copying. Clearly memcpy wasn't good enough for
      them either. Work out how to use their code. It's not a documented
      interface, so some hackishness may be necessary.
 */
static inline void F_mpn_copy(mp_limb_t* dest, const mp_limb_t* src, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
   {
      dest[i] = src[i];
   }
}

static inline void F_mpn_copy_forward(mp_limb_t* dest, const mp_limb_t* src, unsigned long count)
{
   for (long i = 0; i < count; i++)
   {
      dest[i] = src[i];
   }
}


/*
Sets a bunch of limbs to zero.

todo: why does memset have so much overhead????!!?
 */
static inline void F_mpn_clear(mp_limb_t* dest, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = 0;
}

/*
Sets a bunch of limbs to 0xfff....

todo: why does memset have so much overhead????!!?
 */
static inline void F_mpn_set(mp_limb_t* dest, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = (mp_limb_t)(-1L);
}



mp_limb_t F_mpn_divmod_1_preinv(mp_limb_t * qp, mp_limb_t * up, 
             unsigned long un, mp_limb_t d, mp_limb_t dinv, unsigned long norm);
             
mp_limb_t F_mpn_addmul(mp_limb_t * rp, mp_limb_t * s1p, unsigned long s1n, 
                                      mp_limb_t * s2p, unsigned long s2n);
                                      
static inline
void F_mpn_printx(mp_limb_t * mpn, unsigned long count)
{
   if (count) 
      for (unsigned long i = 0; i < count; i++)
         printf("%lx ", mpn[i]);
}
                                      
/* 
   Large integer multiplication
*/

typedef enum {FFT_PRE, KAR_PRE} precomp_type;
 
typedef struct
{
   precomp_type type;
   ZmodF_poly_p poly;
   unsigned long length;
   unsigned long length2;
   unsigned long coeff_limbs;
   unsigned long limbs1;
   unsigned long limbs2;
   unsigned long msl_bits;
} F_mpn_precomp_s;

typedef F_mpn_precomp_s F_mpn_precomp_t[1]; 

void F_mpn_FFT_split_bits(ZmodF_poly_t poly, mp_limb_t * limbs, unsigned long total_limbs,
                               unsigned long bits, unsigned long output_limbs);
                               
void F_mpn_FFT_combine_bits(mp_limb_t * res, ZmodF_poly_t poly, unsigned long bits, 
                             unsigned long output_limbs, unsigned long total_limbs);

mp_limb_t __F_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2, 
                                      unsigned long twk);

mp_limb_t F_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2);
                                      
mp_limb_t F_mpn_mul_trunc(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                        mp_limb_t * data2, unsigned long limbs2, unsigned long trunc);

void F_mpn_mul_precomp_init(F_mpn_precomp_t precomp, mp_limb_t * data1, unsigned long limbs1, unsigned long limbs2);

void F_mpn_mul_precomp_clear(F_mpn_precomp_t precomp);

mp_limb_t F_mpn_mul_precomp(mp_limb_t * res, mp_limb_t * data2, unsigned long limbs2, F_mpn_precomp_t precomp);
      

#endif

