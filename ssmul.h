/****************************************************************************

 ssmul.h

 Copyright (C) 2006, David Harvey and William Hart

****************************************************************************/


#ifndef SSMUL_H
#define SSMUL_H

#include "Zvec.h"

void SSMul(Zvec outpoly, Zvec poly1, Zvec poly2, unsigned long coeff_bits, int sign);



/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*****************************************************************************

Everything below this line is INTERNAL.

Do not rely on it to exist in future versions of FLINT!

It is exported for purposes of test code only.

=============================================================================*/

void ssmul_convert_in_bytes(mp_limb_t* array, mpz_t* coeffs, 
            const unsigned long bundle, const unsigned long coeff_bytes, const unsigned long total_limbs);
            
void ssmul_convert_out_bytes(mpz_t* output, mp_limb_t* array, 
               unsigned long bundle, unsigned long coeff_bytes);

void rotate_right_bits(mp_limb_t* dest, mp_limb_t* src,
                       unsigned long s, unsigned long n);
void reduce_mod_p_exact(mp_limb_t* x, unsigned long n);
void rotate_mod_p_limbs(mp_limb_t* dest, mp_limb_t* src,
                        unsigned long s, unsigned long n);
void rotate_mod_p_bits(mp_limb_t* dest, mp_limb_t* src,
                       unsigned long s, unsigned long n);

void fft_butterfly_limbs(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                         unsigned long s, unsigned long n);
void fft_butterfly_bits(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                        unsigned long s, unsigned long n);
void ifft_butterfly_limbs(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                          unsigned long s, unsigned long n);
void ifft_butterfly_bits(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                         unsigned long s, unsigned long n);

void fft_main(mp_limb_t** start, unsigned long skip,
              unsigned long start_r, unsigned long skip_r,
              unsigned long depth, mp_limb_t** scratch,
              unsigned long n, int first, int crossover);
void ifft_main(mp_limb_t** start, unsigned long skip,
               unsigned long start_r, unsigned long skip_r,
               unsigned long depth, mp_limb_t** scratch,
               unsigned long n, int crossover);

#endif
