/*
   invert.c:  series inversion
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "zn_poly_internal.h"


/*
   Extends an approximation of a power series reciprocal from length n1 to
   length n1 + n2. Must have 1 <= n2 <= n1.
   
   op[0, n1 + n2) represents an input power series f.
   approx[0, n1) should be the first n1 coeffs of the reciprocal of f.
   
   This function computes the next n2 coefficients of the reciprocal of f,
   and stores them at res[0, n2).

   res may not overlap op or approx.
*/
#define zn_array_invert_extend \
    ZNP_zn_array_invert_extend
void
zn_array_invert_extend (ulong* res, const ulong* approx, const ulong* op,
                        size_t n1, size_t n2, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   // The algorithm is basically newton iteration, inspired partly by the
   // algorithm in [HZ04], as follows.

   // Let f be the input series, of length n1 + n2.
   // Let g be the current approximation to 1/f, of length n1.
   
   // By newton iteration, (2*g - g*g*f) is a length n1 + n2 approximation
   // to 1/f. Therefore the output of this function should be terms
   // [n1, n1 + n2) of -g*g*f.
   
   // We have g*f = 1 + h*x^n1 + O(x^(n1 + n2)), where h has length n2,
   // i.e. h consists of terms [n1, n1 + n2) of g*f. Therefore h may be
   // recovered as the middle product of f[1, n1 + n2) and g[0, n1).
   
   // Then g*g*f = g + g*h*x^n1 + O(x^(n1 + n2)). Since g has length
   // n1, the output is (the negative of) the first n2 coefficients of g*h.


   // Compute h, put it in res[0, n2).
   zn_array_mulmid (res, op + 1, n1 + n2 - 1, approx, n1, mod);
   
   // Compute g * h, put it into a scratch buffer.
   ZNP_FASTALLOC (temp, ulong, 6624, n1 + n2 - 1);
   zn_array_mul (temp, approx, n1, res, n2, mod);
   
   // Negate the first n2 coefficients of g * h into the output buffer.
   zn_array_neg (res, temp, n2, mod);
   ZNP_FASTFREE (temp);
}



/*
   Same as zn_array_invert_extend(), but uses Schonhage/Nussbaumer FFT for
   the middle product and product.
   
   Modulus must be odd.

   res may overlap op or approx.
*/
#define zn_array_invert_extend_fft \
    ZNP_zn_array_invert_extend_fft
void
zn_array_invert_extend_fft (ulong* res, const ulong* approx, const ulong* op,
                            size_t n1, size_t n2, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (mod->m & 1);

   // The algorithm here is the same as in zn_array_invert_extend(), except
   // that we work with the FFTs directly. This allows us to save one FFT,
   // since we use the FFT of g in both the middle product step and the
   // product step.

   // Determine FFT parameters for computing h = middle product of
   // f[1, n1 + n2) and g[0, n1). (These parameters will also work for the
   // subsequent product g * h.)
   unsigned lgK, lgM; 
   ulong m1, m2, m3, p;

   mulmid_fft_params (&lgK, &lgM, &m3, &m1, &p, n1 + n2 - 1, n1);
   m2 = m3 - m1 + 1;

   // We now have
   //     m1 = ceil(n1 / (M/2))
   //        = (n1 + p - 1) / (M/2).
   // Therefore
   //     m3 = ceil((n1 + n2 - 1 + p) / (M/2))
   //        = ceil(n2 / (M/2)) + (n1 + p - 1) / (M/2)
   // and
   //     m2 = ceil(n2 / (M/2)) + 1.
   
   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ptrdiff_t skip = M + 1;

   pmfvec_t vec1, vec2;
   pmfvec_init (vec1, lgK, skip, lgM, mod);
   pmfvec_init (vec2, lgK, skip, lgM, mod);

   // Find scaling factor that needs to be applied to both of the products
   // below; takes into account the fudge from the pointwise multiplies, and
   // the division by 2^lgK coming from the FFTs.
   ulong x = pmfvec_mul_fudge (lgM, 0, mod);
   x = zn_mod_mul (x, zn_mod_pow2 (-lgK, mod), mod);

   // Split g[0, n1) into m1 coefficients, apply scaling factor, and compute
   // m3 fourier coefficients, written to vec2.
   fft_split (vec2, approx, n1, 0, x, 0);
   pmfvec_fft (vec2, m3, m1, 0);

   // Split f[1, n1 + n2) into m3 coefficients (in reversed order, with
   // appropriate zero-padding), and compute transposed IFFT of length m3,
   // written to vec1.
   pmfvec_reverse (vec1, m3);
   fft_split (vec1, op + 1, n1 + n2 - 1, p, 1, 0);
   pmfvec_reverse (vec1, m3);
   pmfvec_tpifft (vec1, m3, 0, m3, 0);

   // Pointwise multiply the above FFT and transposed IFFT, into vec1.
   pmfvec_mul (vec1, vec1, vec2, m3, 0);
   
   // Transposed FFT vec1, obtaining m2 coefficients, then reverse and combine.
   pmfvec_tpfft (vec1, m3, m2, 0);
   pmfvec_reverse (vec1, m2);
   fft_combine (res, n2, vec1, m2, 1);
   pmfvec_reverse (vec1, m2);
   
   // At this stage we have obtained the polynomial h in res[0, n2).
   // Now we must compute h * g.
   
   // Split h[0, n2) into m2 - 1 coefficients, and compute m3 - 1 fourier
   // coefficients in vec1. For the splitting step, we set the bias to M,
   // which effectively negates everything, so we're really computing the FFT
   // of -h.
   fft_split (vec1, res, n2, 0, 1, M);
   pmfvec_fft (vec1, m3 - 1, m2 - 1, 0);

   // Pointwise multiply that FFT with the first FFT of g into vec2.
   pmfvec_mul (vec2, vec2, vec1, m3 - 1, 1);
   pmfvec_clear (vec1);

   // IFFT and combine, to obtain the product -h * g. We only need the low n2
   // terms of the product (we throw away the high n1 - 1 terms).
   pmfvec_ifft (vec2, m3 - 1, 0, m3 - 1, 0);
   fft_combine (res, n2, vec2, m3 - 1, 0);
   pmfvec_clear (vec2);
}



void
zn_array_invert (ulong* res, const ulong* op, size_t n, const zn_mod_t mod)
{
   ZNP_ASSERT (n >= 1);
   
   // for now assume input is monic
   ZNP_ASSERT (op[0] == 1);
   
   if (n == 1)
   {
      res[0] = 1;
      return;
   }
   
   size_t half = (n + 1) / 2;    // ceil(n / 2)
   
   // recursively obtain the first half of the output
   zn_array_invert (res, op, half, mod);

   // extend to second half of the output
   if (mod->m & 1)
      zn_array_invert_extend_fft (res + half, res, op, half, n - half, mod);
   else
      zn_array_invert_extend (res + half, res, op, half, n - half, mod);
}


// end of file ****************************************************************
