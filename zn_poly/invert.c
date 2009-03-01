/*
   invert.c:  series inversion
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
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
   Extends an approximation of a power series reciprocal from length len1
   to length len1 + len2. Must have 1 <= len2 <= len1.
   
   op[0, len1 + len2) represents an input power series f.
   approx[0, len1) should be the first len1 coeffs of the reciprocal of f.
   
   This function computes the next len2 coefficients of the reciprocal of f,
   and stores them at res[0, len2).

   _res_ may not overlap _op_ or _approx_.
*/
#define zn_array_invert_extend \
    ZNP_zn_array_invert_extend
void zn_array_invert_extend(
         ulong* res, const ulong* approx, const ulong* op,
         size_t len1, size_t len2, const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   // The algorithm is basically newton iteration, inspired partly by the
   // algorithm in [HZ04], as follows.

   // Let m = len1, k = len2.
   // Let f be the input series, of length m + k.
   // Let g be the current approximation to 1/f, of length m.
   
   // By newton iteration, (2*g - g*g*f) is a length m+k approximation to 1/f.
   // Therefore the output of this function should be terms [m, m+k) of -g*g*f.
   
   // We have g*f = 1 + h*x^m + O(x^(m+k)), where h has length k,
   // i.e. h consists of terms [m, m+k) of g*f. Therefore h may be recovered
   // as the middle product of f[1, m+k) and g[0, m).
   
   // Then g*g*f = g + g*h*x^m + O(x^(m+k)). Since g has length m, the output
   // is (the negative of) the first k coefficients of g*h.


   // Compute h, put it in res[0, k).
   zn_array_midmul(res, op + 1, len1 + len2 - 1, approx, len1, mod);
   
   // Compute g*h, put it into a scratch buffer.
   ZNP_FASTALLOC(temp, ulong, 6624, len1 + len2 - 1);
   zn_array_mul(temp, approx, len1, res, len2, mod);
   
   // Negate the first k coefficients of g*h into the output buffer.
   zn_array_neg(res, temp, len2, mod);
   ZNP_FASTFREE(temp);
}



/*
   Same as zn_array_invert_extend(), but uses Schonhage/Nussbaumer FFT for
   the middle product and product.
   
   Modulus must be odd.

   _res_ may overlap _op_ or _approx_.
*/
#define zn_array_invert_extend_fft \
    ZNP_zn_array_invert_extend_fft
void zn_array_invert_extend_fft(
         ulong* res, const ulong* approx, const ulong* op,
         size_t len1, size_t len2, const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   ZNP_ASSERT(mod->n & 1);

   // The algorithm here is the same as in zn_array_invert_extend(), except
   // that we work with the FFTs directly. This allows us to save one FFT,
   // since we use the FFT of g in both the middle product step and the
   // product step.
   
   // Determine FFT parameters for computing h = middle product of f[1, m+k)
   // and g[0, m). (These parameters will also work for the subsequent
   // product g*h.)
   unsigned lgK, lgM;
   ulong coeffs1, coeffs2, coeffs3, pad;

   midmul_fft_params(&lgK, &lgM, &coeffs3, &coeffs1, &pad,
                     len1 + len2 - 1, len1);

   coeffs2 = coeffs3 - coeffs1 + 1;

   // We now have
   //     coeffs1 = ceil(m / (M/2))
   //             = (m + pad - 1) / (M/2).
   // Therefore
   //     coeffs3 = ceil((m + k - 1 + pad) / (M/2))
   //             = ceil(k / (M/2)) + (m + pad - 1) / (M/2)
   // and
   //     coeffs2 = ceil(k / (M/2)) + 1.
   
   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ptrdiff_t skip = M + 1;

   zn_pmf_vec_t vec1, vec2;
   zn_pmf_vec_init(vec1, lgK, skip, lgM, mod);
   zn_pmf_vec_init(vec2, lgK, skip, lgM, mod);

   // Find scaling factor that needs to be applied to both of the products
   // below; takes into account the fudge from the pointwise multiplies, and
   // the division by 2^lgK coming from the FFTs.
   ulong scale = zn_pmf_vec_mul_get_fudge(lgM, 0, mod);
   scale = zn_mod_mul(scale, zn_mod_pow2(-lgK, mod), mod);

   // Split g[0, m) into coeffs1 coefficients, apply scaling factor, and
   // compute coeffs3 fourier coefficients, written to vec2.
   fft_split(vec2, approx, len1, 0, scale, 0);
   zn_pmf_vec_fft(vec2, coeffs3, coeffs1, 0);

   // Split f[1, m+k) into coeffs3 coefficients (in reversed order, with
   // appropriate zero-padding), and compute transposed IFFT of length
   // coeffs3, written to vec1.
   zn_pmf_vec_reverse(vec1, coeffs3);
   fft_split(vec1, op + 1, len1 + len2 - 1, pad, 1, 0);
   zn_pmf_vec_reverse(vec1, coeffs3);
   zn_pmf_vec_ifft_transposed(vec1, coeffs3, 0, coeffs3, 0);
   
   // Pointwise multiply the above FFT and transposed IFFT, into vec1.
   zn_pmf_vec_mul(vec1, vec1, vec2, coeffs3, 0);
   
   // Transposed FFT vec1, obtaining _coeffs2_ coefficients, then reverse
   // and combine.
   zn_pmf_vec_fft_transposed(vec1, coeffs3, coeffs2, 0);
   zn_pmf_vec_reverse(vec1, coeffs2);
   fft_combine(res, len2, vec1, coeffs2, 1);
   zn_pmf_vec_reverse(vec1, coeffs2);
   
   // At this stage we have obtained the polynomial h in res[0, k).
   // Now we must compute h*g.
   
   // Split h[0, k) into coeffs2 - 1 coefficients, and compute coeffs3 - 1
   // fourier coefficients in vec1. For the splitting step, we set the bias
   // to M, which effectively negates everything, so we're really computing
   // the FFT of -h.
   fft_split(vec1, res, len2, 0, 1, M);
   zn_pmf_vec_fft(vec1, coeffs3 - 1, coeffs2 - 1, 0);

   // Pointwise multiply that FFT with the first FFT of g into vec2.
   zn_pmf_vec_mul(vec2, vec2, vec1, coeffs3 - 1, 1);
   zn_pmf_vec_clear(vec1);

   // IFFT and combine, to obtain product -h*g. We only need the low k terms
   // of the product (we throw away the high m - 1 terms).
   zn_pmf_vec_ifft(vec2, coeffs3 - 1, 0, coeffs3 - 1, 0);
   fft_combine(res, len2, vec2, coeffs3 - 1, 0);
   zn_pmf_vec_clear(vec2);
}



void zn_array_invert(ulong* res, const ulong* op, size_t len,
                     const zn_mod_t mod)
{
   ZNP_ASSERT(len >= 1);
   
   // for now assume input is monic
   ZNP_ASSERT(op[0] == 1);
   
   if (len == 1)
   {
      res[0] = 1;
      return;
   }
   
   size_t half = (len + 1) / 2;    // ceil(length / 2)
   
   // recursively obtain the first half of the output
   zn_array_invert(res, op, half, mod);

   // extend to second half of the output
   if (mod->n & 1)
      zn_array_invert_extend_fft(res + half, res, op, half, len - half, mod);
   else
      zn_array_invert_extend(res + half, res, op, half, len - half, mod);
}


// end of file ****************************************************************
