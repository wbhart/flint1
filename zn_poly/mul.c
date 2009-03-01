/*
   mul.c:  polynomial multiplication
   
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


ulong _zn_array_mul_get_fudge(size_t len1, size_t len2, int squaring,
                              const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   if (!(mod->n & 1))
      // no fudge if the modulus is even.
      return 1;

   tuning_info_t* i = &tuning_info[mod->bits];

   if (!squaring)
   {
      if (len2 < i->mul_KS2_crossover  ||  len2 < i->mul_KS4_crossover ||
          len2 < i->mul_fft_crossover)
         // fudge is -B
         return mod->n - mod->B;
   }
   else
   {
      if (len2 < i->sqr_KS2_crossover  ||  len2 < i->sqr_KS4_crossover ||
          len2 < i->sqr_fft_crossover)
         // fudge is -B
         return mod->n - mod->B;
   }

   // return whatever fudge is used by the fft multiplication code
   return zn_array_mul_fft_get_fudge(len1, len2, squaring, mod);
}


void _zn_array_mul(ulong* res, const ulong* op1, size_t len1,
                   const ulong* op2, size_t len2,
                   int fastred, const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   // we can use REDC reduction if the modulus is odd and the caller is happy
   // to receive the result with a fudge factor
   int odd = (mod->n & 1);
   int redc = fastred && odd;

   if (len2 == 1)
   {
      // special case for 1xN multiplication
      _zn_array_scalar_mul(res, op1, len1, op2[0], redc, mod);
      return;
   }

   tuning_info_t* i = &tuning_info[mod->bits];
   
   if (op1 != op2  ||  len1 != len2)
   {
      // multiplying two distinct inputs
      
      if (len2 < i->mul_KS2_crossover)
         zn_array_mul_KS1(res, op1, len1, op2, len2, redc, mod);

      else if (len2 < i->mul_KS4_crossover)
         zn_array_mul_KS2(res, op1, len1, op2, len2, redc, mod);
         
      else if (!odd  ||  len2 < i->mul_fft_crossover)
         zn_array_mul_KS4(res, op1, len1, op2, len2, redc, mod);
         
      else
      {
         ulong scale = fastred ? 1 :
                       zn_array_mul_fft_get_fudge(len1, len2, 0, mod);
         zn_array_mul_fft(res, op1, len1, op2, len2, scale, mod);
      }
   }
   else
   {
      // squaring a single input

      if (len2 < i->sqr_KS2_crossover)
         zn_array_mul_KS1(res, op1, len1, op2, len2, redc, mod);

      else if (len2 < i->sqr_KS4_crossover)
         zn_array_mul_KS2(res, op1, len1, op2, len2, redc, mod);
         
      else if (!odd  ||  len2 < i->sqr_fft_crossover)
         zn_array_mul_KS4(res, op1, len1, op2, len2, redc, mod);
         
      else
      {
         ulong scale = fastred ? 1 :
                       zn_array_mul_fft_get_fudge(len1, len2, 1, mod);
         zn_array_mul_fft(res, op1, len1, op2, len2, 1, mod);
      }
   }
}


void zn_array_mul(ulong* res, const ulong* op1, size_t len1,
                  const ulong* op2, size_t len2, const zn_mod_t mod)
{
   _zn_array_mul(res, op1, len1, op2, len2, 0, mod);
}


// end of file ****************************************************************
