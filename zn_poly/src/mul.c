/*
   mul.c:  polynomial multiplication
   
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


ulong
_zn_array_mul_fudge (size_t n1, size_t n2, int sqr, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   if (!(mod->m & 1))
      // no fudge if the modulus is even.
      return 1;

   tuning_info_t* i = &tuning_info[mod->bits];

   if (!sqr)
   {
      if (n2 < i->mul_KS2_thresh  ||  n2 < i->mul_KS4_thresh ||
          n2 < i->mul_fft_thresh)
         // fudge is -B
         return mod->m - mod->B;
   }
   else
   {
      if (n2 < i->sqr_KS2_thresh  ||  n2 < i->sqr_KS4_thresh ||
          n2 < i->sqr_fft_thresh)
         // fudge is -B
         return mod->m - mod->B;
   }

   // return whatever fudge is used by the fft multiplication code
   return zn_array_mul_fft_fudge (n1, n2, sqr, mod);
}


void
_zn_array_mul (ulong* res,
               const ulong* op1, size_t n1,
               const ulong* op2, size_t n2,
               int fastred, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   // we can use REDC reduction if the modulus is odd and the caller is happy
   // to receive the result with a fudge factor
   int odd = (mod->m & 1);
   int redc = fastred && odd;

   if (n2 == 1)
   {
      // special case for 1xN multiplication
      _zn_array_scalar_mul (res, op1, n1, op2[0], redc, mod);
      return;
   }

   tuning_info_t* i = &tuning_info[mod->bits];
   
   if (op1 != op2  ||  n1 != n2)
   {
      // multiplying two distinct inputs
      
      if (n2 < i->mul_KS2_thresh)
         zn_array_mul_KS1 (res, op1, n1, op2, n2, redc, mod);

      else if (n2 < i->mul_KS4_thresh)
         zn_array_mul_KS2 (res, op1, n1, op2, n2, redc, mod);
         
      else if (!odd  ||  n2 < i->mul_fft_thresh)
         zn_array_mul_KS4 (res, op1, n1, op2, n2, redc, mod);
         
      else
      {
         ulong x = fastred ? 1 : zn_array_mul_fft_fudge (n1, n2, 0, mod);
         zn_array_mul_fft (res, op1, n1, op2, n2, x, mod);
      }
   }
   else
   {
      // squaring a single input

      if (n2 < i->sqr_KS2_thresh)
         zn_array_mul_KS1 (res, op1, n1, op1, n1, redc, mod);

      else if (n2 < i->sqr_KS4_thresh)
         zn_array_mul_KS2 (res, op1, n1, op1, n1, redc, mod);
         
      else if (!odd  ||  n2 < i->sqr_fft_thresh)
         zn_array_mul_KS4 (res, op1, n1, op1, n1, redc, mod);
         
      else
      {
         ulong x = fastred ? 1 : zn_array_mul_fft_fudge (n1, n1, 1, mod);
         zn_array_mul_fft (res, op1, n1, op1, n1, x, mod);
      }
   }
}


void
zn_array_mul (ulong* res,
              const ulong* op1, size_t n1,
              const ulong* op2, size_t n2,
              const zn_mod_t mod)
{
   _zn_array_mul (res, op1, n1, op2, n2, 0, mod);
}


// end of file ****************************************************************
