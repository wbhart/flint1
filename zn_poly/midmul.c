/*
   midmul.c:  middle products
   
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


ulong zn_array_midmul_fallback_get_fudge(size_t len1, size_t len2,
                                         const zn_mod_t mod)
{
   return _zn_array_mul_get_fudge(len1, len2, 0, mod);
}


void zn_array_midmul_fallback(ulong* res, const ulong* op1, size_t len1,
                              const ulong* op2, size_t len2, int fastred,
                              const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   ZNP_FASTALLOC(temp, ulong, 6624, len1 + len2 - 1);

   // just do full product and extract relevant segment
   _zn_array_mul(temp, op1, len1, op2, len2, fastred, mod);
   zn_array_copy(res, temp + len2 - 1, len1 - len2 + 1);

   ZNP_FASTFREE(temp);
}


void _zn_array_midmul(ulong* res, const ulong* op1, size_t len1,
                      const ulong* op2, size_t len2, int fastred,
                      const zn_mod_t mod)
{
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);

   tuning_info_t* i = &tuning_info[mod->bits];

   if (len2 < i->midmul_fft_crossover  ||  !(mod->n & 1))
      // (can't use FFT algorithm if the modulus is even)
      zn_array_midmul_fallback(res, op1, len1, op2, len2, fastred, mod);
   else
   {
      ulong scale = zn_array_midmul_fft_get_fudge(len1, len2, mod);
      zn_array_midmul_fft(res, op1, len1, op2, len2, scale, mod);
   }
}


void zn_array_midmul(ulong* res, const ulong* op1, size_t len1,
                     const ulong* op2, size_t len2, const zn_mod_t mod)
{
   _zn_array_midmul(res, op1, len1, op2, len2, 0, mod);
}



void zn_array_midmul_precomp1_init(zn_array_midmul_precomp1_t res,
                                   const ulong* op1, size_t len1,
                                   size_t len2, const zn_mod_t mod)
{
   res->len1 = len1;
   res->len2 = len2;
   res->mod = mod;

   int odd = (mod->n & 1);

   // figure out which algorithm to use

   if (!odd)
      // can't use FFT algorithm when modulus is even
      res->algo = ZNP_MIDMUL_ALGO_FALLBACK;
   else
   {
      tuning_info_t* i = &tuning_info[mod->bits];
   
      if (len2 < i->midmul_fft_crossover)
         res->algo = ZNP_MIDMUL_ALGO_FALLBACK;
      else
         res->algo = ZNP_MIDMUL_ALGO_FFT;
   }

   // now perform initialisation for chosen algorithm

   switch (res->algo)
   {
      case ZNP_MIDMUL_ALGO_FALLBACK:
      {
         // Make a copy of op1[0, len1).

         // If modulus is odd, multiply it by the appropriate fudge factor
         // so that we can use faster REDC reduction in the execute() routine.
         res->op1 = (ulong*) malloc(len1 * sizeof(ulong));
         if (odd)
         {
            ulong scale = zn_array_midmul_fallback_get_fudge(len1, len2, mod);
            zn_array_scalar_mul(res->op1, op1, len1, scale, mod);
         }
         else
            zn_array_copy(res->op1, op1, len1);
      }
      break;

      case ZNP_MIDMUL_ALGO_FFT:
      {
         res->precomp_fft = (struct zn_array_midmul_fft_precomp1_struct*)
                               malloc(sizeof(zn_array_midmul_fft_precomp1_t));

         // we do scaling in this init() routine, to avoid doing it during
         // each call to execute()
         ulong scale = zn_array_midmul_fft_precomp1_get_fudge(len1, len2, mod);
         zn_array_midmul_fft_precomp1_init(res->precomp_fft,
                                           op1, len1, len2, scale, mod);
      }
      break;

      default: ZNP_ASSERT(0);
   }
}


void zn_array_midmul_precomp1_clear(zn_array_midmul_precomp1_t op)
{
   // dispatch to appropriate cleanup code
   switch (op->algo)
   {
      case ZNP_MIDMUL_ALGO_FALLBACK:
         free(op->op1);
         break;

      case ZNP_MIDMUL_ALGO_FFT:
         zn_array_midmul_fft_precomp1_clear(op->precomp_fft);
         free(op->precomp_fft);
         break;

      default: ZNP_ASSERT(0);
   }
}



void zn_array_midmul_precomp1_execute(
            ulong* res, const ulong* op2,
            const zn_array_midmul_precomp1_t precomp)
{
   // dispatch to appropriate middle product code
   switch (precomp->algo)
   {
      case ZNP_MIDMUL_ALGO_FALLBACK:
         zn_array_midmul_fallback(res, precomp->op1, precomp->len1,
                                  op2, precomp->len2, precomp->mod->n & 1,
                                  precomp->mod);
         break;

      case ZNP_MIDMUL_ALGO_FFT:
         zn_array_midmul_fft_precomp1_execute(res, op2, 1,
                                              precomp->precomp_fft);
         break;
         
      default: ZNP_ASSERT(0);
   }
}


// end of file ****************************************************************
