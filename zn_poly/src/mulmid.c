/*
   mulmid.c:  middle products
   
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
zn_array_mulmid_fallback_fudge (size_t n1, size_t n2, const zn_mod_t mod)
{
   return _zn_array_mul_fudge (n1, n2, 0, mod);
}


void
zn_array_mulmid_fallback (ulong* res, 
                          const ulong* op1, size_t n1,
                          const ulong* op2, size_t n2,
                          int fastred, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   ZNP_FASTALLOC (temp, ulong, 6624, n1 + n2 - 1);

   // just do full product and extract relevant segment
   _zn_array_mul (temp, op1, n1, op2, n2, fastred, mod);
   zn_array_copy (res, temp + n2 - 1, n1 - n2 + 1);

   ZNP_FASTFREE (temp);
}


ulong
_zn_array_mulmid_fudge (size_t n1, size_t n2, const zn_mod_t mod)
{
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);

   if (!(mod->m & 1))
      // no fudge if the modulus is even.
      return 1;

   tuning_info_t* i = &tuning_info[mod->bits];

   if (n2 < i->mulmid_KS2_thresh  ||  n2 < i->mulmid_KS4_thresh ||
       n2 < i->mulmid_fft_thresh)
      // fudge is -B
      return mod->m - mod->B;

   // return whatever fudge is used by the fft middle product code
   return zn_array_mulmid_fft_fudge (n1, n2, mod);
}


void
_zn_array_mulmid (ulong* res,
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

   tuning_info_t* i = &tuning_info[mod->bits];

   if (n2 < i->mulmid_KS2_thresh)
      zn_array_mulmid_KS1 (res, op1, n1, op2, n2, redc, mod);

   else if (n2 < i->mulmid_KS4_thresh)
      zn_array_mulmid_KS2 (res, op1, n1, op2, n2, redc, mod);
      
   else if (!odd  ||  n2 < i->mulmid_fft_thresh)
      zn_array_mulmid_KS4 (res, op1, n1, op2, n2, redc, mod);
      
   else
   {
      ulong x = fastred ? 1 : zn_array_mulmid_fft_fudge (n1, n2, mod);
      zn_array_mulmid_fft (res, op1, n1, op2, n2, x, mod);
   }
}


void
zn_array_mulmid (ulong* res,
                 const ulong* op1, size_t n1,
                 const ulong* op2, size_t n2,
                 const zn_mod_t mod)
{
   _zn_array_mulmid (res, op1, n1, op2, n2, 0, mod);
}



void
zn_array_mulmid_precomp1_init (zn_array_mulmid_precomp1_t res,
                               const ulong* op1, size_t n1, size_t n2,
                               const zn_mod_t mod)
{
   res->n1 = n1;
   res->n2 = n2;
   res->mod = mod;

   // figure out which algorithm to use

   int odd = (mod->m & 1);

   if (!odd)
      // can't use FFT algorithm when modulus is even
      res->algo = ZNP_MULMID_ALGO_KS;
   else
   {
      tuning_info_t* i = &tuning_info[mod->bits];

      res->algo = (n2 < i->mulmid_fft_thresh) ? ZNP_MULMID_ALGO_KS
                                              : ZNP_MULMID_ALGO_FFT;
   }

   // now perform initialisation for chosen algorithm

   switch (res->algo)
   {
      case ZNP_MULMID_ALGO_KS:
      {
         // Make a copy of op1[0, n1).

         // If modulus is odd, multiply it by the appropriate fudge factor
         // so that we can use faster REDC reduction in the execute() routine.
         res->op1 = (ulong*) malloc (n1 * sizeof (ulong));
         if (odd)
            zn_array_scalar_mul (res->op1, op1, n1, mod->m - mod->B, mod);
         else
            zn_array_copy (res->op1, op1, n1);
      }
      break;

      case ZNP_MULMID_ALGO_FFT:
      {
         res->precomp_fft = (struct zn_array_mulmid_fft_precomp1_struct*)
                              malloc (sizeof (zn_array_mulmid_fft_precomp1_t));

         // we do scaling in this init() routine, to avoid doing it during
         // each call to execute()
         ulong x = zn_array_mulmid_fft_precomp1_fudge (n1, n2, mod);
         zn_array_mulmid_fft_precomp1_init (res->precomp_fft,
                                            op1, n1, n2, x, mod);
      }
      break;

      default: ZNP_ASSERT (0);
   }
}


void
zn_array_mulmid_precomp1_clear (zn_array_mulmid_precomp1_t op)
{
   // dispatch to appropriate cleanup code
   switch (op->algo)
   {
      case ZNP_MULMID_ALGO_KS:
         free (op->op1);
         break;

      case ZNP_MULMID_ALGO_FFT:
         zn_array_mulmid_fft_precomp1_clear (op->precomp_fft);
         free (op->precomp_fft);
         break;

      default: ZNP_ASSERT (0);
   }
}



void
zn_array_mulmid_precomp1_execute (ulong* res, const ulong* op2,
                                  const zn_array_mulmid_precomp1_t precomp)
{
   // dispatch to appropriate middle product code
   switch (precomp->algo)
   {
      case ZNP_MULMID_ALGO_KS:
         _zn_array_mulmid (res, precomp->op1, precomp->n1,
                           op2, precomp->n2, precomp->mod->m & 1,
                           precomp->mod);
         break;

      case ZNP_MULMID_ALGO_FFT:
         zn_array_mulmid_fft_precomp1_execute (res, op2, 1,
                                               precomp->precomp_fft);
         break;
         
      default: ZNP_ASSERT (0);
   }
}


// end of file ****************************************************************
