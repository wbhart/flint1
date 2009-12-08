/*
   mpn_mulmid.c:  middle products of integers
   
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
#include <string.h>


void
ZNP_mpn_smp_basecase (mp_limb_t* res,
                      const mp_limb_t* op1, size_t n1,
                      const mp_limb_t* op2, size_t n2)
{
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n2 >= 1);

#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   mp_limb_t hi0, hi1, hi;
   size_t s, j;
   
   j = n2 - 1;
   s = n1 - j;
   op2 += j;

   hi0 = mpn_mul_1 (res, op1, s, *op2);
   hi1 = 0;
   
   for (op1++, op2--; j; j--, op1++, op2--)
   {
      hi = mpn_addmul_1 (res, op1, s, *op2);
      ZNP_ADD_WIDE (hi1, hi0, hi1, hi0, 0, hi);
   }
   
   res[s] = hi0;
   res[s + 1] = hi1;

#else
#error Not nails-safe yet
#endif
}


/*
   Let x = op1[0, 2*n-1),
       y = op2[0, n),
       z = op3[0, n).
   
   If y >= z, this function computes y - z and the correction term
      SMP(x, y) - SMP(x, z) - SMP(x, y - z)
   and returns 0.
   
   If y < z, it computes z - y and the correction term
      SMP(x, z) - SMP(x, y) - SMP(x, z - y)
   and returns 1.
   
   In both cases abs(y - z) is stored at res[0, n).
   
   The correction term is v - u*B^n, where u is stored at hi[0, 2) and
   v is stored at lo[0, 2).
   
   None of the output buffers are allowed to overlap either each other or
   the input buffers.
*/
#define bilinear2_sub_fixup \
    ZNP_bilinear2_sub_fixup
int
bilinear2_sub_fixup (mp_limb_t* hi, mp_limb_t* lo, mp_limb_t* res,
                     const mp_limb_t* op1, const mp_limb_t* op2,
                     const mp_limb_t* op3, size_t n)
{
   ZNP_ASSERT (n >= 1);

#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS
   
   int sign = 0;
   if (mpn_cmp (op2, op3, n) < 0)
   {
      // swap y and z if necessary
      const mp_limb_t* temp = op2;
      op2 = op3;
      op3 = temp;
      sign = 1;
   }
   // now can assume y >= z
   
   // The correction term is computed as follows. Let
   //
   //    y_0     - z_0                =  u_0     - c_0 B,
   //    y_1     - z_1     - c_0      =  u_1     - c_1 B,
   //    y_2     - z_2     - c_1      =  u_2     - c_2 B,
   //                                ...
   //    y_{n-1} - z_{n-1} - c_{n-2}  =  u_{n-1},
   //
   // i.e. where c_j is the borrow (0 or 1) from the j-th limb of the
   // subtraction y - z, and where u_j is the j-th digit of y - z. Note
   // that c_{-1} = c_{n-1} = 0. By definition we want to compute
   //
   //    \sum_{0 <= i < 2n-1, 0 <= j < n, n-1 <= i+j < 2n-1}
   //                                  (c_{j-1} - c_j B) x_i B^{i+j-(n-1)}
   //
   // After some algebra this collapses down to
   //
   //    \sum_{0 <= i < n-1} c_i (x_{n-2-i} - B^n x_{2n-2-i}).

   // First compute y - z using mpn_sub_n (fast)
   mpn_sub_n (res, op2, op3, n);

   // Now loop through and figure out where the borrows happened
   size_t i;
   mp_limb_t hi0 = 0, hi1 = 0;
   mp_limb_t lo0 = 0, lo1 = 0;

   for (i = n - 1; i; i--, op1++)
   {
      mp_limb_t borrow = res[i] - op2[i] + op3[i];
      ZNP_ADD_WIDE (lo1, lo0, lo1, lo0, 0, borrow & op1[0]);
      ZNP_ADD_WIDE (hi1, hi0, hi1, hi0, 0, borrow & op1[n]);
   }
   
   hi[0] = hi0;
   hi[1] = hi1;
   lo[0] = lo0;
   lo[1] = lo1;

   return sign;
#else
#error Not nails-safe yet
#endif
}


/*
   Let x = op1[0, 2*n-1),
       y = op2[0, 2*n-1),
       z = op3[0, n).
   
   This function computes x + y mod B^(2n-1) and the correction term
      SMP(x, z) + SMP(y, z) - SMP((x + y) mod B^(2n-1), z).

   The value x + y mod B^(2n-1) is stored at res[0, 2n-1).
   
   The correction term is u*B^n - v, where u is stored at hi[0, 2) and
   v is stored at lo[0, 2).

   None of the output buffers are allowed to overlap either each other or
   the input buffers.
*/
#define bilinear1_add_fixup \
    ZNP_bilinear1_add_fixup
void
bilinear1_add_fixup (mp_limb_t* hi, mp_limb_t* lo, mp_limb_t* res,
                     const mp_limb_t* op1, const mp_limb_t* op2,
                     const mp_limb_t* op3, size_t n)
{
   ZNP_ASSERT (n >= 1);
   
#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   // The correction term is computed as follows. Let
   //
   //    x_0      + y_0                  =  u_0      + c_0 B,
   //    x_1      + y_1      + c_0       =  u_1      + c_1 B,
   //    x_2      + y_2      + c_1       =  u_2      + c_2 B,
   //                                   ...
   //    x_{2n-2} + y_{2n-2} + c_{2n-3}  =  u_{2n-2} + c_{2n-1} B,
   //
   // i.e. where c_j is the carry (0 or 1) from the j-th limb of the
   // addition x + y, and u_j is the j-th digit of x + y. Note that
   // c_{-1} = 0. By definition we want to compute
   //
   //    \sum_{0 <= i < 2n-1, 0 <= j < n, n-1 <= i+j < 2n-1}
   //                                  (c_i B - c_{i-1}) z_j B^{i+j-(n-1)}
   //
   // After some algebra this collapses down to
   //
   //     -\sum_{0 <= j < n-1}    c_j z_{n-2-j}  +
   //  B^n \sum_{n-1 <= j < 2n-1} c_j z_{2n-2-j}.

   // First compute x + y using mpn_add_n (fast)
   mp_limb_t last_carry = mpn_add_n (res, op1, op2, 2*n - 1);

   // Now loop through and figure out where the carries happened
   size_t j;
   mp_limb_t fix0 = 0, fix1 = 0;
   op3 += n - 2;
   
   for (j = 0; j < n - 1; j++, op3--)
   {
      // carry = -1 if there was a carry in the j-th limb addition
      mp_limb_t carry = op1[j+1] + op2[j+1] - res[j+1];
      ZNP_ADD_WIDE (fix1, fix0, fix1, fix0, 0, carry & *op3);
   }
   
   lo[0] = fix0;
   lo[1] = fix1;

   fix0 = fix1 = 0;
   op3 += n;
   
   for (; j < 2*n - 2; j++, op3--)
   {
      // carry = -1 if there was a carry in the j-th limb addition
      mp_limb_t carry = op1[j+1] + op2[j+1] - res[j+1];
      ZNP_ADD_WIDE (fix1, fix0, fix1, fix0, 0, carry & *op3);
   }
   
   ZNP_ADD_WIDE (fix1, fix0, fix1, fix0, 0, (-last_carry) & *op3);

   hi[0] = fix0;
   hi[1] = fix1;
#else
#error Not nails-safe yet
#endif
}



void
ZNP_mpn_smp_kara (mp_limb_t* res, const mp_limb_t* op1, const mp_limb_t* op2,
                  size_t n)
{
   ZNP_ASSERT (n >= 2);
   
   if (n & 1)
   {
      // If n is odd, we strip off the bottom row and last diagonal and
      // handle them separately at the end (stuff marked O in the diagram
      // below); the remainder gets handled via karatsuba (stuff marked E):
      
      // EEEEO....
      // .EEEEO...
      // ..EEEEO..
      // ...EEEEO.
      // ....OOOOO
      
      op2++;
   }

   size_t k = n / 2;
   
   ZNP_FASTALLOC (temp, mp_limb_t, 6642, 2 * k + 2);

   mp_limb_t hi[2], lo[2];

   // The following diagram shows the contributions from various regions
   // for k = 3:

   //  AAABBB.....
   //  .AAABBB....
   //  ..AAABBB...
   //  ...CCCDDD..
   //  ....CCCDDD.
   //  .....CCCDDD
   
   // ------------------------------------------------------------------------
   // Step 1: compute contribution from A + contribution from B

   // Let x = op1[0, 2*k-1)
   //     y = op1[k, 3*k-1)
   //     z = op2[k, 2*k).

   // Need to compute SMP(x, z) + SMP(y, z). To do this, we will compute
   // SMP((x + y) mod B^(2k-1), z) and a correction term.
   
   // First compute x + y mod B^(2k-1) and the correction term.
   bilinear1_add_fixup (hi, lo, temp, op1, op1 + k, op2 + k, k);

   // Now compute SMP(x + y mod B^(2k-1), z).
   // Store result in first half of output.
   if (k < ZNP_mpn_smp_kara_thresh)
      ZNP_mpn_smp_basecase (res, temp, 2 * k - 1, op2 + k, k);
   else
      ZNP_mpn_smp_kara (res, temp, op2 + k, k);
   
   // Add in the correction term.
   mpn_sub (res, res, k + 2, lo, 2);
   mpn_add_n (res + k, res + k, hi, 2);
   
   // Save the last two limbs (they're about to get overwritten)
   mp_limb_t saved[2];
   saved[0] = res[k];
   saved[1] = res[k + 1];

   // ------------------------------------------------------------------------
   // Step 2: compute contribution from C + contribution from D

   // Let x = op1[k, 3*k-1)
   //     y = op1[2*k, 4*k-1)
   //     z = op2[0, k).

   // Need to compute SMP(x, z) + SMP(y, z). To do this, we will compute
   // SMP((x + y) mod B^(2k-1), z) and a correction term.
   
   // First compute x + y mod B^(2k-1) and the correction term.
   bilinear1_add_fixup (hi, lo, temp, op1 + k, op1 + 2 * k, op2, k);

   // Now compute SMP(x + y mod B^(2k-1), z).
   // Store result in second half of output.
   if (k < ZNP_mpn_smp_kara_thresh)
      ZNP_mpn_smp_basecase (res + k, temp, 2 * k - 1, op2, k);
   else
      ZNP_mpn_smp_kara (res + k, temp, op2, k);
   
   // Add in the correction term.
   mpn_sub (res + k, res + k, k + 2, lo, 2);
   mpn_add_n (res + 2 * k, res + 2 * k, hi, 2);
   
   // Add back the saved limbs.
   mpn_add (res + k, res + k, k + 2, saved, 2);

   // ------------------------------------------------------------------------
   // Step 3: compute contribution from B - contribution from C

   // Let x = op1[k, 3*k-1)
   //     y = op2[k, 2*k).
   //     z = op2[0, k)
   
   // Need to compute SMP(x, y) - SMP(x, z). To do this, we will compute
   // SMP(x, abs(y - z)), and a correction term.
   
   // First compute abs(y - z) and the correction term.
   int sign = bilinear2_sub_fixup (hi, lo, temp, op1 + k, op2 + k, op2, k);

   // Now compute SMP(x, abs(y - z)).
   // Store it in second half of temp space, in two's complement (mod B^(k+2))
   if (k < ZNP_mpn_smp_kara_thresh)
      ZNP_mpn_smp_basecase (temp + k, op1 + k, 2 * k - 1, temp, k);
   else
      ZNP_mpn_smp_kara (temp + k, op1 + k, temp, k);
   
   // Add in the correction term.
   mpn_add (temp + k, temp + k, k + 2, lo, 2);
   mp_limb_t borrow = mpn_sub_n (temp + 2 * k, temp + 2 * k, hi, 2);

   // ------------------------------------------------------------------------
   // Step 4: put the pieces together

   // First half of output is A + C = t4 - t2
   // Second half of output is B + D = t6 + t2
   if (sign)
   {
      mpn_add (res, res, 2 * k + 2, temp + k, k + 2);
      mpn_sub_1 (res + k + 2, res + k + 2, k, borrow);
      mpn_sub (res + k, res + k, k + 2, temp + k, k + 2);
   }
   else
   {
      mpn_sub (res, res, 2 * k + 2, temp + k, k + 2);
      mpn_add_1 (res + k + 2, res + k + 2, k, borrow);
      mpn_add (res + k, res + k, k + 2, temp + k, k + 2);
   }

   // ------------------------------------------------------------------------
   // Step 5: add in correction if the length was odd

#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS

   if (n & 1)
   {
      op2--;
      
      mp_limb_t hi0 = mpn_addmul_1 (res, op1 + n - 1, n, *op2);
      mp_limb_t hi1 = 0, lo0 = 0, lo1 = 0;

      size_t i;
      for (i = n - 1; i; i--)
      {
         mp_limb_t y0, y1;
         ZNP_MUL_WIDE (y1, y0, op1[2 * n - i - 2], op2[i]);
         ZNP_ADD_WIDE (hi1, hi0, hi1, hi0, 0, y1);
         ZNP_ADD_WIDE (lo1, lo0, lo1, lo0, 0, y0);
      }
      
      res[n + 1] = hi1;
      mpn_add_1 (res + n, res + n, 2, hi0);
      mpn_add_1 (res + n, res + n, 2, lo1);
      mpn_add_1 (res + n - 1, res + n - 1, 3, lo0);
   }
   
   ZNP_FASTFREE (temp);

#else
#error Not nails-safe yet
#endif
}


void
ZNP_mpn_smp_n (mp_limb_t* res, const mp_limb_t* op1, const mp_limb_t* op2,
               size_t n)
{
   if (n < ZNP_mpn_smp_kara_thresh)
      ZNP_mpn_smp_basecase (res, op1, 2*n - 1, op2, n);
   else
      ZNP_mpn_smp_kara (res, op1, op2, n);
}


void
ZNP_mpn_smp (mp_limb_t* res,
             const mp_limb_t* op1, size_t n1,
             const mp_limb_t* op2, size_t n2)
{
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n2 >= 1);

   size_t n3 = n1 - n2 + 1;

   if (n3 < ZNP_mpn_smp_kara_thresh)
   {
      // region is too narrow to make karatsuba worthwhile for any portion
      ZNP_mpn_smp_basecase (res, op1, n1, op2, n2);
      return;
   }
   
   if (n2 > n3)
   {
      // slice region into chunks horizontally, i.e. like this:

      // AA.....
      // .AA....
      // ..BB...
      // ...BB..
      // ....CC.
      // .....CC
      
      // first chunk (marked A in the above diagram)
      op2 += n2 - n3;
      ZNP_mpn_smp_kara (res, op1, op2, n3);
      
      // remaining chunks (B, C, etc)
      ZNP_FASTALLOC (temp, mp_limb_t, 6642, n3 + 2);

      n1 -= n3;
      n2 -= n3;
      
      while (n2 >= n3)
      {
         op1 += n3;
         op2 -= n3;
         ZNP_mpn_smp_kara (temp, op1, op2, n3);
         mpn_add_n (res, res, temp, n3 + 2);
         n1 -= n3;
         n2 -= n3;
      }
      
      if (n2)
      {
         // last remaining chunk
         op1 += n3;
         op2 -= n2;
         ZNP_mpn_smp (temp, op1, n1, op2, n2);
         mpn_add_n (res, res, temp, n3 + 2);
      }
      
      ZNP_FASTFREE (temp);
   }
   else
   {
      mp_limb_t save[2];
      
      // slice region into chunks diagonally, i.e. like this:
      
      // AAABBBCC..
      // .AAABBBCC.
      // ..AAABBBCC

      // first chunk (marked A in the above diagram)
      ZNP_mpn_smp_n (res, op1, op2, n2);
      
      n1 -= n2;
      n3 -= n2;

      // remaining chunks (B, C, etc)
      while (n3 >= n2)
      {
         op1 += n2;
         res += n2;
         
         // save two limbs which are going to be overwritten
         save[0] = res[0];
         save[1] = res[1];
         
         ZNP_mpn_smp_n (res, op1, op2, n2);
         
         // add back saved limbs
         mpn_add (res, res, n2 + 2, save, 2);
         
         n1 -= n2;
         n3 -= n2;
      }

      if (n3)
      {
         // last remaining chunk
         op1 += n2;
         res += n2;

         save[0] = res[0];
         save[1] = res[1];

         ZNP_mpn_smp (res, op1, n1, op2, n2);

         mpn_add (res, res, n3 + 2, save, 2);
      }
   }
}


void
ZNP_mpn_mulmid_fallback (mp_limb_t* res,
                         const mp_limb_t* op1, size_t n1,
                         const mp_limb_t* op2, size_t n2)
{
   if (n1 < n2 + 1)
      return;

   ZNP_FASTALLOC (temp, mp_limb_t, 6642, n1 + n2);
   ZNP_mpn_mul (temp, op1, n1, op2, n2);
   memcpy (res + 2, temp + n2 + 1, sizeof(mp_limb_t) * (n1 - n2 - 1));
   ZNP_FASTFREE (temp);
}


void
ZNP_mpn_mulmid (mp_limb_t* res, const mp_limb_t* op1, size_t n1,
                const mp_limb_t* op2, size_t n2)
{
   ZNP_ASSERT (n1 >= n2);
   ZNP_ASSERT (n2 >= 1);
   
   if (n2 >= ZNP_mpn_mulmid_fallback_thresh)
   {
      ZNP_mpn_mulmid_fallback (res, op1, n1, op2, n2);
      return;
   }

   // try using the simple middle product
   ZNP_mpn_smp (res, op1, n1, op2, n2);

#if GMP_NAIL_BITS == 0  &&  ULONG_BITS == GMP_NUMB_BITS
   
   // If there's a possibility of overflow from lower diagonals, we just give
   // up and do the whole product. (Note: this should happen extremely rarely
   // on uniform random input. However, on data generated by mpn_random2, it
   // seems to happen with non-negligible probability.)
   if (res[1] >= -(mp_limb_t)(n2))
      ZNP_mpn_mulmid_fallback (res, op1, n1, op2, n2);

#else
#error Not nails-safe yet
#endif
}


// end of file ****************************************************************
