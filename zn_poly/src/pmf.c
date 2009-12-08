/*
   pmf.c:  polynomials modulo a fermat polynomial
   
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


#if DEBUG

/* ============================================================================

     debugging stuff

============================================================================ */

#include <stdio.h>
#include "support.h"


void
pmf_rand (pmf_t res, ulong M, const zn_mod_t mod)
{
   ulong i;

   res[0] = random_ulong (2 * M);
   for (i = 1; i <= M; i++)
      res[i] = random_ulong (mod->m);
}


/*
   res := op, with bias set to zero.

   Inplace operation not allowed.
*/
void
pmf_normalise (pmf_t res, pmf_t op, ulong M, const zn_mod_t mod)
{
   ZNP_ASSERT (res != op);
   
   res[0] = 0;

   ulong i, b = op[0] & (2*M - 1);

   op++;
   res++;

   if (b < M)
   {
      for (i = 0; i < M - b; i++)
         res[b + i] = op[i];
      for (i = 0; i < b; i++)
         res[i] = zn_mod_neg (op[i + M - b], mod);
   }
   else
   {
      b -= M;
      for (i = 0; i < M - b; i++)
         res[b + i] = zn_mod_neg (op[i], mod);
      for (i = 0; i < b; i++)
         res[i] = op[i + M - b];
   }
}


int
pmf_cmp (const pmf_t op1, const pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong i;
   ulong b = op2[0] - op1[0];
   ulong* x1 = op1 + 1;
   ulong* x2 = op2 + 1;
   
   if (b & M)
   {
      b &= (M - 1);

      for (i = b; i < M; i++)
         if (x1[i] != zn_mod_neg (x2[i - b], mod))
            return 1;

      for (i = 0; i < b; i++)
         if (x1[i] != x2[i + M - b])
            return 1;
   }
   else
   {
      b &= (M - 1);

      for (i = b; i < M; i++)
         if (x1[i] != x2[i - b])
            return 1;

      for (i = 0; i < b; i++)
         if (x1[i] != zn_mod_neg (x2[i + M - b], mod))
            return 1;
   }
   
   return 0;
}


void
pmf_print (const pmf_t op, ulong M, const zn_mod_t mod)
{
   ZNP_FASTALLOC (buf, ulong, 6624, M + 1);

   pmf_normalise (buf, op, M, mod);
   
   printf ("[%lu", buf[1]);
   ulong i;
   for (i = 1; i < M; i++)
      printf (" %lu", buf[i+1]);
   printf ("]");
   
   ZNP_FASTFREE (buf);
}


void
pmfvec_print (const pmfvec_t op)
{
   printf ("M = %lu, K = %lu\n", op->M, op->K);
   ulong i;
   for (i = 0; i < op->K; i++)
   {
      printf ("%3lu: ", i);
      pmf_print (op->data + i * op->skip, op->M, op->mod);
      printf ("\n");
   }
}


void
pmfvec_print_trunc (const pmfvec_t op, ulong n)
{
   printf ("M = %lu, K = %lu\n", op->M, op->K);
   ulong i;
   for (i = 0; i < n; i++)
   {
      printf ("%3lu: ", i);
      pmf_print (op->data + i * op->skip, op->M, op->mod);
      printf ("\n");
   }
}

#endif



/* ============================================================================

     inplace array butterflies

============================================================================ */


/*
   op1 := op2 + op1
   op2 := op2 - op1
   where op1 and op2 are arrays of length n.

   Inputs must be [0, m); outputs will be in [0, m).
*/
void
zn_array_bfly_inplace (ulong* op1, ulong* op2, ulong n, const zn_mod_t mod)
{
   ulong x, y;
   
   if (zn_mod_is_slim (mod))
   {
      // slim version
      // (unrolled)
      for (; n >= 4; n -= 4)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim (y, x, mod);
         *op2++ = zn_mod_sub_slim (y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim (y, x, mod);
         *op2++ = zn_mod_sub_slim (y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim (y, x, mod);
         *op2++ = zn_mod_sub_slim (y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim (y, x, mod);
         *op2++ = zn_mod_sub_slim (y, x, mod);
      }
      for (; n; n--)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim (y, x, mod);
         *op2++ = zn_mod_sub_slim (y, x, mod);
      }
   }
   else
   {
      // non-slim version
      // (unrolled)
      for (; n >= 4; n -= 4)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add (y, x, mod);
         *op2++ = zn_mod_sub (y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add (y, x, mod);
         *op2++ = zn_mod_sub (y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add (y, x, mod);
         *op2++ = zn_mod_sub (y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add (y, x, mod);
         *op2++ = zn_mod_sub (y, x, mod);
      }
      for (; n; n--)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add (y, x, mod);
         *op2++ = zn_mod_sub (y, x, mod);
      }
   }
}


/*
   op1 := op1 + op2
   where op1 and op2 are arrays of length n.

   Inputs must be [0, m); outputs will be in [0, m).
*/
void
zn_array_add_inplace (ulong* op1, const ulong* op2, ulong n,
                      const zn_mod_t mod)
{
   if (zn_mod_is_slim (mod))
   {
      // slim version
      // (unrolled)
      for (; n >= 4; n -= 4)
      {
         *op1 = zn_mod_add_slim (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add_slim (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add_slim (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add_slim (*op1, *op2, mod);
         op1++; op2++;
      }
      for (; n; n--)
      {
         *op1 = zn_mod_add_slim (*op1, *op2, mod);
         op1++; op2++;
      }
   }
   else
   {
      // non-slim version
      // (unrolled)
      for (; n >= 4; n -= 4)
      {
         *op1 = zn_mod_add (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add (*op1, *op2, mod);
         op1++; op2++;
      }
      for (; n; n--)
      {
         *op1 = zn_mod_add (*op1, *op2, mod);
         op1++; op2++;
      }
   }
}



/*
   op1 := op1 - op2
   where op1 and op2 are arrays of length n.

   Inputs must be [0, m); outputs will be in [0, m).
*/
void
zn_array_sub_inplace (ulong* op1, const ulong* op2, ulong n,
                      const zn_mod_t mod)
{
   if (zn_mod_is_slim (mod))
   {
      // slim version
      // (unrolled)
      for (; n >= 4; n -= 4)
      {
         *op1 = zn_mod_sub_slim (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub_slim (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub_slim (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub_slim (*op1, *op2, mod);
         op1++; op2++;
      }
      for (; n; n--)
      {
         *op1 = zn_mod_sub_slim (*op1, *op2, mod);
         op1++; op2++;
      }
   }
   else
   {
      // non-slim version
      // (unrolled)
      for (; n >= 4; n -= 4)
      {
         *op1 = zn_mod_sub (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub (*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub (*op1, *op2, mod);
         op1++; op2++;
      }
      for (; n; n--)
      {
         *op1 = zn_mod_sub (*op1, *op2, mod);
         op1++; op2++;
      }
   }
}



/* ============================================================================

     inplace pmf_t butterflies

============================================================================ */


/*
   In the following routines, we work with the "relative bias" between the
   two inputs, i.e. b = difference between bias of op1 and op2. This allows
   us to avoid unnecessary normalisation steps; we just add and subtract
   directly into the correct memory locations.
*/


void
pmf_bfly (pmf_t op1, pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong b = op2[0] - op1[0];
   
   if (b & M)
   {
      // bias is in [M, 2M) mod 2M
      b &= (M - 1);
      
      // butterfly on op1[0, b) and op2[M - b, M)
      zn_array_bfly_inplace (op1 + 1, op2 + 1 + M - b, b, mod);
      // butterfly on op1[b, M) and op2[0, M - b)
      zn_array_bfly_inplace (op2 + 1, op1 + 1 + b, M - b, mod);
   }
   else
   {
      // bias is in [0, M) mod 2M
      b &= (M - 1);
      
      // butterfly on op1[0, b) and op2[M - b, M)
      zn_array_bfly_inplace (op2 + 1 + M - b, op1 + 1, b, mod);
      // butterfly on op1[b, M) and op2[0, M - b)
      zn_array_bfly_inplace (op1 + 1 + b, op2 + 1, M - b, mod);
   }
}


void
pmf_add (pmf_t op1, const pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong b = op2[0] - op1[0];

   if (b & M)
   {
      // bias is in [M, 2M) mod 2M
      b &= (M - 1);
      
      // add op2[M - b, M) to op1[0, b)
      zn_array_add_inplace (op1 + 1, op2 + 1 + M - b, b, mod);
      // subtract op2[0, M - b) from op1[b, M)
      zn_array_sub_inplace (op1 + 1 + b, op2 + 1, M - b, mod);
   }
   else
   {
      // bias is in [0, M) mod 2M
      b &= (M - 1);

      // subtract op2[M - b, M) from op1[0, b)
      zn_array_sub_inplace (op1 + 1, op2 + 1 + M - b, b, mod);
      // add op2[0, M - b) to op1[b, M)
      zn_array_add_inplace (op1 + 1 + b, op2 + 1, M - b, mod);
   }
}


void
pmf_sub (pmf_t op1, const pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong b = op2[0] - op1[0];

   if (b & M)
   {
      // bias is in [M, 2M) mod 2M
      b &= (M - 1);
      
      // subtract op2[M - b, M) from op1[0, b)
      zn_array_sub_inplace (op1 + 1, op2 + 1 + M - b, b, mod);
      // add op2[0, M - b) to op1[b, M)
      zn_array_add_inplace (op1 + 1 + b, op2 + 1, M - b, mod);
   }
   else
   {
      // bias is in [0, M) mod 2M
      b &= (M - 1);
      
      // add op2[M - b, M) to op1[0, b)
      zn_array_add_inplace (op1 + 1, op2 + 1 + M - b, b, mod);
      // subtract op2[0, M - b) from op1[b, M)
      zn_array_sub_inplace (op1 + 1 + b, op2 + 1, M - b, mod);
   }
}



/* ============================================================================

     pmfvec stuff

============================================================================ */


void
pmfvec_init (pmfvec_t res, unsigned lgK, ptrdiff_t skip, unsigned lgM,
             const zn_mod_t mod)
{
   ZNP_ASSERT (skip >= (1UL << lgM) + 1);
   
   res->lgK = lgK;
   res->lgM = lgM;
   res->skip = skip;
   res->K = 1UL << lgK;
   res->M = 1UL << lgM;
   res->mod = mod;
   res->data = (ulong*) malloc (sizeof (ulong) * skip * res->K);
}


void
pmfvec_init_nuss (pmfvec_t res, unsigned lgL, const zn_mod_t mod)
{
   unsigned lgK, lgM;
   nuss_params (&lgK, &lgM, lgL);
   pmfvec_init (res, lgK, (1UL << lgM) + 1, lgM, mod);
}


void
pmfvec_clear (pmfvec_t op)
{
   free (op->data);
}


void
pmfvec_set (pmfvec_t res, const pmfvec_t op)
{
   ZNP_ASSERT (pmfvec_compatible (res, op));
   
   ulong i;
   for (i = 0; i < op->K; i++)
      pmf_set (res->data + i * res->skip, op->data + i * op->skip, op->M);
}


void
pmfvec_scalar_mul (pmfvec_t op, ulong n, ulong x)
{
   ZNP_ASSERT (n <= op->K);

   ulong i;
   pmf_t ptr = op->data;
   for (i = 0; i < n; i++, ptr += op->skip)
      pmf_scalar_mul (ptr, op->M, x, op->mod);
}


ulong
pmfvec_mul_fudge (unsigned lgM, int sqr, const zn_mod_t mod)
{
   int use_nuss = (lgM >= (sqr ? tuning_info[mod->bits].nuss_sqr_thresh
                               : tuning_info[mod->bits].nuss_mul_thresh));
                     
   if (use_nuss)
      return nuss_mul_fudge (lgM, sqr, mod);
   else
      return _zn_array_mul_fudge (1UL << lgM, 1UL << lgM, sqr, mod);
}


void
pmfvec_mul (pmfvec_t res, const pmfvec_t op1, const pmfvec_t op2, ulong n,
            int special_first_two)
{
   ZNP_ASSERT (res->mod->m & 1);
   ZNP_ASSERT (pmfvec_compatible (res, op1));
   ZNP_ASSERT (pmfvec_compatible (res, op2));
   ZNP_ASSERT (res->M >= 2 || !special_first_two);

   pmf_const_t p1 = op1->data;
   pmf_const_t p2 = op2->data;
   pmf_t p3 = res->data;

   const zn_mod_struct* mod = res->mod;
   
   ulong M = op1->M;
   unsigned lgM = op1->lgM;
   
   int sqr = (op1 == op2);
   
   // use nussbaumer algorithm if the pointwise mults are large enough
   int use_nuss = (lgM >= (sqr ? tuning_info[mod->bits].nuss_sqr_thresh
                               : tuning_info[mod->bits].nuss_mul_thresh));

   // scratch space for nussbaumer multiplications
   pmfvec_t vec1, vec2;
   unsigned nuss_lgK;

   if (use_nuss)
   {
      pmfvec_init_nuss (vec1, lgM, mod);
      pmfvec_init_nuss (vec2, lgM, mod);
      nuss_lgK = vec1->lgK;
   }

   ulong i = 0;

   if (special_first_two)
   {
      ZNP_FASTALLOC (temp, ulong, 6624, 2 * M);
      
      // need to adjust the fudge factors, so that the fudge factor for these
      // first two special products matches up with the fudge factors for the
      // remaining products
      ulong fudge, fudge1, fudge2;
      fudge2 = use_nuss ? nuss_mul_fudge (lgM, sqr, mod)
                        : _zn_array_mul_fudge (M, M, sqr, mod);
      fudge1 = _zn_array_mul_fudge (M/2, M/2, sqr, mod);
      fudge = (fudge1 == fudge2)
                   ? 1 : zn_mod_mul (fudge1, zn_mod_invert (fudge2, mod), mod);

      // length M/2 multiplications
      for (; i < 2 && i < n;
             i++, p3 += res->skip, p1 += op1->skip, p2 += op2->skip)
      {
         // add biases
         p3[0] = p1[0] + p2[0];

         // do the actual multiplication
         _zn_array_mul (temp, p1 + 1, M/2, p2 + 1, M/2, 1, mod);

         // apply the fudge factor
         zn_array_scalar_mul_or_copy (p3 + 1, temp, M - 1, fudge, mod);
         p3[M] = 0;
      }
      
      ZNP_FASTFREE(temp);
   }
   
   if (use_nuss)
   {
      for (; i < n; i++, p3 += res->skip,
             p1 += op1->skip, p2 += op2->skip)
      {
         // add biases
         p3[0] = p1[0] + p2[0];

         nuss_mul (p3 + 1, p1 + 1, p2 + 1, vec1, vec2);
      }

      pmfvec_clear (vec2);
      pmfvec_clear (vec1);
   }
   else
   {
      // scratch space for KS negacyclic multiplication
      ZNP_FASTALLOC (temp, ulong, 6624, 2*M);
      temp[2*M - 1] = 0;

      for (; i < n; i++, p3 += res->skip, p1 += op1->skip, p2 += op2->skip)
      {
         // add biases
         p3[0] = p1[0] + p2[0];
         
         // ordinary multiplication...
         _zn_array_mul (temp, p1 + 1, M, p2 + 1, M, 1, mod);
         // ... negacyclic reduction
         zn_array_sub (p3 + 1, temp, temp + M, M, mod);
      }

      ZNP_FASTFREE (temp);
   }
}



void
pmfvec_reverse (pmfvec_t op, ulong n)
{
   op->data += op->skip * (n - 1);
   op->skip = -op->skip;
}



// end of file ****************************************************************
