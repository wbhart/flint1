/*
   pmf.c:  polynomials modulo a fermat polynomial
   
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


#if DEBUG

/* ============================================================================

     debugging stuff

============================================================================ */

#include <stdio.h>


/*
   res := op, with bias set to zero.

   Inplace operation not allowed.
*/
void zn_pmf_normalise(zn_pmf_t res, zn_pmf_t op, ulong M, const zn_mod_t mod)
{
   ZNP_ASSERT(res != op);
   
   res[0] = 0;

   ulong i, b = op[0] & (2*M - 1);

   op++;
   res++;

   if (b < M)
   {
      for (i = 0; i < M-b; i++)
         res[b+i] = op[i];
      for (i = 0; i < b; i++)
         res[i] = zn_mod_neg(op[i+M-b], mod);
   }
   else
   {
      b -= M;
      for (i = 0; i < M-b; i++)
         res[b+i] = zn_mod_neg(op[i], mod);
      for (i = 0; i < b; i++)
         res[i] = op[i+M-b];
   }
}


void zn_pmf_print(const zn_pmf_t op, ulong M, const zn_mod_t mod)
{
   ZNP_FASTALLOC(buf, ulong, 6624, M + 1);

   zn_pmf_normalise(buf, op, M, mod);
   
   printf("[%lu", buf[1]);
   ulong i;
   for (i = 1; i < M; i++)
      printf(" %lu", buf[i+1]);
   printf("]");
   
   ZNP_FASTFREE(buf);
}


void zn_pmf_vec_print(const zn_pmf_vec_t op)
{
   printf("M = %lu, K = %lu\n", op->M, op->K);
   ulong i;
   for (i = 0; i < op->K; i++)
   {
      printf("%3lu: ", i);
      zn_pmf_print(op->data + i*op->skip, op->M, op->mod);
      printf("\n");
   }
}


void zn_pmf_vec_print_trunc(const zn_pmf_vec_t op, ulong length)
{
   printf("M = %lu, K = %lu\n", op->M, op->K);
   ulong i;
   for (i = 0; i < length; i++)
   {
      printf("%3lu: ", i);
      zn_pmf_print(op->data + i*op->skip, op->M, op->mod);
      printf("\n");
   }
}

#endif



/* ============================================================================

     inplace array butterflies

============================================================================ */


/*
   op1 := op2 + op1
   op2 := op2 - op1
   where op1 and op2 are arrays of length _len_.

   Inputs must be [0, n); outputs will be in [0, n).
*/
void zn_array_bfly_inplace(ulong* op1, ulong* op2, ulong len,
                           const zn_mod_t mod)
{
   ulong x, y;
   
   if (zn_mod_is_slim(mod))
   {
      // slim version
      // (unrolled)
      for (; len >= 4; len -= 4)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim(y, x, mod);
         *op2++ = zn_mod_sub_slim(y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim(y, x, mod);
         *op2++ = zn_mod_sub_slim(y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim(y, x, mod);
         *op2++ = zn_mod_sub_slim(y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim(y, x, mod);
         *op2++ = zn_mod_sub_slim(y, x, mod);
      }
      for (; len; len--)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add_slim(y, x, mod);
         *op2++ = zn_mod_sub_slim(y, x, mod);
      }
   }
   else
   {
      // non-slim version
      // (unrolled)
      for (; len >= 4; len -= 4)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add(y, x, mod);
         *op2++ = zn_mod_sub(y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add(y, x, mod);
         *op2++ = zn_mod_sub(y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add(y, x, mod);
         *op2++ = zn_mod_sub(y, x, mod);

         x = *op1; y = *op2;
         *op1++ = zn_mod_add(y, x, mod);
         *op2++ = zn_mod_sub(y, x, mod);
      }
      for (; len; len--)
      {
         x = *op1; y = *op2;
         *op1++ = zn_mod_add(y, x, mod);
         *op2++ = zn_mod_sub(y, x, mod);
      }
   }
}


/*
   op1 := op1 + op2
   where op1 and op2 are arrays of length _len_.

   Inputs must be [0, n); outputs will be in [0, n).
*/
void zn_array_add_inplace(ulong* op1, const ulong* op2, ulong len,
                          const zn_mod_t mod)
{
   if (zn_mod_is_slim(mod))
   {
      // slim version
      // (unrolled)
      for (; len >= 4; len -= 4)
      {
         *op1 = zn_mod_add_slim(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add_slim(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add_slim(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add_slim(*op1, *op2, mod);
         op1++; op2++;
      }
      for (; len; len--)
      {
         *op1 = zn_mod_add_slim(*op1, *op2, mod);
         op1++; op2++;
      }
   }
   else
   {
      // non-slim version
      // (unrolled)
      for (; len >= 4; len -= 4)
      {
         *op1 = zn_mod_add(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_add(*op1, *op2, mod);
         op1++; op2++;
      }
      for (; len; len--)
      {
         *op1 = zn_mod_add(*op1, *op2, mod);
         op1++; op2++;
      }
   }
}



/*
   op1 := op1 - op2
   where op1 and op2 are arrays of length _len_.

   Inputs must be [0, n); outputs will be in [0, n).
*/
void zn_array_sub_inplace(ulong* op1, const ulong* op2, ulong len,
                          const zn_mod_t mod)
{
   if (zn_mod_is_slim(mod))
   {
      // slim version
      // (unrolled)
      for (; len >= 4; len -= 4)
      {
         *op1 = zn_mod_sub_slim(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub_slim(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub_slim(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub_slim(*op1, *op2, mod);
         op1++; op2++;
      }
      for (; len; len--)
      {
         *op1 = zn_mod_sub_slim(*op1, *op2, mod);
         op1++; op2++;
      }
   }
   else
   {
      // non-slim version
      // (unrolled)
      for (; len >= 4; len -= 4)
      {
         *op1 = zn_mod_sub(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub(*op1, *op2, mod);
         op1++; op2++;
         *op1 = zn_mod_sub(*op1, *op2, mod);
         op1++; op2++;
      }
      for (; len; len--)
      {
         *op1 = zn_mod_sub(*op1, *op2, mod);
         op1++; op2++;
      }
   }
}



/* ============================================================================

     inplace zn_pmf_t butterflies

============================================================================ */


/*
   In the following routines, we work with the "relative bias" between the
   two inputs, i.e. b = difference between bias of op1 and op2. This allows
   us to avoid unnecessary normalisation steps; we just add and subtract
   directly into the correct memory locations.
*/


void zn_pmf_bfly(zn_pmf_t op1, zn_pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong b = op2[0] - op1[0];
   
   if (b & M)
   {
      // bias is in [M, 2M) mod 2M
      b &= (M - 1);
      
      // butterfly on op1[0, b) and op2[M-b, M)
      zn_array_bfly_inplace(op1+1, op2+1+M-b, b, mod);
      // butterfly on op1[b, M) and op2[0, M-b)
      zn_array_bfly_inplace(op2+1, op1+1+b, M-b, mod);
   }
   else
   {
      // bias is in [0, M) mod 2M
      b &= (M - 1);
      
      // butterfly on op1[0, b) and op2[M-b, M)
      zn_array_bfly_inplace(op2+1+M-b, op1+1, b, mod);
      // butterfly on op1[b, M) and op2[0, M-b)
      zn_array_bfly_inplace(op1+1+b, op2+1, M-b, mod);
   }
}


void zn_pmf_add(zn_pmf_t op1, const zn_pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong b = op2[0] - op1[0];

   if (b & M)
   {
      // bias is in [M, 2M) mod 2M
      b &= (M - 1);
      
      // add op2[M-b, M) to op1[0, b)
      zn_array_add_inplace(op1+1, op2+1+M-b, b, mod);
      // subtract op2[0, M-b) from op1[b, M)
      zn_array_sub_inplace(op1+1+b, op2+1, M-b, mod);
   }
   else
   {
      // bias is in [0, M) mod 2M
      b &= (M - 1);

      // subtract op2[M-b, M) from op1[0, b)
      zn_array_sub_inplace(op1+1, op2+1+M-b, b, mod);
      // add op2[0, M-b) to op1[b, M)
      zn_array_add_inplace(op1+1+b, op2+1, M-b, mod);
   }
}


void zn_pmf_sub(zn_pmf_t op1, const zn_pmf_t op2, ulong M, const zn_mod_t mod)
{
   ulong b = op2[0] - op1[0];

   if (b & M)
   {
      // bias is in [M, 2M) mod 2M
      b &= (M - 1);
      
      // subtract op2[M-b, M) from op1[0, b)
      zn_array_sub_inplace(op1+1, op2+1+M-b, b, mod);
      // add op2[0, M-b) to op1[b, M)
      zn_array_add_inplace(op1+1+b, op2+1, M-b, mod);
   }
   else
   {
      // bias is in [0, M) mod 2M
      b &= (M - 1);
      
      // add op2[M-b, M) to op1[0, b)
      zn_array_add_inplace(op1+1, op2+1+M-b, b, mod);
      // subtract op2[0, M-b) from op1[b, M)
      zn_array_sub_inplace(op1+1+b, op2+1, M-b, mod);
   }
}



/* ============================================================================

     zn_pmf_vec stuff

============================================================================ */


void zn_pmf_vec_init(zn_pmf_vec_t res, unsigned lgK, ptrdiff_t skip,
                     unsigned lgM, const zn_mod_t mod)
{
   ZNP_ASSERT(skip >= (1UL << lgM) + 1);
   
   res->lgK = lgK;
   res->lgM = lgM;
   res->skip = skip;
   res->K = 1UL << lgK;
   res->M = 1UL << lgM;
   res->mod = mod;
   res->data = (ulong*) malloc(sizeof(ulong) * skip * res->K);
}


void zn_pmf_vec_init_nussbaumer(zn_pmf_vec_t res, unsigned lgL,
                                const zn_mod_t mod)
{
   unsigned lgK, lgM;
   nussbaumer_params(&lgK, &lgM, lgL);
   zn_pmf_vec_init(res, lgK, (1UL << lgM) + 1, lgM, mod);
}


void zn_pmf_vec_clear(zn_pmf_vec_t op)
{
   free(op->data);
}


void zn_pmf_vec_set(zn_pmf_vec_t res, const zn_pmf_vec_t op)
{
   ZNP_ASSERT(zn_pmf_vec_compatible(res, op));
   
   ulong i;
   for (i = 0; i < op->K; i++)
      zn_pmf_set(res->data + i * res->skip, op->data + i * op->skip, op->M);
}


void zn_pmf_vec_scalar_mul(zn_pmf_vec_t op, ulong len, ulong x)
{
   ZNP_ASSERT(len <= op->K);

   ulong i;
   zn_pmf_t ptr = op->data;
   for (i = 0; i < len; i++, ptr += op->skip)
      zn_pmf_scalar_mul(ptr, op->M, x, op->mod);
}


ulong zn_pmf_vec_mul_get_fudge(unsigned lgM, int squaring, const zn_mod_t mod)
{
   int use_nussbaumer = (lgM >= (squaring
                     ? tuning_info[mod->bits].nuss_sqr_crossover
                     : tuning_info[mod->bits].nuss_mul_crossover));
                     
   if (use_nussbaumer)
      return nussbaumer_mul_get_fudge(lgM, squaring, mod);
   else
      return _zn_array_mul_get_fudge(1UL << lgM, 1UL << lgM, squaring, mod);
}


void zn_pmf_vec_mul(zn_pmf_vec_t res, const zn_pmf_vec_t op1,
                    const zn_pmf_vec_t op2, ulong length,
                    int special_first_two)
{
   ZNP_ASSERT(res->mod->n & 1);
   ZNP_ASSERT(zn_pmf_vec_compatible(res, op1));
   ZNP_ASSERT(zn_pmf_vec_compatible(res, op2));
   ZNP_ASSERT(res->M >= 2 || !special_first_two);

   zn_pmf_t res_ptr = res->data;
   zn_pmf_const_t op1_ptr = op1->data;
   zn_pmf_const_t op2_ptr = op2->data;

   const zn_mod_struct* mod = res->mod;
   
   ulong M = op1->M;
   unsigned lgM = op1->lgM;
   
   int squaring = (op1 == op2);
   
   // use nussbaumer algorithm if the pointwise mults are large enough
   int use_nussbaumer = (lgM >= (squaring
                     ? tuning_info[mod->bits].nuss_sqr_crossover
                     : tuning_info[mod->bits].nuss_mul_crossover));

   // scratch space for nussbaumer multiplications
   zn_pmf_vec_t temp1, temp2;
   unsigned nuss_lgK;

   if (use_nussbaumer)
   {
      zn_pmf_vec_init_nussbaumer(temp1, lgM, mod);
      zn_pmf_vec_init_nussbaumer(temp2, lgM, mod);
      nuss_lgK = temp1->lgK;
   }

   ulong i = 0;

   if (special_first_two)
   {
      ZNP_FASTALLOC(temp, ulong, 6624, 2*M);
      
      // need to adjust the fudge factors, so that the fudge factor for these
      // first two special products matches up with the fudge factors for the
      // remaining products
      ulong fixup;
      fixup = use_nussbaumer ? nussbaumer_mul_get_fudge(lgM, squaring, mod)
                             : _zn_array_mul_get_fudge(M, M, squaring, mod);
      fixup = zn_mod_mul(_zn_array_mul_get_fudge(M/2, M/2, squaring, mod),
                         zn_mod_invert(fixup, mod), mod);

      // length M/2 multiplications
      for (; i < 2 && i < length; i++, res_ptr += res->skip,
             op1_ptr += op1->skip, op2_ptr += op2->skip)
      {
         // add biases
         res_ptr[0] = op1_ptr[0] + op2_ptr[0];

         // do the actual multiplication
         _zn_array_mul(temp, op1_ptr + 1, M/2, op2_ptr + 1, M/2, 1, mod);

         // apply the fudge factor
         if (fixup == 1)
            zn_array_copy(res_ptr + 1, temp, M - 1);
         else
            zn_array_scalar_mul(res_ptr + 1, temp, M - 1, fixup, mod);
         res_ptr[M] = 0;
      }
      
      ZNP_FASTFREE(temp);
   }
   
   if (use_nussbaumer)
   {
      for (; i < length; i++, res_ptr += res->skip,
             op1_ptr += op1->skip, op2_ptr += op2->skip)
      {
         // add biases
         res_ptr[0] = op1_ptr[0] + op2_ptr[0];

         nussbaumer_mul(res_ptr + 1, op1_ptr + 1, op2_ptr + 1, temp1, temp2);
      }

      zn_pmf_vec_clear(temp2);
      zn_pmf_vec_clear(temp1);
   }
   else
   {
      // scratch space for KS negacyclic multiplication
      ZNP_FASTALLOC(temp, ulong, 6624, 2*M);
      temp[2*M - 1] = 0;

      for (; i < length; i++, res_ptr += res->skip,
             op1_ptr += op1->skip, op2_ptr += op2->skip)
      {
         // add biases
         res_ptr[0] = op1_ptr[0] + op2_ptr[0];
         
         // ordinary multiplication...
         _zn_array_mul(temp, op1_ptr + 1, M, op2_ptr + 1, M, 1, mod);
         // ... negacyclic reduction
         zn_array_sub(res_ptr + 1, temp, temp + M, M, mod);
      }

      ZNP_FASTFREE(temp);
   }
}



void zn_pmf_vec_reverse(zn_pmf_vec_t op, ulong length)
{
   op->data += op->skip * (length - 1);
   op->skip = -op->skip;
}



// end of file ****************************************************************
