/*
   midmul_fft-test.c:  test code for functions in midmul_fft.c
   
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

#include "support.h"
#include "zn_poly_internal.h"


int testcase_zn_pmf_vec_fft_transposed_small(
         unsigned lgK, unsigned lgM, ulong length, ulong nonzero, ulong twist,
         const zn_mod_t mod)
{
   int success = 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i, j, k;

   ptrdiff_t skip = M + 1;

   // ===================================
   // first check linearity, i.e. that ax + by gets mapped to the right thing

   {
      zn_pmf_vec_t X, Y, A, B, TX, TY;

      zn_pmf_vec_init(X, lgK, skip, lgM, mod);
      zn_pmf_vec_init(Y, lgK, skip, lgM, mod);
      zn_pmf_vec_init(A, lgK, skip, lgM, mod);
      zn_pmf_vec_init(B, lgK, skip, lgM, mod);
      zn_pmf_vec_init(TX, lgK, skip, lgM, mod);
      zn_pmf_vec_init(TY, lgK, skip, lgM, mod);

      // generate random X, Y, A, B
      
      for (i = 0; i < K; i++)
      {
         X->data[i*skip] = random_ulong(2*M);
         Y->data[i*skip] = random_ulong(2*M);
         for (j = 1; j <= M; j++)
         {
            X->data[i*skip + j] = random_ulong(mod->n);
            Y->data[i*skip + j] = random_ulong(mod->n);
         }
      }
      
      A->data[0] = random_ulong(2*M);
      B->data[0] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
      {
         A->data[j] = random_ulong(mod->n);
         B->data[j] = random_ulong(mod->n);
      }

      for (i = 1; i < K; i++)
      {
         zn_pmf_set(A->data + A->skip * i, A->data, M);
         zn_pmf_set(B->data + B->skip * i, B->data, M);
      }
      
      // transform X and Y (after throwing in random ignorable crap)
      
      zn_pmf_vec_set(TX, X);
      zn_pmf_vec_set(TY, Y);
      
      for (i = length; i < K; i++)
      {
         TX->data[i*skip] = random_ulong(2*M);
         TY->data[i*skip] = random_ulong(2*M);
         for (j = 1; j <= M; j++)
         {
            TX->data[i*skip + j] = random_ulong(mod->n);
            TY->data[i*skip + j] = random_ulong(mod->n);
         }
      }

      zn_pmf_vec_fft_transposed_small(TX, length, nonzero, twist);
      zn_pmf_vec_fft_transposed_small(TY, length, nonzero, twist);

      // form linear combination of TX and TY

      zn_pmf_vec_mul(TX, TX, A, nonzero, 0);
      zn_pmf_vec_mul(TY, TY, B, nonzero, 0);
      for (i = 0; i < nonzero; i++)
         zn_pmf_add(TX->data + TX->skip * i, TY->data + TY->skip * i, M, mod);

      // form linear combination of X and Y

      zn_pmf_vec_mul(X, X, A, length, 0);
      zn_pmf_vec_mul(Y, Y, B, length, 0);
      for (i = 0; i < length; i++)
         zn_pmf_add(X->data + X->skip * i, Y->data + Y->skip * i, M, mod);

      // transform linear combination of X and Y

      zn_pmf_vec_fft_transposed_small(X, length, nonzero, twist);

      // compare results
      
      for (i = 0; i < nonzero; i++)
      {
         zn_pmf_sub(X->data + X->skip * i, TX->data + TX->skip * i, M, mod);
         for (j = 1; j <= M; j++)
            success = success && (X->data[i * skip + j] == 0);
      }
      
      zn_pmf_vec_clear(X);
      zn_pmf_vec_clear(Y);
      zn_pmf_vec_clear(TX);
      zn_pmf_vec_clear(TY);
      zn_pmf_vec_clear(A);
      zn_pmf_vec_clear(B);
   }
   
   // ===================================
   // now check that the matrix of the transposed FFT is really the transpose
   // of the matrix of the ordinary FFT

   {
      zn_pmf_vec_t* X = (zn_pmf_vec_t*)
                                malloc(nonzero * sizeof(zn_pmf_vec_t));
      for (i = 0; i < nonzero; i++)
         zn_pmf_vec_init(X[i], lgK, skip, lgM, mod);

      zn_pmf_vec_t* Y = (zn_pmf_vec_t*)
                                malloc(length * sizeof(zn_pmf_vec_t));
      for (i = 0; i < length; i++)
         zn_pmf_vec_init(Y[i], lgK, skip, lgM, mod);
         
      // compute images of basis vectors under FFT

      for (i = 0; i < nonzero; i++)
      for (j = 0; j < nonzero; j++)
      {
         for (k = 0; k <= M; k++)
            X[i]->data[j*skip + k] = 0;

         X[i]->data[j*skip + 1] = (i == j);
      }
      
      for (i = 0; i < nonzero; i++)
         zn_pmf_vec_fft(X[i], length, nonzero, twist);

      // compute images of basis vectors under transposed FFT

      for (i = 0; i < length; i++)
      for (j = 0; j < length; j++)
      {
         for (k = 0; k <= M; k++)
            Y[i]->data[j*skip + k] = 0;

         Y[i]->data[j*skip + 1] = (i == j);
      }
      
      for (i = 0; i < length; i++)
         zn_pmf_vec_fft_transposed(Y[i], length, nonzero, twist);
      
      // check that they are transposes of each other

      for (i = 0; i < nonzero; i++)
      for (j = 0; j < length; j++)
      {
         zn_pmf_sub(X[i]->data + j*skip, Y[j]->data + i*skip, M, mod);
         for (k = 1; k <= M; k++)
            success = success && (X[i]->data[j*skip + k] == 0);
      }

      for (i = 0; i < nonzero; i++)
         zn_pmf_vec_clear(X[i]);
      for (i = 0; i < length; i++)
         zn_pmf_vec_clear(Y[i]);
      free(Y);
      free(X);
   }

   return success;
}



int test_zn_pmf_vec_fft_transposed_small()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM;
   ulong nonzero, length, twist;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgM = lgK ? (lgK - 1) : 0; lgM < lgK + 2 && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;
      for (twist = 0; twist < 2*M/K && success; twist++)
      for (length = 0; length <= K && success; length++)
      for (nonzero = 0; nonzero <= K && success; nonzero++)
      {
         bits = test_bitsizes[i];
         
         for (trial = 0; trial < 5; trial++)
         {
            zn_mod_t mod;
            zn_mod_init(mod, random_modulus(bits, 1));
         
            success = success && testcase_zn_pmf_vec_fft_transposed_small(
                         lgK, lgM, length, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



int testcase_zn_pmf_vec_fft_transposed_factor(
         unsigned lgK, unsigned lgM, unsigned lgT,
         ulong length, ulong nonzero, ulong twist, const zn_mod_t mod)
{
   zn_pmf_vec_t A, B;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i, j, k;

   ptrdiff_t skip = M + 1;

   zn_pmf_vec_init(A, lgK, skip, lgM, mod);
   zn_pmf_vec_init(B, lgK, skip, lgM, mod);

   // create random input
   for (i = 0; i < K; i++)
   {
      A->data[i*skip] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
         A->data[i*skip + j] = random_ulong(mod->n);
   }
   
   // make a copy
   zn_pmf_vec_set(B, A);
   
   // put random crap in B to check that factoring algorithm ignores it
   for (i = length; i < K; i++)
   {
      B->data[i*skip] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
         B->data[i*skip + j] = random_ulong(mod->n);
   }
   
   // run transposed FFTs using "factor" and "small" algorithms
   zn_pmf_vec_fft_transposed_small(A, length, nonzero, twist);
   zn_pmf_vec_fft_transposed_factor(B, lgT, length, nonzero, twist);

   // compare results
   int success = 1;
   for (i = 0; i < nonzero; i++)
   {
      zn_pmf_sub(B->data + i * skip, A->data + i * skip, M, mod);
      for (j = 1; j <= M; j++)
         success = success && (B->data[i * skip + j] == 0);
   }
   
   zn_pmf_vec_clear(B);
   zn_pmf_vec_clear(A);

   return success;
}



int test_zn_pmf_vec_fft_transposed_factor()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM, lgT;
   ulong nonzero, length, twist;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = 1; lgT < lgK && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0; lgM < lgK + 2 && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;
      for (twist = 0; twist < 2*M/K && success; twist++)
      for (length = 0; length <= K && success; length++)
      for (nonzero = 0; nonzero <= K && success; nonzero++)
      {
         bits = test_bitsizes[i];
         
         for (trial = 0; trial < 5; trial++)
         {
            zn_mod_t mod;
            zn_mod_init(mod, random_modulus(bits, 1));
         
            success = success && testcase_zn_pmf_vec_fft_transposed_factor(
                         lgK, lgM, lgT, length, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



int testcase_zn_pmf_vec_ifft_transposed_small(
         unsigned lgK, unsigned lgM, ulong length, int forward,
         ulong nonzero, ulong twist, const zn_mod_t mod)
{
   int success = 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i, j, k;

   ptrdiff_t skip = M + 1;

   // ===================================
   // first check linearity, i.e. that ax + by gets mapped to the right thing

   {
      zn_pmf_vec_t X, Y, A, B, TX, TY;

      zn_pmf_vec_init(X, lgK, skip, lgM, mod);
      zn_pmf_vec_init(Y, lgK, skip, lgM, mod);
      zn_pmf_vec_init(A, lgK, skip, lgM, mod);
      zn_pmf_vec_init(B, lgK, skip, lgM, mod);
      zn_pmf_vec_init(TX, lgK, skip, lgM, mod);
      zn_pmf_vec_init(TY, lgK, skip, lgM, mod);

      // generate random X, Y, A, B
      
      for (i = 0; i < K; i++)
      {
         X->data[i*skip] = random_ulong(2*M);
         Y->data[i*skip] = random_ulong(2*M);
         for (j = 1; j <= M; j++)
         {
            X->data[i*skip + j] = random_ulong(mod->n);
            Y->data[i*skip + j] = random_ulong(mod->n);
         }
      }
      
      A->data[0] = random_ulong(2*M);
      B->data[0] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
      {
         A->data[j] = random_ulong(mod->n);
         B->data[j] = random_ulong(mod->n);
      }

      for (i = 1; i < K; i++)
      {
         zn_pmf_set(A->data + A->skip * i, A->data, M);
         zn_pmf_set(B->data + B->skip * i, B->data, M);
      }
      
      // transform X and Y (after throwing in random ignorable crap)
      
      zn_pmf_vec_set(TX, X);
      zn_pmf_vec_set(TY, Y);
      
      for (i = length + forward; i < K; i++)
      {
         TX->data[i*skip] = random_ulong(2*M);
         TY->data[i*skip] = random_ulong(2*M);
         for (j = 1; j <= M; j++)
         {
            TX->data[i*skip + j] = random_ulong(mod->n);
            TY->data[i*skip + j] = random_ulong(mod->n);
         }
      }

      zn_pmf_vec_ifft_transposed_small(TX, length, forward, nonzero, twist);
      zn_pmf_vec_ifft_transposed_small(TY, length, forward, nonzero, twist);

      // form linear combination of TX and TY

      zn_pmf_vec_mul(TX, TX, A, nonzero, 0);
      zn_pmf_vec_mul(TY, TY, B, nonzero, 0);
      for (i = 0; i < nonzero; i++)
         zn_pmf_add(TX->data + TX->skip * i, TY->data + TY->skip * i, M, mod);

      // form linear combination of X and Y

      zn_pmf_vec_mul(X, X, A, length + forward, 0);
      zn_pmf_vec_mul(Y, Y, B, length + forward, 0);
      for (i = 0; i < length + forward; i++)
         zn_pmf_add(X->data + X->skip * i, Y->data + Y->skip * i, M, mod);

      // transform linear combination of X and Y

      zn_pmf_vec_ifft_transposed_small(X, length, forward, nonzero, twist);

      // compare results
      
      for (i = 0; i < nonzero; i++)
      {
         zn_pmf_sub(X->data + X->skip * i, TX->data + TX->skip * i, M, mod);
         for (j = 1; j <= M; j++)
            success = success && (X->data[i * skip + j] == 0);
      }
      
      zn_pmf_vec_clear(X);
      zn_pmf_vec_clear(Y);
      zn_pmf_vec_clear(TX);
      zn_pmf_vec_clear(TY);
      zn_pmf_vec_clear(A);
      zn_pmf_vec_clear(B);
   }

   // ===================================
   // now check that the matrix of the transposed IFFT is really the transpose
   // of the matrix of the ordinary IFFT

   {
      zn_pmf_vec_t* X = (zn_pmf_vec_t*)
                                malloc(nonzero * sizeof(zn_pmf_vec_t));
      for (i = 0; i < nonzero; i++)
         zn_pmf_vec_init(X[i], lgK, skip, lgM, mod);

      zn_pmf_vec_t* Y = (zn_pmf_vec_t*)
                             malloc((length + forward) * sizeof(zn_pmf_vec_t));
      for (i = 0; i < length + forward; i++)
         zn_pmf_vec_init(Y[i], lgK, skip, lgM, mod);
         
      // compute images of basis vectors under FFT

      for (i = 0; i < nonzero; i++)
      for (j = 0; j < nonzero; j++)
      {
         for (k = 0; k <= M; k++)
            X[i]->data[j*skip + k] = 0;

         X[i]->data[j*skip + 1] = (i == j);
      }
      
      for (i = 0; i < nonzero; i++)
         zn_pmf_vec_ifft(X[i], length, forward, nonzero, twist);

      // compute images of basis vectors under transposed FFT

      for (i = 0; i < length + forward; i++)
      for (j = 0; j < length + forward; j++)
      {
         for (k = 0; k <= M; k++)
            Y[i]->data[j*skip + k] = 0;

         Y[i]->data[j*skip + 1] = (i == j);
      }
      
      for (i = 0; i < length + forward; i++)
         zn_pmf_vec_ifft_transposed(Y[i], length, forward, nonzero, twist);
      
      // check that they are transposes of each other

      for (i = 0; i < nonzero; i++)
      for (j = 0; j < length + forward; j++)
      {
         zn_pmf_sub(X[i]->data + j*skip, Y[j]->data + i*skip, M, mod);
         for (k = 1; k <= M; k++)
            success = success && (X[i]->data[j*skip + k] == 0);
      }

      for (i = 0; i < nonzero; i++)
         zn_pmf_vec_clear(X[i]);
      for (i = 0; i < length + forward; i++)
         zn_pmf_vec_clear(Y[i]);
      free(Y);
      free(X);
   }

   return success;
}


int test_zn_pmf_vec_ifft_transposed_small()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM;
   ulong nonzero, length, twist;
   int forward;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgM = lgK ? (lgK - 1) : 0; lgM < lgK + 2 && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;
      for (twist = 0; twist < 2*M/K && success; twist++)
      for (length = 0; length <= K && success; length++)
      for (nonzero = length; nonzero <= K && success; nonzero++)
      for (forward = 0; forward < 2 && success; forward++)
      {
         if (forward + length > K)
            continue;

         bits = test_bitsizes[i];
         
         for (trial = 0; trial < 5; trial++)
         {
            zn_mod_t mod;
            zn_mod_init(mod, random_modulus(bits, 1));

            success = success && testcase_zn_pmf_vec_ifft_transposed_small(
                         lgK, lgM, length, forward, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



int testcase_zn_pmf_vec_ifft_transposed_factor(
         unsigned lgK, unsigned lgM, unsigned lgT,
         ulong length, int forward, ulong nonzero, ulong twist,
         const zn_mod_t mod)
{
   zn_pmf_vec_t A, B;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i, j, k;

   ptrdiff_t skip = M + 1;

   zn_pmf_vec_init(A, lgK, skip, lgM, mod);
   zn_pmf_vec_init(B, lgK, skip, lgM, mod);

   // create random input
   for (i = 0; i < K; i++)
   {
      A->data[i*skip] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
         A->data[i*skip + j] = random_ulong(mod->n);
   }
   
   // make a copy
   zn_pmf_vec_set(B, A);
   
   // put random crap in B to check that factoring algorithm ignores it
   for (i = length + forward; i < K; i++)
   {
      B->data[i*skip] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
         B->data[i*skip + j] = random_ulong(mod->n);
   }
   
   // run transposed IFFTs using "factor" and "small" algorithms
   zn_pmf_vec_ifft_transposed_small(A, length, forward, nonzero, twist);
   zn_pmf_vec_ifft_transposed_factor(B, lgT, length, forward, nonzero, twist);
   
   // compare results
   int success = 1;
   for (i = 0; i < nonzero; i++)
   {
      zn_pmf_sub(B->data + i * skip, A->data + i * skip, M, mod);
      for (j = 1; j <= M; j++)
         success = success && (B->data[i * skip + j] == 0);
   }
   
   zn_pmf_vec_clear(B);
   zn_pmf_vec_clear(A);

   return success;
}



int test_zn_pmf_vec_ifft_transposed_factor()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM, lgT;
   ulong nonzero, length, twist;
   int forward;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = 1; lgT < lgK && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0; lgM < lgK + 3 && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;
      for (twist = 0; twist < 2*M/K && success; twist++)
      for (length = 0; length <= K && success; length++)
      for (nonzero = length; nonzero <= K && success; nonzero++)
      for (forward = 0; forward < 2 && success; forward++)
      {
         if (forward + length > K)
            continue;

         bits = test_bitsizes[i];
         
         for (trial = 0; trial < 3; trial++)
         {
            zn_mod_t mod;
            zn_mod_init(mod, random_modulus(bits, 1));
         
            success = success && testcase_zn_pmf_vec_ifft_transposed_factor(
                       lgK, lgM, lgT, length, forward, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



/*
   Tests zn_array_midmul_fft, for given len1, len2, modulus.

   If scale_extra is set, zn_array_midmul_fft() gets called with a random
   scaling factor, otherwise gets called with scale == 1.

   Returns 1 on success.
*/
int testcase_zn_array_midmul_fft(size_t len1, size_t len2, int scale_extra,
                                 const zn_mod_t mod)
{
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len1);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * len2);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * (len1 - len2 + 1));
   ulong* compare = (ulong*) malloc(sizeof(ulong) * (len1 - len2 + 1));

   // generate random polys
   size_t i;
   for (i = 0; i < len1; i++)
      buf1[i] = random_ulong(mod->n);
   for (i = 0; i < len2; i++)
      buf2[i] = random_ulong(mod->n);
      
   ulong scale = scale_extra ? random_ulong(mod->n) : 1;

   // compare target implementation against reference implementation
   ref_zn_array_midmul(correct, buf1, len1, buf2, len2, mod);
   ref_zn_array_scalar_mul(correct, correct, len1 - len2 + 1, scale, mod);

   zn_array_midmul_fft(compare, buf1, len1, buf2, len2, scale, mod);
   ulong correct_fudge = zn_array_midmul_fft_get_fudge(len1, len2, mod);
   ref_zn_array_scalar_mul(compare, compare, len1 - len2 + 1,
                           correct_fudge, mod);
   
   int success = !zn_array_cmp(correct, compare, len1 - len2 + 1);

   free(compare);
   free(correct);
   free(buf2);
   free(buf1);
   
   return success;
}


/*
   tests zn_array_midmul_fft() on a range of input cases
*/
int test_zn_array_midmul_fft()
{
   int success = 1;
   int i, bits, trial, scale_extra;
   size_t len1, len2;

   // first try a dense range of "small" problems

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (len2 = 1; len2 <= 50 && success; len2++)
   for (len1 = len2; len1 <= 50 && success; len1++)
   for (scale_extra = 0; scale_extra <= 1 && success; scale_extra++)
   for (trial = 0; trial < 3 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(test_bitsizes[i], 1));

      success = success &&
                   testcase_zn_array_midmul_fft(len1, len2, scale_extra, mod);
      
      zn_mod_clear(mod);
   }

   // now try some random larger problems
   // and temporarily change the nussbaumer crossovers so we use that
   // code sometimes
   
   ulong crossover;
   
   for (i = 0; i < num_test_bitsizes && success; i++)
   for (scale_extra = 0; scale_extra <= 1 && success; scale_extra++)
   for (crossover = 2; crossover <= 8; crossover++)
   {
      bits = test_bitsizes[i];

      ulong save_crossover = tuning_info[bits].nuss_mul_crossover;
      tuning_info[bits].nuss_mul_crossover = crossover;
   
      len1 = random_ulong(10000) + 1;
      len2 = random_ulong(10000) + 1;
      
      if (len1 < len2)
      {
         ulong temp = len1;
         len1 = len2;
         len2 = temp;
      }
      
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 1));

      success = success &&
                   testcase_zn_array_midmul_fft(len1, len2, scale_extra, mod);
      
      zn_mod_clear(mod);
      tuning_info[bits].nuss_mul_crossover = save_crossover;
   }

   return success;
}


// end of file ****************************************************************
