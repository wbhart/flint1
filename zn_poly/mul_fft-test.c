/*
   mul_fft-test.c:  test code for functions in mul_fft.c and mul_fft_dft.c
   
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


/*
   If lgT == 0, this tests zn_pmf_vec_fft_small.
   If lgT > 0, it tests zn_pmf_vec_fft_factor.
*/
int testcase_zn_pmf_vec_fft_small_or_factor(
         unsigned lgK, unsigned lgM, unsigned lgT,
         ulong length, ulong nonzero, ulong twist, const zn_mod_t mod)
{
   zn_pmf_vec_t A, B;
   
   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;

   ptrdiff_t skip = M + 1;
   ulong i, j;

   zn_pmf_vec_init(A, lgK, skip, lgM, mod);
   zn_pmf_vec_init(B, lgK, skip, lgM, mod);

   // create random a_i's, with zero padding
   for (i = nonzero; i < K; i++)
      zn_pmf_zero(A->data + i*skip, M);
   for (i = 0; i < nonzero; i++)
      for (j = 1; j <= M; j++)
         A->data[i*skip + j] = random_ulong(mod->n);
   for (i = 0; i < K; i++)
      A->data[i*skip] = random_ulong(2*M);
      
   // run FFT using simple iterative algorithm
   zn_pmf_vec_set(B, A);
   zn_pmf_vec_fft_notrunc_iterative(B, twist);

   // make sure truncated FFT has to deal with random crap
   for (i = nonzero; i < K; i++)
   {
      A->data[i*skip] = random_ulong(2*M);
      for (j = 1; j <= M; j++)
         A->data[i*skip + j] = random_ulong(mod->n);
   }

   // try truncated FFT
   if (lgT > 0)
      zn_pmf_vec_fft_factor(A, lgT, length, nonzero, twist);
   else
      zn_pmf_vec_fft_small(A, length, nonzero, twist);
      
   // compare results
   int success = 1;
   for (i = 0; i < length; i++)
   {
      zn_pmf_sub(B->data + i * skip, A->data + i * skip, M, mod);
      for (j = 1; j <= M; j++)
         success = success && (B->data[i * skip + j] == 0);
   }

   zn_pmf_vec_clear(B);
   zn_pmf_vec_clear(A);
   
   return success;
}



int test_zn_pmf_vec_fft_small()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM;
   ulong nonzero, length, twist;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgM = lgK ? (lgK - 1) : 0; lgM < lgK + 3 && success; lgM++)
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
         
            success = success && testcase_zn_pmf_vec_fft_small_or_factor(
                         lgK, lgM, 0, length, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



int test_zn_pmf_vec_fft_factor()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM, lgT;
   ulong nonzero, length, twist;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = 1; lgT < lgK && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0; lgM < lgK + 3 && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;
      for (twist = 0; twist < 2*M/K && success; twist++)
      for (length = 0; length <= K && success; length++)
      for (nonzero = 0; nonzero <= K && success; nonzero++)
      {
         bits = test_bitsizes[i];
         
         for (trial = 0; trial < 3; trial++)
         {
            zn_mod_t mod;
            zn_mod_init(mod, random_modulus(bits, 1));
         
            success = success && testcase_zn_pmf_vec_fft_small_or_factor(
                       lgK, lgM, lgT, length, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



/*
   If lgT == 0, this tests zn_pmf_vec_ifft_small.
   If lgT > 0, it tests zn_pmf_vec_ifft_factor.
*/
int testcase_zn_pmf_vec_ifft_small_or_factor(
         unsigned lgK, unsigned lgM, unsigned lgT,
         ulong length, int forward, ulong nonzero, ulong twist,
         const zn_mod_t mod)
{
   zn_pmf_vec_t A, B, C;
   
   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong scale = zn_mod_reduce(K, mod);

   ptrdiff_t skip = M + 1;
   ulong i, j;

   zn_pmf_vec_init(A, lgK, skip, lgM, mod);
   zn_pmf_vec_init(B, lgK, skip, lgM, mod);
   zn_pmf_vec_init(C, lgK, skip, lgM, mod);

   // create random a_i's, with zero padding
   for (i = nonzero; i < K; i++)
      zn_pmf_zero(A->data + i*skip, M);
   for (i = 0; i < nonzero; i++)
      for (j = 1; j <= M; j++)
         A->data[i*skip + j] = random_ulong(mod->n);
   for (i = 0; i < K; i++)
      A->data[i*skip] = random_ulong(2*M);
      
   // run FFT
   zn_pmf_vec_set(B, A);
   zn_pmf_vec_fft(B, K, K, twist);
   zn_pmf_vec_set(C, B);
   
   // fill in missing data, plus junk where the implied zeroes should be
   for (i = length; i < nonzero; i++)
   {
      zn_pmf_set(C->data + i * skip, A->data + i * skip, M);
      zn_pmf_scalar_mul(C->data + i * skip, M, scale, mod);
   }
   for (i = nonzero; i < K; i++)
   {
      C->data[i*skip] = -1UL;
      for (j = 1; j <= M; j++)
         C->data[i*skip + j] = random_ulong(mod->n);
   }
   
   // try IFFT
   if (lgT > 0)
      zn_pmf_vec_ifft_factor(C, lgT, length, forward, nonzero, twist);
   else
      zn_pmf_vec_ifft_small(C, length, forward, nonzero, twist);
   
   // compare results
   int success = 1;
   for (i = 0; i < length; i++)
   {
      zn_pmf_scalar_mul(A->data + i * skip, M, scale, mod);
      zn_pmf_sub(C->data + i * skip, A->data + i * skip, M, mod);
   }
   if (forward)
      zn_pmf_sub(C->data + i * skip, B->data + i * skip, M, mod);

   for (i = 0; i < length + forward; i++)
      for (j = 1; j <= M; j++)
         success = success && (C->data[i * skip + j] == 0);

   zn_pmf_vec_clear(C);
   zn_pmf_vec_clear(B);
   zn_pmf_vec_clear(A);
   
   return success;
}


int test_zn_pmf_vec_ifft_small()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgK, lgM;
   ulong nonzero, length, twist;
   int forward;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
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
         
         for (trial = 0; trial < 5; trial++)
         {
            zn_mod_t mod;
            zn_mod_init(mod, random_modulus(bits, 1));
         
            success = success && testcase_zn_pmf_vec_ifft_small_or_factor(
                         lgK, lgM, 0, length, forward, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}



int test_zn_pmf_vec_ifft_factor()
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
         
            success = success && testcase_zn_pmf_vec_ifft_small_or_factor(
                       lgK, lgM, lgT, length, forward, nonzero, twist, mod);
            
            zn_mod_clear(mod);
         }
      }
   }
   
   return success;
}


/*
   Tests zn_array_mul_fft, for given len1, len2, modulus.
   
   If scale_extra is set, zn_array_mul_fft() gets called with a random
   scaling factor, otherwise gets called with scale == 1.
   
   Returns 1 on success.
*/
int testcase_zn_array_mul_fft(size_t len1, size_t len2, int scale_extra,
                              const zn_mod_t mod)
{
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len1);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * len2);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * (len1 + len2 - 1));
   ulong* compare = (ulong*) malloc(sizeof(ulong) * (len1 + len2 - 1));

   // generate random polys
   size_t i;
   for (i = 0; i < len1; i++)
      buf1[i] = random_ulong(mod->n);
   for (i = 0; i < len2; i++)
      buf2[i] = random_ulong(mod->n);

   ulong scale = scale_extra ? random_ulong(mod->n) : 1;

   // compare target implementation against reference implementation
   ref_zn_array_mul(correct, buf1, len1, buf2, len2, mod);
   ref_zn_array_scalar_mul(correct, correct, len1 + len2 - 1, scale, mod);
      
   zn_array_mul_fft(compare, buf1, len1, buf2, len2, scale, mod);
   ulong correct_fudge = zn_array_mul_fft_get_fudge(len1, len2, 0, mod);
   ref_zn_array_scalar_mul(compare, compare, len1 + len2 - 1,
                           correct_fudge, mod);
   
   int success = !zn_array_cmp(correct, compare, len1 + len2 - 1);
   
   free(compare);
   free(correct);
   free(buf2);
   free(buf1);
   
   return success;
}



/*
   tests zn_array_mul_fft() on a range of input cases
*/
int test_zn_array_mul_fft()
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
                  testcase_zn_array_mul_fft(len1, len2, scale_extra, mod);
      
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
                   testcase_zn_array_mul_fft(len1, len2, scale_extra, mod);
      
      zn_mod_clear(mod);
      tuning_info[bits].nuss_mul_crossover = save_crossover;
   }

   return success;
}



/*
   Tests zn_array_mul_fft for squaring, for given len and modulus.
   
   Parameter scale_extra has same meaning as for testcase_zn_array_mul_fft().
   
   Returns 1 on success.
*/
int testcase_zn_array_sqr_fft(size_t len, int scale_extra, const zn_mod_t mod)
{
   ulong* buf = (ulong*) malloc(sizeof(ulong) * len);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * (2*len - 1));
   ulong* compare = (ulong*) malloc(sizeof(ulong) * (2*len - 1));

   // generate random poly
   size_t i;
   for (i = 0; i < len; i++)
      buf[i] = random_ulong(mod->n);
      
   ulong scale = scale_extra ? random_ulong(mod->n) : 1;

   // compare target implementation against reference implementation
   ref_zn_array_mul(correct, buf, len, buf, len, mod);
   ref_zn_array_scalar_mul(correct, correct, 2*len - 1, scale, mod);

   zn_array_mul_fft(compare, buf, len, buf, len, scale, mod);
   ulong correct_fudge = zn_array_mul_fft_get_fudge(len, len, 1, mod);
   ref_zn_array_scalar_mul(compare, compare, 2*len - 1,
                           correct_fudge, mod);

   int success = !zn_array_cmp(correct, compare, 2*len - 1);
   
   free(compare);
   free(correct);
   free(buf);
   
   return success;
}


/*
   tests zn_array_mul_fft() for squaring on a range of input cases
*/
int test_zn_array_sqr_fft()
{
   int success = 1;
   int i, bits, trial, scale_extra;
   size_t len;

   // first try a dense range of "small" problems

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (len = 1; len <= 30 && success; len++)
   for (scale_extra = 0; scale_extra <= 1 && success; scale_extra++)
   for (trial = 0; trial < 10 && success; trial++)
   {
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(test_bitsizes[i], 1));
      
      success = success && testcase_zn_array_sqr_fft(len, scale_extra, mod);
      
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
   
      len = random_ulong(10000) + 1;
      
      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 1));

      success = success && testcase_zn_array_sqr_fft(len, scale_extra, mod);
      
      zn_mod_clear(mod);
      tuning_info[bits].nuss_mul_crossover = save_crossover;
   }

   return success;
}



/*
   Tests zn_array_mul_fft, for given len1, len2, lgT, modulus.
   Returns 1 on success.
*/
int testcase_zn_array_mul_fft_dft(size_t len1, size_t len2, unsigned lgT,
                                  const zn_mod_t mod)
{
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len1);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * len2);
   ulong* correct = (ulong*) malloc(sizeof(ulong) * (len1 + len2 - 1));
   ulong* compare = (ulong*) malloc(sizeof(ulong) * (len1 + len2 - 1));

   // generate random polys
   size_t i;
   for (i = 0; i < len1; i++)
      buf1[i] = random_ulong(mod->n);
   for (i = 0; i < len2; i++)
      buf2[i] = random_ulong(mod->n);
      
   // compare target implementation against reference implementation
   ref_zn_array_mul(correct, buf1, len1, buf2, len2, mod);
   zn_array_mul_fft_dft(compare, buf1, len1, buf2, len2, lgT, mod);
   int success = !zn_array_cmp(correct, compare, len1 + len2 - 1);
   
   free(compare);
   free(correct);
   free(buf2);
   free(buf1);
   
   return success;
}



/*
   tests zn_array_mul_fft_dft() on a range of input cases
*/
int test_zn_array_mul_fft_dft()
{
   int success = 1;
   int i, bits, trial;
   unsigned lgT;
   size_t len1, len2;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (len2 = 1; len2 <= 30 && success; len2++)
   for (len1 = len2; len1 <= 30 && success; len1++)
   for (lgT = 0; lgT < 5 && success; lgT++)
   for (trial = 0; trial < 3 && success; trial++)
   {
      bits = test_bitsizes[i];

      zn_mod_t mod;
      zn_mod_init(mod, random_modulus(bits, 1));
      
      success = success && testcase_zn_array_mul_fft_dft(len1, len2, lgT, mod);
      
      zn_mod_clear(mod);
   }
   
   return success;
}


// end of file ****************************************************************
