/*
   pmfvec_fft-test.c:  test code for functions in pmfvec_fft.c
   
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

#include "support.h"
#include "zn_poly_internal.h"


/*
   If lgT == 0, this tests pmfvec_fft_dc.
   If lgT > 0, it tests pmfvec_fft_huge.
*/
int
testcase_pmfvec_fft_dc_or_huge (unsigned lgK, unsigned lgM, unsigned lgT,
                                ulong n, ulong z, ulong t, const zn_mod_t mod)
{
   pmfvec_t A, B;
   
   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;

   ptrdiff_t skip = M + 1;
   ulong i;

   pmfvec_init (A, lgK, skip, lgM, mod);
   pmfvec_init (B, lgK, skip, lgM, mod);

   // create random a_i's, with zero padding
   for (i = z; i < K; i++)
      pmf_zero (A->data + i * skip, M);
   for (i = z; i < K; i++)
      A->data[i * skip] = random_ulong (2 * M);
   for (i = 0; i < z; i++)
      pmf_rand (A->data + i * skip, M, mod);
      
   // run FFT using simple iterative algorithm
   pmfvec_set (B, A);
   pmfvec_fft_basecase (B, t);

   // make sure truncated FFT has to deal with random crap
   for (i = z; i < K; i++)
      pmf_rand (A->data + i * skip, M, mod);

   // try truncated FFT
   if (lgT > 0)
      pmfvec_fft_huge (A, lgT, n, z, t);
   else
      pmfvec_fft_dc (A, n, z, t);
      
   // compare results
   int success = 1;
   for (i = 0; i < n; i++)
      success = success && !pmf_cmp (A->data + i * skip, B->data + i * skip,
                                     M, mod);

   pmfvec_clear (B);
   pmfvec_clear (A);
   
   return success;
}



/*
   Tests pmfvec_fft_dc (if huge == 0) or pmfvec_fft_huge (if huge == 1)
*/
int
test_pmfvec_fft_dc_or_huge (int huge, int quick)
{
   int success = 1;
   int i;
   unsigned lgK, lgM, lgT;
   ulong z, n, t;
   zn_mod_t mod;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = (huge ? 1 : 0); lgT < (huge ? lgK : 1) && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0;
        lgM < lgK + (quick ? 1 : 3) && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;

      for (t = 0; t < ZNP_MIN (2 * M / K, quick ? 2 : 1000) && success; t++)
      for (n = 1; n <= K && success; n++)
      for (z = 1; z <= K && success; z++)
      {
         zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
         success = success && testcase_pmfvec_fft_dc_or_huge
                                                (lgK, lgM, lgT, n, z, t, mod);
         zn_mod_clear (mod);
      }
   }
   
   return success;
}


int
test_pmfvec_fft_dc (int quick)
{
   return test_pmfvec_fft_dc_or_huge (0, quick);
}


int
test_pmfvec_fft_huge (int quick)
{
   return test_pmfvec_fft_dc_or_huge (1, quick);
}



/*
   If lgT == 0, this tests pmfvec_ifft_dc.
   If lgT > 0, it tests pmfvec_ifft_huge.
*/
int
testcase_pmfvec_ifft_dc_or_huge (unsigned lgK, unsigned lgM, unsigned lgT,
                                 ulong n, int fwd, ulong z, ulong t,
                                 const zn_mod_t mod)
{
   pmfvec_t A, B, C;
   
   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong x = zn_mod_reduce (K, mod);

   ptrdiff_t skip = M + 1;
   ulong i;

   pmfvec_init (A, lgK, skip, lgM, mod);
   pmfvec_init (B, lgK, skip, lgM, mod);
   pmfvec_init (C, lgK, skip, lgM, mod);

   // create random a_i's, with zero padding
   for (i = z; i < K; i++)
      pmf_zero (A->data + i * skip, M);
   for (i = z; i < K; i++)
      A->data[i * skip] = random_ulong (2 * M);
   for (i = 0; i < z; i++)
      pmf_rand (A->data + i * skip, M, mod);
      
   // run FFT
   pmfvec_set (B, A);
   pmfvec_fft (B, K, K, t);
   pmfvec_set (C, B);
   
   // fill in missing data, plus junk where the implied zeroes should be
   for (i = n; i < z; i++)
   {
      pmf_set (C->data + i * skip, A->data + i * skip, M);
      pmf_scalar_mul (C->data + i * skip, M, x, mod);
   }
   for (i = z; i < K; i++)
      pmf_rand (C->data + i * skip, M, mod);
   
   // try IFFT
   if (lgT > 0)
      pmfvec_ifft_huge (C, lgT, n, fwd, z, t);
   else
      pmfvec_ifft_dc (C, n, fwd, z, t);
   
   // compare results
   int success = 1;
   for (i = 0; i < n; i++)
      pmf_scalar_mul (A->data + i * skip, M, x, mod);
   for (i = 0; i < n; i++)
      success = success && !pmf_cmp (C->data + i * skip, A->data + i * skip,
                                     M, mod);
   if (fwd)
      success = success && !pmf_cmp (C->data + i * skip, B->data + i * skip,
                                     M, mod);

   pmfvec_clear (C);
   pmfvec_clear (B);
   pmfvec_clear (A);
   
   return success;
}



/*
   Tests pmfvec_ifft_dc (if huge == 0) or pmfvec_ifft_huge (if huge == 1)
*/
int
test_pmfvec_ifft_dc_or_huge (int huge, int quick)
{
   int success = 1;
   int i;
   unsigned lgK, lgM, lgT;
   ulong z, n, t;
   int fwd;
   zn_mod_t mod;
   
   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = (huge ? 1 : 0); lgT < (huge ? lgK : 1) && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0;
        lgM < lgK + (quick ? 1 : 3) && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;

      for (t = 0; t < ZNP_MIN (2 * M / K, quick ? 2 : 1000) && success; t++)
      for (z = 1; z <= K && success; z++)
      for (fwd = 0; fwd < 2 && success; fwd++)
      for (n = 1 - fwd; n <= K - fwd && n <= z && success; n++)
      {
         zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
         success = success && testcase_pmfvec_ifft_dc_or_huge
                                           (lgK, lgM, lgT, n, fwd, z, t, mod);
         zn_mod_clear (mod);
      }
   }
   
   return success;
}


int
test_pmfvec_ifft_dc (int quick)
{
   return test_pmfvec_ifft_dc_or_huge (0, quick);
}


int
test_pmfvec_ifft_huge (int quick)
{
   return test_pmfvec_ifft_dc_or_huge (1, quick);
}



int
testcase_pmfvec_tpfft_dc (unsigned lgK, unsigned lgM, ulong n, ulong z,
                          ulong t, const zn_mod_t mod)
{
   int success = 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i, j;

   ptrdiff_t skip = M + 1;

   // ===================================
   // first check linearity, i.e. that ax + by gets mapped to the right thing

   {
      pmfvec_t X, Y, A, B, TX, TY;

      pmfvec_init (X, lgK, skip, lgM, mod);
      pmfvec_init (Y, lgK, skip, lgM, mod);
      pmfvec_init (A, lgK, skip, lgM, mod);
      pmfvec_init (B, lgK, skip, lgM, mod);
      pmfvec_init (TX, lgK, skip, lgM, mod);
      pmfvec_init (TY, lgK, skip, lgM, mod);

      // generate random X, Y, A, B
      for (i = 0; i < K; i++)
      {
         pmf_rand (X->data + i * skip, M, mod);
         pmf_rand (Y->data + i * skip, M, mod);
      }

      pmf_rand (A->data, M, mod);
      pmf_rand (B->data, M, mod);
      
      for (i = 1; i < K; i++)
      {
         pmf_set (A->data + i * A->skip, A->data, M);
         pmf_set (B->data + i * B->skip, B->data, M);
      }
      
      // transform X and Y (after throwing in random ignorable crap)
      pmfvec_set (TX, X);
      pmfvec_set (TY, Y);

      for (i = n; i < K; i++)
      {
         pmf_rand (TX->data + i * skip, M, mod);
         pmf_rand (TY->data + i * skip, M, mod);
      }

      pmfvec_tpfft_dc (TX, n, z, t);
      pmfvec_tpfft_dc (TY, n, z, t);

      // form linear combination of TX and TY
      pmfvec_mul (TX, TX, A, z, 0);
      pmfvec_mul (TY, TY, B, z, 0);
      for (i = 0; i < z; i++)
         pmf_add (TX->data + TX->skip * i, TY->data + TY->skip * i, M, mod);

      // form linear combination of X and Y
      pmfvec_mul (X, X, A, n, 0);
      pmfvec_mul (Y, Y, B, n, 0);
      for (i = 0; i < n; i++)
         pmf_add (X->data + X->skip * i, Y->data + Y->skip * i, M, mod);

      // transform linear combination of X and Y
      pmfvec_tpfft_dc (X, n, z, t);

      // compare results
      for (i = 0; i < z; i++)
         success = success && !pmf_cmp (X->data + X->skip * i,
                                        TX->data + TX->skip * i, M, mod);
      
      pmfvec_clear (X);
      pmfvec_clear (Y);
      pmfvec_clear (TX);
      pmfvec_clear (TY);
      pmfvec_clear (A);
      pmfvec_clear (B);
   }
   
   // ===================================
   // now check that the matrix of the transposed FFT is really the transpose
   // of the matrix of the ordinary FFT

   {
      pmfvec_t* X = (pmfvec_t*) malloc (z * sizeof(pmfvec_t));
      for (i = 0; i < z; i++)
         pmfvec_init (X[i], lgK, skip, lgM, mod);

      pmfvec_t* Y = (pmfvec_t*) malloc (n * sizeof(pmfvec_t));
      for (i = 0; i < n; i++)
         pmfvec_init (Y[i], lgK, skip, lgM, mod);
         
      // compute images of basis vectors under FFT
      for (i = 0; i < z; i++)
      for (j = 0; j < z; j++)
      {
         pmf_zero (X[i]->data + j * skip, M);
         X[i]->data[j * skip + 1] = (i == j);
      }
      
      for (i = 0; i < z; i++)
         pmfvec_fft (X[i], n, z, t);

      // compute images of basis vectors under transposed FFT
      for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
      {
         pmf_zero (Y[i]->data + j * skip, M);
         Y[i]->data[j * skip + 1] = (i == j);
      }
      
      for (i = 0; i < n; i++)
         pmfvec_tpfft (Y[i], n, z, t);
      
      // check that they are transposes of each other
      for (i = 0; i < z; i++)
      for (j = 0; j < n; j++)
         success = success && !pmf_cmp (X[i]->data + j * skip,
                                        Y[j]->data + i * skip, M, mod);

      for (i = 0; i < z; i++)
         pmfvec_clear (X[i]);
      for (i = 0; i < n; i++)
         pmfvec_clear (Y[i]);
      free (Y);
      free (X);
   }

   return success;
}



int
testcase_pmfvec_tpfft_huge (unsigned lgK, unsigned lgM, unsigned lgT, ulong n,
                            ulong z, ulong t, const zn_mod_t mod)
{
   pmfvec_t A, B;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i;

   ptrdiff_t skip = M + 1;

   pmfvec_init (A, lgK, skip, lgM, mod);
   pmfvec_init (B, lgK, skip, lgM, mod);

   // create random input
   for (i = 0; i < K; i++)
      pmf_rand (A->data + i * skip, M, mod);
   
   // make a copy
   pmfvec_set (B, A);
   
   // put random crap in B to check that it's ignored
   for (i = n; i < K; i++)
      pmf_rand (B->data + i * skip, M, mod);
   
   // run transposed FFTs using huge and dc algorithms
   pmfvec_tpfft_dc (A, n, z, t);
   pmfvec_tpfft_huge (B, lgT, n, z, t);

   // compare results
   int success = 1;
   for (i = 0; i < z; i++)
      success = success && !pmf_cmp (B->data + i * skip, A->data + i * skip,
                                     M, mod);
   
   pmfvec_clear (B);
   pmfvec_clear (A);

   return success;
}



/*
   Tests pmfvec_tpfft_dc (if huge == 0) or pmfvec_tpfft_huge (if huge == 1)
*/
int
test_pmfvec_tpfft_dc_or_huge (int huge, int quick)
{
   int success = 1;
   int i;
   unsigned lgK, lgM, lgT;
   ulong z, n, t;
   zn_mod_t mod;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = (huge ? 1 : 0); lgT < (huge ? lgK : 1) && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0;
        lgM < lgK + (quick ? 1 : 3) && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;

      for (t = 0; t < ZNP_MIN (2 * M / K, quick ? 2 : 1000) && success; t++)
      for (n = 1; n <= K && success; n++)
      for (z = 1; z <= K && success; z++)
      {
         zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
         success = success && (huge
                ? testcase_pmfvec_tpfft_huge (lgK, lgM, lgT, n, z, t, mod)
                : testcase_pmfvec_tpfft_dc (lgK, lgM, n, z, t, mod));
         zn_mod_clear (mod);
      }
   }
   
   return success;
}


int
test_pmfvec_tpfft_dc (int quick)
{
   return test_pmfvec_tpfft_dc_or_huge (0, quick);
}


int
test_pmfvec_tpfft_huge (int quick)
{
   return test_pmfvec_tpfft_dc_or_huge (1, quick);
}



int
testcase_pmfvec_tpifft_dc (unsigned lgK, unsigned lgM, ulong n, int fwd,
                           ulong z, ulong t, const zn_mod_t mod)
{
   int success = 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i, j;

   ptrdiff_t skip = M + 1;

   // ===================================
   // first check linearity, i.e. that ax + by gets mapped to the right thing

   {
      pmfvec_t X, Y, A, B, TX, TY;

      pmfvec_init (X, lgK, skip, lgM, mod);
      pmfvec_init (Y, lgK, skip, lgM, mod);
      pmfvec_init (A, lgK, skip, lgM, mod);
      pmfvec_init (B, lgK, skip, lgM, mod);
      pmfvec_init (TX, lgK, skip, lgM, mod);
      pmfvec_init (TY, lgK, skip, lgM, mod);

      // generate random X, Y, A, B
      for (i = 0; i < K; i++)
      {
         pmf_rand (X->data + i * skip, M, mod);
         pmf_rand (Y->data + i * skip, M, mod);
      }

      pmf_rand (A->data, M, mod);
      pmf_rand (B->data, M, mod);

      for (i = 1; i < K; i++)
      {
         pmf_set (A->data + i * A->skip, A->data, M);
         pmf_set (B->data + i * B->skip, B->data, M);
      }
      
      // transform X and Y (after throwing in random ignorable crap)
      pmfvec_set (TX, X);
      pmfvec_set (TY, Y);
      
      for (i = n + fwd; i < K; i++)
      {
         pmf_rand (TX->data + i * skip, M, mod);
         pmf_rand (TY->data + i * skip, M, mod);
      }

      pmfvec_tpifft_dc (TX, n, fwd, z, t);
      pmfvec_tpifft_dc (TY, n, fwd, z, t);

      // form linear combination of TX and TY
      pmfvec_mul (TX, TX, A, z, 0);
      pmfvec_mul (TY, TY, B, z, 0);
      for (i = 0; i < z; i++)
         pmf_add (TX->data + TX->skip * i, TY->data + TY->skip * i, M, mod);

      // form linear combination of X and Y
      pmfvec_mul (X, X, A, n + fwd, 0);
      pmfvec_mul (Y, Y, B, n + fwd, 0);
      for (i = 0; i < n + fwd; i++)
         pmf_add (X->data + X->skip * i, Y->data + Y->skip * i, M, mod);

      // transform linear combination of X and Y
      pmfvec_tpifft_dc (X, n, fwd, z, t);

      // compare results
      for (i = 0; i < z; i++)
         success = success && !pmf_cmp (X->data + X->skip * i,
                                        TX->data + TX->skip * i, M, mod);
      
      pmfvec_clear (X);
      pmfvec_clear (Y);
      pmfvec_clear (TX);
      pmfvec_clear (TY);
      pmfvec_clear (A);
      pmfvec_clear (B);
   }

   // ===================================
   // now check that the matrix of the transposed IFFT is really the transpose
   // of the matrix of the ordinary IFFT

   {
      pmfvec_t* X = (pmfvec_t*) malloc (z * sizeof (pmfvec_t));
      for (i = 0; i < z; i++)
         pmfvec_init (X[i], lgK, skip, lgM, mod);

      pmfvec_t* Y = (pmfvec_t*) malloc ((n + fwd) * sizeof (pmfvec_t));
      for (i = 0; i < n + fwd; i++)
         pmfvec_init (Y[i], lgK, skip, lgM, mod);
         
      // compute images of basis vectors under FFT

      for (i = 0; i < z; i++)
      for (j = 0; j < z; j++)
      {
         pmf_zero (X[i]->data + j * skip, M);
         X[i]->data[j * skip + 1] = (i == j);
      }
      
      for (i = 0; i < z; i++)
         pmfvec_ifft (X[i], n, fwd, z, t);

      // compute images of basis vectors under transposed FFT

      for (i = 0; i < n + fwd; i++)
      for (j = 0; j < n + fwd; j++)
      {
         pmf_zero (Y[i]->data + j * skip, M);
         Y[i]->data[j * skip + 1] = (i == j);
      }
      
      for (i = 0; i < n + fwd; i++)
         pmfvec_tpifft (Y[i], n, fwd, z, t);
      
      // check that they are transposes of each other

      for (i = 0; i < z; i++)
      for (j = 0; j < n + fwd; j++)
         success = success && !pmf_cmp (X[i]->data + j * skip,
                                        Y[j]->data + i * skip, M, mod);

      for (i = 0; i < z; i++)
         pmfvec_clear (X[i]);
      for (i = 0; i < n + fwd; i++)
         pmfvec_clear (Y[i]);
      free (Y);
      free (X);
   }

   return success;
}


int
testcase_pmfvec_tpifft_huge (unsigned lgK, unsigned lgM, unsigned lgT, ulong n,
                             int fwd, ulong z, ulong t, const zn_mod_t mod)
{
   pmfvec_t A, B;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ulong i;

   ptrdiff_t skip = M + 1;

   pmfvec_init (A, lgK, skip, lgM, mod);
   pmfvec_init (B, lgK, skip, lgM, mod);

   // create random input
   for (i = 0; i < K; i++)
      pmf_rand (A->data + i * skip, M, mod);
   
   // make a copy
   pmfvec_set (B, A);
   
   // put random crap in B to check that it's ignored
   for (i = n + fwd; i < K; i++)
      pmf_rand (B->data + i * skip, M, mod);
   
   // run transposed IFFTs using huge and dc algorithms
   pmfvec_tpifft_dc (A, n, fwd, z, t);
   pmfvec_tpifft_huge (B, lgT, n, fwd, z, t);
   
   // compare results
   int success = 1;
   for (i = 0; i < z; i++)
      success = success && !pmf_cmp (B->data + i * skip, A->data + i * skip,
                                     M, mod);
   
   pmfvec_clear (B);
   pmfvec_clear (A);

   return success;
}



/*
   Tests pmfvec_tpifft_dc (if huge == 0) or pmfvec_tpifft_huge (if huge == 1)
*/
int
test_pmfvec_tpifft_dc_or_huge (int huge, int quick)
{
   int success = 1;
   int i;
   unsigned lgK, lgM, lgT;
   ulong z, n, t;
   int fwd;
   zn_mod_t mod;

   for (i = 0; i < num_test_bitsizes && success; i++)
   for (lgK = 0; lgK < 5 && success; lgK++)
   for (lgT = (huge ? 1 : 0); lgT < (huge ? lgK : 1) && success; lgT++)
   for (lgM = lgK ? (lgK - 1) : 0;
        lgM < lgK + (quick ? 1 : 3) && success; lgM++)
   {
      ulong K = 1UL << lgK;
      ulong M = 1UL << lgM;

      for (t = 0; t < ZNP_MIN (2 * M / K, quick ? 2 : 1000) && success; t++)
      for (z = 1; z <= K && success; z++)
      for (fwd = 0; fwd < 2 && success; fwd++)
      for (n = 1 - fwd; n <= K - fwd && n <= z && success; n++)
      {
         zn_mod_init (mod, random_modulus (test_bitsizes[i], 1));
         success = success && (huge
               ? testcase_pmfvec_tpifft_huge (lgK, lgM, lgT, n, fwd, z, t, mod)
               : testcase_pmfvec_tpifft_dc (lgK, lgM, n, fwd, z, t, mod));
         zn_mod_clear (mod);
      }
   }
   
   return success;
}


int
test_pmfvec_tpifft_dc (int quick)
{
   return test_pmfvec_tpifft_dc_or_huge (0, quick);
}


int
test_pmfvec_tpifft_huge (int quick)
{
   return test_pmfvec_tpifft_dc_or_huge (1, quick);
}


// end of file ****************************************************************
