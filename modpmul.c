/******************************************************************************

 DEVELOPMENT CODE

 modpmul.c
 Convolutions in (Z/pZ)[x], where p is a word-sized "FFT prime".

 Copyright (C) 2006, David Harvey

 Throughout the documentation in this file, R denotes 2^FLINT_BITS_PER_LIMB.

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "profiler.h"
#include "modpmul.h"
#include "Z.h"
#include "mpn_extras.h"


/******************************************************************************

 REDC precomputations.

******************************************************************************/

/*
Initialises redc_info for the given modulus.

Modulus must be prime and must fit into FLINT_REDC_BITS.

todo: maybe add a flag to skip initialisation of R, RR if they are not needed.
*/
void redc_precomp_init(redc_precomp_t* redc_info, unsigned long p)
{
   redc_info->p = p;

   // First compute -p^(-1) mod R. Do this using Z_invert_long, but we can't
   // quite fit R in a word :-), so first we do it mod R/2, and then get the
   // last bit by trial and error.
   unsigned long half_R = 1L << (FLINT_BITS_PER_LIMB - 1);
   unsigned long guess = Z_invert_long(-p, half_R);
   redc_info->pinv = (guess * p == -1) ? guess : (guess + half_R);
   
   // compute R mod p and R^2 mod p
   unsigned long dummy;
   udiv_qrnnd(dummy, redc_info->R, 1, 0, p);
   redc_info->RR = mod_mul(redc_info->R, redc_info->R, p);
}



/******************************************************************************

Inplace matrix transpose operations, to assist in converting a very long
FFT to a series of smaller FFTs, so that the smaller ones can be done in cache.

The transpose routines recursively divide the matrices into quadrants, and
switch to a simple basecase transpose loop when all the relevant data sits
nicely on a few cache lines in L1.

The "square" versions operate on a portion of a S x S matrix, where
S = 2^log_size, stored in row-major order. Matrix transpose means what you
think it does: the cell with binary index abcdeABCDE gets moved to ABCDEabcde.

The "rectangle" versions operate on a portion of a S x 2S matrix, i.e. S rows
of length 2S. It's not possible to efficiently transpose such a matrix inplace.
Instead we use a "pseudo-transpose", which transposes the left half and the
right half independently, i.e. the cell with binary index abcdeXABCDE gets
moved to the cell ABCDEXabcde. This *can* be done efficiently inplace, and
turns out to work fine for the cache-friendly FFT.

******************************************************************************/

/*
Sets A = transpose(B) and B = transpose(A), where A and B are the submatrices
with top-left corners (x1, y1) and (x2, y2) and side length N.
*/
void transpose_square_swap(unsigned long* data, unsigned long log_size,
                           unsigned long x1, unsigned long y1,
                           unsigned long x2, unsigned long y2,
                           unsigned long N)
{
   if (N <= FLINT_MODP_FFT_TRANSPOSE_THRESHOLD)
   {
      // basecase transpose
      for (unsigned long j = 0; j < N; j++)
         for (unsigned long i = 0; i < N; i++)
         {
            unsigned long index1 = x1 + i + ((y1 + j) << log_size);
            unsigned long index2 = x2 + j + ((y2 + i) << log_size);
            unsigned long temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
         }
   }
   else
   {
      // recurse into four quadrants
      unsigned long half = N >> 1;
      transpose_square_swap(data, log_size, x1, y1, x2, y2, half);
      transpose_square_swap(data, log_size, x1 + half, y1, x2, y2 + half,
                            half);
      transpose_square_swap(data, log_size, x1, y1 + half, x2 + half, y2,
                            half);
      transpose_square_swap(data, log_size, x1 + half, y1 + half,
                            x2 + half, y2 + half, half);
   }
}

/*
Sets A = transpose(A), where A is the submatrix with top-left corner (x1, y1)
and side length N
*/
void transpose_square_inplace(unsigned long* data, unsigned long log_size,
                              unsigned long x, unsigned long y,
                              unsigned long N)
{
   if (N <= FLINT_MODP_FFT_TRANSPOSE_THRESHOLD)
   {
      // basecase transpose
      for (unsigned long j = 0; j < N; j++)
         for (unsigned long i = j + 1; i < N; i++)
         {
            unsigned long index1 = x + i + ((y + j) << log_size);
            unsigned long index2 = x + j + ((y + i) << log_size);
            unsigned long temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
         }
   }
   else
   {
      // recurse into four quadrants
      unsigned long half = N >> 1;
      transpose_square_inplace(data, log_size, x, y, half);
      transpose_square_inplace(data, log_size, x + half, y + half, half);
      transpose_square_swap(data, log_size, x, y + half, x + half, y, half);
   }
}


/*
Sets Ai = transpose(Bi), Bi = transpose(Ai), for i = 1, 2, where A1 and B1
are the submatrices with top-left corners (x1, y1) and (x2, y2) and side length
N, and where A2 and B2 are the corresponding submatrices on the right hand side
(i.e. with co-ordinates (x1 + half, y1) and (x2 + half, y2)).
*/
void transpose_rectangle_swap(unsigned long* data, unsigned long log_size,
                              unsigned long x1, unsigned long y1,
                              unsigned long x2, unsigned long y2,
                              unsigned long N)
{
   if (N <= FLINT_MODP_FFT_TRANSPOSE_THRESHOLD)
   {
      // basecase transpose for left half
      for (unsigned long j = 0; j < N; j++)
         for (unsigned long i = 0; i < N; i++)
         {
            unsigned long index1 = x1 + i + ((y1 + j) << (log_size + 1));
            unsigned long index2 = x2 + j + ((y2 + i) << (log_size + 1));
            unsigned long temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
         }
      // basecase transpose for right half
      for (unsigned long j = 0; j < N; j++)
         for (unsigned long i = 0; i < N; i++)
         {
            unsigned long index1 = x1 + i + ((y1 + j) << (log_size + 1))
                                   + (1 << log_size);
            unsigned long index2 = x2 + j + ((y2 + i) << (log_size + 1))
                                   + (1 << log_size);
            unsigned long temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
         }
   }
   else
   {
      // recurse into four quadrants
      unsigned long half = N >> 1;
      transpose_rectangle_swap(data, log_size, x1, y1, x2, y2, half);
      transpose_rectangle_swap(data, log_size, x1 + half, y1, x2, y2 + half,
                               half);
      transpose_rectangle_swap(data, log_size, x1, y1 + half, x2 + half, y2,
                               half);
      transpose_rectangle_swap(data, log_size, x1 + half, y1 + half,
                               x2 + half, y2 + half, half);
   }
}


/*
Sets Ai = transpose(Ai) for i = 1, 2, where A1 is the submatrix with top-left
corner (x, y) and side length N, and where A2 is the corresponding submatrix
on the right hand side (i.e. with co-ordinates (x + half, y)).
*/
void transpose_rectangle_inplace(unsigned long* data, unsigned long log_size,
                                 unsigned long x, unsigned long y,
                                 unsigned long N)
{
   if (N <= FLINT_MODP_FFT_TRANSPOSE_THRESHOLD)
   {
      // basecase transpose for left half
      for (unsigned long j = 0; j < N; j++)
         for (unsigned long i = j + 1; i < N; i++)
         {
            unsigned long index1 = x + i + ((y + j) << (log_size + 1));
            unsigned long index2 = x + j + ((y + i) << (log_size + 1));
            unsigned long temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
         }
      // basecase transpose for right half
      for (unsigned long j = 0; j < N; j++)
         for (unsigned long i = j + 1; i < N; i++)
         {
            unsigned long index1 = x + i + ((y + j) << (log_size + 1))
                                   + (1 << log_size);
            unsigned long index2 = x + j + ((y + i) << (log_size + 1))
                                   + (1 << log_size);
            unsigned long temp = data[index1];
            data[index1] = data[index2];
            data[index2] = temp;
         }
   }
   else
   {
      // recurse into four quadrants
      unsigned long half = N >> 1;
      transpose_rectangle_inplace(data, log_size, x, y, half);
      transpose_rectangle_inplace(data, log_size, x + half, y + half, half);
      transpose_rectangle_swap(data, log_size, x, y + half, x + half, y, half);
   }
}


/******************************************************************************

 FFT precomputations

******************************************************************************/

/*
Prepares a modp_fft_precomp_t struct with information for performing an
FFT/IFFT for a given combination of transform length, prime modulus, and
primitive root of unity.

INPUTS:
- redc_info is the output of redc_precomp_init for the given prime.
- length of transform is N = 2^logN
- w is a primitive Nth root of unit mod p, in REDC format.
- basecase_threshold is a parameter indicating when to switch over to the
  basecase iterative algorithm. Best value should be
  FLINT_MODP_FFT_BASECASE_THRESHOLD, but different values can be used for
  debugging (e.g. a very large value will force basecase for all transform
  lengths; a small value will give the FFT factoring code a workout.)

OUTPUT:
Precomputed info is stored in info. (This is a recursive data structure;
it should be cleaned out by modp_fft_precomp_clear when not needed any more.)
*/
void modp_fft_precomp_init(modp_fft_precomp_t* info, redc_precomp_t* redc_info,
                           unsigned long logN, unsigned long w,
                           unsigned long basecase_threshold)
{
   info->redc_info = redc_info;
   info->logN = logN;
   unsigned long N = 1 << logN;
   info->N = N;
   info->w = w;

   if (logN <= basecase_threshold)
   {
      // basecase transform
      info->basecase = 1;
      info->roots = (unsigned long*) limb_alloc(N/2, 0);

      // compute roots of unity w, w^2, w^3, ..., w^(N/2 - 1)
      info->roots[1] = w;
      for (unsigned long i = 2; i < N/2; i++)
         info->roots[i] = mul_redc2(info->roots[i-1], w, redc_info);
   }
   else
   {
      // Non-basecase case.
      info->basecase = 0;

      // The plan here is to split logN = logJ + logK (i.e. N = J*K),
      // where logJ and logK are approximately equal. Then in the FFT we
      // will do K transforms of length J followed by J transforms of length K.
      // (or the other way for the IFFT).

      // To make the transposes easier, we ensure that either logK == logJ
      // (square case) or logK == logJ + 1 (rectangular case).
      unsigned long logJ = logN >> 1;
      unsigned long logK = logN - logJ;

      // child1 corresponds to the subtransforms of length J.
      // child2 corresponds to the subtransforms of length K.
      // todo: use smart allocator instead of malloc here:
      info->child1 = (modp_fft_precomp_t*) malloc(sizeof(modp_fft_precomp_t));
      info->child2 = (modp_fft_precomp_t*) malloc(sizeof(modp_fft_precomp_t));

      // Compute inverse of w mod p. (Since w is in REDC format we need to
      // multiply by a correction factor.)
      unsigned long winv = Z_invert_long(w, redc_info->p);
      winv = convert_to_redc2(convert_to_redc2(winv, redc_info), redc_info);

      // Compute w, w^2, w^4, ..., w^(J/2)
      // and w^-1, w^-2, w^-4, ... w^(-J/2),
      // store them in roots array,
      // also compute w_to_J = w^J and w_to_K = w^K
      info->roots = (unsigned long*) limb_alloc(2 * logJ, 0);
      unsigned long w_to_J = w;
      unsigned long winv_to_J = winv;
      for (unsigned long i = 0; i < logJ; i++)
      {
         info->roots[logJ - 1 - i] = w_to_J;
         info->roots[2*logJ - 1 - i] = winv_to_J;
         w_to_J = mul_redc2(w_to_J, w_to_J, redc_info);
         winv_to_J = mul_redc2(winv_to_J, winv_to_J, redc_info);
      }
      unsigned long w_to_K;
      if (logK == logJ)
         w_to_K = w_to_J;
      else
         w_to_K = mul_redc2(w_to_J, w_to_J, redc_info);

      // Recursively build precomputed information for the subtransforms
      // of length J and K
      modp_fft_precomp_init(info->child1, redc_info, logJ,
                            w_to_K, basecase_threshold);
      modp_fft_precomp_init(info->child2, redc_info, logK,
                            w_to_J, basecase_threshold);
   }
}

/*
Recursively deallocates a modp_fft_precomp_t struct that was allocated by
modp_fft_precomp_init.
 */
void modp_fft_precomp_clear(modp_fft_precomp_t* info)
{
   if (info->basecase)
      limb_release();      // release info->roots
   else
   {
      modp_fft_precomp_clear(info->child2);
      modp_fft_precomp_clear(info->child1);
      limb_release();      // release info->roots
      free(info->child2);
      free(info->child1);
   }
}


/******************************************************************************

Basecase FFT routines, which work entirely in cache with precomputed roots
of unity.

******************************************************************************/

/*
Base case iterative FFT.

info must have basecase == 1; it contains precomputed information describing
the transform, including appropriate roots of unity.

data is the input data, an array of length N = 2^logN. All input coefficients
must be in [0, p).

If first == 1, the second half of the input data is assumed to be zero, and
appropriate optimisations are made. (todo: this is not implemented yet)

Output is in-place, in bit-reversed order. So output[bitrev(j)] will hold
\sum_{k=0}^{N-1} input[k] w^{jk}.

Ideally the transform (including the roots of unity) should fit in L1.

*/
void modp_fft_basecase(modp_fft_precomp_t* info, unsigned long* data,
                       int first)
{
   unsigned long p = info->redc_info->p;
   unsigned long pinv = info->redc_info->pinv;
   unsigned long* roots = info->roots;
   unsigned long length = info->N;
   unsigned long half = length >> 1;
   long start;
   unsigned long layer, root_skip = 1, i, j;
   unsigned long sum, diff;

   // todo: special case for first == 1

   // We do the layers in two groups, to ensure that the inner loop is
   // always fairly long, so that we don't waste time keeping track of
   // two counter variables.

   // First group:
   for (layer = info->logN; layer > info->logN / 2;
        layer--, half >>= 1, root_skip <<= 1)
   {
      // loop through blocks
      for (start = length - 2*half; start >= 0; start -= 2*half)
      {
         // loop through roots of unity
         
         // the non-trivial roots
         for (i = half - 1, j = length/2 - root_skip;
              i > 0; i--, j -= root_skip)
         {
            // (a, b) -> (a + b, w*(a - b))
            sum = data[start + i] + data[start + half + i];
            if (sum >= p)
               sum -= p;
            diff = p + data[start + i] - data[start + half + i];
            // don't need to normalise diff because mul_redc doesn't care
            data[start + i] = sum;
            data[start + half + i] = mul_redc(diff, roots[j], p, pinv);
         }

         // special case for root = 1
         // (a, b) -> (a + b, a - b)
         sum = data[start] + data[start + half];
         if (sum >= p)
            sum -= p;
         diff = data[start] - data[start + half];
         if ((long)diff < 0)
            diff += p;
         data[start] = sum;
         data[start + half] = diff;
      }
   }

   // Second group:
   for (; layer > 1; layer--, half >>= 1, root_skip <<= 1)
   {
      // loop through roots of unity

      // special case for root = 1
      for (start = length - 2*half; start >= 0; start -= 2*half)
      {
         // (a, b) -> (a + b, a - b)
         sum = data[start] + data[start + half];
         if (sum >= p)
            sum -= p;
         diff = data[start] - data[start + half];
         if ((long)diff < 0)
            diff += p;
         data[start] = sum;
         data[start + half] = diff;
      }

      // now the non-trivial roots
      for (i = half - 1, j = length/2 - root_skip; i > 0; i--, j -= root_skip)
      {
         unsigned long root = roots[j];

         // loop through blocks
         for (start = length - 2*half; start >= 0; start -= 2*half)
         {
            // (a, b) -> (a + b, w*(a - b))
            sum = data[start + i] + data[start + half + i];
            if (sum >= p)
               sum -= p;
            diff = p + data[start + i] - data[start + half + i];
            // don't need to normalise diff because mul_redc doesn't care
            data[start + i] = sum;
            data[start + half + i] = mul_redc(diff, root, p, pinv);
         }
      }
   }

   // Finally a special case for the last layer (when half == 1)
   for (start = length - 2; start >= 0; start -= 2)
   {
      // (a, b) -> (a + b, a - b)
      sum = data[start] + data[start + 1];
      if (sum >= p)
         sum -= p;
      diff = data[start] - data[start + 1];
      if ((long)diff < 0)
         diff += p;
      data[start] = sum;
      data[start + 1] = diff;
   }
}


/*
Basecase inverse FFT.

info should be the same struct that was passed to the corresponding FFT.

Input data is in bit-reversed order, output is in correct order.

Output will be multiplied by a factor of N = 2^logN from the true inverse
transform.

The input does NOT have to be reduced mod p; it is sufficient for each
coefficient to be in the range [0, 2p). SIMILARLY, the output is only
guaranteed to be in the range [0, 2p), not necessarily reduced into [0, p).

*/
void modp_ifft_basecase(modp_fft_precomp_t* info, unsigned long* data)
{
   unsigned long p = info->redc_info->p;
   unsigned long pinv = info->redc_info->pinv;
   unsigned long* roots = info->roots;
   unsigned long length = info->N;
   long start;
   unsigned long layer, i, j;
   unsigned long sum, diff;

   // We do the layers in two groups, to ensure that the inner loop is
   // always fairly long, so that we don't waste time keeping track of
   // two counter variables.

   // special case for half == 1
   for (start = length - 2; start >= 0; start -= 2)
   {
      // (a + b, a - b) -> (2a, 2b)
      diff = data[start + 1];
      if (diff >= p)
         diff -= p;
      sum = data[start];
      if (sum >= p)
         sum -= p;
      // don't need to normalise outputs
      data[start] = sum + diff;
      data[start + 1] = p + sum - diff;
   }

   unsigned long half = 2;
   unsigned long root_skip = length >> 2;
   
   // First group:
   for (layer = 1; layer < info->logN/2; layer++, half <<= 1, root_skip >>= 1)
   {
      // loop through roots of unity

      // special case for root = 1
      for (start = length - 2*half; start >= 0; start -= 2*half)
      {
         // (a + b, a - b) -> (2a, 2b)
         diff = data[start + half];
         if (diff >= p)
            diff -= p;
         sum = data[start];
         if (sum >= p)
            sum -= p;
         // don't need to normalise outputs
         data[start] = sum + diff;
         data[start + half] = p + sum - diff;
      }

      // now the non-trivial roots
      for (i = half - 1, j = root_skip; i > 0; i--, j += root_skip)
      {
         // (note: here roots[j] = -1/w, not 1/w)
         unsigned long root = roots[j];

         // loop through blocks
         for (start = length - 2*half; start >= 0; start -= 2*half)
         {
            // (a + b, w*(a - b)) -> (2a, 2b)
            // don't need to normalise data[start + half + i] because
            // mul_redc doesn't care
            diff = mul_redc(data[start + half + i], root, p, pinv);
            sum = data[start + i];
            if (sum >= p)
               sum -= p;
            // don't need to normalise outputs
            data[start + i] = p + sum - diff;
            data[start + half + i] = sum + diff;
         }
      }
   }

   // Second group:
   for (; layer < info->logN; layer++, half <<= 1, root_skip >>= 1)
   {
      // loop through blocks
      for (start = length - 2*half; start >= 0; start -= 2*half)
      {
         // loop through roots of unity
         
         // the non-trivial roots
         for (i = half - 1, j = root_skip; i > 0; i--, j += root_skip)
         {
            // (a + b, w*(a - b)) -> (2a, 2b)
            // don't need to normalise data[start + half + i] because
            // mul_redc doesn't care
            // (note: here roots[j] = -1/w, not 1/w)
            diff = mul_redc(data[start + half + i], roots[j], p, pinv);
            sum = data[start + i];
            if (sum >= p)
               sum -= p;
            // don't need to normalise outputs
            data[start + i] = p + sum - diff;
            data[start + half + i] = sum + diff;
         }

         // special case for root = 1
         // (a + b, a - b) -> (2a, 2b)
         diff = data[start + half + i];
         if (diff >= p)
            diff -= p;
         sum = data[start + i];
         if (sum >= p)
            sum -= p;
         // don't need to normalise outputs
         data[start + i] = sum + diff;
         data[start + half + i] = p + sum - diff;
      }
   }
}


/******************************************************************************

Main FFT/IFFT routines.

These routines implement some variant of Bailey's algorithm for cache-friendly
FFTs. Essentially we start with the FFT:

 b_{j_1 + Jk_1} = \sum_{j_2, k_2} w^{(j_1 + Jk_1)(k_2 + Kj_2)} a_{k_2 + Kj_2}

and factor it as

 b_{j_1 + Jk_1} = \sum_{k_2} (w^J)^{k_1 k_2}
                     ( w^{j_1 k_2}  \sum_{j_2} (w^K)^(j_1 j_2) a_{k_2 + Kj_2} )

******************************************************************************/

/*
This code really belongs in modp_fft, but it's in a separate function because
it's most efficient to implement recursively. It basically applies the twiddle
factors w^{j_1 k_2} and then performs the outer FFTs in the above sums.

data is an array of 2^depth blocks of length K, where 0 <= depth <= logJ.
Suppose the first block corresponds to some j = j_1. Then u should be
the root w^j (in REDC). This routine recurses down to all the blocks of
length K, and when it hits the bottom, it applies the twiddles and calls
the subtransforms of length K.

info should be the precomputed information corresponding to the transform
of length N (i.e. the one modp_fft was using).

Note that the j_1's occur in bit-reversed order, because the inner transform
was bit-reversing. That's why this routine is recursive; it avoids explicit
bit-reversal.
 */
void modp_fft_round2(modp_fft_precomp_t* info, unsigned long* data,
                     unsigned long depth, unsigned long u)
{
   if (depth == 0)
   {
      // Apply twiddle factors (except in the case u == R, which means
      // we're on the first block where the twiddles would be 1, 1, 1, ... 1).
      if (u != info->redc_info->R)
      {
         unsigned long p = info->redc_info->p;
         unsigned long pinv = info->redc_info->pinv;

         // i.e. set data[k] = u^k * data[k] for k = 0, ... K-1.
         data[1] = mul_redc(data[1], u, p, pinv);
         unsigned long u_to_k = u;
         for (unsigned long k = 2; k < info->child2->N; k++)
         {
            u_to_k = mul_redc(u_to_k, u, p, pinv);
            data[k] = mul_redc(data[k], u_to_k, p, pinv);
         }
      }

      // Perform the outer transform on this block
      modp_fft(info->child2, data, 0);
   }
   else
   {
      // Recurse into first 2^(depth-1) blocks
      modp_fft_round2(info, data, depth - 1, u);
      // Recurse into second 2^(depth-1) blocks
      u = mul_redc2(u, info->roots[depth-1], info->redc_info);
      modp_fft_round2(info, data + (info->child2->N << (depth - 1)),
                      depth - 1, u);
   }
}

/*
This is essentially the same as modp_fft_round2 but in reverse.
 */
void modp_ifft_round1(modp_fft_precomp_t* info, unsigned long* data,
                      unsigned long depth, unsigned long u)
{
   if (depth == 0)
   {
      modp_ifft(info->child2, data);
      
      if (u != info->redc_info->R)
      {
         unsigned long p = info->redc_info->p;
         unsigned long pinv = info->redc_info->pinv;

         data[1] = mul_redc(data[1], u, p, pinv);
         unsigned long u_to_k = u;
         for (unsigned long k = 2; k < info->child2->N; k++)
         {
            u_to_k = mul_redc(u_to_k, u, p, pinv);
            data[k] = mul_redc(data[k], u_to_k, p, pinv);
         }
      }
   }
   else
   {
      modp_ifft_round1(info, data, depth - 1, u);
      u = mul_redc2(info->roots[info->child1->logN + depth - 1],
                    u, info->redc_info);
      modp_ifft_round1(info, data + (info->child2->N << (depth - 1)),
                       depth - 1, u);
   }
}

/*
Main FFT routine.

info contains precomputed information describing the transform.

data is an array of length N = 2^logN.

If first == 1, the second half of the input is assumed to be zero
(todo: not implemented yet).

Output is inplace, in bit-reversed order.
*/
void modp_fft(modp_fft_precomp_t* info, unsigned long* data, int first)
{
   // todo: implement "first" flag

   if (info->basecase)
   {
      // go straight to basecase if we fit in cache
      modp_fft_basecase(info, data, 0);
      return;
   }

   unsigned long logJ = info->child1->logN;
   unsigned long J = info->child1->N;
   unsigned long K = info->child2->N;
   unsigned long N = info->N;

   // transpose data to make the spaced-out coefficients consecutive
   if (J == K)
      transpose_square_inplace(data, logJ, 0, 0, J);
   else
      transpose_rectangle_inplace(data, logJ, 0, 0, J);

   // recursively FFT K blocks of length J
   for (unsigned long k = 0; k < N; k += J)
      modp_fft(info->child1, data + k, 0);

   // transpose back so we can work on inner blocks
   if (J == K)
      transpose_square_inplace(data, logJ, 0, 0, J);
   else
      transpose_rectangle_inplace(data, logJ, 0, 0, J);

   // apply twiddles, and recursively FFT J blocks of length K
   modp_fft_round2(info, data, logJ, info->redc_info->R);
}


/*
Similar to modp_fft, but performs the inverse transform.

Input in bit-reversed order, output in correct order.

Output is multiplied by a factor of N = 2^logN.
*/
void modp_ifft(modp_fft_precomp_t* info, unsigned long* data)
{
   if (info->basecase)
   {
      // go straight to basecase if we fit in cache
      modp_ifft_basecase(info, data);
      return;
   }
   
   // Basically the plan is to do the FFT backwards!

   unsigned long logJ = info->child1->logN;
   unsigned long J = info->child1->N;
   unsigned long K = info->child2->N;
   unsigned long N = info->N;

   modp_ifft_round1(info, data, logJ, info->redc_info->R);

   if (J == K)
      transpose_square_inplace(data, logJ, 0, 0, J);
   else
      transpose_rectangle_inplace(data, logJ, 0, 0, J);

   for (unsigned long k = 0; k < N; k += J)
      modp_ifft(info->child1, data + k);

   if (J == K)
      transpose_square_inplace(data, logJ, 0, 0, J);
   else
      transpose_rectangle_inplace(data, logJ, 0, 0, J);
}


/******************************************************************************

Convolutions

todo: this code is far from optimised yet. Probably the normalisation
stuff should be left to the caller, since it would be more efficient to roll
it into the reductions, CRT etc. Also surely we want to do both FFTs
simultaneously so that the pointwise multiplications all happen locally.

******************************************************************************/

// data1 and data2 should be length 2^n arrays, coeffs in [0, p),
// w should be a 2^n-th root of unity mod p in REDC format
// convolution is written to data1
// data2 is destroyed
// it's a cyclic convolution mod x^(2^n) - 1
void modp_convolution(redc_precomp_t* redc_info, unsigned long n,
                      unsigned long w,
                      unsigned long* data1, unsigned long* data2)
{
   modp_fft_precomp_t info;

   modp_fft_precomp_init(&info, redc_info, n, w,
                         FLINT_MODP_FFT_BASECASE_THRESHOLD);

   modp_fft(&info, data1, 0);
   modp_fft(&info, data2, 0);

   unsigned long p = redc_info->p;
   unsigned long pinv = redc_info->pinv;
   unsigned long fudge = convert_to_redc2(1L << (FLINT_BITS_PER_LIMB - n),
                                          redc_info);

   for (unsigned long i = 0; i < (1 << n); i++)
      data1[i] = mul_redc(data1[i], data2[i], p, pinv);

   modp_ifft(&info, data1);

   // todo: it would be better to roll this in with the previous loop.
   // The problem is that IFFT output is not normalised into [0, p), so
   // we have to do a normalisation pass anyway.
   // We could fix this by adding a flag to modp_ifft which does an extra
   // normalisation inside the IFFT.
   for (unsigned long i = 0; i < (1 << n); i++)
      data1[i] = mul_redc(data1[i], fudge, p, pinv);

   modp_fft_precomp_clear(&info);
}


// end of file ****************************************************************
