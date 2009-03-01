/*
   mul_fft_dft.c:  multiplication by Schonhage FFT, with a few layers of
                   naive DFT to save memory
   
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


/*
   Returns length n bit reversal of x.
*/
#define bit_reverse \
    ZNP_bit_reverse
ulong bit_reverse(ulong x, unsigned n)
{
   ulong y = 0;
   unsigned i;
   for (i = 0; i < n; i++)
   {
      y <<= 1;
      y += x & 1;
      x >>= 1;
   }
   return y;
}



/*
   Let [a, b) = intersection of [0, len) and [start, start + M/2).
   This functions adds op[a, b) to res.
*/
#define merge_chunk_to_pmf \
    ZNP_merge_chunk_to_pmf
void merge_chunk_to_pmf(zn_pmf_t res, const ulong* op, size_t len,
                        size_t start, ulong M, const zn_mod_t mod)
{
   ZNP_ASSERT((M & 1) == 0);
   
   ulong r = (-res[0]) & (2*M - 1);
   
   size_t end = start + M/2;
   if (end > len)
      end = len;
   if (start >= end)
      // nothing to do
      return;
   
   op += start;
   ulong size = end - start;
   // now we need to handle op[0, size), and we are guaranteed size <= M/2.

   if (r < M)
   {
      if (size <= M - r)
         zn_array_add_inplace(res + 1 + r, op, size, mod);
      else
      {
         zn_array_add_inplace(res + 1 + r, op, M - r, mod);
         // negacyclic wraparound:
         zn_array_sub_inplace(res + 1, op + M - r, size - M + r, mod);
      }
   }
   else
   {
      r -= M;

      if (size <= M - r)
         zn_array_sub_inplace(res + 1 + r, op, size, mod);
      else
      {
         zn_array_sub_inplace(res + 1 + r, op, M - r, mod);
         // negacyclic wraparound:
         zn_array_add_inplace(res + 1, op + M - r, size - M + r, mod);
      }
   }
}


/*
   Adds op into res, starting at index _start_, and not writing beyond
   res + len.
   
   If op == NULL, does no operation.
*/
#define merge_chunk_from_pmf \
    ZNP_merge_chunk_from_pmf
void merge_chunk_from_pmf(ulong* res, size_t len, const zn_pmf_t op,
                          size_t start, ulong M, const zn_mod_t mod)
{
   if (op == NULL)
      return;

   size_t end = start + M;
   if (end > len)
      end = len;
   if (start >= end)
      // nothing to do
      return;
   
   res += start;
   ulong size = end - start;
   // now we need to write to res[0, size), and we are guaranteed size <= M.

   ulong r = op[0] & (2*M - 1);

   if (r < M)
   {
      if (size <= r)
         zn_array_sub_inplace(res, op + 1 + M - r, size, mod);
      else
      {
         zn_array_sub_inplace(res, op + 1 + M - r, r, mod);
         // negacyclic wraparound:
         zn_array_add_inplace(res + r, op + 1, size - r, mod);
      }
   }
   else
   {
      r -= M;
   
      if (size <= r)
         zn_array_add_inplace(res, op + 1 + M - r, size, mod);
      else
      {
         zn_array_add_inplace(res, op + 1 + M - r, r, mod);
         // negacyclic wraparound:
         zn_array_sub_inplace(res + r, op + 1, size - r, mod);
      }
   }
}


/* ============================================================================

     "virtual" zn_pmf_t's

============================================================================ */


/*
   The zn_virtual_pmf_t and zn_virtual_pmf_vec_t are similar to their
   non-virtual counterparts (zn_pmf_t and zn_pmf_vec_t), but the underlying
   representation is are optimised for the case where:
     * many of the coefficients in the vector are zero, and
     * many of the coefficients differ only by multiplication by a root of
       unity.
   
   They are used in the mul_fft_dft routine (below) to perform truncated
   inverse FFTs on vectors containing only a single nonzero entry. This is
   achieved in pretty much a constant (logarithmic?) number of coefficient
   operations, instead of K*lg(K) coefficient operations. (In the *forward*
   FFT case, it's easy to write down a "formula" for the FFT, so we don't even
   need a vector in which to do the computation; but in the inverse case I
   couldn't see a simple formula, so these structs in effect let us figure out
   a formula on the fly.)
*/


#define zn_virtual_pmf_vec_struct \
    ZNP_zn_virtual_pmf_vec_struct
struct zn_virtual_pmf_vec_struct;     // forward declaration


/*
   Virtual version of a zn_pmf_t.
   
   Each zn_virtual_pmf_t belongs to a "parent" zn_virtual_pmf_vec_t.
   
   The _index_ field is:
     * -1 if the value is zero
     * otherwise, an index into parent's _buf_ array, which is where the
       actual coefficient data is stored.
   
   The _bias_ field overrides the bias word in the underlying zn_pmf_t.
   (This lets different zn_virtual_pmf's share memory but still have
   different bias values.)
*/
#define zn_virtual_pmf_struct \
    ZNP_zn_virtual_pmf_struct
struct zn_virtual_pmf_struct
{
   struct zn_virtual_pmf_vec_struct* parent;
   int index;
   ulong bias;
};

#define zn_virtual_pmf_t \
    ZNP_zn_virtual_pmf_t
typedef struct zn_virtual_pmf_struct zn_virtual_pmf_t[1];


/*
   Virtual version of zn_pmf_vec_t.
   
   M, lgM, K, lgK, mod are as for zn_pmf_vec_t.
   
   _data_ is an array of K zn_virtual_pmf_t's (the coefficients in the vector).
   
   The underlying data is managed by three arrays of length max_buffers:
     * buf[i] points to a zn_pmf_t, or NULL if the i-th slot is not
       yet associated to any actual memory.
     * ref_count[i] is a reference count for slot #i, i.e. the number of
       zn_pmf_virtual_pmf_t's pointing at this slot.
     * external[i] is a flag indicating whether the memory belongs to
       someone else (i.e. it's not the responsiblity of the zn_virtual_pmf_vec
       to free it).
*/
struct zn_virtual_pmf_vec_struct
{
   ulong M;
   unsigned lgM;

   ulong K;
   unsigned lgK;
   
   const zn_mod_struct* mod;

   zn_virtual_pmf_t* data;

   unsigned max_buffers;
   ulong** buf;
   unsigned* ref_count;
   int* external;
};

#define zn_virtual_pmf_vec_t \
    ZNP_zn_virtual_pmf_vec_t
typedef struct zn_virtual_pmf_vec_struct zn_virtual_pmf_vec_t[1];



/* ----------------------------------------------------------------------------

   zn_virtual_pmf_t infrastructure

---------------------------------------------------------------------------- */

/*
   Initialises a zn_virtual_pmf_t to zero, with a given parent vector.
*/
#define zn_virtual_pmf_init \
    ZNP_zn_virtual_pmf_init
void zn_virtual_pmf_init(zn_virtual_pmf_t res, zn_virtual_pmf_vec_t parent)
{
   res->index = -1;
   res->parent = parent;
}


/*
   Initialises a zn_virtual_pmf_vec_t to length K, with all zero values.
   
   All slots are initially marked as empty.
*/
#define zn_virtual_pmf_vec_init \
    ZNP_zn_virtual_pmf_vec_init
void zn_virtual_pmf_vec_init(zn_virtual_pmf_vec_t vec, unsigned lgK,
                             unsigned lgM, const zn_mod_t mod)
{
   vec->mod = mod;

   vec->lgM = lgM;
   vec->M = 1UL << lgM;
   
   vec->lgK = lgK;
   vec->K = 1UL << lgK;
   
   vec->data = (zn_virtual_pmf_t*) malloc(vec->K * sizeof(zn_virtual_pmf_t));

   ulong i;
   for (i = 0; i < vec->K; i++)
      zn_virtual_pmf_init(vec->data[i], vec);

   vec->max_buffers = 2 * vec->K;     // should be safe

   vec->buf = (ulong**) malloc(sizeof(ulong*) * vec->max_buffers);
   vec->ref_count = (unsigned*) malloc(sizeof(unsigned) * vec->max_buffers);
   vec->external = (int*) malloc(sizeof(int) * vec->max_buffers);
   
   for (i = 0; i < vec->max_buffers; i++)
   {
      vec->buf[i] = NULL;
      vec->ref_count[i] = 0;
      vec->external[i] = 0;
   }
}


/*
   Sets all values to zero, and detaches externally-allocated memory from
   all slots.
   
   This does *not* free any memory owned by this vector; buffers already
   allocated will get re-used by subsequent operations.
*/
#define zn_virtual_pmf_vec_reset \
    ZNP_zn_virtual_pmf_vec_reset
void zn_virtual_pmf_vec_reset(zn_virtual_pmf_vec_t vec)
{
   ulong i;
   for (i = 0; i < vec->K; i++)
      vec->data[i]->index = -1;

   for (i = 0; i < vec->max_buffers; i++)
   {
      vec->ref_count[i] = 0;
      
      if (vec->external[i])
      {
         vec->buf[i] = NULL;
         vec->external[i] = 0;
      }
   }
}


/*
   Destroys the vector, and frees all memory owned by it.
*/
#define zn_virtual_pmf_vec_clear \
    ZNP_zn_virtual_pmf_vec_clear
void zn_virtual_pmf_vec_clear(zn_virtual_pmf_vec_t vec)
{
   zn_virtual_pmf_vec_reset(vec);

   ulong i;
   for (i = 0; i < vec->max_buffers; i++)
      if (vec->buf[i] && !vec->external[i])
         free(vec->buf[i]);

   free(vec->external);
   free(vec->buf);
   free(vec->ref_count);
   free(vec->data);
}


/*
   Finds a free slot (one not attached to any underlying zn_pmf_t yet),
   and returns its index.
*/
#define zn_virtual_pmf_vec_find_slot \
    ZNP_zn_virtual_pmf_vec_find_slot
unsigned zn_virtual_pmf_vec_find_slot(zn_virtual_pmf_vec_t vec)
{
   unsigned i;
   for (i = 0; i < vec->max_buffers; i++)
      if (!vec->buf[i])
         return i;
   
   // this should never happen; we always should have enough slots
   ZNP_ASSERT(0);
}


/*
   Finds a slot attached to an underlying zn_pmf_t which is not currently
   used by any other zn_virtual_pmf_t's, and returns its index.
   
   If there are no such slots, it allocated more space, attaches a slot to it,
   and returns its index.
   
   In both cases, the reference count of the returned slot will be 1.
*/
#define zn_virtual_pmf_vec_new_buf \
    ZNP_zn_virtual_pmf_vec_new_buf
unsigned zn_virtual_pmf_vec_new_buf(zn_virtual_pmf_vec_t vec)
{
   // first search for an already-allocated buffer that no-one else is using
   unsigned i;
   for (i = 0; i < vec->max_buffers; i++)
      if (vec->buf[i] && !vec->ref_count[i])
         break;
   
   if (i == vec->max_buffers)
   {
      // not found; need to allocate more space
      i = zn_virtual_pmf_vec_find_slot(vec);
      vec->buf[i] = (ulong*) malloc(sizeof(ulong) * (vec->M + 1));
      vec->external[i] = 0;
   }
   
   vec->ref_count[i] = 1;
   return i;
}


/*
   res := 0
*/
#define zn_virtual_pmf_zero \
    ZNP_zn_virtual_pmf_zero
void zn_virtual_pmf_zero(zn_virtual_pmf_t res)
{
   // already zero, nothing to do
   if (res->index == -1)
      return;
   
   // detach from buffer, update refcount
   res->parent->ref_count[res->index]--;
   res->index = -1;
}


/*
   Sets res := op, by attaching a slot in the parent vector to the memory
   occupied by op (no data movement is involved).
   
   Note: this means that subsequent changes to op affect the value of res!
*/
#define zn_virtual_pmf_import \
    ZNP_zn_virtual_pmf_import
void zn_virtual_pmf_import(zn_virtual_pmf_t res, zn_pmf_t op)
{
   zn_virtual_pmf_zero(res);

   res->index = zn_virtual_pmf_vec_find_slot(res->parent);
   res->parent->ref_count[res->index] = 1;
   res->parent->external[res->index] = 1;
   res->parent->buf[res->index] = op;
   res->bias = op[0];
}


/*
   Returns a zn_pmf_t with the value of op.
   
   All this does is overwrite the bias field of the underlying zn_pmf_t
   and return a pointer it. It doesn't copy any data. (This is fragile;
   the returned value should be used immediately, before doing anything else
   with the parent vector.)

   If op is zero, the return value is NULL.
*/
#define zn_virtual_pmf_export \
    ZNP_zn_virtual_pmf_export
zn_pmf_t zn_virtual_pmf_export(zn_virtual_pmf_t op)
{
   if (op->index == -1)
      return NULL;

   zn_pmf_t res = op->parent->buf[op->index];
   res[0] = op->bias;
   return res;
}


/*
   Ensures that op has a reference count of 1, by possibly copying the data
   to a new buffer if necessary. Then it's safe to mutate without affecting
   the value of other zn_virtual_pmf_t's.
   
   If op is zero, this is a no-op.
*/
#define zn_virtual_pmf_isolate \
    ZNP_zn_virtual_pmf_isolate
void zn_virtual_pmf_isolate(zn_virtual_pmf_t op)
{
   if (op->index == -1)
      return;

   struct zn_virtual_pmf_vec_struct* parent = op->parent;

   if (parent->ref_count[op->index] == 1)
      // already has reference count 1
      return;
   
   // detach
   parent->ref_count[op->index]--;
   
   // find new buffer and copy the data
   unsigned index = zn_virtual_pmf_vec_new_buf(parent);
   zn_pmf_set(parent->buf[index], parent->buf[op->index], parent->M);
   op->index = index;
}



/* ----------------------------------------------------------------------------

   zn_virtual_pmf_t coefficient operations
   
   These functions all handle reference counting automatically

---------------------------------------------------------------------------- */


/*
   res := op
*/
#define zn_virtual_pmf_set \
    ZNP_zn_virtual_pmf_set
void zn_virtual_pmf_set(zn_virtual_pmf_t res, zn_virtual_pmf_t op)
{
   if (op == res)
      return;
      
   zn_virtual_pmf_zero(res);
   
   if (op->index == -1)
      return;

   res->bias = op->bias;
   res->index = op->index;
   res->parent->ref_count[op->index]++;
}


/*
   op := Y^r * op
*/
#define zn_virtual_pmf_rotate \
    ZNP_zn_virtual_pmf_rotate
void zn_virtual_pmf_rotate(zn_virtual_pmf_t op, ulong r)
{
   if (op->index != -1)
      op->bias += r;
}


/*
   res += op
*/
#define zn_virtual_pmf_add \
    ZNP_zn_virtual_pmf_add
void zn_virtual_pmf_add(zn_virtual_pmf_t res, zn_virtual_pmf_t op)
{
   ZNP_ASSERT(res->parent == op->parent);
   struct zn_virtual_pmf_vec_struct* parent = res->parent;

   // op == 0
   if (op->index == -1)
      return;
      
   // res == 0
   if (res->index == -1)
   {
      zn_virtual_pmf_set(res, op);
      return;
   }

   zn_virtual_pmf_isolate(res);

   zn_pmf_t res_ptr = parent->buf[res->index];
   zn_pmf_t op_ptr = parent->buf[op->index];

   res_ptr[0] = res->bias;
   op_ptr[0] = op->bias;

   zn_pmf_add(res_ptr, op_ptr, parent->M, parent->mod);
}


/*
   res -= op
*/
#define zn_virtual_pmf_sub \
    ZNP_zn_virtual_pmf_sub
void zn_virtual_pmf_sub(zn_virtual_pmf_t res, zn_virtual_pmf_t op)
{
   ZNP_ASSERT(res->parent == op->parent);
   struct zn_virtual_pmf_vec_struct* parent = res->parent;

   // op == 0
   if (op->index == -1)
      return;
      
   // res == 0
   if (res->index == -1)
   {
      zn_virtual_pmf_set(res, op);
      zn_virtual_pmf_rotate(res, parent->M);
      return;
   }
   
   zn_virtual_pmf_isolate(res);

   zn_pmf_t res_ptr = parent->buf[res->index];
   zn_pmf_t op_ptr = parent->buf[op->index];

   res_ptr[0] = res->bias;
   op_ptr[0] = op->bias;

   zn_pmf_sub(res_ptr, op_ptr, parent->M, parent->mod);
}


/*
   op1 := op1 + op2
   op2 := op2 - op1
*/
#define zn_virtual_pmf_bfly \
    ZNP_zn_virtual_pmf_bfly
void zn_virtual_pmf_bfly(zn_virtual_pmf_t op1, zn_virtual_pmf_t op2)
{
   ZNP_ASSERT(op1->parent == op2->parent);
   struct zn_virtual_pmf_vec_struct* parent = op1->parent;
   
   // op1 == 0
   if (op1->index == -1)
   {
      zn_virtual_pmf_set(op1, op2);
      return;
   }
   
   // op2 == 0
   if (op2->index == -1)
   {
      zn_virtual_pmf_set(op2, op1);
      zn_virtual_pmf_rotate(op2, parent->M);
      return;
   }

   zn_virtual_pmf_isolate(op1);
   zn_virtual_pmf_isolate(op2);

   zn_pmf_t op1_ptr = parent->buf[op1->index];
   zn_pmf_t op2_ptr = parent->buf[op2->index];

   op1_ptr[0] = op1->bias;
   op2_ptr[0] = op2->bias;

   zn_pmf_bfly(op1_ptr, op2_ptr, parent->M, parent->mod);
}


/*
   op := op / 2
*/
#define zn_virtual_pmf_divby2 \
    ZNP_zn_virtual_pmf_divby2
void zn_virtual_pmf_divby2(zn_virtual_pmf_t op)
{
   struct zn_virtual_pmf_vec_struct* parent = op->parent;

   if (op->index == -1)
      return;
   
   zn_virtual_pmf_isolate(op);
   zn_pmf_divby2(parent->buf[op->index], parent->M, parent->mod);
}


/* ----------------------------------------------------------------------------

   virtual IFFT routine

---------------------------------------------------------------------------- */

/*
   Performs truncated IFFT on the subvector vec[start, start + J),
   where J = 2^lgJ.

   The meaning of _length_, _forward_ and _twist_ is the same as for
   zn_pmf_vec_ifft (see mul_fft.c).
   
   The algorithm is essentially the same as zn_pmf_vec_ifft_small, except we
   don't worry about the nonzero parameter (this is handled automatically
   by the underlying optimised representation).
*/
#define zn_virtual_pmf_vec_ifft \
    ZNP_zn_virtual_pmf_vec_ifft
void zn_virtual_pmf_vec_ifft(
            zn_virtual_pmf_vec_t vec, ulong start, unsigned lgJ,
            ulong length, int forward, ulong twist)
{
   ZNP_ASSERT(lgJ <= vec->lgM + 1);
   ZNP_ASSERT((twist << lgJ) < 2*vec->M);
   ZNP_ASSERT(length + forward <= (1UL << lgJ));
   ZNP_ASSERT(lgJ <= vec->lgK);

   if (lgJ == 0)
      return;
   
   const zn_mod_struct* mod = vec->mod;
   ulong J = 1UL << lgJ;
   ulong M = vec->M;
   ulong s, r = M >> (lgJ - 1);     // 2M/J
   long i;

   if (length + forward <= J/2)
   {
      for (i = J/2 - 1; i >= (long)length; i--)
      {
         zn_virtual_pmf_add(vec->data[start + i], vec->data[start + i + J/2]);
         zn_virtual_pmf_divby2(vec->data[start + i]);
      }
      
      zn_virtual_pmf_vec_ifft(vec, start, lgJ - 1,
                              length, forward, twist << 1);
      
      for (; i >= 0; i--)
      {
         zn_virtual_pmf_add(vec->data[start + i], vec->data[start + i]);
         zn_virtual_pmf_sub(vec->data[start + i], vec->data[start + i + J/2]);
      }
   }
   else
   {
      zn_virtual_pmf_vec_ifft(vec, start, lgJ - 1, J/2, 0, twist << 1);

      for (i = J/2 - 1, s = twist + r * i;
           i >= (long)(length - J/2); i--, s -= r)
      {
         zn_virtual_pmf_sub(vec->data[start + i + J/2], vec->data[start + i]);
         zn_virtual_pmf_sub(vec->data[start + i], vec->data[start + i + J/2]);
         zn_virtual_pmf_rotate(vec->data[start + i + J/2], M + s);
      }

      zn_virtual_pmf_vec_ifft(vec, start + J/2, lgJ - 1,
                             length - J/2, forward, twist << 1);
      
      for (; i >= 0; i--, s -= r)
      {
         zn_virtual_pmf_rotate(vec->data[start + i + J/2], M - s);
         zn_virtual_pmf_bfly(vec->data[start + i + J/2], vec->data[start + i]);
      }
   }
}



/* ============================================================================

     main array multiplication routine

============================================================================ */


/*
   The idea of this routine is as follows. We simulate the algorithm used
   in zn_array_mul_fft, in the case that the FFTs are performed via
   zn_pmf_vec_fft_factor() and zn_pmf_vec_ifft_factor(), with K factored into
   T = 2^lgT rows and U = 2^lgU columns. Here we assume that T is quite small
   and U is possibly very large.
   
   However, instead of storing the whole fourier transform, we only work on a
   single row at a time. This means we have to store at most three rows
   simultaneously: the two rows whose transforms are being multiplied together,
   and the "partial" row, which enters the computation right at the
   beginning and is needed until the end. (In the ordinary mul_fft routine,
   we need space for 2T rows simultaneously.)
   
   This also means we cannot use *fast* fourier transforms for the columns,
   since we don't have all the data available. They are done by a naive DFT
   instead.
*/

void zn_array_mul_fft_dft(ulong* res, const ulong* op1, size_t len1,
                          const ulong* op2, size_t len2, unsigned lgT,
                          const zn_mod_t mod)
{
   ZNP_ASSERT(mod->n & 1);
   ZNP_ASSERT(len2 >= 1);
   ZNP_ASSERT(len1 >= len2);
   
   if (lgT == 0)
   {
      // no layers of DFT; just call usual FFT routine
      int squaring = (op1 == op2) && (len1 == len2);
      ulong scale = zn_array_mul_fft_get_fudge(len1, len2, squaring, mod);
      zn_array_mul_fft(res, op1, len1, op2, len2, scale, mod);
      return;
   }

   unsigned lgM, lgK;

   // number of zn_pmf_t coefficients for each input poly
   ulong coeffs1, coeffs2;

   // figure out how big the transform needs to be
   mul_fft_params(&lgK, &lgM, &coeffs1, &coeffs2, len1, len2);

   // number of zn_pmf_t coefficients for output poly
   ulong length = coeffs1 + coeffs2 - 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ptrdiff_t skip = M + 1;

   size_t len_res = len1 + len2 - 1;

   // Split up transform into length K = U * T, i.e. U columns and T rows.
   if (lgT >= lgK)
      lgT = lgK;
   unsigned lgU = lgK - lgT;
   ulong U = 1UL << lgU;
   ulong T = 1UL << lgT;
   
   // space for two input rows, and one partial row
   zn_pmf_vec_t in1, in2, partial;
   zn_pmf_vec_init(in1, lgU, skip, lgM, mod);
   zn_pmf_vec_init(in2, lgU, skip, lgM, mod);
   zn_pmf_vec_init(partial, lgU, skip, lgM, mod);

   // the virtual pmf_vec_t that we use for the column DFTs
   zn_virtual_pmf_vec_t col;
   zn_virtual_pmf_vec_init(col, lgT, lgM, mod);

   // zero the output
   zn_array_zero(res, len_res);
   
   long i, j, k, which;
   int fudge;

   // Write length = U * length_T + length_U, where 0 <= length_U < U
   ulong length_U = length & (U - 1);
   ulong length_T = length >> lgU;

   // for each row (beginning with the last partial row if it exists)....
   for (i = length_T - (length_U == 0); i >= 0; i--)
   {
      ulong i_rev = bit_reverse(i, lgT);
      
      // for each input array....
      for (which = 0; which < 2; which++)
      {
         zn_pmf_vec_struct* in = which ? in2 : in1;
         const ulong* op = which ? op2 : op1;
         size_t len = which ? len2 : len1;
         
         zn_pmf_t in_ptr = in->data;

         for (j = 0; j < U; j++, in_ptr += in->skip)
         {
            // compute the i-th row of the j-th column as it would look after
            // the column FFTs, using naive DFT
            zn_pmf_zero(in_ptr, M);
            ulong r = i_rev << (lgM - lgT + 1);
            
            for (k = 0; k < T; k++)
            {
               merge_chunk_to_pmf(in_ptr, op, len, (k * U + j) * M/2, M, mod);
               zn_pmf_rotate(in_ptr, -r);
            }
            
            zn_pmf_rotate(in_ptr, (i_rev * j) << (lgM - lgK + 1));
         }
         
         // Now we've got the whole row; run FFT on the row
         zn_pmf_vec_fft(in, (i == length_T) ? length_U : U, U, 0);
      }

      if (i == length_T)
      {
         // pointwise multiply the two partial rows
         zn_pmf_vec_mul(partial, in1, in2, length_U, i == 0);
         // remove fudge factor
         zn_pmf_vec_scalar_mul(partial, length_U,
                               zn_pmf_vec_mul_get_fudge(lgM, 0, mod));

         // zero remainder of the partial row; we will subsequently add
         // in contributions from the vertical IFFTs when we process the other
         // rows.
         for (j = length_U; j < U; j++)
            zn_pmf_zero(partial->data + partial->skip * j, M);
      }
      else
      {
         // pointwise multiply the two rows
         zn_pmf_vec_mul(in1, in1, in2, U, i == 0);
         // remove fudge factor
         zn_pmf_vec_scalar_mul(in1, U, zn_pmf_vec_mul_get_fudge(lgM, 0, mod));
         
         // horizontal IFFT this row
         zn_pmf_vec_ifft(in1, U, 0, U, 0);

         // simulate vertical IFFTs with DFTs
         for (j = 0; j < U; j++)
         {
            zn_virtual_pmf_vec_reset(col);
            zn_virtual_pmf_import(col->data[i], in1->data + in1->skip * j);
            zn_virtual_pmf_vec_ifft(col, 0, lgT, length_T + (j < length_U),
                                    (j >= length_U) && length_U,
                                    j << (lgM + 1 - lgK));
            
            if ((j >= length_U) && length_U)
            {
               // add contribution to partial row (only for rightmost columns)
               zn_pmf_t src = zn_virtual_pmf_export(col->data[length_T]);
               if (src)
                  zn_pmf_add(partial->data + partial->skip * j, src, M, mod);
            }

            // add contributions to output
            for (k = 0; k < length_T + (j < length_U); k++)
               merge_chunk_from_pmf(res, len_res,
                                    zn_virtual_pmf_export(col->data[k]),
                                    (k * U + j) * M/2, M, mod);
         }
      }
   }

   // now finish off the partial row
   if (length_U)
   {
      // horizontal IFFT partial row
      zn_pmf_vec_ifft(partial, length_U, 0, U, 0);

      // simulate leftmost vertical IFFTs
      for (j = 0; j < length_U; j++)
      {
         zn_virtual_pmf_vec_reset(col);
         zn_virtual_pmf_import(col->data[length_T],
                               partial->data + partial->skip * j);
         zn_virtual_pmf_vec_ifft(col, 0, lgT, length_T + 1, 0,
                                 j << (lgM + 1 - lgK));
                         
         // add contributions to output
         for (k = 0; k <= length_T; k++)
            merge_chunk_from_pmf(res, len_res,
                                 zn_virtual_pmf_export(col->data[k]),
                                 (k * U + j) * M/2, M, mod);
      }
   }
   
   // normalise result
   zn_array_scalar_mul(res, res, len_res, zn_mod_pow2(-lgK, mod), mod);

   zn_virtual_pmf_vec_clear(col);
   zn_pmf_vec_clear(partial);
   zn_pmf_vec_clear(in2);
   zn_pmf_vec_clear(in1);
}


// end of file ****************************************************************
