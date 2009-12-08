/*
   mul_fft_dft.c:  multiplication by Schonhage/Nussbaumer FFT, with a few
                   layers of naive DFT to save memory
   
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


/*
   Returns length n bit reversal of x.
*/
#define bit_reverse \
    ZNP_bit_reverse
ulong
bit_reverse (ulong x, unsigned n)
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
   Let [a, b) = intersection of [0, n) and [k, k + M/2).
   This functions adds op[a, b) to res.
*/
#define merge_chunk_to_pmf \
    ZNP_merge_chunk_to_pmf
void
merge_chunk_to_pmf (pmf_t res, const ulong* op, size_t n, size_t k, ulong M,
                    const zn_mod_t mod)
{
   ZNP_ASSERT ((M & 1) == 0);
   
   ulong r = (-res[0]) & (2*M - 1);
   
   size_t end = k + M/2;
   if (end > n)
      end = n;
   if (k >= end)
      // nothing to do
      return;
   
   op += k;
   ulong size = end - k;
   // now we need to handle op[0, size), and we are guaranteed size <= M/2.

   if (r < M)
   {
      if (size <= M - r)
         zn_array_add_inplace (res + 1 + r, op, size, mod);
      else
      {
         zn_array_add_inplace (res + 1 + r, op, M - r, mod);
         // negacyclic wraparound:
         zn_array_sub_inplace (res + 1, op + M - r, size - M + r, mod);
      }
   }
   else
   {
      r -= M;

      if (size <= M - r)
         zn_array_sub_inplace (res + 1 + r, op, size, mod);
      else
      {
         zn_array_sub_inplace (res + 1 + r, op, M - r, mod);
         // negacyclic wraparound:
         zn_array_add_inplace (res + 1, op + M - r, size - M + r, mod);
      }
   }
}


/*
   Adds op into res, starting at index k, and not writing beyond res + n.
   
   If op == NULL, does no operation.
*/
#define merge_chunk_from_pmf \
    ZNP_merge_chunk_from_pmf
void
merge_chunk_from_pmf (ulong* res, size_t n, const pmf_t op, size_t k,
                      ulong M, const zn_mod_t mod)
{
   if (op == NULL)
      return;

   size_t end = k + M;
   if (end > n)
      end = n;
   if (k >= end)
      // nothing to do
      return;
   
   res += k;
   ulong size = end - k;
   // now we need to write to res[0, size), and we are guaranteed size <= M.

   ulong r = op[0] & (2*M - 1);

   if (r < M)
   {
      if (size <= r)
         zn_array_sub_inplace (res, op + 1 + M - r, size, mod);
      else
      {
         zn_array_sub_inplace (res, op + 1 + M - r, r, mod);
         // negacyclic wraparound:
         zn_array_add_inplace (res + r, op + 1, size - r, mod);
      }
   }
   else
   {
      r -= M;
   
      if (size <= r)
         zn_array_add_inplace (res, op + 1 + M - r, size, mod);
      else
      {
         zn_array_add_inplace (res, op + 1 + M - r, r, mod);
         // negacyclic wraparound:
         zn_array_sub_inplace (res + r, op + 1, size - r, mod);
      }
   }
}


/* ============================================================================

     "virtual" pmf_t's

============================================================================ */


/*
   The virtual_pmf_t and virtual_pmfvec_t are similar to their non-virtual
   counterparts (pmf_t and pmfvec_t), but the underlying representation is
   optimised for the case where:
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


#define virtual_pmfvec_struct \
    ZNP_virtual_pmfvec_struct
struct virtual_pmfvec_struct;     // forward declaration


/*
   Virtual version of a pmf_t.
   
   Each virtual_pmf_t belongs to a "parent" virtual_pmfvec_t.
   
   The index field is:
     * -1 if the value is zero
     * otherwise, an index into parent's buf, which is where the actual
       coefficient data is stored.
   
   The bias field overrides the bias word in the underlying pmf_t.
   (This lets different virtual_pmf's share memory but still have different
   bias values.)
*/
#define virtual_pmf_struct \
    ZNP_virtual_pmf_struct
struct virtual_pmf_struct
{
   struct virtual_pmfvec_struct* parent;
   int index;
   ulong bias;
};

#define virtual_pmf_t \
    ZNP_virtual_pmf_t
typedef struct virtual_pmf_struct  virtual_pmf_t[1];


/*
   Virtual version of pmfvec_t.
   
   M, lgM, K, lgK, mod are as for pmfvec_t.
   
   data is an array of K virtual_pmf_t's (the coefficients in the vector).
   
   The underlying data is managed by three arrays of length max_buffers:
     * buf[i] points to a pmf_t, or NULL if the i-th slot is not yet
       associated to any actual memory.
     * count[i] is a reference count for slot #i, i.e. the number of
       pmf_virtual_pmf_t's pointing to this slot.
     * external[i] is a flag indicating whether the memory belongs to
       someone else (i.e. it's not the responsiblity of the virtual_pmfvec
       to free it).
*/
struct virtual_pmfvec_struct
{
   ulong M;
   unsigned lgM;

   ulong K;
   unsigned lgK;
   
   const zn_mod_struct* mod;

   virtual_pmf_t* data;

   unsigned max_buffers;
   ulong** buf;
   unsigned* count;
   int* external;
};

#define virtual_pmfvec_t \
    ZNP_virtual_pmfvec_t
typedef struct virtual_pmfvec_struct  virtual_pmfvec_t[1];



/* ----------------------------------------------------------------------------

   virtual_pmf_t infrastructure

---------------------------------------------------------------------------- */

/*
   Initialises a virtual_pmf_t to zero, with a given parent vector.
*/
#define virtual_pmf_init \
    ZNP_virtual_pmf_init
void
virtual_pmf_init (virtual_pmf_t res, virtual_pmfvec_t parent)
{
   res->index = -1;
   res->parent = parent;
}


/*
   Initialises a virtual_pmfvec_t to length K, with all zero values.
   
   All slots are initially marked as empty.
*/
#define virtual_pmfvec_init \
    ZNP_virtual_pmfvec_init
void
virtual_pmfvec_init (virtual_pmfvec_t vec, unsigned lgK, unsigned lgM,
                     const zn_mod_t mod)
{
   vec->mod = mod;

   vec->lgM = lgM;
   vec->M = 1UL << lgM;
   
   vec->lgK = lgK;
   vec->K = 1UL << lgK;
   
   vec->data = (virtual_pmf_t*) malloc (vec->K * sizeof (virtual_pmf_t));

   ulong i;
   for (i = 0; i < vec->K; i++)
      virtual_pmf_init (vec->data[i], vec);

   vec->max_buffers = 2 * vec->K;     // should be safe

   vec->buf = (ulong**) malloc (sizeof (ulong*) * vec->max_buffers);
   vec->count = (unsigned*) malloc (sizeof (unsigned) * vec->max_buffers);
   vec->external = (int*) malloc (sizeof (int) * vec->max_buffers);
   
   for (i = 0; i < vec->max_buffers; i++)
   {
      vec->buf[i] = NULL;
      vec->count[i] = 0;
      vec->external[i] = 0;
   }
}


/*
   Sets all values to zero, and detaches externally-allocated memory from
   all slots.
   
   This does *not* free any memory owned by this vector; buffers already
   allocated will get re-used by subsequent operations.
*/
#define virtual_pmfvec_reset \
    ZNP_virtual_pmfvec_reset
void
virtual_pmfvec_reset (virtual_pmfvec_t vec)
{
   ulong i;
   for (i = 0; i < vec->K; i++)
      vec->data[i]->index = -1;

   for (i = 0; i < vec->max_buffers; i++)
   {
      vec->count[i] = 0;
      
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
#define virtual_pmfvec_clear \
    ZNP_virtual_pmfvec_clear
void
virtual_pmfvec_clear (virtual_pmfvec_t vec)
{
   virtual_pmfvec_reset (vec);

   ulong i;
   for (i = 0; i < vec->max_buffers; i++)
      if (vec->buf[i] && !vec->external[i])
         free (vec->buf[i]);

   free (vec->external);
   free (vec->buf);
   free (vec->count);
   free (vec->data);
}


/*
   Finds a free slot (one not attached to any underlying pmf_t yet),
   and returns its index.
*/
#define virtual_pmfvec_find_slot \
    ZNP_virtual_pmfvec_find_slot
unsigned
virtual_pmfvec_find_slot (virtual_pmfvec_t vec)
{
   unsigned i;
   for (i = 0; i < vec->max_buffers; i++)
      if (!vec->buf[i])
         return i;
   
   // this should never happen; we always should have enough slots
   ZNP_ASSERT (0);
}


/*
   Finds a slot attached to an underlying pmf_t which is not currently
   used by any other virtual_pmf_t's, and returns its index.
   
   If there are no such slots, it allocated more space, attaches a slot to it,
   and returns its index.
   
   In both cases, the reference count of the returned slot will be 1.
*/
#define virtual_pmfvec_new_buf \
    ZNP_virtual_pmfvec_new_buf
unsigned
virtual_pmfvec_new_buf (virtual_pmfvec_t vec)
{
   // first search for an already-allocated buffer that no-one else is using
   unsigned i;
   for (i = 0; i < vec->max_buffers; i++)
      if (vec->buf[i] && !vec->count[i])
         break;
   
   if (i == vec->max_buffers)
   {
      // not found; need to allocate more space
      i = virtual_pmfvec_find_slot (vec);
      vec->buf[i] = (ulong*) malloc (sizeof (ulong) * (vec->M + 1));
      vec->external[i] = 0;
   }
   
   vec->count[i] = 1;
   return i;
}


/*
   res := 0
*/
#define virtual_pmf_zero \
    ZNP_virtual_pmf_zero
void
virtual_pmf_zero (virtual_pmf_t res)
{
   // already zero, nothing to do
   if (res->index == -1)
      return;
   
   // detach from buffer, update refcount
   res->parent->count[res->index]--;
   res->index = -1;
}


/*
   Sets res := op, by attaching a slot in the parent vector to the memory
   occupied by op (no data movement is involved).
   
   Note: this means that subsequent changes to op affect the value of res!
*/
#define virtual_pmf_import \
    ZNP_virtual_pmf_import
void
virtual_pmf_import (virtual_pmf_t res, pmf_t op)
{
   virtual_pmf_zero (res);

   res->index = virtual_pmfvec_find_slot (res->parent);
   res->parent->count[res->index] = 1;
   res->parent->external[res->index] = 1;
   res->parent->buf[res->index] = op;
   res->bias = op[0];
}


/*
   Returns a pmf_t with the value of op.
   
   All this does is overwrite the bias field of the underlying pmf_t
   and return a pointer it. It doesn't copy any data. (This is fragile;
   the returned value should be used immediately, before doing anything else
   with the parent vector.)

   If op is zero, the return value is NULL.
*/
#define virtual_pmf_export \
    ZNP_virtual_pmf_export
pmf_t
virtual_pmf_export (virtual_pmf_t op)
{
   if (op->index == -1)
      return NULL;

   pmf_t res = op->parent->buf[op->index];
   res[0] = op->bias;
   return res;
}


/*
   Ensures that op has a reference count of 1, by possibly copying the data
   to a new buffer if necessary. Then it's safe to mutate without affecting
   the value of other virtual_pmf_t's.
   
   If op is zero, this is a no-op.
*/
#define virtual_pmf_isolate \
    ZNP_virtual_pmf_isolate
void
virtual_pmf_isolate (virtual_pmf_t op)
{
   if (op->index == -1)
      return;

   struct virtual_pmfvec_struct* parent = op->parent;

   if (parent->count[op->index] == 1)
      // already has reference count 1
      return;
   
   // detach
   parent->count[op->index]--;
   
   // find new buffer and copy the data
   unsigned index = virtual_pmfvec_new_buf (parent);
   pmf_set (parent->buf[index], parent->buf[op->index], parent->M);
   op->index = index;
}



/* ----------------------------------------------------------------------------

   virtual_pmf_t coefficient operations
   
   These functions all handle reference counting automatically

---------------------------------------------------------------------------- */


/*
   res := op
*/
#define virtual_pmf_set \
    ZNP_virtual_pmf_set
void
virtual_pmf_set (virtual_pmf_t res, virtual_pmf_t op)
{
   if (op == res)
      return;
      
   virtual_pmf_zero (res);
   
   if (op->index == -1)
      return;

   res->bias = op->bias;
   res->index = op->index;
   res->parent->count[op->index]++;
}


/*
   op := Y^r * op
*/
#define virtual_pmf_rotate \
    ZNP_virtual_pmf_rotate
void
virtual_pmf_rotate (virtual_pmf_t op, ulong r)
{
   if (op->index != -1)
      op->bias += r;
}


/*
   res += op
*/
#define virtual_pmf_add \
    ZNP_virtual_pmf_add
void
virtual_pmf_add (virtual_pmf_t res, virtual_pmf_t op)
{
   ZNP_ASSERT (res->parent == op->parent);
   struct virtual_pmfvec_struct* parent = res->parent;

   // op == 0
   if (op->index == -1)
      return;
      
   // res == 0
   if (res->index == -1)
   {
      virtual_pmf_set (res, op);
      return;
   }

   virtual_pmf_isolate (res);

   pmf_t p2 = parent->buf[res->index];
   pmf_t p1 = parent->buf[op->index];

   p2[0] = res->bias;
   p1[0] = op->bias;

   pmf_add (p2, p1, parent->M, parent->mod);
}


/*
   res -= op
*/
#define virtual_pmf_sub \
    ZNP_virtual_pmf_sub
void
virtual_pmf_sub (virtual_pmf_t res, virtual_pmf_t op)
{
   ZNP_ASSERT (res->parent == op->parent);
   struct virtual_pmfvec_struct* parent = res->parent;

   // op == 0
   if (op->index == -1)
      return;
      
   // res == 0
   if (res->index == -1)
   {
      virtual_pmf_set (res, op);
      virtual_pmf_rotate (res, parent->M);
      return;
   }
   
   virtual_pmf_isolate (res);

   pmf_t p2 = parent->buf[res->index];
   pmf_t p1 = parent->buf[op->index];

   p2[0] = res->bias;
   p1[0] = op->bias;

   pmf_sub (p2, p1, parent->M, parent->mod);
}


/*
   op1 := op1 + op2
   op2 := op2 - op1
*/
#define virtual_pmf_bfly \
    ZNP_virtual_pmf_bfly
void
virtual_pmf_bfly (virtual_pmf_t op1, virtual_pmf_t op2)
{
   ZNP_ASSERT (op1->parent == op2->parent);
   struct virtual_pmfvec_struct* parent = op1->parent;
   
   // op1 == 0
   if (op1->index == -1)
   {
      virtual_pmf_set (op1, op2);
      return;
   }
   
   // op2 == 0
   if (op2->index == -1)
   {
      virtual_pmf_set (op2, op1);
      virtual_pmf_rotate (op2, parent->M);
      return;
   }

   virtual_pmf_isolate (op1);
   virtual_pmf_isolate (op2);

   pmf_t p1 = parent->buf[op1->index];
   pmf_t p2 = parent->buf[op2->index];

   p1[0] = op1->bias;
   p2[0] = op2->bias;

   pmf_bfly (p1, p2, parent->M, parent->mod);
}


/*
   op := op / 2
*/
#define virtual_pmf_divby2 \
    ZNP_virtual_pmf_divby2
void
virtual_pmf_divby2 (virtual_pmf_t op)
{
   struct virtual_pmfvec_struct* parent = op->parent;

   if (op->index == -1)
      return;
   
   virtual_pmf_isolate (op);
   pmf_divby2 (parent->buf[op->index], parent->M, parent->mod);
}


/* ----------------------------------------------------------------------------

   virtual IFFT routine

---------------------------------------------------------------------------- */

/*
   Performs truncated IFFT on vec. The meanings of n, fwd and t are the same
   as for pmfvec_ifft (see pmfvec_fft.c).
   
   The algorithm is essentially the same as pmfvec_ifft_dc(), except that we
   don't worry about the z parameter (zero entries are handled automatically
   by the underlying optimised representation).
*/
#define virtual_pmfvec_ifft \
    ZNP_virtual_pmfvec_ifft
void
virtual_pmfvec_ifft (virtual_pmfvec_t vec, ulong n, int fwd, ulong t)
{
   ZNP_ASSERT (vec->lgK <= vec->lgM + 1);
   ZNP_ASSERT (t * vec->K < 2 * vec->M);
   ZNP_ASSERT (n + fwd <= vec->K);

   if (vec->lgK == 0)
      return;
      
   vec->lgK--;
   vec->K >>= 1;
   
   const zn_mod_struct* mod = vec->mod;
   virtual_pmf_t* data = vec->data;
   ulong M = vec->M;
   ulong K = vec->K;
   ulong s, r = M >> vec->lgK;
   long i;

   if (n + fwd <= K)
   {
      for (i = K - 1; i >= (long) n; i--)
      {
         virtual_pmf_add (data[i], data[i + K]);
         virtual_pmf_divby2 (data[i]);
      }

      virtual_pmfvec_ifft (vec, n, fwd, t << 1);
      
      for (; i >= 0; i--)
      {
         virtual_pmf_add (data[i], data[i]);
         virtual_pmf_sub (data[i], data[i + K]);
      }
   }
   else
   {
      virtual_pmfvec_ifft (vec, K, 0, t << 1);

      for (i = K - 1, s = t + r * i; i >= (long)(n - K); i--, s -= r)
      {
         virtual_pmf_sub (data[i + K], data[i]);
         virtual_pmf_sub (data[i], data[i + K]);
         virtual_pmf_rotate (data[i + K], M + s);
      }

      vec->data += K;
      virtual_pmfvec_ifft (vec, n - K, fwd, t << 1);
      vec->data -= K;
      
      for (; i >= 0; i--, s -= r)
      {
         virtual_pmf_rotate (data[i + K], M - s);
         virtual_pmf_bfly (data[i + K], data[i]);
      }
   }

   vec->K <<= 1;
   vec->lgK++;
}



/* ============================================================================

     main array multiplication routine

============================================================================ */


/*
   The idea of this routine is as follows. We simulate the algorithm used
   in zn_array_mul_fft, in the case that the FFTs are performed via
   pmfvec_fft_huge() and pmfvec_ifft_huge(), with K factored into
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
   instead. The total number of coefficient operations (adds/subs) is
   O(T * U * log(U) + T^2 * U).
*/

void
zn_array_mul_fft_dft (ulong* res,
                      const ulong* op1, size_t n1,
                      const ulong* op2, size_t n2,
                      unsigned lgT, const zn_mod_t mod)
{
   ZNP_ASSERT (mod->m & 1);
   ZNP_ASSERT (n2 >= 1);
   ZNP_ASSERT (n1 >= n2);
   
   if (lgT == 0)
   {
      // no layers of DFT; just call usual FFT routine
      int sqr = (op1 == op2) && (n1 == n2);
      ulong x = zn_array_mul_fft_fudge (n1, n2, sqr, mod);
      zn_array_mul_fft (res, op1, n1, op2, n2, x, mod);
      return;
   }

   unsigned lgM, lgK;

   // number of pmf_t coefficients for each input poly
   ulong m1, m2;

   // figure out how big the transform needs to be
   mul_fft_params (&lgK, &lgM, &m1, &m2, n1, n2);

   // number of pmf_t coefficients for output poly
   ulong m = m1 + m2 - 1;

   ulong M = 1UL << lgM;
   ulong K = 1UL << lgK;
   ptrdiff_t skip = M + 1;

   size_t n3 = n1 + n2 - 1;

   // Split up transform into length K = U * T, i.e. U columns and T rows.
   if (lgT >= lgK)
      lgT = lgK;
   unsigned lgU = lgK - lgT;
   ulong U = 1UL << lgU;
   ulong T = 1UL << lgT;
   
   // space for two input rows, and one partial row
   pmfvec_t in1, in2, part;
   pmfvec_init (in1, lgU, skip, lgM, mod);
   pmfvec_init (in2, lgU, skip, lgM, mod);
   pmfvec_init (part, lgU, skip, lgM, mod);

   // the virtual pmfvec_t that we use for the column DFTs
   virtual_pmfvec_t col;
   virtual_pmfvec_init (col, lgT, lgM, mod);

   // zero the output
   zn_array_zero (res, n3);
   
   long i, j, k;
   int which;

   // Write m = U * mT + mU, where 0 <= mU < U
   ulong mU = m & (U - 1);
   ulong mT = m >> lgU;

   // for each row (beginning with the last partial row if it exists)....
   for (i = mT - (mU == 0); i >= 0; i--)
   {
      ulong i_rev = bit_reverse (i, lgT);
      
      // for each input array....
      for (which = 0; which < 2; which++)
      {
         pmfvec_struct* in = which ? in2 : in1;
         const ulong* op = which ? op2 : op1;
         size_t n = which ? n2 : n1;

         pmf_t p = in->data;

         for (j = 0; j < U; j++, p += in->skip)
         {
            // compute the i-th row of the j-th column as it would look after
            // the column FFTs, using naive DFT
            pmf_zero (p, M);
            ulong r = i_rev << (lgM - lgT + 1);
            
            for (k = 0; k < T; k++)
            {
               merge_chunk_to_pmf (p, op, n, (k * U + j) << (lgM - 1), M, mod);
               pmf_rotate (p, -r);
            }
            
            pmf_rotate (p, (i_rev * j) << (lgM - lgK + 1));
         }
         
         // Now we've got the whole row; run FFT on the row
         pmfvec_fft (in, (i == mT) ? mU : U, U, 0);
      }

      if (i == mT)
      {
         // pointwise multiply the two partial rows
         pmfvec_mul (part, in1, in2, mU, i == 0);
         // remove fudge factor
         pmfvec_scalar_mul (part, mU, pmfvec_mul_fudge (lgM, 0, mod));

         // zero remainder of the partial row; we will subsequently add
         // in contributions from the vertical IFFTs when we process the other
         // rows.
         for (j = mU; j < U; j++)
            pmf_zero (part->data + part->skip * j, M);
      }
      else
      {
         // pointwise multiply the two rows
         pmfvec_mul (in1, in1, in2, U, i == 0);
         // remove fudge factor
         pmfvec_scalar_mul (in1, U, pmfvec_mul_fudge (lgM, 0, mod));
         
         // horizontal IFFT this row
         pmfvec_ifft (in1, U, 0, U, 0);

         // simulate vertical IFFTs with DFTs
         for (j = 0; j < U; j++)
         {
            virtual_pmfvec_reset (col);
            virtual_pmf_import (col->data[i], in1->data + in1->skip * j);
            virtual_pmfvec_ifft (col, mT + (j < mU), (j >= mU) && mU,
                                 j << (lgM + 1 - lgK));
            
            if ((j >= mU) && mU)
            {
               // add contribution to partial row (only for rightmost columns)
               pmf_t src = virtual_pmf_export (col->data[mT]);
               if (src)
                  pmf_add (part->data + part->skip * j, src, M, mod);
            }

            // add contributions to output
            for (k = 0; k < mT + (j < mU); k++)
               merge_chunk_from_pmf (res, n3,
                                     virtual_pmf_export (col->data[k]),
                                     (k * U + j) * M/2, M, mod);
         }
      }
   }

   // now finish off the partial row
   if (mU)
   {
      // horizontal IFFT partial row
      pmfvec_ifft (part, mU, 0, U, 0);

      // simulate leftmost vertical IFFTs
      for (j = 0; j < mU; j++)
      {
         virtual_pmfvec_reset (col);
         virtual_pmf_import (col->data[mT], part->data + part->skip * j);
         virtual_pmfvec_ifft (col, mT + 1, 0, j << (lgM + 1 - lgK));
                         
         // add contributions to output
         for (k = 0; k <= mT; k++)
            merge_chunk_from_pmf (res, n3,
                                  virtual_pmf_export (col->data[k]),
                                  (k * U + j) * M/2, M, mod);
      }
   }
   
   // normalise result
   zn_array_scalar_mul (res, res, n3, zn_mod_pow2 (-lgK, mod), mod);

   virtual_pmfvec_clear (col);
   pmfvec_clear (part);
   pmfvec_clear (in2);
   pmfvec_clear (in1);
}


// end of file ****************************************************************
