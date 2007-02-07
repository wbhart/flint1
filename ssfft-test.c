#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <assert.h>
#include "flint.h"
#include "ssfft.h"

gmp_randstate_t ssfft_test_randstate;


/*
Fills dest with n+1 limbs of random data. It uses mpz_rrandomb to get long
strings of 0's and 1's. The top spare_bits bits are all equal to the sign bit.
spare_bits must be in the range 0 < spare_bits <= FLINT_BITS_PER_LIMB.
*/
void raw_rrandomb(mp_limb_t* dest, unsigned long n, unsigned long spare_bits,
                  gmp_randstate_t state)
{
   mpz_t temp;
   mpz_init(temp);

   unsigned long bits = (n+1)*FLINT_BITS_PER_LIMB - spare_bits + 1;
   mpz_rrandomb(temp, state, bits);
   for (unsigned long i = 0; i <= n; i++)
      dest[i] = 0;
   mpz_export(dest, NULL, -1, sizeof(mp_limb_t), 0, 0, temp);
   
   // sign-extend
   if (mpz_tstbit(temp, bits-1))
      dest[n] |= ~((1UL << (FLINT_BITS_PER_LIMB - spare_bits)) - 1);

   // mpz_rrandomb seems to have a "bug" where the high bit of its output
   // is always 1. We compensate by flipping all the bits with probability 0.5.
   if (gmp_urandomb_ui(state, 1))
      for (unsigned long i = 0; i <= n; i++)
         dest[i] = ~dest[i];
   
   mpz_clear(temp);
}


/*
Counts the number of high bits of x that are equal to the sign bit.
e.g. count_spare_bits(000111001) = 3, count_spare_bits(111010011) = 3
*/
unsigned long count_spare_bits(mp_limb_signed_t x)
{
   if (x < 0)
      x = ~x;
   unsigned long i;
      for (i = 0; x; x >>= 1, i++);
   return FLINT_BITS_PER_LIMB - i;
}


mp_limb_t random_limb()
{
   return gmp_urandomb_ui(ssfft_test_randstate, FLINT_BITS_PER_LIMB);
}


void ssfft_print_buf(mp_limb_t* x, unsigned long n)
{
   for (long i = n; i >= 0; i--)
      printf("%08lx ", x[i]);
}


void ssfft_negate_limbs(mp_limb_t* dest, mp_limb_t* src, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = ~src[i];
   mpn_add_1(dest, dest, count, 1);
}


void ssfft_convert_from_raw(mpz_t output, mp_limb_t* input, unsigned long n)
{
   int negative = ((mp_limb_signed_t) input[n]) < 0;
   
   if (negative)
      ssfft_negate_limbs(input, input, n+1);

   mpz_import(output, n+1, -1, sizeof(mp_limb_t), 0, 0, input);

   if (negative)
   {
      ssfft_negate_limbs(input, input, n+1);
      mpz_neg(output, output);
   }
}


/*
y := x * 2^(s/2)  mod p    (using a very naive algorithm)
y may alias x
*/
void ssfft_naive_rotate(mpz_t y, mpz_t x, unsigned long s,
                        unsigned long n, mpz_t p)
{
   if (s & 1)
   {
      mpz_t temp;
      mpz_init(temp);
      mpz_mul_2exp(y, x, s/2 + n*FLINT_BITS_PER_LIMB/4);
      mpz_mul_2exp(temp, y, n*FLINT_BITS_PER_LIMB/2);
      mpz_sub(y, temp, y);
      mpz_mod(y, y, p);
      mpz_clear(temp);
   }
   else
   {
      mpz_mul_2exp(y, x, s/2);
      mpz_mod(y, y, p);
   }
}


/*
y := x * 2^(-s/2)  mod p    (using a very naive algorithm)
y may alias x
*/
void ssfft_naive_unrotate(mpz_t y, mpz_t x, unsigned long s,
                          unsigned long n, mpz_t p)
{
   ssfft_naive_rotate(y, x, 4*n*FLINT_BITS_PER_LIMB - s, n, p);
}


// ============================================================================
// coefficient operation test case support

/*
An array of K buffers of length n+1.
Handles memory management, buffer overflow checks.
*/
typedef struct buf_array_t
{
   unsigned long K;        // number of buffers
   unsigned long n;        // length of buffers
   mpz_t p;                // p = 2^(Bn) + 1
   mp_limb_t** buffers;    // length K vector of pointers to buffers
   mpz_t* mpz_buffers;     // mpz_t versions of the buffers

   // private:
   mp_limb_t* data;       // block of allocated memory
   mp_limb_t* sentries;   // length 2K vector of sentries (two for each buffer)
} buf_array_t;


/*
Initialises the given array with K buffers of length n+1.
*/
void buf_array_init(buf_array_t* array, unsigned long K, unsigned long n)
{
   array->K = K;
   array->n = n;

   // need n+3 limbs for each buffer (n+1 data plus a sentry at each end)
   array->data = (mp_limb_t*) malloc(sizeof(mp_limb_t) * K * (n+3));
   array->sentries = (mp_limb_t*) malloc(sizeof(mp_limb_t) * K * 2);
   array->buffers = (mp_limb_t**) malloc(sizeof(mp_limb_t*) * K);
   array->mpz_buffers = (mpz_t*) malloc(sizeof(mpz_t) * K);
   
   mpz_init(array->p);
   mpz_set_ui(array->p, 1);
   mpz_mul_2exp(array->p, array->p, n * FLINT_BITS_PER_LIMB);
   mpz_add_ui(array->p, array->p, 1);

   for (unsigned long i = 0; i < K; i++)
   {
      mpz_init(array->mpz_buffers[i]);
      array->buffers[i] = array->data + (n+3)*i + 1;
      array->buffers[i][-1] = array->sentries[2*i] = random_limb();
      array->buffers[i][n+1] = array->sentries[2*i+1] = random_limb();
   }
}


/*
Deallocates a buf_array.
*/
void buf_array_clear(buf_array_t* array)
{
   for (unsigned long i = 0; i < array->K; i++)
      mpz_clear(array->mpz_buffers[i]);
   free(array->mpz_buffers);
   free(array->buffers);
   free(array->data);
   free(array->sentries);
   mpz_clear(array->p);
}


/*
Checks that the sentries haven't been overwritten, and checks that the
list of buffers still point to correct (and distinct) buffers.

Returns true on success.
*/
int buf_array_check(buf_array_t* array)
{
   unsigned long K = array->K;
   unsigned long n = array->n;
   int success = 1;

   // check all the buffers point to valid, distinct buffers
   int* used = malloc(sizeof(int) * K);
   for (unsigned long i = 0; i < K; i++)
      used[i] = 0;
   for (unsigned long i = 0; i < K; i++)
   {
      unsigned offset = array->buffers[i] - array->data - 1;
      if (offset % (n+3))
         success = 0;
      else
      {
         offset /= (n+3);
         if (offset >= K || used[offset])
            success = 0;
         else
            used[offset] = 1;
      }
   }
   free(used);

   // check sentries
   for (unsigned long i = 0; i < K; i++)
   {
      if (array->data[i*(n+3)] != array->sentries[2*i])
         success = 0;
      if (array->data[i*(n+3) + (n+2)] != array->sentries[2*i+1])
         success = 0;
   }

   return success;
}


/*
Converts buffers[index] to mpz_t format.
Result is reduced mod p.
*/
void buf_array_convert_out(buf_array_t* array, mpz_t output,
                           unsigned long index)
{
   ssfft_convert_from_raw(output, array->buffers[index], array->n);
   mpz_mod(output, output, array->p);
}


/*
Converts all raw buffers into mpz_buffers (reduced mod p).
*/
void buf_array_convert_out_all(buf_array_t* array)
{
   for (unsigned i = 0; i < array->K; i++)
      buf_array_convert_out(array, array->mpz_buffers[i], i);
}


/*
Sets buffers[index] to random data, (n + 3/4) limbs long.
Also converts the random number to mpz format, stores it in mpz_buffers[index],
reduced mod p.

It uses mpz_rrandomb to get random data with long strings of 0's and 1's.
*/
void buf_array_randomise(buf_array_t* array, unsigned long index)
{
   raw_rrandomb(array->buffers[index], array->n, FLINT_BITS_PER_LIMB/4,
                ssfft_test_randstate);
   buf_array_convert_out(array, array->mpz_buffers[index], index);
}


/*
Retrieves the (signed) high limb of the index'th buffer.
*/
mp_limb_signed_t
buf_array_get_high_limb(buf_array_t* array, unsigned long index)
{
   return array->buffers[index][array->n];
}


/*
Retrieves the number of spare bits of the high limb of the index'th buffer.
*/
unsigned long buf_array_count_spare_bits(buf_array_t* array,
                                         unsigned long index)
{
   return count_spare_bits(array->buffers[index][array->n]);
}


/* ============================================================================
   testing FFT coefficient operations

No testing for basic_add, basic_sub, basic_rotate_bits since these are
currently just wrappers for corresponding GMP mpn functions.
*/


void test_ssfft_signed_add_1()
{
   printf("testing ssfft_signed_add_1... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);
   
   mp_limb_signed_t limb;

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 1, n);

      for (unsigned long trials = 0; trials < 10000; trials++)
      {
         // initialise random data
         buf_array_randomise(&array, 0);
         raw_rrandomb((mp_limb_t*)(&limb), 0, 1, ssfft_test_randstate);
         
         // compute correct answer
         if (limb >= 0)
            mpz_add_ui(correct, array.mpz_buffers[0], limb);
         else
            mpz_sub_ui(correct, array.mpz_buffers[0], -limb);
         mpz_mod(correct, correct, array.p);
            
         // run target
         ssfft_signed_add_1(array.buffers[0], n+1, limb);

         // check buffer integrity and sentries
         if (!buf_array_check(&array))
            success = 0;
            
         // check correct answer
         buf_array_convert_out_all(&array);
         if (mpz_cmp(array.mpz_buffers[0], correct))
            success = 0;
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_fast_reduce()
{
   printf("testing basic_fast_reduce... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 1, n);

      for (unsigned long trials = 0; trials < 10000; trials++)
      {
         // initialise random data
         buf_array_randomise(&array, 0);
            
         // compute correct answer
         mpz_mod(correct, array.mpz_buffers[0], array.p);
            
         // run target
         basic_fast_reduce(array.buffers[0], n);

         // check buffer integrity and sentries
         if (!buf_array_check(&array))
            success = 0;
            
         // check correct answer
         buf_array_convert_out_all(&array);
         if (mpz_cmp(array.mpz_buffers[0], correct))
            success = 0;
            
         // check output bit guarantees
         mp_limb_signed_t high;
         high = buf_array_get_high_limb(&array, 0);
         if (high < 0 || high > 2)
            success = 0;
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_unrotate_bits()
{
   printf("testing basic_unrotate_bits... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct, temp;
   mpz_init(correct);
   mpz_init(temp);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (int inplace = 0; inplace <= 1; inplace++)
      {
         for (unsigned long s = 1; s < FLINT_BITS_PER_LIMB; s++)
         {
            for (unsigned long trials = 0; trials < 500; trials++)
            {
               // initialise random data
               buf_array_randomise(&array, 0);
               mpz_set(temp, array.mpz_buffers[0]);
               buf_array_randomise(&array, 1);
               
               // compute correct answer
               ssfft_naive_unrotate(correct, array.mpz_buffers[0],
                                    2*s, n, array.p);
               
               // run target
               basic_unrotate_bits(array.buffers[1-inplace], array.buffers[0],
                                   s, n);

               // check buffer integrity and sentries
               if (!buf_array_check(&array))
                  success = 0;

               // check correct answer
               buf_array_convert_out_all(&array);
               if (mpz_cmp(array.mpz_buffers[1-inplace], correct))
                  success = 0;

               // check input hasn't changed mod p
               if (!inplace && mpz_cmp(array.mpz_buffers[0], temp))
                  success = 0;
               
               // check input/output bit guarantees
               mp_limb_signed_t high;
               high = buf_array_get_high_limb(&array, 1-inplace);
               if (high > 1 || high < -1)
                  success = 0;
               if (!inplace)
               {
                  high = buf_array_get_high_limb(&array, 0);
                  if (high > 2 || high < 0)
                     success = 0;
               }
            }
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
   mpz_clear(temp);
}


void test_basic_rotate_limbs()
{
   printf("testing basic_rotate_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (unsigned long s = 1; s < n; s++)
      {
         for (unsigned long trials = 0; trials < 500; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            
            // compute correct answer
            ssfft_naive_rotate(correct, array.mpz_buffers[0],
                               2*FLINT_BITS_PER_LIMB*s, n, array.p);
            
            // run target
            basic_rotate_limbs(array.buffers[1], array.buffers[0], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[1], correct))
               success = 0;
            
            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 1);
            if (high > 0 || high < -2)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_unrotate_limbs()
{
   printf("testing basic_unrotate_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (unsigned long s = 1; s < n; s++)
      {
         for (unsigned long trials = 0; trials < 500; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            
            // compute correct answer
            ssfft_naive_unrotate(correct, array.mpz_buffers[0],
                                 2*FLINT_BITS_PER_LIMB*s, n, array.p);
            
            // run target
            basic_unrotate_limbs(array.buffers[1], array.buffers[0], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[1], correct))
               success = 0;
            
            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 1);
            if (high > 0 || high < -2)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_sub_rotate_limbs()
{
   printf("testing basic_sub_rotate_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 1; s < n; s++)
      {
         for (unsigned long trials = 0; trials < 500; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);
            
            // compute correct answer
            mpz_sub(correct, array.mpz_buffers[0], array.mpz_buffers[1]);
            ssfft_naive_rotate(correct, correct, 2*FLINT_BITS_PER_LIMB*s,
                               n, array.p);
            
            // run target
            basic_sub_rotate_limbs(array.buffers[2], array.buffers[0],
                                   array.buffers[1], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[2], correct))
               success = 0;
            
            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 2);
            if (high > 1 || high < -2)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_rotate_sub_limbs()
{
   printf("testing basic_rotate_sub_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long inplace = 0; inplace <= 2; inplace++)
      {
         // Need to test various combinations of buffer aliasing.
         // inplace == 0 means  c := a - 2^(Bs) b
         // inplace == 1 means  c := a - 2^(Bs) a
         // inplace == 2 means  a := a - 2^(Bs) b
         unsigned long src1, src2, dest;
         switch (inplace)
         {
         case 0: src1 = 0; src2 = 1; dest = 2; break;
         case 1: src1 = 0; src2 = 0; dest = 1; break;
         case 2: src1 = 0; src2 = 1; dest = 0; break;
         }

         for (unsigned long s = 1; s < n; s++)
         {
            for (unsigned long trials = 0; trials < 500; trials++)
            {
               // initialise random data
               buf_array_randomise(&array, 0);
               buf_array_randomise(&array, 1);
               buf_array_randomise(&array, 2);
               
               // compute correct answer
               ssfft_naive_rotate(correct, array.mpz_buffers[src2],
                                  2*FLINT_BITS_PER_LIMB*s, n, array.p);
               mpz_sub(correct, array.mpz_buffers[src1], correct);
               mpz_mod(correct, correct, array.p);
               
               // run target
               mp_limb_signed_t src1_high;
               src1_high = buf_array_get_high_limb(&array, src1);
               basic_rotate_sub_limbs(array.buffers[dest], array.buffers[src1],
                                      array.buffers[src2], s, n);
               
               // check buffer integrity and sentries
               if (!buf_array_check(&array))
                  success = 0;
               
               // check correct answer
               buf_array_convert_out_all(&array);
               if (mpz_cmp(array.mpz_buffers[dest], correct))
                  success = 0;
               
               // check output bit guarantees
               mp_limb_signed_t high;
               high = buf_array_get_high_limb(&array, dest) - src1_high;
               if (high > 1 || high < -2)
                  success = 0;
            }
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_unrotate_sub_limbs()
{
   printf("testing basic_unrotate_sub_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long inplace = 0; inplace <= 2; inplace++)
      {
         // Need to test various combinations of buffer aliasing.
         // inplace == 0 means  c := a - 2^(-Bs) b
         // inplace == 1 means  c := a - 2^(-Bs) a
         // inplace == 2 means  a := a - 2^(-Bs) b
         unsigned long src1, src2, dest;
         switch (inplace)
         {
         case 0: src1 = 0; src2 = 1; dest = 2; break;
         case 1: src1 = 0; src2 = 0; dest = 1; break;
         case 2: src1 = 0; src2 = 1; dest = 0; break;
         }

         for (unsigned long s = 1; s < n; s++)
         {
            for (unsigned long trials = 0; trials < 500; trials++)
            {
               // initialise random data
               buf_array_randomise(&array, 0);
               buf_array_randomise(&array, 1);
               buf_array_randomise(&array, 2);
               
               // compute correct answer
               ssfft_naive_unrotate(correct, array.mpz_buffers[src2],
                                    2*FLINT_BITS_PER_LIMB*s, n, array.p);
               mpz_sub(correct, array.mpz_buffers[src1], correct);
               mpz_mod(correct, correct, array.p);
               
               // run target
               mp_limb_signed_t src1_high;
               src1_high = buf_array_get_high_limb(&array, src1);
               basic_unrotate_sub_limbs(array.buffers[dest],
                                        array.buffers[src1],
                                        array.buffers[src2], s, n);
               
               // check buffer integrity and sentries
               if (!buf_array_check(&array))
                  success = 0;
               
               // check correct answer
               buf_array_convert_out_all(&array);
               if (mpz_cmp(array.mpz_buffers[dest], correct))
                  success = 0;
               
               // check output bit guarantees
               mp_limb_signed_t high;
               high = buf_array_get_high_limb(&array, dest) - src1_high;
               if (high > 2 || high < -1)
                  success = 0;
            }
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_unrotate_add_limbs()
{
   printf("testing basic_unrotate_add_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long inplace = 0; inplace <= 2; inplace++)
      {
         // Need to test various combinations of buffer aliasing.
         // inplace == 0 means  c := a + 2^(-Bs) b
         // inplace == 1 means  c := a + 2^(-Bs) a
         // inplace == 2 means  a := a + 2^(-Bs) b
         unsigned long src1, src2, dest;
         switch (inplace)
         {
         case 0: src1 = 0; src2 = 1; dest = 2; break;
         case 1: src1 = 0; src2 = 0; dest = 1; break;
         case 2: src1 = 0; src2 = 1; dest = 0; break;
         }

         for (unsigned long s = 1; s < n; s++)
         {
            for (unsigned long trials = 0; trials < 500; trials++)
            {
               // initialise random data
               buf_array_randomise(&array, 0);
               buf_array_randomise(&array, 1);
               buf_array_randomise(&array, 2);
               
               // compute correct answer
               ssfft_naive_unrotate(correct, array.mpz_buffers[src2],
                                    2*FLINT_BITS_PER_LIMB*s, n, array.p);
               mpz_add(correct, array.mpz_buffers[src1], correct);
               mpz_mod(correct, correct, array.p);
               
               // run target
               mp_limb_signed_t src1_high;
               src1_high = buf_array_get_high_limb(&array, src1);
               basic_unrotate_add_limbs(array.buffers[dest],
                                        array.buffers[src1],
                                        array.buffers[src2], s, n);
               
               // check buffer integrity and sentries
               if (!buf_array_check(&array))
                  success = 0;
               
               // check correct answer
               buf_array_convert_out_all(&array);
               if (mpz_cmp(array.mpz_buffers[dest], correct))
                  success = 0;
               
               // check output bit guarantees
               mp_limb_signed_t high;
               high = buf_array_get_high_limb(&array, dest) - src1_high;
               if (high > 1 || high < -2)
                  success = 0;
            }
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_copy()
{
   // ok it's kind of stupid to test this, but if someone changes the code
   // and accidentally overflows a buffer, that would be annoying

   printf("testing basic_copy... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (unsigned long trials = 0; trials < 500; trials++)
      {
         // initialise random data
         buf_array_randomise(&array, 0);
         buf_array_randomise(&array, 1);
            
         // compute correct "answer"
         mpz_mod(correct, array.mpz_buffers[0], array.p);
            
         // run target
         basic_copy(array.buffers[1], array.buffers[0], n);

         // check buffer integrity and sentries
         if (!buf_array_check(&array))
            success = 0;
            
         // check correct answer
         buf_array_convert_out_all(&array);
         if (mpz_cmp(array.mpz_buffers[1], correct))
            success = 0;
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_basic_negate()
{
   printf("testing basic_negate... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (int inplace = 0; inplace <= 1; inplace++)
      {
         for (unsigned long trials = 0; trials < 500; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
               
            // compute correct answer
            mpz_neg(correct, array.mpz_buffers[0]);
            mpz_mod(correct, correct, array.p);

            mp_limb_signed_t high_input;
            high_input = buf_array_get_high_limb(&array, 0);
               
            // run target
            basic_negate(array.buffers[1-inplace], array.buffers[0], n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
               
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[1-inplace], correct))
               success = 0;
               
            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 1-inplace) + high_input;
            if (high != -2)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_coeff_forward_simple_butterfly()
{
   printf("testing coeff_forward_simple_butterfly... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct0, correct1;
   mpz_init(correct0);
   mpz_init(correct1);

   for (unsigned long n = 1; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long trials = 0; trials < 500; trials++)
      {
         // initialise random data
         buf_array_randomise(&array, 0);
         buf_array_randomise(&array, 1);
         buf_array_randomise(&array, 2);
         
         unsigned long spare0 = buf_array_count_spare_bits(&array, 0);
         unsigned long spare1 = buf_array_count_spare_bits(&array, 1);
         unsigned long spare_bits = (spare0 < spare1) ? spare0 : spare1;
         
         // compute correct answer
         mpz_add(correct0, array.mpz_buffers[0], array.mpz_buffers[1]);
         mpz_sub(correct1, array.mpz_buffers[0], array.mpz_buffers[1]);
         mpz_mod(correct0, correct0, array.p);
         mpz_mod(correct1, correct1, array.p);
         
         // run target
         coeff_forward_simple_butterfly(&array.buffers[0], &array.buffers[1],
                                        &array.buffers[2], n);

         // check buffer integrity and sentries
         if (!buf_array_check(&array))
            success = 0;
         
         // check correct answer
         buf_array_convert_out_all(&array);
         if (mpz_cmp(array.mpz_buffers[0], correct0))
            success = 0;
         if (mpz_cmp(array.mpz_buffers[1], correct1))
            success = 0;

         // check output bit guarantees
         if (buf_array_count_spare_bits(&array, 0) < spare_bits - 1)
            success = 0;
         if (buf_array_count_spare_bits(&array, 1) < spare_bits - 1)
            success = 0;
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct0);
   mpz_clear(correct1);
}


void test_coeff_forward_butterfly_limbs()
{
   printf("testing coeff_forward_butterfly_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct0, correct1;
   mpz_init(correct0);
   mpz_init(correct1);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 1; s < n; s++)
      {
         for (unsigned long trials = 0; trials < 500; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);
            
            unsigned long spare0 = buf_array_count_spare_bits(&array, 0);
            unsigned long spare1 = buf_array_count_spare_bits(&array, 1);
            unsigned long spare_bits = (spare0 < spare1) ? spare0 : spare1;

            // compute correct answer
            mpz_add(correct0, array.mpz_buffers[0], array.mpz_buffers[1]);
            mpz_mod(correct0, correct0, array.p);
            mpz_sub(correct1, array.mpz_buffers[0], array.mpz_buffers[1]);
            ssfft_naive_rotate(correct1, correct1, 2*FLINT_BITS_PER_LIMB*s,
                               n, array.p);
            
            // run target
            coeff_forward_butterfly_limbs(&array.buffers[0], &array.buffers[1],
                                          &array.buffers[2], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[0], correct0))
               success = 0;
            if (mpz_cmp(array.mpz_buffers[1], correct1))
               success = 0;

            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 1);
            if (high > 1 || high < -2)
               success = 0;
            if (buf_array_count_spare_bits(&array, 0) < spare_bits - 1)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct0);
   mpz_clear(correct1);
}


void test_coeff_inverse_butterfly_limbs()
{
   printf("testing coeff_inverse_butterfly_limbs... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct0, correct1;
   mpz_init(correct0);
   mpz_init(correct1);

   for (unsigned long n = 2; n <= 6; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 1; s < n; s++)
      {
         for (unsigned long trials = 0; trials < 500; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);
            
            mp_limb_signed_t high0 = buf_array_get_high_limb(&array, 0);
            
            // compute correct answer
            ssfft_naive_unrotate(correct0, array.mpz_buffers[1],
                                 2*FLINT_BITS_PER_LIMB*s, n, array.p);
            mpz_sub(correct1, array.mpz_buffers[0], correct0);
            mpz_mod(correct1, correct1, array.p);
            mpz_add(correct0, array.mpz_buffers[0], correct0);
            mpz_mod(correct0, correct0, array.p);
            
            // run target
            coeff_inverse_butterfly_limbs(&array.buffers[0], &array.buffers[1],
                                          &array.buffers[2], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[0], correct0))
               success = 0;
            if (mpz_cmp(array.mpz_buffers[1], correct1))
               success = 0;

            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 0) - high0;
            if (high > 1 || high < -2)
               success = 0;
            high = buf_array_get_high_limb(&array, 1) - high0;
            if (high > 2 || high < -1)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct0);
   mpz_clear(correct1);
}


void test_coeff_sub_rotate_bits()
{
   printf("testing coeff_sub_rotate_bits... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
      {
         for (unsigned long trials = 0; trials < 400; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);

            unsigned long spare0 = buf_array_count_spare_bits(&array, 0);
            unsigned long spare1 = buf_array_count_spare_bits(&array, 1);
            unsigned long spare_bits = (spare0 < spare1) ? spare0 : spare1;
            
            // compute correct answer
            mpz_sub(correct, array.mpz_buffers[0], array.mpz_buffers[1]);
            ssfft_naive_rotate(correct, correct, 2*s, n, array.p);
            
            // run target
            coeff_sub_rotate_bits(&array.buffers[2], &array.buffers[0],
                                  &array.buffers[1], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[2], correct))
               success = 0;
            
            // check output bit guarantees
            if (buf_array_count_spare_bits(&array, 2) < spare_bits - 1)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_coeff_forward_butterfly_bits()
{
   printf("testing coeff_forward_butterfly_bits... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct0, correct1;
   mpz_init(correct0);
   mpz_init(correct1);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
      {
         for (unsigned long trials = 0; trials < 100; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);
            
            unsigned long spare0 = buf_array_count_spare_bits(&array, 0);
            unsigned long spare1 = buf_array_count_spare_bits(&array, 1);
            unsigned long spare_bits = (spare0 < spare1) ? spare0 : spare1;

            // compute correct answer
            mpz_add(correct0, array.mpz_buffers[0], array.mpz_buffers[1]);
            mpz_mod(correct0, correct0, array.p);
            mpz_sub(correct1, array.mpz_buffers[0], array.mpz_buffers[1]);
            ssfft_naive_rotate(correct1, correct1, 2*s, n, array.p);
            
            // run target
            coeff_forward_butterfly_bits(&array.buffers[0], &array.buffers[1],
                                         &array.buffers[2], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[0], correct0))
               success = 0;
            if (mpz_cmp(array.mpz_buffers[1], correct1))
               success = 0;

            // check output bit guarantees
            if (buf_array_count_spare_bits(&array, 0) < spare_bits - 1)
               success = 0;
            if (buf_array_count_spare_bits(&array, 1) < spare_bits - 1)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct0);
   mpz_clear(correct1);
}


void test_coeff_inverse_butterfly_bits()
{
   printf("testing coeff_inverse_butterfly_bits... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct0, correct1;
   mpz_init(correct0);
   mpz_init(correct1);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
      {
         for (unsigned long trials = 0; trials < 100; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);
            
            unsigned long spare0 = buf_array_count_spare_bits(&array, 0);
            unsigned long spare1 = buf_array_count_spare_bits(&array, 1);
            unsigned long spare_bits = (spare0 < spare1) ? spare0 : spare1;

            // compute correct answer
            ssfft_naive_unrotate(correct0, array.mpz_buffers[1],
                                 2*s, n, array.p);
            mpz_sub(correct1, array.mpz_buffers[0], correct0);
            mpz_mod(correct1, correct1, array.p);
            mpz_add(correct0, array.mpz_buffers[0], correct0);
            mpz_mod(correct0, correct0, array.p);
            
            // run target
            coeff_inverse_butterfly_bits(&array.buffers[0], &array.buffers[1],
                                         &array.buffers[2], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[0], correct0))
               success = 0;
            if (mpz_cmp(array.mpz_buffers[1], correct1))
               success = 0;

            // check output bit guarantees
            if (buf_array_count_spare_bits(&array, 0) < spare_bits - 1)
               success = 0;
            if (buf_array_count_spare_bits(&array, 1) < spare_bits - 1)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct0);
   mpz_clear(correct1);
}


void test_coeff_rotate_bits()
{
   printf("testing coeff_rotate_bits... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
      {
         for (unsigned long trials = 0; trials < 100; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            
            unsigned long spare_bits = buf_array_count_spare_bits(&array, 0);

            // compute correct answer
            ssfft_naive_rotate(correct, array.mpz_buffers[0], 2*s, n, array.p);
            
            // run target
            coeff_rotate_bits(&array.buffers[1], &array.buffers[0], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[1], correct))
               success = 0;

            // check output bit guarantees
            if (spare_bits > FLINT_BITS_PER_LIMB - 1)
               spare_bits = FLINT_BITS_PER_LIMB - 1;
            if (buf_array_count_spare_bits(&array, 1) < spare_bits)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_coeff_cross_butterfly_bits()
{
   printf("testing coeff_cross_butterfly_bits... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct0, correct1;
   mpz_init(correct0);
   mpz_init(correct1);

   for (unsigned long n = 1; n <= 3; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (unsigned long s = 0; s < n*FLINT_BITS_PER_LIMB; s++)
      {
         for (unsigned long trials = 0; trials < 100; trials++)
         {
            // initialise random data
            buf_array_randomise(&array, 0);
            buf_array_randomise(&array, 1);
            buf_array_randomise(&array, 2);
            
            unsigned long spare0 = buf_array_count_spare_bits(&array, 0);
            unsigned long spare1 = buf_array_count_spare_bits(&array, 1);
            unsigned long spare_bits = (spare0 < spare1) ? spare0 : spare1;

            // compute correct answer
            mpz_sub(correct1, array.mpz_buffers[0], array.mpz_buffers[1]);
            ssfft_naive_rotate(correct1, correct1, 2*s, n, array.p);
            mpz_add(correct0, array.mpz_buffers[0], array.mpz_buffers[0]);
            mpz_sub(correct0, correct0, array.mpz_buffers[1]);
            mpz_mod(correct0, correct0, array.p);
            
            // run target
            coeff_cross_butterfly_bits(&array.buffers[0], &array.buffers[1],
                                       &array.buffers[2], s, n);

            // check buffer integrity and sentries
            if (!buf_array_check(&array))
               success = 0;
            
            // check correct answer
            buf_array_convert_out_all(&array);
            if (mpz_cmp(array.mpz_buffers[0], correct0))
               success = 0;
            if (mpz_cmp(array.mpz_buffers[1], correct1))
               success = 0;

            // check output bit guarantees
            mp_limb_signed_t high;
            high = buf_array_get_high_limb(&array, 0);
            if (high < 0 || high > 2)
               success = 0;
            if (buf_array_count_spare_bits(&array, 1) < spare_bits - 1)
               success = 0;
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct0);
   mpz_clear(correct1);
}


void test_coeff_sqrt2_helper()
{
   printf("testing coeff_sqrt2_helper... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 16; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 2, n);

      for (unsigned long trials = 0; trials < 100; trials++)
      {
         // initialise random data
         buf_array_randomise(&array, 0);
         buf_array_randomise(&array, 1);
         
         unsigned long spare_bits = buf_array_count_spare_bits(&array, 0);

         // compute correct answer
         mpz_mul_2exp(correct, array.mpz_buffers[0], n*FLINT_BITS_PER_LIMB/2);
         mpz_sub(correct, array.mpz_buffers[0], correct);
         mpz_mod(correct, correct, array.p);
         
         // run target
         coeff_sqrt2_helper(&array.buffers[1], &array.buffers[0], n);

         // check buffer integrity and sentries
         if (!buf_array_check(&array))
            success = 0;
         
         // check correct answer
         buf_array_convert_out_all(&array);
         if (mpz_cmp(array.mpz_buffers[1], correct))
            success = 0;

         // check output bit guarantees
         if (spare_bits > FLINT_BITS_PER_LIMB - 1)
            spare_bits = FLINT_BITS_PER_LIMB - 1;
         if (buf_array_count_spare_bits(&array, 1) < spare_bits - 1)
            success = 0;
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


void test_coeff_rotate_arbitrary()
{
   printf("testing coeff_rotate_arbitrary... ");
   fflush(stdout);
   
   int success = 1;
   
   mpz_t correct;
   mpz_init(correct);

   for (unsigned long n = 1; n <= 8; n++)
   {
      buf_array_t array;
      buf_array_init(&array, 3, n);

      for (int inplace = 0; inplace <= 1; inplace++)
      {
         for (unsigned long s = 0; s < 2*n*FLINT_BITS_PER_LIMB; s++)
         {
            for (unsigned long trials = 0; trials < 30; trials++)
            {
               // initialise random data
               buf_array_randomise(&array, 0);
               buf_array_randomise(&array, 1);
               buf_array_randomise(&array, 2);
               
               unsigned long spare_bits = buf_array_count_spare_bits(&array, 0);

               // compute correct answer
               ssfft_naive_rotate(correct, array.mpz_buffers[0],
                                  s, n, array.p);
               
               // run target
               coeff_rotate_arbitrary(&array.buffers[1], &array.buffers[0],
                                 &array.buffers[2*(1-inplace)], s, n);

               // check buffer integrity and sentries
               if (!buf_array_check(&array))
                  success = 0;
               
               // check correct answer
               buf_array_convert_out_all(&array);
               if (mpz_cmp(array.mpz_buffers[1], correct))
                  success = 0;

               // check output bit guarantees
               if (spare_bits > FLINT_BITS_PER_LIMB - 1)
                  spare_bits = FLINT_BITS_PER_LIMB - 1;
               if (buf_array_count_spare_bits(&array, 1) < spare_bits - 1)
                  success = 0;
            }
         }
      }

      buf_array_clear(&array);
   }
   
   printf("%s\n", success ? "ok" : "FAIL");

   mpz_clear(correct);
}


// ============================================================================
// FFT/IFFT test case support

/*
Stores information used in testing FFTs/IFFTs.
*/
typedef struct fft_testcase_t
{
   unsigned long m, M;         // M = 2^m = transform length
   unsigned long z, g, e, n;   // transform parameters
   unsigned long ru, rU;       // describe roots of unity (power of sqrt2)
   mpz_t p;                    // schonhage-strassen prime = 2^(Bn) + 1
   
   unsigned long spare_bits;  // number of spare bits in input data

   mpz_t* input;           // input data in mpz_t format

   mp_limb_t* buf;         // buffer used by FFT
   mp_limb_t** buf_ptrs;
} fft_testcase_t;


/*
Allocates memory for and prepares random input data for an fft_testcase_t.
*/
void fft_testcase_init(fft_testcase_t* info, unsigned long m,
            unsigned long z, unsigned long g, unsigned long e,
            unsigned long ru, unsigned long rU, unsigned long n)
{
   unsigned long i;

   info->m = m;
   unsigned long M = 1 << m;
   info->M = M;
   info->z = z;
   info->g = g;
   info->e = e;
   info->ru = ru;
   info->rU = rU;
   info->n = n;
   
   mpz_init(info->p);
   mpz_set_ui(info->p, 1);
   mpz_mul_2exp(info->p, info->p, n*FLINT_BITS_PER_LIMB);
   mpz_add_ui(info->p, info->p, 1);

   info->input = (mpz_t*) malloc(M * sizeof(mpz_t));
   info->buf = (mp_limb_t*) malloc((n+1) * (M+1) * sizeof(mp_limb_t));
   info->buf_ptrs = (mp_limb_t**) malloc((M+1) * sizeof(mp_limb_t*));

   // assign buffers in some random order
   for (i = 0; i <= M; i++)
      info->buf_ptrs[i] = info->buf + i*(n+1);
   for (i = 0; i < 4*M; i++)
   {
      // swap a random pair of pointers
      unsigned long x = gmp_urandomm_ui(ssfft_test_randstate, M + 1);
      unsigned long y = gmp_urandomm_ui(ssfft_test_randstate, M + 1);
      mp_limb_t* temp = info->buf_ptrs[x];
      info->buf_ptrs[x] = info->buf_ptrs[y];
      info->buf_ptrs[y] = temp;
   }
      
   // put random data into the buffers
   for (i = 0; i <= M; i++)
   {
      raw_rrandomb(info->buf_ptrs[i], n, 3*FLINT_BITS_PER_LIMB/4,
                   ssfft_test_randstate);
   }
   
   // about half the time, replace top limb with *uniformly* random
   // perturbation, to get better testing for overflow limb bit counts
   if (gmp_urandomb_ui(ssfft_test_randstate, 1))
      for (i = 0; i <= M; i++)
         info->buf_ptrs[i][n] = gmp_urandomm_ui(ssfft_test_randstate, 15) - 8;

   // compute number of spare bits in input
   // (the extra -2 is because sometimes the overflow limb picks up an extra
   // bit at the bottom, regardless of bit buildup)
   info->spare_bits = FLINT_BITS_PER_LIMB - 2;
   for (i = 0; i < M; i++)
   {
      unsigned long temp = count_spare_bits(info->buf_ptrs[i][n]);
      if (temp < info->spare_bits)
         info->spare_bits = temp;
   }
   
   // convert random input to mpz_t format
   for (i = 0; i < M; i++)
   {
      mpz_init(info->input[i]);
      ssfft_convert_from_raw(info->input[i], info->buf_ptrs[i], n);
      mpz_mod(info->input[i], info->input[i], info->p);
   }
}


void fft_testcase_clear(fft_testcase_t* info)
{
   free(info->buf_ptrs);
   free(info->buf);
   for (unsigned long i = 0; i < info->M; i++)
      mpz_clear(info->input[i]);
   free(info->input);
   mpz_clear(info->p);
}


/*
Checks that the buf_ptrs each point to a correct buffer, and that they
are all distinct.

Returns 1 on success.
*/
int check_buf_ptrs(mp_limb_t* buf, mp_limb_t** buf_ptrs,
                   unsigned long M, unsigned long n)
{
   int* used = (int*) malloc((M+1) * sizeof(int));
   int success = 1;
   unsigned long i, diff;

   for (i = 0; i <= M; i++)
      used[i] = 0;

   for (i = 0; i <= M; i++)
   {
      diff = buf_ptrs[i] - buf;
      if (diff % (n+1))
         success = 0;
      diff /= (n+1);
      if (diff > M || used[i])
         success = 0;
      used[i] = 1;
   }
   
   free(used);
   return success;
}


/*
Performs an inplace, nontruncated FFT using simple mpz arithmetic.
*/
void ssfft_naive_fft(mpz_t* data, unsigned long m,
            unsigned long ru, unsigned long rU, unsigned long n, mpz_t p)
{
   unsigned long M = 1 << m;
   unsigned long layer, half, start, i;
   mpz_t temp;
   mpz_init(temp);
   
   for (layer = 0; layer < m; layer++)
   {
      half = 1 << (m - layer - 1);
      for (start = 0; start < M; start += 2*half)
      {
         for (i = 0; i < half; i++)
         {
            mpz_t* x = data + start + i;
            mpz_t* y = data + start + half + i;
            mpz_sub(temp, *x, *y);
            ssfft_naive_rotate(temp, temp, (ru + i*rU) << layer, n, p);
            mpz_add(*x, *x, *y);
            if (mpz_cmp(*x, p) >= 0)
               mpz_sub(*x, *x, p);
            mpz_set(*y, temp);
         }
      }
   }

   mpz_clear(temp);
}


/*
Checks that an FFT was performed correctly:
(1) checks validity of output buffer pointers
(2) checks number of bits used in each overflow limb
(3) performs naive FFT and compares result against data in buffer

Returns 1 on success.
*/
int fft_testcase_check_forward(fft_testcase_t* info)
{
   unsigned long m = info->m;
   unsigned long M = info->M;
   unsigned long n = info->n;
   unsigned long g = info->g;
   unsigned long z = info->z;
   unsigned long ru = info->ru;
   unsigned long rU = info->rU;
   unsigned long layer, i, start, half;

   // check buffer pointers
   if (!check_buf_ptrs(info->buf, info->buf_ptrs, M, n))
      return 0;
  
   // check overflow limbs
   for (i = 0; i < g; i++)
      if (count_spare_bits(info->buf_ptrs[i][n]) < info->spare_bits - m)
         return 0;
   
   // perform naive FFT
   for (i = z; i < M; i++)
      mpz_set_ui(info->input[i], 0);
   ssfft_naive_fft(info->input, m, ru, rU, n, info->p);
   
   // compare naive FFT output to raw output
   mpz_t temp;
   mpz_init(temp);
   int success = 1;
   for (i = 0; i < g; i++)
   {
      ssfft_convert_from_raw(temp, info->buf_ptrs[i], n);
      mpz_mod(temp, temp, info->p);
      if (mpz_cmp(temp, info->input[i]))
         success = 0;
   }
   mpz_clear(temp);
   
   return success;
}


/*
Checks that an IFFT was performed correctly:
(1) checks validity of output buffer pointers
(2) checks number of bits used in each overflow limb
(3) performs naive FFT on raw output and compares result against original data

Returns 1 on success.
*/
int fft_testcase_check_inverse(fft_testcase_t* info)
{
   unsigned long m = info->m;
   unsigned long M = info->M;
   unsigned long n = info->n;
   unsigned long z = info->z;
   unsigned long g = info->g;
   unsigned long e = info->e;
   unsigned long ru = info->ru;
   unsigned long rU = info->rU;
   unsigned long layer, i, start, half;

   // check buffer pointers
   if (!check_buf_ptrs(info->buf, info->buf_ptrs, M, n))
      return 0;

   // check overflow limbs
   for (i = 0; i < g + e; i++)
      if (count_spare_bits(info->buf_ptrs[i][n]) < info->spare_bits - m)
         return 0;

   // construct new array of mpz_t's with the untransformed coefficients
   // (this is pieced together from the raw IFFT output and the original
   // input data), remember this is scaled by 2^m
   mpz_t* array = malloc(M * sizeof(mpz_t));
   for (i = 0; i < M; i++)
      mpz_init(array[i]);
      
   for (i = 0; i < g; i++)
      ssfft_convert_from_raw(array[i], info->buf_ptrs[i], n);
   for (i = g; i < z; i++)
      mpz_set(array[i], info->input[i]);
   for (i = z; i < M; i++)
      mpz_set_ui(array[i], 0);
      
   // run naive FFT on the data
   ssfft_naive_fft(array, m, ru, rU, n, info->p);

   int success = 1;
   
   // compare first g coefficients to original input
   for (i = 0; i < g; i++)
   {
      ssfft_naive_unrotate(array[i], array[i], 2*m, n, info->p);
      if (mpz_cmp(array[i], info->input[i]))
         success = 0;
   }
   
   // if e is set, compare g-th coefficient to IFFT output
   if (e)
   {
      mpz_t temp;
      mpz_init(temp);
      ssfft_convert_from_raw(temp, info->buf_ptrs[g], n);
      mpz_mod(temp, temp, info->p);
      ssfft_naive_unrotate(array[g], array[g], 2*m, n, info->p);
      if (mpz_cmp(temp, array[g]))
         success = 0;
      mpz_clear(temp);
   }
   
   for (i = 0; i < M; i++)
      mpz_clear(array[i]);
   free(array);

   return success;
}


// ============================================================================
// Forward FFT test routines

void test_ssfft_fft_size4_z4g4_limbs()
{
   printf("testing ssfft_fft_size4_z4g4_limbs... ");
   fflush(stdout);

   unsigned long n, trial, ru, rU;
   int success = 1;

   for (n = 2; n < 32; n += 2)
   {
      // fourth root of unity is sqrt2^(nB)
      rU = n * FLINT_BITS_PER_LIMB;
      
      for (ru = 0; ru < rU; ru += 2*FLINT_BITS_PER_LIMB)
      {
         for (trial = 0; trial < 60; trial++)
         {
            fft_testcase_t info;
            fft_testcase_init(&info, 2, 4, 4, 0, ru, rU, n);

            ssfft_fft_size4_z4g4_limbs(info.buf_ptrs, 1,
                  ru / (2*FLINT_BITS_PER_LIMB), rU / (2*FLINT_BITS_PER_LIMB),
                  n, info.buf_ptrs + 4);

            if (!fft_testcase_check_forward(&info))
               success = 0;
            fft_testcase_clear(&info);
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_fft_size8_z8g8_limbs()
{
   printf("testing ssfft_fft_size8_z8g8_limbs... ");
   fflush(stdout);

   unsigned long n, trial, ru, rU;
   int success = 1;

   for (n = 4; n < 64; n += 4)
   {
      // eighth root of unity is sqrt2^(nB/2)
      rU = n * FLINT_BITS_PER_LIMB / 2;
      
      for (ru = 0; ru < rU; ru += 2*FLINT_BITS_PER_LIMB)
      {
         for (trial = 0; trial < 60; trial++)
         {
            fft_testcase_t info;
            fft_testcase_init(&info, 3, 8, 8, 0, ru, rU, n);

            ssfft_fft_size8_z8g8_limbs(info.buf_ptrs, 1,
                  ru / (2*FLINT_BITS_PER_LIMB), rU / (2*FLINT_BITS_PER_LIMB),
                  n, info.buf_ptrs + 8);

            if (!fft_testcase_check_forward(&info))
               success = 0;
            fft_testcase_clear(&info);
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_fft_size4_bits()
{
   printf("testing ssfft_fft_size4_bits... ");
   fflush(stdout);

   int success = 1;

   unsigned long z, g, n, trial, ru, rU;
   
   for (z = 1; z <= 4; z++)
   {
      for (g = 1; g <= 4; g++)
      {
         for (n = 1; n <= 8; n++)
         {
            // fourth root of unity is sqrt2^(nB)
            rU = n * FLINT_BITS_PER_LIMB;

            for (ru = 0; ru < rU; ru += 2)
            {
               for (trial = 0; trial < 10; trial++)
               {
                  fft_testcase_t info;
                  fft_testcase_init(&info, 2, z, g, 0, ru, rU, n);

                  ssfft_fft_size4_bits(info.buf_ptrs, 1, z, g,
                                       ru / 2, rU / 2, n, info.buf_ptrs + 4);

                  if (!fft_testcase_check_forward(&info))
                     success = 0;
                  
                  fft_testcase_clear(&info);
               }
            }
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_fft_size8_bits()
{
   printf("testing ssfft_fft_size8_bits... ");
   fflush(stdout);

   int success = 1;

   unsigned long z, g, n, trial, ru, rU;
   
   for (z = 1; z <= 8; z++)
   {
      for (g = 1; g <= 8; g++)
      {
         for (n = 1; n <= 8; n++)
         {
            // eighth root of unity is sqrt2^(nB/2)
            rU = n * FLINT_BITS_PER_LIMB / 2;

            for (ru = 0; ru < rU; ru += 2)
            {
               for (trial = 0; trial < 10; trial++)
               {
                  fft_testcase_t info;
                  fft_testcase_init(&info, 3, z, g, 0, ru, rU, n);

                  ssfft_fft_size8_bits(info.buf_ptrs, 1, z, g,
                                       ru / 2, rU / 2, n, info.buf_ptrs + 8);

                  if (!fft_testcase_check_forward(&info))
                     success = 0;
                  
                  fft_testcase_clear(&info);
               }
            }
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_fft()
{
   printf("testing ssfft_fft... ");
   fflush(stdout);
   
   // For 1 <= m <= 7 (i.e. transform lengths 2 through 128) we try every
   // legal value of z and g, with n = 16, with ru = 0 and ru = 3, to give the
   // FFT factoring and truncation code a workout.
   
   // Then for 1 <= m <= 12, we try a random selection of ru, z, g.
   
   // (Ideally we'd like to try *everything*, but it just takes too long.)

   int success = 1;

   unsigned long m, M, z, g, n, trial, ru, rU;

   for (m = 1; m <= 7; m++)
   {
      M = 1 << m;
      n = 16;
      // Mth root of unity is sqrt2^(4nB/M)
      rU = 4 * n * FLINT_BITS_PER_LIMB / M;

      for (z = 1; z <= M; z++)
      {
         for (g = 1; g <= M; g++)
         {
            for (ru = 0; ru < rU && ru <= 3; ru += 3)
            {
               fft_testcase_t info;
               fft_testcase_init(&info, m, z, g, 0, ru, rU, n);

               ssfft_fft(info.buf_ptrs, 1, m, z, g, ru, rU,
                         n, info.buf_ptrs + M);

               if (!fft_testcase_check_forward(&info))
                  success = 0;
               
               fft_testcase_clear(&info);
            }
         }
      }
   }

   for (m = 1; m <= 12; m++)
   {
      M = 1 << m;
      n = 64;
      // Mth root of unity is sqrt2^(4nB/M)
      rU = 4 * n * FLINT_BITS_PER_LIMB / M;
      
      for (trial = 0; trial < 5000/M || trial < 25; trial++)
      {
         z = gmp_urandomm_ui(ssfft_test_randstate, M-1) + 1;
         g = gmp_urandomm_ui(ssfft_test_randstate, M-1) + 1;
         ru = gmp_urandomm_ui(ssfft_test_randstate, rU);

         fft_testcase_t info;
         fft_testcase_init(&info, m, z, g, 0, ru, rU, n);

         ssfft_fft(info.buf_ptrs, 1, m, z, g, ru, rU,
                   n, info.buf_ptrs + M);

         if (!fft_testcase_check_forward(&info))
            success = 0;
         
         fft_testcase_clear(&info);
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_fft_threaded()
{
   printf("testing ssfft_fft_threaded... ");
   fflush(stdout);
   
   // For 1 <= m <= 7 (i.e. transform lengths 2 through 128) we try every
   // legal value of z and g, with n = 16, with ru = 0 and ru = 3, to give the
   // FFT factoring and truncation code a workout.
   
   // Then for 1 <= m <= 12, we try a random selection of ru, z, g.
   
   // (Ideally we'd like to try *everything*, but it just takes too long.)

   int success = 1;

   unsigned long m, M, z, g, n, trial, ru, rU;

   for (m = 2; m <= 7; m++)
   {
      M = 1 << m;
      n = 16;
      // Mth root of unity is sqrt2^(4nB/M)
      rU = 4 * n * FLINT_BITS_PER_LIMB / M;

      for (z = 1; z <= M; z++)
      {
         for (g = 1; g <= M; g++)
         {
            for (ru = 0; ru < rU && ru <= 3; ru += 3)
            {
               fft_testcase_t info;
               fft_testcase_init(&info, m, z, g, 0, ru, rU, n);

               ssfft_fft_threaded(info.buf_ptrs, 1, m, z, g, ru, rU, n);

               if (!fft_testcase_check_forward(&info))
                  success = 0;
               
               fft_testcase_clear(&info);
            }
         }
      }
   }

   for (m = 2; m <= 12; m++)
   {
      M = 1 << m;
      n = 64;
      // Mth root of unity is sqrt2^(4nB/M)
      rU = 4 * n * FLINT_BITS_PER_LIMB / M;
      
      for (trial = 0; trial < 5000/M || trial < 25; trial++)
      {
         z = gmp_urandomm_ui(ssfft_test_randstate, M-1) + 1;
         g = gmp_urandomm_ui(ssfft_test_randstate, M-1) + 1;
         ru = gmp_urandomm_ui(ssfft_test_randstate, rU);

         fft_testcase_t info;
         fft_testcase_init(&info, m, z, g, 0, ru, rU, n);

         ssfft_fft_threaded(info.buf_ptrs, 1, m, z, g, ru, rU, n);

         if (!fft_testcase_check_forward(&info))
            success = 0;
         
         fft_testcase_clear(&info);
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}



// ============================================================================
// Inverse FFT test routines


void test_ssfft_ifft_size4_z4g4_limbs()
{
   printf("testing ssfft_ifft_size4_z4g4_limbs... ");
   fflush(stdout);

   unsigned long n, trial, ru, rU;
   int success = 1;

   for (n = 2; n < 32; n += 2)
   {
      // fourth root of unity is sqrt2^(nB)
      rU = n * FLINT_BITS_PER_LIMB;
      
      for (ru = 0; ru < rU; ru += 2*FLINT_BITS_PER_LIMB)
      {
         for (trial = 0; trial < 60; trial++)
         {
            fft_testcase_t info;
            fft_testcase_init(&info, 2, 4, 4, 0, ru, rU, n);

            ssfft_ifft_size4_z4g4_limbs(info.buf_ptrs, 1,
                  ru / (2*FLINT_BITS_PER_LIMB), rU / (2*FLINT_BITS_PER_LIMB),
                  n, info.buf_ptrs + 4);

            if (!fft_testcase_check_inverse(&info))
               success = 0;
            fft_testcase_clear(&info);
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_ifft_size8_z8g8_limbs()
{
   printf("testing ssfft_ifft_size8_z8g8_limbs... ");
   fflush(stdout);

   unsigned long n, trial, ru, rU;
   int success = 1;

   for (n = 4; n < 64; n += 4)
   {
      // eighth root of unity is sqrt2^(nB/2)
      rU = n * FLINT_BITS_PER_LIMB / 2;
      
      for (ru = 0; ru < rU; ru += 2*FLINT_BITS_PER_LIMB)
      {
         for (trial = 0; trial < 60; trial++)
         {
            fft_testcase_t info;
            fft_testcase_init(&info, 3, 8, 8, 0, ru, rU, n);

            ssfft_ifft_size8_z8g8_limbs(info.buf_ptrs, 1,
                  ru / (2*FLINT_BITS_PER_LIMB), rU / (2*FLINT_BITS_PER_LIMB),
                  n, info.buf_ptrs + 8);

            if (!fft_testcase_check_inverse(&info))
               success = 0;
            fft_testcase_clear(&info);
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_ifft_size4_bits()
{
   printf("testing ssfft_ifft_size4_bits... ");
   fflush(stdout);

   int success = 1;

   unsigned long z, g, e, n, trial, ru, rU;
   
   for (z = 1; z <= 4; z++)
   {
      for (g = 0; g <= z; g++)
      {
         for (e = 0; e <= 1; e++)
         {
            // illegal cases:
            if (g == 0 && e == 0)
               continue;
            if (g == 4 && e == 1)
               continue;

            for (n = 1; n <= 8; n++)
            {
               // fourth root of unity is sqrt2^(nB)
               rU = n * FLINT_BITS_PER_LIMB;

               for (ru = 0; ru < rU; ru += 2)
               {
                  for (trial = 0; trial < 10; trial++)
                  {
                     fft_testcase_t info;
                     fft_testcase_init(&info, 2, z, g, e, ru, rU, n);

                     ssfft_ifft_size4_bits(info.buf_ptrs, 1, z, g, e,
                                    ru / 2, rU / 2, n, info.buf_ptrs + 4);

                     if (!fft_testcase_check_inverse(&info))
                        success = 0;
                     
                     fft_testcase_clear(&info);
                  }
               }
            }
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_ifft_size8_bits()
{
   printf("testing ssfft_ifft_size8_bits... ");
   fflush(stdout);

   int success = 1;

   unsigned long z, g, e, n, trial, ru, rU;
   
   for (z = 1; z <= 8; z++)
   {
      for (g = 0; g <= z; g++)
      {
         for (e = 0; e <= 1; e++)
         {
            // illegal cases:
            if (g == 0 && e == 0)
               continue;
            if (g == 8 && e == 1)
               continue;

            for (n = 1; n <= 8; n++)
            {
               // eighth root of unity is sqrt2^(nB/2)
               rU = n * FLINT_BITS_PER_LIMB / 2;

               for (ru = 0; ru < rU; ru += 2)
               {
                  for (trial = 0; trial < 10; trial++)
                  {
                     fft_testcase_t info;
                     fft_testcase_init(&info, 3, z, g, e, ru, rU, n);

                     ssfft_ifft_size8_bits(info.buf_ptrs, 1, z, g, e,
                                    ru / 2, rU / 2, n, info.buf_ptrs + 8);

                     if (!fft_testcase_check_inverse(&info))
                        success = 0;
                     
                     fft_testcase_clear(&info);
                  }
               }
            }
         }
      }
   }
   
   printf("%s\n", success ? "ok" : "FAIL");
}


void test_ssfft_ifft()
{
   printf("testing ssfft_ifft... ");
   fflush(stdout);
   
   // For 1 <= m <= 7 (i.e. transform lengths 2 through 128) we try every
   // legal value of z, g and e, with n = 16, with ru = 0 and ru = 3, to give
   // the FFT factoring and truncation code a workout.
   
   // Then for 1 <= m <= 12, we try a random selection of ru, z, g and e.
   
   // (Ideally we'd like to try *everything*, but it just takes too long.)

   int success = 1;

   unsigned long m, M, z, g, e, n, trial, ru, rU;

   for (m = 1; m <= 7; m++)
   {
      M = 1 << m;
      n = 16;
      // Mth root of unity is sqrt2^(4nB/M)
      rU = 4 * n * FLINT_BITS_PER_LIMB / M;

      for (z = 1; z <= M; z++)
      {
         for (g = 0; g <= z; g++)
         {
            for (e = 0; e <= 1; e++)
            {
               // illegal cases:
               if (g == 0 && e == 0)
                  continue;
               if (g == M && e == 1)
                  continue;

               for (ru = 0; ru < rU && ru <= 3; ru += 3)
               {
                  fft_testcase_t info;
                  fft_testcase_init(&info, m, z, g, e, ru, rU, n);

                  ssfft_ifft(info.buf_ptrs, 1, m, z, g, e, ru, rU,
                            n, info.buf_ptrs + M);

                  if (!fft_testcase_check_inverse(&info))
                     success = 0;
                  
                  fft_testcase_clear(&info);
               }
            }
         }
      }
   }

   for (m = 1; m <= 12; m++)
   {
      M = 1 << m;
      n = 64;
      // Mth root of unity is sqrt2^(4nB/M)
      rU = 4 * n * FLINT_BITS_PER_LIMB / M;
      
      for (trial = 0; trial < 5000/M || trial < 25; trial++)
      {
         z = gmp_urandomm_ui(ssfft_test_randstate, M-1) + 1;
         g = gmp_urandomm_ui(ssfft_test_randstate, z+1);
         if (g == 0)
            e = 1;
         else if (g == M)
            e = 0;
         else
            e = gmp_urandomm_ui(ssfft_test_randstate, 2);

         ru = gmp_urandomm_ui(ssfft_test_randstate, rU);

         fft_testcase_t info;
         fft_testcase_init(&info, m, z, g, e, ru, rU, n);

         ssfft_ifft(info.buf_ptrs, 1, m, z, g, e, ru, rU,
                   n, info.buf_ptrs + M);

         if (!fft_testcase_check_inverse(&info))
            success = 0;
         
         fft_testcase_clear(&info);
      }
   }

   printf("%s\n", success ? "ok" : "FAIL");
}


// ============================================================================
// test function for naive_KS_mul, just here temporarily until it finds a home

void test_naive_KS_mul()
{
   printf("testing naive_KS_mul... ");
   fflush(stdout);
   
   int success = 1;
   
   unsigned long max_degree = 10;
   unsigned long max_bitsize = 10;
   mpz_t* poly1 = malloc((max_degree + 1) * sizeof(mpz_t));
   mpz_t* poly2 = malloc((max_degree + 1) * sizeof(mpz_t));
   mpz_t* poly3 = malloc((max_degree*2 + 1) * sizeof(mpz_t));
   mpz_t* poly4 = malloc((max_degree*2 + 1) * sizeof(mpz_t));
   for (unsigned long i = 0; i <= max_degree; i++)
   {
      mpz_init(poly1[i]);
      mpz_init(poly2[i]);
   }
   for (unsigned long i = 0; i <= 2*max_degree; i++)
   {
      mpz_init(poly3[i]);
      mpz_init(poly4[i]);
   }

   for (unsigned long degree1 = 1; degree1 <= max_degree; degree1++)
      for (unsigned long degree2 = 1; degree2 <= max_degree; degree2++)
         for (unsigned long bitsize1 = 1; bitsize1 <= max_bitsize; bitsize1++)
            for (unsigned long bitsize2 = 1; bitsize2 <= max_bitsize; bitsize2++)
               for (unsigned long trial = 0; trial < 10; trial++)
               {
                  // generate random polys
                  for (unsigned long i = 0; i <= degree1; i++)
                  {
                     mpz_rrandomb(poly1[i], ssfft_test_randstate,
                           gmp_urandomm_ui(ssfft_test_randstate, bitsize1+1));
                     if (gmp_urandomb_ui(ssfft_test_randstate, 1))
                        mpz_neg(poly1[i], poly1[i]);
                  }
                  for (unsigned long i = 0; i <= degree2; i++)
                  {
                     mpz_rrandomb(poly2[i], ssfft_test_randstate,
                           gmp_urandomm_ui(ssfft_test_randstate, bitsize2+1));
                     if (gmp_urandomb_ui(ssfft_test_randstate, 1))
                        mpz_neg(poly2[i], poly2[i]);
                  }
                  
                  // compute naive product
                  unsigned long degree3 = degree1 + degree2;
                  for (unsigned long i = 0; i <= degree3; i++)
                     mpz_set_ui(poly3[i], 0);
                  for (unsigned long i = 0; i <= degree1; i++)
                     for (unsigned long j = 0; j <= degree2; j++)
                        mpz_addmul(poly3[i+j], poly1[i], poly2[j]);
                  
                  // try with naive_KS_mul
                  naive_KS_mul(poly4, poly1, degree1 + 1, poly2, degree2 + 1);
                  
                  // compare results
                  for (unsigned long i = 0; i <= degree3; i++)
                     success = success && !mpz_cmp(poly3[i], poly4[i]);
               }

   for (unsigned long i = 0; i <= 2*max_degree; i++)
   {
      mpz_clear(poly3[i]);
      mpz_clear(poly4[i]);
   }
   for (unsigned long i = 0; i <= max_degree; i++)
   {
      mpz_clear(poly1[i]);
      mpz_clear(poly2[i]);
   }
   free(poly2);
   free(poly1);

   printf("%s\n", success ? "ok" : "FAIL");
}


// ============================================================================
// Main function (calls all test functions)

int main()
{
   gmp_randinit_default(ssfft_test_randstate);

   test_ssfft_fft_threaded();

   test_naive_KS_mul();

   test_ssfft_signed_add_1();
   test_basic_fast_reduce();
   test_basic_unrotate_bits();
   test_basic_rotate_limbs();
   test_basic_unrotate_limbs();
   test_basic_sub_rotate_limbs();
   test_basic_rotate_sub_limbs();
   test_basic_unrotate_sub_limbs();
   test_basic_unrotate_add_limbs();
   test_basic_copy();
   test_basic_negate();
   test_coeff_forward_simple_butterfly();
   test_coeff_forward_butterfly_limbs();
   test_coeff_inverse_butterfly_limbs();
   test_coeff_sub_rotate_bits();
   test_coeff_forward_butterfly_bits();
   test_coeff_inverse_butterfly_bits();
   test_coeff_rotate_bits();
   test_coeff_cross_butterfly_bits();
   test_coeff_sqrt2_helper();
   test_coeff_rotate_arbitrary();

   test_ssfft_fft_size4_z4g4_limbs();
   test_ssfft_fft_size8_z8g8_limbs();
   test_ssfft_fft_size4_bits();
   test_ssfft_fft_size8_bits();
   test_ssfft_fft();

   test_ssfft_ifft_size4_z4g4_limbs();
   test_ssfft_ifft_size8_z8g8_limbs();
   test_ssfft_ifft_size4_bits();
   test_ssfft_ifft_size8_bits();
   test_ssfft_ifft();

   return 0;
}
