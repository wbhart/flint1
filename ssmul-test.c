/*****************************************************************************

 ssmul-test.c
 test routines for ssmul.c

 Copyright (C) 2006, David Harvey

****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gmp.h>
#include "flint.h"
#include "ssmul.h"
#include "Zvec.h"
#include "Zvec.c"
#include "mpn_extras.h"


gmp_randstate_t test_state;


// ============ some auxiliary functions used by test routines ============


mp_limb_t random_limb()
{
   return gmp_urandomb_ui(test_state, FLINT_BITS_PER_LIMB);
}

// Sets p to 2^(nB) + 1.
// Assumes p is already init'd.
void make_p(mpz_t p, unsigned long n)
{
   mpz_set_ui(p, 1);
   mpz_mul_2exp(p, p, n * FLINT_BITS_PER_LIMB);
   mpz_add_ui(p, p, 1);
}

// rotates left by s bits mod p, the naive way
// result is inplace
void naive_rotate_left(mpz_t x, unsigned long s, mpz_t p)
{
   mpz_mul_2exp(x, x, s);
   mpz_mod(x, x, p);
}


// rotates right by s bits mod p, the naive way
// result in inplace
void naive_rotate_right(mpz_t x, unsigned long s, mpz_t p)
{
   mpz_t temp;
   mpz_init(temp);
   mpz_set_ui(temp, 1);
   mpz_mul_2exp(temp, temp, s);
   mpz_invert(temp, temp, p);
   mpz_mul(x, x, temp);
   mpz_mod(x, x, p);
   mpz_clear(temp);
}


/*
This is like convert_raw_to_mpz, except it allows negative input, and it allows
input in the high limb as well.

It's not actually used in ssmul.c, it's only used in the test module.
It's not supposed to be particularly fast.
 */
void convert_signed_raw_to_mpz(mpz_t output, mp_limb_t* input, unsigned long n)
{
   if ((mp_limb_signed_t) input[n] >= 0)
      mpz_import(output, n + 1, -1, sizeof(mp_limb_t), 0, 0, input);
   else
   {
      // input is negative
      negate_limbs(input, input, n + 1);
      mpz_import(output, n + 1, -1, sizeof(mp_limb_t), 0, 0, input);
      mpz_neg(output, output);
      negate_limbs(input, input, n + 1);
   }
}


// ==================== the actual test routines ====================

/*
Tests negate_limbs on random input and a bunch of boundary cases.

todo: this should be moved to an mpn_extras-test.c file...?
 */
void test_negate_limbs()
{
   printf("   testing negate_limbs...\n");

   mp_limb_signed_t cases[][3] = {
      {0, 0, 0},
      {1, 0, 0},
      {2, 0, 0},
      {0, 1, 0},
      {0, 0, 1},
      {-1, 0, 0},
      {-2, 0, 0},
      {1, 1, 1},
      {-1, -1, -1},
      {-2, -1, -1},
      {0, -1, -1}
   };

   for (unsigned i = 0; i < 1000; i++)
   {
      mp_limb_t input[3];
      mp_limb_t output[5];

      mp_limb_t sentry = random_limb();
      output[0] = output[4] = sentry;

      for (unsigned j = 0; j < 3; j++)
         if (i < sizeof(cases) / sizeof(cases[0]))
            // one of the above test cases
            input[j] = cases[i][j];
         else
            // random test case
            input[j] = random_limb();

      // input + negate_limbs(input) should be zero
      negate_limbs(output + 1, input, 3);
      mpn_add_n(output + 1, output + 1, input, 3);
      for (unsigned j = 0; j < 3; j++)
         assert(output[j + 1] == 0);

      assert(output[0] == sentry);
      assert(output[4] == sentry);
   }
}

typedef struct signed_add_test_case
{
   mp_limb_signed_t input[3];
   mp_limb_signed_t limb;
   mp_limb_signed_t output[3];
} signed_add_test_case;

/*
todo: this should be moved to an mpn_extras-test.c file...?
*/
void test_signed_add_1()
{
   printf("   testing signed_add_1...\n");

   signed_add_test_case cases[] = {
      {{0, 0, 0},         0,   {0, 0, 0}},
      {{1, 0, 0},         0,   {1, 0, 0}},
      {{2, 0, 0},         0,   {2, 0, 0}},
      {{-1, 0, 0},        0,   {-1, 0, 0}},
      {{-2, 0, 0},        0,   {-2, 0, 0}},
      {{-1, -1, -1},      0,   {-1, -1, -1}},
      {{-1, -2, -1},      0,   {-1, -2, -1}},
      {{-2, -1, -1},      0,   {-2, -1, -1}},

      {{0, 0, 0},         1,   {1, 0, 0}},
      {{1, 0, 0},         1,   {2, 0, 0}},
      {{2, 0, 0},         1,   {3, 0, 0}},
      {{-1, 0, 0},        1,   {0, 1, 0}},
      {{-2, 0, 0},        1,   {-1, 0, 0}},
      {{-1, -1, -1},      1,   {0, 0, 0}},
      {{-1, -2, -1},      1,   {0, -1, -1}},
      {{-2, -1, -1},      1,   {-1, -1, -1}},

      {{0, 0, 0},         -1,  {-1, -1, -1}},
      {{1, 0, 0},         -1,  {0, 0, 0}},
      {{2, 0, 0},         -1,  {1, 0, 0}},
      {{-1, 0, 0},        -1,  {-2, 0, 0}},
      {{-2, 0, 0},        -1,  {-3, 0, 0}},
      {{-1, -1, -1},      -1,  {-2, -1, -1}},
      {{-1, -2, -1},      -1,  {-2, -2, -1}},
      {{-2, -1, -1},      -1,  {-3, -1, -1}},
   };

   for (unsigned i = 0; i < sizeof(cases) / sizeof(cases[0]); i++)
   {
      mp_limb_t output[5];

      mp_limb_t sentry = random_limb();
      output[0] = output[4] = sentry;

      for (unsigned j = 0; j < 3; j++)
         output[j+1] = (mp_limb_t) cases[i].input[j];

      signed_add_1(output + 1, 3, cases[i].limb);

      for (unsigned j = 0; j < 3; j++)
         assert(output[j + 1] == cases[i].output[j]);

      assert(output[0] == sentry);
      assert(output[4] == sentry);
   }
}

/*void test_ssmul_convert_in_bytes(unsigned long num_trials, unsigned long length, 
          unsigned long coeff_bits, unsigned long pack_bytes)
{
   mpz_t* data1 = (mpz_t*) malloc(length * sizeof(mpz_t));
   mpz_t* data2 = (mpz_t*) malloc(length * sizeof(mpz_t));
   for (unsigned long i = 0; i < length; i++)
   {
      mpz_init(data1[i]);
      mpz_init(data2[i]);
   }
   gmp_randinit_default(test_state);
   
   unsigned long num_limbs = ((length*pack_bytes)*8)/FLINT_BITS_PER_LIMB + 2;
   
   mp_limb_t* array = (mp_limb_t*) limb_alloc(num_limbs,0);
   
   for (unsigned long i = 0; i < num_trials; i++)
   {
      // make up random polys
      for (unsigned long j = 0; j < length; j++)
      {
         mpz_urandomb(data1[j], test_state, coeff_bits);
         if (gmp_urandomm_ui(test_state, 2))
            mpz_neg(data1[j], data1[j]);
         //gmp_printf("%Zx, ",data1[j]);
      }   
      //printf("\n\n");
       
      clear_limbs(array, num_limbs);
      
      ssmul_convert_in_bytes(array, data1, length, pack_bytes, num_limbs);
      
      //for (unsigned long j = 0; j < num_limbs; j++)
      //{
         //printf("%lx ", array[j]);
      //}
      //printf("\n\n");
      
      for (unsigned long i = 0; i < length; i++)
      {
         mpz_set_ui(data2[i],0);
      }
   
      ssmul_convert_out_bytes(data2,array,length,pack_bytes);

      for (unsigned long j = 0; j < length; j++)
      {
         //gmp_printf("%Zx, ",data2[j]);
         assert(!mpz_cmp(data1[j],data2[j]));
      }   
      //printf("\n\n\n\n");
     
   }

   for (unsigned long i = 0; i < length; i++)
   {
      mpz_clear(data1[i]);
      mpz_clear(data2[i]);
   }
   free(data1);
   free(data2);
   
   limb_release();
} */

/*
 Strategy: rotate random sequences of 3 words of various bitlengths, check
answer is correct and also that bitlength doesn't increase.
 */
/*void test_rotate_right_bits()
{
   printf("   testing rotate_right_bits...\n");
   
   mpz_t p, x, y;
   mpz_init(p);
   mpz_init(x);
   mpz_init(y);
   make_p(p, 2);

   for (unsigned trials = 0; trials < 20; trials++)
   {
      for (unsigned bits = 2*FLINT_BITS_PER_LIMB - 4;
           bits < 3*FLINT_BITS_PER_LIMB - 1; bits++)
      {
         for (unsigned s = 1; s < FLINT_BITS_PER_LIMB; s++)
         {
            mpz_urandomb(x, test_state, bits);
            if (gmp_urandomm_ui(test_state, 2))
               mpz_neg(x, x);
            
            mp_limb_t buf[5];
            convert_unsigned_mpz_to_raw(buf + 1, x, 3);
            if (mpz_sgn(x) < 0)
            negate_limbs(buf + 1, buf + 1, 3);
            
            mp_limb_t sentry = random_limb();
            buf[0] = buf[4] = sentry;
            
            rotate_right_bits(buf + 1, buf + 1, s, 2);
            assert(buf[0] == sentry);
            assert(buf[4] == sentry);

            convert_signed_raw_to_mpz(y, buf + 1, 2);
            mpz_mod(y, y, p);
            naive_rotate_right(x, s, p);
            mpz_mod(x, x, p);
            assert(mpz_cmp(x, y) == 0);
         }
      }
   }

   mpz_clear(y);
   mpz_clear(x);
   mpz_clear(p);
}*/

/*
 The strategy is: take p = B^2 + 1 (where B = one machine word). Take all
 possible 3-word integers whose individual words are between -3 and +3
 inclusive. Test reduce_mod_p_exact against naive reduction using GMP.
*/
/*void test_reduce_mod_p_exact()
{
   printf("   testing reduce_mod_p_exact...\n");

   mpz_t p, x;
   mpz_init(p);
   mpz_init(x);
   make_p(p, 2);

   mp_limb_signed_t y[3];
   mp_limb_t buf[5];
   mp_limb_t correct[3];
   mp_limb_t sentry = random_limb();
   buf[0] = buf[4] = sentry;

   for (y[0] = -3; y[0] <= 3; y[0]++)
      for (y[1] = -3; y[1] <= 3; y[1]++)
         for (y[2] = -3; y[2] <= 3; y[2]++)
         {
            buf[1] = y[0];
            buf[2] = y[1];
            buf[3] = y[2];

            convert_signed_raw_to_mpz(x, buf + 1, 2);

            reduce_mod_p_exact(buf + 1, 2);
            assert(buf[0] == sentry);
            assert(buf[4] == sentry);
            
            mpz_mod(x, x, p);
            convert_unsigned_mpz_to_raw(correct, x, 3);
            
            assert(correct[0] == buf[1]);
            assert(correct[1] == buf[2]);
            assert(correct[2] == buf[3]);
         }
   
   mpz_clear(p);
   mpz_clear(x);
}*/

void test_rotate_mod_p_limbs()
{
}

void test_rotate_mod_p_bits()
{
}

void test_fft_butterfly_limbs()
{
}

void test_fft_butterfly_bits()
{
}

void test_ifft_butterfly_limbs()
{
}

void test_ifft_butterfly_bits()
{
}

void test_fft_main()
{
}

void test_ifft_main()
{
}

// NOTE:
// The following block is from the old test.c. I'll be rewriting it somewhat
// soon. So it can live here happily for now. -- david
#if 0

void print_mpz_array(mpz_t* data, unsigned long count)
{
  for (unsigned long i = 0; i < count; i++)
    gmp_printf("%Zx\n", data[i]);
}

void test_rotate_mod_p_bits()
{
    /* We test rotate_mod_p() by taking a small fixed n, trying every allowable
    shift length, on a bunch of random inputs. If we try enough random
    inputs, we expect that every combination of carry/borrow in the
    algorithm will eventually occur (there shouldn't be too many
    combinations).
    
    We check that:
    * if the inputs don't use up more than 3/4 of the last limb, then
      the output don't use up more than 5/8 of the last limb
    * the answers agree modulo p with the naive mpz-based algorithm
    * the code hasn't changes the limb beyond the end of or before the
      beginning of the output buffer (coz that would be easy to do!!!!)  */
    
    const unsigned int n = 6;
    const unsigned int trials = 100;
    unsigned long s, i;

    // compute p = 2^(Bn) + 1
    mpz_t p;
    mpz_init(p);
    mpz_set_ui(p, 1);
    mpz_mul_2exp(p, p, n*FLINT_BITS_PER_LIMB);
    mpz_add_ui(p, p, 1);
    
    gmp_randstate_t test_state;
    gmp_randinit_default(test_state);
    
    mpz_t data[trials];
    for (i = 0; i < trials; i++)
    {
        mpz_init2(data[i], (n+1) * FLINT_BITS_PER_LIMB);
        // make random input take up about 3/4 of the overflow limb
        mpz_urandomb(data[i], test_state, n*FLINT_BITS_PER_LIMB + 3*FLINT_BITS_PER_LIMB/4);
        if (rand() % 2)
            mpz_neg(data[i], data[i]);
    }
    
    // a big power of two that we can use to compute twos complement
    mpz_t big_power;
    mpz_init(big_power);
    mpz_set_ui(big_power, 1);
    mpz_mul_2exp(big_power, big_power, (n+1) * FLINT_BITS_PER_LIMB);
    
    mp_limb_t* buf  = (mp_limb_t*) malloc(sizeof(mp_limb_t) * (n+3));
    mp_limb_t* buf2 = (mp_limb_t*) malloc(sizeof(mp_limb_t) * (n+3));
    
    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);

    for (s = 0; s < n * FLINT_BITS_PER_LIMB; s++)
    {
        for (i = 0; i < trials; i++)
        {
            // ==== Use rotate_mod_p() ==================================
            
            // convert integer to 2's complement form of the correct length
            mpz_add(a, big_power, data[i]);
            // now we have n+1 limbs to play with, plus an extra high limb
            
            // put them in a separate buffer and add sentry limbs at each end
            mpz_export(buf + 1, NULL, -1, sizeof(mp_limb_t), 0, 0, a);
            
            mpz_urandomb(a, test_state, FLINT_BITS_PER_LIMB * 2);
            mp_limb_t sentry1 = mpz_getlimbn(a, 0);
            mp_limb_t sentry2 = mpz_getlimbn(a, 1);
            buf[0] = sentry1;
            buf[n+2] = sentry2;
            
            // make an output buffer with sentry limbs at each end
            mpz_urandomb(a, test_state, FLINT_BITS_PER_LIMB * 2);
            mp_limb_t sentry3 = mpz_getlimbn(a, 0);
            mp_limb_t sentry4 = mpz_getlimbn(a, 1);
            buf2[0] = sentry3;
            buf2[n+2] = sentry4;
            
            // try the rotation
            rotate_mod_p_bits(buf2 + 1, buf + 1, s, n);
            
            // check sentries didn't change
            if (buf[0] != sentry1)
                FAIL("test_rotate_mod_p", "sentry1 overwritten");
            if (buf[n+2] != sentry2)
                FAIL("test_rotate_mod_p", "sentry2 overwritten");
            if (buf2[0] != sentry3)
                FAIL("test_rotate_mod_p", "sentry3 overwritten");
            if (buf2[n+2] != sentry4)
                FAIL("test_rotate_mod_p", "sentry4 overwritten");
            
            // check output doesn't use more than 5/8 of the last limb
            if (buf2[n+1] + (1L << (5*FLINT_BITS_PER_LIMB/8)) >=
                                       (1L << ((5*FLINT_BITS_PER_LIMB/8)+1)))
                FAIL("test_rotate_mod_p", "answer overflowed");

            // convert back to an integer
            mpz_import(a, n+1, -1, sizeof(mp_limb_t), 0, 0, buf2+1);
            if ((mp_limb_signed_t)(buf2[n+1]) < 0)
            {
                mpz_sub(a, big_power, a);
                mpz_neg(a, a);
            }
            mpz_mod(a, a, p);

            // ==== Use naive mpz method ================================
            
            // multiply by 2^s and reduce mod p
            mpz_mul_2exp(b, data[i], s);
            mpz_mod(b, b, p);
            
            if (mpz_cmp(a, b))
            {
                gmp_printf("  input: %Zd\n  input: %Zx\n output: %Zx\ncorrect: %Zx\n", data[i], data[i], a, b);
                char message[100];
                sprintf(message, "case %ld\n", s);
                FAIL("test_rotate_mod_p", message);
            }
        }
    }

    for (i = 0; i < trials; i++)
    
    mpz_clear(data[i]);
    mpz_clear(big_power);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(p);
    free(buf);
    free(buf2);
}

#endif



/* Tests SSMul on random data.

This function does not attempt to test corner cases. It just plugs in
random data and tests that the results match the results of the naive
algorithm.

num_trials = number of trials to perform.
log_size = log of the number of input coefficients.
coeff_bits = number of bits to use in each coefficient. */

void test_SSMul(unsigned long num_trials, unsigned long length,
                unsigned long coeff_bits)
{
   // alloc some space
   mpz_t* data1 = (mpz_t*) malloc(length * sizeof(mpz_t));
   mpz_t* data2 = (mpz_t*) malloc(length * sizeof(mpz_t));
   mpz_t* data3 = (mpz_t*) malloc((2*length-1) * sizeof(mpz_t));
   mpz_t* data4 = (mpz_t*) malloc((2*length-1) * sizeof(mpz_t));
   for (unsigned long i = 0; i < length; i++)
   {
      mpz_init(data1[i]);
      mpz_init(data2[i]);
   }
   for (unsigned long i = 0; i < 2*length-1; i++)
   {
      mpz_init(data3[i]);
      mpz_init(data4[i]);
   }
   
   for (unsigned long i = 0; i < num_trials; i++)
   {
      // make up random polys
      for (unsigned j = 0; j < length; j++)
      {
         mpz_urandomb(data1[j], test_state, coeff_bits);
         if (gmp_urandomm_ui(test_state, 2))
            mpz_neg(data1[j], data1[j]);
         
         mpz_urandomb(data2[j], test_state, coeff_bits);
         if (gmp_urandomb_ui(test_state, 2))
            mpz_neg(data2[j], data2[j]);
      }

      /*for (unsigned j = 0; j < length; j++)
      {
         gmp_printf("%Zd ",data1[j]);
      }
      printf("\n\n");
      for (unsigned j = 0; j < length; j++)
      {
         gmp_printf("%Zd ",data2[j]);
      }
      printf("\n\n");*/
      
      // compute product using naive algorithm
      /*for (unsigned long j = 0; j < 2*length-1; j++)
         mpz_set_ui(data3[j], 0);
      for (unsigned long j = 0; j < length; j++)
         for (unsigned long k = 0; k < length; k++)
           mpz_addmul(data3[j+k], data1[j], data2[k]);*/
      Zvec a;
      Zvec b;
      Zvec d;
      a->coords = data1;
      b->coords = data2;
      a->length = length;
      b->length = length;
      d->coords = data4;
      d->length = 2*length-1;
      // compute product using SSMul
      SSMul(d, a, b, coeff_bits, 1);
      
      /*for (unsigned j = 0; j < 2*length-1; j++)
      {
         gmp_printf("%Zx ",data3[j]);
      }
      printf("\n\n");*/

      // compare results
      for (unsigned j = 0; j < 2*length-1; j++)
      {
         //gmp_printf("%Zx ",data4[j]);
         assert(!mpz_cmp(data3[j], data4[j]));
      }
      //printf("\n\n\n\n");
      
   }
   
   // clean up
   for (unsigned long i = 0; i < length; i++)
   {
      mpz_clear(data1[i]);
      mpz_clear(data2[i]);
   }
   for (unsigned long i = 0; i < 2*length-1; i++)
   {
      mpz_clear(data3[i]);
      mpz_clear(data4[i]);
   }
  
   free(data1);
   free(data2);
   free(data3);
   free(data4);
}

void test_KSMul(unsigned long num_trials, unsigned long lengtha, unsigned lengthb,
                unsigned long coeff_bits)
{
   unsigned long length2 = 1;
   while (length2 < lengtha) length2*=2;
   while (length2 < lengthb) length2*=2;
   // alloc some space
   mpz_t* data1 = (mpz_t*) malloc(length2 * sizeof(mpz_t));
   mpz_t* data2 = (mpz_t*) malloc(length2 * sizeof(mpz_t));
   mpz_t* data3 = (mpz_t*) malloc((2*length2) * sizeof(mpz_t));
   mpz_t* data4 = (mpz_t*) malloc((2*length2) * sizeof(mpz_t));
   for (unsigned long i = 0; i < length2; i++)
   {
      mpz_init(data1[i]);
      mpz_init(data2[i]);
   }
   for (unsigned long i = 0; i < 2*length2; i++)
   {
      mpz_init(data3[i]);
      mpz_init(data4[i]);
   }
      
   Zvec a, b, c, d;
   a->coords = data1;
   b->coords = data2;
   c->coords = data3;
   d->coords = data4;
   
   for (unsigned long i = 0; i < num_trials; i++)
   {
      // make up random polys
      for (unsigned j = 0; j < lengtha; j++)
      {
         mpz_urandomb(data1[j], test_state, coeff_bits);
         /*if (gmp_urandomm_ui(test_state, 2))
            mpz_neg(data1[j], data1[j]);*/
      }   
      for (unsigned j = 0; j < lengthb; j++)
      {
         mpz_urandomb(data2[j], test_state, coeff_bits);
         /*if (gmp_urandomb_ui(test_state, 2))
            mpz_neg(data2[j], data2[j]);*/
      }
      for (unsigned j = lengtha; j < length2; j++) 
      {
          mpz_set_ui(data1[j],0);
      }
      for (unsigned j = lengthb; j < length2; j++) 
      {
          mpz_set_ui(data2[j],0);
      }
      
      for (unsigned j = 0; j < 2*length2; j++) 
      {
          mpz_set_ui(data3[j],0);
          mpz_set_ui(data4[j],0);
      }
      /*for (unsigned j = 0; j < lengtha; j++)
      {
         gmp_printf("%Zd ",data1[j]);
      }
      printf("\n\n");
      for (unsigned j = 0; j < lengthb; j++)
      {
         gmp_printf("%Zd ",data2[j]);
      }
      printf("\n\n");*/
      
      // compute product using naive algorithm
      /*for (unsigned long j = 0; j < 2*length-1; j++)
         mpz_set_ui(data3[j], 0);
      for (unsigned long j = 0; j < length; j++)
         for (unsigned long k = 0; k < length; k++)
           mpz_addmul(data3[j+k], data1[j], data2[k]);*/
      //a.length = length;
      //b.length = length;
      //c.length = 2*length-1;
      //if (a.length < 9) Zvec_karamul(c,a,b,2*coeff_bits+10);
      //else
      //{
       /*  a.length = length2;
         b.length = length2;
         c.length = 2*length2-1;
         Zvec_SSMul(c, a, b, coeff_bits, 0);*/
      //}
      
      a->length = lengtha;
      b->length = lengthb;
      c->length = lengtha+lengthb-1;
      d->length = lengtha+lengthb-1;

      /*unsigned long log_length = 0;
      while ((1<<log_length) < lengtha) log_length++;*/

      Zvec_karamul(c,a,b,coeff_bits);
      // compute product using SSMul
      Zvec_mul(d, a, b);
      //SSMul(data4, data1, data2, length, coeff_bits);
      /*for (unsigned j = 0; j < lengtha+lengthb-1; j++)
      {
         gmp_printf("%Zx ",c.coords[j]);
      }
      printf("\n\n");*/

      // compare results
      for (unsigned j = 0; j < lengtha + lengthb - 1; j++)
      {
         //gmp_printf("%Zx ",data4[j]);
         assert(!mpz_cmp(c->coords[j], d->coords[j]));
      }
      //printf("\n\n\n\n");
      
   }
   
   // clean up
   for (unsigned long i = 0; i < length2; i++)
   {
      mpz_clear(data1[i]);
      mpz_clear(data2[i]);
   }
   for (unsigned long i = 0; i < 2*length2; i++)
   {
      mpz_clear(data3[i]);
      mpz_clear(data4[i]);
   }
  
   free(data1);
   free(data2);
   free(data3);
   free(data4);
}

int main()
{
   initZvec();
    
   printf("Testing ssmul.c...\n");
   gmp_randinit_default(test_state);

   test_negate_limbs();
   test_signed_add_1();
   
  // printf("   testing convert_in_bytes and convert_out_bytes...\n");
   
   //for (unsigned long i = 1; i < 1000; i+=5)
     // for (unsigned long j = 1; j < 25; j+=3)
      //{
         //printf("%ld, %ld\n",i,j);
         //test_ssmul_convert_in_bytes(100, 100, i, i/4+j);  //trials, length, coeff_bits, pack_bytes
      //}
   
   //test_rotate_right_bits();
   //test_reduce_mod_p_exact();
   test_rotate_mod_p_limbs();
   test_rotate_mod_p_bits();
   
   test_fft_butterfly_limbs();
   test_fft_butterfly_bits();
   test_ifft_butterfly_limbs();
   test_ifft_butterfly_bits();
   
   test_fft_main();
   test_ifft_main();
   
   printf("   testing SSMul...\n");

   for (unsigned long i = 1; i <= 1000000; i+=((i/2)+2))
     for (unsigned long j = 1; j <= 10000; j+=((j/2)+2)) 
      for (unsigned long k = 1; k <= j; k+=((k/2)+2))  
      {
        //unsigned long k = 0;
        //while ((1UL<<k) < j) k++;
        //if ((2*i + k + 2 >= 64) && (i > j))
        //{
           FILE* outfile = fopen("testdata15","a");
           fprintf(outfile,"%ld, %ld, %ld\n",i, j, k);
           test_KSMul(10, j, k, i);
           fclose(outfile);
        //}
      }

   printf("All tests passed.\n");
}

// end of file ***************************************************************
