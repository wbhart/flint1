/****************************************************************************

ZmodFpoly-test.c: test module for ZmodFpoly module

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "flint-manager.h"
#include "ZmodFpoly.h"
#include "Zpoly_mpn.h"
#include "Zpoly.h"

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

gmp_randstate_t Zpoly_test_randstate;


/****************************************************************************

   Test code for Conversion Routines
   
****************************************************************************/


unsigned long randint(unsigned long randsup) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(Zpoly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   if (pol->coeffs) Zpoly_clear(pol);
   Zpoly_init3(pol, length, maxbits);
   for (unsigned long i = 0; i < length; i++)
   {
       bits = randint(maxbits);
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, Zpoly_test_randstate, bits);
          if (randint(2)) mpz_neg(temp,temp);
       }
       Zpoly_set_coeff(pol, i, temp);
       
   }
   
   mpz_clear(temp);
} 

void randpoly_unsigned(Zpoly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   if (pol->coeffs) Zpoly_clear(pol);
   Zpoly_init3(pol, length, maxbits);
   for (unsigned long i = 0; i < length; i++)
   {
       bits = randint(maxbits);
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, Zpoly_test_randstate, bits);
       }
       Zpoly_set_coeff(pol, i, temp);
       
   }
   
   mpz_clear(temp);
} 

#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");

int test_ZmodFpoly_convert()
{
   Zpoly_t test_poly, test_poly2;
   Zpoly_mpn_t test_mpn_poly, test_mpn_poly2;
   ZmodFpoly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, depth;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      Zpoly_mpn_init(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      Zpoly_mpn_init(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000);
          depth = 0;
          while ((1<<depth) < length) depth++;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          Zpoly_mpn_realloc(test_mpn_poly2, length);
          Zpoly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits-1); 
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          ZmodFpoly_init(test_modF_poly, depth, (bits-1)/FLINT_BITS_PER_LIMB+1, 0);
          ZmodFpoly_convert_in_mpn(test_modF_poly, test_mpn_poly);
          ZmodFpoly_convert_out_mpn(test_mpn_poly2, test_modF_poly);
          _Zpoly_mpn_convert_out(test_poly2, test_mpn_poly2);
          
          ZmodFpoly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = Zpoly_equal(test_poly, test_poly2);
      }   
      Zpoly_mpn_clear(test_mpn_poly);
      Zpoly_mpn_clear(test_mpn_poly2);
   }
   
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   
   return result;
}

int test_ZmodFpoly_convert_bits()
{
   Zpoly_t test_poly, test_poly2;
   Zpoly_mpn_t test_mpn_poly, test_mpn_poly2;
   ZmodFpoly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, depth, bundle;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,FLINT_BITS_PER_LIMB-2)+ 2;
      
      Zpoly_mpn_init(test_mpn_poly, 1, 1);
      Zpoly_mpn_init(test_mpn_poly2, 1, 10);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;
          bundle = length/5;
          if (bundle == 0) bundle = length;
          depth = 0;
          while ((1<<depth) < (length-1)/bundle + 1) depth++;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          Zpoly_mpn_realloc(test_mpn_poly2, length);
          Zpoly_realloc(test_poly2, length);
          
          randpoly(test_poly, length, bits-1);
          for (unsigned long i = bundle-1; i < length; i+=bundle)
          {
             if (mpz_sgn(test_poly->coeffs[i])<0) // Final coeff in each bundle
                                                  // must be positive
                mpz_neg(test_poly->coeffs[i], test_poly->coeffs[i]);
             if (mpz_sgn(test_poly->coeffs[i]) == 0) 
                mpz_set_ui(test_poly->coeffs[i], 1);
          }
          if (mpz_sgn(test_poly->coeffs[length-1])<0) 
             mpz_neg(test_poly->coeffs[length-1], test_poly->coeffs[length-1]);
          if (mpz_sgn(test_poly->coeffs[length-1]) == 0) 
             mpz_set_ui(test_poly->coeffs[length-1], 1);


#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          ZmodFpoly_init(test_modF_poly, depth, (bits*bundle-1)/FLINT_BITS_PER_LIMB+1, 0);
          
          ZmodFpoly_bit_pack_mpn(test_modF_poly, test_mpn_poly, bundle, bits);
          test_mpn_poly2->length = length;
          
          for (unsigned long i = 0; i < length; i++) // Must clear coeffs in advance
             test_mpn_poly2->coeffs[i*(test_mpn_poly2->limbs+1)] = 0; 
             
          ZmodFpoly_bit_unpack_mpn(test_mpn_poly2, test_modF_poly, bundle, bits);  
          _Zpoly_mpn_convert_out(test_poly2, test_mpn_poly2);
          
          ZmodFpoly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = Zpoly_equal(test_poly, test_poly2);
      }   
      Zpoly_mpn_clear(test_mpn_poly);
      Zpoly_mpn_clear(test_mpn_poly2);
   }
   
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   
   return result;
}

int test_ZmodFpoly_convert_bits_unsigned()
{
   Zpoly_t test_poly, test_poly2;
   Zpoly_mpn_t test_mpn_poly, test_mpn_poly2;
   ZmodFpoly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, depth, bundle;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,FLINT_BITS_PER_LIMB-2)+ 2;
      
      Zpoly_mpn_init(test_mpn_poly, 1, 1);
      Zpoly_mpn_init(test_mpn_poly2, 1, 10);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;
          bundle = length/5;
          if (bundle == 0) bundle = length;
          depth = 0;
          while ((1<<depth) < (length-1)/bundle + 1) depth++;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          Zpoly_mpn_realloc(test_mpn_poly, length);
          Zpoly_mpn_realloc(test_mpn_poly2, length);
          Zpoly_realloc(test_poly2, length);
          
          randpoly_unsigned(test_poly, length, bits);

#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          _Zpoly_mpn_convert_in(test_mpn_poly, test_poly);
          ZmodFpoly_init(test_modF_poly, depth, (bits*bundle-1)/FLINT_BITS_PER_LIMB+1, 0);
          
          ZmodFpoly_bit_pack_mpn(test_modF_poly, test_mpn_poly, bundle, bits);
          test_mpn_poly2->length = length;
          
          for (unsigned long i = 0; i < length; i++) // Must clear coeffs in advance
             test_mpn_poly2->coeffs[i*(test_mpn_poly2->limbs+1)] = 0; 
             
          ZmodFpoly_bit_unpack_unsigned_mpn(test_mpn_poly2, test_modF_poly, bundle, bits);  
          _Zpoly_mpn_convert_out(test_poly2, test_mpn_poly2);
          
          ZmodFpoly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = Zpoly_equal(test_poly, test_poly2);
      }   
      Zpoly_mpn_clear(test_mpn_poly);
      Zpoly_mpn_clear(test_mpn_poly2);
   }
   
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   
   return result;
}


/****************************************************************************

   Test code for Fourier Transform Routines

****************************************************************************/


/*
Prints the ZmodF_t, each limb in a separate block, most significant limb
(i.e. the overflow limb) first.
*/
void ZmodF_print(ZmodF_t x, unsigned long n)
{
   for (long i = n; i >= 0; i--)
#if FLINT_BITS_PER_LIMB == 64
      printf("%016lx ", x[i]);
#else
      printf("%08lx ", x[i]);
#endif
}


/*
Prints each coefficient of the polynomial on a separate line.
*/
void ZmodFpoly_print(ZmodFpoly_t x)
{
   for (unsigned long k = 0; k < (1UL << x->depth); k++)
   {
      ZmodF_print(x->coeffs[k], x->n);
      printf("\n");
   }
}


unsigned long random_ulong(unsigned long max)
{
   return gmp_urandomm_ui(Zpoly_test_randstate, max);
}


/*
Generates a random ZmodFpoly_t with at most overflow_bits used in the
overflow limb for each coefficient.

The ZmodFpoly_t should already be initialised. This function ignores the
"length" attribute.
*/
void ZmodFpoly_random(ZmodFpoly_t x, unsigned long overflow_bits)
{
   unsigned long n = x->n;
   
   mpz_t temp;
   mpz_init(temp);

   for (unsigned long k = 0; k < (1UL << x->depth); k++)
   {
      ZmodF_t y = x->coeffs[k];
   
      ZmodF_zero(y, n);
      mpz_rrandomb(temp, Zpoly_test_randstate, (n+1)*FLINT_BITS_PER_LIMB);
      mpz_export(y, NULL, -1, sizeof(mp_limb_t), 0, 0, temp);

      // GMP has a "bug" where the top bit of the output of mpz_rrandomb
      // is always set. So we flip everything with probability 1/2.
      if (random_ulong(2))
         for (unsigned long i = 0; i <= n; i++)
            y[i] = ~y[i];

      // Copy the sign bit downwards so that only overflow_bits bits are used.
      if ((mp_limb_signed_t) y[n] >= 0)
         y[n] &= (1UL << overflow_bits) - 1;
      else
         y[n] |= ~((1UL << overflow_bits) - 1);
   }

   mpz_clear(temp);
}


mpz_t global_p;
unsigned long global_n = 0;


// Sets:
// global_n := n,
// global_p = 2^(FLINT_BITS_PER_LIMB*n) + 1
void set_global_n(unsigned long n)
{
   if (n != global_n)
   {
      global_n = n;
      mpz_set_ui(global_p, 1);
      mpz_mul_2exp(global_p, global_p, n*FLINT_BITS_PER_LIMB);
      mpz_add_ui(global_p, global_p, 1);
   }
}


/*
Converts given ZmodF_t into mpz_t format, reduced into [0, p) range.
Assumes global_n and global_p are set correctly.
*/
void ZmodF_convert_out(mpz_t output, ZmodF_t input)
{
   int negative = ((mp_limb_signed_t) input[global_n] < 0);
   
   if (negative)
      for (int i = 0; i <= global_n; i++)
         input[i] = ~input[i];
         
   mpz_import(output, global_n+1, -1, sizeof(mp_limb_t), 0, 0, input);
   
   if (negative)
   {
      mpz_add_ui(output, output, 1);
      mpz_neg(output, output);
      for (int i = 0; i <= global_n; i++)
         input[i] = ~input[i];
   }

   mpz_mod(output, output, global_p);
}


/*
Converts input polynomial to Zpoly format. Each output coefficient is
normalised into [0, p). All 2^depth coefficients are converted.

Assumes that output has already initialised.
*/
void ZmodFpoly_convert_out(Zpoly_t output, ZmodFpoly_t input)
{
   unsigned long size = 1UL << input->depth;
   unsigned long n = input->n;
   
   Zpoly_ensure_space(output, size);
   set_global_n(n);
   
   for (unsigned long k = 0; k < size; k++)
      ZmodF_convert_out(output->coeffs[k], input->coeffs[k]);
      
   output->length = size;
}


/*
y := x * 2^(s/2)  mod p    (using a very naive algorithm)
y may alias x
Assumes global_n and global_p are set correctly.
*/
void naive_mul_sqrt2exp(mpz_t y, mpz_t x, unsigned long s)
{
   static mpz_t temp;
   static int init = 0;
   
   if (!init)
   {
      mpz_init(temp);
      init = 1;
   }

   if (s & 1)
   {
      mpz_mul_2exp(y, x, s/2 + global_n*FLINT_BITS_PER_LIMB/4);
      mpz_mul_2exp(temp, y, global_n*FLINT_BITS_PER_LIMB/2);
      mpz_sub(y, temp, y);
      mpz_mod(y, y, global_p);
   }
   else
   {
      mpz_mul_2exp(y, x, s/2);
      mpz_mod(y, y, global_p);
   }
}


// root and twist are powers of sqrt2
void naive_FFT(Zpoly_t x, unsigned long depth, unsigned long root,
               unsigned long twist, unsigned long n)
{
   static mpz_t temp;
   static int init = 0;
   
   if (!init)
   {
      mpz_init(temp);
      init = 1;
   }

   unsigned long size = 1UL << depth;
   
   for (unsigned long d = 0; d < depth; d++)
   {
      unsigned long half = 1UL << (depth - d - 1);
      for (unsigned long start = 0; start < size; start += 2*half)
      {
         for (unsigned long i = 0; i < half; i++)
         {
            mpz_t* a = &x->coeffs[start + i];
            mpz_t* b = &x->coeffs[start + half + i];
            mpz_add(temp, *a, *b);
            mpz_sub(*b, *a, *b);
            naive_mul_sqrt2exp(*b, *b, twist + i*root);
            mpz_mod(*a, temp, global_p);
         }
      }
      root <<= 1;
      twist <<= 1;
   }
}


int test__ZmodFpoly_FFT()
{
   Zpoly_t poly1, poly2;
   Zpoly_init(poly1);
   Zpoly_init(poly2);
   int success = 1;

   for (unsigned long depth = 0; depth <= 10 && success; depth++)
   {
      unsigned long size = 1UL << depth;
   
      // need 4*n*FLINT_BITS_PER_LIMB divisible by 2^depth
      unsigned long n_skip = size / (4*FLINT_BITS_PER_LIMB);
      if (n_skip == 0)
         n_skip = 1;
         
      for (unsigned long n = n_skip; n < 8*n_skip && success; n += n_skip)
      {
         printf("depth = %d, n = %d\n", depth, n);
      
         for (unsigned long trial = 0; trial < 20; trial++)
         {
            ZmodFpoly_t f;
            ZmodFpoly_init(f, depth, n, 1);
            
            set_global_n(n);
            
            ZmodFpoly_random(f, 4);
            ZmodFpoly_convert_out(poly1, f);
            naive_FFT(poly1, depth, 4*n*FLINT_BITS_PER_LIMB / size, 0, n);

            _ZmodFpoly_FFT(f->coeffs, depth, 1, size, size, 0, n, f->scratch);
            ZmodFpoly_convert_out(poly2, f);
            
            if (!Zpoly_equal(poly1, poly2))
               success = 0;
            
            ZmodFpoly_clear(f);
         }
      }
   }

   Zpoly_clear(poly2);
   Zpoly_clear(poly1);
   return success;
}


int test__ZmodFpoly_IFFT()
{
   return 0;
}


int test_ZmodFpoly_convolution()
{
   return 0;
}


int test_ZmodFpoly_negacyclic_convolution()
{
   return 0;
}


/****************************************************************************

   Main test functions

****************************************************************************/


void ZmodFpoly_test_all()
{
   int success, all_success = 1;

#if 0    // just here temporarily so I don't have to run these tests every time
   RUN_TEST(ZmodFpoly_convert);
   RUN_TEST(ZmodFpoly_convert_bits);
   RUN_TEST(ZmodFpoly_convert_bits_unsigned);
#endif
   RUN_TEST(_ZmodFpoly_FFT);
   RUN_TEST(_ZmodFpoly_IFFT);
   RUN_TEST(ZmodFpoly_convolution);
   RUN_TEST(ZmodFpoly_negacyclic_convolution);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   mpz_init(global_p);
   
   ZmodFpoly_test_all();

   return 0;
}



// end of file ****************************************************************
