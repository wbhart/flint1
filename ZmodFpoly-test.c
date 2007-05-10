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
      bits = gmp_urandomm_ui(Zpoly_test_randstate,62)+ 2;
      
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
      bits = gmp_urandomm_ui(Zpoly_test_randstate,62)+ 2;
      
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


int test_ZmodFpoly_FFT()
{
   return 0;
}


int test_ZmodFpoly_IFFT()
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

   RUN_TEST(ZmodFpoly_convert);
   RUN_TEST(ZmodFpoly_convert_bits);
   RUN_TEST(ZmodFpoly_convert_bits_unsigned);
   RUN_TEST(ZmodFpoly_FFT);
   RUN_TEST(ZmodFpoly_IFFT);
   RUN_TEST(ZmodFpoly_convolution);
   RUN_TEST(ZmodFpoly_negacyclic_convolution);

   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}


int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   ZmodFpoly_test_all();

   return 0;
}



// end of file ****************************************************************
