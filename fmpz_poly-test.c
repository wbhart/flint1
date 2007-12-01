/****************************************************************************

fmpz_poly-test.c: Test code for fmpz_poly.c and fmpz_poly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "mpz_poly.h"
#include "fmpz_poly.h"
#include "memory-manager.h"
#include "ZmodF_poly.h"
#include "test-support.h"

#define VARY_BITS 0
#define SIGNS 1
#define SPARSE 1

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

unsigned long randint(unsigned long randsup) 
{
    if (randsup == 0) return 0;
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(mpz_poly_t pol, long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_zero(pol);
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
#if SPARSE
          if (randint(10) == 1) mpz_rrandomb(temp, randstate, bits);
          else mpz_set_ui(temp, 0);
#else
          mpz_rrandomb(temp, randstate, bits);
#endif
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   mpz_clear(temp);
} 

void randpoly_unsigned(mpz_poly_t pol, long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_zero(pol);
   
   for (long i = 0; i < length; i++)
   {
#if VARY_BITS
       bits = randint(maxbits+1);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, randstate, bits);
       }
       mpz_poly_set_coeff(pol, i, temp);
   }
   
   mpz_clear(temp);
} 


#define RUN_TEST(targetfunc) \
   printf("Testing " #targetfunc "()... ");            \
   fflush(stdout);                                     \
   success = test_##targetfunc();                      \
   all_success = all_success && success;               \
   printf(success ? "ok\n" : "FAIL!\n");
   
int test__fmpz_poly_convert()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test__fmpz_poly_truncate()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
          trunc = random_ulong(length+1);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_truncate(test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly);
          mpz_poly_truncate(test_poly, test_poly, trunc);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test__fmpz_poly_max_bits()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   long next_bits, mpz_bits, fmpz_bits, sign;
   
   mpz_poly_init(test_poly); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_bits = 0;
          sign = 1L;
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
             next_bits = mpz_sizeinbase(test_poly->coeffs[i], 2);
             if (next_bits > mpz_bits) mpz_bits = next_bits;
             if (mpz_sgn(test_poly->coeffs[i]) < 0L) sign = -1L;
          }
          mpz_bits = sign*mpz_bits;
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_bits = _fmpz_poly_max_bits(test_fmpz_poly);
          
          result = (mpz_bits == fmpz_bits);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result;
}

int test__fmpz_poly_max_bits1()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   long next_bits, mpz_bits, fmpz_bits, sign;
   
   mpz_poly_init(test_poly); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(64)+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_bits = 0;
          sign = 1L;
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
             next_bits = mpz_sizeinbase(test_poly->coeffs[i], 2);
             if (next_bits > mpz_bits) mpz_bits = next_bits;
             if (mpz_sgn(test_poly->coeffs[i]) < 0L) sign = -1L;
          }
          mpz_bits = mpz_bits*sign;
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_bits = _fmpz_poly_max_bits1(test_fmpz_poly);
          
          result = (mpz_bits == fmpz_bits);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result;
}

int test__fmpz_poly_max_limbs()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length, next_limbs, mpz_limbs, fmpz_limbs;
   
   mpz_poly_init(test_poly); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_limbs = 0;
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
             next_limbs = mpz_size(test_poly->coeffs[i]);
             if (next_limbs > mpz_limbs) mpz_limbs = next_limbs;
          }
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_limbs = _fmpz_poly_max_limbs(test_fmpz_poly);
          
          result = (mpz_limbs == fmpz_limbs);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result;
}

int test__fmpz_poly_attach()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_attach(test_fmpz_poly2, test_fmpz_poly);
          _fmpz_poly_check_normalisation(test_fmpz_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test__fmpz_poly_attach_shift()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, shift;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
          shift = random_ulong(length+1);
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_attach_shift(test_fmpz_poly2, test_fmpz_poly, shift);
          _fmpz_poly_check_normalisation(test_fmpz_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          _fmpz_poly_right_shift(test_fmpz_poly, test_fmpz_poly, shift);
          fmpz_poly_to_mpz_poly(test_poly, test_fmpz_poly);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test__fmpz_poly_attach_truncate()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);
          trunc = random_ulong(length+1);
#if DEBUG
          printf("%ld, %ld, %ld\n",length, bits, trunc);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits);
          if (trunc > test_poly->length) trunc = test_poly->length; 
           
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_attach_truncate(test_fmpz_poly2, test_fmpz_poly, trunc);
          _fmpz_poly_check_normalisation(test_fmpz_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly, trunc);
          fmpz_poly_to_mpz_poly(test_poly, test_fmpz_poly);
          result = mpz_poly_equal(test_poly, test_poly2);
#if DEBUG2
          if (!result) 
          {
             mpz_poly_print(test_poly); printf("\n");
             mpz_poly_print(test_poly2); printf("\n"); 
          }
#endif
          
      }   
          
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test__fmpz_poly_getset_ui()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long coeff, coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
              
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(test_fmpz_poly->length);
              if (test_fmpz_poly->length)
              {
                 _fmpz_poly_set_coeff_ui(test_fmpz_poly, coeff_num, coeff);
                 fmpz_poly_check_normalisation(test_fmpz_poly);
                 result = (_fmpz_poly_get_coeff_ui(test_fmpz_poly, coeff_num) == coeff);
              }
#if DEBUG2
              if (!result) printf("Length = %ld, index = %ld, bits = %ld, coeff = %ld\n", test_fmpz_poly->length, coeff_num, coeff_bits, coeff);
#endif
          }
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_getset_si()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          length = test_fmpz_poly->length;
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(test_fmpz_poly->length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              if (randint(2)) sign = -1L; else sign = 1L;
              coeff = sign*coeff;
              if (test_fmpz_poly->length)
              {
                 _fmpz_poly_set_coeff_si(test_fmpz_poly, coeff_num, coeff);
                 fmpz_poly_check_normalisation(test_fmpz_poly);
                 result = ((_fmpz_poly_get_coeff_si(test_fmpz_poly, coeff_num) == coeff) && (_fmpz_poly_get_coeff_ui(test_fmpz_poly, coeff_num) == sign*coeff));
              }
          }
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_ui()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long coeff, coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
              
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(test_fmpz_poly->length);
              if (test_fmpz_poly->length)
              {
                 fmpz_poly_set_coeff_ui(test_fmpz_poly, coeff_num, coeff);
                 fmpz_poly_check_normalisation(test_fmpz_poly);
                 result = (fmpz_poly_get_coeff_ui(test_fmpz_poly, coeff_num) == coeff);
              }
#if DEBUG2
              if (!result) printf("Length = %ld, index = %ld, bits = %ld, coeff = %ld\n", test_fmpz_poly->length, coeff_num, coeff_bits, coeff);
#endif
          }
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_si()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          length = test_fmpz_poly->length;
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(test_fmpz_poly->length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              if (randint(2)) sign = -1L; else sign = 1L;
              coeff = sign*coeff;
              if (test_fmpz_poly->length)
              {
                 fmpz_poly_set_coeff_si(test_fmpz_poly, coeff_num, coeff);
                 fmpz_poly_check_normalisation(test_fmpz_poly);
                 result = ((fmpz_poly_get_coeff_si(test_fmpz_poly, coeff_num) == coeff) && (_fmpz_poly_get_coeff_ui(test_fmpz_poly, coeff_num) == sign*coeff));
              }
          }
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_get_coeff_ptr()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              if (randint(2)) sign = -1L; else sign = 1L;
              coeff = sign*coeff;
              _fmpz_poly_set_coeff_si(test_fmpz_poly, coeff_num, coeff);
              fmpz_poly_check_normalisation(test_fmpz_poly);
              if (coeff == 0) sign = 0;
              result = (_fmpz_poly_get_coeff_ptr(test_fmpz_poly, coeff_num)[0] == sign);
          }
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_get_coeff_ptr()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          do 
          {
             randpoly(test_poly, length, bits); 
             mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          } while (test_fmpz_poly->length == 0);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(test_fmpz_poly->length);
#if DEBUG
              printf("Length = %ld, index = %ld, bits = %ld, coeff = %ld\n", test_fmpz_poly->length, coeff_num, coeff_bits, coeff);
#endif
              if (randint(2)) sign = -1L; else sign = 1L;
              coeff = sign*coeff;
              _fmpz_poly_set_coeff_si(test_fmpz_poly, coeff_num, coeff);
              fmpz_poly_check_normalisation(test_fmpz_poly);
              if (coeff == 0) sign = 0;
              fmpz_t coeff_ptr = fmpz_poly_get_coeff_ptr(test_fmpz_poly, coeff_num);
              if (coeff_ptr != NULL) result = (coeff_ptr[0] == sign);
          }
          if (test_fmpz_poly->length)
          {
             for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
             {
                 coeff_num = randint(test_fmpz_poly->length)+test_fmpz_poly->length;
                 result = (fmpz_poly_get_coeff_ptr(test_fmpz_poly, coeff_num) == NULL);
             }
          }   
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_normalise()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long nz_coeff;
   long sign;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          nz_coeff = randint(length+1)-1;
          if (randint(2)) sign = -1L; else sign = 1;
          if (nz_coeff != -1L) _fmpz_poly_set_coeff_si(test_fmpz_poly, nz_coeff, sign*1000);
          for (unsigned long i = nz_coeff+1; i < length; i++)
            _fmpz_poly_set_coeff_ui(test_fmpz_poly, i, 0);
              
          _fmpz_poly_normalise(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
#if DEBUG
          printf("length = %ld, nonzero coefficient = %ld\n",_fmpz_poly_length(test_fmpz_poly), nz_coeff);
#endif              
          result = (_fmpz_poly_length(test_fmpz_poly) == nz_coeff+1);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_getset_coeff()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   unsigned long result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mp_limb_t * coeff1 = (mp_limb_t *) calloc(test_fmpz_poly->limbs, sizeof(mp_limb_t));
          mp_limb_t * coeff2 = (mp_limb_t *) calloc(test_fmpz_poly->limbs, sizeof(mp_limb_t));
          
          sign = _fmpz_poly_get_coeff(coeff1, test_fmpz_poly, randint(test_fmpz_poly->length));
          rand_coeff = randint(test_fmpz_poly->length);
          if (test_fmpz_poly->length)
          {
             _fmpz_poly_set_coeff(test_fmpz_poly, rand_coeff, coeff1, sign, test_fmpz_poly->limbs);
             fmpz_poly_check_normalisation(test_fmpz_poly);
             sign2 = _fmpz_poly_get_coeff(coeff2, test_fmpz_poly, rand_coeff);
             
             for (unsigned long i = 0; (i < test_fmpz_poly->limbs) && (result == 1); i++)
             {
                result = (coeff1[i] == coeff2[i]);
             }
                                 
             if (sign != sign2) result = 0;
          }
          
          free(coeff1);
          free(coeff2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_coeff()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   unsigned long result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          do
          {
             randpoly(test_poly, length, bits); 
             mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          } while (test_fmpz_poly->length != length);

          mp_limb_t * coeff1 = (mp_limb_t *) calloc(test_fmpz_poly->limbs, sizeof(mp_limb_t));
          mp_limb_t * coeff2 = (mp_limb_t *) calloc(test_fmpz_poly->limbs, sizeof(mp_limb_t));
          
          sign = fmpz_poly_get_coeff(coeff1, test_fmpz_poly, randint(test_fmpz_poly->length));
          rand_coeff = randint(test_fmpz_poly->length);
          if (test_fmpz_poly->length)
          {
             fmpz_poly_set_coeff(test_fmpz_poly, rand_coeff, coeff1, sign, test_fmpz_poly->limbs);
             fmpz_poly_check_normalisation(test_fmpz_poly);
             sign2 = fmpz_poly_get_coeff(coeff2, test_fmpz_poly, rand_coeff);
             
             for (unsigned long i = 0; (i < test_fmpz_poly->limbs) && (result == 1); i++)
             {
                result = (coeff1[i] == coeff2[i]);
             }
                                 
             if (sign != sign2) result = 0;
          }
          
          free(coeff1);
          free(coeff2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_getset_coeff_fmpz()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   unsigned long result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_t coeff1 = fmpz_init(test_fmpz_poly->limbs);
          fmpz_t coeff2 = fmpz_init(test_fmpz_poly->limbs);
          
          _fmpz_poly_get_coeff_fmpz(coeff1, test_fmpz_poly, randint(test_fmpz_poly->length));
          rand_coeff = randint(test_fmpz_poly->length);
          if (test_fmpz_poly->length)
          {
             _fmpz_poly_set_coeff_fmpz(test_fmpz_poly, rand_coeff, coeff1);
             fmpz_poly_check_normalisation(test_fmpz_poly);
             _fmpz_poly_get_coeff_fmpz(coeff2, test_fmpz_poly, rand_coeff);
             
             result = fmpz_equal(coeff1, coeff2);
          }
          
          fmpz_clear(coeff1);
          fmpz_clear(coeff2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_getset_coeff_mpz()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   unsigned long result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
           do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_t coeff1, coeff2;
          mpz_init(coeff1);
          mpz_init(coeff2);
          
          _fmpz_poly_get_coeff_mpz(coeff1, test_fmpz_poly, randint(test_fmpz_poly->length));
          rand_coeff = randint(test_fmpz_poly->length);
          if (test_fmpz_poly->length)
          {
             _fmpz_poly_set_coeff_mpz(test_fmpz_poly, rand_coeff, coeff1);
             fmpz_poly_check_normalisation(test_fmpz_poly);
             _fmpz_poly_get_coeff_mpz(coeff2, test_fmpz_poly, rand_coeff);
             
             result = (mpz_cmp(coeff1, coeff2) == 0);
          }
          
          mpz_clear(coeff1);
          mpz_clear(coeff2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_coeff_fmpz()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   unsigned long result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_t coeff1 = fmpz_init(test_fmpz_poly->limbs);
          fmpz_t coeff2 = fmpz_init(test_fmpz_poly->limbs);
          
          fmpz_poly_get_coeff_fmpz(coeff1, test_fmpz_poly, randint(test_fmpz_poly->length));
          rand_coeff = randint(test_fmpz_poly->length);
          if (test_fmpz_poly->length)
          {
             fmpz_poly_set_coeff_fmpz(test_fmpz_poly, rand_coeff, coeff1);
             fmpz_poly_check_normalisation(test_fmpz_poly);
             fmpz_poly_get_coeff_fmpz(coeff2, test_fmpz_poly, rand_coeff);
             
             result = fmpz_equal(coeff1, coeff2);
          }
          
          fmpz_clear(coeff1);
          fmpz_clear(coeff2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_coeff_mpz()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly;
   unsigned long result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
           do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_t coeff1, coeff2;
          mpz_init(coeff1);
          mpz_init(coeff2);
          
          fmpz_poly_get_coeff_mpz(coeff1, test_fmpz_poly, randint(test_fmpz_poly->length));
          rand_coeff = randint(test_fmpz_poly->length);
          if (test_fmpz_poly->length)
          {
             fmpz_poly_set_coeff_mpz(test_fmpz_poly, rand_coeff, coeff1);
             fmpz_poly_check_normalisation(test_fmpz_poly);
             fmpz_poly_get_coeff_mpz(coeff2, test_fmpz_poly, rand_coeff);
             
             result = (mpz_cmp(coeff1, coeff2) == 0);
          }
          
          mpz_clear(coeff1);
          mpz_clear(coeff2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_setequal()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long altered_coeff;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1+randint(30));
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly); 
      }
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          altered_coeff = randint(length);
          test_fmpz_poly2->coeffs[altered_coeff*(test_fmpz_poly2->limbs+1)+1]++;
          if (test_fmpz_poly2->coeffs[altered_coeff*(test_fmpz_poly2->limbs+1)] == 0)
             test_fmpz_poly2->coeffs[altered_coeff*(test_fmpz_poly2->limbs+1)] = 1;
          result = !_fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly); 
      }

      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          do
          {
             randpoly(test_poly, length, bits); 
          } while (test_poly->length != length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);          
          altered_coeff = randint(length);
          test_fmpz_poly2->coeffs[altered_coeff*(test_fmpz_poly2->limbs+1)]*=-1L;
          if (test_fmpz_poly2->coeffs[altered_coeff*(test_fmpz_poly2->limbs+1)] == 0)
             test_fmpz_poly2->coeffs[altered_coeff*(test_fmpz_poly2->limbs+1)] = 1;
          
          result = !_fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly); 
      }
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

#if TESTFILE
int test_fmpz_poly_freadprint()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1+randint(30));
      FILE * testfile; 
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          testfile = fopen("testfile", "w");
          fmpz_poly_fprint(test_fmpz_poly, testfile);
          fclose(testfile);
          testfile = fopen("testfile", "r");
          fmpz_poly_fread(test_fmpz_poly2, testfile);
          fclose(testfile);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
           
      }
      remove("testfile");
            
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}
#endif

int test_fmpz_poly_tofromstring()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1+randint(30));
      FILE * testfile; 
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          char * strbuf = fmpz_poly_to_string(test_fmpz_poly);
          int OK = fmpz_poly_from_string(test_fmpz_poly2, strbuf);
          free(strbuf);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly) && OK;
           
      }
            
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_zero_coeffs()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length, zeroes;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1+randint(30));
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
          
          zeroes = randint(test_fmpz_poly->length + 1);
          
          _fmpz_poly_zero_coeffs(test_fmpz_poly2, zeroes);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          unsigned long i;
          fmpz_t coeff1, coeff2;
          
          if (zeroes == test_fmpz_poly->length)
          {
             result = (test_fmpz_poly2->length == 0);
          } else
          {
             for (i = 0; i < zeroes; i++)
             {
                coeff1 = _fmpz_poly_get_coeff_ptr(test_fmpz_poly2, i);
                result &= fmpz_is_zero(coeff1);
             }
             for (i = zeroes; i < test_fmpz_poly->length; i++)
             {
                coeff1 = _fmpz_poly_get_coeff_ptr(test_fmpz_poly, i);
                coeff2 = _fmpz_poly_get_coeff_ptr(test_fmpz_poly2, i);
                result &= fmpz_equal(coeff1, coeff2);
             }
          }
      }
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_zero_coeffs()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length, zeroes;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1+randint(30));
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
          
          zeroes = randint(2*test_fmpz_poly->length);
          
          fmpz_poly_zero_coeffs(test_fmpz_poly2, zeroes);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          unsigned long i;
          fmpz_t coeff1, coeff2;
          
          if (zeroes >= test_fmpz_poly->length)
          {
             result = (test_fmpz_poly2->length == 0);
          } else
          {
             for (i = 0; i < zeroes; i++)
             {
                coeff1 = _fmpz_poly_get_coeff_ptr(test_fmpz_poly2, i);
                result &= fmpz_is_zero(coeff1);
             }
             for (i = zeroes; i < test_fmpz_poly->length; i++)
             {
                coeff1 = _fmpz_poly_get_coeff_ptr(test_fmpz_poly, i);
                coeff2 = _fmpz_poly_get_coeff_ptr(test_fmpz_poly2, i);
                result &= fmpz_equal(coeff1, coeff2);
             }
          }
      }
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_swap()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length, length2;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          fmpz_poly_fit_length(test_fmpz_poly3, length);
          randpoly(test_poly, length2, bits);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly);
          
          _fmpz_poly_set(test_fmpz_poly3, test_fmpz_poly);
          fmpz_poly_swap(test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly3);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test__fmpz_poly_shift()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long shift;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          shift = randint(100);
          fmpz_poly_fit_length(test_fmpz_poly, length+shift);
          fmpz_poly_fit_length(test_fmpz_poly2, length+shift);
          
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
          _fmpz_poly_left_shift(test_fmpz_poly, test_fmpz_poly, shift); 
          fmpz_poly_check_normalisation(test_fmpz_poly);
          _fmpz_poly_right_shift(test_fmpz_poly, test_fmpz_poly, shift);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      }

      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          shift = randint(length);
          fmpz_poly_fit_length(test_fmpz_poly, length+shift);
          fmpz_poly_fit_length(test_fmpz_poly2, length+shift);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          test_fmpz_poly3->limbs = test_fmpz_poly->limbs;
          test_fmpz_poly3->length = test_fmpz_poly->length-shift;
          test_fmpz_poly3->coeffs = test_fmpz_poly->coeffs+shift*(test_fmpz_poly->limbs+1);
          _fmpz_poly_right_shift(test_fmpz_poly2, test_fmpz_poly, shift);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
                
          result = _fmpz_poly_equal(test_fmpz_poly3, test_fmpz_poly2);
      }

      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test__fmpz_poly_neg()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length, check_coeff;
   unsigned long extra_bits1, extra_bits2;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      extra_bits1 = randint(200);
      extra_bits2 = randint(200);
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits+extra_bits1-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits+extra_bits1+extra_bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);      
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          fmpz_poly_fit_length(test_fmpz_poly3, length);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          
          check_coeff = randint(length);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_neg(test_fmpz_poly2, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          _fmpz_poly_neg(test_fmpz_poly3, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          if (length == 0) result = (test_fmpz_poly2->length == 0);
          else result = _fmpz_poly_equal(test_fmpz_poly, test_fmpz_poly3) 
             && (test_fmpz_poly->coeffs[check_coeff*(test_fmpz_poly->limbs+1)] == -test_fmpz_poly2->coeffs[check_coeff*(test_fmpz_poly2->limbs+1)]);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test__fmpz_poly_add()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, max_length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = bits+random_ulong(200);
      bits3 = bits2+random_ulong(200)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits3-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10); 
          length2 = random_ulong(10);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          max_length = (length > length2) ? length : length2;
          fmpz_poly_fit_length(test_fmpz_poly3, max_length);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          mpz_poly_add(test_poly3, test_poly, test_poly2);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          _fmpz_poly_add(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          _fmpz_poly_normalise(test_fmpz_poly3);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);

          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG2
          if (!result)
          {
             mpz_poly_print(test_poly); printf("\n");
             mpz_poly_print(test_poly2); printf("\n");
             mpz_poly_print(test_poly3); printf("\n");
             mpz_poly_print(test_poly4); printf("\n");
          }
#endif    
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }

   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = bits+random_ulong(200) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length + random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          mpz_poly_add(test_poly3, test_poly, test_poly2);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          _fmpz_poly_add(test_fmpz_poly2, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2);
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG2
          if (!result)
          {
             mpz_poly_print(test_poly); printf("\n");
             mpz_poly_print(test_poly2); printf("\n");
             mpz_poly_print(test_poly3); printf("\n");
             mpz_poly_print(test_poly4); printf("\n");
          }
#endif    

          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, bits/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          randpoly(test_poly, length, bits); 
          mpz_poly_add(test_poly3, test_poly, test_poly);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          _fmpz_poly_add(test_fmpz_poly2, test_fmpz_poly, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2);
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG2
          if (!result)
          {
             mpz_poly_print(test_poly); printf("\n");
             mpz_poly_print(test_poly3); printf("\n");
             mpz_poly_print(test_poly4); printf("\n");
          }
#endif    

          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test__fmpz_poly_sub()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, max_length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = bits+random_ulong(200);
      bits3 = bits2+random_ulong(200)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits3-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          max_length = (length > length2) ? length : length2;
          fmpz_poly_fit_length(test_fmpz_poly3, max_length);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          mpz_poly_sub(test_poly3, test_poly, test_poly2);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          _fmpz_poly_sub(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = bits+random_ulong(200) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length + random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          mpz_poly_sub(test_poly3, test_poly2, test_poly);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          _fmpz_poly_sub(test_fmpz_poly2, test_fmpz_poly2, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          randpoly(test_poly, length, bits); 
          mpz_poly_sub(test_poly3, test_poly, test_poly);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          _fmpz_poly_sub(test_fmpz_poly2, test_fmpz_poly, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test_fmpz_poly_add()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, max_length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = bits+random_ulong(200);
      bits3 = bits2+random_ulong(200)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10); 
          length2 = random_ulong(10);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          mpz_poly_add(test_poly3, test_poly, test_poly2);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          fmpz_poly_add(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);

          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG2
          if (!result)
          {
             mpz_poly_print(test_poly); printf("\n");
             mpz_poly_print(test_poly2); printf("\n");
             mpz_poly_print(test_poly3); printf("\n");
             mpz_poly_print(test_poly4); printf("\n");
          }
#endif    
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }

   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = bits+random_ulong(200) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length + random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          mpz_poly_add(test_poly3, test_poly, test_poly2);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          fmpz_poly_add(test_fmpz_poly2, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2);
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG2
          if (!result)
          {
             mpz_poly_print(test_poly); printf("\n");
             mpz_poly_print(test_poly2); printf("\n");
             mpz_poly_print(test_poly3); printf("\n");
             mpz_poly_print(test_poly4); printf("\n");
          }
#endif    

          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test_fmpz_poly_sub()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, max_length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = bits+random_ulong(200);
      bits3 = bits2+random_ulong(200)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          mpz_poly_sub(test_poly3, test_poly, test_poly2);

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          fmpz_poly_sub(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = bits+random_ulong(200) + 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length + random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          mpz_poly_sub(test_poly3, test_poly2, test_poly);
          
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          fmpz_poly_sub(test_fmpz_poly2, test_fmpz_poly2, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test__fmpz_poly_scalar_mul_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739);
          _fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   for (unsigned long count1 = 1; (count1 < 8000) && (result == 1) ; count1++)
   {
      bits = random_ulong(150)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(40);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739);
          _fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_mul_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739);
          fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 8000) && (result == 1) ; count1++)
   {
      bits = random_ulong(150)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = random_ulong(40);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739);
          fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 8000) && (result == 1) ; count1++)
   {
      bits = random_ulong(150)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(40);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739);
          fmpz_poly_scalar_mul_ui(test_fmpz_poly, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test__fmpz_poly_scalar_mul_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length, sign;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = (long) randint(34682739);
          sign = randint(2);
          if (sign) mult = -mult;
          
          _fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_si(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_mul_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length, sign;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = (long) randint(34682739);
          sign = randint(2);
          if (sign) mult = -mult;
          
          fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_si(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }

   for (unsigned long count1 = 1; (count1 < 8000) && (result == 1) ; count1++)
   {
      bits = random_ulong(150)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = random_ulong(40);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = (long) randint(34682739);
          sign = randint(2);
          if (sign) mult = -mult;
          
          fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_si(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 8000) && (result == 1) ; count1++)
   {
      bits = random_ulong(150)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = (long) randint(34682739);
          sign = randint(2);
          if (sign) mult = -mult;
          
          fmpz_poly_scalar_mul_si(test_fmpz_poly, test_fmpz_poly, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_si(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }

   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test__fmpz_poly_scalar_div_exact_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          fmpz_poly_fit_length(test_fmpz_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_fmpz_poly3, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          result = _fmpz_poly_equal(test_fmpz_poly3, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_fmpz_poly2, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits-1)/FLINT_BITS+3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          fmpz_poly_fit_length(test_fmpz_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_fmpz_poly3, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          result = _fmpz_poly_equal(test_fmpz_poly3, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_fmpz_poly2, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_scalar_div_exact_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length;
   long mult;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          fmpz_poly_fit_length(test_fmpz_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_fmpz_poly3, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          result = _fmpz_poly_equal(test_fmpz_poly3, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("2:length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_fmpz_poly2, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_fmpz_poly3, 1, (bits-1)/FLINT_BITS+3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          fmpz_poly_fit_length(test_fmpz_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_fmpz_poly3, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          result = _fmpz_poly_equal(test_fmpz_poly3, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_fmpz_poly2, test_fmpz_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_fmpz_poly2, test_fmpz_poly2, mult);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_scalar_tdiv_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long div;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 20) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      //bits = 64*1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          //length = 10000;       
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          for (unsigned long i = 0; i < 100; i++)
          {
             _fmpz_poly_scalar_tdiv_ui(test_fmpz_poly2, test_fmpz_poly, div);
             fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly2->length; i++)
          {
              mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          for (unsigned long i = test_fmpz_poly2->length; i < test_poly->length; i++)
          {
              mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
          }         
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      //bits = 64*1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 50) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          //length = 10000;       
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          _fmpz_poly_scalar_tdiv_ui(test_fmpz_poly, test_fmpz_poly, div);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly->length; i++)
          {
              mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          for (unsigned long i = test_fmpz_poly->length; i < test_poly->length; i++)
          {
              mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
          }         
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test__fmpz_poly_scalar_div_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long div;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 20) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      //bits = 64*1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          //length = 10000;       
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          for (unsigned long i = 0; i < 100; i++)
          {
             _fmpz_poly_scalar_div_ui(test_fmpz_poly2, test_fmpz_poly, div);
             fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly2->length; i++)
          {
              mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          for (unsigned long i = test_fmpz_poly2->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
          }         
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      //bits = 64*1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          //length = 10000;       
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          _fmpz_poly_scalar_div_ui(test_fmpz_poly, test_fmpz_poly, div);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly->length; i++)
          {
              mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          for (unsigned long i = test_fmpz_poly->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
          }         
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test__fmpz_poly_scalar_tdiv_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   long div;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          if (randint(2)) div = -div;
          _fmpz_poly_scalar_tdiv_si(test_fmpz_poly2, test_fmpz_poly, div);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly2->length; i++)
          {
              if (div < 0)
              {
                 mpz_tdiv_q_ui(temp, test_poly->coeffs[i], -div);
                 mpz_neg(temp, temp);
              } else
                 mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }   
          for (unsigned long i = test_fmpz_poly2->length; i < test_poly->length; i++)
          {
              if (div < 0) mpz_tdiv_q_ui(temp, test_poly->coeffs[i], -div);
              else mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
          }              
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          if (randint(2)) div = -div;
          _fmpz_poly_scalar_tdiv_si(test_fmpz_poly, test_fmpz_poly, div);
          fmpz_poly_check_normalisation(test_fmpz_poly);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly->length; i++)
          {
              if (div < 0)
              {
                 mpz_tdiv_q_ui(temp, test_poly->coeffs[i], -div);
                 mpz_neg(temp, temp);
              } else
                 mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }   
          for (unsigned long i = test_fmpz_poly->length; i < test_poly->length; i++)
          {
              if (div < 0) mpz_tdiv_q_ui(temp, test_poly->coeffs[i], -div);
              else mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
          }              
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test__fmpz_poly_scalar_div_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length;
   long div;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          if (randint(2)) div = -div;
          _fmpz_poly_scalar_div_si(test_fmpz_poly2, test_fmpz_poly, div);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly2->length; i++)
          {
              if (div < 0L)
              {
                 mpz_cdiv_q_ui(temp, test_poly->coeffs[i], -div);
                 mpz_neg(temp, temp);
              } else
                 mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
#if DEBUG2
              if (!result)
              {
                 gmp_printf("%Zd, %ld, %Zd, %Zd\n", test_poly->coeffs[i], div, temp, test_poly2->coeffs[i]);
                 break;
              }
#endif
          }   
          for (unsigned long i = test_fmpz_poly2->length; i < test_poly->length; i++)
          {
              if (div < 0L) mpz_cdiv_q_ui(temp, test_poly->coeffs[i], -div);
              else mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
#if DEBUG2
              if (!result)
              {
                 gmp_printf("%Zd, %ld, %Zd\n", test_poly->coeffs[i], div, temp);
                 break;
              }
#endif
          }              
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000);        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          div = randint(34682739)+1;
          if (randint(2)) div = -div;
          _fmpz_poly_scalar_div_si(test_fmpz_poly, test_fmpz_poly, div);
          fmpz_poly_check_normalisation(test_fmpz_poly);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif    
          for (unsigned long i = 0; i < test_fmpz_poly->length; i++)
          {
              if (div < 0L)
              {
                 mpz_cdiv_q_ui(temp, test_poly->coeffs[i], -div);
                 mpz_neg(temp, temp);
              } else
                 mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
#if DEBUG2
              if (!result)
              {
                 gmp_printf("%Zd, %ld, %Zd, %Zd\n", test_poly->coeffs[i], div, temp, test_poly2->coeffs[i]);
                 break;
              }
#endif
          }   
          for (unsigned long i = test_fmpz_poly->length; i < test_poly->length; i++)
          {
              if (div < 0L) mpz_cdiv_q_ui(temp, test_poly->coeffs[i], -div);
              else mpz_fdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp_ui(temp, 0) == 0);
#if DEBUG2
              if (!result)
              {
                 gmp_printf("%Zd, %ld, %Zd\n", test_poly->coeffs[i], div, temp);
                 break;
              }
#endif
          }              
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test__fmpz_poly_mul_classical()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      
      length = random_ulong(35);   
      unsigned long log_length = 0;
      while ((1<<log_length) < length) log_length++;
      unsigned long output_bits = 2*bits + log_length;
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do randpoly(test_poly, length, bits); 
      while (mpz_poly_length(test_poly) < length);
    
      fmpz_poly_fit_length(test_fmpz_poly, length);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n");
#endif          
      mpz_poly_mul_naive(test_poly3, test_poly, test_poly);
          
      mpz_poly_init(test_poly4);
      if (length) fmpz_poly_init2(test_fmpz_poly3, 2*length-1, 2*test_fmpz_poly->limbs+1); //(output_bits-1)/FLINT_BITS+1);
      else fmpz_poly_init(test_fmpz_poly3);
          
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3); 
#if DEBUG
          mpz_poly_print(test_poly3);printf("\n");
          mpz_poly_print(test_poly4);printf("\n");
#endif          
          result = mpz_poly_equal(test_poly4, test_poly3);
      }
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_fit_length(test_fmpz_poly, length+length2-1);
          fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_fit_length(test_fmpz_poly2, length+length2-1);
          fmpz_poly_fit_limbs(test_fmpz_poly2, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly2, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_classical_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100); 
          if (length+length2 == 0) trunc = 0;
          else trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_classical_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          if (length == 0) trunc = 0;
          else trunc = random_ulong(2*length);       
#if DEBUG
          printf("length = %ld, trunc = %ld, bits = %ld\n", length, trunc, bits);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_fmpz_poly3, 2*length-1, (2*bits-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (2*bits-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_classical_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100); 
          if (length+length2 == 0) trunc = 0;
          else trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_fit_length(test_fmpz_poly, trunc);
          fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_classical_trunc(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_classical_trunc_left()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100); 
          length2 = random_ulong(100); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2 < 2)
          {
             fmpz_poly_init2(test_fmpz_poly3, 1, (bits+bits2-1)/FLINT_BITS+2);
             fmpz_poly_init2(test_fmpz_poly4, 1, (bits+bits2-1)/FLINT_BITS+2);
          } else
          {
             fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
             fmpz_poly_init2(test_fmpz_poly4, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          }
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, FLINT_MIN(trunc, test_fmpz_poly3->length));
          
          _fmpz_poly_mul_classical_trunc_left(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100); 
          if (length < 1) trunc = 0;
          else trunc = random_ulong(2*length-1)+1;       
#if DEBUG
          printf("length = %ld, trunc = %ld, bits = %ld\n", length, trunc, bits);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_poly_init(test_poly4);
          if (length < 1)
          {
             fmpz_poly_init2(test_fmpz_poly3, 1, (2*bits-1)/FLINT_BITS+2);
             fmpz_poly_init2(test_fmpz_poly4, 1, (2*bits-1)/FLINT_BITS+2);
          } else
          {
             fmpz_poly_init2(test_fmpz_poly3, 2*length-1, (2*bits-1)/FLINT_BITS+2);
             fmpz_poly_init2(test_fmpz_poly4, 2*length-1, (2*bits-1)/FLINT_BITS+2);
          }
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, FLINT_MIN(trunc, test_fmpz_poly3->length));
          
          _fmpz_poly_mul_classical_trunc_left(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100); 
          length2 = random_ulong(100); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2 < 2)
          {
             fmpz_poly_init2(test_fmpz_poly3, 1, (bits+bits2-1)/FLINT_BITS+2);
             fmpz_poly_fit_length(test_fmpz_poly, 1);
             fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          } else
          {
             fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
             fmpz_poly_fit_length(test_fmpz_poly, length+length2-1);
             fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          }
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, FLINT_MIN(trunc, test_fmpz_poly3->length));
          
          _fmpz_poly_mul_classical_trunc_left(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_karatsuba()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2, max_length, log_length, output_bits;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(35); 
      length = random_ulong(35);   
      max_length = FLINT_MIN(length, length2);
      log_length = 0;
      while ((1<<log_length) < max_length) log_length++;
      output_bits = bits + bits2 + log_length;
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do randpoly(test_poly, length, bits); 
      while (mpz_poly_length(test_poly) < length);

      do randpoly(test_poly2, length2, bits2); 
      while (mpz_poly_length(test_poly2) < length2);
          
      fmpz_poly_fit_length(test_fmpz_poly, length);
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n");
      mpz_poly_print(test_poly2);printf("\n");
#endif          
      mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1); //(output_bits-1)/FLINT_BITS+1);
      else fmpz_poly_init(test_fmpz_poly3);
          
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_karatsuba(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3); 
#if DEBUG
          mpz_poly_print(test_poly3);printf("\n");
          mpz_poly_print(test_poly4);printf("\n");
#endif          
          result = mpz_poly_equal(test_poly4, test_poly3);
      }
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      
      length = random_ulong(35);   
      log_length = 0;
      while ((1<<log_length) < length) log_length++;
      output_bits = 2*bits + log_length;
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do randpoly(test_poly, length, bits); 
      while (mpz_poly_length(test_poly) < length);
    
      fmpz_poly_fit_length(test_fmpz_poly, length);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n");
#endif          
      mpz_poly_mul_naive(test_poly3, test_poly, test_poly);
          
      mpz_poly_init(test_poly4);
      if (length) fmpz_poly_init2(test_fmpz_poly3, 2*length-1, 2*test_fmpz_poly->limbs+1); //(output_bits-1)/FLINT_BITS+1);
      else fmpz_poly_init(test_fmpz_poly3);
          
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_karatsuba(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3); 
#if DEBUG
          mpz_poly_print(test_poly3);printf("\n");
          mpz_poly_print(test_poly4);printf("\n");
#endif          
          result = mpz_poly_equal(test_poly4, test_poly3);
      }
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(35); 
      length = random_ulong(35);   
      max_length = FLINT_MIN(length, length2);
      log_length = 0;
      while ((1<<log_length) < max_length) log_length++;
      output_bits = bits + bits2 + log_length;
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do randpoly(test_poly, length, bits); 
      while (mpz_poly_length(test_poly) < length);

      do randpoly(test_poly2, length2, bits2); 
      while (mpz_poly_length(test_poly2) < length2);
          
      fmpz_poly_fit_length(test_fmpz_poly, length);
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n");
      mpz_poly_print(test_poly2);printf("\n");
#endif          
      mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      if (length + length2) 
      {
         fmpz_poly_fit_length(test_fmpz_poly, length+length2-1);
         fmpz_poly_fit_limbs(test_fmpz_poly, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1); //(output_bits-1)/FLINT_BITS+1);
      }
          
      _fmpz_poly_mul_karatsuba(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2);
      fmpz_poly_check_normalisation(test_fmpz_poly);
          
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 
#if DEBUG
      mpz_poly_print(test_poly3);printf("\n");
      mpz_poly_print(test_poly4);printf("\n");
#endif          
      result = mpz_poly_equal(test_poly4, test_poly3);
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(35); 
      length = random_ulong(35);   
      max_length = FLINT_MIN(length, length2);
      log_length = 0;
      while ((1<<log_length) < max_length) log_length++;
      output_bits = bits + bits2 + log_length;
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do randpoly(test_poly, length, bits); 
      while (mpz_poly_length(test_poly) < length);

      do randpoly(test_poly2, length2, bits2); 
      while (mpz_poly_length(test_poly2) < length2);
          
      fmpz_poly_fit_length(test_fmpz_poly, length);
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n");
      mpz_poly_print(test_poly2);printf("\n");
#endif          
      mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      if (length + length2) 
      {
         fmpz_poly_fit_length(test_fmpz_poly2, length+length2-1);
         fmpz_poly_fit_limbs(test_fmpz_poly2, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1); //(output_bits-1)/FLINT_BITS+1);
      }
          
      _fmpz_poly_mul_karatsuba(test_fmpz_poly2, test_fmpz_poly, test_fmpz_poly2);
      fmpz_poly_check_normalisation(test_fmpz_poly2);
          
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly2); 
#if DEBUG
      mpz_poly_print(test_poly3);printf("\n");
      mpz_poly_print(test_poly4);printf("\n");
#endif          
      result = mpz_poly_equal(test_poly4, test_poly3);
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test__fmpz_poly_mul_karatsuba_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length2 = random_ulong(35); 
          length = random_ulong(35);   
          if (length+length2 == 0) trunc = 0;
          else trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          
          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_karatsuba_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          result = _fmpz_poly_equal(test_fmpz_poly4, test_fmpz_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(35);   
          if (length == 0) trunc = 0;
          else trunc = random_ulong(length*2);       
#if DEBUG
          printf("length = %ld, trunc = %ld, bits = %ld\n", length, trunc, bits);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          if (length) fmpz_poly_init2(test_fmpz_poly3, 2*length-1, (2*bits-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (2*bits-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_karatsuba_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          result = _fmpz_poly_equal(test_fmpz_poly4, test_fmpz_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length2 = random_ulong(35); 
          length = random_ulong(35);   
          if (length+length2 == 0) trunc = 0;
          else trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          
          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_fit_length(test_fmpz_poly, trunc);
          fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_karatsuba_trunc(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          result = _fmpz_poly_equal(test_fmpz_poly, test_fmpz_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          fmpz_poly_clear(test_fmpz_poly3); 
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result; 
}

int test__fmpz_poly_mul_karatsuba_trunc_left()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100); 
          length2 = random_ulong(100); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) 
          {
             fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
             fmpz_poly_init2(test_fmpz_poly4, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          } else
          {
              fmpz_poly_init(test_fmpz_poly3);
              fmpz_poly_init(test_fmpz_poly4);
          }
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, FLINT_MIN(trunc, test_fmpz_poly3->length));
          
          _fmpz_poly_mul_karatsuba_trunc_left(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100); 
          if (length < 1) trunc = 0;
          else trunc = random_ulong(2*length-1)+1;       
#if DEBUG
          printf("length = %ld, trunc = %ld, bits = %ld\n", length, trunc, bits);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_poly_init(test_poly4);
          if (length) 
          {
             fmpz_poly_init2(test_fmpz_poly3, 2*length-1, (2*bits-1)/FLINT_BITS+2);
             fmpz_poly_init2(test_fmpz_poly4, 2*length-1, (2*bits-1)/FLINT_BITS+2);
          } else
          {
              fmpz_poly_init(test_fmpz_poly3);
              fmpz_poly_init(test_fmpz_poly4);
          }
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, FLINT_MIN(trunc, test_fmpz_poly3->length));
          
          _fmpz_poly_mul_karatsuba_trunc_left(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100); 
          length2 = random_ulong(100); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) 
          {
             fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
             fmpz_poly_fit_length(test_fmpz_poly, length+length2-1);
             fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          } else
          {
              fmpz_poly_init(test_fmpz_poly3);
          }
          
          _fmpz_poly_mul_classical(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, FLINT_MIN(trunc, test_fmpz_poly3->length));
          
          _fmpz_poly_mul_karatsuba_trunc_left(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_KS()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   ZmodF_poly_t test_modF_poly;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   mpz_poly_init(test_poly4); 
      
   for (unsigned long count1 = 0; count1 < 60; count1++)
   {
      bits = random_ulong(100) + 1;
      bits2 = random_ulong(100) + 1;
      //bits = 4;
      //bits2 = 4;
      length2 = random_ulong(100); 
      length = random_ulong(100);   
      //length = 32000;
      //length2 = 32000;   
   
      _fmpz_poly_stack_init(test_fmpz_poly, length, (bits-1)/FLINT_BITS+1);
      _fmpz_poly_stack_init(test_fmpz_poly2, length2, (bits2-1)/FLINT_BITS+1);
      if (length + length2) _fmpz_poly_stack_init(test_fmpz_poly3, length + length2 - 1, test_fmpz_poly->limbs + test_fmpz_poly2->limbs + 1);
      else fmpz_poly_init(test_fmpz_poly3);

#if DEBUG
      printf("%ld, %ld, %ld, %ld\n", length, length2, bits, bits2);
#endif
      mpz_poly_realloc(test_poly, length);
      mpz_poly_realloc(test_poly2, length2);
      if (length + length2)
      {
         mpz_poly_realloc(test_poly3, length + length2 - 1);
         mpz_poly_realloc(test_poly4, length + length2 - 1);
      }    
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
         do randpoly(test_poly, length, bits);
         while (mpz_poly_length(test_poly) < length);

         do randpoly(test_poly2, length2, bits2);
         while (mpz_poly_length(test_poly2) < length2);
          
#if DEBUG
         if (bits2 == 64)
         {
            printf("Input poly 1:\n");
            for (unsigned j = 0; j < test_poly->length; j++)
               gmp_printf("%Zx, ",test_poly->coeffs[j]);
            printf("\n\n");
            printf("Input poly 2:\n");
            for (unsigned j = 0; j < test_poly2->length; j++)
               gmp_printf("%Zx, ",test_poly2->coeffs[j]);
            printf("\n\n");
         }
#endif
         mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
         mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
         for (unsigned long i = 0; i < 10; i++)
            _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
         fmpz_poly_check_normalisation(test_fmpz_poly3);
         fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
                  
         result = mpz_poly_equal(test_poly3, test_poly4);
         
#if DEBUG
         if (!result)
         {
            printf("Output poly correct\n");
            for (unsigned j = 0; j < test_poly3->length; j++)
               gmp_printf("%Zx, ",test_poly3->coeffs[j]);
            printf("\n\n");
            printf("Output poly incorrect\n");
            for (unsigned j = 0; j < test_poly4->length; j++)
               gmp_printf("%Zx, ",test_poly4->coeffs[j]);
            printf("\n\n");
         }
#endif         

      }   
   
      _fmpz_poly_stack_clear(test_fmpz_poly3);
      _fmpz_poly_stack_clear(test_fmpz_poly2);
      _fmpz_poly_stack_clear(test_fmpz_poly);
   
   }
   
   for (unsigned long count1 = 0; count1 < 60; count1++)
   {
      bits = random_ulong(100) + 1;

      length = random_ulong(100);   
   
      _fmpz_poly_stack_init(test_fmpz_poly, length, (bits-1)/FLINT_BITS+1);
      if (length) _fmpz_poly_stack_init(test_fmpz_poly3, 2*length - 1, 2*test_fmpz_poly->limbs + 1);
      else fmpz_poly_init(test_fmpz_poly3);

#if DEBUG
      printf("%ld, %ld\n", length, bits);
#endif
      mpz_poly_realloc(test_poly, length);
      if (length)
      {
         mpz_poly_realloc(test_poly3, 2*length - 1);
         mpz_poly_realloc(test_poly4, 2*length - 1);
      }    
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
         do randpoly(test_poly, length, bits);
         while (mpz_poly_length(test_poly) < length);
          
#if DEBUG
         if (bits2 == 64)
         {
            printf("Input poly 1:\n");
            for (unsigned j = 0; j < test_poly->length; j++)
               gmp_printf("%Zx, ",test_poly->coeffs[j]);
            printf("\n\n");
         }
#endif
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
         mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly);
          
         _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly3);
         fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
                  
         result = mpz_poly_equal(test_poly3, test_poly4);
         
#if DEBUG
         if (!result)
         {
            printf("Output poly correct\n");
            for (unsigned j = 0; j < test_poly3->length; j++)
               gmp_printf("%Zx, ",test_poly3->coeffs[j]);
            printf("\n\n");
            printf("Output poly incorrect\n");
            for (unsigned j = 0; j < test_poly4->length; j++)
               gmp_printf("%Zx, ",test_poly4->coeffs[j]);
            printf("\n\n");
         }
#endif         

      }   
   
      _fmpz_poly_stack_clear(test_fmpz_poly3);
      _fmpz_poly_stack_clear(test_fmpz_poly);
   
   }
   
   for (unsigned long count1 = 0; count1 < 60; count1++)
   {
      bits = random_ulong(100) + 1;
      bits2 = random_ulong(100) + 1;
      //bits = 4;
      //bits2 = 4;
      length2 = random_ulong(100); 
      length = random_ulong(100);   
      //length = 32000;
      //length2 = 32000;   
   
      fmpz_poly_init2(test_fmpz_poly, length, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, length2, (bits2-1)/FLINT_BITS+1);
      if (length + length2) 
      {
         fmpz_poly_fit_length(test_fmpz_poly, length + length2 - 1);
         fmpz_poly_fit_limbs(test_fmpz_poly, test_fmpz_poly->limbs + test_fmpz_poly2->limbs + 1);
      }
      
#if DEBUG
      printf("%ld, %ld, %ld, %ld\n", length, length2, bits, bits2);
#endif
      mpz_poly_realloc(test_poly, length);
      mpz_poly_realloc(test_poly2, length2);
      if (length + length2)
      {
         mpz_poly_realloc(test_poly3, length + length2 - 1);
         mpz_poly_realloc(test_poly4, length + length2 - 1);
      }    
      
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
         do randpoly(test_poly, length, bits);
         while (mpz_poly_length(test_poly) < length);

         do randpoly(test_poly2, length2, bits2);
         while (mpz_poly_length(test_poly2) < length2);
          
#if DEBUG
         if (bits2 == 64)
         {
            printf("Input poly 1:\n");
            for (unsigned j = 0; j < test_poly->length; j++)
               gmp_printf("%Zx, ",test_poly->coeffs[j]);
            printf("\n\n");
            printf("Input poly 2:\n");
            for (unsigned j = 0; j < test_poly2->length; j++)
               gmp_printf("%Zx, ",test_poly2->coeffs[j]);
            printf("\n\n");
         }
#endif
         mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
         mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
         _fmpz_poly_mul_KS(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2);
         fmpz_poly_check_normalisation(test_fmpz_poly);
         fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly);
                  
         result = mpz_poly_equal(test_poly3, test_poly4);
         
#if DEBUG
         if (!result)
         {
            printf("Output poly correct\n");
            for (unsigned j = 0; j < test_poly3->length; j++)
               gmp_printf("%Zx, ",test_poly3->coeffs[j]);
            printf("\n\n");
            printf("Output poly incorrect\n");
            for (unsigned j = 0; j < test_poly4->length; j++)
               gmp_printf("%Zx, ",test_poly4->coeffs[j]);
            printf("\n\n");
         }
#endif         

      }   
   
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly);
   
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   mpz_poly_clear(test_poly4);
   
   return result;
}

int test__fmpz_poly_mul_KS_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      bits2 = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(2000); 
          length2 = random_ulong(2000); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_KS_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(2000); 
          if (length < 1) trunc = 0;
          else trunc = random_ulong(2*length-1)+1;       
#if DEBUG
          printf("length = %ld, trunc = %ld, bits = %ld\n", length, trunc, bits);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_poly_init(test_poly4);
          if (length) fmpz_poly_init2(test_fmpz_poly3, 2*length-1, (2*bits-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (2*bits-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_KS_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      bits2 = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(2000); 
          length2 = random_ulong(2000); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_fit_length(test_fmpz_poly, trunc);
          fmpz_poly_fit_limbs(test_fmpz_poly, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_KS_trunc(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
#if 1
   for (unsigned long count1 = 1; (count1 < 10) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(500); 
          length2 = random_ulong(500); 
          if (length+length2 < 2) trunc = 0;
          else trunc = random_ulong(length+length2-1)+1;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_KS_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
#endif

   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test__fmpz_poly_mul_SS()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(35); 
      length = random_ulong(35);   
      //length = length2 = 256;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do randpoly(test_poly, length, bits);
      while (mpz_poly_length(test_poly) < length);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
          
      fmpz_poly_fit_length(test_fmpz_poly, length);
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
      mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
      else fmpz_poly_init(test_fmpz_poly3);
          
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_SS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
      }
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
           
      result = mpz_poly_equal(test_poly4, test_poly3);
      if (!result)
      {
#if DEBUG
          mpz_poly_print(test_poly);printf("\n\n");
          mpz_poly_print(test_poly2);printf("\n\n");
          mpz_poly_print(test_poly3);printf("\n\n");
          mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      }
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      
      length = random_ulong(35);   
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do randpoly(test_poly, length, bits);
      while (mpz_poly_length(test_poly) < length);
     
      fmpz_poly_fit_length(test_fmpz_poly, length);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
      mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly);
          
      mpz_poly_init(test_poly4);
      if (length) fmpz_poly_init2(test_fmpz_poly3, 2*length-1, 2*test_fmpz_poly->limbs+1);
      else fmpz_poly_init(test_fmpz_poly3);
          
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_SS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
      }
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
           
      result = mpz_poly_equal(test_poly4, test_poly3);
      if (!result)
      {
#if DEBUG
          mpz_poly_print(test_poly);printf("\n\n");
          mpz_poly_print(test_poly3);printf("\n\n");
          mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      }
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(35); 
      length = random_ulong(35);   
      //length = length2 = 256;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do randpoly(test_poly, length, bits);
      while (mpz_poly_length(test_poly) < length);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
          
      fmpz_poly_fit_length(test_fmpz_poly, length);
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
      mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      if (length + length2) 
      {
         fmpz_poly_fit_length(test_fmpz_poly, length+length2-1);
         fmpz_poly_fit_limbs(test_fmpz_poly, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
      }    
      
      _fmpz_poly_mul_SS(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2);
      fmpz_poly_check_normalisation(test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly);
           
      result = mpz_poly_equal(test_poly4, test_poly3);
      if (!result)
      {
#if DEBUG
          mpz_poly_print(test_poly);printf("\n\n");
          mpz_poly_print(test_poly2);printf("\n\n");
          mpz_poly_print(test_poly3);printf("\n\n");
          mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      }
      
      mpz_poly_clear(test_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_SS_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 30) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = random_ulong(1000); 
          if (length+length2 == 0) trunc = 0;
          else trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_SS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_SS_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 

          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 30) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          if (length == 0) trunc = 0;
          else trunc = random_ulong(2*length);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          mpz_poly_init(test_poly4);
          if (length) fmpz_poly_init2(test_fmpz_poly3, 2*length-1, (2*bits-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (2*bits-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_SS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_SS_trunc(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 

          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 30) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = random_ulong(1000); 
          if (length+length2 == 0) trunc = 0;
          else trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) 
          {
             fmpz_poly_fit_length(test_fmpz_poly, length+length2-1);
             fmpz_poly_fit_limbs(test_fmpz_poly, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
             fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          } else fmpz_poly_init(test_fmpz_poly3);
          
          _fmpz_poly_mul_SS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_SS_trunc(test_fmpz_poly, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly); 

          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3); 
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_trunc_n()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 25) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length; 
          trunc = length;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init2(test_fmpz_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          _fmpz_poly_mul_trunc_n(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_trunc_n()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 25) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length; 
          trunc = length;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init(test_fmpz_poly4);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_truncate(test_fmpz_poly3, trunc);
          _fmpz_poly_normalise(test_fmpz_poly3);
          
          fmpz_poly_mul_trunc_n(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul_trunc_left_n()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 25) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length; 
          trunc = length;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly4, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly4);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, trunc);
          
          _fmpz_poly_mul_trunc_left_n(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          _fmpz_poly_zero_coeffs(test_fmpz_poly4, trunc);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_trunc_left_n()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 25) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000); 
          length2 = length; 
          trunc = length;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length2);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          if (length + length2) fmpz_poly_init2(test_fmpz_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          else fmpz_poly_init(test_fmpz_poly3);
          fmpz_poly_init(test_fmpz_poly4);
          
          _fmpz_poly_mul_KS(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
          _fmpz_poly_zero_coeffs(test_fmpz_poly3, trunc);
          
          fmpz_poly_mul_trunc_left_n(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly2, trunc);
          fmpz_poly_check_normalisation(test_fmpz_poly4);
          _fmpz_poly_zero_coeffs(test_fmpz_poly4, trunc);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_fmpz_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_fmpz_poly3);
          fmpz_poly_clear(test_fmpz_poly4);
      }
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test__fmpz_poly_mul()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   ZmodF_poly_t test_modF_poly;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   mpz_poly_init(test_poly4); 
      
   for (unsigned long count1 = 0; count1 < 20; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = random_ulong(1000) + 1;
      //bits = 4;
      //bits2 = 4;
      length2 = random_ulong(1000); 
      length = random_ulong(1000);   
      //length = 32000;
      //length2 = 32000;   
   
      _fmpz_poly_stack_init(test_fmpz_poly, length, (bits-1)/FLINT_BITS+1);
      _fmpz_poly_stack_init(test_fmpz_poly2, length2, (bits2-1)/FLINT_BITS+1);
      if (length + length2) _fmpz_poly_stack_init(test_fmpz_poly3, length + length2 - 1, test_fmpz_poly->limbs + test_fmpz_poly2->limbs + 1);
      else fmpz_poly_init(test_fmpz_poly3);

#if DEBUG
      printf("%ld, %ld, %ld, %ld\n", length, length2, bits, bits2);
#endif
      mpz_poly_realloc(test_poly, length);
      mpz_poly_realloc(test_poly2, length2);
      if (length + length2)
      {
         mpz_poly_realloc(test_poly3, length + length2 - 1);
         mpz_poly_realloc(test_poly4, length + length2 - 1);
      }    
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      { 
         do randpoly(test_poly, length, bits);
         while (mpz_poly_length(test_poly) < length);

         do randpoly(test_poly2, length2, bits2);
         while (mpz_poly_length(test_poly2) < length2);
          
         mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
         mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
         for (unsigned long i = 0; i < 2; i++)
            _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
         fmpz_poly_check_normalisation(test_fmpz_poly3);
         fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
                  
         result = mpz_poly_equal(test_poly3, test_poly4);
         
#if DEBUG
         if (!result)
         {
            printf("Output poly correct\n");
            for (unsigned j = 0; j < test_poly3->length; j++)
               gmp_printf("%Zx, ",test_poly3->coeffs[j]);
            printf("\n\n");
            printf("Output poly incorrect\n");
            for (unsigned j = 0; j < test_poly4->length; j++)
               gmp_printf("%Zx, ",test_poly4->coeffs[j]);
            printf("\n\n");
         }
#endif         

      }   
   
      _fmpz_poly_stack_clear(test_fmpz_poly3);
      _fmpz_poly_stack_clear(test_fmpz_poly2);
      _fmpz_poly_stack_clear(test_fmpz_poly);
   
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   mpz_poly_clear(test_poly4);
   
   return result;
}

int test_fmpz_poly_mul()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   ZmodF_poly_t test_modF_poly;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   mpz_poly_init(test_poly4); 
      
   for (unsigned long count1 = 0; count1 < 20; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = random_ulong(1000) + 1;
      //bits = 4;
      //bits2 = 4;
      length2 = random_ulong(1000); 
      length = random_ulong(1000);   
      //length = 32000;
      //length2 = 32000;   
   
      _fmpz_poly_stack_init(test_fmpz_poly, length, (bits-1)/FLINT_BITS+1);
      _fmpz_poly_stack_init(test_fmpz_poly2, length2, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly3);

#if DEBUG
      printf("%ld, %ld, %ld, %ld\n", length, length2, bits, bits2);
#endif
      mpz_poly_realloc(test_poly, length);
      mpz_poly_realloc(test_poly2, length2);
      if (length + length2)
      {
         mpz_poly_realloc(test_poly3, length + length2 - 1);
         mpz_poly_realloc(test_poly4, length + length2 - 1);
      }    
      
      for (unsigned long count2 = 0; (count2 < 20) && (result == 1); count2++)
      { 
         do randpoly(test_poly, length, bits);
         while (mpz_poly_length(test_poly) < length);

         do randpoly(test_poly2, length2, bits2);
         while (mpz_poly_length(test_poly2) < length2);
          
         mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
         mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
         for (unsigned long i = 0; i < 2; i++)
            fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
         fmpz_poly_check_normalisation(test_fmpz_poly3);
         fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
                  
         result = mpz_poly_equal(test_poly3, test_poly4);
         
#if DEBUG
         if (!result)
         {
            printf("Output poly correct\n");
            for (unsigned j = 0; j < test_poly3->length; j++)
               gmp_printf("%Zx, ",test_poly3->coeffs[j]);
            printf("\n\n");
            printf("Output poly incorrect\n");
            for (unsigned j = 0; j < test_poly4->length; j++)
               gmp_printf("%Zx, ",test_poly4->coeffs[j]);
            printf("\n\n");
         }
#endif         

      }   
   
      fmpz_poly_clear(test_fmpz_poly3);
      _fmpz_poly_stack_clear(test_fmpz_poly2);
      _fmpz_poly_stack_clear(test_fmpz_poly);
   
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   mpz_poly_clear(test_poly4);
   
   return result;
}

int test__fmpz_poly_scalar_mul()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, bits2, limbs2, length;
   mpz_t temp, x_mpz;
   mpz_init(temp);
   mpz_init(x_mpz);
   mp_limb_t * x;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 7) && (result == 1) ; count1++)
   {
      bits = randint(100000) + 150000;
      bits2 = randint(10000) + 150000;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+(bits2-1)/FLINT_BITS+2);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 3) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 1; j++)
          {
            _fmpz_poly_scalar_mul(test_fmpz_poly2, test_fmpz_poly, x);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = randint(1000) + 1500;
      bits2 = randint(1000) + 1500;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+(bits2-1)/FLINT_BITS+2);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            _fmpz_poly_scalar_mul(test_fmpz_poly2, test_fmpz_poly, x);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   mpz_clear(x_mpz);
   
   return result; 
}

int test__fmpz_poly_scalar_div()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, bits2, limbs, limbs2, length;
   mpz_t temp, x_mpz;
   mpz_init(temp);
   mpz_init(x_mpz);
   mp_limb_t * x;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(300)+1;
      bits2 = randint(300)+1;
      limbs = (bits-1)/FLINT_BITS+1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      if (limbs >= limbs2) fmpz_poly_init2(test_fmpz_poly2, 1, limbs-limbs2+1);
      else fmpz_poly_init2(test_fmpz_poly2, 1, 1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          while (!x[limbs2]) mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            _fmpz_poly_scalar_div(test_fmpz_poly2, test_fmpz_poly, x);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (long i = 0; i < test_poly2->length; i++)
          {
              mpz_fdiv_q(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          } 
          for (long i = test_poly2->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q(temp, test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp_ui(temp, 0L) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld not zero\n", i, test_poly->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 gmp_printf("%Zd, %Zd, %Zd\n", temp, test_poly->coeffs[i], x_mpz);
                 break;
              }
          } 
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(300)+1;
      bits2 = randint(300)+1;
      limbs = (bits-1)/FLINT_BITS+1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          while (!x[limbs2]) mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          _fmpz_poly_scalar_div(test_fmpz_poly, test_fmpz_poly, x);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (long i = 0; i < test_poly2->length; i++)
          {
              mpz_fdiv_q(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          } 
          for (long i = test_poly2->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q(temp, test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp_ui(temp, 0L) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld not zero\n", i, test_poly->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 gmp_printf("%Zd, %Zd, %Zd\n", temp, test_poly->coeffs[i], x_mpz);
                 break;
              }
          } 
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }

   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   mpz_clear(x_mpz);
   
   return result; 
}

int test_fmpz_poly_scalar_div()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, bits2, limbs, limbs2, length;
   mpz_t temp, x_mpz;
   mpz_init(temp);
   mpz_init(x_mpz);
   mp_limb_t * x;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(300)+1;
      bits2 = randint(300)+1;
      limbs = (bits-1)/FLINT_BITS+1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly2);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          while (!x[limbs2]) mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            fmpz_poly_scalar_div(test_fmpz_poly2, test_fmpz_poly, x);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (long i = 0; i < test_poly2->length; i++)
          {
              mpz_fdiv_q(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          } 
          for (long i = test_poly2->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q(temp, test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp_ui(temp, 0L) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld not zero\n", i, test_poly->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 gmp_printf("%Zd, %Zd, %Zd\n", temp, test_poly->coeffs[i], x_mpz);
                 break;
              }
          } 
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(300)+1;
      bits2 = randint(300)+1;
      limbs = (bits-1)/FLINT_BITS+1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 30) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          while (!x[limbs2]) mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_poly_scalar_div(test_fmpz_poly, test_fmpz_poly, x);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (long i = 0; i < test_poly2->length; i++)
          {
              mpz_fdiv_q(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          } 
          for (long i = test_poly2->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q(temp, test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp_ui(temp, 0L) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld not zero\n", i, test_poly->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 gmp_printf("%Zd, %Zd, %Zd\n", temp, test_poly->coeffs[i], x_mpz);
                 break;
              }
          } 
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }

   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   mpz_clear(x_mpz);
   
   return result; 
}

int test_fmpz_poly_scalar_div_mpz()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, bits2, limbs, limbs2, length;
   mpz_t temp, x_mpz;
   mpz_init(temp);
   mpz_init(x_mpz);
   mp_limb_t * x;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(300)+1;
      bits2 = randint(300)+1;
      limbs = (bits-1)/FLINT_BITS+1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly2);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          while (!x[limbs2]) mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            fmpz_poly_scalar_div_mpz(test_fmpz_poly2, test_fmpz_poly, x_mpz);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (long i = 0; i < test_poly2->length; i++)
          {
              mpz_fdiv_q(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          } 
          for (long i = test_poly2->length; i < test_poly->length; i++)
          {
              mpz_fdiv_q(temp, test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp_ui(temp, 0L) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld not zero\n", i, test_poly->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 gmp_printf("%Zd, %Zd, %Zd\n", temp, test_poly->coeffs[i], x_mpz);
                 break;
              }
          } 
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   mpz_clear(x_mpz);
   
   return result; 
}

int test_fmpz_poly_scalar_mul()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, bits2, limbs2, length;
   mpz_t temp, x_mpz;
   mpz_init(temp);
   mpz_init(x_mpz);
   mp_limb_t * x;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 7) && (result == 1) ; count1++)
   {
      bits = randint(100000) + 150000;
      bits2 = randint(10000) + 150000;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 3) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 1; j++)
          {
             fmpz_poly_scalar_mul(test_fmpz_poly2, test_fmpz_poly, x);
             fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = randint(1000) + 1500;
      bits2 = randint(1000) + 1500;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            fmpz_poly_scalar_mul(test_fmpz_poly2, test_fmpz_poly, x);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }

   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(150) + 1;
      bits2 = randint(150) + 1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = randint(40);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 2; j++)
          {
            fmpz_poly_scalar_mul(test_fmpz_poly2, test_fmpz_poly, x);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 printf("bits 2 actually equals %ld\n",mpz_sizeinbase(x_mpz,2));
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(10) + 1;
      bits2 = randint(10) + 1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_poly_scalar_mul(test_fmpz_poly, test_fmpz_poly, x);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 printf("bits 2 actually equals %ld\n", mpz_sizeinbase(x_mpz, 2));
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = randint(6400) + 1;
      bits2 = randint(6400) + 1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_poly_scalar_mul(test_fmpz_poly, test_fmpz_poly, x);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 printf("bits 2 actually equals %ld\n",mpz_sizeinbase(x_mpz,2));
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 0; (count1 < 20) && (result == 1) ; count1++)
   {
      bits = randint(128000) + 1;
      bits2 = randint(128000) + 1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          fmpz_poly_scalar_mul(test_fmpz_poly, test_fmpz_poly, x);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 printf("bits 2 actually equals %ld\n",mpz_sizeinbase(x_mpz,2));
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   mpz_clear(x_mpz);
   
   return result; 
}

int test_fmpz_poly_scalar_mul_mpz()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, bits2, limbs2, length;
   mpz_t temp, x_mpz;
   mpz_init(temp);
   mpz_init(x_mpz);
   mp_limb_t * x;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 7) && (result == 1) ; count1++)
   {
      bits = randint(100000) + 150000;
      bits2 = randint(10000) + 150000;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 3) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 1; j++)
          {
             fmpz_poly_scalar_mul_mpz(test_fmpz_poly2, test_fmpz_poly, x_mpz);
             fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = randint(1000) + 1500;
      bits2 = randint(1000) + 1500;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = randint(100);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            fmpz_poly_scalar_mul_mpz(test_fmpz_poly2, test_fmpz_poly, x_mpz);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }

   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = randint(150) + 1;
      bits2 = randint(150) + 1;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          fmpz_poly_init(test_fmpz_poly2);
          length = randint(40);        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          
          F_mpn_clear(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          mpz_import(x_mpz, ABS(x[0]), -1, sizeof(mp_limb_t), 0, 0, x+1);
          if ((long) x[0] < 0)
             mpz_neg(x_mpz, x_mpz);
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          
          for (unsigned long j = 0; j < 2; j++)
          {
            fmpz_poly_scalar_mul_mpz(test_fmpz_poly2, test_fmpz_poly, x_mpz);
            fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_fmpz_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul(test_poly->coeffs[i], test_poly->coeffs[i], x_mpz);
              result &= (mpz_cmp(test_poly->coeffs[i], test_poly2->coeffs[i]) == 0);
              if (!result) 
              {
                 printf("Coefficient %ld of %ld incorrect\n", i, test_poly2->length); 
                 printf("bits = %ld, bits2 = %ld, length1 = %ld, length2 = %ld\n", bits, bits2, test_poly->length, test_poly2->length);
                 printf("bits 2 actually equals %ld\n",mpz_sizeinbase(x_mpz,2));
              }
          }  
#if DEBUG
          if (!result)
          {
             mpz_poly_print(test_poly);printf("\n\n");
             mpz_poly_print(test_poly2);printf("\n\n");
          }
#endif
         
          mpz_poly_clear(test_poly2);
          fmpz_poly_clear(test_fmpz_poly2);
      }
      free(x);
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   mpz_clear(x_mpz);
   
   return result; 
}

int test_fmpz_poly_div_classical()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 1000000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(100); 
      length = random_ulong(100)+1; 
      //length = length2 = 20;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_div_classical(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_div_classical(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
      _fmpz_poly_normalise(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_divrem_classical()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 1000000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(100); 
      length = random_ulong(100)+1; 
      //length = length2 = 20;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_divrem_classical(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_divrem_classical(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_div_divconquer_recursive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 0; i < 10; i++)
      {
         fmpz_poly_div_divconquer_recursive(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_div_divconquer_recursive(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_divrem_divconquer()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 0; i < 10; i++)
      {
         fmpz_poly_divrem_divconquer(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_divrem_divconquer(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_divrem()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 0; i < 10; i++)
      {
         fmpz_poly_divrem(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_divrem(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   for (unsigned long count1 = 0; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 2;
      bits2 = random_ulong(200)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      fmpz_poly_divrem(test_fmpz_poly3, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly3);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly5);
   }

   for (unsigned long count1 = 0; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 2;
      bits2 = random_ulong(200)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      fmpz_poly_divrem(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }

   for (unsigned long count1 = 0; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 2;
      bits2 = random_ulong(200)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      fmpz_poly_divrem(test_fmpz_poly4, test_fmpz_poly, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }

   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_div()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 3000) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 0; i < 10; i++)
      {
         fmpz_poly_div(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_div(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   for (unsigned long count1 = 0; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 2;
      bits2 = random_ulong(200)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      fmpz_poly_div(test_fmpz_poly3, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly3);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }

   for (unsigned long count1 = 0; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(200)+ 2;
      bits2 = random_ulong(200)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      fmpz_poly_div(test_fmpz_poly, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }

   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_div_divconquer()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      //bits = bits2 = 100;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(256)+1; 
      length = random_ulong(256)+1;
      //length = 1000;
      //length2 = 1000;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 1; i < 10; i++)
      {
         fmpz_poly_div_divconquer(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         //fmpz_poly_clear(test_fmpz_poly4);
         //fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_div_divconquer(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_div_mulders()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      bits2 = random_ulong(100)+ 1;
      //bits = bits2 = 10000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 1000;
      //length2 = 1000;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 1; i < 10; i++)
      {
         fmpz_poly_div_mulders(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly2);
         //fmpz_poly_clear(test_fmpz_poly4);
         //fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_div_mulders(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_newton_invert_basecase()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length, n;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 20000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 2;
      //bits = 100000;
      
      fmpz_poly_init(test_fmpz_poly);
      fmpz_poly_init(test_fmpz_poly2);
      fmpz_poly_init(test_fmpz_poly3);
      
      
      length = random_ulong(64)+1;
      //length = 12;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      fmpz_poly_set_coeff_ui(test_fmpz_poly, test_fmpz_poly->length - 1, 1L);
      
      n = randint(test_fmpz_poly->length) + 1;
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      fmpz_poly_newton_invert_basecase(test_fmpz_poly2, test_fmpz_poly, n);
      fmpz_poly_check_normalisation(test_fmpz_poly2);
      
      fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
           
      for (unsigned long i = 0; i < n - 1; i++)
      {
          result &= (test_fmpz_poly3->coeffs[(i+test_fmpz_poly3->length-n)*(test_fmpz_poly3->limbs+1)] == 0);
      }
      result &= (test_fmpz_poly3->coeffs[(test_fmpz_poly3->length-1)*(test_fmpz_poly3->limbs+1)] == 1);
      result &= (test_fmpz_poly3->coeffs[(test_fmpz_poly3->length-1)*(test_fmpz_poly3->limbs+1)+1] == 1);
      
#if DEBUG
      if (!result)
      {
         fmpz_poly_print(test_fmpz_poly); printf("\n");
         fmpz_poly_print(test_fmpz_poly2); printf("\n");
         fmpz_poly_print(test_fmpz_poly3); printf("\n");
      }
#endif
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test__fmpz_poly_reverse()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   int result = 1;
   unsigned long bits, length, length2;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);   
      
      length = random_ulong(100);
      length2 = length + randint(200);
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
#endif

      randpoly(test_poly, length, bits); 
      fmpz_poly_fit_length(test_fmpz_poly, length);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      _fmpz_poly_normalise(test_fmpz_poly);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      _fmpz_poly_reverse(test_fmpz_poly2, test_fmpz_poly, length2);
      fmpz_poly_check_normalisation(test_fmpz_poly2);
      _fmpz_poly_reverse(test_fmpz_poly2, test_fmpz_poly2, length2);
      fmpz_poly_check_normalisation(test_fmpz_poly2);
           
      result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);   
      
      length = random_ulong(100);
      length2 = length + randint(200);
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
#endif

      randpoly(test_poly, length, bits); 
      fmpz_poly_fit_length(test_fmpz_poly, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      _fmpz_poly_normalise(test_fmpz_poly);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length);
      
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      _fmpz_poly_set(test_fmpz_poly2, test_fmpz_poly);
      _fmpz_poly_reverse(test_fmpz_poly, test_fmpz_poly, length2);
      fmpz_poly_check_normalisation(test_fmpz_poly);
      _fmpz_poly_reverse(test_fmpz_poly, test_fmpz_poly, length2);
      fmpz_poly_check_normalisation(test_fmpz_poly);
           
      result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }

   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_newton_invert()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly2);
      fmpz_poly_init(test_fmpz_poly3);
            
      length = random_ulong(250)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      _fmpz_poly_set_coeff_ui(test_fmpz_poly, test_fmpz_poly->length - 1, 1);
      length = test_fmpz_poly->length;
      
      _fmpz_poly_reverse(test_fmpz_poly, test_fmpz_poly, length);
            
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      fmpz_poly_newton_invert(test_fmpz_poly2, test_fmpz_poly, length);
      fmpz_poly_check_normalisation(test_fmpz_poly2);
      
      fmpz_poly_mul_trunc_n(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2, length);
      
      _fmpz_poly_normalise(test_fmpz_poly3);
            
      result &= (test_fmpz_poly3->length == 1);
      result &= (test_fmpz_poly3->coeffs[0] == 1);
      result &= (test_fmpz_poly3->coeffs[1] == 1);
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_div_series()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 500) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_fmpz_poly3);
      fmpz_poly_init(test_fmpz_poly4);
            
      length = random_ulong(200)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);
      
      _fmpz_poly_set_coeff_ui(test_fmpz_poly, test_fmpz_poly->length - 1, 1);
      length = test_fmpz_poly->length;
      
      randpoly(test_poly, length, bits); 
      fmpz_poly_fit_length(test_fmpz_poly2, length);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly);

      _fmpz_poly_reverse(test_fmpz_poly, test_fmpz_poly, length);
            
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      fmpz_poly_div_series(test_fmpz_poly3, test_fmpz_poly2, test_fmpz_poly, length);
      fmpz_poly_check_normalisation(test_fmpz_poly3);
      
      fmpz_poly_mul_trunc_n(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly, length);
      
      _fmpz_poly_normalise(test_fmpz_poly4);
            
      result = _fmpz_poly_equal(test_fmpz_poly4, test_fmpz_poly2);
      
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_div_newton()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(10)+ 1;
      bits2 = random_ulong(10)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128); 
      length = random_ulong(128)+1;
      //length = 100000;
      //length2 = 100000;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      randpoly(test_poly, length, bits); 
      mpz_poly_set_coeff_ui(test_poly, length - 1, 1);
      
      fmpz_poly_fit_length(test_fmpz_poly, length);
      mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      _fmpz_poly_normalise(test_fmpz_poly2);
      
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 1; i < 10; i++)
      {
         fmpz_poly_div_newton(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_div_newton(test_fmpz_poly4, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
      
#if DEBUG2
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_power()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, temp;
   int result = 1;
   unsigned long bits, length, exp;
   
   mpz_poly_init(test_poly); 
   
   fmpz_poly_init(temp);
   fmpz_poly_init(test_fmpz_poly2);
   fmpz_poly_init(test_fmpz_poly3);
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10); 
          exp = random_ulong(20);
#if DEBUG
          printf("length = %ld, bits = %ld, exp = %ld\n", length, bits, exp);
#endif
          do 
          {
            randpoly(test_poly, length, bits);
            mpz_poly_normalise(test_poly); 
          }
          while (test_poly->length != length);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_normalise(test_fmpz_poly);
            
          fmpz_poly_fit_length(test_fmpz_poly2, 1);
          fmpz_poly_fit_limbs(test_fmpz_poly2, 1);
          fmpz_poly_set_coeff_ui(test_fmpz_poly2, 0, 1);
          test_fmpz_poly2->length = 1;
          
          for (unsigned long i = 0; i < exp; i++)
          {
             fmpz_poly_mul(temp, test_fmpz_poly2, test_fmpz_poly);
             fmpz_poly_fit_length(test_fmpz_poly2, temp->length);
             fmpz_poly_fit_limbs(test_fmpz_poly2, temp->limbs);
             _fmpz_poly_set(test_fmpz_poly2, temp);
          }
          
          fmpz_poly_power(test_fmpz_poly3, test_fmpz_poly, exp);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly3);
                    
#if DEBUG2
          if (!result)
          {
             fmpz_poly_print(test_fmpz_poly); printf("\n");
             fmpz_poly_print(test_fmpz_poly2); printf("\n");
             fmpz_poly_print(test_fmpz_poly3); printf("\n");
          }
#endif
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10); 
          exp = random_ulong(20);
#if DEBUG
          printf("length = %ld, bits = %ld, exp = %ld\n", length, bits, exp);
#endif
          do 
          {
            randpoly(test_poly, length, bits);
            mpz_poly_normalise(test_poly); 
          }
          while (test_poly->length != length);
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_normalise(test_fmpz_poly);
            
          fmpz_poly_fit_length(test_fmpz_poly2, 1);
          fmpz_poly_fit_limbs(test_fmpz_poly2, 1);
          fmpz_poly_set_coeff_ui(test_fmpz_poly2, 0, 1);
          test_fmpz_poly2->length = 1;
          
          for (unsigned long i = 0; i < exp; i++)
          {
             fmpz_poly_mul(temp, test_fmpz_poly2, test_fmpz_poly);
             fmpz_poly_fit_length(test_fmpz_poly2, temp->length);
             fmpz_poly_fit_limbs(test_fmpz_poly2, temp->limbs);
             _fmpz_poly_set(test_fmpz_poly2, temp);
          }
          
          fmpz_poly_power(test_fmpz_poly, test_fmpz_poly, exp);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
                    
#if DEBUG2
          if (!result)
          {
             printf("Exp = %ld\n", exp);
             fmpz_poly_print(test_fmpz_poly); printf("\n");
             fmpz_poly_print(test_fmpz_poly2); printf("\n");
          }
#endif
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   fmpz_poly_clear(temp);
   fmpz_poly_clear(test_fmpz_poly2);
   fmpz_poly_clear(test_fmpz_poly3);
   
   return result; 
}

int test_fmpz_poly_power_trunc_n()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, temp;
   int result = 1;
   unsigned long bits, length, exp, n;
   
   mpz_poly_init(test_poly); 
   
   fmpz_poly_init(temp);
   fmpz_poly_init(test_fmpz_poly2);
   fmpz_poly_init(test_fmpz_poly3);
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10); 
          exp = random_ulong(20);
          n = random_ulong(20);
#if DEBUG
          printf("length = %ld, bits = %ld, exp = %ld\n", length, bits, exp);
#endif
          randpoly(test_poly, length, bits); 
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_normalise(test_fmpz_poly);
            
          fmpz_poly_fit_length(test_fmpz_poly2, 1);
          fmpz_poly_fit_limbs(test_fmpz_poly2, 1);
          fmpz_poly_set_coeff_ui(test_fmpz_poly2, 0, 1);
          test_fmpz_poly2->length = 1;
          
          for (unsigned long i = 0; i < exp; i++)
          {
             fmpz_poly_mul(temp, test_fmpz_poly2, test_fmpz_poly);
             fmpz_poly_fit_length(test_fmpz_poly2, temp->length);
             fmpz_poly_fit_limbs(test_fmpz_poly2, temp->limbs);
             _fmpz_poly_set(test_fmpz_poly2, temp);
          }
          _fmpz_poly_truncate(test_fmpz_poly2, n);
          if (test_fmpz_poly->length == 0) _fmpz_poly_zero(test_fmpz_poly2);
          
          fmpz_poly_power_trunc_n(test_fmpz_poly3, test_fmpz_poly, exp, n);
          fmpz_poly_check_normalisation(test_fmpz_poly3);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly3);
                    
#if DEBUG2
          if (!result)
          {
             fmpz_poly_print(test_fmpz_poly); printf("\n");
             fmpz_poly_print(test_fmpz_poly2); printf("\n");
             fmpz_poly_print(test_fmpz_poly3); printf("\n");
          }
#endif
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10); 
          exp = random_ulong(20);
          n = random_ulong(20);
#if DEBUG
          printf("length = %ld, bits = %ld, exp = %ld\n", length, bits, exp);
#endif
          randpoly(test_poly, length, bits); 
          
          fmpz_poly_fit_length(test_fmpz_poly, length);
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          _fmpz_poly_normalise(test_fmpz_poly);
            
          fmpz_poly_fit_length(test_fmpz_poly2, 1);
          fmpz_poly_fit_limbs(test_fmpz_poly2, 1);
          fmpz_poly_set_coeff_ui(test_fmpz_poly2, 0, 1);
          test_fmpz_poly2->length = 1;
          
          for (unsigned long i = 0; i < exp; i++)
          {
             fmpz_poly_mul(temp, test_fmpz_poly2, test_fmpz_poly);
             fmpz_poly_fit_length(test_fmpz_poly2, temp->length);
             fmpz_poly_fit_limbs(test_fmpz_poly2, temp->limbs);
             _fmpz_poly_set(test_fmpz_poly2, temp);
          }
          _fmpz_poly_truncate(test_fmpz_poly2, n);
          if (test_fmpz_poly->length == 0) _fmpz_poly_zero(test_fmpz_poly2);
          
          fmpz_poly_power_trunc_n(test_fmpz_poly, test_fmpz_poly, exp, n);
          fmpz_poly_check_normalisation(test_fmpz_poly);
          
          result = _fmpz_poly_equal(test_fmpz_poly2, test_fmpz_poly);
                    
#if DEBUG2
          if (!result)
          {
             fmpz_poly_print(test_fmpz_poly); printf("\n");
             fmpz_poly_print(test_fmpz_poly2); printf("\n");
          }
#endif
      }
      fmpz_poly_clear(test_fmpz_poly);
   }
   
   mpz_poly_clear(test_poly);
   fmpz_poly_clear(temp);
   fmpz_poly_clear(test_fmpz_poly2);
   fmpz_poly_clear(test_fmpz_poly3);
   
   return result; 
}

int test_fmpz_poly_power2()
{
    fmpz_poly_t poly, power;
    fmpz_poly_init(power);
    fmpz_poly_init(poly);
    fmpz_poly_set_coeff_ui(poly, 0, 743);
    fmpz_poly_set_coeff_ui(poly, 1, 423);
    fmpz_poly_power(power, poly, 2000);//(1UL<<13));
    fmpz_poly_check_normalisation(power);
#if DEBUG
    fmpz_poly_print(power); printf("\n");
#endif
    
    return 1;
}

int test_fmpz_poly_pseudo_divrem_cohen()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(300)+ 2;
      bits2 = random_ulong(300)+ 1;
      //bits = bits2 = 1000000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(40); 
      length = random_ulong(40)+1; 
      //length = length2 = 20;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_pseudo_divrem_cohen(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_pseudo_divrem_cohen(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = (test_fmpz_poly5->length == 0);//mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_divrem_shoup()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(300)+ 2;
      bits2 = random_ulong(300)+ 1;
      //bits = bits2 = 1000000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(40); 
      length = random_ulong(40)+1; 
      //length = length2 = 20;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_pseudo_divrem_shoup(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_pseudo_divrem_shoup(test_fmpz_poly4, test_fmpz_poly5, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = (test_fmpz_poly5->length == 0);//mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_divrem_basecase()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2, d;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 6000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 2;
      bits2 = random_ulong(100)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(100); 
      length = random_ulong(100)+1; 
      //length = 100;
      //length2 = 199;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_pseudo_divrem_basecase(test_fmpz_poly4, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_pseudo_divrem_basecase(test_fmpz_poly4, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_div_basecase()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, d;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 6000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 2;
      bits2 = random_ulong(100)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(100); 
      length = random_ulong(100)+1; 
      //length = 100;
      //length2 = 199;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_pseudo_div_basecase(test_fmpz_poly4, &d, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_pseudo_div_basecase(test_fmpz_poly4, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_divrem_recursive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2, d;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
      //length = 100;
      //length2 = 199;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_pseudo_divrem_recursive(test_fmpz_poly4, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly5);
         fmpz_poly_init(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly5);
      }
      fmpz_poly_pseudo_divrem_recursive(test_fmpz_poly4, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_divrem()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4, test_fmpz_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2, d;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_divrem(test_fmpz_poly4, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_divrem(test_fmpz_poly3, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly3);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly5);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_divrem(test_fmpz_poly, test_fmpz_poly5, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly5);
   }
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_divrem(test_fmpz_poly4, test_fmpz_poly3, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly3);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_divrem(test_fmpz_poly4, test_fmpz_poly, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_div_recursive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, d;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
      //length = 100;
      //length2 = 199;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_pseudo_div_recursive(test_fmpz_poly4, &d, test_fmpz_poly3, test_fmpz_poly);
         fmpz_poly_check_normalisation(test_fmpz_poly4);
         fmpz_poly_clear(test_fmpz_poly4);
         fmpz_poly_init(test_fmpz_poly4);
      }
      fmpz_poly_pseudo_div_recursive(test_fmpz_poly4, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_pseudo_div()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2, test_fmpz_poly3, test_fmpz_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, d;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      fmpz_poly_init(test_fmpz_poly4);
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_div(test_fmpz_poly4, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
      fmpz_poly_clear(test_fmpz_poly4);
   }
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_div(test_fmpz_poly3, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly3);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly3);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   for (unsigned long count1 = 0; (count1 < 600) && (result == 1) ; count1++)
   {
      bits = random_ulong(20)+ 2;
      bits2 = random_ulong(20)+ 1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length = random_ulong(300)+1; 
      length2 = random_ulong(300); 
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_fit_length(test_fmpz_poly, length);
         mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
         _fmpz_poly_normalise(test_fmpz_poly);
      } while (test_fmpz_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_fit_length(test_fmpz_poly2, length2);
      mpz_poly_to_fmpz_poly(test_fmpz_poly2, test_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_fmpz_poly3, length+length2-1, test_fmpz_poly->limbs+test_fmpz_poly2->limbs+1);
          
      _fmpz_poly_mul(test_fmpz_poly3, test_fmpz_poly, test_fmpz_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      fmpz_poly_pseudo_div(test_fmpz_poly, &d, test_fmpz_poly3, test_fmpz_poly);
      fmpz_poly_check_normalisation(test_fmpz_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_fmpz_poly);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
      fmpz_poly_clear(test_fmpz_poly3);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_to_ZmodF_poly()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, depth;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
          depth = 0;
          while ((1<<depth) < length) depth++;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly(test_poly, length, bits-1); 
          while (mpz_poly_length(test_poly) < length);
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, depth, (bits-1)/FLINT_BITS+1, 0);
          fmpz_poly_to_ZmodF_poly(test_modF_poly, test_fmpz_poly, length);
          ZmodF_poly_to_fmpz_poly(test_fmpz_poly2, test_modF_poly, 1);
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      } 
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_bit_pack()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, depth, bundle;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(FLINT_BITS-2)+ 2;
      
      fmpz_poly_init2(test_fmpz_poly, 1, 1);
      fmpz_poly_init2(test_fmpz_poly2, 1, 10);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
          bundle = length/5;
          if (bundle == 0) bundle = length;
          depth = 0;
          while ((1<<depth) < (length-1)/bundle + 1) depth++;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly(test_poly, length, bits-1);
          while (mpz_poly_length(test_poly) < length);
          
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
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, depth, (bits*bundle-1)/FLINT_BITS+1, 0);
          
          fmpz_poly_bit_pack(test_modF_poly, test_fmpz_poly, bundle, -bits, length, 1L);
          test_fmpz_poly2->length = length;
          
          for (unsigned long i = 0; i < length; i++) // Must clear coeffs in advance
             test_fmpz_poly2->coeffs[i*(test_fmpz_poly2->limbs+1)] = 0; 
             
          fmpz_poly_bit_unpack(test_fmpz_poly2, test_modF_poly, bundle, bits); 
          fmpz_poly_check_normalisation(test_fmpz_poly2); 
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_bit_pack_unsigned()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, depth, bundle;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(FLINT_BITS-2)+ 2;
      
      fmpz_poly_init2(test_fmpz_poly, 1, 1);
      fmpz_poly_init2(test_fmpz_poly2, 1, 10);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
          bundle = length/5;
          if (bundle == 0) bundle = length;
          depth = 0;
          while ((1<<depth) < (length-1)/bundle + 1) depth++;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly_unsigned(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, depth, (bits*bundle-1)/FLINT_BITS+1, 0);
          
          fmpz_poly_bit_pack(test_modF_poly, test_fmpz_poly, bundle, bits, length, 1L);
          test_fmpz_poly2->length = length;
          
          for (unsigned long i = 0; i < length; i++) // Must clear coeffs in advance
             test_fmpz_poly2->coeffs[i*(test_fmpz_poly2->limbs+1)] = 0; 
             
          fmpz_poly_bit_unpack_unsigned(test_fmpz_poly2, test_modF_poly, bundle, bits);
          fmpz_poly_check_normalisation(test_fmpz_poly2);  
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_limb_pack_unsigned()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   int result = 1;
   unsigned long bits, length, length2, depth, bundle, limbs;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      limbs = (bits-1)/FLINT_BITS + 1;
      fmpz_poly_init2(test_fmpz_poly, 1, limbs);
      fmpz_poly_init2(test_fmpz_poly2, 1, limbs);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
      
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly_unsigned(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, 0, length*limbs, 0);
          
          fmpz_poly_limb_pack(test_modF_poly, test_fmpz_poly, length, limbs);
                  
          fmpz_poly_limb_unpack_unsigned(test_fmpz_poly2, test_modF_poly, length, limbs);
          fmpz_poly_check_normalisation(test_fmpz_poly2);  
          test_fmpz_poly2->length = length;
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_limb_pack()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   int result = 1;
   unsigned long bits, length, length2, depth, bundle, limbs;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      limbs = (bits-1)/FLINT_BITS + 1;
      fmpz_poly_init2(test_fmpz_poly, 1, limbs);
      fmpz_poly_init2(test_fmpz_poly2, 1, limbs);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
      
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly(test_poly, length, bits-1);
          while (mpz_poly_length(test_poly) < length);

#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, 0, length*limbs, 0);
          
          fmpz_poly_limb_pack(test_modF_poly, test_fmpz_poly, length, limbs);
                  
          fmpz_poly_limb_unpack(test_fmpz_poly2, test_modF_poly, length, limbs);  
          fmpz_poly_check_normalisation(test_fmpz_poly2);
          test_fmpz_poly2->length = length;
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_byte_pack_unsigned()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, bytes;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 64;
      bytes = ((bits-1)>>3)+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly_unsigned(test_poly, length, bits/2);
          while (mpz_poly_length(test_poly) < length);

#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, 0, ((bytes*length-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);
          
          for (unsigned long j = 0; j < 100; j++)
          {
             fmpz_poly_byte_pack(test_modF_poly, test_fmpz_poly, length, bytes, length, 1L);
             test_fmpz_poly2->length = length;
          
             for (unsigned long i = 0; i < length; i++) // Must clear coeffs in advance
                test_fmpz_poly2->coeffs[i*(test_fmpz_poly2->limbs+1)] = 0; 
             
             fmpz_poly_byte_unpack_unsigned(test_fmpz_poly2, test_modF_poly->coeffs[0], length, bytes);  
             fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_byte_pack()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_fmpz_poly, test_fmpz_poly2;
   ZmodF_poly_t test_modF_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length, bytes;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 5) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 130;
      bytes = ((bits-1)>>3)+1;
      
      fmpz_poly_init2(test_fmpz_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_fmpz_poly2, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_fit_length(test_fmpz_poly, length);
          fmpz_poly_fit_length(test_fmpz_poly2, length);
          mpz_poly_realloc(test_poly2, length);
          
          do randpoly(test_poly, length, bits/2);
          while (mpz_poly_length(test_poly) < length);

#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_fmpz_poly, test_poly);
          ZmodF_poly_init(test_modF_poly, 0, ((bytes*length-1)>>FLINT_LG_BYTES_PER_LIMB)+1, 0);
          
          for (unsigned long j = 0; j < 100; j++)
          {
             fmpz_poly_byte_pack(test_modF_poly, test_fmpz_poly, length, bytes, length, 1L);
             test_fmpz_poly2->length = length;
          
             for (unsigned long i = 0; i < length; i++) // Must clear coeffs in advance
                test_fmpz_poly2->coeffs[i*(test_fmpz_poly2->limbs+1)] = 0; 
             
             fmpz_poly_byte_unpack(test_fmpz_poly2, test_modF_poly->coeffs[0], length, bytes);  
             fmpz_poly_check_normalisation(test_fmpz_poly2);
          }
          fmpz_poly_to_mpz_poly(test_poly2, test_fmpz_poly2);
          
          ZmodF_poly_clear(test_modF_poly);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
          if (!result) 
          {
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zx, ",test_poly->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zx, ",test_poly2->coeffs[j]);
          printf("\n\n");
          }
      }   
      fmpz_poly_clear(test_fmpz_poly);
      fmpz_poly_clear(test_fmpz_poly2);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

#if TESTFILE
   RUN_TEST(fmpz_poly_freadprint); 
#endif
   RUN_TEST(fmpz_poly_tofromstring); 
   RUN_TEST(fmpz_poly_to_ZmodF_poly); 
   RUN_TEST(fmpz_poly_bit_pack); 
   RUN_TEST(fmpz_poly_bit_pack_unsigned); 
   RUN_TEST(fmpz_poly_byte_pack_unsigned); 
   RUN_TEST(fmpz_poly_byte_pack); 
   RUN_TEST(fmpz_poly_limb_pack_unsigned); 
   RUN_TEST(fmpz_poly_limb_pack); 
   RUN_TEST(_fmpz_poly_attach);
   RUN_TEST(_fmpz_poly_attach_shift);
   RUN_TEST(_fmpz_poly_attach_truncate);
   RUN_TEST(_fmpz_poly_truncate);
   RUN_TEST(_fmpz_poly_max_bits);
   RUN_TEST(_fmpz_poly_max_bits1);
   RUN_TEST(_fmpz_poly_max_limbs);
   RUN_TEST(_fmpz_poly_convert);
   RUN_TEST(_fmpz_poly_getset_ui);
   RUN_TEST(fmpz_poly_getset_ui);
   RUN_TEST(_fmpz_poly_getset_si);
   RUN_TEST(fmpz_poly_getset_si);
   RUN_TEST(_fmpz_poly_get_coeff_ptr);
   RUN_TEST(fmpz_poly_get_coeff_ptr);
   RUN_TEST(_fmpz_poly_normalise);
   RUN_TEST(_fmpz_poly_getset_coeff);
   RUN_TEST(fmpz_poly_getset_coeff);
   RUN_TEST(_fmpz_poly_getset_coeff_fmpz);
   RUN_TEST(fmpz_poly_getset_coeff_fmpz);
   RUN_TEST(_fmpz_poly_getset_coeff_mpz);
   RUN_TEST(fmpz_poly_getset_coeff_mpz);
   RUN_TEST(_fmpz_poly_setequal);
   RUN_TEST(_fmpz_poly_zero_coeffs);
   RUN_TEST(fmpz_poly_zero_coeffs);
   RUN_TEST(fmpz_poly_swap);
   RUN_TEST(_fmpz_poly_reverse);
   RUN_TEST(_fmpz_poly_neg);
   RUN_TEST(_fmpz_poly_shift);
   RUN_TEST(_fmpz_poly_add);
   RUN_TEST(fmpz_poly_add);
   RUN_TEST(_fmpz_poly_sub);
   RUN_TEST(fmpz_poly_sub);
   RUN_TEST(_fmpz_poly_scalar_mul_ui);
   RUN_TEST(fmpz_poly_scalar_mul_ui);
   RUN_TEST(_fmpz_poly_scalar_mul_si);
   RUN_TEST(fmpz_poly_scalar_mul_si);
   RUN_TEST(_fmpz_poly_scalar_mul); 
   RUN_TEST(fmpz_poly_scalar_mul); 
   RUN_TEST(fmpz_poly_scalar_mul_mpz); 
   RUN_TEST(_fmpz_poly_scalar_div_exact_ui);
   RUN_TEST(_fmpz_poly_scalar_div_exact_si);
   RUN_TEST(_fmpz_poly_scalar_div_ui);
   RUN_TEST(_fmpz_poly_scalar_tdiv_ui);
   RUN_TEST(_fmpz_poly_scalar_div_si);
   RUN_TEST(_fmpz_poly_scalar_tdiv_si);
   RUN_TEST(_fmpz_poly_scalar_div); 
   RUN_TEST(fmpz_poly_scalar_div); 
   RUN_TEST(fmpz_poly_scalar_div_mpz); 
   RUN_TEST(_fmpz_poly_mul_classical);
   RUN_TEST(_fmpz_poly_mul_classical_trunc);
   RUN_TEST(_fmpz_poly_mul_classical_trunc_left);
   RUN_TEST(_fmpz_poly_mul_karatsuba);
   RUN_TEST(_fmpz_poly_mul_karatsuba_trunc);
   RUN_TEST(_fmpz_poly_mul_karatsuba_trunc_left);
   RUN_TEST(_fmpz_poly_mul_KS);
   RUN_TEST(_fmpz_poly_mul_KS_trunc);
   RUN_TEST(_fmpz_poly_mul_SS);
   RUN_TEST(_fmpz_poly_mul_SS_trunc);
   RUN_TEST(_fmpz_poly_mul);
   RUN_TEST(fmpz_poly_mul);
   RUN_TEST(_fmpz_poly_mul_trunc_n);
   RUN_TEST(fmpz_poly_mul_trunc_n);
   RUN_TEST(_fmpz_poly_mul_trunc_left_n);
   RUN_TEST(fmpz_poly_mul_trunc_left_n);
   RUN_TEST(fmpz_poly_div_classical);
   RUN_TEST(fmpz_poly_divrem_classical);
   RUN_TEST(fmpz_poly_div_divconquer_recursive);
   RUN_TEST(fmpz_poly_divrem_divconquer);
   RUN_TEST(fmpz_poly_div_divconquer);
   RUN_TEST(fmpz_poly_newton_invert_basecase);
   RUN_TEST(fmpz_poly_newton_invert);
   RUN_TEST(fmpz_poly_div_series);
   RUN_TEST(fmpz_poly_div_newton);
   RUN_TEST(fmpz_poly_div_mulders);
   RUN_TEST(fmpz_poly_divrem);
   RUN_TEST(fmpz_poly_div);
   RUN_TEST(fmpz_poly_pseudo_divrem_recursive); 
   RUN_TEST(fmpz_poly_pseudo_divrem_basecase); 
   RUN_TEST(fmpz_poly_pseudo_div_basecase); 
   RUN_TEST(fmpz_poly_pseudo_div_recursive); 
   RUN_TEST(fmpz_poly_pseudo_divrem_cohen); 
   RUN_TEST(fmpz_poly_pseudo_divrem_shoup); 
   RUN_TEST(fmpz_poly_pseudo_divrem); 
   RUN_TEST(fmpz_poly_pseudo_div); 
   RUN_TEST(fmpz_poly_power);
   RUN_TEST(fmpz_poly_power_trunc_n);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   fmpz_poly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();

   return 0;
}


