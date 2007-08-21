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

#define VARY_BITS 1
#define SIGNS 1

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

unsigned long randint(unsigned long randsup) 
{
    static unsigned long randval = 4035456057U;
    randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
    
    return (unsigned long)randval%randsup;
}

void randpoly(mpz_poly_t pol, unsigned long length, unsigned long maxbits)
{
   unsigned long bits;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_zero(pol);
   
   for (unsigned long i = 0; i < length; i++)
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
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
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
   
int test_fmpz_poly_convert()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 300) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          mpz_poly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits); 
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = mpz_poly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_clear(temp);
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_getset_ui()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long coeff, coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(randstate, coeff_bits);
              coeff_num = randint(length);
#if DEBUG
              printf("Index = %ld, bits = %ld, coeff = %ld\n", coeff_num, coeff_bits, coeff);
#endif
              _fmpz_poly_set_coeff_ui(test_mpn_poly, coeff_num, coeff);
              result = (_fmpz_poly_get_coeff_ui(test_mpn_poly, coeff_num) == coeff);
          }
      }
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_si()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
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
              _fmpz_poly_set_coeff_si(test_mpn_poly, coeff_num, coeff);
              result = ((_fmpz_poly_get_coeff_si(test_mpn_poly, coeff_num) == coeff) && (_fmpz_poly_get_coeff_ui(test_mpn_poly, coeff_num) == sign*coeff));
          }
      }
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_get_coeff_ptr()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
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
              _fmpz_poly_set_coeff_si(test_mpn_poly, coeff_num, coeff);
              if (coeff == 0) sign = 0;
              result = (_fmpz_poly_get_coeff_ptr(test_mpn_poly, coeff_num)[0] == sign);
          }
      }
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_normalise()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long nz_coeff;
   long sign;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          nz_coeff = randint(length+1)-1;
          if (randint(2)) sign = -1L; else sign = 1;
          if (nz_coeff != -1L) _fmpz_poly_set_coeff_si(test_mpn_poly, nz_coeff, sign*1000);
          for (unsigned long i = nz_coeff+1; i < length; i++)
            _fmpz_poly_set_coeff_ui(test_mpn_poly, i, 0);
              
          _fmpz_poly_normalise(test_mpn_poly);
#if DEBUG
          printf("length = %ld, nonzero coefficient = %ld\n",_fmpz_poly_length(test_mpn_poly), nz_coeff);
#endif              
          result = (_fmpz_poly_length(test_mpn_poly) == nz_coeff+1);
      }
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test_fmpz_poly_getset_coeff()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   mpz_poly_init(test_poly); 
           
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mp_limb_t * coeff1 = (mp_limb_t *) calloc(test_mpn_poly->limbs, sizeof(mp_limb_t));
          mp_limb_t * coeff2 = (mp_limb_t *) calloc(test_mpn_poly->limbs, sizeof(mp_limb_t));
          
          sign = _fmpz_poly_get_coeff(coeff1, test_mpn_poly, randint(test_mpn_poly->length));
          rand_coeff = randint(test_mpn_poly->length);
          _fmpz_poly_set_coeff(test_mpn_poly, rand_coeff, coeff1, sign, test_mpn_poly->limbs);
          sign2 = _fmpz_poly_get_coeff(coeff2, test_mpn_poly, rand_coeff);
          
          for (unsigned long i = 0; (i < test_mpn_poly->limbs) && (result == 1); i++)
          {
              result = (coeff1[i] == coeff2[i]);
          }
          
          free(coeff1);
          free(coeff2);
                    
          if (sign != sign2) result = 0;
      }
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_setequal()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long altered_coeff, extra_zeroes;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+1+randint(30));
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+extra_zeroes);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          for (unsigned long i = mpz_poly_length(test_poly); i < length + extra_zeroes; i++)
          {
             _fmpz_poly_set_coeff_ui(test_mpn_poly, i, 0);
          }
          test_mpn_poly->length = length;
          _fmpz_poly_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly); 
      }

      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+extra_zeroes);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          for (unsigned long i = mpz_poly_length(test_poly); i < length + extra_zeroes; i++)
          {
             _fmpz_poly_set_coeff_ui(test_mpn_poly, i, 0);
          }
          test_mpn_poly->length = length;
          _fmpz_poly_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          altered_coeff = randint(length);
          test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)+1]++;
          if (test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] == 0)
             test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] = 1;
          result = !_fmpz_poly_equal(test_mpn_poly2, test_mpn_poly); 
      }

      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+extra_zeroes);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          for (unsigned long i = mpz_poly_length(test_poly); i < length + extra_zeroes; i++)
          {
             _fmpz_poly_set_coeff_ui(test_mpn_poly, i, 0);
          }
          test_mpn_poly->length = length;
          _fmpz_poly_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          altered_coeff = randint(length);
          test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)]*=-1L;
          if (test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] == 0)
             test_mpn_poly2->coeffs[altered_coeff*(test_mpn_poly2->limbs+1)] = 1;
          
          result = !_fmpz_poly_equal(test_mpn_poly2, test_mpn_poly); 
      }
      
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);         
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_swap()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length, length2;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          length2 = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          fmpz_poly_realloc(test_mpn_poly2, length2);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length2, bits);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly);
          
          _fmpz_poly_set(test_mpn_poly3, test_mpn_poly);
          _fmpz_poly_swap(test_mpn_poly, test_mpn_poly2);
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly3);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test_fmpz_poly_shift()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long shift;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+1);
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          shift = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+shift);
          fmpz_poly_realloc(test_mpn_poly2, length+shift);
          
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          _fmpz_poly_set(test_mpn_poly2, test_mpn_poly);
          _fmpz_poly_left_shift(test_mpn_poly, test_mpn_poly, shift); 
          _fmpz_poly_right_shift(test_mpn_poly, test_mpn_poly, shift);
          
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly);
      }

      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          shift = randint(length);
          fmpz_poly_realloc(test_mpn_poly, length+shift);
          fmpz_poly_realloc(test_mpn_poly2, length+shift);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          test_mpn_poly3->limbs = test_mpn_poly->limbs;
          test_mpn_poly3->length = test_mpn_poly->length-shift;
          test_mpn_poly3->coeffs = test_mpn_poly->coeffs+shift*(test_mpn_poly->limbs+1);
          _fmpz_poly_right_shift(test_mpn_poly2, test_mpn_poly, shift);
                
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly2);
      }

      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test_fmpz_poly_neg()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length, check_coeff;
   unsigned long extra_bits1, extra_bits2;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      extra_bits1 = randint(200);
      extra_bits2 = randint(200);
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits+extra_bits1-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits+extra_bits1+extra_bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;      
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          
          check_coeff = randint(length);

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          _fmpz_poly_neg(test_mpn_poly2, test_mpn_poly);
          _fmpz_poly_neg(test_mpn_poly3, test_mpn_poly2);
          
          result = _fmpz_poly_equal(test_mpn_poly, test_mpn_poly3) 
             && (test_mpn_poly->coeffs[check_coeff*(test_mpn_poly->limbs+1)] == -test_mpn_poly2->coeffs[check_coeff*(test_mpn_poly2->limbs+1)]);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}


int test_fmpz_poly_add()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
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
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits3-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(10)+1; 
          length2 = random_ulong(10)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          max_length = (length > length2) ? length : length2;
          fmpz_poly_realloc(test_mpn_poly3, max_length);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          mpz_poly_add(test_poly3, test_poly, test_poly2);

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_add(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_normalise(test_mpn_poly3);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly3);

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
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }

   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = bits+random_ulong(200) + 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          length2 = length + random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          mpz_poly_add(test_poly3, test_poly, test_poly2);
          
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_add(test_mpn_poly2, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_normalise(test_mpn_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly2);
          
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
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test_fmpz_poly_sub()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
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
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits3-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          length2 = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          max_length = (length > length2) ? length : length2;
          fmpz_poly_realloc(test_mpn_poly3, max_length);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          mpz_poly_sub(test_poly3, test_poly, test_poly2);

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_sub(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_normalise(test_mpn_poly3);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly3);
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
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000) + 1;
      bits2 = bits+random_ulong(200) + 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          length2 = length + random_ulong(1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          mpz_poly_sub(test_poly3, test_poly, test_poly2);
          
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_sub(test_mpn_poly2, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_normalise(test_mpn_poly2);
          mpz_poly_init(test_poly4);
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly2);
          
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
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test_fmpz_poly_scalar_mul_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_mul_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length, sign;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = (long) randint(34682739);
          sign = randint(2);
          if (sign) mult = -mult;
          
          _fmpz_poly_scalar_mul_si(test_mpn_poly2, test_mpn_poly, mult);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_si(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}


int test_fmpz_poly_scalar_div_exact_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_mpn_poly3, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_mpn_poly2, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS+3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739) + 1;
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_mpn_poly3, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_scalar_div_exact_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   long mult;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_mpn_poly3, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_mpn_poly2, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS+3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          mult = randint(34682739) + 1;
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_mpn_poly3, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_scalar_div_ui()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
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
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          //length = 10000;       
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          div = randint(34682739)+1;
          for (unsigned long i = 0; i < 100; i++)
          {
             _fmpz_poly_scalar_div_ui(test_mpn_poly2, test_mpn_poly, div);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_div_si()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   long div;
   mpz_t temp;
   mpz_init(temp);
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 400) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          div = randint(34682739)+1;
          if (randint(2)) div = -div;
          _fmpz_poly_scalar_div_si(test_mpn_poly2, test_mpn_poly, div);
         
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              if (div < 0)
              {
                 mpz_tdiv_q_ui(temp, test_poly->coeffs[i], -div);
                 mpz_neg(temp, temp);
              } else
                 mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          mpz_poly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_mul_naive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_naive(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly3); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_naive_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(100)+1; 
          length2 = random_ulong(100)+1; 
          trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_mpn_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_naive(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_truncate(test_mpn_poly3, trunc);
          
          _fmpz_poly_mul_naive_trunc(test_mpn_poly4, test_mpn_poly, test_mpn_poly2, trunc);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_mpn_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3);
          fmpz_poly_clear(test_mpn_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_karatsuba()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2, max_length, log_length, output_bits;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 10000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(35)+1; 
      length = random_ulong(35)+1;   
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
          
      fmpz_poly_realloc(test_mpn_poly, length);
      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n");
      mpz_poly_print(test_poly2);printf("\n");
#endif          
      mpz_poly_mul_naive(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1); //(output_bits-1)/FLINT_BITS+1);
          
      for (unsigned long count2 = 0; (count2 < 1) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_karatsuba(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly3); 
#if DEBUG
          mpz_poly_print(test_poly3);printf("\n");
          mpz_poly_print(test_poly4);printf("\n");
#endif          
          result = mpz_poly_equal(test_poly4, test_poly3);
      }
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test_fmpz_poly_mul_karatsuba_trunc()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   
   for (unsigned long count1 = 1; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(2000)+ 1;
      bits2 = random_ulong(2000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(35)+1; 
          length2 = random_ulong(35)+1; 
          trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          
          do randpoly(test_poly2, length2, bits2); 
          while (mpz_poly_length(test_poly2) < length2);

          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_mpn_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_naive(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_truncate(test_mpn_poly3, trunc);
          
          _fmpz_poly_mul_karatsuba_trunc(test_mpn_poly4, test_mpn_poly, test_mpn_poly2, trunc);
          
          result = _fmpz_poly_equal(test_mpn_poly4, test_mpn_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          fmpz_poly_clear(test_mpn_poly3); 
          fmpz_poly_clear(test_mpn_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   
   return result; 
}


int test_fmpz_poly_mul_KS()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
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
      length = random_ulong(100)+1;
      length2 = random_ulong(100)+1;
      //length = 32000;
      //length2 = 32000;   
   
      _fmpz_poly_stack_init(test_mpn_poly, length, (bits-1)/FLINT_BITS+1);
      _fmpz_poly_stack_init(test_mpn_poly2, length2, (bits2-1)/FLINT_BITS+1);
      _fmpz_poly_stack_init(test_mpn_poly3, length + length2 - 1, test_mpn_poly->limbs + test_mpn_poly2->limbs + 1);
#if DEBUG
      printf("%ld, %ld, %ld, %ld\n", length, length2, bits, bits2);
#endif
      mpz_poly_realloc(test_poly, length);
      mpz_poly_realloc(test_poly2, length2);
      mpz_poly_realloc(test_poly3, length + length2 - 1);
      mpz_poly_realloc(test_poly4, length + length2 - 1);
          
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
         mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
         mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
         for (unsigned long i = 0; i < 10; i++)
            _fmpz_poly_mul_KS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
         fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly3);
                  
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
   
      _fmpz_poly_stack_clear(test_mpn_poly3);
      _fmpz_poly_stack_clear(test_mpn_poly2);
      _fmpz_poly_stack_clear(test_mpn_poly);
   
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   mpz_poly_clear(test_poly4);
   
   return result;
}

int test_fmpz_poly_mul_KS_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      bits2 = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(2000)+1; 
          length2 = random_ulong(2000)+1; 
          trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_mpn_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_truncate(test_mpn_poly3, trunc);
          _fmpz_poly_normalise(test_mpn_poly3);
          
          _fmpz_poly_mul_KS_trunc(test_mpn_poly4, test_mpn_poly, test_mpn_poly2, trunc);
          _fmpz_poly_normalise(test_mpn_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_mpn_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3); 
          fmpz_poly_clear(test_mpn_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
#if 1
   for (unsigned long count1 = 1; (count1 < 10) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(500)+1; 
          length2 = random_ulong(500)+1; 
          trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_mpn_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_truncate(test_mpn_poly3, trunc);
          _fmpz_poly_normalise(test_mpn_poly3);
          
          _fmpz_poly_mul_KS_trunc(test_mpn_poly4, test_mpn_poly, test_mpn_poly2, trunc);
          _fmpz_poly_normalise(test_mpn_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_mpn_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3); 
          fmpz_poly_clear(test_mpn_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
#endif

   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}


int test_fmpz_poly_mul_SS()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
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
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(1000)+1; 
      length = random_ulong(1000)+1; 
      //length = length2 = 256;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do randpoly(test_poly, length, bits);
      while (mpz_poly_length(test_poly) < length);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
          
      fmpz_poly_realloc(test_mpn_poly, length);
      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
      mpz_poly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_SS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
      }
      
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly3);
           
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
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_SS_trunc()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 30) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          length2 = random_ulong(1000)+1; 
          trunc = random_ulong(length+length2);       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_mpn_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_SS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_truncate(test_mpn_poly3, trunc);
          _fmpz_poly_normalise(test_mpn_poly3);
          
          _fmpz_poly_mul_SS_trunc(test_mpn_poly4, test_mpn_poly, test_mpn_poly2, trunc);
          _fmpz_poly_normalise(test_mpn_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_mpn_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4); 

          result = mpz_poly_equal(test_poly4, test_poly3);
#if DEBUG          
          if (!result)
          {
             mpz_poly_print(test_poly3); printf("\n\n");
             mpz_poly_print(test_poly4);
          }
#endif          
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3); 
          fmpz_poly_clear(test_mpn_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_trunc_n()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2, trunc;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 25) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 1;
      bits2 = random_ulong(1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = random_ulong(1000)+1; 
          length2 = length; 
          trunc = length;       
#if DEBUG
          printf("length = %ld, length2 = %ld, trunc = %ld, bits = %ld, bits2 = %ld\n", length, length2, trunc, bits, bits2);
#endif
          do randpoly(test_poly, length, bits);
          while (mpz_poly_length(test_poly) < length);

          do randpoly(test_poly2, length2, bits2);
          while (mpz_poly_length(test_poly2) < length2);
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
          
          mpz_poly_init(test_poly4);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS+2);
          fmpz_poly_init2(test_mpn_poly4, trunc, (bits+bits2-1)/FLINT_BITS+2);
          
          _fmpz_poly_mul_KS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          _fmpz_poly_truncate(test_mpn_poly3, trunc);
          _fmpz_poly_normalise(test_mpn_poly3);
          
          _fmpz_poly_mul_trunc_n(test_mpn_poly4, test_mpn_poly, test_mpn_poly2, trunc);
          _fmpz_poly_normalise(test_mpn_poly4);
          
          fmpz_poly_to_mpz_poly(test_poly3, test_mpn_poly3); 
          fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4); 
          
          result = mpz_poly_equal(test_poly4, test_poly3);
          mpz_poly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3);
          fmpz_poly_clear(test_mpn_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_scalar_mul()
{
   mpz_poly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
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
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+(bits2-1)/FLINT_BITS+2);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 3) && (result == 1); count2++)
      { 
          length = randint(100)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          
          clear_limbs(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          do randpoly(test_poly, length, bits); 
          while (mpz_poly_length(test_poly) < length);
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          for (unsigned long j = 0; j < 1; j++)
          {
            _fmpz_poly_scalar_mul(test_mpn_poly2, test_mpn_poly, x);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
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
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = randint(1000) + 1500;
      bits2 = randint(1000) + 1500;
      limbs2 = (bits2-1)/FLINT_BITS+1;
      bits2 = limbs2*FLINT_BITS;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+(bits2-1)/FLINT_BITS+2);
      x = (mp_limb_t*) malloc(sizeof(mp_limb_t)*(limbs2+1));
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = randint(100)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld, bits2 = %ld\n",length, bits, bits2);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          
          clear_limbs(x, limbs2+1);
          mpn_random2(x+1, limbs2);
          if (randint(2)) 
              x[0] = limbs2;
          else x[0] = -limbs2;
          
          randpoly(test_poly, length, bits); 
          mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
          
          for (unsigned long j = 0; j < 5; j++)
          {
            _fmpz_poly_scalar_mul(test_mpn_poly2, test_mpn_poly, x);
          }
          
          mpz_poly_init(test_poly2);
          fmpz_poly_to_mpz_poly(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
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
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }

   mpz_poly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_div_naive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 1000000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(100)+1; 
      length = random_ulong(100)+1; 
      //length = length2 = 20;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      fmpz_poly_init(test_mpn_poly4);
      _fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_div_naive(test_mpn_poly4, test_mpn_poly3, test_mpn_poly);
         fmpz_poly_clear(test_mpn_poly4);
         fmpz_poly_init(test_mpn_poly4);
      }
      fmpz_poly_div_naive(test_mpn_poly4, test_mpn_poly3, test_mpn_poly);
      _fmpz_poly_normalise(test_mpn_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_divrem_naive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4, test_mpn_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 1000) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 1000000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(100)+1; 
      length = random_ulong(100)+1; 
      //length = length2 = 20;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);

      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      fmpz_poly_init(test_mpn_poly4);
      fmpz_poly_init(test_mpn_poly5);
      _fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      
      for (unsigned long i = 1; i < 5; i++)
      {
         fmpz_poly_divrem_naive(test_mpn_poly4, test_mpn_poly5, test_mpn_poly3, test_mpn_poly);
         fmpz_poly_clear(test_mpn_poly4);
         fmpz_poly_clear(test_mpn_poly5);
         fmpz_poly_init(test_mpn_poly4);
         fmpz_poly_init(test_mpn_poly5);
      }
      fmpz_poly_divrem_naive(test_mpn_poly4, test_mpn_poly5, test_mpn_poly3, test_mpn_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly4);printf("\n\n");
#endif               
      
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
      fmpz_poly_clear(test_mpn_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_div_karatsuba_recursive()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4, test_mpn_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 150) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128)+1; 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);

      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      fmpz_poly_init(test_mpn_poly4);
      fmpz_poly_init(test_mpn_poly5);
      _fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 0; i < 10; i++)
      {
         fmpz_poly_div_karatsuba_recursive(test_mpn_poly4, test_mpn_poly5, test_mpn_poly3, test_mpn_poly);
         fmpz_poly_clear(test_mpn_poly4);
         fmpz_poly_clear(test_mpn_poly5);
         fmpz_poly_init(test_mpn_poly4);
         fmpz_poly_init(test_mpn_poly5);
      }
      fmpz_poly_div_karatsuba_recursive(test_mpn_poly4, test_mpn_poly5, test_mpn_poly3, test_mpn_poly);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
      fmpz_poly_clear(test_mpn_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_divrem_karatsuba()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4, test_mpn_poly5;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 150) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128)+1; 
      length = random_ulong(128)+1;
      //length = 12;
      //length2 = 5;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      do randpoly(test_poly2, length2, bits2);
      while (mpz_poly_length(test_poly2) < length2);
      
      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      fmpz_poly_init(test_mpn_poly4);
      fmpz_poly_init(test_mpn_poly5);
      _fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 0; i < 10; i++)
      {
         fmpz_poly_divrem_karatsuba(test_mpn_poly4, test_mpn_poly5, test_mpn_poly3, test_mpn_poly);
         fmpz_poly_clear(test_mpn_poly4);
         fmpz_poly_clear(test_mpn_poly5);
         fmpz_poly_init(test_mpn_poly4);
         fmpz_poly_init(test_mpn_poly5);
      }
      fmpz_poly_divrem_karatsuba(test_mpn_poly4, test_mpn_poly5, test_mpn_poly3, test_mpn_poly);
      _fmpz_poly_normalise(test_mpn_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
      fmpz_poly_clear(test_mpn_poly5);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_div_karatsuba()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 150) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      bits2 = random_ulong(1000)+ 1;
      //bits = bits2 = 10;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128)+1; 
      length = random_ulong(128)+1;
      //length = 1000;
      //length2 = 1000;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      fmpz_poly_init(test_mpn_poly4);
      _fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 1; i < 10; i++)
      {
         fmpz_poly_div_karatsuba(test_mpn_poly4, test_mpn_poly3, test_mpn_poly);
         fmpz_poly_clear(test_mpn_poly4);
         fmpz_poly_init(test_mpn_poly4);
      }
      fmpz_poly_div_karatsuba(test_mpn_poly4, test_mpn_poly3, test_mpn_poly);
      _fmpz_poly_normalise(test_mpn_poly4);
      
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_newton_invert_basecase()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length, n;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = random_ulong(1000)+ 2;
      //bits = 100000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_mpn_poly2);
      fmpz_poly_init(test_mpn_poly3);
      
      
      length = random_ulong(128)+1;
      //length = 12;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      _fmpz_poly_set_coeff_ui(test_mpn_poly, test_mpn_poly->length - 1, 1);
      
      n = randint(test_mpn_poly->length) + 1;
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      fmpz_poly_newton_invert_basecase(test_mpn_poly2, test_mpn_poly, n);
      
      fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
           
      for (unsigned long i = 0; i < n - 1; i++)
      {
          result &= (test_mpn_poly3->coeffs[(i+test_mpn_poly3->length-n)*(test_mpn_poly3->limbs+1)] == 0);
      }
      result &= (test_mpn_poly3->coeffs[(test_mpn_poly3->length-1)*(test_mpn_poly3->limbs+1)] == 1);
      result &= (test_mpn_poly3->coeffs[(test_mpn_poly3->length-1)*(test_mpn_poly3->limbs+1)+1] == 1);
      
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_reverse()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length, length2;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 5000) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+1);   
      
      length = random_ulong(100)+1;
      length2 = length + randint(200);
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld\n", length, length2, bits);
#endif

      randpoly(test_poly, length, bits); 
      fmpz_poly_realloc(test_mpn_poly, length);
      mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
      _fmpz_poly_normalise(test_mpn_poly);
      
      fmpz_poly_realloc(test_mpn_poly2, length2);
      
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      _fmpz_poly_reverse(test_mpn_poly2, test_mpn_poly, length2);
      _fmpz_poly_reverse(test_mpn_poly2, test_mpn_poly2, length2);
           
      result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly);
      
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_newton_invert()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_mpn_poly2);
      fmpz_poly_init(test_mpn_poly3);
            
      length = random_ulong(250)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      _fmpz_poly_set_coeff_ui(test_mpn_poly, test_mpn_poly->length - 1, 1);
      length = test_mpn_poly->length;
      
      _fmpz_poly_reverse(test_mpn_poly, test_mpn_poly, length);
            
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      fmpz_poly_newton_invert(test_mpn_poly2, test_mpn_poly, length);
      
      fmpz_poly_mul_trunc_n(test_mpn_poly3, test_mpn_poly, test_mpn_poly2, length);
      
      _fmpz_poly_normalise(test_mpn_poly3);
            
      result &= (test_mpn_poly3->length == 1);
      result &= (test_mpn_poly3->coeffs[0] == 1);
      result &= (test_mpn_poly3->coeffs[1] == 1);
      
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_div_series()
{
   mpz_poly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, length;
   
   mpz_poly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 50) && (result == 1) ; count1++)
   {
      bits = random_ulong(100)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init(test_mpn_poly3);
      fmpz_poly_init(test_mpn_poly4);
            
      length = random_ulong(200)+1;
       
#if DEBUG
      printf("length = %ld, bits = %ld\n", length, bits);
#endif

      do {
         randpoly(test_poly, length, bits); 
         fmpz_poly_realloc(test_mpn_poly, length);
         mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
         _fmpz_poly_normalise(test_mpn_poly);
      } while (test_mpn_poly->length == 0);
      
      _fmpz_poly_set_coeff_ui(test_mpn_poly, test_mpn_poly->length - 1, 1);
      length = test_mpn_poly->length;
      
      randpoly(test_poly, length, bits); 
      fmpz_poly_realloc(test_mpn_poly2, length);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly);

      _fmpz_poly_reverse(test_mpn_poly, test_mpn_poly, length);
            
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
#endif          
          
      fmpz_poly_div_series(test_mpn_poly3, test_mpn_poly2, test_mpn_poly, length);
      
      fmpz_poly_mul_trunc_n(test_mpn_poly4, test_mpn_poly3, test_mpn_poly, length);
      
      _fmpz_poly_normalise(test_mpn_poly4);
            
      result = _fmpz_poly_equal(test_mpn_poly4, test_mpn_poly2);
      
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
   }
   
   mpz_poly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_div_newton()
{
   mpz_poly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3, test_mpn_poly4;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   mpz_poly_init(test_poly); 
   mpz_poly_init(test_poly2); 
   mpz_poly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = random_ulong(10)+ 1;
      bits2 = random_ulong(10)+ 1;
      //bits = bits2 = 100000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS+1);
      
      length2 = random_ulong(128)+1; 
      length = random_ulong(128)+1;
      //length = 100000;
      //length2 = 100000;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      randpoly(test_poly, length, bits); 
      mpz_poly_set_coeff_ui(test_poly, length - 1, 1);
      
      fmpz_poly_realloc(test_mpn_poly, length);
      mpz_poly_to_fmpz_poly(test_mpn_poly, test_poly);
      
      randpoly(test_poly2, length2, bits2);     
      fmpz_poly_realloc(test_mpn_poly2, length2);
      mpz_poly_to_fmpz_poly(test_mpn_poly2, test_poly2);
      _fmpz_poly_normalise(test_mpn_poly2);
      
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
#endif          
          
      mpz_poly_init(test_poly4);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      fmpz_poly_init(test_mpn_poly4);
      _fmpz_poly_mul(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
      
#if DEBUG
      mpz_poly_print(test_poly);printf("\n\n");
      mpz_poly_print(test_poly2);printf("\n\n");
      mpz_poly_print(test_poly3);printf("\n\n");
#endif               
      for (unsigned long i = 1; i < 10; i++)
      {
         fmpz_poly_div_newton(test_mpn_poly4, test_mpn_poly3, test_mpn_poly);
         fmpz_poly_clear(test_mpn_poly4);
         fmpz_poly_init(test_mpn_poly4);
      }
      fmpz_poly_div_newton(test_mpn_poly4, test_mpn_poly3, test_mpn_poly);
      _fmpz_poly_check(test_mpn_poly4);
      
      _fmpz_poly_normalise(test_mpn_poly4);
      fmpz_poly_to_mpz_poly(test_poly4, test_mpn_poly4);
           
      result = mpz_poly_equal(test_poly4, test_poly2);
      
#if DEBUG
      if (!result) 
      {
         mpz_poly_print(test_poly2);printf("\n\n");
         mpz_poly_print(test_poly4);printf("\n\n");
      }
#endif               
            
      mpz_poly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
      fmpz_poly_clear(test_mpn_poly4);
   }
   
   mpz_poly_clear(test_poly);
   mpz_poly_clear(test_poly2);
   mpz_poly_clear(test_poly3);
   
   return result; 
}

void fmpz_poly_test_all()
{
   int success, all_success = 1;

   RUN_TEST(fmpz_poly_convert);
   RUN_TEST(fmpz_poly_getset_ui);
   RUN_TEST(fmpz_poly_getset_si);
   RUN_TEST(fmpz_poly_get_coeff_ptr);
   RUN_TEST(fmpz_poly_normalise);
   RUN_TEST(fmpz_poly_getset_coeff);
   RUN_TEST(fmpz_poly_setequal);
   RUN_TEST(fmpz_poly_swap);
   RUN_TEST(fmpz_poly_neg);
   RUN_TEST(fmpz_poly_shift);
   RUN_TEST(fmpz_poly_add);
   RUN_TEST(fmpz_poly_sub);
   RUN_TEST(fmpz_poly_scalar_mul_ui);
   RUN_TEST(fmpz_poly_scalar_mul_si);
   RUN_TEST(fmpz_poly_scalar_div_exact_ui);
   RUN_TEST(fmpz_poly_scalar_div_exact_si);
   RUN_TEST(fmpz_poly_scalar_div_ui);
   RUN_TEST(fmpz_poly_scalar_div_si);
   RUN_TEST(fmpz_poly_mul_naive);
   RUN_TEST(fmpz_poly_mul_naive_trunc);
   RUN_TEST(fmpz_poly_mul_karatsuba);
   RUN_TEST(fmpz_poly_mul_karatsuba_trunc);
   RUN_TEST(fmpz_poly_mul_KS);
   RUN_TEST(fmpz_poly_mul_KS_trunc);
   RUN_TEST(fmpz_poly_mul_SS);
   RUN_TEST(fmpz_poly_mul_SS_trunc);
   RUN_TEST(fmpz_poly_mul_trunc_n);
   RUN_TEST(fmpz_poly_scalar_mul);
   RUN_TEST(fmpz_poly_div_naive);
   RUN_TEST(fmpz_poly_divrem_naive);
   RUN_TEST(fmpz_poly_div_karatsuba_recursive);
   RUN_TEST(fmpz_poly_divrem_karatsuba);
   RUN_TEST(fmpz_poly_div_karatsuba);
   RUN_TEST(fmpz_poly_newton_invert_basecase);
   RUN_TEST(fmpz_poly_reverse);
   RUN_TEST(fmpz_poly_newton_invert);
   RUN_TEST(fmpz_poly_div_series);
   RUN_TEST(fmpz_poly_div_newton);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   fmpz_poly_test_all();
   test_support_cleanup();

   return 0;
}


