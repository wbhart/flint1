/****************************************************************************

fmpz_poly-test.c: Test code for fmpz_poly.c and fmpz_poly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz_poly.h"
#include "Zpoly.h"
#include "memory-manager.h"
#include "ZmodFpoly.h"

#define VARY_BITS 1
#define SIGNS 1

#define DEBUG 0 // prints debug information
#define DEBUG2 1 

gmp_randstate_t Zpoly_test_randstate;

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
#if VARY_BITS
       bits = randint(maxbits);
#else
       bits = maxbits;
#endif
       if (bits == 0) mpz_set_ui(temp,0);
       else 
       {
          mpz_rrandomb(temp, Zpoly_test_randstate, bits);
#if SIGNS
          if (randint(2)) mpz_neg(temp,temp);
#endif
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
   
int test_fmpz_poly_convert()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly;
   mpz_t temp;
   mpz_init(temp);
   int result = 1;
   unsigned long bits, length;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   for (unsigned long count1 = 1; (count1 < 100) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;
#if DEBUG
          printf("%ld, %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          Zpoly_realloc(test_poly2, length);
          randpoly(test_poly, length, bits); 
#if DEBUG
          for (unsigned j = 0; j < test_poly->length; j++)
             gmp_printf("%Zd, ",test_poly->coeffs[j]);
          printf("\n\n");
#endif
          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_convert_out(test_poly2, test_mpn_poly);
#if DEBUG
          for (unsigned j = 0; j < test_poly2->length; j++)
             gmp_printf("%Zd, ",test_poly2->coeffs[j]);
          printf("\n\n");
#endif
          
          result = Zpoly_equal(test_poly, test_poly2);
      }   
      fmpz_poly_clear(test_mpn_poly);
   }
   
   mpz_clear(temp);
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   
   return result;
}

int test_fmpz_poly_getset_ui()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long coeff, coeff_bits, coeff_num;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS_PER_LIMB);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(Zpoly_test_randstate, coeff_bits);
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
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_si()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS_PER_LIMB-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(Zpoly_test_randstate, coeff_bits);
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
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_get_coeff_ptr()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   long coeff, sign;
   unsigned long coeff_bits, coeff_num;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          for (unsigned long count3 = 1; (count3 < 1000) && result == 1; count3++)
          {
              coeff_bits = randint(FLINT_BITS_PER_LIMB-1);
              if (coeff_bits == 0) coeff = 0;
              else coeff = gmp_urandomb_ui(Zpoly_test_randstate, coeff_bits);
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
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_normalise()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length;
   unsigned long nz_coeff;
   long sign;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
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
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_getset_coeff()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly;
   int result = 1;
   unsigned long bits, length, rand_coeff;
   long sign, sign2;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mp_limb_t * coeff1 = (mp_limb_t *) calloc(test_mpn_poly->limbs,sizeof(mp_limb_t));
          mp_limb_t * coeff2 = (mp_limb_t *) calloc(test_mpn_poly->limbs,sizeof(mp_limb_t));
          
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
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_setequal()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long altered_coeff, extra_zeroes;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+1+randint(30));
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+extra_zeroes);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          for (unsigned long i = 0; i < extra_zeroes; i++)
          {
             _fmpz_poly_set_coeff_ui(test_mpn_poly, length+i, 0);
          }
          test_mpn_poly->length = length;
          _fmpz_poly_set(test_mpn_poly2, test_mpn_poly);
          test_mpn_poly->length = length+extra_zeroes;
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly); 
      }
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+extra_zeroes);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          for (unsigned long i = 0; i < extra_zeroes; i++)
          {
             _fmpz_poly_set_coeff_ui(test_mpn_poly, length+i, 0);
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
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          extra_zeroes = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+extra_zeroes);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          for (unsigned long i = 0; i < extra_zeroes; i++)
          {
             _fmpz_poly_set_coeff_ui(test_mpn_poly, length+i, 0);
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
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_swap()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length, length2;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
          length2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length2, bits); 

          _fmpz_poly_convert_in(test_mpn_poly2, test_poly);
          
          _fmpz_poly_set(test_mpn_poly3, test_mpn_poly);
          _fmpz_poly_swap(test_mpn_poly, test_mpn_poly2);
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly3);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_shift()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long shift;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          shift = randint(100);
          fmpz_poly_realloc(test_mpn_poly, length+shift);
          fmpz_poly_realloc(test_mpn_poly2, length+shift);
          
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_set(test_mpn_poly2, test_mpn_poly);
          _fmpz_poly_left_shift(test_mpn_poly, test_mpn_poly, shift); 
          _fmpz_poly_right_shift(test_mpn_poly, test_mpn_poly, shift);
          
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly);
      }
      
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          shift = randint(length);
          fmpz_poly_realloc(test_mpn_poly, length+shift);
          fmpz_poly_realloc(test_mpn_poly2, length+shift);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          test_mpn_poly3->limbs = test_mpn_poly->limbs;
          test_mpn_poly3->length = test_mpn_poly->length-shift;
          test_mpn_poly3->coeffs = test_mpn_poly->coeffs+shift*(test_mpn_poly->limbs+1);
          _fmpz_poly_right_shift(test_mpn_poly2, test_mpn_poly, shift);
                
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly2);
      }
     
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_negate()
{
   Zpoly_t test_poly;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length, check_coeff;
   unsigned long extra_bits1, extra_bits2;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      extra_bits1 = randint(200);
      extra_bits2 = randint(200);
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits+extra_bits1-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits+extra_bits1+extra_bits2-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;      
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 
          
          check_coeff = randint(length);

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_negate(test_mpn_poly2, test_mpn_poly);
          _fmpz_poly_negate(test_mpn_poly3, test_mpn_poly2);
          
          result = _fmpz_poly_equal(test_mpn_poly, test_mpn_poly3) 
             && (test_mpn_poly->coeffs[check_coeff*(test_mpn_poly->limbs+1)] == -test_mpn_poly2->coeffs[check_coeff*(test_mpn_poly2->limbs+1)]);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_add()
{
   Zpoly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, max_length;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   Zpoly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      bits2 = bits+gmp_urandomm_ui(Zpoly_test_randstate,200);
      bits3 = bits2+gmp_urandomm_ui(Zpoly_test_randstate,200)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits3-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
          length2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          max_length = (length > length2) ? length : length2;
          fmpz_poly_realloc(test_mpn_poly3, max_length);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          Zpoly_add(test_poly3, test_poly, test_poly2);

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_add(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          Zpoly_init3(test_poly4, max_length, bits3+1);
          _fmpz_poly_convert_out(test_poly4, test_mpn_poly3);
#if DEBUG
          Zpoly_print(stdout, test_poly); printf("\n");
          Zpoly_print(stdout, test_poly2); printf("\n");
          Zpoly_print(stdout, test_poly3); printf("\n");
          Zpoly_print(stdout, test_poly4); printf("\n");
#endif    
          result = Zpoly_equal(test_poly4, test_poly3);
          Zpoly_clear(test_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000) + 1;
      bits2 = bits+gmp_urandomm_ui(Zpoly_test_randstate,200) + 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
          length2 = length + gmp_urandomm_ui(Zpoly_test_randstate,1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          Zpoly_add(test_poly3, test_poly, test_poly2);
          
          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_add(test_mpn_poly2, test_mpn_poly, test_mpn_poly2);
          Zpoly_init3(test_poly4, length2, bits2);
          _fmpz_poly_convert_out(test_poly4, test_mpn_poly2);
          
#if DEBUG
          Zpoly_print(stdout, test_poly); printf("\n");
          Zpoly_print(stdout, test_poly2); printf("\n");
          Zpoly_print(stdout, test_poly3); printf("\n");
          Zpoly_print(stdout, test_poly4); printf("\n");
#endif    
          result = Zpoly_equal(test_poly4, test_poly3);
          Zpoly_clear(test_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   Zpoly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_sub()
{
   Zpoly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, bits3, length, length2, max_length;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   Zpoly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      bits2 = bits+gmp_urandomm_ui(Zpoly_test_randstate,200);
      bits3 = bits2+gmp_urandomm_ui(Zpoly_test_randstate,200)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits3-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
          length2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld\n",length, length2, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          max_length = (length > length2) ? length : length2;
          fmpz_poly_realloc(test_mpn_poly3, max_length);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          Zpoly_sub(test_poly3, test_poly, test_poly2);

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_sub(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          Zpoly_init3(test_poly4, max_length, bits3+1);
          _fmpz_poly_convert_out(test_poly4, test_mpn_poly3);
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = Zpoly_equal(test_poly4, test_poly3);
          Zpoly_clear(test_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000) + 1;
      bits2 = bits+gmp_urandomm_ui(Zpoly_test_randstate,200) + 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
          length2 = length + gmp_urandomm_ui(Zpoly_test_randstate,1000);        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n",length, length2, bits, bits2);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2-1); 
          Zpoly_sub(test_poly3, test_poly, test_poly2);
          
          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
          
          _fmpz_poly_sub(test_mpn_poly2, test_mpn_poly, test_mpn_poly2);
          Zpoly_init3(test_poly4, length2, bits2);
          _fmpz_poly_convert_out(test_poly4, test_mpn_poly2);
          
#if DEBUG
          for (unsigned j = 0; j < test_poly3->length; j++)
             gmp_printf("%Zd, ",test_poly3->coeffs[j]);
          printf("\n\n");
          for (unsigned j = 0; j < test_poly4->length; j++)
             gmp_printf("%Zd, ",test_poly4->coeffs[j]);
          printf("\n\n");
#endif    
          result = Zpoly_equal(test_poly4, test_poly3);
          Zpoly_clear(test_poly4);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   Zpoly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_scalar_mul_ui()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   mpz_t temp;
   mpz_init(temp);
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
         
          Zpoly_init3(test_poly2, length, bits+FLINT_BITS_PER_LIMB);
          _fmpz_poly_convert_out(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif              
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_ui(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          Zpoly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_mul_si()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length, sign;
   long mult;
   mpz_t temp;
   mpz_init(temp);
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = (long) randint(34682739);
          sign = randint(2);
          if (sign) mult = -mult;
          
          _fmpz_poly_scalar_mul_si(test_mpn_poly2, test_mpn_poly, mult);
         
          Zpoly_init3(test_poly2, length, bits+FLINT_BITS_PER_LIMB);
          _fmpz_poly_convert_out(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_mul_si(temp, test_poly->coeffs[i], mult);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }          
          Zpoly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_div_exact_ui()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   unsigned long mult;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
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
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_mpn_poly2, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly2, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS_PER_LIMB+3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
          _fmpz_poly_scalar_mul_ui(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_ui(test_mpn_poly3, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_scalar_div_exact_si()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, length;
   long mult;
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
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
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
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
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      fmpz_poly_init2(test_mpn_poly3, 1, (bits-1)/FLINT_BITS_PER_LIMB+3);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          fmpz_poly_realloc(test_mpn_poly3, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          mult = randint(34682739);
          if (randint(2)) mult = -mult;
          _fmpz_poly_scalar_mul_si(test_mpn_poly2, test_mpn_poly, mult);
          _fmpz_poly_scalar_div_exact_si(test_mpn_poly3, test_mpn_poly2, mult);
          
          result = _fmpz_poly_equal(test_mpn_poly3, test_mpn_poly);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   Zpoly_clear(test_poly);
   
   return result; 
}

int test_fmpz_poly_scalar_div_ui()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   unsigned long div;
   mpz_t temp;
   mpz_init(temp);
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 0; (count1 < 10) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      //bits = 64*1000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
          //length = 10000;       
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          div = randint(34682739)+1;
          for (unsigned long i = 0; i < 100; i++)
          {
             _fmpz_poly_scalar_div_ui(test_mpn_poly2, test_mpn_poly, div);
          }
          
          Zpoly_init3(test_poly2, length, bits);
          _fmpz_poly_convert_out(test_poly2, test_mpn_poly2); 
          
#if DEBUG
          printf("length = %ld\n",_fmpz_poly_length(test_mpn_poly));
#endif    
          for (unsigned long i = 0; i < test_poly->length; i++)
          {
              mpz_tdiv_q_ui(temp, test_poly->coeffs[i], div);
              result &= (mpz_cmp(temp, test_poly2->coeffs[i]) == 0);
          }         
          Zpoly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_scalar_div_si()
{
   Zpoly_t test_poly, test_poly2;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2;
   int result = 1;
   unsigned long bits, length;
   long div;
   mpz_t temp;
   mpz_init(temp);
   
   Zpoly_init(test_poly); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits-1)/FLINT_BITS_PER_LIMB+2);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;        
#if DEBUG
          printf("length = %ld, bits = %ld\n",length, bits);
#endif
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length);
          randpoly(test_poly, length, bits); 

          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
          div = randint(34682739)+1;
          if (randint(2)) div = -div;
          _fmpz_poly_scalar_div_si(test_mpn_poly2, test_mpn_poly, div);
         
          Zpoly_init3(test_poly2, length, bits);
          _fmpz_poly_convert_out(test_poly2, test_mpn_poly2); 
          
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
          Zpoly_clear(test_poly2);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   mpz_clear(temp);
   
   return result; 
}

int test_fmpz_poly_mul_naieve()
{
   Zpoly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   Zpoly_init(test_poly3); 
   
   for (unsigned long count1 = 1; (count1 < 200) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      bits2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
          length = gmp_urandomm_ui(Zpoly_test_randstate,200)+1; 
          length2 = gmp_urandomm_ui(Zpoly_test_randstate,200)+1;        
#if DEBUG
          printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif
          randpoly(test_poly, length, bits); 
          randpoly(test_poly2, length2, bits2); 
          
          fmpz_poly_realloc(test_mpn_poly, length);
          fmpz_poly_realloc(test_mpn_poly2, length2);
          _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
          
          Zpoly_mul_naive(test_poly3, test_poly, test_poly2);
          
          Zpoly_init3(test_poly4, length+length2-1, bits+bits2+FLINT_BITS_PER_LIMB);
          fmpz_poly_init2(test_mpn_poly3, length+length2-1, (bits+bits2-1)/FLINT_BITS_PER_LIMB+2);
          
          _fmpz_poly_mul_naive(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          
          _fmpz_poly_convert_out(test_poly4, test_mpn_poly3); 
          
          result = _Zpoly_equal(test_poly4, test_poly3);
          Zpoly_clear(test_poly4);
          fmpz_poly_clear(test_mpn_poly3);
      }
      fmpz_poly_clear(test_mpn_poly);
      fmpz_poly_clear(test_mpn_poly2);
   }
   
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   Zpoly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_karatsuba()
{
   Zpoly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   Zpoly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 20) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      bits2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      
      length2 = gmp_urandomm_ui(Zpoly_test_randstate,200)+1; 
      length = gmp_urandomm_ui(Zpoly_test_randstate,200)+1;   
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif
      randpoly(test_poly, length, bits); 
      randpoly(test_poly2, length2, bits2); 
          
      fmpz_poly_realloc(test_mpn_poly, length);
      fmpz_poly_realloc(test_mpn_poly2, length2);
      _fmpz_poly_convert_in(test_mpn_poly, test_poly);
      _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
#if DEBUG
      Zpoly_print(stdout, test_poly);printf("\n");
      Zpoly_print(stdout, test_poly2);printf("\n");
#endif          
      Zpoly_mul_naive(test_poly3, test_poly, test_poly2);
          
      Zpoly_init3(test_poly4, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+2);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+2);
          
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_karatsuba(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
          
          _fmpz_poly_convert_out(test_poly4, test_mpn_poly3); 
#if DEBUG
          Zpoly_print(stdout, test_poly3);printf("\n");
          Zpoly_print(stdout, test_poly4);printf("\n");
#endif          
          result = _Zpoly_equal(test_poly4, test_poly3);
          
      }
      
      Zpoly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   Zpoly_clear(test_poly3);
   
   return result; 
}

int test_fmpz_poly_mul_KS()
{
   Zpoly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   ZmodFpoly_t test_modF_poly;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   Zpoly_init(test_poly3); 
   Zpoly_init(test_poly4); 
      
   for (unsigned long count1 = 0; count1 < 10; count1++)
   {
      
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000) + 1;
      bits2 = gmp_urandomm_ui(Zpoly_test_randstate,1000) + 1;
      //bits = 4;
      //bits2 = 4;
      length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;
      length2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1;
      //length = 32000;
      //length2 = 32000;   
   
      _fmpz_poly_stack_init(test_mpn_poly, length, (bits-1)/FLINT_BITS_PER_LIMB+1);
      _fmpz_poly_stack_init(test_mpn_poly2, length2, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      _fmpz_poly_stack_init(test_mpn_poly3, length + length2 - 1, test_mpn_poly->limbs + test_mpn_poly2->limbs + 1);
#if DEBUG
      printf("%ld, %ld, %ld, %ld\n", length, length2, bits, bits2);
#endif
      Zpoly_realloc(test_poly, length);
      Zpoly_realloc(test_poly2, length2);
      Zpoly_realloc(test_poly3, length + length2 - 1);
      Zpoly_realloc(test_poly4, length + length2 - 1);
          
      for (unsigned long count2 = 0; (count2 < 10) && (result == 1); count2++)
      { 
         randpoly(test_poly, length, bits);
         randpoly(test_poly2, length2, bits2);      
          
#if DEBUG
         for (unsigned j = 0; j < test_poly->length; j++)
            gmp_printf("%Zx, ",test_poly->coeffs[j]);
         printf("\n\n");
         for (unsigned j = 0; j < test_poly2->length; j++)
            gmp_printf("%Zx, ",test_poly2->coeffs[j]);
         printf("\n\n");
#endif
         _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
         _fmpz_poly_convert_in(test_mpn_poly, test_poly);
          
         Zpoly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
         for (unsigned long i = 0; i < 10; i++)
            _fmpz_poly_mul_KS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
         _fmpz_poly_convert_out(test_poly4, test_mpn_poly3);
                  
#if DEBUG
         for (unsigned j = 0; j < test_poly3->length; j++)
            gmp_printf("%Zx, ",test_poly3->coeffs[j]);
         printf("\n\n");
         for (unsigned j = 0; j < test_poly4->length; j++)
            gmp_printf("%Zx, ",test_poly4->coeffs[j]);
         printf("\n\n");
#endif         
         result = Zpoly_equal(test_poly3, test_poly4);
      }   
   
      _fmpz_poly_stack_clear(test_mpn_poly3);
      _fmpz_poly_stack_clear(test_mpn_poly2);
      _fmpz_poly_stack_clear(test_mpn_poly);
   
   }
   
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   Zpoly_clear(test_poly3);
   Zpoly_clear(test_poly4);
   
   return result;
}

int test_fmpz_poly_mul_SS()
{
   Zpoly_t test_poly, test_poly2, test_poly3, test_poly4;
   fmpz_poly_t test_mpn_poly, test_mpn_poly2, test_mpn_poly3;
   int result = 1;
   unsigned long bits, bits2, length, length2;
   
   Zpoly_init(test_poly); 
   Zpoly_init(test_poly2); 
   Zpoly_init(test_poly3); 
   
   for (unsigned long count1 = 0; (count1 < 10) && (result == 1) ; count1++)
   {
      bits = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      bits2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+ 1;
      //bits = bits2 = 1000;
      
      fmpz_poly_init2(test_mpn_poly, 1, (bits-1)/FLINT_BITS_PER_LIMB+1);
      fmpz_poly_init2(test_mpn_poly2, 1, (bits2-1)/FLINT_BITS_PER_LIMB+1);
      
      length2 = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
      length = gmp_urandomm_ui(Zpoly_test_randstate,1000)+1; 
      //length = length2 = 256;
       
#if DEBUG
      printf("length = %ld, length2 = %ld, bits = %ld, bits2 = %ld\n", length, length2, bits, bits2);
#endif

      randpoly(test_poly, length, bits); 
      randpoly(test_poly2, length2, bits2); 
          
      fmpz_poly_realloc(test_mpn_poly, length);
      fmpz_poly_realloc(test_mpn_poly2, length2);
      _fmpz_poly_convert_in(test_mpn_poly, test_poly);
      _fmpz_poly_convert_in(test_mpn_poly2, test_poly2);
#if DEBUG
      Zpoly_print(stdout, test_poly);printf("\n\n");
      Zpoly_print(stdout, test_poly2);printf("\n\n");
#endif          
      Zpoly_mul_naive_KS(test_poly3, test_poly, test_poly2);
          
      Zpoly_init3(test_poly4, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
      fmpz_poly_init2(test_mpn_poly3, length+length2-1, test_mpn_poly->limbs+test_mpn_poly2->limbs+1);
          
      for (unsigned long count2 = 0; (count2 < 100) && (result == 1); count2++)
      { 
          _fmpz_poly_mul_SS(test_mpn_poly3, test_mpn_poly, test_mpn_poly2);
      }
      
      _fmpz_poly_convert_out(test_poly4, test_mpn_poly3);
           
      result = _Zpoly_equal(test_poly4, test_poly3);
      if (!result)
      {
#if DEBUG
          Zpoly_print(stdout, test_poly);printf("\n\n");
          Zpoly_print(stdout, test_poly2);printf("\n\n");
          Zpoly_print(stdout, test_poly3);printf("\n\n");
          Zpoly_print(stdout, test_poly4);printf("\n\n");
#endif               
      }
      
      Zpoly_clear(test_poly4);
      fmpz_poly_clear(test_mpn_poly3);
   }
   
   Zpoly_clear(test_poly);
   Zpoly_clear(test_poly2);
   Zpoly_clear(test_poly3);
   
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
   RUN_TEST(fmpz_poly_negate);
   RUN_TEST(fmpz_poly_shift);
   RUN_TEST(fmpz_poly_add);
   RUN_TEST(fmpz_poly_sub);
   RUN_TEST(fmpz_poly_scalar_mul_ui);
   RUN_TEST(fmpz_poly_scalar_mul_si);
   RUN_TEST(fmpz_poly_scalar_div_exact_ui);
   RUN_TEST(fmpz_poly_scalar_div_exact_si);
   RUN_TEST(fmpz_poly_scalar_div_ui);
   RUN_TEST(fmpz_poly_scalar_div_si);
   RUN_TEST(fmpz_poly_mul_naieve);
   RUN_TEST(fmpz_poly_mul_karatsuba);
   RUN_TEST(fmpz_poly_mul_KS);
   RUN_TEST(fmpz_poly_mul_SS);
   
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   gmp_randinit_default(Zpoly_test_randstate);
   fmpz_poly_test_all();

   return 0;
}

