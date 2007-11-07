/*
   zmod_poly-test.c: Test functions for z_modp_poly-test.c
   
   Copyright (C) 2007, David Howden
*/

#include <stdlib.h>
#include "zmod_poly.h"
#include "long_extras.h"

void tests1();
void tests2();
unsigned long tests3(unsigned long poly_size);
void tests4();

int main (int argc, char const *argv[])
{
   unsigned long p;
   unsigned long tests;

   if (argc > 1)
   {
      p = z_nextprime(atoi(argv[1]));
      
      if (argc > 2)
      {
         tests = atoi(argv[2]);
      }
      else
      {
         tests = 1;
      }
   }
   else
   {
      p = 2;
      tests = 2000;
   }
   
   unsigned int i = 0;
   unsigned long failures = 0;

   while (i < tests)
   {
      p = z_nextprime(p);
      failures += tests3(p);
      i++;
   }
   
   if (failures)
   {
      printf("%d Failures Occurred in %d tests:-(\n", failures, tests);
   }
   else
   {
      printf("%d tests completed successfully :-)\n", tests);
   }
   
   return 0;
}

void tests4()
{
   zmod_poly_t p1, p2;
}

unsigned long tests3(unsigned long poly_size)
{
   zmod_poly_t poly1, poly2;
   zmod_poly_init2(poly1, poly_size, poly_size);
   
   srand(time(NULL));
   
   for (unsigned long i = 0; i < poly_size; i++)
   {
      zmod_poly_set_coeff(poly1, i, rand());
   }
   
   zmod_poly_set_coeff(poly1, poly_size - 1, 1);
   
   unsigned long limbs;
   
   unsigned long failures = 0;

   mp_limb_t * mpn;
   
   unsigned long bits = FLINT_BIT_COUNT(poly_size);
   
   //for (; bits <= (2*FLINT_BITS); bits++)
   //{
      zmod_poly_init2(poly2, poly_size, poly_size);
      if (bits != FLINT_BITS)
      {
         limbs = bits * poly_size / FLINT_BITS + 1;
      }
      else
      {
         limbs = poly_size;
      }
      
      mpn = (mp_limb_t*) malloc(limbs*sizeof(mp_limb_t));

      zmod_poly_bit_pack_mpn(mpn, poly1, bits);
      
      // if(bits == 48) {
      // for( unsigned long i = 0; i < limbs; i++)
      // {
      //    print_binary(mpn[i], FLINT_BITS);
      //    printf("\n");
      // }
      // }

      zmod_poly_bit_unpack_mpn(poly2, mpn, poly_size, bits);

      poly2->length = poly_size;
      
      if (zmod_poly_equal(poly1, poly2))
      {
      }
      else
      {
         failures++;
         printf("FAILURE:  poly_size=%d\t bits=%d\n", poly_size, bits);
         printf("%s\n%s\n", zmod_poly_to_string(poly1), zmod_poly_to_string(poly2));
         // zmod_poly_print(poly2);
         //printf(" Result on poly_size %d with %d bits\n", poly_size, bits);
         // unsigned long space_bit = bits;
         // for( unsigned long i = 0; i < limbs; i++)
         // {
         //    print_binary(mpn[i], FLINT_BITS);
         //    printf("\n");
         // }
      }

      free(mpn);
      zmod_poly_clear(poly2);
   //}
   zmod_poly_clear(poly1);
   
   
   //printf("Tests Completed, %d OK, %d Failed.\n", count, FLINT_BITS - count);
   return failures;
}


void tests2()
{
   zmod_poly_t poly1, poly2;
   zmod_poly_init2(poly1, 7, 2);
   zmod_poly_init2(poly2, 7, 2);
   for (unsigned long i = 0; i < 2; i++)
   {
      zmod_poly_set_coeff(poly1, i, i+2);
      zmod_poly_set_coeff(poly2, i, i+2);
   }
   
   zmod_poly_t ans, ans2;
   zmod_poly_init(ans, 7);
   zmod_poly_init(ans2, 7);
   
   printf("poly1 = %s\n", zmod_poly_to_string(poly1));
   printf("poly2 = %s\n", zmod_poly_to_string(poly2));
   
   zmod_poly_mul_KS(ans, poly1, poly2, 44);
   zmod_poly_mul(ans2, poly1, poly2);
   printf("ans  = %s\n", zmod_poly_to_string(ans));
   printf("ans2 = %s\n", zmod_poly_to_string(ans2));
}

void tests1()
{
   zmod_poly_t poly, poly2;
   zmod_poly_init2(poly, 7, 10);
   zmod_poly_init2(poly2, 7, 10);
   for(unsigned long i = 0; i < 10; i++)
   {
      zmod_poly_set_coeff(poly, i, i);
      zmod_poly_set_coeff(poly2, i, i);
   }
   
   printf("poly = %s\n", zmod_poly_to_string(poly));
   
   zmod_poly_lshift(poly, poly, 5);
   
   printf("poly = %s\n", zmod_poly_to_string(poly));
   
   zmod_poly_rshift(poly, poly, 5);
   
   printf("poly  = %s\n", zmod_poly_to_string(poly));
   printf("poly2 = %s\n", zmod_poly_to_string(poly2));
   
   if(zmod_poly_equal(poly, poly2))
   {
      printf("The polys are equal...\n");
   }
   else
   {
      printf("The polys are not equal...\n");
   }
   
   zmod_poly_t sum, diff;
   zmod_poly_init2(sum, 7, 10);
   zmod_poly_init2(diff, 7, 10);
   
   zmod_poly_add(sum, poly, poly2);
   
   printf("sum  = %s\n", zmod_poly_to_string(sum));
   
   zmod_poly_t temp;
   zmod_poly_init(temp, 7);
   zmod_poly_sub(temp, sum, poly);
   
   printf("temp = %s\n", zmod_poly_to_string(temp));
   

   
   zmod_poly_t temp2;
   zmod_poly_init2(temp2, 7, 2);
   if(!zmod_poly_from_string(temp2, "10  3  1 1 1 1 0 2 2 2 2 1 "))
   {
      printf("Error in creating temp2");
   }
   zmod_poly_print(temp2);
   
   printf("\n\nPlease enter a polynomial: ");
   zmod_poly_read(temp2);
   printf("\n\n");
   zmod_poly_print(temp2);
   
   zmod_poly_t product;
   zmod_poly_init(product, 7);
   
   printf("poly = %s\n\n", zmod_poly_to_string(poly));
   
   zmod_poly_mul(product, poly, poly2);
   
   printf("(%s) * (%s) = %s\n\n", zmod_poly_to_string(poly), zmod_poly_to_string(poly2), zmod_poly_to_string(product));
   
   zmod_poly_clear(product);
   zmod_poly_clear(temp2);
   zmod_poly_clear(temp);
   zmod_poly_clear(sum);
   zmod_poly_clear(diff);
   zmod_poly_clear(poly);
   zmod_poly_clear(poly2);
}
