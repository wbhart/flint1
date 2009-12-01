/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

F_mpz_mpoly-test.c: Test code for F_mpz_mpoly.c and F_mpz_mpoly.h

Copyright (C) 2009, William Hart

*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <gmp.h>
#include <time.h>
#include "flint.h"
#include "F_mpz_mpoly.h"
#include "memory-manager.h"
#include "test-support.h"
#include "long_extras.h"

#define VARY_BITS 1 // random coefficients have random number of bits up to the limit given
#define SIGNS 1 // random coefficients will be randomly signed
#define SPARSE 1 // polynomials are spares (triggers more corner cases)
#define ITER 1 // if you want all tests to run longer, increase this

#define TESTFILE 0 // Set this to test polynomial reading and writing to a file in the current dir

#define DEBUG 0 // allows easy switching of debugging code on and off when debugging (if inserted)
#define DEBUG2 1 

void rand_mpoly(F_mpz_mpoly_t poly, ulong length, ulong vars, ulong expmax, ulong bits)
{
   ulong sumlast = 0;
	ulong i = 0, k = 0;

	if (z_randint(2))
	{
		F_mpz_mpoly_set_coeff_ui(poly, 0, z_randint(1L << bits));
		i++;
		k++;
	}

	ulong summax = 0;
	
	for ( ; (i < length) && (summax <= expmax); i++, k++)
	{
		summax = (expmax*k)/length;
		//printf("%ld\n", summax);
		if ((summax > sumlast) && (summax <= expmax))
		{
			F_mpz_mpoly_set_coeff_ui(poly, i, z_randint(1L << bits));
		
			ulong sum = 0;
		   ulong exp;
			long j;
			for (j = 0; j < vars - 1; j++)
		   {
			   exp = z_randint(2*summax/vars);
			   if (sum + exp > summax) exp = summax - sum;
			   if (exp) F_mpz_mpoly_set_var_exp(poly, i, j, exp);
				sum += exp;
		   }  
         exp = summax - sum;
         if (exp) F_mpz_mpoly_set_var_exp(poly, i, j, exp);
		} else i--;
		sumlast = summax;
	}
}

int test_F_mpz_mpoly_add_inplace()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t poly2;

	char * syms[5];
	syms[0] = "s";
	syms[1] = "t";
	syms[2] = "x";
	syms[3] = "y";
   syms[4] = "z";
	
	ulong i;
	for (i = 0; i < 100; i++)
	{
		F_mpz_mpoly_init2(poly, 10, 5, GRLEX);
      F_mpz_mpoly_init2(poly2, 10, 5, GRLEX);

	   ulong length = z_randint(20)+1;
      ulong length2 = z_randint(20)+1;

	   rand_mpoly(poly, length, 5, 30, 8);
	   rand_mpoly(poly2, length2, 5, 30, 8);
	
	   printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	   printf("B = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n");
   
	   F_mpz_mpoly_add_inplace(poly, 0, poly2);

	   printf("A + B = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n\n");
	
	   F_mpz_mpoly_clear(poly);
	   F_mpz_mpoly_clear(poly2);
	}
   
   return result;
}

int test__F_mpz_mpoly_mul_mxn()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t poly2;
   F_mpz_mpoly_t res;

	char * syms[5];
	syms[0] = "s";
	syms[1] = "t";
	syms[2] = "x";
	syms[3] = "y";
   syms[4] = "z";
	
	// 2 x 2
	F_mpz_mpoly_init2(poly, 2, 5, GRLEX);
   F_mpz_mpoly_init2(poly2, 2, 5, GRLEX);
   F_mpz_mpoly_init2(res, 4, 5, GRLEX);

	rand_mpoly(poly, 2, 5, 30, 8);
	rand_mpoly(poly2, 2, 5, 30, 8);
	
	printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	printf("B = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n");
   
	_F_mpz_mpoly_mul_small_mxn(res, poly, 0, poly2, 0);

	printf("A * B = "); F_mpz_mpoly_print_pretty(res, syms); printf("\n\n");
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(res);
   
   // 3 x 3
	
	F_mpz_mpoly_init2(poly, 3, 5, GRLEX);
   F_mpz_mpoly_init2(poly2, 3, 5, GRLEX);
   F_mpz_mpoly_init2(res, 9, 5, GRLEX);

	rand_mpoly(poly, 3, 5, 30, 8);
	rand_mpoly(poly2, 3, 5, 30, 8);
	
	printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	printf("B = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n");
   
	_F_mpz_mpoly_mul_small_mxn(res, poly, 0, poly2, 0);

	printf("A * B = "); F_mpz_mpoly_print_pretty(res, syms); printf("\n\n");
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(res);
   
   // 4 x 4
	F_mpz_mpoly_init2(poly, 4, 5, GRLEX);
   F_mpz_mpoly_init2(poly2, 4, 5, GRLEX);
   F_mpz_mpoly_init2(res, 16, 5, GRLEX);

	rand_mpoly(poly, 4, 5, 30, 8);
	rand_mpoly(poly2, 4, 5, 30, 8);
	
	printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	printf("B = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n");
   
	ulong i;
	for (i = 0; i < 1000000; i++) _F_mpz_mpoly_mul_small_mxn(res, poly, 0, poly2, 0);

	printf("A * B = "); F_mpz_mpoly_print_pretty(res, syms); printf("\n\n");
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(res);
   
   return result;
}

int test__F_mpz_mpoly_mul_Mxn()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t poly2;
   F_mpz_mpoly_t res;

	char * syms[5];
	syms[0] = "s";
	syms[1] = "t";
	syms[2] = "x";
	syms[3] = "y";
   syms[4] = "z";
	
	// 32 x 4
	F_mpz_mpoly_init2(poly, 32, 5, GRLEX);
   F_mpz_mpoly_init2(poly2, 4, 5, GRLEX);
   F_mpz_mpoly_init2(res, 128, 5, GRLEX);

	rand_mpoly(poly, 32, 5, 20, 8);
	rand_mpoly(poly2, 4, 5, 20, 8);
	
	printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	printf("B = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n");
   
	_F_mpz_mpoly_mul_small_Mxn(res, poly, poly2, 0, poly2->length);

	printf("A * B = "); F_mpz_mpoly_print_pretty(res, syms); printf("\n\n");
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(res);
      
   return result;
}

int test_F_mpz_mpoly_mul_small_recursive()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t poly2;
   F_mpz_mpoly_t res;

	char * syms[5];
	syms[0] = "s";
	syms[1] = "t";
	syms[2] = "x";
	syms[3] = "y";
   syms[4] = "z";
	
	// 32 x 4
	F_mpz_mpoly_init2(poly, 32, 5, GRLEX);
   F_mpz_mpoly_init2(poly2, 32, 5, GRLEX);
   F_mpz_mpoly_init2(res, 0, 5, GRLEX);

	rand_mpoly(poly, 32, 5, 20, 8);
	rand_mpoly(poly2, 32, 5, 20, 8);
	
	printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	printf("B = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n");
   
	F_mpz_mpoly_mul_small_recursive(res, poly, poly2, 0, poly2->length);

	printf("A * B = "); F_mpz_mpoly_print_pretty(res, syms); printf("\n\n");
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(res);
      
   return result;
}

int test_F_mpz_mpoly_mul_fateman_heap()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t poly2;
   F_mpz_mpoly_t poly4;
   F_mpz_mpoly_t poly8;
   F_mpz_mpoly_t poly10;
   F_mpz_mpoly_t poly20;
   F_mpz_mpoly_t poly40;

	char * syms[4];
	syms[0] = "t";
	syms[1] = "x";
	syms[2] = "y";
   syms[3] = "z";
	
	F_mpz_mpoly_init2(poly, 5, 4, GRLEX);
   F_mpz_mpoly_init2(poly2, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly4, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly8, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly10, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly20, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly40, 0, 4, GRLEX);

	F_mpz_mpoly_set_coeff_ui(poly, 0, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 1, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 2, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 3, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 4, 1);
	F_mpz_mpoly_set_var_exp(poly, 1, 3, 1);
	F_mpz_mpoly_set_var_exp(poly, 2, 2, 1);
	F_mpz_mpoly_set_var_exp(poly, 3, 1, 1);
	F_mpz_mpoly_set_var_exp(poly, 4, 0, 1);

   printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	
	F_mpz_mpoly_mul_small_heap(poly2, poly, poly);

	printf("A^2 = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_heap(poly4, poly2, poly2);

	printf("A^4 = "); F_mpz_mpoly_print_pretty(poly4, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_heap(poly8, poly4, poly4);

	//printf("A^8 = "); F_mpz_mpoly_print_pretty(poly8, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_heap(poly10, poly8, poly2);

	printf("A^10 = "); F_mpz_mpoly_print_pretty(poly10, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_heap(poly20, poly10, poly10);

	printf("A^20 has length %ld\n", poly20->length);
	
	F_mpz_mpoly_mul_small_heap(poly40, poly20, poly20);

	printf("A^40 has length %ld\n", poly40->length);
	
	//printf("A^20 = "); F_mpz_mpoly_print_pretty(poly20, syms); printf("\n\n");
	
	//printf("A^40 = "); F_mpz_mpoly_print_pretty(poly40, syms); printf("\n\n");

	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(poly4);
   F_mpz_mpoly_clear(poly8);
   F_mpz_mpoly_clear(poly10);
   F_mpz_mpoly_clear(poly20);
   F_mpz_mpoly_clear(poly40);
      
   return result;
}

int test_F_mpz_mpoly_mul_5sparse_heap()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t polyb;
	F_mpz_mpoly_t poly2;
	F_mpz_mpoly_t poly3;
   F_mpz_mpoly_t poly6;
   F_mpz_mpoly_t poly12;
   F_mpz_mpoly_t poly24;
   F_mpz_mpoly_t poly12b;
   
	char * syms[5];
	syms[0] = "x";
	syms[1] = "y";
	syms[2] = "z";
   syms[3] = "t";
	syms[4] = "u";
	
	F_mpz_mpoly_init2(poly, 6, 5, GRLEX);
   F_mpz_mpoly_init2(polyb, 6, 5, GRLEX);
   F_mpz_mpoly_init2(poly2, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly3, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly6, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly12, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly24, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly12b, 0, 5, GRLEX);
   
	F_mpz_mpoly_set_coeff_ui(poly, 0, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 1, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 2, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 3, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 4, 1);
	F_mpz_mpoly_set_coeff_ui(poly, 5, 1);
	F_mpz_mpoly_set_var_exp(poly, 1, 0, 1);
	F_mpz_mpoly_set_var_exp(poly, 2, 1, 2);
	F_mpz_mpoly_set_var_exp(poly, 3, 2, 3);
	F_mpz_mpoly_set_var_exp(poly, 4, 3, 5);
   F_mpz_mpoly_set_var_exp(poly, 5, 4, 7);

   F_mpz_mpoly_set_coeff_ui(polyb, 0, 1);
   F_mpz_mpoly_set_coeff_ui(polyb, 1, 1);
   F_mpz_mpoly_set_coeff_ui(polyb, 2, 1);
   F_mpz_mpoly_set_coeff_ui(polyb, 3, 1);
   F_mpz_mpoly_set_coeff_ui(polyb, 4, 1);
	F_mpz_mpoly_set_coeff_ui(polyb, 5, 1);
	F_mpz_mpoly_set_var_exp(polyb, 1, 4, 1);
	F_mpz_mpoly_set_var_exp(polyb, 2, 3, 2);
	F_mpz_mpoly_set_var_exp(polyb, 3, 2, 3);
	F_mpz_mpoly_set_var_exp(polyb, 4, 1, 5);
   F_mpz_mpoly_set_var_exp(polyb, 5, 0, 7);

   printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	
	F_mpz_mpoly_mul_small_heap(poly2, poly, poly);
   F_mpz_mpoly_mul_small_heap(poly3, poly2, poly);

	printf("A^3 = "); F_mpz_mpoly_print_pretty(poly3, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_heap(poly6, poly3, poly3);

	printf("A^6 = "); F_mpz_mpoly_print_pretty(poly6, syms); printf("\n\n");
	
	printf("A^6 has length %ld\n", poly6->length);
	
	F_mpz_mpoly_mul_small_heap(poly12, poly6, poly6);

	printf("A^12 has length %ld\n", poly12->length);
	
	printf("B = "); F_mpz_mpoly_print_pretty(polyb, syms); printf("\n");
	
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(poly3);
	F_mpz_mpoly_clear(poly6);
	F_mpz_mpoly_init2(poly2, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly3, 0, 5, GRLEX);
   F_mpz_mpoly_init2(poly6, 0, 5, GRLEX);
   
	F_mpz_mpoly_mul_small_heap(poly2, polyb, polyb);
   F_mpz_mpoly_mul_small_heap(poly3, poly2, polyb);

	printf("B^3 = "); F_mpz_mpoly_print_pretty(poly3, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_heap(poly6, poly3, poly3);

	//printf("B^6 = "); F_mpz_mpoly_print_pretty(poly6, syms); printf("\n\n");
	
	printf("B^6 has length %ld\n", poly6->length);
	
	F_mpz_mpoly_mul_small_heap(poly12b, poly6, poly6);

	printf("B^12 has length %ld\n", poly12b->length);

	F_mpz_mpoly_mul_small_heap(poly24, poly12, poly12b);

	printf("A^12 * B^12 has length %ld\n", poly24->length);
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(polyb);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(poly3);
   F_mpz_mpoly_clear(poly6);
   F_mpz_mpoly_clear(poly12);
   F_mpz_mpoly_clear(poly24);
   F_mpz_mpoly_clear(poly12b);
      
   return result;
}

int test_F_mpz_mpoly_mul_fateman()
{
   int result = 1;

	F_mpz_mpoly_t poly;
	F_mpz_mpoly_t poly2;
   F_mpz_mpoly_t poly4;
   F_mpz_mpoly_t poly8;
   F_mpz_mpoly_t poly10;
   F_mpz_mpoly_t poly20;
   F_mpz_mpoly_t poly40;

	char * syms[4];
	syms[0] = "t";
	syms[1] = "x";
	syms[2] = "y";
   syms[3] = "z";
	
	F_mpz_mpoly_init2(poly, 5, 4, GRLEX);
   F_mpz_mpoly_init2(poly2, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly4, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly8, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly10, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly20, 0, 4, GRLEX);
   F_mpz_mpoly_init2(poly40, 0, 4, GRLEX);

	F_mpz_mpoly_set_coeff_ui(poly, 0, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 1, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 2, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 3, 1);
   F_mpz_mpoly_set_coeff_ui(poly, 4, 1);
	F_mpz_mpoly_set_var_exp(poly, 1, 3, 1);
	F_mpz_mpoly_set_var_exp(poly, 2, 2, 1);
	F_mpz_mpoly_set_var_exp(poly, 3, 1, 1);
	F_mpz_mpoly_set_var_exp(poly, 4, 0, 1);

   printf("A = "); F_mpz_mpoly_print_pretty(poly, syms); printf("\n");
	
	F_mpz_mpoly_mul_small_recursive(poly2, poly, poly, 0, poly->length);

	printf("A^2 = "); F_mpz_mpoly_print_pretty(poly2, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_recursive(poly4, poly2, poly2, 0, poly2->length);

	//printf("A^4 = "); F_mpz_mpoly_print_pretty(poly4, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_recursive(poly8, poly4, poly4, 0, poly4->length);

	//printf("A^8 = "); F_mpz_mpoly_print_pretty(poly8, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_recursive(poly10, poly8, poly2, 0, poly2->length);

	//printf("A^10 = "); F_mpz_mpoly_print_pretty(poly10, syms); printf("\n\n");
	
	F_mpz_mpoly_mul_small_recursive(poly20, poly10, poly10, 0, poly10->length);

	printf("A^20 has length %ld\n", poly20->length);
	
	//F_mpz_mpoly_mul_small_recursive(poly40, poly20, poly20, 0, poly20->length);

	//printf("A^40 has length %ld\n", poly40->length);
	
	//printf("A^20 = "); F_mpz_mpoly_print_pretty(poly20, syms); printf("\n\n");
	
	//printf("A^40 = "); F_mpz_mpoly_print_pretty(poly40, syms); printf("\n\n");
	
	F_mpz_mpoly_clear(poly);
	F_mpz_mpoly_clear(poly2);
	F_mpz_mpoly_clear(poly4);
   F_mpz_mpoly_clear(poly8);
   F_mpz_mpoly_clear(poly10);
   F_mpz_mpoly_clear(poly20);
   F_mpz_mpoly_clear(poly40);
      
   return result;
}

void F_mpz_mpoly_test_all()
{
   int success, all_success = 1;
   printf("FLINT_BITS = %ld\n", FLINT_BITS);

#if TESTFILE
#endif
   /*RUN_TEST(F_mpz_mpoly_add_inplace); 
	RUN_TEST(_F_mpz_mpoly_mul_mxn); 
	RUN_TEST(_F_mpz_mpoly_mul_Mxn); 
	RUN_TEST(F_mpz_mpoly_mul_small_recursive); */
	//RUN_TEST(F_mpz_mpoly_mul_fateman); 
	//RUN_TEST(F_mpz_mpoly_mul_fateman_heap); 
	RUN_TEST(F_mpz_mpoly_mul_5sparse_heap); 
	
   printf(all_success ? "\nAll tests passed\n" :
                        "\nAt least one test FAILED!\n");
}

int main()
{
   test_support_init();
   F_mpz_mpoly_test_all();
   test_support_cleanup();
   
   flint_stack_cleanup();
	_F_mpz_cleanup();
   
   return 0;
}
