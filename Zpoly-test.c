/****************************************************************************

Zpoly-test.c: Test code for Zpoly.c and Zpoly.h

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "Zpoly.h"
#include <stdio.h>


int main()
{
   char *s;

   Zpoly_mpz_t poly1, poly2, poly3;
   Zpoly_mpz_init(poly1);
   Zpoly_mpz_init(poly2);
   Zpoly_mpz_init(poly3);

   Zpoly_mpz_set_from_string(poly1, " -1 2 0  34   -123123123123123123 34 ");
   Zpoly_mpz_set_from_string(poly2, "0 49 28");
   Zpoly_mpz_add(poly3, poly1, poly2);

   s = malloc(Zpoly_mpz_get_string_size(poly1));
   Zpoly_mpz_get_as_string(s, poly1);
   printf(" first poly is: %s\n", s);
   free(s);
   
   s = malloc(Zpoly_mpz_get_string_size(poly2));
   Zpoly_mpz_get_as_string(s, poly2);
   printf("second poly is: %s\n", s);
   free(s);
   
   s = malloc(Zpoly_mpz_get_string_size(poly3));
   Zpoly_mpz_get_as_string(s, poly3);
   printf("  their sum is: %s\n", s);
   free(s);
   
   Zpoly_mpz_clear(poly1);
   Zpoly_mpz_clear(poly2);
   Zpoly_mpz_clear(poly3);

   return 0;
}


// *************** end of file
