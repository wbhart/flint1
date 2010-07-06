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
/******************************************************************************

mpq_mat.c: Matrices over Q, implemented as an array of mpq_t's
           Not intended to be efficient

Copyright (C) 2010 Andy Novocin, Max Flander

******************************************************************************/

#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "mpq_mat.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpq_mat_init(mpq_mat_t mat, ulong r, ulong c)
{
   if ((r) && (c)) mat->entries = (mpq_t*) flint_heap_alloc_bytes(r*c*sizeof(mpq_t));
	else mat->entries = NULL;

   long i;
   for (i = 0; i < r*c; i++)
	{
	   mpq_init(mat->entries[i]);
	}
	mat->r = r;
	mat->c = c;
}

void mpq_mat_clear(mpq_mat_t mat)
{
   long i;
   for (i = 0; i < mat->r*mat->c; i++)
      mpq_clear(mat->entries[i]);

   if (mat->entries) flint_heap_free(mat->entries);
	mat->entries = NULL;
	mat->r = 0;
	mat->c = 0;
}

void mpz_mat_to_mpq_mat(mpq_mat_t res, mpz_mat_t mat){

   if ((res->r != mat->r) || (res->c != mat->c)){
      printf("FLINT exception: dimensions don't match\n");
      abort();
   }

   long i;
   for (i = 0; i < mat->r*mat->c; i++)
      mpq_set_z(res->entries[i], mat->entries[i]);

   return;
}

void mpq_mat_print(mpq_mat_t mat){

   printf("%ld %ld  ", mat->r, mat->c);

   long i;
   for (i = 0; i < mat->r*mat->c; i++){
      mpq_out_str(NULL, 10, mat->entries[i]); printf(" ");
   }

   return;
}

void mpq_mat_print_pretty(mpq_mat_t mat){

   printf("[\n[");
   long i;
   for (i = 0; i < mat->r*mat->c; i++){
      mpq_out_str(NULL, 10, mat->entries[i]);
      if (((i+1) % mat->c) != 0) printf(" ");
      else if ((i != 0) && ((i+1)!= mat->r * mat->c) ) printf("]\n[");
   }
   printf("]\n]\n");

   return;
}

void mpq_mat_row_scalar_product(mpq_mat_t mat1, ulong r1, mpq_mat_t mat2, ulong r2, mpq_t scalar){

   if (mat2->c != mat1->c){
      printf("FLINT exception: dimensions don't match\n");
      abort();
   }

   long i;
   for (i = 0; i< mat1->c; i++)
      mpq_mul(mat2->entries[r2*mat2->c + i],mat1->entries[r1*mat1->c + i], scalar);

   return;

}

void mpq_mat_row_inner_product(mpq_t res, mpq_mat_t mat1, ulong r1, mpq_mat_t mat2, ulong r2){

   if (mat2->c != mat1->c){
      printf("FLINT exception: dimensions don't match\n");
      abort();
   }

   mpq_set_ui(res, 0L, 1L);
   mpq_t temp;
   mpq_init(temp);

   long i;
   for (i = 0; i< mat1->c; i++){
      mpq_mul(temp, mat2->entries[r2*mat2->c + i],  mat1->entries[r1*mat1->c + i]);
      mpq_add(res, res, temp);
   }

   mpq_clear(temp);

   return;
}

void mpq_mat_row_add(mpq_mat_t mat1, ulong r1, mpq_mat_t mat2, ulong r2){

   if (mat2->c != mat1->c){
      printf("FLINT exception: dimensions don't match\n");
      abort();
   }

   long i;
   for (i = 0; i< mat1->c; i++)
      mpq_add(mat1->entries[r1*mat1->c + i], mat2->entries[r2*mat2->c + i], mat1->entries[r1*mat1->c + i]);

   return;

}

void mpq_mat_GS(mpq_mat_t mu, mpq_mat_t GS, mpq_mat_t mat){

   if ( ( GS->c != mat->c ) || ( GS->r != mat->r ) ){
      printf("FLINT exception: dimensions don't match\n");
      abort();
   }

   if ( ( mu->r != mu->c) || (mu->r != mat->r) ){
      printf("FLINT exception: mu dimensions don't match\n");
      abort();
   }

  //I'm going to use mu[0] to store <GS[i], GS[i]> until the end

   mpq_t temp;
   mpq_init(temp);

   long i, j, k;
   //setting GS[0] := mat[0]
   for (i = 0; i < mat->c; i++)
      mpq_set(GS->entries[i], mat->entries[i]);

   mpq_mat_row_inner_product(mu->entries[0], GS, 0, GS, 0);

   //mu[i,i] := 1
   for (i = 1; i < mu->r; i++)
      mpq_set_ui(mu->entries[i*mu->c + i], 1L, 1L);

   for (i = 1; i < mat->r; i++){
      //in this loop we want to find GS[i] := mat[i] - sum (mu[i,j]*mat[j])
      //start with GS[i] = mat[i] then for each j < i compute mu and subtract
      for (k = 0; k < mat->c; k++)
         mpq_set(GS->entries[i*mat->c + k], mat->entries[i*mat->c + k]);

      for (j = 0; j < i; j++){
         //temp will be the numerator of mu[i,j] which is <mat[i], GS[j]>
         mpq_mat_row_inner_product(temp, mat, i, GS, j);

         if (mpq_sgn(mu->entries[j]) != 0)
            mpq_div(mu->entries[i*mu->c + j], temp, mu->entries[j]);
         else
            mpq_set_ui(mu->entries[i*mu->c + j], 0L, 1L);

         //now need GS[i] := GS[i] - mu[i,j]*GS[j]
         for (k=0; k < mat->c; k++){
            //temp = mu[i,j] * GS[j,k]
            mpq_mul(temp, mu->entries[i*mu->c + j], GS->entries[j*GS->c + k]);
            mpq_neg(temp, temp);
            mpq_add(GS->entries[i*GS->c + k], GS->entries[i*GS->c + k], temp);
         }
      }
      if (i+1 < mu->c) 
         mpq_mat_row_inner_product(mu->entries[i], GS, i, GS, i);
   }

   mpq_set_ui(mu->entries[0], 1L, 1L);
   for (k = 1; k < mu->c; k++)
      mpq_set_ui(mu->entries[k], 0L, 1L);

   mpq_clear(temp);
   return;
}

int mpq_mat_is_reduced(mpq_mat_t mu, mpq_mat_t GS, double delta, double eta){

   //want to return 1 if this data could come from a reduced matrix 0 otherwise

   mpq_mat_t gs_len;
   mpq_mat_init( gs_len, 1, GS->c);

   mpq_t temp, temp1, temp2;
   mpq_init(temp);
   mpq_init(temp1);
   mpq_init(temp2);

   long i, j;
   int result = 1;
   for (i = 0; (i < GS->r) && (result == 1); i++){
      mpq_mat_row_inner_product(gs_len->entries[i], GS, i, GS, i);
      if (i > 0){
         mpq_div(temp, gs_len->entries[i], gs_len->entries[i-1]);
         mpq_mul(temp1, mu->entries[i*mu->r + i-1], mu->entries[i*mu->r + i-1]); 
         mpq_add(temp, temp, temp1);
         mpq_set_d(temp1, delta - eta*eta);
         if (mpq_cmp(temp, temp1) < 0){
            result = 0;
         }
         else{
            mpq_set_d(temp2, eta);
            for( j = 0 ; (j < i) && (result == 1); j++)
               if ( mpq_cmp( mu->entries[i*mu->r + j], temp2) >= 0){
                  result = 0;
               }
         }
// if temp < (3/4 or 1/(delta - eta^2))==temp1 then not reduced...
      }
   }

   mpq_clear(temp);
   mpq_clear(temp1);
   mpq_clear(temp2);
   mpq_mat_clear(gs_len);
   return result;
}


// *************** end of file
