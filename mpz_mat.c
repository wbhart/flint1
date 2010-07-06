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

mpz_mat.c: Matrices over Z, implemented as an array of mpz_t's
           Not intended to be efficient

Copyright (C) 2009, Andy Novocin
Copyright (C) 2008, William Hart

******************************************************************************/

#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "mpz_mat.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/


void mpz_mat_init(mpz_mat_t mat, ulong r, ulong c)
{
   if ((r) && (c)) mat->entries = (mpz_t*) flint_heap_alloc_bytes(r*c*sizeof(mpz_t));
	else mat->entries = NULL;

   long i;
   for (i = 0; i < r*c; i++)
	{
	   mpz_init(mat->entries[i]);
	}
	mat->r = r;
	mat->c = c;
}

void mpz_mat_clear(mpz_mat_t mat)
{
   long i;
   for (i = 0; i < mat->r*mat->c; i++)
      mpz_clear(mat->entries[i]);

   if (mat->entries) flint_heap_free(mat->entries);
	mat->entries = NULL;
	mat->r = 0;
	mat->c = 0;
}

/****************************************************************************

   Add/sub

****************************************************************************/

void mpz_mat_add(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2)
{
	long i;
	for (i = 0; i < mat1->r*mat1->c; i++)
	   mpz_add(res->entries[i], mat1->entries[i], mat2->entries[i]);
}

void mpz_mat_sub(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2)
{
	long i;
	for (i = 0; i < mat1->r*mat1->c; i++)
	   mpz_sub(res->entries[i], mat1->entries[i], mat2->entries[i]);
}

/****************************************************************************

   I/O 

****************************************************************************/

int mpz_mat_from_string(mpz_mat_t mat, const char *s)
{

   const char* whitespace = " \t\n\r";

   //read mat->rows
   unsigned long r;
   if (!sscanf(s, "%ld", &r))
      return 0;

   // jump to next whitespace
   s += strcspn(s, whitespace);

   // skip whitespace
   s += strspn(s, whitespace);

   //read mat->columns
   unsigned long c;
   if (!sscanf(s, "%ld", &c))
      return 0;

   // jump to next whitespace
   s += strcspn(s, whitespace);

   // skip 1 whitespace
   s += strspn(s, whitespace);


   mpz_mat_clear(mat);
   mpz_mat_init(mat,r,c);

   unsigned long i;
   for (i = 0; i < r*c; i++)
   {

      // skip whitespace
      s += strspn(s, whitespace);

      if (!gmp_sscanf(s, "%Zd", mat->entries[i]))
         return 0;

      // jump to next whitespace
      s += strcspn(s, whitespace);

   }

   return 1;
}

char* mpz_mat_to_string(mpz_mat_t mat)
{
   // estimate the size of the string
   // 41 = enough room for null terminator and space and row and column info
   unsigned long size = 41;
   unsigned long i;
   for (i = 0; i < mat->r * mat->c; i++)
      // +2 is for the sign and a space
      size += mpz_sizeinbase(mat->entries[i], 10) + 2;

   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "%ld %ld  ", mat->r, mat->c);
   for (i = 0; i < mat->r * mat->c; i++)
   {
      mpz_get_str(ptr, 10, mat->entries[i]);
      ptr += strlen(ptr);
      *ptr = ' ';
      ptr++;
   }
   
   ptr--;
   *ptr = 0;
   
   return buf;
}

int mpz_mat_from_string_pretty(mpz_mat_t mat, char *s)
{

   char* pnt;

   unsigned long r = 0;
   unsigned long c = 0;

   pnt = s;
//calculates the number of rows by counting the ']'s
   while (pnt != NULL)
   {
      pnt++;
      pnt = strchr(pnt, ']');
      r++;
   }

   r = r - 2;
//reset the pointer and count the number of columns by the number of numbers then later divides by r (not optimal)
   pnt = s;

   while ( pnt != NULL)
   {
      pnt += strspn(pnt,"-0123456789");
      pnt = strpbrk(pnt, "-0123456789");
      if ( pnt != NULL)
         c++;
   }

   pnt = s + strcspn(s, "[")+1;


   if (r == 0){
      mpz_mat_clear(mat);
      mpz_mat_init(mat,0,0);
      return 1;
   }

   if (c == 0){
      mpz_mat_clear(mat);
      mpz_mat_init(mat,r,c);
      return 1;
   }

   c = c/r;

   mpz_mat_clear(mat);
   mpz_mat_init(mat,r,c);


   ulong i;
   for (i = 0; i < r*c; i++){
//searches for the next digit of - then calls gmp's mpz scanner
         pnt = strpbrk(pnt,"-0123456789");
         if (!gmp_sscanf(pnt, "%Zd", mat->entries[i]))
            return 0;
//skips the big number
         pnt += strspn(pnt,"-0123456789");
   }
   
   return 1;

}

char* mpz_mat_to_string_pretty(mpz_mat_t mat)
{

   // estimate the size of the string
   // 4 + 3*r = enough room for null terminator, [,],\n and []\n per row
   unsigned long size = 4 + 3*mat->r;
   unsigned long i;
   for (i = 0; i < mat->r * mat->c; i++)
      // +2 is for the sign and a space
      size += mpz_sizeinbase(mat->entries[i], 10) + 2;

   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "[");
   for (i = 0; i < mat->r; i++)
   {
      *ptr = '[';
      ptr++;
      unsigned long j;
      for (j = 0; j < mat->c; j++)
      {
         mpz_get_str(ptr, 10, mat->entries[i*mat->c + j]);
         ptr += strlen(ptr);
	      if (j < mat->c - 1)
            {
            *ptr = ' ';
            ptr++;
            }
      }
      if (i != mat->r - 1)
         {
            *ptr = ']';
            ptr++;
            *ptr = '\n';
            ptr++;
         }
   }
   *ptr = ']';
   ptr++;
   *ptr = ']';
   ptr++;
   *ptr = '\n';
   ptr++;

   
   ptr--;
   *ptr = 0;
   
   return buf;

}

void mpz_mat_fprint(mpz_mat_t mat, FILE* f)
{
   char* s = mpz_mat_to_string(mat);
   fputs(s, f);
   free(s);
}

void mpz_mat_fprint_pretty(mpz_mat_t mat, FILE* f)
{
   char* s = mpz_mat_to_string_pretty(mat);
   fputs(s, f);
   free(s);
}

int mpz_mat_fread(mpz_mat_t mat, FILE* f)
{

   //read mat->rows
   unsigned long r;
   unsigned long c;

   if (!fscanf(f, "%ld %ld  ", &r, &c))
      return 0;

   mpz_mat_clear(mat);
   mpz_mat_init(mat,r,c);

   unsigned long i;
   for (i = 0; i < r*c; i++)
   {
      if (!mpz_inp_str(mat->entries[i], f, 10))
         return 0;
   }
   return 1;
}

int mpz_mat_fread_pretty(mpz_mat_t mat, FILE* f)
{

   unsigned long f_size;
   unsigned long lof=0;
   int ok;
   char* s;
   char c = ' ';

   fseek(f,0,SEEK_END);
   f_size = ftell(f);
   rewind(f);
   s = (char *)malloc(sizeof(char)*f_size+5);   


   c = fgetc(f);
   while(!feof(f))
   {
      s[lof] = c;
      c = fgetc(f);
      lof++;
   }
   s[lof] = 0;

   ok = mpz_mat_from_string_pretty(mat, s);

   return ok;
}

void mpz_mat_mul_classical(mpz_mat_t res, mpz_mat_t mat1, mpz_mat_t mat2)
{
   for (ulong i = 0; i < mat1->r; i++)
      for (ulong j = 0; j < mat2->c; j++)
      {
         mpz_set_ui(res->entries[i*res->c + j], 0);
         for (ulong k = 0; k < mat1->c; k++)
            mpz_addmul(res->entries[i*res->c + j], mat1->entries[i*mat1->c + k], mat2->entries[k*mat2->c + j]);
      }
}

// *************** end of file
