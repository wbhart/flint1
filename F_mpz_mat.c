/*============================================================================

    F_mpz_mat.c: Matrices over Z (FLINT 2.0 matrices)

    Copyright (C) 2008, William Hart 

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

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "flint.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"
#include "mpz_mat.h"

/*===============================================================================

	Memory management

================================================================================*/

void F_mpz_mat_init(F_mpz_mat_t mat, const ulong r, const ulong c)
{
   if ((r) && (c)) // allocate space for r*c small entries
   {
      mat->entries = (mp_limb_t *) flint_heap_alloc(r*c);
		F_mpn_clear(mat->entries, r*c); // zero all entries
		mat->rows = (F_mpz **) flint_heap_alloc(r); // initialise rows
		for (ulong i = 0; i < r; i++)
		   mat->rows[i] = mat->entries + i*c;
   } else mat->entries = NULL;
   
	mat->r = r;
	mat->c = c;
	mat->r_alloc = r;
	mat->c_alloc = c;
}

void F_mpz_mat_init_identity(F_mpz_mat_t mat, const ulong n)
{
   if (n) // allocate space for n*n small entries
   {
      mat->entries = (mp_limb_t *) flint_heap_alloc(n*n);
		F_mpn_clear(mat->entries, n*n); // zero all entries
		mat->rows = (F_mpz **) flint_heap_alloc(n); // initialise rows
		for (ulong i = 0; i < n; i++)
		{
			mat->rows[i] = mat->entries + i*n;
			F_mpz_set_ui(mat->rows[i] + i, 1); // set diagonal entries to 1
	   }
   } else mat->entries = NULL;
   
	mat->r = n;
	mat->c = n;
	mat->r_alloc = n;
	mat->c_alloc = n;
}

void F_mpz_mat_clear(F_mpz_mat_t mat)
{
   if (mat->entries) 
	{
		for (ulong i = 0; i < mat->r * mat->c; i++) 
			F_mpz_clear(mat->entries + i); // Clear all coefficients
		flint_heap_free(mat->entries); // clean up array of entries
		flint_heap_free(mat->rows); // clean up row array
	}
}

// Todo : add r_alloc_new which should be 8 more than r if new rows need to be allocated
// allocate that many rows instead of r rows

void F_mpz_mat_resize(F_mpz_mat_t mat, const ulong r, const ulong c)
{
	if (!mat->entries) // matrix has zero rows or columns, so we can just init
	{
		F_mpz_mat_init(mat, r, c);
		return;
	}
	
	if ((r == 0) || (c == 0)) // clear everything
	{
		for (ulong i = 0; i < mat->r; i++)
			for (ulong j = 0; j < mat->c; j++)
				F_mpz_zero(mat->rows[i] + j);

		mat->r = 0;
		mat->c = 0;

		return;
	}

	if (c <= mat->c_alloc) // don't need to alloc new columns
	{
		if (c < mat->c) // need to clear some columns
		{
			for (ulong i = 0; i < mat->r; i++)
            for (ulong j = c; j < mat->c; j++)
				   F_mpz_zero(mat->rows[i] + j);
		} 

		mat->c = c;
		
		if (r <= mat->r_alloc) // don't need to realloc rows
	   {
	      if (r < mat->r) // need to clear rows, but not realloc
		   {
			   for (long i = mat->r - 1; i >= r; i--)
				{
               F_mpz * start = mat->entries + i*mat->c_alloc; // Start of physically last row in memory
					
					if (mat->rows[i] == start) // row to be removed is physically last in memory
					{
					   for (ulong j = 0; j < mat->c; j++) // zero row to be removed
						   F_mpz_zero(mat->rows[i] + j);
					} else
					{
						ulong k;
						for (k = 0; k < mat->r - 1; k++) // find which row is physically last in memory
							if (mat->rows[k] == start) break;

                  for (ulong j = 0; j < mat->c; j++) // move row which is physically last with row i
						   F_mpz_swap(mat->rows[k] + j, mat->rows[i] + j);
						
						F_mpz * temp = mat->rows[k]; // swap pointers in row array
						mat->rows[k] = mat->rows[i];
						mat->rows[i] = temp;

                  for (ulong j = 0; j < mat->c; j++) // zero row to be removed which is now physically last
						   F_mpz_zero(start + j);
					}
				}
		   } 			
		} else // need to alloc new rows
		{
         F_mpz * old_entries = mat->entries; // save pointer to old entries

			mat->entries = (mp_limb_t *) flint_heap_realloc(mat->entries, r*mat->c_alloc);
			if (r > mat->r_alloc) mat->rows = (F_mpz **) flint_heap_realloc(mat->rows, r);
			mat->r_alloc = r;

			long diff = mat->entries - old_entries;
			for (ulong i = 0; i < mat->r; i++) // update row pointers
			   mat->rows[i] += diff;
 
			for (ulong i = mat->r; i < r; i++) // clear new rows
			{
				mat->rows[i] = mat->entries + mat->c_alloc*i; // add pointers to new rows
				for (ulong j = 0; j < mat->c_alloc; j++)
		         mat->rows[i][j] = 0L;
			}
		} 
	   
		mat->r = r;
	} else // need to alloc new columns
	{
		if (r*c > mat->r_alloc*mat->c_alloc) // need to realloc
		{
         F_mpz * old_entries = mat->entries; // save pointer to old entries

			mat->entries = (mp_limb_t *) flint_heap_realloc(mat->entries, r*c);
			if (r > mat->r_alloc) 
			   mat->rows = (F_mpz **) flint_heap_realloc(mat->rows, r);

			for (ulong i = mat->r_alloc*mat->c_alloc; i < r*c; i++) // clear new entries
			   mat->entries[i] = 0L;

			long diff = mat->entries - old_entries;
			for (ulong i = 0; i < mat->r; i++) // update row pointers
			   mat->rows[i] += diff;
		}

	   if (r < mat->r) // need to clear some existing rows
		{
			for (ulong i = mat->r - 1; i >= r; i--)
		   {
            F_mpz * start = mat->entries + i*mat->c_alloc; // Start of physically last row in memory
					
			   if (mat->rows[i] == start) // row to be removed is physically last in memory
				{
					for (ulong j = 0; j < mat->c; j++) // zero row to be removed
						F_mpz_zero(mat->rows[i] + j);
				} else
				{
					ulong k;
					for (k = 0; k < mat->r - 1; k++) // find which row is physically last in memory
						if (mat->rows[k] == start) break;

               for (ulong j = 0; j < mat->c; j++) // move row which is physically last with row i
						F_mpz_swap(mat->rows[k] + j, mat->rows[i] + j);
						
					F_mpz * temp = mat->rows[k]; // swap pointers in row array
					mat->rows[k] = mat->rows[i];
					mat->rows[i] = temp;

               for (ulong j = 0; j < mat->c; j++) // zero row to be removed which is now physically last
						F_mpz_zero(start + j);
				}
			}

			mat->r = r; // need to set mat->r so we know how many rows to move below
		}

		for (long i = mat->r - 1; i > 0L; i--) // move rows into new positions, 
			                                    // starting from physically last in memory 
															// first row stays where it is
		{
         F_mpz * start_old = mat->entries + i*mat->c_alloc;
			F_mpz * start_new = mat->entries + i*c;
			for (long j = mat->c - 1; j >= 0L; j--)
			{
				start_new[j] = start_old[j]; // copy entry data forwards
				start_old[j] = 0L; // zero old entry data
			}
		}

		// update old row pointers
		for (ulong i = 0; i < mat->r; i++)
		{ 
			ulong phys_row = (mat->rows[i] - mat->entries)/mat->c_alloc;
         mat->rows[i] = mat->entries + phys_row*c;
		}

		// put in new row pointers (if any)
		for (ulong i = mat->r; i < r; i++)
		   mat->rows[i] = mat->entries + c*i;

		mat->c_alloc = c;
		mat->r_alloc = r;
		mat->r = r;
		mat->c = c;
	}
}

/* ==============================================================================

   Input/Output

===============================================================================*/


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

   for (unsigned long i = 0; i < r*c; i++)
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
   for (unsigned long i = 0; i < mat->r * mat->c; i++)
      // +2 is for the sign and a space
      size += mpz_sizeinbase(mat->entries[i], 10) + 2;

   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "%ld %ld  ", mat->r, mat->c);
   for (unsigned long i = 0; i < mat->r * mat->c; i++)
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


   for(ulong i = 0; i < r*c; i++){
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
   for (unsigned long i = 0; i < mat->r * mat->c; i++)
      // +2 is for the sign and a space
      size += mpz_sizeinbase(mat->entries[i], 10) + 2;

   // write the string
   char* buf = (char*) malloc(size);
   char* ptr = buf + sprintf(buf, "[");
   for (unsigned long i = 0; i < mat->r; i++)
   {
      *ptr = '[';
      ptr++;
      for (unsigned long j = 0; j < mat->c; j++)
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

   for (unsigned long i = 0; i < r*c; i++)
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

int F_mpz_mat_from_string(F_mpz_mat_t mat, const char *s)
{

   int ok;
   
   mpz_mat_t m;
   mpz_mat_init(m,0,0);
   ok = mpz_mat_from_string(m, s);
   if (ok)
   {
      F_mpz_mat_clear(mat);
      F_mpz_mat_init(mat,m->r,m->c);
      mpz_mat_to_F_mpz_mat(mat, m);
   }
   mpz_mat_clear(m);
   
   return ok;
}

char* F_mpz_mat_to_string(F_mpz_mat_t mat)
{
   char* buf;
   mpz_mat_t m;
   mpz_mat_init(m,mat->r,mat->c);
   F_mpz_mat_to_mpz_mat(m, mat);
   buf = mpz_mat_to_string(m);
   mpz_mat_clear(m);
   return buf;
}

int F_mpz_mat_from_string_pretty(F_mpz_mat_t mat, char *s)
{

   int ok;
   
   mpz_mat_t m;
   mpz_mat_init(m,0,0);
   ok = mpz_mat_from_string_pretty(m, s);
   if (ok)
   {
      F_mpz_mat_clear(mat);
      F_mpz_mat_init(mat,m->r,m->c);
      mpz_mat_to_F_mpz_mat(mat, m);
   }
   mpz_mat_clear(m);
   
   return ok;
}

char* F_mpz_mat_to_string_pretty(F_mpz_mat_t mat)
{
   char* buf;
   mpz_mat_t m;
   mpz_mat_init(m, mat->r, mat->c);
   F_mpz_mat_to_mpz_mat(m, mat);
   buf = mpz_mat_to_string_pretty(m);
   mpz_mat_clear(m);
   return buf;
}

void F_mpz_mat_print(F_mpz_mat_t mat) 
{
   ulong i, j; 
   ulong r = mat->r;
   ulong c = mat->c;
	
   printf("%ld %ld  ", r, c);
   for (i = 0; i < r; i++) 
   {
      for (j = 0; j < c; j++) 
	   { 
	      F_mpz_print(mat->rows[i] + j); 
	      if (j < c - 1) printf(" "); 
	   }
      if (i != r - 1) printf(" "); 
   }  
   printf("\n"); 
}

void F_mpz_mat_print_pretty(F_mpz_mat_t mat) 
{
   ulong i, j; 
   ulong r = mat->r;
   ulong c = mat->c;
	
   printf("[");
   for (i = 0; i < r; i++) 
   {
      printf("[");
      for (j = 0; j < c; j++) 
	   { 
	      F_mpz_print(mat->rows[i] + j); 
	      if (j < c - 1) printf(" "); 
	   }
      if (i != r - 1) printf("]\n"); 
   }  
   printf("]]\n"); 
}

void F_mpz_mat_fprint(F_mpz_mat_t mat, FILE* f)
{
   char* s = F_mpz_mat_to_string(mat);
   fputs(s, f);
   free(s);
}

void F_mpz_mat_fprint_pretty(F_mpz_mat_t mat, FILE* f)
{
   char* s = F_mpz_mat_to_string_pretty(mat);
   fputs(s, f);
   free(s);
}

int F_mpz_mat_fread(F_mpz_mat_t mat, FILE* f)
{

   int ok;

   mpz_mat_t m;
   mpz_mat_init(m,0,0);
   
   ok = mpz_mat_fread(m, f);

   F_mpz_mat_clear(mat);
   F_mpz_mat_init(mat,m->r,m->c);

   mpz_mat_to_F_mpz_mat(mat,m);
   mpz_mat_clear(m);

   return ok;

}

/*===============================================================================

	Conversions

================================================================================*/

void mpz_mat_to_F_mpz_mat(F_mpz_mat_t F_mat, const mpz_mat_t m_mat)
{
	for (ulong r = 0; r < m_mat->r; r++)
	{
		ulong row = r*m_mat->c;
		for (ulong c = 0; c < m_mat->c; c++)
		   F_mpz_set_mpz(F_mat->rows[r] + c, m_mat->entries[row+c]);
	}
}

void F_mpz_mat_to_mpz_mat(mpz_mat_t m_mat, const F_mpz_mat_t F_mat)
{
	for (ulong r = 0; r < m_mat->r; r++)
	{
		ulong row = r*m_mat->c;
		for (ulong c = 0; c < m_mat->c; c++)
	      F_mpz_get_mpz(m_mat->entries[row+c], F_mat->rows[r] + c);
	}
}

long F_mpz_mat_set_line_d(double * appv, const F_mpz_mat_t mat, const ulong r, const int n)
{
   long * exp, i, maxexp = 0L;
   exp = (long *) malloc(n * sizeof(long)); 
  
   for (i = 0; i < n; i++)
   {
      appv[i] = F_mpz_get_d_2exp(&exp[i], mat->rows[r] + i);
      if (exp[i] > maxexp) maxexp = exp[i];
   }

   for (i = 0; i < n; i++) appv[i] = ldexp(appv[i], exp[i] - maxexp);

   free(exp);
   return maxexp;
}

/*===============================================================================

	Assignment/swap

================================================================================*/

void F_mpz_mat_set(F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{	
	if (mat1 != mat2) // trivial if aliased
	{
   	ulong r = mat1->r;
	   ulong c = mat1->c;
		
		for (ulong i = 0; i < r; i++)
			for (ulong j = 0; j < c; j++)
				F_mpz_set(mat1->rows[i] + j, mat2->rows[i] + j);
	}
}

/*===============================================================================

	Comparison

================================================================================*/

int F_mpz_mat_equal(const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   if (mat1 == mat2) return 1; // same matrix
	
	ulong r = mat1->r;
	ulong c = mat1->c;

	for (ulong i = 0; i < r; i++) // check if entries the same
		for (ulong j = 0; j < c; j++)
			if (!F_mpz_equal(mat1->rows[i] + j, mat2->rows[i] + j)) 
			   return 0;

	return 1;
}

/*===============================================================================

	Negation

================================================================================*/

void F_mpz_mat_neg(F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
	ulong r = mat1->r;
	ulong c = mat1->c;
	
	for (ulong i = 0; i < r; i++)
		for (ulong j = 0; j < c; j++)
		   F_mpz_neg(mat1->rows[i] + j, mat2->rows[i] + j);
}

/*===============================================================================

	Addition/subtraction

================================================================================*/

void F_mpz_mat_row_add(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, 
							  const ulong r1, const F_mpz_mat_t mat2, const ulong r2, 
							                         const ulong start, const ulong n)
{
   for (ulong i = start; i < start + n; i++) 
		F_mpz_add(res->rows[r3] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);   
}

void F_mpz_mat_row_sub(F_mpz_mat_t res, const ulong r3, const F_mpz_mat_t mat1, 
							  const ulong r1, const F_mpz_mat_t mat2, const ulong r2, 
							                          const ulong start, const ulong n)
{
   for (ulong i = start; i < start + n; i++) 
		F_mpz_sub(res->rows[r3] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);   
}

void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   ulong r = mat1->r;
	ulong c = mat1->c;
		
	for (ulong i = 0; i < r; i++) // add up to the length of the shorter mat
		for (ulong j = 0; j < c; j++)
			F_mpz_add(res->rows[i] + j, mat1->rows[i] + j, mat2->rows[i] + j);   
}

void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   ulong r = mat1->r;
	ulong c = mat1->c;
		
	for (ulong i = 0; i < r; i++) // add up to the length of the shorter mat
		for (ulong j = 0; j < c; j++)
			F_mpz_sub(res->rows[i] + j, mat1->rows[i] + j, mat2->rows[i] + j);   
}

/*===============================================================================

	Scalar multiplication

================================================================================*/

void F_mpz_mat_row_mul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)
{
	// either scalar of input mat is zero
	if (x == 0L)
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_zero(mat1->rows[r1] + i);

		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_set(mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_mul_ui(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

void F_mpz_mat_row_mul_si(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, long x)
{
	// either scalar of input mat is zero
	if (x == 0L)
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_zero(mat1->rows[r1] + i);

		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_set(mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	if (x == -1L) 
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_neg(mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_mul_si(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

void F_mpz_mat_row_mul_F_mpz(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, F_mpz_t x)
{
	// either scalar or input mat is zero
	if (*x == 0)
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_zero(mat1->rows[r1] + i);

		return;
	}
	
	// special cases, muliply by +/- 1
	if (*x == -1L)
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_neg(mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	if (*x == 1L)
	{
		for (ulong i = start; i < start + n; i++)
		   F_mpz_set(mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_mul2(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

/*===============================================================================

	Scalar addmul/submul

================================================================================*/

void F_mpz_mat_row_addmul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)
{
	// scalar is zero
	if (x == 0L)
		return;
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
			F_mpz_add(mat1->rows[r1] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_addmul_ui(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

void F_mpz_mat_row_addmul(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                             ulong start, ulong n, F_mpz_t x)
{
	// scalar is zero
	if ((*x) == 0L)
		return;
	
	// special case, multiply by 1
	if ((*x) == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
			F_mpz_add(mat1->rows[r1] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	// special case, multiply by -1
	if ((*x) == -1L) 
	{
	   for (ulong i = start; i < start + n; i++)
			F_mpz_sub(mat1->rows[r1] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_addmul(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

void F_mpz_mat_row_submul_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                               ulong start, ulong n, ulong x)
{
	// scalar is zero
	if (x == 0L)
		return;
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
			F_mpz_sub(mat1->rows[r1] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_submul_ui(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

void F_mpz_mat_row_submul(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								                             ulong start, ulong n, F_mpz_t x)
{
	// scalar is zero
	if ((*x) == 0L)
		return;
	
	// special case, multiply by 1
	if ((*x) == 1L) 
	{
	   for (ulong i = start; i < start + n; i++)
			F_mpz_sub(mat1->rows[r1] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	// special case, multiply by -1
	if ((*x) == -1L) 
	{
	   for (ulong i = start; i < start + n; i++)
			F_mpz_add(mat1->rows[r1] + i, mat1->rows[r1] + i, mat2->rows[r2] + i);

		return;
	}
	
	for (ulong i = start; i < start + n; i++)
		F_mpz_submul(mat1->rows[r1] + i, mat2->rows[r2] + i, x);
}

void F_mpz_mat_row_addmul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp)
{
	// scalar is zero, nothing to add
	if (c == 0)
	{
	   return;
	}
	
	// scalar is 1, just add 2^exp times the entry
	if (c == 1)
	{
	   F_mpz_t temp;
		F_mpz_init(temp);
		
		for (ulong i = start; i < start + n; i++)
		{
			F_mpz_mul_2exp(temp, mat2->rows[r2] + i, exp);
         F_mpz_add(mat1->rows[r1] + i, mat1->rows[r1] + i, temp);
		}

		F_mpz_clear(temp);

		return;
	}
	
	// exp is 1, just do addmul
	if (exp == 0)
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_addmul_ui(mat1->rows[r1] + i, mat2->rows[r2] + i, c);

		return;
	}
	
	F_mpz_t temp;
   F_mpz_init(temp);
		
	for (ulong i = start; i < start + n; i++)
	{
		F_mpz_mul_2exp(temp, mat2->rows[r2] + i, exp);
	   F_mpz_addmul_ui(mat1->rows[r1] + i, temp, c);
	}

	F_mpz_clear(temp);
}

void F_mpz_mat_row_submul_2exp_ui(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, ulong r2, 
								               ulong start, ulong n, ulong c, ulong exp)
{
	// scalar is zero, nothing to subtract
	if (c == 0)
	{
	   return;
	}
	
	// scalar is 1, just subtract 2^exp times the entry
	if (c == 1)
	{
	   F_mpz_t temp;
		F_mpz_init(temp);
		
		for (ulong i = start; i < start + n; i++)
		{
			F_mpz_mul_2exp(temp, mat2->rows[r2] + i, exp);
         F_mpz_sub(mat1->rows[r1] + i, mat1->rows[r1] + i, temp);
		}

		F_mpz_clear(temp);

		return;
	}
	
	// exp is 1, just do submul
	if (exp == 0)
	{
	   for (ulong i = start; i < start + n; i++)
		   F_mpz_submul_ui(mat1->rows[r1] + i, mat2->rows[r2] + i, c);

		return;
	}
	
	F_mpz_t temp;
   F_mpz_init(temp);
		
	for (ulong i = start; i < start + n; i++)
	{
		F_mpz_mul_2exp(temp, mat2->rows[r2] + i, exp);
	   F_mpz_submul_ui(mat1->rows[r1] + i, temp, c);
	}

	F_mpz_clear(temp);
}

void F_mpz_mat_row_swap(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, 
								                 ulong r2, ulong start, ulong n)
{
	if ((mat1 == mat2) && (n == mat1->c)) // matrices are the same, 
		                                   // just swap row pointers
	{
		if (r1 == r2) return; // rows are the same, nothing to do

		F_mpz * temp = mat1->rows[r1]; // swap rows via temporary
		mat1->rows[r1] = mat1->rows[r2];
		mat2->rows[r2] = temp;
	} else // swap entries in rows
	{
		for (ulong i = start; i < start + n; i++)
			F_mpz_swap(mat1->rows[r1] + i, mat2->rows[r2] + i);
	}
}

void F_mpz_mat_row_neg(F_mpz_mat_t mat1, ulong r1, F_mpz_mat_t mat2, 
								                 ulong r2, ulong start, ulong n)
{
	for (ulong i = start; i < start + n; i++)
		F_mpz_neg(mat1->rows[r1] + i, mat2->rows[r2] + i);
}
