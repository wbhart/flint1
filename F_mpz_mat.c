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
#include <mpfr.h>

#include "flint.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL_wrapper.h"
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
		ulong i;
		for (i = 0; i < r; i++)
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
		ulong i;
		for (i = 0; i < n; i++)
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
		ulong i;
		for (i = 0; i < mat->r * mat->c; i++) 
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
		ulong i, j;
		for (i = 0; i < mat->r; i++)
			for (j = 0; j < mat->c; j++)
				F_mpz_zero(mat->rows[i] + j);

		mat->r = 0;
		mat->c = 0;

		return;
	}

	if (c <= mat->c_alloc) // don't need to alloc new columns
	{
		if (c < mat->c) // need to clear some columns
		{
			ulong i, j;
			for (i = 0; i < mat->r; i++)
            for (j = c; j < mat->c; j++)
				   F_mpz_zero(mat->rows[i] + j);
		} 

		mat->c = c;
		
		if (r <= mat->r_alloc) // don't need to realloc rows
	   {
	      if (r < mat->r) // need to clear rows, but not realloc
		   {
			   long i;
			   for (i = mat->r - 1; i >= r; i--)
				{
               F_mpz * start = mat->entries + i*mat->c_alloc; // Start of physically last row in memory
					
					if (mat->rows[i] == start) // row to be removed is physically last in memory
					{
					   ulong j;
					   for (j = 0; j < mat->c; j++) // zero row to be removed
						   F_mpz_zero(mat->rows[i] + j);
					} else
					{
						ulong k;
						for (k = 0; k < mat->r - 1; k++) // find which row is physically last in memory
							if (mat->rows[k] == start) break;

                  ulong j;
                  for (j = 0; j < mat->c; j++) // move row which is physically last with row i
						   F_mpz_swap(mat->rows[k] + j, mat->rows[i] + j);
						
						F_mpz * temp = mat->rows[k]; // swap pointers in row array
						mat->rows[k] = mat->rows[i];
						mat->rows[i] = temp;

                  for (j = 0; j < mat->c; j++) // zero row to be removed which is now physically last
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
			ulong i;
			for (i = 0; i < mat->r; i++) // update row pointers
			   mat->rows[i] += diff;
 
			for (i = mat->r; i < r; i++) // clear new rows
			{
				mat->rows[i] = mat->entries + mat->c_alloc*i; // add pointers to new rows
				ulong j;
				for (j = 0; j < mat->c_alloc; j++)
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

			ulong i;
			for (i = mat->r_alloc*mat->c_alloc; i < r*c; i++) // clear new entries
			   mat->entries[i] = 0L;

			long diff = mat->entries - old_entries;
			for (i = 0; i < mat->r; i++) // update row pointers
			   mat->rows[i] += diff;
		}

	   if (r < mat->r) // need to clear some existing rows
		{
			ulong i;
			for (i = mat->r - 1; i >= r; i--)
		   {
            F_mpz * start = mat->entries + i*mat->c_alloc; // Start of physically last row in memory
					
			   if (mat->rows[i] == start) // row to be removed is physically last in memory
				{
					ulong j;
					for (j = 0; j < mat->c; j++) // zero row to be removed
						F_mpz_zero(mat->rows[i] + j);
				} else
				{
					ulong k;
					for (k = 0; k < mat->r - 1; k++) // find which row is physically last in memory
						if (mat->rows[k] == start) break;

               ulong j;
               for (j = 0; j < mat->c; j++) // move row which is physically last with row i
						F_mpz_swap(mat->rows[k] + j, mat->rows[i] + j);
						
					F_mpz * temp = mat->rows[k]; // swap pointers in row array
					mat->rows[k] = mat->rows[i];
					mat->rows[i] = temp;

               for (j = 0; j < mat->c; j++) // zero row to be removed which is now physically last
						F_mpz_zero(start + j);
				}
			}

			mat->r = r; // need to set mat->r so we know how many rows to move below
		}

		long i;
		for (i = mat->r - 1; i > 0L; i--) // move rows into new positions, 
			                                    // starting from physically last in memory 
															// first row stays where it is
		{
         F_mpz * start_old = mat->entries + i*mat->c_alloc;
			F_mpz * start_new = mat->entries + i*c;
			long j;
			for (j = mat->c - 1; j >= 0L; j--)
			{
				start_new[j] = start_old[j]; // copy entry data forwards
				start_old[j] = 0L; // zero old entry data
			}
		}

		// update old row pointers
		for (i = 0; i < mat->r; i++)
		{ 
			ulong phys_row = (mat->rows[i] - mat->entries)/mat->c_alloc;
         mat->rows[i] = mat->entries + phys_row*c;
		}

		// put in new row pointers (if any)
		for (i = mat->r; i < r; i++)
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

int F_mpz_mat_fread_pretty(F_mpz_mat_t mat, FILE* f)
{

   int ok;

   mpz_mat_t m;
   mpz_mat_init(m,0,0);
   
   ok = mpz_mat_fread_pretty(m, f);

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
	ulong r;
	for (r = 0; r < m_mat->r; r++)
	{
		ulong row = r*m_mat->c;
		ulong c;
		for (c = 0; c < m_mat->c; c++)
		   F_mpz_set_mpz(F_mat->rows[r] + c, m_mat->entries[row+c]);
	}
}

void F_mpz_mat_to_mpz_mat(mpz_mat_t m_mat, const F_mpz_mat_t F_mat)
{
	ulong r;
	for (r = 0; r < m_mat->r; r++)
	{
		ulong row = r*m_mat->c;
		ulong c;
		for (c = 0; c < m_mat->c; c++)
	      F_mpz_get_mpz(m_mat->entries[row+c], F_mat->rows[r] + c);
	}
}

long _F_mpz_vec_to_d_vec_2exp(double * appv, const F_mpz * vec, const ulong n)
{
   long * exp, i, maxexp = 0L;
   exp = (long *) malloc(n * sizeof(long)); 
  
   for (i = 0; i < n; i++)
   {
      appv[i] = F_mpz_get_d_2exp(&exp[i], vec + i);
      if (exp[i] > maxexp) maxexp = exp[i];
   }

   for (i = 0; i < n; i++) appv[i] = ldexp(appv[i], exp[i] - maxexp);

   free(exp);
   return maxexp;
}

void _F_mpz_vec_to_mpfr_vec(__mpfr_struct * appv, const F_mpz * vec, const ulong n)
{
   ulong i;
   for (i = 0; i < n; i++)
	  F_mpz_get_mpfr(appv + i, vec + i);

   return;
}

void _F_mpz_vec_2exp_to_mpfr_vec(__mpfr_struct * appv, const F_mpz * vec, const ulong n, int * cexpo)
{
   ulong i;
   for (i = 0; i < n; i++)
   {
      F_mpz_get_mpfr(appv + i, vec + i);
      mpfr_mul_2si(appv + i, appv + i, cexpo[i], GMP_RNDN);
   }

   return;
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
		
		ulong i, j;
		for (i = 0; i < r; i++)
			for (j = 0; j < c; j++)
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

	ulong i, j;
	for (i = 0; i < r; i++) // check if entries the same
		for (j = 0; j < c; j++)
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
	
	ulong i, j;
	for (i = 0; i < r; i++)
		for (j = 0; j < c; j++)
		   F_mpz_neg(mat1->rows[i] + j, mat2->rows[i] + j);
}

/*===============================================================================

	Addition/subtraction

================================================================================*/

void _F_mpz_vec_add(F_mpz * res, const F_mpz * vec1, 
					                     const F_mpz * vec2, const ulong n)
{
   ulong i;
   for (i = 0; i < n; i++) 
		F_mpz_add(res + i, vec1 + i, vec2 + i);   
}

void _F_mpz_vec_sub(F_mpz * res, const F_mpz * vec1, 
					                     const F_mpz * vec2, const ulong n)
{
   ulong i;
   for (i = 0; i < n; i++) 
		F_mpz_sub(res + i, vec1 + i, vec2 + i);   
}

void F_mpz_mat_add(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   ulong r = mat1->r;
	ulong c = mat1->c;
		
	ulong i, j;
	for (i = 0; i < r; i++) // add up to the length of the shorter mat
		for (j = 0; j < c; j++)
			F_mpz_add(res->rows[i] + j, mat1->rows[i] + j, mat2->rows[i] + j);   
}

void F_mpz_mat_sub(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   ulong r = mat1->r;
	ulong c = mat1->c;
		
	ulong i, j;
	for (i = 0; i < r; i++) // add up to the length of the shorter mat
		for (j = 0; j < c; j++)
			F_mpz_sub(res->rows[i] + j, mat1->rows[i] + j, mat2->rows[i] + j);   
}

/*===============================================================================

	Scalar multiplication

================================================================================*/

void _F_mpz_vec_mul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x)
{
	// either scalar of input mat is zero
	if (x == 0L)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_zero(vec1 + i);

		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_set(vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_mul_ui(vec1 + i, vec2 + i, x);
}

void _F_mpz_vec_mul_si(F_mpz * vec1, F_mpz * vec2, ulong n, long x)
{
	// either scalar of input mat is zero
	if (x == 0L)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_zero(vec1 + i);

		return;
	}
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_set(vec1 + i, vec2 + i);

		return;
	}
	
	if (x == -1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_neg(vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_mul_si(vec1 + i, vec2 + i, x);
}

void _F_mpz_vec_mul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x)
{
	// either scalar or input mat is zero
	if (*x == 0)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_zero(vec1 + i);

		return;
	}
	
	// special cases, muliply by +/- 1
	if (*x == -1L)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_neg(vec1 + i, vec2 + i);

		return;
	}
	
	if (*x == 1L)
	{
		ulong i;
		for (i = 0; i < n; i++)
		   F_mpz_set(vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_mul2(vec1 + i, vec2 + i, x);
}

/*===============================================================================

	Scalar addmul/submul

================================================================================*/

void _F_mpz_vec_addmul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x)
{
	// scalar is zero
	if (x == 0L)
		return;
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
			F_mpz_add(vec1 + i, vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_addmul_ui(vec1 + i, vec2 + i, x);
}

void _F_mpz_vec_addmul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x)
{
	// scalar is zero
	if ((*x) == 0L)
		return;
	
	// special case, multiply by 1
	if ((*x) == 1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
			F_mpz_add(vec1 + i, vec1 + i, vec2 + i);

		return;
	}
	
	// special case, multiply by -1
	if ((*x) == -1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
			F_mpz_sub(vec1 + i, vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_addmul(vec1 + i, vec2 + i, x);
}

void _F_mpz_vec_submul_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong x)
{
	// scalar is zero
	if (x == 0L)
		return;
	
	// special case, multiply by 1
	if (x == 1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
			F_mpz_sub(vec1 + i, vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_submul_ui(vec1 + i, vec2 + i, x);
}

void _F_mpz_vec_submul_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t x)
{
	// scalar is zero
	if ((*x) == 0L)
		return;
	
	// special case, multiply by 1
	if ((*x) == 1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
			F_mpz_sub(vec1 + i, vec1 + i, vec2 + i);

		return;
	}
	
	// special case, multiply by -1
	if ((*x) == -1L) 
	{
	   ulong i;
	   for (i = 0; i < n; i++)
			F_mpz_add(vec1 + i, vec1 + i, vec2 + i);

		return;
	}
	
	ulong i;
	for (i = 0; i < n; i++)
		F_mpz_submul(vec1 + i, vec2 + i, x);
}

void _F_mpz_vec_addmul_2exp_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong c, ulong exp)
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
		
		ulong i;
		for (i = 0; i < n; i++)
		{
			F_mpz_mul_2exp(temp, vec2 + i, exp);
         F_mpz_add(vec1 + i, vec1 + i, temp);
		}

		F_mpz_clear(temp);

		return;
	}
	
	// exp is 1, just do addmul
	if (exp == 0)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_addmul_ui(vec1 + i, vec2 + i, c);

		return;
	}
	
	F_mpz_t temp;
   F_mpz_init(temp);
		
	ulong i;
	for (i = 0; i < n; i++)
	{
		F_mpz_mul_2exp(temp, vec2 + i, exp);
	   F_mpz_addmul_ui(vec1 + i, temp, c);
	}

	F_mpz_clear(temp);
}

void _F_mpz_vec_submul_2exp_ui(F_mpz * vec1, F_mpz * vec2, ulong n, ulong c, ulong exp)
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
		
		ulong i;
		for (i = 0; i < n; i++)
		{
			F_mpz_mul_2exp(temp, vec2 + i, exp);
         F_mpz_sub(vec1 + i, vec1 + i, temp);
		}

		F_mpz_clear(temp);

		return;
	}
	
	// exp is 1, just do submul
	if (exp == 0)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_submul_ui(vec1 + i, vec2 + i, c);

		return;
	}
	
	F_mpz_t temp;
   F_mpz_init(temp);
		
	ulong i;
	for (i = 0; i < n; i++)
	{
		F_mpz_mul_2exp(temp, vec2 + i, exp);
	   F_mpz_submul_ui(vec1 + i, temp, c);
	}

	F_mpz_clear(temp);
}

void _F_mpz_vec_submul_2exp_F_mpz(F_mpz * vec1, F_mpz * vec2, ulong n, F_mpz_t c, ulong exp)
{
	// scalar is zero, nothing to subtract
	if (F_mpz_is_zero(c))
	{
	   return;
	}
	
	// scalar is 1, just subtract 2^exp times the entry
	if (F_mpz_is_one(c) || F_mpz_is_m1(c))
	{
	   F_mpz_t temp;
		F_mpz_init(temp);
		
		if (F_mpz_sgn(c) > 0)
      {
         ulong i;
         for (i = 0; i < n; i++)
		   {
			   F_mpz_mul_2exp(temp, vec2 + i, exp);
            F_mpz_sub(vec1 + i, vec1 + i, temp);
		   }
      } else
      {
         ulong i;
         for (i = 0; i < n; i++)
		   {
			   F_mpz_mul_2exp(temp, vec2 + i, exp);
            F_mpz_add(vec1 + i, vec1 + i, temp);
		   }
      }

		F_mpz_clear(temp);

		return;
	}
	
	// exp is 1, just do submul
	if (exp == 0)
	{
	   ulong i;
	   for (i = 0; i < n; i++)
		   F_mpz_submul(vec1 + i, vec2 + i, c);

		return;
	}
	
	F_mpz_t temp;
   F_mpz_init(temp);
		
	ulong i;
	for (i = 0; i < n; i++)
	{
		F_mpz_mul_2exp(temp, vec2 + i, exp);
	   F_mpz_submul(vec1 + i, temp, c);
	}

	F_mpz_clear(temp);
}

/* ======================================================================================================

 Classical Multiplication

=========================================================================================================*/

void _F_mpz_mat_mul_classical(F_mpz_mat_t res, const F_mpz_mat_t mat1, const F_mpz_mat_t mat2)
{
   //res=mat1*mat2

   ulong r1 = mat1->r;
   ulong c1 = mat1->c;
   ulong r2 = mat2->r;
   ulong c2 = mat2->c;

   if (c1!=r2)
      return; //dimensions don't match up

   F_mpz_mat_clear(res);
   F_mpz_mat_init(res,r1,c2);

   ulong i, j, c;
   for (i = 0; i < r1; i++) // add up to the length of the shorter mat
      for (j = 0; j < c2; j++)
         for (c=0; c<c1;c++)
            F_mpz_addmul(res->rows[i] + j, mat1->rows[i]+c, mat2->rows[c]+j);   

}

/*============================================================================

   assorted F_mpz_mat functions

=============================================================================*/

#define MAT_POS_UPDATE \
   do { \
      if (row != NULL) \
      { \
         (*row) = i; \
		 (*col) = j; \
      } \
   } while (0)

long F_mpz_mat_max_bits2(ulong * row, ulong * col, const F_mpz_mat_t M)
{
   int sign = 0;
   ulong max = 0;
   ulong bits = 0;
   ulong max_limbs = 1;
   ulong size;
   ulong i, j;
   F_mpz c;

   if (row != NULL)
   {
      *row = 0;
	  *col = 0;
   }

   // search until we find an mpz_t coefficient or one of at least FLINT_BITS - 2 bits
   for (i = 0; i < M->r; i++) 
   {
      for (j = 0; j < M->c; j++)
	  {
		 c = M->rows[i][j];
		 if (COEFF_IS_MPZ(c)) goto check_mpzs; // found an mpz_t coeff
         if (c < 0L) 
	     {
		    sign = 1;
            bits = FLINT_BIT_COUNT(-c);
	     } else bits = FLINT_BIT_COUNT(c);
	     if (bits > max) 
	     {
		    max = bits;
			MAT_POS_UPDATE;
	        if (max >= FLINT_BITS - 2) goto check_mpzs; // coeff is at least FLINT_BITS - 2 bits
	     }
	  }
	}

check_mpzs:

    // search through mpz coefficients for largest size in bits
	
	for ( ; i < M->r; i++)
   	{
	   for ( ; j < M->c; j++)
   	   {
		  c =  M->rows[i][j];
          if (COEFF_IS_MPZ(c))
	   	  {
	   		 __mpz_struct * mpz_ptr = F_mpz_ptr_mpz(c);
	   		 if (mpz_sgn(mpz_ptr) < 0) sign = 1;
	   		 size = mpz_size(mpz_ptr);
	   		 if (size > max_limbs)
	   		 {
	   		    max_limbs = size;
	   			mp_limb_t * data = mpz_ptr->_mp_d;
	   		    bits = FLINT_BIT_COUNT(data[max_limbs - 1]);
				MAT_POS_UPDATE;
	   			max = bits;
	   		 } else if (size == max_limbs)
	   		 {
	   			mp_limb_t * data = mpz_ptr->_mp_d;
	   		    bits = FLINT_BIT_COUNT(data[max_limbs - 1]);
	   		    if (bits > max) 
				{
				   max = bits;
				   MAT_POS_UPDATE;
				}
	   		 }
	   	  } else if ((long) c < 0L) sign = 1; // still need to check the sign of small coefficients
	   }
	   j = 0;
	}
	
	if (sign) return -(max + FLINT_BITS*(max_limbs - 1));
	else return max + FLINT_BITS*(max_limbs - 1);
}

void F_mpz_mat_mul_2exp(F_mpz_mat_t res, F_mpz_mat_t M, ulong n)
{
   if (n == 0)
   {
      F_mpz_mat_set(res, M);
      return;
   }

   if (n == 1)
   {
      F_mpz_mat_add(res, M, M);
      return;
   }

   ulong i, j;
   for (i = 0; i < M->r; i++)
      for (j = 0; j < M->c; j++)
         F_mpz_mul_2exp(res->rows[i] + j, M->rows[i] + j, n);   
}

void F_mpz_mat_div_2exp(F_mpz_mat_t res, F_mpz_mat_t M, ulong n)
{
   if (n == 0)
   {
      F_mpz_mat_set(res, M);
      return;
   }

   ulong i, j;
   for (i = 0; i < M->r; i++)
      for (j = 0; j < M->c; j++)
         F_mpz_div_2exp(res->rows[i] + j, M->rows[i] + j, n);

   return;
}

void _F_mpz_vec_scalar_product(F_mpz_t sp, F_mpz * vec1, F_mpz * vec2, ulong n)
{
	ulong i;
   
   F_mpz_mul2(sp, vec1, vec2);
   
   for (i = 1; i < n; i++)
      F_mpz_addmul(sp, vec1 + i, vec2 + i);

   return;
}

int F_mpz_mat_col_equal(F_mpz_mat_t M, ulong a, ulong b)
{
   ulong i;
   
   for (i = 0; i < M->r; i++)
      if (!F_mpz_equal(M->rows[i] + a, M->rows[i] + b))
	     return 0;

   return 1;
}

void F_mpz_mat_col_copy(F_mpz_mat_t M, ulong a, ulong b)
{
   ulong i;

   for (i = 0; i < M->r; i++)
      F_mpz_set(M->rows[i] + a, M->rows[i] + b);
}

typedef struct
{
   ulong col;
   ulong hash;
} col_hash_t;

void F_mpz_mat_col_hash(col_hash_t * colh, F_mpz_mat_t M)
{
   ulong i, j, hash;
   for (i = 0; i < M->c; i++)
   {
      colh[i].col = i;
      hash = 0;
	  for (j = 0; j < M->r; j++)
	     hash += ((j+1)*F_mpz_get_ui(M->rows[j] + i));
	  colh[i].hash = hash;
   }
}

void F_mpz_mat_hash_bubble_sort(col_hash_t * colh, ulong cols)
{
   int swapped;
   ulong i;
   col_hash_t t;

   do
   {
      swapped = 0;
      for (i = 1; i < cols; i++)
	  {
	     if (colh[i-1].hash < colh[i].hash) // need to swap
		 {
		    t = colh[i-1]; // swap
			colh[i-1] = colh[i];
			colh[i] = t;
			swapped = 1;
		 }
	  }
   } while (swapped);
}

int F_mpz_mat_col_partition(ulong * part, F_mpz_mat_t M)
{
   ulong start = 0, upto = 1, p = 0;
   
   col_hash_t * colh = malloc(sizeof(col_hash_t)*M->c);
   F_mpz_mat_col_hash(colh, M);
   F_mpz_mat_hash_bubble_sort(colh, M->c);
   
   while (start < M->c)
   {
	  p++;
	  if (p > M->r) return 0; // already too many partitions

      part[colh[start].col] = p;
	  for (upto = start + 1; upto < M->c; upto++)
	  {
		  if (!F_mpz_mat_col_equal(M, colh[start].col, colh[upto].col)) break;
		  part[colh[upto].col] = p;
	  }

      start = upto;
   }

   free(colh);
   return p;
}

void F_mpz_mat_window_init(F_mpz_mat_t U, F_mpz_mat_t M, 
                           ulong r0, ulong c0, ulong rows, ulong cols)
{
   ulong i;

   U->r = rows;
   U->c = cols;
   if (rows)
      U->rows = malloc(sizeof(F_mpz *)*rows);
   
   for (i = 0; i < rows; i++)
      U->rows[i] = M->rows[r0 + i] + c0;
}

void F_mpz_mat_smod(F_mpz_mat_t res, F_mpz_mat_t M, F_mpz_t P)
{
   ulong i, j;
   
   if (F_mpz_is_zero(P))
   {
      printf("FLINT Exception: Division by zero\n");
      abort();
   }

   if (F_mpz_is_one(P))
   {
      F_mpz_mat_clear(res);
      F_mpz_mat_init(res, M->r, M->c);
      return;
   }

   mp_ptr Pinv = F_mpz_precompute_inverse(P);
   
   for (i = 0; i < M->r; i++)
      for (j = 0; j < M->c; j++)
         F_mpz_smod_preinv(res->rows[i] + j, M->rows[i] + j, P, Pinv);

   F_mpz_preinv_clear(Pinv);
}

int _trunc_col_test(F_mpz_mat_t temp_col, F_mpz_t trunc_P, long max_bits)
{
   F_mpz_t sum;
   F_mpz_init(sum);

   F_mpz_set_ui(sum, 0L);

   long i;
   for (i = 0; i < temp_col->r; i++)
      F_mpz_add(sum, sum, temp_col->rows[i]);

   F_mpz_smod(sum, sum, trunc_P);

   if (F_mpz_bits(sum) > max_bits)
   {
      F_mpz_clear(sum);
      return 0;
   }

   return 1;
}

/*
   Goal here is to take a matrix M, get U, multiply U by col look at max bits of U*col and P 
   subtract exp and decide if it's worth calling LLL
   if is not return 0
   U_exp is the assumed power of the scalar multiple of U (so 2^U_exp is the scalar weight of U) 
   this should match the virtual precision.
   if it is then return weight of last column (check that it should not be zero... ??) and 
   augment M with the new column and new row
   This new column should be truncated to the correct amount before re-multiplying by U
   first make sure there are enough bits to even bother
*/
int _F_mpz_mat_next_col(F_mpz_mat_t M, F_mpz_t P, F_mpz_mat_t col, long exp, long U_exp)
{   
   ulong r = col->r;
   ulong B = r + 2;
   ulong bit_r = FLINT_MAX(r, 20);
   
   long ISD = F_mpz_bits(P) - bit_r - bit_r/2;
   if (ISD < exp)
      return 0;

   F_mpz_mat_t U;
   F_mpz_mat_window_init(U, M, 0, 0, M->r, r);
   
   F_mpz_mat_t temp_col;
   F_mpz_mat_init(temp_col, M->r, 1);

   // full precision column for deciding truncation levels
   F_mpz_mat_mul_classical(temp_col, U, col);

   if (U_exp >= 0)
      F_mpz_mat_div_2exp(temp_col, temp_col, (ulong) U_exp);
   else
      F_mpz_mat_mul_2exp(temp_col, temp_col, (ulong) (-1*U_exp));
    
   F_mpz_mat_smod(temp_col, temp_col, P);

   long mbts = FLINT_ABS(F_mpz_mat_max_bits(temp_col));

   // bare minimum of data above the bound
   if (mbts < (long) (0.973 * (double) bit_r - 0.1 + (double) exp))
   {
      F_mpz_mat_clear(temp_col);
      F_mpz_mat_window_clear(U);      
      return 0;
   }

   // Now do a test to see if we can just skip the p^a vector
   // double version of no vector test:
   double S = B * (ldexp(1.51, M->r) * (1/0.51) - (1/0.51)- 1);
   
   F_mpz_t temp;
   F_mpz_init(temp);
   
   F_mpz_set_d_2exp(temp, (double) B, ISD);
   F_mpz_sub(temp, P, temp);

   // this is the bits version of K*B*(2*(3/2)^(s-1)-2) < P - B*2^ISD check
   double no_vec_check = (double) F_mpz_bits(temp) - (1 + log2(S) + mbts); 
   F_mpz_clear(temp);
   
   int no_vec = (no_vec_check > 0.0);
   
   ulong prec = U_exp; // taken as an argument to match the scaling on U;
   
   long take_away;
   int virt_exp;
   
   F_mpz_mat_t trunc_col;
   F_mpz_mat_init(trunc_col, col->r, col->c);

   // rare for the first time, frequent for repeated scalings
   // In here we're going to scale to make the new entries use their 2*r bits
   // Now decide the scaling based on mbts... should be mbits - 0.973r
   // maybe have mbits be more precise... 
   if (no_vec)
      ISD = mbts - (long)( .973 * (double) bit_r - .1);

   take_away = ISD - prec;
   virt_exp = -prec;

   if (take_away >= 0)
      F_mpz_mat_div_2exp(trunc_col, col, (ulong) take_away);
   else
      F_mpz_mat_mul_2exp(trunc_col, col, (ulong) (-1*take_away)); 

   // Here we'll run with some extra bits with the new vector, 
   // just for complexity's sake
   // This means using ISD as the scale down.  Since we want to truncate we will 
   // take away ISD - s bits then return -s as the virtual exponent
   F_mpz_mat_mul_classical(temp_col, U, trunc_col);

   if (U_exp >= 0)
      F_mpz_mat_div_2exp(temp_col, temp_col, (ulong) U_exp);
   else
      F_mpz_mat_mul_2exp(temp_col, temp_col, (ulong) (-1*U_exp));

   F_mpz_t trunc_P;
   F_mpz_init(trunc_P);
      
   ulong new_r = M->r + !(no_vec);
   ulong old_r = M->r;
   long i;
   
   F_mpz_mat_resize(M, new_r, M->c + 1); 
   for (i = old_r - 1; i >= 0; i--)
      F_mpz_mat_swap_rows(M, i + new_r - old_r, i);

   if (take_away >= 0)
      F_mpz_div_2exp(trunc_P, P, (ulong) take_away);
   else
      F_mpz_mul_2exp(trunc_P, P, (ulong) (-1)*take_away);         
   
   if (!F_mpz_is_zero(trunc_P))
      F_mpz_mat_smod(temp_col, temp_col, trunc_P);

   FLINT_ASSERT(_trunc_col_test(col, P, exp));

   if (!no_vec)
      F_mpz_set(M->rows[0] + M->c - 1, trunc_P);

   for (ulong j = !(no_vec); j < M->r; j++)
      F_mpz_set(M->rows[j] + M->c - 1, temp_col->rows[j - !(no_vec)]);
   
   F_mpz_clear(trunc_P);
   F_mpz_mat_clear(trunc_col);
   F_mpz_mat_clear(temp_col);
   F_mpz_mat_window_clear(U);

   return virt_exp;
}


/*
   We are given the same column we just ran with.
   Last time we used 1.5*col->r bits. 
   Here we want to add more data to see if we can get enough to make 
   the later Gram-Schmidt lengths large without upsetting the stability 
   of the G-S computation.
   If we can do this cheaply its beter than running LLL at full size.  
   We'll return the new dimension with our findings, or -1 which means 
   that there's not enough new data.
*/
int F_mpz_mat_check_rest(F_mpz_mat_t M, F_mpz_t P, F_mpz_mat_t col, long exp, long U_exp, F_mpz_t B)
{

   ulong r = col->r;
   ulong bit_r = FLINT_MAX(r, 20);
   
   F_mpz_mat_t U;
   F_mpz_mat_window_init(U, M, 0, 0, M->r, r);
   
   F_mpz_mat_t temp_col;
   F_mpz_mat_init(temp_col, M->r, 1);

   // full precision column for deciding truncation levels
   F_mpz_mat_mul_classical(temp_col, U, col);

/*   if (U_exp >= 0)
      F_mpz_mat_div_2exp(temp_col, temp_col, (ulong) U_exp);
   else
      F_mpz_mat_mul_2exp(temp_col, temp_col, (ulong) (-1*U_exp));
*/    
   F_mpz_mat_smod(temp_col, temp_col, P);

   long mbts = FLINT_ABS(F_mpz_mat_max_bits(temp_col));

   // bare minimum of data above the bound
   if (mbts < exp)
   {
      F_mpz_mat_clear(temp_col);
      F_mpz_mat_window_clear(U);      
      return -1;
   }

   F_mpz_mat_t new_col;
   F_mpz_mat_init(new_col, M->r, 1);

   long i;
   for( i = 0; i < M->r; i++)
   {
      F_mpz_set(temp_col->rows[i],  M->rows[i] + M->c - 1);
   }

   if (U_exp >= 0)
      F_mpz_mat_mul_2exp(new_col, temp_col, (ulong) U_exp);
   else
      F_mpz_mat_div_2exp(new_col, temp_col, (ulong) (-1*U_exp));

   for( i = 0; i < M->r; i++)
   {
      F_mpz_set(M->rows[i] + M->c - 1, new_col->rows[i]);
   }

   ulong newd;
   newd = F_mpz_mat_gs_d(M, B);

   for (i = 0; i < M->r; i++)
      F_mpz_set(M->rows[i] + M->c - 1, temp_col->rows[i]);

   F_mpz_mat_clear(new_col);
   F_mpz_mat_clear(temp_col);
   F_mpz_mat_window_clear(U); 

   return newd;
}

