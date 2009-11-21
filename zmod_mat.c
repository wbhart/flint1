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
/*****************************************************************************

   zmod_mat.c: Matrices over (unsigned) long mod p, for p prime.
   
   Copyright (C) 2008, 2009 William Hart.
   Copyright (C) 2008, Richard Howell-Peak
   Copyright (C) 2008, Martin Albrecht <M.R.Albrecht@rhu.ac.uk> (Some Strassen code)

*****************************************************************************/

#include "zmod_mat.h"
#include "zmod_poly.h"
#include "long_extras.h"
#include "flint.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void zmod_mat_init(zmod_mat_t mat, ulong p, ulong rows, ulong cols)
{
   zmod_mat_init_precomp(mat, p, z_precompute_inverse(p), rows, cols);
}

void zmod_mat_init_precomp(zmod_mat_t mat, ulong p, double p_inv, 
						                             ulong rows, ulong cols)
{
   mat->arr = (ulong **) flint_heap_alloc(rows);
   mat->ptr = (ulong *) flint_heap_alloc(rows*cols);
   
   // Set up the rows
   ulong * ptr = mat->ptr;
   for (ulong i = 0; i < rows; i++, ptr += cols)
      mat->arr[i] = ptr;

   mat->p = p;
   mat->p_inv = p_inv;
   mat->rows = rows;
   mat->cols = cols;
}

void zmod_mat_clear(zmod_mat_t mat)
{
   flint_heap_free(mat->ptr);
   flint_heap_free(mat->arr);
}

/*******************************************************************************************

   Matrix Windows

*******************************************************************************************/

void zmod_mat_window_init(zmod_mat_t out, zmod_mat_t mat, ulong r1, ulong c1, ulong r2, ulong c2)
{
   out->arr = (ulong **) flint_heap_alloc(r2 - r1);

   for (ulong i = 0; i < r2 - r1; i++)
      out->arr[i] = mat->arr[r1 + i] + c1;

   out->p = mat->p;
   out->p_inv = mat->p_inv;
   out->rows = r2 - r1;
   out->cols = c2 - c1;
}

void zmod_mat_window_clear(zmod_mat_t mat)
{
   flint_heap_free(mat->arr);
}

/*******************************************************************************************

   Conversions

*******************************************************************************************/

/*
   Set a row to the coefficients of a polynomial, starting with the constant coefficient
   Assumes that poly->length <= mat->cols
*/
void zmod_poly_to_zmod_mat_row(zmod_mat_t mat, ulong row, zmod_poly_t poly)
{
   ulong * r1 = mat->arr[row];
   ulong * coeffs = poly->coeffs;
   ulong cols = mat->cols;

   long i;
   
   for (i = 0; i < poly->length; i++)
      r1[i] = coeffs[i];
   
   for ( ; i < cols; i++)
      r1[i] = 0L;
}
/*
   Set a column to the coefficients of a polynomial, starting with the constant coefficient
   Assumes that poly->length <= mat->rows
*/

void zmod_poly_to_zmod_mat_col(zmod_mat_t mat, ulong col, zmod_poly_t poly)
{
   ulong * r1; 
   ulong * coeffs = poly->coeffs;
   ulong rows = mat->rows;

   long i;
   
   for (i = 0; i < poly->length; i++)
   {
	  r1 = mat->arr[i];
	  r1[col] = coeffs[i];
   }
   
   for ( ; i < rows; i++)
   {
	  r1 = mat->arr[i];
	  r1[col] = 0L;
   }
}

/*
   Set a zmod_poly's coefficients to the entries in a column, starting with the constant coefficient
*/
void zmod_mat_col_to_zmod_poly(zmod_poly_t poly, zmod_mat_t mat, ulong col)
{
   ulong rows = mat->rows;
   ulong * ptr;
   
   zmod_poly_fit_length(poly, rows);
   for (ulong i = 0; i < rows; i++)
   {  
	  ptr = mat->arr[i];
      poly->coeffs[i] = ptr[col];
   }

   poly->length = rows;
   __zmod_poly_normalise(poly);
}

/*
   Set a zmod_poly's coefficients to the entries in a column, starting with the constant coefficient
   but shifting along by one for every non-zero entry in shift
*/
void zmod_mat_col_to_zmod_poly_shifted(zmod_poly_t poly, zmod_mat_t mat, ulong col, ulong * shift)
{
   ulong rows = mat->rows;
   ulong * ptr;
   
   zmod_poly_fit_length(poly, rows);
   for (ulong i = 0, j = 0; j < rows; j++)
   {  
	  if (shift[j]) poly->coeffs[j] = 0L;
	  else
	  {
		 ptr = mat->arr[i];
         poly->coeffs[j] = ptr[col];
	     i++;
	  }
   }

   poly->length = rows;
   __zmod_poly_normalise(poly);
}

/*******************************************************************************************

   Elementary row operations

*******************************************************************************************/

/*
   row1 = row1 + u * row2
   only columns [start..mat->cols) are affected
   assumes u is reduced mod mat->p
*/
void zmod_mat_row_scalar_addmul_right(zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start)
{
   ulong * r1 = mat->arr[row1];
   ulong * r2 = mat->arr[row2];
   ulong p = mat->p;
   double p_inv = mat->p_inv;
   ulong prod;
   ulong cols = mat->cols;

#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) >= FLINT_D_BITS)
   {
	  for (ulong i = start; i < cols; i++)
      {
	     prod = z_mulmod2_precomp(r2[i], u, p, p_inv);
	     r1[i] = z_addmod(r1[i], prod, p);
      }
   } else
#endif
   {
	  for (ulong i = start; i < cols; i++)
      {
	     prod = z_mulmod_precomp(r2[i], u, p, p_inv);
	     r1[i] = z_addmod(r1[i], prod, p);
      }
   }
}

/*
   row1 = row1 - u * row2
   only columns [start..mat->cols) are affected
   assumes u is reduced mod mat->p
*/
void zmod_mat_row_scalar_submul_right(zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start)
{
   ulong * r1 = mat->arr[row1];
   ulong * r2 = mat->arr[row2];
   ulong p = mat->p;
   double p_inv = mat->p_inv;
   ulong prod;
   ulong cols = mat->cols;

#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) >= FLINT_D_BITS)
   {
	  for (ulong i = start; i < cols; i++)
      {
	     prod = z_mulmod2_precomp(r2[i], u, p, p_inv);
	     r1[i] = z_submod(r1[i], prod, p);
      }
   } else
#endif
   {
	  for (ulong i = start; i < cols; i++)
      {
	     prod = z_mulmod_precomp(r2[i], u, p, p_inv);
	     r1[i] = z_submod(r1[i], prod, p);
      }
   }
}

/*
   row = row * u
   only columns [start..mat->cols) are affected
   assumes u is reduced mod mat->p
*/
void zmod_mat_row_scalar_mul_right(zmod_mat_t mat, ulong row, ulong u, ulong start)
{
   ulong * r1 = mat->arr[row];
   ulong cols = mat->cols;
   ulong p = mat->p;
   double p_inv = mat->p_inv;

#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) >= FLINT_D_BITS)
   {
	  for (ulong i = start; i < cols; i++)
	     r1[i] = z_mulmod2_precomp(r1[i], u, p, p_inv);
   } else
#endif
   {
	  for (ulong i = start; i < cols; i++)
	     r1[i] = z_mulmod_precomp(r1[i], u, p, p_inv);
   }
}

/*
   r1 = r2 * u
   only columns [start..end) are affected
   assumes u is reduced mod p
*/
void zmod_vec_scalar_mul_range(ulong * r1, ulong * r2, ulong u, ulong p, double p_inv, ulong start, ulong end)
{

#if FLINT_BITS == 64
   if (FLINT_BIT_COUNT(p) >= FLINT_D_BITS)
   {
	  for (ulong i = start; i < end; i++)
	     r1[i] = z_mulmod2_precomp(r2[i], u, p, p_inv);
   } else
#endif
   {
	  for (ulong i = start; i < end; i++)
	     r1[i] = z_mulmod_precomp(r2[i], u, p, p_inv);
   }
}

/*
   r1 = r1 - r2
   only columns [start..end) are affected
*/
void zmod_vec_sub_range(ulong * r1, ulong * r2, ulong p, ulong start, ulong end)
{
   for (ulong i = start; i < end; i++)
   {
      r1[i] = z_submod(r1[i], r2[i], p);
   }
}

/*******************************************************************************************

   Comparison

*******************************************************************************************/

int zmod_mat_equal(zmod_mat_t mat1, zmod_mat_t mat2)
{
   for (ulong i = 0; i < mat1->rows; i++)
      for (ulong j = 0; j < mat1->cols; j++)
         if (mat1->arr[i][j] != mat2->arr[i][j]) return 0;

   return 1;
}

/*******************************************************************************************

   Input/output

*******************************************************************************************/

void zmod_mat_print(zmod_mat_t mat)
{
   if (mat->rows == 0) 
   {
      printf("[[]]");
	  return;
   }

   if (mat->cols == 0)
   {
      printf("[");
	  for (ulong i = 0; i < mat->rows - 1; i++)
	  {
	     printf("[]\n");
	  }
	  printf("[]]");
	  return;
   }
   
   ulong * ptr;
   printf("[");
   for (ulong i = 0; i < mat->rows - 1; i++)
   {
      ptr = mat->arr[i];
	  printf("[");
	  for (ulong j = 0; j < mat->cols - 1; j++)
	  {
	     printf("%ld, ", ptr[j]);
	  }
	  printf("%ld]\n", ptr[mat->cols - 1]);
   }
   ptr = mat->arr[mat->rows - 1];
   printf("[");
   for (ulong j = 0; j < mat->cols - 1; j++)
   {
	  printf("%ld, ", ptr[j]);
   }
   printf("%ld]]", ptr[mat->cols - 1]);
}

/*******************************************************************************************

   Reduction

*******************************************************************************************/

/*
   Gaussian Elimination
   After the algorithm has run, A will be in upper-triangular form
   Returns the rank of the matrix
   Assumes the modulus is a very small prime
 */
ulong zmod_mat_row_reduce_gauss_small(zmod_mat_t mat)
{
	ulong i = 0, j = 0, k;
    ulong rows = mat->rows;
	ulong cols = mat->cols;
	ulong coeff;
	ulong p = mat->p;
	double p_inv = mat->p_inv;
	ulong * temp = (ulong *) flint_heap_alloc(cols);

	while ((i < rows) && (j < cols))
	{
		for(k = i; k < rows; k++)
		   if (zmod_mat_get_coeff_ui(mat, k, j)) break;
		
		if (k < rows)
		{
			if (k != i) zmod_mat_swap_rows(mat, i, k);
			ulong n = zmod_mat_get_coeff_ui(mat, i, j);
			ulong n_inv = z_invert(n, p);
			zmod_mat_row_scalar_mul_right(mat, i, n_inv, j);
			for (ulong lead = 1; lead < p; lead++)
			{
			   zmod_vec_scalar_mul_range(temp, mat->arr[i], lead, p, p_inv, j, cols);
			   for(ulong u = i + 1; u < rows; u++)
			   {
			      if (zmod_mat_get_coeff_ui(mat, u, j) == lead)
			          zmod_vec_sub_range(mat->arr[u], temp, p, j, cols);
			   }
			}
			i++;
		}
		j++;
	}

	flint_heap_free(temp);
	return i;
}

/*
   Gaussian Elimination
   After the algorithm has run, A will be in upper-triangular form
   Returns the rank of the matrix
 */
ulong zmod_mat_row_reduce_gauss(zmod_mat_t mat)
{
	ulong p = mat->p;
    if (p < 16)
	{
	   return zmod_mat_row_reduce_gauss_small(mat); 
	}
    ulong i = 0, j = 0, k;
    ulong rows = mat->rows;
	ulong cols = mat->cols;
	ulong coeff;
	
	while ((i < rows) && (j < cols))
	{
		for(k = i; k < rows; k++)
		   if (zmod_mat_get_coeff_ui(mat, k, j)) break;
		
		if (k < rows)
		{
			if (k != i) zmod_mat_swap_rows(mat, i, k);
			ulong n = zmod_mat_get_coeff_ui(mat, i, j);
			ulong n_inv = z_invert(n, p);
			zmod_mat_row_scalar_mul_right(mat, i, n_inv, j);
			for(ulong u = i + 1; u < rows; u++)
			{
			   coeff = zmod_mat_get_coeff_ui(mat, u, j);
			   if (coeff) zmod_mat_row_scalar_submul_right(mat, u, i, coeff, j);
			}
			i++;
		}
		j++;
	}
	return i;
}

/*
   Gauss-Jordan Elimination
   After the algorithm has run, A will be in upper-triangular form in reduced row echelon format
   Returns the rank of the matrix
 */
ulong zmod_mat_row_reduce_gauss_jordan(zmod_mat_t mat)
{
	ulong i = 0, j = 0, k;
    ulong rows = mat->rows;
	ulong cols = mat->cols;
	ulong coeff;
    ulong p = mat->p;

	while ((i < rows) && (j < cols))
	{
		for(k = i; k < rows; k++)
		   if (zmod_mat_get_coeff_ui(mat, k, j)) break;
		
		if (k < rows)
		{
			if (k != i) zmod_mat_swap_rows(mat, i, k);
			ulong n = zmod_mat_get_coeff_ui(mat, i, j);
			ulong n_inv = z_invert(n, p);
			zmod_mat_row_scalar_mul_right(mat, i, n_inv, j);
			for(ulong u = 0; u < i; u++)
			{
			   coeff = zmod_mat_get_coeff_ui(mat, u, j);
			   if (coeff) zmod_mat_row_scalar_submul_right(mat, u, i, coeff, j);
			}
			for(ulong u = i + 1; u < rows; u++)
			{
			   coeff = zmod_mat_get_coeff_ui(mat, u, j);
			   if (coeff) zmod_mat_row_scalar_submul_right(mat, u, i, coeff, j);
			}
			i++;
		}
		j++;
	}
	return i;
}

/*******************************************************************************************

   Addition/subtraction

*******************************************************************************************/

void zmod_mat_add(zmod_mat_t C, zmod_mat_t A, zmod_mat_t B)
{
   ulong p = A->p;
   
   ulong r = A->rows;
   ulong c = A->cols;
   
   for (ulong i = 0; i < r; i++)
      for (ulong j = 0; j < c; j++)
         C->arr[i][j] = z_addmod(A->arr[i][j], B->arr[i][j], p);
}

void zmod_mat_sub(zmod_mat_t C, zmod_mat_t A, zmod_mat_t B)
{
   ulong p = A->p;
   
   ulong r = A->rows;
   ulong c = A->cols;
   
   for (ulong i = 0; i < r; i++)
      for (ulong j = 0; j < c; j++)
         C->arr[i][j] = z_submod(A->arr[i][j], B->arr[i][j], p);
}

/*******************************************************************************************

   Multiplication

*******************************************************************************************/

ulong zmod_mat_scalar_mul(ulong * r, ulong ** arr, ulong c, ulong n, ulong p, double p_inv)
{
   ulong res = 0;

#if FLINT_BITS == 64
   ulong bits = FLINT_BIT_COUNT(p);

   if (bits > FLINT_D_BITS)
   {
      for (ulong i = 0; i < n; i++)
      {
         res = z_addmod(res, z_mulmod2_precomp(r[i], arr[i][c], p, p_inv), p);
      }
   } else
   {
#endif
      for (ulong i = 0; i < n; i++)
      {
         res = z_addmod(res, z_mulmod_precomp(r[i], arr[i][c], p, p_inv), p);
      }
#if FLINT_BITS == 64
   }
#endif

   return res;
}

void zmod_mat_mul_classical(zmod_mat_t prod, zmod_mat_t A, zmod_mat_t B)
{
   ulong p = A->p;
   double p_inv = A->p_inv;

   ulong r1 = A->rows;
   ulong c1 = A->cols;
   ulong r2 = B->rows;
   ulong c2 = B->cols;

   for (ulong i = 0; i < r1; i++)
      for (ulong j = 0; j < c2; j++)
         prod->arr[i][j] = zmod_mat_scalar_mul(A->arr[i], B->arr, j, c1, p, p_inv);
}

void zmod_mat_addmul_classical(zmod_mat_t prod, zmod_mat_t A, zmod_mat_t B)
{
   ulong p = A->p;
   double p_inv = A->p_inv;

   ulong r1 = A->rows;
   ulong c1 = A->cols;
   ulong r2 = B->rows;
   ulong c2 = B->cols;

   for (ulong i = 0; i < r1; i++)
      for (ulong j = 0; j < c2; j++)
         prod->arr[i][j] = z_addmod(prod->arr[i][j], zmod_mat_scalar_mul(A->arr[i], B->arr, j, c1, p, p_inv), p);
}

#define CUTOFF 8

void zmod_mat_mul_strassen(zmod_mat_t C, zmod_mat_t A, zmod_mat_t B)
{
   ulong a, b, c;
   ulong anr, anc, bnr, bnc;

   a = A->rows;
   b = A->cols;
   c = B->cols;

   if (a <= CUTOFF || b <= CUTOFF || c <= CUTOFF)
   {
      zmod_mat_mul_classical(C, A, B);
      return;
   }

   anr = a/2;
   anc = b/2;
   bnr = anc;
   bnc = c/2;

   zmod_mat_t A11; zmod_mat_window_init(A11, A, 0, 0, anr, anc);
   zmod_mat_t A12; zmod_mat_window_init(A12, A, 0, anc, anr, 2*anc);
   zmod_mat_t A21; zmod_mat_window_init(A21, A, anr, 0, 2*anr, anc);
   zmod_mat_t A22; zmod_mat_window_init(A22, A, anr, anc, 2*anr, 2*anc);

   zmod_mat_t B11; zmod_mat_window_init(B11, B, 0, 0, bnr, bnc);
   zmod_mat_t B12; zmod_mat_window_init(B12, B, 0, bnc, bnr, 2*bnc);
   zmod_mat_t B21; zmod_mat_window_init(B21, B, bnr, 0, 2*bnr, bnc);
   zmod_mat_t B22; zmod_mat_window_init(B22, B, bnr, bnc, 2*bnr, 2*bnc);

   zmod_mat_t C11; zmod_mat_window_init(C11, C, 0, 0, anr, bnc);
   zmod_mat_t C12; zmod_mat_window_init(C12, C, 0, bnc, anr, 2*bnc);
   zmod_mat_t C21; zmod_mat_window_init(C21, C, anr, 0, 2*anr, bnc);
   zmod_mat_t C22; zmod_mat_window_init(C22, C, anr, bnc, 2*anr, 2*bnc);

   zmod_mat_t X1; zmod_mat_init(X1, A->p, anr, FLINT_MAX(bnc, anc));
   zmod_mat_t X2; zmod_mat_init(X2, A->p, anc, bnc);

   X1->cols = anc;

   /**
   * \note See Jean-Guillaume Dumas, Clement Pernet, Wei Zhou; "Memory
   * efficient scheduling of Strassen-Winograd's matrix multiplication
   * algorithm"; http://arxiv.org/pdf/0707.2347v3 for reference on the
   * used operation scheduling.
   */
  
   zmod_mat_sub(X1, A11, A21);
   zmod_mat_sub(X2, B22, B12);
   zmod_mat_mul_strassen(C21, X1, X2);
   
   zmod_mat_add(X1, A21, A22);
   zmod_mat_sub(X2, B12, B11);
   zmod_mat_mul_strassen(C22, X1, X2);
   
   zmod_mat_sub(X1, X1, A11);
   zmod_mat_sub(X2, B22, X2);
   zmod_mat_mul_strassen(C12, X1, X2);

   zmod_mat_sub(X1, A12, X1);
   zmod_mat_mul_strassen(C11, X1, B22);

   X1->cols = bnc;
   zmod_mat_mul_strassen(X1, A11, B11);

   zmod_mat_add(C12, X1, C12);
   zmod_mat_add(C21, C12, C21);
   zmod_mat_add(C12, C12, C22);
   zmod_mat_add(C22, C21, C22);
   zmod_mat_add(C12, C12, C11);
   zmod_mat_sub(X2, X2, B21);
   zmod_mat_mul_strassen(C11, A22, X2);

   zmod_mat_clear(X2);

   zmod_mat_sub(C21, C21, C11);
   zmod_mat_mul_strassen(C11, A12, B21);

   zmod_mat_add(C11, X1, C11);

   zmod_mat_clear(X1);

   zmod_mat_window_clear(A11);
   zmod_mat_window_clear(A12);
   zmod_mat_window_clear(A21);
   zmod_mat_window_clear(A22);

   zmod_mat_window_clear(B11);
   zmod_mat_window_clear(B12);
   zmod_mat_window_clear(B21);
   zmod_mat_window_clear(B22);

   zmod_mat_window_clear(C11);
   zmod_mat_window_clear(C12);
   zmod_mat_window_clear(C21);
   zmod_mat_window_clear(C22);

   if (c > 2*bnc) // A by last col of B -> last col of C
   {
      zmod_mat_t Bc; zmod_mat_window_init(Bc, B, 0, 2*bnc, b, c);
      zmod_mat_t Cc; zmod_mat_window_init(Cc, C, 0, 2*bnc, a, c);
      zmod_mat_mul_classical(Cc, A, Bc);
      zmod_mat_window_clear(Bc);
      zmod_mat_window_clear(Cc);
   }

   if (a > 2*anr) // last row of A by B -> last row of C
   {
      zmod_mat_t Ar; zmod_mat_window_init(Ar, A, 2*anr, 0, a, b);
      zmod_mat_t Cr; zmod_mat_window_init(Cr, C, 2*anr, 0, a, c);
      zmod_mat_mul_classical(Cr, Ar, B);
      zmod_mat_window_clear(Ar);
      zmod_mat_window_clear(Cr);
   }

   if (b > 2*anc) // last col of A by last row of B -> C
   {
      zmod_mat_t Ac; zmod_mat_window_init(Ac, A, 0, 2*anc, 2*anr, b);
      zmod_mat_t Br; zmod_mat_window_init(Br, B, 2*bnr, 0, b, 2*bnc);
      zmod_mat_t Cb; zmod_mat_window_init(Cb, C, 0, 0, 2*anr, 2*bnc);
      zmod_mat_addmul_classical(Cb, Ac, Br);
      zmod_mat_window_clear(Ac);
      zmod_mat_window_clear(Br);
      zmod_mat_window_clear(Cb);
   }
}
