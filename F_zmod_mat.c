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

   zmod_mat.c: Matrices over (unsigned) long mod p, for p prime with packed
	            representation (using packed_vec) (FLINT 2.0).

   Copyright (C) 2008, William Hart.
   Copyright (C) 2008, Richard Howell-Peak
   
*****************************************************************************/

#include "F_zmod_mat.h"
#include "zmod_poly.h"
#include "long_extras.h"
#include "flint.h"
#include "packed_vec.h"

/****************************************************************************

   Initialisation and memory management

****************************************************************************/

void F_zmod_mat_init(F_zmod_mat_t mat, ulong p, ulong rows, ulong cols)
{
   F_zmod_mat_init_precomp(mat, p, z_precompute_inverse(p), rows, cols);
}

void F_zmod_mat_init_precomp(F_zmod_mat_t mat, ulong p, double p_inv, 
						                             ulong rows, ulong cols)
{
   pv_init(&mat->arr, rows*cols, pv_bit_fit(FLINT_BIT_COUNT(p)));
   mat->rows = (ulong *) flint_heap_alloc(rows);
   
   // Set up the rows
   ulong offset = 0;
	for (ulong i = 0; i < rows; i++, offset += cols)
      mat->rows[i] = offset;

   mat->p = p;
   mat->p_inv = p_inv;
   mat->r = rows;
   mat->c = cols;
}

void F_zmod_mat_clear(F_zmod_mat_t mat)
{
   flint_heap_free(mat->rows);
   pv_clear(&mat->arr);
}

/*******************************************************************************************

   Arithmetic

*******************************************************************************************/

void F_zmod_mat_add(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2)
{
	ulong p = mat1->p;
	
	for (ulong i = 0; i < mat1->r; i++)
	{
		ulong m1, m2;
		pv_iter_s i1, i2, i3;
		PV_ITER_INIT(i1, mat1->arr, mat1->rows[i]); // row i for mat1
      PV_ITER_INIT(i2, mat2->arr, mat2->rows[i]); // row i for mat2
      PV_ITER_INIT(i3, res->arr, res->rows[i]); // row i for res
      for (ulong j = 0; j < mat1->c; j++)
		{
			PV_GET_NEXT(m1, i1);
			PV_GET_NEXT(m2, i2);
			PV_SET_NEXT(i3, z_addmod(m1, m2, p));
		}
	}
}

void F_zmod_mat_sub(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2)
{
	ulong p = mat1->p;
	
	for (ulong i = 0; i < mat1->r; i++)
	{
		ulong m1, m2;
		pv_iter_s i1, i2, i3;
		PV_ITER_INIT(i1, mat1->arr, mat1->rows[i]); // row i for mat1
      PV_ITER_INIT(i2, mat2->arr, mat2->rows[i]); // row i for mat2
      PV_ITER_INIT(i3, res->arr, res->rows[i]); // row i for res
      for (ulong j = 0; j < mat1->c; j++)
		{
			PV_GET_NEXT(m1, i1);
			PV_GET_NEXT(m2, i2);
			PV_SET_NEXT(i3, z_submod(m1, m2, p));
		}
	}
}

void F_zmod_mat_neg(F_zmod_mat_t res, F_zmod_mat_t mat1)
{
	ulong p = mat1->p;
	
	for (ulong i = 0; i < mat1->r; i++)
	{
		ulong m1;
		pv_iter_s i1, i2;
		PV_ITER_INIT(i1, mat1->arr, mat1->rows[i]); // row i for mat1
      PV_ITER_INIT(i2, res->arr, res->rows[i]); // row i for res
      for (ulong j = 0; j < mat1->c; j++)
		{
			PV_GET_NEXT(m1, i1);
			PV_SET_NEXT(i2, z_negmod(m1, p));
		}
	}
}

void F_zmod_mat_mul_classical(F_zmod_mat_t res, F_zmod_mat_t mat1, F_zmod_mat_t mat2)
{
   ulong c1 = mat1->c;
	ulong r2 = mat2->r;
   
	if ((c1 != r2) || (c1 == 0))
	{
		printf("FLINT exception : invalid matrix multiplication!\n");
		abort();
	}

	ulong r1 = mat1->r;
   ulong c2 = mat2->c;
	ulong p = mat1->p;
   double pinv = mat1->p_inv;

	mp_limb_t * temp = (mp_limb_t *) flint_stack_alloc(2*c2);
   mp_limb_t * temp2 = (mp_limb_t *) flint_stack_alloc(2*c2);

   ulong spare_bits = 2*FLINT_BITS - 2*(FLINT_BIT_COUNT(p));
	if (spare_bits >= FLINT_BITS) spare_bits = FLINT_BITS - 1; 
	ulong red_max = (1L << spare_bits); // number of loops before reduction must be done

	for (ulong i = 0; i < r1; i++) // for each row of mat1
	{
		pv_iter_s iter1;
		PV_ITER_INIT(iter1, mat1->arr, mat1->rows[i]); // iterate along row i of mat1
		
		ulong c;
		PV_GET_NEXT(c, iter1); // get first coefficient of row i of mat1

		pv_iter_s iter2;
		PV_ITER_INIT(iter2, mat2->arr, mat2->rows[0]); 
		ulong j = 0, d;
		for (ulong k = 0; k < 2*c2; k+=2) // do initial scalar product of row 1 of mat2 by c
		{
         PV_GET_NEXT(d, iter2);
			umul_ppmm(temp[k+1], temp[k], d, c);
		}

		for (j = 1; j < c1; j+= (red_max - 1)) // add up to red_max - 1 
			                                         // scalar products at a time
		{
         for (ulong s = 0; s < FLINT_MIN(red_max - 1, c1 - j); s++) // don't exceed c1 rows
			{
				PV_GET_NEXT(c, iter1); // get next coefficient of row i of mat1

			   PV_ITER_INIT(iter2, mat2->arr, mat2->rows[j+s]); // iterate along row j+s of mat2
            for (ulong k = 0; k < 2*c2; k+=2) // do scalar product of row j+s of mat2 by c
		      {
               PV_GET_NEXT(d, iter2);
			      umul_ppmm(temp2[k+1], temp2[k], d, c);
				   add_ssaaaa(temp[k+1], temp[k], temp[k+1], temp[k], temp2[k+1], temp2[k]); //add
				                                                              // to existing sum
		      }
			}

			for (ulong k = 0; k < 2*c2; k+=2) // do a reduction
			{
				temp[k] = z_ll_mod_precomp(temp[k+1], temp[k], p, pinv); // reduce mod p
				temp[k+1] = 0L;
			}
		}

		pv_iter_s iter3;
		PV_ITER_INIT(iter3, res->arr, res->rows[i]); // iterate along row i of res
      for (ulong k = 0; k < 2*c2; k+=2) // store row i of res
			PV_SET_NEXT(iter3, temp[k]);
	}
		
	flint_stack_release(); // temp2
   flint_stack_release(); // temp
}

/*******************************************************************************************

   Conversions

*******************************************************************************************/

/*
   Set a row to the coefficients of a polynomial, starting with the constant coefficient
   Assumes that poly->length <= mat->cols
*/
/*void zmod_poly_to_zmod_mat_row(zmod_mat_t mat, ulong row, zmod_poly_t poly)
{
   ulong * r1 = mat->arr[row];
   ulong * coeffs = poly->coeffs;
   ulong cols = mat->cols;

   long i;
   
   for (i = 0; i < poly->length; i++)
      r1[i] = coeffs[i];
   
   for ( ; i < cols; i++)
      r1[i] = 0L;
}*/

/*
   Set a column to the coefficients of a polynomial, starting with the constant coefficient
   Assumes that poly->length <= mat->rows
*/

/*void zmod_poly_to_zmod_mat_col(zmod_mat_t mat, ulong col, zmod_poly_t poly)
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
}*/

/*
   Set a zmod_poly's coefficients to the entries in a column, starting with the constant coefficient
*/

/*void zmod_mat_col_to_zmod_poly(zmod_poly_t poly, zmod_mat_t mat, ulong col)
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
}*/

/*
   Set a zmod_poly's coefficients to the entries in a column, starting with the constant coefficient
   but shifting along by one for every non-zero entry in shift
*/

/*void zmod_mat_col_to_zmod_poly_shifted(zmod_poly_t poly, zmod_mat_t mat, ulong col, ulong * shift)
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
}*/

/*******************************************************************************************

   Elementary row operations

*******************************************************************************************/

/*
   row1 = row1 + u * row2
   only columns [start..mat->cols) are affected
   assumes u is reduced mod mat->p
*/

/*void zmod_mat_row_scalar_addmul_right(zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start)
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
}*/

/*
   row1 = row1 - u * row2
   only columns [start..mat->cols) are affected
   assumes u is reduced mod mat->p
*/

/*void zmod_mat_row_scalar_submul_right(zmod_mat_t mat, ulong row1, ulong row2, ulong u, ulong start)
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
}*/

/*
   row = row * u
   only columns [start..mat->cols) are affected
   assumes u is reduced mod mat->p
*/

/*void zmod_mat_row_scalar_mul_right(zmod_mat_t mat, ulong row, ulong u, ulong start)
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
}*/

/*
   r1 = r2 * u
   only columns [start..end) are affected
   assumes u is reduced mod p
*/

/*void zmod_vec_scalar_mul_range(ulong * r1, ulong * r2, ulong u, ulong p, double p_inv, ulong start, ulong end)
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
}*/

/*
   r1 = r1 - r2
   only columns [start..end) are affected
*/

/*void zmod_vec_sub_range(ulong * r1, ulong * r2, ulong p, ulong start, ulong end)
{
   for (ulong i = start; i < end; i++)
   {
      r1[i] = z_submod(r1[i], r2[i], p);
   }
}*/

/*******************************************************************************************

   Input/output

*******************************************************************************************/

/*void zmod_mat_print(zmod_mat_t mat)
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
}*/

/*
   Gaussian Elimination
   After the algorithm has run, A will be in upper-triangular form
   Returns the rank of the matrix
   Assumes the modulus is a very small prime
 */
/*ulong zmod_mat_row_reduce_gauss_small(zmod_mat_t mat)
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
}*/

/*
   Gaussian Elimination
   After the algorithm has run, A will be in upper-triangular form
   Returns the rank of the matrix
 */

/*ulong zmod_mat_row_reduce_gauss(zmod_mat_t mat)
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
}*/

/*
   Gauss-Jordan Elimination
   After the algorithm has run, A will be in upper-triangular form in reduced row echelon format
   Returns the rank of the matrix
 */

/*ulong zmod_mat_row_reduce_gauss_jordan(zmod_mat_t mat)
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
}*/

