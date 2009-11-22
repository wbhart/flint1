/*
   Copyright 2005, 2006 Damien Stehlé.
   Copyright 2008 William Hart

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2 of the License, or (at your
   option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along
   with this program; see the file COPYING.  If not, write to the Free
   Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
   02111-1307, USA.

   This program implements ideas from the paper "Floating-point LLL Revisited", 
   by Phong Nguyen and Damien Stehlé, in the Proceedings of Eurocrypt'2005, 
   Springer-Verlag; and was partly inspired by Shoup's NTL library: 
   http://www.shoup.net/ntl/

*/

/****************************************************************************

   F_mpz_LLL_fast_d.h: Lattice reduction of multiprecision integer matrices using doubles
	                    for storing approximate Gram Schmidt Orthogonalisations

*****************************************************************************/

#ifndef FLINT_FMPZ_LLL_FAST_D_H
#define FLINT_FMPZ_LLL_FAST_D_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "flint.h"
#include "F_mpz_mat.h"

#ifndef NAN
#define NAN (0.0/0.0)
#endif 

double halfplus, onedothalfplus, ctt;

#if FLINT_BITS == 32
#define CPU_SIZE_1 31
#define MAX_LONG 0x1p31
#else
#define CPU_SIZE_1 53
#define MAX_LONG 0x1p53
#endif

#ifndef ETA
#define ETA 0.51
#endif
#ifndef DELTA
#define DELTA 0.99
#endif

ulong getShift(F_mpz_mat_t B);

void Babai (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
                            double **appB, int *expo, double **appSP, 
                         int a, int zeros, int kappamax, int n);

void Babai_heuristic_d_2exp(int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
                            double **appB, int *expo, double **appSP, 
                         int a, int zeros, int kappamax, int n, int *cexpo);

void Babai_heuristic_d (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n);
                         
void LLL (F_mpz_mat_t B);

void LLL_heuristic_d_2exp (F_mpz_mat_t B, int *cexpo);

int LLL_heuristic_d_2exp_with_removal(F_mpz_mat_t B, int *cexpo, F_mpz_t gs_B);

int LLL_heuristic_d_with_removal(F_mpz_mat_t B, F_mpz_t gs_B);
       
#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
