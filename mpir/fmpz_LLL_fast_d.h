/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   fmpz_mat.h: Matrices over the "flat" multi-precision integer format

   Copyright (C) 2007, William Hart

*****************************************************************************/

#ifndef MPIR_FMPZ_LLL_FAST_D_H
#define MPIR_FMPZ_LLL_FAST_D_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "mpir.h"
#include "fmpz.h"

double halfplus, onedothalfplus, ctt;

#if MPIR_BITS == 32
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

void Babai (int kappa, fmpz_mat_t B, double **mu, double **r, double *s, 
                            double **appB, int *expo, double **appSP, 
                         int a, int zeros, int kappamax, int n, mpz_t ztmp);
                         
void LLL (fmpz_mat_t B);
       
#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
