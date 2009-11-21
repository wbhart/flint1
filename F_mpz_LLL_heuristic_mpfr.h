/*
   Copyright 2005, 2006 Damien Stehl�.
   Copyright 2008 William Hart, 
   Copyright 2009 William Hart, Andy Novocin

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
   by Phong Nguyen and Damien Stehl�, in the Proceedings of Eurocrypt'2005, 
   Springer-Verlag; and was partly inspired by Shoup's NTL library: 
   http://www.shoup.net/ntl/

*/

/****************************************************************************

   F_mpz_LLL_heuristic_mpfr.h: Lattice reduction of multiprecision integer matrices using mpfrs
	                    for storing approximate Gram Schmidt Orthogonalisations and full precision scalar products occasionally.

*****************************************************************************/

#ifndef FLINT_FMPZ_LLL_FAST_D_H
#define FLINT_FMPZ_LLL_FAST_D_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include <mpfr.h>
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

void Babai_heuristic(int kappa, F_mpz_mat_t B, mpfr_t **mu, mpfr_t **r, mpfr_t *s, 
       mpfr_t **appB, mpfr_t **appSP, int a, int zeros, int kappamax, 
       int n, mpfr_t tmp, mpfr_t rtmp);
                         
void LLL_heuristic(F_mpz_mat_t B);
       
#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
