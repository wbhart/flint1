/*
   Copyright 2005, 2006 Damien Stehlé.
   Copyright 2009, 2010 William Hart, Andy Novocin

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

   F_mpz_LLL_helper.h

*****************************************************************************/

#ifndef FLINT_FMPZ_LLL_HELPER_H
#define FLINT_FMPZ_LLL_HELPER_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <gmp.h>
#include "flint.h"
#include "F_mpz_mat.h"

double heuristic_scalar_product(double * vec1, double * vec2, ulong n, 
								F_mpz_mat_t B, ulong k, ulong j, long exp_adj);


#ifdef __cplusplus
 }
#endif
 
#endif

// *************** end of file
