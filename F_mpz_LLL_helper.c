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

   F_mpz_LLL_helper.c: Helper functions for F_mpz versions of LLL

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <ctype.h>
#include <limits.h>
#include "gmp.h"
#include "flint.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL_helper.h"
#include "mpfr.h"
#include "mpfr_mat.h"
#include "d_mat.h"

/*
   Computes the scalar product of two vectors of doubles vec1 and vec2, which are 
   respectively double approximations (up to scaling by a power of 2) to rows k and
   j in the exact integer matrix B. If massive cancellation is detected an exact
   computation is made.

   The exact computation is scaled by 2^-exp_adj, where exp_adj = r2 + r1 where
   r2 is the exponent for row j and r1 is the exponent for row k (i.e. row j is 
   notionally thought of as being multiplied by 2^r2, etc).

   The final scalar product computed by this function is then notionally the return
   value times 2^exp_adj.
*/
double heuristic_scalar_product(double * vec1, double * vec2, ulong n, 
								F_mpz_mat_t B, ulong k, ulong j, long exp_adj)
{
  double sum = _d_vec_scalar_product(vec1, vec2, n);
  double tmp = _d_vec_norm(vec1, n);
  double tmp2 = _d_vec_norm(vec2, n);

  tmp = ldexp(tmp*tmp2, -70);
  tmp2 = sum*sum;

  if (tmp2 <= tmp)
  {
     F_mpz_t sp;
     F_mpz_init(sp);
     ulong exp;
     _F_mpz_vec_scalar_product(sp, B->rows[k], B->rows[j], n);
     sum = F_mpz_get_d_2exp(&exp, sp);
     sum = ldexp(sum, exp - exp_adj);
     F_mpz_clear(sp);
  }

  return sum;
} 



