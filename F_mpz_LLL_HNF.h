/*============================================================================

    F_mpz_mat.h: Header file for F_mpz_LLL_HNF.c

    Copyright (C) 2009, William Hart 

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

==============================================================================*/

#ifndef FLINT_F_MPZ_LLL_HNF_H
#define FLINT_F_MPZ_LLL_HNF_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint.h"
#include "F_mpz_mat.h"
#include "d_mat.h"

void F_mpz_mat_HNF_LLL(F_mpz_mat_t A, F_mpz_mat_t B, F_mpz_mat_t G);

#ifdef __cplusplus
 }
#endif
 
#endif

 // *************** end of file
