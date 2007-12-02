/*============================================================================

    Copyright (C) 2007, William Hart

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
#ifndef FLINT_F_MPN_MUL_TUNING_H
#define FLINT_F_MPN_MUL_TUNING_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#define FLINT_FFT_LIMBS_CROSSOVER 2300

#define MUL_TWK_SMALL_CUTOFF 2000
#define MUL_TWK_SMALL_DEFAULT 64
#define MUL_TWK_LARGE_CUTOFF 8350000
#define MUL_TWK_LARGE_DEFAULT 1
#define MUL_TWK_COUNT 20
   
#define SQR_TWK_SMALL_CUTOFF 1564
#define SQR_TWK_SMALL_DEFAULT 16
#define SQR_TWK_LARGE_CUTOFF 8350000
#define SQR_TWK_LARGE_DEFAULT 1
#define SQR_TWK_COUNT 13

#ifdef __cplusplus
 }
#endif
 
#endif
