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
/****************************************************************************

   theta.h: theta functions

   Copyright (C) 2008, Gonzalo Tornaria and William Hart

*****************************************************************************/

#ifndef THETA_H
#define THETA_H

#ifdef __cplusplus
 extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"

void theta_print(long *coeff, ulong start, ulong len);

void theta_fprint(long *coeff, ulong start, ulong len, FILE * f);

void theta(long * out, ulong start, ulong len);

void theta_1d(unsigned long a, unsigned long b, unsigned long c,
				                                          long * out, ulong start, ulong len);

void theta_1d_0(unsigned long a, unsigned long b, unsigned long c,
				                                          long * out, ulong start, ulong len);

void theta_1d_quadchar(long * character, unsigned long a, unsigned long b, unsigned long c,
				                                          long * out, ulong start, ulong len);

void theta_1d_quadchar_0(long * character, unsigned long a, unsigned long b, unsigned long c,
				                        unsigned long m, long * out, ulong start, ulong len);

void theta_1d_quadchar_2(long * character, unsigned long a, unsigned long b, unsigned long c,
				                        unsigned long m, long * out, ulong start, ulong len);

void theta_1d_B(long * out, ulong start, ulong len);

void theta_2d(long *out, ulong start, ulong len);

void theta_2d_A1(long *out, ulong start, ulong len);

void theta_2d_A2(long *out, ulong start, ulong len);

void theta_2d_B(long *out, ulong start, ulong len);

void theta_2d_C(long *out, ulong start, ulong len);

#ifdef __cplusplus
 }
#endif

#endif

// *************** end of file
