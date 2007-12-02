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
/*
   ZmodF_mul-tuning.h
   
   Tuning values for ZmodF_mul.c
   
   (C) 2007 David Harvey

NOTE: the tuning values in this file are for sage.math only.
TODO: write an automatic tuning utility!!

*/

#ifndef FLINT_ZMODF_MUL_TUNING_H
#define FLINT_ZMODF_MUL_TUNING_H

#ifdef __cplusplus
 extern "C" {
#endif
 
// for ZmodF_mul_fft_table[] and ZmodF_sqr_fft_table[],
// first value is crossover n from depth 3 to depth 4,
// then crossover from depth 4 to depth 5, etc.

extern unsigned long ZmodF_mul_plain_threeway_threshold;
extern unsigned long ZmodF_mul_plain_fft_threshold;
extern unsigned long ZmodF_mul_threeway_fft_threshold;
extern unsigned long ZmodF_mul_fft_table[];

extern unsigned long ZmodF_sqr_plain_threeway_threshold;
extern unsigned long ZmodF_sqr_plain_fft_threshold;
extern unsigned long ZmodF_sqr_threeway_fft_threshold;
extern unsigned long ZmodF_sqr_fft_table[];

#ifdef __cplusplus
 }
#endif

#endif
// end of file ****************************************************************
