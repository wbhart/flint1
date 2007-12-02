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
   mpz_poly_tuning.h
   
   Tuning values for mpz_poly.c
   
   (C) 2007 David Harvey

NOTE: the tuning values in this file are for sage.math only.
TODO: write an automatic tuning utility!!

*/

#ifndef FLINT_MPZ_POLY_TUNING_H
#define FLINT_MPZ_POLY_TUNING_H

#ifdef __cplusplus
 extern "C" {
#endif

/*
   mpz_poly_kara_crossover_table[k] is the smallest length for which
   karatsuba should be used when the coefficients have k+1 limbs.
   The number of entries in the table is mpz_poly_kara_crossover_table_size.
*/
extern unsigned long mpz_poly_kara_crossover_table[];
extern unsigned long mpz_poly_kara_crossover_table_size;

#ifdef __cplusplus
 extern "C" {
#endif

#endif
// end of file ****************************************************************
