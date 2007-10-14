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
