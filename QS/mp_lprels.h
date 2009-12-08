/*============================================================================
    Copyright 2006 William Hart    

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

#ifndef LPRELS_H
#define LPRELS_H

#include "block_lanczos.h"
#include "common.h"
#include "mp_poly.h"
#include "mp_linear_algebra.h"

#define MPQS_STRING_LENGTH (4 * 1024UL)

typedef struct {
  long q;
  char Y[MPQS_STRING_LENGTH];
  char E[MPQS_STRING_LENGTH];
} mpqs_lp_entry;

char * get_filename(char *dir, char *s);

int mpqs_relations_cmp(const void *a, const void *b);

void flint_fputs(char *s, FILE *file);

void flint_remove(char * name);

long sort_lp_file(char *filename);

long append_file(FILE *fp, FILE *fp1);

long mpqs_mergesort_lp_file_internal(FILE *LPREL, FILE *LPNEW, FILE *COMB, FILE *TMP);

long mergesort_lp_file(char *REL_str, char *NEW_str, char *TMP_str, FILE *COMB);

void add_factor(char **last, unsigned long ei, unsigned long pi);

void add_0(char **last);

void set_exponents(unsigned long *ei, char *r);

void set_lp_entry(mpqs_lp_entry *e, char *buf);

unsigned long combine_large_primes(QS_t * qs_inf, linalg_t * la_inf, poly_t * poly_inf, 
                                                       FILE * COMB, mpz_t factor);
                                                       
FILE * flint_fopen(char * name, char * mode);

char * unique_filename(char *s);

#endif
