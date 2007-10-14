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
