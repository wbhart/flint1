/*
   ZmodF_mul-tuning.h
   
   Tuning values for ZmodF_mul.c
   
   (C) 2007 David Harvey

NOTE: the tuning values in this file are for sage.math only.
TODO: write an automatic tuning utility!!

*/

#ifndef FLINT_ZMODF_MUL_TUNING_H
#define FLINT_ZMODF_MUL_TUNING_H


// ===========================================
//   multiplication


// crossover from plain to threeway algorithm, assuming 3 divides n
#define ZMODF_MUL_PLAIN_THREEWAY_THRESHOLD 36

// crossover from threeway to negacyclic algorithm, assuming 3 divides n
#define ZMODF_MUL_THREEWAY_NEGACYCLIC_THRESHOLD 750

// crossover from plain to negacyclic algorithm, if 3 does not divide n
#define ZMODF_MUL_PLAIN_NEGACYCLIC_THRESHOLD 320


// the minimum depth used for negacyclic transforms
#define ZMODF_MUL_MIN_NEGACYCLIC_DEPTH 4

// the i-th entry in this table is the value of n for which we switch from
// depth K+i to depth K+i+1, where K = ZMODF_MUL_MIN_NEGACYCLIC_DEPTH
// (the first three values are tuned for sage.math, I made the rest up)
unsigned long ZMODF_MUL_NEGACYCLIC_THRESHOLD[] =
   {500, 1400, 3000, 6000, 12000, 24000, 48000, 0};


// ===========================================
//   squaring


// same as above, but for squaring
// todo: for now these are the same as the multiplication values
#define ZMODF_SQR_PLAIN_THREEWAY_THRESHOLD 36
#define ZMODF_SQR_THREEWAY_NEGACYCLIC_THRESHOLD 750
#define ZMODF_SQR_PLAIN_NEGACYCLIC_THRESHOLD 320
#define ZMODF_SQR_MIN_NEGACYCLIC_DEPTH 4
unsigned long ZMODF_SQR_NEGACYCLIC_THRESHOLD[] =
   {500, 1400, 3000, 6000, 12000, 24000, 48000, 0};


#endif
