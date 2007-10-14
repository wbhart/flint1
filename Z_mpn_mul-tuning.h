#ifndef FLINT_ZMUL_TUNING_H
#define FLINT_ZMUL_TUNING_H

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
