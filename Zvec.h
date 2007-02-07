#ifndef ZVEC_H
#define ZVEC_H
//#include "Z.h"

typedef struct Zvec
{
   mpz_t * coords;
   unsigned long length;
} Zvec;

void initZvec(void);
void Zvec_init1(Zvec&, unsigned long, unsigned long);
void Zvec_init3(Zvec&, unsigned long, unsigned long);
void Zvec_alloc(Zvec&, unsigned long, unsigned long, int);
void Zvec_release(Zvec&);
static inline void LeftRotate(mpz_t, mpz_t, unsigned long, mpz_t, unsigned long);
static inline void Rotate(mpz_t, mpz_t, long, mpz_t, unsigned long);
void Zvec_mul(Zvec&, Zvec&, Zvec&);
void Zvec_mul_naieve(Zvec&, Zvec&, Zvec&);
void Zvec_karamul(Zvec&, Zvec&, Zvec&, unsigned long);
void Zvec_radixMul(Zvec&, Zvec&, Zvec&, unsigned long);
unsigned long Zvec_max_bits(Zvec);
unsigned long Zvec_max_length(Zvec, Zvec);
void Zvec_clear(Zvec&);
static inline void Zvec_free2(Zvec&);
static inline void Truncate(mpz_t, mpz_t, unsigned long);
void Zvec_SSMul(Zvec&, Zvec&, Zvec&, unsigned long, int);
void Zvec_mul(Zvec& res, Zvec& a, Zvec& b);

#endif
