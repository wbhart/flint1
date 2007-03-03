#ifndef ZVEC_H
#define ZVEC_H

typedef struct Zvec_t
{
   mpz_t * coords;
   unsigned long length;
} Zvec_t;

typedef Zvec_t Zvec[1];

void initZvec(void);
void Zvec_init1(Zvec, unsigned long, unsigned long);
void Zvec_init3(Zvec, unsigned long, unsigned long);
void Zvec_alloc(Zvec, unsigned long, unsigned long, int);
void Zvec_release(Zvec);
void Zvec_mul(Zvec, Zvec, Zvec);
void Zvec_mul_naieve(Zvec, Zvec, Zvec);
void Zvec_karamul(Zvec, Zvec, Zvec, unsigned long);
void Zvec_radixMul(Zvec, Zvec, Zvec, unsigned long);
unsigned long Zvec_max_bits(Zvec);
unsigned long Zvec_max_length(Zvec, Zvec);
void Zvec_clear(Zvec);
void Zvec_SSMul(Zvec, Zvec, Zvec, unsigned long, int);
void Zvec_mul(Zvec res, Zvec a, Zvec b);

#endif
