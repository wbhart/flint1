#ifndef ZSSMUL_H
#define ZSSMUL_H
void Z_SSMul(mp_limb_t* res, mp_limb_t* data1, mp_limb_t* data2,
              unsigned long limbs, unsigned long limbs2);//, unsigned long tweak);
void Z_fast_mul(mpz_t, mpz_t, mpz_t);//, unsigned long tweak);
#endif
