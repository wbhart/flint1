/****************************************************************************

Z_mpn.h: Z arithmetic

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_Z_MPN_H
#define FLINT_Z_MPN_H

#include <gmp.h>

/****************************************************************************/

void Z_mpn_mul(mp_limb_t * res, mp_limb_t * data1, unsigned long limbs1, 
                                      mp_limb_t * data2, unsigned long limbs2);
                                      
void Z_mul(mpz_t res, mpz_t a, mpz_t b);
                                      
#endif
