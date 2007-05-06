/****************************************************************************

ZmodFermatPoly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "flint-manager.h"
#include "ZmodFermatPoly.h"


/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFPoly_init(ZmodFPoly_t poly, unsigned long m, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->m = m;
   poly->scratch_count = scratch_count;
   poly->M = 1 << m;
   poly->length = 0;
   
   unsigned long bufs = poly->M + scratch_count;
   
   poly->storage = (mp_limb_t*) flint_malloc_limbs(bufs * (n+1));

   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) flint_malloc_bytes(bufs * sizeof(ZmodF_t*));
   poly->scratch = poly->coeffs + poly->M;
   
   poly->coeffs[0] = poly->storage;
   for (unsigned long i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}


void ZmodFPoly_clear(ZmodFPoly_t poly)
{
   flint_free(poly->coeffs);
   flint_free(poly->storage);
}



/****************************************************************************

   Conversion Routines
   
****************************************************************************/

void ZmodFPoly_convert_in_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   abort();
}

void ZmodFPoly_convert_out_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   abort();
}

void ZmodFPoly_bit_pack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
                            unsigned long bundle, unsigned long bits)
{
   abort();
}

void ZmodFPoly_bit_unpack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
                              unsigned long bundle, unsigned long bits)
{
   abort();
}
     
void ZmodFPoly_byte_pack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
                             unsigned long bundle, unsigned long bytes)
{
   abort();
}
     
void ZmodFPoly_byte_unpack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
                               unsigned long bundle, unsigned long bytes)
{
   abort();
}

     
void ZmodFPoly_split_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
                         unsigned long bundle, unsigned long limbs)
{
   abort();
}
     
void ZmodFPoly_unsplit_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
                           unsigned long bundle, unsigned long limbs)
{
   abort();
}
     


// end of file ****************************************************************
