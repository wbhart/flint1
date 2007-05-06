/****************************************************************************

ZmodFpoly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include "flint.h"
#include "flint-manager.h"
#include "ZmodFpoly.h"


/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFPoly_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;
   
   poly->storage = (mp_limb_t*) flint_malloc_limbs(bufs * (n+1));

   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) flint_malloc_bytes(bufs * sizeof(ZmodF_t*));
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   for (unsigned long i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}


void ZmodFPoly_clear(ZmodFpoly_t poly)
{
   flint_free(poly->coeffs);
   flint_free(poly->storage);
}



/****************************************************************************

   Conversion Routines
   
****************************************************************************/

void ZmodFPoly_convert_in_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   abort();
}

void ZmodFPoly_convert_out_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   abort();
}

void ZmodFPoly_bit_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                            unsigned long bundle, unsigned long bits)
{
   abort();
}

void ZmodFPoly_bit_unpack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                              unsigned long bundle, unsigned long bits)
{
   abort();
}
     
void ZmodFPoly_byte_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                             unsigned long bundle, unsigned long bytes)
{
   abort();
}
     
void ZmodFPoly_byte_unpack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                               unsigned long bundle, unsigned long bytes)
{
   abort();
}

     
void ZmodFPoly_split_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                         unsigned long bundle, unsigned long limbs)
{
   abort();
}
     
void ZmodFPoly_unsplit_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                           unsigned long bundle, unsigned long limbs)
{
   abort();
}
     


// end of file ****************************************************************
