/****************************************************************************

ZmodFermatPoly.c: Polynomials over Z mod a Fermat number p = 2^n+1 where 
                  n = r*B*2^l where B is the number of bits per limb

Copyright (C) 2007, William Hart and David Harvey

(Development code)

*****************************************************************************/

#include "flint.h"
#include "ZmodFermatPoly.h"
#include "flint-manager.h"

/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFPoly_init(ZmodFPoly_t poly, unsigned long r, unsigned long l, 
                                                         unsigned long length)
{
   unsigned long limbs = r*(1<<l);
   f_coeff_t data = (f_coeff_t) flint_malloc_limbs(length*(limbs+1));
   
   poly->coeffs = (f_coeff_t *) flint_malloc_limbs(length); 
   for (unsigned long i = 0; i < length; i++)
   {
      poly->coeffs[i] = data;
      data += (limbs+1);
   }
   
   poly->r = r;
   poly->l = l;
   poly->n = FLINT_BITS_PER_LIMB*limbs;
}
                                              
void ZmodFPoly_clear(ZmodFPoly_t poly)
{
   flint_free(poly->coeffs[0]);
   flint_free(poly->coeffs);
}

/****************************************************************************

   Conversion Routines
   
****************************************************************************/

/* 
   Converts Zpoly_mpn_t "poly_mpn" to a ZmodFPoly. 
   
   Each coefficient of poly_mpn is assumed to fit into a coefficient 
   of poly_f. 
*/

void ZmodFPoly_convert_in_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn);
{
   abort();
}

/* 
   Converts ZmodFPoly "poly_f" to a Zpoly_mpn_t. 
   
   Each coefficient of poly_f is assumed to fit into a coefficient 
   of poly_mpn. 
*/

void ZmodFPoly_convert_out_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn);
{
   abort();
}

/*
   Packs poly_mpn down to the bit into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it. Each of the original
   coefficients is packed into a bitfield "bits" bits wide. 
   
   "bits" is assumed to be less than FLINT_BITS_PER_LIMB.
*/ 
   
void ZmodFPoly_bit_pack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long bits);
{
   abort();
}

/*
   Unpacks poly_f into poly_mpn. This is the inverse of ZmodFPoly_bitpack_mpn.
   Each coeff of poly_f is assumed to contain "bundle" coefficients, each stored 
   in a bitfield "bits" bits wide. 
   
   The total number of coefficients to be unpacked is given by the length of 
   poly_mpn.
   
   "bits" is assumed to be less than FLINT_BITS_PER_LIMB.
*/ 
   
void ZmodFPoly_bit_unpack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long bits);
{
   abort();
}
     
/*
   Packs poly_mpn down to the byte into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it, each packed into a field
   "bytes" bytes wide.
   
   "bytes" is assumed to be at least FLINT_BITS_PER_LIMB/8, i.e. the coefficients
   are assumed to be at least a limb wide.
*/ 
   
void ZmodFPoly_byte_pack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long bytes);
{
   abort();
}
     
/*
   Unpacks poly_f into poly_mpn. Each coefficient of poly_f will have "bundle" 
   coefficients, each packed into a field "bytes" bytes wide.
   
   The total number of coefficients to be unpacked is given by the length of 
   poly_mpn.
   
   "bytes" is assumed to be at least FLINT_BITS_PER_LIMB/8, i.e. the coefficients
   are assumed to be at least a limb wide.
*/ 
   
void ZmodFPoly_byte_unpack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long bytes);
{
   abort();
}

     
/*
   Splits each coefficient of poly_mpn into pieces "limbs" limbs long and 
   stores each piece into bundle coefficients of poly_f. 
*/
 
void ZmodFPoly_split_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long limbs);
{
   abort();
}
     
/*
   Combines each "bundle" coefficients of poly_f, each taken to be "limbs" 
   limbs long, into a coefficient of poly_mpn. 
   
   This function is used for testing purposed only, and is the exact inverse
   of ZmodFPoly_split_mpn.
   
   The number of coefficients extracted is given by the length of poly_mpn.
*/
 
void ZmodFPoly_unsplit_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long limbs)
{
   abort();
}
     


