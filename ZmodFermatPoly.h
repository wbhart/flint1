/****************************************************************************

ZmodFermatPoly_mpn.h: Polynomials over Polynomials over Z mod a Fermat 
                      number p = 2^n+1 where n = r*B*2^l where B is the 
                      number of bits per limb

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/
#ifndef FLINT_ZMODFPOLY_H
#define FLINT_ZMODFPOLY_H

#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "flint-manager.h"
#include "mpn_extras.h"
#include "Zpoly_mpn.h"

/****************************************************************************

   ZmodFermatPoly_t
   -----------

/*
   Polynomials modulo p = 2^n+1 where n = r*B*2^l are stored as arrays of
   pointers into an array of limbs. Each coefficient has a leading limb
   which allows one to perform a certain number of operations on the 
   coefficient before having to reduce mod p.
*/
 
typedef mp_limb_t * f_coeff_t;  // Each coefficient is an array of limbs

typedef struct
{
   f_coeff_t * coeffs;
   unsigned long r;
   unsigned long l;
   unsigned long n; // n = r*B*2^l
   unsigned long length;
} ZmodFPoly_struct;

// Zpoly_mpn_t allows reference-like semantics for Zpoly_mpn_struct:
typedef ZmodFPoly_struct ZmodFPoly_t[1];

/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFPoly_init(ZmodFPoly_t poly, unsigned long r, unsigned long l, 
                                                         unsigned long length);
                                              
void ZmodFPoly_clear(ZmodFPoly_t poly);

/****************************************************************************

   Conversion Routines
   
****************************************************************************/

/* 
   Converts Zpoly_mpn_t "poly_mpn" to a ZmodFPoly. 
   
   Each coefficient of poly_mpn is assumed to fit into a coefficient 
   of poly_f. 
*/

void ZmodFPoly_convert_in_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn);

/* 
   Converts ZmodFPoly "poly_f" to a Zpoly_mpn_t. 
   
   Each coefficient of poly_f is assumed to fit into a coefficient 
   of poly_mpn. 
*/

void ZmodFPoly_convert_out_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn);

/*
   Packs poly_mpn down to the bit into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it. Each of the original
   coefficients is packed into a bitfield "bits" bits wide. 
   
   "bits" is assumed to be less than FLINT_BITS_PER_LIMB.
*/ 
   
void ZmodFPoly_bit_pack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long bits);

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
     
/*
   Packs poly_mpn down to the byte into poly_f. Each coefficient of poly_f
   will have "bundle" coefficients packed into it, each packed into a field
   "bytes" bytes wide.
   
   "bytes" is assumed to be at least FLINT_BITS_PER_LIMB/8, i.e. the coefficients
   are assumed to be at least a limb wide.
*/ 
   
void ZmodFPoly_byte_pack_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long bytes);
     
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

     
/*
   Splits each coefficient of poly_mpn into pieces "limbs" limbs long and 
   stores each piece into bundle coefficients of poly_f. 
*/
 
void ZmodFPoly_split_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long limbs);
     
/*
   Combines each "bundle" coefficients of poly_f, each taken to be "limbs" 
   limbs long, into a coefficient of poly_mpn. 
   
   This function is used for testing purposed only, and is the exact inverse
   of ZmodFPoly_split_mpn.
   
   The number of coefficients extracted is given by the length of poly_mpn.
*/
 
void ZmodFPoly_unsplit_mpn(ZmodFPoly_t poly_f, Zpoly_mpn_t poly_mpn,
     unsigned long bundle, unsigned long limbs);

#endif
