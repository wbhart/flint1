/******************************************************************************

 mpn_extras.c
 
 Extra functions for manipulating mpn's and limbs.

 Copyright (C) 2006, William Hart

 mp_limb_t mpn_divmod_1_preinv was adapted from GMP, (C) Free Software Foundation

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "long_extras.h"

/*=======================================================================================*/

/* 
    Performs division by a limb d and places the quotient in qp and returns the 
    remainder. Requires a single limb approximation to 1/d as input. If the most
    significant bit of d is not 1 it expects d to be shifted left (by norm bits)
    until the most significant bit is 1 before the inverse is computed. However 
    the original d should be supplied to the function, not the shifted d. 
    
    This code has been adapted from code found in the GMP package version 4.2.1
    (divrem_1.c) (C) Free Software Foundation
*/
mp_limb_t mpn_divmod_1_preinv(mp_limb_t * qp, mp_limb_t * up, 
                                  unsigned long un, mp_limb_t d, mp_limb_t dinv, unsigned long norm)
{
  mp_size_t  n;
  mp_size_t  i;
  mp_limb_t  n1, n0;
  mp_limb_t  r = 0;

  n = un;
  if (n == 0)
    return 0;
  
  qp += (n - 1);   /* Make qp point at most significant quotient limb */

  if ((d & (1L<<(FLINT_BITS-1))) != 0)
  {
     if (un != 0)
     {
        /* High quotient limb is 0 or 1, skip a divide step. */
	    mp_limb_t q;
	    r = up[un - 1];
	    q = (r >= d);
	    *qp-- = q;
	    r -= (d & -q);
	    n--;
	    un--;
	 }

     /* Multiply-by-inverse, divisor already normalized. */
     for (i = un - 1; i >= 0; i--)
     {
        n0 = up[i];
        udiv_qrnnd_preinv (*qp, r, r, n0, d, dinv);
        qp--;
     }
     return r;
  } else
  {
     /* Most significant bit of divisor == 0.  */
     
     /* Skip a division if high < divisor (high quotient 0).  Testing here
	 before normalizing will still skip as often as possible.  */
     if (un != 0)
	 {
	    n1 = up[un - 1];
	    if (n1 < d)
        {
           r = n1;
	       *qp-- = 0;
	       n--;
	       if (n == 0) return r;
	       un--;
        }
	 }  

     d <<= norm;
     r <<= norm;

     if (un != 0)
     {
        n1 = up[un - 1];
        r |= (n1 >> (FLINT_BITS - norm));
        for (i = un - 2; i >= 0; i--)
		{
		  n0 = up[i];
		  udiv_qrnnd_preinv (*qp, r, r, 
				     ((n1 << norm) | (n0 >> (FLINT_BITS - norm))), d, dinv);
		  qp--;
		  n1 = n0;
		}
        udiv_qrnnd_preinv (*qp, r, r, n1 << norm, d, dinv);
        qp--;
     }
     
     return r >> norm;
  }
}

mp_limb_t mpn_addmul(mp_limb_t * rp, mp_limb_t * s1p, unsigned long s1n, 
                                      mp_limb_t * s2p, unsigned long s2n)
{
   if (s2n == 0) return 0;
   
   mp_limb_t carry;
   
   carry = mpn_addmul_1(rp, s1p, s1n, s2p[0]);
   for (unsigned long i = 1; i < s2n; i++)
   {
      carry = mpn_add_1(rp+i+s1n-1, rp+i+s1n-1, 1, carry); 
      if (s2p[i]) carry += mpn_addmul_1(rp+i, s1p, s1n, s2p[i]);
   }
   return carry;
}

