/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

test-support.c: Support code for test modules

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <string.h>
#include "mpir.h"
#include "test_support.h"

gmp_randstate_t randstate;

/*
   Generate a random integer with up to the given limit. 
   If the limit is 0, a random limb is generated.
*/

unsigned long randint(unsigned long limit) 
{
#if MPIR_BITS == 32
    static uint64_t randval = 4035456057U;
    randval = ((uint64_t)randval*(uint64_t)1025416097U+(uint64_t)286824430U)%(uint64_t)4294967311U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)randval%limit;
#else
    static unsigned long randval = 4035456057U;
    static unsigned long randval2 = 6748392731U;
    randval = ((unsigned long)randval*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    randval2 = ((unsigned long)randval2*(unsigned long)1647637699U+(unsigned long)286824428U)%(unsigned long)4294967357U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)(randval+(randval2<<32))%limit;
#endif
}

/*
   Generate a random integer with up to the given number of bits [0, MPIR_BITS]
*/

unsigned long randbits(unsigned long bits)
{
   return randint(MPIR_LSHIFT(1L, bits));
}

char randbyte() 
{
    static ulong randval = 4035456057U;
    randval = ((ulong)randval*1025416097U+286824428U)%(ulong)4294967291U;
    
    return (char)randval;
}

// end of file ****************************************************************
