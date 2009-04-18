/*
   array-profile.c:  routines for profiling simple array operations
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.9).
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) version 3 of the License.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/


#include <math.h>
#include "support.h"
#include "profiler.h"
#include "zn_poly_internal.h"
#include "zn_poly.h"


typedef void  (*bfly_func)(ulong*, ulong*, ulong, const zn_mod_t);


/*
   Profiles one of the butterfly routines.

   arg points to an array of ulongs:
      * First is 0 for add, 1 for subtract, 2 for inplace butterfly.
      * Second is 0 for safe version, 1 for slim version.
   
   Returns total cycle count for _count_ calls to butterfly of length 1000.
*/
double
profile_bfly (void* arg, unsigned long count)
{
   ulong type = ((ulong*) arg)[0];
   ulong speed = ((ulong*) arg)[1];
   ulong m = 123 + (1UL << (ULONG_BITS - (speed ? 2 : 1)));
   
   zn_mod_t mod;
   zn_mod_init (mod, m);
   
   const ulong n = 1000;
   
   ulong* buf1 = (ulong*) malloc (sizeof (ulong) * n);
   ulong* buf2 = (ulong*) malloc (sizeof (ulong) * n);

   // generate random inputs
   size_t i;
   for (i = 0; i < n; i++)
      buf1[i] = random_ulong (m);
   for (i = 0; i < n; i++)
      buf2[i] = random_ulong (m);

   bfly_func target;

   if (type == 0)
      target = (bfly_func) zn_array_add_inplace;
   else if (type == 1)
      target = (bfly_func) zn_array_sub_inplace;
   else   // type == 2
      target = zn_array_bfly_inplace;

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      target (buf1, buf2, n, mod);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      target(buf1, buf2, n, mod);

   cycle_count_t t1 = get_cycle_counter ();
   
   free (buf2);
   free (buf1);
   
   zn_mod_clear (mod);

   return cycle_diff (t0, t1);
}



/*
   Profiles mpn_add_n or mpn_sub_n.

   arg points to a single ulong: 0 for mpn_add_n, 1 for mpn_sub_n.
   
   Returns total cycle count for _count_ calls to length 1000 call.
*/
double
profile_mpn_aors (void* arg, unsigned long count)
{
   ulong type = ((ulong*) arg)[0];

   const ulong n = 1000;
   
   mp_limb_t* buf1 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n);
   mp_limb_t* buf2 = (mp_limb_t*) malloc (sizeof (mp_limb_t) * n);

   mp_limb_t (*target)(mp_limb_t*, const mp_limb_t*, const mp_limb_t*,
                       mp_size_t n);
   
   target = type ? mpn_sub_n : mpn_add_n;

   // generate random inputs
   size_t i;
   for (i = 0; i < n; i++)
      buf1[i] = random_ulong (1UL << (ULONG_BITS - 1));
   for (i = 0; i < n; i++)
      buf2[i] = random_ulong (1UL << (ULONG_BITS - 1));

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      target (buf1, buf1, buf2, n);

   // do the actual profile
   cycle_count_t t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      target (buf1, buf1, buf2, n);

   cycle_count_t t1 = get_cycle_counter ();

   free (buf2);
   free (buf1);

   return cycle_diff (t0, t1);
}


/*
   Profiles scalar multiplication.

   arg points to an array of ulongs:
      * First is modulus size in bits.
      * Second is 0 for regular multiply, 1 for REDC multiply
   
   Returns total cycle count for _count_ calls to zn_array_scalar_mul
   of length 1000.
*/
double
profile_scalar_mul (void* arg, unsigned long count)
{
   int bits = ((ulong*) arg)[0];
   int algo = ((ulong*) arg)[1];
   
   zn_mod_t mod;
   ulong m = random_modulus (bits, 1);
   zn_mod_init (mod, m);
   
   ulong x = random_ulong (m);
   const ulong n = 1000;

   // generate random input
   ulong* buf = (ulong*) malloc (sizeof (ulong) * n);
   size_t i;
   for (i = 0; i < n; i++)
      buf[i] = random_ulong (m);

   cycle_count_t t0, t1;

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      _zn_array_scalar_mul (buf, buf, n, x, algo, mod);

   // do the actual profile
   t0 = get_cycle_counter ();

   for (j = 0; j < count; j++)
      _zn_array_scalar_mul (buf, buf, n, x, algo, mod);

   t1 = get_cycle_counter ();

   free (buf);
   zn_mod_clear (mod);

   return cycle_diff (t0, t1);
}



// end of file ****************************************************************
