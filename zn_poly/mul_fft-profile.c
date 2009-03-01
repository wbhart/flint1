/*
   mul_fft-profile.c:  profiling for mul_fft.c module
   
   Copyright (C) 2007, David Harvey
   
   This file is part of the zn_poly library (version 0.4.1).
   
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
#include "zn_poly_internal.h"
#include "profiler.h"


/*
   Profiles the main mul_fft() routine.

   arg points to an array of ulongs.
   First ulong is poly length.
   Second is modulus n.
   
   Returns total cycle count for _count_ calls.
*/
double profile_mul_fft(void* arg, unsigned long count)
{
   size_t len = ((ulong*) arg)[0];
   ulong n = ((ulong*) arg)[1];
   
   zn_mod_t mod;
   zn_mod_init(mod, n);
   
   ulong* buf1 = (ulong*) malloc(sizeof(ulong) * len);
   ulong* buf2 = (ulong*) malloc(sizeof(ulong) * len);
   ulong* buf3 = (ulong*) malloc(sizeof(ulong) * 2 * len);
   
   size_t i;
   for (i = 0; i < len; i++)
      buf1[i] = random_ulong(n);
   for (i = 0; i < len; i++)
      buf2[i] = random_ulong(n);

   // warm up
   ulong j;
   for (j = 0; j < count; j++)
      zn_array_mul_fft(buf3, buf1, len, buf2, len, mod);
      
   cycle_count_t t0 = get_cycle_counter();

   for (j = 0; j < count; j++)
      zn_array_mul_fft(buf3, buf1, len, buf2, len, mod);

   cycle_count_t t1 = get_cycle_counter();
   
   free(buf3);
   free(buf2);
   free(buf1);
   
   zn_mod_clear(mod);
   
   return cycle_diff(t0, t1);
}



void prof_main(int argc, char* argv[])
{
   double result;
   unsigned long arg[2];
   unsigned bitsizes[9] = {4, 8, 16, 24, 32, 40, 48, 56, 64};

   int i, j;
   size_t len;
   unsigned bits;

   // loop over bitsizes in above table
   for (i = 0; i < sizeof(bitsizes) / sizeof(bitsizes[0]); i++)
   {
      bits = bitsizes[i];

      // loop over lengths, spaced out logarithmically
      for (j = 0; j < 120; j++)
      {
         size_t new_len = (size_t) floor(pow(1.1, (double) j));
         if (new_len == len)
            continue;
         len = new_len;

         arg[0] = len;

         // choose an odd modulus exactly _bits_ bits long, and not equal to
         // a power of two (so that the residues are _bits_ bits long too),
         ulong n = (1UL << (bits - 1))
                      + 2 * random_ulong((1UL << (bits - 2))) + 1;
         arg[1] = n;

         result = profile(profile_mul_fft, arg);
         
         printf("len = %5lu, bits = %2u, mul_fft = %.2le\n",
                len, bits, result);
      }
      printf("-------------------------------------------\n");
   }
}

// end of file ****************************************************************
