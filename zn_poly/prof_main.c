/*
   prof_main.c:  main() routine for profiling programs
   
   Copyright (C) 2007, 2008, David Harvey
   
   This file is part of the zn_poly library (version 0.8).
   
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

#include "support.h"
#include "profiler.h"


/*
   Profiling programs link against this file, and implement a function
   void prof_main(argc, argv)
*/

int main(int argc, char* argv[])
{
   calibrate_cycle_scale_factor();

#if !ZNP_HAVE_CYCLE_COUNTER
   printf("Cannot run profiles; no cycle counter on this system!\n");
#else
   gmp_randinit_default(randstate);
   prof_main(argc, argv);
   gmp_randclear(randstate);
#endif
   
   return 0;
}


// end of file ****************************************************************
