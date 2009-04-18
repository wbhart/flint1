/*
   ntl-profile-dummy.c:  dummy routines replacing NTL profiling routines when
                         no NTL support is compiled in
   
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

#include <stdlib.h>
#include <stdio.h>


void
no_ntl_support ()
{
   printf ("\n\n");
   printf ("no NTL profiling support compiled in!\n");
   abort ();
}


double
profile_mul_ntl (void* arg, unsigned long count)
{
   no_ntl_support ();
}

double
profile_invert_ntl (void* arg, unsigned long count)
{
   no_ntl_support ();
}


// end of file ****************************************************************
