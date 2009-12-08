/*
   misc.c:  various random things that don't belong anywhere else
   
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

#include "zn_poly_internal.h"


char* ZNP_version_string = "0.9";


const char*
zn_poly_version_string ()
{
   return ZNP_version_string;
}


int
floor_lg (ulong x)
{
   int result = -1;
   while (x)
   {
      x >>= 1;
      result++;
   }
   
   return result;
}


int
ceil_lg (ulong x)
{
   ZNP_ASSERT (x >= 1);

   return floor_lg (x - 1) + 1;
}


// end of file ****************************************************************
