/*============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

   theta.c: theta functions

   Copyright (C) 2007, Gonzalo Tornaria and William Hart

*****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "flint.h"
#include "theta.h"
#include "long_extras.h"

/*
   Print theta series as q expansion to standard out
*/

void theta_print(long * coeff, ulong start, ulong len) 
{
    theta_fprint(coeff, start, len, stdout);
}

/*
   Print theta series as q expansion -- it aims to be equal to gp
   printout of power series
*/

void theta_fprint(long * coeff, ulong start, ulong len, FILE * f)
{
    // flag first term
    int first = 1;

    for(ulong i=0; i<len; i++) {
        ulong n = i + start;
        long c = coeff[i];

        if(!c) continue;

        // print operator
        if (first) {
            first = 0;
            if(c<0) fprintf(f, "-");
        } else {
            if(c<0) fprintf(f, " - ");
            else fprintf(f, " + ");
        }

        if (c<0) c=-c;

        // print coefficient
        if (c != 1 || n == 0)
            fprintf(f, "%ld", c);

        if (c != 1 && n != 0)
            fprintf(f, "*");

        if (n == 1)
            fprintf(f, "q", n);

        if (n > 1)
            fprintf(f, "q^%ld", n);
    }
    fprintf(f, " + O(q^%ld)\n", start+len);
}

/*
   Compute Theta(z)
   i.e. theta series of  x^2
*/

void theta(long * out, ulong start, ulong len)
{

    // zero all the coefficients
    for(ulong i=0; i<len; i++)
        out[i] = 0;

	ulong x = (start == 0 ? 0 : z_intsqrt(start-1)+1); // start <= x^2
    ulong i = x * x - start; // x^2 - start

    while(i < len) {
        out[i] = 2;

        // iterate x++ and update i = x^2 - start
        i += x;
        x += 1;
        i += x;
    }

    if(start == 0)
        out[0] = 1; // the constant term is 1 not 2
}

/*
   Compute theta with a quasi-character modulo 6
   i.e. theta series of  x^2, summed only over non-negative exponents
	note the character starts from the starting position (start), not from
	zero.
*/

void theta_mod6_char(long * out, long * character, ulong start, ulong len)
{

    // zero all the coefficients
    for(ulong i = 0; i < len; i++)
        out[i] = 0;

	ulong x = (start == 0 ? 0 : z_intsqrt(start-1)+1); // start <= x^2
   ulong i = x * x - start; // x^2 - start

    while(i < len) {
        out[i] = character[i%6];

        // iterate x++ and update i = x^2 - start
        i += x;
        x += 1;
        i += x;
    }
}

/*
   Compute Sum q^{ax^2+bx+c} 
   i.e. theta series Sum q^{ax^2+bx+c}
   Assumes that start is either 0 or quite large compared to a and b
   Also assumes a, b and c are positive
*/

void theta_1d(unsigned long a, unsigned long b, unsigned long c,
				                                          long * out, ulong start, ulong len)
{
   // zero all the coefficients
   for(ulong i=0; i<len; i++)
      out[i] = 0;

   ulong i, x;
   
   if (start == 0)
   {
      long disc = b*b - 4*a*c;

	  if (disc < 0L) // no intersection with x-axis
	  {
	     x = (b/a);
         
		 i = a*x*x - b*x + c;

		 while (i < len)
         {
            out[i] = 1;
	  
	        i += (a*x - b);
	        x += 1;
			i += (a*x);
         }
         
		 x = (b/a);
		 
		 i = a*x*x - b*x + c;

		 while ((x > 0) && (i < len))
		 {
		 	i -= (a*x - b);
	        x -= 1;
			i -= (a*x);
            
			out[i] += 1;
		 }

		 x = 1;

		 i = a*x*x + b*x + c;

		 while (i < len)
         {
            out[i] += 1;
	  
	        i += (a*x + b);
	        x += 1;
		    i += (a*x);
         }  
	  } else // discriminant >= 0
	  {
	     x = (b/a);
         
		 long j = i = a*x*x - b*x + c;

		 while (j < (long)len)
         {
            if (j >= 0L) out[j] = 1;
	  
	        j += (a*x - b);
	        x += 1;
			j += (a*x);
         }
         
		 x = (b/a);
		 
		 j = a*x*x - b*x + c;

		 while ((x > 0) && (j < (long) len))
		 {
		 	j -= (a*x - b);
	        x -= 1;
			j -= (a*x);
            
			if (j >= 0L) out[j] += 1;
		 }

		 x = 1;

		 j = a*x*x + b*x + c;

		 while (j < (long) len)
         {
            if (j >= 0L) out[j] += 1;
	  
	        j += (a*x + b);
	        x += 1;
		    j += (a*x);
         }
	  }
   } else
   {
      x = z_intsqrt(start/a-1)+1; // start <= ax^2

      while (a*x*x + b*x + c >= start) x--; 
      x++; // start <= ax^2+bx+c

      i = a*x*x + b*x + c - start;

      while (i < len)
      {
         out[i] = 1;
	  
	     i += (a*x + b);
	     x += 1;
		 i += (a*x);
      }

      x = z_intsqrt(start/a-1)+1; // start <= ax^2
   
      while (a*x*x - b*x + c < start) x++; // start <= ax^2-bx+c (deals with negative x)

      i = a*x*x - b*x + c - start;

      while (i < len)
      {
         out[i] += 1;
	  
	     i += (a*x - b);
	     x += 1;
		 i += (a*x);
      }
   }
}


/*
   Compute Sum q^{ax^2+bx+c} with a quadratic character mod 4
   i.e. theta series Sum a_{x mod 4} q^{ax^2+bx+c} where a_i = quad[i], i = 0..3 are in -1, 0, 1
   Assumes that start is either 0 or quite large compared to a an b
   Also assumes a, b and c are positive
*/

void theta_1d_quadchar(long * character, unsigned long a, unsigned long b, unsigned long c,
				                                          long * out, ulong start, ulong len)
{
   // zero all the coefficients
   for(ulong i=0; i<len; i++)
      out[i] = 0;

   ulong i, x;
   
   if (start == 0)
   {
      long disc = b*b - 4*a*c;

	  if (disc < 0L) // no intersection with x-axis
	  {
	     x = (b/a);
         
		 i = a*x*x - b*x + c;

		 while (i < len)
         {
            out[i] = character[(-x)&3];
	  
	        i += (a*x - b);
	        x += 1;
			i += (a*x);
         }
         
		 x = (b/a);
		 
		 i = a*x*x - b*x + c;

		 while ((x > 0) && (i < len))
		 {
		 	i -= (a*x - b);
	        x -= 1;
			i -= (a*x);
            
			out[i] += character[(-x)&3];
		 }

		 x = 1;

		 i = a*x*x + b*x + c;

		 while (i < len)
         {
            out[i] += character[x&3];
	  
	        i += (a*x + b);
	        x += 1;
		    i += (a*x);
         }  
	  } else // discriminant >= 0
	  {
	     x = (b/a);
         
		 long j = a*x*x - b*x + c;

		 while (j < (long)len)
         {
            if (j >= 0L) out[j] = character[(-x)&3];
	  
	        j += (a*x - b);
	        x += 1;
			j += (a*x);
         }
         
		 x = (b/a);
		 
		 j = a*x*x - b*x + c;

		 while ((x > 0) && (j < (long) len))
		 {
		 	j -= (a*x - b);
	        x -= 1;
			j -= (a*x);
            
			if (j >= 0L) out[j] += character[(-x)&3];
		 }

		 x = 1;

		 j = a*x*x + b*x + c;

		 while (j < (long) len)
         {
            if (j >= 0L) out[j] += character[x&3];
	  
	        j += (a*x + b);
	        x += 1;
		    j += (a*x);
         }
	  }
   } else
   {
      x = z_intsqrt(start/a-1)+1; // start <= ax^2

      while (a*x*x + b*x + c >= start) x--; 
      x++; // start <= ax^2+bx+c

      i = a*x*x + b*x + c - start;

      while (i < len)
      {
         out[i] = character[x&3];
	  
	     i += (a*x + b);
	     x += 1;
		 i += (a*x);
      }

      x = z_intsqrt(start/a-1)+1; // start <= ax^2
   
      while (a*x*x - b*x + c < start) x++; // start <= ax^2-bx+c (deals with negative x)

      i = a*x*x - b*x + c - start;

      while (i < len)
      {
         out[i] += character[(-x)&3];
	  
	     i += (a*x - b);
	     x += 1;
		 i += (a*x);
      }
   }
}

/*
   Compute Sum_{i=0,infty} (2*a*x+b)/m * q^{ax^2+bx+c} with a quadratic character mod 4
   i.e. theta series Sum_{i=0, infty} a_{x mod 4} * (2*a*x+b)/m * q^{ax^2+bx+c} where a_i = quad[i], i = 0..3 
   are in -1, 0, 1
   Assumes that start is either 0 or quite large compared to a an b
   Also assumes a, b and c are positive
*/

void theta_1d_quadchar_0(long * character, unsigned long a, unsigned long b, unsigned long c,
				                        unsigned long m, long * out, ulong start, ulong len)
{
   unsigned long a1 = (2*a)/m;
   unsigned long b1 = b/m;

   // zero all the coefficients
   for(ulong i=0; i<len; i++)
      out[i] = 0;

   ulong i, x;
   long fac;
   
   if (start == 0)
   {
      long disc = b*b - 4*a*c;

	  if (disc < 0L) // no intersection with x-axis
	  {
		 x = 0;

		 fac = b1;

		 i = c;

		 while (i < len)
         {
            out[i] += fac*character[x&3];
	        
	        i += (a*x + b);
	        fac += a1;
			x += 1;
		    i += (a*x);
         }  
	  } else // discriminant >= 0
	  {
		 x = 0;

		 fac = b1;

		 long j = c;

		 while (j < (long) len)
         {
            if (j >= 0L) out[j] += fac*character[x&3];
			
	        j += (a*x + b);
	        fac += a1;
			x += 1;
		    j += (a*x);
         }
	  }
   } else
   {
      x = z_intsqrt(start/a-1)+1; // start <= ax^2

      while (a*x*x + b*x + c >= start) x--; 
      x++; // start <= ax^2+bx+c

      fac = a1*x + b1;

	  i = a*x*x + b*x + c - start;

      while (i < len)
      {
         out[i] = fac*character[x&3];
	     
	     i += (a*x + b);
	     fac += a1;
		 x += 1;
		 i += (a*x);
      }
   }
}


/*
   Compute Sum (2*a*x+b)/m * q^{ax^2+bx+c} with a quadratic character mod 4
   i.e. theta series Sum a_{x mod 4} * (2*a*x+b)/m * q^{ax^2+bx+c} where a_i = quad[i], i = 0..3 
   are in -1, 0, 1
   Assumes that start is either 0 or quite large compared to a an b
   Also assumes a, b and c are positive
*/

void theta_1d_quadchar_2(long * character, unsigned long a, unsigned long b, unsigned long c,
				                        unsigned long m, long * out, ulong start, ulong len)
{
   unsigned long a1 = (2*a)/m;
   unsigned long b1 = b/m;

   // zero all the coefficients
   for(ulong i=0; i<len; i++)
      out[i] = 0;

   ulong i, x;
   long fac;
   
   if (start == 0)
   {
      long disc = b*b - 4*a*c;

	  if (disc < 0L) // no intersection with x-axis
	  {
	     x = (b/a);
         
		 i = a*x*x - b*x + c;
         
		 fac = -a1*x + b1;

		 while (i < len)
         {
            out[i] = fac*character[(-x)&3];
	        
	        i += (a*x - b);
			fac -= a1;
	        x += 1;
			i += (a*x);
         }
         
		 x = (b/a);
		 
		 fac = -a1*x + b1;

		 i = a*x*x - b*x + c;

		 while ((x > 0) && (i < len))
		 {
		 	i -= (a*x - b);
			fac -= a1;
	        x -= 1;
			i -= (a*x);
            
			out[i] += fac*character[(-x)&3];
		 }

		 x = 1;

		 fac = a1*x + b1;

		 i = a*x*x + b*x + c;

		 while (i < len)
         {
            out[i] += fac*character[x&3];
	        
	        i += (a*x + b);
	        fac += a1;
			x += 1;
		    i += (a*x);
         }  
	  } else // discriminant >= 0
	  {
	     x = (b/a);
         
		 fac = a1*x + b1;

		 long j = a*x*x - b*x + c;

		 while (j < (long)len)
         {
            if (j >= 0L) out[j] = fac*character[(-x)&3];
	        
	        j += (a*x - b);
			fac += a1;
	        x += 1;
			j += (a*x);
         }
         
		 x = (b/a);
		 
		 fac = -a1*x + b1;

		 j = a*x*x - b*x + c;

		 while ((x > 0) && (j < (long) len))
		 {
		 	j -= (a*x - b);
	        fac -= a1;
			x -= 1;
			j -= (a*x);
            
			if (j >= 0L) out[j] += fac*character[(-x)&3];
		 }

		 x = 1;

		 fac = a1*x + b1;

		 j = a*x*x + b*x + c;

		 while (j < (long) len)
         {
            if (j >= 0L) out[j] += fac*character[x&3];
			
	        j += (a*x + b);
	        fac += a1;
			x += 1;
		    j += (a*x);
         }
	  }
   } else
   {
      x = z_intsqrt(start/a-1)+1; // start <= ax^2

      while (a*x*x + b*x + c >= start) x--; 
      x++; // start <= ax^2+bx+c

      fac = a1*x + b1;

	  i = a*x*x + b*x + c - start;

      while (i < len)
      {
         out[i] = fac*character[x&3];
	     
	     i += (a*x + b);
	     fac += a1;
		 x += 1;
		 i += (a*x);
      }

      x = z_intsqrt(start/a-1)+1; // start <= ax^2
   
      while (a*x*x - b*x + c < start) x++; // start <= ax^2-bx+c (deals with negative x)

      fac = -a1*x + b1;

	  i = a*x*x - b*x + c - start;

      while (i < len)
      {
         out[i] += fac*character[(-x)&3];
	     
	     i += (a*x - b);
	     fac -= a1;
		 x += 1;
		 i += (a*x);
      }
   }
}

/*
   Compute Theta(8z) - 2*Theta(32z)
   i.e. theta series of  8x^2 with coefficient +2 for odd x, -2 for even x
*/

void theta_1d_B(long * out, ulong start, ulong len)
{ 
	// for now I assume start = 0 ; shift = 1

    // zero all the coefficients
    for(ulong i=0; i<len; i++)
        out[i] = 0;

    ulong x = (start < 8 ? 0 : z_intsqrt(start/8-1)+1); // start <= 8x^2
	if (8*x*x < start) x++;

    ulong x_8 = 8 * x; // 8x
    ulong i = x_8 * x - start; // 8x^2 - start
    long val = (start+i) & 8 ? 2 : -2;

    while(i < len) {
        out[i] = val;

        // iterate x++ and update i = 8x^2 - start
        i += x_8;
        x_8 += 8;
        i += x_8;

        val = -val;
    }

    if(start == 0)
        out[0] = -1; // the constant term is -1 not -2
}

/*
   Compute Theta(z) ^ 2
   i.e. theta series of  x^2+y^2
*/

void theta_2d(long *out, ulong start, ulong len)
{
    for(ulong i=0; i<len; i++)
        out[i] = 0;

    long A = start;
    long B = start+len-1;

    long y = 1;
    long N0 = 1; // N0 = y^2

    long x = z_intsqrt(B-N0); // largest x such that x^2 + y^2 <= B

    // First loop, where x won't reach 0
    while(N0 < A) {

        long N = N0 + x * x; // x^2 + y^2

        while(N > B) {
            // iterate x-- and update N = x^2 + y^2
            N -= x;
            x -= 1;
            N -= x;
        }

        while(N >= A) {
            out[N-A] += 4;

            // iterate x-- and update N = x^2 + y^2
            N -= x;
            x -= 1;
            N -= x;
        }

        // iterate y++ and update N0 = y^2
        N0 += y;
        y += 1;
        N0 += y;
    }

    // Second loop, where x starts from 0
    while(N0 <= B) {

        // count x = 0
        out[N0-A] += 4;
        
        // start from x=1
        long x = 1;
        long N = N0 + x*x; // x^2 + y^2

        while(N <= B) {
            out[N-A] += 4;

            // iterate x+=1 and update N = x^2 + y^2
            N += x;
            x += 1;
            N += x;
        }

        // iterate y++ and update N0 = y^2
        N0 += y;
        y += 1;
        N0 += y;
    }

	if (start == 0) out[0] = 1;
}

/*
   Compute Theta(2z) * (Theta(z) - Theta(4z))
   i.e. theta series of  2x^2+y^2  with odd y
*/

void theta_2d_A1(long *out, ulong start, ulong len)
{
    for(ulong i=0; i<len; i++)
        out[i] = 0;

    long A = start;
    long B = start+len-1;

    long y_2 = 2; // 2 * y
    long N0 = 1; // N0 = y^2

    long x = z_intsqrt((B-1)/2); // largest x such that 2x^2 + y^2 <= B

    // First loop, where x won't reach 0
    while(N0 < A) {

        long x_2 = 2 * x;  // 2x
        long N = N0 + x_2 * x; // 2x^2 + y^2

        while(N > B) {
            // iterate x-- and update N = 2x^2 + y^2
            N -= x_2;
            x_2 -= 2;
            N -= x_2;
        }

        // remember x for next iteration of y
        x = x_2 / 2;

        while(N >= A) {
            out[N-A] += 2;

            // iterate x-- and update N = 2x^2 + y^2
            N -= x_2;
            x_2 -= 2;
            N -= x_2;
        }

        // iterate y+=2 and update N0 = y^2
        N0 += y_2;
        y_2 += 4;
        N0 += y_2;
    }

    // Second loop, where x starts from 0
    while(N0 <= B) {

        out[N0-A] += 1;
        
        long x_2 = 2;  // 2x
        long N = N0 + 2; // 2x^2 + y^2

        while(N <= B) {
            out[N-A] += 2;

            // iterate x++ and update N = 2x^2 + y^2
            N += x_2;
            x_2 += 2;
            N += x_2;
        }

        // iterate y+=2 and update N0 = y^2
        N0 += y_2;
        y_2 += 4;
        N0 += y_2;
    }
}

/*
   Compute Theta(4z) * (Theta(z) - Theta(4z))
   i.e. theta series of  4x^2+y^2  with odd y
*/

void theta_2d_A2(long *out, ulong start, ulong len)
{
    for(ulong i=0; i<len; i++)
        out[i] = 0;

    long A = start;
    long B = start+len-1;

    long y_2 = 2; // 2y
    long N0 = 1; // N0 = y^2

    long x = z_intsqrt((B-1)/2); // largest x such that 2x^2 + y^2 <= B

    // First loop, where x won't reach 0
    while(N0 < A) {

        long x_4 = 4 * x;  // 4x
        long N = N0 + x_4 * x; // 4x^2 + y^2

        while(N > B) {
            // iterate x-- and update N = 4x^2 + y^2
            N -= x_4;
            x_4 -= 4;
            N -= x_4;
        }

        // remember x for next iteration of y
        x = x_4 / 4;

        while(N >= A) {
            out[N-A] += 2;

            // iterate x-- and update N = 4x^2 + y^2
            N -= x_4;
            x_4 -= 4;
            N -= x_4;
        }

        // iterate y+=2 and update N0 = y^2
        N0 += y_2;
        y_2 += 4;
        N0 += y_2;
    }

    // Second loop, where x starts from 0
    while(N0 <= B) {

        out[N0-A] += 1;
        
        long x_4 = 4;  // 4x
        long N = N0 + 4; // 4x^2 + y^2

        while(N <= B) {
            out[N-A] += 2;

            // iterate x++ and update N = 4x^2 + y^2
            N += x_4;
            x_4 += 4;
            N += x_4;
        }

        // iterate y+=2 and update N0 = y^2
        N0 += y_2;
        y_2 += 4;
        N0 += y_2;
    }
}

/*
   Compute Sum  s * y * q ^(x^2 + y^2)
   for x even, y odd, s = +1 for x+y = 1 mod 4, s = -1 for x+y = 3 mod 4
*/

void theta_2d_B(long *out, ulong start, ulong len)
{
    for(ulong i=0; i<len; i++)
        out[i] = 0;

    long A = start;
    long B = start+len-1;

    long y = 1;
    long y_2 = 2 * y; // 2 * y
    long N0 = 1; // N0 = y^2

    long x = z_intsqrt((B-1))/2 * 2; // largest even x such that x^2 + y^2 <= B

    // First loop, where x won't reach 0
    while(N0 < A) {

        long x_2 = 2 * x;  // 2x
        long N = N0 + x * x; // x^2 + y^2

        while(N > B) {
            // iterate x-=2 and update N = x^2 + y^2
            N -= x_2;
            x_2 -= 4;
            N -= x_2;
        }

        // remember x for next iteration of y
        x = x_2 / 2;

        long val = (x+y)%4 == 1 ? y : -y;
        while(N >= A) {
            out[N-A] += 2*val;

            // iterate x-=2 and update N = x^2 + y^2
            N -= x_2;
            x_2 -= 4;
            N -= x_2;

            val = -val;
        }

        // iterate y+=2 and update N0 = y^2
        y += 2;
        N0 += y_2;
        y_2 += 4;
        N0 += y_2;
    }

    // Second loop, where x starts from 0
    while(N0 <= B) {

        long val = y%4 == 1 ? y : -y;

        // count x = 0
        out[N0-A] += val;
        
        // start from x=2
        long x_2 = 4;  // 2x
        long N = N0 + 4; // x^2 + y^2

        while(N <= B) {
            val = -val;
            out[N-A] += 2*val;

            // iterate x+=2 and update N = x^2 + y^2
            N += x_2;
            x_2 += 4;
            N += x_2;
        }

        // iterate y+=2 and update N0 = y^2
        y += 2;
        N0 += y_2;
        y_2 += 4;
        N0 += y_2;
    }
}

/*
   Compute theta series Sum (-1)^(x+y) x^2+y^2
*/

void theta_2d_C(long *out, ulong start, ulong len)
{
    for(ulong i=0; i<len; i++)
        out[i] = 0;

    long A = start;
    long B = start+len-1;

    long y = 1;
    long N0 = 1; // N0 = y^2

    long x = z_intsqrt(B-N0); // largest x such that x^2 + y^2 <= B

    // First loop, where x won't reach 0
    while(N0 < A) {

        long N = N0 + x * x; // x^2 + y^2

        while(N > B) {
            // iterate x-- and update N = x^2 + y^2
            N -= x;
            x -= 1;
            N -= x;
        }

        while(N >= A) {
            if (N & 1) out[N-A] -= 4;
			else out[N-A] += 4;

            // iterate x-- and update N = x^2 + y^2
            N -= x;
            x -= 1;
            N -= x;
        }

        // iterate y++ and update N0 = y^2
        N0 += y;
        y += 1;
        N0 += y;
    }

    // Second loop, where x starts from 0
    while(N0 <= B) {

        // count x = 0
        if (N0 & 1) out[N0-A] -= 4;
		else out[N0-A] += 4;
        
        // start from x=1
        long x = 1;
        long N = N0 + x*x; // x^2 + y^2

        while(N <= B) {
            if (N & 1) out[N-A] -= 4;
			else out[N-A] += 4;

            // iterate x+=1 and update N = x^2 + y^2
            N += x;
            x += 1;
            N += x;
        }

        // iterate y++ and update N0 = y^2
        N0 += y;
        y += 1;
        N0 += y;
    }

	if (start == 0) out[0] = 1;
}
