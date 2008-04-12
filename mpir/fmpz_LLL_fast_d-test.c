/*
  Copyright 2005, 2006 Damien Stehlé.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctype.h>
#include <limits.h>
#include <gmp.h>

#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpz_LLL_fast_d.h"
#include "test_support.h"

/* 
   Intrel generates a random lattice of the kind:  (speaking with columns)
   x1   x2 ....   xd
   1    0         0
   0    1         0
   ...  
   0    0         1 
   
   with x_i <= 2^bits, a random integer.
*/

void intrel(fmpz_mat_t B, int bits, int d)
{
   long i, j;
   
   for (i = 0; i < d; i++)
   {
      fmpz_t * row = B->row_arr[i];
      fmpz_random(row, bits);
      for (j = 1; j < i+1; j++)
	     fmpz_set_ui(row + j, 0L);
      fmpz_set_ui(row + (i+1), 1L);
      for (j = i+2; j < d+1; j++)
	     fmpz_set_ui(row + j, 0L);
   }
}

/* 
   Simdioph generates a random lattice of the kind:  (speaking with rows)
   2^bits2  x1    x2 ....  xd
   0      2^bits  0  ....   0
   0        0   2^bits      0
   ...
   0        0     0       2^bits 
   
   with x_i <= 2^bits, a random integer.
*/

void simdioph(fmpz_mat_t B, int bits, int bits2, int d)
{
   int i, j;

   fmpz_t * row = B->row_arr[0];
   fmpz_set_ui(row, 1L);
   fmpz_mul_2exp(row, row, bits2);
   for (i = 1; i < d+1; i++) fmpz_random(row + i, bits);
   for (i = 1; i < d+1; i++)
   {
      row = B->row_arr[i];
      for (j = 1; j < i; j++) fmpz_set_ui(row + j, 0L);
      fmpz_set_ui(row + i, 1L);
      fmpz_mul_2exp(row + i, row + i, bits);
      for (j = i + 1; j < d + 1; j++) fmpz_set_ui(row + j, 0L);
   }
}


/* 
   Generates "random" lattices according to Goldstein-Mayer's algorithm:
   On the equidistribution of Hecke points, Forum Mathematicum, 15:165-189, 2003
*/

void goldstein_mayer(fmpz_mat_t B, int bits, int d)
{
   int i, j;
   fmpz_t * row;
   
   row = B->row_arr[0];
   do fmpz_random(row, bits);
   while ((fmpz_probab_prime_p(row, 10)) == 0);
   for (i = 1; i < d; i++) 
      fmpz_set_ui(row + i, 0L);

   for (i = 1; i < d; i++)
   {
      row = B->row_arr[i];
      fmpz_randomm(row, B->row_arr[0]);
      for (j = 1; j < i; j++)
	     fmpz_set_ui(row + j, 0L);
      fmpz_set_ui(row + i, 1L);
      for (j = i + 1; j < d; j++)
	     fmpz_set_ui(row + j, 0L);
   }
}

void uniform(fmpz_mat_t B, int bits, int d)
{
   int i, j;

   for (i=0; i < d; i++)
   {
     fmpz_t * row = B->row_arr[i];
     for (j = 0; j < d; j++)
         fmpz_random(row +j, bits);
   }
}

/*void ntrulike(fmpz_mat_t B, int bits, int d, int q)
{
   int i, j, k;
   Z_NR * h;

   h = init_vec (d+1);
  
   for (i=0; i<d; i++)
      ZRANDB (h[i], state, bits);
  
   for (i=1; i<=d; i++)
   {
      for (j=1; j<i; j++)
	     ZSET_UI (B->coeff[i][j], 0);
      ZSET_UI (B->coeff[i][i], 1);
      for (j=i+1; j<=d; j++)
	     ZSET_UI (B->coeff[i][j], 0);
   }

   for (i=d+1; i<=2*d; i++)
      for (j=1; j<=d; j++)
         ZSET_UI (B->coeff[i][j], 0);

   for (i=d+1; i<=2*d; i++)
   {
      for (j=d+1; j<i; j++)
	     ZSET_UI (B->coeff[i][j], 0);
      ZSET_UI (B->coeff[i][i], q);
      for (j=i+1; j<=2*d; j++)
	     ZSET_UI (B->coeff[i][j], 0);
   }

   for (i=1; i<=d; i++)
      for (j=d+1; j<=2*d; j++)
      { 
	     k = j+i-2;
	     while (k>=d)
	        k -=d;
         ZSET (B->coeff[i][j], h[k]);
      }

   clear_vec(h, d+1);
}*/

/* 
   The same as ntrulike, but with the q-vectors coming first, and written
   as a lower-triangular matrix 
*/

/*void ntrulike2(fmpz_mat_t B, int bits, int d, int q)
{
   int i, j, k;
   Z_NR * h;

   h = init_vec (d+1);
  
   for (i=0; i<d; i++)
      ZRANDB (h[i], state, bits);
  
   for (i=1; i<=d; i++)
      for (j=1; j<=2*d; j++)
         ZSET_UI (B->coeff[i][j], 0);

   for (i=1; i<=d; i++)
      ZSET_UI (B->coeff[i][i], q);


   for (i=d+1; i<=2*d; i++)
      for (j=d+1; j<=2*d; j++)
         ZSET_UI (B->coeff[i][i], 0);
      
   for (i=d+1; i<=2*d; i++)
      ZSET_UI (B->coeff[i][i], 1);

   for (i=d+1; i<=2*d; i++)
      for (j=1; j<=d; j++)
      { 
	     k = i+j-2;
	     while (k>=d)
	        k -=d;
	     ZSET (B->coeff[i][j], h[k]);
      }

   clear_vec (h, d+1);
}*/

/*
   Ajtai generates a random lattice of the kind: (speaking with rows)
   a1     00000000000000000000000000
   a1/2   a2    00000000000000000000
   a1/2   a2/2   a3    0000000000000 
   ...
   a1/2   a2/2   a3/2  ......     ad
   
   with ai = random( 2^((2d-i+1)^alpha))   
   a1/2 is random(a1/2).
*/

/*void ajtai (fmpz_mat_t B, int d, double alpha)
{
   int i, j, bits;
   Z_NR ztmp;

   ZINIT (ztmp);
   for (i=1; i<=d; i++)
   {
      bits = (int) pow((double) (2*d-i+1), alpha);
      ZSET_UI (ztmp, 1);
      ZMUL_2EXP (ztmp, ztmp, bits);	  
      ZRANDM (B->coeff[i][i], state, ztmp);
      ZADD_UI ( B->coeff[i][i], B->coeff[i][i], 1);
      ZDIV_2EXP (ztmp, B->coeff[i][i], 1);
      for (j=i+1; j<=d; j++)
	  {
	     ZRANDM (B->coeff[j][i], state, ztmp);
	     ZSET_UI (B->coeff[i][j], 0);
      }
   }
   ZCLEAR (ztmp);
}*/

/* 
   Generates a random unimodular matrix U by performing
   random elementary operations Li <- Li \pm Lj, starting from
   the identity matrix, and stopping when the mean of the absolute
   values of the entries of U has more than bits bits. 
*/

/*void unimodular(fmpz_mat_t B, int d, int bits)
{
   int i, j, k, bbits;

   for (i=1; i<=d; i++)
      for (j=1; j<=d; j++)
      {
	     if (i==j) {ZSET_UI(B->coeff[i][i], 1);}
	     else {ZSET_UI(B->coeff[i][j], 0);}
      }

   do
   {
      i = rand()%d+1;
      j = rand()%d+1;
      if (j!=i)
	  {
	     k = rand();
	     if (k%2==0)
	     {
            for (k=1; k<=d;k++)
		       ZADD (B->coeff[i][k], B->coeff[i][k], B->coeff[j][k]);
	     } else
	     {
	        for (k=1; k<=d;k++)
	           ZSUB (B->coeff[i][k], B->coeff[i][k], B->coeff[j][k]);
	     }
	  }

      bbits=0;
      for (i=1; i<=d; i++)
	     for (j=1; j<=d; j++)
	        bbits+=ZSIZEINBASE2 (B->coeff[i][j]);
   } while (bbits<d*d*bits);
}*/

/* ********************** */
/*  MAIN **************** */
/* ********************** */

int main(int argc, char *argv[])
{
   gmp_randinit_default(state);
   
   int d=0, i, j, n, bits, bits2=0, decal;
   double alpha;
   fmpz_mat_t B;
   char c = 0;
   argc = argc;

   decal = 1;

   ctt = DELTA;
   halfplus = ETA;

   if (strcmp(argv[1],"-delta")==0) 
   {
      decal+=2;
      ctt = atof(argv[2]);
   }

   if (strcmp(argv[decal],"-eta")==0)
   {
      halfplus = atof(argv[decal+1]);
      decal+=2;
   }

   if (strcmp(argv[decal],"-delta")==0)
   {
      ctt = atof(argv[decal+1]);
      decal+=2;
   }

   onedothalfplus = 1.0+halfplus;

   if ((ctt >=1.0 ) || (ctt <=0.25) || (halfplus <=0.5) || (halfplus>=sqrt(ctt)))
   {
      fprintf (stderr, "Incorrect parameters! You must choose\ndelta in (0.25,1) and eta in (0.5, sqrt(delta))\n"); 
      abort();
   }

   if ( strcmp(argv[decal], "h")==0) 
   {
      printf("Please have a look at the README file.\n");
      abort();
   }
  
   if (strcmp(argv[decal], "m")==0)
   {
      d = atoi(argv[decal+1]);
      n = atoi(argv[decal+2]);
      fmpz_mat_init(B, d, n);

      /* read initial '[' */
      while (isspace(c = getchar ()) || c == '\n');
      if (c != '[')
      {
         fprintf (stderr, "Error: '[' expected\n");
         abort();
      }

      for (i = 0; i < d; i++)
	  {
	     while (isspace(c = getchar ()) || c == '\n');
	     if (c != '[')
	     {
	        fprintf (stderr, "Error at row %d: '[' expected instead of %c\n", i, c);
            abort();
	     }
	     for (j = 0; j < n; j++) fmpz_read(B->row_arr[i] + j);

	     while (isspace(c = getchar ()) || c=='\n');
	     if (c != ']')
         {
	        fprintf (stderr, "Error: ']' expected at line %u\n", i);
	        abort();
         }
      }
   }

   if (strcmp(argv[decal], "r")==0)
   {
      d = atoi(argv[decal+1]);
      bits = atoi(argv[decal+2]);
      fmpz_mat_init(B, d, d+1);
      intrel(B, bits, d); 
   }      

   if (strcmp(argv[decal], "s")==0)
   {
      d = atoi(argv[decal+1]);
      bits = atoi(argv[decal+2]);
      bits2 = atoi(argv[decal+3]);
      fmpz_mat_init(B, d+1, d+1);
      simdioph(B, bits, bits2, d);
   }   


   /*if (strcmp(argv[decal], "a")==0)
   {
      d = atoi(argv[decal+1]);
      alpha = atof(argv[decal+2]);
      fmpz_mat_init(B, d, d);
      ajtai(B, d, alpha);
   } */ 

   if (strcmp(argv[decal], "u")==0)
   {
      d = atoi(argv[decal+1]);
      bits = atoi(argv[decal+2]);
      fmpz_mat_init(B, d, d);
      uniform(B, bits, d); 
   }  

   if (strcmp(argv[decal], "g")==0)
   {
      d = atoi(argv[decal+1]);
      bits = atoi(argv[decal+2]);
      fmpz_mat_init(B, d, d);
      goldstein_mayer(B, bits, d); 
   }  

   /*if (strcmp(argv[decal], "n")==0)
   {
      d = atoi(argv[decal+1]);
      bits = atoi(argv[decal+2]);
      i = atoi(argv[decal+3]);
      fmpz_mat_init(B, 2*d, 2*d);
      ntrulike(B, bits, d, i); 
   }  

   if (strcmp(argv[decal], "n2")==0)
   {
      d = atoi(argv[decal+1]);
      bits = atoi(argv[decal+2]);
      i = atoi(argv[decal+3]);
      fmpz_mat_init(B, 2*d, 2*d);
      ntrulike2(B, bits, d, i); 
   }  */
    
   //fmpz_mat_print(B, B->rows, B->cols); printf("\n");
   LLL(B);
   //fmpz_mat_print(B, B->rows, B->cols); printf("\n");
   
   fmpz_mat_clear(B);
   
   gmp_randclear(state);

   return 0;
}
