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

  This program implements ideas from the paper "Floating-point LLL Revisited", 
  by Phong Nguyen and Damien Stehlé, in the Proceedings of Eurocrypt'2005, 
  Springer-Verlag; and was partly inspired by Shoup's NTL library: 
  http://www.shoup.net/ntl/

*/

/*
   Modifications: This code has been modified to work with mpir's fmpz "flat" 
   multiprecision integer package. 
   
   Copyright 2008 William Hart.
   
   Modifications distributed under the GPL.
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
#include "gmp.h"
#include "mpir.h"
#include "fmpz_mat.h"
#include "fmpz_LLL_fast_d.h"
#include "d_mat.h"

#define LOOPS_BABAI 10

/* computes the prec MSBs of the inner product of the approximate 
   vec1 and vec2. */

double fpScalarProduct(double *vec1, double *vec2, int n)
{
  int i;
  double sum;

  sum = vec1[0] * vec2[0];
  for (i = 1; i < n; i++)
     sum += vec1[i] * vec2[i];

  return sum;
} 

double fpNorm(double *vec, int n)
{
  int i;
  double sum;

  sum = vec[0] * vec[0];
  for (i = 1 ; i < n ; i++)
     sum += vec[i] * vec[i];

  return sum;

} 

ulong getShift(fmpz_mat_t B)
{
   ulong n = B->cols;
   ulong shift = 0;
   for (ulong i = 0; i < B->rows; i++)
   {
      fmpz_t * row = B->row_arr[i];

      ulong j;
      for (j = n-1; j >= 0 && fmpz_size(row + j) == 0L; j--);  
      
      if (shift < j-i) shift = j-i;
      
   }

   return shift;
}

/***********************************/
/* Babai's Nearest Plane algorithm */
/***********************************/

/* 
   Size-reduces b_kappa using mu_ij and r_ij for j<=i <kappa
   updates B (kappa)
   computes mu_kappaj, r_kappaj for j<=kappa, and s(kappa) 
   The algorithm is the iterative Babai algorithm of the paper
*/

void Babai (int kappa, fmpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n, fmpz_t * ztmp)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros+1;
  
   do
   {
      test=0;
      
#ifdef DEBUG
      if (loops++ > LOOPS_BABAI) 
      {
	     printf("INFINITE LOOP?\n"); 
	     abort();
      }
#endif 
      
      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j=aa; j<kappa; j++)
	  {	  
	     if (appSP[kappa][j]!=appSP[kappa][j])
	     {
	        appSP[kappa][j] = fpScalarProduct (appB[kappa], appB[j], n);
	     }
	  	  
         if (j > zeros+2)
	     {
	        tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	        rtmp = appSP[kappa][j] - tmp;
	        for (k=zeros+2; k<j-1; k++)
		    {
		       tmp = mu[j][k] * r[kappa][k];
		       rtmp = rtmp - tmp;
		    }
	        tmp = mu[j][j-1] * r[kappa][j-1];
	        r[kappa][j] = rtmp - tmp;
         } else if (j==zeros+2)
	     {
	        tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	        r[kappa][j] = appSP[kappa][j] - tmp;
	     } else r[kappa][j] = appSP[kappa][j];

	     mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j=kappa-1; j>zeros; j--)
	  {

	     /* test of the relaxed size-reduction condition */
	     tmp = fabs (mu[kappa][j]);
	     tmp = ldexp (tmp, expo[kappa]-expo[j]);
	  
	     if (tmp > halfplus) 
	     {
	        test = 1; 
	        exponent = expo[j] - expo[kappa];
	      
	        /* we consider separately the cases X = +-1 */     
	        if (tmp <= onedothalfplus)   
		    {		  
		       if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
               {
		          for (k=zeros+1; k<j; k++)
			      {
			         tmp = ldexp (mu[j][k], exponent);
			         mu[kappa][k] =  mu[kappa][k] - tmp; 
			      }
		      
		          fmpz_mat_row_sub(B, kappa, j, n);
		  
		       } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			      {
			         tmp = ldexp (mu[j][k], exponent);
			         mu[kappa][k] = mu[kappa][k] + tmp;
			      }
		      
                  fmpz_mat_row_add(B, kappa, j, n); 
               }
		    } else   /* we must have |X| >= 2 */
		    {
		       tmp = ldexp (mu[kappa][j] , -exponent);

	           if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		       {
		          tmp = rint (tmp);
		      
		          for (k=zeros+1; k<j; k++)
			      {
			         rtmp = tmp * mu[j][k];
			         rtmp = ldexp (rtmp, exponent);
			         mu[kappa][k] = mu[kappa][k] - rtmp;
			      }

		          xx = (signed long int) tmp;
		      
                  if (xx > 0L)
                  { 
                     fmpz_mat_row_submul_ui(B, kappa, j, (ulong) xx, n);  
                  } else
                  {
                     fmpz_mat_row_addmul_ui(B, kappa, j, (ulong) -xx, n);  
                  } 
               } else
		       {
		          tmp = frexp(mu[kappa][j], &exponent); 
		          tmp = tmp * MAX_LONG;
		          xx = (signed long int) tmp;
		          exponent += expo[kappa]-expo[j] - CPU_SIZE_1 ;

		          /* This case is extremely rare: never happened for me */
		          if (exponent <= 0) 
			      {
			         xx = xx << -exponent;
			         exponent = 0;
			  
			         if (xx > 0)
                     {
                        fmpz_mat_row_submul_ui(B, kappa, j, xx, n);  
                     } else
                     {
                        fmpz_mat_row_addmul_ui(B, kappa, j, -xx, n);  
                     }
              			    
			         for (k=zeros+1; k<j; k++)
			         {
                        rtmp = ((double) xx) * mu[j][k];
			            rtmp = ldexp (rtmp, expo[j]-expo[kappa]);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }
			      } else
			      {
			         if (xx > 0)
                     {
                        fmpz_mat_row_submul_2exp_ui(B, kappa, j, (ulong) xx, exponent, n, ztmp);  
                     } else
                     {
                        fmpz_mat_row_addmul_2exp_ui(B, kappa, j, (ulong) -xx, exponent, n, ztmp);  
                     }
			         for (k=zeros+1; k<j; k++)
			         {
			            rtmp = ((double) xx) * mu[j][k];
			            rtmp = ldexp (rtmp, exponent+expo[j]-expo[kappa]);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }
			      }
               }		  
		    }
	     }
      }
      
      if (test)   /* Anything happened? */
	  {
	     expo[kappa] = set_line(appB[kappa], B->row_arr[kappa], n);
	     aa = zeros+1;
	     for (i=zeros+1; i<=kappa; i++) 
	        appSP[kappa][i] = NAN;//0.0/0.0;
	     for (i=kappa+1; i<=kappamax; i++) 
	        appSP[i][kappa] = NAN;//0.0/0.0;
	  }
   } while (test);



   if (appSP[kappa][kappa]!=appSP[kappa][kappa]) 
   {
      appSP[kappa][kappa] = fpNorm (appB[kappa], n);
   }
   s[zeros+1] = appSP[kappa][kappa];
  
   for (k=zeros+1; k<kappa-1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }
}

/* ****************** */
/* The LLL Algorithm  */
/* ****************** */

/* LLL-reduces the integer matrix B "in place" */

void LLL (fmpz_mat_t B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double **mu, **r, **appB, **appSP;
   double *s, *mutmp, *appBtmp, *appSPtmp;
   double tmp=0.0;
   int *expo, *alpha;
   fmpz_t * ztmp = fmpz_init();
   fmpz_t * Btmp;
   
   n = B->cols;
   d = B->rows;
   ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = (double **) malloc (d * sizeof (double*));
   for (i=0; i<d; i++)
      mu[i] = (double*) malloc (d * sizeof (double));

   r = (double **) malloc (d * sizeof (double*));
   for (i=0; i<d; i++)
      r[i] = (double*) malloc (d * sizeof (double));

   appB = (double **) malloc (d * sizeof (double*));
   for (i=0; i<d; i++)
      appB[i] = (double*) malloc (n * sizeof (double));

   appSP = (double **) malloc (d * sizeof (double*));
   for (i=0; i<d; i++)
      appSP[i] = (double*) malloc (d * sizeof (double));

   s = (double *) malloc (d * sizeof (double));
   appSPtmp = (double *) malloc (d * sizeof (double));

   for (i=0; i<d; i++)
      for (j=0; j<d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i=0; i<d; i++)
      expo[i] = set_line(appB[i], B->row_arr[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = fpNorm (appB[i], n); 
   while ((appSP[i][i] <=0.0) && (++i < d));
   zeros = i-1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i+1;
  
   if (zeros < d-1)  r[i][i] = appSP[i][i];

   for (i=zeros+1; i<d; i++)
      alpha[i]=0;
    
   while (kappa < d)
   {      
      if (kappa>kappamax) kappamax++;

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai (kappa, B, mu, r, s, appB, expo, appSP, 
                                      alpha[kappa], zeros, kappamax, MPIR_MIN(kappamax+1+shift, n), ztmp); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1]-expo[kappa]));

      if (tmp <= s[kappa-1]) 
	  {
	     alpha[kappa] = kappa;
	     tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	     r[kappa][kappa] = s[kappa-1] - tmp;
	     kappa++;
	  } else
	  {

	     /* ******************************************* */
	     /* Step5: Find the right insertion index kappa */
         /* kappa2 remains the initial kappa            */
	     /* ******************************************* */  


	     kappa2 = kappa;
	     do
	     {
	        kappa--;
	        if (kappa>zeros+1) 
		    {
		       tmp = r[kappa-1][kappa-1] * ctt;
	           tmp = ldexp (tmp, 2*(expo[kappa-1]-expo[kappa2]));
	        }
         } while ( (kappa >= zeros+2) && (s[kappa-1] <= tmp));

         for (i=kappa; i<kappa2; i++)
	        if (kappa <= alpha[i]) alpha[i] = kappa;

	     for (i=kappa2; i>kappa; i--) alpha[i] = alpha[i-1];

	     for (i=kappa2+1; i<=kappamax; i++)
	        if (kappa < alpha[i]) alpha[i] = kappa;
	  
	     alpha[kappa] = kappa;

	     /* ****************************** */
	     /* Step6: Update the mu's and r's */
	     /* ****************************** */  
	  
	     mutmp = mu[kappa2];
	     for (i=kappa2; i>kappa; i--) mu[i] = mu[i-1];
	     mu[kappa] = mutmp;
	  
	     mutmp = r[kappa2];
	     for (i=kappa2; i>kappa; i--) r[i] = r[i-1];
	     r[kappa] = mutmp;

	     r[kappa][kappa] = s[kappa];
	  
	     /* ************************ */
	     /* Step7: Update B and appB */
	     /* ************************ */  	  
	  
	     Btmp = B->row_arr[kappa2];
         for (i=kappa2; i>kappa; i--) B->row_arr[i] = B->row_arr[i-1];
         B->row_arr[kappa] = Btmp;
      
	     appBtmp = appB[kappa2];
	     for (i=kappa2; i>kappa; i--) appB[i] = appB[i-1];
	     appB[kappa] = appBtmp;

	     j = expo[kappa2];
	     for (i=kappa2; i>kappa; i--) expo[i] = expo[i-1];
	     expo[kappa] = j;

	     /* *************************** */
	     /* Step8: Update appSP: tricky */
	     /* *************************** */  	 
	  
	     for (i=0; i<=kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	     for (i=kappa2+1; i<=kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	     for (i=kappa2; i>kappa; i--)
	     {
	        for (j=0; j<kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	        appSP[i][kappa] = appSPtmp[i-1];
	      
	        for (j=kappa+1; j<=i; j++) appSP[i][j] = appSP[i-1][j-1];

	        for (j=kappa2+1; j<=kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	     }
	  
	     for (i=0; i<kappa; i++) appSP[kappa][i] = appSPtmp[i];
	     appSP[kappa][kappa] = appSPtmp[kappa2];

	     for (i=kappa2+1; i<=kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	     if (r[kappa][kappa] <= 0.0)
	     {
	        zeros++;
	        kappa++;
	        appSP[kappa][kappa] = fpNorm (appB[kappa], n);
	        r[kappa][kappa] = appSP[kappa][kappa];
	     }
	  
	     kappa++;
	  }
   } 
  
   free (alpha);
   free (expo);
   fmpz_clear(ztmp);
   clear_matrixf (mu, d);
   clear_matrixf (r, d);
   clear_matrixf (appB, d);
   clear_matrixf (appSP, d);
   free (s);
   free (appSPtmp);
}
