/*
   Copyright 2005, 2006 Damien Stehlé.
   Copyright 2008 William Hart

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

/****************************************************************************

   F_mpz_LLL_fast_d.c: Lattice reduction of multiprecision integer matrices using doubles
	                    for storing approximate Gram Schmidt Orthogonalisations

*****************************************************************************/

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
#include "flint.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL_helper.h"
#include "F_mpz_LLL_fast_d.h"
#include "d_mat.h"

#define LOOPS_BABAI 10

/* Computes the largest number of non-zero entries after the diagonal. */

ulong getShift(F_mpz_mat_t B)
{
   ulong n = B->c;
   ulong shift = 0;
   ulong i;
   for (i = 0; i < B->r; i++)
   {
      ulong j;
      for (j = n - 1; j >= i + shift + 1 && F_mpz_size(B->rows[i] + j) == 0L; j--);  
      
      if (shift < j - i) shift = j - i;
      
   }

   return shift;
}

/***********************************/
/* Babai's Nearest Plane algorithm */
/***********************************/

/* 
   Computes mu[kappa][j] and r[kappa][j] for j < kappa.
	
	Size-reduces b_kappa using r[i][j] for j <= i < kappa
   and mu[i][j] for j < i < kappa.
	
	Compute s(kappa).
   
	Updates B(kappa).
   
	The algorithm is the iterative Babai algorithm of the paper.
*/

void Babai (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;


   do
   {
      test = 0;
      
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
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
	         appSP[kappa][j] = _d_vec_scalar_product(appB[kappa], appB[j], n);
	      }
	  	  
         if (j > zeros + 2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         rtmp = appSP[kappa][j] - tmp;
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		      }
	         tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
         } else if (j == zeros+2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		      {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
		      
		            F_mpz_mat_row_sub(B, kappa, B, kappa, B, j, 0, n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                  F_mpz_mat_row_add(B, kappa, B, kappa, B, j, 0, n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
		         tmp = ldexp (mu[kappa][j] , -exponent);

	            if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		            tmp = rint (tmp); 
		      
		            for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
			            rtmp = ldexp (rtmp, exponent);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		            xx = (long) tmp;
		      
                  if (xx > 0L)
                  { 
                     F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, (ulong) xx);  
                  } else
                  {
                     F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, (ulong) -xx);  
                  } 
               } else
		         {
		            tmp = frexp(mu[kappa][j], &exponent); 
		            tmp = tmp * MAX_LONG;
		            xx = (long) tmp;
		            exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		            /* This case is extremely rare: never happened for me */
		            if (exponent <= 0) 
			         {
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, xx);  
                     } else
                     {
                        F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, -xx);  
                     }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                        rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_2exp_ui(B, kappa, B, j, 0, n, (ulong) xx, exponent);  
                     } else
                     {
                        F_mpz_mat_row_addmul_2exp_ui(B, kappa, B, j, 0, n, (ulong) -xx, exponent);  
                     }
			            for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					      }
				      }	    
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_ldexp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
   }
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }
}

/***********************************/
/* Babai's Nearest Plane algorithm */
/***********************************/

/* 
   Computes mu[kappa][j] and r[kappa][j] for j < kappa.
	
	Size-reduces b_kappa using r[i][j] for j <= i < kappa
   and mu[i][j] for j < i < kappa.
	
	Compute s(kappa).
   
	Updates B(kappa).
   
	The algorithm is the iterative Babai algorithm of the paper.
*/

//### This is different -------
void Babai_heuristic_d (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
//-----------------------------
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;


   do
   {
      test = 0;
      
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
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
//### This is different -----
            appSP[kappa][j] = heuristic_scalar_product(appB[kappa], appB[j], n, B, kappa, j, expo[kappa]+expo[j]);
//---------------------------
         }
	  	  
         if (j > zeros + 2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         rtmp = appSP[kappa][j] - tmp;
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		      }
	         tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
         } else if (j == zeros+2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		      {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
		      
		            F_mpz_mat_row_sub(B, kappa, B, kappa, B, j, 0, n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                  F_mpz_mat_row_add(B, kappa, B, kappa, B, j, 0, n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
		         tmp = ldexp (mu[kappa][j] , -exponent);

	            if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		            tmp = rint (tmp); 
		      
		            for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
			            rtmp = ldexp (rtmp, exponent);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		            xx = (long) tmp;
		      
                  if (xx > 0L)
                  { 
                     F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, (ulong) xx);  
                  } else
                  {
                     F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, (ulong) -xx);  
                  } 
               } else
		         {
		            tmp = frexp(mu[kappa][j], &exponent); 
		            tmp = tmp * MAX_LONG;
		            xx = (long) tmp;
		            exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		            /* This case is extremely rare: never happened for me */
		            if (exponent <= 0) 
			         {
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, xx);  
                     } else
                     {
                        F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, -xx);  
                     }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                        rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_2exp_ui(B, kappa, B, j, 0, n, (ulong) xx, exponent);  
                     } else
                     {
                        F_mpz_mat_row_addmul_2exp_ui(B, kappa, B, j, 0, n, (ulong) -xx, exponent);  
                     }
			            for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					      }
				      }	    
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_ldexp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
//### This is different -------
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
//-----------------------------
   }
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }
}

#if 0
//### This is different -------
void Babai_heuristic_d_2exp (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n, int *cexpo)
//-----------------------------
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;


   do
   {
      test = 0;
      
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
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
//### This is different -----
            appSP[kappa][j] = d_2exp_vec_scalar_product(appB[kappa], appB[j], n, cexpo, B, kappa, j);
//---------------------------
         }
	  	  
         if (j > zeros + 2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         rtmp = appSP[kappa][j] - tmp;
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         tmp = mu[j][k] * r[kappa][k];
		         rtmp = rtmp - tmp;
		      }
	         tmp = mu[j][j-1] * r[kappa][j-1];
	         r[kappa][j] = rtmp - tmp;
         } else if (j == zeros+2)
	      {
	         tmp = mu[j][zeros+1] * r[kappa][zeros+1];
	         r[kappa][j] = appSP[kappa][j] - tmp;
	      } else r[kappa][j] = appSP[kappa][j];

	      mu[kappa][j] = r[kappa][j] / r[j][j];
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
	  
	      if (tmp > halfplus) 
	      {
	         test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		      {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
		      
		            F_mpz_mat_row_sub(B, kappa, B, kappa, B, j, 0, n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                  F_mpz_mat_row_add(B, kappa, B, kappa, B, j, 0, n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
		         tmp = ldexp (mu[kappa][j] , -exponent);

	            if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		            tmp = rint (tmp); 
		      
		            for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
			            rtmp = ldexp (rtmp, exponent);
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		            xx = (long) tmp;
		      
                  if (xx > 0L)
                  { 
                     F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, (ulong) xx);  
                  } else
                  {
                     F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, (ulong) -xx);  
                  } 
               } else
		         {
		            tmp = frexp(mu[kappa][j], &exponent); 
		            tmp = tmp * MAX_LONG;
		            xx = (long) tmp;
		            exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		            /* This case is extremely rare: never happened for me */
		            if (exponent <= 0) 
			         {
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, xx);  
                     } else
                     {
                        F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, -xx);  
                     }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                        rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_2exp_ui(B, kappa, B, j, 0, n, (ulong) xx, exponent);  
                     } else
                     {
                        F_mpz_mat_row_addmul_2exp_ui(B, kappa, B, j, 0, n, (ulong) -xx, exponent);  
                     }
			            for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					      }
				      }	    
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_ldexp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
//### This is different -------
      appSP[kappa][kappa] = d_2exp_vec_norm(appB[kappa], n, cexpo);
//-----------------------------
   }
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }
}
#endif

/* ****************** */
/* The LLL Algorithm  */
/* ****************** */

/* LLL-reduces the integer matrix B "in place" */

void LLL(F_mpz_mat_t B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_ldexp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n)); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

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
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
         } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

         for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
	         r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);
}

/* LLL-reduces the integer matrix B "in place" 
   uses an array of virtual weights for each column (cexpo), allows powers of 2
   so cexpo = [0,...,0] is normal LLL 
   and cexpo = [2,0,...,0] will weigh the first column as 4 times the importance of the others
*/

#if 0
//### This is different ------
void LLL_heuristic_d_2exp(F_mpz_mat_t B, int *cexpo)
//----------------------------
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_ldexp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
//### This is different -----
      appSP[i][i] = d_2exp_vec_norm(appB[i], n, cexpo); 
//---------------------------
      while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

//### This is different -----
      Babai_heuristic_d_2exp(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n), cexpo); 
//---------------------------

      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

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
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
         } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

         for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)
	      {
	         zeros++;
	         kappa++;
//### This is different ------
            appSP[kappa][kappa] = d_2exp_vec_norm(appB[kappa], n, cexpo);
//----------------------------
            r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);
}
#endif

/* 
   LLL-reduces the integer matrix B "in place"
   uses a virtual weight for each column stored as a power of 2 in the array cexpo (0,...,0) would be normal LLL
   also returns the number of rows who's G-S lengths are guaranteed to be <= gs_B 
*/

#if 0
//### This is different ------
int LLL_heuristic_d_2exp_with_removal(F_mpz_mat_t B, int *cexpo, F_mpz_t gs_B)
//----------------------------
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_ldexp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = d_2exp_vec_norm(appB[i], n, cexpo); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai_heuristic_d_2exp(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n), cexpo); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

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
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
         } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

         for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)//fixme should be EPS
	      {
	         zeros++;
	         kappa++;
	         appSP[kappa][kappa] = d_2exp_vec_norm(appB[kappa], n, cexpo);
	         r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }

//### This is different --------
   F_mpz_t tmp_gs;
   F_mpz_init(tmp_gs);
 
   int ok = 1;
   int newd = d;
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
      //tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      F_mpz_set_d_2exp(tmp_gs, appSP[i][i], 2*expo[i] - 1);
      ok = F_mpz_cmpabs(tmp_gs, gs_B);
      if (ok > 0) newd--;
   }

   F_mpz_clear(tmp_gs);
//-------------------------------

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);
   return newd;
}
#endif

/* 
   LLL-reduces the integer matrix B "in place"
   uses a virtual weight for each column stored as a power of 2 in the array cexpo. (0,...,0) would be normal LLL
   also returns the number of rows who's G-S lengths are guaranteed to be <= gs_B 
*/

//### This is different ------
int LLL_heuristic_d_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
//----------------------------
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_ldexp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
//### This is different -----
   appSP[i][i] = _d_vec_norm(appB[i], n); 
//---------------------------
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n)); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */
      
      tmp = r[kappa-1][kappa-1] * ctt;
      tmp = ldexp (tmp, 2*(expo[kappa-1] - expo[kappa]));

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
	         if (kappa > zeros + 1) 
		      {
		         tmp = r[kappa-1][kappa-1] * ctt;
	            tmp = ldexp(tmp, 2*(expo[kappa-1] - expo[kappa2]));
	         }
         } while ((kappa >= zeros + 2) && (s[kappa-1] <= tmp));

         for (i = kappa; i < kappa2; i++)
	         if (kappa <= alpha[i]) alpha[i] = kappa;

	      for (i = kappa2; i > kappa; i--) alpha[i] = alpha[i-1];

	      for (i = kappa2 + 1; i <= kappamax; i++)
	         if (kappa < alpha[i]) alpha[i] = kappa;
	  
	      alpha[kappa] = kappa;

	      /* ****************************** */
	      /* Step6: Update the mu's and r's */
	      /* ****************************** */  
	  
	      mutmp = mu[kappa2];
	      for (i = kappa2; i > kappa; i--) mu[i] = mu[i-1];
	      mu[kappa] = mutmp;
	  
	      mutmp = r[kappa2];
	      for (i = kappa2; i > kappa; i--) r[i] = r[i-1];
	      r[kappa] = mutmp;

	      r[kappa][kappa] = s[kappa];
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;

	      j = expo[kappa2];
	      for (i = kappa2; i > kappa; i--) expo[i] = expo[i-1];
	      expo[kappa] = j;

	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) appSPtmp[i] = appSP[kappa2][i];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSPtmp[i] = appSP[i][kappa2];
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) appSP[i][j] = appSP[i-1][j];	      
	         appSP[i][kappa] = appSPtmp[i-1];
	      
	         for (j = kappa + 1; j <= i; j++) appSP[i][j] = appSP[i-1][j-1];

	         for (j = kappa2 + 1; j <= kappamax; j++) appSP[j][i] = appSP[j][i-1];     
	      }
	  
	      for (i = 0; i < kappa; i++) appSP[kappa][i] = appSPtmp[i];
	      appSP[kappa][kappa] = appSPtmp[kappa2];

	      for (i = kappa2 + 1; i <= kappamax; i++) appSP[i][kappa] = appSPtmp[i];
	  
	      if (r[kappa][kappa] <= 0.0)//fixme should be EPS
	      {
	         zeros++;
	         kappa++;
//### This is different -----
	         appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
//---------------------------
            r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }

//### This is different --------
   F_mpz_t tmp_gs;
   F_mpz_init(tmp_gs);
 
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
//tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
//F_mpz_set_d_2exp(tmp_gs, r[i][i], 2*expo[i] - 1);
      d_rii = ldexp(r[i][i], 2*expo[i] - 1);
      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newd--;
   }

/*   for (i = d-1; (i >= 0); i--)
   {
      //tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      F_mpz_set_d_2exp(tmp_gs, r[i][i], 2*expo[i] - 1);
      F_mpz_print(tmp_gs); printf(" is r[i][i] for %d\n", i);
   }*/

   F_mpz_clear(tmp_gs);
//-------------------------------

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);
   return newd;
}
