/*
   Copyright 2005, 2006 Damien Stehlé.
   Copyright 2008  William Hart 
   Copyright 2009, William Hart, Andy Novocin

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

   F_mpz_LLL_heuristic_mpfr.c: Lattice reduction of multiprecision integer matrices using mpfrs
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
#include <mpfr.h>
#include "flint.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL_fast_d.h"
#include "F_mpz_LLL_heuristic_mpfr.h"
#include "mpfr_mat.h"
#include "gmp.h"

#define LOOPS_BABAI 10



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

void Babai_heuristic(int kappa, F_mpz_mat_t B, mpfr_t **mu, mpfr_t **r, mpfr_t *s, 
       mpfr_t **appB, mpfr_t **appSP, 
       int a, int zeros, int kappamax, int n, mpfr_t tmp, mpfr_t rtmp)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   F_mpz_t ztmp, X;
   F_mpz_init(ztmp);
   F_mpz_init(X);
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;

   long loops = 0;

   do
   {
      test = 0;
      loops++;
      if (loops > 2000)
         abort();

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if ( mpfr_nan_p(appSP[kappa][j]) ) // if appSP[kappa][j] == NAN
	      {
            if (!(mpfr_vec_scalar_product(appSP[kappa][j], appB[kappa], appB[j], n) ) ){
//In this case a heuristic told us that some cancelation probably happened so we just compute the scalar product at full precision
               F_mpz_mat_row_scalar_product(ztmp, B, kappa, B, j, 0, n);
               F_mpz_get_mpfr(appSP[kappa][j], ztmp);
            }
	      }
         if (j > zeros + 2)
	      {
	         mpfr_mul(tmp, mu[j][zeros+1], r[kappa][zeros+1], GMP_RNDN);
	         mpfr_sub(rtmp, appSP[kappa][j], tmp, GMP_RNDN);
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         mpfr_mul(tmp, mu[j][k], r[kappa][k], GMP_RNDN);
		         mpfr_sub(rtmp, rtmp, tmp, GMP_RNDN);
		      }
	         mpfr_mul(tmp, mu[j][j-1], r[kappa][j-1], GMP_RNDN);
	         mpfr_sub(r[kappa][j], rtmp, tmp, GMP_RNDN);
         } else if (j == zeros+2)
	      {
	         mpfr_mul(tmp, mu[j][zeros+1], r[kappa][zeros+1], GMP_RNDN);
	         mpfr_sub(r[kappa][j], appSP[kappa][j], tmp, GMP_RNDN);
	      } else mpfr_set(r[kappa][j], appSP[kappa][j], GMP_RNDN);

	      mpfr_div(mu[kappa][j], r[kappa][j], r[j][j], GMP_RNDN);
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {
	      /* test of the relaxed size-reduction condition */
         mpfr_abs(tmp, mu[kappa][j], GMP_RNDN); 
	  
	      if ( mpfr_cmp_d(tmp, halfplus) > 0) 
	      {
	         test = 1; 
	      
	         /* we consider separately the cases X = +-1 */     
	         if (mpfr_cmp_d(tmp, onedothalfplus) <= 0)   
		      {
               int sgn = mpfr_sgn(mu[kappa][j]);		  
		         if (sgn >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
                     mpfr_sub(mu[kappa][k], mu[kappa][k], mu[j][k], GMP_RNDN);
			         }
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            mpfr_add(mu[kappa][k], mu[kappa][k], mu[j][k], GMP_RNDN);
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
               mpfr_round(tmp, mu[kappa][j]);
		         for (k = zeros + 1; k < j; k++)
			      {
			         mpfr_mul(rtmp, tmp, mu[j][k], GMP_RNDN);
			         mpfr_sub(mu[kappa][k], mu[kappa][k], rtmp, GMP_RNDN);
			      }

	            if (mpfr_get_exp(tmp) < CPU_SIZE_1 - 2)  
		         {
                  /* X is stored in a long signed int */
                  xx = mpfr_get_si(tmp, GMP_RNDN);		      
                  if (xx > 0L)
                  { 
                     F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, (ulong) xx);  
                  } else
                  {
                     F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, (ulong) -xx);  
                  } 
               } else
		         {
                  exponent = F_mpz_set_mpfr_2exp(ztmp, tmp);
                  if (exponent <= 0){
                     F_mpz_div_2exp(ztmp, ztmp, -exponent);
                     F_mpz_mat_row_submul(B, kappa, B, j, 0, n, ztmp);
                  }
                  else{
                     F_mpz_mat_row_submul_2exp_F_mpz(B, kappa, B, j, 0, n, ztmp, exponent);
                  }
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      _F_mpz_vec_ldexp_mpfr(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         mpfr_set_nan(appSP[kappa][i]);//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         mpfr_set_nan(appSP[i][kappa]);//0.0/0.0;
	   }
   } while (test);




   if (mpfr_nan_p(appSP[kappa][kappa])) 
   {
      mpfr_vec_norm(appSP[kappa][kappa], appB[kappa], n);
   }

   mpfr_set(s[zeros + 1], appSP[kappa][kappa], GMP_RNDN);

   for (k = zeros + 1; k < kappa - 1; k++)
   {
      mpfr_mul( tmp, mu[kappa][k], r[kappa][k], GMP_RNDN);
      mpfr_sub( s[k+1], s[k], tmp, GMP_RNDN);
   }
   mpfr_set(r[kappa][kappa], s[kappa - 1], GMP_RNDN);

   F_mpz_clear(ztmp);
   F_mpz_clear(X);
}

/* ****************** */
/* The LLL Algorithm  */
/* ****************** */

/* LLL-reduces the integer matrix B "in place" */

void LLL_heuristic(F_mpz_mat_t B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   mpfr_t ** mu, ** r, ** appB, ** appSP;
   mpfr_t * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init(rtmp);
   mpfr_init(tmp);

   mu = mpfr_mat_init(d, d);
   r = mpfr_mat_init(d, d);
   appB = mpfr_mat_init(d, n);
   appSP = mpfr_mat_init(d, d);

   s = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));
   appSPtmp = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));

   for (i = 0; i < d+1; i++){
      mpfr_init(s[i]);
      mpfr_init(appSPtmp[i]);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i][j]);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_ldexp_mpfr(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      mpfr_vec_norm(appSP[i][i], appB[i], n); 
   while ( (mpfr_sgn(appSP[i][i]) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) mpfr_set(r[i][i], appSP[i][i], GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n),  tmp, rtmp); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1][kappa - 1], ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s[kappa -1]) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa][kappa-1], r[kappa][kappa-1], GMP_RNDN);
         mpfr_sub(r[kappa][kappa], s[kappa - 1], tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1][kappa-1], ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s[kappa-1],tmp) <= 0) );

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

	      mpfr_set(r[kappa][kappa], s[kappa], GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;


	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp[i], appSP[kappa2][i], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp[i], appSP[i][kappa2], GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i][j], appSP[i-1][j], GMP_RNDN);	      
	         mpfr_set(appSP[i][kappa], appSPtmp[i-1], GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i][j], appSP[i-1][j-1], GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j][i], appSP[j][i-1], GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa][i], appSPtmp[i], GMP_RNDN);
	      mpfr_set(appSP[kappa][kappa], appSPtmp[kappa2], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i][kappa], appSPtmp[i], GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa][kappa]) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         mpfr_vec_norm(appSP[kappa][kappa], appB[kappa], n);
	         mpfr_set(r[kappa][kappa], appSP[kappa][kappa], GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < d+1; i++){
      mpfr_clear(s[i]);
      mpfr_clear(appSPtmp[i]);
   }


   mpfr_mat_clear(mu, d, d);
   mpfr_mat_clear(r, d, d);
   mpfr_mat_clear(appB, d, n);
   mpfr_mat_clear(appSP, d, d);
   free(s);
   free(appSPtmp);
}


/* ****************** */
/* The LLL Algorithm  */
/* ****************** */

/* LLL-reduces the integer matrix B "in place" */

long LLL_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   mpfr_t ** mu, ** r, ** appB, ** appSP;
   mpfr_t * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init(rtmp);
   mpfr_init(tmp);

   mu = mpfr_mat_init(d, d);
   r = mpfr_mat_init(d, d);
   appB = mpfr_mat_init(d, n);
   appSP = mpfr_mat_init(d, d);

   s = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));
   appSPtmp = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));

   for (i = 0; i < d+1; i++){
      mpfr_init(s[i]);
      mpfr_init(appSPtmp[i]);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i][j]);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_ldexp_mpfr(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      mpfr_vec_norm(appSP[i][i], appB[i], n); 
   while ( (mpfr_sgn(appSP[i][i]) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) mpfr_set(r[i][i], appSP[i][i], GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n),  tmp, rtmp); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1][kappa - 1], ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s[kappa -1]) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa][kappa-1], r[kappa][kappa-1], GMP_RNDN);
         mpfr_sub(r[kappa][kappa], s[kappa - 1], tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1][kappa-1], ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s[kappa-1],tmp) <= 0) );

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

	      mpfr_set(r[kappa][kappa], s[kappa], GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;


	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp[i], appSP[kappa2][i], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp[i], appSP[i][kappa2], GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i][j], appSP[i-1][j], GMP_RNDN);	      
	         mpfr_set(appSP[i][kappa], appSPtmp[i-1], GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i][j], appSP[i-1][j-1], GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j][i], appSP[j][i-1], GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa][i], appSPtmp[i], GMP_RNDN);
	      mpfr_set(appSP[kappa][kappa], appSPtmp[kappa2], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i][kappa], appSPtmp[i], GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa][kappa]) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         mpfr_vec_norm(appSP[kappa][kappa], appB[kappa], n);
	         mpfr_set(r[kappa][kappa], appSP[kappa][kappa], GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 

   int ok = 1;
   long newd = d;

   F_mpz_get_mpfr(tmp, gs_B);

   for (i = d-1; (i >= 0) && (ok > 0); i--){
//tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      mpfr_set(rtmp, appSP[i][i], GMP_RNDN);
      mpfr_div_d(rtmp, rtmp, 2.0, GMP_RNDN);
//      mpfr_div_2ui(rtmp, rtmp, 1UL, GMP_RNDN);
      mpfr_printf(" gs length[%d] = %.10Rf, tmp = %.10Rf \n", i, rtmp, tmp);
      ok = mpfr_cmp(rtmp, tmp);
      if (ok > 0){
         newd--;
      }
   }
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < d+1; i++){
      mpfr_clear(s[i]);
      mpfr_clear(appSPtmp[i]);
   }


   mpfr_mat_clear(mu, d, d);
   mpfr_mat_clear(r, d, d);
   mpfr_mat_clear(appB, d, n);
   mpfr_mat_clear(appSP, d, d);
   free(s);
   free(appSPtmp);

   return newd;
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

void Babai_heuristic_2exp(int kappa, F_mpz_mat_t B, mpfr_t **mu, mpfr_t **r, mpfr_t *s, 
       mpfr_t **appB, mpfr_t **appSP, 
       int a, int zeros, int kappamax, int n, mpfr_t tmp, mpfr_t rtmp, int * cexpo)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   F_mpz_t ztmp, X;
   F_mpz_init(ztmp);
   F_mpz_init(X);
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;

   do
   {
      test = 0;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if ( mpfr_nan_p(appSP[kappa][j]) ) // if appSP[kappa][j] == NAN
	      {
            if (!(mpfr_vec_scalar_product(appSP[kappa][j], appB[kappa], appB[j], n) ) ){
//In this case a heuristic told us that some cancelation probably happened so we just compute the scalar product at full precision
               long exp = F_mpz_mat_row_scalar_product_2exp(ztmp, B, kappa, B, j, 0, n, cexpo);
               F_mpz_2exp_get_mpfr(appSP[kappa][j], ztmp, exp);
            }
	      }
         if (j > zeros + 2)
	      {
	         mpfr_mul(tmp, mu[j][zeros+1], r[kappa][zeros+1], GMP_RNDN);
	         mpfr_sub(rtmp, appSP[kappa][j], tmp, GMP_RNDN);
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         mpfr_mul(tmp, mu[j][k], r[kappa][k], GMP_RNDN);
		         mpfr_sub(rtmp, rtmp, tmp, GMP_RNDN);
		      }
	         mpfr_mul(tmp, mu[j][j-1], r[kappa][j-1], GMP_RNDN);
	         mpfr_sub(r[kappa][j], rtmp, tmp, GMP_RNDN);
         } else if (j == zeros+2)
	      {
	         mpfr_mul(tmp, mu[j][zeros+1], r[kappa][zeros+1], GMP_RNDN);
	         mpfr_sub(r[kappa][j], appSP[kappa][j], tmp, GMP_RNDN);
	      } else mpfr_set(r[kappa][j], appSP[kappa][j], GMP_RNDN);

	      mpfr_div(mu[kappa][j], r[kappa][j], r[j][j], GMP_RNDN);
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {
	      /* test of the relaxed size-reduction condition */
         mpfr_abs(tmp, mu[kappa][j], GMP_RNDN); 
	  
	      if ( mpfr_cmp_d(tmp, halfplus) > 0) 
	      {
	         test = 1; 
	      
	         /* we consider separately the cases X = +-1 */     
	         if (mpfr_cmp_d(tmp, onedothalfplus) <= 0)   
		      {
               int sgn = mpfr_sgn(mu[kappa][j]);		  
		         if (sgn >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
                     mpfr_sub(mu[kappa][k], mu[kappa][k], mu[j][k], GMP_RNDN);
			         }
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            mpfr_add(mu[kappa][k], mu[kappa][k], mu[j][k], GMP_RNDN);
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
               mpfr_round(tmp, mu[kappa][j]);
		         for (k = zeros + 1; k < j; k++)
			      {
			         mpfr_mul(rtmp, tmp, mu[j][k], GMP_RNDN);
			         mpfr_sub(mu[kappa][k], mu[kappa][k], rtmp, GMP_RNDN);
			      }

	            if (mpfr_get_exp(tmp) < CPU_SIZE_1 - 2)  
		         {
                  /* X is stored in a long signed int */
                  xx = mpfr_get_si(tmp, GMP_RNDN);		      
                  if (xx > 0L)
                  { 
                     F_mpz_mat_row_submul_ui(B, kappa, B, j, 0, n, (ulong) xx);  
                  } else
                  {
                     F_mpz_mat_row_addmul_ui(B, kappa, B, j, 0, n, (ulong) -xx);  
                  } 
               } else
		         {
                  exponent = F_mpz_set_mpfr_2exp(ztmp, tmp);
                  if (exponent <= 0){
                     F_mpz_div_2exp(ztmp, ztmp, -exponent);
                     F_mpz_mat_row_submul(B, kappa, B, j, 0, n, ztmp);
                  }
                  else{
                     F_mpz_mat_row_submul_2exp_F_mpz(B, kappa, B, j, 0, n, ztmp, exponent);
                  }
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      _F_mpz_vec_ldexp_mpfr_2exp(appB[kappa], B->rows[kappa], n, cexpo);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         mpfr_set_nan(appSP[kappa][i]);//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         mpfr_set_nan(appSP[i][kappa]);//0.0/0.0;
	   }
   } while (test);




   if (mpfr_nan_p(appSP[kappa][kappa])) 
   {
      mpfr_vec_norm(appSP[kappa][kappa], appB[kappa], n);
   }

   mpfr_set(s[zeros + 1], appSP[kappa][kappa], GMP_RNDN);

   for (k = zeros + 1; k < kappa - 1; k++)
   {
      mpfr_mul( tmp, mu[kappa][k], r[kappa][k], GMP_RNDN);
      mpfr_sub( s[k+1], s[k], tmp, GMP_RNDN);
   }
   mpfr_set(r[kappa][kappa], s[kappa - 1], GMP_RNDN);

   F_mpz_clear(ztmp);
   F_mpz_clear(X);
}

/* ****************** */
/* The LLL Algorithm  */
/* ****************** */

/* LLL-reduces the integer matrix B "in place" */

void LLL_heuristic_2exp(F_mpz_mat_t B, int * cexpo)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   mpfr_t ** mu, ** r, ** appB, ** appSP;
   mpfr_t * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init(rtmp);
   mpfr_init(tmp);

   mu = mpfr_mat_init(d, d);
   r = mpfr_mat_init(d, d);
   appB = mpfr_mat_init(d, n);
   appSP = mpfr_mat_init(d, d);

   s = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));
   appSPtmp = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));

   for (i = 0; i < d+1; i++){
      mpfr_init(s[i]);
      mpfr_init(appSPtmp[i]);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i][j]);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_ldexp_mpfr_2exp(appB[i], B->rows[i], n, cexpo);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      mpfr_vec_norm(appSP[i][i], appB[i], n); 
   while ( (mpfr_sgn(appSP[i][i]) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) mpfr_set(r[i][i], appSP[i][i], GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai_heuristic_2exp(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n),  tmp, rtmp, cexpo); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1][kappa - 1], ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s[kappa -1]) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa][kappa-1], r[kappa][kappa-1], GMP_RNDN);
         mpfr_sub(r[kappa][kappa], s[kappa - 1], tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1][kappa-1], ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s[kappa-1],tmp) <= 0) );

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

	      mpfr_set(r[kappa][kappa], s[kappa], GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;


	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp[i], appSP[kappa2][i], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp[i], appSP[i][kappa2], GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i][j], appSP[i-1][j], GMP_RNDN);	      
	         mpfr_set(appSP[i][kappa], appSPtmp[i-1], GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i][j], appSP[i-1][j-1], GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j][i], appSP[j][i-1], GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa][i], appSPtmp[i], GMP_RNDN);
	      mpfr_set(appSP[kappa][kappa], appSPtmp[kappa2], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i][kappa], appSPtmp[i], GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa][kappa]) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         mpfr_vec_norm(appSP[kappa][kappa], appB[kappa], n);
	         mpfr_set(r[kappa][kappa], appSP[kappa][kappa], GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < d+1; i++){
      mpfr_clear(s[i]);
      mpfr_clear(appSPtmp[i]);
   }


   mpfr_mat_clear(mu, d, d);
   mpfr_mat_clear(r, d, d);
   mpfr_mat_clear(appB, d, n);
   mpfr_mat_clear(appSP, d, d);
   free(s);
   free(appSPtmp);
}


/* ****************** */
/* The LLL Algorithm  */
/* ****************** */

/* LLL-reduces the integer matrix B "in place" */

long LLL_heuristic_2exp_with_removal(F_mpz_mat_t B, int * cexpo, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   mpfr_t ** mu, ** r, ** appB, ** appSP;
   mpfr_t * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init(rtmp);
   mpfr_init(tmp);

   mu = mpfr_mat_init(d, d);
   r = mpfr_mat_init(d, d);
   appB = mpfr_mat_init(d, n);
   appSP = mpfr_mat_init(d, d);

   s = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));
   appSPtmp = (mpfr_t *) malloc ((d + 1) * sizeof(mpfr_t));

   for (i = 0; i < d+1; i++){
      mpfr_init(s[i]);
      mpfr_init(appSPtmp[i]);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i][j]);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_ldexp_mpfr_2exp(appB[i], B->rows[i], n, cexpo);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      mpfr_vec_norm(appSP[i][i], appB[i], n); 
   while ( (mpfr_sgn(appSP[i][i]) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
  
   if (zeros < d - 1) mpfr_set(r[i][i], appSP[i][i], GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      Babai_heuristic_2exp(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n),  tmp, rtmp, cexpo); 
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1][kappa - 1], ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s[kappa -1]) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa][kappa-1], r[kappa][kappa-1], GMP_RNDN);
         mpfr_sub(r[kappa][kappa], s[kappa - 1], tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1][kappa-1], ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s[kappa-1],tmp) <= 0) );

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

	      mpfr_set(r[kappa][kappa], s[kappa], GMP_RNDN);
	  
	      /* ************************ */
	      /* Step7: Update B and appB */
	      /* ************************ */  	  
	  
	      Btmp = B->rows[kappa2];
         for (i = kappa2; i > kappa; i--) B->rows[i] = B->rows[i-1];
         B->rows[kappa] = Btmp;
      
	      appBtmp = appB[kappa2];
	      for (i = kappa2; i > kappa; i--) appB[i] = appB[i-1];
	      appB[kappa] = appBtmp;


	      /* *************************** */
	      /* Step8: Update appSP: tricky */
	      /* *************************** */  	 
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp[i], appSP[kappa2][i], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp[i], appSP[i][kappa2], GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i][j], appSP[i-1][j], GMP_RNDN);	      
	         mpfr_set(appSP[i][kappa], appSPtmp[i-1], GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i][j], appSP[i-1][j-1], GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j][i], appSP[j][i-1], GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa][i], appSPtmp[i], GMP_RNDN);
	      mpfr_set(appSP[kappa][kappa], appSPtmp[kappa2], GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i][kappa], appSPtmp[i], GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa][kappa]) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         mpfr_vec_norm(appSP[kappa][kappa], appB[kappa], n);
	         mpfr_set(r[kappa][kappa], appSP[kappa][kappa], GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 
 
   int ok = 1;
   int gok = 1;
   int gap = 1;
   long newd = d;

   F_mpz_get_mpfr(tmp, gs_B);

   for (i = d-1; (i >= 0) && (ok > 0); i--){
//tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      mpfr_set(rtmp, appSP[i][i], GMP_RNDN);
      mpfr_div_d(rtmp, rtmp, 2.0, GMP_RNDN);
//mpfr_div_2ui(rtmp, rtmp, 1UL, GMP_RNDN);
      mpfr_printf(" gs length[%d] = %.10Rf, tmp = %.10Rf \n", i, rtmp, tmp);
      ok = mpfr_cmp(rtmp, tmp);
      if (ok > 0){
         newd--;
      }
   }
   
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < d+1; i++){
      mpfr_clear(s[i]);
      mpfr_clear(appSPtmp[i]);
   }


   mpfr_mat_clear(mu, d, d);
   mpfr_mat_clear(r, d, d);
   mpfr_mat_clear(appB, d, n);
   mpfr_mat_clear(appSP, d, d);
   free(s);
   free(appSPtmp);

   return newd;
}


