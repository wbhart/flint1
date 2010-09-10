/*
   Copyright 2009 William Hart, Andy Novocin

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
/****************************************************************************

   F_mpz_LLL_wrapper.c: A program for navigating the various versions of LLL

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <ctype.h>
#include <limits.h>
#include "gmp.h"
#include "flint.h"
#include "profiler.h"
#include "F_mpz_mat.h"
#include "F_mpz_LLL_helper.h"
#include "F_mpz_LLL_wrapper.h"
#include "F_mpz_LLL_fast_d.h"
#include "F_mpz_LLL_heuristic_mpfr.h"
#include "mpfr.h"
#include "mpfr_mat.h"
#include "d_mat.h"


/****************************************************************************

   The various Babai's: check_Babai, check_Babai_heuristic_d, check_Babai_heuristic

****************************************************************************/

#define PROFILE 1

#if PROFILE
//LLL profiles are on, checking babai total, time spent updating B-- but not including the time in advanced babai which is totaled
   double babai_start, babai_stop;
   double babai_total = 0;

   double adv_babai_start, adv_babai_stop;
   double adv_babai_total = 0;

   double update_start, update_stop;
   double update_total = 0.0;

   double convert_start, convert_stop;
   double convert_total = 0.0;

   double inner_start, inner_stop;
   double inner_total = 0.0;

   double ldexp_start, ldexp_stop;
   double ldexp_total = 0.0;

   double hbabai_start, hbabai_stop;
   double hbabai_total = 0;

   double hadv_babai_start, hadv_babai_stop;
   double hadv_babai_total = 0;

   double hupdate_start, hupdate_stop;
   double hupdate_total = 0.0;

   double hconvert_start, hconvert_stop;
   double hconvert_total = 0.0;

   double hinner_start, hinner_stop;
   double hinner_total = 0.0;

   double hldexp_start, hldexp_stop;
   double hldexp_total = 0.0;

#endif


int check_Babai (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
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

   long loops = 0;

   do
   {
      test = 0;

      loops++;
      if (loops > 2){
         return -1;
      }
            
      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
#if PROFILE
   inner_start = get_cycle_counter();
#endif
	         appSP[kappa][j] = _d_vec_scalar_product(appB[kappa], appB[j], n);
#if PROFILE
   inner_stop = get_cycle_counter();
   inner_total = inner_total + inner_stop - inner_start;
#endif
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
#if PROFILE
   ldexp_start = get_cycle_counter();
#endif
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
#if PROFILE
   ldexp_stop = get_cycle_counter();
   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif	  
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
#if PROFILE
   update_start = get_cycle_counter();
#endif
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
#if PROFILE
   update_start = get_cycle_counter();
#endif
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
               }
		      } else   /* we must have |X| >= 2 */
		      {
#if PROFILE
   ldexp_start = get_cycle_counter();
#endif
		         tmp = ldexp (mu[kappa][j] , -exponent);
#if PROFILE
   ldexp_stop = get_cycle_counter();
   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
	            if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		            tmp = rint (tmp); 
		      
		            for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
#if PROFILE
   ldexp_start = get_cycle_counter();
#endif
			            rtmp = ldexp (rtmp, exponent);
#if PROFILE
   ldexp_stop = get_cycle_counter();
   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		            xx = (long) tmp;
		      
                  if (xx > 0L)
                  { 
#if PROFILE
   update_start = get_cycle_counter();
#endif
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
                  } else
                  {
#if PROFILE
   update_start = get_cycle_counter();
#endif
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
                  } 
               } else
		         {
#if PROFILE
   ldexp_start = get_cycle_counter();
#endif
		            tmp = frexp(mu[kappa][j], &exponent); 
#if PROFILE
   ldexp_stop = get_cycle_counter();
   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
		            tmp = tmp * MAX_LONG;
		            xx = (long) tmp;
		            exponent += expo[kappa] - expo[j] - CPU_SIZE_1;

		            /* This case is extremely rare: never happened for me */
		            if (exponent <= 0) 
			         {
//                     printf("rare case kappa = %d, j = %d ******************************************************\n", kappa, j);
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                     {
#if PROFILE
   update_start = get_cycle_counter();
#endif
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
                     } else
                     {
#if PROFILE
   update_start = get_cycle_counter();
#endif
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
                     }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                        rtmp = ((double) xx) * mu[j][k];
#if PROFILE
   ldexp_start = get_cycle_counter();
#endif
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
#if PROFILE
   ldexp_stop = get_cycle_counter();
   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                     {
#if PROFILE
   update_start = get_cycle_counter();
#endif
                        _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
                     } else
                     {
#if PROFILE
   update_start = get_cycle_counter();
#endif
                        _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
#if PROFILE
   update_stop = get_cycle_counter();
   update_total = update_total + update_stop - update_start;
#endif
                     }
			            for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
#if PROFILE
   ldexp_start = get_cycle_counter();
#endif
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
#if PROFILE
   ldexp_stop = get_cycle_counter();
   ldexp_total = ldexp_total + ldexp_stop - ldexp_start;
#endif
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					      }
				      }	    
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
#if PROFILE
   convert_start = get_cycle_counter();
#endif
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
#if PROFILE
   convert_stop = get_cycle_counter();
   convert_total = convert_total + convert_stop - convert_start;
#endif
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
#if PROFILE
   inner_start = get_cycle_counter();
#endif
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
#if PROFILE
   inner_stop = get_cycle_counter();
   inner_total = inner_total + inner_stop - inner_start;
#endif
   }
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }

   return 0;
}

int check_Babai_heuristic_d (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
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


   int loops = 0;

   do
   {
      test = 0;
            
      loops++;
      if (loops > 10)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
#if PROFILE
   hinner_start = get_cycle_counter();
#endif
//### This is different -----
            appSP[kappa][j] = heuristic_scalar_product(appB[kappa], appB[j], n, B, kappa, j, expo[kappa]+expo[j]);
//---------------------------
#if PROFILE
   hinner_stop = get_cycle_counter();
   hinner_total = hinner_total + hinner_stop - hinner_start;
#endif
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
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
	  
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
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
			            tmp = ldexp (mu[j][k], exponent);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
			            tmp = ldexp (mu[j][k], exponent);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }

#if PROFILE
   hupdate_start = get_cycle_counter();
#endif		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
               }
		      } else   /* we must have |X| >= 2 */
		      {
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
		         tmp = ldexp (mu[kappa][j] , -exponent);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif

	            if ((tmp < (double) MAX_LONG) &&(tmp > (double) -MAX_LONG))  
		         {
		            tmp = rint (tmp); 
		      
		            for (k = zeros + 1; k < j; k++)
			         {
			            rtmp = tmp * mu[j][k];
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
			            rtmp = ldexp (rtmp, exponent);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			            mu[kappa][k] = mu[kappa][k] - rtmp;
			         }

		            xx = (long) tmp;
		      
                  if (xx > 0L)
                  { 
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                  } else
                  {
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                  } 
               } else
		         {
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
		            tmp = frexp(mu[kappa][j], &exponent); 
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
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
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                     } else
                     {
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                     }
              			    
			            for (k = zeros + 1; k < j; k++)
			            {
                        rtmp = ((double) xx) * mu[j][k];
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
			               rtmp = ldexp (rtmp, expo[j] - expo[kappa]);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			               mu[kappa][k] = mu[kappa][k] - rtmp;
			            }
			         } else
			         {
			            if (xx > 0)
                     {
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif
                        _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                     } else
                     {
#if PROFILE
   hupdate_start = get_cycle_counter();
#endif
                        _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
#if PROFILE
   hupdate_stop = get_cycle_counter();
   hupdate_total = hupdate_total + hupdate_stop - hupdate_start;
#endif
                     }
			            for (k = zeros + 1; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
#if PROFILE
   hldexp_start = get_cycle_counter();
#endif
			               rtmp = ldexp (rtmp, exponent + expo[j] - expo[kappa]);
#if PROFILE
   hldexp_stop = get_cycle_counter();
   hldexp_total = hldexp_total + hldexp_stop - hldexp_start;
#endif
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					      }
				      }	    
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
#if PROFILE
   hconvert_start = get_cycle_counter();
#endif
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
#if PROFILE
   hconvert_stop = get_cycle_counter();
   hconvert_total = hconvert_total + hconvert_stop - hconvert_start;
#endif
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
#if PROFILE
   hinner_start = get_cycle_counter();
#endif
      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
#if PROFILE
   hinner_stop = get_cycle_counter();
   hinner_total = hinner_total + hinner_stop - hinner_start;
#endif
//-----------------------------
   }
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }

   return 0;
}

int check_Babai_heuristic(int kappa, F_mpz_mat_t B, __mpfr_struct **mu, __mpfr_struct **r, __mpfr_struct *s, 
       __mpfr_struct **appB, __mpfr_struct **appSP, 
       int a, int zeros, int kappamax, int n, mpfr_t tmp, mpfr_t rtmp, mp_prec_t prec)
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
      if (loops > 3)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < kappa; j++)
	   {	  
	      if ( mpfr_nan_p(appSP[kappa] + j) ) // if appSP[kappa] + j == NAN
	      {
            if (!(_mpfr_vec_scalar_product2(appSP[kappa] + j, appB[kappa], appB[j], n, prec) ) ){
//In this case a heuristic told us that some cancelation probably happened so we just compute the scalar product at full precision
               _F_mpz_vec_scalar_product(ztmp, B->rows[kappa], B->rows[j], n);
               F_mpz_get_mpfr(appSP[kappa] + j, ztmp);
            }
	      }
         if (j > zeros + 2)
	      {
	         mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	         mpfr_sub(rtmp, appSP[kappa] + j, tmp, GMP_RNDN);
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         mpfr_mul(tmp, mu[j] + k, r[kappa] + k, GMP_RNDN);
		         mpfr_sub(rtmp, rtmp, tmp, GMP_RNDN);
		      }
	         mpfr_mul(tmp, mu[j] + j - 1, r[kappa] + j - 1, GMP_RNDN);
	         mpfr_sub(r[kappa] + j, rtmp, tmp, GMP_RNDN);
         } else if (j == zeros+2)
	      {
	         mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	         mpfr_sub(r[kappa] + j, appSP[kappa] + j, tmp, GMP_RNDN);
	      } else mpfr_set(r[kappa] + j, appSP[kappa] + j, GMP_RNDN);

	      mpfr_div(mu[kappa] + j, r[kappa] + j, r[j] + j, GMP_RNDN);
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {
	      /* test of the relaxed size-reduction condition */
         mpfr_abs(tmp, mu[kappa] + j, GMP_RNDN); 
	  
	      if ( mpfr_cmp_d(tmp, halfplus) > 0) 
	      {
	         test = 1; 
	      
	         /* we consider separately the cases X = +-1 */     
	         if (mpfr_cmp_d(tmp, onedothalfplus) <= 0)   
		      {
               int sgn = mpfr_sgn(mu[kappa] + j);		  
		         if (sgn >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
                     mpfr_sub(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			         }
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            mpfr_add(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
               mpfr_round(tmp, mu[kappa] + j);
		         for (k = zeros + 1; k < j; k++)
			      {
			         mpfr_mul(rtmp, tmp, mu[j] + k, GMP_RNDN);
			         mpfr_sub(mu[kappa] + k, mu[kappa] + k, rtmp, GMP_RNDN);
			      }

	            if (mpfr_get_exp(tmp) < CPU_SIZE_1 - 2)  
		         {
                  /* X is stored in a long signed int */
                  xx = mpfr_get_si(tmp, GMP_RNDN);		      
                  if (xx > 0L)
                  { 
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                  } else
                  {
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                  } 
               } else
		         {
                  exponent = F_mpz_set_mpfr_2exp(ztmp, tmp);
                  if (exponent <= 0){
                     F_mpz_div_2exp(ztmp, ztmp, -exponent);
                     _F_mpz_vec_submul_F_mpz(B->rows[kappa], B->rows[j], n, ztmp);
                  }
                  else{
                     _F_mpz_vec_submul_2exp_F_mpz(B->rows[kappa], B->rows[j], n, ztmp, exponent);
                  }
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      _F_mpz_vec_to_mpfr_vec(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         mpfr_set_nan(appSP[kappa] + i);//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         mpfr_set_nan(appSP[i] + kappa);//0.0/0.0;
	   }
   } while (test);




   if (mpfr_nan_p(appSP[kappa] + kappa)) 
   {
      _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
   }

   mpfr_set(s + zeros + 1, appSP[kappa] + kappa, GMP_RNDN);

   for (k = zeros + 1; k < kappa - 1; k++)
   {
      mpfr_mul( tmp, mu[kappa] + k, r[kappa] + k, GMP_RNDN);
      mpfr_sub( s + k + 1, s + k, tmp, GMP_RNDN);
   }
   mpfr_set(r[kappa] + kappa, s + kappa - 1, GMP_RNDN);

   F_mpz_clear(ztmp);
   F_mpz_clear(X);

   return 0;
}

int advance_check_Babai (int cur_kappa, int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   ulong temp_expo;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;

   long loops = 0;

   do
   {
      test = 0;

      loops++;
      if (loops > 5){
         return -1;
      }
            
      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < cur_kappa; j++)
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
      for (j = cur_kappa - 1; j > zeros; j--)
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
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
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
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                  } else
                  {
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
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
//                     printf("rare case kappa = %d, j = %d ******************************************************\n", kappa, j);
			            xx = xx << -exponent;
			            exponent = 0;
			  
			            if (xx > 0)
                     {
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
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
                        _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
                     } else
                     {
                        _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
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

      if (test == 1)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= cur_kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = cur_kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
      else
      {
         for (i = zeros + 1; i <= cur_kappa; i++)
            appSP[kappa][i] = NAN;
      }
   } while (test == 1);

//   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
//   {
//      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
//   }
//   s[zeros + 1] = appSP[kappa][kappa];
  
//   for (k = zeros + 1; k < kappa - 1; k++)
//   {
//      tmp = mu[kappa][k] * r[kappa][k];
//      s[k+1] = s[k] - tmp;
//   }
   if (test == 0)
      return 0;
   else
      return -2;
}

int advance_check_Babai_heuristic_d (int cur_kappa, int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
//-----------------------------
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   ulong temp_expo;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;


   int loops = 0;

   do
   {
      test = 0;
            
      loops++;
      if (loops > 3)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < cur_kappa; j++)
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
      
      for (j = cur_kappa - 1; j > zeros; j--)
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
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
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
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                  } else
                  {
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
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
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
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
                        _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
                     } else
                     {
                        _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
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

      if (test == 1)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= cur_kappa; i++) 
	         appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = cur_kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
	   }
      else
      {
         for (i = zeros + 1; i <= cur_kappa; i++)
            appSP[kappa][i] = NAN;
      }
   } while (test == 1);

//   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
//   {
//### This is different -------
//      appSP[kappa][kappa] = _d_vec_norm(appB[kappa], n);
//-----------------------------
//   }
/*   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }
*/
   if (test == 0)
      return 0;
   else
      return -2;
}

int advance_check_Babai_heuristic(int cur_kappa, int kappa, F_mpz_mat_t B, __mpfr_struct **mu, __mpfr_struct **r, __mpfr_struct *s, 
       __mpfr_struct **appB, __mpfr_struct **appSP, 
       int a, int zeros, int kappamax, int n, mpfr_t tmp, mpfr_t rtmp, mp_prec_t prec)
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
      if (loops > 20)
         return -1;

      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      for (j = aa; j < cur_kappa; j++)
	   {	  
	      if ( mpfr_nan_p(appSP[kappa] + j) ) // if appSP[kappa] + j == NAN
	      {
            if (!(_mpfr_vec_scalar_product2(appSP[kappa] + j, appB[kappa], appB[j], n, prec) ) ){
//In this case a heuristic told us that some cancelation probably happened so we just compute the scalar product at full precision
               _F_mpz_vec_scalar_product(ztmp, B->rows[kappa], B->rows[j], n);
               F_mpz_get_mpfr(appSP[kappa] + j, ztmp);
            }
	      }
         if (j > zeros + 2)
	      {
	         mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	         mpfr_sub(rtmp, appSP[kappa] + j, tmp, GMP_RNDN);
	         for (k = zeros + 2; k < j - 1; k++)
		      {
		         mpfr_mul(tmp, mu[j] + k, r[kappa] + k, GMP_RNDN);
		         mpfr_sub(rtmp, rtmp, tmp, GMP_RNDN);
		      }
	         mpfr_mul(tmp, mu[j] + j - 1, r[kappa] + j - 1, GMP_RNDN);
	         mpfr_sub(r[kappa] + j, rtmp, tmp, GMP_RNDN);
         } else if (j == zeros+2)
	      {
	         mpfr_mul(tmp, mu[j] + zeros + 1, r[kappa] + zeros + 1, GMP_RNDN);
	         mpfr_sub(r[kappa] + j, appSP[kappa] + j, tmp, GMP_RNDN);
	      } else mpfr_set(r[kappa] + j, appSP[kappa] + j, GMP_RNDN);

	      mpfr_div(mu[kappa] + j, r[kappa] + j, r[j] + j, GMP_RNDN);
      }
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = cur_kappa - 1; j > zeros; j--)
	   {
	      /* test of the relaxed size-reduction condition */
         mpfr_abs(tmp, mu[kappa] + j, GMP_RNDN); 
	  
	      if ( mpfr_cmp_d(tmp, halfplus) > 0) 
	      {
	         test = 1; 
	      
	         /* we consider separately the cases X = +-1 */     
	         if (mpfr_cmp_d(tmp, onedothalfplus) <= 0)   
		      {
               int sgn = mpfr_sgn(mu[kappa] + j);		  
		         if (sgn >= 0)   /* in this case, X is 1 */
               {
		            for (k = zeros + 1; k < j; k++)
			         {
                     mpfr_sub(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			         }
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            mpfr_add(mu[kappa] + k, mu[kappa] + k, mu[j] + k, GMP_RNDN);
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
               }
		      } else   /* we must have |X| >= 2 */
		      {
               mpfr_round(tmp, mu[kappa] + j);
		         for (k = zeros + 1; k < j; k++)
			      {
			         mpfr_mul(rtmp, tmp, mu[j] + k, GMP_RNDN);
			         mpfr_sub(mu[kappa] + k, mu[kappa] + k, rtmp, GMP_RNDN);
			      }

	            if (mpfr_get_exp(tmp) < CPU_SIZE_1 - 2)  
		         {
                  /* X is stored in a long signed int */
                  xx = mpfr_get_si(tmp, GMP_RNDN);		      
                  if (xx > 0L)
                  { 
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                  } else
                  {
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
                  } 
               } else
		         {
                  exponent = F_mpz_set_mpfr_2exp(ztmp, tmp);
                  if (exponent <= 0){
                     F_mpz_div_2exp(ztmp, ztmp, -exponent);
                     _F_mpz_vec_submul_F_mpz(B->rows[kappa], B->rows[j], n, ztmp);
                  }
                  else{
                     _F_mpz_vec_submul_2exp_F_mpz(B->rows[kappa], B->rows[j], n, ztmp, exponent);
                  }
			      }
		      }
		   }
	   }

      if (test)   /* Anything happened? */
	   {
	      _F_mpz_vec_to_mpfr_vec(appB[kappa], B->rows[kappa], n);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
	         mpfr_set_nan(appSP[kappa] + i);//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         mpfr_set_nan(appSP[i] + kappa);//0.0/0.0;
	   }
   } while (test);




   if (mpfr_nan_p(appSP[kappa] + kappa)) 
   {
      _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
   }

/*   mpfr_set(s + zeros + 1, appSP[kappa] + kappa, GMP_RNDN);

   for (k = zeros + 1; k < kappa - 1; k++)
   {
      mpfr_mul( tmp, mu[kappa] + k, r[kappa] + k, GMP_RNDN);
      mpfr_sub( s + k + 1, s + k, tmp, GMP_RNDN);
   }
   mpfr_set(r[kappa] + kappa, s + kappa - 1, GMP_RNDN);
*/
   F_mpz_clear(ztmp);
   F_mpz_clear(X);

   return 0;
}

int check_Babai_heuristic_d_zero_vec (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
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


   int loops = 0;

   do
   {
      test = 0;
            
      loops++;
      if (loops > 20)
      {
         F_mpz_mat_print_pretty(B);
         d_mat_print(mu, expo, kappa+1, kappa+1);
         printf("kappa = %d \n", kappa);
         return -1;
      }

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
		      
		            _F_mpz_vec_sub(B->rows[kappa], B->rows[kappa], B->rows[j], n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  for (k=zeros+1; k<j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] = mu[kappa][k] + tmp;
			         }
		      
                  _F_mpz_vec_add(B->rows[kappa], B->rows[kappa], B->rows[j], n); 
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
                     _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, (ulong) xx);  
                  } else
                  {
                     _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx);  
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
                        _F_mpz_vec_submul_ui(B->rows[kappa], B->rows[j], n, xx);  
                     } else
                     {
                        _F_mpz_vec_addmul_ui(B->rows[kappa], B->rows[j], n, -xx);  
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
                        _F_mpz_vec_submul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) xx, exponent);  
                     } else
                     {
                        _F_mpz_vec_addmul_2exp_ui(B->rows[kappa], B->rows[j], n, (ulong) -xx, exponent);  
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

      if (test == 1)   /* Anything happened? */
	   {
	      expo[kappa] = _F_mpz_vec_to_d_vec_2exp(appB[kappa], B->rows[kappa], n);
         if (expo[kappa] != 0)
         {
   	      aa = zeros + 1;
	         for (i = zeros + 1; i <= kappa; i++) 
	            appSP[kappa][i] = NAN;//0.0/0.0;
	         for (i = kappa + 1; i <= kappamax; i++) 
	            appSP[i][kappa] = NAN;//0.0/0.0;
         }
         else
            test = 10;
	   }
   } while (test == 1);

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

   if (test == 0)
      return 0;
   else
      return 10;
}

// This is a mildly greedy version, tries the fast version unless that fails then 
// switches to heuristic version for only one loop and right back to fast... 

int LLL_d(F_mpz_mat_t B)
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
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int num_failed_fast = 0;
   int babai_ok = 0;
   int heuristic_fail = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   
      if (num_failed_fast < 50)
      {
         babai_ok = check_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n)); 
      }
      else{
         babai_ok = -1;
      }

      if (babai_ok == -1)
      {
         num_failed_fast++;
         heuristic_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, FLINT_MIN(kappamax + 1 + shift, n)); 
      }

      if (heuristic_fail == -1)
      {
//Got rough in there gotta switch to mpfr
         return -1;
      }

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

   return 0;
}

int LLL_d_heuristic(F_mpz_mat_t B)
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
	
//	ulong shift = getShift(B);

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
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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

   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax++; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n);//FLINT_MIN(kappamax + 1 + shift, n)); 
      if (babai_fail == -1)
         return -1;

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

   free(alpha);
   free(expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   free(s);
   free(appSPtmp);
   return 0;
}
//mpfr_init2 not mpfr_init and set_prec not set_default_prec

int LLL_mpfr2(F_mpz_mat_t B, mp_prec_t prec)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   __mpfr_struct ** mu, ** r, ** appB, ** appSP;
   __mpfr_struct * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
//	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init2(rtmp, prec);
   mpfr_init2(tmp, prec);

   mu = mpfr_mat_init2(d, d, prec);
   r = mpfr_mat_init2(d, d, prec);
   appB = mpfr_mat_init2(d, n, prec);
   appSP = mpfr_mat_init2(d, d, prec);

   s = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));
   appSPtmp = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));

   for (i = 0; i < d+1; i++){
      mpfr_init2(s + i, prec);
      mpfr_init2(appSPtmp + i, prec);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i] + j);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_to_mpfr_vec(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      _mpfr_vec_norm2(appSP[i] + i, appB[i], n, prec); 
   while ( (mpfr_sgn(appSP[i] + i) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) mpfr_set(r[i] + i, appSP[i] + i, GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int babai_fail = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      babai_fail = check_Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, n,  tmp, rtmp, prec);

      if (babai_fail == -1)
         return -1;
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1] + kappa - 1, ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s + kappa - 1) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa] + kappa - 1, r[kappa] + kappa - 1, GMP_RNDN);
         mpfr_sub(r[kappa] + kappa, s + kappa - 1, tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1] + kappa - 1, ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s + kappa - 1,tmp) <= 0) );

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

	      mpfr_set(r[kappa] + kappa, s + kappa, GMP_RNDN);
	  
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
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp + i, appSP[kappa2] + i, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp + i, appSP[i] + kappa2, GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j, GMP_RNDN);	      
	         mpfr_set(appSP[i] + kappa, appSPtmp + i - 1, GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j - 1, GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j] + i, appSP[j] + i - 1, GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa] + i, appSPtmp + i, GMP_RNDN);
	      mpfr_set(appSP[kappa] + kappa, appSPtmp + kappa2, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i] + kappa, appSPtmp + i, GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa] + kappa) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
	         mpfr_set(r[kappa] + kappa, appSP[kappa] + kappa, GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < d+1; i++){
      mpfr_clear(s + i);
      mpfr_clear(appSPtmp + i);
   }


   mpfr_mat_clear(mu, d, d);
   mpfr_mat_clear(r, d, d);
   mpfr_mat_clear(appB, d, n);
   mpfr_mat_clear(appSP, d, d);
   free(s);
   free(appSPtmp);

   return 0;
}

int LLL_mpfr(F_mpz_mat_t B)
{

   mp_prec_t prec;
//prec was 53
   prec = 53;

   int result = -1;
   int num_loops = 1;
   while ((result == -1) && (prec < MPFR_PREC_MAX)){
      result = LLL_mpfr2(B, prec);
      printf("called LLL_mpfr with prec = %ld\n", prec);
      if (result == -1){
         if (num_loops < 20)
            prec = prec + 53;
         else
            prec = prec*2;
         num_loops++;
      }
   }
   if (result >= 0)
      return result;
   else
      return -1;
}

int LLL_wrapper(F_mpz_mat_t B){

   int res = LLL_d(B);
   if (res >= 0){ //hooray worked first time
      printf("first time through, doubles are enough\n");
      return res;
   }
   else if (res == -1){ //just in case the fast/heuristic switch has any impact
      res = LLL_d_heuristic(B);
      printf("finished heuristic\n");
   }

   if (res == -1){ //Now try the mpfr version
      printf("third time through, mpfr is called\n");
      res = LLL_mpfr(B);
   }

   if (res >= 0){ //finally worked
      printf("second time through, doubles with heuristic was enough unless you saw mpfr\n");
      return res;
   }
   else //we've got big problems if this is the exit...
      return -1;
}

// This is a mildly greedy version, tries the fast version unless that fails then 
// switches to heuristic version for only one loop and right back to fast... 

int LLL_d_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
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
	
//	ulong shift = getShift(B);

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
	   expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int num_failed_fast = 0;
   int babai_ok = 0;
   int heuristic_fail = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   
      if (num_failed_fast < 500)
      {
         babai_ok = check_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); //FLINT_MIN(kappamax + 1 + shift, n)); 
      }
      else{
         babai_ok = -1;
      }

      if (babai_ok == -1)
      {
         num_failed_fast++;
         heuristic_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, kappamax, n); 
      }

      if (heuristic_fail == -1)
      {
//Got bad in there, have to switch to mpfr...
         return -1;
      }

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

//Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
//ldexp might not be the right choice as we move on... should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i],        2*expo[i] - 1);
//      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }


  
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

int LLL_d_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
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
	
//	ulong shift = getShift(B);

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
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); //FLINT_MIN(kappamax + 1 + shift, n)); 
      if (babai_fail == -1)
         return -1;

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

//newd stuff here... 
//Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
//ldexp might not be the right choice as we move on... should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i],        2*expo[i] - 1);
//      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }

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
//mpfr_init2 not mpfr_init and set_prec not set_default_prec

int LLL_mpfr2_with_removal(F_mpz_mat_t B, mp_prec_t prec, F_mpz_t gs_B)
{
   int kappa, kappa2, d, D, n, i, j, zeros, kappamax;
   __mpfr_struct ** mu, ** r, ** appB, ** appSP;
   __mpfr_struct * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;
   D = d;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
//	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init2(rtmp, prec);
   mpfr_init2(tmp, prec);

   mu = mpfr_mat_init2(d, d, prec);
   r = mpfr_mat_init2(d, d, prec);
   appB = mpfr_mat_init2(d, n, prec);
   appSP = mpfr_mat_init2(d, d, prec);

   s = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));
   appSPtmp = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));

   for (i = 0; i < d+1; i++){
      mpfr_init2(s + i, prec);
      mpfr_init2(appSPtmp + i, prec);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i] + j);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      _F_mpz_vec_to_mpfr_vec(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      _mpfr_vec_norm2(appSP[i] + i, appB[i], n, prec); 
   while ( (mpfr_sgn(appSP[i] + i) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;

   if (zeros < d - 1) mpfr_set(r[i] + i, appSP[i] + i, GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int babai_fail = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      babai_fail = check_Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, n,  tmp, rtmp, prec);

      if (babai_fail == -1)
         return -1;
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1] + kappa - 1, ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s + kappa - 1) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa] + kappa - 1, r[kappa] + kappa - 1, GMP_RNDN);
         mpfr_sub(r[kappa] + kappa, s + kappa - 1, tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1] + kappa - 1, ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s + kappa - 1,tmp) <= 0) );

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

	      mpfr_set(r[kappa] + kappa, s + kappa, GMP_RNDN);
	  
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
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp + i, appSP[kappa2] + i, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp + i, appSP[i] + kappa2, GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j, GMP_RNDN);	      
	         mpfr_set(appSP[i] + kappa, appSPtmp + i - 1, GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j - 1, GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j] + i, appSP[j] + i - 1, GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa] + i, appSPtmp + i, GMP_RNDN);
	      mpfr_set(appSP[kappa] + kappa, appSPtmp + kappa2, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i] + kappa, appSPtmp + i, GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa] + kappa) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
	         mpfr_set(r[kappa] + kappa, appSP[kappa] + kappa, GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 

//newd stuff goes here...
   int ok = 1;
   long newd = d;

   F_mpz_get_mpfr(tmp, gs_B);

   for (i = d-1; (i >= 0) && (ok > 0); i--){
//tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      mpfr_set(rtmp, r[i] + i, GMP_RNDN);
      mpfr_div_d(rtmp, rtmp, 8.0, GMP_RNDN);
//      mpfr_div_2ui(rtmp, rtmp, 1UL, GMP_RNDN);
      ok = mpfr_cmp(rtmp, tmp);
      if (ok > 0){
         newd--;
      }
   }
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < D+1; i++){
      mpfr_clear(s + i);
      mpfr_clear(appSPtmp + i);
   }


   mpfr_mat_clear(mu, D, D);
   mpfr_mat_clear(r, D, D);
   mpfr_mat_clear(appB, D, n);
   mpfr_mat_clear(appSP, D, D);
   free(s);
   free(appSPtmp);

   return newd;
}

int LLL_mpfr_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{

   mp_prec_t prec;
   prec = 53;

   int result = -1;
   int num_loops = 1;
   while ((result == -1) && (prec < MPFR_PREC_MAX)){
      printf("mpfr LLL with prec = %ld\n", prec);
      result = LLL_mpfr2_with_removal(B, prec, gs_B);
      if (result == -1){
         if (num_loops < 20)
            prec = prec + 53;
         else
            prec = prec*2;
         num_loops++;
      }
   }
   if (result >= 0)
      return result;
   else
      return -1;
}

int LLL_wrapper_with_removal(F_mpz_mat_t B, F_mpz_t gs_B){

   int res = LLL_d_with_removal(B, gs_B);
   if (res >= 0){ //hooray worked first time
      printf("first time through, doubles are enough\n");
      return res;
   }
   else if (res == -1) //just in case the fast/heuristic switch has any impact
      res = LLL_d_heuristic_with_removal(B, gs_B);

   if (res == -1){ //Now try the mpfr version
      res = LLL_mpfr_with_removal(B, gs_B);
   }

   if (res >= 0) //finally worked
      return res;
   else //we've got big problems if this is the exit...
      return -1;
}

int knapsack_LLL_wrapper_with_removal(F_mpz_mat_t B, F_mpz_t gs_B){

   int res = knapsack_LLL_d_with_removal(B, gs_B);
   if (res >= 0){ //hooray worked first time
      printf("first time through, doubles are enough\n");
      return res;
   }
   else if (res == -1) //just in case the fast/heuristic switch has any impact
      res = LLL_d_heuristic_with_removal(B, gs_B);

   if (res == -1){ //Now try the mpfr version
      printf("called mpfr!!\n");
      res = LLL_mpfr_with_removal(B, gs_B);
   }

   if (res >= 0) //finally worked
      return res;
   else //we've got big problems if this is the exit...
      return -1;
}

// This is a mildly greedy version, tries the fast version unless that fails then 
// switches to heuristic version for only one loop and right back to fast... 

int knapsack_LLL_d_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, D, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;

   int copy_kappa, copy_kappamax;
   double ** copy_mu, ** copy_r, ** copy_appB, ** copy_appSP;
   double * copy_s;
   int * copy_expo, * copy_alpha;

   int ok = 1;
   int newd = d;
   int newnewd;
   double d_rii;
   double d_gs_B;
   ulong exp;

   n = B->c;
   d = B->r;
   D = d;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
//	ulong shift = getShift(B);

   alpha = (int *) malloc(d * sizeof(int)); 
   expo = (int *) malloc(d * sizeof(int)); 

   copy_alpha = (int *) malloc(d * sizeof(int)); 
   copy_expo = (int *) malloc(d * sizeof(int)); 

   mu = d_mat_init(d, d);
   r = d_mat_init(d, d);
   appB = d_mat_init(d, n);
   appSP = d_mat_init(d, d);

   copy_mu = d_mat_init(d, d);
   copy_r = d_mat_init(d, d);
   copy_appB = d_mat_init(d, n);
   copy_appSP = d_mat_init(d, d);


   s = (double *) malloc (d * sizeof(double));
   copy_s = (double *) malloc (d * sizeof(double));
   appSPtmp = (double *) malloc (d * sizeof(double));

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         appSP[i][j] = NAN;//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int num_failed_fast = 0;
   int babai_ok = 0;
   int heuristic_fail = 0;
   long new_kappa, newvec, newvec_max;


      newvec = 0;
      newvec_max = 1;

   while (kappa < newd)
   {

      new_kappa = 0;
      if (kappa > kappamax)
      {

//In the first time we hit a new kappa we're going to size-reduce in advance...
         kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++
//         if (kappa >= 5)
         newvec++;

         if (newvec > newvec_max){
            newvec_max = newvec_max * 2;
            newvec = 0;
            new_kappa = 1;
         }
      }
      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */ 
  
      if (num_failed_fast < 150)
      {
#if PROFILE
   babai_start = get_cycle_counter();
#endif
         babai_ok = check_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
#if PROFILE
   babai_stop = get_cycle_counter();
   babai_total = babai_total + babai_stop - babai_start;
#endif
      }
      else
      {
         babai_ok = -1;
      }

      if (babai_ok == -1)
      {
#if PROFILE
   hbabai_start = get_cycle_counter();
#endif
         num_failed_fast++;
         heuristic_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, kappamax, n); 
#if PROFILE
   hbabai_stop = get_cycle_counter();
   hbabai_total = hbabai_total + hbabai_stop - hbabai_start;
#endif
      }


      if (heuristic_fail == -1)
      {
   free(alpha);
   free(expo);
   free(copy_alpha);
   free(copy_expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   d_mat_clear(copy_mu);
   d_mat_clear(copy_r);
   d_mat_clear(copy_appB);
   d_mat_clear(copy_appSP);
   free(s);
   free(copy_s);
   free(appSPtmp);
//Got bad in there, have to switch to mpfr...
         return -1;
      }
//End of the real Babai part...
      if (new_kappa == 1)
      {
//Perhaps a bit naive but just making a copy of everything but B, running ahead
//to kappa = d, without upsetting LLL... we'll see what happens.
         copy_kappa = kappa + 1;
         copy_kappamax = copy_kappa;
/*         for (i = 0; i < kappa; i++)
            for (j = 0; j < i; j++)
            {
               copy_mu[i][j] = mu[i][j];
               copy_r[i][j] = r[i][j];
               copy_appSP[i][j] = appSP[i][j];               
            }
*/
/*
         for (i = 0; i < d; i++)
            for (j = 0; j < n; j++)
            {
               copy_appB[i][j] = appB[i][j];
            }
*/

/*         for (i = 0; i < d; i++)
         {
            copy_s[i] = s[i];
//            copy_expo[i] = expo[i];
//            copy_alpha[i] = alpha[i];
         }
*/

         for (copy_kappa = newd-1; copy_kappa >  kappa; copy_kappa--)
         {
#if PROFILE
   adv_babai_start = get_cycle_counter();
#endif
            babai_ok = advance_check_Babai(kappa, copy_kappa, B, mu, r, copy_s, appB, expo, appSP, alpha[copy_kappa], zeros, copy_kappamax, n);
#if PROFILE
   adv_babai_stop = get_cycle_counter();
   adv_babai_total = adv_babai_total + adv_babai_stop - adv_babai_start;
#endif

            heuristic_fail = 0;
            if (babai_ok == -1)
            {
#if PROFILE
   hadv_babai_start = get_cycle_counter();
#endif
//               printf("heur_fail_advance ");printf("copy_kappa == %d\n", copy_kappa);
               heuristic_fail = advance_check_Babai_heuristic_d(kappa, copy_kappa, B, mu, r, copy_s, appB, expo, appSP, alpha[copy_kappa], zeros, copy_kappamax, n);
#if PROFILE
   hadv_babai_stop = get_cycle_counter();
   hadv_babai_total = hadv_babai_total + hadv_babai_stop - hadv_babai_start;
#endif

            }
         }


/*
  This isn't stable...
            if ( (heuristic_fail >= 0) && (copy_kappa == d-1) )
            {
//Use the newd stuff here...
//ldexp might not be the right choice as we move on... should make a straight d_2exp comparison
               ok = 1;
               newd = d;
               d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
               d_gs_B = ldexp( d_gs_B, exp);
               for (i = d-1; (i >= 0) && (ok > 0); i--)
               {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
                  d_rii = ldexp(copy_r[i][i], 2*copy_expo[i] - 3);
            //      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
                  if (d_rii > d_gs_B) newd--;
      else (ok = 0);
               }
               d = newd;
            }
*/
      }
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

         if (kappa >= newd - 1){
            printf("kappa2 = %ld, newd = %ld \n", kappa, newd);
            ok = 1;
            newnewd = newd;
            d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
            d_gs_B = ldexp( d_gs_B, exp);
            for (i = newd-1; (i >= 0) && (ok > 0); i--)
            {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
               d_rii = ldexp(r[i][i],        2*expo[i] - 1);
//      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
               if (d_rii > d_gs_B) newnewd--;
               else (ok = 0);
            }
            if (newnewd < newd)
               printf(" ooh something could be done %ld\n", newnewd);

         }

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

//Use the newd stuff here...
//ldexp might not be the right choice as we move on... should make a straight d_2exp comparison
   ok = 1;
   newnewd = newd;
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = newd-1; (i >= 0) && (ok > 0); i--)
   {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i],        2*expo[i] - 1);
//      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newnewd--;
      else (ok = 0);
   }
  
   free(alpha);
   free(expo);
   free(copy_alpha);
   free(copy_expo);
   d_mat_clear(mu);
   d_mat_clear(r);
   d_mat_clear(appB);
   d_mat_clear(appSP);
   d_mat_clear(copy_mu);
   d_mat_clear(copy_r);
   d_mat_clear(copy_appB);
   d_mat_clear(copy_appSP);
   free(s);
   free(copy_s);
   free(appSPtmp);

   return newnewd;
}

int knapsack_LLL_d_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
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
	
//	ulong shift = getShift(B);

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
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
      if (babai_fail == -1)
         return -1;

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

//newd stuff here... 
//Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
//ldexp might not be the right choice as we move on... should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i],        2*expo[i] - 1);
//      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }

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
//mpfr_init2 not mpfr_init and set_prec not set_default_prec

int knapsack_LLL_mpfr2_with_removal(F_mpz_mat_t B, mp_prec_t prec, F_mpz_t gs_B)
{
   int kappa, kappa2, d, D, n, i, j, zeros, kappamax;
   __mpfr_struct ** mu, ** r, ** appB, ** appSP;
   __mpfr_struct * s, * mutmp, * appBtmp, * appSPtmp;
   mpfr_t tmp, rtmp;
   F_mpz_t ztmp;
   int * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;
   D = d;

   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;
	
//	ulong shift = getShift(B);

   alpha = (int *) malloc((d + 1) * sizeof(int)); 

   F_mpz_init(ztmp);
   mpfr_init2(rtmp, prec);
   mpfr_init2(tmp, prec);

   mu = mpfr_mat_init2(d, d, prec);
   r = mpfr_mat_init2(d, d, prec);
   appB = mpfr_mat_init2(d, n, prec);
   appSP = mpfr_mat_init2(d, d, prec);

   s = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));
   appSPtmp = (__mpfr_struct *) malloc ((d + 1) * sizeof(__mpfr_struct));

   for (i = 0; i < d+1; i++){
      mpfr_init2(s + i, prec);
      mpfr_init2(appSPtmp + i, prec);
   }

   for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
         mpfr_set_nan(appSP[i] + j);//0.0/0.0;
  
   /* ************************** */
   /* Step1: Initialization Step */
   /* ************************** */     
    
   for (i = 0; i < d; i++)
	   _F_mpz_vec_to_mpfr_vec(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      _mpfr_vec_norm2(appSP[i] + i, appB[i], n, prec); 
   while ( (mpfr_sgn(appSP[i] + i) == 0) && (++i < d));

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) mpfr_set(r[i] + i, appSP[i] + i, GMP_RNDN);

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;

   int babai_fail = 0;
    
   while (kappa < d)
   {      


      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      babai_fail = check_Babai_heuristic(kappa, B, mu, r, s, appB, appSP, alpha[kappa], zeros, 
			                        kappamax, n,  tmp, rtmp, prec);

      if (babai_fail == -1)
         return -1;
      
      /* ************************************ */
      /* Step4: Success of Lovasz's condition */
      /* ************************************ */  
      /* ctt * r.coeff[kappa-1][kappa-1] <= s[kappa-2] ?? */

      mpfr_mul_d( tmp, r[kappa - 1] + kappa - 1, ctt, GMP_RNDN);
      if ( mpfr_cmp(tmp, s + kappa - 1) <= 0) 
	   {
	      alpha[kappa] = kappa;
         mpfr_mul(tmp, mu[kappa] + kappa - 1, r[kappa] + kappa - 1, GMP_RNDN);
         mpfr_sub(r[kappa] + kappa, s + kappa - 1, tmp, GMP_RNDN);
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
               mpfr_mul_d(tmp, r[kappa-1] + kappa - 1, ctt, GMP_RNDN);
	         }
         } while ( (kappa >= zeros + 2) && (mpfr_cmp(s + kappa - 1,tmp) <= 0) );

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

	      mpfr_set(r[kappa] + kappa, s + kappa, GMP_RNDN);
	  
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
	  
	      for (i = 0; i <= kappa2; i++) mpfr_set(appSPtmp + i, appSP[kappa2] + i, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSPtmp + i, appSP[i] + kappa2, GMP_RNDN);
	  
	      for (i = kappa2; i > kappa; i--)
	      {
	         for (j = 0; j < kappa; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j, GMP_RNDN);	      
	         mpfr_set(appSP[i] + kappa, appSPtmp + i - 1, GMP_RNDN);
	      
	         for (j = kappa + 1; j <= i; j++) mpfr_set(appSP[i] + j, appSP[i-1] + j - 1, GMP_RNDN);

	         for (j = kappa2 + 1; j <= kappamax; j++) mpfr_set(appSP[j] + i, appSP[j] + i - 1, GMP_RNDN);     
	      }
	  
	      for (i = 0; i < kappa; i++) mpfr_set(appSP[kappa] + i, appSPtmp + i, GMP_RNDN);
	      mpfr_set(appSP[kappa] + kappa, appSPtmp + kappa2, GMP_RNDN);

	      for (i = kappa2 + 1; i <= kappamax; i++) mpfr_set(appSP[i] + kappa, appSPtmp + i, GMP_RNDN);
	  
	      if ( mpfr_sgn(r[kappa] + kappa) <= 0.0)
	      {
	         zeros++;
	         kappa++;
	         _mpfr_vec_norm2(appSP[kappa] + kappa, appB[kappa], n, prec);
	         mpfr_set(r[kappa] + kappa, appSP[kappa] + kappa, GMP_RNDN);
	      }
	  
	      kappa++;
	   }
   } 

//newd stuff goes here...
   int ok = 1;
   long newd = d;

   F_mpz_get_mpfr(tmp, gs_B);

   for (i = d-1; (i >= 0) && (ok > 0); i--){
//tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      mpfr_set(rtmp, r[i] + i, GMP_RNDN);
      mpfr_div_d(rtmp, rtmp, 8.0, GMP_RNDN);
//      mpfr_div_2ui(rtmp, rtmp, 1UL, GMP_RNDN);
      ok = mpfr_cmp(rtmp, tmp);
      if (ok > 0){
         newd--;
      }
   }
  
   free(alpha);

   F_mpz_clear(ztmp);
   mpfr_clear(rtmp);
   mpfr_clear(tmp);

   for (i = 0; i < D+1; i++){
      mpfr_clear(s + i);
      mpfr_clear(appSPtmp + i);
   }


   mpfr_mat_clear(mu, D, D);
   mpfr_mat_clear(r, D, D);
   mpfr_mat_clear(appB, D, n);
   mpfr_mat_clear(appSP, D, D);
   free(s);
   free(appSPtmp);

   return newd;
}

int knapsack_LLL_mpfr_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
{

   mp_prec_t prec;
   prec = 53;

   int result = -1;
   int num_loops = 1;
   while ((result == -1) && (prec < MPFR_PREC_MAX)){
      result = knapsack_LLL_mpfr2_with_removal(B, prec, gs_B);
      if (result == -1){
         if (num_loops < 20)
            prec = prec + 53;
         else
            prec = prec*2;
         num_loops++;
      }
   }
   if (result >= 0)
      return result;
   else
      return -1;
}

int knapsack_LLL_with_removal(F_mpz_mat_t B, F_mpz_t gs_B){

   int res = knapsack_LLL_d_with_removal(B, gs_B);
   if (res >= 0) //hooray worked first time
      return res;
   else if (res == -1) //just in case the fast/heuristic switch has any impact
      res = knapsack_LLL_d_heuristic_with_removal(B, gs_B);

   if (res == -1) //Now try the mpfr version
      res = knapsack_LLL_mpfr_with_removal(B, gs_B);

   if (res >= 0) //finally worked
      return res;
   else //we've got big problems if this is the exit...
      return -1;
}


int LLL_d_zero_vec_heuristic_with_removal(F_mpz_mat_t B, F_mpz_t gs_B)
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
	
//	ulong shift = getShift(B);

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
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
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
   kappamax = kappa;

   printf("zeros = %d \n", zeros);
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;
    
   while (kappa < d)
   {      
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      /* ********************************** */
      /* Step3: Call to the Babai algorithm */
      /* ********************************** */   

      int babai_fail = check_Babai_heuristic_d_zero_vec(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
      if (babai_fail == -1)
         return -1;
      if (babai_fail == 10)
      {
//This means that the kappa^th vector is a 0-vector... 
         printf(" found a 0 vector!\n");

      }

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

//newd stuff here... 
//Use the newd stuff here...
   int ok = 1;
   int newd = d;
   double d_rii;
   double d_gs_B;
   ulong exp;
//ldexp might not be the right choice as we move on... should make a straight d_2exp comparison
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i],        2*expo[i] - 1);
//      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }

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

int LLL_wrapper_zero_vec_with_removal(F_mpz_mat_t B, F_mpz_t gs_B){

   int res;

   res = LLL_d_zero_vec_heuristic_with_removal(B, gs_B);

   if (res == -1)
   { //Now try the mpfr version
      printf("trying mpfr version... bad\n");
      abort();
      res = LLL_mpfr_with_removal(B, gs_B);
   }

   if (res >= 0) //finally worked
      return res;
   else //we've got big problems if this is the exit...
      return -1;
}

/**

   Stripped down LLL_d to just calculate G-S lengths.

**/

//### This is different -------
void gs_Babai(int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n)
//-----------------------------
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
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

      for (j = kappa - 1; j > zeros; j--)
	   {
	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      tmp = ldexp(tmp, expo[kappa] - expo[j]);

      }
      
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

ulong F_mpz_mat_gs_d( F_mpz_mat_t B, F_mpz_t gs_B)
{
   int kappa, kappa2, d, n, i, j, zeros, kappamax;
   double ** mu, ** r, ** appB, ** appSP;
   double * s, * mutmp, * appBtmp, * appSPtmp;
   double tmp = 0.0;
   int * expo, * alpha;
   mp_limb_t * Btmp;
   
   n = B->c;
   d = B->r;
	
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
      expo[i] = _F_mpz_vec_to_d_vec_2exp(appB[i], B->rows[i], n);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  

   i = 0; 
  
   do
      appSP[i][i] = _d_vec_norm(appB[i], n); 
   while ((appSP[i][i] <= 0.0) && (++i < d)); // Fixme : should this be EPS not 0.0

   zeros = i - 1; /* all vectors B[i] with i <= zeros are zero vectors */
   kappa = i + 1;
   kappamax = kappa;
  
   if (zeros < d - 1) r[i][i] = appSP[i][i];

   for (i = zeros + 1; i < d; i++)
      alpha[i] = 0;   

   while (kappa < d)
   {
      if (kappa > kappamax) kappamax = kappa; // Fixme : should this be kappamax = kappa instead of kappamax++

      gs_Babai(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
			                        kappamax, n); 
	   alpha[kappa] = kappa;
	   tmp = mu[kappa][kappa-1] * r[kappa][kappa-1];
	   r[kappa][kappa] = s[kappa-1] - tmp;
	   kappa++;
   }

   ulong ok, newd, exp;
   double d_gs_B, d_rii; 
   ok = 1;
   newd = d;
   d_gs_B = F_mpz_get_d_2exp(&exp, gs_B);
   d_gs_B = ldexp( d_gs_B, exp);
   for (i = d-1; (i >= 0) && (ok > 0); i--)
   {
//d_rii is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      d_rii = ldexp(r[i][i], expo[i]);
      printf("%5f r[%d] and gs_B = %5f\n", d_rii, i, d_gs_B);
      if (d_rii > d_gs_B) newd--;
      else (ok = 0);
   }


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

int U_LLL_with_removal(F_mpz_mat_t FM, long new_size, F_mpz_t gs_B){

   long r, c, bits, i, j;
   int full_prec = 1;
   int done = 0;
   clock_t lll_start, lll_stop, lll_total, sum_start, sum_stop;
   int is_U_I;

   lll_total = 0;
   sum_start = clock();

   r = FM->r;
   c = FM->c;
   bits = FLINT_ABS(F_mpz_mat_max_bits(FM));

   F_mpz_mat_t U;
   //F_mpz_mat_init(U, r, r);

   F_mpz_mat_t I;
   F_mpz_mat_init_identity(I, r);

   F_mpz_mat_t full_U;
   F_mpz_mat_init_identity(full_U, r);

   F_mpz_mat_t big_FM;
   F_mpz_mat_init(big_FM, r, c + r);

   F_mpz_mat_t full_data;
   F_mpz_mat_init(full_data, r, c);

   F_mpz_mat_t trunc_data;
   F_mpz_mat_init(trunc_data, r, c);

   long mbits;

   int k = 1;

   int newd;
   long prev_mbits = bits;

   if (bits > new_size){
      full_prec = 0;
//do some truncating
      for ( i = 0; i < r; i++)
         for ( j = 0; j < c; j++)
            F_mpz_set(full_data->rows[i]+j, FM->rows[i]+j);

      mbits = FLINT_ABS(F_mpz_mat_max_bits(full_data));

      if ((mbits - new_size) > 0){
         F_mpz_mat_resize(trunc_data, full_data->r, full_data->c);
		 F_mpz_mat_div_2exp(trunc_data, full_data, (ulong) (mbits - new_size));
//Make this iterate over i and j, make a LARGE lattice which has identity in one corner and FM in the other
         for ( i = 0; i < r; i++){
            for (j = 0; j < i; j++)
               F_mpz_set_ui(big_FM->rows[i]+j, 0L);
            F_mpz_set_ui(big_FM->rows[i]+i, 1L);
            for (j = i+1; j < r; j++)
               F_mpz_set_ui(big_FM->rows[i]+j, 0L);
            for (j = r; j < r+c; j++)
               F_mpz_set(big_FM->rows[i]+j, trunc_data->rows[i] + j-r);
         }
      }
      else{
         printf("something odd here\n");
         full_prec = 1;
      }
   }

   while( done == 0){
      k++;
      if (full_prec == 0){
         lll_start = clock();
         knapsack_LLL_wrapper_with_removal(big_FM, gs_B);
         lll_stop = clock();
         printf("was big_FM\n");
      }
      else{
         lll_start = clock();
         newd = knapsack_LLL_wrapper_with_removal(FM, gs_B);
         lll_stop = clock();
         printf("was FM\n");
      }

      lll_total = lll_total + lll_stop - lll_start;

      if (full_prec == 1)
         done = 1;
      else {
//add more bits

         F_mpz_mat_window_init(U, big_FM, 0, 0, big_FM->r, r);

         printf("U bits == %ld\n", FLINT_ABS(F_mpz_mat_max_bits(U)));

         is_U_I = F_mpz_mat_equal(U, I);

//do some truncating
         F_mpz_mat_mul_classical(full_data, U, full_data);

         mbits = FLINT_ABS(F_mpz_mat_max_bits(full_data));
//make this condition better?
         if ( ( (mbits - new_size) > 0) &&  ( mbits <= prev_mbits - (long)(new_size/4) ) && (is_U_I == 0)){
            F_mpz_mat_div_2exp(trunc_data, full_data, (ulong) (mbits - new_size));
          //  F_mpz_mat_mul_2exp(trunc_data, trunc_data, (ulong) (new_size*2));
         }
         else{
     //       printf(" the mbits business == %d\n", ( mbits <= prev_mbits - (long)(new_size/4)));
            full_prec = 1;
         }

         prev_mbits = mbits;

         if (full_prec == 1){
//can switch to FM, no need for a new identity
            for ( i = 0; i < r; i++){
               for (j = 0; j < c; j++)
                  F_mpz_set(FM->rows[i]+j, full_data->rows[i] + j);
            }
         }
         else{
//keep with the big_FM concept
            for ( i = 0; i < r; i++){
               for (j = 0; j < i; j++)
                  F_mpz_set_ui(big_FM->rows[i]+j, 0L);
               F_mpz_set_ui(big_FM->rows[i]+i, 1L);
               for (j = i+1; j < r; j++)
                  F_mpz_set_ui(big_FM->rows[i]+j, 0L);
               for (j = r; j < r+c; j++)
                  F_mpz_set(big_FM->rows[i]+j, trunc_data->rows[i] + j-r);
            }
         }
         F_mpz_mat_window_clear(U);
      }

   }

   sum_stop = clock();

#if PROFILE
   printf(" spent a total of %3f seconds on regular Babai\n", (double) babai_total / 2.4E9);
   printf(" of which %3f cycles spent doing ldexps\n", (double) ldexp_total / 2.4E9);
   printf(" of which %3f cycles spent updating full precision B\n", (double) update_total / 2.4E9);
   printf(" of which %3f cycles spent converting full precision B\n", (double) convert_total / 2.4E9);
   printf(" of which %3f cycles spent computing inner products\n", (double) inner_total /2.4E9);
   printf(" spent a total of %3f seconds on advanced Babai\n", (double) adv_babai_total /2.4E9);

   printf(" spent a total of %3f seconds on h regular Babai\n", (double) hbabai_total / 2.4E9);
   printf(" of which %3f cycles spent doing h ldexps\n", (double) hldexp_total / 2.4E9);
   printf(" of which %3f cycles spent updating h full precision B\n", (double) hupdate_total / 2.4E9);
   printf(" of which %3f cycles spent converting h full precision B\n", (double) hconvert_total / 2.4E9);
   printf(" of which %3f cycles spent computing h inner products\n", (double) hinner_total /2.4E9);
   printf(" spent a total of %3f seconds on h advanced Babai\n", (double) hadv_babai_total /2.4E9);
#endif

   printf(" spent a total of %f seconds on inner LLL\n", (double) lll_total / (double)CLOCKS_PER_SEC);

   printf(" spent a total of %f seconds in ULLL\n", (double) (sum_stop - sum_start) / (double)CLOCKS_PER_SEC);

   F_mpz_mat_clear(full_data);
   F_mpz_mat_clear(trunc_data);
   F_mpz_mat_clear(big_FM);
   F_mpz_mat_clear(I);
   F_mpz_mat_clear(full_U);

   return newd;
}

