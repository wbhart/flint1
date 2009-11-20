#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "flint.h"
#include "F_mpz_poly.h"
#include "F_mpz_LLL_fast_d.h"
#include "d_mat.h"

ulong getShift(F_mpz_mat_t B)
{
   ulong n = B->c;
   ulong shift = 0;
   for (ulong i = 0; i < B->r; i++)
   {
      ulong j;
      for (j = n - 1; j >= i + shift + 1 && F_mpz_size(B->rows[i] + j) == 0L; j--);  
      
      if (shift < j - i) shift = j - i;
      
   }

   return shift;
}

long F_mpz_mat_set_line_d_2exp(double * appv, const F_mpz_mat_t mat, const ulong r, const int n, int * cexpo)
{
   long * exp, i, maxexp = 0L;
   exp = (long *) malloc(n * sizeof(long)); 
  
   for (i = 0; i < n; i++)
   {
      appv[i] = F_mpz_get_d_2exp(&exp[i], mat->rows[r] + i);
      //exp[i] = exp[i] + cexpo[i];
      if (exp[i] > maxexp) maxexp = exp[i];
   }

   for (i = 0; i < n; i++) appv[i] = ldexp(appv[i], exp[i] - maxexp);

   free(exp);
   return maxexp;
}


void Babai_2exp (int kappa, F_mpz_mat_t B, double **mu, double **r, double *s, 
       double **appB, int *expo, double **appSP, 
       int a, int zeros, int kappamax, int n, int *cexpo)
{
   int i, j, k, test, aa, exponent;
   signed long xx;
   double tmp, rtmp;
   
   aa = (a > zeros) ? a : zeros + 1;
  
   ctt = DELTA;
   halfplus = ETA;
   onedothalfplus = 1.0+halfplus;

   int loops = 0;
   int count_trig = 0;
   double try_temp[kappa];

   printf("kappa = %d\n", kappa);
   
   do
   {
      test = 0;
      
      loops++;
      
      /* ************************************** */
      /* Step2: compute the GSO for stage kappa */
      /* ************************************** */
      
      printf(" aa = %d\n", aa);
      for (j = aa; j < kappa; j++)
	   {	  
	      printf("appSP[kappa][j] = %le\n", appSP[kappa][j]);
         if (appSP[kappa][j] != appSP[kappa][j]) // if appSP[kappa][j] == NAN
	      {
	         appSP[kappa][j] = d_vec_scalar_product(appB[kappa], appB[j], n);
            printf("worked %le\n", appSP[kappa][j]);
            
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

	      printf("r[kappa][j] = %le, r[j][j] = %le\n", r[kappa][j], r[j][j]);
         mu[kappa][j] = r[kappa][j] / r[j][j];
      }

      d_mat_print(mu, expo, kappa, kappa);
      d_mat_print(r, expo, kappa, kappa);
      
      /* **************************** */
      /* Step3--5: compute the X_j's  */
      /* **************************** */
      
      for (j = kappa - 1; j > zeros; j--)
	   {

	      /* test of the relaxed size-reduction condition */
	      tmp = fabs(mu[kappa][j]);
	      printf("******************************** %lu, %lu\n", expo[kappa], expo[j]);
         tmp = ldexp(tmp, expo[kappa] - expo[j]);

         if (count_trig > 0)
            printf("j=%d, fabs(mu)=%f150,exp=%d\n", j, fabs(mu[kappa][j]),expo[kappa]-expo[j]);

         try_temp[j] = tmp;
	  
	      printf("tmp = %lf, halfplus = %lf\n", tmp, halfplus);
         if (tmp > halfplus) 
	      {
	         printf("j = %d\n", j);
            test = 1; 
	         exponent = expo[j] - expo[kappa];
	      
	         /* we consider separately the cases X = +-1 */     
	         if (tmp <= onedothalfplus)   
		      {		  
		         if (mu[kappa][j] >= 0)   /* in this case, X is 1 */
               {
		            printf("X is 1 ...................................\n");
                  for (k = zeros + 1; k < j; k++)
			         {
			            tmp = ldexp (mu[j][k], exponent);
			            mu[kappa][k] =  mu[kappa][k] - tmp; 
			         }
		      
		            F_mpz_mat_row_sub(B, kappa, B, kappa, B, j, 0, n);
		  
		         } else          /* otherwise X is -1 */ 
               {
                  printf("X is minus 1 ...................................\n");
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

		            printf("|X| is < MAX_LONG ...................................\n");
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
		            printf("mu[kappa][j] = %le\n", mu[kappa][j]);
                  tmp = frexp(mu[kappa][j], &exponent);
                  printf("exponent = %d, tmp = %le\n", exponent, tmp);
		            tmp = tmp * MAX_LONG;
		            xx = (long) tmp;
		            printf("xx = %ld, expo[kappa] = %ld, expo[j] = %ld\n", xx, expo[kappa], expo[j]);
                  exponent += (expo[kappa] - expo[j] - CPU_SIZE_1);

		            /* This case is extremely rare: never happened for me */
		            if (exponent <= 0) 
			         {
                     printf("X is in rare case ...................................\n");
                  
			            xx <<= -exponent;
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

                     printf("bigger than 2 notrare case xx=%ld, exp = %ld, %ld\n", xx, exponent, CPU_SIZE_1);

			            if (xx > 0)
                     {
                        F_mpz_mat_row_submul_2exp_ui(B, kappa, B, j, 0, n, (ulong) xx, exponent);  
                     } else
                     {
                        F_mpz_mat_row_addmul_2exp_ui(B, kappa, B, j, 0, n, (ulong) -xx, exponent);  
                     }
			            for (k = zeros + 1 ; k < j; k++)
			            {
			               rtmp = ((double) xx) * mu[j][k];
			               rtmp = ldexp(rtmp, exponent + expo[j] - expo[kappa]);
			               mu[kappa][k] = mu[kappa][k] - rtmp;
					      }
				      }	    
			      }
		      }
		   }
	   }

      if( loops > 10000){
         printf("infinite loop at kappa = %d\n", kappa);
         F_mpz_mat_print_pretty(B);
         d_mat_print(appB, expo, B->r, B->c);
         d_mat_print(appSP, expo, B->r, B->r);
         d_mat_print(mu, expo, B->r, B->r);
         d_mat_print(r, expo, B->r, B->r);
         count_trig++;
         for (j = zeros; j < kappa; j++)
            printf("tmpj=%f,",try_temp[j]);
         printf("\n");

         if (count_trig == 5)
            abort();
      }

      if (test)   /* Anything happened? */
	   {
	      expo[kappa] = F_mpz_mat_set_line_d_2exp(appB[kappa], B, kappa, n, cexpo);
	      aa = zeros + 1;
	      for (i = zeros + 1; i <= kappa; i++) 
            appSP[kappa][i] = NAN;//0.0/0.0;
	      for (i = kappa + 1; i <= kappamax; i++) 
	         appSP[i][kappa] = NAN;//0.0/0.0;
      }
   } while (test);

   if (appSP[kappa][kappa] != appSP[kappa][kappa]) 
   {
      appSP[kappa][kappa] = d_vec_norm(appB[kappa], n);
   }
   s[zeros + 1] = appSP[kappa][kappa];
  
   for (k = zeros + 1; k < kappa - 1; k++)
   {
      tmp = mu[kappa][k] * r[kappa][k];
      s[k+1] = s[k] - tmp;
   }
}


int LLL_2exp_with_removal(F_mpz_mat_t B, int *cexpo, F_mpz_t gs_B)
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
      expo[i] = F_mpz_mat_set_line_d_2exp(appB[i], B, i, n, cexpo);  
  
   /* ********************************* */
   /* Step2: Initializing the main loop */
   /* ********************************* */   
  
   kappamax = 0;
   i = 0; 
  
   do
      appSP[i][i] = d_vec_norm(appB[i], n); 
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
      
      Babai_2exp(kappa, B, mu, r, s, appB, expo, appSP, alpha[kappa], zeros, 
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
	         appSP[kappa][kappa] = d_vec_norm(appB[kappa], n);
	         r[kappa][kappa] = appSP[kappa][kappa];
	      }
	  
	      kappa++;
	   }
   }
   F_mpz_t tmp_gs;
   F_mpz_init(tmp_gs);
 
   int ok = 1;
   int newd = d;
   double rtmp;
   for (i = d-1; (i >= 0) && (ok > 0); i--){
//tmp_gs is the G-S length of ith vector divided by 2 (we shouldn't make a mistake and remove something valuable)
      rtmp = sqrt(r[i][i]);
      if (rtmp < 0.00001)
         rtmp = 0.00001;
      F_mpz_set_d_2exp(tmp_gs, rtmp, expo[i] - 1);
      printf("sqrt(rii)=%f,expo=%d,", rtmp, expo[i]-1);
      ok = F_mpz_cmpabs(tmp_gs, gs_B);
      printf("i=%d,tmp_gs=", i); F_mpz_print(tmp_gs); printf(",ok=%d:::",ok);
      if (ok > 0) newd--;
   }
   printf("\n");
   F_mpz_clear(tmp_gs);
   
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

int main(int argc, char * argv[]){

   F_mpz_mat_t M;
   F_mpz_mat_init(M, 0, 0);

   F_mpz_mat_fread(M, stdin);

//   printf("trying to run LLL...\n");

   int cexpo[135];

   for( int i = 0; i < 135; i++){
      cexpo[i]=0;
   }

   cexpo[128] = 0;
   cexpo[129] = 0;
   cexpo[130] = 0;
   cexpo[131] = 0;

   F_mpz_t B;
   F_mpz_init(B);

   F_mpz_set_ui(B, 129);

//   F_mpz_mat_print_pretty(M);

//   F_mpz_mat_resize(M, 3, M->c);

   F_mpz_mat_print_pretty(M);

   LLL_2exp_with_removal(M, cexpo, B);

//   LLL(M);

   F_mpz_mat_print_pretty(M);

   F_mpz_clear(B);

   F_mpz_mat_clear(M);

   return 0;
}
