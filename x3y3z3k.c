/****************************************************************************

   x3y3z3k.c: program to find solutions to x^3 +y^3 +z^3 = k 
              based on the algorithm of Elkies and the Magma 
              implementation of Samir Siksek

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <sys/types.h>

#include "vector.h"
#include "matrix.h"
#include "vec3d.h"
#include "mat3d.h"
#include "memory-manager.h"
#include "qd/c_dd.h"

mpz_t tmp, tmp2;

FILE * outfile;
char * done;

#define DEBUG 0
#define DEBUG2 1

#define COORD(A,x) ((A)[x-1])
#define SQ(A,x) (((A)[x-1])*((A)[x-1]))
#define ROW(R, x) ((R)[x-1])
#define COL(R, x) ((R)[x-1])
#define MAT_R(R, x, y) ((R)[x-1][y-1])
#define MAT_C(R, x, y) ((R)[y-1][x-1])

void elkies(double a, long N, double eps)
{
   static long told;
   long t1, t2, t3, k, l1, l2, l3;
   double b = -powl(1.0+a*a*a, 1.0/3.0);
   double alpha = -(5.0*pow(a,6.0)+6.0*a*a*a*b*b*b+powl(b,6.0))/(6.0*b*b*b*pow(a,4.0)+6.0*a*pow(b,6.0));
   double c = -a*a/(b*b);
   double d = 1.0/(-6.0*pow(b,5.0))*pow(a,4.0) - 5.0/(6.0*b*b)*a;
   double Ne = (double)N*eps;
   double Ne2 = Ne*eps;
   d_vec3d v;
   d_mat3dc_t L, M, R;
   z_mat3dc_t Z;
   d_mat3dc_stack_init(&M);
   d_mat3dc_stack_init(&L);
   d_mat3dc_stack_init(&R);
   z_mat3dc_stack_init(&Z);
   d_vec3d_stack_init(&v);
   MAT_C(M, 1, 1) = (double)N*a;
   MAT_C(M, 2, 1) = (double)N*b;
   MAT_C(M, 3, 1) = (double)N;
   MAT_C(M, 1, 2) = Ne;
   MAT_C(M, 2, 2) = Ne*c;
   MAT_C(M, 3, 2) = 0.0;
   MAT_C(M, 1, 3) = alpha*Ne2;
   MAT_C(M, 2, 3) = d*Ne2;
   MAT_C(M, 3, 3) = 0.0;
#if DEBUG2
   printf("M = "); d_mat3dc_printf(M); printf("\n");
#endif
   d_mat3dc_invert(L, M);
#if DEBUG2
   printf("M^(-1) = "); d_mat3dc_printf(L);printf("\n");
#endif
   if (d_mat3dc_LLL(Z, R, L, 0.999999))
   {
#if DEBUG2
   printf("Z = "); z_mat3dc_printf(Z);printf("\n\n");
#endif
#if DEBUG2
   printf("R = "); d_mat3dc_printf(R);printf("\n\n");
#endif
   l1 = (long)ceill(3.0/d_vec3d_norm(COL(R, 1)));
   l2 = (long)ceill(3.0/d_vec3d_norm(COL(R, 2)));
   l3 = (long)ceill(3.0/d_vec3d_norm(COL(R, 3)));
   if ((l1 > 100) || (l1 < 0)) l1 = 100;
   if ((l2 > 100)  || (l2 < 0)) l2 = 100;
   if ((l3 > 100)  || (l3 < 0)) l3 = 100;
   for (long i1 = -l1; i1 <= l1; i1++)
   {
   for (long i2 = -l2; i2 <= l2; i2++)
   {
   for (long i3 = -l3; i3 <= l3; i3++)
   {
      if (!i1) i1++;
      if (!i2) i2++;
      if (!i3) i3++;
      d_vec3d_mul_scalar(v, COL(R,1), (double) i1);
      d_vec3d_add_scalar_mul(v, v, COL(R,2), (double) i2);
      d_vec3d_add_scalar_mul(v, v, COL(R,3), (double) i3);
      if (d_vec3d_norm(v) < 3.0) 
      {
      t1 = (long) round(COORD(v, 1)*MAT_C(M, 1, 1) + COORD(v, 2)*MAT_C(M, 1, 2) + COORD(v,3)*MAT_C(M, 1, 3));
      if ((t1 != told) && (t1 != - told) && ((t1 > 100) || (t1 < -100)))
      {
         t2 = (long) round(COORD(v, 1)*MAT_C(M, 2, 1) + COORD(v, 2)*MAT_C(M, 2, 2) + COORD(v, 3)*MAT_C(M, 3, 3));
         t3 = (long) round(COORD(v, 1)*MAT_C(M, 3, 1) + COORD(v, 2)*MAT_C(M, 3, 2) + COORD(v, 3)*MAT_C(M, 3, 3));
         mpz_set_si(tmp, t1);
         mpz_pow_ui(tmp, tmp, 3);
         mpz_set_si(tmp2, t2);
         mpz_pow_ui(tmp2, tmp2, 3);
         mpz_add(tmp, tmp, tmp2);
         mpz_set_si(tmp2, t3);
         mpz_pow_ui(tmp2, tmp2, 3);
         mpz_add(tmp, tmp, tmp2);
         if (mpz_size(tmp) == 1) k = mpz_get_ui(tmp);
         else k = 0;
         if ((k >= 2) && (k < 10000) && ((fabs(t1) > 100000) || (fabs(t2) > 100000) || (fabs(t3) > 100000))) 
         {
            if (!done[k])
            {
               fprintf(outfile,"\n%ld^3 + %ld^3 + %ld^3  = %ld\n\n", t1, t2, t3, k);
               fflush(outfile);
            }
            done[k] = 1;
            told = t1;
         }
      }
      }
   }
   }
   }
   } 
   d_vec3d_stack_clear();
   z_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
   d_mat3dc_stack_clear();
}

int main(void)
{
   long N;
   double eps, a;
   char buf[1000];
   char* ptr;
   
   mpz_init(tmp);
   mpz_init(tmp2);
   
   done = (char *) flint_stack_alloc_bytes(10000); 
   for (unsigned long i = 0; i < 10000; i++) done[i] = 1;
   FILE * unknown = fopen("unknown", "r");
   while (fgets(buf, 1000, unknown))
   {
      ptr = strtok(buf, ", \n\r");
      done[atoi(ptr)] = 0;
      while (ptr = strtok(NULL, ", \n\r"))
      {
         done[atoi(ptr)] = 0;
      }
   }
   
   outfile = fopen("test-comp.out", "w");
   N = 100000000;
   eps = ((double)1.0)/((double)N);
   for (long k = -20000000; k > -25000005; k--) //k > -N+20000000; k--) 
   {
      if ((k&0x7ffff) == 0) 
      {
         fprintf(outfile, "Passed k = %ld\n", k);
         fflush(outfile);
      }
      a = ((double)k)/((double) N);
      elkies(a, N, eps);
   }
   /*for (long k = 20000000; k < N-20000000; k++) 
   {
      if ((k&0x7ffff) == 0) 
      {
         fprintf(outfile, "Passed k = %ld\n", k);
         fflush(outfile);
      }
      a = k*eps;
      elkies(a, N, eps);
   }*/ 
   
   fclose(unknown);
   flint_stack_release();
   fclose(outfile);
   mpz_clear(tmp);
   mpz_clear(tmp2);
}
