/****************************************************************************

   x3y3z3k.c: program to find solutions to x^3 +y^3 +z^3 = k 
              based on the algorithm of Elkies and the Magma 
              implementation of Samir Siksek

   Copyright (C) 2007, William Hart

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>

#include "vector.h"
#include "matrix.h"
#include "vec3d.h"
#include "mat3d.h"

mpz_t tmp, tmp2;

FILE * outfile;

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
   double b = -pow(1.0+a*a*a, 1.0/3.0);
   double alpha = -(5*pow(a,6)+6*a*a*a*b*b*b+pow(b,6))/(6*b*b*b*pow(a,4)+6*a*pow(b,6));
   double c = -a*a/(b*b);
   double d = 1/(-6*pow(b,5))*pow(a,4) - 5/(6*b*b)*a;
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
   MAT_C(M, 1, 1) = N*a;
   MAT_C(M, 2, 1) = N*b;
   MAT_C(M, 3, 1) = (double)N;
   MAT_C(M, 1, 2) = Ne;
   MAT_C(M, 2, 2) = Ne*c;
   MAT_C(M, 3, 2) = 0;
   MAT_C(M, 1, 3) = alpha*Ne2;
   MAT_C(M, 2, 3) = d*Ne2;
   MAT_C(M, 3, 3) = 0;
#if DEBUG
   printf("M = "); d_mat3dc_printf(M); printf("\n");
#endif
   d_mat3dc_invert(L, M);
#if DEBUG
   d_mat3dc_printf(L);printf("\n");
#endif
   d_mat3dc_LLL(Z, L, 0.75);
#if DEBUG
   z_mat3dc_printf(Z);printf("\n\n");
#endif
   d_mat3dc_mul_z_mat3dc(R, L, Z);
   l1 = ceil(48.0/d_vec3d_norm(COL(R, 1)));
   l2 = ceil(48.0/d_vec3d_norm(COL(R, 2)));
   l3 = ceil(48.0/d_vec3d_norm(COL(R, 3)));
   for (long i1 = -l1; i1 <= l1; i1++)
   {
   for (long i2 = -l2; i2 <= l2; i2++)
   {
   for (long i3 = -l3; i3 <= l3; i3++)
   {
      d_vec3d_mul_scalar(v, COL(R,1), (double) i1);
      d_vec3d_add_scalar_mul(v, v, COL(R,2), (double) i2);
      d_vec3d_add_scalar_mul(v, v, COL(R,3), (double) i3);
      t1 = round(COORD(v, 1)*MAT_C(M, 1, 1) + COORD(v, 2)*MAT_C(M, 1, 2) + COORD(v,3)*MAT_C(M, 1, 3));
      if ((t1 != told) && (t1 != - told) && ((t1 > 100) || (t1 < -100)))
      {
         t2 = round(COORD(v, 1)*MAT_C(M, 2, 1) + COORD(v, 2)*MAT_C(M, 2, 2) + COORD(v, 3)*MAT_C(M, 3, 3));
         t3 = round(COORD(v, 1)*MAT_C(M, 3, 1) + COORD(v, 2)*MAT_C(M, 3, 2) + COORD(v, 3)*MAT_C(M, 3, 3));
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
         if ((k >= 2) && (k < 1000)) 
         {
            fprintf(outfile,"\n%ld^3 + %ld^3 + %ld^3  = %ld\n\n", t1, t2, t3, k);
            fflush(outfile);
            told = t1;
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
   mpz_init(tmp);
   mpz_init(tmp2);
   
   outfile = fopen("k3.out", "w");
   N = 100000000000;
   eps = (double)1.0/(double)N;
   for (long k = -10000000000; k > -N+10000000000; k--) 
   {
      if ((k&0xfffff) == 0) 
      {
         fprintf(outfile, "Passed k = %ld\n", k);
         fflush(outfile);
      }
      a = (double)k*eps;
      elkies(a, N, eps);
   }
   for (long k = 10000000000; k < N-10000000000; k++) 
   {
      if ((k&0xfffff) == 0) 
      {
         fprintf(outfile, "Passed k = %ld\n", k);
         fflush(outfile);
      }
      a = k*eps;
      elkies(a, N, eps);
   } 
   
   fclose(outfile);
   mpz_clear(tmp);
   mpz_clear(tmp2);
}
