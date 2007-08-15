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

#define D_COORD(A, x) ((A) + 2*(x-1))
#define D_MAT_C(R, x, y) D_COORD(((R)[y-1]),x)
#define D_COL(A, y) ((A)[y-1])

void elkies(double * a, long N, double * eps)
{
   static long told;
   int test;
   long t1, t2, t3, k, l1, l2, l3;
   
   double b[2];
   double c[2];
   double d[2];
   double alpha[2];
   double Ne[2];
   double Ne2[2];
   
   double temp[2];
   double temp2[2];
   double temp3[2];
   double temp4[2];
   
   c_dd_npwr(a, 3, temp);
   c_dd_add_dd_d(temp, 1.0, temp);
   c_dd_nroot(temp, 3, temp);
   c_dd_neg(temp, b);
   
   c_dd_npwr(b, 6, temp);
   c_dd_npwr(b, 3, temp2);
   c_dd_npwr(a, 3, temp3);
   c_dd_mul(temp2, temp3, temp2);
   c_dd_mul_dd_d(temp2, 6.0, temp2);
   c_dd_add(temp, temp2, temp);
   c_dd_npwr(a, 6, temp2);
   c_dd_mul_dd_d(temp2, 5.0, temp2);
   c_dd_add(temp, temp2, temp);
   c_dd_npwr(b, 6, temp2);
   c_dd_mul(temp2, a, temp2);
   c_dd_mul_dd_d(temp2, 6.0, temp2);
   c_dd_npwr(a, 4, temp3);
   c_dd_npwr(b, 3, temp4);
   c_dd_mul(temp3, temp4, temp3);
   c_dd_mul_dd_d(temp3, 6.0, temp3);
   c_dd_add(temp3, temp2, temp2);
   c_dd_div(temp, temp2, alpha);
   c_dd_neg(alpha, alpha);
   
   c_dd_mul(a, a, temp);
   c_dd_mul(b, b, temp2);
   c_dd_div(temp, temp2, c);
   c_dd_neg(c, c);
   
   c_dd_mul(b, b, temp);
   c_dd_mul_dd_d(temp, 6.0, temp);
   c_dd_div(a, temp, temp);
   c_dd_mul_dd_d(temp, 5.0, temp);
   c_dd_npwr(b, 5, temp2);
   c_dd_mul_dd_d(temp2, 6.0, temp2);
   c_dd_npwr(a, 4, temp3);
   c_dd_div(temp3, temp2, temp2);
   c_dd_add(temp, temp2, d);
   c_dd_neg(d, d);
   
   c_dd_mul_dd_d(eps, N, Ne);
   
   c_dd_mul(Ne, eps, Ne2);
   
   dd_vec3d v;
   dd_mat3dc_t L, M, R;
   z_mat3dc_t Z;
   dd_mat3dc_stack_init(&M);
   dd_mat3dc_stack_init(&L);
   dd_mat3dc_stack_init(&R);
   z_mat3dc_stack_init(&Z);
   dd_vec3d_stack_init(&v);
   
   c_dd_mul_dd_d(a, N, D_MAT_C(M, 1, 1));
   c_dd_mul_dd_d(b, N, D_MAT_C(M, 2, 1));
   c_dd_copy_d(N, D_MAT_C(M, 3, 1));
   c_dd_copy(Ne, D_MAT_C(M, 1, 2));
   c_dd_mul(Ne, c, D_MAT_C(M, 2, 2));
   c_dd_copy_d(0.0, D_MAT_C(M, 3, 2));
   c_dd_mul(alpha, Ne2, D_MAT_C(M, 1, 3));
   c_dd_mul(d, Ne2, D_MAT_C(M, 2, 3));
   c_dd_copy_d(0.0, D_MAT_C(M, 3, 3));
   
#if DEBUG
   printf("M = "); dd_mat3dc_printf(M); printf("\n");
#endif
   dd_mat3dc_invert(L, M);
#if DEBUG
   printf("M^(-1) = "); dd_mat3dc_printf(L);printf("\n");
#endif
   if (dd_mat3dc_LLL(Z, R, L, 0.9999))
   {
#if DEBUG
   printf("Z = "); z_mat3dc_printf(Z);printf("\n\n");
#endif
#if DEBUG
   printf("R = "); dd_mat3dc_printf(R);printf("\n\n");
#endif
   dd_vec3d_norm(temp, D_COL(R, 1));
   c_dd_div_d_dd(3.0, temp, temp);
   c_dd_ceil(temp, temp);
   l1 = (long) temp[0];
   dd_vec3d_norm(temp, D_COL(R, 2));
   c_dd_div_d_dd(3.0, temp, temp);
   c_dd_ceil(temp, temp);
   l2 = (long) temp[0];
   dd_vec3d_norm(temp, D_COL(R, 3));
   c_dd_div_d_dd(3.0, temp, temp);
   c_dd_ceil(temp, temp);
   l3 = (long) temp[0];
   if ((l1 > 100) || (l1 < 0)) l1 = 100;
   if ((l2 > 100)  || (l2 < 0)) l2 = 100;
   if ((l3 > 100)  || (l3 < 0)) l3 = 100;
   for (long i1 = -l1; i1 <= l1; i1++)
   {
   for (long i2 = -l2; i2 <= l2; i2++)
   {
   for (long i3 = -l3; i3 <= l3; i3++)
   {
      if ((!i2) && (!i3)) i3++;
      dd_vec3d_mul_scalar_d(v, D_COL(R,1), (double) i1);
      dd_vec3d_add_scalar_mul_d(v, v, D_COL(R,2), (double) i2);
      dd_vec3d_add_scalar_mul_d(v, v, D_COL(R,3), (double) i3);
      dd_vec3d_norm(temp, v);
      c_dd_comp_dd_d(temp, 3.0, &test);
      if (test < 0) 
      {
      c_dd_mul(D_COORD(v, 1), D_MAT_C(M, 1, 1), temp);
      c_dd_mul(D_COORD(v, 2), D_MAT_C(M, 1, 2), temp2);
      c_dd_add(temp, temp2, temp);
      c_dd_mul(D_COORD(v, 3), D_MAT_C(M, 1, 3), temp2);
      c_dd_add(temp, temp2, temp);
      c_dd_nint(temp, temp2);
      t1 = (long) temp2[0];
      if ((t1 != told) && (t1 != - told) && ((t1 > 100) || (t1 < -100)))
      {
         c_dd_mul(D_COORD(v, 1), D_MAT_C(M, 2, 1), temp);
         c_dd_mul(D_COORD(v, 2), D_MAT_C(M, 2, 2), temp2);
         c_dd_add(temp, temp2, temp);
         c_dd_mul(D_COORD(v, 3), D_MAT_C(M, 2, 3), temp2);
         c_dd_add(temp, temp2, temp);
         c_dd_nint(temp, temp2);
         t2 = (long) temp2[0];
         
         c_dd_mul(D_COORD(v, 1), D_MAT_C(M, 3, 1), temp);
         c_dd_mul(D_COORD(v, 2), D_MAT_C(M, 3, 2), temp2);
         c_dd_add(temp, temp2, temp);
         c_dd_mul(D_COORD(v, 3), D_MAT_C(M, 3, 3), temp2);
         c_dd_add(temp, temp2, temp);
         c_dd_nint(temp, temp2);
         t3 = (long) temp2[0];
         
         mpz_set_si(tmp, t1);
         mpz_pow_ui(tmp, tmp, 3);
         mpz_set_si(tmp2, t2);
         mpz_pow_ui(tmp2, tmp2, 3);
         mpz_add(tmp, tmp, tmp2);
         mpz_set_si(tmp2, t3);
         mpz_pow_ui(tmp2, tmp2, 3);
         mpz_add(tmp, tmp, tmp2);
         mpz_abs(tmp, tmp);
         if (mpz_size(tmp) == 1) k = mpz_get_ui(tmp);
         else k = 0;
         if ((k >= 2) && (k < 10000) && ((fabs(t1) > 1000) || (fabs(t2) > 1000) || (fabs(t3) > 1000))) 
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
   dd_vec3d_stack_clear();
   z_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
   dd_mat3dc_stack_clear();
}

int main(void)
{
   long N;
   double eps[2];
   double a[2];
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
   
   outfile = fopen("x3y3z3-test.out", "w");
   N = 100000000000;
   c_dd_copy_d(1.0, eps);
   c_dd_div_dd_d(eps, N, eps); // eps = 1/N
   /*for (long k = -50000000000; k > -N+20000000000; k--) 
   {
      if ((k&0x7ffff) == 0) 
      {
         fprintf(outfile, "Passed k = %ld\n", k);
         fflush(outfile);
      }
      c_dd_copy_d(k, a);
      c_dd_div_dd_d(a, N, a); // a = k/N
      elkies(a, N, eps);
   }*/
   for (long k = 50000000000; k < N-20000000000; k++) 
   {
      if ((k&0x7ffff) == 0) 
      {
         fprintf(outfile, "Passed k = %ld\n", k);
         fflush(outfile);
      }
      c_dd_copy_d(k, a);
      c_dd_div_dd_d(a, N, a); // a = k/N
      elkies(a, N, eps);
   } 
   
   fclose(unknown);
   flint_stack_release();
   fclose(outfile);
   mpz_clear(tmp);
   mpz_clear(tmp2);
}
