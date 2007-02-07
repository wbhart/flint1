#include <stdio.h>
#include <gmp.h>
#include "Z.h"
#include "Zvec.h"

#define DEGREE 56 
#define LGDEG 6
#define BITS 135
#define ITERATIONS 10

int main(void)
{
    unsigned long num1, num2, temp2;
    
    gmp_randstate_t state;
    gmp_randinit_default(state);
    
    mpz_t random;
    mpz_init(random);
    
    mpz_t n1, n, temp;
    mpz_init(n1); mpz_init(n); mpz_init(temp);
    
    initZvec();
    
    
    //printf("Enter value 1: ");scanf("%ld",&num1);getchar();
    //printf("Enter value 2: ");scanf("%ld",&num2);getchar();
    //printf("gcd(%ld,%ld) = %ld\n\n",num1,num2,gcd(num1,num2));
    
    //for (int i = 0; i<10000; i++)
    //{
    //    randomprime(random,100);
        //gmp_printf("%Zd\n",random);
    //}
    
    /*for (long int i = 0; i<10000000; i++)
    {
        mpz_urandomb(random,state,30);
        num1 = mpz_get_ui(random);
        mpz_urandomb(random,state,30);
        num2 = mpz_get_ui(random);
        //printf("%ld\n",gcd(num1,num2));
    }*/
    
    /*mpz_set_ui(n1,39);
    mpz_ui_pow_ui(n1,10,mpz_get_ui(n1));
    mpz_nextprime(n,n1);
    mpz_set_ui(n1,41);
    mpz_ui_pow_ui(n1,10,mpz_get_ui(n1));
    mpz_nextprime(n1,n1);
    mpz_mul(n,n,n1);
    
    for (long int i = 0; i<1000000; i++)
    {
        randomprime(temp,20);
        mpz_fdiv_r(n1,n,temp);
        sqrtmod(n1,n1,temp);
    }*/
    
    /*for (long int i = 0; i<10000000; i++)
    {
        temp2 = i%100;
        power(n,temp2,5);
        power(n,temp2,7);
        power(n,temp2,19);
        power(n,temp2,34);
        power(n,temp2,42);
    }*/
    
    /*mpz_set_ui(temp,2);
    
    for (long i=0; i<1000000; i++)
    {
       mpz_pow_ui(n,temp,100);
       mpz_set_ui(temp,2);
       mpz_set(n1,n);
       while (mpz_cmp_ui(n1,1)!=0)
       {
          if (mpz_divisible_p(n1,temp)) mpz_remove(n1,n1,temp);
          mpz_nextprime(temp,temp);
       }
       randomprime(temp,10);
    }*/
    
    Zvec aVec;
    Zvec_init3(aVec,DEGREE+1,BITS/64+1);
    Zvec bVec;
    Zvec_init3(bVec,DEGREE+1,BITS/64+1);
    Zvec cVec;
    Zvec_init3(cVec,DEGREE+DEGREE+1,(2+LGDEG+2*BITS)/64+1);
    
    /*printf("Start Radix\n");
    for (int j =0; j<ITERATIONS; j++)
    {
       if (j%20==0)
       {
         for (unsigned long i = 0; i<DEGREE+1; i++)
         {
          mpz_urandomb(aVec.coords[i],state,BITS);
          mpz_urandomb(bVec.coords[i],state,BITS); 
         }*/
         
         /*printf("[");
         for (unsigned long i = 0; i<DEGREE; i++)
         {
             gmp_printf("%Zd, ",aVec.coords[i]);
         }
         gmp_printf("%Zd]\n",aVec.coords[DEGREE]);
         printf("[");
         for (unsigned long i = 0; i<DEGREE; i++)
         {
             gmp_printf("%Zd, ",bVec.coords[i]);
         }
         gmp_printf("%Zd]\n",bVec.coords[DEGREE]);*/    
       //}
       //Zvec_clear(cVec);
        //Zvec_karamul(cVec,aVec,bVec,24); 
           
           //bound = 2 + NumBits(min(da, db)+1) + MaxBits(a) + MaxBits(b);
       //Zvec_radixMul(cVec,aVec,bVec,2+LGDEG+2*BITS);
       
       /*printf("[");
       for (unsigned long i = 0; i<2*DEGREE; i++)
       {
             gmp_printf("%Zd, ",cVec.coords[i]);
       }
       gmp_printf("%Zd]\n",cVec.coords[2*DEGREE]);*/
    //}
    
    printf("Start SSMul\n");
    
    for (int j =0; j<ITERATIONS; j++)
    {
       if (j%20==0)
       {
         for (unsigned long i = 0; i<DEGREE+1; i++)
         {
          mpz_urandomb(aVec.coords[i],state,BITS);
          mpz_urandomb(bVec.coords[i],state,BITS); 
         }
         /*printf("[");
         for (unsigned long i = 0; i<DEGREE; i++)
         {
             gmp_printf("%Zd, ",aVec.coords[i]);
         }
         gmp_printf("%Zd]\n",aVec.coords[DEGREE]);
         printf("[");
         for (unsigned long i = 0; i<DEGREE; i++)
         {
             gmp_printf("%Zd, ",bVec.coords[i]);
         }
         gmp_printf("%Zd]\n",bVec.coords[DEGREE]);*/
       }
       //Zvec_radixMul(cVec,aVec,bVec,175);
       //Zvec_karamul(cVec,aVec,bVec,4688);
       Zvec_SSMul(cVec,aVec,bVec,BITS,1);
       //Zvec_KSMul(cVec,aVec,bVec,BITS);
       /*printf("[");
       for (unsigned long i = 0; i<2*DEGREE; i++)
       {
             gmp_printf("%Zd, ",cVec.coords[i]);
       }
       gmp_printf("%Zd]\n\n",cVec.coords[2*DEGREE]); */
    }
      
    printf("Stop\n");
      
    return 0;
}

