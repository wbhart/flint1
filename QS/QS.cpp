/*============================================================================
    Copyright 2006 William Hart    

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

Description:

This is a relatively fast implementation of the self-initialising quadratic sieve.
If you manage to improve the code, the author would like to hear about it.

Contact: hart_wb {at-thingy} yahoo.com
================================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <string.h>
#include <sys/times.h>
#include <limits.h>

#include "TonelliShanks.h"
#include "ModuloArith.h"
#include "F2matrix.h"
#include "lanczos.h"

//===========================================================================
//Uncomment these for various pieces of debugging information

#define COUNT    // Shows the number of relations generated and curves used during sieving
//#define RELS     // Shows the actual factorizations of the relations
//#define ERRORS   // Error if relation should be divisible by a prime but isn't 
//#define POLS     // Shows the polynomials being used by the sieve
//#define ADETAILS // Prints some details about the factors of the A coefficients of the polys
//#define LARGESTP // Prints the size of the largest factorbase prime
//#define CURPARTS // Prints the number of curves used and number of partial relations
//#define TIMING //displays some relative timings, if feature is available
//#define REPORT //report sieve size, multiplier and number of primes used

//===========================================================================
//Architecture dependent fudge factors

#if ULONG_MAX == 4294967295U
#define SIEVEMASK 0xC0C0C0C0U
#define MIDPRIME 1500
#define SIEVEDIV 1
#elif ULONG_MAX == 18446744073709551615U
#define SIEVEMASK 0xC0C0C0C0C0C0C0C0U
#define MIDPRIME       6000 
#define SIEVEDIV 2
#endif

#define CACHEBLOCKSIZE 64000 //Should be a little less than the L1/L2 cache size
                             //and a multiple of 64000
#define MEDIUMPRIME    900   
#define SECONDPRIME    6000 //This should be lower for slower machines
#define FUDGE          0.15 //Every program needs a mysterious fudge factor

//===========================================================================
                      
#define LARGEPRIME 6000000 //This is of course irrelevant since no partial relations are collected 

#define MINDIG 40 //Will not factor numbers with less than this number of decimal digits

#define PREFETCH(addr,n) __builtin_prefetch((unsigned long*)addr+n,0,1)

//===========================================================================
//Knuth-Schroeppel multipliers and a macro to count them

static const u_int32_t multipliers[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 
                                                23, 29, 31, 37, 41, 43};

#define NUMMULTS (sizeof(multipliers)/sizeof(u_int32_t))

//============================================================================
// Number of primes to use in factor base, given the number of decimal digits specified
u_int32_t primesNo[] = 
{
     1500, 1500, 1600, 1700, 1700, 1800, 1900, 2000, 2000, 2100, //40-49
     2200, 2300, 2500, 2700, 2900, 3200, 3500, 3800, 4100, 4400, //50-59
     4700, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, //60-69
     9500, 10000, 11500, 13000, 15000, //70-74
     17000, 24000, 24000, 30000, 37000, //75-79
     45000, 55000, 56000, 57000, 58000,  //80-84
     59000, 60000, 64000, 68000, 72000,  //85-89
     76000, 80000 //90-91
};

//============================================================================
// First prime actually sieved for
u_int32_t firstPrimes[] = 
{
     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, //40-49
     9, 9, 9, 9, 9, 9, 10, 10, 10, 10, //50-59
     10, 10, 11, 11, 12, 12, 13, 13, 14, 15, //60-69
     16, 17, 18, 19, 20, //70-74
     21, 22, 22, 23, 24, //75-79
     24, 25, 26, 27, 28, //80-84
     29, 30, 31, 32, 33, //85-89
     34, 35 //90-91
};

//============================================================================
// Logs of primes are rounded and errors accumulate; this specifies how great an error to allow
u_int32_t errorAmounts[] = 
{
     12, 13, 13, 13, 13, 14, 14, 14, 14, 14, //40-49
     14, 15, 15, 15, 16, 16, 16, 17, 17, 17, //50-59
     17, 17, 17, 18, 18, 18, 19, 18, 19, 19, //60-69
     19, 19, 20, 20, 20, //70-74
     21, 21, 22, 22, 22, //75-79
     22, 22, 23, 23, 23, //80-84
     24, 24, 24, 25, 25, //85-89
     25, 25 //90-91
};

//============================================================================
// This is the threshold the sieve value must exceed in order to be considered for smoothness
u_int32_t thresholds[] = 
{
     66, 67, 68, 69, 70, 70, 71, 72, 73, 74, //40-49
     75, 76, 77, 78, 79, 79, 80, 81, 82, 83, //50-59
     84, 84, 85, 86, 87, 88, 89, 90, 90, 91, //60-69
     91, 92, 92, 93, 94, //70-74
     95, 96, 97, 98, 98, //75-79
     99, 100, 100, 101, 102, //80-84
     103, 104, 104, 105, 106, //85-89
     107, 108 //90-91  
};

//============================================================================
// Usually less than the full number of relations is required. This specifies how many fewer are needed
u_int32_t lessrels[] = 
{
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, //40-49
     50, 50, 50, 50, 50, 50, 50, 50, 50, 50, //50-59
     150, 150, 150, 150, 150, 150, 150, 150, 150, 150, //60-69
     400, 400, 400, 400, 400, //70-74
     600, 600, 700, 800, 9000, //75-79
     1400, 2000, 2000, 2100, 2100, //80-81
     2200, 2200, 2300, 2400, 2500, //85-89
     2600, 2700 //90-91
};

//============================================================================
// Size of sieve to use divided by 2, given the number of decimal digits specified
//N.B: probably optimal if chosen to be a multiple of 32000, though other sizes are supported
u_int32_t sieveSize[] = 
{
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //40-49
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //50-59
     64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, //60-69
     64000, 64000, 64000, 64000, 64000, //70-74
     128000, 128000, 128000, 128000, 128000, //75-79
     160000, 160000, 160000, 160000, 160000, //80-84
     192000, 192000, 192000, 192000, 192000, //85-89
     192000, 192000 //90-91
};

//============================================================================
int32_t decdigits; //number of decimal digits of n
u_int32_t secondprime; //min(numprimes, SECONDPRIME) = cutoff for using flags when sieving
u_int32_t firstprime;  //first prime actually sieved with
unsigned char errorbits;  //first prime actually sieved with
unsigned char threshold;  //sieve threshold cutoff for smooth relations
u_int32_t less;  //number of relations less than numPrimes that will be needed
u_int32_t midprime;

u_int32_t * factorBase; //array of factor base primes
u_int32_t numPrimes; //number of primes in factor base
u_int32_t relSought; //number of relations sought, i.e. a "few" more than numPrimes
unsigned char * primeSizes; //array of sizes in bits, of the factor base primes
unsigned char * sieve; //actual array where sieving takes place
unsigned char * * offsets; //offsets for each prime to use in sieve 
unsigned char * * offsets2; //offsets for each prime to use in sieve (we switch between these)
u_int32_t relsFound =0; //number of relations found so far
unsigned char * flags; //flags used for speeding up sieving for large primes
u_int32_t partials = 0; //number of partial relations
u_int32_t Mdiv2; //size of sieving interval divide 2 
u_int32_t mat2off; //offset of second square block in matrix

mpz_t * sqrts; //square roots of n modulo each prime in the factor base

mpz_t n; //number to be factored 
mpz_t res; //smooth values which are trial factored

mpz_t temp, temp2, temp3; //temporary variables
mpz_t q,r; //quotient and remainder

//Variables used for keeping time

u_int32_t clockstart;
u_int32_t clocktotal = 0;  

//Variables used by the modular inversion macro function
int32_t u1, u3;
int32_t v1, v3;
int32_t t1, t3, quot;

//Variable used for random function
u_int32_t randval = 2994439072;

//==========================================================================================
//Timing: provides some relative timings on X86 machines running gcc

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))

#ifdef TIMING
#define TIMES
#endif

double counterfirst[4];
double countertotal[4] = {0,0,0,0};

static unsigned counthi = 0;
static unsigned countlo = 0;

void counterasm(unsigned *hi, unsigned *lo)
{
 asm("rdtsc; movl %%edx,%0; movl %%eax,%1" 
 : "=r" (*hi), "=r" (*lo) 
 : 
 : "%edx", "%eax");
}

double getcounter()
{
   double total;

   counterasm(&counthi, &countlo);

   total = (double) counthi * (1 << 30) * 4 + countlo;
   return total;
}

#endif
   
/*========================================================================
   Modular Inversion:

   Function: GMP has a modular inverse function, but believe it or not, 
             this clumsy implementation is apparently quite a bit faster. 
             It inverts the value a, modulo the prime p, using the extended 
             gcd algorithm.

========================================================================*/

inline u_int32_t modinverse(u_int32_t a, u_int32_t p)
{
   u1=1; u3=a;
   v1=0; v3=p;
   t1=0; t3=0;
   while (v3)
   {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               t1 = u1; u1 = v1; v1 = t1;
               t3 = u3; u3 = v3; v3 = t3;
            } else 
            {
               t1 = u1 - v1; u1 = v1; v1 = t1;
               t3 = u3 - v3; u3 = v3; v3 = t3;
            }
         } else if (quot < (v3<<1))
         {  
            t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
            t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
         } else
         {
            t1 = u1 - v1*3; u1 = v1; v1 = t1;
            t3 = u3 - v3*3; u3 = v3; v3 = t3;
         }
      } else
      {
         quot=u3/v3;
         t1 = u1 - v1*quot; u1 = v1; v1 = t1;
         t3 = u3 - v3*quot; u3 = v3; v3 = t3;
      }
   } 
   
   if (u1<0) u1+=p;
   
   return u1;
}

/*=========================================================================
   Knuth_Schroeppel Multiplier:
 
   Function: Find the best multiplier to use (allows 2 as a multiplier).
             The general idea is to find a multiplier k such that kn will
             be faster to factor. This is achieved by making kn a square 
             modulo lots of small primes. These primes will then be factor
             base primes, and the more small factor base primes, the faster
             relations will accumulate, since they hit the sieving interval
             more often. 
 
==========================================================================*/
u_int32_t knuthSchroeppel(mpz_t n)
{
    float bestFactor = -10.0f;
    u_int32_t multiplier = 1;
    u_int32_t nmod8;
    float factors[NUMMULTS];
    float logpdivp;
    mpz_t prime, r, mult;
    int32_t kron, multindex;
    
    mpz_init(prime);
    mpz_init(r);
    mpz_init(mult);
    
    nmod8 = mpz_fdiv_r_ui(r,n,8);
    
    for (multindex = 0; multindex < NUMMULTS; multindex++)
    {
       int32_t mod = nmod8*multipliers[multindex]%8;
       factors[multindex] = 0.34657359; // ln2/2 
       if (mod == 1) factors[multindex] *= 4.0;   
       if (mod == 5) factors[multindex] *= 2.0;   
       factors[multindex] -= (log((float) multipliers[multindex]) / 2.0);
    }
    
    mpz_set_ui(prime,3);
    while (mpz_cmp_ui(prime,10000)<0)
    {
          logpdivp = log((float)mpz_get_ui(prime)) / mpz_get_ui(prime);
          kron = mpz_kronecker(n,prime);
          for (multindex = 0; multindex < NUMMULTS; multindex++)
          {
              mpz_set_ui(mult,multipliers[multindex]);
              switch (kron*mpz_kronecker(mult,prime))
              {
                 case 0:
                 {
                      factors[multindex] += logpdivp;
                 } break;
                 case 1:
                 {
                      factors[multindex] += 2.0*logpdivp;
                 } break;
                 default: break;
              }
          }
          
          mpz_nextprime(prime,prime);
    }
    
    for (multindex=0; multindex<NUMMULTS; multindex++)
    {
      if (factors[multindex] > bestFactor)
      { 
        bestFactor = factors[multindex];
        multiplier = multipliers[multindex];
      }
    } 
    
    mpz_clear(prime);
    mpz_clear(r);
    mpz_clear(mult);
    
    return multiplier;
}



/*========================================================================
   Initialize Quadratic Sieve:
  
   Function: Initialises the global gmp variables.

========================================================================*/
void initSieve(void)
{
    mpz_init(n);
    mpz_init(temp); 
    mpz_init(temp2);
    mpz_init(temp3);
    mpz_init(res);
    mpz_init(q);
    mpz_init(r);
    
    return;
}

/*========================================================================
   Compute Factor Base:
 
   Function: Computes primes p up to B for which n is a square mod p,  
   allocates memory and stores them in an array pointed to by factorBase
   Returns: number of primes actually in the factor base

========================================================================*/
void computeFactorBase(mpz_t n, u_int32_t B,u_int32_t multiplier)
{
     mpz_t currentPrime;
     u_int32_t primesinbase = 0;
     
     factorBase = (u_int32_t *) calloc(sizeof(u_int32_t),B); 
     
     factorBase[primesinbase] = multiplier;
     primesinbase++;
     if (multiplier!=2)
     {
        factorBase[primesinbase] = 2;
        primesinbase++;
     }
     mpz_init_set_ui(currentPrime,3);
     while (primesinbase < B)
     {
          if (mpz_kronecker(n,currentPrime)==1)
          {
              factorBase[primesinbase] = mpz_get_ui(currentPrime);
              primesinbase++;
          } 
          mpz_nextprime(currentPrime,currentPrime);
     }
#ifdef LARGESTP
     gmp_printf("Largest prime less than %Zd\n",currentPrime);
#endif
      
     mpz_clear(currentPrime);
     return;
}

/*===========================================================================
   Compute Prime Sizes:
 
   Function: Computes the size in bits of each prime in the factor base
     allocates memory for an array, primeSizes, to store the sizes
     stores the size for each of the numPrimes primes in the array 
 
===========================================================================*/
void computeSizes(u_int32_t numPrimes)
{
     primeSizes = (unsigned char *) calloc(sizeof(unsigned char),numPrimes);
     for (u_int32_t i = 0; i<numPrimes; i++)
     {
         primeSizes[i]=(unsigned char)floor(log(factorBase[i])/log(2.0)-FUDGE+0.5);
     }
     
     return;
}

/*===========================================================================
   Tonelli-Shanks:

   Function: Performs Tonelli-Shanks on n mod every prime in the factor base
      allocates memory for the results to be stored in the array sqrts

===========================================================================*/
void tonelliShanks(u_int32_t numPrimes,mpz_t n)
{
     sqrts = (mpz_t *) calloc(sizeof(mpz_t),numPrimes); 
     mpz_array_init(sqrts[0],numPrimes,8*sizeof(u_int32_t));
     
     for (u_int32_t i = 1; i<numPrimes; i++) 
     {
         mpz_set_ui(temp,factorBase[i]);
         sqrtmod(sqrts[i],n,temp);
     }
     
     return;
}

/*==========================================================================
   insertColEntry:

   Function: insert an entry into a column of the matrix, 
   reallocating the space for the column if necessary

===========================================================================*/
inline void insertColEntry(la_col_t* colarray, u_int32_t colNum, u_int32_t entry)
{
   u_int32_t* temp;
   
       if (((colarray[colNum].weight>>4)<<4)==colarray[colNum].weight) //need more space
       {
           temp = colarray[colNum].data;
           colarray[colNum].data = (u_int32_t*)malloc((colarray[colNum].weight+16)*sizeof(u_int32_t));
           for (u_int32_t i = 0; i<colarray[colNum].weight; i++)
           {
               colarray[colNum].data[i] = temp[i];
           }
           if (colarray[colNum].weight!=0) free(temp);
       }
   
   colarray[colNum].data[colarray[colNum].weight] = entry;
   colarray[colNum].weight++;
   colarray[colNum].orig = colNum;
}

/*==========================================================================
   xorColEntry:

   Function: xor entry corresponding to a prime dividing A, which will be 
   either the last entry in the column, or not there at all, so we either
   add it in or take it out

===========================================================================*/
inline void xorColEntry(la_col_t* colarray, u_int32_t colNum, u_int32_t entry)
{
   if ((colarray[colNum].weight !=0) && (colarray[colNum].data[colarray[colNum].weight-1] == entry)) colarray[colNum].weight--;
   else insertColEntry(colarray,colNum,entry);
}

/*==========================================================================
   clearCol:

   Function: clear a column
   
===========================================================================*/
inline void clearCol(la_col_t* colarray, u_int32_t colNum)
{
   colarray[colNum].weight =0;
}

/*==========================================================================
   evaluateSieve:

   Function: searches sieve for relations and sticks them into a matrix, then
             sticks their X and Y values into two arrays XArr and YArr

===========================================================================*/
void evaluateSieve(u_int32_t ** relations, u_int32_t ctimesreps, u_int32_t M, unsigned char * sieve, mpz_t A, mpz_t B, mpz_t C, u_int32_t * soln1, u_int32_t * soln2, int32_t polyadd, u_int32_t * polycorr, matrix m, mpz_t * XArr, u_int32_t * aind, int32_t min, int32_t s,u_int32_t multiplier, int32_t * exponents, la_col_t* colarray)
{
     long i,j;
     register u_int32_t k;
     u_int32_t exponent, vv;
     unsigned char extra;
     register u_int32_t modp;
     unsigned long * sieve2;
     unsigned char bits;
     int32_t numfactors;
     
     i = 0;
     j=0;
     sieve2 = (unsigned long *) sieve;
#ifdef POLS
     gmp_printf("%Zdx^2%+Zdx\n%+Zd\n",A,B,C);
#endif
     
     while (j<M/sizeof(unsigned long))
     {
        do
        {
           while (!(sieve2[j] & SIEVEMASK)) j++;
           i=j*sizeof(unsigned long);
           j++;
           while ((i<j*sizeof(unsigned long))&&(sieve[i] < threshold)) i++;
        } while (sieve[i] < threshold);
           
        if ((i<M) && (relsFound < relSought)) 
        {
           mpz_set_ui(temp,i+ctimesreps);
           mpz_sub_ui(temp,temp,Mdiv2);
              
           mpz_set(temp3,B);  //B
           mpz_addmul(temp3,A,temp);  //AX+B
           mpz_add(temp2,temp3,B);  //AX+2B
           mpz_mul(temp2,temp2,temp);  //AX^2+2BX
           mpz_add(res,temp2,C);  //AX^2+2BX+C
              
           bits=mpz_sizeinbase(res,2);
           bits-=errorbits;
              
           numfactors=0;
              
           extra = 0;
           if (factorBase[0]!=1)
           {
              mpz_set_ui(temp,factorBase[0]);
              exponent = mpz_remove(res,res,temp);
              if (exponent) 
              { 
                 extra+=primeSizes[0];
                 for (int32_t i = 0; i<exponent; i++) relations[relsFound][++numfactors] = 0;
              }
              if (exponent & 1) 
              {
                 insertEntry(m,0,0);
                 //insertColEntry(colarray,relsFound,0); 
              }
           }
             
           mpz_set_ui(temp,factorBase[1]);
           exponent = mpz_remove(res,res,temp);
           if (exponent) 
           { 
              extra+=primeSizes[1];
              for (int32_t i = 0; i<exponent; i++) relations[relsFound][++numfactors] = 1;
           }
           if (exponent & 1) 
           {
              insertEntry(m,0,1);
              //insertColEntry(colarray,relsFound,1);  
           }
                
           for (k = 2; k<firstprime; k++)
           {
              modp=(i+ctimesreps)%factorBase[k];
                
              if (soln2[k]!=0xFFFFFFFFl)
              {
                 if ((modp==soln1[k]) || (modp==soln2[k]))
                 {
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
             
#ifdef ERRORS
                    if (exponent==0) printf("Error!\n");
#endif
                    extra+=primeSizes[k];
#ifdef RELS
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%d",exponent);
#endif
                    exponents[k] = exponent;
                    if (exponent&1) 
                    {
                       insertEntry(m,0,k);
                       //insertColEntry(colarray,relsFound,k); 
                    }
                 } else exponents[k] = 0;
              } else
              {
                 mpz_set_ui(temp,factorBase[k]);
                 exponent = mpz_remove(res,res,temp);
                 if (exponent) extra+=primeSizes[k];
#ifdef RELS
                 if (exponent > 0) gmp_printf(" %Zd",factorBase[k]);
                 if (exponent > 1) printf("^%d",exponent);
#endif
                 exponents[k] = exponent;
                 if (exponent &1) 
                 {
                    insertEntry(m,0,k);
                    //insertColEntry(colarray,relsFound,k);
                 } 
              }  
           }  
           sieve[i]+=extra;
           if (sieve[i] >= bits)
           {
              vv=((unsigned char)1<<(i&7));
              for (k = firstprime; (k<secondprime)&&(extra<sieve[i]); k++)
              {
                 modp=(i+ctimesreps)%factorBase[k];
                 if (soln2[k]!=0xFFFFFFFFl)
                 {
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       extra+=primeSizes[k];
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
              
#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif
                        
#ifdef RELS
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%d",exponent);
#endif
                       if (exponent) 
                       { 
                          for (int32_t i = 0; i<exponent; i++) relations[relsFound][++numfactors] = k;
                       }
                       if (exponent&1) 
                       {
                          insertEntry(m,0,k);
                          //insertColEntry(colarray,relsFound,k); 
                       }
                    }  
                 } else
                 {
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
                    if (exponent) extra+=primeSizes[k];
                        
#ifdef RELS
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%d",exponent);
#endif
                    if (exponent) 
                    { 
                       for (int32_t i = 0; i<exponent; i++) relations[relsFound][++numfactors] = k;
                    }
                    if (exponent &1) 
                    {
                       insertEntry(m,0,k); 
                       //insertColEntry(colarray,relsFound,k); 
                    }
                 }  
              }  
              
              
              for (k = secondprime; (k<numPrimes)&&(extra<sieve[i]); k++)
              {
                 if (flags[k]&vv)
                 {
                    modp=(i+ctimesreps)%factorBase[k];
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif
                       extra+=primeSizes[k];
#ifdef RELS
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%d",exponent);
#endif
                       if (exponent) 
                       { 
                          for (int32_t i = 0; i<exponent; i++) relations[relsFound][++numfactors] = k;
                       }
                       if (exponent&1) 
                       {
                          insertEntry(m,0,k);
                          //insertColEntry(colarray,relsFound,k); 
                       }
                    }  
                 }  
              }  
              
              for (int32_t i =0; i<s; i++)
              {
                 xorEntry(m,0,aind[i]+min);
                 relations[relsFound][++numfactors] = aind[i]+min;
              }
              
              for (u_int32_t i=0; i<numPrimes; i++)
              {
                  if (getEntry(m,0,i)) insertColEntry(colarray,relsFound,i);
              }
              
              if (mpz_cmp_ui(res,1000)>0)
              {
                 if (mpz_cmp_ui(res,LARGEPRIME)<0) 
                 {
                    partials++;
                 }
                 clearRow(m,numPrimes,0);
                 clearCol(colarray,relsFound);
#ifdef RELS
                 gmp_printf(" %Zd\n",res);
#endif
              } else
              { 
                 mpz_neg(res,res);
                 if (mpz_cmp_ui(res,1000)>0)
                 {
                    if (mpz_cmp_ui(res,LARGEPRIME)<0) 
                    {
                       partials++;
                    }
                    clearRow(m,numPrimes,0);
                    clearCol(colarray,relsFound);
#ifdef RELS
                    gmp_printf(" %Zd\n",res);
#endif
                 } else 
                 {
#ifdef RELS
                    printf("....R\n");
#endif
                    for (int32_t i = 1; i<firstprime; i++) 
                    {
                       for (int32_t j = 0; j < exponents[i]; j++) relations[relsFound][++numfactors] = i;
                    }
                    relations[relsFound][0] = numfactors;
                       
                    mpz_init(XArr[relsFound]);
                    mpz_set(XArr[relsFound],temp3); //(AX+B)
                       
                    relsFound++;
                    clearRow(m,numPrimes,0);
#ifdef COUNT
                    if (relsFound%100==0) fprintf(stderr,"%ld relations, %ld partials.\n",(long)relsFound,(long)partials);
#endif
                 }  
              }  
           } else 
           {
              clearRow(m,numPrimes,0);
              clearCol(colarray,relsFound);
#ifdef RELS
              printf("\r                                                                    \r");
#endif
                 
           }   
              i++;
              
        } else if (relsFound >= relSought) i++;      
     }  
     
     return;
}


/*=============================================================================
   Sieve:

   Function: Allocates space for a sieve of M integers and sieves the interval
             starting at start

=============================================================================*/
void sieveInterval(u_int32_t M, u_int32_t numPrimes, unsigned char * sieve, int32_t last, int32_t first, int32_t polyadd, u_int32_t * soln1, u_int32_t * soln2, u_int32_t * polycorr, unsigned char * * offsets, unsigned char * * offsets2)
{
     register unsigned char currentprimesize; 
     register u_int32_t currentprime;
     unsigned char * position2;
     register unsigned char * position;
     register long diff;
     unsigned char * end;
     u_int32_t ptimes4;
     int32_t correction;
     
     end = sieve+M;
     
     if (first)
     {
        for (u_int32_t prime=1; prime<firstprime; prime++) 
        {
            if (soln2[prime] == 0xFFFFFFFF) continue;
            currentprime = factorBase[prime];
            correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
            soln1[prime]+=correction;
            while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
            soln2[prime]+=correction;
            while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
        }  
     }
     
     for (u_int32_t prime=firstprime; prime<MEDIUMPRIME; prime++) 
     {
        if (soln2[prime] == 0xFFFFFFFF) continue;
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        if (first)
        {
           correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
           soln1[prime]+=correction;
           while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
           soln2[prime]+=correction;
           while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;

           position = sieve+soln1[prime];
           position2 = sieve+soln2[prime];
        } else
        {
           position = offsets[prime];
           position2 = offsets2[prime];
        }
        diff=position2-position;
        
        ptimes4 = currentprime*4;
        register unsigned char * bound=end-ptimes4;
        while (bound - position > 0)  
        {  
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
        }
        while ((end - position > 0)&&(end - position - diff > 0))
        { 
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
              
        }
        position2 = position+diff;
        if (end - position2 > 0)
        { 
              (* position2)+=currentprimesize, position2+=currentprime;
        }
        if (end - position > 0)
        { 
              (* position)+=currentprimesize, position+=currentprime;
        }
        
        if (!last)
        {
           offsets[prime] = position;
           offsets2[prime] = position2;
        }
     } 
     
      for (u_int32_t prime=MEDIUMPRIME; prime<midprime; prime++) 
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        if (first)
        {
           correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
           soln1[prime]+=correction;
           while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
           soln2[prime]+=correction;
           while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
           
           position = sieve+soln1[prime];
           position2 = sieve+soln2[prime];
        } else
        {
           position = offsets[prime];
           position2 = offsets2[prime];
        }
        diff=position2-position;
           
        ptimes4 = 2*currentprime;
        register unsigned char * bound=end-ptimes4;
        while (bound - position > 0)  
        {  
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
        }
        position2 = position+diff;
        while ((end - position > 0)&&(end - position2 > 0))
        { 
              (* position)+=currentprimesize, position+=currentprime, (* position2)+=currentprimesize, position2+=currentprime;
              
        }
        
        if (end - position2 > 0)
        { 
              (* position2)+=currentprimesize, position2+=currentprime;
        }
        if (end - position > 0)
        { 
              (* position)+=currentprimesize, position+=currentprime;
        }
        if (!last)
        {
           offsets[prime] = position;
           offsets2[prime] = position2;
        }
     } 
     
     return;
}

/*===========================================================================
   Sieve 2:

   Function: Second sieve for larger primes

=========================================================================== */
void sieve2(u_int32_t M, u_int32_t numPrimes, unsigned char * sieve, int32_t last, int32_t first, int32_t polyadd, u_int32_t * soln1, u_int32_t * soln2, u_int32_t * polycorr, unsigned char * * offsets, unsigned char * * offsets2)
{
     register unsigned char currentprimesize; 
     register u_int32_t currentprime;
     register unsigned char * position2;
     register unsigned char * position;
     unsigned char * end;
     int32_t correction;
     
     memset(sieve,0,M*sizeof(unsigned char));
     memset(flags,0,numPrimes*sizeof(unsigned char));
     end = sieve+M;
     *end = 255; //sentinel to speed up sieve evaluators inner loop

     for (u_int32_t prime=midprime; prime<secondprime; prime++) 
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
        soln1[prime]+=correction;
        while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
        soln2[prime]+=correction;
        while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
           
        position = sieve+soln1[prime];
        position2 = sieve+soln2[prime];
           
        while ((end - position > 0)&&(end - position2 > 0))
        { 
             (* position)+=currentprimesize, position+=currentprime, (* position2)+=currentprimesize, position2+=currentprime;
        }
        
        if (end - position2 > 0)
        { 
              (* position2)+=currentprimesize;
        }
        if (end - position > 0)
        { 
              (* position)+=currentprimesize;
        }        
     }
     
     for (u_int32_t prime=secondprime; prime<numPrimes; prime++) 
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];
        
        correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
        soln1[prime]+=correction;
        while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
        soln2[prime]+=correction;
        while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
           
        position = sieve+soln1[prime];
        position2 = sieve+soln2[prime];
           
        while (end - position > 0)
        { 
              flags[prime]|=((unsigned char)1<<((position-sieve)&7)), (* position)+=currentprimesize, position+=currentprime;
        }
        
        while (end - position2 > 0)
        { 
              flags[prime]|=((unsigned char)1<<((position2-sieve)&7)), (* position2)+=currentprimesize, position2+=currentprime;
        }
     }
     
     return;
}

/*============================================================================

   random: 
           
   Function: Generates a pseudo-random integer between 0 and n-1 inclusive
   
============================================================================*/
u_int32_t random(u_int32_t upto)
{
   randval = ((u_int64_t)randval*1025416097+286824428)%(u_int64_t)4294967291;
   //randval = randval<<16+randval>>16;
   return randval%upto;
}


/*============================================================================
   mainRoutine:

   Function: Generates the polynomials, initialises and calls the sieve, 
             implementing cache blocking (breaking the sieve interval into
             small blocks for the small primes.

============================================================================*/
void mainRoutine(u_int32_t Mdiv2, mpz_t n, u_int32_t multiplier)
{
    mpz_t A; mpz_init(A);
    mpz_t B; mpz_init(B);
    mpz_t C; mpz_init(C);
    mpz_t D; mpz_init(D);
    mpz_t temp; mpz_init(temp);
    mpz_t temp2; mpz_init(temp2);
    mpz_t q; mpz_init(q);
    mpz_t r; mpz_init(r);
    mpz_t Bdivp2; mpz_init(Bdivp2);
     
    //gmp_randstate_t state;
    //gmp_randinit_default(state);
     
    u_int32_t u1;
     
    int32_t s, fact, span, min;
    u_int32_t p;     
    u_int32_t reps;
     
    u_int32_t curves = 0; 
     
    u_int32_t ** relations;
    int32_t * primecount;

    int32_t * exponents = (int32_t *) calloc(firstprime,sizeof(int32_t));
    if (exponents==NULL) 
    {
       printf("Unable to allocate memory!\n");
       exit(1);
    }
 
#ifdef TIMES
    counterfirst[2] = getcounter();
#endif
    s = mpz_sizeinbase(n,2)/28+1;
     
    u_int32_t * aind = (u_int32_t*) calloc(sizeof(u_int32_t),s);  
    u_int32_t * amodp = (u_int32_t*) calloc(sizeof(u_int32_t),s); 
    u_int32_t * Ainv = (u_int32_t*) calloc(sizeof(u_int32_t),numPrimes); 
    u_int32_t * soln1 = (u_int32_t*) calloc(sizeof(u_int32_t),numPrimes); 
    u_int32_t * soln2 = (u_int32_t*) calloc(sizeof(u_int32_t),numPrimes); 
    u_int32_t ** Ainv2B = (u_int32_t**) calloc(sizeof(u_int32_t*),s);
    if (Ainv2B==NULL) 
    {
       printf("Unable to allocate memory!\n");
       exit(1);
    }
    for (int32_t i=0; i<s; i++)
    {
       Ainv2B[i] = (u_int32_t *) calloc(sizeof(u_int32_t),numPrimes);
       if (Ainv2B[i]==NULL) 
       {
          printf("Unable to allocate memory!\n");
          exit(1);
       }
    } 
     
    mpz_t * Bterms = (mpz_t *)calloc(sizeof(mpz_t),s);
    mpz_array_init(*Bterms,s,mpz_sizeinbase(n,2)); 

    row m[1];
    m[0] = (row)calloc(numPrimes/32+1,sizeof(u_int32_t));
    //m = constructMat(numPrimes, relSought);
    
    la_col_t* colarray = (la_col_t*)calloc(relSought,sizeof(la_col_t)); //initialise a zeroed array of column structures
    
    relsFound = 0;
     
    mpz_t XArr[relSought]; 
     
    sieve = (unsigned char *) calloc(sizeof(unsigned char),Mdiv2*2+4); //one dword extra for sentinel to speed up sieve evaluation loop
    if (sieve==NULL) 
    {
       printf("Unable to allocate memory for sieve!\n");
       exit(1);
    }                
 
     
    flags = (unsigned char*) calloc(sizeof(unsigned char),numPrimes);
     
    offsets = (unsigned char * *)calloc(numPrimes,sizeof(unsigned char *));
    offsets2 = (unsigned char * *)calloc(numPrimes,sizeof(unsigned char *));
     
    relations = (u_int32_t * *) calloc(numPrimes,sizeof(u_int32_t *));
    for (int32_t i = 0; i < numPrimes; i++)
    {
       relations[i] = (u_int32_t *) calloc(50, sizeof(u_int32_t));
    }
     
    primecount = (int32_t *) calloc(numPrimes, sizeof(int32_t));
 
//Compute min A_prime and A_span
     
    mpz_mul_ui(temp,n,2);
    mpz_sqrt(temp,temp);
    mpz_div_ui(temp2,temp,Mdiv2);
    mpz_root(temp,temp2,s);
    for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++); 
    span = numPrimes/s/s/2;
    min=fact-span/2;
    
#ifdef ADETAILS
    printf("s = %d, fact = %d, min = %d, span = %d\n",s,fact,min,span);
#endif
     
//Compute first polynomial and adjustments
     
    while (relsFound < relSought)
    {
        int32_t i;
        mpz_set_ui(A,1);
        for (i=0; i<s-1; )
        {
           int32_t j,ran;
           //mpz_set_ui(temp,span);
           //mpz_urandomm(temp,state,temp);
           //ran = mpz_get_ui(temp);
           ran = random(span);
           for (j=0;((j<i)&&(aind[j]!=ran));j++);
           if (j==i) 
           {
              aind[i] = ran;
              mpz_mul_ui(A,A,factorBase[ran+min]);
              i++;
           }
        }  
        mpz_div(temp,temp2,A);
        for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++); 
        fact-=min;
        int32_t j;
        do
        {
           for (j=0;((j<i)&&(aind[j]!=fact));j++);
           fact++;
        } while (j!=i);
        fact--;
        aind[i] = fact;
        mpz_mul_ui(A,A,factorBase[fact+min]); 
       
        for (int32_t i=0; i<s; i++)
        {
           p = factorBase[aind[i]+min];
           mpz_div_ui(temp,A,p);
           //mpz_fdiv_r_ui(temp,temp,p*p);
	   amodp[i] = mpz_fdiv_r_ui(temp,temp,p);
                          
           //mpz_set_ui(temp3,p);
	   //mpz_invert(temp,temp,temp3);
	   mpz_set_ui(temp,modinverse(mpz_get_ui(temp),p));
           mpz_mul(temp,temp,sqrts[aind[i]+min]);
           mpz_fdiv_r_ui(temp,temp,p);
           if (mpz_cmp_ui(temp,p/2)>0)
           {
              mpz_sub_ui(temp,temp,p);
              mpz_neg(temp,temp);
           }
           mpz_mul(temp,temp,A);
           mpz_div_ui(Bterms[i],temp,p);
        }  
         
        mpz_set(B,Bterms[0]);
        for (int32_t i=1; i<s; i++)
        {
           mpz_add(B,B,Bterms[i]);
        }
         
        for (int32_t i=0; i<numPrimes; i++)
        {
           p = factorBase[i];
           //mpz_set_ui(temp3,p);
	   //mpz_fdiv_r_ui(temp,A,(u_int64_t)p*(u_int64_t)p);
	   //mpz_fdiv_r_ui(temp,temp,p);
	   //mpz_fdiv_r(temp,A,temp3);
	   //mpz_invert(temp3,temp,temp3);
	   //Ainv[i] = mpz_get_ui(temp3);
	   Ainv[i] = modinverse(mpz_fdiv_r_ui(temp,A,p),p);
             
           for (int32_t j=0; j<s; j++)
           {
              mpz_fdiv_r_ui(temp,Bterms[j],p);
	      //mpz_fdiv_r_ui(temp,temp,p);
              mpz_mul_ui(temp,temp,2*Ainv[i]);
              Ainv2B[j][i] = mpz_fdiv_r_ui(temp,temp,p);
           }
             
           mpz_fdiv_r_ui(temp,B,p);
	   //mpz_fdiv_r_ui(temp,temp,p);
           mpz_sub(temp,sqrts[i],temp);
           mpz_add_ui(temp,temp,p);
           mpz_mul_ui(temp,temp,Ainv[i]);
           mpz_add_ui(temp,temp,Mdiv2);
           soln1[i] = mpz_fdiv_r_ui(temp,temp,p);
           //mpz_set(temp,sqrts[i]);
           mpz_sub_ui(temp,sqrts[i],p);
           mpz_neg(temp,temp);
           mpz_mul_ui(temp,temp,2*Ainv[i]);
           soln2[i] = mpz_fdiv_r_ui(temp,temp,p)+soln1[i];
        }  
         
        for (int32_t polyindex=1; polyindex<(1<<(s-1))-1; polyindex++)
        {
           int32_t j,polyadd;
           u_int32_t * polycorr;
           for (j=0; j<s; j++)
           {
              if (((polyindex>>j)&1)!=0) break;
           }
           if ((polyadd = (((polyindex>>j)&2)!=0)))
           {
              mpz_add(B,B,Bterms[j]);
              mpz_add(B,B,Bterms[j]);
           } else
           {
              mpz_sub(B,B,Bterms[j]); 
              mpz_sub(B,B,Bterms[j]); 
           }
           polycorr = Ainv2B[j];
             
           int32_t index;
           for (int32_t j=0; j<s; j++)
           {
              index = aind[j]+min;
              p = factorBase[index];
              //mpz_set_ui(temp,p*p);
              //mpz_mul(temp,temp,temp);
              mpz_fdiv_r_ui(D,n,p*p);
              mpz_fdiv_r_ui(Bdivp2,B,p*p);
              mpz_mul_ui(temp,Bdivp2,amodp[j]);
              mpz_realloc2(temp3,64);
	      //mpz_set_ui(temp3,p);
	      //mpz_fdiv_r_ui(temp,temp,p*p);
	      mpz_fdiv_r_ui(temp,temp,p);
	      //mpz_invert(temp3,temp,temp3);
	      //u1=mpz_get_ui(temp3);
	      u1 = modinverse(mpz_fdiv_r_ui(temp,temp,p),p);        
              mpz_mul(temp,Bdivp2,Bdivp2);
              mpz_sub(temp,temp,D);
              mpz_neg(temp,temp);
              mpz_div_ui(temp,temp,p);
              mpz_mul_ui(temp,temp,u1);
              mpz_add_ui(temp,temp,Mdiv2);
              mpz_add_ui(temp,temp,p);
              //mpz_fdiv_r_ui(temp,temp,p*p);
	      soln1[index]=mpz_fdiv_r_ui(temp,temp,p);
              //soln1[index] = mpz_get_ui(temp);
              soln2[index]=0xFFFFFFFFl;
           }
           
// Count the number of polynomial curves used so far and compute the C coefficient of our polynomial
     
           curves++;
             
           mpz_mul(C,B,B);
           mpz_sub(C,C,n);
           mpz_divexact(C,C,A);
           
// Do the sieving and relation collection
     
           mpz_set_ui(temp,Mdiv2*2);
           mpz_fdiv_qr_ui(q,r,temp,CACHEBLOCKSIZE);
#ifdef TIMES
           counterfirst[3] = getcounter();
#endif
           sieve2(mpz_get_ui(temp),numPrimes,sieve,1,1,polyadd,soln1,soln2,polycorr,offsets,offsets2);
#ifdef TIMES
           countertotal[3]+=(getcounter()-counterfirst[3]);
           counterfirst[0] = getcounter();
#endif
           sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve,0,1,polyadd,soln1,soln2,polycorr,offsets,offsets2);
           if (mpz_cmp_ui(q,1)>0)
           {
              for (reps = 1;reps < mpz_get_ui(q)-1; reps++)
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,0,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              }
              if (mpz_cmp_ui(r,0)==0)
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,1,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              } else
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,0,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
                 reps++;
                 sieveInterval(mpz_get_ui(r),numPrimes,sieve+CACHEBLOCKSIZE*reps,1,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              }
           }
           
#ifdef TIMES  
           countertotal[0]+=(getcounter()-counterfirst[0]);
           counterfirst[1] = getcounter();
#endif
           evaluateSieve(relations,0,mpz_get_ui(temp),sieve,A,B,C,soln1, soln2, polyadd, polycorr, m,XArr,aind,min,s,multiplier,exponents,colarray);
#ifdef TIMES
           countertotal[1]+=(getcounter()-counterfirst[1]);
#endif
    }  
     
#ifdef COUNT
        if (curves%20==0) printf("%ld curves.\n",(long)curves);
#endif
    }
     
#ifdef CURPARTS
    printf("%ld curves, %ld partials.\n",(long)curves,(long)partials);
#endif
             
#ifdef REPORT
    printf("Done with sieving!\n");
#endif
     
u_int32_t ncols = relSought;
u_int32_t nrows = numPrimes;

reduce_matrix(&nrows, &ncols, colarray);
     
#ifdef REPORT
     //u_int32_t kernel = gaussReduce(m, numPrimes, relSought, 0);
     u_int64_t* nullrows;
     do {
         nullrows = block_lanczos(nrows, 0, ncols, colarray);
     } while (nullrows == NULL);
     //printf("%d relations in kernel.\n",kernel); 
#else
     //gaussReduce(m, numPrimes, relSought, 0);
     u_int64_t* nullrows;
     do {
         nullrows = block_lanczos(nrows, 0, ncols, colarray);
     } while (nullrows == NULL);
     
#endif
     

#ifdef TIMES
    countertotal[2]+=(getcounter()-counterfirst[2]);
    printf("Total time = %.0f\n",countertotal[2]);
    printf("Polynomial generation time = %.0f\n",countertotal[2]-countertotal[0]-countertotal[1]-countertotal[3]);
    printf("Small prime sieve time = %.0f\n",countertotal[0]);
    printf("Large prime sieve time = %.0f\n",countertotal[3]);
    printf("Evaluate sieve time = %.0f\n",countertotal[1]); 
#endif
    
// We want factors of n, not kn, so divide out by the multiplier
     
    mpz_div_ui(n,n,multiplier);
    
// Now do the "square root" and GCD steps hopefully obtaining factors of n
    for (int32_t l = 0;l<64;l++)
    {
        mpz_set_ui(temp,1);
        mpz_set_ui(temp2,1);
        memset(primecount,0,numPrimes*sizeof(int32_t));
        for (int32_t i = 0; i<ncols; i++)
        {
           if (getNullEntry(nullrows,i,l)) 
           {
              mpz_mul(temp2,temp2,XArr[colarray[i].orig]);
              for (int32_t j=1; j<=relations[colarray[i].orig][0]; j++)
              {
                 primecount[relations[colarray[i].orig][j]]++;
              }
           }
           if (i%30==0) mpz_mod(temp2,temp2,n);
        }
        for (int32_t j = 0; j<numPrimes; j++) 
        {
           mpz_set_ui(temp3,factorBase[j]);
           mpz_pow_ui(temp3,temp3,primecount[j]/2);
           mpz_mul(temp,temp,temp3);
           if (j%30==0) mpz_mod(temp,temp,n);
        }
        mpz_sub(temp,temp2,temp);
        mpz_gcd(temp,temp,n);
        if ((mpz_cmp(temp,n)!=0)&&(mpz_cmp_ui(temp,1)!=0)) //print only nontrivial factors
        {
           gmp_printf("%Zd\n",temp);
        }
    }
     
    return;
}

/*===========================================================================
   Main Program:

   Function: Factors a user specified number using a quadratic sieve

===========================================================================*/
int main(int argc, unsigned char *argv[])
{
    u_int32_t multiplier;

    initSieve(); 
    
    printf("Input number to factor [ >=40 decimal digits]: "); 
    gmp_scanf("%Zd",n);getchar();
    
    decdigits = mpz_sizeinbase(n,10);
    if (decdigits < 40) 
    {
       printf("Error in input or number has too few digits.\n");
       exit(1);
    }
    
    multiplier = knuthSchroeppel(n);
    mpz_mul_ui(n,n,multiplier);
    
  if (decdigits<=91) 
  {
    numPrimes=primesNo[decdigits-MINDIG];
    
    Mdiv2 = sieveSize[decdigits-MINDIG]/SIEVEDIV;
    if (Mdiv2*2 < CACHEBLOCKSIZE) Mdiv2 = CACHEBLOCKSIZE/2;

#ifdef REPORT
    printf("Using multiplier: %ld\n",(long)multiplier);
    printf("%ld primes in factor base.\n",(long)numPrimes);
    printf("Sieving interval M = %ld\n",(long)Mdiv2*2);
#endif
    
    if (numPrimes < SECONDPRIME) secondprime = numPrimes;
    else secondprime = SECONDPRIME;
    if (numPrimes < MIDPRIME) midprime = numPrimes;
    else midprime = MIDPRIME;
    
    firstprime = firstPrimes[decdigits-MINDIG];
    errorbits = errorAmounts[decdigits-MINDIG];
    threshold = thresholds[decdigits-MINDIG];
    less = lessrels[decdigits-MINDIG];
    
  } else //all bets are off
  {
     numPrimes = 64000;
     Mdiv2 = 192000/SIEVEDIV;

#ifdef REPORT
    printf("Using multiplier: %ld\n",(long)multiplier);
    printf("%ld primes in factor base.\n",(long)numPrimes);
    printf("Sieving interval M = %ld\n",(long)Mdiv2*2);
#endif
    
    secondprime = SECONDPRIME;
    midprime = MIDPRIME;
    firstprime = 30;
    errorbits = decdigits/4 + 2;
    threshold = 43+(7*decdigits)/10;
    less = 2000;
    
  }
    relSought = numPrimes-less; 
    computeFactorBase(n, numPrimes, multiplier);

    computeSizes(numPrimes);
    
    TonelliInit();
    tonelliShanks(numPrimes,n);
    
    mainRoutine(Mdiv2, n,multiplier);
    
    getchar();     
    return 0;
}

