#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include "Z.h"
#include "flint.h"

#define DEBUG2 1
#define DEBUG 0

gmp_randstate_t state_Z;
int rand_init = 0;

/*
   Memory manager to allocate a single Z_t. It returns a pointer to the Z_t. 
   Z_t's should be released in the order they were allocated.
*/

#define RESALLOC 100 //allocate this many Z_t's at once to save on overheads

Z_t** reservoir; // Array of pointers to Z_t's in the reservoir
unsigned long rescount=0; //Next available Z_t in reservoir
unsigned long currentalloc=0; //total number of Z_t's in reservoir

Z_t* Z_alloc(void)
{
   static int initialised = 0;
   static Z_t** tempres;
   Z_t* alloc_d;
   
   //allocate another block of Z_t's if none are currently allocated, or the reservoir is depleted
   if (rescount==currentalloc) // need more limb_memp_t's
   {
      if (!initialised) 
      {
         reservoir = (Z_t**)malloc(RESALLOC*sizeof(Z_t*)); //allocate space for the array of pointers
         reservoir[0] = (Z_t*)malloc(RESALLOC*sizeof(Z_t)); //allocate space for the Z_t's
         for (unsigned long i=0; i<RESALLOC-1; i++)
         {
             reservoir[i+1]=reservoir[i]+1; //initialise the array
             mpz_init(*reservoir[i]); //initialise the Z_t's
         }
         mpz_init(*reservoir[RESALLOC-1]);
         rescount=0;
         initialised = 1;
         currentalloc = RESALLOC;
      } else
      {
         //copy old reservoir into larger one
         tempres = reservoir;
         reservoir = (Z_t**)malloc((currentalloc+RESALLOC)*sizeof(Z_t*));
         reservoir[currentalloc] = (Z_t*)malloc(RESALLOC*sizeof(Z_t));  
         memcpy(reservoir,tempres,currentalloc*sizeof(Z_t*)); 
         for (unsigned long i=currentalloc; i<RESALLOC+currentalloc-1; i++)
         {
             reservoir[i+1]=reservoir[i]+1; //initialise the array
             mpz_init(*reservoir[i]); //initialise the Z_t's
         }
         mpz_init(*reservoir[currentalloc+RESALLOC-1]);
         
         currentalloc+=RESALLOC;
         //free old reservoir
         free(tempres);  
      }       
   }
   
   alloc_d = reservoir[rescount];
   rescount++;
   return alloc_d;
}

/*
    Release a Z_t back into the reservoir
*/

void Z_release(void)
{   
    rescount--;
}

/*
    returns a^exp
*/

unsigned long Z_pow_long(unsigned long a, unsigned long exp)
{
   if (exp == 0) return 1;
   if (a == 1) return 1;
   
   unsigned long power = a;

   for (unsigned long i = 1; i < exp; i++)
      power *= a;

   return power;
}

/* 
    sets res to a*b modulo p
    assumes res is not p
*/

void Z_mulmod(Z_t res, Z_t a, Z_t b, Z_t p)
{
     Z_t* temp = Z_alloc();
     
     mpz_fdiv_r(*temp,a,p);
     mpz_fdiv_r(res,b,p);
     mpz_mul(res,*temp,res);
     mpz_fdiv_r(res,res,p);
     
     Z_release();
     
     return;
}


/* 
     sets res to a*b modulo p
     DIRTY for 64 bit due to long long
*/

unsigned long Z_mulmod_ui(Z_t res, Z_t a, Z_t b, unsigned long p)
{
     unsigned long result;
     
     result = ((unsigned long long)mpz_fdiv_r_ui(res,a,p)*(unsigned long long)mpz_fdiv_r_ui(res,b,p))%p;
     mpz_set_ui(res,result);
     
     return result;
}

/*
         returns a^exp modulo p
         assumes a is not (much) bigger than p
         DIRTY for 64 bit due to long long
*/

long Z_powm_long(long a, long exp, long p)
{
   long x, y;

   unsigned long e;

   if (exp < 0)
      e = (unsigned long) -exp;
   else
      e = exp;

   x = 1;
   y = a;
   while (e) {
      if (e & 1) x = ((long long)x * (long long)y)% p;
      y = ((long long)y * (long long)y)% p;
      e = e >> 1;
   }

   if (exp < 0) x = Z_invert_long(x, p);

   return x;
}

/* 
    returns the next prime after n (does not check if the result is too big)
*/

unsigned long Z_nextprime_long(unsigned long n)
{
   Z_t* temp = Z_alloc();
     
   mpz_set_ui(*temp,n);
   mpz_nextprime(*temp,*temp);
   
   Z_release();
   
   return mpz_get_ui(*temp);
}

/*
    returns a random prime up to the specified number of bits
*/

void Z_randomprime(Z_t res, unsigned long bitsize)
{
   static int rand_running=0;
   
   if (!rand_init) //required for thread safety
   {
      gmp_randinit_default(state_Z);
      rand_init = 1;
   }
   
   while (rand_running){} //required for thread safety to prevent simultaneous access to state_Z
   
   rand_running=1; //required for thread safety
   
   do { 
      mpz_urandomb(res,state_Z,bitsize);
      mpz_setbit(res,0);
      while (!mpz_probab_prime_p(res,5)) mpz_add_ui(res,res,2);
   } while (mpz_sizeinbase(res,2)>bitsize); 
   
   rand_running=0; //required for thread safety
}

/* 
    returns the inverse of a modulo p
*/

unsigned long Z_invert_long(unsigned long a, unsigned long p)
{
   long u1=1, u3=a;
   long v1=0, v3=p;
   long t1=0, t3=0;
   long quot;
   
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

/* 
     returns gcd(x, y) = a*x + b*y.
*/

long Z_extgcd_long(long* a, long* b, long x, long y)
{
   long u1=1, v1=0;
   long u2=0, v2=1;
   long t1, t2;
   long u3, v3;
   long quot, rem;
   long xsign = 0;
   long ysign = 0;
   
   if (x < 0) {
      x = -x;
      xsign = 1;
   }

   if (y < 0) {
      y = -y;
      ysign = 1;
   }
   
   u3 = x, v3 = y;

   while (v3) {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               rem = u3;
               t1 = u2; u2 = u1; u1 = t1; u3 = v3;
               t2 = v2; v2 = v1; v1 = t2; v3 = rem;
            } else 
            {
               t1 = u2; u2 = u1 - u2; u1 = t1; u3 = v3;
               t2 = v2; v2 = v1 - v2; v1 = t2; v3 = quot;
            }
         } else if (quot < (v3<<1))
         {  
            t1 = u2; u2 = u1 - (u2<<1); u1 = t1; u3 = v3;
            t2 = v2; v2 = v1 - (v2<<1); v1 = t2; v3 = quot-u3;
         } else
         {
            t1 = u2; u2 = u1 - 3*u2; u1 = t1; u3 = v3;
            t2 = v2; v2 = v1 - 3*v2; v1 = t2; v3 = quot-(u3<<1);
         }
      } else
      {
         quot=u3/v3;
         rem = u3 % v3;
         t1 = u2; u2 = u1 - quot*u2; u1 = t1; u3 = v3;
         t2 = v2; v2 = v1 - quot*v2; v1 = t2; v3 = rem;
      }
   }

   if (xsign)
      u1 = -u1;

   if (ysign)
      v1 = -v1;

   *a = u1;
   *b = v1;
   
   return u3;
}

/* 
     returns gcd(x, y)
*/

unsigned long Z_gcd_long(long x, long y)
{
   if (x < 0) {
      x = -x;
   }

   if (y < 0) {
      y = -y;
   }

   long u3 = x, v3 = y;
   long quot, rem;

   while (v3) {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            { 
               rem = u3;
               u3 = v3;
               v3 = rem;
            } else 
            {
               u3 = v3;
               v3 = quot;
            }
         } else if (quot < (v3<<1))
         {  
            u3 = v3;
            v3 = quot-u3;
         } else
         {
            u3 = v3;
            v3 = quot-(u3<<1);
         }
      } else
      {
         rem = u3 % v3;
         u3 = v3;
         v3 = rem;
      }
   }

   return u3;
}

/* 
     sets res to the square root of a modulo p for a prime p
     returns 0 if a is not a square modulo p
*/

int Z_sqrtmod(Z_t res, Z_t a, Z_t p) 
{
     int r,k,m;
     Z_t* p1 = Z_alloc();
     Z_t* two = Z_alloc();
     Z_t* b = Z_alloc();
     Z_t* g = Z_alloc();
     Z_t* bpow = Z_alloc();
     Z_t* gpow = Z_alloc();
     Z_t* mk = Z_alloc();
     
     if (mpz_kronecker(a,p)!=1) 
     {
         mpz_set_ui(res,0);
         return 0;   //return 0 if a is not a square mod p
     }
     
     if ((mpz_cmp_ui(a,0)==0)||(mpz_cmp_ui(a,1)==0)) 
     {
        mpz_set(res,a);
        return 1;
     }
     
     if ((Z_tstbit(p,0)==1)||(Z_tstbit(p,1)==1))
     {
        mpz_add_ui(*p1,p,1);
        mpz_fdiv_q_2exp(*p1,*p1,2);
        mpz_powm(res,a,*p1,p);
        return 1;
     }
     
     mpz_set_ui(*two,2);
     
     mpz_sub_ui(*p1,p,1);
     r = mpz_remove(*p1,*p1,*two);
     mpz_powm(*b,a,*p1,p);
     for (k=2; ;k++)
     {
         if (mpz_ui_kronecker(k,p) == -1) break;
     }
     mpz_set_ui(*mk,k);
     mpz_powm(*g,*mk,*p1,p);
     mpz_add_ui(*p1,*p1,1);
     mpz_divexact_ui(*p1,*p1,2);
     mpz_powm(res,a,*p1,p);
     if (!mpz_cmp_ui(*b,1)) return 1;
     
     while (mpz_cmp_ui(*b,1))
     {
           mpz_set(*bpow,*b);
           for (m=1; (m<=r-1) && (mpz_cmp_ui(*bpow,1));m++)
           {
               mpz_powm_ui(*bpow,*bpow,2,p);
           }
           mpz_set(*gpow,*g);
           for (int i = 1;i<= r-m-1;i++)
           {
               mpz_powm_ui(*gpow,*gpow,2,p);
           };
           Z_mulmod(res,res,*gpow,p);
           mpz_powm_ui(*gpow,*gpow,2,p);
           Z_mulmod(*b,*b,*gpow,p);
           mpz_set(*gpow,*g);
           r = m;
     }
          
     Z_release();Z_release();Z_release();Z_release();
     Z_release();Z_release();Z_release();
     
     return 1;
}

/* 
     computes the square root of a modulo p^k when given z, the square root mod p^(k-1)
*/

inline void sqrtmodpow(Z_t res, Z_t z, Z_t a, Z_t pk, Z_t tempsqpow, Z_t inv)
{
     mpz_mul_ui(inv,z,2);
     mpz_invert(inv,inv,pk);
     mpz_set(tempsqpow,a);
     mpz_submul(tempsqpow,z,z);
     mpz_fdiv_r(tempsqpow,tempsqpow,pk);
     Z_mulmod(tempsqpow,tempsqpow,inv,pk);
     mpz_add(tempsqpow,tempsqpow,z);
     mpz_set(res,tempsqpow);
     
     return;
} 

/* 
     computes the square root of a modulo p^k when given z, the square root mod p^{k-1}
*/

void Z_sqrtmodpklift(Z_t res, Z_t z, Z_t a, Z_t pk)
{
     Z_t* tempsqpow = Z_alloc();
     Z_t* inv = Z_alloc();
     
     sqrtmodpow( res, z, a, pk, *tempsqpow, *inv);
     
     Z_release();Z_release();
}

/* 
     computes the square root of a modulo p^k when given the square root modulo p
*/

void Z_sqrtmodptopk(Z_t res, Z_t sqrt, Z_t a, Z_t p, int k)
{
     Z_t* tempsqpow = Z_alloc();
     Z_t* inv = Z_alloc();
     Z_t* pk = Z_alloc();
     
     mpz_set(res,sqrt);
     mpz_set(*pk,p);
     for (int i = 2; i<=k; i++)
     {
            mpz_mul(*pk,*pk,p);
            sqrtmodpow(res,res,a,*pk, *tempsqpow, *inv);
     }
     
     Z_release();Z_release();Z_release();
}

/* 
     computes the square root of a modulo p^k 
*/

int Z_sqrtmodpk(Z_t res, Z_t a, Z_t p, int k)
{
     int exists;
     Z_t* sqrtmodp = Z_alloc();
     
     exists = Z_sqrtmod(*sqrtmodp,a,p);
     Z_sqrtmodptopk(res,*sqrtmodp,a,p,k);
     
     Z_release();
     
     return exists;
}


/* 
     Find res mod n=n1*n2 such that res = x1 mod n1 and res = x2 mod n2
*/
 
void Z_CRT(Z_t res, Z_t n, Z_t x1, Z_t x2, Z_t n1, Z_t n2)
{
     Z_t* ntemp = Z_alloc();
     Z_t* restemp = Z_alloc();
     Z_t* chtemp = Z_alloc();
     
     mpz_mul(*ntemp,n1,n2);
     mpz_invert(*restemp,n2,n1);
     Z_mulmod(*restemp,res,n2,*ntemp);
     Z_mulmod(*restemp,*restemp,x1,*ntemp);
     
     mpz_invert(*chtemp,n1,n2);
     Z_mulmod(*chtemp,*chtemp,n1,*ntemp);
     Z_mulmod(*chtemp,*chtemp,x2,*ntemp);
     
     mpz_add(res,*restemp,*chtemp);
     mpz_fdiv_r(res,res,*ntemp);
     mpz_set(n,*ntemp);
     
     Z_release();Z_release();Z_release();
     
     return;
}

void Z_urandomb(Z_t res, unsigned long n)
{
   static int rand_running=0;
   
   if (!rand_init) //required for thread safety
   {
      gmp_randinit_default(state_Z);
      rand_init = 1;
   }
   
   while (rand_running){} //required for thread safety to prevent simultaneous access to state_Z
   
   rand_running=1; //required for thread safety
   
   mpz_urandomb(res,state_Z,n);
   
   rand_running=0; //required for thread safety
}

void Z_urandomm(Z_t res, Z_t n)
{
   static int rand_running=0;
   
   if (!rand_init) //required for thread safety
   {
      gmp_randinit_default(state_Z);
      rand_init = 1;
   }
   
   while (rand_running){} //required for thread safety to prevent simultaneous access to state_Z
   
   rand_running=1; //required for thread safety
   
   mpz_urandomm(res,state_Z,n);
   
   rand_running=0; //required for thread safety
}

void Z_rrandomb(Z_t res, unsigned long n)
{
   static int rand_running=0;
   
   if (!rand_init) //required for thread safety
   {
      gmp_randinit_default(state_Z);
      rand_init = 1;
   }
   
   while (rand_running){} //required for thread safety to prevent simultaneous access to state_Z
   
   rand_running=1; //required for thread safety
   
   mpz_rrandomb(res,state_Z,n);
   
   rand_running=0; //required for thread safety
}

/*
    Compute the Montgomery reduced form of a mod m
    Returns n such that m < 2^n (n will be divisible by FLINT_BITS)
    Assumes a is already reduced mod m
*/

unsigned long Z_mont_red(mpz_t res, mpz_t a, mpz_t m)
{
   unsigned long n = mpz_size(m)*FLINT_BITS;
   
   mpz_mul_2exp(res, a, n);
   mpz_mod(res, res, m);
   
   return n;
}

/* 
    Compute the Montgomery multiplication r = a*b mod m assuming a and b are in 
    Montgomery form with respect to 2^n where m < 2^n and R is -m mod 2^n
*/

void Z_mont_mul(mpz_t res, mpz_t a, mpz_t b, mpz_t m, mpz_t R, unsigned long n)
{
   mpz_t x, s;
   mpz_init(x);
   mpz_init(s);
   
   mpz_mul(x, a, b);
   mpz_fdiv_r_2exp(s, x, n);
   mpz_mul(s, s, R);
   mpz_fdiv_r_2exp(s, s, n);
   mpz_mul(res, s, m);
   mpz_add(res, res, x);
   mpz_fdiv_q_2exp(res, res, n);
   
   if (mpz_cmp(res, m) >= 0) mpz_sub(res, res, m);
    
   mpz_clear(x);
   mpz_clear(s);
}

/* 
    Compute a^exp mod m using Montgomery reduction
    Requires that m is odd and positive and that exp is positive
*/

void Z_expmod_mont(mpz_t res, mpz_t a, mpz_t exp, mpz_t m)
{
   unsigned long n;
   unsigned long bits = mpz_sizeinbase(exp, 2);
   mpz_t aRED;
   mpz_t powRED;
   mpz_t R;
   mpz_t temp;
   int flag = 0;
   
   mpz_init(aRED);
   mpz_init(powRED);
   mpz_init(R);
   mpz_init(temp);
   
   n = Z_mont_red(aRED, a, m);
   
   mpz_set_ui(temp, 1);
   mpz_mul_2exp(temp, temp, n);
   mpz_invert(R, m, temp);
   mpz_sub(R, temp, R);
   if (mpz_cmp(R, temp) == 0) mpz_sub(R, R, temp);
   
   mpz_set(powRED, aRED);
#ifdef DEBUG
   gmp_printf("powRED = %Zd\n", powRED);
#endif
   
   for (unsigned long i = 0; i < bits - 1; i++)
   {
      if (mpz_tstbit(exp, i))
      {
         if (flag) Z_mont_mul(res, res, powRED, m, R, n);
         else 
         {
            mpz_set(res, powRED);
            flag = 1;
         }
      }
      Z_mont_mul(powRED, powRED, powRED, m, R, n);
#ifdef DEBUG
      gmp_printf("powRED = %Zd\n", powRED);
#endif
   }
   
   if (flag) Z_mont_mul(res, res, powRED, m, R, n);
   else mpz_set(res, powRED);
   
   mpz_set_ui(temp, 1);
   Z_mont_mul(res, res, temp, m, R, n);
   
   mpz_clear(temp);
   mpz_clear(R);
   mpz_clear(powRED);
   mpz_clear(aRED);
}

void Z_divrem_jebelean(mpz_t Q, mpz_t R, mpz_t A, mpz_t B)
{
   unsigned long n = mpz_size(B);
   unsigned long m = mpz_size(A) - n;
   
   if ((long) m < 0)
   {
      mpz_set_ui(Q, 0);
      mpz_set(R, A);
      return;
   }
   
   if (m < 64) 
   {
      mpz_fdiv_qr(Q, R, A, B);
      return;
   }
   
   unsigned long k = m/2;
   
   Z_t * B0 = Z_alloc();
   Z_t * B1 = Z_alloc();
   Z_t * A0 = Z_alloc();
   Z_t * A1 = Z_alloc();
   Z_t * Q0 = Z_alloc();
   Z_t * Q1 = Z_alloc();
   Z_t * R0 = Z_alloc();
   Z_t * R1 = Z_alloc();
   Z_t * temp = Z_alloc();
   Z_t * temp2 = Z_alloc();
   Z_t * temp3 = Z_alloc();
   
   
   mpz_fdiv_q_2exp(*B1, B, FLINT_BITS*k);
   mpz_fdiv_q_2exp(*A1, A, FLINT_BITS*2*k);
   
   Z_divrem_jebelean(*Q1, *R1, *A1, *B1);
   mpz_fdiv_r_2exp(*B0, B, FLINT_BITS*k);
   mpz_fdiv_r_2exp(*A0, A, FLINT_BITS*2*k);
   mpz_mul_2exp(*temp, *R1, FLINT_BITS*2*k);
   mpz_add(*temp, *temp, *A0);
   mpz_mul_2exp(*temp2, *Q1, FLINT_BITS*k);
   mpz_mul(*temp2, *temp2, *B0);
   mpz_sub(*temp, *temp, *temp2);
   mpz_mul_2exp(*temp2, B, FLINT_BITS*k);
   
   while (mpz_cmp_ui(*temp, 0) < 0) 
   {
      mpz_sub_ui(*Q1, *Q1, 1);
      mpz_add(*temp, *temp, *temp2);
   }
   mpz_fdiv_q_2exp(*temp2, *temp, FLINT_BITS*k); 
   Z_divrem_jebelean(*Q0, *R0, *temp2, *B1);
   
   mpz_fdiv_r_2exp(*temp2, *temp, FLINT_BITS*k);
   mpz_mul_2exp(R, *R0, FLINT_BITS*k);
   mpz_add(R, R, *temp2);
   mpz_submul(R, *Q0, *B0);
   while (mpz_cmp_ui(R, 0) < 0) 
   {
      mpz_sub_ui(*Q0, *Q0, 1);
      mpz_add(R, R, B);
   }
   mpz_mul_2exp(Q, *Q1, FLINT_BITS*k);
   mpz_add(Q, Q, *Q0);
   
   Z_release(); Z_release(); Z_release(); Z_release();
   Z_release(); Z_release(); Z_release(); Z_release();
   Z_release(); Z_release(); Z_release();
}

void Z_rem_jebelean(mpz_t R, mpz_t A, mpz_t B)
{
   unsigned long n = mpz_size(B);
   unsigned long m = mpz_size(A) - n;
   
   if ((long) m < 0)
   {
      mpz_set(R, A);
      return;
   }
   
   if (m < 64) 
   {
      mpz_fdiv_r(R, A, B);
      return;
   }
   
   unsigned long k = m/2;
   
   Z_t * B0 = Z_alloc();
   Z_t * B1 = Z_alloc();
   Z_t * A0 = Z_alloc();
   Z_t * A1 = Z_alloc();
   Z_t * Q0 = Z_alloc();
   Z_t * Q1 = Z_alloc();
   Z_t * R0 = Z_alloc();
   Z_t * R1 = Z_alloc();
   Z_t * temp = Z_alloc();
   Z_t * temp2 = Z_alloc();
   Z_t * temp3 = Z_alloc();
   
   
   mpz_fdiv_q_2exp(*B1, B, FLINT_BITS*k);
   mpz_fdiv_q_2exp(*A1, A, FLINT_BITS*2*k);
   
   Z_divrem_jebelean(*Q1, *R1, *A1, *B1);
   mpz_fdiv_r_2exp(*B0, B, FLINT_BITS*k);
   mpz_fdiv_r_2exp(*A0, A, FLINT_BITS*2*k);
   mpz_mul_2exp(*temp, *R1, FLINT_BITS*2*k);
   mpz_add(*temp, *temp, *A0);
   mpz_mul_2exp(*temp2, *Q1, FLINT_BITS*k);
   mpz_mul(*temp2, *temp2, *B0);
   mpz_sub(*temp, *temp, *temp2);
   mpz_mul_2exp(*temp2, B, FLINT_BITS*k);
   
   while (mpz_cmp_ui(*temp, 0) < 0) 
   {
      mpz_sub_ui(*Q1, *Q1, 1);
      mpz_add(*temp, *temp, *temp2);
   }
   mpz_fdiv_q_2exp(*temp2, *temp, FLINT_BITS*k); 
   Z_divrem_jebelean(*Q0, *R0, *temp2, *B1);
   
   mpz_fdiv_r_2exp(*temp2, *temp, FLINT_BITS*k);
   mpz_mul_2exp(R, *R0, FLINT_BITS*k);
   mpz_add(R, R, *temp2);
   mpz_submul(R, *Q0, *B0);
   while (mpz_cmp_ui(R, 0) < 0) 
   {
      mpz_add(R, R, B);
   }
   
   Z_release(); Z_release(); Z_release(); Z_release();
   Z_release(); Z_release(); Z_release(); Z_release();
   Z_release(); Z_release(); Z_release();
}

void Z_mulmod_jebelean(mpz_t res, mpz_t a, mpz_t b, mpz_t m)
{
   Z_t * temp = Z_alloc();
   
   mpz_mul(*temp, a, b);
   Z_rem_jebelean(res, *temp, m);
   
   Z_release();
}


void Z_expmod_jebelean(mpz_t res, mpz_t a, mpz_t exp, mpz_t m)
{
   unsigned long n;
   unsigned long bits = mpz_sizeinbase(exp, 2);
   mpz_t aRED;
   mpz_t powRED;
   mpz_t temp;
   int flag = 0;
   
   mpz_init(aRED);
   mpz_init(powRED);
   mpz_init(temp);
   
   mpz_set(powRED, a);
#if DEBUG
   gmp_printf("powRED = %Zd\n", powRED);
#endif
   
   for (unsigned long i = 0; i < bits - 1; i++)
   {
      if (mpz_tstbit(exp, i))
      {
         if (flag) Z_mulmod_jebelean(res, res, powRED, m);
         else 
         {
            mpz_set(res, powRED);
            flag = 1;
         }
      }
      Z_mulmod_jebelean(powRED, powRED, powRED, m);
#if DEBUG
      gmp_printf("powRED = %Zd\n", powRED);
#endif
   }
   
   if (flag) Z_mulmod_jebelean(res, res, powRED, m);
   else mpz_set(res, powRED);
   
   mpz_clear(temp);
   mpz_clear(powRED);
   mpz_clear(aRED);
}
