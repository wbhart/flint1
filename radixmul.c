#include <gmp.h>
#include <stdlib.h>
#include "Z.h"
#include "flint.h"
#include "mpn_extras.h"

#define UWtype mp_limb_t
#define UHWtype mp_limb_t
#define UDWtype mp_limb_t 
#define W_TYPE_SIZE FLINT_BITS_PER_LIMB
#define SItype u_int32_t
#define USItype int32_t
#define DItype int64_t
#define UDItype u_int64_t

#include "longlong.h"

//todo: move *_mod* into Zmodpoly.c

/* returns hi*R+lo modulo p assuming p is of the form 2^63-a where R=2^64
   and a is no more than 9 bits                                            */
static inline mp_limb_t modp(mp_limb_t hi, mp_limb_t lo, u_int32_t a)
{
    mp_limb_t temp;
    mp_limb_t hi2,lo2;
    
    //first reduction
    temp=(u_int32_t)((hi&18446741874686296064U)>>41)*a; //result is at most 32 bits
    hi2=(temp>>22);
    lo2=(temp<<44);
    add_ssaaaa(hi2,lo2,hi2,lo2,hi&2199023255551,lo); //hi2 has at most 42 bits
    
    //second reduction
    temp=(u_int32_t)((hi2&4398045986816U)>>19)*a; //result is at most 32 bits
    temp<<=22;
    add_ssaaaa(temp,lo2,hi2&524287U,lo2,0,temp); //temp has at most 20 bits
    
    //third reduction
    temp=(u_int32_t)((temp<<3)+(lo2>>61))*a; //result is at most 32 bits
    temp<<=4;
    temp=((lo2&2305843009213693951U)+temp);  //temp is at most 62 bits

    if ((temp+a)&9223372036854775808U) return (temp+a)-9223372036854775808U;
    else return temp;  
}


static inline void add_mod(mp_limb_t* asum, mp_limb_t* a1, mp_limb_t* a2, unsigned long size, mp_limb_t p)
{
    for (unsigned long i=0; i<size; i++)
    {
       asum[i] = a1[i]+a2[i];
       if (asum[i]>=p) asum[i]-=p;
    }
}

static inline void sub_mod(mp_limb_t* adiff, mp_limb_t* a1, mp_limb_t* a2, unsigned long size, mp_limb_t p)
{
    for (unsigned long i=0; i<size; i++)
    {
       adiff[i] = a1[i]+p;
       adiff[i]-=a2[i];
       if (adiff[i]>=p) adiff[i]-=p;
    }
}

// assumes x and y are both in the range [0,p)
static inline mp_limb_t mul_mod(mp_limb_t x, mp_limb_t y, mp_limb_t p)
{
    mp_limb_t hi,lo,quo,rem;
    
    umul_ppmm(hi,lo,x,y);
    udiv_qrnnd(quo,rem,hi,lo,p);
    
    return rem;
}

static inline mp_limb_t addmul_mod(mp_limb_t x, mp_limb_t y, mp_limb_t p)
{
    mp_limb_t hi,lo,quo,rem;
    
    umul_ppmm(hi,lo,x,y);
    add_ssaaaa(hi,lo,hi,lo,0,x);
    udiv_qrnnd(quo,rem,hi,lo,p);
    
    return rem;
}

/* this can do at most a 4x4 multiplication due to p being 62 bits */
static inline void mul_mod_naieve(mp_limb_t* res, mp_limb_t* a,  mp_limb_t* b, mp_limb_t* scr, unsigned long size, mp_limb_t p)
{
     mp_limb_t hi1,lo1,hi2,lo2,quo,rem;
     unsigned long i,d;
     
     if (size == 1) 
     {
        res[0] = mul_mod(a[0],b[0],p);
        return;
     }
     
     if (size == 2) {
        res[0] = mul_mod(a[0], b[0],p);
        umul_ppmm(hi1,lo1,a[0], b[1]);
        res[2] = mul_mod(a[1], b[1],p);
        umul_ppmm(hi2,lo2,a[1], b[0]);
        add_ssaaaa(hi1,lo1,hi1,lo1,hi2,lo2);
        udiv_qrnnd(quo,res[1],hi1,lo1,p);
        return;
     } 
     
     if (size == 4)
     {
        d=9223372036854775808U-p;
        /* do 4x4 multiplication */ 
        umul_ppmm(scr[1],scr[0],a[0], b[0]); //0x0
        umul_ppmm(scr[3],scr[2],a[0], b[1]); //0x1
        umul_ppmm(scr[5],scr[4],a[0], b[2]); //0x2
        umul_ppmm(scr[7],scr[6],a[0], b[3]); //0x3
        umul_ppmm(scr[9],scr[8],a[1], b[3]); //1x3
        umul_ppmm(scr[11],scr[10],a[2], b[3]); //2x3
        umul_ppmm(scr[13],scr[12],a[3], b[3]); //3x3
     
        umul_ppmm(hi1,lo1,a[1], b[0]); //1x0
        add_ssaaaa(scr[3],scr[2],scr[3],scr[2],hi1,lo1);
        umul_ppmm(hi2,lo2,a[1], b[1]); //1x1
        add_ssaaaa(scr[5],scr[4],scr[5],scr[4],hi2,lo2);
        umul_ppmm(hi1,lo1,a[1], b[2]); //1x2
        add_ssaaaa(scr[7],scr[6],scr[7],scr[6],hi1,lo1);
        umul_ppmm(hi2,lo2,a[2], b[0]); //2x0
        add_ssaaaa(scr[5],scr[4],scr[5],scr[4],hi2,lo2);
        umul_ppmm(hi1,lo1,a[2], b[1]); //2x1
        add_ssaaaa(scr[7],scr[6],scr[7],scr[6],hi1,lo1);
        umul_ppmm(hi2,lo2,a[2], b[2]); //2x2
        add_ssaaaa(scr[9],scr[8],scr[9],scr[8],hi2,lo2);
        umul_ppmm(hi1,lo1,a[3], b[0]); //3x0
        add_ssaaaa(scr[7],scr[6],scr[7],scr[6],hi1,lo1);
        umul_ppmm(hi2,lo2,a[3], b[1]); //3x1
        add_ssaaaa(scr[9],scr[8],scr[9],scr[8],hi2,lo2);
        umul_ppmm(hi1,lo1,a[3], b[2]); //3x2
        add_ssaaaa(scr[11],scr[10],scr[11],scr[10],hi1,lo1);
     
        for (i=0; i<7; i++) 
        {
           res[i]=modp(scr[2*i+1],scr[2*i],d);
           //udiv_qrnnd(quo,res[i],scr[2*i+1],scr[2*i],p);
        }
     }       
}
       
//todo: move karamul_* into Zmodpoly.c

/* 
       Sets res = a*b using iterated Karatsuba multiplication
       assumes length a = length b at the top level of the recursion.
*/
static inline void karamul_modp_recursive(mp_limb_t* res, mp_limb_t* a, mp_limb_t* b, mp_limb_t* scratch, unsigned long size, mp_limb_t p)
{
     mp_limb_t* temp1;
     mp_limb_t* temp;
     
     /*
        As karamul is called recursively we need to ensure that dirty data is not 
        present in the output Zvec res. In fact only one mpz_t needs cleaning.
     */
     res[size-1]=0;
     
     /*
        It's faster if the bottom levels of the recursion are done by dedicated code
        rather than handled as instances of the general case. */
     
     if (size<=4) {
        mul_mod_naieve(res,a,b,scratch,size,p);
        return;
     }
     
      /*
         Set up some fake vecs with half the lengths of the input polynomials
         and point them to the appropriate parts of the input polynomials. 
         a = a1a2, b=b1b2.
     */
     
     mp_limb_t *a1,*a2,*b1,*b2;
     
     a1 = a;
     a2 = a+(size>>1);
     b1 = b;
     b2 = b+(size>>1);
     
     
     /* Allocate space for variables from scratch space
        asum has length size/2
        bsum has length size/2
        prodsum has length size-1                           */
     
     mp_limb_t *asum, *bsum, *prodsum;
     asum = scratch;
     bsum = scratch+(size>>1);
     prodsum = scratch+size;
     
     /*
        (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
     */
     
     karamul_modp_recursive(res,a1,b1,scratch,size>>1,p);
     temp = res+size;
     karamul_modp_recursive(temp,a2,b2,scratch,size>>1,p);
     add_mod(asum,a1,a2,size>>1,p);
     add_mod(bsum,b1,b2,size>>1,p);
     karamul_modp_recursive(prodsum,asum,bsum,scratch+(size<<1)-1,size>>1,p);
     
     sub_mod(prodsum,prodsum,res,size-1,p);
     
     temp = res+size;
     sub_mod(prodsum,prodsum,temp,size-1,p);
     
     temp = res+(size>>1);
     add_mod(temp,temp,prodsum,size-1,p);
}

/*  data1 contains an array of 2^(log_size) mpz_t's
    data1 contains an array of 2^(log_size) mpz_t's  
    we do radixmul with the the given no_primes which
    are stored in the array "primes"                  */

void RadixMul(mpz_t* out, mpz_t* data1, mpz_t* data2, unsigned long log_size,
           unsigned long no_primes, mp_limb_t * primes)
{
   unsigned long i, j;
   unsigned long size = 1 << log_size;
   unsigned long p;
   
   // We need the following memory allocated:
   //   Space for size coefficients each 1 limb in length
   //   for each of the input polys, i.e. 2*size limbs
   
   mp_limb_t* array_all = (mp_limb_t*) limb_alloc(4*size, 0);
   
   mp_limb_t* array1;
   mp_limb_t* array2;
   mp_limb_t* mul_out = array_all+2*size;
   Z_t* temp = Z_alloc();
   Z_t* temp1 = Z_alloc();
   Z_t* prod = Z_alloc();
   Z_t* coeff = Z_alloc();
   
   mp_limb_t* scratch = (mp_limb_t*)limb_alloc(4*size,0);
       
   unsigned long temp2;
   
   Z_set_ui(*prod,1);
   for (i=0; i<no_primes; i++)
   {
       Z_mul_ui(*prod,*prod,primes[i]);
   }
   
   for (i=0; i<no_primes; i++)
   {
       /* p is the current prime modulus*/
       p = primes[i];
       
       array1 = array_all;
       array2 = array_all+size;
       for (j = 0; j < size; j++)
       {
           if (!(j&7)) 
           {
              FLINT_PREFETCH(array1, j+8);
              FLINT_PREFETCH(array2, j+8);
           }
           array1[j] = (mp_limb_t)mpz_fdiv_ui(data1[j],p);
           array2[j] = (mp_limb_t)mpz_fdiv_ui(data2[j],p);
       }
       
       /* precompute coeff = prod/p * (p/prod mod p) */
       mpz_divexact_ui(*temp1,*prod,p);
       temp2 = mpz_fdiv_r_ui(*temp,*temp1,p);
       temp2 = Z_invert_long(temp2,p);
       mpz_mul_ui(*coeff,*temp1,temp2);
       
       /* 
          Compute product of the reduced polynomials and reduce mod p
          TODO, naieve multiplication may be quicker over some regimes here */
       
       karamul_modp_recursive(mul_out,array1,array2,scratch,size,p);
       
       /*
           CRT is done on the fly by multiplying by the precomputed coefficients
           and adding the results to c. */
       
       for (unsigned long j =0; j<2*size-1; j++)
       {
           mpz_mul_ui(*temp1,*coeff,mul_out[j]);
           mpz_add(out[j],out[j],*temp1);
       }
   }
   
   /* Compute and store the final result */
   mpz_fdiv_q_2exp(*temp,*prod,1);
   for (unsigned long j=0; j<2*size-1; j++)
   {
      mpz_mod(out[j],out[j],*temp);
   }
   
   limb_release();
   Z_release();Z_release();Z_release();Z_release();
   limb_release();
}
