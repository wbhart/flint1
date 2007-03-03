/****************************************************************************

Zvec.cpp: Routines for polynomials, implemented as vectors of mpz_t's.

Copyright 2006 William Hart and David Harvey

*****************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include "flint.h"
#include "Zvec.h"
#include "Z.h"
#include "Z-ssmul.h"
#include "ssmul.h"

mpz_t tempZvec1;
mpz_t tempZvec2;
mpz_t tempZvec3;
mpz_t tempZvec4;

unsigned long numPrimes = 0;

/*
    Initializes some temporary mpz_t's to be used by functions in this module.
*/
void initZvec(void)
{
     mpz_init(tempZvec1);
     mpz_init(tempZvec2);
     mpz_init(tempZvec3);
     mpz_init(tempZvec4);
}

/*==========================================================================================

Memory Allocation Routines

============================================================================================*/

#define EXPIRE_AFTER 2 // 1/5-th number of allocs before expiring unused allocated space
#define RESALLOC 100 //allocate this many Zvec_memp_t's at once to save on overheads

typedef struct Zvec_mem_t //record for managing allocations of mpz_t's of a certain bitsize
{
   unsigned long size; //size of coefficients in bits
   unsigned long remaining; //number of remaining mpz_t's which are available for allocation
   unsigned long length; //total length of array of mpz_t's this record manages
   mpz_t* point; //pointer to the next available mpz_t
   int expire; //how long till we expire
   int allocated; //1 if some of the space is allocated, 0 otherwise (will not expire if allocated)
   struct Zvec_mem_t* next; //next record
   struct Zvec_mem_t* prev; //previous record
} Zvec_mem_t;

typedef struct Zvec_memp_t //record for managing a particular allocation of memory for a Zvec
{
   Zvec_mem_t* point;
} Zvec_memp_t;

Zvec_mem_t* head_Zvec = NULL; //start of doubly linked list of records
Zvec_mem_t* last_Zvec = NULL; //last record in linked list 
Zvec_memp_t* top_Zvec = NULL; //top of stack of Zvec_memp_t's
Zvec_memp_t* reservoir_Zvec; //array of preallocated Zvec_memp_t's
unsigned int rescount_Zvec=0; //counter for which Zvec_memp_t we are upto in the reservoir_Zvec

/*
    Same as Zvec_init1 above, but does not set the global pointer 'point', thus this
    function should not be used for initialising a stack of mpz_t's but for initialising
    a new Zvec. Size is set to the number of machine words used to store each coefficient.
*/

void Zvec_init3(Zvec vect, unsigned long length, unsigned long size)
{
     vect->length = length;
     vect->coords = (mpz_t *) calloc(sizeof(mpz_t),length);
     for (unsigned long i=0; i<length; i++)
     {
         mpz_init2(vect->coords[i],size);
     }
}

/*
     Frees memory allocated for a Zvec by Zvec_init2. Note this really frees the
     memory, and doesn't just release it to the stack. So this function should only
     be called when you are completely finished with a Zvec. (See Zvec_init3.)
*/
static inline void Zvec_free2(Zvec vect)
{
     
     for (unsigned long i=0; i<vect->length; i++)
     {
         mpz_clear(vect->coords[i]);
     }
     free(vect->coords);
}

// todo: deal with possible out of memory situation when allocating

void Zvec_alloc(Zvec vect, unsigned long length, unsigned long size, int reserved)
{
   static Zvec_mem_t* curr;
   static Zvec_mem_t* temp;
   static int initialised = 0; //has Zvec_alloc been initialised
   static unsigned int currentalloc = 0; //total number of Zvec_memp_t's in reservoir_Zvec
   static Zvec_memp_t* tempres;
   static int check = 0;
   
   //allocate another block of Zvec_memp_t's if none are currently allocated, or the reservoir_Zvec is depleted
   if (rescount_Zvec==currentalloc) // need more Zvec_memp_t's
   {
      if (!initialised) 
      {
         reservoir_Zvec = (Zvec_memp_t*)malloc(RESALLOC*sizeof(Zvec_memp_t));
         rescount_Zvec=0;
         initialised = 1;
         currentalloc = RESALLOC;
      } else
      {
         //copy old reservoir_Zvec into larger one
         tempres = reservoir_Zvec;
         reservoir_Zvec = (Zvec_memp_t*)malloc((currentalloc+RESALLOC)*sizeof(Zvec_memp_t));  
         memcpy(reservoir_Zvec,tempres,currentalloc*sizeof(Zvec_memp_t)); 
         currentalloc+=RESALLOC;
         //free old reservoir_Zvec
         free(tempres);  
      }       
   }
   
   curr = head_Zvec;
   if (curr != NULL)
   {
      do 
      {
         //search for data of requested type
         if ((curr->size >= size)&&(curr->size <= size+(size>>2))&&(curr->remaining >= length))
         {
             vect->length = length;
             vect->coords = curr->point;
             curr->point+=length;
             curr->remaining-=length;
             curr->allocated=1;
             top_Zvec=&reservoir_Zvec[rescount_Zvec];
             top_Zvec->point=curr;
             //check if any remaining nodes have expired, expire them and return
             if ((check&3)==0)
             {
               do 
               {
                 if (!curr->allocated) 
                 {
                    curr->expire--;
                    if (curr->expire == 0) 
                    {
                       for (unsigned long i = 0; i< curr->remaining; i++) mpz_clear(curr->point[i]);
                       free(curr->point);
                       if (curr==last_Zvec) last_Zvec = curr->prev;
                       else curr->next->prev = curr->prev;
                       if (curr==head_Zvec) head_Zvec=curr->next;
                       else (curr->prev->next = curr->next);
                       temp=curr;
                       curr = curr->next;
                       free(temp);
                    } else curr = curr->next;
                 } else curr = curr->next;
               } while (curr != NULL);
             }
             rescount_Zvec++;
             return;
         }
         //update expiry information for curr if necessary and possibly expire space
         if (((check&3)==0)&&(!curr->allocated)) 
         {
            curr->expire--;
            if (curr->expire == 0) 
            {
               for (unsigned long i = 0; i< curr->remaining; i++) mpz_clear(curr->point[i]);
               free(curr->point);
               if (curr==last_Zvec) last_Zvec = curr->prev;
               else curr->next->prev = curr->prev;
               if (curr==head_Zvec) head_Zvec=curr->next;
               else (curr->prev->next = curr->next);
               temp=curr;
               curr = curr->next;
               free(temp);
            } else curr = curr->next;
         } else curr = curr->next;
      } while (curr != NULL);
      //Nothing suitable found, so initialise data of the requested type
      //attach to last_Zvec->next
      vect->length = length;
      vect->coords = (mpz_t *) malloc(sizeof(mpz_t)*length);  //### maybe should be calloc
      for (unsigned long i=0; i<length; i++)
      {
         mpz_init2(vect->coords[i],size);
      }
      //set up the new record, and last_Zvec to point to this new record
      last_Zvec->next = (Zvec_mem_t*) malloc(sizeof(Zvec_mem_t));
      last_Zvec->next->prev = last_Zvec;
      last_Zvec=last_Zvec->next;
      last_Zvec->size=size;
      last_Zvec->point = vect->coords+length;
      last_Zvec->next = NULL;
      last_Zvec->remaining=0;
      last_Zvec->allocated=1;
      last_Zvec->length=length;
      top_Zvec=&reservoir_Zvec[rescount_Zvec];
      top_Zvec->point=last_Zvec;
      rescount_Zvec++;
      return;
   } 
   /*first time anything has been allocated
     so do the actual allocation of mpz_t's*/
   vect->length = length;
   vect->coords = (mpz_t *) malloc(sizeof(mpz_t)*length);  //### maybe should be calloc
   for (unsigned long i=0; i<length; i++)
   {
      mpz_init2(vect->coords[i],size);
   }
   //set up the new record, and head_Zvec and last_Zvec to point to this single new record
   head_Zvec = (Zvec_mem_t*) malloc(sizeof(Zvec_mem_t));
   head_Zvec->size=size;
   head_Zvec->point = vect->coords+length;
   head_Zvec->next = NULL;
   head_Zvec->prev = NULL;
   head_Zvec->remaining=0;
   head_Zvec->allocated=1;
   head_Zvec->length=length;
   last_Zvec = head_Zvec;
   top_Zvec=&reservoir_Zvec[rescount_Zvec];
   top_Zvec->point=head_Zvec;
   rescount_Zvec++;
}

void Zvec_release(Zvec vec)
{
    //adjust record to reflect the fact that the mpz_t's have been released back to the stack
    top_Zvec->point->point-=vec->length;
    top_Zvec->point->remaining+=vec->length;
    //if no mpz_t's of that size are allocated any more set them to expire if not eventually used
    if (top_Zvec->point->remaining == top_Zvec->point->length) 
    {
       top_Zvec->point->allocated=0;
       top_Zvec->point->expire=EXPIRE_AFTER;
    }
    //release Zvec_memp_t back into reservoir_Zvec
    top_Zvec--;
    rescount_Zvec--;
}

/*==========================================================================================*/

/* 
     Copy the contents of one Zvec into another, i.e. sets Zvec res = a 
*/
static inline void Zvec_copy(Zvec res, Zvec a)
{
     for (unsigned long i = 0; i<a->length; i++)
     {
             mpz_set(res->coords[i],a->coords[i]);
     }
}

/* 
     Sets res = a+b 
     assumes length b<=length a
*/
static inline void Zvec_add(Zvec res, Zvec a, Zvec b)
{
     unsigned long i;
     for (i = 0; i<b->length; i++)
     {
             mpz_add(res->coords[i],a->coords[i],b->coords[i]);
     }
     for (; i<a->length; i++)
     {
             mpz_set(res->coords[i],a->coords[i]);
     }
}

/*
       Sets res = a-b 
       assumes length b<=length a
*/
static inline void Zvec_sub1(Zvec res, Zvec a, Zvec b)
{
     unsigned long i;
     for (i = 0; i<b->length; i++)
     {
             mpz_sub(res->coords[i],a->coords[i],b->coords[i]);
     }
     for (; i<a->length; i++)
     {
             mpz_set(res->coords[i],a->coords[i]);
     }
}

/* 
     Sets res = a*b using the naieve multiplication algorithm
*/
void Zvec_mul_naieve(Zvec res, Zvec a, Zvec b)
{
     for (unsigned long i = 0; i<a->length+b->length-1;i++) mpz_set_ui(res->coords[i],0);
     for (unsigned long i = 0; i<a->length; i++)
     {
         for (unsigned long j = 0; j<b->length; j++)
         {
             mpz_addmul(res->coords[i+j],a->coords[i],b->coords[j]);
         }
     }
     res->length = a->length + b->length -1;
}

/* 
     Sets res = a*b where b is a poly of degree zero
*/
static inline void Zvec_mul1(Zvec res, Zvec a, Zvec b)
{
     for (unsigned long i = 0; i<a->length; i++)
     {
             Z_fast_mul(res->coords[i],a->coords[i],b->coords[0]);
     }
     res->length = a->length + b->length -1;
}

/* 
       Sets res = a*b using iterated Karatsuba multiplication
       assumes length a = length b at the top level of the recursion.
*/
void karamul_recursive(Zvec res, Zvec a, Zvec b, mpz_t* scratch)
{
     Zvec temp1;
     Zvec temp;
     
     /*
        If a or b is a constant simply multiply through by it and return. 
     */
     
     if (a->length == 1) 
     {
        if (b->length ==1)
        {
           Z_fast_mul(res->coords[0],a->coords[0],b->coords[0]);
        } else Zvec_mul1(res,b,a);
        res->length = a->length + b->length - 1;
        return;
     }
     
     if (b->length == 1) 
     {
        Zvec_mul1(res,a,b);
        res->length = a->length + b->length - 1;
        return;
     }
     
     /*
        As karamul is called recursively we need to ensure that dirty data is not 
        present in the output Zvec res. In fact only two mpz_t's need cleaning.
     */

     for (unsigned long i = 0; i < a->length + b->length - 1; i++) mpz_set_ui(res->coords[i],0);
     //mpz_set_ui(res->coords[a->length-1],0);
     //mpz_set_ui(res->coords[a->length],0);
     

     /*
        It's faster if the bottom levels of the recursion are done by dedicated code
        rather than handled as instances of the general case. */
     
     
      if (a->length ==2 && b->length == 2) {
         Z_fast_mul(res->coords[0], a->coords[0], b->coords[0]);
         mpz_add(scratch[0], a->coords[0], a->coords[1]);
         Z_fast_mul(res->coords[2], a->coords[1], b->coords[1]);
         mpz_add(scratch[1], b->coords[0], b->coords[1]);
         Z_fast_mul(res->coords[1], scratch[0], scratch[1]);
         mpz_sub(res->coords[1], res->coords[1], res->coords[0]);
         mpz_sub(res->coords[1], res->coords[1], res->coords[2]);
         res->length = a->length + b->length - 1;
      
         return;
      }
     
      /* This code is faster than the above when coefficients 
         are smaller. Crossover points have not been determined yet. TODO, 
         have karamul decide whether to execute this or the above code for the
         bottom levels of the recursion. */
      
      
      /*if (a->length == 2 && b->length == 2) {
        Z_fast_mul(res->coords[0], a->coords[0], b->coords[0]);
        Z_fast_mul(res->coords[1],a->coords[0],b->coords[1]);
        Z_fast_mul(res->coords[2], a->coords[1], b->coords[1]);
        mpz_addmul(res->coords[1],a->coords[1],b->coords[0]);
        res->length = a->length + b->length - 1;
        return;
      }
      
      if ((a->length < 6) || (b->length < 6)) {
        Zvec_mul_naieve(res,a,b);
        res->length = a->length + b->length - 1;
        return;
      }*/
     
     /*
         Set up some fake Zvecs with half the lengths of the input polynomials
         and point them to the appropriate parts of the input polynomials. 
         a = a1a2, b=b1b2.
     */
     
     Zvec a1,a2,b1,b2;
     
     a1->length = (a->length+1)/2;
     a2->length = a->length-a1->length;
     a1->coords = a->coords;
     a2->coords = a->coords+a1->length;
     
   if (a1->length < b->length) //ordinary case
   {
     b1->length = a1->length;
     b2->length = b->length - b1->length;
     b1->coords = b->coords;
     b2->coords = b->coords + b1->length;   
     
     /* Allocate space for temporary Zvecs from the stack.*/
     
     Zvec asum, bsum, prodsum;
     asum->length = a1->length;
     asum->coords = scratch;
     bsum->length = a1->length;
     bsum->coords = scratch+a1->length;
     prodsum->length = (a1->length<<1)-1;
     prodsum->coords = scratch+(a1->length<<1);
     
     /*
        (a1+a2*x)*(b1+b2*x) = a1*b1 + a2*b2*x^2 + (a1+a2)*(b1+b2)*x-a1*b1*x-a2*b2*x;
     */
     
     // res_lo = a1*b1
     karamul_recursive(res,a1,b1,scratch+(a1->length<<2)-1);
     // res_hi = a2*b2
     temp->coords = res->coords+(a1->length<<1);
     karamul_recursive(temp,a2,b2,scratch+(a1->length<<2)-1);
     // asum = a1+a2
     Zvec_add(asum,a1,a2);
     // bsum = b1+b2
     Zvec_add(bsum,b1,b2);
     // prodsum = asum*bsum
     karamul_recursive(prodsum,asum,bsum,scratch+(a1->length<<2)-1);
     
     // prodsum = prodsum - res_lo
     temp->coords = res->coords;
     temp->length = (a1->length<<1)-1;
     Zvec_sub1(prodsum,prodsum,temp);
     
     // prodsum = prodsum - res_hi
     temp->coords = res->coords+(a1->length<<1);
     temp->length = a2->length+b2->length-1;
     Zvec_sub1(prodsum,prodsum,temp);
     
     // res_mid += prodsum
     temp->coords = res->coords+a1->length;
     temp->length=prodsum->length;
     Zvec_add(temp,temp,prodsum);
     
     res->length = a->length + b->length - 1;
     
   } else //degenerate case
   {
      // res_lo = a1*b
      karamul_recursive(res,a1,b,scratch);
      
      //temp = a2*b
      temp->coords = scratch;
      temp->length = a2->length + b->length - 1;
      if (b->length <= a2->length) karamul_recursive(temp,a2,b,scratch+temp->length);
      else karamul_recursive(temp,b,a2,scratch+temp->length);
      
      // res_mid += temp
      temp1->coords = res->coords+a1->length;
      temp1->length = temp->length;
      Zvec_add(temp1,temp1,temp);
      
      res->length = a->length + b->length - 1;
   }
}

/* Set res to the polynomial a*b
   Size should be the maximum output coefficient bit size 
   Assumes length a >= length b
*/

void Zvec_karamul(Zvec res, Zvec a, Zvec b, unsigned long bits)
{
     unsigned long log_length = 0;
     while ((1<<log_length) < a->length) log_length++;
     
     Zvec scratch;
     Zvec_alloc(scratch,5*a->length,2*bits+log_length,0);
     
     karamul_recursive(res,a,b,scratch->coords);
     
     Zvec_release(scratch);
}
     

/*
     Set vec = 0, i.e. zero out all coefficients
*/ 
void Zvec_clear(Zvec vec)
{
     for (unsigned long i = 0; i<vec->length; i++)
     {
         mpz_set_ui(vec->coords[i],0);
     }
}

/*
Returns the maximum size in bits of the coefficients of the Zvec a and if
sign is not NULL, it sets sign to -1 if the coefficients are signed, and 
sets it to 1 if they are all non-negative.
*/
unsigned long Zvec_max_size(int* sign, Zvec a, Zvec b)
{
   unsigned long max_size = 0;
   unsigned long new_size;
   int int_sign = 0;
   
   for (unsigned long i = 0; i < a->length; i++)
   {
       if (mpz_sgn(a->coords[i]) < 0) int_sign = 1;
       new_size = mpz_size(a->coords[i]);
       if (new_size >= max_size)
       {
          max_size = new_size;
       }
   }   
   for (unsigned long i = 0; i < b->length; i++)
   {
       if (mpz_sgn(b->coords[i]) < 0) int_sign = 1;
       new_size = mpz_size(b->coords[i]);
       if (new_size >= max_size)
       {
          max_size = new_size;
       }
   }   
   
   if (sign != NULL) *sign = int_sign;
   
   return max_size;   
}

unsigned long Zvec_max_length(Zvec a, Zvec b)
{
   return (a->length > b->length) ? a->length : b->length;
}

/*
     Set res = a*b using radix multiplication algorithm. Similar to NTL's HomMul. 
     Size must be set to the maximum size of the output coefficients (in bits).
*/
void Zvec_radixMul(Zvec res, Zvec a, Zvec b, unsigned long size)
{
   static unsigned long primes[2000]; // Primes that can be used in by RadixMul
   
   mpz_t prod; mpz_init(prod);
   mpz_t temp1; mpz_init(temp1);
   mpz_t temp; mpz_init(temp);
   mpz_t coeff; mpz_init(coeff);
   
   unsigned long temp2;
    
   Zvec A, B, C; // Zvec's to hold a, b and the product modulo a prime p
   
   /* 
      Generate primes to be used. Note they are only computed once, not every time 
      radixmul is called.
   */
   
   unsigned long prime = (65536UL*32768UL*65536UL*65536UL); //Primes start here, i.e. they are larger than this value
   while (numPrimes < 2000)
   {
      prime = Z_nextprime_long(prime);
      primes[numPrimes] = prime;
      numPrimes++;
   } 
   
   /* Find primes whose product is bigger than the size of the output coefficients */
   
   unsigned long primecount=0; // Actual number of primes that will be used this time.
   mpz_set_ui(prod,1);
   while (mpz_sizeinbase(prod,2)<size)
   {
         mpz_mul_ui(prod,prod,primes[primecount]);
         primecount++;
   }
   
   /* 
      Allocate a temporary Zvec capable of holding the increments for the CRT, which is 
      basically computed as you go with this algorithm.
   */
   
   Zvec c;
   c->length = a->length + b->length - 1;
   Zvec_alloc(c,c->length,size,0);
   Zvec_clear(c);
   
   /* Initialize the Zvecs which will hold a, b and c modulo a prime p */
   
   Zvec_init3(A,a->length,1);
   Zvec_init3(B,b->length,1);
   Zvec_init3(C,c->length,1);
   
   for (unsigned long i=0; i<primecount; i++) // for each prime
   {
       /* p is the current prime modulus*/
       unsigned long p = primes[i];
       
       /* precompute coeff = prod/p * (p/prod mod p) */
       mpz_divexact_ui(temp1,prod,p);
       temp2 = mpz_fdiv_r_ui(temp,temp1,p);
       temp2 = Z_invert_long(temp2,p);
       mpz_mul_ui(coeff,temp1,temp2);
       
       /* do reductions of a and b modulo p */
       for (unsigned long j = 0; j<a->length; j++)
       {
           mpz_mod_ui(A->coords[j],a->coords[j],p);
       }
       
       for (unsigned long j = 0; j<b->length; j++)
       {
           mpz_mod_ui(B->coords[j],b->coords[j],p);
       }
       
       /* 
          Compute product of the reduced polynomials and reduce mod p
          TODO, naieve multiplication may be quicker over some regimes here
       */
       Zvec_karamul(C,A,B,size);
       for (unsigned long j=0; j<C->length; j++) mpz_fdiv_r_ui(C->coords[j],C->coords[j],p);
       
       /*
           CRT is done on the fly by multiplying by by the precomputed coefficients
           and adding the results to c.
       */
       for (unsigned long j =0; j<C->length; j++)
       {
           Z_fast_mul(temp1,coeff,C->coords[j]);
           mpz_add(c->coords[j],c->coords[j],temp1);
       }
   }
       
   /* make sure the length of the final result is set correctly */
   res->length = a->length + b->length - 1;
   
   /* Compute and store the final result */
   mpz_fdiv_q_2exp(temp,prod,1);
   for (unsigned long j=0; j<a->length+b->length-1; j++)
   {
      mpz_mod(res->coords[j],c->coords[j],temp);
   }
    
   /* Free temporarily allocated space */
   Zvec_free2(A); Zvec_free2(B); Zvec_free2(C);
   Zvec_release(c);
   mpz_clear(prod); mpz_clear(temp1); mpz_clear(temp);  
}

/*
     a = b mod 2^n
*/
static inline void Truncate(mpz_t a, mpz_t b, unsigned long n)
{
     mpz_fdiv_r_2exp(a,b,n);
}
          
/*
    a = b * 2^l mod p where p = 2^n+1
    We assume 0 <= l < p
*/
static inline void LeftRotate(mpz_t a, mpz_t b, unsigned long l, mpz_t p, unsigned long n)
{
     if (l==0)
     {
        if (a!=b)
        {
           mpz_set(a,b);
        }
        return;
     }
     
     mpz_div_2exp(tempZvec4,b,n-l);
     Truncate(a,b,n-l);
     mpz_mul_2exp(a,a,l);
     mpz_sub(a,a,tempZvec4);
     if (mpz_sgn(a)<0) mpz_add(a,a,p);
}

/*
    a = b * 2^l mod p where p = 2^n+1
    We make no assumptions about l. In particular it may be negative!!
*/
static inline void Rotate(mpz_t a, mpz_t b, long l, mpz_t p, unsigned long n)
{  
   if (mpz_cmp_ui(b,0) == 0) 
   {
      mpz_set_ui(a,0);
      return;
   }
   
   if (l >= 0) l %= (n << 1); 
   else l = (n << 1) - 1 - (-(l + 1) % (n << 1));

   if (l < (long) n) LeftRotate(a, b, l, p, n);
   else 
   {
     LeftRotate(a, b, l - n, p, n);
     mpz_sub(a, p, a);
   }
}

void Zvec_fft(Zvec a, unsigned long r, unsigned long l, mpz_t p, unsigned long n)
{
     unsigned long round;
     long i;
     unsigned long off, j, e;
     unsigned long halfsize;
     
     for (round = 0; round < l; round++, r <<= 1) 
     {
         halfsize =  1UL << (l - 1 - round);
         for (i = (1L << round) - 1, off = 0; i >= 0; i--, off += halfsize) 
         {
             for (j = 0, e = 0; j < halfsize; j++, off++, e+=r) 
             {
                 mpz_sub(tempZvec1, a->coords[off], a->coords[off + halfsize]);
                 if (mpz_sgn(tempZvec1) < 0) mpz_add(tempZvec1, tempZvec1, p);
	             mpz_add(a->coords[off], a->coords[off], a->coords[off + halfsize]);
	             mpz_sub(tempZvec2, a->coords[off], p);
	             if (mpz_sgn(tempZvec2) >= 0) mpz_set(a->coords[off],tempZvec2);
	             Rotate(a->coords[off + halfsize], tempZvec1, e, p, n);
             }
         }    
     } 
}

void Zvec_fft3(mpz_t* a, unsigned long l, unsigned long m, mpz_t p, unsigned long n, unsigned long r)
{
     unsigned long mmulinc=m/l;
     mmulinc*=r;
     
     if (l<=1) return;
     
     if (l==2)
     {
        mpz_sub (tempZvec1, a[0], a[1]); 
        mpz_add (a[0], a[0], a[1]);  
        
        if (mpz_cmp(a[0],p)>=0) mpz_sub(a[0],a[0],p);
        //mpz_mod(a[1], tempZvec1, p);    
        if (mpz_sgn(tempZvec1) < 0) mpz_add(a[1], tempZvec1, p);
        else mpz_set(a[1], tempZvec1);
        
        //mpz_mod(a[0], a[0], p); 
        
        return;  
     }
     
     l/=2;
     mpz_sub (tempZvec1, a[0], a[l]); 
     mpz_add (a[0], a[0], a[l]);  
     
     //mpz_mod(a[l], tempZvec1, p);    
     if (mpz_sgn(tempZvec1) < 0) mpz_add(a[l], tempZvec1, p);
     else mpz_set(a[l], tempZvec1);
     
     //mpz_mod(a[0], a[0], p);
     if (mpz_cmp(a[0],p)>=0) mpz_sub(a[0],a[0],p);
        
     for (unsigned long j=1, mmul=mmulinc; j < l; j++, mmul+=mmulinc) {
        mpz_sub(tempZvec1,a[j],a[j+l]);
        mpz_add(a[j],a[j],a[j+l]);
        
        if (mpz_cmp(a[j],p)>=0) mpz_sub(a[j],a[j],p); // this appears to be faster
        Rotate(a[j+l],tempZvec1,mmul,p,n);
        //mpz_mul_2exp(a[j+l],tempZvec1,r*mmul);
        //mpz_mod(a[j+l],a[j+l],p);
        
        //mpz_mod(a[j],a[j],p);
        
        //mpz_sub(tempZvec1,a[j],p);
        //if (mpz_sgn(tempZvec1)>=0) mpz_set(a[j],tempZvec1);
            //TODO check that it is always faster!!
     }
     Zvec_fft3(a+l,l,m,p,n,r);
     Zvec_fft3(a,l,m,p,n,r);
     
}

void Zvec_ifft(Zvec a, unsigned long r, unsigned long l, mpz_t p, unsigned long n)
{
  long round;
  long i, e;
  unsigned long j, off;
  unsigned long halfsize;
  
  for (round = l - 1, r <<= l - 1; round >= 0; round--, r >>= 1) {
    halfsize = 1UL << (l - 1 - round);
    for (i = (1L << round) - 1, off = 0; i >= 0; i--, off += halfsize) {
      for (j = 0, e = 0; j < halfsize; j++, off++, e+=r) {
         Rotate(a->coords[off + halfsize], a->coords[off + halfsize], -e, p, n);
         mpz_sub(tempZvec2, a->coords[off], a->coords[off + halfsize]);
         mpz_add(a->coords[off], a->coords[off], a->coords[off + halfsize]);
         mpz_sub(tempZvec1, a->coords[off], p);
         if (mpz_sgn(tempZvec1) >= 0) mpz_set(a->coords[off],tempZvec1);
	     if (mpz_sgn(tempZvec2) < 0) mpz_add(a->coords[off+halfsize], tempZvec2, p);
	     else mpz_set(a->coords[off+halfsize],tempZvec2);
      }
    }
  }
}

void Zvec_ifft3(mpz_t* a, unsigned long l, unsigned long m, mpz_t p, unsigned long n, unsigned long r)
{
     unsigned long mmulinc=m/l;
     mmulinc*=r;
     
     if (l<=1) return;
     
     if (l==2)
     {
        mpz_sub (tempZvec1, a[0], a[1]); 
        mpz_add (a[0], a[0], a[1]);  
        
        //mpz_mod(a[1], tempZvec1, p); 
        if (mpz_sgn(tempZvec1) < 0) mpz_add(a[1], tempZvec1, p);
        else mpz_set(a[1], tempZvec1);   
        
        //mpz_mod(a[0], a[0], p);   
        if (mpz_cmp(a[0],p)>=0) mpz_sub(a[0],a[0],p);
        
        return;
     }
     
     l/=2;
     
     Zvec_ifft3(a+l,l,m,p,n,r);
     Zvec_ifft3(a,l,m,p,n,r);
     
     mpz_sub (tempZvec1, a[0], a[l]); 
     mpz_add (a[0], a[0], a[l]);  
     
     if (mpz_cmp(a[0],p)>=0) mpz_sub(a[0],a[0],p);
     //mpz_mod(a[l], tempZvec1, p);    
     if (mpz_sgn(tempZvec1) < 0) mpz_add(a[l], tempZvec1, p);
     else mpz_set(a[l], tempZvec1);   
     
     //mpz_mod(a[0], a[0], p);
     
     for (unsigned long j=1, mmul=r*m-mmulinc; j < l; j++, mmul-=mmulinc) {
        
        //mpz_mul_2exp(a[j+l],a[j+l],r*mmul);
        //mpz_mod(a[j+l],a[j+l],p);
        Rotate(a[j+l],a[j+l],mmul,p,n);
        
        mpz_sub(tempZvec1,a[j],a[j+l]);
        mpz_add(a[j],a[j],a[j+l]);
        
        if (mpz_cmp(a[j],p)>=0) mpz_sub(a[j],a[j],p); //this appears to be faster
        //mpz_mod(a[j+l],tempZvec1,p);
        if (mpz_sgn(tempZvec1)<0) mpz_add(a[j+l],tempZvec1,p);
        else mpz_set(a[j+l],tempZvec1);
             
        //mpz_mod(a[j],a[j],p);
        //mpz_sub(tempZvec1,a[j],p);
        //if (mpz_sgn(tempZvec1)>=0) mpz_set(a[j],tempZvec1);
     }
}


/*
     Set res = a*b using Schoenhage-Strassen. This works modulo a prime p = 2^mr+1 where
     m is 2^l. This allows a convolution of size 2^(l+1) since 2 is a 2^(l+1)-th
     root of unity modulo p. The output coefficients must be smaller than p and the degree
     of the output poly (which is assumed to be the same - they can be padded if needs be)
     must be smaller than 2^(l+1).
     
     Type can be 0 = iterative FFT's or 1 = recursive FFT's
     TODO see if iterative FFT is ever faster.
     
     Size is the size in bits of the input coefficients.
*/
void Zvec_SSMul(Zvec res, Zvec a, Zvec b, unsigned long size, int type)
{
     /* 
        Compute appropriate l, m2 = 2^(l+1), bound = size of output coefficients, r and mr.
        Note m = 2^l.
     */
     mpz_set_ui(tempZvec1,a->length+b->length-1);
     unsigned long l = mpz_sizeinbase(tempZvec1,2)-1;
     unsigned long m2 = 1L << (l+1);
     unsigned long templ1 = ((a->length < b->length) ? a->length - 1 : b->length - 1);
     mpz_set_ui(tempZvec1,templ1);
     templ1 = mpz_sizeinbase(tempZvec1,2);
     unsigned long bound = 2 + templ1 + 2*size;
     unsigned long r = (bound >> l) + 1;
     unsigned long mr = r << l;
     
     /* Allocate space for and compute p = 2^mr */
     mpz_t p;
     mpz_init(p);
     mpz_set_ui(p,1);
     mpz_mul_2exp(p,p,mr);
     mpz_add_ui(p,p,1);
     
     /* 
        Allocate temporary space to hold the result of the FFT's. We don't want to destroy 
        the input values as the FFT's are computed in place.
     */
     
     Zvec aa, bb;
     Zvec_alloc(aa,m2,size,0);
     Zvec_alloc(bb,m2,size,0);
     /* Adjust coefficients to deal with the case where they are negative. The FFT wants
        them all positive.
     */
     unsigned long i;
     for (i=0; i<b->length; i++)
     {
         if (mpz_cmp_ui(b->coords[i],0) < 0) mpz_add(bb->coords[i],b->coords[i],p);
         else mpz_set(bb->coords[i],b->coords[i]);
     }
     for (i=b->length;i<bb->length;i++) mpz_set_ui(bb->coords[i],0); // Zero possibly dirty mpz_t's
     for (i=0; i<a->length; i++)
     {
         if (mpz_cmp_ui(a->coords[i],0) < 0) mpz_add(aa->coords[i],a->coords[i],p);
         else mpz_set(aa->coords[i],a->coords[i]);
     }
     for (i=a->length;i<aa->length;i++) mpz_set_ui(aa->coords[i],0); // Zero possibly dirty mpz_t's
     
     /* perform forward FFT */
     if (type==1)
     {
        Zvec_fft3(aa->coords, m2, m2, p, mr, r);
        Zvec_fft3(bb->coords, m2, m2, p, mr, r);
     } else
     {
        Zvec_fft(aa, r, l+1, p, mr);
        Zvec_fft(bb, r, l+1, p, mr);
     }
     /* perform pointwise mults */
     for (i=0; i<m2; i++)
     {
         Z_fast_mul(tempZvec2,aa->coords[i],bb->coords[i]);
         if (mpz_sizeinbase(tempZvec2,2) >= mr)
         {
            mpz_div_2exp(tempZvec1,tempZvec2,mr);
            Truncate(tempZvec2,tempZvec2,mr);
            mpz_sub(tempZvec2,tempZvec2,tempZvec1);
            
            //mpz_mod(tempZvec2,tempZvec2,p);
            if (mpz_sgn(tempZvec2)<0) mpz_add(tempZvec2,tempZvec2,p);
         }
         mpz_set(aa->coords[i],tempZvec2);
     }
     
     /* perform inverse FFT */
     
     if (type==1) 
     {
        Zvec_ifft3(aa->coords, m2, m2, p, mr, r);
     } else 
     {
        Zvec_ifft(aa, r, l + 1, p, mr);
     }
     
     /* Scale and adjust so coefficients have the right sign */
     for (i=0; i< a->length+b->length-1; i++)
     {
         mpz_set(tempZvec2,aa->coords[i]);
         if(mpz_sgn(tempZvec2)!=0)
         {
            LeftRotate(tempZvec2, tempZvec2, mr - l - 1, p, mr);
            mpz_sub(tempZvec1, p, tempZvec2);
            if (mpz_sizeinbase(tempZvec1,2) >= mr) {
	           mpz_neg(res->coords[i], tempZvec2); 
            } else mpz_set(res->coords[i],tempZvec1);
         } else mpz_set_ui(res->coords[i],0);
     }
     
     /* Free temporary space */
     
     Zvec_release(bb);
     Zvec_release(aa);
     mpz_clear(p);   
}

void Zvec_mul(Zvec res, Zvec a, Zvec b)
{
   int sign;
   unsigned long max_bits;
   unsigned long max_size;
   mpz_t* acoords = a->coords;
   mpz_t* bcoords = b->coords;
   unsigned long i;
   unsigned long log_length = 0;
   while ((1<<log_length) < a->length) log_length++;
   while ((1<<log_length) < b->length) log_length++;
          
   if ((a->length < 16384) && (b->length < 16384))
   {
      max_size = Zvec_max_size(&sign,a,b);
      if (max_size == 1)
      {
         unsigned long mask = -2UL;
         max_bits = 1;
         {
            for (i = 0; i < a->length; i++)
            {
               while (*(acoords[i]->_mp_d) & mask) 
               {
                  mask <<= 1;
                  max_bits++;
               }
            }
            for (i = 0; i < b->length; i++)
            {
               while (*(bcoords[i]->_mp_d) & mask) 
               {
                  mask <<= 1;
                  max_bits++;
               }
            }
         }
      } else
      {
         unsigned long new_bits;
         max_bits = FLINT_BITS_PER_LIMB+1;
         for (i = 0; i < a->length; i++)
         {
            new_bits = mpz_sizeinbase(acoords[i],2);
            if (new_bits > max_bits) max_bits = new_bits; 
         }
         for (i = 0; i < b->length; i++)
         {
            new_bits = mpz_sizeinbase(bcoords[i],2);
            if (new_bits > max_bits) max_bits = new_bits; 
         }
      }
   } else
   {
      sign = 0;
      unsigned long new_bits;
      max_bits = 1;
      for (i = 0; i < a->length; i++)
      {
         if (mpz_sgn(acoords[i])<0) sign = 1;
         new_bits = mpz_sizeinbase(acoords[i],2);
         if (new_bits > max_bits) max_bits = new_bits; 
      }
      for (i = 0; i < b->length; i++)
      {
         if (mpz_sgn(acoords[i])<0) sign = 1;
         new_bits = mpz_sizeinbase(bcoords[i],2);
         if (new_bits > max_bits) max_bits = new_bits; 
      }
   }
   
   if ((Zvec_max_length(a,b) < 10) && ((max_bits<<1)+log_length+(sign<<1) >= FLINT_BITS_PER_LIMB))
   {
      if (a->length >= b->length) Zvec_karamul(res, a, b, max_bits);
      else Zvec_karamul(res, b, a, max_bits);
   } else
   {
      if (a->length >= b->length) SSMul(res, a, b, max_bits, sign);
      else SSMul(res, b, a, max_bits, sign);
   }
}

