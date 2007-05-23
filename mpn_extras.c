/******************************************************************************

 mpn_extras.c
 
 Extra functions for manipulating mpn's and limbs.

 Copyright (C) 2006, William Hart
 
 mp_limb_t mpn_divmod_1_preinv was adapted from GMP.

******************************************************************************/

#include <gmp.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "flint.h"
#include "extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"

/*
   Memory manager to allocate an array of limbs of the given length. It returns a (void*)
   which needs to be typecast to the required object, e.g. mp_limb_t* or mp_limb_t**, etc.
   
   Limbs must be released in the reverse order to that in which they were allocated. 
*/

#define EXPIRE_AFTER 3 // 1/4-th number of allocs before expiring unused allocated space
#define RESALLOC 100 //allocate this many mp_limb_t's at once to save on overheads

inline unsigned long max(unsigned long a, unsigned long b)
{
   return (a<b) ? b : a;
}

typedef struct limb_mem_t //record for managing all allocations of limbs of a certain bitsize
{
   unsigned long remaining; //number of remaining limbs which are available for allocation
   unsigned long length; //total length of array of limbs this record manages
   mp_limb_t* point; //pointer to the next available mp_limb_t
   int expire; //how long till we expire
   int allocated; //1 if some of the space is allocated, 0 otherwise (will not expire if allocated)
   struct limb_mem_t* next; //next record
   struct limb_mem_t* prev; //previous record
} limb_mem_t;

typedef struct limb_memp_t //record for managing a particular allocation of memory for a limb array
{
   limb_mem_t* point; //which record controls this allocation
   unsigned long length; //how many limbs allocated
} limb_memp_t;

limb_mem_t* head_mpn = NULL; //start of doubly linked list of records
limb_mem_t* last_mpn = NULL; //last record in linked list 
limb_memp_t* top_mpn = NULL; //top of stack of limb_memp_t's
limb_memp_t* reservoir_mpn; //array of preallocated limb_memp_t's
unsigned int rescount_mpn=0; //counter for which limb_memp_t we are upto in the reservoir_mpn

// todo: deal with possible out of memory situation when allocating

void* limb_alloc(unsigned long length, int reserved)
{
   printf("%ld limbs allocated\n", length);
   static limb_mem_t* curr;
   static limb_mem_t* temp;
   static int initialised = 0; //has limb_alloc been initialised
   static unsigned int currentalloc = 0; //total number of limb_memp_t's in reservoir_mpn
   static limb_memp_t* tempres;
   static int check=0;
   void* alloc_d;
   
   check++;
   //allocate another block of limb_memp_t's if none are currently allocated, or the reservoir_mpn is depleted
   if (rescount_mpn==currentalloc) // need more limb_memp_t's
   {
      if (!initialised) 
      {
         reservoir_mpn = (limb_memp_t*)malloc(RESALLOC*sizeof(limb_memp_t));
         rescount_mpn=0;
         initialised = 1;
         currentalloc = RESALLOC;
      } else
      {
         //copy old reservoir_mpn into larger one
         tempres = reservoir_mpn;
         reservoir_mpn = (limb_memp_t*)malloc((currentalloc+RESALLOC)*sizeof(limb_memp_t));  
         memcpy(reservoir_mpn,tempres,currentalloc*sizeof(limb_memp_t)); 
         currentalloc+=RESALLOC;
         //free old reservoir_mpn
         free(tempres);  
      }       
   }
   
   curr = head_mpn;
   if (curr != NULL)
   {
      do 
      {
         //search for data of requested size
         if ((curr->remaining >= length) && (curr->remaining < 2*length))//appropriate unallocated space found, so allocate it
         {
             alloc_d = (void*)curr->point;
             curr->point+=length;
             curr->remaining-=length;
             curr->allocated=1;
             top_mpn=&reservoir_mpn[rescount_mpn];
             top_mpn->point=curr;
             top_mpn->length=length;
             
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
                       free(curr->point);
                       if (curr==last_mpn) last_mpn = curr->prev;
                       else curr->next->prev = curr->prev;
                       if (curr==head_mpn) head_mpn=curr->next;
                       else (curr->prev->next = curr->next);
                       temp=curr;
                       curr = curr->next;
                       free(temp);
                    } else curr = curr->next;
                 } else curr = curr->next;
               } while (curr != NULL);
             }
             rescount_mpn++;
             return alloc_d;
         }
         //update expiry information for curr if necessary and possibly expire space
         if (((check&3)==0)&&(!curr->allocated))
         {
            curr->expire--;
            if (curr->expire == 0) 
            {
               free(curr->point);
               if (curr==last_mpn) last_mpn = curr->prev;
               else curr->next->prev = curr->prev;
               if (curr==head_mpn) head_mpn=curr->next;
               else (curr->prev->next = curr->next);
               temp=curr;
               curr = curr->next;
               free(temp);
            } else curr = curr->next;
         } else curr = curr->next;
      } while (curr != NULL);
      //Nothing suitable found, so initialise data of the requested type
      //attach to last_mpn->next
      alloc_d = malloc(length*sizeof(mp_limb_t));
      //set up the new record, and last_mpn to point to this new record
      last_mpn->next = (limb_mem_t*) malloc(sizeof(limb_mem_t));
      last_mpn->next->prev = last_mpn;
      last_mpn=last_mpn->next;
      printf("Here---1\n");
      last_mpn->point = (mp_limb_t*)alloc_d+length;
      last_mpn->next = NULL;
      last_mpn->remaining=0;
      last_mpn->allocated=1;
      last_mpn->length=length;
      top_mpn=&reservoir_mpn[rescount_mpn];
      top_mpn->point=last_mpn;
      top_mpn->length=length;
      rescount_mpn++;
      return alloc_d;
   } 
   /*first time anything has been allocated
     so do the actual allocation of limbs*/
   //set up the new record, and head_mpn and last_mpn to point to this single new record
   alloc_d = malloc(length*sizeof(mp_limb_t));
   head_mpn = (limb_mem_t*) malloc(sizeof(limb_mem_t));
   printf("Here---2\n");
   head_mpn->point = (mp_limb_t*)alloc_d+length;
   head_mpn->next = NULL;
   head_mpn->prev = NULL;
   head_mpn->remaining=0;
   head_mpn->allocated=1;
   head_mpn->length=length;
   last_mpn = head_mpn;
   top_mpn=&reservoir_mpn[rescount_mpn];
   top_mpn->point=head_mpn;
   top_mpn->length=length;
   rescount_mpn++;
   return alloc_d;
}

void limb_release()
{
    unsigned long length = top_mpn->length;
    
    //adjust record to reflect the fact that the limbs have been released back to the stack
    top_mpn->point->point-=length;
    top_mpn->point->remaining+=length;
    //if no limbs of that size are allocated any more set them to expire if not eventually used
    if (top_mpn->point->remaining == top_mpn->point->length) 
    {
       top_mpn->point->allocated=0;
       top_mpn->point->expire=EXPIRE_AFTER;
    }
    //release limb_memp_t back into reservoir_mpn
    top_mpn--;
    rescount_mpn--;
}

/*=======================================================================================*/

/* 
    Performs division by a limb d and places the quotient in qp and returns the 
    remainder. Requires a single limb approximation to 1/d as input. If the most
    significant bit of d is not 1 it expects d to be shifted left (by norm bits)
    until the most significant bit is 1 before the inverse is computed. However 
    the original d should be supplied to the function, not the shifted d. 
    
    This code has been adapted from code found in the GMP package version 4.2.1
    (divrem_1.c)
*/
mp_limb_t mpn_divmod_1_preinv(mp_limb_t * qp, mp_limb_t * up, 
                                  unsigned long un, mp_limb_t d, mp_limb_t dinv, unsigned long norm)
{
  mp_size_t  n;
  mp_size_t  i;
  mp_limb_t  n1, n0;
  mp_limb_t  r = 0;

  n = un;
  if (n == 0)
    return 0;
  
  qp += (n - 1);   /* Make qp point at most significant quotient limb */

  if ((d & (1L<<(FLINT_BITS_PER_LIMB-1))) != 0)
  {
     if (un != 0)
     {
        /* High quotient limb is 0 or 1, skip a divide step. */
	    mp_limb_t q;
	    r = up[un - 1];
	    q = (r >= d);
	    *qp-- = q;
	    r -= (d & -q);
	    n--;
	    un--;
	 }

     /* Multiply-by-inverse, divisor already normalized. */
     for (i = un - 1; i >= 0; i--)
     {
        n0 = up[i];
        udiv_qrnnd_preinv (*qp, r, r, n0, d, dinv);
        qp--;
     }
     return r;
  } else
  {
     /* Most significant bit of divisor == 0.  */
     
     /* Skip a division if high < divisor (high quotient 0).  Testing here
	 before normalizing will still skip as often as possible.  */
     if (un != 0)
	 {
	    n1 = up[un - 1];
	    if (n1 < d)
        {
           r = n1;
	       *qp-- = 0;
	       n--;
	       if (n == 0) return r;
	       un--;
        }
	 }  

     d <<= norm;
     r <<= norm;

     if (un != 0)
     {
        n1 = up[un - 1];
        r |= (n1 >> (FLINT_BITS_PER_LIMB - norm));
        for (i = un - 2; i >= 0; i--)
		{
		  n0 = up[i];
		  udiv_qrnnd_preinv (*qp, r, r, 
				     ((n1 << norm) | (n0 >> (FLINT_BITS_PER_LIMB - norm))), d, dinv);
		  qp--;
		  n1 = n0;
		}
        udiv_qrnnd_preinv (*qp, r, r, n1 << norm, d, dinv);
        qp--;
     }
     
     return r >> norm;
  }
}

mp_limb_t mpn_addmul(mp_limb_t * rp, mp_limb_t * s1p, unsigned long s1n, 
                                      mp_limb_t * s2p, unsigned long s2n)
{
   if (s2n == 0) return 0;
   
   mp_limb_t carry;
   
   carry = mpn_addmul_1(rp, s1p, s1n, s2p[0]);
   for (unsigned long i = 1; i < s2n; i++)
   {
      carry = mpn_add_1(rp+i+s1n-1, rp+i+s1n-1, 1, carry); 
      if (s2p[i]) carry += mpn_addmul_1(rp+i, s1p, s1n, s2p[i]);
   }
   return carry;
}

