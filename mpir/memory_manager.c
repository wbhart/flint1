/*============================================================================

    This file is part of MPIR.

    MPIR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    MPIR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MPIR; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/
/****************************************************************************

memory-manager.c: MPIR-wide memory management

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "mpir.h"
#include "memory_manager.h"

#define DEBUG 0 // Switches to debugging memory allocators
#define DEBUG_PRINT 0 // Prints info about all allocations and releases

#if DEBUG

#define OVERRUN 100
#define CHECK_CHAR 0x55

/*-----------------------------------------------------------------------------------------------
 
    Debug version of stack based memory managers
    
------------------------------------------------------------------------------------------------*/

void * mempts[200000];
ulong upto = 0;

void * mpir_stack_alloc(ulong bytes)
{
#if DEBUG_PRINT
   printf("Allocating %ld bytes\n",  bytes);
#endif
   if (upto == 200000) 
   {
      printf("Error: no free stack nodes in mpir_stack_alloc\n");
      abort();         
   }
   char * block = malloc(bytes+OVERRUN+sizeof(ulong));
   if (block == NULL)
   {
      printf("Error: unable to allocate memory in mpir_stack_alloc\n");
      abort();         
   }
   mempts[upto] = (void *) block;
   upto++;
   ulong * ptr = (ulong *) block;
   ptr[0] = bytes;
   for (ulong i = 0; i < 100; i++)
      block[bytes+i+sizeof(ulong)] = CHECK_CHAR;
   return (void*) (block+sizeof(ulong));
}

void mpir_stack_release()
{
   if (upto == 0) 
   {
      printf("Error: attempt to free unallocated block in mpir_stack_release\n");
      abort();         
   }
   upto--;
   char * block = mempts[upto];
   ulong * ptr = (ulong *) block;
   ulong bytes = (ulong) ptr[0];
#if DEBUG_PRINT
   printf("Releasing %ld bytes\n", bytes);
#endif
   for (ulong i = 0; i < 100; i++)
      if (block[bytes+i+sizeof(ulong)] != CHECK_CHAR) 
      {
         printf("Error: Block overrun detected by stack memory allocator!!\n");
         abort();
      }
   free(mempts[upto]);
}

void * mpir_stack_alloc_small(ulong bytes)
{
   return mpir_stack_alloc(bytes);
}

void mpir_stack_release_small(void)
{
   mpir_stack_release();
}

void mpir_stack_cleanup(void)
{
   if (upto) printf("Error: stack allocator detected mismatch on cleanup!\n"); 
}

#else

/*
   Stack based memory manager to allocate an array of bytes of the given length. 
   It returns a (void*) which needs to be typecast to the required object, 
   e.g. mp_limb_t* or char**, etc.
   
   Bytes must be released in the reverse order to that in which they were allocated. 
*/

#define EXPIRE_AFTER 3 // 1/4-th number of allocs before expiring unused allocated space
#define RESALLOC 100 // allocate this many bytes at once to save on overheads

typedef struct byte_mem_t // record for managing all allocations of bytes of a certain size
{
   ulong remaining; // number of remaining bytes which are available for allocation
   ulong bytes; // total length of array of bytes this record manages
   char * point; // pointer to the next available byte
   int expire; // how long till we expire
   int allocated; // 1 if some of the space is allocated, 0 otherwise (will not expire if allocated)
   struct byte_mem_t* next; //next record
   struct byte_mem_t* prev; //previous record
} byte_mem_t;

typedef struct byte_memp_t //record for managing a particular allocation of memory for a byte array
{
   byte_mem_t* point; // which record controls this allocation
   ulong bytes; // how many bytes allocated
} byte_memp_t;

byte_mem_t* head = NULL; // start of doubly linked list of records
byte_mem_t* last = NULL; // last record in linked list 
byte_memp_t* top = NULL; // top of stack of byte_memp_t's
byte_memp_t* reservoir; // array of preallocated byte_memp_t's
ulong rescount = 0; // counter for which byte_memp_t we are upto in the reservoir

void* mpir_stack_alloc(ulong bytes)
{
   static byte_mem_t* curr;
   static byte_mem_t* temp;
   static int initialised = 0; // has stack_alloc been initialised
   static ulong currentalloc = 0; // total number of byte_memp_t's in reservoir
   static byte_memp_t* tempres;
   static int check = 0;
   void* alloc_d;
   
   check++;
   // allocate another block of byte_memp_t's if none are currently allocated, or the reservoir is depleted
   if (rescount == currentalloc) // need more byte_memp_t's
   {
      if (!initialised) 
      {
         reservoir = (byte_memp_t*) malloc(RESALLOC*sizeof(byte_memp_t));
         rescount = 0;
         initialised = 1;
         currentalloc = RESALLOC;
      } else
      {
         // copy old reservoir into larger one
         tempres = reservoir;
         reservoir = (byte_memp_t*) malloc((currentalloc+RESALLOC)*sizeof(byte_memp_t));  
         memcpy(reservoir, tempres, currentalloc*sizeof(byte_memp_t)); 
         currentalloc+=RESALLOC;
         // free old reservoir
         free(tempres);  
      }       
   }
   
   curr = head;
   if (curr != NULL)
   {
      do 
      {
         // search for data of requested size
         if ((curr->remaining >= bytes) && (curr->remaining < 2*bytes)) // appropriate unallocated space found, so allocate it
         {
             alloc_d = (void*)curr->point;
             curr->point+=bytes;
             curr->remaining-=bytes;
             curr->allocated = 1;
             top = reservoir + rescount;
             top->point = curr;
             top->bytes = bytes;
             
             // check if any remaining nodes have expired, expire them and return
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
                       if (curr==last) last = curr->prev;
                       else curr->next->prev = curr->prev;
                       if (curr==head) head=curr->next;
                       else (curr->prev->next = curr->next);
                       temp=curr;
                       curr = curr->next;
                       free(temp);
                    } else curr = curr->next;
                 } else curr = curr->next;
               } while (curr != NULL);
             }
             rescount++;
             return alloc_d;
         }
         // update expiry information for curr if necessary and possibly expire space
         if (((check&3)==0)&&(!curr->allocated))
         {
            curr->expire--;
            if (curr->expire == 0) 
            {
               free(curr->point);
               if (curr==last) last = curr->prev;
               else curr->next->prev = curr->prev;
               if (curr==head) head=curr->next;
               else (curr->prev->next = curr->next);
               temp=curr;
               curr = curr->next;
               free(temp);
            } else curr = curr->next;
         } else curr = curr->next;
      } while (curr != NULL);
      // Nothing suitable found, so initialise data of the requested type
      // attach to last->next
      alloc_d = malloc(bytes);
      // set up the new record, and last to point to this new record
      last->next = (byte_mem_t*) malloc(sizeof(byte_mem_t));
      last->next->prev = last;
      last = last->next;
      last->point = (char*)alloc_d + bytes;
      last->next = NULL;
      last->remaining = 0;
      last->allocated = 1;
      last->bytes = bytes;
      top = reservoir + rescount;
      top->point = last;
      top->bytes = bytes;
      rescount++;
      return alloc_d;
   } 
   /* first time anything has been allocated
      so do the actual allocation of bytes */
   // set up the new record, and head and last to point to this single new record
   alloc_d = malloc(bytes);
   head = (byte_mem_t*) malloc(sizeof(byte_mem_t));
   head->point = (char*)alloc_d + bytes;
   head->next = NULL;
   head->prev = NULL;
   head->remaining = 0;
   head->allocated = 1;
   head->bytes = bytes;
   last = head;
   top = reservoir + rescount;
   top->point = head;
   top->bytes = bytes;
   rescount++;
   return alloc_d;
}

void mpir_stack_release()
{
    ulong bytes = top->bytes;
    
    // adjust record to reflect the fact that the bytes have been released back to the stack
    top->point->point-=bytes;
    top->point->remaining+=bytes;
    // if no blocks of that size are allocated any more set them to expire if not eventually used
    if (top->point->remaining == top->point->bytes) 
    {
       top->point->allocated = 0;
       top->point->expire = EXPIRE_AFTER;
    }
    // release byte_memp_t back into reservoir
    top--;
    rescount--;
}


/*-----------------------------------------------------------------------------------------------*/

#define MPIR_SMALL_BLOCK_SIZE 10000L

char * block_ptr = NULL;
ulong block_left = 0;
   
void * mpir_stack_alloc_small(ulong bytes)
{
   if (bytes + sizeof(ulong) > block_left) // not enough space left, allocate a new block
   {
      if (bytes + 2*sizeof(ulong) + sizeof(char *) >= MPIR_SMALL_BLOCK_SIZE)
      {
         printf("Error: attempt to allocate %ld bytes in small stack memory manager!\n", bytes);
         abort();
      }
      if (block_ptr == NULL)
      {
         block_ptr = (char *) mpir_alloc(MPIR_SMALL_BLOCK_SIZE);
         block_left = MPIR_SMALL_BLOCK_SIZE - 2*sizeof(ulong);
         ulong * ptr = (ulong *) block_ptr;
         ptr[0] = 0UL;
         char ** ptr2 = (char **) (ptr + 1);
         ptr2[0] = NULL;
         block_ptr += sizeof(ulong) + sizeof(char *);
      } else
      {
         char * temp = (char *) mpir_alloc(MPIR_SMALL_BLOCK_SIZE);
         ulong * ptr = (ulong *) temp;
         ptr[0] = block_left;
         char ** ptr2 = (char **) (ptr + 1);
         ptr2[0] = (char *) block_ptr; 
         block_ptr = temp + sizeof(ulong) + sizeof(char *);
         block_left = MPIR_SMALL_BLOCK_SIZE - sizeof(ulong) - sizeof(char *);
      }
   }
   
   ulong * ptr = (ulong *) (block_ptr + bytes);
   ptr[0] = bytes;
   block_ptr += (bytes+sizeof(ulong)); 
   block_left -= (bytes+sizeof(ulong));
   return (void *) (block_ptr - (bytes + sizeof(ulong))); 
}

void mpir_stack_release_small(void)
{
   if (block_left == MPIR_SMALL_BLOCK_SIZE - sizeof(ulong) - sizeof(char *))
   {
      block_ptr -= (sizeof(ulong) - sizeof(char *));
      ulong * ptr = (ulong *) block_ptr;
      block_left = ptr[0];
      char * temp = block_ptr;
      ulong * ptr2 = (ulong *) (ptr + 1);
      block_ptr = (char *) ptr2[0]; 
      mpir_free(temp);           
   } 
   
   ulong * ptr = (ulong *) block_ptr;
   ptr--;
   ulong temp = (*ptr);
   block_left += (temp+sizeof(ulong));
   block_ptr -= (temp + sizeof(ulong));
}

void mpir_stack_cleanup()
{
   byte_mem_t* curr = head;
   byte_mem_t* temp;
   
   if (curr != NULL)
   {
      do 
      {
         if (curr->allocated) 
         {
            printf("Warning: MPIR stack memory allocation cleanup detected mismatched allocation/releases\n"); 
         } 
         free(curr->point);
         if (curr==last) last = curr->prev;
         else curr->next->prev = curr->prev;
         if (curr==head) head=curr->next;
         else (curr->prev->next = curr->next);
         temp=curr;
         curr = curr->next;
         free(temp);
      } while (curr != NULL);
      free(reservoir);
   }
   
   if (block_ptr != NULL)
   {
      if (block_left != MPIR_SMALL_BLOCK_SIZE - sizeof(ulong) - sizeof(char *))
      {
         printf("Warning: MPIR small stack memory allocator detected mismatched alloc/release\n");
         while (block_left != MPIR_SMALL_BLOCK_SIZE - sizeof(ulong) - sizeof(char *)) mpir_stack_release_small();
      }  
      
      block_ptr -= (sizeof(ulong) + sizeof(char *));
      mpir_free(block_ptr);           
   } 
}

#endif

/*-----------------------------------------------------------------------------------------------*/

void mpir_memory_failure(void)
{
   printf("Error: unable to alloc/realloc memory\n");
   abort();
}
  
#if DEBUG

void* mpir_alloc(ulong bytes)
{
#if DEBUG_PRINT
   printf("Allocating %ld bytes on heap\n", bytes);
#endif
   char * buf = malloc(OVERRUN+bytes+sizeof(ulong));
   if (buf == NULL)
      mpir_memory_failure();
   ulong * ptr = (ulong *) buf;
   ptr[0] = bytes;
   for (ulong i = 0; i < 100; i++)
      buf[bytes+i+sizeof(ulong)] = CHECK_CHAR;
   return (void*) (buf+sizeof(ulong));
}

void* mpir_realloc(void * block_void, ulong bytes)
{
   char * block = (char *) block_void;
   block-=sizeof(ulong);
   ulong * ptr = (ulong *) block;
   ulong old_bytes = ptr[0];
   for (ulong i = 0; i < 100; i++)
      if (block[old_bytes+i+sizeof(ulong)] != CHECK_CHAR) 
      {
         printf("Error: Block overrun detected by heap memory (re)allocator!!\n");
         abort();
      }
#if DEBUG_PRINT
   printf("Reallocing from %ld to %ld bytes on heap\n", old_bytes, bytes);
#endif
   char * buf = realloc(block, bytes+OVERRUN+sizeof(ulong));
   if (buf == NULL)
      mpir_memory_failure();
   ptr = (ulong *) buf;
   ptr[0] = bytes;
   for (ulong i = 0; i < 100; i++)
      buf[bytes+i+sizeof(ulong)] = CHECK_CHAR;
   return (void*) (buf+sizeof(ulong));
}

void mpir_free(void * block_void)
{
   char * block = (char *) block_void;
   block-=sizeof(ulong);
   ulong * ptr = (ulong *) block;
   ulong bytes = ptr[0];
#if DEBUG_PRINT
   printf("Releasing %ld bytes from heap\n", bytes);
#endif
   for (ulong i = 0; i < 100; i++)
      if (block[bytes+i+sizeof(ulong)] != CHECK_CHAR) 
      {
         printf("Error: Block overrun detected by heap memory allocator!!\n");
         abort();
      }
   free(block);
}

#else

void* mpir_alloc(ulong bytes)
{
   void* buf = malloc(bytes);
   if (buf == NULL)
      mpir_memory_failure();
   return buf;
}

void* mpir_realloc(void* block, ulong bytes)
{
   void* buf = realloc(block, bytes);
   if (buf == NULL)
      mpir_memory_failure();
   return buf;
}

void mpir_free(void* block)
{
   free(block);
}

#endif

void* mpir_aligned_alloc(ulong bytes)
{
   void* block = malloc(sizeof(void *) + sizeof(ulong) + MPIR_ALIGN + bytes - 1);
   void ** ptr = (void **) block;
   ptr++; 
   ulong * ptr_u = (ulong *) ptr;
   ptr_u++;
   char* ptr2 = (char *) ptr_u;
   ptr2 += ((MPIR_ALIGN - ((ulong)ptr2 & (MPIR_ALIGN - 1UL))) & (MPIR_ALIGN - 1UL)); 
   ptr = (void **) ptr2;
   ptr--; 
   ptr[0] = block;
   ulong * ptr3 = (ulong *) ptr;
   ptr3--;
   ptr3[0] = bytes;
   return (void*) ptr2;  
}

void mpir_aligned_free(void * block)
{
   void ** ptr = (void **) block;
   ptr--;
   free(ptr[0]); 
}

void* mpir_aligned_realloc(void * block, ulong bytes)
{
   void ** ptr = (void **) block;
   ptr--;
   ulong * ptr2 = (ulong *) ptr;
   ptr2--;
   void * new_block = mpir_aligned_alloc(bytes);
   
   ulong * old_b = (ulong *) block;
   ulong * new_b = (ulong *) new_block;
   ulong i;
   for (i = 0; i < (ptr2[0]>>MPIR_LG_BYTES); i++)
   {
      new_b[i] = old_b[i];
   }
   char * old_e = (char *) (old_b + i); 
   char * new_e = (char *) (new_b + i); 
   for (ulong j = 0; j < (ptr2[0] & (MPIR_BYTES-1L)); j++)
   {
      new_e[i] = old_e[i];      
   }

   free(ptr[0]);

   return new_block;  
}
