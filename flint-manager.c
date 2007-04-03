#include <stdlib.h>
#include <gmp.h>
#include "flint-manager.h"

//////////////////////////////////////////////////////////
/*
Temporary baby implementations of memory management functions
*/


/*
This function gets called if malloc fails (i.e. returns NULL).
For now we just die via abort(); if we change this to do nothing, then
anyone who calls flint_malloc() etc will get a NULL pointer returned, and
perhaps get a chance to clean up before completely dying...
*/
void flint_memory_failure()
{
   abort();
}

void* flint_malloc_limbs(unsigned long limbs)
{
   void* buf = malloc(limbs * sizeof(mp_limb_t));
   if (!buf)
      flint_memory_failure();
   return buf;
}

void* flint_malloc(unsigned long bytes)
{
   void* buf = malloc(bytes);
   if (!buf)
      flint_memory_failure();
   return buf;
}

void* flint_realloc_limbs(void* block, unsigned long limbs)
{
   void* buf = realloc(block, limbs * sizeof(mp_limb_t));
   if (!buf)
      flint_memory_failure();
   return buf;
}

void* flint_realloc(void* block, unsigned long bytes)
{
   void* buf = realloc(block, bytes);
   if (!buf)
      flint_memory_failure();
   return buf;
}

void flint_free(void* block)
{
   free(block);
}

//////////////////////////////////////////////////////////
