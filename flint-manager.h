#ifndef FLINT_MANAGER_H
#define FLINT_MANAGER_H
//////////////////////////////////////////////////////////
/*
flint_malloc_limbs and flint_malloc use the same memory manager, they just
measure the size in different units

At first they would just be a wrapper for malloc/free, but eventually we
probably want an extra layer in between flint_malloc and malloc, which finds
free memory faster, possibly at the expense of greater fragmentation
*/

void* flint_malloc_limbs(unsigned long limbs);
void* flint_malloc(unsigned long bytes);
void* flint_realloc_limbs(void* block, unsigned long limbs);
void* flint_realloc(void* block, unsigned long bytes);
void flint_free(void* block);
//////////////////////////////////////////////////////////

#endif
