/****************************************************************************

memory-manager.h: FLINT-wide memory management

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef FLINT_MANAGER_H
#define FLINT_MANAGER_H

void* flint_stack_alloc(unsigned long length);

void* flint_stack_alloc_bytes(unsigned long bytes);

void* flint_stack_alloc_small(unsigned long length);

void flint_stack_release();

void flint_stack_release_small();

void* flint_heap_alloc(unsigned long limbs);

void* flint_heap_alloc_bytes(unsigned long bytes);

void* flint_heap_realloc(void* block, unsigned long limbs);

void* flint_heap_realloc_bytes(void* block, unsigned long limbs);

void flint_heap_free(void* block);

#endif
