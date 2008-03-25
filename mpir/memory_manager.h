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

memory-manager.h: MPIR-wide memory management

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#ifndef MPIR_MANAGER_H
#define MPIR_MANAGER_H

#ifdef __cplusplus
 extern "C" {
#endif

#include "mpir.h"
 
void* mpir_stack_alloc(ulong bytes);

void* mpir_stack_alloc_small(ulong bytes);

void mpir_stack_release();

void mpir_stack_release_small();

void mpir_stack_cleanup();

void* mpir_alloc(ulong bytes);

void* mpir_aligned_alloc(ulong bytes);

void* mpir_realloc(void* block, ulong bytes);

void* mpir_aligned_realloc(void* block, ulong bytes);

void mpir_free(void* block);

void mpir_aligned_free(void* block);

#ifdef __cplusplus
 }
#endif
 
#endif
