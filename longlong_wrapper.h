/*============================================================================

    Copyright (C) 2007, William Hart, David Harvey

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

===============================================================================*/
#ifndef LONGLONGWRAP_H
#define LONGLONGWRAP_H

#ifdef __cplusplus
 extern "C" {
#endif
 
#include <stdint.h>

// todo: I think perhaps UDWtype is not quite right. It needs to be
// twice the length. But how to do this on a 64-bit machine?
#if FLINT_BITS == 32
#define HAVE_HOST_CPU_pentiummmx 1 
#endif

#define UWtype mp_limb_t
#define UHWtype mp_limb_t
#define UDWtype mp_limb_t 
#define W_TYPE_SIZE FLINT_BITS
#define SItype int32_t
#define USItype uint32_t
#define DItype int64_t
#define UDItype uint64_t

#define LONGLONG_STANDALONE

// todo: longlong.h requires there to be an ASSERT macro around.
// For now I'm killing it, but we should hook this up to our own assertion
// code.

#define ASSERT(condition)

#ifdef __cplusplus
 }
#endif
 
#endif
