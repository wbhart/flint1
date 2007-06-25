#ifndef LONGLONGWRAP_H
#define LONGLONGWRAP_H

#include <stdint.h>

// todo: I think perhaps UDWtype is not quite right. It needs to be
// twice the length. But how to do this on a 64-bit machine?

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

#endif
