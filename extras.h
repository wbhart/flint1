#ifndef EXTRAS_H
#define EXTRAS_H

#include <gmp.h>
#include "flint.h"
#include "longlong_wrapper.h"

static inline unsigned long r_shift(unsigned long in, unsigned long shift)
{
   if (shift == FLINT_BITS_PER_LIMB) return 0L;
   return (in>>shift);
}

static inline unsigned long l_shift(unsigned long in, unsigned long shift)
{
   if (shift == FLINT_BITS_PER_LIMB) return 0L;
   return (in<<shift);
}

#define invert_limb(invxl,xl)                   \
  do {                                          \
    mp_limb_t dummy;                            \
    udiv_qrnnd (invxl, dummy, ~(xl), ~(0L), xl);  \
  } while (0)

#define LIMB_HIGHBIT_TO_MASK(n)                                 \
  (((mp_limb_signed_t) -1 >> 1) < 0                             \
   ? (mp_limb_signed_t) (n) >> (FLINT_BITS_PER_LIMB - 1)              \
   : (n) & (1L<<(FLINT_BITS_PER_LIMB-1)) ? (~ (mp_limb_t) 0L) : (0L))

#define udiv_qrnnd_preinv(q, r, nh, nl, d, di)				\
  do {									\
    mp_limb_t _n2, _n10, _nmask, _nadj, _q1;				\
    mp_limb_t _xh, _xl;							\
    _n2 = (nh);								\
    _n10 = (nl);							\
    _nmask = LIMB_HIGHBIT_TO_MASK (_n10);				\
    _nadj = _n10 + (_nmask & (d));					\
    umul_ppmm (_xh, _xl, di, _n2 - _nmask);				\
    add_ssaaaa (_xh, _xl, _xh, _xl, _n2, _nadj);			\
    _q1 = ~_xh;								\
    umul_ppmm (_xh, _xl, _q1, d);					\
    add_ssaaaa (_xh, _xl, _xh, _xl, nh, nl);				\
    _xh -= (d);					/* xh = 0 or -1 */	\
    (r) = _xl + ((d) & _xh);						\
    (q) = _xh - _q1;							\
  } while (0)

#endif
