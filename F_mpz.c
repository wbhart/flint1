/*============================================================================

    F_mpz.c: The FLINT integer format (FLINT 2.0)

    Copyright (C) 2008, William Hart 

	 This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================*/

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <gmp.h>

#include "flint.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz.h"
#include "mpz_extras.h"

/*===============================================================================

	mpz_t memory management

================================================================================*/

// The array of mpz's used by the F_mpz type
__mpz_struct * F_mpz_arr;

// Total number of mpz's initialised in F_mpz_arr;
ulong F_mpz_allocated;

// An array of indices of mpz's which are not being used presently. These are stored with the second
// most significant bit set so that they are a valid index as per the F_mpz_t type below.
F_mpz * F_mpz_unused_arr;

// The number of mpz's not being used presently
ulong F_mpz_num_unused = 0;

F_mpz _F_mpz_new_mpz(void)
{
	if (!F_mpz_num_unused) // time to allocate MPZ_BLOCK more mpz_t's
	{
	   if (F_mpz_allocated) // realloc mpz_t's and unused array
		{
			F_mpz_arr = (__mpz_struct *) flint_heap_realloc_bytes(F_mpz_arr, (F_mpz_allocated + MPZ_BLOCK)*sizeof(__mpz_struct));
			F_mpz_unused_arr = (F_mpz *) flint_heap_realloc_bytes(F_mpz_unused_arr, (F_mpz_allocated + MPZ_BLOCK)*sizeof(F_mpz));
		} else // first time alloc of mpz_t's and unused array
		{
			F_mpz_arr = (__mpz_struct*) flint_heap_alloc_bytes(MPZ_BLOCK*sizeof(__mpz_struct));	
			F_mpz_unused_arr = (F_mpz *) flint_heap_alloc_bytes(MPZ_BLOCK*sizeof(F_mpz));
		}
		
		// initialise the new mpz_t's and unused array
		for (ulong i = 0; i < MPZ_BLOCK; i++)
		{
			mpz_init(F_mpz_arr + F_mpz_allocated + i);
			F_mpz_unused_arr[F_mpz_num_unused] = OFF_TO_COEFF(F_mpz_allocated + i);
         F_mpz_num_unused++;
		}

		F_mpz_num_unused--;
		F_mpz_allocated += MPZ_BLOCK;

		return F_mpz_unused_arr[F_mpz_num_unused];
	} else // unused mpz's are available
	{
		F_mpz_num_unused--;
		return F_mpz_unused_arr[F_mpz_num_unused];
	}
}

void _F_mpz_clear_mpz(F_mpz f)
{
   F_mpz_unused_arr[F_mpz_num_unused] = f;
   F_mpz_num_unused++;	
}

void _F_mpz_cleanup(void)
{
	for (ulong i = 0; i < F_mpz_num_unused; i++)
	{
		mpz_clear(F_mpz_arr + COEFF_TO_OFF(F_mpz_unused_arr[i]));
	}
	
	if (F_mpz_num_unused) free(F_mpz_unused_arr);
	if (F_mpz_allocated) free(F_mpz_arr);
}

/*===============================================================================

	Promotion/Demotion

================================================================================*/

__mpz_struct * _F_mpz_promote(F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) *f = _F_mpz_new_mpz(); // f is small so promote it first
	// if f is large already, just return the pointer
      
   return F_mpz_arr + COEFF_TO_OFF(*f);
}

__mpz_struct * _F_mpz_promote_val(F_mpz_t f)
{
   F_mpz c = *f;
	if (!COEFF_IS_MPZ(c)) // f is small so promote it
	{
	   *f = _F_mpz_new_mpz();
	   __mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(*f);
		mpz_set_si(mpz_ptr, c);
		return mpz_ptr;
	} else // f is large already, just return the pointer
      return F_mpz_arr + COEFF_TO_OFF(*f);
}

void _F_mpz_demote_val(F_mpz_t f)
{
   __mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(*f);

	long size = mpz_ptr->_mp_size;
	
	if (size == 0L) // value is zero
	{
		_F_mpz_clear_mpz(*f);
		*f = 0;
	} else if (size == 1L) // value is positive and 1 limb
	{
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= (ulong) COEFF_MAX) 
		{
			_F_mpz_clear_mpz(*f);
			*f = (F_mpz) uval;
		}
	} else if (size == -1L) // value is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(mpz_ptr);
		if (uval <= (ulong) COEFF_MAX) 
		{
			_F_mpz_clear_mpz(*f);
			*f = (F_mpz) -uval;
		}
	}
	// don't do anything if value has to be multi precision
}

/*===============================================================================

	F_mpz_t memory management

================================================================================*/

void F_mpz_init2(F_mpz_t f, ulong limbs)
{
	if (limbs) 
	{
		*f = _F_mpz_new_mpz();
		_mpz_realloc(F_mpz_arr + COEFF_TO_OFF(*f), limbs);
		return;
	} else 
	{
		*f = 0;
		return;
	}
}

/* ==============================================================================

   Random generation

===============================================================================*/

/*
  These are not serious random generators, they are just here for testing 
  purposes at this stage

  We require bits to be non-zero
*/

#define NORM(xxx, coeffxxx) \
do { \
   if ((long) xxx->_mp_size < 0L) \
   { \
      while ((xxx->_mp_size) && (!(coeffxxx)[-xxx->_mp_size - 1])) xxx->_mp_size++; \
   } else if ((long) xxx->_mp_size > 0L) \
   { \
      while ((xxx->_mp_size) && (!(coeffxxx)[xxx->_mp_size - 1])) xxx->_mp_size--; \
   } \
} while (0);

void F_mpz_random(F_mpz_t f, const ulong bits)
{
	if (bits <= FLINT_BITS - 2)
   {
      ulong temp;
      mpn_random(&temp, 1L);
      ulong mask = ((1L<<bits)-1L);
      *f = temp & mask;
      return;
   }
   
	ulong limbs = ((bits-1)>>FLINT_LG_BITS_PER_LIMB)+1;
   ulong rem = (bits & (FLINT_BITS - 1));
   
   __mpz_struct * mpz_ptr = _F_mpz_promote(f);
   mpz_realloc2(mpz_ptr, bits);
	
	mp_limb_t * fp = mpz_ptr->_mp_d;
   mpz_ptr->_mp_size = limbs;
   mpn_random(fp, limbs);
   if (rem)
   {
      ulong mask = ((1L<<rem)-1L);
      fp[limbs-1] &= mask;
   }
   NORM(mpz_ptr, fp);
	_F_mpz_demote_val(f);
}

void F_mpz_randomm(F_mpz_t f, const mpz_t in)
{
   if (mpz_size(in) > 1) 
   {
      __mpz_struct * mpz_ptr = _F_mpz_promote(f);
		mpz_urandomm(mpz_ptr, F_mpz_state, in);
		_F_mpz_demote_val(f);
   } else
   {
      ulong val = mpz_get_ui(in);
		ulong rnd = (val == 0 ? 0L : z_randint(val));
		F_mpz_set_ui(f, rnd);
   }
}

/*===============================================================================

	Get/set

================================================================================*/

void F_mpz_set_si(F_mpz_t f, const long val)
{
   if (FLINT_ABS(val) > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(f);
		mpz_set_si(mpz_coeff, val);
	} else 
	{
		_F_mpz_demote(f);
		*f = val; // val is small
	}
}

void F_mpz_set_ui(F_mpz_t f, const ulong val)
{
   if (val > COEFF_MAX) // val is large
	{
		__mpz_struct * mpz_coeff = _F_mpz_promote(f);
		mpz_set_ui(mpz_coeff, val);
	} else 
	{
		_F_mpz_demote(f);
		*f = val; // val is small
	}
}

long F_mpz_get_si(const F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) return *f; // value is small
	return mpz_get_si(F_mpz_arr + COEFF_TO_OFF(*f)); // value is large
}

ulong F_mpz_get_ui(const F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) 
	{ 
		if (*f < 0L) return -*f; // value is small
		else return *f;
	}
	return mpz_get_ui(F_mpz_arr + COEFF_TO_OFF(*f)); // value is large
}

void F_mpz_get_mpz(mpz_t x, const F_mpz_t f)
{
	if (!COEFF_IS_MPZ(*f)) mpz_set_si(x, *f); // set x to small value
	else mpz_set(x, F_mpz_arr + COEFF_TO_OFF(*f)); // set x to large value
}

extern double __gmpn_get_d(mp_limb_t *, size_t, size_t, long);

double F_mpz_get_d_2exp(long * exp, const F_mpz_t f)
{
   F_mpz d = *f;

	if (!COEFF_IS_MPZ(d))
   {
      if (d == 0L) 
      {
         (*exp) = 0L;
         return 0.0;
      }
      ulong d_abs = FLINT_ABS(d);
      (*exp) = FLINT_BIT_COUNT(d_abs);
      if (d < 0L) return __gmpn_get_d(&d_abs, 1L, -1L, -*exp);
      else return __gmpn_get_d(&d, 1L, 1L, -*exp);
   } else 
	   return mpz_get_d_2exp(exp, F_mpz_arr + COEFF_TO_OFF(d));
}

void F_mpz_set_mpz(F_mpz_t f, const mpz_t x)
{
   long size = (long) x->_mp_size;
	
	if (size == 0L) // x is zero
	{
		F_mpz_zero(f);
	} else if (size == 1L) // x is positive and 1 limb
	{
	   F_mpz_set_ui(f, mpz_get_ui(x));
	} else if (size == -1L) // x is negative and 1 limb
   {
	   ulong uval = mpz_get_ui(x);
		if (uval <= COEFF_MAX) // x is small
		{
		   _F_mpz_demote(f);
		   *f = -uval;
		} else // x is large but one limb
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			mpz_set_ui(mpz_ptr, uval);
			mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // x is more than one limb
	{
		__mpz_struct * mpz_ptr = _F_mpz_promote(f);
		mpz_set(mpz_ptr, x);
	}			
}

void F_mpz_set_limbs(F_mpz_t f, const mp_limb_t * x, const ulong limbs)
{
	if (limbs == 0L) // x is zero
	{
		F_mpz_zero(f);
	} else if (limbs == 1L) // x is 1 limb
	{
	   F_mpz_set_ui(f, x[0]);
	} else // x is more than one limb
	{
		__mpz_struct * mpz_ptr = _F_mpz_promote(f);
		// read limbs, least significant first, native endianness, no nails
		mpz_import(mpz_ptr, limbs, -1, sizeof(mp_limb_t), 0, 0, x);
	}			
}

ulong F_mpz_get_limbs(mp_limb_t * x, const F_mpz_t f)
{
	ulong limbs = F_mpz_size(f);
	
	if (limbs == 0L) return 0; // f is zero, no limbs to get

	if (limbs == 1L) // f is 1 limb
	{
	   x[0] = F_mpz_get_ui(f);
	} else // f is more than one limb
	{
		__mpz_struct * mpz_ptr = F_mpz_ptr_mpz(*f);
		// discard count, read least significant first, native endianness, no nails
		mpz_export(x, NULL, -1, sizeof(mp_limb_t), 0, 0, mpz_ptr);
	}	

	return limbs;
}

void F_mpz_set(F_mpz_t f, const F_mpz_t g)
{
	if (f == g) return; // aliased inputs
	
	if (!COEFF_IS_MPZ(*g)) // g is small
	{
		_F_mpz_demote(f);
		*f = *g;
	} else // g is large
	{
	   __mpz_struct * mpz_ptr = _F_mpz_promote(f);
		mpz_set(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*g));
	}
}

void F_mpz_swap(F_mpz_t f, F_mpz_t g)
{
	if (f == g) return; // swapping not required
	
	if (!COEFF_IS_MPZ(*f))
	{
		if (!COEFF_IS_MPZ(*g)) // both values are small
		{
	      F_mpz t = *f;
			*f = *g;
	      *g = t;
		} else // f is small, g is large
		{
			F_mpz t = *f;
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			mpz_set(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*g));
			_F_mpz_demote(g);
         *g = t;
		}
	} else
	{
      if (!COEFF_IS_MPZ(*g)) // g is small, f is large
		{
         F_mpz t = *g;
			__mpz_struct * mpz_ptr = _F_mpz_promote(g);
			mpz_set(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*f));
			_F_mpz_demote(f);
         *f = t;
		} else // both values are large
		{
			mpz_swap(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g));
		}
	}
}

/*===============================================================================

	Comparison

================================================================================*/

int F_mpz_equal(const F_mpz_t f, const F_mpz_t g)
{
	if (f == g) return 1; // aliased inputs
	
	if (!COEFF_IS_MPZ(*f)) return (*f == *g); // if f is large it can't be equal to g
	else if (!COEFF_IS_MPZ(*g)) return 0; // f is large, so if g isn't....
	else return (mpz_cmp(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)) == 0); 
}

int F_mpz_cmpabs(const F_mpz_t f, const F_mpz_t g)
{
	if (f == g) return 0; // aliased inputs
	
	if (!COEFF_IS_MPZ(*f)) 
	{
		if (!COEFF_IS_MPZ(*g)) 
		{
         ulong uf = FLINT_ABS(*f);
         ulong ug = FLINT_ABS(*g);
         if (uf < ug) return -1;
			else return (uf > ug);
		} else return -1;
	} else if (!COEFF_IS_MPZ(*g)) return 1; // f is large, so if g isn't....
	else return mpz_cmpabs(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)); 
}

int F_mpz_cmp(const F_mpz_t f, const F_mpz_t g)
{
	if (f == g) return 0; // aliased inputs
	
	if (!COEFF_IS_MPZ(*f)) 
	{
		if (!COEFF_IS_MPZ(*g)) // both coeffs small
		{
         if (*f < *g) return -1;
			else return (*f > *g);
		} else // f is small, g is large 
		{
			if (mpz_sgn(F_mpz_arr + COEFF_TO_OFF(*g)) < 0) return 1; // g is a large negative 
			else return -1; // g is a large positive
		}
	} else if (!COEFF_IS_MPZ(*g)) 
	{
		if (mpz_sgn(F_mpz_arr + COEFF_TO_OFF(*f)) < 0) return -1; // f is large negative
		else return 1; // f is large positive
	} else // both f and g are large 
		return mpz_cmp(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)); 
}

/*===============================================================================

	Properties

================================================================================*/

ulong F_mpz_size(const F_mpz_t f)
{
	F_mpz d = *f;

	if (d == 0) return 0;
   if (!COEFF_IS_MPZ(d)) // c1 is small
	{
		return 1;
	}

	return mpz_size(F_mpz_arr + COEFF_TO_OFF(d));
}

int F_mpz_sgn(const F_mpz_t f)
{
	F_mpz d = *f;

	if (d == 0) return 0;
   if (!COEFF_IS_MPZ(d)) // c1 is small
	{
		if (d > 0L) return 1;
		else return -1;
	}

	return mpz_sgn(F_mpz_arr + COEFF_TO_OFF(d));
}

ulong F_mpz_bits(const F_mpz_t f)
{
	F_mpz d = *f;

	if (!COEFF_IS_MPZ(d)) // c1 is small
	{
		return FLINT_BIT_COUNT(FLINT_ABS(d));
	}

	return mpz_sizeinbase(F_mpz_arr + COEFF_TO_OFF(d), 2);
}

__mpz_struct * F_mpz_ptr_mpz(F_mpz f)
{
	return F_mpz_arr + COEFF_TO_OFF(f);
}

/*===============================================================================

	Input/output

================================================================================*/

void F_mpz_read(F_mpz_t f)
{
	mpz_t temp;
	mpz_init(temp);

	mpz_inp_str(temp, stdin, 10);
	F_mpz_set_mpz(f, temp);

	mpz_clear(temp);
}

/*===============================================================================

	Arithmetic

================================================================================*/

void F_mpz_neg(F_mpz_t f1, const F_mpz_t f2)
{
   if (!COEFF_IS_MPZ(*f2)) // coeff is small
	{
		F_mpz t = -*f2; // Need to save value in case of aliasing
		_F_mpz_demote(f1);
		*f1 = t;
	} else // coeff is large
	{
	   // No need to retain value in promotion, as if aliased, both already large
		__mpz_struct * mpz_ptr = _F_mpz_promote(f1);
		mpz_neg(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*f2));
	}
}

void F_mpz_add(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // both inputs are small
		{
			F_mpz_set_si(f, c1 + c2);
		} else // g is small, h is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // g is saved and h is large
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			if (c1 < 0L) mpz_sub_ui(mpz3, mpz2, -c1);	
		   else mpz_add_ui(mpz3, mpz2, c1);
			_F_mpz_demote_val(f); // may have cancelled
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // g is large, h is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // h is saved and g is large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			if (c2 < 0L) mpz_sub_ui(mpz3, mpz1, -c2);	
			else mpz_add_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(f); // may have cancelled
		} else // g and h are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // aliasing means f is already large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			mpz_add(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(f); // may have cancelled
		}
	}
}

void F_mpz_add_mpz(F_mpz_t f, const F_mpz_t g, mpz_t h)
{
	F_mpz c1 = *g;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   __mpz_struct * mpz3 = _F_mpz_promote(f); // g is saved and h is large
		if (c1 < 0L) mpz_sub_ui(mpz3, h, -c1);	
		else mpz_add_ui(mpz3, h, c1);
		_F_mpz_demote_val(f); // may have cancelled
	} else // g is large
	{
		__mpz_struct * mpz3 = _F_mpz_promote(f); // aliasing means f is already large
		__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
		mpz_add(mpz3, mpz1, h);
		_F_mpz_demote_val(f); // may have cancelled
	}
}

void F_mpz_sub(F_mpz_t f, const F_mpz_t g, F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // both inputs are small
		{
			F_mpz_set_si(f, c1 - c2);
		} else // g is small, h is large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // g is saved and h is large
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			if (c1 < 0L) 
			{
				mpz_add_ui(mpz3, mpz2, -c1);
				mpz_neg(mpz3, mpz3);
			} else mpz_ui_sub(mpz3, c1, mpz2);
			_F_mpz_demote_val(f); // may have cancelled
		}
	} else
	{
		if (!COEFF_IS_MPZ(c2)) // g is large, h is small
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // h is saved and g is large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			if (c2 < 0L) mpz_add_ui(mpz3, mpz1, -c2);	
			else mpz_sub_ui(mpz3, mpz1, c2);
			_F_mpz_demote_val(f); // may have cancelled
		} else // g and h are large
		{
         __mpz_struct * mpz3 = _F_mpz_promote(f); // aliasing means f is already large
			__mpz_struct * mpz1 = F_mpz_arr + COEFF_TO_OFF(c1);
			__mpz_struct * mpz2 = F_mpz_arr + COEFF_TO_OFF(c2);
			mpz_sub(mpz3, mpz1, mpz2);
			_F_mpz_demote_val(f); // may have cancelled
		}
	}
}

void F_mpz_mul_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
{
	F_mpz c2 = *g;
	
	if (x == 0)
	{
		F_mpz_zero(f);
		return;
	} else if (!COEFF_IS_MPZ(c2)) // coeff2 is small
	{
		mp_limb_t prod[2];
		ulong uc2 = FLINT_ABS(c2);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], uc2, x); 
		if (!prod[1]) // result fits in one limb
		{
			F_mpz_set_ui(f, prod[0]);
			if (c2 < 0L) F_mpz_neg(f, f);
		} else // result takes two limbs
		{
		   __mpz_struct * mpz_ptr = _F_mpz_promote(f);
			// two limbs, least significant first, native endian, no nails, stored in prod
         mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if (c2 < 0L) mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // coeff2 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f); // promote without val as if aliased both are large
      mpz_mul_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c2), x);
	}
}

void F_mpz_mul_si(F_mpz_t f, const F_mpz_t g, const long x)
{
	F_mpz c2 = *g;

	if (x == 0)
	{
		F_mpz_zero(f);
		return;
	} else if (!COEFF_IS_MPZ(c2)) // coeff2 is small
	{
		mp_limb_t prod[2];
		ulong uc2 = FLINT_ABS(c2);
		ulong ux = FLINT_ABS(x);
		
		// unsigned limb by limb multiply (assembly for most CPU's)
		umul_ppmm(prod[1], prod[0], uc2, ux); 
		if (!prod[1]) // result fits in one limb
		{
			F_mpz_set_ui(f, prod[0]);
			if ((c2 ^ x) < 0L) F_mpz_neg(f, f);
		} else // result takes two limbs
		{
		   __mpz_struct * mpz_ptr = _F_mpz_promote(f);
         // two limbs, least significant first, native endian, no nails, stored in prod
			mpz_import(mpz_ptr, 2, -1, sizeof(mp_limb_t), 0, 0, prod);
			if ((c2 ^ x) < 0L) mpz_neg(mpz_ptr, mpz_ptr);
		}
	} else // coeff2 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f); // ok without val as if aliased both are large
      mpz_mul_si(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c2), x);
	}
}

void F_mpz_mul2(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
   
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
		F_mpz_mul_si(f, h, c1);
      return;
	}
	
	F_mpz c2 = *h; // save h in case it is aliased with f
   
	if (!c2) // special case, h = 0 
	{
		F_mpz_zero(f);
		return;
	}

   __mpz_struct * mpz_ptr = _F_mpz_promote(f); // h is saved, g is already large
		
	if (!COEFF_IS_MPZ(c2)) // g is large, h is small
	   mpz_mul_si(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), c2);
   else // c1 and c2 are large
	   F_mpz_mul(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
}

void F_mpz_mul_2exp(F_mpz_t f, const F_mpz_t g, const ulong exp)
{
	F_mpz d = *g;

	if (!COEFF_IS_MPZ(d)) // g is small
	{
		ulong dabs = FLINT_ABS(d);
		ulong bits = FLINT_BIT_COUNT(dabs);
		if (bits + exp <= FLINT_BITS - 2) // result will fit in small
		{
			F_mpz_set_si(f, d<<exp);
		} else // result is large
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(f); // g is saved
         mpz_set_si(mpz_ptr, d); 
	      mpz_mul_2exp(mpz_ptr, mpz_ptr, exp);
		}
	} else // g is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f); // g is already large
      mpz_mul_2exp(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(d), exp);   
	}
}

void F_mpz_div_2exp(F_mpz_t f, const F_mpz_t g, const ulong exp)
{
	F_mpz d = *g;

	if (!COEFF_IS_MPZ(d)) // g is small
	{
		F_mpz_set_si(f, d>>exp);
	} else // g is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f); // g is already large
		mpz_div_2exp(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(d), exp);   
		_F_mpz_demote_val(f); // division may make value small
	}
}

void F_mpz_add_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
{
	F_mpz c = *g;

	if (!COEFF_IS_MPZ(c)) // g is small
	{
      mp_limb_t sum[2];
		if (c >= 0L) // both operands non-negative
		{
			add_ssaaaa(sum[1], sum[0], 0, c, 0, x);
			if (sum[1] == 0) F_mpz_set_ui(f, sum[0]); // result fits in 1 limb
			else // result takes two limbs
			{
				mpz_t temp;
				temp->_mp_d = sum;
				temp->_mp_size = 2; // result is sum of two non-negative numbers and is hence non-negative
				__mpz_struct * mpz_ptr = _F_mpz_promote(f); // g has already been read
				mpz_set(mpz_ptr, temp);
			}
		} else // coeff is negative, x positive
		{
			if (-c > x) F_mpz_set_si(f, x + c); // can't overflow as g is small and x smaller
			else F_mpz_set_ui(f, x + c); // won't be negative and has to be less than x
		}
	} else
	{
		__mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(c);
		__mpz_struct * mpz_ptr2 = _F_mpz_promote(f); // g is already large
		mpz_add_ui(mpz_ptr2, mpz_ptr, x);
		_F_mpz_demote_val(f); // cancellation may have occurred
	}
}

void F_mpz_sub_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
{
	F_mpz c = *g;

	if (!COEFF_IS_MPZ(c)) // coeff is small
	{
      mp_limb_t sum[2];
		if (c < 0L) // g negative, x positive, so difference is negative
		{
			add_ssaaaa(sum[1], sum[0], 0, -c, 0, x);
			if (sum[1] == 0) 
			{   
				F_mpz_set_ui(f, sum[0]); // result fits in 1 limb
				F_mpz_neg(f, f);
			} else // result takes two limbs
			{
				mpz_t temp;
				temp->_mp_d = sum;
				temp->_mp_size = -2; // result is negative number minus negative number, hence negative
				__mpz_struct * mpz_ptr = _F_mpz_promote(f); // g has already been read
				mpz_set(mpz_ptr, temp);
			}
		} else // coeff is non-negative, x non-negative
		{
			if (x < c) F_mpz_set_ui(f, c - x); // won't be negative and is smaller than c
			else 
			{
				F_mpz_set_ui(f, x - c); // positive or zero
				F_mpz_neg(f, f);
			}
		}
	} else
	{
		__mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(c);
		__mpz_struct * mpz_ptr2 = _F_mpz_promote(f); // g is already large
		mpz_sub_ui(mpz_ptr2, mpz_ptr, x);
		_F_mpz_demote_val(f); // cancellation may have occurred
	}
}

void F_mpz_addmul_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
{
	F_mpz c1 = *g;
   if ((x == 0) || (c1 == 0)) return; // product is zero
   
	F_mpz r = *f;
	if (r == 0) 
	{
		F_mpz_mul_ui(f, g, x); // we are adding product to 0
		return;
	}

	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
      mp_limb_t prod[2];
	   ulong uc1 = FLINT_ABS(c1);
      
		umul_ppmm(prod[1], prod[0], uc1, x); // compute product

		if (prod[1] == 0) // product fits in one limb
		{
			if (c1 < 0L) 
			   F_mpz_sub_ui(f, f, prod[0]);
			else F_mpz_add_ui(f, f, prod[0]);
			return;
		} else if ((prod[1] == 1) && (!COEFF_IS_MPZ(r)) && ((r ^ c1) < 0L))
		{
			// only chance at cancellation is if product is one bit past a limb
			// and res is small and opposite sign to this product
			
			ulong ur = FLINT_ABS(r);
			if (ur > prod[0]) // cancellation will occur
			{
				F_mpz_set_ui(f, prod[0] - ur);
				if (r > 0L) F_mpz_neg(f, f);
				return;
			} 
		}
		
		// in all remaining cases res is either big already, or will be big in the end
	   __mpz_struct * mpz_ptr = _F_mpz_promote_val(f); 
		
		mpz_t temp; // set up a temporary, cheap mpz_t to contain prod
	   temp->_mp_d = prod;
	   temp->_mp_size = (c1 < 0L ? -2 : 2);
	   mpz_add(mpz_ptr, mpz_ptr, temp);
		_F_mpz_demote_val(f); // cancellation may have occurred
	} else // c1 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote_val(f);
		
      mpz_addmul_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), x);
		_F_mpz_demote_val(f); // cancellation may have occurred
	}
}

void F_mpz_submul_ui(F_mpz_t f, const F_mpz_t g, const ulong x)
{
	F_mpz c1 = *g;
   if ((x == 0) || (c1 == 0)) return; // product is zero
   
	F_mpz r = *f;
	if (r == 0) 
	{
		F_mpz_mul_ui(f, g, x); // we are subtracting product from 0
		F_mpz_neg(f, f);
		return;
	}

	if (!COEFF_IS_MPZ(c1)) // c1 is small
	{
      mp_limb_t prod[2];
	   ulong uc1 = FLINT_ABS(c1);
      
		umul_ppmm(prod[1], prod[0], uc1, x); // compute product

		if (prod[1] == 0) // product fits in one limb
		{
			if (c1 < 0L) F_mpz_add_ui(f, f, prod[0]);
			else F_mpz_sub_ui(f, f, prod[0]);
			return;
		} else if ((prod[1] == 1) && (!COEFF_IS_MPZ(r)) && ((r ^ c1) >= 0L))
		{
			// only chance at cancellation is if product is one bit past a limb
			// and f is small and same sign as this product
			
			ulong ur = FLINT_ABS(r);
			if (ur > prod[0]) // cancellation will occur
			{
				F_mpz_set_ui(f, prod[0] - ur);
				if (r > 0L) F_mpz_neg(f, f);
				return;
			} 
		}
		
		// in all remaining cases res is either big already, or will be big in the end
	   __mpz_struct * mpz_ptr = _F_mpz_promote_val(f); 
		
		mpz_t temp; // set up a temporary, cheap mpz_t to contain prod
	   temp->_mp_d = prod;
	   temp->_mp_size = (c1 < 0L ? -2 : 2);
	   mpz_sub(mpz_ptr, mpz_ptr, temp);
		_F_mpz_demote_val(f); // cancellation may have occurred
	} else // c1 is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote_val(f);
		
      mpz_submul_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), x);
		_F_mpz_demote_val(f); // cancellation may have occurred
	}
}

void F_mpz_addmul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
		if (c1 < 0L) F_mpz_submul_ui(f, h, -c1);
		else F_mpz_addmul_ui(f, h, c1);
		return;
	} 

	F_mpz c2 = *h;
   
	if (!COEFF_IS_MPZ(c2)) // h is small
	{
		if (c2 < 0L) F_mpz_submul_ui(f, g, -c2);
		else F_mpz_addmul_ui(f, g, c2);
		return;
	} 

	// both g and h are large
   __mpz_struct * mpz_ptr = _F_mpz_promote_val(f);
	
   mpz_addmul(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
	_F_mpz_demote_val(f); // cancellation may have occurred
}

void F_mpz_submul(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	
	if (!COEFF_IS_MPZ(c1)) // g is small
	{
		if (c1 < 0L) F_mpz_addmul_ui(f, h, -c1);
		else F_mpz_submul_ui(f, h, c1);
		return;
	} 

	F_mpz c2 = *h;
   
	if (!COEFF_IS_MPZ(c2)) // h is small
	{
		if (c2 < 0L) F_mpz_addmul_ui(f, g, -c2);
		else F_mpz_submul_ui(f, g, c2);
		return;
	} 

	// both g and h are large
   __mpz_struct * mpz_ptr = _F_mpz_promote_val(f);
	
   mpz_submul(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
	_F_mpz_demote_val(f); // cancellation may have occurred
}

void F_mpz_mod(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
   if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // h is also small
		   F_mpz_set_si(f, c1 % c2);
		else // h is large and g is small
			F_mpz_set_si(f, c1);
	} else // g is large
	{
      if (!COEFF_IS_MPZ(c2)) // h is small
		   F_mpz_set_ui(f, mpz_fdiv_ui(F_mpz_arr + COEFF_TO_OFF(c1), c2));
		else // both are large
		{
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			mpz_mod(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // reduction mod h may result in small value
		}	
	}
}

void F_mpz_divexact(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
   if (!COEFF_IS_MPZ(c1)) // g is small, h must be also or division isn't exact
	{
	   F_mpz_set_si(f, c1 / c2);
	} else // g is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f);
		
		if (!COEFF_IS_MPZ(c2)) // h is small
		{
		   if (c2 > 0) // h > 0
			{
            mpz_divexact_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), c2);
			   _F_mpz_demote_val(f); // division by h may result in small value
			} else
			{
            mpz_divexact_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), -c2);
			   _F_mpz_demote_val(f); // division by h may result in small value
				F_mpz_neg(f, f);
			}
		} else // both are large
		{
			mpz_divexact(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // division by h may result in small value
		}	
	}
}

void F_mpz_cdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
   if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // h is also small
		{
			F_mpz q = c1 / c2; // compute C quotient
			F_mpz r = c1 - c2*q; // compute remainder
			if (r > 0L) q++; // q cannot overflow as remainder implies |c2| != 1
			F_mpz_set_si(f, q);
		} else // h is large and g is small
		{
			if (c1 == 0L) F_mpz_set_ui(f, 0L); // g is zero
			else if (((c1 < 0L) && (F_mpz_sgn(h) < 0)) || 
				      ((c1 > 0L) && (F_mpz_sgn(h) > 0))) // signs are the same
				F_mpz_set_ui(f, 1); // quotient is positive, round up to one
			else F_mpz_zero(f); // quotient is negative, round up to zero
            
		}
	} else // g is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f);
		
		if (!COEFF_IS_MPZ(c2)) // h is small
		{
		   if (c2 > 0) // h > 0
			{
            mpz_cdiv_q_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), c2);
			   _F_mpz_demote_val(f); // division by h may result in small value
			} else
			{
            mpz_fdiv_q_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), -c2);
			   _F_mpz_demote_val(f); // division by h may result in small value
				F_mpz_neg(f, f);
			}
		} else // both are large
		{
			mpz_cdiv_q(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // division by h may result in small value
		}	
	}
}

void F_mpz_fdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	
   if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // h is also small
		{
			F_mpz q = c1 / c2; // compute C quotient
			F_mpz r = c1 - c2*q; // compute remainder
			if (r < 0L) q--; // q cannot overflow as remainder implies |c2| != 1
			F_mpz_set_si(f, q);
		} else // h is large and g is small
		{
			if (c1 == 0L) F_mpz_set_ui(f, 0L); // g is zero
			else if (((c1 < 0L) && (F_mpz_sgn(h) < 0)) || 
				      ((c1 > 0L) && (F_mpz_sgn(h) > 0))) // signs are the same
				F_mpz_zero(f); // quotient is positive, round down to zero
			else F_mpz_set_si(f, -1L); // quotient is negative, round down to minus one
            
		}
	} else // g is large
	{
      __mpz_struct * mpz_ptr = _F_mpz_promote(f);
		
		if (!COEFF_IS_MPZ(c2)) // h is small
		{
		   if (c2 > 0) // h > 0
			{
            mpz_fdiv_q_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), c2);
			   _F_mpz_demote_val(f); // division by h may result in small value
			} else
			{
            mpz_cdiv_q_ui(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), -c2);
			   _F_mpz_demote_val(f); // division by h may result in small value
				F_mpz_neg(f, f);
			}
		} else // both are large
		{
			mpz_fdiv_q(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // division by h may result in small value
		}	
	}
}

void F_mpz_rdiv_q(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz_t h2;  
	F_mpz_init(h2);

	F_mpz_div_2exp(h2, h, 1); // h2 = h >> 1
	F_mpz_add(h2, h2, g); // h2 = g + (h >> 1) 
   F_mpz_fdiv_q(f, h2, h); // f = floor((g + (h >> 1) / h)
   
	F_mpz_clear(h2);
}
