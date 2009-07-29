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
#include <pthread.h>

#include "flint.h"
#include "mpn_extras.h"
#include "longlong_wrapper.h"
#include "longlong.h"
#include "memory-manager.h"
#include "long_extras.h"
#include "F_mpz.h"
#include "mpz_extras.h"
#include "zn_poly/src/zn_poly.h"

/*===============================================================================

	mpz_t memory management

================================================================================*/

// The array of mpz's used by the F_mpz type
__mpz_struct * F_mpz_arr;

// Total number of mpz's initialised in F_mpz_arr;
ulong F_mpz_allocated;

// An array of indices of mpz's which are not being used presently. These are stored with the second
// most significant bit set so that they are a valid index as per the F_mpz_t type below.
__thread F_mpz * F_mpz_unused_arr;

// The number of mpz's not being used presently
__thread ulong F_mpz_num_unused = 0;

#define MAX_SEM 16 

int sem;
pthread_mutex_t sem_mutex;
pthread_mutex_t sem_acquire_mutex;

void semaphore_init(void)
{
   pthread_mutex_init(&sem_mutex, NULL);
   pthread_mutex_init(&sem_acquire_mutex, NULL);
   pthread_mutex_init(&F_mpz_random_mutex, NULL);
	sem = 0;
}

void semaphore_up(void)
{
	int done = 0;
	while (!done)
	{
      pthread_mutex_lock(&sem_mutex);
	   if (sem < MAX_SEM) 
		{
			sem++;
			done = 1;
		}
		pthread_mutex_unlock(&sem_mutex);   
	}
}

void semaphore_down(void)
{
	pthread_mutex_lock(&sem_mutex);
	sem--;
	pthread_mutex_unlock(&sem_mutex);   
}

void semaphore_acquire(void)
{
	pthread_mutex_lock(&sem_acquire_mutex);
	for (int i = 0; i < MAX_SEM; i++)
	{
		semaphore_up();
	}
}

void semaphore_release(void)
{
	pthread_mutex_lock(&sem_mutex);
	sem = 1; // our thread still has one semaphore unit to release
	pthread_mutex_unlock(&sem_mutex);  
	pthread_mutex_unlock(&sem_acquire_mutex);
}

F_mpz _F_mpz_new_mpz(void)
{
	if (!F_mpz_num_unused) // time to allocate MPZ_BLOCK more mpz_t's
	{
	   semaphore_down();
	   semaphore_acquire();
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
	   semaphore_release();
      
		F_mpz_num_unused--;
		F_mpz_allocated += MPZ_BLOCK;
	} else // unused mpz's are available
	{
		F_mpz_num_unused--;
	}
	F_mpz ret = F_mpz_unused_arr[F_mpz_num_unused];
	
	return ret;
}

void _F_mpz_clear_mpz(F_mpz f)
{
   F_mpz_unused_arr[F_mpz_num_unused] = f;
   F_mpz_num_unused++;	
}

void _F_mpz_cleanup(void)
{
	if (sem != 0) printf("Warning: unbalanced semaphores, sem = %d\n", sem);
	for (ulong i = 0; i < F_mpz_num_unused; i++)
	{
		mpz_clear(F_mpz_arr + COEFF_TO_OFF(F_mpz_unused_arr[i]));
	}
	
	if (F_mpz_num_unused) free(F_mpz_unused_arr);
	if (F_mpz_allocated) free(F_mpz_arr);
}

void _F_mpz_cleanup2(void)
{
	for (ulong i = 0; i < F_mpz_num_unused; i++)
	{
		mpz_clear(F_mpz_arr + COEFF_TO_OFF(F_mpz_unused_arr[i]));
	}
	F_mpz_num_unused = 0;
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
   pthread_mutex_lock(&F_mpz_random_mutex);
	mpn_random(fp, limbs);
   pthread_mutex_unlock(&F_mpz_random_mutex);
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
		pthread_mutex_lock(&F_mpz_random_mutex);
	   mpz_urandomm(mpz_ptr, F_mpz_state, in);
		pthread_mutex_unlock(&F_mpz_random_mutex);
	   _F_mpz_demote_val(f);
		
   } else
   {
      ulong val = mpz_get_ui(in);
		pthread_mutex_lock(&F_mpz_random_mutex);
	   ulong rnd = (val == 0 ? 0L : z_randint(val));
		pthread_mutex_unlock(&F_mpz_random_mutex);
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
	
	long ret = mpz_get_si(F_mpz_arr + COEFF_TO_OFF(*f)); // value is large
	
	return ret;
}

ulong F_mpz_get_ui(const F_mpz_t f)
{
   if (!COEFF_IS_MPZ(*f)) 
	{ 
		if (*f < 0L) return -*f; // value is small
		else return *f;
	}
	
	ulong ret = mpz_get_ui(F_mpz_arr + COEFF_TO_OFF(*f)); // value is large
	
	return ret;
}

void F_mpz_get_mpz(mpz_t x, const F_mpz_t f)
{
	if (!COEFF_IS_MPZ(*f)) mpz_set_si(x, *f); // set x to small value
	else 
	{
		
		mpz_set(x, F_mpz_arr + COEFF_TO_OFF(*f)); // set x to large value
		
	}	
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
	{
		
		double ret = mpz_get_d_2exp(exp, F_mpz_arr + COEFF_TO_OFF(d));
		
		return ret;
	}
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
	else 
	{
		
		int ret = (mpz_cmp(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)) == 0); 
		
		return ret;
	}
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
	else 
	{
		
		int ret = mpz_cmpabs(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)); 
		
		return ret;
	}
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
			int ret = -1;
			
		   if (mpz_sgn(F_mpz_arr + COEFF_TO_OFF(*g)) < 0) ret = 1; // g is a large negative 
			
			return ret; // g is a large positive
		}
	} else if (!COEFF_IS_MPZ(*g)) 
	{
		int ret = 1;
		
		if (mpz_sgn(F_mpz_arr + COEFF_TO_OFF(*f)) < 0) ret = -1; // f is large negative
		
		return ret; // f is large positive
	} else // both f and g are large 
	{
		
		int ret = mpz_cmp(F_mpz_arr + COEFF_TO_OFF(*f), F_mpz_arr + COEFF_TO_OFF(*g)); 
		
		return ret;
	}
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

	
   ulong ret = mpz_size(F_mpz_arr + COEFF_TO_OFF(d));
	
	return ret;
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

	
   int ret = mpz_sgn(F_mpz_arr + COEFF_TO_OFF(d));
	
	return ret;
}

ulong F_mpz_bits(const F_mpz_t f)
{
	F_mpz d = *f;

	if (!COEFF_IS_MPZ(d)) // c1 is small
	{
		return FLINT_BIT_COUNT(FLINT_ABS(d));
	}

	
   ulong ret = mpz_sizeinbase(F_mpz_arr + COEFF_TO_OFF(d), 2);
	
	return ret;
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

void F_mpz_abs(F_mpz_t f1, const F_mpz_t f2)
{
   if (!COEFF_IS_MPZ(*f2)) // coeff is small
	{
		F_mpz t = FLINT_ABS(*f2); // Need to save value in case of aliasing
		
		_F_mpz_demote(f1);
		
		*f1 = t;
	} else // coeff is large
	{
	   // No need to retain value in promotion, as if aliased, both already large
		
		__mpz_struct * mpz_ptr = _F_mpz_promote(f1);
		mpz_abs(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(*f2));
		
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
		
		__mpz_struct * mpz_ptr2 = _F_mpz_promote(f); // g is already large
		__mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(c);
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
		
		__mpz_struct * mpz_ptr2 = _F_mpz_promote(f); // g is already large
		__mpz_struct * mpz_ptr = F_mpz_arr + COEFF_TO_OFF(c);
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

ulong F_mpz_mod_ui(F_mpz_t f, const F_mpz_t g, const ulong h)
{
	F_mpz c1 = *g;
	ulong r;
	
   if (!COEFF_IS_MPZ(c1)) // g is small
	{
		if (c1 < 0L) 
		{
			r = h - (-c1 % h); // C doesn't correctly handle negative mods
			if (r == h) r = 0;
		} else r = c1 % h;
		
		F_mpz_set_ui(f, r);
		return r;
	} else // g is large
	{
		
		r = mpz_fdiv_ui(F_mpz_arr + COEFF_TO_OFF(c1), h);
		
		F_mpz_set_ui(f, r);
		return r;
	}
}

void F_mpz_mod(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	ulong r;
	
   if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // h is also small
		{
			if (c2 < 0L) c2 = -c2;
			if (c1 < 0L) 
		   {
			   r = c2 - (-c1 % c2); // C doesn't correctly handle negative mods
			   if (r == c2) r = 0;
		   } else r = c1 % c2;
			
		   F_mpz_set_si(f, r);
		}
		else // h is large and g is small
		{
			if (c1 < 0L) 
			{
				F_mpz_abs(f, h);
			   F_mpz_sub_ui(f, f, -c1);			
			} else F_mpz_set_ui(f, c1);
		}
	} else // g is large
	{
      if (!COEFF_IS_MPZ(c2)) // h is small
		{
			if (c2 < 0L) F_mpz_set_si(f, mpz_fdiv_ui(F_mpz_arr + COEFF_TO_OFF(c1), -c2));
			else 
			{
				
				ulong r = mpz_fdiv_ui(F_mpz_arr + COEFF_TO_OFF(c1), c2);
				
				F_mpz_set_ui(f, r);
			}
		} else // both are large
		{
			
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			mpz_mod(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // reduction mod h may result in small value
			
		}	
	}
}

int F_mpz_invert(F_mpz_t f, const F_mpz_t g, const F_mpz_t h)
{
	F_mpz c1 = *g;
	F_mpz c2 = *h;
	int val;

   if (!COEFF_IS_MPZ(c1)) // g is small
	{
	   if (!COEFF_IS_MPZ(c2)) // h is also small
		{
			long inv;
			if (c2 < 0L) c2 = -c2;
			if (c2 == 1L)  return 0; // special case not handled by z_gcd_invert
			long gcd = z_gcd_invert(&inv, c1, c2);
			if (gcd == 1L) 
			{
				F_mpz_set_si(f, inv); // check gcd is 1
				return 1;
			} else return 0;
		} else // h is large and g is small
		{
			__mpz_struct temp; // put g into a temporary mpz_t
			if (c1 < 0L) 
			{
				c1 = -c1;
				temp._mp_d = &c1;
			   temp._mp_size = -1;
			} else if (c1 == 0L) temp._mp_size = 0;
			else 
			{
				temp._mp_d = &c1;
				temp._mp_size = 1;
			}
			
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			val = mpz_invert(mpz_ptr, &temp, F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // inverse mod h may result in small value
			
			return val;
		}
	} else // g is large
	{
      if (!COEFF_IS_MPZ(c2)) // h is small
		{
			long inv;
			if (c2 < 0L) c2 = -c2;
			if (c2 == 1L) return 0; // special case not handled by z_gcd_invert
			// reduce g mod h first
			
			ulong r = mpz_fdiv_ui(F_mpz_arr + COEFF_TO_OFF(c1), c2);
			
			long gcd = z_gcd_invert(&inv, r, c2);
			if (gcd == 1L) 
			{
				F_mpz_set_si(f, inv); // check gcd is 1
				return 1;
			} else return 0;
		} else // both are large
		{
			
			__mpz_struct * mpz_ptr = _F_mpz_promote(f);
			val = mpz_invert(mpz_ptr, F_mpz_arr + COEFF_TO_OFF(c1), F_mpz_arr + COEFF_TO_OFF(c2));
			_F_mpz_demote_val(f); // reduction mod h may result in small value
			
			return val;
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

void F_mpz_comb_init(F_mpz_comb_t comb, ulong * primes, ulong num_primes)
{
   ulong i, j, k;
	
	comb->primes = primes;
   comb->num_primes = num_primes;

   ulong n = ceil_log2(num_primes);
   comb->n = n;
	ulong num;

   // create zn_poly modulus information
	comb->mod = (zn_mod_t *) flint_heap_alloc_bytes(sizeof(zn_mod_t)*num_primes);
   for (ulong i = 0; i < num_primes; i++) 
      zn_mod_init(comb->mod[i], primes[i]);

	if (n == 0) return; // nothing to do

	// allocate space for comb and res
   comb->comb = (F_mpz **) flint_heap_alloc(n);
   comb->res = (F_mpz **) flint_heap_alloc(n);
   
   j = (1L << (n - 1)); // size of top level
   
	for (i = 0; i < n; i++) // initialise arrays at each level
   {
      comb->comb[i] = (F_mpz *) flint_heap_alloc(j);
      comb->res[i] = (F_mpz *) flint_heap_alloc(j);
		
		for (ulong k = 0; k < j; k++)
		{
			F_mpz_init(comb->comb[i] + k);
		   F_mpz_init(comb->res[i] + k);
		}

      j/=2;
   }

	// compute products of pairs of primes and place in comb
	for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
      F_mpz_set_ui(comb->comb[0] + j, primes[i]);
      F_mpz_mul_ui(comb->comb[0] + j, comb->comb[0] + j, primes[i+1]);
   }
   if (i < num_primes) // in case number of primes is odd
   {
      F_mpz_set_ui(comb->comb[0] + j, primes[i]);
	   i+=2;
	   j++;
	}
   num = (1L<<n); // set the rest of the entries on that row of the comb to 1
	for (; i < num; i += 2, j++)
   {
      F_mpz_set_ui(comb->comb[0] + j, 1L);
	}

	// compute rest of comb by multiplying in pairs
	ulong log_comb = 1;
	num /= 2;
   while (num >= 2)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         F_mpz_mul2(comb->comb[log_comb] + j, comb->comb[log_comb-1] + i, comb->comb[log_comb-1] + i + 1);
      }
      log_comb++;
      num /= 2;
   }

	// compute inverses from pairs of primes
	F_mpz_t temp, temp2;
	F_mpz_init(temp);
	F_mpz_init(temp2);
	
	for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
       F_mpz_set_ui(temp, primes[i]);
       F_mpz_set_ui(temp2, primes[i+1]);
       F_mpz_invert(comb->res[0] + j, temp, temp2);
	}
   
	F_mpz_clear(temp);
	F_mpz_clear(temp2);
   
   // compute remaining inverses, each level combining pairs from the level below
	ulong log_res = 1;
   num = (1L << (n - 1));

   while (log_res < n)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         F_mpz_invert(comb->res[log_res] + j, comb->comb[log_res-1] + i, comb->comb[log_res-1] + i + 1);
      }
      log_res++;
      num /= 2;
   }  
}

void F_mpz_comb_clear(F_mpz_comb_t comb)
{
   unsigned long n = comb->n;

   ulong j = (1L << (n - 1)); // size of top level
   
	for (ulong i = 0; i < n; i++) // initialise arrays at each level
   {
      for (ulong k = 0; k < j; k++)
		{
			F_mpz_clear(comb->comb[i] + k);
		   F_mpz_clear(comb->res[i] + k);
		}
      
		flint_heap_free(comb->comb[i]);
      flint_heap_free(comb->res[i]);

      j/=2;
   }
	
	if (n)
	{
		flint_heap_free(comb->res);
      flint_heap_free(comb->comb);
	}

   flint_heap_free(comb->mod);
}

void F_mpz_multi_mod_ui_basecase(ulong * out, F_mpz_t in, 
                               ulong * primes, ulong num_primes)
{
   F_mpz_t temp;
	F_mpz_init(temp);
	
	for (ulong i = 0; i < num_primes; i++)
   {
      out[i] = F_mpz_mod_ui(temp, in, primes[i]);
   }

	F_mpz_clear(temp);
}

void F_mpz_multi_mod_ui(ulong * out, F_mpz_t in, F_mpz_comb_t comb)
{
   ulong i, j, k;
   ulong n = comb->n;
   long log_comb;
	ulong num;
	ulong num_primes = comb->num_primes;

   if (num_primes == 1) // we are reducing modulo a single prime which is assumed to be big enough
   {
	  if ((*in) > 0L) out[0] = (*in);
	  else if ((*in) < 0L) out[0] = comb->primes[0] + (*in);
	  else out[0] = 0L;
	  return;
   }

   log_comb = n - 1;
   
	// allocate space for temp
	F_mpz ** temp = (F_mpz **) flint_heap_alloc(n);
   j = (1L << (n - 1));
   
   for (i = 0; i < n; i++)
   {
      temp[i] = (F_mpz *) flint_heap_alloc(j);

      for (k = 0; k < j; k++)
      {
         F_mpz_init(temp[i] + k);
      }

      j/=2;
   }

   // find level in comb with entries bigger than the input integer
	log_comb = 0;
	if ((*in) < 0L)
	  while ((F_mpz_bits(in) >= F_mpz_bits(comb->comb[log_comb]) - 1) && (log_comb < comb->n - 1)) log_comb++;
   else
		while (F_mpz_cmpabs(in, comb->comb[log_comb]) >= 0) log_comb++;
   num = (1L << (n - log_comb - 1));

	// set each entry of this level of temp to the input integer
	for (i = 0; i < num; i++)
   {
      F_mpz_set(temp[log_comb] + i, in);
   }
   log_comb--;
   num *= 2;
   
   // fill in other entries of temp by taking entries of temp at higher level mod pairs from comb
	while (log_comb > FLINT_F_MPZ_LOG_MULTI_MOD_CUTOFF) // keep going until we reach the basecase
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         F_mpz_mod(temp[log_comb] + i, temp[log_comb + 1] + j, comb->comb[log_comb] + i);
         F_mpz_mod(temp[log_comb] + i + 1, temp[log_comb + 1] + j, comb->comb[log_comb] + i + 1);
      }
      num *= 2;
      log_comb--;
   }
   
   // do basecase
	num /= 2;
   log_comb++;
   ulong stride = (1L << (log_comb + 1));
   for (i = 0, j = 0; j < num_primes; i++, j += stride)
   {
	   F_mpz_multi_mod_ui_basecase(out + j, temp[log_comb] + i, comb->primes + j, FLINT_MIN(stride, num_primes - j));
   }

	// free temp
	j = (1L << (n - 1));
   
   for (i = 0; i < n; i++)
   {
      for (k = 0; k < j; k++)
      {
         F_mpz_clear(temp[i] + k);
      }
      
		flint_heap_free(temp[i]);
      
		j/=2;
   }
	
	flint_heap_free(temp);
}

void __F_mpz_multi_CRT_ui_sign(F_mpz_t output, F_mpz_t input, F_mpz_comb_t comb)
{
   ulong n = comb->n;
   if (n == 0L) 
   {
      if (F_mpz_is_zero(input)) 
	   {
	      F_mpz_set_ui(output, 0L);
	      return;
	   }

	   long p = comb->primes[0];
		if ((p - (*input)) < (*input)) F_mpz_set_si(output, (long) ((*input) - p));
	   else F_mpz_set_ui(output, (*input));
	   return;
   }

   F_mpz_t temp;
	F_mpz_init(temp);
	
   F_mpz_sub(temp, input, comb->comb[comb->n - 1]);
	
   if (F_mpz_cmpabs(temp, input) <= 0) F_mpz_set(output, temp);
   else F_mpz_set(output, input);

   F_mpz_clear(temp);
   return;
}

void __F_mpz_multi_CRT_ui(F_mpz_t output, ulong * residues, F_mpz_comb_t comb, int sign)
{
   ulong i, j, k;

   ulong n = comb->n;
   ulong num;
   ulong log_res;
	ulong num_primes = comb->num_primes;

   if (num_primes == 1) // the output is less than a single prime, so just output the result
   {
      if (sign)
		{
		   ulong p = comb->primes[0];
	      if ((p - residues[0]) < residues[0]) F_mpz_set_si(output, residues[0] - p);
	      else F_mpz_set_ui(output, residues[0]);
	   } else F_mpz_set_ui(output, residues[0]);
	   
		return;
	}

   // allocate space for comb_temp
	F_mpz ** comb_temp = (F_mpz **) flint_heap_alloc(n);
   j = (1L << (n - 1));
   
	for (i = 0; i < n; i++)
   {
      comb_temp[i] = (F_mpz *) flint_heap_alloc(j);

      for (k = 0; k < j; k++)
      {
         F_mpz_init(comb_temp[i] + k);
      }

      j/=2;
   }

	// first layer of reconstruction
	num = (1L << n);
   F_mpz_t temp, temp2;
	F_mpz_init(temp);
	F_mpz_init(temp2);
	
	for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
   {
      F_mpz_set_ui(temp, residues[i]);
      F_mpz_mod_ui(temp2, temp, comb->primes[i+1]);
      F_mpz_sub_ui(temp2, temp2, residues[i + 1]);
      F_mpz_neg(temp2, temp2);
		F_mpz_mul2(temp, temp2, comb->res[0] + j);
      F_mpz_mod_ui(temp2, temp, comb->primes[i+1]);
      F_mpz_mul_ui(temp, temp2, comb->primes[i]); 
		F_mpz_add_ui(comb_temp[0] + j, temp, residues[i]);
   }

   if (i < num_primes) F_mpz_set_ui(comb_temp[0] + j, residues[i]);
    
   // compute other layers of reconstruction
	num /= 2;
   log_res = 1;
   
	while (log_res < n)
   {
      for (i = 0, j = 0; i < num; i += 2, j++)
      {
         if (F_mpz_is_one(comb->comb[log_res-1] + i + 1))
		   {
		      if (!F_mpz_is_one(comb->comb[log_res-1] + i)) 
					F_mpz_set(comb_temp[log_res] + j, comb_temp[log_res-1] + i);
		   } else
		   {
			   F_mpz_mod(temp2, comb_temp[log_res-1] + i, comb->comb[log_res-1] + i + 1);
            F_mpz_sub(temp, comb_temp[log_res-1] + i + 1, temp2);
            F_mpz_mul2(temp2, temp, comb->res[log_res] + j);
            F_mpz_mod(temp, temp2, comb->comb[log_res-1] + i + 1);
            F_mpz_mul2(temp2, temp, comb->comb[log_res-1] + i);
            F_mpz_add(comb_temp[log_res] + j, temp2, comb_temp[log_res-1] + i);
		   }
      }     
      log_res++;
      num /= 2; 
   }

	F_mpz_clear(temp);
   F_mpz_clear(temp2);

   // write out the output
	if (sign) __F_mpz_multi_CRT_ui_sign(output, comb_temp[log_res - 1], comb);
	else F_mpz_set(output, comb_temp[log_res - 1]);

	// free comb_temp
	j = (1L << (n - 1));
   
	for (i = 0; i < n; i++)
   {
      for (k = 0; k < j; k++)
      {
         F_mpz_clear(comb_temp[i] + k);
      }
      
		flint_heap_free(comb_temp[i]);

      j/=2;
   }
   
	flint_heap_free(comb_temp);
}

void F_mpz_multi_CRT_ui_unsigned(F_mpz_t output, ulong * residues, F_mpz_comb_t comb)
{
	__F_mpz_multi_CRT_ui(output, residues, comb, 0);
}

void F_mpz_multi_CRT_ui(F_mpz_t output, ulong * residues, F_mpz_comb_t comb)
{
	__F_mpz_multi_CRT_ui(output, residues, comb, 1);
}
