#ifndef MPN_EXTRAS_H
#define MPN_EXTRAS_H

void* limb_alloc(unsigned long, int);
void limb_release();



/*============================================================================

"mpn-wannabe" code.

These are functions that I wish were in GMP's mpn layer.

=============================================================================*/

/*
Computes the negation of a multiple-precision integer in 2's complement.

Input is count limbs stored at src. Output is stored at dest.

src and dest can be the same buffer. If they're not, they should be disjoint.

todo: currently this code will make only 1 pass over the data, EXCEPT in the
      case where all limbs are zero, in which case it will make two passes.
      FIX THIS!

todo: try writing another version that makes a block of zeroes and then
      uses mpn_sub_n repeatedly. This could be faster, if GMP's assembler
      is better than what gcc can come up with.

todo: consider writing this in assembly

todo: write test function for this

todo: consider using GMP's mpn_com_n (undocumented)

*/

static inline void negate_limbs(mp_limb_t* dest, mp_limb_t* src, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = ~src[i];
   mpn_add_1(dest, dest, count, 1);
}

/*
Copies a bunch of limbs from one buffer to another.

Input is count limbs stored at src. Output is stored at dest.

src and dest can be the same buffer. If they're not, they should be disjoint.

todo: it's completely insane that we're not using memcpy. But memcpy seems
      to have crazy overhead and is slow!! Why is this?

todo: GMP has code to do limb copying. Clearly memcpy wasn't good enough for
      them either. Work out how to use their code. It's not a documented
      interface, so some hackishness may be necessary.
 */
static inline void copy_limbs(mp_limb_t* dest, mp_limb_t* src, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
   {
      dest[i] = src[i];
   }
}

static inline void forward_copy_limbs(mp_limb_t* dest, mp_limb_t* src, unsigned long count)
{
   for (long i = 0; i < count; i++)
   {
      dest[i] = src[i];
   }
}


/*
Sets a bunch of limbs to zero.

todo: why does memset have so much overhead????!!?
 */
static inline void clear_limbs(mp_limb_t* dest, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = 0;
}

/*
Sets a bunch of limbs to 0xfff....

todo: why does memset have so much overhead????!!?
 */
static inline void set_limbs(mp_limb_t* dest, unsigned long count)
{
   for (long i = count - 1; i >= 0; i--)
      dest[i] = (mp_limb_t)(-1L);
}


/*
   Add the given *signed* limb to the buffer [x, x+count), much like
   mpn_add_1 and mpn_sub_1 (except it's always inplace).

   PRECONDITIONS:
      count >= 1
   
   NOTE:
      The branch predictability of this function is optimised for the case that
      abs(limb) is relatively small and that the first limb of x is randomly
      distributed, which should be the normal usage in the FFT routines.
*/
static inline
void signed_add_1(mp_limb_t* x, unsigned long count, mp_limb_signed_t limb)
{
   FLINT_ASSERT(count >= 1);
   
   // If the high bit of x[0] doesn't change when we add "limb" to it,
   // then there's no possibility of overflow.
   mp_limb_t temp = x[0] + limb;
   if ((mp_limb_signed_t)(temp ^ x[0]) >= 0)
      // the likely case
      x[0] = temp;
   else
   {
      // the unlikely case; here we need to branch based on the sign of
      // the limb being added
      if (limb >= 0)
         mpn_add_1(x, x, count, limb);
      else
         mpn_sub_1(x, x, count, -limb);
   }
}


mp_limb_t mpn_divmod_1_preinv(mp_limb_t * qp, mp_limb_t * up, 
             unsigned long un, mp_limb_t d, mp_limb_t dinv, unsigned long norm);
             
mp_limb_t mpn_addmul(mp_limb_t * rp, mp_limb_t * s1p, unsigned long s1n, 
                                      mp_limb_t * s2p, unsigned long s2n);

#endif

