/******************************************************************************

Code for Schoenhage-Strassen and Kronecker-Schoenhage polynomial 
         multiplication over Z.

 (C) 2006 William Hart and David Harvey

This implementation of SS/KS multiplication uses the following optimisations.

* Coefficients are represented as raw arrays of limbs. Arithmetic is performed
mainly with GMP's mpn layer. We don't use mpz_t's; they are too slow. We make 
certain assumptions about the underlying structure of GMP's mpz_t's, allowing 
us to convert data in and out to raw limbs quickly. GMP_COMPATIBLE mode (see
flint.h for the relevant #define) will not make such assumptions for 
guaranteed compatibility with future releases of GMP, however it has not been 
implemented for ssmul.c yet.

* We reduce mod p as rarely as possible. Instead, we let bits build up in the
last limb then reduce only when necessary.

* Butterflies are accomplished using as few passes over the data as possible.
In particular, for those operations that can't be performed in-place, we swap
pointers rather than making unnecessary copies of the data.

* We use a cache-friendly algorithm for large transform lengths. It
decomposes the transform into many much smaller transforms which utilise the
cache much better. We also supply some cache hints (for compilers/processors
that understand them).

* We make use of our SS code to do fast integer multiplication and use this 
for the pointwise multiplications when it is faster than GMP.

* For very large coefficients sizes, we split the coefficients up into smaller
ones. This technique was developed by William Hart.

* In the region where KS is used (coefficient bit size < degree) we do not pile
all the coefficients up into two large integers and multiply, but we pile
sufficiently many up so that the SS code can handle it optimally. This technique
was independently suggested by Paul Zimmerman and David Harvey.

* For very small coefficient sizes (where the output coefficients will fit in 
64 bits) we use an efficient bit packing routine and use ordinary KS. The routine
reads and writes limbs only, but shifts the data through a 64 bit integer as 
needed. These techniques were developed/coded by William Hart.

* TODO: Use a negacyclic convolution to reduce the time taken for pointwise
multiplication modulo p.

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <pthread.h>
#include "flint.h"
#include "Z-ssmul.h"
#include "mpn_extras.h"
#include "Zvec.h"
#include "ssmul.h"
#include "ssfft.h"

// Whether we should use a multi-threaded version of ssmul 
#define USE_THREADS 0
// Log_2 of number of threads (currently only up to 16 threads are allowed)
#define LG_THREADS 1

#define THREADS (1<<LG_THREADS)

// Whether we should use the new truncated FFT where possible
// Note: if using threads this should be 0 for now
#define USE_TRUNCATED_FFT 1

#if FLINT_BITS_PER_LIMB == 32
#define half_ulong u_int16_t
#define half_long int16_t
#define HALF_FLINT_BITS_PER_LIMB 16
#else
#define half_ulong u_int32_t
#define half_long int32_t
#define HALF_FLINT_BITS_PER_LIMB 32
#endif

typedef struct fft_ss
{
   mp_limb_t** start;
   unsigned long length;
   unsigned long skip;
   unsigned long next_skip;
   unsigned long start_r;
   unsigned long next_skip_r; 
   unsigned long skip_r;
   unsigned long depth;
   mp_limb_t** scratch;
   unsigned long n;
   int first;
   unsigned long * thread;
} fft_t;

typedef struct fft_s
{
   mp_limb_t** start;
   unsigned long length;
   mp_limb_t** scratch;
   unsigned long r; 
   unsigned long n;
} fft_s;

typedef struct point_t
{
   unsigned long length;
   mp_limb_t** array1;
   mp_limb_t** array2;
   mp_limb_t* array3;
   unsigned long n;
} point_t;

/*============================================================================

For some strange reason, if one shifts by FLINT_BITS_PER_LIMB, the result is 
not zero as one would expect. Therefore we provide these shift functions
which take care of this case as required by the convert in/out routines below.

=============================================================================*/

inline unsigned long shift_r(unsigned long n, unsigned long r)
{
   if (r == FLINT_BITS_PER_LIMB) return 0;
   return (n>>r);
}

inline unsigned long shift_l(unsigned long n, unsigned long l)
{
   if (l == FLINT_BITS_PER_LIMB) return 0;
   return (n<<l);
}
 
/*============================================================================

Routines for converting data in from mpz_t's to arrays of limbs and for
converting out from arrays of limbs to mpz_t's. 

=============================================================================*/
  
/*
Packs bundle input coefficients down to the bit into an array. It assumes 
each coefficient is less than FLINT_BITS_PER_LIMB bits and that they are being packed into 
bitfields of the given number of bits wide. Assumes bits < FLINT_BITS_PER_LIMB. Assumes n 
is set to the length of the array in limbs and that coeffs_per_limb is the
maximum number of whole coefficients that can be packed into a limb given
bitfields of the given number of bits wide.

Unlike the signed version below, we rarely deal explicitly with the case 
where k > FLINT_BITS_PER_LIMB. After each section, k is guaranteed to be less than HALF_FLINT_BITS_PER_LIMB. When 
the next coefficient is packed in, all the non-zero bits of the coefficient 
(which cannot number more than HALF_FLINT_BITS_PER_LIMB bits) can fit inside a FLINT_BITS_PER_LIMB bit integer 
even after shifted left by k bits. Thus although k may now be greater than
FLINT_BITS_PER_LIMB, the bits which overhang are all just zero padding in the bitfield and 
when a right shift is done, zeroes are shifted in. Thus we never need to 
explicitly keep track of the overhanging bits.
*/
void ssmul_convert_in_bits(mp_limb_t* array, const mpz_t* data, 
          const unsigned long n, const unsigned long bundle, 
          const unsigned long bits, const unsigned long coeffs_per_limb)
{   
   unsigned long k, l, skip;

   unsigned long temp = 0;
   half_ulong lower;
   
   // pile up input coefficients into single large coefficient
   
    k=0; skip=0; l = 0;
    while (l+coeffs_per_limb+1 < bundle)
    {
       // k is guaranteed to be less than FLINT_BITS_PER_LIMB at this point
       while (k<HALF_FLINT_BITS_PER_LIMB)
       {
          temp+=((data[l]->_mp_d[0])<<k); 
          l++;
          k+=bits;
       }
       // k may now exceed FLINT_BITS_PER_LIMB and is guaranteed to be >= HALF_FLINT_BITS_PER_LIMB

       if (k>FLINT_BITS_PER_LIMB)
       {
          // if k exceeds FLINT_BITS_PER_LIMB, write out a whole limb and reduce k by FLINT_BITS_PER_LIMB
          array[skip] = temp;
          skip++;
          temp=((data[l-1]->_mp_d[0])>>(bits+FLINT_BITS_PER_LIMB-k));
          k=(k-FLINT_BITS_PER_LIMB);
          // k is now guaranteed to be less than HALF_FLINT_BITS_PER_LIMB
       } else
       {
          // k < FLINT_BITS_PER_LIMB
          lower = (half_ulong)temp;
          k-=HALF_FLINT_BITS_PER_LIMB;
          temp>>=HALF_FLINT_BITS_PER_LIMB;
          // k is now guaranteed to be less than HALF_FLINT_BITS_PER_LIMB

          while (k<HALF_FLINT_BITS_PER_LIMB)
          {
             temp+=((data[l]->_mp_d[0])<<k);  
             l++;
             k+=bits;
          }
          // k may now exceed FLINT_BITS_PER_LIMB
          k-=HALF_FLINT_BITS_PER_LIMB;
          array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
          temp>>=HALF_FLINT_BITS_PER_LIMB;
          skip++;
          // k is now guaranteed to be less than FLINT_BITS_PER_LIMB
      }
   }
   // deal with remainder of coefficients

   while (l< bundle)
   {
       // k is guaranteed to be less than FLINT_BITS_PER_LIMB at this point
       while ((k<HALF_FLINT_BITS_PER_LIMB)&(l<bundle))
       {
          temp+=((data[l]->_mp_d[0])<<k);  
          l++;
          k+=bits;
       }
       // k may now exceed FLINT_BITS_PER_LIMB and is guaranteed to be >= HALF_FLINT_BITS_PER_LIMB

       if (k>FLINT_BITS_PER_LIMB)
       {
          // if k exceeds FLINT_BITS_PER_LIMB, write out a whole limb and reduce k by FLINT_BITS_PER_LIMB
          array[skip] = temp;
          skip++;
          temp=((data[l-1]->_mp_d[0])>>(bits+FLINT_BITS_PER_LIMB-k));
          k=(k-FLINT_BITS_PER_LIMB);
          // k is now guaranteed to be less than HALF_FLINT_BITS_PER_LIMB
       } else
       {
          // k < FLINT_BITS_PER_LIMB
          lower = (half_ulong)temp;
          k-=HALF_FLINT_BITS_PER_LIMB;
          temp>>=HALF_FLINT_BITS_PER_LIMB;
          // k is now guaranteed to be less than HALF_FLINT_BITS_PER_LIMB

          while ((k<HALF_FLINT_BITS_PER_LIMB)&(l<bundle))
          {
             temp+=((data[l]->_mp_d[0])<<k);  
             l++;
             k+=bits;
          }
          // k may now exceed FLINT_BITS_PER_LIMB
          k-=HALF_FLINT_BITS_PER_LIMB;
          array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
          temp>>=HALF_FLINT_BITS_PER_LIMB;
          skip++;
          // k is now guaranteed to be less than FLINT_BITS_PER_LIMB
      }
   }
   if (skip < n) 
   {
     array[skip] = temp;
   } 
}

/*
Packs bundle signed input coefficients down to the bit into an array. It 
assumes each coefficient is less than HALF_FLINT_BITS_PER_LIMB bits and that they are being 
packed into bitfields of the given number of bits wide. Assumes bits < FLINT_BITS_PER_LIMB. 
Assumes n is set to the length of the array in limbs and that coeffs_per_limb 
is the maximum number of whole coefficients that can be packed into a limb 
given bitfields of the given number of bits wide.
*/
void ssmul_convert_in_bits_signed(mp_limb_t* array, const mpz_t* data, 
          const unsigned long n, const unsigned long bundle, 
          const unsigned long bits, const unsigned long coeffs_per_limb, int sign_extend)
{   
    unsigned long k, l, m, skip;

    unsigned long temp = 0;
    half_ulong lower;
    long coeff;
    long borrow = 0UL;

    const unsigned long mask = (1UL<<bits)-1;
   
    k=0; skip=0; l = 0;
    while (l+coeffs_per_limb+1 < bundle)
    {
       // k is guaranteed to be < FLINT_BITS_PER_LIMB at this point
       while (k<HALF_FLINT_BITS_PER_LIMB)
       {
          if (mpz_sgn(data[l])>=0) coeff = data[l]->_mp_d[0] - borrow;
          else coeff = (-data[l]->_mp_d[0] - borrow);
          borrow = 0UL;
          if (coeff < 0) borrow = 1UL;
          coeff&=mask;
          temp+=(coeff<<k);
          l++;
          k+=bits;
       }
       // k may now exceed FLINT_BITS_PER_LIMB but is less than 96
       if (k>FLINT_BITS_PER_LIMB)
       {
          //if k > FLINT_BITS_PER_LIMB write out whole FLINT_BITS_PER_LIMB bit limb and load the remaining bits from coeff
          array[skip] = temp;
          skip++;
          temp=(coeff>>(bits+FLINT_BITS_PER_LIMB-k));
          k=(k-FLINT_BITS_PER_LIMB);
          // k is now guaranteed to be < HALF_FLINT_BITS_PER_LIMB
       } else
       {
          // k <= FLINT_BITS_PER_LIMB
          lower = (half_ulong)temp;
          k-=HALF_FLINT_BITS_PER_LIMB;
          temp>>=HALF_FLINT_BITS_PER_LIMB;
          // k <= HALF_FLINT_BITS_PER_LIMB
          while (k<HALF_FLINT_BITS_PER_LIMB)
          {
             if (mpz_sgn(data[l])>=0) coeff = data[l]->_mp_d[0] - borrow;
             else coeff = (-data[l]->_mp_d[0] - borrow);
             borrow = 0UL;
             if (coeff < 0) borrow = 1UL;
             coeff&=mask;
             temp+=(coeff<<k);
             l++;
             k+=bits;
          }
          // k may now exceed FLINT_BITS_PER_LIMB but is less than 96
          if (k>FLINT_BITS_PER_LIMB)
          {
             // write out bottom HALF_FLINT_BITS_PER_LIMB bits of temp (along with HALF_FLINT_BITS_PER_LIMB bits from lower)
             // load remaining bits of coeff and reduce k by HALF_FLINT_BITS_PER_LIMB
             array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+(unsigned long)lower;
             skip++;
             temp>>=HALF_FLINT_BITS_PER_LIMB;
             temp+=((coeff>>(bits+FLINT_BITS_PER_LIMB-k))<<HALF_FLINT_BITS_PER_LIMB);
             k=(k-HALF_FLINT_BITS_PER_LIMB);
             // k is now less than FLINT_BITS_PER_LIMB and we are ready to read in the next coefficient
          } else
          {
             // k <= FLINT_BITS_PER_LIMB
             k-=HALF_FLINT_BITS_PER_LIMB;
             array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
             temp>>=HALF_FLINT_BITS_PER_LIMB;
             skip++;
             // k <= HALF_FLINT_BITS_PER_LIMB and we are ready to read in the next coefficient
          }
      }
   }
   // deal with remainder of coefficients

   while (l< bundle)
   {
       // k is guaranteed to be less than FLINT_BITS_PER_LIMB at this point
       while ((k<HALF_FLINT_BITS_PER_LIMB)&(l<bundle))
       {
          if (mpz_sgn(data[l])>=0) coeff = data[l]->_mp_d[0] - borrow;
          else coeff = (-data[l]->_mp_d[0] - borrow);
          borrow = 0UL;
          if (coeff < 0) borrow = 1UL;
          coeff&=mask;
          temp+=(coeff<<k);
          l++;
          k+=bits;
       }
       // k may exceed FLINT_BITS_PER_LIMB at this point but is less than 96

       if (k>FLINT_BITS_PER_LIMB)
       {
          // if k > FLINT_BITS_PER_LIMB write out a whole limb and read in remaining bits of coeff
          array[skip] = temp;
          skip++;
          temp=(coeff>>(bits+FLINT_BITS_PER_LIMB-k));
          k=(k-FLINT_BITS_PER_LIMB);
          // k < HALF_FLINT_BITS_PER_LIMB
       } else
       {
          // k <= FLINT_BITS_PER_LIMB
          if (k >= HALF_FLINT_BITS_PER_LIMB)
          {
             // if k >= HALF_FLINT_BITS_PER_LIMB store bottom HALF_FLINT_BITS_PER_LIMB bits
             lower = (half_ulong)temp;
             k-=HALF_FLINT_BITS_PER_LIMB;
             temp>>=HALF_FLINT_BITS_PER_LIMB;
             // k is now <= HALF_FLINT_BITS_PER_LIMB

             while ((k<HALF_FLINT_BITS_PER_LIMB)&(l<bundle))
             {
                if (mpz_sgn(data[l])>=0) coeff = data[l]->_mp_d[0] - borrow;
                else coeff = (-data[l]->_mp_d[0] - borrow);
                borrow = 0UL;
                if (coeff < 0) borrow = 1UL;
                coeff&=mask;
                temp+=(coeff<<k);
                l++;
                k+=bits;
             }
             // k may again exceed FLINT_BITS_PER_LIMB bits but is less than 96
             if (k>FLINT_BITS_PER_LIMB)
             {
                // if k > FLINT_BITS_PER_LIMB, write out bottom HALF_FLINT_BITS_PER_LIMB bits (along with HALF_FLINT_BITS_PER_LIMB bits from lower)
                // read remaining bits from coeff and reduce k by HALF_FLINT_BITS_PER_LIMB
                array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+(unsigned long)lower;
                skip++;
                temp>>=HALF_FLINT_BITS_PER_LIMB;
                temp+=((coeff>>(bits+FLINT_BITS_PER_LIMB-k))<<HALF_FLINT_BITS_PER_LIMB);
                k=(k-HALF_FLINT_BITS_PER_LIMB);
                // k < FLINT_BITS_PER_LIMB and we are ready to read next coefficient if there is one
             } else if (k >= HALF_FLINT_BITS_PER_LIMB) 
             {
                // k <= FLINT_BITS_PER_LIMB
                // if k >= HALF_FLINT_BITS_PER_LIMB write out bottom HALF_FLINT_BITS_PER_LIMB bits (along with lower)
                // and reduce k by HALF_FLINT_BITS_PER_LIMB
                k-=HALF_FLINT_BITS_PER_LIMB;
                array[skip] = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
                temp>>=HALF_FLINT_BITS_PER_LIMB;
                skip++;
                // k is now less than or equal to HALF_FLINT_BITS_PER_LIMB and we are now ready to read 
                // the next coefficient if there is one
             } else
             {
                // k < HALF_FLINT_BITS_PER_LIMB
                // there isn't enough to write out a whole FLINT_BITS_PER_LIMB bits, so put it all 
                // together in temp
                temp = (temp<<HALF_FLINT_BITS_PER_LIMB)+lower;
                k+=HALF_FLINT_BITS_PER_LIMB;
                // k is now guaranteed to be less than FLINT_BITS_PER_LIMB and we are ready for the
                // next coefficient if there is one
             }
          }
      }
   }

   // sign extend the last FLINT_BITS_PER_LIMB bits we write out
   if (skip < n)
   {
     if (borrow) temp+=(-1UL << k);
     array[skip] = temp;
     skip++;
   } 
   // if sign_extend was set by caller, sign extend the remainder of the array,
   // reducing modulo p 
   if ((sign_extend) && (borrow))
   {
      while (skip < n) 
      {
         array[skip] = -1UL;
         skip++;
      }
      mpn_add_1(array, array, n+1, 1UL);
   }
}

/*
Convert a bundle of bit packed coefficients out to mpz_t's. Assumes bits < FLINT_BITS_PER_LIMB
*/
void ssmul_convert_out_bits(mpz_t* res, const mp_limb_t* array, 
            const unsigned long bundle, const unsigned long bits, 
            const unsigned long coeffs_per_limb)
{
    unsigned long k, l, skip;

    unsigned long temp = 0;
    unsigned long full_limb;
    half_ulong lower;
    
    const unsigned long mask = (1UL<<bits)-1;

    unsigned long s;
    k=0; skip=0; l = 0;
    while (l+coeffs_per_limb+1 < bundle)
    {
       // read in a FLINT_BITS_PER_LIMB bit limb 
       full_limb = array[skip];
       temp += shift_l(full_limb,k);
       s=FLINT_BITS_PER_LIMB-k;
       k+=s;
       while (k>=bits)
       {
          (res[l]->_mp_d)[0]+=(temp&mask);
          res[l]->_mp_size = 1;
          if (!(res[l]->_mp_d)[0]) res[l]->_mp_size = 0;
          temp>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       // read in the remainder of full_limb 
       temp += shift_l(shift_r(full_limb,s),k);
       k+=(FLINT_BITS_PER_LIMB-s);
       
       while (k>=bits)
       {
          (res[l]->_mp_d)[0]+=(temp&mask);
          res[l]->_mp_size = 1;
          if (!(res[l]->_mp_d)[0]) res[l]->_mp_size = 0;
          temp>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       skip++;
    }
   
    //deal with any remaining coefficients
    while (l< bundle)
    {
       // read in a FLINT_BITS_PER_LIMB bit limb
       full_limb = array[skip];
       temp += shift_l(full_limb,k);
       s=FLINT_BITS_PER_LIMB-k;
       k+=s;
       while ((k>=bits)&&(l<bundle))
       {
          (res[l]->_mp_d)[0]+=(temp&mask);
          res[l]->_mp_size = 1;
          if (!(res[l]->_mp_d)[0]) res[l]->_mp_size = 0;
          temp>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       // read in the remainder of full_limb 
       temp += shift_l(shift_r(full_limb,s),k);
       k+=(FLINT_BITS_PER_LIMB-s);
       
       while ((k>=bits)&&(l<bundle))
       {
          (res[l]->_mp_d)[0]+=(temp&mask);
          res[l]->_mp_size = 1;
          if (!(res[l]->_mp_d)[0]) res[l]->_mp_size = 0;
          temp>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       skip++;
    }
}

/*
Convert a bundle of signed bit packed coefficients out to mpz_t's. Assumes bits < FLINT_BITS_PER_LIMB
*/
void ssmul_convert_out_bits_signed(mpz_t* res, const mp_limb_t* array, 
            const unsigned long bundle, const unsigned long bits, 
            const unsigned long coeffs_per_limb)
{
    unsigned long k, l, skip;

    unsigned long temp2 = 0;
    unsigned long temp;
    unsigned long full_limb;
    half_ulong lower;
    unsigned long carry = 0UL;
    
    const unsigned long mask = (1UL<<bits)-1;
    const unsigned long sign_mask = (1UL<<(bits-1));

    unsigned long s;
    k=0; skip=0; l = 0;
    while (l+coeffs_per_limb+1 < bundle)
    {
       // read in a FLINT_BITS_PER_LIMB bit limb
       full_limb = array[skip];
       temp2 += shift_l(full_limb,k);
       s=FLINT_BITS_PER_LIMB-k;
       k+=s;
       while (k>=bits)
       {
          if (!(temp2&sign_mask)) 
          {
             mpz_add_ui(res[l],res[l],(temp2&mask)+carry);
             carry = 0UL;
          }
          else
          {
             temp = ((-temp2)&mask)-carry;
             mpz_sub_ui(res[l],res[l],temp);
             carry = 1UL;
          }
          temp2>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       // read in remainder of full_limb
       temp2 += shift_l(shift_r(full_limb,s),k);
       k+=(FLINT_BITS_PER_LIMB-s);
       
       while (k>=bits)
       {
          if (!(temp2&sign_mask)) 
          {
             mpz_add_ui(res[l],res[l],(temp2&mask)+carry);
             carry = 0UL;
          }
          else
          {
             temp = ((-temp2)&mask)-carry;
             mpz_sub_ui(res[l],res[l],temp);
             carry = 1UL;
          }
          temp2>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       skip++;
    }
   
    //deal with any remaining coefficients
    while (l< bundle)
    {
       // read in a full limb
       full_limb = array[skip];
       temp2 += shift_l(full_limb,k);
       s=FLINT_BITS_PER_LIMB-k;
       k+=s;
       while ((k>=bits)&&(l<bundle))
       {
          if (!(temp2&sign_mask)) 
          {
             mpz_add_ui(res[l],res[l],(temp2&mask)+carry);
             carry = 0UL;
          }  
          else
          {
             temp = ((-temp2)&mask)-carry;
             mpz_sub_ui(res[l],res[l],temp);
             carry = 1UL;
          }
          temp2>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       // read in remainder of full_limb
       temp2 += shift_l(shift_r(full_limb,s),k);
       k+=(FLINT_BITS_PER_LIMB-s);
       
       while ((k>=bits)&&(l<bundle))
       {
          if (!(temp2&sign_mask)) 
          {
             mpz_add_ui(res[l],res[l],(temp2&mask)+carry);
             carry = 0UL;
          }
          else
          {
             temp = ((-temp2)&mask)-carry;
             mpz_sub_ui(res[l],res[l],temp);
             carry = 1UL;
          }
          temp2>>=bits;
          l++;
          k-=bits;
       }
       // k is now less than bits
       skip++;
    }
}

/*
Convert a bundle of input coefficients from mpz_t format to an array of 
limbs, representing a single FFT coefficient. The coefficients are packed 
to the nearest byte. Each coefficient is to take up coeff_bytes. 

** IMPORTANT ** This code assumes the input coefficient size is at least a FLINT_BITS_PER_LIMB bit
limb. This is no restriction, since it wouldn't be efficient for smaller
coefficients anyway.   
*/
void ssmul_convert_in_bytes(mp_limb_t* array, mpz_t* coeffs, 
            const unsigned long bundle, const unsigned long coeff_bytes, const unsigned long total_limbs)
{
    const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
    // Start limb of the current coefficient within array
    unsigned long coeff_limb = 0;
    // Additional offset in bytes of current coefficient within array
    unsigned long coeff_byte = 0;
    
    // Where we are up to in the current coefficient: limbs + bytes
    unsigned long offset_limb = coeff_limb;
    
    // Next limb to be written to bytes
    unsigned long next_limb;
    unsigned long temp = 0;

    unsigned long shift_1, shift_2;
    
    // when a coefficient is negative, we need to borrow from the next coefficient
    int borrow = 0;
    int borrowed = 0;
    
    unsigned long i, j;
    
    for (i = 0; i < bundle; i++)
    {
        // compute shifts to be used
        shift_1 = coeff_byte<<3;
        shift_2 = FLINT_BITS_PER_LIMB-shift_1;
        
        borrowed = borrow;
        if (borrow) mpz_sub_ui(coeffs[i],coeffs[i],1);
        
        //negative coefficient
        if (mpz_sgn(coeffs[i]) < 0) 
        {
           // mpz_t's store the absolute value only, so add 1 then complement
           mpz_add_ui(coeffs[i],coeffs[i],1);
           // deal with first limb of coefficient
           next_limb = ~(coeffs[i]->_mp_d)[0];
           if (limbs_per_coeff == 0) 
           {
              if (i == bundle-1) 
              {
                 temp += shift_l(next_limb,shift_1);
                 array[offset_limb] = temp + ((shift_l(1UL,shift_1)-1)&array[offset_limb]);
                 offset_limb++;
                 temp = shift_r(next_limb,shift_2);
                 temp += shift_l(-1UL,shift_1);
                 array[offset_limb] = temp;
                 offset_limb++;
                 while (offset_limb < total_limbs)
                 {
                    array[offset_limb] = -1UL;
                    offset_limb++;
                 }
              } else 
              {
                 next_limb &= ((1UL<<(extra_bytes_per_coeff<<3))-1);
                 temp += shift_l(next_limb,shift_1);
                 array[offset_limb] = temp + ((shift_l(1UL,shift_1)-1)&array[offset_limb]);
                 offset_limb++;
                 temp = shift_r(next_limb,shift_2);
                 array[offset_limb] = temp;
              }
           } else
           {
              temp += shift_l(next_limb,shift_1);
              array[offset_limb] = temp + ((shift_l(1UL,shift_1)-1)&array[offset_limb]);
              offset_limb++;
              temp = shift_r(next_limb,shift_2);
              // deal with remaining limbs
              for (j = 1; j < mpz_size(coeffs[i]); j++, offset_limb++)
              {
                 next_limb = ~(coeffs[i]->_mp_d)[j];
                 temp += shift_l(next_limb,shift_1);
                 array[offset_limb] = temp;
                 temp = shift_r(next_limb,shift_2);
              }
              // write remaining part of coefficient and fill 
              // remainder of coeff_bytes with binary 1's
              if ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                  limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte) 
              {
                 temp += shift_l(-1UL,shift_1);
                 array[offset_limb] = temp;
                 offset_limb++;
              }
              for (; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
              {
                 array[offset_limb] = -1UL;
              }
              while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                  limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
              {
                 array[offset_limb] = -1UL;
                 offset_limb++;
              }
              if (i == bundle-1)
              {
                 while (offset_limb < total_limbs)
                 {
                    array[offset_limb] = -1UL;
                    offset_limb++;
                 }
              }
           }
           temp = 0;
           // restore the original coefficient
           mpz_sub_ui(coeffs[i],coeffs[i],1);
           borrow = 1;
        }
        
        //positive coefficient
        else if (mpz_sgn(coeffs[i]) > 0)
        {
           // deal with first limb of coefficient
           next_limb = (coeffs[i]->_mp_d)[0];
           temp = shift_l(next_limb,shift_1);
           array[offset_limb] = temp + ((shift_l(1UL,shift_1)-1)&array[offset_limb]);
           offset_limb++;
           temp = shift_r(next_limb,shift_2);
           if (shift_2 == FLINT_BITS_PER_LIMB) temp = 0;
           // deal with remaining limbs
           for (j = 1; j < mpz_size(coeffs[i]); j++)
           {
              next_limb = (coeffs[i]->_mp_d)[j];
              temp += shift_l(next_limb,shift_1);
              array[offset_limb] = temp;
              offset_limb++;
              temp = shift_r(next_limb,shift_2);
           }
           // write remaining part of coefficient
           array[offset_limb] = temp;
           temp = 0;
           borrow = 0;
        }

        else
        {
           array[offset_limb] = ((shift_l(1UL,shift_1)-1)&array[offset_limb]);
           temp = 0;
           borrow = 0;
        }
        // update information for next coefficient
        coeff_limb += limbs_per_coeff;
        coeff_byte += extra_bytes_per_coeff;
        if (coeff_byte > FLINT_BYTES_PER_LIMB) 
        {
           coeff_byte -= FLINT_BYTES_PER_LIMB;
           coeff_limb++;
        }
        offset_limb = coeff_limb;
        // set coefficient back to what it was before borrow
        if (borrowed) mpz_add_ui(coeffs[i],coeffs[i],1);
    }
}

/*
Used for unpacking a single coefficient which has been byte packed into a larger one.
This is used by ssmul_convert_out_bytes below. limb_start gives the offset in limbs
of the first limb containing bytes of the coefficient of interest. byte_start gives 
the offset in bytes into that first limb that gives the starting position of the 
coefficient to be extracted. num_bytes is the number of bytes that the coefficient
takes up. 
*/
inline unsigned long ssmul_unpack_signed_bytes(mp_limb_t* output, mp_limb_t* array, 
            const unsigned long limb_start, const unsigned long byte_start, 
            const unsigned long num_bytes)
{
    const unsigned long limbs_to_extract = (num_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const unsigned long extra_bytes_to_extract = num_bytes 
                                    - (limbs_to_extract<<FLINT_LG_BYTES_PER_LIMB);
                                    
    unsigned long next_limb;
    unsigned long temp = 0;
    
    // the limb we are up to in the array and output respectively
    unsigned long coeff_limb = limb_start;
    unsigned long output_limb = 0;

    unsigned long shift_1, shift_2;
    
    unsigned long i, j;
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS_PER_LIMB - shift_1;
    
    unsigned long sign;

    if (byte_start + extra_bytes_to_extract > FLINT_BYTES_PER_LIMB)
    {
       sign = array[limb_start+limbs_to_extract+1]&(1UL<<(((byte_start 
            + extra_bytes_to_extract - FLINT_BYTES_PER_LIMB)<<3)-1));
    } else if (byte_start + extra_bytes_to_extract == FLINT_BYTES_PER_LIMB)
    {
       sign = array[limb_start+limbs_to_extract]&(1UL<<(FLINT_BITS_PER_LIMB-1));
    } else if (byte_start + extra_bytes_to_extract == 0)
    {
       sign = array[limb_start+limbs_to_extract-1]&(1UL<<(FLINT_BITS_PER_LIMB-1));
    } else
    {
       sign = array[limb_start+limbs_to_extract]&(1UL<<(((byte_start 
            + extra_bytes_to_extract)<<3)-1));
    }
    
    if (sign)
    {
        temp = ~array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract)
        {
           next_limb = shift_r(temp,shift_1);
           temp = ~array[coeff_limb];
           coeff_limb++;
           next_limb += shift_l(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
        {
           next_limb = shift_r(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else
        {
           next_limb = shift_r(temp,shift_1);
           temp = ~array[coeff_limb];
           next_limb += shift_l(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    else
    {
        temp = array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract)
        {
           next_limb = shift_r(temp,shift_1);
           temp = array[coeff_limb];
           coeff_limb++;
           next_limb += shift_l(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
        {
           next_limb = shift_r(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else
        {
           next_limb = shift_r(temp,shift_1);
           temp = array[coeff_limb];
           next_limb += shift_l(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    return sign;
}   

/*
Used for unpacking a single signed coefficient which has been byte packed into a 
larger one. This is used by ssmul_convert_out_bytes_signed below. limb_start gives 
the offset in limbs of the first limb containing bytes of the coefficient of 
interest. byte_start gives the offset in bytes into that first limb that gives the 
starting position of the coefficient to be extracted. num_bytes is the number of 
bytes that the coefficient takes up. 
*/
inline void ssmul_unpack_bytes(mp_limb_t* output, mp_limb_t* array, 
            const unsigned long limb_start, const unsigned long byte_start, 
            const unsigned long num_bytes)
{
    const unsigned long limbs_to_extract = (num_bytes>>FLINT_LG_BYTES_PER_LIMB);
    const unsigned long extra_bytes_to_extract = num_bytes 
                                    - (limbs_to_extract<<FLINT_LG_BYTES_PER_LIMB);
                                    
    unsigned long next_limb;
    unsigned long temp = 0;
    
    // the limb we are up to in the array and output respectively
    unsigned long coeff_limb = limb_start;
    unsigned long output_limb = 0;

    unsigned long shift_1, shift_2;
    
    unsigned long i, j;
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS_PER_LIMB - shift_1;
    
    temp = array[coeff_limb];
    coeff_limb++;
    while (output_limb < limbs_to_extract)
    {
       next_limb = shift_r(temp,shift_1);
       temp = array[coeff_limb];
       coeff_limb++;
       next_limb += shift_l(temp,shift_2);
       output[output_limb] = next_limb;
       output_limb++;
    }
    if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
    {
       next_limb = shift_r(temp,shift_1);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    } else
    {
       next_limb = shift_r(temp,shift_1);
       temp = array[coeff_limb];
       next_limb += shift_l(temp,shift_2);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    }
}

/*
Assumes array points to a large coefficient which contains bundle smaller signed 
coefficients, each coeff_bytes long, which are to be converted out to mpz_t's. temp
is an mpz_t provided by the caller, to be used for scratch space. It need not be 
set to zero.
*/
void ssmul_convert_out_bytes_signed(mpz_t* output, mp_limb_t* array, 
               unsigned long bundle, unsigned long coeff_bytes, mpz_t temp)
{
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   mp_limb_t* temp_array = (mp_limb_t*) limb_alloc(((coeff_bytes)>>FLINT_LG_BYTES_PER_LIMB)+1,0);
   
   unsigned long limb_upto = 0;
   unsigned long byte_offset = 0;
   
   unsigned long sign;
   unsigned long borrow = 0;
   
   for (unsigned long i = 0; i < bundle; i++)
   {
       clear_limbs(temp_array,((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB)+1);
       sign = ssmul_unpack_signed_bytes(temp_array, array, limb_upto, 
                                             byte_offset, coeff_bytes);
       mpz_import(temp,((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB)+1, -1, sizeof(mp_limb_t), 0, 0, temp_array);
       if (sign)
       {
          mpz_neg(temp,temp);
          mpz_sub_ui(temp,temp,1);
       }
       if (borrow) mpz_add_ui(temp,temp,1);
       mpz_add(output[i],output[i],temp);
       borrow = 0;
       if (sign) borrow = 1;
       limb_upto += limbs_per_coeff;
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
   }                        
   limb_release();
}

/*
Assumes array points to a large coefficient which contains bundle smaller 
coefficients, each coeff_bytes long, which are to be converted out to mpz_t's. temp
is an mpz_t provided by the caller, to be used for scratch space. It need not be 
set to zero.
*/
void ssmul_convert_out_bytes(mpz_t* output, mp_limb_t* array, 
               unsigned long bundle, unsigned long coeff_bytes, mpz_t temp)
{
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   mp_limb_t* temp_array = (mp_limb_t*) limb_alloc(((coeff_bytes)>>FLINT_LG_BYTES_PER_LIMB)+1,0);
   
   unsigned long limb_upto = 0;
   unsigned long byte_offset = 0;
   
   for (unsigned long i = 0; i < bundle; i++)
   {
       clear_limbs(temp_array,((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB)+1);
       ssmul_unpack_bytes(temp_array, array, limb_upto, 
                                             byte_offset, coeff_bytes);
       mpz_import(temp,((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB)+1, -1, sizeof(mp_limb_t), 0, 0, temp_array);
       if (mpz_sgn(output[i])) mpz_add(output[i],output[i],temp);
       else mpz_set(output[i],temp);
       
       limb_upto += limbs_per_coeff;
       
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
   }                        
   limb_release();
}

/*
Converts an mpz_t to raw coefficient format.

Only the absolute value of the input is used (the sign is ignored), and
it must lie in the range 0 <= abs(input) < 2^((n+1)B - 1).

*/
void convert_unsigned_mpz_to_raw(mp_limb_t* output, mpz_t input,
                                 unsigned long n)
{
   unsigned long size = mpz_size(input);

   copy_limbs(output, input->_mp_d, size);
   clear_limbs(output + size, n + 1 - size);
}


/*
Converts a possibly signed mpz_t to our raw coefficient format, taking into 
account the sign. Input must fit into n+1 limbs.
*/
void ssmul_convert_in_limbs(mp_limb_t* output, mpz_t input, unsigned long n)
{
   // non-portable because we access _mp_d directly
   unsigned long size = mpz_size(input);
   if (mpz_sgn(input) >= 0)
   {
      copy_limbs(output, input->_mp_d, size);
      clear_limbs(output + size, n + 1 - size);
   }
   else
   {
      negate_limbs(output, input->_mp_d, size);
      set_limbs(output + size, n + 1 - size);
   }
}


/*
Converts our raw coefficient format to an mpz_t. Input must be no longer 
than n limbs. The input should not be negative.
*/
inline void convert_raw_to_mpz(mpz_t output, mp_limb_t* input, unsigned long n)
{
   mpz_import(output, n, -1, sizeof(mp_limb_t), 0, 0, input);
}

/*
This is similar to convert_raw_to_mpz. However, the input MUST be in the
normalised range 0 <= x < p, and if it is >= (p-1)/2, it will be converted
to the appropriate negative integer (i.e. congruent mod p to the input).

NOTE: the data in the input buffer will be destroyed.
 */
void ssmul_convert_out_limbs_signed(mpz_t output, mp_limb_t* input, unsigned long n)
{
   if (input[n] || (input[n-1] >> (FLINT_BITS_PER_LIMB - 1)))
   {
      // If the coefficient is at least (p-1)/2, then it actually represents
      // a negative number so compute p - coefficient
      negate_limbs(input, input, n);
      //mpn_add_1(input, input, n, 1);
      convert_raw_to_mpz(output, input, n);
      mpz_neg(output, output);
   }
   else
   {
      // positive case
      convert_raw_to_mpz(output, input, n);
   }
}

/*
An unsigned version of ssmul_convert_out_limbs_signed.
*/
inline void ssmul_convert_out_limbs(mpz_t output, mp_limb_t* input, unsigned long n)
{
   convert_raw_to_mpz(output, input, n);
}

/*
Splits a large coefficient into split pieces and split-1 zeroes. The chunks
are each allocated n+1 limbs of space apiece, are stored in array. Each 
chunk of the original coefficient is taken to be input_limbs long. 
*/
inline void ssmul_convert_in_split(mp_limb_t** array, mpz_t data, 
            unsigned long input_limbs, unsigned long split, unsigned long n)
{
   unsigned long k, l, m;
   const unsigned long size = mpz_size(data);
   
   for (k = 0, l = 0; k < split-1, l + input_limbs < size; k++, l+= input_limbs)
   {
       clear_limbs(array[k],n+1);
       for (m = 0; m < n; m += 8) FLINT_PREFETCH(array[k+1], m);
       copy_limbs(array[k],(data->_mp_d) + l, input_limbs);       
   }
   clear_limbs(array[k], n + 1);
   copy_limbs(array[k], (data->_mp_d) + l, size - l);
   k++;
   for (;k < 2*split - 1; k++) clear_limbs(array[k], n + 1);
}           

/*
Splits a large signed coefficient into split pieces and split-1 zeroes. The 
chunks are each allocated n+1 limbs of space apiece, are stored in array. 
Each chunk of the original coefficient is taken to be input_limbs long. 
*/
inline void ssmul_convert_in_split_signed(mp_limb_t** array, mpz_t data, 
            unsigned long input_limbs, unsigned long split, unsigned long n)
{
   unsigned long k, l, m;
   const unsigned long size = mpz_size(data);
   
   if (mpz_sgn(data) >= 0) // positive coefficient
   {
      for (k = 0, l = 0; k < split-1, l + input_limbs < size; k++, l+= input_limbs)
      {
         clear_limbs(array[k], n + 1);
         for (m = 0; m < n; m += 8) FLINT_PREFETCH(array[k+1], m);
         copy_limbs(array[k], (data->_mp_d) + l, input_limbs);       
      }
      clear_limbs(array[k], n + 1);
      copy_limbs(array[k], (data->_mp_d) + l, size - l);
      k++;
      for (;k < 2*split - 1; k++) clear_limbs(array[k], n + 1);
   } else // negative coefficient
   {
      for (k = 0, l = 0; k < split - 1, l + input_limbs < size; k++, l+= input_limbs)
      {
         for (m = 0; m < n; m += 8) FLINT_PREFETCH(array[k+1], m);
         negate_limbs(array[k], (data->_mp_d) + l, input_limbs);  
         set_limbs(array[k] + input_limbs, n + 1 - input_limbs);     
      }
      negate_limbs(array[k], (data->_mp_d) + l, size - l);
      set_limbs(array[k] + size - l, n + 1 - (size - l));
      k++;
      for (;k < 2*split - 1; k++) clear_limbs(array[k], n + 1);
   }
}           

inline void ssmul_convert_in(mp_limb_t** array, mpz_t* data, unsigned long length,
             int sign, int bit_pack, unsigned long split, unsigned long bundle,
             unsigned long orig_length, unsigned long n, unsigned long input_limbs,
             unsigned long skip_limbs, unsigned long orig_output_bits, 
             unsigned long coeffs_per_limb)
{
   unsigned long i, j, k, m;
   
   // convert data for first FFT   
   if (bit_pack)
   {
      for (k = 0, i = 0; i < length/2; i++, k+=bundle)
      {
         clear_limbs(array[i],n+1);
         // prefetch entire bundled coefficient
         for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j);
         // convert a bundle of coefficients
         if (k + bundle <= orig_length)
         {
            if (!sign) ssmul_convert_in_bits(array[i], data+k, n, bundle, orig_output_bits, coeffs_per_limb);
            else ssmul_convert_in_bits_signed(array[i], data+k, n, bundle, orig_output_bits, coeffs_per_limb, 1);
         } else if (k < orig_length)
         {
            if (!sign) ssmul_convert_in_bits(array[i], data+k, n, orig_length - k, orig_output_bits, coeffs_per_limb);
            else ssmul_convert_in_bits_signed(array[i], data+k, n, orig_length - k, orig_output_bits, coeffs_per_limb, 1);
         }
      }       
   } else if (split != 1) // we need to split coefficients
   {
      // for each orig coefficient
      for (i = 0, j = 0; i < orig_length; i++, j+=(2*split-1))
      {
          if (!sign) ssmul_convert_in_split(array+j, data[i], input_limbs, split, n);
          else ssmul_convert_in_split_signed(array+j, data[i], input_limbs, split, n);      
      }
      // clear remainder of fft input coeffs
      for (;j<length/2;j++)
      {
         for (m = 0; m < n; m += 8) FLINT_PREFETCH(array[j+1], m);
         clear_limbs(array[j],n+1);
      }
   } else if (bundle == 1)
   {
      for (i = 0; i < orig_length; i++)
      {
         for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j);
         ssmul_convert_in_limbs(array[i], data[i], n);
      }
      for (; i < length/2; i++)
      {
         clear_limbs(array[i],n+1);
      }
   } else
   {
      for (k = 0, i = 0; i < length/2; i++, k+=bundle)
      {
         clear_limbs(array[i],n+1);
         // prefetch entire bundled coefficient
         for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j);
         // convert a bundle of coefficients
         if (k + bundle <= orig_length)
         {
            ssmul_convert_in_bytes(array[i], data+k, bundle, skip_limbs, n+1);
         } else if (k < orig_length)
         {
            ssmul_convert_in_bytes(array[i], data+k, orig_length-k, skip_limbs, n+1);
         }
      }
   }
}

/*
Divide by appropriate normalising power of 2 (I'm assuming
here that the transform length will always be less than 2^B)
and return 1 if the coefficient is negative
*/
inline int scale_coefficient(mp_limb_t** array, unsigned long i, 
            unsigned long log_length, unsigned long n, int sign)
{
      unsigned long j;
      
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j);
      rotate_right_bits(array[i], array[i], log_length+1, n);
      reduce_mod_p_exact(array[i], n);
      if (sign)
      {
         if (array[i][n] || (array[i][n-1] >> (FLINT_BITS_PER_LIMB - 1)))
         {
            // If the coefficient is at least (p-1)/2, then it actually represents
            // a negative number. Convert to a two's complement negative number
            mpn_sub_1(array[i], array[i], n, 1);
            return 1;
         } else return 0;
      }
      return 0;
}

inline void ssmul_convert_out(mpz_t* res, mp_limb_t** array, 
            unsigned long orig_length, unsigned long length2, unsigned long length, 
            unsigned long log_length, unsigned long n)
{
   unsigned long i, j;
   
   for (i = 0; i < length-1; i++)
   { 
      scale_coefficient(array, i, log_length, n, 0);
      
      // prefetch entire next coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j);
      
      if (i < orig_length + length2 - 1) 
      {
         ssmul_convert_out_limbs(res[i], array[i], n);
      }
   }   
}

inline void ssmul_convert_out_signed(mpz_t* res, mp_limb_t** array, 
            unsigned long orig_length, unsigned long length2,  unsigned long length, 
            unsigned long log_length, unsigned long n)
{
   unsigned long i, j;
   
   for (i = 0; i < length-1; i++)
   { 
      scale_coefficient(array, i, log_length, n, 1);
      
      // prefetch entire next coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j);
      
      if (i < orig_length + length2 - 1) 
      {
         ssmul_convert_out_limbs_signed(res[i], array[i], n);
      }
   }   
}

inline void ssmul_bundle_convert_out(mpz_t* res, mp_limb_t** array, 
            unsigned long orig_length, unsigned long length2, unsigned long length, 
            unsigned long log_length, unsigned long n, unsigned long skip_limbs, 
            unsigned long bundle)
{
   unsigned long i, j, k;
   
   mpz_t temp;
   mpz_init(temp);
   
   for (k=0, i = 0; i < length-1; i++, k+=bundle)
   { 
      scale_coefficient(array, i, log_length, n, 0);
      
      // prefetch entire next coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j); 
      
      if (k + 2*bundle <= orig_length + length2)
      {
         ssmul_convert_out_bytes(res + k, array[i], 2*bundle-1, skip_limbs, temp);
      } else if (k < orig_length + length2 - 1)
      {
         ssmul_convert_out_bytes(res + k, array[i], orig_length + length2 - 1 - k, skip_limbs, temp);
      }
   }
}

inline void ssmul_bundle_convert_out_signed(mpz_t* res, mp_limb_t** array, 
            unsigned long orig_length, unsigned long length2, unsigned long length, 
            unsigned long log_length, unsigned long n, unsigned long skip_limbs, 
            unsigned long bundle)
{
   unsigned long i, j, k;
   
   mpz_t temp;
   mpz_init(temp);
   
   for (k=0, i = 0; i < length-1; i++, k+=bundle)
   { 
      scale_coefficient(array, i, log_length, n, 1);
      
      // prefetch entire next coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j); 
      
      if (k + 2*bundle <= orig_length + length2)
      {
         ssmul_convert_out_bytes_signed(res + k, array[i], 2*bundle-1, skip_limbs, temp);
      } else if (k < orig_length + length2 - 1)
      {
         ssmul_convert_out_bytes_signed(res + k, array[i], orig_length + length2 - 1 - k, skip_limbs, temp);
      }
   }
   
   mpz_clear(temp);
}

inline void ssmul_split_convert_out(mpz_t* res, mp_limb_t** array, 
            unsigned long input_limbs, unsigned long orig_length, unsigned long length2, 
            unsigned long length, unsigned long log_length, 
            unsigned long n, unsigned long orig_limbs, unsigned long split)
{
   
   unsigned long output_limbs;
   
   mp_limb_t* temp_coeff;
   
   unsigned long h, i, j, l, m;
   
   // each output coefficient can't have more than 2*orig_limbs+2 limbs
   // even on a 32 bit machine
   // but the final addition has n+1 limbs and might be adding to the final 
   // limb of this, so just to be on the safe side....
   output_limbs = 2*orig_limbs+n+2;
   // allocate space for a full output coefficient
   temp_coeff = (mp_limb_t*) limb_alloc(output_limbs, 0);
   clear_limbs(temp_coeff,output_limbs);
   // l+m gives which partial output (fft) coefficient we are up to
   l = 0;
   // m gives which partial output coefficient we are up to within a full output coeff
   m = 0;
   // which full output coefficient we are up to
   h = 0; 
         
   for (i = 0; i < length-1; i++)
   { 
      // divide by appropriate normalising power of 2 (I'm assuming
      // here that the transform length will always be less than 2^B)
      // and set coeff_sign if the coefficient is negative
      
      scale_coefficient(array, i, log_length, n, 0);
      
      // prefetch entire next coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j); 
      
      // convert all the coefficients out to mpz_t's, adding coefficients 
      // together where necessitated by the splitting trick.
      if (i < (orig_length + length2 - 1)*(2*split-1))
      {
         if (l == 2*split - 1)
         {
            mpz_import(res[h], output_limbs, -1, sizeof(mp_limb_t), 0, 0, temp_coeff);
               
            h++;
            l = 0;
            clear_limbs(temp_coeff, output_limbs);
            m+=(2*split-1);
         } 
            
         mpn_add(temp_coeff+input_limbs*l, temp_coeff+input_limbs*l, output_limbs-input_limbs*l, array[m+l], n); 
         l++;
      }
   } 
   
   // write out final coefficient
   mpz_import(res[h], output_limbs, -1, sizeof(mp_limb_t), 0, 0, temp_coeff);
       
   // free allocated space
   limb_release(); 
}

inline void ssmul_split_convert_out_signed(mpz_t* res, mp_limb_t** array, 
            unsigned long input_limbs, unsigned long orig_length, unsigned long length2, 
            unsigned long length, unsigned long log_length, 
            unsigned long n, unsigned long orig_limbs, unsigned long split)
{
   
   unsigned long output_limbs;
   int coeff_sign;
   
   mp_limb_t* temp_coeff;
   
   unsigned long h, i, j, l, m;
   
   // each output coefficient can't have more than 2*orig_limbs+2 limbs
   // even on a 32 bit machine
   // but the final addition has n+1 limbs and might be adding to the final 
   // limb of this, so just to be on the safe side....
   output_limbs = 2*orig_limbs+n+2;
   // allocate space for a full output coefficient
   temp_coeff = (mp_limb_t*) limb_alloc(output_limbs, 0);
   clear_limbs(temp_coeff,output_limbs);
   // l+m gives which partial output (fft) coefficient we are up to
   l = 0;
   // m gives which partial output coefficient we are up to within a full output coeff
   m = 0;
   // which full output coefficient we are up to
   h = 0; 
         
   for (i = 0; i < length-1; i++)
   { 
      // divide by appropriate normalising power of 2 (I'm assuming
      // here that the transform length will always be less than 2^B)
      // and set coeff_sign if the coefficient is negative
      
      coeff_sign = scale_coefficient(array, i, log_length, n, 1);
      
      // prefetch entire next coefficient
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array[i+1], j); 
      
      // convert all the coefficients out to mpz_t's, adding coefficients 
      // together where necessitated by the splitting trick.
      if (i < (orig_length + length2 - 1)*(2*split-1))
      {
         if (l == 2*split - 1)
         {
            if (temp_coeff[output_limbs-1]>>(FLINT_BITS_PER_LIMB-1)) 
            {
               negate_limbs(temp_coeff,temp_coeff,output_limbs);
               mpz_import(res[h], output_limbs, -1, sizeof(mp_limb_t), 0, 0, temp_coeff);
               mpz_neg(res[h],res[h]);
            } else
            {
               mpz_import(res[h], output_limbs, -1, sizeof(mp_limb_t), 0, 0, temp_coeff);
            }
               
            h++;
            l = 0;
            clear_limbs(temp_coeff,output_limbs);
            m+=(2*split-1);
         } 
            
         if (coeff_sign) 
         {
            mpn_add(temp_coeff+input_limbs*l, temp_coeff+input_limbs*l, output_limbs-input_limbs*l, array[m+l], n); 
            mpn_sub_1(temp_coeff+input_limbs*l+n,temp_coeff+input_limbs*l+n,output_limbs-input_limbs*l-n,1);
         } else
         {
            mpn_add(temp_coeff+input_limbs*l, temp_coeff+input_limbs*l, output_limbs-input_limbs*l, array[m+l], n); 
         }
         l++;
      }
   } 
   
   // write out final coefficient
   if (temp_coeff[output_limbs-1]>>(FLINT_BITS_PER_LIMB-1)) 
   {
      negate_limbs(temp_coeff,temp_coeff,output_limbs);
      mpz_import(res[h], output_limbs, -1, sizeof(mp_limb_t), 0, 0, temp_coeff);
      mpz_neg(res[h],res[h]);
   } else
   {
      mpz_import(res[h], output_limbs, -1, sizeof(mp_limb_t), 0, 0, temp_coeff);
   }
       
   // free allocated space
   limb_release(); 
}

/*==============================================================================

Arithmetic Support Functions for the FFT and IFFT. 

Everywhere we are working modulo a fermat number p = 2^(nB) + 1, where n >= 1.

All coefficients are stored in arrays of limbs of length n+1 (least significant
limb first). The limbs represent a single signed integer in 2's complement;
the high bit of x[n] is the sign bit.

In the FFT, we try to do reductions mod p as rarely as possible; we just let
bits build up in the extra limb as long as we can get away with it. To ensure
correctness, each function's documentation describes exactly how much of the
last limb can get chewed up.

===============================================================================*/

/*
Divides a coefficient by 2^s mod p. The shift s MUST be in the range 0 < s < B;
i.e. this is a rotation by a fraction of a limb.

src and dest can be the same buffer. If they're not, they must be disjoint.

Total bitlength of coefficient is not increased. That is, if
    0 <= abs(src) < 2^(nB + k),
then also
    0 <= abs(dest) < 2^(nB + k).

*/
void rotate_right_bits(mp_limb_t* dest, mp_limb_t* src,
                       unsigned long s, unsigned long n)
{
   mp_limb_t high = src[n];
   mp_limb_t overflow = mpn_rshift(dest, src, n + 1, s);

   if ((mp_limb_signed_t)(high) < 0)
      // If src was negative, we have to add extra 1's at the top to
      // make the shifted version negative (because mpn_rshift doesn't
      // respect sign)
      dest[n] += (mp_limb_signed_t)(-1L) << (FLINT_BITS_PER_LIMB - s);
   
   // rotate the overflow back around to the top limb
   mpn_sub_1(dest + n - 1, dest + n - 1, 2, overflow);
}


/*
Reduces the input into the canonical range 0 <= x < p. Result is inplace.

In most cases, this involves simply distributing the top limb down to the
bottom.

In rare cases, it will need to propagate a carry to the next limb. (How rare
depends on how full the top limb is. The probability will be lower for
smaller FFT transform lengths.)

In extremely rare cases it will need to propagate the carry further, possibly
doing a whole pass over the data.

In pathologically rare cases it will perform TWO passes over the data. For
example if the input is exactly p-1, it will first subtract p to obtain -1,
and then add p back to get p-1.

*/
void reduce_mod_p_exact(mp_limb_t* x, unsigned long n)
{
   mp_limb_t hi = x[n];

   if ((mp_limb_signed_t) hi < 0)
   {
      // If top limb (hi) is negative, we add -hi multiples of p
      x[n] = 0;
      mpn_add_1(x, x, n + 1, -hi);

      // If the result is >= p (very unlikely)...
      if (x[n] && x[0])
      {
         // ... need to subtract off p.
         x[n] = 0;
         x[0]--;
      }
   }
   else
   {
      // If top limb (hi) is non-negative, we subtract hi multiples of p
      x[n] = 0;
      mpn_sub_1(x, x, n + 1, hi);

      // If the result is negative (very unlikely)...
      if (x[n])
      {
         // ... need to add back p.
         x[n] = 0;
         mpn_add_1(x, x, n + 1, 1);
      }
   }
}


/*
Computes dest = src * 2^(Bs) mod p.

i.e. shifts left by a whole number of limbs.

MUST have 0 <= s < n.

dest and src may NOT overlap.

The top limb is distributed down, so the result is guaranteed to satisfy
    -2^B + 1 <= dest <= 2^(nB) + 2^B - 2
In other words it doesn't use more than one bit of the top limb.

*/
void rotate_mod_p_limbs(mp_limb_t* dest, mp_limb_t* src,
                        unsigned long s, unsigned long n)
{
   // put negative of high limbs of src into low limbs of dest
   negate_limbs(dest, src + n - s, s + 1);
   mp_limb_t carry = dest[s];

   // put low limbs of src into high limbs of dest
   copy_limbs(dest + s, src, n - s);
   dest[n] = 0;

   // propagate carry
   signed_add_1(dest + s, n - s + 1, carry);
}


/*
Computes dest = src * 2^s mod p.

i.e. shifts left by s bits.

MUST have 0 <= s < Bn.

dest and src may NOT overlap.

The top limb is distributed down, so the result is guaranteed to satisfy
    -2^B + 1 <= dest <= 2^(nB) + 2^B - 2
In other words it doesn't use more than one bit of the top limb.

*/
void rotate_mod_p_bits(mp_limb_t* dest, mp_limb_t* src,
                       unsigned long s, unsigned long n)
{
   // split rotation into a whole number of limbs and leftover number of bits

   // dirty: word size must be a power of two
   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
    
   if (bits == 0)
   {
      // need special case here for bits == 0 because rotate_right_bits
      // does not allow shift count to be zero
      rotate_mod_p_limbs(dest, src, s, n);
      return;
   }

   // Instead of shifting left by s limbs and bits bits, we will
   // shift left by s+1 limbs and RIGHT by (B - bits) bits.
   s++;

   // put the negative of the high limbs of src into the low limbs of dest
   negate_limbs(dest, src + n - s, s + 1);
   mp_limb_t carry = dest[s];
   
   // put the low limbs of src into the high limbs of dest
   if (n != s)
      copy_limbs(dest + s, src, n - s);
   dest[n] = 0;
   
   // propagate carry
   signed_add_1(dest + s, n - s + 1, carry);
   
   // do the fractional limb rotation
   rotate_right_bits(dest, dest, FLINT_BITS_PER_LIMB - bits, n);
}


/*============================================================================

FFT butterflies.

We have separate functions for doing butterflies with rotations by an
arbitrary number of bits and by a whole number of limbs; this saves on some
arithmetic and branching in the deeper recursion levels.

All of these functions operate on raw coefficient buffers (n+1 limbs as
described above). However, they accept POINTERS to mp_limb_t*. The reason is
that it might not be efficient to store the output in the same buffers as
the input, so the functions are permitted to store one of the results in the
scratch buffer, and then permute the pointers so that *a is the first output
and *b is the second output (and *scratch points to whatever the third
remaining buffer is).

============================================================================*/

/* 
Swaps arrays of limbs efficiently

*/
inline void swap_limb_ptrs(mp_limb_t** x, mp_limb_t** y)
{
    mp_limb_t* temp;
    temp = *x;
    *x = *y;
    *y = temp;
}


/*
Computes  (a, b) -> (a + b, 2^(Bs) (a - b)).

i.e. performs a forward FFT butterfly with a rotation by a whole number
of limbs.

MUST have 0 < s < n.

Suppose the inputs satisfy
   0 <= abs(a) < 2^(nB + k),
   0 <= abs(b) < 2^(nB + k),
for some integer 0 <= k < B-1.

Then the outputs have the following guarantees:
   0 <= abs(a) < 2^(nB + k + 1),
   0 <= abs(b) < 2^(nB + 1).

In other words, a is at most one bit bigger than either of the inputs, and b
does not use more than one bit of the last limb.

*/
void fft_butterfly_limbs(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                         unsigned long s, unsigned long n)
{
   // get low limbs of a - b into high limbs of output
   (*scratch)[n] = -mpn_sub_n(*scratch + s, *a, *b, n - s);
    
   // get high limbs of b - a into low limbs of output
   mp_limb_t overflow =
         (*b)[n] - (*a)[n] - mpn_sub_n(*scratch, *b + n - s, *a + n - s, s);
   signed_add_1(*scratch + s, n + 1 - s, overflow);
    
   // a = a + b
   mpn_add_n(*a, *a, *b, n+1);
    
   // swap rotated version of a - b back into b
   swap_limb_ptrs(b, scratch);
}


/*
Computes  (a, b) -> (a + b, 2^s (a - b)).

i.e. performs a forward FFT butterfly with a rotation by an arbitrary number
of bits.

MUST have 0 < s < Bn.

Output guarantees are the same as for fft_butterfly_limbs.

*/
void fft_butterfly_bits(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                        unsigned long s, unsigned long n)
{
   // split rotation into a whole number of limbs and leftover number of bits

   // dirty: word size must be a power of two
   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
    
   if (bits == 0)
   {
      // need special case here for bits == 0 because rotate_right_bits
      // does not allow shift count to be zero
      fft_butterfly_limbs(a, b, scratch, s, n);
      return;
   }
   
   // Instead of shifting left by s limbs and bits bits, we will
   // shift left by s+1 limbs and RIGHT by (B - bits) bits.
   s++;
   
   // First do the subtraction and whole-limb rotation in a single pass.
   // We use pretty much the same code sequence as in fft_butterfly_limbs,
   // except we also need to cover the case s == n.
   if (s == n)
      (*scratch)[n] = 0;
   else
      (*scratch)[n] = -mpn_sub_n(*scratch + s, *a, *b, n - s);
    
   mp_limb_t overflow =
           (*b)[n] - (*a)[n] - mpn_sub_n(*scratch, *b + n - s, *a + n - s, s);
   signed_add_1(*scratch + s, n + 1 - s, overflow);
   
   // a = a + b
   mpn_add_n(*a, *a, *b, n + 1);
   
   // Now do the fractional-limb rotation
   rotate_right_bits(*b, *scratch, FLINT_BITS_PER_LIMB - bits, n);
}


/*
Computes  (a, b) -> (a - 2^(Bs) b, a + 2^(Bs) b).

i.e. performs an inverse FFT butterfly with a rotation by a whole number
of limbs.

MUST have 0 < s < n.

The output guarantees for this function are very tight: the high limbs of
both outputs differ by at most 2 from the high limb of the input a.

*/
void ifft_butterfly_limbs(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                          unsigned long s, unsigned long n)
{
   // add low limbs of b to high limbs of a
   (*scratch)[n] = (*a)[n] + mpn_add_n(*scratch + s, *a + s, *b, n - s);
    
   // subtract high limbs of b from low limbs of a
   mp_limb_t overflow = -(*b)[n] - mpn_sub_n(*scratch, *a, *b + n - s, s);
   signed_add_1(*scratch + s, n + 1 - s, overflow);
   
   // subtract low limbs of b from high limbs of a
   (*a)[n] -= mpn_sub_n(*a + s, *a + s, *b, n - s);
   
   // add low limbs of a to high limbs of b
   overflow = (*b)[n] + mpn_add_n(*a, *a, *b + n - s, s);
   signed_add_1(*a + s, n + 1 - s, overflow);
   
   // swap answer into b
   swap_limb_ptrs(b, scratch);
}


/*
Computes  (a, b) -> (a - 2^s b, a + 2^s b).

i.e. performs an inverse FFT butterfly with a rotation by an arbitrary number
of bits.

MUST have 0 < s < nB.

Output guarantees are the same as for ifft_butterfly_limbs.

*/
void ifft_butterfly_bits(mp_limb_t** a, mp_limb_t** b, mp_limb_t** scratch,
                         unsigned long s, unsigned long n)
{
   // split shift into a whole number of limbs and leftover number of bits

   // dirty: word size must be a power of two
   unsigned long bits = s & (FLINT_BITS_PER_LIMB - 1);
   s /= FLINT_BITS_PER_LIMB;
   
   if (bits == 0)
   {
      // need special case here for bits == 0 because rotate_right_bits
      // does not allow shift count to be zero
      ifft_butterfly_limbs(a, b, scratch, s, n);
      return;
   }
   
   // Instead of shifting left by s limbs and bits bits, we will
   // shift left by s+1 limbs and RIGHT by (B - bits) bits.
   s++;
   
   // First do the fractional-limb rotation.
   rotate_right_bits(*scratch, *b, FLINT_BITS_PER_LIMB - bits, n);
   
   // Now basically use the same code sequence as in ifft_butterfly_limbs,
   // except we also need to cover the case s == n.
   
   if (n == s)
      (*b)[n] = (*a)[n];
   else
      (*b)[n] = (*a)[n] + mpn_add_n(*b + s, *a + s, *scratch, n - s);
   
   mp_limb_t overflow =
               -(*scratch)[n] - mpn_sub_n(*b, *a, *scratch + n - s, s);
   signed_add_1(*b + s, n + 1 - s, overflow);
   
   if (n != s)
      (*a)[n] -= mpn_sub_n(*a + s, *a + s, *scratch, n - s);
   
   overflow = (*scratch)[n] + mpn_add_n(*a, *a, *scratch + n - s, s);
   signed_add_1(*a + s, n + 1 - s, overflow);
}

/*============================================================================

Main FFT functions. These functions call the butterflies in the right order to
perform FFTs and IFFTs.

=============================================================================*/

/*
Basecase of the FFT. Once the FFT has been broken up into small enough pieces
this function actually performs the FFT recursively for each of those pieces.

*/

void fft_recursive(mp_limb_t** start, unsigned long skip,
                   unsigned long start_r, unsigned long skip_r,
                   unsigned long depth, mp_limb_t** scratch,
                   unsigned long n, int first)
{
   int fractional_limb_shift = (skip_r | start_r) & (FLINT_BITS_PER_LIMB - 1);

   unsigned long half = skip << (depth - 1);
   mp_limb_t** middle = start + half;
   unsigned long offset, r;

   if (first)
   {
      if (fractional_limb_shift)
      {
         for (offset = 0, r = start_r; offset < half;
              offset += skip, r += skip_r)
            rotate_mod_p_bits(middle[offset], start[offset], r, n);
      }
      else
      {
         unsigned long skip_r_limbs = skip_r / FLINT_BITS_PER_LIMB;

         for (offset = 0, r = start_r / FLINT_BITS_PER_LIMB; offset < half;
              offset += skip, r += skip_r_limbs)
         {
            rotate_mod_p_limbs(middle[offset], start[offset], r, n);
         }
      }

      if (depth == 1)
         return;
   }
   else
   {
      if (start_r == 0)
      {
         // scratch = a - b
         mpn_sub_n(*scratch, *start, *middle, n + 1);
         
         // a += b
         mpn_add_n(*start, *start, *middle, n + 1);
         
         swap_limb_ptrs(middle, scratch);
      }
      else
         fft_butterfly_bits(start, middle, scratch, start_r, n);

      if (depth == 1)
         return;
      
      if (fractional_limb_shift)
      {
         for (offset = skip, r = start_r + skip_r; offset < half;
              offset += skip, r += skip_r)
         {
             fft_butterfly_bits(start + offset, middle + offset, scratch, r, n);
         } 

      }
      else
      {
         unsigned long skip_r_limbs = skip_r / FLINT_BITS_PER_LIMB;
         unsigned long start_r_limbs = start_r / FLINT_BITS_PER_LIMB;
         
         for (offset = skip, r = start_r_limbs + skip_r_limbs; offset < half;
              offset += skip, r += skip_r_limbs)
         {
             fft_butterfly_limbs(start + offset, middle + offset,
                                 scratch, r, n);
         }
      }
   }

   fft_recursive(start, skip, start_r << 1, skip_r << 1, depth - 1,
                 scratch, n, 0);
   fft_recursive(middle, skip, start_r << 1, skip_r << 1, depth - 1,
                 scratch, n, 0);
}

void* fft_loop(void* fft_p)
{
   unsigned long i;
   
   fft_t fft_params = *((fft_t*) fft_p);
   
   unsigned long length = fft_params.length;
   
   mp_limb_t** start = fft_params.start;
   unsigned long skip = fft_params.skip;
   unsigned long next_skip = fft_params.next_skip;
   
   unsigned long start_r = fft_params.start_r; 
   unsigned long next_skip_r = fft_params.next_skip_r;
   
   unsigned long depth = fft_params.depth;
   unsigned long n = fft_params.n;
   mp_limb_t** scratch = fft_params.scratch+*fft_params.thread;
   
   int first = fft_params.first;
   
   for (i = 0; i < length; i++, start += next_skip)
   {
      fft_recursive(start, skip, start_r, next_skip_r, depth, scratch, n, first);
   }
}

void* fft_loop2(void* fft_p)
{
   unsigned long i;
   
   fft_t fft_params = *((fft_t*) fft_p);
   
   unsigned long length = fft_params.length;
   
   mp_limb_t** start = fft_params.start;
   unsigned long skip = fft_params.skip;
   unsigned long next_skip = fft_params.next_skip;
   
   unsigned long start_r = fft_params.start_r; 
   unsigned long next_skip_r = fft_params.next_skip_r;
   unsigned long skip_r = fft_params.skip_r;
   
   unsigned long depth = fft_params.depth;
   
   unsigned long n = fft_params.n;
   mp_limb_t** scratch = fft_params.scratch+*fft_params.thread;
   
   int first = fft_params.first;

   for (i = 0; i < length; i++, start += skip, start_r += skip_r)
   {
      fft_recursive(start, next_skip, start_r, next_skip_r, depth, scratch, n, first);
   }
}

/*
transform of length 2^depth

coefficients are at start, start + skip, ... start + (length-1)*skip

corresponding values of r are:
 start_r, start_r + skip_r, ... start_r + (length/2 - 1)*skip_r

the idea of this routine is to get the transforms down to fit in cache
cache if possible, then switch to a simple recursive algorithm
 
*/
void fft_main(mp_limb_t** start, unsigned long skip,
              unsigned long start_r, unsigned long skip_r,
              unsigned long depth, mp_limb_t** scratch,
              unsigned long n, int first, int crossover = -1)
{
   pthread_t thread_arr[THREADS];
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
   
   unsigned long thread_num[16] = {0,16,32,48,64,80,96,112,128,144,160,176,192,208,224,240};
   unsigned long i;
   
   unsigned long threads;
   unsigned long lg_threads;
   
   fft_t fft_params[THREADS];
   
   if (crossover == -1)
   {
      // work out crossover = optimal depth to switch over to plain
      // recursive algorithm
      unsigned long cache_length =
                        2*FLINT_CACHE_SIZE / ((n+1) * sizeof(mp_limb_t)) - 1;
      for (crossover = 0; cache_length > 0; cache_length >>= 1, crossover++);
      if (crossover == 0)
         crossover = 1;
         
      if (depth <= crossover)
      {
         fft_recursive(start, skip, start_r, skip_r, depth, scratch, n, first);
         return;
      }
   }

   // Factor FFT into two pieces.
   // Write depth = depth1 + depth2.
   // First piece does the first depth1 layers.
   // Second piece does the last depth2 layers.
   
   // In other words we do
   // length2 transforms of length length1, followed by
   // length1 transforms of length length2,
   // where length1 * length2 = 2^depth.

   unsigned long depth2 = crossover;
   unsigned long depth1 = depth - depth2;

   unsigned long length1 = 1 << depth1;
   unsigned long length2 = 1 << depth2;

   mp_limb_t** next_start = start;
   unsigned long next_skip = skip << depth2;
   unsigned long next_start_r = start_r;
   unsigned long next_skip_r = skip_r << depth2;

   if (depth1 <= crossover)
   {
#if USE_THREADS
      threads = THREADS;
      lg_threads = LG_THREADS;
      while (threads > length2) 
      {
         lg_threads--;
         threads >>= 1;
      }
      
      unsigned long new_length2 = ((length2-1)>>lg_threads)+1;
      for (i = 0; i < threads; i++)
      {
         fft_params[i].next_skip = next_skip;
         fft_params[i].skip = skip;
         fft_params[i].next_skip_r = next_skip_r;
         fft_params[i].skip_r = skip_r;
         fft_params[i].depth = depth1;
         fft_params[i].scratch = scratch;
         fft_params[i].n = n;
         fft_params[i].first = 0;
      }
      for (i = 0; i < threads-1; i++)
      {
         fft_params[i].length = new_length2;
         fft_params[i].start = next_start + i*new_length2*skip;
         fft_params[i].start_r = next_start_r + i*new_length2*skip_r;
      }
      fft_params[i].length = length2 - new_length2*i;
      fft_params[i].start = next_start + i*new_length2*skip;
      fft_params[i].start_r = next_start_r + i*new_length2*skip_r;
      
      for (i = 0; i < threads; i++)
      {
         fft_params[i].thread = thread_num+i;
         pthread_create(thread_arr+i, &attr, fft_loop2, (void*) (fft_params+i));
             
      }
      for (unsigned long i = 0; i < threads; i++)
      {
         pthread_join(thread_arr[i], NULL);
      } 
      
#else
      for (unsigned long i = 0; i < length2; i++, next_start += skip, 
              next_start_r += skip_r)
      {
         fft_recursive(next_start, skip, next_start_r, skip_r, depth1, scratch, n, 0);
      }    
#endif
   } else
   {
      for (unsigned long i = 0; i < length2; i++, next_start += skip,
          next_start_r += skip_r)
      {
          fft_main(next_start, next_skip, next_start_r, next_skip_r,
               depth1, scratch, n, first, crossover);
      }
   }
   
   next_start = start;
   next_start_r = start_r << depth1;
   next_skip_r = skip_r << depth1;
   
   if (depth2 <= crossover)
   {
#if USE_THREADS
      threads = THREADS;
      lg_threads = LG_THREADS;
      while (threads > length1) 
      {
         lg_threads--;
         threads >>= 1;
      }
      
      unsigned long new_length1 = ((length1-1)>>lg_threads)+1;
      for (i = 0; i < threads; i++)
      {
         fft_params[i].next_skip = next_skip;
         fft_params[i].skip = skip;
         fft_params[i].next_skip_r = next_skip_r;
         fft_params[i].start_r = next_start_r;
         fft_params[i].depth = depth2;
         fft_params[i].scratch = scratch;
         fft_params[i].n = n;
         fft_params[i].first = 0;
      }
      for (i = 0; i < threads-1; i++)
      {
         fft_params[i].length = new_length1;
         fft_params[i].start = next_start + i*new_length1*next_skip;
      }
      fft_params[i].length = length1 - new_length1*i;
      fft_params[i].start = next_start + i*new_length1*next_skip;
      
      for (i = 0; i < threads; i++)
      {
         fft_params[i].thread = thread_num+i;
         pthread_create(thread_arr+i, &attr, fft_loop, (void*) (fft_params+i));
             
      }
      for (unsigned long i = 0; i < threads; i++)
      {
         pthread_join(thread_arr[i], NULL);
      } 
      
#else
      for (unsigned long i = 0; i < length1; i++, next_start += next_skip)
      {
         fft_recursive(next_start, skip, next_start_r, skip_r, depth2, scratch, n, 0);
      }    
#endif  
   } else
   {
      for (unsigned long i = 0; i < length1; i++, next_start += next_skip)
      {
         fft_main(next_start, skip, next_start_r, next_skip_r,
               depth2, scratch, n, 0, crossover);
      }
   }
   pthread_attr_destroy(&attr);
}

/*
Basecase for the IFFT. Just the inverse of fft_recursive.

*/
void ifft_recursive(mp_limb_t** start, unsigned long skip,
                    unsigned long start_r, unsigned long skip_r,
                    unsigned long depth, mp_limb_t** scratch,
                    unsigned long n)
{
   int fractional_limb_shift = (skip_r | start_r) & (FLINT_BITS_PER_LIMB - 1);

   unsigned long half = skip << (depth - 1);
   mp_limb_t** middle = start + half;
   unsigned long offset, r;

   if (depth > 1)
   {
      ifft_recursive(start, skip, start_r << 1, skip_r << 1, depth - 1,
                     scratch, n);
      ifft_recursive(middle, skip, start_r << 1, skip_r << 1, depth - 1,
                     scratch, n);
   }

   if (start_r == 0)
   {
      // scratch = a - b
      mpn_sub_n(*scratch, *start, *middle, n + 1);
      
      // a += b
      mpn_add_n(*start, *start, *middle, n + 1);
      
      swap_limb_ptrs(middle, scratch);
   }
   else
      ifft_butterfly_bits(start, middle, scratch, 
           n*FLINT_BITS_PER_LIMB - start_r, n);
   
   if (depth == 1)
      return;
   
   if (fractional_limb_shift)
   {
      for (offset = skip, r = n*FLINT_BITS_PER_LIMB - start_r - 
           skip_r; offset < half; offset += skip, r -= skip_r)
      {
          ifft_butterfly_bits(start + offset, middle + offset, scratch, r, n);
      }
   }
   else
   {
      unsigned long skip_r_limbs = skip_r / FLINT_BITS_PER_LIMB;
      unsigned long start_r_limbs = start_r / FLINT_BITS_PER_LIMB;
      
      for (offset = skip, r = n - start_r_limbs - skip_r_limbs;
             offset < half; offset += skip, r -= skip_r_limbs)
      {
         ifft_butterfly_limbs(start + offset, middle + offset,
                              scratch, r, n);
      }
   }
}

void* ifft_loop2(void* fft_p)
{
   unsigned long i;
   
   fft_t fft_params = *((fft_t*) fft_p);
   
   unsigned long length = fft_params.length;
   
   mp_limb_t** start = fft_params.start;
   unsigned long skip = fft_params.skip;
   unsigned long next_skip = fft_params.next_skip;
   
   unsigned long start_r = fft_params.start_r; 
   unsigned long next_skip_r = fft_params.next_skip_r;
   
   unsigned long depth = fft_params.depth;
   unsigned long n = fft_params.n;
   mp_limb_t** scratch = fft_params.scratch+*fft_params.thread;
   
   for (i = 0; i < length; i++, start += next_skip)
   {
      ifft_recursive(start, skip, start_r, next_skip_r, depth, scratch, n);
   }
}

void* ifft_loop(void* fft_p)
{
   unsigned long i;
   
   fft_t fft_params = *((fft_t*) fft_p);
   
   unsigned long length = fft_params.length;
   
   mp_limb_t** start = fft_params.start;
   unsigned long skip = fft_params.skip;
   unsigned long next_skip = fft_params.next_skip;
   
   unsigned long start_r = fft_params.start_r; 
   unsigned long next_skip_r = fft_params.next_skip_r;
   unsigned long skip_r = fft_params.skip_r;
   
   unsigned long depth = fft_params.depth;
   
   unsigned long n = fft_params.n;
   mp_limb_t** scratch = fft_params.scratch+*fft_params.thread;
   
   for (i = 0; i < length; i++, start += skip, start_r += skip_r)
   {
      ifft_recursive(start, next_skip, start_r, next_skip_r, depth, scratch, n);
   }
}

/*
Main IFFT routine. It breaks the IFFT up into pieces and calls 
ifft_recursive on the pieces. It is the inverse of fft_main.

*/
void ifft_main(mp_limb_t** start, unsigned long skip,
               unsigned long start_r, unsigned long skip_r,
               unsigned long depth, mp_limb_t** scratch,
               unsigned long n, int crossover = -1)
{
   pthread_t thread_arr[THREADS];
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
   
   unsigned long thread_num[16] = {0,16,32,48,64,80,96,112,128,144,160,176,192,208,224,240};
   unsigned long i;
   
   unsigned long lg_threads;
   unsigned long threads;
   
   fft_t fft_params[THREADS];
   
   if (crossover == -1)
   {
      // work out crossover = optimal depth to switch over to plain
      // recursive algorithm
      unsigned long cache_length =
                        2*FLINT_CACHE_SIZE / ((n+1) * sizeof(mp_limb_t)) - 1;
      for (crossover = 0; cache_length > 0; cache_length >>= 1, crossover++);
      if (crossover == 0)
         crossover = 1;
      
      if (depth <= crossover)
      {
         ifft_recursive(start, skip, start_r, skip_r, depth, scratch, n);
         return;
      }
   }

   // do the forward FFT backwards

   unsigned long depth2 = crossover;
   unsigned long depth1 = depth - depth2;

   unsigned long length1 = 1 << depth1;
   unsigned long length2 = 1 << depth2;

   mp_limb_t** next_start = start;
   unsigned long next_skip = skip << depth2;
   unsigned long next_start_r = start_r << depth1;
   unsigned long next_skip_r = skip_r << depth1;

   if (depth2 <= crossover)
   {
#if USE_THREADS
      threads = THREADS;
      lg_threads = LG_THREADS;
      while (threads > length1) 
      {
         lg_threads--;
         threads >>= 1;
      }
      
      unsigned long new_length1 = ((length1-1)>>lg_threads)+1;
      for (i = 0; i < threads; i++)
      {
         fft_params[i].next_skip = next_skip;
         fft_params[i].skip = skip;
         fft_params[i].next_skip_r = next_skip_r;
         fft_params[i].start_r = next_start_r;
         fft_params[i].depth = depth2;
         fft_params[i].scratch = scratch;
         fft_params[i].n = n;
      }
      for (i = 0; i < threads-1; i++)
      {
         fft_params[i].length = new_length1;
         fft_params[i].start = next_start + i*new_length1*next_skip;
      }
      fft_params[i].length = length1 - new_length1*i;
      fft_params[i].start = next_start + i*new_length1*next_skip;
      
      for (i = 0; i < threads; i++)
      {
         fft_params[i].thread = thread_num+i;
         pthread_create(thread_arr+i, &attr, ifft_loop2, (void*) (fft_params+i));
             
      }
      for (unsigned long i = 0; i < threads; i++)
      {
         pthread_join(thread_arr[i], NULL);
      } 
      
#else
      for (unsigned long i = 0; i < length1; i++, next_start += next_skip)
      {
         ifft_recursive(next_start, skip, next_start_r, skip_r, depth2, scratch, n);
      }    
#endif
   } else
   {
      for (unsigned long i = 0; i < length1; i++, next_start += next_skip)
         ifft_main(next_start, skip, next_start_r, next_skip_r,
                depth2, scratch, n, crossover);
   }
   next_start = start;
   next_start_r = start_r;
   next_skip_r = skip_r << depth2;

   if (depth1 <= crossover)
   {
#if USE_THREADS
      threads = THREADS;
      lg_threads = LG_THREADS;
      while (threads > length2) 
      {
         lg_threads--;
         threads >>= 1;
      }
      
      unsigned long new_length2 = ((length2-1)>>lg_threads)+1;
      for (i = 0; i < threads; i++)
      {
         fft_params[i].next_skip = next_skip;
         fft_params[i].skip = skip;
         fft_params[i].next_skip_r = next_skip_r;
         fft_params[i].skip_r = skip_r;
         fft_params[i].depth = depth1;
         fft_params[i].scratch = scratch;
         fft_params[i].n = n;
      }
      for (i = 0; i < threads-1; i++)
      {
         fft_params[i].length = new_length2;
         fft_params[i].start = next_start + i*new_length2*skip;
         fft_params[i].start_r = next_start_r + i*new_length2*skip_r;
      }
      fft_params[i].length = length2 - new_length2*i;
      fft_params[i].start = next_start + i*new_length2*skip;
      fft_params[i].start_r = next_start_r + i*new_length2*skip_r;
      
      for (i = 0; i < threads; i++)
      {
         fft_params[i].thread = thread_num+i;
         pthread_create(thread_arr+i, &attr, ifft_loop, (void*) (fft_params+i));
             
      }
      for (unsigned long i = 0; i < threads; i++)
      {
         pthread_join(thread_arr[i], NULL);
      } 
      
#else
      for (unsigned long i = 0; i < length2; i++, next_start += skip, 
              next_start_r += skip_r)
      {
         ifft_recursive(next_start, skip, next_start_r, skip_r, depth1, scratch, n);
      }    
#endif  
   } else
   {
       for (unsigned long i = 0; i < length2; i++, next_start += skip,
        next_start_r += skip_r)
           ifft_main(next_start, next_skip, next_start_r, next_skip_r,
                depth1, scratch, n, crossover);
   }
   pthread_attr_destroy(&attr);
}

/*=============================================================================

Old recursive FFT code, but parallelised

==============================================================================*/

void* fft_inner(void* fft_p)
{
   fft_s* fft_params = (fft_s*) fft_p;
   
   mp_limb_t** start = fft_params->start;
   unsigned long length = fft_params->length;
   mp_limb_t** scratch = fft_params->scratch;
   unsigned long r = fft_params->r; 
   unsigned long n = fft_params->n;
   
   unsigned long half = length >> 1;
   mp_limb_t** middle = start + half;
   unsigned long offset, r_mult;
    
   
      // scratch = a - b
      mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
      
      // a += b
      mpn_add_n(start[0], start[0], middle[0], n + 1);
      
      swap_limb_ptrs(middle, scratch);
      
      if (half == 1)
         return NULL;
      
      if (r & (FLINT_BITS_PER_LIMB - 1))
      {
         // not shifting by a multiple of the limb length
         for (offset = 1, r_mult = r; offset < half; offset++, r_mult += r)
            fft_butterfly_bits(start + offset, middle + offset,
                               scratch, r_mult, n);
      }
      else
      {
         // shifting by a multiple of the limb length
         unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
         for (offset = 1, r_mult = rdivB; offset < half;
              offset++, r_mult += rdivB)
            fft_butterfly_limbs(start + offset, middle + offset,
                                scratch, r_mult, n);
      }
   
   fft_s fft_params1;
   
   fft_params1.start = middle;
   fft_params1.length = (length >> 1);
   fft_params1.scratch = scratch;
   fft_params1.r = (r << 1); 
   fft_params1.n = n;
   
   fft_inner((void*)&fft_params1);
   
   fft_params1.start = start;
   
   fft_inner((void*)&fft_params1);
}

void fft(mp_limb_t** start, unsigned long length, mp_limb_t** scratch,
         unsigned long r, unsigned long n)
{
   pthread_t thread1;
   pthread_t thread2;
   
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
    
   unsigned long half = length >> 1;
   mp_limb_t** middle = start + half;
   unsigned long offset, r_mult;
    
   
      copy_limbs(middle[0], start[0], n + 1);
      
      if (r & (FLINT_BITS_PER_LIMB - 1))
      {
         // not shifting by a multiple of the limb length
         for (offset = 1, r_mult = r; offset < half; offset++, r_mult += r)
            rotate_mod_p_bits(middle[offset], start[offset], r_mult, n);
      }
      else
      {
         // shifting by a multiple of the limb length
         unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
         for (offset = 1, r_mult = rdivB; offset < half;
              offset++, r_mult += rdivB)
            rotate_mod_p_limbs(middle[offset], start[offset], r_mult, n);
      }
      
   fft_s fft_params1;
   fft_s fft_params2;
   
   fft_params1.start = middle;
   fft_params1.length = (length >> 1);
   fft_params1.scratch = scratch;
   fft_params1.r = (r << 1); 
   fft_params1.n = n;
   
   fft_params2.start = start;
   fft_params2.length = (length >> 1);
   fft_params2.scratch = scratch+16;
   fft_params2.r = (r << 1); 
   fft_params2.n = n;
   
#if USE_THREADS
   pthread_create(&thread1, &attr, fft_inner, (void*) &fft_params1);
   pthread_create(&thread2, &attr, fft_inner, (void*) &fft_params2);
   
   pthread_join(thread1, NULL);
   pthread_join(thread2, NULL);
    
   pthread_attr_destroy(&attr);         
#else
   fft_inner((void*) &fft_params1);
   fft_inner((void*) &fft_params2);
#endif
}

void * ifft_inner(void* ifft_params)
{
   fft_s * ifft_p = (fft_s*) ifft_params;
    
   mp_limb_t** start = ifft_p->start;
   unsigned long length = ifft_p->length;
   mp_limb_t** scratch = ifft_p->scratch;
   unsigned long r = ifft_p->r;
   unsigned long n = ifft_p->n;  
     
   unsigned long half = length >> 1;
   mp_limb_t** middle = start + half;
   unsigned long offset, r_mult = 0;
   unsigned long nB = FLINT_BITS_PER_LIMB * n;
   
   if (half == 1)
   {
      // scratch = a - b
      mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
      
      // a += b
      mpn_add_n(start[0], start[0], middle[0], n + 1);
      
      swap_limb_ptrs(middle, scratch);
      
      return NULL;
   }
   
   fft_s ifft_params1;
   
   ifft_params1.start = middle;
   ifft_params1.length = half;
   ifft_params1.scratch = scratch;
   ifft_params1.r = (r << 1); 
   ifft_params1.n = n;
   
   ifft_inner((void*)&ifft_params1);
   
   ifft_params1.start = start;
   
   ifft_inner((void*)&ifft_params1);
   
   // scratch = a - b
   mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
   
   // a += b
   mpn_add_n(start[0], start[0], middle[0], n + 1);
   
   swap_limb_ptrs(middle, scratch);
   
   if (r & (FLINT_BITS_PER_LIMB - 1))
   {
      // not shifting by a multiple of the limb length
      for (offset = 1, r_mult = nB - r; offset < half; offset++, r_mult -= r)
      {
         ifft_butterfly_bits(start + offset, middle + offset,
                             scratch, r_mult, n);
      }
   }
   else
   {
      // shifting by a multiple of the limb length
      unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
      for (offset = 1, r_mult = n - rdivB; offset < half;
           offset++, r_mult -= rdivB)
         ifft_butterfly_limbs(start + offset, middle + offset,
                              scratch, r_mult, n);
   }
}

void ifft(mp_limb_t** start, unsigned long length, mp_limb_t** scratch,
          unsigned long r, unsigned long n)
{
   pthread_t thread1;
   pthread_t thread2;
   
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
   
   unsigned long half = length >> 1;
   mp_limb_t** middle = start + half;
   unsigned long offset, r_mult = 0;
   unsigned long nB = FLINT_BITS_PER_LIMB * n;
   
   if (half == 1)
   {
      // scratch = a - b
      mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
      
      // a += b
      mpn_add_n(start[0], start[0], middle[0], n + 1);
      
      swap_limb_ptrs(middle, scratch);
      
      return;
   }
   
   fft_s ifft_params1;
   fft_s ifft_params2;
   
   ifft_params1.start = middle;
   ifft_params1.length = half;
   ifft_params1.scratch = scratch;
   ifft_params1.r = (r << 1); 
   ifft_params1.n = n;
   
   ifft_params2.start = start;
   ifft_params2.length = half;
   ifft_params2.scratch = scratch+16;
   ifft_params2.r = (r << 1); 
   ifft_params2.n = n;
   
#if USE_THREADS
   pthread_create(&thread1, &attr, ifft_inner, (void*) &ifft_params1);
   pthread_create(&thread2, &attr, ifft_inner, (void*) &ifft_params2);
   
   pthread_join(thread1, NULL);
   pthread_join(thread2, NULL);
    
   pthread_attr_destroy(&attr);         
#else
   ifft_inner((void*) &ifft_params1);
   ifft_inner((void*) &ifft_params2);
#endif
   
   // scratch = a - b
   mpn_sub_n(scratch[0], start[0], middle[0], n + 1);
   
   // a += b
   mpn_add_n(start[0], start[0], middle[0], n + 1);
   
   swap_limb_ptrs(middle, scratch);
   
   if (r & (FLINT_BITS_PER_LIMB - 1))
   {
      // not shifting by a multiple of the limb length
      for (offset = 1, r_mult = nB - r; offset < half; offset++, r_mult -= r)
      {
         ifft_butterfly_bits(start + offset, middle + offset,
                             scratch, r_mult, n);
      }
   }
   else
   {
      // shifting by a multiple of the limb length
      unsigned long rdivB = r / FLINT_BITS_PER_LIMB;
      for (offset = 1, r_mult = n - rdivB; offset < half;
           offset++, r_mult -= rdivB)
         ifft_butterfly_limbs(start + offset, middle + offset,
                              scratch, r_mult, n);
   }
}

/*============================================================================

Kronecker-Schonhage routines 

=============================================================================*/

void KSMul_bits(mpz_t* res, const mpz_t* data1, const mpz_t* data2, 
           const unsigned long orig_length, const unsigned long length2,
           const unsigned long coeff_bits, const unsigned long log_length)
{
   unsigned long output_bits = 2*coeff_bits+log_length;
   const unsigned long array_bits = output_bits*orig_length;
   const unsigned long array_limbs = (array_bits-1)/FLINT_BITS_PER_LIMB+1;

   const unsigned long coeffs_per_limb = FLINT_BITS_PER_LIMB/output_bits;
   
   mp_limb_t* array1 = (mp_limb_t*) limb_alloc(array_limbs, 0);
   mp_limb_t* array2 = (mp_limb_t*) limb_alloc(array_limbs, 0);
   mp_limb_t* output = (mp_limb_t*) limb_alloc(2*array_limbs, 0);
   
   unsigned long k;

   clear_limbs(array1, array_limbs);
   clear_limbs(array2, array_limbs);
   
   // pile up first lot of input coefficients into single large integer
   ssmul_convert_in_bits(array1, data1, array_limbs, 
                                   orig_length, output_bits, coeffs_per_limb);

   // pile up second lot of input coefficients into single large integer
   ssmul_convert_in_bits(array2, data2, array_limbs, 
                                       length2, output_bits, coeffs_per_limb);

   // do multiplication of large integers
   
   if (array_bits > 164000)
   {   
      Z_SSMul(output, array1, array2, array_limbs);  
   } else 
   {
      mpn_mul_n(output, array1, array2, array_limbs);
   }
   
   // convert cofficients out
   for (k = 0; k < orig_length + length2 - 1; k++) mpz_set_ui(res[k],0);
   
   ssmul_convert_out_bits(res, output, 
                       orig_length + length2 - 1, output_bits, coeffs_per_limb);
                
   limb_release();limb_release();limb_release();
   
}

void KSMul_bits_signed(mpz_t* res, mpz_t* data1, mpz_t* data2, 
           const unsigned long orig_length, const unsigned long length2,
           const unsigned long coeff_bits, const unsigned long log_length)
{
   int sign1 = 1;
   int sign2 = 2;

   unsigned long k = orig_length-1;
   while ((mpz_sgn(data1[k]) == 0)&&(k>0)) k--;

   if ((k == 0) && (mpz_sgn(data1[0]) == 0))
   {
      for (k = 0; k < orig_length + length2 - 1; k++) mpz_set_ui(res[k],0);
      return;
   }

   if (mpz_sgn(data1[k])<0)
   {
      for (k = 0; k < orig_length; k++) mpz_neg(data1[k],data1[k]);
      sign1 = -1;
   }

   k = length2 - 1;
   while ((mpz_sgn(data2[k]) == 0)&&(k>0)) k--;

   if ((k == 0) && (mpz_sgn(data2[0]) == 0))
   {
      if (sign1<0)
      {
         for (k = 0; k < orig_length; k++) mpz_neg(data1[k],data1[k]);
      }
      for (k = 0; k < orig_length + length2 - 1; k++) mpz_set_ui(res[k],0);
      return;
   }

   if (mpz_sgn(data2[k])<0)
   {
      for (k = 0; k < length2; k++) mpz_neg(data2[k],data2[k]);
      sign2 = -1;
   }

   const unsigned long output_bits = 2*coeff_bits+log_length+2;
   const unsigned long array_bits = output_bits*orig_length;
   const unsigned long array_limbs = (array_bits-1)/FLINT_BITS_PER_LIMB+1;

   const unsigned long coeffs_per_limb = FLINT_BITS_PER_LIMB/output_bits;

   mp_limb_t* array1 = (mp_limb_t*) limb_alloc(array_limbs, 0);
   mp_limb_t* array2 = (mp_limb_t*) limb_alloc(array_limbs, 0);
   mp_limb_t* output = (mp_limb_t*) limb_alloc(2*array_limbs, 0);
   
   clear_limbs(array1, array_limbs);
   clear_limbs(array2, array_limbs);
   ssmul_convert_in_bits_signed(array1, data1, array_limbs, 
                                orig_length, output_bits, coeffs_per_limb, 0);
   
   ssmul_convert_in_bits_signed(array2, data2, array_limbs, 
                                    length2, output_bits, coeffs_per_limb, 0);
   // do multiplication of large integers

   if (array_bits > 164000)
   {
      Z_SSMul(output, array1, array2, array_limbs);
   } else
   {
      mpn_mul_n(output, array1, array2, array_limbs);
   }

   for (k = 0; k < orig_length + length2 - 1; k++) mpz_set_ui(res[k],0);
   // convert cofficients out
   ssmul_convert_out_bits_signed(res, output, orig_length + length2 - 1, output_bits, coeffs_per_limb);

   if (sign1*sign2 < 0)
   {
      for (k = 0; k < orig_length + length2 - 1; k++) mpz_neg(res[k],res[k]);
   }
   if (sign1<0)
   {
      for (k = 0; k < orig_length; k++) mpz_neg(data1[k],data1[k]);
   }
   if (sign2<0)
   {
      for (k = 0; k < length2; k++) mpz_neg(data2[k],data2[k]);
   }

   limb_release();limb_release();limb_release();

}

void KSMul_bytes(mpz_t* res, mpz_t* data1, mpz_t* data2, 
           const unsigned long orig_length, const unsigned long length2,
           const unsigned long coeff_bits, const unsigned long log_length)
{
   const unsigned long output_bits = 2*coeff_bits+log_length;
   const unsigned long coeff_bytes = ((output_bits-1)>>3)+1;
   const unsigned long array_bits = (coeff_bytes<<3)*orig_length;
   const unsigned long array_limbs = array_bits/FLINT_BITS_PER_LIMB+1;

   mp_limb_t* array1 = (mp_limb_t*) limb_alloc(array_limbs,0);
   mp_limb_t* array2 = (mp_limb_t*) limb_alloc(array_limbs,0);
   mp_limb_t* output = (mp_limb_t*) limb_alloc(2*array_limbs,0);
   
   unsigned long i;
   
   mpz_t temp;
   mpz_init(temp);
   
   clear_limbs(array1,array_limbs);
   clear_limbs(array2,array_limbs);
   
   ssmul_convert_in_bytes(array1, data1, orig_length, coeff_bytes, array_limbs); 
   ssmul_convert_in_bytes(array2, data2, length2, coeff_bytes, array_limbs); 
   
   if (array_bits > 164000)
   {   
      Z_SSMul(output, array1, array2, array_limbs);  
   } else 
   {
      mpn_mul_n(output, array1, array2, array_limbs);
   }
   
   for (i = 0; i < orig_length + length2 - 1; i++) mpz_set_ui(res[i],0);
   
   ssmul_convert_out_bytes(res, output, orig_length + length2 - 1, coeff_bytes, temp);
   
   mpz_clear(temp);
   limb_release();limb_release();limb_release();
   
}

void KSMul_bytes_signed(mpz_t* res, mpz_t* data1, mpz_t* data2, 
           const unsigned long orig_length, const unsigned long length2,
           const unsigned long coeff_bits, const unsigned long log_length)
{
   int sign1 = 1;
   int sign2 = 2;

   unsigned long k = orig_length-1;
   while ((mpz_sgn(data1[k]) == 0)&&(k>0)) k--;

   if ((k == 0) && (mpz_sgn(data1[0]) == 0))
   {
      for (k = 0; k < orig_length + length2 - 1; k++) mpz_set_ui(res[k],0);
      return;
   }

   if (mpz_sgn(data1[k])<0)
   {
      for (k = 0; k < orig_length; k++) mpz_neg(data1[k],data1[k]);
      sign1 = -1;
   }

   k = length2 - 1;
   while ((mpz_sgn(data2[k]) == 0)&&(k>0)) k--;

   if ((k == 0) && (mpz_sgn(data2[0]) == 0))
   {
      if (sign1<0)
      {
         for (k = 0; k < orig_length; k++) mpz_neg(data1[k],data1[k]);
      }
      for (k = 0; k < orig_length + length2 - 1; k++) mpz_set_ui(res[k],0);
      return;
   }

   if (mpz_sgn(data2[k])<0)
   {
      for (k = 0; k < length2; k++) mpz_neg(data2[k],data2[k]);
      sign2 = -1;
   }

   const unsigned long output_bits = 2*coeff_bits+log_length+2;
   const unsigned long coeff_bytes = ((output_bits-1)>>3)+1;
   const unsigned long array_bits = (coeff_bytes<<3)*orig_length;
   const unsigned long array_limbs = array_bits/FLINT_BITS_PER_LIMB+1;

   mp_limb_t* array1 = (mp_limb_t*) limb_alloc(array_limbs,0);
   mp_limb_t* array2 = (mp_limb_t*) limb_alloc(array_limbs,0);
   mp_limb_t* output = (mp_limb_t*) limb_alloc(2*array_limbs,0);
   
   unsigned long i;
   
   mpz_t temp;
   mpz_init(temp);
   
   clear_limbs(array1,array_limbs);
   clear_limbs(array2,array_limbs);
   
   // We don't need a special function to convert in signed bytes
   ssmul_convert_in_bytes(array1, data1, orig_length, coeff_bytes, array_limbs); 
   ssmul_convert_in_bytes(array2, data2, length2, coeff_bytes, array_limbs); 
   
   if (array_bits > 164000)
   {   
      Z_SSMul(output, array1, array2, array_limbs);  
   } else 
   {
      mpn_mul_n(output, array1, array2, array_limbs);
   }
   
   for (i = 0; i < orig_length + length2 - 1; i++) mpz_set_ui(res[i],0);
   
   ssmul_convert_out_bytes_signed(res, output, orig_length + length2 - 1, coeff_bytes, temp);
 
   if (sign1*sign2 < 0)
   {
      for (k = 0; k < orig_length + length2 - 1; k++) mpz_neg(res[k],res[k]);
   }
   if (sign1<0)
   {
      for (k = 0; k < orig_length; k++) mpz_neg(data1[k],data1[k]);
   }
   if (sign2<0)
   {
      for (k = 0; k < length2; k++) mpz_neg(data2[k],data2[k]);
   }
  
   mpz_clear(temp);
   limb_release();limb_release();limb_release();
   
}

void* pointwise_mult(void* point_params)
{
   point_t * point_p = (point_t*) point_params;
   unsigned long length = point_p->length;
   mp_limb_t** array1 = point_p->array1;
   mp_limb_t** array2 = point_p->array2;
   mp_limb_t* array3 = point_p->array3;
   unsigned long n = point_p->n;
   
   unsigned long i, j;
   
   // pointwise multiplies
   for (i = 0; i < length; i++)
   {  
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array2[i+1], j);

      reduce_mod_p_exact(array1[i], n);
      reduce_mod_p_exact(array2[i], n);

      // We've reduced the coefficients to the range 0 <= x < p.
      // If they are < p-1, then they fit into n limbs. Otherwise we have
      // special cases as dealt with in the following code.

      if (array1[i][n])
      {
         // special case: array1[i] == -1 (mod p),
         // so answer is (-1) * array2[i]
         negate_limbs(array1[i], array2[i], n + 1);
         reduce_mod_p_exact(array1[i], n);
      }
      else if (array2[i][n])
      {
         // special case: array2[i] == -1 (mod p),
         // so answer is (-1) * array1[i]
         
         negate_limbs(array1[i], array1[i], n + 1);
         reduce_mod_p_exact(array1[i], n);
      }
      else
      {
#if !USE_THREADS
         // usual case, just do the multiplication
         if (n > 164000/FLINT_BITS_PER_LIMB)
         {   
             Z_SSMul(array3, array1[i], array2[i], n);  
         } else 
         {
            mpn_mul_n(array3, array1[i], array2[i], n);
         }
#else
         // Z_SSMul is not threadsafe due to calls to limb_alloc, so for now:
         mpn_mul_n(array3, array1[i], array2[i], n);
         // reduce the result mod p
#endif    
         array1[i][n] = -mpn_sub_n(array1[i], array3, array3 + n, n);
      }
   }
}

/*
The core of the ssmul routine. It performs the forward FFT's, pointwise
mults and IFFT. It accepts two arrays containing pointers to 
2^log_length
coefficients, each coefficient being packed with n+1 limbs space. Only 
the
first length1 (resp. length2) coefficients of array1 (resp. array2) are
assumed to contain actual input data, the remainder are treated as 
zeroes.
length1 and length2 must be positive. scratch contains n+1 limbs of
scratch space for the FFTs and IFFT.

The output is stored in array1; only the first (length1 + length2 - 1)
coefficients will contain actual data, the remainder is garbage. The 
output
coefficients are not reduced mod p. log_length must be large enough so 
that
the output will fit (i.e. length1+length2-1 <= 2^log_length).
*/
void ssmul_main(mp_limb_t** array1, unsigned long length1,
                 mp_limb_t** array2, unsigned long length2,
                 mp_limb_t** scratch, unsigned long log_length,
                 unsigned long n)
{
    unsigned long i, j;

    // root of unity is sqrt2^r
    // todo: add assertion here to check divisibility by 2^log_length
    // todo: change this to shift by FLINT_LG_BITS_PER_LIMB
    unsigned long r = (4*n*FLINT_BITS_PER_LIMB) >> log_length;

    // number of fourier coefficients required (i.e. length of output)
    // todo: add assertion here to check output_length <= 2^log_length
    unsigned long output_length = length1 + length2 - 1;
    
#if USE_THREADS
    unsigned long threads = THREADS;
    unsigned long lg_threads = LG_THREADS;
    while (threads > output_length-1)
    {
       lg_threads--;
       threads >>= 1;
    } 
    lg_threads++;
    threads <<= 1;
    if (output_length < 2)
    {
       threads = 1;
       lg_threads = 0;
    }
    pthread_t thread[threads];
   
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
   
#else
    unsigned long threads = 1;
    unsigned long lg_threads = 0;
#endif

    point_t point[threads];

    // do FFTs

    ssfft_fft(array1, 1, log_length, length1, output_length, 0, r, n, scratch);
    ssfft_fft(array2, 1, log_length, length2, output_length, 0, r, n, scratch);

    // pointwise multiplies

    mp_limb_t* mul_scratch = (mp_limb_t*) limb_alloc(2*n*threads, 0);
    
   for (i = 0; i < threads-1; i++)
   {
      point[i].length = (output_length >> lg_threads);
      point[i].array1 = array1 + i*(output_length >> lg_threads);
      point[i].array2 = array2 + i*(output_length >> lg_threads);
      point[i].array3 = mul_scratch + 2*n*i;
      point[i].n = n;
   }
   point[i].length = (output_length - (output_length >> lg_threads)*i);
   point[i].array1 = array1 + i*(output_length >> lg_threads);
   point[i].array2 = array2 + i*(output_length >> lg_threads);
   point[i].array3 = mul_scratch + 2*n*i;
   point[i].n = n;

#if USE_THREADS   
   for (unsigned long i = 0; i < threads; i++)
   {
      pthread_create(&thread[i], &attr, pointwise_mult, (void*) &point[i]);
   }
   
   for (unsigned long i = 0; i < threads; i++)
   {
      pthread_join(thread[i], NULL);
   }
    
   pthread_attr_destroy(&attr);  
#else
   pointwise_mult((void*) point);
#endif

    limb_release();

    // do inverse FFT
    ssfft_ifft(array1, 1, log_length, output_length, output_length, 0,
               0, r, n, scratch);
}

/* 
The core of the old ssmul routine, but parallelised. It performs the forward 
FFT's, pointwise mults and IFFT. It accepts two arrays containing pointers 
to length coefficients, each coefficient being packed with n+1 limbs space. 
log_length is the ceiling of the log_2 of length. array3 contains scratch 
space to contain the result of a pointwise mult. scratch contains n+1 
limbs of scratch space for the FFT's and IFFT.
*/
void ssmul_main_old(mp_limb_t** array1, mp_limb_t** array2, mp_limb_t* array3,
           unsigned long length, unsigned long log_length, unsigned long n,
           unsigned long r, mp_limb_t** scratch)
{
#if USE_THREADS
   unsigned long threads = THREADS;
   unsigned long lg_threads = LG_THREADS;
   while (threads > length)
   {
      lg_threads--;
      threads >>= 1;
   } 
   pthread_t thread[threads];
   
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setscope(&attr,PTHREAD_SCOPE_SYSTEM);
   
#else
   unsigned long threads = 1;
   unsigned long lg_threads = 0;
#endif

   point_t point[threads];
   
   unsigned long i, j;
   
   // do FFT's

   fft_main(array1, 1, 0, r, log_length+1, scratch, n, 1);
   fft_main(array2, 1, 0, r, log_length+1, scratch, n, 1);

   for (unsigned long i = 0; i < threads; i++)
   {
      point[i].length = (length >> lg_threads);
      point[i].array1 = array1 + i*(length >> lg_threads);
      point[i].array2 = array2 + i*(length >> lg_threads);
      point[i].array3 = array3 + 2*n*i;
      point[i].n = n;
   }

#if USE_THREADS   
   for (unsigned long i = 0; i < threads; i++)
   {
      pthread_create(&thread[i], &attr, pointwise_mult, (void*) &point[i]);
   }
   
   for (unsigned long i = 0; i < threads; i++)
   {
      pthread_join(thread[i], NULL);
   }
    
   pthread_attr_destroy(&attr);  
#else
   pointwise_mult((void*) point);
#endif

   // do inverse FFT
   ifft_main(array1, 1, 0, r, log_length+1, scratch, n);  
}

/* 
For some reason, when applying the Kronecker-Schoenhage trick, with the
Harvey-Zimmerman improvements, for less than FLINT_BITS_PER_LIMB output bits, the code needs
to all be in the one place, or too much time is lost jumping around everywhere.
So we repeat all the necessary things here.
*/
void KSHZ_mul(mpz_t* res, mpz_t* data1, mpz_t* data2, 
           const unsigned long orig_length, const unsigned long length2,
           const unsigned long coeff_bits, unsigned long log_length, int sign)
{
   unsigned long output_bits = 2*coeff_bits + log_length;
   if (sign) output_bits+=2;
   unsigned long orig_output_bits = output_bits;
   unsigned long coeffs_per_limb = FLINT_BITS_PER_LIMB/orig_output_bits;    
   
   unsigned long i, j, k; 
   unsigned long length;
   
   unsigned long bundle = 1; 
   
   // approximate value for bundle (should be right to within 3 or so)
   while ((2*bundle-1)*bundle*coeff_bits < orig_length) bundle++;
      
   length = (orig_length-1)/bundle + 1;
   while ((1<<log_length) > length) log_length--;
      
   //recompute a more optimal bundle size
   length = 1<<log_length;
   bundle = (orig_length-1)/length+1;
      
   unsigned long bundle_bits;
   
   if (!sign) bundle_bits = (bundle-1)*output_bits+coeff_bits;
   else bundle_bits = bundle*output_bits;
   //bundle_bits already incorporates extra space for sign bits
   output_bits = 2*bundle_bits+log_length;
      
   // compute parameters for fft
   output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
   unsigned long n = (output_bits - 1) / FLINT_BITS_PER_LIMB + 1;
   length = 1 << (log_length + 1);
   unsigned long r = n*FLINT_BITS_PER_LIMB * 2 / length;

   mp_limb_t* array = (mp_limb_t*) limb_alloc((2*length + 16) * (n+1) + 2*n*THREADS,0);
   mp_limb_t** pointarr = (mp_limb_t**) limb_alloc(2*length + 16*16,0);

   // partition up array into coefficient buffers of length n+1...
   
   pointarr[0] = array;
   for (i = 0; i < 2*length-1; i++)
      pointarr[i+1] = pointarr[i] + (n+1);
   for (i = 0; i < 16; i++)
      pointarr[2*length+16*i] = array+(2*length+i)*(n+1);
      
   // ...and pointwise multiplication working space
   
   mp_limb_t* array3 = array + (2*length+16)*(n+1);   

   mp_limb_t** array1 = pointarr;
   mp_limb_t** array2 = pointarr + length;
   mp_limb_t** scratch = pointarr + 2*length;
   
   ssmul_convert_in(array1, data1, length, sign, 1, 1, bundle, orig_length, 
                          n, 0, 0, orig_output_bits, coeffs_per_limb);
    
   ssmul_convert_in(array2, data2, length, sign, 1, 1, bundle, length2, 
                          n, 0, 0, orig_output_bits, coeffs_per_limb);
   
   ssmul_main_old(array1, array2, array3, length, log_length, n, r, scratch);
   
   // We have to clear the output polynomial, since we have to add to its
   // coefficients rather than just copy into them
   
   for (i = 0; i < orig_length+length2-1; i++) mpz_set_ui(res[i],0);
            
   for (k=0, i = 0; i < length-1; i++, k+=bundle)
   { 
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
      rotate_right_bits(array1[i], array1[i], log_length+1, n);
      reduce_mod_p_exact(array1[i], n);
      if (sign)
      {
         if (array1[i][n] || (array1[i][n-1] >> (FLINT_BITS_PER_LIMB - 1)))
         {
            mpn_sub_1(array1[i], array1[i], n, 1);
         }
      }
      
      for (j = 0; j < n; j += 8) FLINT_PREFETCH(array1[i+1], j);
      if (k + 2*bundle - 1 <= orig_length + length2 - 1)
      {
            if (!sign) ssmul_convert_out_bits(res + k, array1[i], 2*bundle - 1, orig_output_bits, coeffs_per_limb);
            else ssmul_convert_out_bits_signed(res + k, array1[i], 2*bundle - 1, orig_output_bits, coeffs_per_limb);
      } else if (k < orig_length + length2 - 1)
      {
            if (!sign) ssmul_convert_out_bits(res + k, array1[i], orig_length + length2 - 1 - k, orig_output_bits, coeffs_per_limb);
            else ssmul_convert_out_bits_signed(res + k, array1[i], orig_length +length2 - 1 - k, orig_output_bits, coeffs_per_limb);
      }    
   } 
   
   limb_release();
   limb_release();
}


/*============================================================================

SSMul tuning code for split, bundle and bitpack tricks

=============================================================================*/

/*
Tuning parameters for split coefficients. Returns the optimal number of smaller
coefficients to split each original coefficient into (which may be 1)
*/
inline unsigned long get_split_tuning_parameters(unsigned long * input_limbs, 
                    unsigned long * orig_limbs, unsigned long * output_bits, 
                    unsigned long * log_length, unsigned long coeff_bits, 
                    unsigned long orig_length, int sign)
{     
   unsigned long split = 1;
   
   // the number of limbs in each piece of the (possibly) split coefficients
   *input_limbs = (coeff_bits-1)/FLINT_BITS_PER_LIMB+1;
      
   // the original number of limbs per coefficient, before splitting
   *orig_limbs = *input_limbs;
   
   // if we should split coefficients
   if ((4*orig_length < *output_bits) && (coeff_bits > 32000)) 
   {
      while (((1<<(*log_length + 2)) <= *output_bits)&&(*input_limbs > 1)) 
      {
         (*log_length)++;
         // we will need nearly twice as many limbs as usual so full 
         // output coefficients don't overlap
         *input_limbs = (*orig_limbs-1)/((1<<(*log_length))/(2*orig_length))+1;
         if (sign) *output_bits = 2*FLINT_BITS_PER_LIMB*(*input_limbs)+*log_length + 2;
         else *output_bits = 2*FLINT_BITS_PER_LIMB*(*input_limbs)+*log_length;
      }
      // check if we can fit in a smaller convolution
      split = (*orig_limbs-1)/(*input_limbs)+1;
      if (orig_length*(2*split-1) <= 1<<(*log_length-1)) (*log_length)--;
   }
    
   return split;
}

/*
Tuning parameters for bundling of parameters (KS-Harvey/Zimmerman). Returns
the optimal number of coefficients to bundle into a single coefficient.
*/
inline unsigned long get_bundle_tuning_parameters(unsigned long * skip_limbs, 
                   unsigned long * output_bits, unsigned long * input_limbs, 
                   unsigned long * log_length, unsigned long * length, 
                   unsigned long orig_length, unsigned long coeff_bits, int sign)
{
   unsigned long bundle = 1;
   
   // if we should use KS (bundle coefficients)
   if ((4*coeff_bits < orig_length) || ((coeff_bits < orig_length) && (coeff_bits < 256)))
   {
      // the coefficients will be put in bundles with each coefficient
      // starting this many limbs from the start of the previous one
      // We break into 8 bit limbs
      *skip_limbs = (*output_bits-1)/8+1;
         
      // approximate value for bundle (should be right to within 3 or so)
      while ((2*bundle-1)*bundle*coeff_bits < orig_length) bundle++;
        
      // set new parameters for input length, *input* coefficient bits, and output_bits
      // based on the computed optimal bundle size
      *length = (orig_length-1)/bundle + 1;
      while (1<<(*log_length) > *length) (*log_length)--;
      
      if (*log_length > 3)
      {
         if (coeff_bits <= 16) 
         {
            if (orig_length <= 60) (*log_length)-=3;
         } else if (coeff_bits <= 100) 
         {
            if (orig_length <= 720) (*log_length)-=2;
            else if (orig_length <= 8192) (*log_length)--;
         } else if (coeff_bits > 512) (*log_length)++;
      }
       
      *input_limbs = (coeff_bits-1)/8+1;
        
      // output_bits will be a multiple of the limb size
      *output_bits = *skip_limbs*8;
      
      //recompute a more optimal bundle size
      *length = 1<<(*log_length);
      bundle = (orig_length-1)/(*length)+1;
      
      if (sign)
      {
         unsigned long new_coeff_bits = bundle*(*output_bits);
         // new_coeff_bits already allows space for signs
         *output_bits = 2*new_coeff_bits+(*log_length);
      } else
      {
         unsigned long new_coeff_bits = (bundle-1)*(*output_bits)+coeff_bits;
         *output_bits = 2*new_coeff_bits+(*log_length);
      }
   }
   
   return bundle;
}

/*============================================================================

Main Schonhage-Strassen routine exported from this module.

=============================================================================*/

/* 
Main polynomial multiplication routine. It performs a Schoenhage-Strassen
algorithm, calling the above FFT and IFFT functions. res is set to the 
polynomial data1*data2.

At approximately the
point when the degree of the polynomials exceeds the input coefficient bit
size, it starts bundling coefficients together using the Kronecker-Schoenhage
technique as refined by Paul Zimmerman and David Harvey.

log_length is the ceiling of the log_2 of the length of the *input* 
polynomial (length = deg+1).  coeff_bits is the number of bits in the 
maximimum *input* coefficient. We assume that the length of the output
polynomial will always be less than 2^B where B = FLINT_BITS_PER_LIMB

**IMPORTANT** Expects the length of poly1 to be >= length of poly2. 
orig_length is the greater of the two.

todo: write an SSSqr for squaring instead of multiplying

*/
void SSMul(Zvec outpoly, Zvec poly1, Zvec poly2, unsigned long coeff_bits, int sign)
{
   mpz_t* data1 = poly1.coords;
   mpz_t* data2 = poly2.coords;
   unsigned long orig_length = poly1.length;;
   unsigned long length2 = poly2.length;;
   unsigned long trunc_length;

   mpz_t* res = outpoly.coords;
   
   unsigned long log_length=0;
   while ((1<<log_length) < orig_length) log_length++;
             
   unsigned long output_bits;
   if (sign) output_bits = 2*coeff_bits + log_length + 2;
   else output_bits = 2*coeff_bits + log_length;

   unsigned long orig_output_bits = output_bits;
   unsigned long coeffs_per_limb = FLINT_BITS_PER_LIMB/orig_output_bits;    
   
   // Check if we want to use special high speed KS routine if output 
   // coefficients fit into FLINT_BITS_PER_LIMB bits
   
   if ((output_bits < FLINT_BITS_PER_LIMB) && (orig_length < 4096))
   {
      if (!sign) KSMul_bits(res, data1, data2, orig_length, length2, coeff_bits, log_length);
      else KSMul_bits_signed(res, data1, data2, orig_length, length2, coeff_bits, log_length);
      return;
   } 

   if (((coeff_bits*orig_length < 55000UL) && (coeff_bits < 360) && (orig_length < 4096)) || ((orig_length == 9) && (coeff_bits < 32000)))
   {
      if (!sign) KSMul_bytes(res, data1, data2, orig_length, length2, coeff_bits, log_length);
      else KSMul_bytes_signed(res, data1, data2, orig_length, length2, coeff_bits, log_length);
      return;
   }
   
   //----------------------------------------------------------------------------
   // Start of general SS/KS code
   
   unsigned long h, i, j, k, l, m, skip, skip2; 
   unsigned long length, skip_limbs, input_limbs, orig_limbs, coeff_sign;
   
   // the Schoenhage-Strassen technique as implemented below, handles
   // the situtation where the output coefficient size is *bigger* than the 
   // input poly length, but when it is much smaller, the code will use a 
   // non-optimal coefficient size and essentially zero-pad the coefficients.
   // To get more optimal input coefficients to SS, we can use the Kronecker-
   // Schoenhage trick, which as we implement it, stacks more than one 
   // polynomial coefficient into the Schoenhage-Strassen input coefficients 
   // at once. This was suggested to us by Paul Zimmerman and independently 
   // derived by David Harvey.

   
   // First we determine how many coefficients should be bundled together
   // we start with a default of 1 and then compute what is optimal
   unsigned long bundle; 
   
   // number of smaller coefficients into which each original coefficient will be split
   unsigned long split;
   
   // Set to 1 if we should pack coefficients in bundles, down to the bit
   int bit_pack = 0; 
   
   // If we should bit pack
   if (output_bits < FLINT_BITS_PER_LIMB)
   {
      KSHZ_mul(res, data1, data2, orig_length, length2, coeff_bits, log_length, sign);
      return;
   } 
   
   split = get_split_tuning_parameters(&input_limbs, &orig_limbs, &output_bits,
                    &log_length, coeff_bits, orig_length, sign);
      
   // if we should not split, check if we should bundle
   if (split == 1)
   bundle = get_bundle_tuning_parameters(&skip_limbs, &output_bits, 
                        &input_limbs, &log_length, &length, orig_length, 
                        coeff_bits, sign);

   // if we aren't bit packing, splitting or bundling, put everything back
   // how it should be, assuming one input coefficient to each fft coefficient
   if ((bundle == 1) && (split == 1))
   {
      length = orig_length;
      log_length = 0;
      while ((1<<log_length) < orig_length) log_length++;
      if (sign) output_bits = 2*coeff_bits + log_length + 2;
      else output_bits = 2*coeff_bits + log_length;
      // Round up to a multiple of the minimum coefficient size supported by a 
      // transform able to deal with that length input adjusting for sqrt(2) trick
#if USE_TRUNCATED_FFT
      if (2*coeff_bits <= orig_length)
      output_bits = (((output_bits - 1) >> (log_length-1)) + 1) << (log_length-1);
      else 
#endif
      output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
   
   } else
   // Round up to a multiple of the minimum coefficient size supported by a 
   // transform able to deal with that length input adjusting for sqrt(2) trick
   output_bits = (((output_bits - 1) >> log_length) + 1) << log_length;
   
   // make nB larger than the output coefficient size (B = FLINT_BITS_PER_LIMB)
   // this effectively rounds the output coefficients up to a multiple of the 
   // limb length
   unsigned long n = (output_bits - 1) / FLINT_BITS_PER_LIMB + 1;
   
   // now we set length to the *output* polynomial length (deg+1)
   // since the input polynomial length is not needed further
   length = 1 << (log_length + 1);
   
#if USE_TRUNCATED_FFT
   if ((bundle == 1) && (split == 1)) trunc_length = orig_length; 
   if (bundle != 1) trunc_length = (orig_length - 1)/bundle + 1; 
   if (split != 1) trunc_length = split*orig_length; 
#endif
   
   // We want p = 2^rm+1 to be just larger than the largest output
   // coefficient and where m=2^k is bigger than B = FLINT_BITS_PER_LIMB
   // and in such a way that p will have 2l roots of unity where l is the
   // length of the input polynomials (i.e. l = deg+1 of input polys), so that
   // we can do a convolution of length 2l. In other words: k >= l and 2^k > B. 
   // But the output coefficients are smaller than 2^nB+1 by definition, 
   // and we already ensured nB was a multiple of l, so these two 
   // conditions are satisfied if we set r = nB / (length/2).
   // (Recall l = length/2).
   unsigned long r = n*FLINT_BITS_PER_LIMB * 2 / length;

   // We need the following memory allocated:
   //
   // * Space for 2*length + 1 coefficients of size n+1 limbs
   //   (the extra 20 is for the scratch buffer)
   //
   // * Space for 2*length + 1 pointers to coefficients
   //
   // * Working space for two pointwise multiplications of length 2n
   //
   //   We use n+1 limbs instead of n so that we can allow
   //   carries to build up without reducing mod p all the time
   
   mp_limb_t* array = (mp_limb_t*) limb_alloc((2*length + 16) * (n+1) + 2*n*THREADS,0);
   mp_limb_t** pointarr = (mp_limb_t**) limb_alloc(2*length + 16*16,0);

   // partition up array into coefficient buffers of length n+1...
   
   pointarr[0] = array;
   for (i = 0; i < 2*length-1; i++)
      pointarr[i+1] = pointarr[i] + (n+1);
   for (i = 0; i < 16; i++)
      pointarr[2*length+16*i] = array+(2*length+i)*(n+1);
      
   // ...and pointwise multiplication working space
   
   mp_limb_t* array3 = array + (2*length+16)*(n+1);   

   mp_limb_t** array1 = pointarr;
   mp_limb_t** array2 = pointarr + length;
   mp_limb_t** scratch = pointarr + 2*length;
   
   // output_bits: number of bits required for each fft output coefficient (multiple of B)
   // input_limbs: number of limbs per fft coeff
   // orig_limbs: orig no. limbs per full orig coefficient
   // orig_length: length of orig input poly
   // log_length: log of the new total fft input length
   // length: 2^(log_length+1)
   // n: number of limbs per fft output coefficient (and spacing for input coeffs)
   // split: number of pieces we split each coeff into
   
   // Convert data in from mpz_t's
   ssmul_convert_in(array1, data1, length, sign, bit_pack, split, bundle, orig_length, 
                          n, input_limbs, skip_limbs, orig_output_bits, coeffs_per_limb);
          
   ssmul_convert_in(array2, data2, length, sign, bit_pack, split, bundle, length2,
                          n, input_limbs, skip_limbs, orig_output_bits, coeffs_per_limb);
             
   // Do main ssmul (FFT's, pointwise mults, IFFT)
   
#if USE_TRUNCATED_FFT
   ssmul_main(array1, trunc_length, array2, trunc_length, scratch, log_length+1, n);
#else
   ssmul_main_old(array1, array2, array3, length, log_length, n, r, scratch);
#endif

   // We have to clear the output polynomial, since we have to add to its
   // coefficients rather than just copy into them
   
   for (i = 0; i < orig_length + length2 - 1; i++) mpz_set_ui(res[i],0);
    
   // Convert data out to mpz_t's
   if (!sign)
   {
      if (split != 1) ssmul_split_convert_out(res, array1, input_limbs, 
                          orig_length, length2, length, log_length, n, orig_limbs, split);
      else if (bundle == 1) ssmul_convert_out(res, array1, orig_length, length2, length,  
                                                                  log_length, n);
      else ssmul_bundle_convert_out(res, array1, orig_length, length2, length, log_length, 
                                                          n, skip_limbs, bundle);
   } else
   {
      if (split != 1) ssmul_split_convert_out_signed(res, array1, input_limbs, 
                            orig_length, length2, length, log_length, n, orig_limbs, split);
      else if (bundle == 1) ssmul_convert_out_signed(res, array1, orig_length, length2, 
                                                            length, log_length, n);
      else ssmul_bundle_convert_out_signed(res, array1, orig_length, length2, length, 
                                                log_length, n, skip_limbs, bundle);
   }
   limb_release();
   limb_release();
}

// end of file ****************************************************************
