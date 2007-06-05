/****************************************************************************

ZmodFpoly.c

Polynomials over Z/pZ, where p = the Fermat number B^n + 1, where
B = 2^FLINT_BITS_PER_LIMB. Routines for truncated Schoenhage-Strassen FFTs
and convolutions.

Copyright (C) 2007, William Hart and David Harvey

*****************************************************************************/

#include <sys/types.h>
#include "flint.h"
#include "memory-manager.h"
#include "ZmodFpoly.h"
#include "ZmodF_mul.h"
#include "Zpoly_mpn.h"
#include "extras.h"
#include "mpn_extras.h"

/****************************************************************************

   Memory Management Routines
   
****************************************************************************/

void ZmodFpoly_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;

   poly->storage = (mp_limb_t*) flint_heap_alloc(bufs * (n+1));
     
   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) flint_heap_alloc_bytes(bufs*sizeof(ZmodF_t));
   
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   for (unsigned long i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}


void ZmodFpoly_clear(ZmodFpoly_t poly)
{
   flint_heap_free(poly->coeffs);
   flint_heap_free(poly->storage);
}

void ZmodFpoly_stack_init(ZmodFpoly_t poly, unsigned long depth, unsigned long n,
                    unsigned long scratch_count)
{
   poly->n = n;
   poly->depth = depth;
   poly->scratch_count = scratch_count;
   poly->length = 0;
   
   unsigned long bufs = (1 << depth) + scratch_count;

   poly->storage = (mp_limb_t*) flint_stack_alloc(bufs * (n+1));
     
   // put scratch array immediately after coeffs array
   poly->coeffs = (ZmodF_t*) flint_stack_alloc_bytes(bufs*sizeof(ZmodF_t));
   
   poly->scratch = poly->coeffs + (1 << depth);
   
   poly->coeffs[0] = poly->storage;
   for (unsigned long i = 1; i < bufs; i++)
      poly->coeffs[i] = poly->coeffs[i-1] + (n+1);
}

void ZmodFpoly_stack_clear(ZmodFpoly_t poly)
{
   flint_stack_release();
   flint_stack_release();
}



/****************************************************************************

   Conversion Routines
   
****************************************************************************/

long ZmodFpoly_convert_in_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn)
{
   unsigned long size_f = poly_f->n + 1;
   unsigned long size_m = poly_mpn->limbs+1;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;
   
   unsigned long mask = -1L;
   long bits = 0;
   long limbs = 0;
   long sign = 1;
   
   long size_j;
   
   for (unsigned long i = 0, j = 0; i < poly_mpn->length; i++, j += size_m)
   {
      size_j = coeffs_m[j];
      if ((long) size_j < 0) sign = -1L;
      if (ABS(size_j) > limbs + 1)
      {
         limbs = ABS(size_j) - 1;
         bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]); 
         if (bits == FLINT_BITS_PER_LIMB) mask = 0L;
         else mask = -1L - ((1L<<bits)-1);  
      } else if (ABS(size_j) == limbs + 1)
      {
         if (coeffs_m[j+ABS(size_j)] & mask)
         {
            bits = FLINT_BIT_COUNT(coeffs_m[j+ABS(size_j)]);   
            if (bits == FLINT_BITS_PER_LIMB) mask = 0L;
            else mask = -1L - ((1L<<bits)-1);
         }
      }
      if (size_j < 0)
      {
         negate_limbs(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         set_limbs(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      } else
      {
         copy_limbs(coeffs_f[i], coeffs_m + j + 1, ABS(size_j)); 
         clear_limbs(coeffs_f[i] + ABS(size_j), size_f - ABS(size_j)); 
      }
   }
   poly_f->length = poly_mpn->length; 
   
   return sign*(FLINT_BITS_PER_LIMB*limbs+bits);  
}

void ZmodFpoly_convert_out_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, long sign)
{
   unsigned long n = poly_f->n;
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long limbs = FLINT_MIN(n, size_m-1);
   
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   ZmodF_t * coeffs_f = poly_f->coeffs;

   if (sign)
   {
      for (unsigned long i = 0, j = 0; i < poly_f->length; i++, j += size_m)
      {
         ZmodF_normalise(coeffs_f[i], n);
         if (coeffs_f[i][n-1]>>(FLINT_BITS_PER_LIMB-1) || coeffs_f[i][n])
         {
            negate_limbs(coeffs_m + j + 1, coeffs_f[i], limbs);
            mpn_add_1(coeffs_m + j + 1, coeffs_m + j + 1, limbs, 1L);
            coeffs_m[j] = -limbs;
            NORM(coeffs_m + j);
         } else
         {
            copy_limbs(coeffs_m + j + 1, coeffs_f[i], limbs);
            coeffs_m[j] = limbs;
            NORM(coeffs_m + j);         
         }
      }
   } else 
   {
      for (unsigned long i = 0, j = 0; i < poly_f->length; i++, j += size_m)
      {
         ZmodF_normalise(coeffs_f[i], n);
         copy_limbs(coeffs_m + j + 1, coeffs_f[i], limbs);
         coeffs_m[j] = limbs;
         NORM(coeffs_m + j);         
      }
   }
   
   poly_mpn->length = poly_f->length;   

}

static inline long __get_next_coeff(mp_limb_t * coeff_m, long * borrow, long * coeff, long mask)
{ 
   if ((long) coeff_m[0] == 0) *coeff = -*borrow;
   else if ((long) coeff_m[0] > 0) *coeff = coeff_m[1] - *borrow;
   else *coeff = (-coeff_m[1] - *borrow);
   *borrow = 0UL;
   if (*coeff < 0) 
   {
      *borrow = 1UL;
   }
   *coeff&=mask;
   
   return *coeff;
}

static inline long __get_next_coeff_unsigned(mp_limb_t * coeff_m, long * coeff)
{ 
   if ((long) coeff_m[0] == 0) *coeff = 0;
   else *coeff = coeff_m[1];
   
   return *coeff;
}

void ZmodFpoly_bit_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                            unsigned long bundle, long bits)
{   
   unsigned long i, k, skip;
   unsigned long n = poly_f->n;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * array, * next_point;
    
   unsigned long temp;
   half_ulong lower;
   long coeff;
   long borrow;
   mp_limb_t extend;
   
   int sign = (bits < 0);
   if (sign) bits = ABS(bits);
   
   unsigned long coeffs_per_limb = FLINT_BITS_PER_LIMB/bits;

   const unsigned long mask = (1UL<<bits)-1;
      
   poly_f->length = 0;
   i=0;
      
   while (coeff_m < poly_mpn->coeffs + 2*poly_mpn->length)
   {
      k=0; skip=0;
      coeff = 0; borrow = 0L; temp = 0;
      array = poly_f->coeffs[i];
      i++;
   
      next_point = coeff_m + 2*bundle;
      if (next_point >= poly_mpn->coeffs + 2*poly_mpn->length) next_point = poly_mpn->coeffs + 2*poly_mpn->length;
      else for (unsigned long j = 0; j < n; j += 8) FLINT_PREFETCH(poly_f->coeffs[i+1], j);
         
      while (coeff_m < next_point)
      {
         if ((unsigned long)coeff_m&7 == 0) FLINT_PREFETCH(coeff_m,64);
         // k is guaranteed to be less than FLINT_BITS_PER_LIMB at this point
         while ((k<HALF_FLINT_BITS_PER_LIMB)&&(coeff_m < next_point))
         {
            if (sign) temp+=(__get_next_coeff(coeff_m, &borrow, &coeff, mask) << k);
            else temp+=(__get_next_coeff_unsigned(coeff_m, &coeff) << k);
            coeff_m+=2; k+=bits;
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

               while ((k<HALF_FLINT_BITS_PER_LIMB)&&(coeff_m < next_point))
               {
                  if (sign) temp+=(__get_next_coeff(coeff_m, &borrow, &coeff, mask) << k);
                  else temp+=(__get_next_coeff_unsigned(coeff_m, &coeff) << k);
                  coeff_m+=2; k+=bits;
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
            } // if
         } // else
         poly_f->length++;
      } // while

      // sign extend the last FLINT_BITS_PER_LIMB bits we write out
      if (skip < n)
      {
        if (borrow) temp+=(-1UL << k);
        array[skip] = temp;
        skip++;
      } 
      // sign extend the remainder of the array, reducing modulo p 
      extend = 0;
      if (borrow) extend = -1L;
      while (skip < n+1) 
      {
         array[skip] = extend;
         skip++;
      }
      
   } // while
#if DEBUG
   for (unsigned long i = 0; i < n+1; i++)
   {
       printf("%lx ",poly_f->coeffs[0][i]);
   }
   printf("\n");
#endif

}


void ZmodFpoly_bit_unpack_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, 
                              unsigned long bundle, unsigned long bits)
{
   unsigned long k, skip;

   unsigned long temp2;
   unsigned long temp;
   unsigned long full_limb;
   unsigned long carry;
    
   mp_limb_t* array;
    
   const unsigned long mask = (1UL<<bits)-1;
   const unsigned long sign_mask = (1UL<<(bits-1));

   unsigned long s;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * next_point;
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long n = poly_f->n;
   
#if DEBUG
   for (unsigned long i = 0; i < n+1; i++)
   {
       printf("%lx ",poly_f->coeffs[0][i]);
   }
   printf("\n");
#endif
    
   for (unsigned long i = 0; coeff_m < poly_mpn->coeffs + poly_mpn->length*size_m; i++)
   {
      array = poly_f->coeffs[i];
      ZmodF_normalise(array, n);

      k=0; skip=0; carry = 0UL; temp2 = 0;
      next_point = coeff_m + size_m*bundle;
      if (next_point >= poly_mpn->coeffs + poly_mpn->length*size_m) next_point = poly_mpn->coeffs + poly_mpn->length*size_m;
      else for (unsigned long j = 0; j < n; j += 8) FLINT_PREFETCH(poly_f->coeffs[i+1], j);
      
      while (coeff_m < next_point)
      {
         // read in a full limb
         full_limb = array[skip];
         temp2 += l_shift(full_limb,k);
         s=FLINT_BITS_PER_LIMB-k;
         k+=s;
         while ((k >= bits)&&(coeff_m < next_point))
         {
            if (!(temp2&sign_mask)) 
            {
               __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask)+carry);
               carry = 0UL;
            }  
            else
            {
               temp = ((-temp2)&mask)-carry;
               __Zpoly_mpn_sub_coeff_ui(coeff_m, temp);
               carry = 1UL;
            }
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         // read in remainder of full_limb
         temp2 += l_shift(r_shift(full_limb,s),k);
         k+=(FLINT_BITS_PER_LIMB-s);
       
         while ((k >= bits)&&(coeff_m < next_point))
         {
            if (!(temp2&sign_mask)) 
            {
               __Zpoly_mpn_add_coeff_ui(coeff_m, (temp2&mask)+carry);
               carry = 0UL;
            }
            else
            {
               temp = ((-temp2)&mask)-carry;
               __Zpoly_mpn_sub_coeff_ui(coeff_m, temp);
               carry = 1UL;
            }
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         skip++;
      }
   }
}
 
void ZmodFpoly_bit_unpack_unsigned_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, 
                              unsigned long bundle, unsigned long bits)
{
   unsigned long k, l, skip;

   unsigned long temp2;
   unsigned long temp;
   unsigned long full_limb;
    
   mp_limb_t* array;
    
   const unsigned long mask = (1UL<<bits)-1;

   unsigned long s;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * next_point;
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long n = poly_f->n;
       
   for (unsigned long i = 0; coeff_m < poly_mpn->coeffs + poly_mpn->length*size_m; i++)
   {
      array = poly_f->coeffs[i];
      
      ZmodF_normalise(array, n);

      k=0; skip=0; temp2 = 0;
      next_point = coeff_m + size_m*bundle;
      if (next_point >= poly_mpn->coeffs + poly_mpn->length*size_m) next_point = poly_mpn->coeffs + poly_mpn->length*size_m;
      else for (unsigned long j = 0; j < n; j += 8) FLINT_PREFETCH(poly_f->coeffs[i+1], j);
      
      while (coeff_m < next_point)
      {
         if (skip&7 == 0) FLINT_PREFETCH(array+skip,64);
         // read in a full limb
         full_limb = array[skip];
         temp2 += l_shift(full_limb,k);
         s=FLINT_BITS_PER_LIMB-k;
         k+=s;
         while ((k >= bits)&&(coeff_m < next_point))
         {
            __Zpoly_mpn_add_coeff2_ui(coeff_m, (temp2&mask));
            coeff_m += size_m;
            temp2>>=bits;
            k-=bits;
         }
         // k is now less than bits
         // read in remainder of full_limb
         temp2 += l_shift(r_shift(full_limb,s),k);
         k+=(FLINT_BITS_PER_LIMB-s);
       
         while ((k >= bits)&&(coeff_m < next_point))
         {
            __Zpoly_mpn_add_coeff2_ui(coeff_m, temp2&mask);
            coeff_m += size_m;
            temp2>>=bits;
            l++;
            k-=bits;
         }
         // k is now less than bits
         skip++;
      }
   }
}

void ZmodFpoly_limb_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                                           unsigned long bundle, long limbs)
{
   unsigned long size_m = poly_mpn->limbs + 1;
   long size_j;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   mp_limb_t * coeffs_f = poly_f->coeffs[0];
   long carry = 0;
   
   for (unsigned long i = 0, j = 0, k = 0; i < bundle; i++, j += size_m, k += limbs)
   {
      size_j = (long) coeffs_m[j];
      if (size_j < 0)
      {
         negate_limbs(coeffs_f + k, coeffs_m + j + 1, ABS(size_j)); 
         set_limbs(coeffs_f + k + ABS(size_j), limbs - ABS(size_j));
         if (carry) mpn_sub_1(coeffs_f + k, coeffs_f + k, limbs, 1L);
         carry = 1L;
      } else if (size_j > 0)
      {
         copy_limbs(coeffs_f + k, coeffs_m + j + 1, ABS(size_j)); 
         clear_limbs(coeffs_f + k + ABS(size_j), limbs - ABS(size_j)); 
         if (carry) mpn_sub_1(coeffs_f + k, coeffs_f + k, limbs, 1L);
         carry = 0L;
      } else
      {
         if (carry) 
         {
            set_limbs(coeffs_f + k, limbs);
            carry = 1L;
         } else 
         {
            clear_limbs(coeffs_f + k, limbs);
            carry = 0L;
         }
      }
   }
}

void ZmodFpoly_limb_unpack_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, 
                                  unsigned long bundle, unsigned long limbs)
{
   unsigned long size_m = poly_mpn->limbs + 1;
   unsigned long n = poly_f->n;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   mp_limb_t * coeffs_f = poly_f->coeffs[0];
   unsigned long carry = 0L;
   
   for (unsigned long i = 0, j = 0, k = 0; i < bundle; i++, j += size_m, k += limbs)
   {
      if (carry) mpn_add_1(coeffs_f + k, coeffs_f + k, n - k, 1L);
      if (coeffs_f[k+limbs-1]>>(FLINT_BITS_PER_LIMB-1))
      {
         negate_limbs(coeffs_m + j + 1, coeffs_f + k, limbs);
         coeffs_m[j] = -limbs;
         NORM(coeffs_m + j);
         carry = 1L;
      } else
      {
         copy_limbs(coeffs_m + j + 1, coeffs_f + k, limbs);
         coeffs_m[j] = limbs;
         NORM(coeffs_m + j); 
         carry = 0L;        
      }
   }
}

void ZmodFpoly_limb_unpack_unsigned_mpn(Zpoly_mpn_t poly_mpn, ZmodFpoly_t poly_f, 
                                  unsigned long bundle, unsigned long limbs)
{
   unsigned long size_m = poly_mpn->limbs + 1;
   unsigned long n = poly_f->n;
   mp_limb_t * coeffs_m = poly_mpn->coeffs;
   mp_limb_t * coeffs_f = poly_f->coeffs[0];
   
   for (unsigned long i = 0, j = 0, k = 0; i < bundle; i++, j += size_m, k += limbs)
   {
         copy_limbs(coeffs_m + j + 1, coeffs_f+k, limbs);
         coeffs_m[j] = limbs;
         NORM(coeffs_m + j);          
   }
   
}

void __ZmodFpoly_write_next_limb(mp_limb_t * array, unsigned long * temp, unsigned long * offset_limb, 
             unsigned long next_limb, unsigned long shift_1, unsigned long shift_2)
{
   *temp += l_shift(next_limb, shift_1);
   array[*offset_limb] = *temp + ((l_shift(1UL,shift_1)-1)&array[*offset_limb]);
   (*offset_limb)++;
   *temp = r_shift(next_limb, shift_2);
}

void __ZmodFpoly_write_whole_limb(mp_limb_t * array, unsigned long * temp, unsigned long * offset_limb, 
             unsigned long next_limb, unsigned long shift_1, unsigned long shift_2)
{
   *temp += l_shift(next_limb,shift_1);
   array[*offset_limb] = *temp;
   (*offset_limb)++;
   *temp = r_shift(next_limb,shift_2);
}

void ZmodFpoly_byte_pack_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                             unsigned long bundle, unsigned long coeff_bytes)
{
   unsigned long size_m = poly_mpn->limbs+1;
   unsigned long total_limbs = poly_f->n + 1;
   mp_limb_t * coeff_m = poly_mpn->coeffs;
   mp_limb_t * array;
    
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   // Start limb of the current coefficient within array
   unsigned long coeff_limb;
   // Additional offset in bytes of current coefficient within array
   unsigned long coeff_byte;
    
   // Where we are up to in the current coefficient: limbs + bytes
   unsigned long offset_limb;
    
   // Next limb to be written to bytes
   unsigned long next_limb;
   unsigned long temp;
   unsigned long extend;

   unsigned long shift_1, shift_2;
    
   // when a coefficient is negative, we need to borrow from the next coefficient
   int borrow;
   int borrowed;
    
   unsigned long i, j;
   
   mp_limb_t * end = poly_mpn->coeffs + size_m*poly_mpn->length;
    
   i = 0;
   poly_f->length = 0;
   
   while (coeff_m < end)
   {
      coeff_limb = 0;
      coeff_byte = 0;
      offset_limb = 0;
      temp = 0;
      borrow = 0;
      borrowed = 0;
            
      array = poly_f->coeffs[i];
      i++;

      mp_limb_t * next_point = coeff_m + size_m*bundle;
      if (next_point > end) next_point = end;
         
      while (coeff_m < next_point)
      {
          // compute shifts to be used
          shift_1 = coeff_byte<<3;
          shift_2 = FLINT_BITS_PER_LIMB-shift_1;
        
          borrowed = borrow;
             
          if (borrow) __Zpoly_mpn_sub_coeff_ui(coeff_m, 1);
          
          /* Coefficient is negative after borrow */
          if ((long) coeff_m[0] < 0) 
          {
             // mpz_t's store the absolute value only, so add 1 then complement
             __Zpoly_mpn_add_coeff_ui(coeff_m, 1);
             
             // deal with first limb of coefficient
             next_limb = ~coeff_m[1];
             if (limbs_per_coeff == 0) 
             {
                if (coeff_m == next_point-size_m) 
                {
                   __ZmodFpoly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                   temp += l_shift(-1UL,shift_1);
                   array[offset_limb] = temp;
                   offset_limb++;
                   extend = -1L;
                } else 
                {
                   next_limb &= ((1UL<<(extra_bytes_per_coeff<<3))-1);
                   __ZmodFpoly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                   array[offset_limb] = temp;
                }
             } else
             {
                __ZmodFpoly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                // deal with remaining limbs
                for (j = 1; j < ABS(coeff_m[0]); j++)
                {
                   next_limb = ~coeff_m[j+1];
                   __ZmodFpoly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
                }
                // write remaining part of coefficient and fill 
                // remainder of coeff_bytes with binary 1's
                if ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                    limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte) 
                {
                   temp += l_shift(-1UL,shift_1);
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
                extend = -1L;
             }
             temp = 0;
             // restore the original coefficient
             __Zpoly_mpn_sub_coeff_ui(coeff_m, 1);
             borrow = 1;
          }
          
          /* Coefficient is positive after borrow */
          else if ((long) coeff_m[0] > 0)
          {
             // deal with first limb of coefficient
             next_limb = coeff_m[1];
             __ZmodFpoly_write_next_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             if (shift_2 == FLINT_BITS_PER_LIMB) temp = 0;
             // deal with remaining limbs
             for (j = 1; j < ABS(coeff_m[0]); j++)
             {
                next_limb = coeff_m[j+1];
                __ZmodFpoly_write_whole_limb(array, &temp, &offset_limb, next_limb, shift_1, shift_2);
             }
             // write remaining part of coefficient
             array[offset_limb] = temp;
             offset_limb++;
             for (; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = 0UL;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                    limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = 0UL;
                offset_limb++;
             }
             extend = 0L;
             temp = 0;
             borrow = 0;
          }
          /* Coefficient is zero after borrow */
          else
          {
             array[offset_limb] = ((l_shift(1UL,shift_1)-1)&array[offset_limb]);
             offset_limb++;
             for (; offset_limb < coeff_limb + limbs_per_coeff; offset_limb++)
             {
                array[offset_limb] = 0UL;
             }
             while ((offset_limb<<FLINT_LG_BYTES_PER_LIMB) < ((coeff_limb +
                    limbs_per_coeff)<<FLINT_LG_BYTES_PER_LIMB)+extra_bytes_per_coeff + coeff_byte)
             {
                array[offset_limb] = 0UL;
                offset_limb++;
             }

             extend = 0L;
             temp = 0;
             borrow = 0;
          }
          
          // Extend sign of final input coefficient to end of output coefficient 
          if (coeff_m == next_point-size_m)
          {
             while (offset_limb < total_limbs)
             {
                array[offset_limb] = extend;
                offset_limb++;
             } 
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
          if (borrowed) __Zpoly_mpn_add_coeff_ui(coeff_m, 1);
          coeff_m += size_m;
      }
      poly_f->length++;
   }
}
     
static inline void __ZmodFpoly_unpack_bytes(mp_limb_t* output, mp_limb_t* array, 
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
    
    
    shift_1 = (byte_start<<3);
    shift_2 = FLINT_BITS_PER_LIMB - shift_1;
    
    temp = array[coeff_limb];
    coeff_limb++;
    while (output_limb < limbs_to_extract)
    {
       next_limb = r_shift(temp,shift_1);
       temp = array[coeff_limb];
       coeff_limb++;
       next_limb += l_shift(temp,shift_2);
       output[output_limb] = next_limb;
       output_limb++;
    }
    if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
    {
       next_limb = r_shift(temp,shift_1);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    } else
    {
       next_limb = r_shift(temp,shift_1);
       temp = array[coeff_limb];
       next_limb += l_shift(temp,shift_2);
       output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
    }
}

static inline unsigned long __ZmodFpoly_unpack_signed_bytes(mp_limb_t* output, mp_limb_t* array, 
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
           next_limb = r_shift(temp,shift_1);
           temp = ~array[coeff_limb];
           coeff_limb++;
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
        {
           next_limb = r_shift(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else
        {
           next_limb = r_shift(temp,shift_1);
           temp = ~array[coeff_limb];
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    else
    {
        temp = array[coeff_limb];
        coeff_limb++;
        while (output_limb < limbs_to_extract)
        {
           next_limb = r_shift(temp,shift_1);
           temp = array[coeff_limb];
           coeff_limb++;
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb;
           output_limb++;
        }
        if (extra_bytes_to_extract <= FLINT_BYTES_PER_LIMB - byte_start)
        {
           next_limb = r_shift(temp,shift_1);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        } else
        {
           next_limb = r_shift(temp,shift_1);
           temp = array[coeff_limb];
           next_limb += l_shift(temp,shift_2);
           output[output_limb] = next_limb&((1UL<<(extra_bytes_to_extract<<3))-1);
        }
    }
    return sign;
}   

void ZmodFpoly_byte_unpack_unsigned_mpn(Zpoly_mpn_t poly_m, mp_limb_t* array,
                               unsigned long bundle, unsigned long coeff_bytes)
{
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   const unsigned long limbs = ((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
   
   mp_limb_t* temp = (mp_limb_t*) flint_stack_alloc(limbs+2);
   
   unsigned long limb_upto = 0;
   unsigned long byte_offset = 0;
   
   mp_limb_t * coeff_m = poly_m->coeffs;
   unsigned long size_m = poly_m->limbs+1;
   poly_m->length = bundle;
   
   for (unsigned long i = 0; i < bundle; i++)
   {
       clear_limbs(temp, limbs+2);
       __ZmodFpoly_unpack_bytes(temp + 1, array, limb_upto, 
                                             byte_offset, coeff_bytes);
       temp[0] = limbs;
       NORM(temp);
       
       __Zpoly_mpn_add_coeffs(coeff_m, coeff_m, temp);
       
       limb_upto += limbs_per_coeff;
       
       if (byte_offset + extra_bytes_per_coeff >= FLINT_BYTES_PER_LIMB)
       {
          limb_upto++;
          byte_offset = byte_offset + extra_bytes_per_coeff - FLINT_BYTES_PER_LIMB;
       } else
       {
          byte_offset = byte_offset + extra_bytes_per_coeff;
       }
       coeff_m += size_m;
   }                        
   flint_stack_release();
}

void ZmodFpoly_byte_unpack_mpn(Zpoly_mpn_t poly_m, mp_limb_t* array,
                               unsigned long bundle, unsigned long coeff_bytes)
{
   const unsigned long limbs_per_coeff = (coeff_bytes>>FLINT_LG_BYTES_PER_LIMB);
   const unsigned long extra_bytes_per_coeff = coeff_bytes 
                                    - (limbs_per_coeff<<FLINT_LG_BYTES_PER_LIMB);
                                    
   const unsigned long limbs = ((coeff_bytes-1)>>FLINT_LG_BYTES_PER_LIMB) + 1;
   
   mp_limb_t* temp = (mp_limb_t*) flint_stack_alloc(limbs+2);
   
   unsigned long limb_upto = 0;
   unsigned long byte_offset = 0;
   
   unsigned long sign;
   unsigned long borrow = 0;
   
   mp_limb_t * coeff_m = poly_m->coeffs;
   unsigned long size_m = poly_m->limbs+1;
   poly_m->length = bundle;
   
   for (unsigned long i = 0; i < bundle; i++)
   {
       clear_limbs(temp,limbs+2);
       sign = __ZmodFpoly_unpack_signed_bytes(temp + 1, array, limb_upto, 
                                             byte_offset, coeff_bytes);
       if (sign) temp[0] = -limbs;
       else temp[0] = limbs;
       NORM(temp);
       
       if (sign)
       {
          __Zpoly_mpn_sub_coeff_ui(temp, 1);
       }
       if (borrow) __Zpoly_mpn_add_coeff_ui(temp, 1);
       
       __Zpoly_mpn_add_coeffs(coeff_m, coeff_m, temp);
       
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
       coeff_m += size_m;
   }                        
   flint_stack_release();
}


     
void ZmodFpoly_split_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                         unsigned long bundle, unsigned long limbs)
{
   abort();
}
     
void ZmodFpoly_unsplit_mpn(ZmodFpoly_t poly_f, Zpoly_mpn_t poly_mpn,
                           unsigned long bundle, unsigned long limbs)
{
   abort();
}
     

/****************************************************************************

   Basic Arithmetic Routines
   
****************************************************************************/


void ZmodFpoly_set(ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->n == y->n);

   for (unsigned long i = 0; i < y->length; i++)
      ZmodF_set(x->coeffs[i], y->coeffs[i], x->n);

   x->length = y->length;
}


void ZmodFpoly_pointwise_mul(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   unsigned long j;

   ZmodF_mul_precomp_t info;
   ZmodF_mul_precomp_init(info, x->n, x == y);
   
   if (x != y)
      for (unsigned long i = 0; i < x->length; i++)
      {
         for (j = 0; j < x->n; j += 8) FLINT_PREFETCH(x->coeffs[i+8], j);
         for (j = 0; j < y->n; j += 8) FLINT_PREFETCH(y->coeffs[i+8], j);
         ZmodF_mul_precomp(info, res->coeffs[i], x->coeffs[i], y->coeffs[i]);
      }
   else
      for (unsigned long i = 0; i < x->length; i++)
      {
         for (j = 0; j < x->n; j += 8) FLINT_PREFETCH(x->coeffs[i+8], j);
         ZmodF_sqr_precomp(info, res->coeffs[i], x->coeffs[i]);
      }

   ZmodF_mul_precomp_clear(info);

   res->length = x->length;
}


void ZmodFpoly_add(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   for (unsigned long i = 0; i < x->length; i++)
      ZmodF_add(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodFpoly_sub(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);
   FLINT_ASSERT(x->length == y->length);
   
   for (unsigned long i = 0; i < x->length; i++)
      ZmodF_sub(res->coeffs[i], x->coeffs[i], y->coeffs[i], x->n);

   res->length = x->length;
}


void ZmodFpoly_normalise(ZmodFpoly_t poly)
{
   for (unsigned long i = 0; i < poly->length; i++)
      ZmodF_normalise(poly->coeffs[i], poly->n);
}


void ZmodFpoly_rescale(ZmodFpoly_t poly)
{
   if (poly->depth == 0)
      return;

   for (unsigned long i = 0; i < poly->length; i++)
      ZmodF_short_div_2exp(poly->coeffs[i], poly->coeffs[i],
                           poly->depth, poly->n);
}


/****************************************************************************

   Forward fourier transforms (internal code)

****************************************************************************/


void _ZmodFpoly_FFT_iterative(
            ZmodF_t* x, unsigned long depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   FLINT_ASSERT(depth >= 1);

   unsigned long i, s, start;
   ZmodF_t* y, * z;

   // root is the (2^depth)-th root of unity for the current layer,
   // measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);

   // half = half the current block length
   unsigned long half = 1UL << (depth - 1);
   unsigned long half_skip = half * skip;

   unsigned long layer;

   // =========================================================================
   // Special case for first layer, if root and/or twist involve sqrt2
   // rotations.

   if ((root | twist) & 1)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length <= half)
      {
         // Only need first output for each butterfly,
         // i.e. (a, b) -> (a + b, ?)
         if (nonzero > half)
            for (i = 0, y = x; i < nonzero - half; i++, y += skip)
               ZmodF_add(y[0], y[0], y[half_skip], n);
      }
      else
      {
         // Need both outputs for each butterfly.
         if (nonzero <= half)
         {
            // The second half of each butterfly input are zeroes, so we just
            // computing (a, 0) -> (a, ra), where r is the appropriate root
            // of unity.
            for (i = 0, s = twist, y = x; i < nonzero;
                 i++, s += root, y += skip)
            {
               ZmodF_mul_sqrt2exp(y[half_skip], y[0], s, n);
            }
         }
         else
         {
            // If nonzero > half, then we need some full butterflies...
            for (i = 0, s = twist, y = x; i < nonzero - half;
                 i++, s += root, y += skip)
            {
               ZmodF_forward_butterfly_sqrt2exp(y, y + half_skip, 
                                                scratch, s, n);
            }
            // and also some partial butterflies (a, 0) -> (a, ra).
            for (; i < half; i++, s += root, y += skip)
               ZmodF_mul_sqrt2exp(y[half_skip], y[0], s, n);
         }
      }

      // Here we switch to measuring roots as powers of 2, but we also need
      // to update to the next layer's roots, and these two actions cancel
      // each other out :-)
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
      if (nonzero > 2*half)
         nonzero = 2*half;
   
      layer = 1;
   }
   else
   {
      // no special case for first layer
      layer = 0;

      // switch to measuring roots as powers of 2
      root >>= 1;
      twist >>= 1;
   }

   // =========================================================================
   // This section handles the layers where there are still zero coefficients
   // to take advantage of. In most cases, this will only happen for the
   // first layer or two, so we don't bother with specialised limbshift-only
   // code for these layers.

   // Note: from here on there are no sqrt2 rotations, and we measure all
   // roots as powers of 2.
   
   for (; (layer < depth) && (nonzero < 2*half); layer++)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length_remainder > half)
      {
         // If length overhangs by more than half the block, then we need to
         // perform full butterflies on the last block (i.e. the last block
         // doesn't get any special treatment).
         length_quantised += 2*half;
      }
      else if (length_remainder)
      {
         // If length overhangs the block by at most half the block size,
         // then we only need to compute the first output of each butterfly
         // for this block, i.e. (a, b) -> (a + b)
         if (nonzero > half)
         {
            y = x + skip * length_quantised;
            for (i = 0; i < nonzero - half; i++, y += skip)
               ZmodF_add(y[0], y[0], y[half_skip], n);
         }
      }

      if (nonzero <= half)
      {
         // If nonzero <= half, then the second half of each butterfly input
         // are zeroes, so we just computing (a, 0) -> (a, ra), where r is the
         // appropriate root of unity.
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, s = twist, z = y; i < nonzero;
                 i++, s += root, z += skip)
            {
               ZmodF_mul_2exp(z[half_skip], z[0], s, n);
            }
         }
      }
      else
      {
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            // If nonzero > half, then we need some full butterflies...
            for (i = 0, s = twist, z = y; i < nonzero - half;
                 i++, s += root, z += skip)
            {
               ZmodF_forward_butterfly_2exp(z, z + half_skip, scratch, s, n);
            }
            // and also some partial butterflies (a, 0) -> (a, ra).
            for (; i < half; i++, s += root, z += skip)
               ZmodF_mul_2exp(z[half_skip], z[0], s, n);
         }
      }

      // Update roots of unity
      twist <<= 1;
      root <<= 1;
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
      
      if (nonzero > 2*half)
         // no more zero coefficients to take advantage of:
         nonzero = 2*half;
   }

   // =========================================================================
   // Now we may assume there are no more zero coefficients.

   for (; layer < depth; layer++)
   {
      // Let length = multiple of block size plus a remainder.
      unsigned long length_quantised = length & (-2*half);
      unsigned long length_remainder = length - length_quantised;

      if (length_remainder > half)
      {
         // If length overhangs by more than half the block, then we need to
         // perform full butterflies on the last block (i.e. the last block
         // doesn't get any special treatment).
         length_quantised += 2*half;
      }
      else if (length_remainder)
      {
         // If length overhangs the block by at most half the block size,
         // then we only need to compute the first output of each butterfly
         // for this block, i.e. (a, b) -> (a + b)
         y = x + skip * length_quantised;
         for (i = 0; i < half; i++, y += skip)
            ZmodF_add(y[0], y[0], y[half_skip], n);
      }
      
      // To keep the inner loops long, we have two versions of the next loop.
      if (layer < depth/2)
      {
         // Version 1: only a few relatively long blocks.
         
         for (start = 0, y = x; start < length_quantised;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, s = twist, z = y; i < half; i++, s += root, z += skip)
               ZmodF_forward_butterfly_2exp(z, z + half_skip, scratch, s, n);
         }
      }
      else
      {
         // Version 2: lots of short blocks.
         
         // Two sub-versions, depending on whether the rotations are all by
         // a whole number of limbs.
         if ((root | twist) & (FLINT_BITS_PER_LIMB - 1))
         {
            // Version 2a: rotations still involve bitshifts.
            for (i = 0, s = twist, y = x; i < half; i++, s += root, y += skip)
               for (start = 0, z = y; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_forward_butterfly_2exp(z, z + half_skip,
                                               scratch, s, n);
               }
         }
         else
         {
            // Version 2b: rotations involve only limbshifts.
            unsigned long root_limbs = root >> FLINT_LG_BITS_PER_LIMB;

            if (twist == 0)
            {
               // special case, since ZmodF_forward_butterfly_Bexp doesn't
               // allow zero rotation count
               for (start = 0, z = x; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_simple_butterfly(z, z + half_skip, scratch, n);
               }
               i = 1;
               y = x + skip;
               s = root_limbs;
            }
            else
            {
               i = 0;
               y = x;
               s = twist >> FLINT_LG_BITS_PER_LIMB;
            }
            
            for (; i < half; i++, s += root_limbs, y += skip)
               for (start = 0, z = y; start < length_quantised;
                    start += 2*half, z += 2*half_skip)
               {
                  ZmodF_forward_butterfly_Bexp(z, z + half_skip,
                                               scratch, s, n);
               }
         }
      }

      // Update roots of unity
      twist <<= 1;
      root <<= 1;
      
      // Update block length.
      half >>= 1;
      half_skip >>= 1;
   }
}


/*
Factors FFT of length 2^depth into length 2^rows_depth and length 2^cols_depth
transforms
*/
void _ZmodFpoly_FFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(rows_depth >= 1);
   FLINT_ASSERT(cols_depth >= 1);
   
   unsigned long depth = rows_depth + cols_depth;
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);

   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long length_whole_rows = length_cols ?
                                     (length_rows + 1) : length_rows;
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // column transforms
   for (i = 0, y = x, j = twist; i < nonzero_cols; i++, y += skip, j += root)
      _ZmodFpoly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                     length_whole_rows, j, n, scratch);

   if (nonzero_rows)
   {
      for (; i < cols; i++, y += skip, j += root)
         _ZmodFpoly_FFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                        length_whole_rows, j, n, scratch);
      nonzero_cols = cols;
   }
   
   // row transforms
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, cols,
                     twist << rows_depth, n, scratch);

   if (length_cols)
      // The relevant portion of the last row:
      _ZmodFpoly_FFT(y, cols_depth, skip, nonzero_cols, length_cols,
                     twist << rows_depth, n, scratch);
}



/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of buffers assumed to be nonzero
length = number of fourier coefficients requested
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_FFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                    unsigned long nonzero, unsigned long length,
                    unsigned long twist, unsigned long n,
                    ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1 << depth));
   FLINT_ASSERT(length >= 1 && length <= (1 << depth));
   FLINT_ASSERT(depth >= 1);

   // If the data fits in L1 (2^depth coefficients of length n+1, plus a
   // scratch buffer), then use the iterative transform. Otherwise factor the
   // FFT into two chunks.
   if (depth == 1 ||
       ((1 << depth) + 1) * (n+1) <= ZMODFPOLY_FFT_FACTOR_THRESHOLD)
   {
      _ZmodFpoly_FFT_iterative(x, depth, skip, nonzero, length,
                               twist, n, scratch);
   }
   else
   {
      unsigned long rows_depth = depth >> 1;
      unsigned long cols_depth = depth - rows_depth;
      _ZmodFpoly_FFT_factor(x, rows_depth, cols_depth, skip, nonzero, length,
                            twist, n, scratch);
   }
}



/****************************************************************************

   Inverse fourier transforms (internal code)

****************************************************************************/


/*
This one is for when there is no truncation.
*/
void _ZmodFpoly_IFFT_iterative(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(depth >= 1);

   // root is the (2^(layer+1))-th root unity for each layer,
   // measured as a power of sqrt2
   long root = 2*n*FLINT_BITS_PER_LIMB;
   twist <<= (depth - 1);
   FLINT_ASSERT(twist < root);

   unsigned long half = 1;
   unsigned long half_skip = skip;
   unsigned long size = 1UL << depth;
   unsigned long layer, start, i, s;
   ZmodF_t* y, * z;
   
   // First group of layers; lots of small blocks.

   for (layer = 0; layer < depth/2; layer++)
   {
      // no sqrt2 should be involved here
      FLINT_ASSERT(!((twist | root) & 1));

      // change roots to be measured as powers of 2
      // (also this updates for the next layer in advance)
      root >>= 1;
      twist >>= 1;
      
      if ((root | twist) & (FLINT_BITS_PER_LIMB-1))
      {
         // This version allows bitshifts
         for (i = 0, y = x, s = twist; i < half; i++, s += root, y += skip)
            for (start = 0, z = y; start < size;
                 start += 2*half, z += 2*half_skip)
            {
               ZmodF_inverse_butterfly_2exp(z, z + half_skip, scratch, s, n);
            }
      }
      else
      {
         // This version is limbshifts only
         unsigned long root_limbs = root >> FLINT_LG_BITS_PER_LIMB;

         if (twist == 0)
         {
            // special case since ZmodF_inverse_butterfly_Bexp doesn't allow
            // zero rotation count
            for (start = 0, z = x; start < size;
                 start += 2*half, z += 2*half_skip)
            {
               ZmodF_simple_butterfly(z, z + half_skip, scratch, n);
            }
         
            i = 1;
            s = root_limbs;
            y = x + skip;
         }
         else
         {
            i = 0;
            s = twist >> FLINT_LG_BITS_PER_LIMB;
            y = x;
         }
         
         for (; i < half; i++, s += root_limbs, y += skip)
            for (start = 0, z = y; start < size;
                 start += 2*half, z += 2*half_skip)
            {
               ZmodF_inverse_butterfly_Bexp(z, z + half_skip, scratch, s, n);
            }
      }
      
      half <<= 1;
      half_skip <<= 1;
   }


   // Second group of layers; just a few large blocks.
   
   for (; layer < depth; layer++)
   {
      if ((root | twist) & 1)
      {
         // sqrt2 is involved. This had better be the last layer.
         FLINT_ASSERT(layer == depth - 1);
         
         for (i = 0, z = x, s = twist; i < half; i++, s += root, z += skip)
            ZmodF_inverse_butterfly_sqrt2exp(z, z + half_skip, scratch, s, n);
         
         return;
      }
      else
      {
         // Only bitshifts.

         // change roots to be measured as powers of 2
         // (also this updates for the next layer in advance)
         twist >>= 1;
         root >>= 1;
         
         for (start = 0, y = x; start < size;
              start += 2*half, y += 2*half_skip)
         {
            for (i = 0, z = y, s = twist; i < half; i++, s += root, z += skip)
               ZmodF_inverse_butterfly_2exp(z, z + half_skip, scratch, s, n);
         }
      }
   
      half <<= 1;
      half_skip <<= 1;
   }
}



/*
This one's for working in L1 when truncation is involved. It splits into
two halves.
*/
void _ZmodFpoly_IFFT_recursive(
               ZmodF_t* x, unsigned long depth, unsigned long skip,
               unsigned long nonzero, unsigned long length, int extra,
               unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   FLINT_ASSERT(depth >= 1);

   long size = 1UL << depth;

   if (length == size)
   {
      // no truncation necessary
      _ZmodFpoly_IFFT_iterative(x, depth, skip, twist, n, scratch);
      return;
   }

   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);

   long cols = size >> 1;
   long half = skip << (depth - 1);

   // symbols in the following diagrams:
   // A = fully untransformed coefficient
   // a = fully untransformed coefficient (implied zero)
   // B = intermediate coefficient
   // b = intermediate coefficient (implied zero)
   // C = fully transformed coefficient
   // c = fully transformed coefficient (implied zero)
   // ? = garbage that we don't care about
   // * = the extra C coefficient, or "?" if no extra coefficient requested
   
   // the horizontal transforms convert between B and C
   // the vertical butterflies convert between A and B

   if ((length < cols) || (length == cols && !extra))
   {
      // The input could look like one of the following:
      // CCCCAAAA      CCCCAAAA      CCCCAAaa      CCCCaaaa
      // AAAAAAaa  or  AAaaaaaa  or  aaaaaaaa  or  aaaaaaaa

      long i, last_zero_forward_butterfly, last_zero_cross_butterfly;

      if (nonzero <= cols)
      {
         i = nonzero - 1;
         last_zero_forward_butterfly = length;
         last_zero_cross_butterfly = 0;
      }
      else
      {
         i = cols - 1;
         if (nonzero > length + cols)
         {
            last_zero_forward_butterfly = nonzero - cols;
            last_zero_cross_butterfly = length;
         }
         else
         {
            last_zero_forward_butterfly = length;
            last_zero_cross_butterfly = nonzero - cols;
         }
      }
      
      ZmodF_t* y = x + skip*i;

      // First some forward butterflies ("Aa" => "B?") to make them look like:
      // CCCCAABB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAAAA??  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= last_zero_forward_butterfly; i--, y -= skip)
      {
         // (2*a0, ?) -> (a0, ?)   = (b0, ?)
         ZmodF_short_div_2exp(y[0], y[0], 1, n);
      }

      // Then some forward butterflies ("AA" => "B?") to make them look like:
      // CCCCBBBB      CCCCBBBB      CCCCBBaa      CCCCaaaa
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      for (; i >= (long)length; i--, y -= skip)
      {
         // (2*a0, 2*a1) -> (a0 + a1, ?)   = (b0, ?)
         ZmodF_add(y[0], y[0], y[half], n);
         ZmodF_short_div_2exp(y[0], y[0], 1, n);
      }

      // Transform the first row to make them look like:
      // BBBB*???      BBBB*???      BBBB*???      BBBB*???
      // AAAA????  or  AAaa????  or  aaaa??aa  or  aaaaaaaa
      if (depth > 1)
         _ZmodFpoly_IFFT_recursive(x, depth - 1, skip,
                                   (nonzero < cols) ? nonzero : cols,
                                   length, extra, twist << 1, n, scratch);
      
      // Cross butterflies ("Ba" => "A?") to make them look like:
      // BBBB*???      BBAA*???      AAAA*???      AAAA*???
      // AAAA????  or  AA??????  or  ??????aa  or  ????aaaa
      for (; i >= last_zero_cross_butterfly; i--, y -= skip)
      {
         // (b0, ?) -> (2*b0, ?)    = (2*a0, ?)
         ZmodF_add(y[0], y[0], y[0], n);
      }
         
      // Cross butterflies ("BA" => "A?") to make them look like:
      // AAAA*???      AAAA*???      AAAA*???      AAAA*???
      // ????????  or  ????????  or  ??????aa  or  ????aaaa
      for (; i >= 0; i--, y -= skip)
      {
         // (b0, 2*a1) -> (2*b0 - 2*a1, ?)     = (2*a0, ?)
         ZmodF_add(y[0], y[0], y[0], n);
         ZmodF_sub(y[0], y[0], y[half], n);
      }
   }
   else
   {
      // The input looks like one of these:
      // CCCCCCCC                   CCCCCCCC
      // AAAAaaaa (extra == 1)  or  CCCAAAaa
   
      // Transform first row (no truncation necessary) to make them look like:
      // BBBBBBBB                   BBBBBBBB
      // AAAAaaaa (extra == 1)  or  CCCAAAaa
      if (depth > 1)
         _ZmodFpoly_IFFT_iterative(x, depth - 1, skip, twist << 1, n, scratch);

      long i = cols - 1;
      unsigned long s = twist + root*i;
      ZmodF_t* y = x + skip*i;
      
      long last_zero_cross_butterfly = nonzero - cols;
      long last_cross_butterfly = length - cols;
   
      // Cross butterflies ("Ba" => "AB") to make them look like:
      // BBBBAAAA                   BBBBBBAA
      // AAAABBBB (extra == 1)  or  CCCAAABB
      for (; i >= last_zero_cross_butterfly; i--, s -= root, y -= skip)
      {
         // (b0, ?) -> (2*b0, w*b0)     = (2*a0, b1)
         ZmodF_mul_sqrt2exp(y[half], y[0], s, n);
         ZmodF_add(y[0], y[0], y[0], n);
      }
         
      // Cross butterflies ("BA" => "AB") to make them look like:
      // AAAAAAAA                   BBBAAAAA
      // BBBBBBBB (extra == 1)  or  CCCBBBBB
      for (; i >= last_cross_butterfly; i--, s -= root, y -= skip)
      {
         // (b0, 2*a1) -> (2*(b0-a1), w*(b0-2*a1))    = (2*a0, b1)
         ZmodF_sub(scratch[0], y[0], y[half], n);
         ZmodF_add(y[0], y[0], scratch[0], n);
         ZmodF_mul_sqrt2exp(y[half], scratch[0], s, n);
      }
      
      // Transform second row to make them look like:
      // AAAAAAAA                   BBBAAAAA
      // *??????? (extra == 1)  or  BBB*????
      if (depth > 1)
         _ZmodFpoly_IFFT_recursive(x + skip*cols, depth - 1, skip, cols,
                                   length - cols, extra, twist << 1, n,
                                   scratch);

      // Inverse butterflies ("BB" => "AA") to make them look like:
      // AAAAAAAA                   AAAAAAAA
      // *??????? (extra == 1)  or  AAA*????
      for (; i >= 0; i--, s -= root, y -= skip)
      {
         // (b0, b1) -> (b0 + w*b1, b0 - w*b1)    = (2*a0, 2*a1)
         ZmodF_inverse_butterfly_sqrt2exp(y, y + half, scratch, s, n);
      }
   }
}



/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_IFFT_factor(
            ZmodF_t* x, unsigned long rows_depth, unsigned long cols_depth,
            unsigned long skip, unsigned long nonzero, unsigned long length,
            int extra, unsigned long twist, unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(rows_depth >= 1);
   FLINT_ASSERT(cols_depth >= 1);

   unsigned long depth = rows_depth + cols_depth;
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   
   // root is the (2^depth)-th root unity, measured as a power of sqrt2
   unsigned long root = (4*n*FLINT_BITS_PER_LIMB) >> depth;
   FLINT_ASSERT(twist < root);
   
   unsigned long rows = 1UL << rows_depth;
   unsigned long cols = 1UL << cols_depth;

   unsigned long length_rows = length >> cols_depth;
   unsigned long length_cols = length & (cols-1);
   unsigned long nonzero_rows = nonzero >> cols_depth;
   unsigned long nonzero_cols = nonzero & (cols-1);

   unsigned long i, j;
   ZmodF_t* y;

   // row transforms for the rows where we have all fourier coefficients
   for (i = 0, y = x; i < length_rows; i++, y += (skip << cols_depth))
      _ZmodFpoly_IFFT(y, cols_depth, skip, cols, cols, 0,
                      twist << rows_depth, n, scratch);

   // column transforms where we have enough information
   for (i = length_cols, y = x + (skip * length_cols),
        j = twist + (root*length_cols);
        i < nonzero_cols; i++, y += skip, j += root)
   {
      _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                      length_rows, length_cols ? 1 : extra, j, n, scratch);
   }
   if (nonzero_rows)
      for (; i < cols; i++, y += skip, j += root)
         _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                         length_rows, length_cols ? 1 : extra, j, n, scratch);

   if (length_cols)
   {
      // a single switcheroo row transform
      _ZmodFpoly_IFFT(x + length_rows * (skip << cols_depth), cols_depth,
                      skip, (nonzero_rows ? cols : nonzero_cols),
                      length_cols, extra, twist << rows_depth, n, scratch);

      // remaining column transforms
      for (i = 0, y = x, j = twist; i < length_cols && i < nonzero_cols;
           i++, y += skip, j += root)
      {
         _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows + 1,
                         length_rows + 1, 0, j, n, scratch);
      }
      if (nonzero_rows)
      {
         for (; i < length_cols; i++, y += skip, j += root)
            _ZmodFpoly_IFFT(y, rows_depth, skip << cols_depth, nonzero_rows,
                            length_rows + 1, 0, j, n, scratch);
      }
   }
   else if (extra)
   {
      // need one extra trivial fourier coefficient
      x += length_rows * (skip << cols_depth);
      for (i = 1, y = x + skip; i < (nonzero_rows ? cols : nonzero_cols);
           i++, y += skip)
      {
         ZmodF_add(x[0], x[0], y[0], n);
      }
      ZmodF_short_div_2exp(x[0], x[0], cols_depth, n);
   }
}


/*
This is an internal function. It's just a temporary implementation so that
we can get started on higher level code. It is not optimised particularly
well yet.

x = array of buffers to operate on
skip = distance between buffers
depth = log2(number of buffers)
nonzero = number of *output* buffers assumed to be nonzero
length = number of untransformed coefficients requested
extra = indicates whether an extra *forward* coefficient should be computed
twist = twisting power of sqrt2
n = coefficient length
scratch = a scratch buffer
*/
void _ZmodFpoly_IFFT(ZmodF_t* x, unsigned long depth, unsigned long skip,
                     unsigned long nonzero, unsigned long length, int extra,
                     unsigned long twist, unsigned long n,
                     ZmodF_t* scratch)
{
   FLINT_ASSERT((4*n*FLINT_BITS_PER_LIMB) % (1 << depth) == 0);
   FLINT_ASSERT(skip >= 1);
   FLINT_ASSERT(n >= 1);
   FLINT_ASSERT(nonzero >= 1 && nonzero <= (1UL << depth));
   FLINT_ASSERT(length <= nonzero);
   FLINT_ASSERT((length == 0 && extra) ||
                (length == (1UL << depth) && !extra) ||
                (length > 0 && length < (1UL << depth)));
   FLINT_ASSERT(depth >= 1);

   // If the data fits in L1 (2^depth coefficients of length n+1, plus a
   // scratch buffer), then use the iterative transform. Otherwise factor the
   // FFT into two chunks.
   if (depth == 1 ||
       ((1 << depth) + 1) * (n+1) <= ZMODFPOLY_FFT_FACTOR_THRESHOLD)
   {
      _ZmodFpoly_IFFT_recursive(x, depth, skip, nonzero, length, extra,
                                twist, n, scratch);
   }
   else
   {
      unsigned long rows_depth = depth >> 1;
      unsigned long cols_depth = depth - rows_depth;
      _ZmodFpoly_IFFT_factor(x, rows_depth, cols_depth, skip, nonzero, length,
                             extra, twist, n, scratch);
   }
}


/****************************************************************************

   Forward "dual" fourier transforms (internal code)

(twists are applied *before* the transform instead of afterwards, so these
are used for e.g. negacyclic transforms)

****************************************************************************/

/*
Let M = 2^depth
2^root = Mth root of unity
input is assumed to be mod x^M - a^M, where a = 2^twist

assumes twist nonzero
*/
void _ZmodFpoly_FFT_dual_recursive(
            ZmodF_t* x, unsigned long depth,
            unsigned long twist, unsigned long root,
            unsigned long n, ZmodF_t* scratch)
{
   FLINT_ASSERT(twist);
   FLINT_ASSERT(twist < root);

   // =========================================================================
   // special cases for length <= 4

   if (depth == 2)
   {
      // length == 4
      
      // ----------------------------------------------------------------------
      // Do the outer layer of two butterflies first. This is basically an
      // unrolled version of the length >= 8 case below.

      unsigned long bits = (2*twist) & (FLINT_BITS_PER_LIMB-1);
      unsigned long limbs = n - (twist >> (FLINT_LG_BITS_PER_LIMB-1));
      
      if (bits)
      {
         // each butterfly needs a bitshift
         bits = FLINT_BITS_PER_LIMB - bits;
         if (--limbs)
         {
            ZmodF_short_div_2exp(*scratch, x[2], bits, n);
            ZmodF_div_Bexp_add(x[2], x[0], *scratch, limbs, n);
            ZmodF_div_Bexp_sub(x[0], x[0], *scratch, limbs, n);

            ZmodF_short_div_2exp(*scratch, x[3], bits, n);
            ZmodF_div_Bexp_add(x[3], x[1], *scratch, limbs, n);
            ZmodF_div_Bexp_sub(x[1], x[1], *scratch, limbs, n);
         }
         else
         {
            ZmodF_short_div_2exp(*scratch, x[2], bits, n);
            ZmodF_add(x[2], x[0], *scratch, n);
            ZmodF_sub(x[0], x[0], *scratch, n);

            ZmodF_short_div_2exp(*scratch, x[3], bits, n);
            ZmodF_add(x[3], x[1], *scratch, n);
            ZmodF_sub(x[1], x[1], *scratch, n);
         }
      }
      else
      {
         // no bitshifts needed
         ZmodF_div_Bexp_add(*scratch, x[0], x[2], limbs, n);
         ZmodF_swap(scratch, x+2);
         ZmodF_div_Bexp_sub(x[0], x[0], *scratch, limbs, n);

         ZmodF_div_Bexp_add(*scratch, x[1], x[3], limbs, n);
         ZmodF_swap(scratch, x+3);
         ZmodF_div_Bexp_sub(x[1], x[1], *scratch, limbs, n);
      }

      // ----------------------------------------------------------------------
      // Now do the bottom layer, two "blocks" of one butterfly each.

      twist = n*FLINT_BITS_PER_LIMB - twist;
      ZmodF_inverse_butterfly_2exp(x, x+1, scratch, twist, n);
      ZmodF_swap(x, x+1);
      ZmodF_inverse_butterfly_2exp(x+2, x+3, scratch, twist - root, n);
      ZmodF_swap(x+2, x+3);

      return;
   }

   if (depth <= 1)
   {
      // length == 1 or 2
      if (depth == 1)
      {
         ZmodF_inverse_butterfly_2exp(x, x+1, scratch,
                                      n*FLINT_BITS_PER_LIMB - twist, n);
         ZmodF_swap(x, x+1);
      }
      return;
   }
   
   // =========================================================================
   // general case for length >= 8

   // butterflies (a, b) -> (a + w*b, a - w*b), where w = 2^(amount).
   unsigned long half = 1 << (depth - 1);
   ZmodF_t* y = x + half;
   unsigned long amount = twist << (depth - 1);
   unsigned long bits = amount & (FLINT_BITS_PER_LIMB-1);
   unsigned long limbs = n - (amount >> FLINT_LG_BITS_PER_LIMB);
   
   if (bits)
   {
      // each butterfly needs a bitshift
      bits = FLINT_BITS_PER_LIMB - bits;
      if (--limbs)
      {
         for (unsigned long i = 0; i < half; i++)
         {
            ZmodF_short_div_2exp(*scratch, y[i], bits, n);
            ZmodF_div_Bexp_add(y[i], x[i], *scratch, limbs, n);
            ZmodF_div_Bexp_sub(x[i], x[i], *scratch, limbs, n);
         }
      }
      else
      {
         for (unsigned long i = 0; i < half; i++)
         {
            ZmodF_short_div_2exp(*scratch, y[i], bits, n);
            ZmodF_add(y[i], x[i], *scratch, n);
            ZmodF_sub(x[i], x[i], *scratch, n);
         }
      }
   }
   else
   {
      // all butterflies are limbshifts only
      for (unsigned long i = 0; i < half; i++)
      {
         ZmodF_div_Bexp_add(*scratch, x[i], y[i], limbs, n);
         ZmodF_swap(scratch, y+i);
         ZmodF_div_Bexp_sub(x[i], x[i], *scratch, limbs, n);
      }
   }
   
   // =========================================================================
   // recurse into two halves

   _ZmodFpoly_FFT_dual_recursive(x, depth-1, twist, root << 1, n, scratch);
   _ZmodFpoly_FFT_dual_recursive(x + half, depth-1, twist + root, root << 1,
                                 n, scratch);
}



void _ZmodFpoly_IFFT_dual_recursive(
            ZmodF_t* x, unsigned long depth,
            unsigned long twist, unsigned long root,
            unsigned long n, ZmodF_t* scratch)
{
   // =========================================================================
   // special cases for length <= 4
   
   if (depth == 2)
   {
      // ----------------------------------------------------------------------
      // Do the inner layer of two "blocks" of one butterfly each.

      unsigned long temp = n*FLINT_BITS_PER_LIMB - twist;
      ZmodF_forward_butterfly_2exp(x+3, x+2, scratch, temp - root, n);
      ZmodF_swap(x+2, x+3);
      ZmodF_forward_butterfly_2exp(x+1, x, scratch, temp, n);
      ZmodF_swap(x, x+1);

      // ----------------------------------------------------------------------
      // Now do the outer layer of two butterflies. This is basically an
      // unrolled version of the length >= 8 case below.

      unsigned long amount = 2*twist;
      unsigned long bits = amount & (FLINT_BITS_PER_LIMB-1);
      unsigned long limbs = n - (amount >> FLINT_LG_BITS_PER_LIMB);

      if (bits)
      {
         // each butterfly needs a bitshift
         if (limbs != n)
         {
            ZmodF_sub_mul_Bexp(*scratch, x[2], x[0], limbs, n);
            ZmodF_add(x[0], x[0], x[2], n);
            ZmodF_short_div_2exp(x[2], *scratch, bits, n);

            ZmodF_sub_mul_Bexp(*scratch, x[3], x[1], limbs, n);
            ZmodF_add(x[1], x[1], x[3], n);
            ZmodF_short_div_2exp(x[3], *scratch, bits, n);
         }
         else
         {
            ZmodF_sub(*scratch, x[0], x[2], n);
            ZmodF_add(x[0], x[0], x[2], n);
            ZmodF_short_div_2exp(x[2], *scratch, bits, n);

            ZmodF_sub(*scratch, x[1], x[3], n);
            ZmodF_add(x[1], x[1], x[3], n);
            ZmodF_short_div_2exp(x[3], *scratch, bits, n);
         }
      }
      else
      {
         // no bitshifts required
         ZmodF_sub_mul_Bexp(*scratch, x[2], x[0], limbs, n);
         ZmodF_add(x[0], x[0], x[2], n);
         ZmodF_swap(x+2, scratch);

         ZmodF_sub_mul_Bexp(*scratch, x[3], x[1], limbs, n);
         ZmodF_add(x[1], x[1], x[3], n);
         ZmodF_swap(x+3, scratch);
      }

      return;
   }

   if (depth <= 1)
   {
      if (depth == 1)
      {
         ZmodF_forward_butterfly_2exp(x+1, x, scratch, twist, n);
         ZmodF_swap(x, x+1);
      }
      return;
   }

   unsigned long half = 1 << (depth - 1);
   
   // =========================================================================
   // recurse into two halves

   _ZmodFpoly_IFFT_dual_recursive(x, depth-1, twist, root << 1, n, scratch);
   _ZmodFpoly_IFFT_dual_recursive(x + half, depth-1, twist + root, root << 1,
                                  n, scratch);

   // =========================================================================
   // general case for length >= 8

   // butterflies (a, b) -> (a + b, w*(a - b)), where w = 2^(-amount).
   ZmodF_t* y = x + half;
   unsigned long amount = twist << (depth - 1);
   unsigned long bits = amount & (FLINT_BITS_PER_LIMB-1);
   unsigned long limbs = n - (amount >> FLINT_LG_BITS_PER_LIMB);

   if (bits)
   {
      // each butterfly needs a bitshift
      if (limbs != n)
      {
         for (unsigned long i = 0; i < half; i++)
         {
            ZmodF_sub_mul_Bexp(*scratch, y[i], x[i], limbs, n);
            ZmodF_add(x[i], x[i], y[i], n);
            ZmodF_short_div_2exp(y[i], *scratch, bits, n);
         }
      }
      else
      {
         for (unsigned long i = 0; i < half; i++)
         {
            ZmodF_sub(*scratch, x[i], y[i], n);
            ZmodF_add(x[i], x[i], y[i], n);
            ZmodF_short_div_2exp(y[i], *scratch, bits, n);
         }
      }
   }
   else
   {
      // all butterflies are limbshifts only
      for (unsigned long i = 0; i < half; i++)
      {
         ZmodF_sub_mul_Bexp(*scratch, y[i], x[i], limbs, n);
         ZmodF_add(x[i], x[i], y[i], n);
         ZmodF_swap(y+i, scratch);
      }
   }
}


/****************************************************************************

   Fourier Transform Routines

****************************************************************************/


void ZmodFpoly_FFT(ZmodFpoly_t poly, unsigned long length)
{
   FLINT_ASSERT(length <= (1UL << poly->depth));
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (length != 0)
   {
      if (poly->length == 0)
      {
         // input is zero, so output is zero too
         for (unsigned long i = 0; i < length; i++)
            ZmodF_zero(poly->coeffs[i], poly->n);
      }
      else
      {
         if (poly->depth >= 1)
            _ZmodFpoly_FFT(poly->coeffs, poly->depth, 1, poly->length,
                           length, 0, poly->n, poly->scratch);
      }
   }

   poly->length = length;
}


void ZmodFpoly_IFFT(ZmodFpoly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((4 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   if (poly->length && poly->depth)
      _ZmodFpoly_IFFT(poly->coeffs, poly->depth, 1, poly->length,
                      poly->length, 0, 0, poly->n, poly->scratch);
}


// res may alias x or y
// x and y may alias each other
void ZmodFpoly_convolution(ZmodFpoly_t res, ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long length = x->length + y->length - 1;
   unsigned long size = 1UL << res->depth;
   if (length > size)
      length = size;
   
   ZmodFpoly_FFT(x, length);
   if (x != y)    // take care of aliasing
      ZmodFpoly_FFT(y, length);
      
   ZmodFpoly_pointwise_mul(res, x, y);
   ZmodFpoly_IFFT(res);
   ZmodFpoly_rescale(res);
}


/****************************************************************************

   Negacyclic Fourier Transform Routines
   
****************************************************************************/


/*
ignores length of poly
*/
void ZmodFpoly_negacyclic_FFT(ZmodFpoly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   unsigned long twist = (poly->n * FLINT_BITS_PER_LIMB) >> poly->depth;

   _ZmodFpoly_FFT_dual_recursive(poly->coeffs, poly->depth, twist, 2*twist, poly->n, poly->scratch);
   poly->length = 1 << poly->depth;
}


void ZmodFpoly_negacyclic_IFFT(ZmodFpoly_t poly)
{
   // check the right roots of unity are available
   FLINT_ASSERT((2 * poly->n * FLINT_BITS_PER_LIMB) % (1 << poly->depth) == 0);
   FLINT_ASSERT(poly->scratch_count >= 1);

   unsigned long twist = (poly->n * FLINT_BITS_PER_LIMB) >> poly->depth;
   _ZmodFpoly_IFFT_dual_recursive(poly->coeffs, poly->depth, twist, 2*twist, poly->n, poly->scratch);
   poly->length = 1 << poly->depth;
}


void ZmodFpoly_negacyclic_convolution(ZmodFpoly_t res,
                                      ZmodFpoly_t x, ZmodFpoly_t y)
{
   FLINT_ASSERT(x->depth == y->depth);
   FLINT_ASSERT(x->depth == res->depth);
   FLINT_ASSERT(x->n == y->n);
   FLINT_ASSERT(x->n == res->n);

   unsigned long size = 1UL << res->depth;
   
   ZmodFpoly_negacyclic_FFT(x);
   if (x != y)    // take care of aliasing
      ZmodFpoly_negacyclic_FFT(y);
      
   ZmodFpoly_pointwise_mul(res, x, y);
   ZmodFpoly_negacyclic_IFFT(res);
   ZmodFpoly_rescale(res);
   res->length = size;
}


// end of file ****************************************************************
