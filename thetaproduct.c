/*============================================================================

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

/*
   Demo FLINT program for computing products of theta functions.
   
   (C) 2008 William Hart
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include "flint.h"
#include "F_mpz.h"
#include "F_mpz_poly.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "theta.h"
#include "profiler.h"
#include "mpn_extras.h"

#define MUL_MOD_TRACE 1

#define FILES1 250
#define FILES2 250

#define LIMIT 62500000000L 
#define BLOCK  100000L
#define BUNDLE 1000L
#define BYTES 2L

#define MOD 16L
#define K 10L
#define SCALE 2

#define COUNT (LIMIT/BLOCK)

void mul_modular_trunc_file(ulong len1, ulong len2,
        zn_mod_t * mod, ulong n, ulong numprimes, const ulong trunc)
{
    ulong page_size = (ulong) sysconf (_SC_PAGESIZE);

	 // Multimodular reduction of poly1, place result into block1
    ulong block1block = (len1+FILES1-1)/FILES1;
	 ulong block1tail = len1 - block1block*(FILES1-1);
    ulong filesize1 = (block1block << n)*sizeof(ulong);
	 filesize1 = ((filesize1+page_size-1)/page_size)*page_size;
    ulong * block1map;
	 char suff[30];

	 int filer;
	 int filer2;
	 char filename[30];

	 if (MUL_MOD_TRACE) printf("len1 = %ld\n", len1);

	 ulong blocklen;

	 ulong * block0map;

	 // Multimodular reduction of poly2, place result into block2
    ulong block2block = (len2+FILES1-1)/FILES1;
	 ulong block2tail = len2 - block2block*(FILES1-1);
    ulong filesize2 = (block2block << n)*sizeof(ulong);
	 filesize2 = ((filesize2+page_size-1)/page_size)*page_size;
    ulong * block2map;
    
    // initialise space for transposed output coefficients
	 ulong len_out = len1 + len2 - 1;

    // Inputs for zn_array_mul (we deal with 16 primes at a time)
    ulong * in1 = flint_heap_alloc(len1*16);
    ulong * in2 = flint_heap_alloc(len2*16);

    // Output for zn_array_mul
    ulong * out = flint_heap_alloc(len_out*16);

	 ulong block3block;
	 ulong block3tail;
    ulong filesize3;
    ulong * block3map;

    ulong i = 0;
	
	for(i = 0; i + 16 <= numprimes; i+= 16) 
	{
 
	 if (MUL_MOD_TRACE) printf("%ld/%ld primes...\n", i, numprimes);
	 // in1 := poly1 % comb->primes[i]
		
	 ulong k;
	 for (k = 0; k < FILES1; k++)
	 {
		 if (k == FILES1 - 1) blocklen = block1tail;
		 else blocklen = block1block;

		 snprintf(suff, 20, "%ld", k);
		 strcpy(filename, "/storage/block1file");
		 filer = open(strcat(filename, suff), O_RDONLY);
		 ulong size = 16*block1block*sizeof(ulong);
		 ulong wanted_offset = i*block1block*sizeof(ulong);
		 ulong offset = (wanted_offset/page_size)*page_size;
       ulong diff = wanted_offset - offset;
		 block1map = mmap(0, size + diff, PROT_READ, MAP_SHARED, filer, offset);
		 diff = diff/sizeof(ulong);
		 block1map += diff;

#pragma omp parallel
	   {		  
#pragma omp for
          for(long j = 0; j < blocklen; j++) 
	      {
             ulong s;
             for (s = 0; s < 16; s++)
				 {
		          ulong l = k*block1block;
					 in1[(j + l) + s*len1] = block1map[j + s*block1block];
				 }
          }
	   }
      
		block1map -= diff;
		diff = diff*sizeof(ulong);
		munmap(block1map, size + diff);
	   close(filer);

	 }
	   
	   // in2 := poly2 % comb->primes[i]

	 ulong k;
	 for (k = 0; k < FILES1; k++)
	 {
		 if (k == FILES1 - 1) blocklen = block2tail;
		 else blocklen = block2block;

		 snprintf(suff, 20, "%ld", k);
		 strcpy(filename, "/storage/block2file");
		 filer = open(strcat(filename, suff), O_RDONLY);
       ulong size = 16*block2block*sizeof(ulong);
		 ulong wanted_offset = i*block2block*sizeof(ulong);
		 ulong offset = (wanted_offset/page_size)*page_size;
       ulong diff = wanted_offset - offset;
		 block2map = mmap(0, size + diff, PROT_READ, MAP_SHARED, filer, offset);
       diff = diff/sizeof(ulong);
		 block2map += diff;

#pragma omp parallel
	   {		  
#pragma omp for
          for(long j = 0; j < blocklen; j++) 
	      {
             ulong s;
             for (s = 0; s < 16; s++)
				 {
		          ulong l = k*block2block;
					 in2[(j + l) + s*len2] = block2map[j + s*block2block];
				 }
          }
	   }

		block2map -= diff;
		diff = diff*sizeof(ulong);
		munmap(block2map, size + diff);
	   close(filer);

	 }
    
       // multiply using zn_poly (requires len1>=len2>=1)
#pragma omp parallel
	   {		  
#pragma omp for
	      long s;
	      for (s = 0; s < 16; s++)
		  {
		     if(len1 >= len2)
                zn_array_mul(out + s*len_out, in1 + s*len1, len1, in2 + s*len2, len2, mod[i + s]);
             else
                zn_array_mul(out + s*len_out, in2 + s*len2, len2, in1 + s*len1, len1, mod[i + s]);
		  }
	   }
		 
		block3block = (trunc+FILES2-1)/FILES2;
	   block3tail = trunc - block3block*(FILES2-1);
      filesize3 = (block3block << n)*sizeof(ulong);
		filesize3 = ((filesize3+page_size-1)/page_size)*page_size;
      
	   ulong k;
	   for (k = 0; k < FILES2; k++)
	   {
		   if (k == FILES2 - 1) blocklen = block3tail;
		   else blocklen = block3block;
			
			snprintf(suff, 20, "%ld", k);
		   strcpy(filename, "/storage/block1file");
		   filer = open(strcat(filename, suff), O_RDWR, (mode_t)0600);
			ulong size = 16*sizeof(ulong)*block3block;
			ulong wanted_offset = i*sizeof(ulong)*block3block;
		   ulong offset = (wanted_offset/page_size)*page_size;
         ulong diff = wanted_offset - offset;
		   block3map = mmap(0, size + diff, PROT_READ | PROT_WRITE, MAP_SHARED, filer, offset);
         diff = diff/sizeof(ulong);
			block3map += diff;

			ulong l = k*block3block;
			   
	   // place result in block_out
#pragma omp parallel
	   {		  
#pragma omp for
         long j;
         for (j = 0; j < blocklen; j++) 
		   {
            ulong s;
            for (s = 0; s < 16; s++)
			   block3map[j + s*block3block] = out[l + j + s*len_out];
         }
	   }   

		   block3map -= diff;
			diff = diff*sizeof(ulong);
		   munmap(block3map, size+diff);
	      close(filer);
		}

   }

   ulong k;
   for (k = 0; k < FILES1; k++)
   {
      snprintf(suff, 20, "%ld", k);
		strcpy(filename, "/storage/block2file");
		remove(strcat(filename, suff));  
   }

	flint_heap_free(in1);
   flint_heap_free(in2);
   flint_heap_free(out);

	// transpose disk files (currently each file has block3block across and 2^ceil_log2(num_primes) down)

	if (MUL_MOD_TRACE) printf("Multimodular multiplications done, trunc = %ld\n", trunc);
}

int main(void)
{
   ulong page_size = (ulong) sysconf (_SC_PAGESIZE);
  
   //--------------------------------------------------------------
   
   fmpz_poly_t theta_1, theta_2, p1, p2, out;
 	
#if FLINT_BITS == 32
   ulong p0 = z_nextprime(1UL << 30, 0);
    
	// primes_per_limb = 32/log2(p0)
   double primes_per_limb = 1.067;
#else 
   ulong p0 = z_nextprime(1L << 62, 0);
    
	// primes_per_limb = 64/log2(p0)
   double primes_per_limb = 1.0323;
#endif

   ulong length = LIMIT/BUNDLE;
   ulong output_bits = BUNDLE*BYTES*8*2 + ceil_log2(length) + 1;

   ulong numprimes = (output_bits * primes_per_limb)/FLINT_BITS + 1;
	numprimes = (((numprimes + 15)>>4)<<4);
	
	// allocate space for primes
	ulong * primes = flint_heap_alloc(numprimes);

   // compute primes
	ulong p = p0;
   ulong i;
   for (i = 0; i < numprimes; i++) 
	{
      primes[i] = p;
      p = z_nextprime(p, 0);
   }

   // precompute comb
   fmpz_comb_t comb;
   fmpz_comb_init(comb, primes, numprimes);

   if (MUL_MOD_TRACE) printf("comb initialised\n");

   long * array1 = (long *) flint_heap_alloc(16*BLOCK);
   
   ulong COEFF_BLOCK = LIMIT/FILES1;
   
   fmpz_poly_init(theta_1);
   fmpz_poly_fit_limbs(theta_1, 1);
   fmpz_poly_fit_length(theta_1, COEFF_BLOCK);

   fmpz_poly_init(p1);

   // Multimodular reduction of poly1, place result into block1
   ulong len1 = LIMIT/BUNDLE;
	ulong block1block = (len1+FILES1-1)/FILES1;
	ulong block1tail = len1 - block1block*(FILES1-1);
   ulong filesize1 = (block1block << comb->n)*sizeof(ulong);
	filesize1 = ((filesize1+page_size-1)/page_size)*page_size;
   ulong * block0map = flint_heap_alloc(block1block<<comb->n);
	ulong * block1map;
	char suff[30];

	int filer;
	int filer2;
	char filename[30];

	ulong blocklen;

	ulong j;
	for (j = 0; j < FILES1; j++)
	{
      theta_1->length = COEFF_BLOCK;
   
#pragma omp parallel
	{
#pragma omp for
	   long start;
	   for (start = 0; start < COEFF_BLOCK*MOD; start += SCALE*BLOCK)
      {   
	      ulong s = omp_get_thread_num();
	      theta_2d_A2(array1 + s*BLOCK, (start + j*COEFF_BLOCK*MOD)/SCALE, BLOCK);
   
         ulong start2 = start/MOD;
	  
	      ulong i;
	      for (i = 0; i < (SCALE*BLOCK)/MOD; i++)
         {
            fmpz_set_si(theta_1->coeffs + (start2 + i)*2, array1[s*BLOCK + (MOD*i+K)/SCALE]);
         }
      }
	}
   
      //_F_mpz_poly_normalise(theta_1);

      fmpz_poly_pack_bytes(p1, theta_1, BUNDLE, BYTES);
      
		if (j == FILES1 - 1) blocklen = block1tail;
		else blocklen = block1block;

      ulong size_p1 = p1->limbs + 1;

      fmpz_t ** comb_temp[16];
      ulong i;
      for (i = 0; i < 16; i++)
      {
         comb_temp[i] = fmpz_comb_temp_init(comb);
      }

#pragma omp parallel
	{		  
#pragma omp for
      for(long i = 0; i < blocklen; i++) 
	   {
         ulong v = omp_get_thread_num();
         fmpz_multi_mod_ui(block0map + (i << comb->n), p1->coeffs + i*size_p1, comb, comb_temp[v]);
		}
	}

      ulong i;
      for (i = 0; i < 16; i++)
      {
         fmpz_comb_temp_clear(comb_temp[i], comb);
      }

	   snprintf(suff, 20, "%ld", j);
		strcpy(filename, "/storage/block1file");
		filer2 = open(strcat(filename, suff), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
		lseek(filer2, filesize1-1, SEEK_SET);
		int res = write(filer2, "", 1);
		block1map = mmap(0, filesize1, PROT_READ | PROT_WRITE, MAP_SHARED, filer2, 0);

#pragma omp parallel
	{		  
#pragma omp for
      long i;
      for (i = 0; i < blocklen; i++)
		{
			ulong k;
			for (k = 0; k < numprimes; k++)
				block1map[i + k*block1block] = block0map[k + (i<<comb->n)];
		}
	}

		munmap(block1map, filesize1);
		close(filer2);

	}

   flint_heap_free(block0map);

	if (MUL_MOD_TRACE) printf("Packing and multimodular reduction and transpose of poly1 done\n");
   
   flint_heap_free(array1);

   fmpz_poly_clear(theta_1);
   
	fmpz_poly_clear(p1);
     
   //-------------------------------------------------------------
   
   array1 = (long *) flint_heap_alloc(16*BLOCK);
	
   fmpz_poly_init(theta_2);
   fmpz_poly_fit_limbs(theta_2, 1);
   fmpz_poly_fit_length(theta_2, COEFF_BLOCK); 
   
   fmpz_poly_init(p2);


   ulong len2 = LIMIT/BUNDLE;
   ulong block2block = (len2+FILES1-1)/FILES1;
	ulong block2tail = len2 - block2block*(FILES1-1);
   ulong filesize2 = (block2block << comb->n)*sizeof(ulong);
	filesize2 = ((filesize2+page_size-1)/page_size)*page_size;
   ulong * block2map;

   block0map = flint_heap_alloc(block2block<<comb->n);
	 
	ulong j;
	for (j = 0; j < FILES1; j++)
	{

      theta_2->length = COEFF_BLOCK;
   
#pragma omp parallel
	{
#pragma omp for
	   long start;
	   for (start = 0; start < COEFF_BLOCK*MOD; start += SCALE*BLOCK)
      {   
	      ulong s = omp_get_thread_num();
	      theta_2d_B(array1 + s*BLOCK, (start + j*COEFF_BLOCK*MOD)/SCALE, BLOCK);
   
         ulong start2 = start/MOD;
	  
	      ulong i;
	      for (i = 0; i < (SCALE*BLOCK)/MOD; i++)
         {
            fmpz_set_si(theta_2->coeffs + (start2 + i)*2, array1[s*BLOCK + (MOD*i)/SCALE]);
         }
      }
	}
   
      //_F_mpz_poly_normalise(theta_2);

   
      fmpz_poly_pack_bytes(p2, theta_2, BUNDLE, BYTES);

		if (j == FILES1 - 1) blocklen = block2tail;
		else blocklen = block2block;

      ulong size_p2 = p2->limbs + 1;

      fmpz_t ** comb_temp[16];
      ulong i;
      for (i = 0; i < 16; i++)
      {
         comb_temp[i] = fmpz_comb_temp_init(comb);
      }

#pragma omp parallel
	{		  
#pragma omp for
      for(long i = 0; i < blocklen; i++) 
	   {
         ulong v = omp_get_thread_num();
         fmpz_multi_mod_ui(block0map + (i << comb->n), p2->coeffs + i*size_p2, comb, comb_temp[v]);
		}
	}
       
	   ulong i;
	   for (i = 0; i < 16; i++)
      {
         fmpz_comb_temp_clear(comb_temp[i], comb);
      }

      snprintf(suff, 20, "%ld", j);
		strcpy(filename, "/storage/block2file");
		filer2 = open(strcat(filename, suff), O_RDWR | O_CREAT | O_TRUNC, (mode_t)0600);
		lseek(filer2, filesize2-1, SEEK_SET);
		int res = write(filer2, "", 1);
		block2map = mmap(0, filesize2, PROT_READ | PROT_WRITE, MAP_SHARED, filer2, 0);

#pragma omp parallel
	{		  
#pragma omp for
      long i;
      for (i = 0; i < blocklen; i++)
		{
			ulong k;
			for (k = 0; k < numprimes; k++)
				block2map[i + k*block2block] = block0map[k + (i<<comb->n)];
		}
	}

		munmap(block2map, filesize2);
		close(filer2);

	}

   flint_heap_free(block0map);
    
	if (MUL_MOD_TRACE) printf("Packing and multimodular reduction and transpose of poly2 done\n");

   flint_heap_free(array1);

   fmpz_poly_clear(theta_2);

   fmpz_poly_clear(p2);
   
   //----------------------------------------------------------------------

   fmpz_poly_init(out);

   ulong trunc = LIMIT/BUNDLE;
   
   mul_modular_trunc_file(LIMIT/BUNDLE, LIMIT/BUNDLE, comb->mod, comb->n, comb->num_primes, trunc);
	
   printf("First product computed\n");

   //F_mpz_poly_unpack_bytes(theta_prod, out, BUNDLE, 2);

	FILE * myfile, * myfile2;
   
   ulong limbs = ((2*BUNDLE)*BYTES*8 - 1)/FLINT_BITS + 2; // max number of limbs of each large coeff
	ulong extras = (FLINT_BYTES_PER_LIMB - 1)/BYTES + 1;

	ulong length_max = BUNDLE*trunc + BUNDLE + extras;

	mp_limb_t * arr1 = flint_heap_alloc(limbs);
	mp_limb_t * arr2 = flint_heap_alloc(limbs);
	mp_limb_t * arr3 = flint_heap_alloc(limbs);

	int neg1 = 0;
	int neg2 = 0;
	int neg3 = 0;

   mp_limb_t * temp;
	int temp_n;
	long limb;
	long carry = 0L;
	int borrow = 0;
   int borrow2 = 0;

	F_mpn_clear(arr1, limbs);
	F_mpn_clear(arr2, limbs);
	F_mpn_clear(arr3, limbs);

	short int * s1 = (short int *) arr1;
	short int * s2 = (short int *) arr2;
	short int * s3 = (short int *) arr3;
	
   ulong limbs2;

   ulong block3block;
	ulong block3tail;
   ulong filesize3;
   ulong * block3map;
   ulong * block4map;

   block3block = (trunc+FILES2-1)/FILES2;
	block3tail = trunc - block3block*(FILES2-1);
   filesize3 = (block3block << comb->n)*sizeof(ulong);
	filesize3 = ((filesize3+page_size-1)/page_size)*page_size;

   short int * theta_orig = malloc(block3block*BUNDLE*sizeof(short int));
   short int * theta_prod = theta_orig;
  
   #define OFFSET (1L<<15)
   #define LEN (1L<<16)

	// zero coeffs, as we will be adding to them
	ulong i;
	for (i = 0; i < block3block*BUNDLE; i++)
		theta_prod[i] = 0;

   ulong off_length;

   unsigned long * arr = flint_heap_alloc(LEN);
   for(long i = 0; i < LEN; i++)
      arr[i] = 0L;

   unsigned long zerocount = 0L;
   long maxneg = 0L;
   long maxpos = 0L;
   long coeff;

   ulong k = 0;
   blocklen = block3block;

   snprintf(suff, 20, "%ld", k);
   strcpy(filename, "/storage/block1file");
	filer = open(strcat(filename, suff), O_RDONLY);
   block3map = mmap(0, filesize3, PROT_READ, MAP_SHARED, filer, 0);

   block4map = flint_heap_alloc(block3block<<comb->n);

#pragma omp parallel
	{		  
#pragma omp for
   long i;
   for (i = 0; i < blocklen; i++)
   {
	   ulong j;
	   for (j = 0; j < numprimes; j++)
		   block4map[j + (i<<comb->n)] = block3map[i + j*block3block];
   }
	}

   munmap(block3map, filesize3);
   close(filer);
   
   strcpy(filename, "/storage/block1file");
	remove(strcat(filename, suff));
   
   k++;
       
   fmpz_poly_fit_length(out, blocklen);
   fmpz_poly_fit_limbs(out, limbs);
   ulong size_out = limbs + 1;

   fmpz_t ** comb_temp[16];
   ulong i;
   for (i = 0; i < 16; i++)
   {
      comb_temp[i] = fmpz_comb_temp_init(comb);
   }

#pragma omp parallel
	{		  
#pragma omp for
      for(long i = 0; i < blocklen; i++) 
	   {
	      ulong v = omp_get_thread_num();
         fmpz_multi_CRT_ui(out->coeffs + i*size_out, block4map + (i << comb->n), comb, comb_temp[v]);
	   }
   }

   ulong psize = z_intsqrt(LIMIT*MOD);
   ulong * prime_arr = malloc(sizeof(ulong)*psize/5+1000);
   ulong prime_num = 0;
   ulong pcurr = 2L;
   while (pcurr <= psize)
   {
      prime_arr[prime_num] = pcurr;
      pcurr = z_nextprime(pcurr, 1);
      prime_num++;
   }

   //fmpz_get_limbs(arr1, out->coeffs); // initialise with first big coefficient
   limbs2 = FLINT_ABS((long)(out->coeffs[0]));
   if (limbs2) F_mpn_copy(arr1, out->coeffs + 1, limbs2);
   if (fmpz_sgn(out->coeffs) < 0) neg1 = 1;
   if (limbs2 > limbs) printf("0 %ld, %ld\n", limbs, limbs2);
	ulong i;
	for (i = 0; i < length_max; i+=BUNDLE)
	{
		ulong j;

	   for (j = 0; (j < extras) && (i + j < length_max); j++)
	   {
         if (neg1) 
			{
				borrow2 = -(s1[j] < 0);
			   carry -= (long) s1[j];
			} else 
			{
				borrow2 = (s1[j] < 0);
			   carry += (long) s1[j];
			}

			if (neg2) 
			{
				borrow2 -= (s2[BUNDLE + j] < 0);
			   carry -= (long) s2[BUNDLE + j];
			} else 
			{
				borrow2 += (s2[BUNDLE + j] < 0);
			   carry += (long) s2[BUNDLE + j];
			}

			if (neg3) 
			{
				borrow2 -= (s3[2*BUNDLE + j] < 0);
			   carry -= (long) s3[2*BUNDLE + j];
			} else 
			{
				borrow2 += (s3[2*BUNDLE + j] < 0);
			   carry += (long) s3[2*BUNDLE + j];
			}

			theta_prod[i + j] = ((short int) carry);
			carry -= (long) theta_prod[i + j];
			carry >>= (BYTES*8);
         
			theta_prod[i + j] += (short int) borrow;

			borrow = borrow2;
		}
      
		for ( ; (j < BUNDLE) && (i + j < length_max); j++)
		{
         if (neg1) 
			{
				borrow2 = -(s1[j] < 0);
			   carry -= (long) s1[j];
			} else 
			{
				borrow2 = (s1[j] < 0);
			   carry += (long) s1[j];
			}

			if (neg2) 
			{
				borrow2 -= (s2[BUNDLE + j] < 0);
			   carry -= (long) s2[BUNDLE + j];
			} else 
			{
				borrow2 += (s2[BUNDLE + j] < 0);
			   carry += (long) s2[BUNDLE + j];
			}

			theta_prod[i + j] = ((short int) carry);
			carry -= (long) theta_prod[i + j];
			carry >>= (BYTES*8);
         
			theta_prod[i + j] += (short int) borrow;
			
			borrow = borrow2;
		}
      
		temp = arr3;
		arr3 = arr2;
		arr2 = arr1;
		arr1 = temp;

		s1 = (short int *) arr1;
	   s2 = (short int *) arr2;
	   s3 = (short int *) arr3;
		
		neg3 = neg2;
		neg2 = neg1;

		F_mpn_clear(arr1, limbs);
	   ulong s = i/BUNDLE + 1;
		ulong t = (s%block3block);
      if (s < trunc)
		{
			if (t == 0UL)
         {
            if (k == FILES2 - 1) blocklen = block3tail;
		      else blocklen = block3block;

            snprintf(suff, 20, "%ld", k);
		      strcpy(filename, "/storage/block1file");
		      filer = open(strcat(filename, suff), O_RDONLY);
            block3map = mmap(0, filesize3, PROT_READ, MAP_SHARED, filer, 0);

#pragma omp parallel
	   {		  
#pragma omp for
            long i;
            for (i = 0; i < blocklen; i++)
		      {
			      ulong j;
			      for (j = 0; j < numprimes; j++)
				      block4map[j + (i<<comb->n)] = block3map[i + j*block3block];
		      }
		}

		      munmap(block3map, filesize3);
            close(filer);
		      
            strcpy(filename, "/storage/block1file");
		      remove(strcat(filename, suff));
            
#pragma omp parallel
	{		  
#pragma omp for
            for(long r = 0; r < blocklen; r++) 
	         {
	            ulong v = omp_get_thread_num();
               fmpz_multi_CRT_ui(out->coeffs + r*size_out, block4map + (r << comb->n), comb, comb_temp[v]);
	         }
    }

         }

         //limbs2 = fmpz_get_limbs(arr1, out->coeffs + t); 
		   limbs2 = FLINT_ABS((long)(out->coeffs[t*size_out]));
         if (limbs2) F_mpn_copy(arr1, out->coeffs + t*size_out + 1, limbs2);
    
			if (fmpz_sgn(out->coeffs + t*size_out) < 0) neg1 = 1;
		   else neg1 = 0;
         if (limbs2 > limbs) printf("%ld %ld, %ld\n", s, limbs, limbs2);
	   }

      if (t == 0UL)
      {
         ulong startoff = (s-block3block)*BUNDLE*MOD + K;
         off_length = block3block*BUNDLE;
         if (k == FILES2 - 1)
         {
            while (!theta_orig[off_length - 1]) off_length--;
         }
         ulong endoff = startoff+off_length*MOD;
         ulong ip;
         for (ip = 0; ip < prime_num; ip++)
         {
            ulong p = prime_arr[ip];
            ulong p2 = p*p;
            ulong off = ((startoff+p2-1)/p2)*p2;
            while (off < endoff)
            {
               if (((off - K) % MOD) == 0L)
               {
                  theta_orig[(off - startoff)/MOD] = -OFFSET;
               }
               off += p2;
            }
         }
         
         snprintf(suff, 20, "%ld", k);
		   strcpy(filename, "zeros2mod16-");
		   myfile = fopen(strcat(filename, suff), "w");
         
		   for(unsigned long j = 0; j < off_length; j ++)
         {
            coeff = theta_orig[j];

            // skip non-squarefree coefficients
            if(coeff == -OFFSET) 
            {
                //printf("arr[%ld] : non-squarefree\n", j);
                continue;
            }

            //printf("arr[%ld] = %ld\n", j, coeff);

            arr[OFFSET+coeff]++;
            if (coeff > maxpos) maxpos = coeff;
            if (coeff < maxneg) maxneg = coeff;
            if (!coeff)
            {
               fprintf(myfile, "%ld ", startoff+MOD*j);
               zerocount++;
            }
         }
         
         fprintf(myfile, "\n");
         fclose(myfile);
         
         k++;
         theta_prod -= (block3block*BUNDLE);
      }
	}

   flint_heap_free(block4map);

   ulong i;
   for (i = 0; i < 16; i++)
   {
      fmpz_comb_temp_clear(comb_temp[i], comb);
   }
     
   flint_heap_free(arr1);
	flint_heap_free(arr2);
	flint_heap_free(arr3);
      
   fmpz_poly_clear(out);

   free(prime_arr);
   
   /*ulong theta_prod_length = LIMIT;
   while (theta_prod[theta_prod_length - 1] == 0) theta_prod_length--;*/

   fmpz_comb_clear(comb);

   // free allocated space
   flint_heap_free(primes);

	printf("Final transpose and unpacking and sieving computed, theta_prod has length %ld\n", ((FILES2 - 1)*block3block*BUNDLE+off_length)*MOD);

   myfile2 = fopen("stats2mod16", "w"); 
   
   for(long i = maxneg; i <= maxpos; i++) {
      fprintf(myfile2, "VALUE = %ld, count = %ld\n", i, arr[i+OFFSET]);
   }
   
   fclose(myfile2);
   flint_heap_free(arr);

   printf("\n\nmaxneg = %ld, maxpos = %ld\nnumzeros = %ld\n\n", maxneg, maxpos, zerocount);

   free(theta_orig);
   
   return 0;
}
