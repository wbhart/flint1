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
/****************************************************************************

   zmod_poly-profile.c : Profiling code for zmod_poly

   Copyright (C) 2007, David Howden

*****************************************************************************/

#include "profiler-main.h"
#include "zmod_poly.h"
#include "long_extras.h"
#include "flint.h"
#include <string.h>
#include <math.h>

// ============================================================================

unsigned long randint(unsigned long limit) 
{
#if FLINT_BITS == 32
    static uint64_t randval = 4035456057U;
    randval = ((uint64_t)randval*(uint64_t)1025416097U+(uint64_t)286824430U)%(uint64_t)4294967311U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)randval%limit;
#else
    static unsigned long randval = 4035456057U;
    static unsigned long randval2 = 6748392731U;
    randval = ((unsigned long)randval*(unsigned long)1025416097U+(unsigned long)286824428U)%(unsigned long)4294967311U;
    randval2 = ((unsigned long)randval2*(unsigned long)1647637699U+(unsigned long)286824428U)%(unsigned long)4294967357U;
    
    if (limit == 0L) return (unsigned long) randval;
    
    return (unsigned long)(randval+(randval2<<32))%limit;
#endif
}

/*
   Generate a random integer with up to the given number of bits [0, FLINT_BITS]
*/

unsigned long randbits(unsigned long bits)
{
   return randint(l_shift(1L, bits));
}

/* Generate a random zmod polynomial with the modulus n of the given length with 
   normalised coefficients */

void randpoly(zmod_poly_t poly, long length, unsigned long n)
{
   if (length == 0) 
   {
      zmod_poly_fit_length(poly, 1);
      poly->length = 0;
      return;
   }
              
   zmod_poly_fit_length(poly, length);
   
   for (unsigned long i = 0; i < length; i++)
      poly->coeffs[i] = randint(n);
   poly->length = length;
      
   __zmod_poly_normalise(poly);
} 

void sample_zmod_poly_mul_KS(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, pol2, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(pol2, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(pol2);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(pol2, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
         randpoly(pol2, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_mul_KS(res1, pol1, pol2, 0);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(pol2);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_mul_KS(char* params)
{
   return
   "zmod_poly_mul_KS over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_mul_KS()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_mul_KS(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_mul_KS);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==================================================================================

void sample_zmod_poly_mul_KS_trunc(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, pol2, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(pol2, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(pol2);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(pol2, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
         randpoly(pol2, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_mul_KS_trunc(res1, pol1, pol2, 0, length);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(pol2);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_mul_KS_trunc(char* params)
{
   return
   "zmod_poly_mul_KS_trunc over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_mul_KS_trunc()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_mul_KS_trunc(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_mul_KS_trunc);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==================================================================================

void sample_zmod_poly_mul_classical(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, pol2, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(pol2, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(pol2);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(pol2, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
         randpoly(pol2, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_mul_classical(res1, pol1, pol2);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(pol2);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_mul_classical(char* params)
{
   return
   "zmod_poly_mul_classical over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_mul_classical()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_mul_classical(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_mul_classical);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==============================================================================

void sample_zmod_poly_mul_classical_trunc(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, pol2, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(pol2, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(pol2);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(pol2, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
         randpoly(pol2, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_mul_classical_trunc(res1, pol1, pol2, length);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(pol2);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_mul_classical_trunc(char* params)
{
   return
   "zmod_poly_mul_classical_trunc over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_mul_classical_trunc()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_mul_classical_trunc(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_mul_classical_trunc);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==============================================================================

void sample_zmod_poly_mul_classical_trunc_left(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, pol2, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(pol2, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(pol2);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(pol2, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
         randpoly(pol2, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_mul_classical_trunc_left(res1, pol1, pol2, length);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(pol2);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_mul_classical_trunc_left(char* params)
{
   return
   "zmod_poly_mul_classical_trunc_left over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_mul_classical_trunc_left()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_mul_classical_trunc_left(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_mul_classical_trunc_left);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==============================================================================

void sample_zmod_poly_sqr_KS(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_mul_KS(res1, pol1, pol1, 0);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_sqr_KS(char* params)
{
   return
   "zmod_poly_mul_KS squaring over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_sqr_KS()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_sqr_KS(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_sqr_KS);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==================================================================================

void sample_zmod_poly_sqr_classical(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, res1;
   unsigned long modulus;
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(res1, modulus);
    
   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = randbits(bits);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(res1);  
         
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_sqr_classical(res1, pol1);
      prof_stop();
      
   }
      
   zmod_poly_clear(pol1);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_sqr_classical(char* params)
{
   return
   "zmod_poly_sqr_classical over various lengths and various bit sizes.\n"
   "Parameters: n_min, n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_sqr_classical()
{
   return "1 1000 1.2";
}


void profDriver_zmod_poly_sqr_classical(char* params)
{
   unsigned long n, n_min, n_max;
   double n_ratio;
   sscanf(params, "%ld %ld %lf", &n_min, &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_sqr_classical);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
   int min_iter = (int) ceil(log((double) n_min) / log(n_ratio));
      
   for (unsigned long i = min_iter; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (unsigned long bits = 2; bits < 64; bits++)
         {
             unsigned long log_length = 0;
             while ((1L<<log_length)<n) log_length++;
             if (2*bits + log_length <= 2*FLINT_BITS)
                prof2d_sample(n, bits, NULL);
         }
      }
   }
}

//==================================================================================

void sample_zmod_poly_factor(unsigned long length, unsigned long bits, void* arg, unsigned long count)
{
   zmod_poly_t pol1, pol2, res1;
	zmod_poly_factor_t factors;
   ulong modulus = z_nextprime(z_randbits(bits), 0);
      
   zmod_poly_init(pol1, modulus);
   zmod_poly_init(pol2, modulus);
   zmod_poly_init(res1, modulus);
    
	zmod_poly_factor_init(factors);

   unsigned long r_count;    // how often to generate new random data
   
   if (count >= 1000) r_count = 100;
   else if (count >= 100) r_count = 10;
   else if (count >= 20) r_count = 5;
   else if (count >= 8) r_count = 2;
   else r_count = 1;
     
   for (unsigned long count2 = 0; count2 < count; count2++)
   {     
                
      if (count2 % r_count == 0)
      {
         do {modulus = z_nextprime(z_randbits(bits), 0);} while (modulus < 2);
         
         zmod_poly_clear(pol1);
         zmod_poly_clear(pol2);
         zmod_poly_clear(res1);  
         zmod_poly_init(pol1, modulus);
         zmod_poly_init(pol2, modulus);
         zmod_poly_init(res1, modulus);

         randpoly(pol1, length, modulus);
         randpoly(pol2, length, modulus);

			zmod_poly_mul(res1, pol1, pol2);
      }
        
#if DEBUG
      printf("bits = %ld, length = %ld, modulus = %ld\n", bits, length, modulus);
#endif
   
      prof_start();
      zmod_poly_factor_berlekamp(factors, res1);
      prof_stop();
      zmod_poly_factor_clear(factors);
	
      zmod_poly_factor_init(factors);
	   
   }
      
   zmod_poly_factor_clear(factors);
	zmod_poly_clear(pol1);
   zmod_poly_clear(pol2);
   zmod_poly_clear(res1);  
}


char* profDriverString_zmod_poly_factor(char* params)
{
   return
   "zmod_poly_factor over various lengths and various bit sizes.\n"
   "Parameters: n_max, n_ratio.\n";
}


char* profDriverDefaultParams_zmod_poly_factor()
{
   return "10000 1.2";
}


void profDriver_zmod_poly_factor(char* params)
{
   unsigned long n, n_max, j, old_bits;
   double n_ratio;
   sscanf(params, "%ld %lf", &n_max, &n_ratio);
   unsigned long last_n = 0;
   
   prof2d_set_sampler(sample_zmod_poly_factor);
   
   int max_iter = (int) ceil(log((double) n_max) / log(n_ratio));
      
   for (unsigned long i = 0; i < max_iter; i++)
   {
      n = (unsigned long) floor(pow(n_ratio, i));
      if (n != last_n)
      {
         last_n = n;
         for (ulong bits = 1, j = 0; bits < FLINT_BITS; j++)
         {      
				 prof2d_sample(n, bits, NULL);
				 old_bits = bits;
				 while (old_bits == bits) {bits = (ulong) floor(pow(n_ratio, j)); j++;}
				 j--;
         }
      }
   }
}


// end of file ****************************************************************
