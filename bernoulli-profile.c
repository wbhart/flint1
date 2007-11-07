/****************************************************************************

   zmod_poly-profile.c : Profiling code for zmod_poly

   Copyright (C) 2007, David Howden

*****************************************************************************/

#include "profiler-main.h"
#include "zmod_poly.h"
#include "mpz_poly.h"
#include "fmpz_poly.h"
#include "long_extras.h"
#include "flint.h"
#include <string.h>
#include <math.h>


/*
   Computes the bernoulli numbers B_0, B_2, ..., B_{p-3}
   for prime p
   
   Requires that res be allocated for (p-1)/2 unsigned longs
   which will hold the result.
   
   If returns 0, then the factoring of p has failed, otherwise
   will always return 1.
*/

int bernoulli_mod_p_mpz(unsigned long *res, unsigned long p)
{
   FLINT_ASSERT(p > 2);
   FLINT_ASSERT(z_isprime(p) == 1);
   
   unsigned long g, g_inv, g_sqr, g_sqr_inv;
   double p_inv = z_precompute_inverse(p);
   g = z_primitive_root_precomp(p, p_inv);
   
   if(!g)
   {
      return 0;
   }
   
   g_inv = z_invert(g, p);
   g_sqr = z_mulmod_precomp(g, g, p, p_inv);
   g_sqr_inv = z_mulmod_precomp(g_inv, g_inv, p, p_inv);
   
   unsigned long poly_size = (p-1)/2;
   
   int is_odd = poly_size % 2;
   
   unsigned long g_power, g_power_inv;
   g_power = g_inv;
   g_power_inv = 1;
   
   // constant is (g-1)/2 mod p
   unsigned long constant;
   if(g % 2)
   {
      constant = (g-1)/2;
   }
   else
   {
      constant = (g+p-1)/2;
   }
   
   // fudge holds g^{i^2}, fudge_inv holds g^{-i^2}
   unsigned long fudge, fudge_inv;
   fudge = fudge_inv = 1;
   
   // compute the polynomials F(X) and G(X)
   mpz_poly_t F, G;
   
   mpz_poly_init2(F, poly_size);
   mpz_poly_init2(G, poly_size);
   
   unsigned long i, temp, h;
   
   for(i = 0; i < poly_size; i++)
   {  
      // compute h(g^i)/g^i (h(x) is as in latex notes)
      temp = g * g_power;
            
      h = z_mulmod_precomp(p + constant - (temp / p), g_power_inv, p, p_inv);
      
      g_power = z_mod_precomp(temp, p, p_inv);
      g_power_inv = z_mulmod_precomp(g_power_inv, g_inv, p, p_inv);
      
      // store coefficient g^{i^2} h(g^i)/g^i
      mpz_poly_set_coeff_ui(G, i, z_mulmod_precomp(h, fudge, p, p_inv));
      mpz_poly_set_coeff_ui(F, i, fudge_inv);
      
      // update fudge and fudge_inv
      fudge = z_mulmod_precomp(z_mulmod_precomp(fudge, g_power, p, p_inv), z_mulmod_precomp(g_power, g, p, p_inv), p, p_inv);
      fudge_inv = z_mulmod_precomp(z_mulmod_precomp(fudge_inv, g_power_inv, p, p_inv), z_mulmod_precomp(g_power_inv, g, p, p_inv), p, p_inv);
   }
   
   mpz_poly_set_coeff_ui(F, 0, 0);
   
   // step 2: multiply the polynomials...
   mpz_poly_t product;
   mpz_poly_init(product);
   mpz_poly_mul(product, G, F);
   
   // step 3: assemble the result...   
   unsigned long g_sqr_power, value;
   g_sqr_power = g_sqr;
   fudge = g;

   res[0] = 1;
   
   mpz_t value_coeff;
   mpz_init(value_coeff);
   
   unsigned long value_coeff_ui;

   for(i = 1; i < poly_size; i++)
   {
      mpz_poly_get_coeff(value_coeff, product, i + poly_size);
      value = mpz_fdiv_ui(value_coeff, p);
      
      value = z_mod_precomp(mpz_poly_get_coeff_ui(product, i + poly_size), p, p_inv);
      
      mpz_poly_get_coeff(value_coeff, product, i);
      if(is_odd)
      {
         value = z_mod_precomp(mpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + p - value, p, p_inv);
      }
      else
      {
         value = z_mod_precomp(mpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + value, p, p_inv);
      }
      
      value = z_mulmod_precomp(z_mulmod_precomp(z_mulmod_precomp(4, i, p, p_inv), fudge, p, p_inv), value, p, p_inv);
      value = z_mulmod_precomp(value, z_invert(p+1-g_sqr_power, p), p, p_inv);

      res[i] = value;
      
      g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
      fudge = z_mulmod_precomp(fudge, g_sqr_power, p, p_inv);
      g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
   }
   
   mpz_poly_clear(product);
   mpz_poly_clear(F);
   mpz_poly_clear(G);
   mpz_clear(value_coeff);
   
   return 1;
}

int bernoulli_mod_p_fmpz(unsigned long *res, unsigned long p)
{
   FLINT_ASSERT(p > 2);
   FLINT_ASSERT(z_isprime(p) == 1);
   
   unsigned long g, g_inv, g_sqr, g_sqr_inv;
   double p_inv = z_precompute_inverse(p);
   g = z_primitive_root_precomp(p, p_inv);
   
   if(!g)
   {
      return 0;
   }
   
   g_inv = z_invert(g, p);
   g_sqr = z_mulmod_precomp(g, g, p, p_inv);
   g_sqr_inv = z_mulmod_precomp(g_inv, g_inv, p, p_inv);
   
   unsigned long poly_size = (p-1)/2;
   
   int is_odd = poly_size % 2;
   
   unsigned long g_power, g_power_inv;
   g_power = g_inv;
   g_power_inv = 1;
   
   // constant is (g-1)/2 mod p
   unsigned long constant;
   if(g % 2)
   {
      constant = (g-1)/2;
   }
   else
   {
      constant = (g+p-1)/2;
   }
   
   // fudge holds g^{i^2}, fudge_inv holds g^{-i^2}
   unsigned long fudge, fudge_inv;
   fudge = fudge_inv = 1;
   
   // compute the polynomials F(X) and G(X)
   fmpz_poly_t F, G;
   
   fmpz_poly_init2(F, poly_size, 2);
   fmpz_poly_init2(G, poly_size, 2);
   
   unsigned long i, temp, h;
   
   for(i = 0; i < poly_size; i++)
   {  
      // compute h(g^i)/g^i (h(x) is as in latex notes)
      temp = g * g_power;
            
      h = z_mulmod_precomp(p + constant - (temp / p), g_power_inv, p, p_inv);
      
      g_power = z_mod_precomp(temp, p, p_inv);
      g_power_inv = z_mulmod_precomp(g_power_inv, g_inv, p, p_inv);
      
      // store coefficient g^{i^2} h(g^i)/g^i
      fmpz_poly_set_coeff_ui(G, i, z_mulmod_precomp(h, fudge, p, p_inv));
      fmpz_poly_set_coeff_ui(F, i, fudge_inv);
      
      // update fudge and fudge_inv
      fudge = z_mulmod_precomp(z_mulmod_precomp(fudge, g_power, p, p_inv), z_mulmod_precomp(g_power, g, p, p_inv), p, p_inv);
      fudge_inv = z_mulmod_precomp(z_mulmod_precomp(fudge_inv, g_power_inv, p, p_inv), z_mulmod_precomp(g_power_inv, g, p, p_inv), p, p_inv);
   }
   
   fmpz_poly_set_coeff_ui(F, 0, 0);
   
   // step 2: multiply the polynomials...
   
   fmpz_poly_t product;
   fmpz_poly_init(product);
   fmpz_poly_mul(product, G, F);

   // step 3: assemble the result...   
   unsigned long g_sqr_power, value;
   g_sqr_power = g_sqr;
   fudge = g;

   res[0] = 1;
   
   mpz_t value_coeff;
   mpz_init(value_coeff);
   
   unsigned long value_coeff_ui;
   unsigned long value2;

   // know that there are either 1 limbs per coeff or 2 limbs per coeff (since we have a limit on p)
   //if(_fmpz_poly_limbs(product) == 1)
   //{
      for(i = 1; i < poly_size; i++)
      {
         fmpz_poly_get_coeff_mpz(value_coeff, product, i + poly_size);

         //value = z_mod_precomp(_fmpz_poly_get_coeff_ui(product, i+poly_size), p, p_inv);
         //value_coeff_ui = z_mod_precomp(_fmpz_poly_get_coeff_ui(product, i), p, p_inv);

         value = mpz_fdiv_ui(value_coeff, p);
         //if(value != value2) printf("ERROR!!!! %d != %d\n", value, value2);


         fmpz_poly_get_coeff_mpz(value_coeff, product, i);
         if(is_odd)
         {
            value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + p - value, p, p_inv);
            //value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + value_coeff_ui + p - value, p, p_inv);
         }
         else
         {
            value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + value, p, p_inv);
            //value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + value_coeff_ui + value, p, p_inv);
         }

         value = z_mulmod_precomp(z_mulmod_precomp(z_mulmod_precomp(4, i, p, p_inv), fudge, p, p_inv), value, p, p_inv);

         value = z_mulmod_precomp(value, z_invert(p+1-g_sqr_power, p), p, p_inv);

         res[i] = value;

         g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
         fudge = z_mulmod_precomp(fudge, g_sqr_power, p, p_inv);
         g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
      }
   // }
   //    else
   //    {
   //       for(i = 1; i < poly_size; i++)
   //       {
   //          //fmpz_poly_get_coeff_mpz(value_coeff, product, i + poly_size);
   // 
   //          value2 = z_ll_mod_precomp(product->coeffs[i+poly_size], product->coeffs[i+poly_size + 1], p, p_inv);
   //          value_coeff_ui = z_ll_mod_precomp(product->coeffs[i], product->coeffs[i+1], p, p_inv);
   // 
   //          //value = mpz_fdiv_ui(value_coeff, p);
   //          //if(value != value2) printf("ERROR!!!! %d != %d\n", value, value2);
   // 
   // 
   //          //fmpz_poly_get_coeff_mpz(value_coeff, product, i);
   //          if(is_odd)
   //          {
   //             //value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + p - value, p, p_inv);
   //             value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + value_coeff_ui + p - value, p, p_inv);
   //          }
   //          else
   //          {
   //             //value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + mpz_fdiv_ui(value_coeff, p) + value, p, p_inv);
   //             value = z_mod_precomp(_fmpz_poly_get_coeff_ui(G, i) + value_coeff_ui + value, p, p_inv);
   //          }
   // 
   //          value = z_mulmod_precomp(z_mulmod_precomp(z_mulmod_precomp(4, i, p, p_inv), fudge, p, p_inv), value, p, p_inv);
   // 
   //          value = z_mulmod_precomp(value, z_invert(p+1-g_sqr_power, p), p, p_inv);
   // 
   //          res[i] = value;
   // 
   //          g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
   //          fudge = z_mulmod_precomp(fudge, g_sqr_power, p, p_inv);
   //          g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
   //       }
   // 
   //    }
      fmpz_poly_clear(product);
      fmpz_poly_clear(F);
      fmpz_poly_clear(G);
   return 1;
}

int bernoulli_mod_p_zmod(unsigned long *res, unsigned long p)
{
   FLINT_ASSERT(p > 2);
   FLINT_ASSERT(z_isprime(p) == 1);
   
   unsigned long g, g_inv, g_sqr, g_sqr_inv;
   double p_inv = z_precompute_inverse(p);
   g = z_primitive_root_precomp(p, p_inv);
   
   if(!g)
   {
      return 0;
   }
   
   g_inv = z_invert(g, p);
   g_sqr = z_mulmod_precomp(g, g, p, p_inv);
   g_sqr_inv = z_mulmod_precomp(g_inv, g_inv, p, p_inv);
   
   unsigned long poly_size = (p-1)/2;
   
   int is_odd = poly_size % 2;
   
   unsigned long g_power, g_power_inv;
   g_power = g_inv;
   g_power_inv = 1;
   
   // constant is (g-1)/2 mod p
   unsigned long constant;
   if(g % 2)
   {
      constant = (g-1)/2;
   }
   else
   {
      constant = (g+p-1)/2;
   }
   
   // fudge holds g^{i^2}, fudge_inv holds g^{-i^2}
   unsigned long fudge, fudge_inv;
   fudge = fudge_inv = 1;
   
   // compute the polynomials F(X) and G(X)   
   zmod_poly_t F, G;
   zmod_poly_init2(F, p, poly_size);
   zmod_poly_init2(G, p, poly_size);
   
   unsigned long i, temp, h;
   
   for(i = 0; i < poly_size; i++)
   {  
      // compute h(g^i)/g^i (h(x) is as in latex notes)
      temp = g * g_power;
            
      h = z_mulmod_precomp(p + constant - (temp / p), g_power_inv, p, p_inv);
      
      g_power = z_mod_precomp(temp, p, p_inv);
      g_power_inv = z_mulmod_precomp(g_power_inv, g_inv, p, p_inv);
      
      // store coefficient g^{i^2} h(g^i)/g^i
      zmod_poly_set_coeff(G, i, z_mulmod_precomp(h, fudge, p, p_inv));
      zmod_poly_set_coeff(F, i, fudge_inv);
      
      // update fudge and fudge_inv
      fudge = z_mulmod_precomp(z_mulmod_precomp(fudge, g_power, p, p_inv), z_mulmod_precomp(g_power, g, p, p_inv), p, p_inv);
      fudge_inv = z_mulmod_precomp(z_mulmod_precomp(fudge_inv, g_power_inv, p, p_inv), z_mulmod_precomp(g_power_inv, g, p, p_inv), p, p_inv);
   }
   
   zmod_poly_set_coeff(F, 0, 0);
   
   // step 2: multiply the polynomials...
   
   zmod_poly_t product;
   zmod_poly_init(product, p);
   zmod_poly_mul_KS(product, G, F, 0);

   // step 3: assemble the result...   
   unsigned long g_sqr_power, value;
   g_sqr_power = g_sqr;
   fudge = g;

   res[0] = 1;
   
   unsigned long value_coeff_ui;

   for(i = 1; i < poly_size; i++)
   {
      value = zmod_poly_get_coeff(product, i + poly_size);
      
      if(is_odd)
      {
         value = z_mod_precomp(zmod_poly_get_coeff(G, i) + zmod_poly_get_coeff(product, i) + p - value, p, p_inv);
      }
      else
      {
         value = z_mod_precomp(zmod_poly_get_coeff(G, i) + zmod_poly_get_coeff(product, i) + value, p, p_inv);
      }
      
      value = z_mulmod_precomp(z_mulmod_precomp(z_mulmod_precomp(4, i, p, p_inv), fudge, p, p_inv), value, p, p_inv);
      value = z_mulmod_precomp(value, z_invert(p+1-g_sqr_power, p), p, p_inv);

      res[i] = value;
      
      g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
      fudge = z_mulmod_precomp(fudge, g_sqr_power, p, p_inv);
      g_sqr_power = z_mulmod_precomp(g_sqr_power, g, p, p_inv);
   }
   
   zmod_poly_clear(product);
   zmod_poly_clear(F);
   zmod_poly_clear(G);
   
   return 1;
}

/*
   Verifies that the ouput of bernoulli_mod_p above is correct.
   
   Takes the result from bernoulli_mod_p (res - an array of (p-1)/2
   unsigned longs), and the prime p.
   
   Returns 0 if res is incorrect, 1 if res is correct.
*/

int verify_bernoulli_mod_p(unsigned long *res, unsigned long p)
{
   unsigned long N, i, product, sum, value, element;
   double p_inv;
   N = (p-1)/2;
   product = 1;
   sum = 0;
   
   p_inv = z_precompute_inverse(p);
   
   for(i = 0; i < N; i++)
   {
      element = res[i];
      // if((signed long)element < 0)
      // {
      //    printf("NEGATIVE NUMBER!!!!!\n");
      // }
      // if(element > p)
      // {
      //    printf("OVERFLOW!!!!!\n");
      // }
      value = z_mulmod_precomp(z_mulmod_precomp(product, 2*i+1, p, p_inv), element, p, p_inv);
      sum = z_mod_precomp(sum + value, p, p_inv);
      product = z_mulmod_precomp(product, 4, p, p_inv);
   }
   
   if(z_mod_precomp(sum + 2,  p, p_inv))
   {   
      return 0;
   }
   
   return 1;
}

/*
This is a helper function used by the other sampler functions below.
*/
// void sample_ZmodF_mul_helper(ZmodF_mul_info_t info, unsigned long n,
//                              unsigned long count)
// {
//    mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
//    mp_limb_t* x2 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
//    mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
// 
//    profiler_random_limbs(x1, n);
//    x1[n] = 0;
//    profiler_random_limbs(x2, n);
//    x2[n] = 0;
//    
//    prof_start();
// 
//    for (unsigned long i = 0; i < count; i++)
//       ZmodF_mul_info_mul(info, x3, x1, x2);
// 
//    prof_stop();
// 
//    free(x3);
//    free(x2);
//    free(x1);
// }



/*
This is a helper function used by the other sampler functions below.
*/
// void sample_ZmodF_sqr_helper(ZmodF_mul_info_t info, unsigned long n,
//                              unsigned long count)
// {
//    mp_limb_t* x1 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
//    mp_limb_t* x3 = (mp_limb_t*) malloc((n+1) * sizeof(mp_limb_t));
// 
//    profiler_random_limbs(x1, n);
//    x1[n] = 0;
//    
//    prof_start();
// 
//    for (unsigned long i = 0; i < count; i++)
//       ZmodF_mul_info_mul(info, x3, x1, x1);
// 
//    prof_stop();
// 
//    free(x3);
//    free(x1);
// }


// ============================================================================

void sample_bernoulli_mpz(unsigned long n, void* arg, unsigned long count)
{
   unsigned long *res = (unsigned long*) malloc(sizeof(unsigned long)*((n-1)/2));
   prof_start();
   while (count)
   {
      if(!bernoulli_mod_p_mpz(res, n))
      {
         printf("Could not factor p = %d\n", n);
      }
      count--;
   }
   prof_stop();
   free(res);
}


char* profDriverString_bernoulli_mpz(char* params)
{
   return
   "Bernoulli mpz implementation.\n"
   "Parameters: n (number of primes to test above 2).\n";
}


char* profDriverDefaultParams_bernoulli_mpz()
{
   return "2 1000";
}


void profDriver_bernoulli_mpz(char* params)
{
   unsigned long p;
   unsigned long n;
   sscanf(params, "%ld %ld", &p, &n);

   prof1d_set_sampler(sample_bernoulli_mpz);
   
   for (unsigned long i = 0; i < n; i++)
   {
      p = z_nextprime(p);
      prof1d_sample(p, NULL);
   }
}


// ============================================================================

void sample_bernoulli_fmpz(unsigned long n, void* arg, unsigned long count)
{
   unsigned long *res = (unsigned long*) malloc(sizeof(unsigned long)*((n-1)/2));
   prof_start();
   while (count)
   {
      if(!bernoulli_mod_p_fmpz(res, n))
      {
         printf("Could not factor p = %d\n", n);
      }
      count--;
   }
   prof_stop();
   free(res);
}


char* profDriverString_bernoulli_fmpz(char* params)
{
   return
   "Bernoulli fmpz implementation.\n"
   "Parameters: n (number of primes to test above 2).\n";
}


char* profDriverDefaultParams_bernoulli_fmpz()
{
   return "2 1000";
}


void profDriver_bernoulli_fmpz(char* params)
{
   unsigned long p;
   unsigned long n;
   sscanf(params, "%ld %ld", &p, &n);

   prof1d_set_sampler(sample_bernoulli_fmpz);
   
   for (unsigned long i = 0; i < n; i++)
   {
      p = z_nextprime(p);
      prof1d_sample(p, NULL);
   }
}


// ============================================================================

void sample_bernoulli_zmod(unsigned long n, void* arg, unsigned long count)
{
   unsigned long *res = (unsigned long*) malloc(sizeof(unsigned long)*((n-1)/2));
   prof_start();
   while (count)
   {
      if(!bernoulli_mod_p_zmod(res, n))
      {
         printf("Could not factor p = %d\n", n);
      }
      count--;
   }
   prof_stop();
   free(res);
}


char* profDriverString_bernoulli_zmod(char* params)
{
   return
   "Bernoulli zmod implementation.\n"
   "Parameters: n (number of primes to test above 2).\n";
}


char* profDriverDefaultParams_bernoulli_zmod()
{
   return "2 1000";
}


void profDriver_bernoulli_zmod(char* params)
{
   unsigned long p;
   unsigned long n;
   sscanf(params, "%ld %ld", &p, &n);

   prof1d_set_sampler(sample_bernoulli_zmod);
   
   for (unsigned long i = 0; i < n; i++)
   {
      p = z_nextprime(p);
      prof1d_sample(p, NULL);
   }  
}


// ============================================================================


// void sample_ZmodF_mul_threeway(unsigned long n, void* arg, unsigned long count)
// {
//    ZmodF_mul_info_t info;
//    ZmodF_mul_info_init_threeway(info, n, 0);
//    sample_ZmodF_mul_helper(info, n, count);
//    ZmodF_mul_info_clear(info);
// }
// 
// 
// char* profDriverString_ZmodF_mul_threeway(char* params)
// {
//    return
//    "ZmodF_mul using threeway algorithm.\n"
//    "Parameters: n_min, n_max, n_skip.\n"
//    "Note: n not divisible by 3 are skipped.\n";
// }
// 
// 
// char* profDriverDefaultParams_ZmodF_mul_threeway()
// {
//    return "1 1000 1";
// }
// 
// 
// void profDriver_ZmodF_mul_threeway(char* params)
// {
//    unsigned long n_min, n_max, n_skip;
//    sscanf(params, "%ld %ld %ld", &n_min, &n_max, &n_skip);
// 
//    prof1d_set_sampler(sample_ZmodF_mul_threeway);
//    
//    // round up n_min so we start on a permissible value
//    while (n_min % 3)
//       n_min++;
//    
//    for (unsigned long n = n_min; n <= n_max; n += n_skip)
//    {
//       if (n % 3 == 0)
//          prof1d_sample(n, NULL);
//    }
// }


// ============================================================================


// void sample_ZmodF_mul_auto(unsigned long n, void* arg, unsigned long count)
// {
//    ZmodF_mul_info_t info;
//    ZmodF_mul_info_init(info, n, 0);
//    sample_ZmodF_mul_helper(info, n, count);
//    ZmodF_mul_info_clear(info);
// }
// 
// 
// char* profDriverString_ZmodF_mul_auto(char* params)
// {
//    return
//    "ZmodF_mul using automatically selected algorithm.\n"
//    "Parameters: n_min, n_max, n_skip.\n";
// }
// 
// 
// char* profDriverDefaultParams_ZmodF_mul_auto()
// {
//    return "1 1000 1";
// }
// 
// 
// void profDriver_ZmodF_mul_auto(char* params)
// {
//    unsigned long n_min, n_max, n_skip;
//    sscanf(params, "%ld %ld %ld", &n_min, &n_max, &n_skip);
// 
//    prof1d_set_sampler(sample_ZmodF_mul_auto);
//    
//    for (unsigned long n = n_min; n <= n_max; n += n_skip)
//       prof1d_sample(n, NULL);
// }


// ============================================================================


// void sample_ZmodF_sqr_plain(unsigned long n, void* arg, unsigned long count)
// {
//    ZmodF_mul_info_t info;
//    ZmodF_mul_info_init_plain(info, n, 1);
//    sample_ZmodF_sqr_helper(info, n, count);
//    ZmodF_mul_info_clear(info);
// }
// 
// 
// char* profDriverString_ZmodF_sqr_plain(char* params)
// {
//    return
//    "ZmodF_sqr using plain algorithm.\n"
//    "Parameters: n_min, n_max, n_skip.\n";
// }
// 
// 
// char* profDriverDefaultParams_ZmodF_sqr_plain()
// {
//    return "1 1000 1";
// }
// 
// 
// void profDriver_ZmodF_sqr_plain(char* params)
// {
//    unsigned long n_min, n_max, n_skip;
//    sscanf(params, "%ld %ld %ld", &n_min, &n_max, &n_skip);
// 
//    prof1d_set_sampler(sample_ZmodF_sqr_plain);
//    
//    for (unsigned long n = n_min; n <= n_max; n += n_skip)
//       prof1d_sample(n, NULL);
// }



// ============================================================================


// void sample_ZmodF_sqr_threeway(unsigned long n, void* arg, unsigned long count)
// {
//    ZmodF_mul_info_t info;
//    ZmodF_mul_info_init_threeway(info, n, 1);
//    sample_ZmodF_sqr_helper(info, n, count);
//    ZmodF_mul_info_clear(info);
// }
// 
// 
// char* profDriverString_ZmodF_sqr_threeway(char* params)
// {
//    return
//    "ZmodF_sqr using threeway algorithm.\n"
//    "Parameters: n_min, n_max, n_skip.\n"
//    "Note: n not divisible by 3 are skipped.\n";
// }
// 
// 
// char* profDriverDefaultParams_ZmodF_sqr_threeway()
// {
//    return "1 1000 1";
// }
// 
// 
// void profDriver_ZmodF_sqr_threeway(char* params)
// {
//    unsigned long n_min, n_max, n_skip;
//    sscanf(params, "%ld %ld %ld", &n_min, &n_max, &n_skip);
// 
//    prof1d_set_sampler(sample_ZmodF_sqr_threeway);
//    
//    // round up n_min so we start on a permissible value
//    while (n_min % 3)
//       n_min++;
//    
//    for (unsigned long n = n_min; n <= n_max; n += n_skip)
//    {
//       if (n % 3 == 0)
//          prof1d_sample(n, NULL);
//    }
// }


// ============================================================================


// void sample_ZmodF_sqr_auto(unsigned long n, void* arg, unsigned long count)
// {
//    ZmodF_mul_info_t info;
//    ZmodF_mul_info_init(info, n, 1);
//    sample_ZmodF_sqr_helper(info, n, count);
//    ZmodF_mul_info_clear(info);
// }
// 
// 
// char* profDriverString_ZmodF_sqr_auto(char* params)
// {
//    return
//    "ZmodF_sqr using automatically selected algorithm.\n"
//    "Parameters: n_min, n_max, n_skip.\n";
// }
// 
// 
// char* profDriverDefaultParams_ZmodF_sqr_auto()
// {
//    return "1 1000 1";
// }
// 
// 
// void profDriver_ZmodF_sqr_auto(char* params)
// {
//    unsigned long n_min, n_max, n_skip;
//    sscanf(params, "%ld %ld %ld", &n_min, &n_max, &n_skip);
// 
//    prof1d_set_sampler(sample_ZmodF_sqr_auto);
//    
//    for (unsigned long n = n_min; n <= n_max; n += n_skip)
//       prof1d_sample(n, NULL);
// }


// end of file ****************************************************************
