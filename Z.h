#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <gmp.h>

#ifndef FLINT_Z_H
#define FLINT_Z_H

#define Z_t mpz_t

/*-----------------------------------------------------------------------------

    Memory Management Functions
    
-----------------------------------------------------------------------------*/

Z_t* Z_alloc(void);

void Z_release(void);

/*-----------------------------------------------------------------------------

    Movement and Assignment Functions
    
-----------------------------------------------------------------------------*/

static inline void Z_set(Z_t res, Z_t op)
{
   mpz_set(res, op);
}

static inline void Z_set_ui(Z_t res, unsigned long op)
{
   mpz_set_ui(res, op);
}

static inline void Z_set_si(Z_t res, signed long op)
{
   mpz_set_si(res, op);
}

static inline void Z_set_d(Z_t res, double op)
{
   mpz_set_d(res, op);
}

static inline void Z_set_q(Z_t res, mpq_t op)
{
   mpz_set_q(res, op);
}

static inline void Z_set_f(Z_t res, mpf_t op)
{
   mpz_set_f(res, op);
}

static inline int Z_set_str(Z_t res, char *str, int base)
{
   return mpz_set_str(res, str, base);
}

static inline void Z_swap(Z_t op1, Z_t op2)
{
   mpz_swap(op1, op2);
}

/*-----------------------------------------------------------------------------

    Conversion Functions
    
-----------------------------------------------------------------------------*/

static inline unsigned long Z_get_ui(Z_t op)
{
   return mpz_get_ui(op);
}

static inline signed long Z_get_si(Z_t op)
{
   return mpz_get_si(op);
}

static inline double Z_get_d(Z_t op)
{
   return mpz_get_d(op);
}

static inline double Z_get_d_2exp(signed long *exp, Z_t op)
{
   return mpz_get_d_2exp(exp, op);
}

static inline char * Z_get_str(char *str, int base, Z_t op)
{
   return mpz_get_str(str, base, op);
}
     
/*-----------------------------------------------------------------------------

    Arithmetic Functions
    
-----------------------------------------------------------------------------*/

static inline void Z_add(Z_t res, Z_t a, Z_t b)
{
   mpz_add(res, a, b);
}

static inline void Z_add_ui(Z_t res, Z_t a, unsigned long b)
{
   mpz_add_ui(res, a, b);
}

static inline void Z_sub(Z_t res, Z_t a, Z_t b)
{
   mpz_sub(res, a, b);
}

static inline void Z_sub_ui(Z_t res, Z_t a, unsigned long b)
{
   mpz_sub_ui(res, a, b);
}

static inline void Z_ui_sub(Z_t res, unsigned long a, Z_t b)
{
   mpz_ui_sub(res, a, b);
}

static inline void Z_mul(Z_t res, Z_t a, Z_t b)
{
   mpz_mul(res, a, b);
}

static inline void Z_mul_si(Z_t res, Z_t a, long int b)
{
   mpz_mul_si(res, a, b);
}

static inline void Z_mul_ui(Z_t res, Z_t a, unsigned long b)
{
   mpz_mul_ui(res, a, b);
}

static inline void Z_addmul(Z_t res, Z_t a, Z_t b)
{
   mpz_addmul(res, a, b);     
}

static inline void Z_addmul_ui(Z_t res, Z_t a, unsigned long b)
{
   mpz_addmul_ui(res, a, b);  
}

static inline void Z_submul(Z_t res, Z_t a, Z_t b)
{
   mpz_submul(res, a, b); 
}

static inline void Z_submul_ui(Z_t res, Z_t a, unsigned long b)
{
   mpz_submul_ui(res, a, b); 
}

static inline void Z_mul_2exp(Z_t res, Z_t a, unsigned long exp)
{
   mpz_mul_2exp(res, a, exp); 
}

static inline void Z_neg(Z_t res, Z_t n)
{
   mpz_neg(res, n); 
}

static inline void Z_abs(Z_t res, Z_t n)
{
   mpz_abs(res, n);
}

/*-----------------------------------------------------------------------------

    Division Functions
    
-----------------------------------------------------------------------------*/

static inline void Z_cdiv_q(Z_t q, Z_t n, Z_t d)
{
   mpz_cdiv_q(q, n, d);
}

static inline void Z_cdiv_r(Z_t r, Z_t n, Z_t d)
{
   mpz_cdiv_r(r, n, d);
}
 
static inline void Z_cdiv_qr(Z_t q, Z_t r, Z_t n, Z_t d)
{
   mpz_cdiv_qr(q, r, n, d);
}
 
static inline unsigned long Z_cdiv_q_ui(Z_t q, Z_t n, unsigned long d)
{
   return mpz_cdiv_q_ui(q, n, d);
}
 
static inline unsigned long Z_cdiv_r_ui(Z_t r, Z_t n, unsigned long d)
{
   return mpz_cdiv_r_ui(r, n, d);
}
 
static inline unsigned long Z_cdiv_qr_ui(Z_t q, Z_t r, Z_t n, unsigned long d)
{
   return mpz_cdiv_qr_ui(q, r, n, d);
}
 
static inline unsigned long Z_cdiv_ui(Z_t n, unsigned long d)
{
   return mpz_cdiv_ui(n, d);
}
 
static inline void Z_cdiv_q_2exp(Z_t q, Z_t n, unsigned long exp)
{
   mpz_cdiv_q_2exp(q, n, exp);
}
 
static inline void Z_cdiv_r_2exp(Z_t r, Z_t n, unsigned long exp)
{
   mpz_cdiv_r_2exp(r, n, exp);
}

static inline void Z_fdiv_q(Z_t q, Z_t n, Z_t d)
{
   mpz_fdiv_q(q, n, d);
}

static inline void Z_fdiv_r(Z_t r, Z_t n, Z_t d)
{
   mpz_fdiv_r(r, n, d);
}
 
static inline void Z_fdiv_qr(Z_t q, Z_t r, Z_t n, Z_t d)
{
   mpz_fdiv_qr(q, r, n, d);
}
 
static inline unsigned long Z_fdiv_q_ui(Z_t q, Z_t n, unsigned long d)
{
   return mpz_fdiv_q_ui(q, n, d);
}
 
static inline unsigned long Z_fdiv_r_ui(Z_t r, Z_t n, unsigned long d)
{
   return mpz_fdiv_r_ui(r, n, d);
}
 
static inline unsigned long Z_fdiv_qr_ui(Z_t q, Z_t r, Z_t n, unsigned long d)
{
   return mpz_fdiv_qr_ui(q, r, n, d);
}
 
static inline unsigned long Z_fdiv_ui(Z_t n, unsigned long d)
{
   return mpz_fdiv_ui(n, d);
}

static inline void Z_fdiv_q_2exp(Z_t q, Z_t n, unsigned long exp)
{
   mpz_fdiv_q_2exp(q, n, exp);
}

static inline void Z_fdiv_r_2exp(Z_t r, Z_t n, unsigned long exp)
{
   mpz_fdiv_r_2exp(r, n, exp);
}

static inline void Z_tdiv_q(Z_t q, Z_t n, Z_t d)
{
   mpz_tdiv_q(q, n, d);
}

static inline void Z_tdiv_r(Z_t r, Z_t n, Z_t d)
{
   mpz_tdiv_r(r, n, d);
}
 
static inline void Z_tdiv_qr(Z_t q, Z_t r, Z_t n, Z_t d)
{
   mpz_tdiv_qr(q, r, n, d);
}
 
static inline unsigned long Z_tdiv_q_ui(Z_t q, Z_t n, unsigned long d)
{
   return mpz_tdiv_q_ui(q, n, d);
}
 
static inline unsigned long Z_tdiv_r_ui(Z_t r, Z_t n, unsigned long d)
{
   return mpz_tdiv_r_ui(r, n, d);
}
 
static inline unsigned long Z_tdiv_qr_ui(Z_t q, Z_t r, Z_t n, unsigned long d)
{
   return mpz_tdiv_qr_ui(q, r, n, d);
}
 
static inline unsigned long Z_tdiv_ui(Z_t n, unsigned long d)
{
   return mpz_tdiv_ui(n, d);
}
 
static inline void Z_tdiv_q_2exp(Z_t q, Z_t n, unsigned long exp)
{
   mpz_tdiv_q_2exp(q, n, exp);
}
 
static inline void Z_tdiv_r_2exp(Z_t r, Z_t n, unsigned long exp)
{
   mpz_tdiv_r_2exp(r, n, exp);
}

static inline void Z_divexact(Z_t q, Z_t n, Z_t d)
{
   mpz_divexact(q, n, d);
}

static inline void Z_divexact_ui(Z_t q, Z_t n, unsigned long d)
{
   mpz_divexact_ui(q, n, d);
}

/*-----------------------------------------------------------------------------

    Modular Arithmetic
    
-----------------------------------------------------------------------------*/

static inline void Z_mod(Z_t res, Z_t n, Z_t d)
{
   mpz_mod(res, n, d); 
}

static inline unsigned long Z_mod_ui(Z_t res, Z_t n, unsigned long d)
{
   return mpz_mod_ui(res, n, d);
}

void Z_mulmod(Z_t, Z_t, Z_t, Z_t);

/*
  sets res to a*b modulo p
  assumes a and b are not (much) bigger than p and that res is not p       
*/
static inline void mulmod2(mpz_t res, mpz_t a, mpz_t b, mpz_t p)
{
   mpz_mul(res,a,b);
   mpz_fdiv_r(res,res,p);
}

unsigned long Z_mulmod_ui(Z_t, Z_t, Z_t, unsigned long);

static inline int Z_congruent_p(Z_t n, Z_t c, Z_t d)
{
   return mpz_congruent_p(n, c, d);
}

static inline int Z_congruent_ui_p(Z_t n, unsigned long c, unsigned long d)
{
   return mpz_congruent_ui_p(n, c, d);
}

static inline int Z_congruent_2exp_p(Z_t n, Z_t c, unsigned long exp)
{
   return mpz_congruent_2exp_p(n, c, exp);
}

static inline int Z_divisible_p(Z_t n, Z_t d)
{
   return mpz_divisible_p(n, d);
}

static inline int Z_divisible_ui_p(Z_t n, unsigned long d)
{
   return mpz_divisible_ui_p(n, d);
}

static inline int Z_divisible_2exp_p(Z_t n, unsigned long exp)
{
   return mpz_divisible_2exp_p(n, exp);
}

long Z_powm_long(long, long, long);

static inline void Z_powm(Z_t res, Z_t base, Z_t exp, Z_t mod)
{
   mpz_powm(res, base, exp, mod);
}

static inline void Z_powm_ui(Z_t res, Z_t base, unsigned long exp, Z_t mod)
{
   mpz_powm_ui(res, base, exp, mod);
}

unsigned long Z_invert_long(unsigned long, unsigned long);

static inline int Z_invert(Z_t res, Z_t a, Z_t p)
{
   return mpz_invert(res, a, p);
}

int Z_sqrtmod(Z_t, Z_t, Z_t);

void Z_sqrtmodpklift(Z_t, Z_t, Z_t, Z_t);

void Z_sqrtmodptopk(Z_t, Z_t, Z_t, Z_t, int);

int Z_sqrtmodpk(Z_t, Z_t, Z_t, int);

/*-----------------------------------------------------------------------------

    Exponentiation
    
-----------------------------------------------------------------------------*/

unsigned long Z_pow_long(unsigned long, unsigned long);

static inline void Z_pow_ui(Z_t res, Z_t base, unsigned long exp)
{
   mpz_pow_ui(res, base, exp);
}

static inline void Z_ui_pow_ui(Z_t res, unsigned long base, unsigned long exp)
{
   mpz_ui_pow_ui(res, base, exp);
}

/*-----------------------------------------------------------------------------

    Root Extraction
    
-----------------------------------------------------------------------------*/

static inline int Z_root(Z_t res, Z_t op, unsigned long n)
{
   return mpz_root(res, op, n);
}

/*static inline void Z_rootrem(Z_t root, Z_t rem, Z_t op, unsigned long n)
{
   mpz_rootrem(root, rem, op, n);
}*/

static inline void Z_sqrt(Z_t rop, Z_t op)
{
   mpz_sqrt(rop, op);
}

static inline void Z_sqrtrem(Z_t rop, Z_t rem, Z_t op)
{
   mpz_sqrtrem(rop, rem, op);
}

static inline int Z_perfect_power_p(Z_t op)
{
   return mpz_perfect_power_p(op);
}

static inline int Z_perfect_square_p(Z_t op)
{
   return mpz_perfect_square_p(op);
}

/*-----------------------------------------------------------------------------

    Number Theoretic
    
-----------------------------------------------------------------------------*/

unsigned long Z_nextprime_long(unsigned long);

static inline void Z_nextprime(Z_t res, Z_t op)
{
   mpz_nextprime(res, op);
}

void Z_randomprime(Z_t, unsigned long);

static inline int Z_probab_prime_p(Z_t n, int reps)
{
   return mpz_probab_prime_p(n, reps);
}

long Z_extgcd_long(long*, long*, long, long);

unsigned long Z_gcd_long(long, long);

static inline void Z_gcd(Z_t res, Z_t op1, Z_t op2)
{
   mpz_gcd(res, op1, op2);
}

static inline unsigned long Z_gcd_ui(Z_t res, Z_t op1, unsigned long op2)
{
   return mpz_gcd_ui(res, op1, op2);
}

static inline void Z_gcdext(Z_t g, Z_t s, Z_t t, Z_t a, Z_t b)
{
   mpz_gcdext(g, s, t, a, b);
}

static inline void Z_lcm(Z_t res, Z_t op1, Z_t op2)
{
   mpz_lcm(res, op1, op2);
}

static inline void Z_lcm_ui(Z_t res, Z_t op1, unsigned long op2)
{
   mpz_lcm_ui(res, op1, op2);
}

static inline int Z_jacobi(Z_t a, Z_t b)
{
   return mpz_jacobi(a, b);
}

static inline int Z_legendre(Z_t a, Z_t p)
{
   return mpz_legendre(a, p);
}

static inline int Z_kronecker(Z_t a, Z_t b)
{
   return mpz_kronecker(a, b);
}

static inline int Z_kronecker_si(Z_t a, long b)
{
   return mpz_kronecker_si(a, b);
}

static inline int Z_kronecker_ui(Z_t a, unsigned long b)
{
   return mpz_kronecker_ui(a, b);
}

static inline int Z_si_kronecker(long a, Z_t b)
{
   return mpz_si_kronecker(a, b);
}

static inline int Z_ui_kronecker(unsigned long a, Z_t b)
{
   return mpz_ui_kronecker(a, b);
}

static inline unsigned long Z_remove(Z_t res, Z_t op, Z_t n)
{
   return mpz_remove(res, op, n);
}

static inline void Z_fac_ui(Z_t res, unsigned long op)
{
   mpz_fac_ui(res, op);
}

static inline void Z_bin_ui(Z_t res, Z_t n, unsigned long k)
{
   mpz_bin_ui(res, n, k);
}

static inline void Z_bin_uiui(Z_t res, unsigned long n, unsigned long k)
{
   mpz_bin_uiui(res, n, k);
}

static inline void Z_fib_ui(Z_t fn, unsigned long n)
{
   mpz_fib_ui(fn, n);
}

static inline void Z_fib2_ui(Z_t fn, Z_t fnsub1, unsigned long n)
{
   mpz_fib2_ui(fn, fnsub1, n);
}

static inline void Z_lucnum_ui(Z_t ln, unsigned long n)
{
   mpz_lucnum_ui(ln, n);
}

static inline void Z_lucnum2_ui(Z_t ln, Z_t lnsub1, unsigned long n)
{
   mpz_lucnum2_ui(ln, lnsub1, n);
}

void Z_CRT(Z_t, Z_t, Z_t, Z_t, Z_t, Z_t);

/*-----------------------------------------------------------------------------

    Comparison Functions
    
-----------------------------------------------------------------------------*/

static inline int Z_cmp(Z_t op1, Z_t op2)
{
   return mpz_cmp(op1, op2);
}

static inline int Z_cmp_d(Z_t op1, double op2)
{
   return mpz_cmp_d(op1, op2);
}

static inline int Z_cmp_si(Z_t op1, signed long op2)
{
   return mpz_cmp_si(op1, op2);
}

static inline int Z_cmp_ui(Z_t op1, unsigned long op2)
{
  return  mpz_cmp_ui(op1, op2);
}

static inline int Z_cmpabs(Z_t op1, Z_t op2)
{
   return mpz_cmpabs(op1, op2);
}

static inline int Z_cmpabs_d(Z_t op1, double op2)
{
   return mpz_cmpabs_d(op1, op2);
}

static inline int Z_cmpabs_ui(Z_t op1, unsigned long op2)
{
   return mpz_cmpabs_ui(op1, op2);
}

static inline int Z_sgn(Z_t op)
{
   return mpz_sgn(op);
}

/*-----------------------------------------------------------------------------

    Logical and Bitwise Manipulation Functions
    
-----------------------------------------------------------------------------*/

static inline void Z_and(Z_t res, Z_t op1, Z_t op2)
{
   mpz_and(res, op1, op2);
}

static inline void Z_ior(Z_t res, Z_t op1, Z_t op2)
{
   mpz_ior(res, op1, op2);
}

static inline void Z_xor(Z_t res, Z_t op1, Z_t op2)
{
   mpz_xor(res, op1, op2);
}

static inline void Z_com(Z_t res, Z_t op)
{
   mpz_com(res, op);
}

static inline unsigned long Z_popcount(Z_t op)
{
   return mpz_popcount(op);
}

static inline unsigned long Z_hamdist(Z_t op1, Z_t op2)
{
   return mpz_hamdist(op1, op2);
}

static inline unsigned long Z_scan0(Z_t op, unsigned long starting_bit)
{
   return mpz_scan0(op, starting_bit);
}

static inline unsigned long Z_scan1(Z_t op, unsigned long starting_bit)
{
   return mpz_scan1(op, starting_bit);
}

static inline void Z_setbit(Z_t res, unsigned long bit_index)
{
   mpz_setbit(res, bit_index);
}

static inline void Z_clrbit(Z_t res, unsigned long bit_index)
{
   mpz_clrbit(res, bit_index);
}

/*static inline void Z_combit(Z_t res, unsigned long bit_index)
{
   mpz_combit(res, bit_index);
}*/

static inline int Z_tstbit(Z_t op, unsigned long bit_index)
{
   return mpz_tstbit(op, bit_index);
}

/*-----------------------------------------------------------------------------

    Input and Output Functions
    
-----------------------------------------------------------------------------*/

/*static inline size_t Z_out_str(FILE *stream, int base, Z_t op)
{
   return mpz_out_str(stream, base, op);
}

static inline size_t Z_inp_str(Z_t res, FILE *stream, int base)
{
   return mpz_inp_str(res, stream, base);
}

static inline size_t Z_out_raw(FILE *stream, Z_t op)
{
   return mpz_out_raw(stream, op);
}

static inline size_t Z_inp_raw(Z_t res, FILE *stream)
{
   return mpz_inp_raw(res, stream);
}*/

/*-----------------------------------------------------------------------------

    Random Number Generation
    
-----------------------------------------------------------------------------*/

void Z_urandomb(Z_t, unsigned long);

void Z_urandomm(Z_t, Z_t);

void Z_rrandomb(Z_t, unsigned long);

/*-----------------------------------------------------------------------------

    Import/Export Functions
    
-----------------------------------------------------------------------------*/

static inline void Z_import(Z_t rop, size_t count, int order, int size,
                     int endian, size_t nails, const void *op)
{
   mpz_import(rop, count, order, size, endian, nails, op);
}

static inline void * Z_export(void *rop, size_t *countp, int order, int size,
                       int endian, size_t nails, Z_t op)
{
   return mpz_export(rop, countp, order, size, endian, nails, op);
}

/*-----------------------------------------------------------------------------

    Miscellaneous Functions
    
-----------------------------------------------------------------------------*/

static inline int Z_fits_ulong_p(Z_t op)
{
   return mpz_fits_ulong_p(op);
}

static inline int Z_fits_slong_p(Z_t op)
{
   return mpz_fits_slong_p(op);
}

static inline int Z_fits_uint_p(Z_t op)
{
   return mpz_fits_uint_p(op);
}

static inline int Z_fits_sint_p(Z_t op)
{
   return mpz_fits_sint_p(op);
}

static inline int Z_fits_ushort_p(Z_t op)
{
   return mpz_fits_ushort_p(op);
}

static inline int Z_fits_sshort_p(Z_t op)
{
   return mpz_fits_sshort_p(op);
}

static inline int Z_odd_p(Z_t op)
{
   return mpz_odd_p(op);
}

static inline int Z_even_p(Z_t op)
{
   return mpz_even_p(op);
}

static inline size_t Z_sizeinbase(Z_t op, int base)
{
   return mpz_sizeinbase(op, base);
}

static inline size_t Z_size(Z_t op)
{
   return mpz_size(op);
}

/*===================================================================================

Montgomery routines

====================================================================================*/

unsigned long Z_mont_red(mpz_t res, mpz_t a, mpz_t m);

void Z_mont_mul(mpz_t res, mpz_t a, mpz_t b, mpz_t m, mpz_t R, unsigned long n);

void Z_expmod_mont(mpz_t res, mpz_t a, mpz_t exp, mpz_t m);

void Z_divrem_jebelean(mpz_t Q, mpz_t R, mpz_t A, mpz_t B);

void Z_rem_jebelean(mpz_t R, mpz_t A, mpz_t B);

void Z_mulmod_jebelean(mpz_t res, mpz_t a, mpz_t b, mpz_t m);

void Z_expmod_jebelean(mpz_t res, mpz_t a, mpz_t exp, mpz_t m);

#endif
