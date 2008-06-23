#ifndef FMPZ_MONTGOMERY_H
#define FMPZ_MONTGOMERY_H

typedef struct
{
    fmpz_t M; // Holds preconditioned modulus M
    fmpz_t R; // Holds 1/M  mod  2^#
    fmpz_t t; // Holds 2^#
    fmpz_t b; // Holds extra value for mulmod, divmod, mod
} fmpz_montgomery_t[1];


void fmpz_montgomery_init(fmpz_montgomery_t mont, fmpz_t m);
void fmpz_montgomery_clear(fmpz_montgomery_t mont);
void fmpz_montgomery_redc(fmpz_t res, fmpz_t x, fmpz_montgomery_t mont);
void fmpz_montgomery_mulmod_init(fmpz_montgomery_t mont, fmpz_t b, fmpz_t m);
void fmpz_montgomery_divmod_init(fmpz_montgomery_t mont, fmpz_t b, fmpz_t m);
void fmpz_montgomery_mod_init(fmpz_montgomery_t mont, fmpz_t m);
void fmpz_montgomery_mulmod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont);

// montgomery precomputed divmod (same as mulmod)
static inline
void fmpz_montgomery_divmod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont)
{
    fmpz_montgomery_mulmod(res, a, mont);
}

// montgomery precomputed mod (same as mulmod)
static inline
void fmpz_montgomery_mod(fmpz_t res, fmpz_t a, fmpz_montgomery_t mont)
{
    fmpz_montgomery_mulmod(res, a, mont);
}

#endif
