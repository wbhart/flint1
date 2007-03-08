/****************************************************************************

Zpoly.h: Polynomials over Z

Copyright (C) 2007, William Hart and David Harvey

!!!!! only PROPOSED data formats and interface

There are two entirely separate data formats for polynomials over Z:
  -- Zpoly_mpz_t uses an array of mpz_t's
  -- Zpoly_mpn_t uses a single block of memory with each coefficient occupying
     the same number of limbs.

*****************************************************************************/


//////////////////////////////////////////////////////////
/*
These allocation functions should eventually go in another file....

flint_malloc_limbs and flint_malloc use the same memory manager, they just
measure the size in different units

At first they would just be a wrapper for malloc/free, but eventually we
probably want an extra layer in between flint_malloc and malloc, which finds
free memory faster, possibly at the expense of greater fragmentation
*/

void* flint_malloc_limbs(unsigned long limbs);
void* flint_malloc(unsigned long bytes);
void* flint_realloc_limbs(void* block, unsigned long limbs);
void* flint_realloc(void* block, unsigned long bytes);
void flint_free(void* block);
//////////////////////////////////////////////////////////


/****************************************************************************

   Zpoly_mpz_t
   -----------

Zpoly_mpz_t represents a dense polynomial in Z[x] using a single mpz_t for each
coefficient.

"coeffs" is an array of mpz_t's of length "alloc". They are all mpz_init'd.
coeffs[n] is the coefficient of x^n.

Only the first "length" coefficients actually represent coefficients of the
polynomial; i.e. it's a polynomial of degree at most length-1. There is no
requirement for coeff[length-1] to be nonzero. If length == 0, this is the
zero polynomial. Obviously always alloc >= length.

There are two classes of functions operating on Zpoly_mpz_t:

-- The Zpoly_mpz_raw_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough mpz_t's for
   the result.

-- The Zpoly_mpz_* functions ASSUME that "coeffs" was allocated via
   flint_malloc, and they MAY free or reallocate "coeffs" using flint_realloc,
   flint_free etc, whenever they feel the need. Furthermore they assume that
   always alloc >= 1.

 */
typedef struct
{
   mpz_t* coeffs;
   unsigned long alloc;
   unsigned long length;
} Zpoly_mpz_struct;

// Zpoly_mpz_t allows reference-like semantics for Zpoly_mpz_struct:
typedef Zpoly_mpz_struct Zpoly_mpz_t[1];


// ============================================================================
// functions in Zpoly_mpz_raw_* layer

// returns x^n coefficient with no bounds checking
static inline
mpz_t* Zpoly_mpz_raw_get_coeff_ptr(Zpoly_mpz_t poly, unsigned long n)
{
    return &poly->coeffs[n];
}

// copies out x^n coefficient (via mpz_set) with no bounds checking
static inline
void Zpoly_mpz_raw_get_coeff(mpz_t output, Zpoly_mpz_t poly,
                             unsigned long n)
{
    mpz_set(output, poly->coeffs[n]);
}

// copies out x^n coefficient (via mpz_get_ui) with no bounds checking
// the coefficient is assumed to fit into an unsigned long
static inline
unsigned long Zpoly_mpz_raw_get_coeff_ui(Zpoly_mpz_t poly, unsigned long n)
{
    // todo: check mpz_get_ui is a macro; if not, write our own non-GMP_COMPLIANT version
    return mpz_get_ui(poly->coeffs[n]);
}

// todo: do we want a signed version of the above? i.e. Zpoly_mpz_raw_get_coeff_si

// sets x^n coefficient (via mpz_set) with no bounds checking
static inline
void Zpoly_mpz_raw_set_coeff(Zpoly_mpz_t poly, unsigned long n, mpz_t x)
{
    mpz_set(poly->coeffs[n], x);
}

static inline
void Zpoly_mpz_raw_set_coeff_ui(Zpoly_mpz_t poly, unsigned long n,
                                unsigned long x)
{
    mpz_set_ui(poly->coeffs[n], x);
}

static inline
void Zpoly_mpz_raw_set_coeff_si(Zpoly_mpz_t poly, unsigned long n, long x)
{
    mpz_set_si(poly->coeffs[n], x);
}


// Decreases poly.length to point at the last non-zero coefficient
// (i.e. so that the degree really is length-1)
void Zpoly_mpz_raw_normalise(Zpoly_mpz_t poly);

// output = input
// assumes output.alloc >= input.length
void Zpoly_mpz_raw_set(Zpoly_mpz_t output, Zpoly_mpz_t input);

// swaps coeffs, alloc, length
static inline
void Zpoly_mpz_raw_swap(Zpoly_mpz_t x, Zpoly_mpz_t y)
{
    mpz_t* temp1;
    unsigned long temp2;

    temp1 = y->coeffs;
    y->coeffs = x->coeffs;
    x->coeffs = coeffs_temp;
    
    temp2 = x->alloc;
    x->alloc = y->alloc;
    y->alloc = temp2;

    temp2 = x->length;
    x->length = y->length;
    y->length = temp2;
}


// output = input1 + input2
// Assumes output.alloc >= max(input1.length, input2.length)
// all combinations of parameter aliasing are allowed
// output.length is set to the maximum of the two lengths, and no normalisation
// is performed (so the output may have leading zeroes)
void Zpoly_mpz_raw_add(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);
void Zpoly_mpz_raw_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);

// output = -input
// assumes output.alloc >= input.length
// output may alias input
void Zpoly_mpz_raw_negate(Zpoly_mpz_t output, Zpoly_mpz_t input);


// scalar multiplication by x
void Zpoly_mpz_raw_scalar_mul(Zpoly_mpz_t poly, mpz_t x);
void Zpoly_mpz_raw_scalar_mul_ui(Zpoly_mpz_t poly, unsigned long x);
// scalar division by x
// todo: what about all the variations... floor, ceiling, truncating,
// exact division, etc?
void Zpoly_mpz_raw_scalar_div(Zpoly_mpz_t poly, mpz_t x);
void Zpoly_mpz_raw_scalar_div_ui(Zpoly_mpz_t poly, unsigned long x);


// output = input1 * input2
// assumes output.alloc >= input1.length + input2.length - 1
// output may NOT alias either input
void Zpoly_mpz_raw_mul(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);
// as above, but always uses naive multiplication algorithm
void Zpoly_mpz_raw_mul_naive(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                             Zpoly_mpz_t input2);
// as above, but always uses karatsuba
void Zpoly_mpz_raw_mul_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                                 Zpoly_mpz_t input2);

// output = input * input
// output may NOT alias input
void Zpoly_mpz_raw_sqr(Zpoly_mpz_t output, Zpoly_mpz_t input);
void Zpoly_mpz_raw_sqr_naive(Zpoly_mpz_t output, Zpoly_mpz_t input);
void Zpoly_mpz_raw_sqr_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input);


// output = x^n * input
// assumes output.alloc >= input.length + n
// output may alias input
void Zpoly_mpz_raw_left_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                              unsigned long n);
// output = input / x^n (truncating division)
// assumes output.alloc >= input.length - n >= 0
// output may alias input
void Zpoly_mpz_raw_right_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                               unsigned long n);


// todo: figure out better division interface once we understand the
// underlying algorithms a bit better.... in particular I'm not even sure to
// what extent these things make sense when the divisor is not monic.
// Perhaps we want specialised versions for the monic case.
// Also perhaps a power series inversion function should go here somewhere;
// i.e. computes inverse modulo a given x^n. And for that matter direct
// power series division, which might be faster than "invert & multiply".

// quotient = input1 / input2 (throw away remainder)
// asumes quotient.alloc >= max(0, input1.length - input2.length)
void Zpoly_mpz_raw_div(Zpoly_mpz_t quotient, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);
// remainder = input1 % input2 (throw away quotient)
// assumes remainder.alloc >= input2.length
void Zpoly_mpz_raw_rem(Zpoly_mpz_t remainder, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);
void Zpoly_mpz_raw_div_rem(Zpoly_mpz_t quotient, Zpoly_mpz_t remainder,
                           Zpoly_mpz_t input1, Zpoly_mpz_t input2);


// output = gcd(input1, input2)
// assumes output.alloc >= max(input1.length, input2.length)
// (or perhaps only requires output.alloc >= length of gcd?)
void Zpoly_mpz_raw_gcd(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);
// also sets a, b so that a*input1 + b*input2 = output
// (is this even always possible in Z[x]?)
void Zpoly_mpz_raw_xgcd(Zpoly_mpz_t a, Zpoly_mpz_t b, Zpoly_mpz_t output,
                        Zpoly_mpz_t input1, Zpoly_mpz_t input2);



// ============================================================================
// Zpoly_mpz_* layer ("non-raw" functions)

// note: the implementation of many of these functions will be essentially:
// check for available space in output, reallocate if necessary, then call
// the corresponding raw version


// allocate coeffs to length 1, sets length = 0 (i.e. zero polynomial)
void Zpoly_mpz_init(Zpoly_mpz_t poly);
// allocate coeffs to given length, sets length = 0 (i.e. zero polynomial)
void Zpoly_mpz_init2(Zpoly_mpz_t poly, unsigned long alloc);
// allocate coeffs to given length, with space for coeff_size bits in each
// coefficient, sets length = 0 (i.e. zero polynomial)
void Zpoly_mpz_init3(Zpoly_mpz_t poly, unsigned long alloc,
                     unsigned long coeff_bits);


// Changes allocated space to be alloc.
// Current value is preserved if it fits, or truncated if it doesn't fit.
void Zpoly_mpz_realloc(Zpoly_mpz_t poly, unsigned long alloc);

// todo: should have available a version of Zpoly_mpz_realloc() that also
// allows the user to specify the space allocated for new mpz's


// This is the actual implementation that's called for Zpoly_mpz_ensure_space()
// (see below) if a reallocation is required
void Zpoly_mpz_ensure_space_IMPL(Zpoly_mpz_t poly, unsigned long alloc);

// Ensures that the polynomial has at least alloc coefficients allocated.
// If the polynomial already has enough space allocated, nothing happens.
// If more space is required, then a realloc occurs, and the space allocated
// will at least *double* from its current amount. This strategy ensures that
// repeated calls have amortised constant cost.
// (NOTE: I'm making only the initial comparison inline, since this will
// happen very frequently; the actual reallocation will be less frequent,
// and will chew up too many bytes of code if I make it inline.)
static inline
void Zpoly_mpz_ensure_space(Zpoly_mpz_t poly, unsigned long alloc)
{
   if (poly->alloc < alloc)
      Zpoly_mpz_ensure_space_IMPL(poly, alloc);
}


// deallocates space (poly becomes uninitialised)
void Zpoly_mpz_clear(Zpoly_mpz_t poly);


// returns x^n coefficient, or NULL if out of range
mpz_t* Zpoly_mpz_get_coeff_ptr(Zpoly_mpz_t poly, unsigned long n);

// copies out x^n coefficient (via mpz_set), or retrieves zero if out of range
void Zpoly_mpz_get_coeff(mpz_t output, Zpoly_mpz_t poly,
                         unsigned long n);
unsigned long Zpoly_mpz_get_coeff_ui(Zpoly_mpz_t poly, unsigned long n);

// and the rest of them are just like the mpz_raw versions, but reallocate
// on demand

void Zpoly_mpz_set_coeff(Zpoly_mpz_t poly, unsigned long n, mpz_t x);
void Zpoly_mpz_set_coeff_ui(Zpoly_mpz_t poly, unsigned long n, unsigned long x);
void Zpoly_mpz_set_coeff_si(Zpoly_mpz_t poly, unsigned long n, long x);


// sets poly from a string, where the string is just a sequence of
// coefficients in decimal notation, separated by whitespace.
// An example string: "0 -3 45" represents 45x^2 - 3x
void Zpoly_mpz_set_from_string(Zpoly_mpz_t output, char* s);



void Zpoly_mpz_normalise(Zpoly_mpz_t poly);
void Zpoly_mpz_set(Zpoly_mpz_t output, Zpoly_mpz_t input);
void Zpoly_mpz_swap(Zpoly_mpz_t x, Zpoly_mpz_t y);
void Zpoly_mpz_add(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2);
void Zpoly_mpz_sub(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2);
void Zpoly_mpz_negate(Zpoly_mpz_t output, Zpoly_mpz_t input);
void Zpoly_mpz_scalar_mul(Zpoly_mpz_t poly, mpz_t x);
void Zpoly_mpz_scalar_mul_ui(Zpoly_mpz_t poly, unsigned long x);
void Zpoly_mpz_scalar_div(Zpoly_mpz_t poly, mpz_t x);
void Zpoly_mpz_scalar_div_ui(Zpoly_mpz_t poly, unsigned long x);
void Zpoly_mpz_mul(Zpoly_mpz_t output, Zpoly_mpz_t input1, Zpoly_mpz_t input2);
void Zpoly_mpz_mul_naive(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                         Zpoly_mpz_t input2);
void Zpoly_mpz_mul_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                             Zpoly_mpz_t input2);
void Zpoly_mpz_sqr(Zpoly_mpz_t output, Zpoly_mpz_t input);
void Zpoly_mpz_sqr_naive(Zpoly_mpz_t output, Zpoly_mpz_t input);
void Zpoly_mpz_sqr_karatsuba(Zpoly_mpz_t output, Zpoly_mpz_t input);

void Zpoly_mpz_left_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                          unsigned long n);
void Zpoly_mpz_right_shift(Zpoly_mpz_t output, Zpoly_mpz_t input,
                           unsigned long n);

void Zpoly_mpz_div(Zpoly_mpz_t quotient, Zpoly_mpz_t input1,
                   Zpoly_mpz_t input2);
void Zpoly_mpz_rem(Zpoly_mpz_t remainder, Zpoly_mpz_t input1,
                   Zpoly_mpz_t input2);
void Zpoly_mpz_div_rem(Zpoly_mpz_t quotient, Zpoly_mpz_t remainder,
                       Zpoly_mpz_t input1, Zpoly_mpz_t input2);

void Zpoly_mpz_gcd(Zpoly_mpz_t output, Zpoly_mpz_t input1,
                       Zpoly_mpz_t input2);
void Zpoly_mpz_xgcd(Zpoly_mpz_t a, Zpoly_mpz_t b, Zpoly_mpz_t output,
                    Zpoly_mpz_t input1, Zpoly_mpz_t input2);



/****************************************************************************

   Zpoly_mpn_t
   -----------

Zpoly_mpn_t represents a dense polynomial in Z[x] using a single block of
memory to hold all the coefficients.

This type is better suited to handling very dense polynomials with relatively
small coefficients, where the memory management overhead of Zpoly_mpz_t would
be too expensive.

"coeffs" is an array of limbs of length (alloc * (coeff_size+1)). Each
coefficient uses coeff_size+1 limbs. For each coefficient, the first limb is
a sign limb: 0 means positive and 1 means negative. (Zero may be stored as
either positive or negative.) The remaining "coeff_size" limbs represent the
absolute value of the coefficient, stored in GMP's mpn format.

Only the first "length" coefficients actually represent coefficients of the
polynomial; i.e. it's a polynomial of degree at most length-1. There is no
requirement for coeff[length-1] to be nonzero. If length == 0, this is the
zero polynomial. Obviously always alloc >= length.

There are two classes of functions operating on Zpoly_mpn_t:

-- The Zpoly_mpn_raw_* functions NEVER free or reallocate "coeffs", so they
   don't care how "coeffs" was allocated, and they never even look at the
   "alloc" attribute. They always assume the output has enough space for
   the result. They also NEVER modify the coeff_size attribute (since this
   would screw up the block size).

-- The Zpoly_mpn_* functions ASSUME that "coeffs" was allocated via
   flint_malloc, and they MAY free or reallocate "coeffs" using flint_realloc,
   flint_free etc, whenever they feel the need. Furthermore they assume that
   always alloc >= 1.

 */
typedef struct
{
   mp_limb_t* coeffs;
   unsigned long alloc;
   unsigned long length;
   unsigned long coeff_size;
} Zpoly_mpn_struct;

// Zpoly_mpn_t allows reference-like semantics for Zpoly_mpn_struct:
typedef Zpoly_mpn_struct Zpoly_mpn_t[1];


// ============================================================================
// functions in Zpoly_mpn_raw_* layer


// ... hmmmmmm ...


// ============================================================================
// functions in Zpoly_mpn_* layer (non-raw stuff)


// ... hmmmmmm ...


// *************** end of file
