/*============================================================================

    Copyright (C) 2007, William Hart
 
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
#include <stdio.h>

#include <gmp.h>
#include "mpz_extras.h"

/*int main(void)
{
   mpz_t a, exp, m, res;
   
   mpz_init(a);
   mpz_init(exp);
   mpz_init(m);
   mpz_init(res);
   
   printf("Enter a: ");
   gmp_scanf("%Zd", a);
   printf("Enter exponent: ");
   gmp_scanf("%Zd", exp);
   printf("Enter modulus: ");
   gmp_scanf("%Zd", m);
   
   F_mpz_expmod_mont(res, a, exp, m);
   
   gmp_printf("a^exp mod m is %Zd\n", res);
   
   mpz_clear(a);
   mpz_clear(exp);
   mpz_clear(m);
   mpz_clear(res);
}*/

int main(void)
{
   mpz_t a, exp, temp, p, res;
   
   mpz_init(a);
   mpz_init(exp);
   mpz_init(p);
   mpz_init(temp);
   mpz_init(res);
   
   mpz_set_ui(p, 1);
   mpz_mul_2exp(p, p, 29440);
   mpz_set_ui(temp, 1);
   mpz_mul_2exp(temp, temp, 27392);
   mpz_sub(p, p, temp);
   mpz_add_ui(p, p, 1);
   
   mpz_set_ui(a, 1);
   mpz_mul_2exp(a, a, 32);
   
   mpz_sub_ui(exp, p, 1);
   mpz_fdiv_q_2exp(exp, exp, 10);
   
   F_mpz_expmod_BZ(res, a, exp, p);
   //F_mpz_expmod_mont(res, a, exp, p);
   //mpz_powm(res, a, exp, p);
   
   mpz_clear(a);
   mpz_clear(exp);
   mpz_clear(p);
   mpz_init(temp);
   mpz_clear(res);   
}
