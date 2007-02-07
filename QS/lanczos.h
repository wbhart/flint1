/*============================================================================
    Copyright 2006 William Hart    

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

============================================================================*/

#include <stdlib.h>

typedef struct {
	u_int32_t *data;		/* The list of occupied rows in this column */
	u_int32_t weight;		/* Number of nonzero entries in this column */
	u_int32_t orig;         /* Original relation number */
} la_col_t;

u_int64_t getNullEntry(u_int64_t *, int32_t, int32_t);
void reduce_matrix(u_int32_t *, u_int32_t *, la_col_t *);
u_int64_t * block_lanczos(u_int32_t, u_int32_t, u_int32_t, la_col_t*);
