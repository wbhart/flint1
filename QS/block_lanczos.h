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
#ifndef LANCZOS_H
#define LANCZOS_H

#define DISPLAY 0 // Display some info about the linear algebra phase

#include <stdlib.h>
#include <stdint.h>

#include "../memory-manager.h"

typedef struct {
	unsigned long *data;		/* The list of occupied rows in this column */
	unsigned long weight;		/* Number of nonzero entries in this column */
	unsigned long orig;         /* Original relation number */
} la_col_t;

uint64_t get_null_entry(uint64_t *, long, long);

void reduce_matrix(unsigned long *, unsigned long *, la_col_t *);

uint64_t * block_lanczos(unsigned long, unsigned long, unsigned long, la_col_t*);

/*==========================================================================
   insert_col_entry:

   Function: insert an entry into a column of the matrix, 
   reallocating the space for the column if necessary

===========================================================================*/
static inline void insert_col_entry(la_col_t* col, unsigned long entry)
{
   unsigned long* temp;
   
   if (((col->weight >> 4) << 4) == col->weight) //need more space
   {
      if (col->weight != 0) col->data = (unsigned long *) flint_heap_realloc(col->data, col->weight+16);
      else col->data = (unsigned long*) flint_heap_alloc(16);
   }
   
   col->data[col->weight] = entry;
   col->weight++;
}

/*==========================================================================
   copy_col:

   Function: clear a column
   
===========================================================================*/

static inline void copy_col(la_col_t* col2, la_col_t* col1)
{
   col2->weight = col1->weight;
   col2->data = col1->data;
   col2->orig = col1->orig;
}

/*==========================================================================
   swap_cols:

   Function: swap two columns
   
===========================================================================*/

static inline void swap_cols(la_col_t* col2, la_col_t* col1)
{
   la_col_t temp;
   
   temp.weight = col1->weight;
   temp.data = col1->data;
   temp.orig = col1->orig;
   
   col1->weight = col2->weight;
   col1->data = col2->data;
   col1->orig = col2->orig;
   
   col2->weight = temp.weight;
   col2->data = temp.data;
   col2->orig = temp.orig;
}

/*==========================================================================
   clear_col:

   Function: clear a column
   
===========================================================================*/

static inline void clear_col(la_col_t* col)
{
   col->weight = 0;
}

/*==========================================================================
   free_col:

   Function: free the memory used by a column
   
===========================================================================*/

static inline void free_col(la_col_t* col)
{
   if (col->weight) flint_heap_free(col->data);
}

#endif
