/* Copyright (C) 2010 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>

#include <adaptive-sampling/exception.h>
#include <adaptive-sampling/linalg.h>

#include "interface.h"
#include "bayesian-binning.h"

// export functions for interface.py
vector_t * _alloc_vector(int size)              { return alloc_vector(size); }
void       _free_vector(vector_t *v)            { free_vector(v); }
matrix_t * _alloc_matrix(int rows, int columns) { return alloc_matrix(rows, columns); }
void       _free_matrix(matrix_t *m)            { free_matrix(m); }
void       _free(void *ptr)                     { free(ptr); }

void _init_(double epsilon)
{
        __init__(epsilon);
}

void _free_()
{
        __free__();
}
