/* Copyright (C) 2011 Philipp Benner
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

#ifndef MGS_H
#define MGS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stddef.h>

#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/logarithmetic.h>

typedef struct {
        size_t from;
        size_t to;
} bin_t;

typedef struct {
        size_t n;        /* length of breaks[] */
        size_t n_breaks; /* number of possible breaks */
        size_t n_bins;   /* current number of bins */
        uint32_t* breaks;
} multibin_t;

multibin_t* new_multibin(size_t L);
multibin_t* clone_multibin(multibin_t* multibin);
void free_multibin(multibin_t* multibin);
void insert_break(multibin_t* multibin, size_t i);
void remove_break(multibin_t* multibin, size_t i);
void switch_break(multibin_t* multibin, size_t i);
void print_mutlibin(multibin_t* mutlibin);
void get_bins(multibin_t* multibin, bin_t* bins);
void get_breaks(multibin_t* multibin, size_t *breaks);

#endif /* MGS_H */
