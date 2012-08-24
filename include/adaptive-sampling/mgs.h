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

#ifndef ADAPTIVE_SAMPLING_MGS_H
#define ADAPTIVE_SAMPLING_MGS_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/probtype.h>

void mgs(prob_t *result, prob_t *g, prob_t (*f)(int, int, void*), size_t L, void *data);
void mgs_init(size_t R, size_t N, prob_t *g, prob_t (*f)(int, int, void*), size_t L, void *data);
void mgs_free();
size_t * mgs_get_counts();
void mgs_get_bprob(vector_t *bprob, size_t L);

#endif /* ADAPTIVE_SAMPLING_MGS_H */
