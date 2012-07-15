/* Copyright (C) 2010, 2011 Philipp Benner
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

#ifndef MODEL_H
#define MODEL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <datatypes.h>
#include <adaptive-sampling/datatypes.h>

void __init_model__();
void __free_model__();

prob_t mbeta_log(prob_t *p, binProblem *bp);
prob_t iec_log(int kk, int k, binProblem *bp);

void hmm_forward(vector_t *result, binProblem* bp);
void hmm_backward(vector_t *result, binProblem* bp);
void hmm_fb(vector_t *result, vector_t *forward, vector_t *backward, prob_t (*f)(int, int, binProblem*), binProblem* bp);

#endif /* MODEL_H */
