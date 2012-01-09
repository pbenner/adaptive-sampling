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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/logarithmetic.h>

typedef struct {
        size_t max_breaks;
        prob_t *result;
        prob_t *g;
        prob_t (*f)(int, int, void*);
        size_t L;
        void *data;
} approximate_multibin_t;

static
void approximate_multibin(
        size_t from,
        size_t pos,
        size_t n_breaks,
        prob_t partial_result,
        approximate_multibin_t* am)
{
        if (n_breaks > am->max_breaks+1) {
                return;
        }

        if (pos == am->L-1) {
                /* we arrived at a leaf */
                partial_result += (*am->f)(from, pos, am->data);
                partial_result += am->g[n_breaks];

                am->result[n_breaks] =
                        logadd(am->result[n_breaks], partial_result);
        }
        else {
                /* no break */
                approximate_multibin(from,  pos+1, n_breaks,   partial_result, am);
                /* break */
                partial_result += (*am->f)(from, pos, am->data);
                approximate_multibin(pos+1, pos+1, n_breaks+1, partial_result, am);
        }
}

void prombs_tree(
        prob_t *result,
        prob_t *g,
        prob_t (*f)(int, int, void*),
        size_t L,
        size_t m,
        void *data)
{
        approximate_multibin_t am = {m, result, g, f, L, data};
        size_t i;

        for (i = 0; i < L; i++) {
                result[i] = -HUGE_VAL;
        }

        approximate_multibin(0, 0, 0, 0, &am);
}
