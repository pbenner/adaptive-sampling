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
#include <math.h>

#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/logarithmetic.h>

/* Algorithm from D. Barry and J. A. Hartigan (1992)
 *
 * prombs_rec(bd.L, execPrombs_f, (void *)bp)
 */

static
prob_t prombs_rec_(
        prob_t *result,
        size_t j,
        prob_t (*f)(int, int, void*),
        void *data)
{
        if (j == 0) {
                return (*f)(0, 0, data);
        }
        else {
                size_t i;
                prob_t sum = -HUGE_VAL;

                for (i = 0; i < j; i++) {
                        if (result[i] == -HUGE_VAL) {
                                result[i] = prombs_rec_(result, i, f, data);
                        }
                        sum = logadd(sum, result[i] + f(i+1, j, data));
                }
                return logadd(sum, f(0, j, data));
        }
}

prob_t prombs_rec(
        size_t L,
        prob_t (*f)(int, int, void*),
        void *data)
{
        prob_t result[L];
        size_t i;

        for (i = 0; i < L; i++) {
                result[i] = -HUGE_VAL;
        }

        return prombs_rec_(result, L-1, f, data);
}
