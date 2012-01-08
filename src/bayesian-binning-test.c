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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>
#include <sys/time.h>

#include <adaptive-sampling/exception.h>
#include <adaptive-sampling/logarithmetic.h>
#include <adaptive-sampling/mgs.h>
#include <adaptive-sampling/prombs.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/uthash.h>

#include <datatypes.h>
#include <model.h>
#include <utility.h>

////////////////////////////////////////////////////////////////////////////////
// Prombs test functions
////////////////////////////////////////////////////////////////////////////////

static
prob_t prombsTest_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(i, j, bp);
}

static
prob_t prombsTest_h(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        // - Log[f(b)]
        return -iec_log(i, j, bp);
}

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////

void prombsTest(binData *bd)
{
        MET_INIT;
        prob_t result1[bd->L];
        prob_t result2[bd->L];
        prob_t sum;
        size_t i;
        matrix_t *ak = alloc_matrix(bd->L, bd->L);

        // set prior to 1
        for (i = 0; i < bd->L; i++) {
                bd->prior_log[i] = 0;
        }
        MET("Testing prombs",
            prombs   (result1, ak, bd->prior_log, &prombsTest_f, bd->L, bd->L-1, NULL));
        MET("Testing prombsExt",
            prombsExt(result2, ak, bd->prior_log, &prombsTest_f, &prombsTest_h, bd->L, bd->L-1, NULL));

        sum = -HUGE_VAL;
        for (i = 0; i < bd->L; i++) {
                (void)printf("prombs[%02du]: %.10f\n", (unsigned int)i, (double)result1[i]);
                sum = logadd(sum, result1[i]);
        }
        (void)printf("prombs: %.10f\n", (double)sum);

        sum = -HUGE_VAL;
        for (i = 0; i < bd->L; i++) {
                (void)printf("prombsExt[%02du]: %.10f\n", (unsigned int)i, (double)result2[i]);
                sum = logadd(sum, result2[i]);
        }
        (void)printf("prombsExt: %.10f\n", (double)sum);

        free_matrix(ak);
}
