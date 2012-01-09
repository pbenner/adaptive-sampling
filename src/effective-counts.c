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

/******************************************************************************
 * Prombs effective counts functions
 ******************************************************************************/

static
prob_t effectiveCounts_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        if (i <= bp->counts_pos && bp->counts_pos <= j) {
                int k;
                prob_t n = 0;
                for (k = 0; k < bp->bd->events; k++) {
                        n += countStatistic(k, i, j, bp) + countAlpha(k, i, j, bp);
                }
                return log(n) + iec_log(i, j, bp);
        }
        else {
                return iec_log(i, j, bp);
        }
}
static
prob_t effectiveCounts(size_t pos, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];

        bp->counts_pos = pos;
        prombs(ev_log, bp->ak, bp->bd->prior_log, &effectiveCounts_f, bp->bd->L, minM(bp), (void *)bp);

        return EXP(sumModels(ev_log, bp) - evidence_ref);
}

/******************************************************************************
 * Main
 ******************************************************************************/

void computeEffectiveCountsUtility(
        vector_t *result,
        prob_t evidence_ref,
        binData* bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        size_t i;

        for (i = 0; i < bd->L; i++) {
                notice(NONE, "Computing effective counts... %.1f%%", (float)100*(i+1)/bd->L);
                result->content[i] = -effectiveCounts(i, evidence_ref, &bp);
        }

        binProblemFree(&bp);
}
