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
#include <tools.h>

/******************************************************************************
 * Prombs effective counts functions
 ******************************************************************************/

static
prob_t effectivePosteriorCounts_f(int i, int j, void *data)
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
prob_t effectivePosteriorCounts(size_t pos, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];

        bp->counts_pos = pos;
        prombs(ev_log, bp->ak, bp->bd->prior_log, &effectivePosteriorCounts_f, bp->bd->L, minM(bp), (void *)bp);

        return EXP(sumModels(ev_log, bp) - evidence_ref);
}

static
prob_t effectiveCounts_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        if (i <= bp->counts_pos && bp->counts_pos <= j) {
                int k;
                prob_t n = 0;
                for (k = 0; k < bp->bd->events; k++) {
                        n += countStatistic(k, i, j, bp);
                }
                if (n == 0) {
                        return -HUGE_VAL;
                }
                else {
                        return log(n) + iec_log(i, j, bp);
                }
        }
        else {
                return iec_log(i, j, bp);
        }
}
static
prob_t effectiveCounts(size_t pos, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];
        prob_t sum;

        bp->counts_pos = pos;
        prombs(ev_log, bp->ak, bp->bd->prior_log, &effectiveCounts_f, bp->bd->L, minM(bp), (void *)bp);

        sum = sumModels(ev_log, bp);
        if (sum == -HUGE_VAL) {
                return 0;
        }
        else {
                return EXP(sum - evidence_ref);
        }
}

/******************************************************************************
 * Expectation
 ******************************************************************************/

static
prob_t expectation(
        size_t i,
        size_t j,
        prob_t evidence_ref,
        binProblem *bp)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bp->bd->L];

        bp->add_event.n     = 1;
        bp->add_event.which = i;
        bp->add_event.pos   = j;
        evidence_log        = evidence(evidence_log_tmp, bp);
        bp->add_event.pos   = -1;
        bp->add_event.n     = 0;

        return EXP(evidence_log - evidence_ref);
}

static
void computeExpectation(utility_t* result, size_t i, prob_t evidence_ref, binProblem *bp)
{
        size_t j;

        for (j = 0; j < bp->bd->events; j++) {
                result->expectation->content[j][i] = expectation(j, i, evidence_ref, bp);
        }
}

/******************************************************************************
 * Main
 ******************************************************************************/

void computeEffectiveCountsUtility(
        utility_t *result,
        prob_t evidence_ref,
        binData* bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        size_t i;

        for (i = 0; i < bd->L; i++) {
                notice(NONE, "Computing effective counts... %.1f%%", (float)100*(i+1)/bd->L);

                /* for each event compute its expectation */
                computeExpectation(result, i, evidence_ref, &bp);

                /* compute utilities */
                result->utility->content[i] = -effectiveCounts(i, evidence_ref, &bp);
        }

        binProblemFree(&bp);
}

void computeEffectivePosteriorCountsUtility(
        utility_t *result,
        prob_t evidence_ref,
        binData* bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        size_t i;

        for (i = 0; i < bd->L; i++) {
                notice(NONE, "Computing posterior effective counts... %.1f%%", (float)100*(i+1)/bd->L);

                /* for each event compute its expectation */
                computeExpectation(result, i, evidence_ref, &bp);

                /* compute utilities */
                result->utility->content[i] = -effectivePosteriorCounts(i, evidence_ref, &bp);
        }

        binProblemFree(&bp);
}
