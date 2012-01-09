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

#include <gsl/gsl_sf_psi.h>

#include <datatypes.h>
#include <model.h>
#include <utility.h>

/******************************************************************************
 * Prombs entropy functions
 ******************************************************************************/

static
prob_t singlebinEntropy(int i, int j, binProblem *bp)
{
        size_t k;
        prob_t n = 0;
        prob_t c[bp->bd->events];
        prob_t sum = 0;

        for (k = 0; k < bp->bd->events; k++) {
                c[k] = countStatistic(k, i, j, bp) + countAlpha(k, i, j, bp);
                sum += (c[k] - 1.0)*gsl_sf_psi(c[k]);
                n   +=  c[k];
        }

        return mbeta_log(c, bp) + (n - bp->bd->events)*gsl_sf_psi(n) - sum;
}

static
prob_t differentialEntropy_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(i, j, bp);
}
static
prob_t differentialEntropy_h(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return -singlebinEntropy(i, j, bp);
}

static
prob_t differentialEntropy(prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];
        prob_t sum;

        prombsExt(ev_log, bp->ak, bp->bd->prior_log, &differentialEntropy_f, &differentialEntropy_h, bp->bd->L, minM(bp), (void *)bp);

        sum = sumModels(ev_log, bp);
        if (sum == -HUGE_VAL) {
                return 0.0;
        }
        else {
                return -EXP(sum - evidence_ref);
        }
}

static prob_t multibinEntropy_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(i, j, bp);
}
static prob_t multibinEntropy_h(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        /* - Log[f(b)] */
        return -iec_log(i, j, bp);
}

static
prob_t multibinEntropy(prob_t evidence_ref, binProblem *bp)
{
        size_t i;
        prob_t g[bp->bd->L];
        prob_t result1[bp->bd->L];
        prob_t result2[bp->bd->L];
        prob_t result[bp->bd->L];
        prob_t sum;

        for (i = 0; i < bp->bd->L; i++) {
                g[i] = LOG(bp->bd->prior_log[i] - evidence_ref) + bp->bd->prior_log[i];
        }
        prombsExt(result1, bp->ak, bp->bd->prior_log, &multibinEntropy_f, &multibinEntropy_h, bp->bd->L, minM(bp), (void *)bp);
        prombs(result2, bp->ak, g, &multibinEntropy_f, bp->bd->L, minM(bp), (void *)bp);

        for (i = 0; i < minM(bp)+1; i++) {
                result[i] = logsub(result1[i], result2[i]);
        }

        sum = sumModels(result, bp);
        if (sum == -HUGE_VAL) {
                return 0.0;
        }
        else {
                return EXP(sum - evidence_ref);
        }
}

/******************************************************************************
 * Main
 ******************************************************************************/

prob_t computeDifferentialEntropy(prob_t evidence_ref, binData* bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        prob_t differential_entropy;

        differential_entropy = differentialEntropy(evidence_ref, &bp);
        binProblemFree(&bp);

        return differential_entropy;
}

prob_t computeMultibinEntropy(prob_t evidence_ref, binData* bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        prob_t multibin_entropy;

        multibin_entropy = multibinEntropy(evidence_ref, &bp);
        binProblemFree(&bp);

        return multibin_entropy;
}

prob_t computeEntropy(prob_t evidence_ref, binData* bd)
{
        prob_t entropy = 0;

        if (bd->options->differential_entropy) {
                entropy += computeDifferentialEntropy(evidence_ref, bd);
        }
        if (bd->options->multibin_entropy) {
                entropy += computeMultibinEntropy(evidence_ref, bd);
        }

        return entropy;
}

void computeEntropicUtility(
        vector_t *result,
        prob_t evidence_ref,
        binData* bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        size_t i, j;
        prob_t differential_entropy;
        prob_t multibin_entropy;
        prob_t expectation;
        prob_t evidence_log;
        prob_t evidence_log_tmp[bd->L];

        bp.add_event.n = 1;
        for (i = 0; i < bd->L; i++) {
                /* for all items */
                notice(NONE, "Computing utilities... %.1f%%", (float)100*(i+1)/bd->L);
                bp.add_event.pos = i;
                result->content[i] = 0;
                for (j = 0; j < bd->events; j++) {
                        /* for all events */
                        bp.add_event.which = j;
                        /* recompute the evidence */
                        evidence_log = evidence(evidence_log_tmp, &bp);
                        expectation  = EXP(evidence_log - evidence_ref);
                        /* compute entropies */
                        if (bd->options->differential_entropy) {
                                differential_entropy = differentialEntropy(evidence_log, &bp);
                                result->content[i] += expectation*differential_entropy;
                        }
                        if (bd->options->multibin_entropy) {
                                multibin_entropy = multibinEntropy(evidence_log, &bp);
                                result->content[i] += expectation*multibin_entropy;
                        }
                }
                /* entropy -> utility */
                result->content[i] = -result->content[i];
        }

        binProblemFree(&bp);
}
