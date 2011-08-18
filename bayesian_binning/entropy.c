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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>
#include <sys/time.h>

#include <bayes/exception.h>
#include <bayes/logarithmetic.h>
#include <bayes/mgs.h>
#include <bayes/prombs.h>
#include <bayes/datatypes.h>
#include <bayes/uthash.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_psi.h>

#include <datatypes.h>
#include <model.h>
#include <utility.h>

////////////////////////////////////////////////////////////////////////////////
// Prombs entropy functions
////////////////////////////////////////////////////////////////////////////////

static
prob_t singlebinEntropy(binProblem *bp, int i, int j)
{
        unsigned int k;
        prob_t n = 0;
        prob_t c[bd.events];
        prob_t sum = 0;

        for (k = 0; k < bd.events; k++) {
                c[k] = countStatistic(bp, k, i, j) + countAlpha(k, i, j);
                sum += (c[k] - 1.0)*gsl_sf_psi(c[k]);
                n   +=  c[k];
        }

        return mbeta_log(c) + (n - bd.events)*gsl_sf_psi(n) - sum;
}

static
prob_t differentialEntropy_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(bp, i, j);
}
static
prob_t differentialEntropy_h(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return -singlebinEntropy(bp, i, j);
}

static
prob_t differentialEntropy(binProblem *bp, prob_t evidence_ref)
{
        prob_t ev_log[bd.L];
        prob_t sum;

        prombsExt(ev_log, bp->ak, bd.prior_log, &differentialEntropy_f, &differentialEntropy_h, bd.L, minM(), (void *)bp);

        sum = sumModels(ev_log);
        if (sum == -HUGE_VAL) {
                return 0.0;
        }
        else {
                return -expl(sum - evidence_ref);
        }
}

static prob_t multibinEntropy_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(bp, i, j);
}
static prob_t multibinEntropy_h(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        // - Log[f(b)]
        return -iec_log(bp, i, j);
}

static
prob_t multibinEntropy(binProblem *bp, prob_t evidence_ref)
{
        size_t i;
        prob_t g[bd.L];
        prob_t result1[bd.L];
        prob_t result2[bd.L];
        prob_t result[bd.L];
        prob_t sum;

        for (i = 0; i < bd.L; i++) {
                g[i] = logl(bd.prior_log[i] - evidence_ref) + bd.prior_log[i];
        }
        prombsExt(result1, bp->ak, bd.prior_log, &multibinEntropy_f, &multibinEntropy_h, bd.L, minM(), (void *)bp);
        prombs(result2, bp->ak, g, &multibinEntropy_f, bd.L, minM(), (void *)bp);

        for (i = 0; i < bd.L; i++) {
                result[i] = logsub(result1[i], result2[i]);
        }

        sum = sumModels(result);
        if (sum == -HUGE_VAL) {
                return 0.0;
        }
        else {
                return expl(sum - evidence_ref);
        }
}

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////

void computeEntropicUtility(
        prob_t *result,
        prob_t evidence_ref)
{
        binProblem bp; binProblemInit(&bp);
        unsigned int i, j;
        prob_t differential_entropy;
        prob_t multibin_entropy;
        prob_t expectation;
        prob_t evidence_log;
        prob_t evidence_log_tmp[bd.L];

        bp.add_event.n = 1;
        for (i = 0; i < bd.L; i++) {
                // for all items
                notice(NONE, "Computing utilities... %.1f%%", (float)100*(i+1)/bd.L);
                bp.add_event.pos = i;
                result[i] = 0;
                for (j = 0; j < bd.events; j++) {
                        // for all events
                        bp.add_event.which = j;
                        // recompute the evidence
                        evidence_log = evidence(&bp, evidence_log_tmp);
                        expectation  = expl(evidence_log - evidence_ref);
                        // compute entropies
                        differential_entropy = differentialEntropy(&bp, evidence_log);
                        multibin_entropy     = multibinEntropy(&bp, evidence_log);
                        // initialize sum
//                        printf("differential entropy: %f\n", (float)differential_entropy);
//                        printf("multibin     entropy: %f\n", (float)multibin_entropy);
//                        result[i] += expectation*(differential_entropy - multibin_entropy);
//                        result[i] += expectation*differential_entropy;
                        result[i] += expectation*multibin_entropy;
                }
                result[i] = -result[i];
        }

        binProblemFree(&bp);
}
