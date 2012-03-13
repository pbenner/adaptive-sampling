/* Copyright (C) 2010, 2011, 2012 Philipp Benner
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
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

#include <adaptive-sampling/exception.h>
#include <adaptive-sampling/logarithmetic.h>
#include <adaptive-sampling/mgs.h>
#include <adaptive-sampling/prombs.h>
#include <adaptive-sampling/datatypes.h>

#include <gsl/gsl_sf_gamma.h>

#include <bayesian-binning-test.h>
#include <break-probabilities.h>
#include <datatypes.h>
#include <effective-counts.h>
#include <marginal.h>
#include <model.h>
#include <model-posterior.h>
#include <moment.h>
#include <utility.h>
#include <tools.h>

/******************************************************************************
 * Main binning function
 ******************************************************************************/

static
void computeModelPrior(binData* bd)
{
        size_t m_b;
        for (m_b = 0; m_b < bd->L; m_b++) {
                if (bd->beta->content[m_b] == 0) {
                        bd->prior_log[m_b] = -HUGE_VAL;
                }
                else {
                        bd->prior_log[m_b] = -gsl_sf_lnchoose(bd->L-1, m_b) +
                                LOG(bd->beta->content[m_b]);
                }
        }
}

static
void computeUtility(vector_t *utility, prob_t evidence_ref, binData* bd)
{
        /* compute kl-divergence */
        if (bd->options->kl_component || bd->options->kl_multibin) {
                computeKLUtility(utility, evidence_ref, bd);
        }
        /* compute effective counts */
        else if (bd->options->effective_counts) {
                computeEffectiveCountsUtility(utility, evidence_ref, bd);
        }
        else if (bd->options->effective_posterior_counts) {
                computeEffectivePosteriorCountsUtility(utility, evidence_ref, bd);
        }
}

static
void computeBinning(
        BinningResult* result,
        binData *bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        prob_t evidence_ref;
        prob_t evidence_log_tmp[bd->L];

        /* init sampler */
        if (bd->options->algorithm == 2) {
                mgs_init(bd->options->samples[0], bd->options->samples[1],
                         bd->prior_log, &execPrombs_f, bd->L, (void *)&bp);
        }
        /* compute evidence P(D) */
        evidence_ref = evidence(evidence_log_tmp, &bp);
        /* compute model posteriors P(m_B|D) */
        if (bd->options->model_posterior) {
                computeModelPosteriors(evidence_log_tmp, result->mpost, evidence_ref, bd);
        }
        /* compute marginal */
        if (bd->options->marginal) {
                computeMarginal(result->marginals, evidence_ref, bd);
        }
        /* compute break probability */
        if (bd->options->bprob) {
                computeBreakProbabilities(result->bprob, evidence_ref, bd);
        }
        /* compute the first n moments */
        if (bd->options->n_moments > 0) {
                computeMoments(result->moments, evidence_ref, bd);
        }
        /* compute sampling utility */
        if (bd->options->utility && bd->options->algorithm == 0) {
                computeUtility(result->utility, evidence_ref, bd);
        }

        if (bd->options->algorithm == 2) {
                mgs_free();
        }
        binProblemFree(&bp);
}

/******************************************************************************
 * Initialization of common data structures
 ******************************************************************************/

void __init_rand__() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

void __init__(double epsilon)
{
        __init_rand__();
        __init_model__();
        __init_prombs__(epsilon);
}

void __free__() {
        __free_model__();
}

/******************************************************************************
 * Library entry point
 ******************************************************************************/

void bin_init(
        size_t events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        Options* options,
        binData* bd)
{
        size_t L = counts[0]->columns;

        assert(options->which < events);

        verbose       = options->verbose;
        bd->options   = options;
        bd->L         = L;
        bd->events    = events;
        bd->counts    = counts;
        bd->alpha     = alpha;
        bd->beta      = beta;
        bd->gamma     = gamma;
        bd->prior_log = (prob_t *)malloc(L*sizeof(prob_t));

        /* compute the model prior once for all computations */
        computeModelPrior(bd);
}

void bin_free(binData* bd)
{
        free(bd->prior_log);
}

BinningResult *
binning(
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        Options *options)
{
        binData bd;
        size_t L = counts[0]->columns;

        BinningResult *result = (BinningResult *)malloc(sizeof(BinningResult));
        result->moments   = (options->n_moments       ? alloc_matrix(options->n_moments, L)   : NULL);
        result->marginals = (options->marginal        ? alloc_matrix(L, options->n_marginals) : NULL);
        result->bprob     = (options->bprob           ? alloc_vector(L)                       : NULL);
        result->mpost     = (options->model_posterior ? alloc_vector(L)                       : NULL);
        result->utility   = (options->utility         ? alloc_vector(L)                       : NULL);

        bin_init(events, counts, alpha, beta, gamma, options, &bd);

        if (options->prombsTest) {
                prombsTest(&bd);
        }
        computeBinning(result, &bd);

        bin_free(&bd);

        return result;
}
