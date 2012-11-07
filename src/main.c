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

#include <break-probabilities.h>
#include <datatypes.h>
#include <density.h>
#include <effective-counts.h>
#include <main-test.h>
#include <model.h>
#include <model-posterior.h>
#include <moment.h>
#include <utility.h>
#include <tools.h>

/******************************************************************************
 * Main binning function
 ******************************************************************************/

static
void copyModelPrior(binData* bd)
{
        size_t m_b;
        for (m_b = 0; m_b < bd->L; m_b++) {
                if (bd->beta->content[m_b] <= -HUGE_VAL) {
                        bd->prior_log[m_b] = (prob_t)-HUGE_VAL;
                }
                else {
                        bd->prior_log[m_b] = (prob_t)bd->beta->content[m_b];
                }
        }
}

static
void computeUtility(utility_t *result, prob_t evidence_ref, binData* bd)
{
        /* compute kl-divergence */
        if (bd->options->kl_psi || bd->options->kl_multibin) {
                computeKLUtility(result, evidence_ref, bd);
        }
        /* compute effective counts */
        else if (bd->options->effective_counts) {
                computeEffectiveCountsUtility(result, evidence_ref, bd);
        }
        /* compute effective posterior counts */
        else if (bd->options->effective_posterior_counts) {
                computeEffectivePosteriorCountsUtility(result, evidence_ref, bd);
        }
}

static
void computeBinning(
        marginal_t* result,
        binData *bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        prob_t evidence_ref;
        prob_t evidence_log_tmp[bd->L];

        /* init sampler */
        if (bd->options->algorithm == 1) {
                mgs_init(bd->options->samples[0], bd->options->samples[1],
                         bd->prior_log, &execPrombs_f, bd->L, (void *)&bp);
        }
        /* compute evidence P(D) */
        evidence_ref = evidence(evidence_log_tmp, &bp);
        /* compute model posteriors P(m_B|D) */
        if (bd->options->model_posterior) {
                computeModelPosteriors(evidence_log_tmp, result->mpost, evidence_ref, bd);
        }
        /* compute density */
        if (bd->options->density) {
                computeDensity(result->density, evidence_ref, bd);
        }
        /* compute break probability */
        if (bd->options->bprob) {
                computeBreakProbabilities(result->bprob, evidence_ref, bd);
        }
        /* compute the first n moments */
        if (bd->options->n_moments > 0) {
                computeMoments(result->moments, evidence_ref, bd);
        }

        if (bd->options->algorithm == 1) {
                mgs_free();
        }
        binProblemFree(&bp);
}

static
void computeHMM(
        marginal_t* result,
        binData *bd)
{
        binProblem bp; binProblemInit(&bp, bd);

        prob_t forward [bd->L];
        prob_t backward[bd->L];

        hmm_forward (forward,  &bp);
        hmm_backward(backward, &bp);

        /* compute the first n moments */
        if (bd->options->n_moments > 0) {
                hmm_computeMoments(result->moments, forward, backward, &bp);
        }
        /* compute density */
        if (bd->options->density) {
                hmm_computeDensity(result->density, forward, backward, &bp);
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
 * Binning init and free
 ******************************************************************************/

void bin_init(
        size_t events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        options_t* options,
        binData* bd)
{
        size_t L = counts[0]->columns;

        assert(options->which < events);

        verbose         = options->verbose;
        bd->options     = options;
        bd->L           = L;
        bd->events      = events;
        bd->counts      = counts;
        bd->alpha       = alpha;
        bd->beta        = beta;
        bd->gamma       = gamma;
        bd->prior_log   = (prob_t *)malloc(L*sizeof(prob_t));

        /* compute the model prior once for all computations */
        copyModelPrior(bd);
}

void bin_free(binData* bd)
{
        free(bd->prior_log);
}

/******************************************************************************
 * Library entry point
 ******************************************************************************/

marginal_t *
posterior(
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        options_t *options)
{
        binData bd;
        size_t L = counts[0]->columns;

        marginal_t *result = (marginal_t *)malloc(sizeof(marginal_t));
        result->moments    = (options->n_moments       ? alloc_matrix(options->n_moments, L)   : NULL);
        result->density    = (options->density         ? alloc_matrix(L, options->n_density)   : NULL);
        result->bprob      = (options->bprob           ? alloc_vector(L)                       : NULL);
        result->mpost      = (options->model_posterior ? alloc_vector(L)                       : NULL);

        bin_init(events, counts, alpha, beta, gamma, options, &bd);

        if (options->prombsTest) {
                prombsTest(&bd);
        }
        if (options->hmm) {
                computeHMM(result, &bd);
        }
        else {
                computeBinning(result, &bd);
        }
        bin_free(&bd);

        return result;
}

/*
 * Compute the expected utility N steps ahead
 */
utility_t*
utility(
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        options_t *options)
{
        binData bd;
        bin_init(events, counts, alpha, beta, gamma, options, &bd);
        binProblem bp; binProblemInit(&bp, &bd);

        utility_t *result   = (utility_t *)malloc(sizeof(utility_t));
        result->expectation = alloc_matrix(events, bp.bd->L);
        result->utility     = alloc_vector(bp.bd->L);

        if (options->hmm) {
                prob_t forward [bd.L];
                prob_t backward[bd.L];

                hmm_forward (forward,  &bp);
                hmm_backward(backward, &bp);
                hmm_computeUtility(result, forward, backward, &bp);
        }
        else {
                prob_t evidence_log_tmp[bd.L];
                prob_t evidence_ref = evidence(evidence_log_tmp, &bp);

                computeUtility(result, evidence_ref, &bd);
        }
        bin_free(&bd);

        return result;
}

/*
 * Compute the expected utility for sampling at position j
 * with a one-step look-ahead
 */
vector_t*
utilityAt(
        int pos,
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        options_t *options)
{
        binData bd;
        bin_init(events, counts, alpha, beta, gamma, options, &bd);
        binProblem bp; binProblemInit(&bp, &bd);

        /* result stores (expectation_1, ..., expectation_n, utility) */
        vector_t *result = alloc_vector(events+1);

        if (options->hmm) {
                prob_t forward [bd.L];
                prob_t backward[bd.L];

                hmm_forward (forward,  &bp);
                hmm_backward(backward, &bp);
                hmm_computeUtilityAt(pos, result, forward, backward, &bp);
        }
        else {
                warn(NONE, "This function is only implemented the hidden Markov model.");
        }
        bin_free(&bd);

        return result;
}

/*
 * Compute the Kullback-Leibler distance between the distribution
 * given by the counts and the one by adding the event (x,y)
 */
double
distance(
        int x,
        int y,
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        options_t *options)
{
        binData bd;
        bin_init(events, counts, alpha, beta, gamma, options, &bd);
        binProblem bp; binProblemInit(&bp, &bd);

        prob_t result = 0.0;

        if (options->hmm) {
                prob_t forward [bd.L];
                prob_t backward[bd.L];

                hmm_forward (forward,  &bp);
                hmm_backward(backward, &bp);
                result = hmm_computeDistance(x, y, forward, backward, &bp);
        }
        else {
                warn(NONE, "This function is only implemented the hidden Markov model.");
        }
        bin_free(&bd);

        return result;
}
