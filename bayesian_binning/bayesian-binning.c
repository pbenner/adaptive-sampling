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

#include <bayesian-binning-test.h>
#include <break-probabilities.h>
#include <datatypes.h>
#include <effective-counts.h>
#include <entropy.h>
#include <marginal.h>
#include <model.h>
#include <model-posterior.h>
#include <moment.h>
#include <utility.h>

////////////////////////////////////////////////////////////////////////////////
// Global variables
////////////////////////////////////////////////////////////////////////////////

binData bd;

////////////////////////////////////////////////////////////////////////////////
// Main binning function
////////////////////////////////////////////////////////////////////////////////

static
void computeModelPrior()
{
        unsigned int m_b;
        for (m_b = 0; m_b < bd.L; m_b++) {
                if (gsl_vector_get(bd.beta, m_b) == 0) {
                        bd.prior_log[m_b] = -HUGE_VAL;
                }
                else {
                        bd.prior_log[m_b] = -gsl_sf_lnchoose(bd.L-1, m_b) +
                                logl(gsl_vector_get(bd.beta, m_b));
                }
        }
}

static
void computeUtility(prob_t *utility, prob_t evidence_ref)
{
        if (bd.options->differential_entropy || bd.options->multibin_entropy) {
                computeEntropicUtility(utility, evidence_ref);
        }
        // compute effective counts
        else if (bd.options->effective_counts) {
                computeEffectiveCountsUtility(utility, evidence_ref);
        }
}

static
void computeBinning(
        prob_t **moments,
        prob_t **marginals,
        prob_t  *bprob,
        prob_t  *mpost,
        prob_t  *utility)
{
        binProblem bp; binProblemInit(&bp);
        prob_t evidence_ref;
        prob_t evidence_log_tmp[bd.L];

        // init sampler
        if (bd.options->algorithm == 2) {
                mgs_init(bd.options->samples[0], bd.options->samples[1],
                         bd.prior_log, &execPrombs_f, bd.L, (void *)&bp);
        }
        // compute evidence P(D)
        evidence_ref = evidence(&bp, evidence_log_tmp);
        // compute model posteriors P(m_B|D)
        if (bd.options->model_posterior) {
                computeModelPosteriors(evidence_log_tmp, mpost, evidence_ref);
        }
        // compute moments
        if (bd.options->marginal) {
                computeMarginal(marginals, evidence_ref);
        }
        // compute break probability
        if (bd.options->bprob) {
                computeBreakProbabilities(bprob, evidence_ref);
        }
        // compute the first n moments
        if (bd.options->n_moments > 0) {
                computeMoments(moments, evidence_ref);
        }
        // compute sampling utility
        if (bd.options->utility) {
                computeUtility(utility, evidence_ref);
        }

        if (bd.options->algorithm == 2) {
                mgs_free();
        }
        binProblemFree(&bp);
}

////////////////////////////////////////////////////////////////////////////////
// Initialization of common data structures
////////////////////////////////////////////////////////////////////////////////

void __init_rand__() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
}

void __init__(
        size_t events,
        gsl_matrix **counts,
        gsl_matrix **alpha,
        gsl_vector  *beta,
        gsl_matrix  *gamma,
        Options* options)
{
        size_t L = counts[0]->size2;

        __init_rand__();
        __init_model__();

        prombs_init(options->epsilon);

        verbose      = options->verbose;
        bd.options   = options;
        bd.L         = L;
        bd.events    = events;
        bd.counts    = counts;
        bd.alpha     = alpha;
        bd.beta      = beta;
        bd.gamma     = gamma;
        bd.prior_log = (prob_t *)malloc(L*sizeof(prob_t));

        // compute the model prior once for all computations
        computeModelPrior();
}

void __free__() {
        __free_model__();

        free(bd.prior_log);
}

////////////////////////////////////////////////////////////////////////////////
// Library entry point
////////////////////////////////////////////////////////////////////////////////

prob_t
bin_entropy(
        size_t events,
        gsl_matrix **counts,
        gsl_matrix **alpha,
        gsl_vector  *beta,
        gsl_matrix  *gamma,
        Options *options)
{
        __init__(events, counts, alpha, beta, gamma, options);

        prob_t entropy;
        prob_t evidence_ref;
        prob_t evidence_ref_tmp[bd.L];

        binProblem bp; binProblemInit(&bp);

        evidence_ref = evidence(&bp, evidence_ref_tmp);
        entropy = computeEntropy(evidence_ref);

        binProblemFree(&bp);

        __free__();

        return entropy;
}

BinningResultGSL *
bin_log(
        size_t events,
        gsl_matrix **counts,
        gsl_matrix **alpha,
        gsl_vector  *beta,
        gsl_matrix  *gamma,
        Options *options)
{
        size_t L = counts[0]->size2;
        BinningResultGSL *result = (BinningResultGSL *)malloc(sizeof(BinningResultGSL));
        result->moments   = (options->n_moments ?
                             gsl_matrix_alloc(options->n_moments, L)   : NULL);
        result->marginals = (options->marginal  ?
                             gsl_matrix_alloc(L, options->n_marginals) : NULL);
        result->bprob     = gsl_vector_alloc(L);
        result->mpost     = gsl_vector_alloc(L);
        result->utility   = gsl_vector_alloc(L);
        prob_t * moments[options->n_moments];
        prob_t * marginals[L];
        prob_t bprob[L];
        prob_t mpost[L];
        prob_t utility[L];
        unsigned int i, j;

        for (i = 0; i < options->n_moments; i++) {
                moments[i]   = (prob_t *)calloc(L, sizeof(prob_t));
        }
        for (i = 0; i < L; i++) {
                marginals[i] = (prob_t *)calloc(options->n_marginals, sizeof(prob_t));
        }
        bzero(bprob,   L*sizeof(prob_t));
        bzero(mpost,   L*sizeof(prob_t));
        bzero(utility, L*sizeof(prob_t));

        __init__(events, counts, alpha, beta, gamma, options);

        if (options->prombsTest) {
                prombsTest();
        }
        computeBinning(moments, marginals, bprob, mpost, utility);
        for (i = 0; i <= bd.L-1; i++) {
                for (j = 0; j < options->n_moments; j++) {
                        gsl_matrix_set(result->moments, j, i, moments[j][i]);
                }
                if (options->marginal) {
                        for (j = 0; j < options->n_marginals; j++) {
                                gsl_matrix_set(result->marginals, i, j, marginals[i][j]);
                        }
                }
                gsl_vector_set(result->bprob,   i, bprob[i]);
                gsl_vector_set(result->mpost,   i, mpost[i]);
                gsl_vector_set(result->utility, i, utility[i]);
        }

        for (i = 0; i < options->n_moments; i++) {
                free(moments[i]);
        }
        for (i = 0; i < L; i++) {
                free(marginals[i]);
        }

        __free__();

        return result;
}
