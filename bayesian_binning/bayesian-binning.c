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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_psi.h>

////////////////////////////////////////////////////////////////////////////////
// Data structures
////////////////////////////////////////////////////////////////////////////////

// data that has to be immutable
typedef struct {
        Options *options;
        // number of timesteps
        unsigned int L;
        unsigned int events;
        prob_t *prior_log;     // P(p,B|m_B)
        // counts and parameters
        gsl_matrix **counts;
        gsl_matrix **alpha;
        gsl_vector  *beta;     // P(m_B)
        gsl_matrix  *gamma;
} binData;

// mutable data, local to each thread
typedef struct {
        // temporary memory for prombs
        Matrix *ak;
        // break probability
        int bprob_pos;
        // effective counts
        int counts_pos;
        // moments
        struct {
                int pos;
                int n;
                int which;
        } add_event;
        // marginals
        struct {
                int pos;
                prob_t val;
                int which;
        } fix_prob;
} binProblem;

static binData bd;

////////////////////////////////////////////////////////////////////////////////
// Count statistics
////////////////////////////////////////////////////////////////////////////////

static
unsigned int countStatistic(binProblem *bp, unsigned int event, int ks, int ke)
{
        if (bp != NULL && bp->add_event.which == event &&
            ks <= bp->add_event.pos && bp->add_event.pos <= ke) {
                return gsl_matrix_get(bd.counts[event], ks, ke) +
                        bp->add_event.n;
        }
        else if (ks <= ke) {
                return gsl_matrix_get(bd.counts[event], ks, ke);
        }
        else {
                return 0;
        }
}

static
prob_t countAlpha(unsigned int event, int ks, int ke)
{
        if (ks <= ke) {
                return gsl_matrix_get(bd.alpha[event], ks, ke);
        }
        else {
                return 0;
        }
}

////////////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////////////

static inline
prob_t sumModels(prob_t *ev_log)
{
        // sum up the vector returned by prombs
        prob_t sum = -HUGE_VAL;
        int i;

        for (i = 0; i < bd.L; i++) {
                if (gsl_vector_get(bd.beta, i) > 0) {
                        sum = logadd(sum, ev_log[i]);
                }
        }

        return sum;
}

static
prob_t mbeta_log(prob_t *p)
{
        unsigned int i;
        prob_t sum1, sum2;

        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < bd.events; i++) {
                sum1 += p[i];
                sum2 += gsl_sf_lngamma(p[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

static /* P(E|B) */
prob_t iec_log(binProblem *bp, int kk, int k)
{
        unsigned int i;
        prob_t c[bd.events];
        prob_t alpha[bd.events];
        prob_t gamma = gsl_matrix_get(bd.gamma, kk, k);
        if (gamma == 0) {
                return -HUGE_VAL;
        }
        for (i = 0; i < bd.events; i++) {
                c[i]     = countStatistic(bp, i, kk, k) + countAlpha(i, kk, k);
                alpha[i] = countAlpha(i, kk, k);
        }
        if (bp != NULL && kk <= bp->fix_prob.pos && bp->fix_prob.pos <= k) {
                // compute marginals
                // TODO: extend to multinomial case
                if (bp->fix_prob.which == 0) {
                        return logl(gamma) + (c[0]-1)*logl(bp->fix_prob.val)
                                + (c[1]-1)*logl(1-bp->fix_prob.val)
                                - mbeta_log(alpha);
                }
                else {
                        return logl(gamma) + (c[0]-1)*log(1-bp->fix_prob.val)
                                + (c[1]-1)*log(bp->fix_prob.val)
                                - mbeta_log(alpha);
                }
        }
        else {
                return logl(gamma) + (mbeta_log(c) - mbeta_log(alpha));
        }
}

/* Find the smallest m_B for which the prior P(m_B) is nonzero. */
static
int minM()
{
        int i;
        for (i = bd.L-1; i > 0; i--) {
                if (gsl_vector_get(bd.beta, i) > 0) {
                        return i;
                }
        }
        return i;
}

static
void binProblemInit(binProblem *bp)
{
        if (bd.options->algorithm == 0) {
                bp->ak      = allocMatrix(bd.L, bd.L);
        }
        else {
                bp->ak      = NULL;
        }
        bp->bprob_pos       = -1;
        bp->counts_pos      = -1;
        bp->add_event.pos   = -1;
        bp->add_event.n     = 0;
        bp->add_event.which = bd.options->which;
        bp->fix_prob.pos    = -1;
        bp->fix_prob.val    = 0;
        bp->fix_prob.which  = bd.options->which;
}

static
void binProblemFree(binProblem *bp)
{
        if (bp->ak) {
                freeMatrix(bp->ak);
        }
}

////////////////////////////////////////////////////////////////////////////////
// Prombs
////////////////////////////////////////////////////////////////////////////////

static
void callBinningAlgorithm(
        binProblem *bp,
        prob_t (*f)(int, int, void*),
        prob_t *ev_log)
{
        switch (bd.options->algorithm) {
        default:
        case 0:
                prombs(ev_log, bp->ak, bd.prior_log, f, bd.L, minM(), (void *)bp);
                break;
        case 1:
                prombs_tree(ev_log, bd.prior_log, f, bd.L, minM(), (void *)bp);
                break;
        case 2:
                mgs(ev_log, bd.prior_log, f, bd.L, (void *)bp);
                break;
        }
}

static
prob_t execPrombs_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(bp, i, j);
}

static
prob_t evidence(binProblem *bp, prob_t *ev_log)
{
        callBinningAlgorithm(bp, execPrombs_f, ev_log);

        return sumModels(ev_log);
}

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
prob_t differentialEntropy(binProblem *bp, int n, int i, int j, prob_t evidence_ref)
{
        prob_t ev_log[bd.L];
        prob_t sum;

        bp->add_event.n     = n;
        bp->add_event.pos   = i;
        bp->add_event.which = j;

        prombsExt(ev_log, bp->ak, bd.prior_log, &differentialEntropy_f, &differentialEntropy_h, bd.L, minM(), (void *)bp);

        sum = sumModels(ev_log);
        if (sum == -HUGE_VAL) {
                return 0.0;
        }
        else {
                return -expl(sumModels(ev_log) - evidence_ref);
        }
}

static
prob_t effectiveCounts_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        if (i <= bp->counts_pos && bp->counts_pos <= j) {
                int k;
                prob_t n = 0;
                for (k = 0; k < bd.events; k++) {
                        n += countStatistic(NULL, k, i, j) + countAlpha(k, i, j);
                }
                return log(n) + iec_log(NULL, i, j);
        }
        else {
                return iec_log(NULL, i, j);
        }
}
static
prob_t effectiveCounts(binProblem *bp, unsigned int pos, prob_t evidence_ref)
{
        prob_t ev_log[bd.L];

        bp->counts_pos = pos;
        prombs(ev_log, bp->ak, bd.prior_log, &effectiveCounts_f, bd.L, minM(), (void *)bp);

        return expl(sumModels(ev_log) - evidence_ref);
}

static
prob_t breakProb_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        // include only those bins that don't cover
        // position pos, which means, that only multi-bins
        // are included, that have a break at position
        // pos
        if (i < bp->bprob_pos && bp->bprob_pos < j) {
                return -HUGE_VAL;
        }
        return iec_log(NULL, i, j);
}
static
prob_t breakProb(binProblem *bp, unsigned int pos, prob_t evidence_ref)
{
        prob_t ev_log[bd.L];

        bp->bprob_pos = pos;
        callBinningAlgorithm(bp, &breakProb_f, ev_log);

        return expl(sumModels(ev_log) - evidence_ref);
}

////////////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////////////

static
prob_t moment(
        binProblem *bp,
        unsigned int nth,
        unsigned int pos,
        prob_t evidence_ref)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bd.L];

        bp->add_event.pos = pos;
        bp->add_event.n   = nth;
        evidence_log      = evidence(bp, evidence_log_tmp);
        bp->add_event.pos = -1;
        bp->add_event.n   = 0;

        return expl(evidence_log - evidence_ref);
}

static
prob_t marginal(
        binProblem *bp,
        int pos,
        prob_t val,
        int which,
        prob_t evidence_ref)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bd.L];

        bp->fix_prob.pos   = pos;
        bp->fix_prob.val   = val;
        bp->fix_prob.which = which;
        evidence_log       = evidence(bp, evidence_log_tmp);
        bp->fix_prob.pos   = -1;
        bp->fix_prob.val   =  0;
        bp->fix_prob.which =  0;

        return expl(evidence_log - evidence_ref);
}

////////////////////////////////////////////////////////////////////////////////
// Loop through all X
////////////////////////////////////////////////////////////////////////////////

static
void computeEffectiveCounts(
        prob_t *result,
        prob_t evidence_ref)
{
        binProblem bp; binProblemInit(&bp);
        unsigned int i;

        for (i = 0; i < bd.L; i++) {
                notice(NONE, "Computing effective counts... %.1f%%", (float)100*(i+1)/bd.L);
                result[i] = effectiveCounts(&bp, i, evidence_ref);
        }

        binProblemFree(&bp);
}

static
void computeModelPosteriors(
        prob_t *ev_log,
        prob_t *mpost,
        prob_t evidence_ref)
{
        unsigned int j;

        if (bd.options->algorithm == 2) {
                size_t *counts = mgs_get_counts();
                for (j = 0; j < bd.L; j++) {
                        if (gsl_vector_get(bd.beta, j) > 0) {
                                mpost[j] = (prob_t)counts[j]/bd.options->samples[1];
                        }
                        else {
                                mpost[j] = 0;
                        }
                }
        }
        else {
                for (j = 0; j < bd.L; j++) {
                        if (gsl_vector_get(bd.beta, j) > 0) {
                                mpost[j] = expl(ev_log[j] - evidence_ref);
                        }
                        else {
                                mpost[j] = 0;
                        }
                }
        }
}

typedef struct {
        binProblem *bp;
        int i;
        prob_t *bprob;
        prob_t evidence_ref;
} pthread_data_bprob;

static
void * computeBreakProbabilities_thread(void* data_)
{
        pthread_data_bprob *data  = (pthread_data_bprob *)data_;
        binProblem *bp = data->bp;
        int i = data->i;
        prob_t *bprob = data->bprob;
        prob_t evidence_ref = data->evidence_ref;

        bprob[i] = breakProb(bp, i, evidence_ref);
        return NULL;
}

static
void computeBreakProbabilities(
        prob_t *bprob,
        prob_t evidence_ref)
{
        unsigned int i, j, rc;

        binProblem bp[bd.options->threads];
        pthread_t threads[bd.options->threads];
        pthread_data_bprob data[bd.options->threads];
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        if (bd.options->stacksize < PTHREAD_STACK_MIN) {
                if (pthread_attr_setstacksize (&attr, PTHREAD_STACK_MIN) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }
        else {
                if (pthread_attr_setstacksize (&attr, (size_t)bd.options->stacksize) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }

        for (j = 0; j < bd.options->threads && j < bd.L; j++) {
                binProblemInit(&bp[j]);
                data[j].bp = &bp[j];
                data[j].bprob = bprob;
                data[j].evidence_ref = evidence_ref;
        }
        for (i = 0; i < bd.L; i += bd.options->threads) {
                for (j = 0; j < bd.options->threads && i+j < bd.L; j++) {
                        notice(NONE, "Computing break probabilities: %.1f%%", (float)100*(i+j+1)/bd.L);
                        data[j].i = i+j;
                        rc = pthread_create(&threads[j], &attr, computeBreakProbabilities_thread, (void *)&data[j]);
                        if (rc) {
                                std_err(NONE, "Couldn't create thread.");
                        }
                }
                for (j = 0; j < bd.options->threads && i+j < bd.L; j++) {
                        rc = pthread_join(threads[j], NULL);
                        if (rc) {
                                std_err(NONE, "Couldn't join thread.");
                        }
                }
        }
        for (j = 0; j < bd.options->threads && j < bd.L; j++) {
                binProblemFree(&bp[j]);
        }
}

static
void computeDifferentialUtility(
        prob_t *result,
        prob_t evidence_ref)
{
        binProblem bp; binProblemInit(&bp);
        unsigned int i, j;
        prob_t expected_entropy;
        prob_t entropy;
        prob_t evidence_log;
        prob_t evidence_log_tmp[bd.L];

        entropy = differentialEntropy(&bp, 0, -1, bd.options->which, evidence_ref);

        for (i = 0; i < bd.L; i++) {
                notice(NONE, "Computing utilities... %.1f%%", (float)100*(i+1)/bd.L);
                // recompute the evidence
                result[i] = 0;
                for (j = 0; j < bd.events; j++) {
                        evidence_log = evidence(&bp, evidence_log_tmp);
                        // expected entropy for event j
                        expected_entropy = differentialEntropy(&bp, 1, i, j, evidence_log);
                        // initialize sum
                        result[i] += expl(evidence_log - evidence_ref)*expected_entropy;
                }
                result[i] = entropy - result[i];
        }

        binProblemFree(&bp);
}

typedef struct {
        binProblem *bp;
        int i;
        prob_t **moments;
        prob_t evidence_ref;
} pthread_data_moments;

static
void * computeMoments_thread(void* data_)
{
        pthread_data_moments *data  = (pthread_data_moments *)data_;
        binProblem *bp = data->bp;
        int i = data->i, j;
        prob_t **moments = data->moments;
        prob_t evidence_ref = data->evidence_ref;

        // Moments
        for (j = 0; j < bd.options->n_moments; j++) {
                moments[j][i] = moment(bp, j+1, i, evidence_ref);
        }
        return NULL;
}

static
void computeMoments(
        prob_t **moments,
        prob_t evidence_ref)
{
        int i, j, rc;

        binProblem bp[bd.options->threads];
        pthread_t threads[bd.options->threads];
        pthread_data_moments data[bd.options->threads];
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        if (bd.options->stacksize < PTHREAD_STACK_MIN) {
                if (pthread_attr_setstacksize (&attr, PTHREAD_STACK_MIN) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }
        else {
                if (pthread_attr_setstacksize (&attr, (size_t)bd.options->stacksize) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }

        for (j = 0; j < bd.options->threads && j < bd.L; j++) {
                binProblemInit(&bp[j]);
                data[j].bp = &bp[j];
                data[j].moments = moments;
                data[j].evidence_ref = evidence_ref;
        }
        for (i = 0; i < bd.L; i += bd.options->threads) {
                for (j = 0; j < bd.options->threads && i+j < bd.L; j++) {
                        notice(NONE, "Computing moments... %.1f%%", (float)100*(i+j+1)/bd.L);
                        data[j].i = i+j;
                        rc = pthread_create(&threads[j], &attr, computeMoments_thread, (void *)&data[j]);
                        if (rc) {
                                std_err(NONE, "Couldn't create thread.");
                        }
                }
                for (j = 0; j < bd.options->threads && i+j < bd.L; j++) {
                        rc = pthread_join(threads[j], NULL);
                        if (rc) {
                                std_err(NONE, "Couldn't join thread.");
                        }
                }
        }
        for (j = 0; j < bd.options->threads && j < bd.L; j++) {
                binProblemFree(&bp[j]);
        }
}

static
void computeMarginal(
        prob_t **marginals,
        prob_t evidence_ref)
{
        binProblem bp; binProblemInit(&bp);
        int i, j;

        for (i = 0; i < bd.L; i++) {
                notice(NONE, "Computing marginals... %.1f%%", (float)100*(i+1)/bd.L);
                marginals[i][0] = 0;
                for (j = 1; j < bd.options->n_marginals; j++) {
                        prob_t p = j*bd.options->marginal_step;
                        if (bd.options->marginal_range.from <= p &&
                            bd.options->marginal_range.to   >= p) {
                                marginals[i][j] = marginal(&bp, i, p, bd.options->which, evidence_ref);
                        }
                        else {
                                marginals[i][j] = 0;
                        }
                }
        }

        binProblemFree(&bp);
}

////////////////////////////////////////////////////////////////////////////////
// Simple prombs test
////////////////////////////////////////////////////////////////////////////////

static
prob_t prombsTest_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(bp, i, j);
}
static
prob_t prombsTest_h(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        // - Log[f(b)]
        return -iec_log(bp, i, j);
}
static
void prombsTest()
{
        MET_INIT;
        prob_t result1[bd.L];
        prob_t result2[bd.L];
        prob_t sum;
        unsigned int i;
        Matrix *ak = allocMatrix(bd.L, bd.L);

        // set prior to 1
        for (i = 0; i < bd.L; i++) {
                bd.prior_log[i] = 0;
        }
        MET("Testing prombs",
            prombs   (result1, ak, bd.prior_log, &prombsTest_f, bd.L, minM(), NULL));
        MET("Testing prombsExt",
            prombsExt(result2, ak, bd.prior_log, &prombsTest_f, &prombsTest_h, bd.L, minM(), NULL));

        sum = -HUGE_VAL;
        for (i = 0; i < bd.L; i++) {
                (void)printf("prombs[%02d]: %.10f\n", i, (double)result1[i]);
                sum = logadd(sum, result1[i]);
        }
        (void)printf("prombs: %.10f\n", (double)sum);

        sum = -HUGE_VAL;
        for (i = 0; i < bd.L; i++) {
                (void)printf("prombsExt[%02d]: %.10f\n", i, (double)result2[i]);
                sum = logadd(sum, result2[i]);
        }
        (void)printf("prombsExt: %.10f\n", (double)sum);

        freeMatrix(ak);
}

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
void computeBinning(
        prob_t **moments,
        prob_t **marginals,
        prob_t  *bprob,
        prob_t  *mpost,
        prob_t  *differential_gain,
        prob_t  *effective_counts)
{
        binProblem bp; binProblemInit(&bp);
        prob_t evidence_ref;
        prob_t evidence_log_tmp[bd.L];

        // compute the model prior once for all computations
        computeModelPrior();
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
        // compute the differential entropy
        if (bd.options->differential_gain) {
                computeDifferentialUtility(differential_gain, evidence_ref);
        }
        // compute effective counts
        if (bd.options->effective_counts) {
                computeEffectiveCounts(effective_counts, evidence_ref);
        }

        if (bd.options->algorithm == 2) {
                mgs_free();
        }
        binProblemFree(&bp);
}

////////////////////////////////////////////////////////////////////////////////
// Library entry point, initializes data structures
////////////////////////////////////////////////////////////////////////////////

void __init_rand__() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;

        srand(seed);
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
        size_t K = counts[0]->size2;
        BinningResultGSL *result = (BinningResultGSL *)malloc(sizeof(BinningResultGSL));
        result->moments   = (options->n_moments ?
                             gsl_matrix_alloc(options->n_moments, K)   : NULL);
        result->marginals = (options->marginal  ?
                             gsl_matrix_alloc(K, options->n_marginals) : NULL);
        result->bprob             = gsl_vector_alloc(K);
        result->mpost             = gsl_vector_alloc(K);
        result->differential_gain = gsl_vector_alloc(K);
        result->effective_counts  = gsl_vector_alloc(K);
        prob_t * moments[options->n_moments];
        prob_t * marginals[K];
        prob_t bprob[K];
        prob_t mpost[K];
        prob_t differential_gain[K];
        prob_t effective_counts[K];
        prob_t prior_log[K];
        unsigned int i, j;

        __init_rand__();

        for (i = 0; i < options->n_moments; i++) {
                moments[i]   = (prob_t *)calloc(K, sizeof(prob_t));
        }
        for (i = 0; i < K; i++) {
                marginals[i] = (prob_t *)calloc(options->n_marginals, sizeof(prob_t));
        }
        bzero(bprob,             K*sizeof(prob_t));
        bzero(mpost,             K*sizeof(prob_t));
        bzero(prior_log,         K*sizeof(prob_t));
        bzero(differential_gain, K*sizeof(prob_t));
        bzero(effective_counts,  K*sizeof(prob_t));

        prombs_init(options->epsilon);

        verbose            = options->verbose;
        bd.options         = options;
        bd.beta            = beta;
        bd.L               = K;
        bd.prior_log       = prior_log;
        bd.events          = events;
        bd.counts          = counts;
        bd.alpha           = alpha;
        bd.beta            = beta;
        bd.gamma           = gamma;
        if (options->prombsTest) {
                prombsTest();
        }
        computeBinning(moments, marginals, bprob, mpost, differential_gain,
                       effective_counts);
        for (i = 0; i <= bd.L-1; i++) {
                for (j = 0; j < options->n_moments; j++) {
                        gsl_matrix_set(result->moments, j, i, moments[j][i]);
                }
                if (options->marginal) {
                        for (j = 0; j < options->n_marginals; j++) {
                                gsl_matrix_set(result->marginals, i, j, marginals[i][j]);
                        }
                }
                gsl_vector_set(result->bprob, i, bprob[i]);
                gsl_vector_set(result->mpost, i, mpost[i]);
                gsl_vector_set(result->differential_gain, i, differential_gain[i]);
                gsl_vector_set(result->effective_counts,  i, effective_counts[i]);
        }

        for (i = 0; i < options->n_moments; i++) {
                free(moments[i]);
        }
        for (i = 0; i < K; i++) {
                free(marginals[i]);
        }

        return result;
}
