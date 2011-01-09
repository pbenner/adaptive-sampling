/* Copyright (C) 2010 Philipp Benner
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

#include <bayes_exception.h>
#include <bayes_logarithmetic.h>
#include <bayes_prombs.h>
#include <bayes_datatypes.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

typedef struct {
        // precision for entropy estimates
        prob_t epsilon;
        // number of timesteps
        unsigned int T;
        unsigned int events;
        gsl_matrix  *counts;
        gsl_vector  *mprior;   // P(m_B)
        prob_t *prior_log;     // P(p,B|m_B)
        // hyperparameters
        unsigned int *alpha;
        // internal data
        gsl_matrix **counts_m;
        struct {
                int pos;
                int n;
                int which;
        } add_event;
} binProblem;

static
unsigned int countStatistic(binProblem *bp, unsigned int event, int ks, int ke)
{
        if (bp->add_event.which == event &&
            ks <= bp->add_event.pos && bp->add_event.pos <= ke) {
                return gsl_matrix_get(bp->counts_m[event], ks, ke) +
                        bp->add_event.n;
        }
        else if (ks <= ke) {
                return gsl_matrix_get(bp->counts_m[event], ks, ke);
        }
        else {
                return 0;
        }
}

static
prob_t mbeta_log(binProblem *bp, unsigned int *p)
{
        unsigned int i;
        prob_t sum1, sum2;

        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < bp->events; i++) {
                sum1 += p[i];
                sum2 += gsl_sf_lngamma(p[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
}

static
prob_t mbetaInv_log(binProblem *bp, unsigned int *p)
{
        unsigned int i;
        prob_t sum1, sum2;

        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < bp->events; i++) {
                sum1 += p[i];
                sum2 += gsl_sf_lngamma(p[i]);
        }

        return gsl_sf_lngamma(sum1) - sum2;
}

static
void computeCountStatistics(binProblem *bp)
{
        int ks, ke, i, j;
        unsigned int c[bp->events];

        for (ks = 0; ks < bp->T; ks++) {
                for (ke = ks; ke < bp->T; ke++) {
                        // set all counts to zero
                        for (j = 0; j < bp->events; j++) {
                                c[j] = 0;
                        }
                        // count
                        for (i = ks; i <= ke; i++) {
                                for (j = 0; j < bp->events; j++) {
                                        c[j] += gsl_matrix_get(bp->counts, j, i);
                                }
                        }
                        // save result
                        for (j = 0; j < bp->events; j++) {
                                gsl_matrix_set(bp->counts_m[j], ks, ke, c[j]);
                        }
                }
        }
}

static
void computePrior_log(binProblem *bp)
{
        unsigned int m;
        prob_t b;
        b = mbetaInv_log(bp, bp->alpha);
        for (m = 0; m < bp->T; m++) {
                bp->prior_log[m] = (m+1)*b - gsl_sf_lnchoose(bp->T-1, m);
        }
}

/* P(D|B,m_B)/P(p|m_B) */
static
prob_t iec_log(binProblem *bp, int kk, int k)
{
        unsigned int i;
        unsigned int c[bp->events];
        for (i = 0; i < bp->events; i++) {
                c[i] = countStatistic(bp, i, kk, k)+bp->alpha[i];
        }
        return mbeta_log(bp, c);
}

/* Find the smallest m_B for which the prior P(m_B) is nonzero. */
static
int minM(binProblem *bp)
{
        int i;
        for (i = bp->T-1; i>0; i--) {
                if (gsl_vector_get(bp->mprior, i) > 0) {
                        return i;
                }
        }
        return i;
}

static
void execPrombs(binProblem *bp, prob_t *ev_log, int pos)
{
        prob_t f(int i, int j)
        {
                // include only those bins that don't cover
                // position pos, which means, that only multi-bins
                // are included, that have a break at position
                // pos
                if (i < pos && pos < j) {
                        return -HUGE_VAL;
                }
                return iec_log(bp, i, j);
        }

        prombs(ev_log, bp->prior_log, &f, bp->T, minM(bp));
}

static
void computeEntropy(binProblem *bp, prob_t *result, prob_t evidence)
{
        size_t i;
        prob_t resulta[bp->T];
        prob_t resultb[bp->T];
        prob_t epsilona = bp->epsilon;
        prob_t epsilonb = bp->epsilon*2;
        prob_t result2[bp->T];
        prob_t g[bp->T];
        prob_t f(int i, int j)
        {
                return iec_log(bp, i, j);
        }
        prob_t h(int i, int j)
        {
                // - Log[f(b)]
                return -iec_log(bp, i, j);
        }

        for (i = 0; i < bp->T; i++) {
                g[i] = logl(bp->prior_log[i] - evidence) + bp->prior_log[i];
        }
        prombsExt(resulta, bp->prior_log, &f, &h, epsilona, bp->T, bp->T-1);
        prombsExt(resultb, bp->prior_log, &f, &h, epsilonb, bp->T, bp->T-1);
        prombs(result2, g, &f, bp->T, bp->T-1);

        for (i = 0; i < bp->T; i++) {
                prob_t pa = expl(logsub(resulta[i], result2[i]) - evidence);
                prob_t pb = expl(logsub(resultb[i], result2[i]) - evidence);
                result[i] = pa - (pb-pa)/(epsilonb-epsilona)*epsilona;
        }
}

static
prob_t computeEvidence(binProblem *bp, prob_t *ev_log)
{
        prob_t sum;
        unsigned int j;

        execPrombs(bp, ev_log, -1);
        sum = -HUGE_VAL;
        for (j=0; j<bp->T; j++) {
                if (gsl_vector_get(bp->mprior, j) > 0) {
                        prob_t mprior = gsl_vector_get(bp->mprior, j);
                        sum = logadd(sum, ev_log[j] + logl(mprior));
                }
        }
        return sum;
}

static
void computeModelPosteriors(binProblem *bp, prob_t *ev_log, prob_t *mpost, prob_t P_D)
{
        unsigned int j;

        for (j=0; j<bp->T; j++) {
                if (gsl_vector_get(bp->mprior, j) > 0) {
                        prob_t mprior = gsl_vector_get(bp->mprior, j);
                        mpost[j] = expl(ev_log[j] + logl(mprior) - P_D);
                }
                else {
                        mpost[j] = 0;
                }
        }
}

static
void computeBreakProbabilities(binProblem *bp, prob_t *bprob, prob_t P_D)
{
        prob_t ev_log[bp->T];
        prob_t sum;
        unsigned int i, j;

        for (i = 0; i < bp->T; i++) {
                notice(NONE, "break prob.: %.1f%%", (float)100*i/bp->T);
                execPrombs(bp, ev_log, i);
                sum = -HUGE_VAL;
                for (j=0; j<bp->T; j++) {
                        if (gsl_vector_get(bp->mprior, j) > 0) {
                                prob_t mprior = gsl_vector_get(bp->mprior, j);
                                sum = logadd(sum, ev_log[j] + logl(mprior));
                        }
                }
                bprob[i] = expl(sum - P_D);
        }
}

static
void prombsTest(binProblem *bp)
{
        MET_INIT;
        prob_t result1[bp->T];
        prob_t result2[bp->T];
        prob_t sum;
        unsigned int i;

        prob_t f(int i, int j)
        {
                return iec_log(bp, i, j);
        }
        prob_t h(int i, int j)
        {
                // - Log[f(b)]
                return -iec_log(bp, i, j);
        }
        // set prior to 1
        for (i = 0; i < bp->T; i++) {
                bp->prior_log[i] = 0;
        }
        MET("Testing prombs",
            prombs   (result1, bp->prior_log, &f, bp->T, bp->T-1));
        MET("Testing prombsExt",
            prombsExt(result2, bp->prior_log, &f, &h, bp->epsilon, bp->T, bp->T-1));

        sum = -HUGE_VAL;
        for (i = 0; i < bp->T; i++) {
                (void)printf("prombs[%02d]: %.10f\n", i, (double)result1[i]);
                sum = logadd(sum, result1[i]);
        }
        (void)printf("prombs: %.10f\n", (double)sum);

        sum = -HUGE_VAL;
        for (i = 0; i < bp->T; i++) {
                (void)printf("prombsExt[%02d]: %.10f\n", i, (double)result2[i]);
                sum = logadd(sum, result2[i]);
        }
        (void)printf("prombsExt: %.10f\n", (double)sum);
}

static
prob_t binomial_transform(prob_t *raw_moments, unsigned int n)
{
        unsigned int k;
        prob_t sum;

        sum = pow(-1, n)*pow(raw_moments[0], n);
        for (k = 1; k <= n; k++) {
                sum = pow(-1, n-k)*gsl_sf_choose(n, k)*raw_moments[k-1]*pow(raw_moments[0],n-k);
        }

        return sum;
}

static
void computeBinning(
        binProblem *bp,
        prob_t **moments,
        prob_t *bprob,
        prob_t *mpost,
        prob_t *entropy,
        Options *options)
{
        prob_t  ev_log[options->n_moments+1][bp->T];
        prob_t sum_log[options->n_moments+1];
        prob_t raw_moments[options->n_moments];
        unsigned int i, j, k;

        // compute evidence P(D)
        computePrior_log(bp);
        sum_log[0] = computeEvidence(bp, ev_log[0]);
        // compute model posteriors P(m_B|D)
        computeModelPosteriors(bp, ev_log[0], mpost, sum_log[0]);
        // break probability
        if (options->bprob) {
                computeBreakProbabilities(bp, bprob, sum_log[0]);
        }
        // compute the multibin entropy
        if (options->entropy) {
                computeEntropy(bp, entropy, sum_log[0]);
        }
        // for each timestep compute the first three moments
        for (i=0; i<bp->T; i++) {
                notice(NONE, "Computing moments... %.1f%%", (float)100*(i+1)/bp->T);
                bp->add_event.pos = i;
                // recompute the evidence 
                for (k = 1; k <= options->n_moments; k++) {
                        bp->add_event.n   = k;
                        computePrior_log(bp);
                        execPrombs(bp, ev_log[k], -1);
                        // initialize sum
                        sum_log[k] = -HUGE_VAL;
                }

                for (j=0; j<bp->T; j++) {
                        if (gsl_vector_get(bp->mprior, j) > 0) {
                                prob_t mprior = gsl_vector_get(bp->mprior, j);
                                for (k = 1; k <= options->n_moments; k++) {
                                        sum_log[k] = logadd(sum_log[k], ev_log[k][j] + logl(mprior));
                                }
                        }
                }
                // Moments
                for (k = 0; k < options->n_moments; k++) {
                        raw_moments[k] = expl(sum_log[k+1] - sum_log[0]);
                }
                moments[0][i] = raw_moments[0];
                for (k = 1; k < options->n_moments; k++) {
                        moments[k][i] = binomial_transform(raw_moments, k+1);
                }
/*                 // Expectation */
/*                 moments[0][i] = raw_moments[0]; */
/*                 // Variance */
/*                 moments[1][i] = raw_moments[1] - raw_moments[0]*raw_moments[0]; */
/*                 // Skewness */
/*                 moments[2][i] = raw_moments[2] - 3*raw_moments[0]*raw_moments[1] + 2*raw_moments[0]*raw_moments[0]*raw_moments[0]; */
        }
}

BinningResultGSL *
bin_log(
        gsl_matrix *counts,
        gsl_vector *alpha,
        gsl_vector *mprior,
        Options *options)
{
        size_t K = counts->size2;
        BinningResultGSL *result = (BinningResultGSL *)malloc(sizeof(BinningResultGSL));
        result->moments = gsl_matrix_alloc(options->n_moments, K);
        result->bprob   = gsl_vector_alloc(K);
        result->mpost   = gsl_vector_alloc(K);
        result->entropy = gsl_vector_alloc(K);
        prob_t * moments[options->n_moments];
        prob_t bprob[K];
        prob_t mpost[K];
        prob_t entropy[K];
        prob_t prior_log[K];
        unsigned int i, j;
        binProblem bp;

        for (i = 0; i < options->n_moments; i++) {
                moments[i] = (prob_t *)calloc(K, sizeof(prob_t));
        }
        bzero(bprob,     K*sizeof(prob_t));
        bzero(mpost,     K*sizeof(prob_t));
        bzero(entropy,   K*sizeof(prob_t));
        bzero(prior_log, K*sizeof(prob_t));

        verbose       = options->verbose;

        bp.epsilon    = options->epsilon;
        bp.mprior     = mprior;
        bp.T          = K;
        bp.prior_log  = prior_log;
        bp.add_event.pos   = 0;
        bp.add_event.n     = 0;
        bp.add_event.which = options->which;
        bp.events     = counts->size1;
        bp.counts     = counts;
        bp.alpha      = (unsigned int *)malloc(bp.events*sizeof(unsigned int));
        bp.counts_m   = (gsl_matrix  **)malloc(bp.events*sizeof(gsl_matrix *));
        for (i = 0; i < bp.events; i++) {
                bp.counts_m[i] = gsl_matrix_alloc(K, K);
                bp.alpha[i]    = gsl_vector_get(alpha, i);
        }
        computeCountStatistics(&bp);
        if (options->prombsTest) {
                prombsTest(&bp);
        }
        computeBinning(&bp, moments, bprob, mpost, entropy, options);
        for (i = 0; i <= bp.T-1; i++) {
                for (j = 0; j < options->n_moments; j++) {
                        gsl_matrix_set(result->moments, j, i, moments[j][i]);
                }
                gsl_vector_set(result->bprob,   i, bprob[i]);
                gsl_vector_set(result->mpost,   i, mpost[i]);
                gsl_vector_set(result->entropy, i, entropy[i]);
        }

        for (i = 0; i < bp.events; i++) {
                gsl_matrix_free(bp.counts_m[i]);
        }
        for (i = 0; i < options->n_moments; i++) {
                free(moments[i]);
        }
        free(bp.alpha);
        free(bp.counts_m);

        return result;
}
