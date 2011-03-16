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
#include <gsl/gsl_sf_psi.h>

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
        gsl_matrix *alpha;
        // internal data
        gsl_matrix **counts_m;
        gsl_matrix **alpha_m;
        struct {
                int pos;
                int n;
                int which;
        } add_event;
        struct {
                int pos;
                prob_t val;
                int which;
        } fix_prob;
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
prob_t countAlpha(binProblem *bp, unsigned int event, int ks, int ke)
{
        if (ks <= ke) {
                return gsl_matrix_get(bp->alpha_m[event], ks, ke);
        }
        else {
                return 0;
        }
}

static
prob_t mbeta_log(binProblem *bp, prob_t *p)
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
void computeCountStatistics(binProblem *bp)
{
        int ks, ke, i, j;
        unsigned int c[bp->events];
        prob_t       d[bp->events];

        for (ks = 0; ks < bp->T; ks++) {
                for (ke = ks; ke < bp->T; ke++) {
                        // set all counts to zero
                        for (j = 0; j < bp->events; j++) {
                                c[j] = 0;
                                d[j] = 0;
                        }
                        // count
                        for (i = ks; i <= ke; i++) {
                                for (j = 0; j < bp->events; j++) {
                                        c[j] += gsl_matrix_get(bp->counts, j, i);
                                        d[j] += gsl_matrix_get(bp->alpha,  j, i)/(ke - ks + 1);
                                }
                        }
                        // save result
                        for (j = 0; j < bp->events; j++) {
                                gsl_matrix_set(bp->counts_m[j], ks, ke, c[j]);
                                gsl_matrix_set(bp->alpha_m[j],  ks, ke, d[j]);
                        }
                }
        }
}

static
void computePrior_log(binProblem *bp)
{
        unsigned int m_b;
        for (m_b = 0; m_b < bp->T; m_b++) {
                bp->prior_log[m_b] = -gsl_sf_lnchoose(bp->T-1, m_b);
        }
}

static /* P(E|B) */
prob_t iec_log(binProblem *bp, int kk, int k)
{
        unsigned int i;
        prob_t c[bp->events];
        prob_t alpha[bp->events];
        for (i = 0; i < bp->events; i++) {
                c[i]     = countStatistic(bp, i, kk, k) + countAlpha(bp, i, kk, k);
                alpha[i] = countAlpha(bp, i, kk, k);
        }
        if (kk <= bp->fix_prob.pos && bp->fix_prob.pos <= k) {
                // compute marginals
                // TODO: extend to multinomial case
                if (bp->fix_prob.which == 0) {
                        return (c[0]-1)*log(bp->fix_prob.val)
                                + (c[1]-1)*log(1-bp->fix_prob.val)
                                - mbeta_log(bp, alpha);
                }
                else {
                        return (c[0]-1)*log(1-bp->fix_prob.val)
                                + (c[1]-1)*log(bp->fix_prob.val)
                                - mbeta_log(bp, alpha);
                }
        }
        else {
                return mbeta_log(bp, c) - mbeta_log(bp, alpha);
        }
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
prob_t computeEvidence(binProblem *bp, prob_t *ev_log)
{
        prob_t sum;
        unsigned int j;

        execPrombs(bp, ev_log, -1);
        sum = -HUGE_VAL;
        for (j = 0; j < bp->T; j++) {
                if (gsl_vector_get(bp->mprior, j) > 0) {
                        prob_t mprior = gsl_vector_get(bp->mprior, j);
                        sum = logadd(sum, ev_log[j] + logl(mprior));
                }
        }
        return sum;
}

static
prob_t computeMoment(
        binProblem *bp,
        unsigned int nth,
        unsigned int pos,
        prob_t evidence_ref)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bp->T];

        bp->add_event.pos = pos;
        bp->add_event.n   = nth;
        evidence_log      = computeEvidence(bp, evidence_log_tmp);
        bp->add_event.pos = -1;
        bp->add_event.n   = 0;

        return expl(evidence_log - evidence_ref);
}

static
prob_t computeMarginal(
        binProblem *bp,
        int pos,
        prob_t val,
        int which,
        prob_t evidence_ref)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bp->T];

        bp->fix_prob.pos   = pos;
        bp->fix_prob.val   = val;
        bp->fix_prob.which = which;
        evidence_log       = computeEvidence(bp, evidence_log_tmp);
        bp->fix_prob.pos   = -1;
        bp->fix_prob.val   =  0;
        bp->fix_prob.which =  0;

        return expl(evidence_log - evidence_ref);
}

static
prob_t singlebinEntropy(binProblem *bp, int i, int j)
{
        unsigned int k;
        prob_t n = 0;
        prob_t c[bp->events];
        prob_t sum = 0;
        for (k = 0; k < bp->events; k++) {
                c[k] = countStatistic(bp, k, i, j) + countAlpha(bp, k, i, j);
                sum += (c[k] - 1.0)*gsl_sf_psi(c[k]);
                n   +=  c[k];
        }
        return mbeta_log(bp, c) + (n - bp->events)*gsl_sf_psi(n) - sum;
}

static
prob_t differentialEntropy(binProblem *bp, prob_t evidence)
{
        unsigned int i;
        prob_t ev_log[bp->T];
        prob_t sum;
        prob_t epsilon = bp->epsilon;
        prob_t f(int i, int j)
        {
                return iec_log(bp, i, j);
        }
        prob_t h(int i, int j)
        {
                return -singlebinEntropy(bp, i, j);
        }

        prombsExt(ev_log, bp->prior_log, &f, &h, epsilon, bp->T, bp->T-1);
        sum = -HUGE_VAL;
        for (i = 0; i < bp->T; i++) {
                if (gsl_vector_get(bp->mprior, i) > 0) {
                        prob_t mprior = gsl_vector_get(bp->mprior, i);
                        if (!isfinite(ev_log[i])) {
                                return 0;
                        }
                        sum = logadd(sum, ev_log[i] + logl(mprior));
                }
        }
        return -expl(sum - evidence);
}

static
void differentialUtility(binProblem *bp, prob_t *result, prob_t evidence_ref, Options *options)
{
        unsigned int i, j;
        prob_t expected_entropy;
        prob_t entropy;
        prob_t evidence;
        prob_t evidence_log_tmp[bp->T];

        computePrior_log(bp);
        entropy = differentialEntropy(bp, evidence_ref);

        bp->add_event.n = 1;
        for (i = 0; i < bp->T; i++) {
                notice(NONE, "Computing utilities... %.1f%%", (float)100*(i+1)/bp->T);
                bp->add_event.pos = i;
                // recompute the evidence
                result[i] = 0;
                for (j = 0; j < bp->events; j++) {
                        bp->add_event.which = j;
                        evidence = computeEvidence(bp, evidence_log_tmp);
                        // expected entropy for event j
                        expected_entropy = differentialEntropy(bp, evidence);
                        // initialize sum
                        result[i] += expl(evidence - evidence_ref)*expected_entropy;
                }
                result[i] = entropy - result[i];
        }
        bp->add_event.n     =  0;
        bp->add_event.pos   = -1;
        bp->add_event.which = options->which;
}

static
prob_t effectiveCounts(binProblem *bp, unsigned int pos, prob_t evidence)
{
        unsigned int i;
        prob_t ev_log[bp->T];
        prob_t sum;
        prob_t f(int i, int j)
        {
                if (i <= pos && pos <= j) {
                        int k;
                        prob_t n = 0;
                        for (k = 0; k < bp->events; k++) {
                                n += countStatistic(bp, k, i, j) + countAlpha(bp, k, i, j);
                        }
                        return log(n) + iec_log(bp, i, j);
                }
                else {
                        return iec_log(bp, i, j);
                }
        }

        prombs(ev_log, bp->prior_log, &f, bp->T, bp->T-1);
        sum = -HUGE_VAL;
        for (i = 0; i < bp->T; i++) {
                if (gsl_vector_get(bp->mprior, i) > 0) {
                        prob_t mprior = gsl_vector_get(bp->mprior, i);
                        if (!isfinite(ev_log[i])) {
                                return 0;
                        }
                        sum = logadd(sum, ev_log[i] + logl(mprior));
                }
        }
        return expl(sum - evidence);
}

static
void computeEffectiveCounts(binProblem *bp, prob_t *result, prob_t evidence_ref, Options *options)
{
        unsigned int i;

        computePrior_log(bp);
        for (i = 0; i < bp->T; i++) {
                notice(NONE, "Computing effective counts... %.1f%%", (float)100*(i+1)/bp->T);
                result[i] = effectiveCounts(bp, i, evidence_ref);
        }
}

static
void computeMultibinEntropy(binProblem *bp, prob_t *result, prob_t evidence)
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
void computeModelPosteriors(binProblem *bp, prob_t *ev_log, prob_t *mpost, prob_t P_D)
{
        unsigned int j;

        for (j = 0; j < bp->T; j++) {
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
                notice(NONE, "Computing break probabilities: %.1f%%", (float)100*i/bp->T);
                execPrombs(bp, ev_log, i);
                sum = -HUGE_VAL;
                for (j = 0; j < bp->T; j++) {
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
void computeBinning(
        binProblem *bp,
        prob_t **moments,
        prob_t **marginals,
        prob_t *bprob,
        prob_t *mpost,
        prob_t *differential_gain,
        prob_t *effective_counts,
        prob_t *multibin_entropy,
        Options *options)
{
        prob_t evidence_ref;
        prob_t evidence_log_tmp[bp->T];
        unsigned int i, j;

        // compute evidence P(D)
        computePrior_log(bp);
        evidence_ref = computeEvidence(bp, evidence_log_tmp);
        // compute model posteriors P(m_B|D)
        if (options->model_posterior) {
                computeModelPosteriors(bp, evidence_log_tmp, mpost, evidence_ref);
        }
        if (options->marginal) {
                for (i = 0; i < bp->T; i++) {
                        notice(NONE, "Computing marginals... %.1f%%", (float)100*(i+1)/bp->T);
                        marginals[i][0] = 0;
                        for (j = 1; j < options->n_marginals; j++) {
                                prob_t p = j*options->marginal_step;
                                if (options->marginal_range.from <= p &&
                                    options->marginal_range.to   >= p) {
                                        marginals[i][j] = computeMarginal(bp, i, p, options->which, evidence_ref);
                                }
                                else {
                                        marginals[i][j] = 0;
                                }
                        }
                }
        }
        // break probability
        if (options->bprob) {
                computeBreakProbabilities(bp, bprob, evidence_ref);
        }
        // compute the multibin entropy
        if (options->multibin_entropy) {
                computeMultibinEntropy(bp, multibin_entropy, evidence_ref);
        }
        // for each timestep compute the first n moments
        if (options->n_moments > 0) {
                for (i = 0; i < bp->T; i++) {
                        notice(NONE, "Computing moments... %.1f%%", (float)100*(i+1)/bp->T);
                        // Moments
                        for (j = 0; j < options->n_moments; j++) {
                                moments[j][i] = computeMoment(bp, j+1, i, evidence_ref);
                        }
                }
        }
        // compute the differential entropy
        if (options->differential_gain) {
                differentialUtility(bp, differential_gain, evidence_ref, options);
        }
        // compute effective counts
        if (options->effective_counts) {
                computeEffectiveCounts(bp, effective_counts, evidence_ref, options);
        }
}

BinningResultGSL *
bin_log(
        gsl_matrix *counts,
        gsl_matrix *alpha,
        gsl_vector *mprior,
        Options *options)
{
        size_t K = counts->size2;
        BinningResultGSL *result = (BinningResultGSL *)malloc(sizeof(BinningResultGSL));
        if (options->n_moments > 0) {
                result->moments = gsl_matrix_alloc(options->n_moments, K);
        }
        else {
                result->moments = NULL;
        }
        if (options->marginal) {
                result->marginals = gsl_matrix_alloc(K, options->n_marginals);
        }
        else {
                result->marginals = NULL;
        }
        result->bprob   = gsl_vector_alloc(K);
        result->mpost   = gsl_vector_alloc(K);
        result->differential_gain = gsl_vector_alloc(K);
        result->effective_counts  = gsl_vector_alloc(K);
        result->multibin_entropy  = gsl_vector_alloc(K);
        prob_t * moments[options->n_moments];
        prob_t * marginals[K];
        prob_t bprob[K];
        prob_t mpost[K];
        prob_t differential_gain[K];
        prob_t effective_counts[K];
        prob_t multibin_entropy[K];
        prob_t prior_log[K];
        unsigned int i, j;
        binProblem bp;

        for (i = 0; i < options->n_moments; i++) {
                moments[i] = (prob_t *)calloc(K, sizeof(prob_t));
        }
        for (i = 0; i < K; i++) {
                marginals[i] = (prob_t *)calloc(options->n_marginals, sizeof(prob_t));
        }
        bzero(bprob,     K*sizeof(prob_t));
        bzero(mpost,     K*sizeof(prob_t));
        bzero(prior_log, K*sizeof(prob_t));
        bzero(differential_gain, K*sizeof(prob_t));
        bzero(effective_counts,  K*sizeof(prob_t));
        bzero(multibin_entropy,  K*sizeof(prob_t));

        verbose       = options->verbose;

        bp.epsilon    = options->epsilon;
        bp.mprior     = mprior;
        bp.T          = K;
        bp.prior_log  = prior_log;
        bp.add_event.pos   = -1;
        bp.add_event.n     = 0;
        bp.add_event.which = options->which;
        bp.fix_prob.pos    = -1;
        bp.fix_prob.val    = 0;
        bp.fix_prob.which  = options->which;
        bp.events          = counts->size1;
        bp.counts          = counts;
        bp.alpha           = alpha;
        bp.counts_m        = (gsl_matrix  **)malloc(bp.events*sizeof(gsl_matrix *));
        bp.alpha_m         = (gsl_matrix  **)malloc(bp.events*sizeof(gsl_matrix *));
        for (i = 0; i < bp.events; i++) {
                bp.counts_m[i] = gsl_matrix_alloc(K, K);
                bp.alpha_m[i]  = gsl_matrix_alloc(K, K);
        }
        computeCountStatistics(&bp);
        if (options->prombsTest) {
                prombsTest(&bp);
        }
        computeBinning(&bp, moments, marginals, bprob, mpost, differential_gain,
                       effective_counts, multibin_entropy, options);
        for (i = 0; i <= bp.T-1; i++) {
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
                gsl_vector_set(result->multibin_entropy,  i, multibin_entropy[i]);
        }

        for (i = 0; i < bp.events; i++) {
                gsl_matrix_free(bp.counts_m[i]);
                gsl_matrix_free(bp.alpha_m[i]);
        }
        for (i = 0; i < options->n_moments; i++) {
                free(moments[i]);
        }
        for (i = 0; i < K; i++) {
                free(marginals[i]);
        }
        free(bp.counts_m);
        free(bp.alpha_m);

        return result;
}
