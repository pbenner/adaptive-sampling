/* Copyright (C) 2012 Philipp Benner
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

#include <gsl/gsl_sf_psi.h>

#include <datatypes.h>
#include <model.h>
#include <utility.h>
#include <threading.h>
#include <tools.h>

/******************************************************************************
 * Tools
 ******************************************************************************/

static
prob_t predictive_f(int kk, int k, binProblem *bp)
{
        prob_t sum = 0;
        int i;

        for (i = 0; i < bp->bd->events; i++) {
                sum += bp->bd->counts[i]->content[kk][k] + bp->bd->alpha[i]->content[kk][k];
        }

        return LOG(  bp->bd->counts[bp->add_event.which]->content[kk][k]
                   + bp->bd->alpha [bp->add_event.which]->content[kk][k])
              -LOG(sum);
}

/******************************************************************************
 * KL Divergence
 ******************************************************************************/

static
prob_t KLPsiUtility_f(int kk, int k, void *data)
{
        binProblem *bp = (binProblem *)data;
        size_t i;
        prob_t count[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t gamma = bp->bd->gamma->content[kk][k];

        if (gamma == 0) {
                return -HUGE_VAL;
        }
        else {
                for (i = 0; i < bp->bd->events; i++) {
                        count[i] = countStatistic(i, kk, k, bp) + countAlpha(i, kk, k, bp);
                        alpha[i] = countAlpha(i, kk, k, bp);
                }
                if (kk <= bp->add_event.pos && bp->add_event.pos <= k) {
                        gamma   *= -predictive_f(kk, k, bp);
                }
                return LOG(gamma) + (mbeta_log(count, bp) - mbeta_log(alpha, bp));
        }
}

static
prob_t KLPsiUtility_g(int kk, int k, void *data)
{
        binProblem *bp = (binProblem *)data;
        size_t i;
        prob_t count[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t gamma = bp->bd->gamma->content[kk][k];
        prob_t sum   = 0;

        if (gamma == 0) {
                return -HUGE_VAL;
        }
        else {
                for (i = 0; i < bp->bd->events; i++) {
                        count[i] = countStatistic(i, kk, k, bp) + countAlpha(i, kk, k, bp);
                        alpha[i] = countAlpha(i, kk, k, bp);
                        sum     += count[i];
                }
                if (kk <= bp->add_event.pos && bp->add_event.pos <= k) {
                        /* negate this to get a positive term */
                        gamma *= -(gsl_sf_psi(count[bp->add_event.which])-gsl_sf_psi(sum));
                }
                return LOG(gamma) + (mbeta_log(count, bp) - mbeta_log(alpha, bp));
        }
}

static
void KLPsiUtility(utility_t* result, size_t i, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];
        prob_t sum1, sum2;
        size_t j;

        for (j = 0; j < bp->bd->events; j++) {
                bp->add_event.n     = 1;
                bp->add_event.pos   = i;
                bp->add_event.which = j;

                /* we use two functions for the utility to keep track
                 * of negative terms */

                /* localUtility_f */
                prombs(ev_log, bp->ak, bp->bd->prior_log, KLPsiUtility_f, bp->bd->L, minM(bp), (void *)bp);
                sum1    = sumModels(ev_log, bp);
                /* localUtility_g */
                prombs(ev_log, bp->ak, bp->bd->prior_log, KLPsiUtility_g, bp->bd->L, minM(bp), (void *)bp);
                sum2    = sumModels(ev_log, bp);

                result->utility->content[i] += +EXP(sum1 - evidence_ref);
                result->utility->content[i] += -EXP(sum2 - evidence_ref);
        }
}

/******************************************************************************
 * KL Multibin Divergence D_KL (P(B|x,y) || P(B))
 ******************************************************************************/

static
prob_t KLMultibinUtility_f(int kk, int k, void *data)
{
        binProblem *bp = (binProblem *)data;
        size_t i;
        prob_t count[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t gamma = bp->bd->gamma->content[kk][k];

        if (gamma == 0) {
                return -HUGE_VAL;
        }
        else {
                for (i = 0; i < bp->bd->events; i++) {
                        count[i] = countStatistic(i, kk, k, bp) + countAlpha(i, kk, k, bp);
                        alpha[i] = countAlpha(i, kk, k, bp);
                }
                if (kk <= bp->add_event.pos && bp->add_event.pos <= k) {
                        gamma   *= -predictive_f(kk, k, bp);
                }
                return LOG(gamma) + (mbeta_log(count, bp) - mbeta_log(alpha, bp));
        }
}

static
void KLMultibinUtility(utility_t* result, size_t i, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];
        prob_t sum;
        size_t j;

        for (j = 0; j < bp->bd->events; j++) {

                bp->add_event.n     = 1;
                bp->add_event.pos   = i;
                bp->add_event.which = j;

                /* localUtility_f */
                prombs(ev_log, bp->ak, bp->bd->prior_log, KLMultibinUtility_f, bp->bd->L, minM(bp), (void *)bp);
                sum = sumModels(ev_log, bp);

                result->utility->content[i] += -EXP(sum - evidence_ref);
                result->utility->content[i] += -result->expectation->content[j][i]*LOG(result->expectation->content[j][i]);
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
 * Threading
 ******************************************************************************/

static
void * computeKLUtility_thread(void* data_)
{
        pthread_data_t *data  = (pthread_data_t *)data_;
        binProblem *bp = data->bp;
        int i = data->i;
        utility_t *result    = (utility_t *)data->result;
        prob_t evidence_ref = data->evidence_ref;

        /* for each event compute its expectation */
        computeExpectation(result, i, evidence_ref, bp);

        /* compute utilities */
        if (bp->bd->options->kl_psi) {
                KLPsiUtility(result, i, evidence_ref, bp);
        }
        if (bp->bd->options->kl_multibin) {
                KLMultibinUtility(result, i, evidence_ref, bp);
        }

        return NULL;
}

/******************************************************************************
 * HMM utility
 ******************************************************************************/

static
prob_t hmm_he(int from, int to, binProblem* bp)
{
        size_t i;
        prob_t c1[bp->bd->events];
        prob_t c2[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t result = 0;

        /* alpha */
        for (i = 0; i < bp->bd->events; i++) {
                alpha[i] = countAlpha(i, from, from, bp); 
        }
        result -= mbeta_log(alpha, bp);
        /* counts */
        for (i = 0; i < bp->bd->events; i++) {
                c1[i] = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
                c2[i] = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
        }
        c2[bp->add_event.which] += bp->add_event.n;
        /* marginal */
        result += mbeta_log(c1, bp);
        /* expectation */
        result += mbeta_log(c2, bp) - mbeta_log(c1, bp);

        return result;
}

static
prob_t hmm_hu(int from, int to, binProblem* bp)
{
        size_t i;
        prob_t c[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t result = 0.0;
        prob_t sum    = 0.0;

        /* alpha */
        for (i = 0; i < bp->bd->events; i++) {
                alpha[i] = countAlpha(i, from, from, bp);
        }
        result -= mbeta_log(alpha, bp);
        /* counts */
        for (i = 0; i < bp->bd->events; i++) {
                c[i]  = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
        }
        if (bp->add_event.n) {
                c[bp->add_event.which] += bp->add_event.n;
        }
        result += mbeta_log(c, bp);
        /* psi */
        for (i = 0; i < bp->bd->events; i++) {
                sum += c[i];
        }
        result += LOG(-(gsl_sf_psi(c[bp->add_event.which]) - gsl_sf_psi(sum)));

        return result;
}

void hmm_computeUtility(
        utility_t *result,
        prob_t *forward,
        prob_t *backward,
        binProblem *bp)
{
        size_t i, j;
        prob_t tmp[bp->bd->L];

        for (i = 0; i < bp->bd->events; i++) {
                bp->add_event.n     = 1;
                bp->add_event.which = i;
                /* compute expectation */
                hmm_fb(tmp, forward, backward, &hmm_he, bp);
                for (j = 0; j < bp->bd->L; j++) {
                        result->expectation->content[i][j] = +EXP(tmp[j]);
                }
                /* compute utility */
                hmm_fb(tmp, forward, backward, &hmm_hu, bp);
                for (j = 0; j < bp->bd->L; j++) {
                        /* predictive entropy */
                        result->utility->content[j] += -result->expectation->content[i][j]*LOG(result->expectation->content[i][j]);
                        /* parameter entropy */
                        result->utility->content[j] += -EXP(tmp[j]);
                }
                bp->add_event.pos   = -1;
                bp->add_event.n     =  0;
        }
}

/******************************************************************************
 * HMM N-Step utility
 ******************************************************************************/

static
prob_t hmm_forward_rec(prob_t *result, size_t j, size_t to, prob_t (*f)(int, int, binProblem*), binProblem* bp)
{
        size_t k;

        prob_t tmp;

        tmp = to*LOG(bp->bd->options->rho) + f(0, to, bp);
        for (k = 0; k < j; k++) {
                tmp = logadd(tmp, (to-k-1)*LOG(bp->bd->options->rho) + LOG(1.0-bp->bd->options->rho) + result[k] + f(k+1, to, bp));
        }
        return tmp;
}

static
prob_t hmm_fb_rec(prob_t *forward, prob_t *backward, size_t j, prob_t (*f)(int, int, binProblem*), binProblem* bp)
{
        size_t k;
        prob_t tmp;

        tmp = hmm_forward_rec(forward, j, bp->bd->L-1, f, bp);
        for (k = j+1; k < bp->bd->L; k++) {
                tmp = logadd(tmp, LOG(1-bp->bd->options->rho) + hmm_forward_rec(forward, j, k-1, f, bp) + backward[k]);
        }
        return tmp;
}

void hmm_computeUtilityAt(
        size_t pos,
        vector_t* result,
        prob_t *forward,
        prob_t *backward,
        binProblem *bp)
{
        size_t i;
        prob_t utility = 0;
        prob_t tmp;

        /* compute predictive entropy */
        bp->add_event.n     =  0;
        bp->add_event.which = -1;
        for (i = 0; i < bp->bd->events; i++) {
                bp->add_event.n     = 1;
                bp->add_event.which = i;

                tmp = hmm_fb_rec(forward, backward, pos, &hmm_he, bp) - forward[bp->bd->L-1];

                utility += -EXP(tmp)*tmp;
                result->content[i] = EXP(tmp);
        }
        /* compute parameter entropy */
        for (i = 0; i < bp->bd->events; i++) {
                bp->add_event.n     = 1;
                bp->add_event.which = i;

                tmp = hmm_fb_rec(forward, backward, pos, &hmm_hu, bp) - forward[bp->bd->L-1];

                utility += -EXP(tmp);
        }
        bp->add_event.n     =  0;
        bp->add_event.which = -1;

        result->content[bp->bd->events] = utility;
}

prob_t hmm_computeDistance(
        size_t x,
        size_t y,
        prob_t *forward,
        prob_t *backward,
        binProblem *bp)
{
        prob_t expectation;
        prob_t result;
        prob_t tmp;

        /* add one event y */
        bp->add_event.n     = 1;
        bp->add_event.which = y;

        /* compute predictive entropy */
        tmp = hmm_fb_rec(forward, backward, x, &hmm_he, bp) - forward[bp->bd->L-1];

        expectation = EXP(tmp);
        /* result gets -log(expectation) */
        result      = -tmp;

        /* compute parameter entropy */
        tmp = hmm_fb_rec(forward, backward, x, &hmm_hu, bp) - forward[bp->bd->L-1];

        result += -EXP(tmp)/expectation;

        /* clean up */
        bp->add_event.n     =  0;
        bp->add_event.which = -1;

        return result;
}

/******************************************************************************
 * Main
 ******************************************************************************/

void computeKLUtility(
        utility_t *result,
        prob_t evidence_ref,
        binData* bd)
{
        /* compute utilities */
        threaded_computation((void *)result, evidence_ref, bd, computeKLUtility_thread,
                             "Computing utility: %.1f%%");
}
