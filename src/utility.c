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
prob_t KLComponentUtility_f(int kk, int k, void *data)
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
prob_t KLComponentUtility_g(int kk, int k, void *data)
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
                        gamma *= -(gsl_sf_psi(count[bp->add_event.which])-gsl_sf_psi(sum));
                }
                return LOG(gamma) + (mbeta_log(count, bp) - mbeta_log(alpha, bp));
        }
}

static
prob_t KLComponentUtility(size_t i, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];
        prob_t sum1, sum2, result = 0;
        size_t j;

        for (j = 0; j < bp->bd->events; j++) {

                bp->add_event.n     = 1;
                bp->add_event.pos   = i;
                bp->add_event.which = j;

                /* localUtility_f */
                prombs(ev_log, bp->ak, bp->bd->prior_log, KLComponentUtility_f, bp->bd->L, minM(bp), (void *)bp);
                sum1    = sumModels(ev_log, bp);
                /* localUtility_g */
                prombs(ev_log, bp->ak, bp->bd->prior_log, KLComponentUtility_g, bp->bd->L, minM(bp), (void *)bp);
                sum2    = sumModels(ev_log, bp);

                result += +EXP(sum1 - evidence_ref);
                result += -EXP(sum2 - evidence_ref);
        }
        return result;
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
prob_t KLMultibinUtility(size_t i, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];
        prob_t sum, result = 0;
        prob_t tmp;
        size_t j;

        for (j = 0; j < bp->bd->events; j++) {

                tmp = expectation(j, i, evidence_ref, bp);

                bp->add_event.n     = 1;
                bp->add_event.pos   = i;
                bp->add_event.which = j;

                /* localUtility_f */
                prombs(ev_log, bp->ak, bp->bd->prior_log, KLMultibinUtility_f, bp->bd->L, minM(bp), (void *)bp);
                sum     = sumModels(ev_log, bp);

                result += -EXP(sum - evidence_ref);
                result += -tmp*LOG(tmp);
        }
        return result;
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
        vector_t *result    = (vector_t *)data->result;
        prob_t evidence_ref = data->evidence_ref;

        result->content[i] = 0;
        if (bp->bd->options->kl_component) {
                result->content[i] += KLComponentUtility(i, evidence_ref, bp);
        }
        if (bp->bd->options->kl_multibin) {
                result->content[i] += KLMultibinUtility(i, evidence_ref, bp);
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
        prob_t ca[bp->bd->events];
        prob_t cb[bp->bd->events];

        for (i = 0; i < bp->bd->events; i++) {
                ca[i] = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
                cb[i] = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
        }
        ca[bp->add_event.which] += bp->add_event.n;

        return mbeta_log(ca, bp) - mbeta_log(cb, bp);
}

static
prob_t hmm_hu(int from, int to, binProblem* bp)
{
        size_t i;
        prob_t sum1, sum2 = 0;

        sum1 = countAlpha(bp->add_event.which, from, from, bp) + countStatistic(bp->add_event.which, from, to, bp);
        for (i = 0; i < bp->bd->events; i++) {
                sum2 += countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
        }
        return LOG(-(gsl_sf_psi(sum1+1) - gsl_sf_psi(sum2+1))*sum1/sum2);
}

void hmm_computeUtility(
        vector_t *utility,
        vector_t *forward,
        vector_t *backward,
        binProblem *bp)
{
        size_t i, j;
        vector_t* tmp = alloc_vector(bp->bd->L);
        matrix_t* expectation = alloc_matrix(bp->bd->events, bp->bd->L);
        matrix_t* utility_tmp = alloc_matrix(bp->bd->events, bp->bd->L);

        for (i = 0; i < bp->bd->events; i++) {
                bp->add_event.n     = 1;
                bp->add_event.which = i;
                // compute expectation
                hmm_fb(tmp, forward, backward, &hmm_he, bp);
                for (j = 0; j < bp->bd->L; j++) {
                        expectation->content[i][j] = +EXP(tmp->content[j]);
                }
                // compute utility
                hmm_fb(tmp, forward, backward, &hmm_hu, bp);
                for (j = 0; j < bp->bd->L; j++) {
                        utility_tmp->content[i][j] = -LOG(expectation->content[i][j])-EXP(tmp->content[j])/expectation->content[i][j];
                }
                bp->add_event.pos   = -1;
                bp->add_event.n     = 0;
        }
        // compute average
        for (j = 0; j < bp->bd->L; j++) {
                utility->content[j] = 0;
                for (i = 0; i < bp->bd->events; i++) {
                        utility->content[j] += expectation->content[i][j]*utility_tmp->content[i][j];
                }
        }
        free_vector(tmp);
        free_matrix(expectation);
        free_matrix(utility_tmp);
}


/******************************************************************************
 * Main
 ******************************************************************************/

void computeKLUtility(
        vector_t *result,
        prob_t evidence_ref,
        binData* bd)
{
        /* compute utilities */
        threaded_computation((void *)result, evidence_ref, bd, computeKLUtility_thread,
                             "Computing utility: %.1f%%");
}
