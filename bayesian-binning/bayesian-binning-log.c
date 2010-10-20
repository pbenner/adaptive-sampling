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
#include <math.h>

#include <exception.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "datatypes.h"
#include "logadd.h"
#include "prombs.h"

typedef struct {
        // number of timesteps
        unsigned int T;
        gsl_vector  *successes;
        gsl_vector  *failures;
        gsl_vector  *mprior;     // P(m_B)
        long double *prior_log;  // P(p,B|m_B)
        // type of the likelihood
        int likelihood;
        // hyperparameters
        long double gamma;
        long double sigma;
        // internal data
        gsl_matrix *success_m;
        gsl_matrix *failure_m;
        unsigned int add_success[2];
} binProblem;

static
unsigned int successes(binProblem *bp, int ks, int ke)
{
        if (ks <= bp->add_success[0] && bp->add_success[0] <= ke) {
                return gsl_matrix_get(bp->success_m, ks, ke) +
                        bp->add_success[1];
        }
        else if (ks <= ke) {
                return gsl_matrix_get(bp->success_m, ks, ke);
        }
        else {
                return 0;
        }
}

static
unsigned int failures(binProblem *bp, int ks, int ke)
{
        if (ks <= ke) {
                return gsl_matrix_get(bp->failure_m, ks, ke);
        }
        else {
                return 0;
        }
}

static
long double beta_log(unsigned int p, unsigned int q)
{
        return gsl_sf_lngamma(p+q) - gsl_sf_lngamma(p) - gsl_sf_lngamma(q);
}

static
long double betaInv_log(unsigned int p, unsigned int q)
{
         return gsl_sf_lngamma(p) + gsl_sf_lngamma(q) - gsl_sf_lngamma(p+q);
}

static
unsigned int computeSuccesses(binProblem *bp)
{
        int ks, ke, i;
        unsigned int s, f;

        for (ks = 0; ks < bp->T; ks++) {
                for (ke = ks; ke < bp->T; ke++) {
                        s = 0; f = 0;
                        for (i = ks; i <= ke; i++) {
                                s += (unsigned int)gsl_vector_get(bp->successes, i);
                                f += (unsigned int)gsl_vector_get(bp->failures,  i);
                        }
                        gsl_matrix_set(bp->success_m, ks, ke, s);
                        gsl_matrix_set(bp->failure_m, ks, ke, f);
                }
        }
        return s;
}

static
void computePrior_log(binProblem *bp)
{
        unsigned int m;
        // multinomial
        if (bp->likelihood == 1) {
                for (m = 0; m < bp->T; m++) {
                        unsigned int N = successes(bp, 0, bp->T-1);

                        bp->prior_log[m] = gsl_sf_lnfact(bp->T-1-m) + 2*gsl_sf_lnfact(m)
                                - gsl_sf_lnfact(bp->T-1) - gsl_sf_lnfact(N+m);
                }
        }
        // binomial
        if (bp->likelihood == 2) {
                for (m = 0; m < bp->T; m++) {
                        bp->prior_log[m] = (m+1)*beta_log(bp->sigma, bp->gamma)
                                - gsl_sf_lnchoose(bp->T-1, m);
                }
        }
}

/* P(D|B,m_B)/P(p|m_B) */
static
long double iec_log(binProblem *bp, int kk, int k)
{
        // multinomial
        if (bp->likelihood == 1) {
                unsigned int n = successes(bp, kk, k);
                return gsl_sf_lnfact(n) - n*logl(k-kk+1);
        }
        // binomial
        if (bp->likelihood == 2) {
                unsigned int s = successes(bp, kk, k);
                unsigned int f = failures (bp, kk, k);
                return betaInv_log(s+bp->sigma, f+bp->gamma);
        }
        err(NONE, "Unknown likelihood function.");
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
void execPrombs(binProblem *bp, long double *ev_log, int pos)
{
        long double f(int i, int j)
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

        prombs(ev_log, &f, bp->prior_log, bp->T, minM(bp));
}

static
long double computeEvidence(binProblem *bp, long double *ev_log)
{
        long double sum;
        unsigned int j;

        execPrombs(bp, ev_log, -1);
        sum = -HUGE_VAL;
        for (j=0; j<bp->T; j++) {
                if (gsl_vector_get(bp->mprior, j) > 0) {
                        long double mprior = gsl_vector_get(bp->mprior, j);
                        sum = logadd(sum, ev_log[j] + logl(mprior));
                }
        }
        return sum;
}

static
void computeModelPosteriors(binProblem *bp, long double *ev_log, long double *mpost, long double P_D)
{
        unsigned int j;

        for (j=0; j<bp->T; j++) {
                if (gsl_vector_get(bp->mprior, j) > 0) {
                        long double mprior = gsl_vector_get(bp->mprior, j);
                        mpost[j] = expl(ev_log[j] + logl(mprior) - P_D);
                }
                else {
                        mpost[j] = 0;
                }
        }
}

static
void computeBreakProbabilities(binProblem *bp, long double *bprob, long double P_D)
{
        long double ev_log[bp->T];
        long double sum;
        unsigned int i, j;

        for (i=0; i<bp->T; i++) {
                notice(NONE, "break prob.: %.1f%%", (float)100*i/bp->T);
                execPrombs(bp, ev_log, i);
                sum = -HUGE_VAL;
                for (j=0; j<bp->T; j++) {
                        if (gsl_vector_get(bp->mprior, j) > 0) {
                                long double mprior = gsl_vector_get(bp->mprior, j);
                                sum = logadd(sum, ev_log[j] + logl(mprior));
                        }
                }
                bprob[i] = expl(sum - P_D);
        }
}

static
void computeBinning(
        binProblem *bp, long double *pdf, long double *var, long double *bprob, long double *mpost,
        Options *options)
{
        long double ev1_log[bp->T], ev2_log[bp->T], ev3_log[bp->T];
        long double sum1; // P(D) = sum_{m_B \in M} P(D|m_B)P(m_B)
        long double sum2; // E[p|D]
        long double sum3; // Var[p|D]
        unsigned int i, j;

        // compute evidence P(D)
        computePrior_log(bp);
        sum1 = computeEvidence(bp, ev1_log);
        // compute model posteriors P(m_B|D)
        computeModelPosteriors(bp, ev1_log, mpost, sum1);
        // break probability
        if (options->bprob) {
                computeBreakProbabilities(bp, bprob, sum1);
        }
        // for each timestep compute expectation and variance
        // from the model average
        notice(NONE, "T: %d", bp->T);
        for (i=0; i<bp->T; i++) {
                notice(NONE, "exp./var.: %.1f%%", (float)100*i/bp->T);
                // expectation
                bp->add_success[0] = i;
                bp->add_success[1] = 1;
                computePrior_log(bp);
                execPrombs(bp, ev2_log, -1);

                // variance
                bp->add_success[1] = 2;
                computePrior_log(bp);
                execPrombs(bp, ev3_log, -1);

                sum2 = -HUGE_VAL;
                sum3 = -HUGE_VAL;
                for (j=0; j<bp->T; j++) {
                        if (gsl_vector_get(bp->mprior, j) > 0) {
                                long double mprior = gsl_vector_get(bp->mprior, j);
                                sum2 = logadd(sum2, ev2_log[j] + logl(mprior));
                                sum3 = logadd(sum3, ev3_log[j] + logl(mprior));
                        }
                }
                // P(D' |M)/P(D|M)
                pdf[i] = expl(sum2 - sum1);
                // P(D''|M)/P(D|M) - (P(D' |M)/P(D|M))^2
                var[i] = expl(sum3 - sum1) - expl(sum2 - sum1)*expl(sum2 - sum1);
        }
}

gsl_matrix * bin_log(
        gsl_vector *successes_v,
        gsl_vector *failures_v,
        gsl_vector *mprior,
        Options *options)
{
        size_t K = successes_v->size;
        gsl_matrix *m = gsl_matrix_alloc(4, K);
        long double pdf[K];
        long double var[K];
        long double bprob[K];
        long double mpost[K];
        long double prior_log[K];
        unsigned int i;
        binProblem bp;

        verbose       = options->verbose;

        bp.successes  = successes_v;
        bp.failures   = failures_v;
        bp.mprior     = mprior;
        bp.T          = successes_v->size;
        bp.sigma      = options->sigma;
        bp.gamma      = options->gamma;
        bp.likelihood = options->likelihood;
        bp.prior_log  = prior_log;
        bp.success_m  = gsl_matrix_alloc(K, K);
        bp.failure_m  = gsl_matrix_alloc(K, K);
        bp.add_success[0] = 0; // number of the bin
        bp.add_success[1] = 0; // how many successes

        computeSuccesses(&bp);
        computeBinning(&bp, pdf, var, bprob, mpost, options);
        for (i = 0; i <= bp.T-1; i++) {
                gsl_matrix_set(m, 0, i, pdf[i]);
                gsl_matrix_set(m, 1, i, var[i]);
                gsl_matrix_set(m, 2, i, bprob[i]);
                gsl_matrix_set(m, 3, i, mpost[i]);
                notice(NONE, "pdf[%03d]=%Lf var[%03d]=%Lf", i, pdf[i], i, var[i]);
        }

        gsl_matrix_free(bp.success_m);
        gsl_matrix_free(bp.failure_m);

        return m;
}
