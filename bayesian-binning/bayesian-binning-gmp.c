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

#include <gmp.h>

#include "datatypes.h"

typedef struct {
        // number of timesteps
        unsigned int T;
        gsl_vector *successes;
        gsl_vector *failures;
        gsl_vector *prior;
        mpf_t      *mprior;
        // type of the likelihood
        int likelihood;
        // hyperparameters
        double gamma;
        double sigma;
        // internal data
        gsl_matrix *success_m;
        gsl_matrix *failure_m;
        unsigned int add_success[2];
} binProblem;

mpz_t tmp1, tmp2;
mpf_t tmp3, tmp4;

static
mpf_t * allocMPFArray(size_t size)
{
        size_t i;
        mpf_t *a = (mpf_t *)malloc((size+1)*sizeof(mpf_t));

        for (i = 0; i<size; i++) {
                mpf_init(a[i]);
        }

        return a;
}

static
void freeMPFArray(mpf_t *a, size_t size)
{
        size_t i;

        for (i = 0; i<size; i++) {
                mpf_clear(a[i]);
        }
        free(a);
}

static
unsigned int successes(binProblem *bp, int ks, int ke)
{
        if (ks+1 <= bp->add_success[0] && bp->add_success[0] <= ke) {
                return gsl_matrix_get(bp->success_m, ks+1, ke) +
                        bp->add_success[1];
        }
        else if (ks+1 <= ke) {
                return gsl_matrix_get(bp->success_m, ks+1, ke);
        }
        else {
                return 0;
        }
}


static
unsigned int failures(binProblem *bp, int ks, int ke)
{
        if (ks+1 <= ke) {
                return gsl_matrix_get(bp->failure_m, ks+1, ke);
        }
        else {
                return 0;
        }
}

static
mpf_t * beta(unsigned int p, unsigned int q)
{
        mpz_fac_ui(tmp1, p-1);     // tmp1 = (p-1)!
        mpz_fac_ui(tmp2, q-1);     // tmp2 = (q-1)!
        mpz_mul(tmp1, tmp1, tmp2); // tmp1 = (p-1)!(q-1)!
        mpz_fac_ui(tmp2, p+q-1);   // tmp2 = (p+q-1)!
        
        mpf_set_z(tmp3, tmp1);     // tmp3 = (p-1)!(q-1)!
        mpf_set_z(tmp4, tmp2);     // tmp4 = (p+q-1)!
        mpf_div(tmp4, tmp4, tmp3); // tmp4 = (p+q-1)! / (p-1)!(q-1)!

        return &tmp4;
}

static
mpf_t * betaInv(unsigned int p, unsigned int q)
{
        mpz_fac_ui(tmp1, p-1);     // tmp1 = (p-1)!
        mpz_fac_ui(tmp2, q-1);     // tmp2 = (q-1)!
        mpz_mul(tmp1, tmp1, tmp2); // tmp1 = (p-1)!(q-1)!
        mpz_fac_ui(tmp2, p+q-1);   // tmp2 = (p+q-1)!
        
        mpf_set_z(tmp3, tmp1);     // tmp3 = (p-1)!(q-1)!
        mpf_set_z(tmp4, tmp2);     // tmp4 = (p+q-1)!
        mpf_div(tmp4, tmp3, tmp4); // tmp4 = (p-1)!(q-1)! / (p+q-1)!

        return &tmp4;
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
void computeMPrior(binProblem *bp)
{
        unsigned int m;
        // multinomial
        if (bp->likelihood == 1) {
                for (m = 0; m < bp->T; m++) {
                        unsigned int N = successes(bp, -1, bp->T-1);
                        mpz_fac_ui(tmp1, bp->T-1-m);
                        mpz_fac_ui(tmp2, m);
                        mpz_mul(tmp2, tmp2, tmp2);
                        mpz_mul(tmp2, tmp1, tmp2);
                        mpf_set_z (tmp3, tmp2);
                        mpz_fac_ui(tmp1, bp->T-1);
                        mpz_fac_ui(tmp2, N+m);
                        mpz_mul(tmp2, tmp1, tmp2);
                        mpf_set_z (tmp4, tmp2);
                        mpf_div(bp->mprior[m], tmp3, tmp4);
                }
        }
        // binomial
        if (bp->likelihood == 2) {
                for (m = 0; m < bp->T; m++) {
                        mpf_pow_ui(tmp4, *beta(bp->sigma, bp->gamma), m+1);
                        mpz_fac_ui(tmp1, m);          // tmp1 = M!
                        mpz_fac_ui(tmp2, bp->T-1-m);  // tmp2 = (K-M-1)!
                        mpz_mul(tmp2, tmp1, tmp2);    // tmp2 = (K-M-1)!M!
                        mpf_set_z(tmp3, tmp2);        // tmp3 = (K-M-1)!M!
                        mpf_mul(tmp3, tmp3, tmp4);    // tmp3 = (K-M-1)!M!Beta

                        mpz_fac_ui(tmp1, bp->T-1);    // tmp1 = (K-1)!
                        mpf_set_z(tmp4, tmp1);        // tmp4 = (K-1)!
                        mpf_div(bp->mprior[m], tmp3, tmp4);
                }
        }
}

static
mpf_t * iec(binProblem *bp, int kk, int k)
{
        // multinomial
        if (bp->likelihood == 1) {
                unsigned int n;
                n = successes(bp, kk, k);
                mpz_fac_ui(tmp1, n);
                mpf_set_ui(tmp3, k-kk);
                mpf_pow_ui(tmp4, tmp3, n);
                mpf_set_z (tmp3, tmp1);
                mpf_div(tmp4, tmp3, tmp4);
        }
        // binomial
        if (bp->likelihood == 2) {
                unsigned int s = successes(bp, kk, k);
                unsigned int f = failures (bp, kk, k);
                mpf_set(tmp4, *betaInv(s+bp->sigma, f+bp->gamma));
        }
        return &tmp4;
}

static
int minM(binProblem *bp)
{
        int i;
        for (i = bp->T-1; i>0; i--) {
                if (gsl_vector_get(bp->prior, i) > 0) {
                        return i;
                }
        }
        return i;
}

static
void evidences(binProblem *bp, mpf_t *ev)
{
        unsigned int m, lb;
        unsigned int M = minM(bp);
        mpf_t *a = allocMPFArray(bp->T);
        int k, kk;
        for (k = 0; k <= bp->T-1; k++) {
                mpf_set(a[k], *iec(bp, -1, k));
        }
        mpf_mul(ev[0], a[bp->T-1], bp->mprior[0]);

        for (m = 1; m <= M; m++) {
                if (m==M) { lb = bp->T-1; }
                else      { lb = m; }

                for (k = bp->T-1; k >= lb; k--) {
                        mpf_set_ui(a[k], 0);
                        for (kk = m-1; kk <= k-1; kk++) {
                                mpf_mul(tmp3, a[kk], *iec(bp, kk, k));
                                mpf_add(a[k], a[k],  tmp3);
                        }
                }
                mpf_mul(ev[m], a[bp->T-1], bp->mprior[m]);
        }

        freeMPFArray(a, bp->T);
}

static
void pdensity(binProblem *bp, double *pdf, double *var, double *mpost)
{
        unsigned int i, j;
        mpf_t *ev1 = allocMPFArray(bp->T+1);
        mpf_t *ev2 = allocMPFArray(bp->T+1);
        mpf_t *ev3 = allocMPFArray(bp->T+1);
        mpf_t sum1, sum2, sum3;
        mpf_t tmp1, tmp2;

        mpf_init(sum1);
        mpf_init(sum2);
        mpf_init(sum3);
        mpf_init(tmp1);
        mpf_init(tmp2);

        // compute evidence
        computeMPrior(bp);
        evidences(bp, ev1);
        // compute model posteriors
        mpf_set_ui(sum1, 0);
        for (j=0; j<bp->T; j++) {
                mpf_set_d(tmp2, gsl_vector_get(bp->prior, j));
                mpf_mul(tmp1, ev1[j], tmp2); // tmp1 = P(D|m)P(m)
                mpf_add(sum1, sum1, tmp1);   // sum1 = Sum_{m \in M} P(D|m)P(m)
        }
        for (j=0; j<bp->T; j++) {
                if (gsl_vector_get(bp->prior, j) > 0) {
                        mpf_set_d(tmp2, gsl_vector_get(bp->prior, j));
                        mpf_mul(tmp1, ev1[j], tmp2);
                        mpf_div(tmp1, tmp1, sum1);
                        mpost[j] = mpf_get_d(tmp1);
                }
                else {
                        mpost[j] = 0;
                }
        }
        // for each timestep compute expectation and variance
        // from the model average
        notice(NONE, "T: %d", bp->T);
        for (i=0; i<bp->T; i++) {
                notice(NONE, "%.1f%%", (float)100*i/bp->T);
                // expectation
                bp->add_success[0] = i;
                bp->add_success[1] = 1;
                computeMPrior(bp);
                evidences(bp, ev2);

                // variance
                bp->add_success[1] = 2;
                computeMPrior(bp);
                evidences(bp, ev3);

                mpf_set_ui(sum2, 0);
                mpf_set_ui(sum3, 0);
                for (j=0; j<bp->T; j++) {
                        if (gsl_vector_get(bp->prior, j) > 0) {
                                mpf_set_d(tmp2, gsl_vector_get(bp->prior, j));
                                mpf_mul(tmp1, ev2[j], tmp2); // mult. with prior
                                mpf_add(sum2, sum2, tmp1);
                                mpf_mul(tmp1, ev3[j], tmp2); // mult. with prior
                                mpf_add(sum3, sum3, tmp1);
                        }
                }
                mpf_div(tmp1, sum2, sum1); // tmp1 = P(D' |M)/P(D|M)
                mpf_div(tmp2, sum3, sum1); // tmp2 = P(D''|M)/P(D|M)
                pdf[i] = mpf_get_d(tmp1);

                mpf_pow_ui(tmp1, tmp1, 2); // tmp1 = (P(D' |M)/P(D|M))^2
                mpf_sub(tmp2, tmp2, tmp1); // tmp2 =  P(D''|M)/P(D|M) - (P(D' |M)/P(D|M))^2
                var[i] = mpf_get_d(tmp2);
        }

        mpf_clear(sum1);
        mpf_clear(sum2);
        mpf_clear(sum3);
        mpf_clear(tmp1);
        mpf_clear(tmp2);
        freeMPFArray(ev1, bp->T+1);
        freeMPFArray(ev2, bp->T+1);
        freeMPFArray(ev3, bp->T+1);
}

gsl_matrix * bin_gmp(
        gsl_vector *successes_v,
        gsl_vector *failures_v,
        gsl_vector *prior,
        Options *options)
{
        size_t K = successes_v->size;
        gsl_matrix *m = gsl_matrix_alloc(3, K);
        double pdf[K];
        double var[K];
        double mpost[K];
        unsigned int i;
        binProblem bp;

        verbose       = options->verbose;

        bp.successes  = successes_v;
        bp.failures   = failures_v;
        bp.prior      = prior;
        bp.T          = successes_v->size;
        bp.sigma      = options->sigma;
        bp.gamma      = options->gamma;
        bp.likelihood = options->likelihood;
        bp.mprior     = allocMPFArray(K);
        bp.success_m  = gsl_matrix_alloc(K, K);
        bp.failure_m  = gsl_matrix_alloc(K, K);
        bp.add_success[0] = 0; // number of the bin
        bp.add_success[1] = 0; // how many successes

        mpz_init(tmp1);
        mpz_init(tmp2);
        mpf_init(tmp3);
        mpf_init(tmp4);

        computeSuccesses(&bp);

        pdensity(&bp, pdf, var, mpost);
        for (i = 0; i <= bp.T-1; i++) {
                gsl_matrix_set(m, 0, i, pdf[i]);
                gsl_matrix_set(m, 1, i, var[i]);
                gsl_matrix_set(m, 2, i, mpost[i]);
                notice(NONE, "pdf[%03d]=%f var[%03d]=%f", i, pdf[i], i, var[i]);
        }

        mpz_clear(tmp1);
        mpz_clear(tmp2);
        mpf_clear(tmp3);
        mpf_clear(tmp4);

        freeMPFArray(bp.mprior, K);
        gsl_matrix_free(bp.success_m);
        gsl_matrix_free(bp.failure_m);

        return m;
}
