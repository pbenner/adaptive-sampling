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
        int T; // number of timesteps
        unsigned int trials;
        gsl_vector *counts;
        gsl_vector *prior;
        // hyperparameters
        double gamma;
        double sigma;
} binProblem;

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
unsigned int spikes(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int spikes = 0, i;

        for (i = ks; i<ke; i++) {
                spikes += (unsigned int)gsl_vector_get(bp->counts, i);
        }
        return spikes;
}

static
unsigned int gaps(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int gaps = 0, i;

        for (i = ks; i<ke; i++) {
                gaps += bp->trials - (unsigned int)gsl_vector_get(bp->counts, i);
        }
        return gaps;
}

static
unsigned int getCount(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int spikes = 0, i;

        for (i = ks+1; i <= ke; i++) {
                spikes += (unsigned int)gsl_vector_get(bp->counts, i);
        }
        return spikes;
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
        unsigned int k, kk, n, m, lb;
        unsigned int M = minM(bp);
        unsigned int N = getCount(bp, -1, bp->T-1);
        mpf_t *a = allocMPFArray(bp->T);
        mpz_t tmp1, tmp2;
        mpf_t tmp3, tmp4;

        mpz_init(tmp1);
        mpz_init(tmp2);
        mpf_init(tmp3);
        mpf_init(tmp4);

        for (k = 0; k <= bp->T-1; k++) {
                n = getCount(bp, -1, k);
                mpz_fac_ui(tmp1, n);
                mpf_set_ui(tmp3, k+1);
                mpf_pow_ui(tmp4, tmp3, n);
                mpf_set_z (tmp3, tmp1);
                mpf_div(a[k], tmp3, tmp4);
        }
        mpz_fac_ui(tmp1, N);
        mpf_set_z (tmp3, tmp1);
        mpf_div(ev[0], a[M], tmp3);

        for (m = 1; m <= M; m++) {
                if (m==M) { lb = bp->T-1; }
                else      { lb = m; }

                for (k = bp->T-1; k >= lb; k--) {
                        mpf_set_ui(a[k], 0);
                        for (kk = m-1; kk <= k-1; kk++) {
                                n = getCount(bp, kk, k);
                                mpz_fac_ui(tmp1, n);
                                mpf_set_ui(tmp3, k-kk);
                                mpf_pow_ui(tmp4, tmp3, n);
                                mpf_set_z (tmp3, tmp1);
                                mpf_div(tmp4, tmp3, tmp4);
                                mpf_mul(tmp3, a[kk],tmp4);
                                mpf_add(a[k], a[k], tmp3);
                        }
                }
                mpz_fac_ui(tmp1, bp->T-1-m);
                mpz_fac_ui(tmp2, m);
                mpz_mul(tmp2, tmp2, tmp2);
                mpz_mul(tmp2, tmp1, tmp2);
                mpf_set_z (tmp3, tmp2);
                mpz_fac_ui(tmp1, bp->T-1);
                mpz_fac_ui(tmp2, N+m);
                mpz_mul(tmp2, tmp1, tmp2);
                mpf_set_z (tmp4, tmp2);
                mpf_div(tmp4, tmp3, tmp4);
                mpf_mul(ev[m], a[bp->T-1], tmp4);
        }

        mpz_clear(tmp1);
        mpz_clear(tmp2);
        mpf_clear(tmp3);
        mpf_clear(tmp4);
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
        evidences(bp, ev1);
        // compute model posteriors
        mpf_set_ui(sum1, 0);
        for (j=0; j<bp->T; j++) {
                mpf_set_d(tmp2, gsl_vector_get(bp->prior, j));
                mpf_mul(tmp1, ev1[j], tmp2);
                mpf_add(sum1, sum1, tmp1);
        }
        for (j=0; j<bp->T; j++) {
                mpf_set_d(tmp2, gsl_vector_get(bp->prior, j));
                mpf_mul(tmp1, ev1[j], tmp2);
                mpf_div(tmp1, tmp1, sum1);
                mpost[j] = mpf_get_d(tmp1);
        }
        // for each timestep compute expectation and variance
        // from the model average
        notice(NONE, "T: %d", bp->T);
        for (i=0; i<bp->T; i++) {
                notice(NONE, "%.1f%%", (float)100*i/bp->T);
                // expectation
                gsl_vector_set(
                        bp->counts, i,
                        gsl_vector_get(bp->counts, i)+1);
                evidences(bp, ev2);
                // variance
                gsl_vector_set(
                        bp->counts, i,
                        gsl_vector_get(bp->counts, i)+1);
                evidences(bp, ev3);
                // reset data
                gsl_vector_set(
                        bp->counts, i,
                        gsl_vector_get(bp->counts, i)-2);

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
                mpf_div(tmp1, sum2, sum1);
                mpf_div(tmp2, sum3, sum1);
                pdf[i] = mpf_get_d(tmp1);
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

gsl_matrix * bin(
        gsl_vector *counts,
        unsigned int trials,
        gsl_vector *prior,
        Options *options)
{
        gsl_matrix *m = gsl_matrix_alloc(3, counts->size);
        double pdf[counts->size];
        double var[counts->size];
        double mpost[counts->size];
        unsigned int i;
        binProblem bp;
        verbose = options->verbose;

        bp.trials = trials;
        bp.counts = counts;
        bp.prior  = prior;
        bp.T      = counts->size;
        bp.gamma  = 32;
        bp.sigma  = 1;

        pdensity(&bp, pdf, var, mpost);
        for (i = 0; i<=bp.T-1; i++) {
                gsl_matrix_set(m, 0, i, pdf[i]);
                gsl_matrix_set(m, 1, i, var[i]);
                gsl_matrix_set(m, 2, i, mpost[i]);
        }

        return m;
}
