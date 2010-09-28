
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

typedef struct {
        int T; // number of timesteps
        unsigned int trials;
        gsl_vector *counts;
        // hyperparameters
        double gamma;
        double sigma;
} binProblem;

unsigned int spikes(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int spikes = 0, i;

        for (i = ks; i<ke; i++) {
                spikes += (unsigned int)gsl_vector_get(bp->counts, i);
        }
        return spikes;
}

unsigned int gaps(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int gaps = 0, i;

        for (i = ks; i<ke; i++) {
                gaps += bp->trials - (unsigned int)gsl_vector_get(bp->counts, i);
        }
        return gaps;
}

double prior(binProblem *bp, unsigned int M)
{
//        double tmp = gsl_sf_gamma(bp->sigma + bp->gamma)/
//                (gsl_sf_gamma(bp->sigma)*gsl_sf_gamma(bp->gamma));
        double tmp = gsl_ran_beta_pdf(1, bp->sigma, bp->gamma);

        printf("gsl_ran_beta_pdf: %.20f\n", gsl_ran_beta_pdf(1,1,1));
        printf("gsl_ran_beta_pdf: %.20f\n", tmp);

        return (M+1)*tmp/gsl_sf_choose(bp->T-1,M);
}

double getIEC(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int s = spikes(bp, ks, ke);
        unsigned int g = gaps(bp, ks, ke);
        double r;

//        r = (gsl_sf_gamma(s+bp->sigma)*gsl_sf_gamma(g+bp->gamma))/
//                gsl_sf_gamma(s+bp->sigma+g+bp->gamma);
        r = gsl_ran_beta_pdf(1, s+bp->sigma, g+bp->gamma);

        return r;
}

void reset(binProblem *bp, size_t k[], size_t m, unsigned int M)
{
        size_t i;

             if (m == 0) k[0] = 0;
        else if (m <  M) k[m] = k[m-1]+1;

        for (i = m+1; i<M-1; i++) {
                k[i] = k[i-1]+1;
        }
        k[M-1] = bp->T-1;
}

int increment(binProblem *bp, size_t k[], unsigned int M)
{
        size_t i, m;

        for (i = 1; i<=M; i++) {
                m=M-i;
                if (k[m] < k[m+1]-1) {
                        k[m] += 1;
                        reset(bp, k, m+1, M);
                        return 1;
                }
        }
        return 0;
}

void printk(size_t k[], unsigned int M)
{
        int m;

        for (m = 0; m<M; m++) {
                printf("%ld ", k[m]);
        }
        printf("\n");
}

double evidence(binProblem *bp, unsigned int M)
{
        size_t k[M], m;
        double tmp1 = 0, tmp2;

        reset(bp, k, 0, M);

        do {
                tmp2 = 1;
                for (m = 0; m < M-1; m++) {
                        tmp2 *= getIEC(bp, k[m]+1, k[m+1])*getIEC(bp, 0, k[0]);
                }
                tmp1 += tmp2;
        } while(increment(bp, k, M));

        return prior(bp, M)*tmp1;
}

unsigned int getCount(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int spikes = 0, i;

        for (i = ks+1; i <= ke; i++) {
                spikes += (unsigned int)gsl_vector_get(bp->counts, i);
        }
        return spikes;
}

void evidences(binProblem *bp, double *ev)
{
        unsigned int k, kk, n, m, lb;
        unsigned int M = bp->T-1;
        unsigned int N = getCount(bp, -1, bp->T-1);
        double a[bp->T];

        for (k = 0; k <= bp->T-1; k++) {
                n = getCount(bp, -1, k);
//                a[k] = gsl_sf_fact(n)/gsl_sf_pow_int(k+1, n);
                a[k] = gsl_sf_lnfact(n) - n*gsl_sf_log(k+1);
//                a[k] = gsl_sf_exp(gsl_sf_lnfact(n) - n*gsl_sf_log(k+1));
//                printf("a[%d]=%f\n", k, a[k]);
        }
//        ev[0] = a[M]/gsl_sf_fact(N);
        ev[0] = a[M] - gsl_sf_lnfact(N);
//        ev[0] = gsl_sf_log(a[M]) - gsl_sf_lnfact(N);

        for (m = 1; m <= M; m++) {
                if (m==M) { lb = bp->T-1; }
                else      { lb = m; }

                for (k = bp->T-1; k >= lb; k--) {
                        a[k] = 0;
                        for (kk = m-1; kk <= k-1; kk++) {
                                double tmp;
                                n = getCount(bp, kk, k);
//                                a[k] += a[kk]*
//                                        gsl_sf_fact(n)/gsl_sf_pow_int(k-kk, n);
//                                gsl_sf_log(k-kk);
                                tmp = a[kk] + gsl_sf_lnfact(n) - n*gsl_sf_log(k-kk);
//                                tmp = gsl_sf_log(a[kk]) + gsl_sf_lnfact(n) - n*gsl_sf_log(k-kk);
                                a[k] = a[k] + gsl_sf_log(1 + gsl_sf_exp(tmp-a[k]));
//                                a[k] = a[k] + gsl_sf_exp(tmp);
//                                printf("a[%d]=%f\n", k, a[k]);
                        }
                }
//                ev[m] = a[bp->T-1]*
//                        gsl_sf_fact(bp->T-1-m)*gsl_sf_fact(m)*gsl_sf_fact(m)/
//                        (gsl_sf_fact(bp->T-1)*gsl_sf_fact(N+m));
                ev[m] = a[bp->T-1] +
                        gsl_sf_lnfact(bp->T-1-m) + 2*gsl_sf_lnfact(m) -
                        gsl_sf_lnfact(bp->T-1)   -   gsl_sf_lnfact(N+m);
        }
}

void pdensity(binProblem *bp, double *pdf)
{
        unsigned int i, j;
        double *ev1 = (double *)malloc((bp->T+1)*sizeof(double));
        double *ev2 = (double *)malloc((bp->T+1)*sizeof(double));
        double sum1, sum2;

        for (i=0; i<bp->T; i++) {
                evidences(bp, ev1);
                gsl_vector_set(
                        bp->counts, i,
                        gsl_vector_get(bp->counts, i)+1);
                evidences(bp, ev2);
                gsl_vector_set(
                        bp->counts, i,
                        gsl_vector_get(bp->counts, i)-1);

                sum1 = 0; sum2 = 0;
                for (j=0; j<bp->T; j++) {
                        sum1 = ev1[j];
                }
        }

        free(ev1);
        free(ev2);
}

gsl_matrix * bin(gsl_vector *counts, unsigned int trials)
{
        unsigned int i;
        gsl_matrix *m = gsl_matrix_alloc(2,2);
        binProblem bp;

        bp.trials = trials;
        bp.counts = counts;
        bp.T      = counts->size;
        bp.gamma  = 32;
        bp.sigma  = 1;

        double *ev = (double *)malloc((bp.T+1)*sizeof(double));
        evidences(&bp, ev);
        for (i = 0; i<=bp.T-1; i++) {
                printf("%.15f ", ev[i]);
        }
        printf("\n");
        free(ev);

        return m;
}
