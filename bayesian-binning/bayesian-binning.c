
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_gamma.h>

typedef struct {
        int M; // number of bin boundaries
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
        double tmp = gsl_sf_gamma(bp->sigma + bp->gamma)/
                (gsl_sf_gamma(bp->sigma)*gsl_sf_gamma(bp->gamma));

        return (M+1)*tmp/gsl_sf_choose(bp->T-1,M);
}

double getIEC(binProblem *bp, size_t ks, size_t ke)
{
        unsigned int s = spikes(bp, ks, ke);
        unsigned int g = gaps(bp, ks, ke);
        double r;

        r = (gsl_sf_gamma(s+bp->sigma)*gsl_sf_gamma(g+bp->gamma))/
                gsl_sf_gamma(s+bp->sigma+g+bp->gamma);

        return r;
}

double evidence(binProblem *bp, unsigned int M)
{

        return 1.0;
}

gsl_matrix * bin(gsl_vector *counts, unsigned int trials)
{
        gsl_matrix *m = gsl_matrix_alloc(2,2);
        binProblem bp;

        bp.trials = trials;
        bp.counts = counts;
        bp.T      = counts->size;
        bp.M      = 0; // so far unknown
        bp.gamma  = 32;
        bp.sigma  = 1;

        printf("%.20f\n", getIEC(&bp, 2, 3));
//        spikes(&bp, 3, 6);

        gsl_matrix_set(m, 0, 0, 1);
        gsl_matrix_set(m, 0, 1, 2);
        gsl_matrix_set(m, 1, 0, 3);
        gsl_matrix_set(m, 1, 1, 4);

        return m;
}
