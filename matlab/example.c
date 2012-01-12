
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include <mex.h>

#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/interface.h>

/* Simple and naive implementation of
 * ../adaptive_sampling/config.py:generate_alpha()
 */
static
void generate_alpha(vector_t* in, matrix_t* out) {
        int i, j, k;
        for (i = 0; i < in->size; i++) {
                for (j = i; j < in->size; j++) {
                        double sum = 0;
                        for (k = i; k <= j; k++) {
                                sum += in->content[k];
                        }
                        out->content[i][j] = sum/(double)(j-i+1);
                }
        }
}

/* Simple and naive implementation of
 * ../adaptive_sampling/statistics.py:countStatistic()
 */
static
void count_statistic(vector_t* in, matrix_t* out) {
        int i, j, k;
        for (i = 0; i < in->size; i++) {
                for (j = i; j < in->size; j++) {
                        double sum = 0;
                        for (k = i; k <= j; k++) {
                                sum += in->content[k];
                        }
                        out->content[i][j] = sum;
                }
        }
}

static
void print_vector(vector_t* v) {
        int i;
        for (i = 0; i < v->size; i++) {
                printf("%f ", v->content[i]);
        }
        printf("\n");
}

static
void print_matrix(matrix_t* m) {
        int i, j;
        for (i = 0; i < m->rows; i++) {
                for (j = 0; j < m->columns; j++) {
                        printf("%f ", m->content[i][j]);
                }
                printf("\n");
        }
}

void example(double x[]) {
        int i, j;
        /* in this example, we have L=6 possible stimuly */
        int L = 6;
        /* and two responses */
        int K = 2;

        /***********************************************************************
         * we have counts for two possible events: {success, failure} */
        vector_t* counts[K];
        counts[0] = alloc_vector(L); /* for each we have L stimuly */
        counts[1] = alloc_vector(L);

        /* set all counts to one */
        for (i = 0; i < L; i++ ) {
                counts[0]->content[i] = 1.0;
                counts[1]->content[i] = 1.0;
        }
        /* and compute the full matrix for the sampling library... */
        matrix_t* counts_m[K];
        counts_m[0] = alloc_matrix(L, L);
        counts_m[1] = alloc_matrix(L, L);
        count_statistic(counts[0], counts_m[0]);
        count_statistic(counts[1], counts_m[1]);

        /***********************************************************************
         * now the very same procedure for the pseudo counts alpha */
        vector_t* alpha[K];
        alpha[0] = alloc_vector(L);
        alpha[1] = alloc_vector(L);

        /* set all alpha to one */
        for (i = 0; i < L; i++ ) {
                alpha[0]->content[i] = 1.0;
                alpha[1]->content[i] = 1.0;
        }
        /* and one value to two */
        alpha[0]->content[4] = 2.0;
        /* and compute the full matrix for the sampling library... */
        matrix_t* alpha_m[K];
        alpha_m[0] = alloc_matrix(L, L);
        alpha_m[1] = alloc_matrix(L, L);
        generate_alpha(alpha[0], alpha_m[0]);
        generate_alpha(alpha[1], alpha_m[1]);

        /* show the result */
        printf("alpha[0]:\n");
        print_vector(alpha[0]);
        printf("alpha[1]:\n");
        print_vector(alpha[1]);
        printf("alpha_m[0]:\n");
        print_matrix(alpha_m[0]);
        printf("alpha_m[1]:\n");
        print_matrix(alpha_m[1]);

        /***********************************************************************
         * define beta */
        vector_t* beta = alloc_vector(L);
        for (i = 0; i < L; i++ ) {
                beta->content[i] = 1.0;
        }
        /* how simple... :)
         * but there is some more elaborate version in
         * ../adaptive_sampling/config.py:generate_beta() */

        /***********************************************************************
         * define gamma */
        matrix_t* gamma = alloc_matrix(L, L);
        for (i = 0; i < gamma->rows; i++) {
                for (j = 0; j < gamma->columns; j++) {
                        gamma->content[i][j] = 1.0;
                }
        }
        /* also this should be considered:
         * ../adaptive_sampling/config.py:generate_gamma() */

        /***********************************************************************
         * goodness... we need to define the options as well... */
        Options options;

        options.epsilon = 0.00001;    /* precision for the extended
                                       * prombs */
        options.verbose = 0;          /* do I use this?! */
        options.prombsTest = 0;
        options.bprob = 1;            /* we want break probabilities */
        options.threads = 1;          /* no threads for the moment */
        options.stacksize = 256*1024; /* some memory for prombs */
        options.algorithm = 0;        /* 0: prombs, 1: no idea,
                                       * 2: mgs */
        options.utility = 1;          /* compute utility gain */
        options.differential_entropy = 1;
        options.multibin_entropy = 0;
        options.effective_counts = 0;
        options.which = 0;            /* compute everything for the
                                       * first event (success) */
        options.samples[0] = 0;       /* mgs burn in */
        options.samples[1] = 0;       /* mgs samples */
        options.marginal = 1;         /* we also want the marginal */
        options.marginal_step = 0.01; /* step size at which the
                                       * density is evaluated */
        options.marginal_range.from = 0.0;
        options.marginal_range.to   = 1.0;
        options.n_moments = 2;        /* two moments are fine */
        options.n_marginals = (int)floor(1.0/options.marginal_step) + 1;
                                      /* why is not computed in the
                                       * library? */
        options.model_posterior = 1;  /* and we want a model
                                       * posterior! */

        /***********************************************************************
         * let's run prombs! */
        printf("\nRunning prombs!\n\n");

        BinningResult * result = binning(K, counts_m, alpha_m, beta, gamma, &options);

        /***********************************************************************
         * print the result */

        printf("moments:\n");
        print_matrix(result->moments);
        printf("break probability:\n");
        print_vector(result->bprob);
        printf("utility:\n");
        print_vector(result->utility);
        printf("model posterior:\n");
        print_vector(result->mpost);

        /***********************************************************************
         * free everything...  puhh... */
        free_vector(alpha[0]);
        free_vector(alpha[1]);
        free_matrix(alpha_m[0]);
        free_matrix(alpha_m[1]);
        free_vector(counts[0]);
        free_vector(counts[1]);
        free_matrix(counts_m[0]);
        free_matrix(counts_m[1]);
        free_vector(beta);
        free_matrix(gamma);

        free_matrix(result->moments);
        free_matrix(result->marginals);
        free_vector(result->bprob);
        free_vector(result->mpost);
        free_vector(result->utility);
        free(result);
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
        double *x,*y;
        mwSize mrows,ncols;
  
        /* Check for proper number of arguments. */
        if(nrhs!=1) {
                mexErrMsgTxt("One input required.");
        } else if(nlhs>1) {
                mexErrMsgTxt("Too many output arguments.");
        }
  
        /* The input must be a noncomplex scalar double.*/
        mrows = mxGetM(prhs[0]);
        ncols = mxGetN(prhs[0]);
        if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
            !(mrows==1 && ncols==1) ) {
                mexErrMsgTxt("Input must be a noncomplex scalar double.");
        }
  
        /* Create matrix for the return argument. */
        plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
  
        /* Assign pointers to each input and output. */
        x = mxGetPr(prhs[0]);
        y = mxGetPr(plhs[0]);
        y[0] = 1.;
        /* Call the timestwo subroutine. */
        example(x);
}

