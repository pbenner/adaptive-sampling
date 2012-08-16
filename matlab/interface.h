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
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <mex.h>

#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/interface.h>

static
void invalidOptions(const char* name) {
        mexPrintf("Field `%s' does not exist.\n", name);
        mexErrMsgTxt("Invalid options.");
}

static
double getScalar(const mxArray* array, const char* name) {
        mxArray* tmp = mxGetField(array, 0, name);

        if (tmp == 0) invalidOptions(name);

        return *mxGetPr(tmp);
}

static
mxArray* copyMatrixToMatlab(matrix_t* in) {
        mxArray *out = mxCreateDoubleMatrix(in->rows, in->columns, mxREAL);
        mwSize  nsubs = mxGetNumberOfDimensions(out);
        mwIndex* subs = mxCalloc(nsubs,sizeof(mwIndex));
        double* m = mxGetPr(out);
        mwIndex index;
        size_t i, j;

        for (i = 0; i < in->rows; i++) {
                for (j = 0; j < in->columns; j++) {
                        subs[0] = i;
                        subs[1] = j;
                        index = mxCalcSingleSubscript(out, nsubs, subs);
                        m[index] = in->content[i][j];
                }
        }
        mxFree(subs);

        return out;
}

static
mxArray* copyVectorToMatlab(vector_t* in) {
        mxArray *out = mxCreateDoubleMatrix(1, in->size, mxREAL);
        mwSize  nsubs = mxGetNumberOfDimensions(out);
        mwIndex* subs = mxCalloc(nsubs,sizeof(mwIndex));
        double* m = mxGetPr(out);
        mwIndex index;
        size_t i;

        for (i = 0; i < in->size; i++) {
                subs[1] = i;
                index = mxCalcSingleSubscript(out, nsubs, subs);
                m[index] = in->content[i];
        }
        mxFree(subs);

        return out;
}

static
options_t* getOptions(const mxArray *array)
{
        options_t* options = (options_t*)malloc(sizeof(options_t));
        mxArray* tmp;
        double*  ptr;

        options->model_posterior = getScalar(array, "model_posterior");
        options->kl_psi = getScalar(array, "kl_psi");
        options->kl_multibin  = getScalar(array, "kl_multibin");
        options->effective_counts = getScalar(array, "effective_counts");
        options->effective_posterior_counts = getScalar(array, "effective_posterior_counts");
        options->bprob = getScalar(array, "bprob");

        options->density = getScalar(array, "density");
        options->density_step = getScalar(array, "density_step");

        tmp = mxGetField(array, 0, "density_range");
        if (tmp == 0) invalidOptions("density_range");
        ptr = mxGetPr(tmp);
        options->density_range.from = ptr[0];
        options->density_range.to   = ptr[1];
        options->n_moments = getScalar(array, "n_moments");
        options->n_density = (int)floor(1.0/options->density_step) + 1;

        options->epsilon = getScalar(array, "epsilon");
        options->verbose = 0;
        options->prombsTest = 0;
        options->threads = getScalar(array, "threads");
        options->stacksize = getScalar(array, "stacksize");
        options->algorithm = getScalar(array, "algorithm");
        options->which = getScalar(array, "which");
        options->hmm = getScalar(array, "hmm");
        options->rho = getScalar(array, "rho");

        tmp = mxGetField(array, 0, "samples");
        if (tmp == 0) invalidOptions("samples");
        ptr = mxGetPr(tmp);
        options->samples[0] = ptr[0];
        options->samples[1] = ptr[1];

        return options;
}

static
void copyMatrix(matrix_t* to, const mxArray* from, size_t offset) {
        const size_t R = to->rows;
        const size_t C = to->columns;
        double* m = mxGetPr(from);
        size_t i, j;
        mwSize  nsubs = mxGetNumberOfDimensions(from);
        mwIndex* subs = mxCalloc(nsubs,sizeof(mwIndex));
        mwIndex index;

        for (i = 0; i < R; i++) {
                for (j = 0; j < C; j++) {
                        subs[0] = i+offset;
                        subs[1] = j;
                        index = mxCalcSingleSubscript(from, nsubs, subs);
                        to->content[i][j] = m[index];
                }
        }
        mxFree(subs);
}

static
void copyVector(vector_t* to, const mxArray* from) {
        double* v = mxGetPr(from);
        size_t i;
        mwSize  nsubs = mxGetNumberOfDimensions(from);
        mwIndex* subs = mxCalloc(nsubs,sizeof(mwIndex));
        mwIndex index;

        for (i = 0; i < to->size; i++) {
                subs[0] = i;
                index = mxCalcSingleSubscript(from, nsubs, subs);
                to->content[i] = v[index];
        }
        mxFree(subs);
}

static
matrix_t** getCounts(const mxArray* array) {
        size_t L = mxGetN(array);
        size_t K = mxGetM(array)/L;
        size_t k;
        matrix_t** counts = (matrix_t**)malloc((K+1)*sizeof(matrix_t*));

        for (k = 0; k < K; k++) {
                counts[k] = alloc_matrix(L, L);
                copyMatrix(counts[k], array, k*L);
        }
        counts[k] = NULL;

        return counts;
}

static
void freeCounts(matrix_t** counts) {
        size_t i;

        for (i = 0; counts[i]; i++) {
                free_matrix(counts[i]);
        }
        free(counts);
}

static
matrix_t** getAlpha(const mxArray* array, size_t K, size_t L) {
        if (mxGetN(array) != L || mxGetM(array)/L != K) {
                mexErrMsgTxt("Invalid dimension of alpha matrix.");
        }
        size_t k;
        matrix_t** alpha = (matrix_t**)malloc((K+1)*sizeof(matrix_t*));

        for (k = 0; k < K; k++) {
                alpha[k] = alloc_matrix(L, L);
                copyMatrix(alpha[k], array, k*L);
        }
        alpha[k] = NULL;

        return alpha;
}

static
void freeAlpha(matrix_t** alpha) {
        size_t i;

        for (i = 0; alpha[i]; i++) {
                free_matrix(alpha[i]);
        }
        free(alpha);
}

static
vector_t* getBeta(const mxArray* array, size_t L) {
        if (mxGetN(array) != L) {
                mexErrMsgTxt("Invalid dimension of beta vector.");
        }
        vector_t* beta = alloc_vector(L);

        copyVector(beta, array);

        return beta;
}

static
matrix_t* getGamma(const mxArray* array, size_t L) {
        if (mxGetN(array) != L || mxGetM(array) != L) {
                mexErrMsgTxt("Invalid dimension of gamma matrix.");
        }
        matrix_t* gamma = alloc_matrix(L, L);

        copyMatrix(gamma, array, 0);

        return gamma;
}
