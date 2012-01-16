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
Options* getOptions(const mxArray *array)
{
        Options* options = (Options*)malloc(sizeof(Options));
        mxArray* tmp;
        double*  ptr;

        options->model_posterior = getScalar(array, "model_posterior");
        options->utility = getScalar(array, "utility");
        options->differential_entropy = getScalar(array, "differential_entropy");
        options->multibin_entropy = getScalar(array, "multibin_entropy");
        options->effective_counts = getScalar(array, "effective_counts");
        options->bprob = getScalar(array, "bprob");

        options->marginal = getScalar(array, "marginal");
        options->marginal_step = getScalar(array, "marginal_step");

        tmp = mxGetField(array, 0, "marginal_range");
        if (tmp == 0) invalidOptions("marginal_range");
        ptr = mxGetPr(tmp);
        options->marginal_range.from = ptr[0];
        options->marginal_range.to   = ptr[1];
        options->n_moments = getScalar(array, "n_moments");
        options->n_marginals = (int)floor(1.0/options->marginal_step) + 1;

        options->epsilon = getScalar(array, "epsilon");
        options->verbose = 0;
        options->prombsTest = 0;
        options->threads = getScalar(array, "threads");
        options->stacksize = getScalar(array, "stacksize");
        options->algorithm = getScalar(array, "algorithm");
        options->which = getScalar(array, "which");

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
        size_t  L = mxGetN(array);
        size_t  K = mxGetM(array)/L;
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

static
void freeResult(BinningResult * result) {
        if (result->moments) {
                free_matrix(result->moments);
        }
        if (result->marginals) {
                free_matrix(result->marginals);
        }
        if (result->bprob) {
                free_vector(result->bprob);
        }
        if (result->mpost) {
                free_vector(result->mpost);
        }
        if (result->utility) {
                free_vector(result->utility);
        }
        free(result);
}

static
BinningResult * callBinning(const mxArray *prhs[]) {

        size_t  L = mxGetN(prhs[0]);
        size_t  K = mxGetM(prhs[0])/L;
        matrix_t** counts;
        matrix_t** alpha;
        vector_t* beta;
        matrix_t* gamma;
        Options* options;
        BinningResult * result;

        counts  = getCounts(prhs[0]);
        alpha   = getAlpha(prhs[1], K, L);
        beta    = getBeta(prhs[2], L);
        gamma   = getGamma(prhs[3], L);
        options = getOptions(prhs[4]);

        result = binning(K, counts, alpha, beta, gamma, options);

        freeCounts(counts);
        freeAlpha(alpha);
        free_vector(beta);
        free_matrix(gamma);
        free(options);

        return result;
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
void copyResult(BinningResult* result, mxArray *plhs[]) {
        const char **fnames;
        const int nfields = 5;

        /* allocate memory  for storing pointers */
        fnames = mxCalloc(nfields, sizeof(*fnames));
        fnames[0] = "moments";
        fnames[1] = "marginals";
        fnames[2] = "bprob";
        fnames[3] = "mpost";
        fnames[4] = "utility";

        plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames);
        mxFree((void *)fnames);

        if (result->moments) {
                mxSetField(plhs[0], 0, "moments", copyMatrixToMatlab(result->moments));
        }
        if (result->marginals) {
                mxSetField(plhs[0], 0, "marginals", copyMatrixToMatlab(result->marginals));
        }
        if (result->bprob) {
                mxSetField(plhs[0], 0, "bprob", copyVectorToMatlab(result->bprob));
        }
        if (result->mpost) {
                mxSetField(plhs[0], 0, "mpost", copyVectorToMatlab(result->mpost));
        }
        if (result->utility) {
                mxSetField(plhs[0], 0, "utility", copyVectorToMatlab(result->utility));
        }
}

/*  the gateway routine.  */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
        BinningResult * result;

        /* check proper input and output */
        if (nrhs != 5) {
                mexErrMsgTxt("Usage: binning(counts, alpha, beta, gamma, options).");
        }
        if (nlhs > 1) {
                mexErrMsgTxt("Too many output arguments.");
        }
        if (!mxIsClass(prhs[0], "double") || !mxIsClass(prhs[1], "double") ||
            !mxIsClass(prhs[2], "double") || !mxIsClass(prhs[3], "double")) {
                mexErrMsgTxt("All matrices must be of type double.");
        }
        if (!mxIsStruct(prhs[4])) {
                mexErrMsgTxt("Input must be a structure.");
        }

        result = callBinning(prhs);
        copyResult(result, plhs);
        freeResult(result);

        return;
}
