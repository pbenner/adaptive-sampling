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

#include "interface.h"

static
marginal_t * callBinning(const mxArray *prhs[]) {
        if (mxGetNumberOfDimensions(prhs[0]) != 3) {
                mexErrMsgTxt("Invalid dimension of counts matrix.");
        }

        size_t K = mxGetDimensions(prhs[0])[0];
        size_t L = mxGetDimensions(prhs[0])[1];

        matrix_t** counts;
        matrix_t** alpha;
        vector_t* beta;
        matrix_t* gamma;
        options_t* options;
        marginal_t * result;

        counts  = getCounts(prhs[0], K, L);
        alpha   = getAlpha(prhs[1], K, L);
        beta    = getBeta(prhs[2], L);
        gamma   = getGamma(prhs[3], L);
        options = getOptions(prhs[4]);

        result  = posterior(K, counts, alpha, beta, gamma, options);

        freeCounts(counts);
        freeAlpha(alpha);
        free_vector(beta);
        free_matrix(gamma);
        free(options);

        return result;
}

static
void copyResult(marginal_t* result, mxArray *plhs[]) {
        const char **fnames;
        const int nfields = 4;

        /* allocate memory  for storing pointers */
        fnames = mxCalloc(nfields, sizeof(*fnames));
        fnames[0] = "moments";
        fnames[1] = "density";
        fnames[2] = "bprob";
        fnames[3] = "mpost";

        plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames);
        mxFree((void *)fnames);

        if (result->moments) {
                mxSetField(plhs[0], 0, "moments", copyMatrixToMatlab(result->moments));
        }
        if (result->density) {
                mxSetField(plhs[0], 0, "density", copyMatrixToMatlab(result->density));
        }
        if (result->bprob) {
                mxSetField(plhs[0], 0, "bprob", copyVectorToMatlab(result->bprob));
        }
        if (result->mpost) {
                mxSetField(plhs[0], 0, "mpost", copyVectorToMatlab(result->mpost));
        }
}

static
void freeResult(marginal_t * result) {
        if (result->moments) {
                free_matrix(result->moments);
        }
        if (result->density) {
                free_matrix(result->density);
        }
        if (result->bprob) {
                free_vector(result->bprob);
        }
        if (result->mpost) {
                free_vector(result->mpost);
        }
        free(result);
}

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
        marginal_t * result;

        /* check proper input and output */
        if (nrhs != 5) {
                mexErrMsgTxt("Usage: binningPosterior(counts, alpha, beta, gamma, options).");
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
