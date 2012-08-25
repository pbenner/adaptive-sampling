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
utility_t * callUtility(const mxArray *prhs[]) {
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
        utility_t * result;

        counts  = getCounts(prhs[0], K, L);
        alpha   = getAlpha(prhs[1], K, L);
        beta    = getBeta(prhs[2], L);
        gamma   = getGamma(prhs[3], L);
        options = getOptions(prhs[4]);

        result  = utility(K, counts, alpha, beta, gamma, options);

        freeCounts(counts);
        freeAlpha(alpha);
        free_vector(beta);
        free_matrix(gamma);
        free(options);

        return result;
}

static
void copyResult(utility_t* result, mxArray *plhs[]) {
        const char **fnames;
        const int nfields = 2;

        /* allocate memory  for storing pointers */
        fnames = mxCalloc(nfields, sizeof(*fnames));
        fnames[0] = "expectation";
        fnames[1] = "utility";

        plhs[0] = mxCreateStructMatrix(1, 1, nfields, fnames);
        mxFree((void *)fnames);

        if (result->expectation) {
                mxSetField(plhs[0], 0, "expectation", copyMatrixToMatlab(result->expectation));
        }
        if (result->utility) {
                mxSetField(plhs[0], 0, "utility", copyVectorToMatlab(result->utility));
        }
}

static
void freeResult(utility_t * result) {
        if (result->expectation) {
                free_matrix(result->expectation);
        }
        if (result->utility) {
                free_vector(result->utility);
        }
        free(result);
}

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
        utility_t * result;

        /* check proper input and output */
        if (nrhs != 5) {
                mexErrMsgTxt("Usage: samplingUtility(counts, alpha, beta, gamma, options).");
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

        result = callUtility(prhs);
        copyResult(result, plhs);
        freeResult(result);

        return;
}
