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

#include <adaptive-sampling/prombs.h>

static __inline__
prob_t** alloc_prombs_matrix(size_t L) {
        size_t i;
        prob_t** m = (prob_t **)malloc(L*sizeof(prob_t*));
        for (i = 0; i < L; i++) {
                m[i] = (prob_t *)malloc(L*sizeof(prob_t));
        }
        return m;
}

static __inline__
void free_prombs_matrix(prob_t** m, size_t L) {
        size_t i;
        for (i = 0; i < L; i++) {
                free(m[i]);
        }
        free(m);
}

typedef struct _data_t_ {
        const mxArray* f;
        const mxArray* h;
        mwSize  nsubs;
        mwIndex* subs;
} data_t;

prob_t prombs_f(int i, int j, void* _data)
{
        data_t* data = (data_t *)_data;

        data->subs[0] = i;
        data->subs[1] = j;
        return mxGetPr(data->f)[mxCalcSingleSubscript(data->f, data->nsubs, data->subs)];
}

prob_t prombs_h(int i, int j, void* _data)
{
        data_t* data = (data_t *)_data;

        data->subs[0] = i;
        data->subs[1] = j;
        return mxGetPr(data->h)[mxCalcSingleSubscript(data->h, data->nsubs, data->subs)];
}

static
mxArray* callPrombs(const mxArray *prhs[], size_t L, size_t m)
{
        prob_t** ak = alloc_prombs_matrix(L);
        prob_t g[L];
        prob_t result[L];
        mwSize nsubs = mxGetNumberOfDimensions(prhs[1]);
        data_t data  = { prhs[1], prhs[2], nsubs, mxCalloc(nsubs, sizeof(mwIndex)) };

        copyArray (g, prhs[0], L);

        prombsExt(result, ak, g, prombs_f, prombs_h, L, m, &data);

        free_prombs_matrix(ak, L);

        return copyArrayToMatlab(result, L);
}

static
void checkInput(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
        /* check proper input and output */
        if (nrhs != 4 && nrhs != 5) {
                mexErrMsgTxt("Usage: prombsExtended(g, f, h, epsilon) OR prombs(g, f, h, epsilon, m).");
        }
        if (nlhs > 1) {
                mexErrMsgTxt("Too many output arguments.");
        }
        if (!mxIsClass(prhs[0], "double") || !mxIsClass(prhs[1], "double") || !mxIsClass(prhs[2], "double")) {
                mexErrMsgTxt("All vectors and matrices must be of type double.");
        }
        if (nrhs == 5 && ((mxGetM(prhs[4]) > 1) || (mxGetN(prhs[4]) > 1))) {
                mexErrMsgTxt("Fifth input argument is not a scalar.");
        }
        if (mxGetM(prhs[0]) > 1) {
                mexErrMsgTxt("First input is not a row vector.");
        }
        if (mxGetM(prhs[1]) != mxGetN(prhs[0]) || mxGetN(prhs[1]) != mxGetN(prhs[0])) {
                mexErrMsgTxt("Second input is not an LxL matrix.");
        }
        if (mxGetM(prhs[2]) != mxGetN(prhs[0]) || mxGetN(prhs[2]) != mxGetN(prhs[0])) {
                mexErrMsgTxt("Third input is not an LxL matrix.");
        }
        if (*mxGetPr(prhs[3]) == 0.0 || *mxGetPr(prhs[3]) > 1.0) {
                mexErrMsgTxt("Invalid epsilon.");
        }
}

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
        size_t L;
        size_t m;

        checkInput(nlhs, plhs, nrhs, prhs);

        L = mxGetN(prhs[0]);

        /* if we have four arguments the obtain m */
        if (nrhs == 5) {
                m = (size_t)*mxGetPr(prhs[4]) > 0 ? (size_t)*mxGetPr(prhs[4]) - 1 : 0;
                if (m > L-1) {
                        m = L-1;
                }
        }
        else {
                m = L-1;
        }
        __init_prombs__(*mxGetPr(prhs[3]));

        plhs[0] = callPrombs(prhs, L, m);

        return;
}
