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

static
mxArray* callPrombs(const mxArray *prhs[], size_t L, size_t m)
{
        matrix_t* f = alloc_matrix(L, L);
        prob_t g[L];
        prob_t result[L];

        copyArray (g, prhs[0], L);
        copyMatrix(f, prhs[1], 0);

        prombs(result, f, g, NULL, L, m, NULL);

        free_matrix(f);

        return copyArrayToMatlab(result, L);
}

static
void checkInput(int nlhs, mxArray *plhs[],
                int nrhs, const mxArray *prhs[])
{
        /* check proper input and output */
        if (nrhs != 2 && nrhs != 3) {
                mexErrMsgTxt("Usage: prombs(g, f) OR prombs(g, f, m).");
        }
        if (nlhs > 1) {
                mexErrMsgTxt("Too many output arguments.");
        }
        if (!mxIsClass(prhs[0], "double") || !mxIsClass(prhs[1], "double")) {
                mexErrMsgTxt("All matrices must be of type double.");
        }
        if (nrhs == 3 && ((mxGetM(prhs[2]) > 1) || (mxGetN(prhs[2]) > 1))) {
                mexErrMsgTxt("Third input argument is not a scalar.");
        }
        if (mxGetM(prhs[0]) > 1) {
                mexErrMsgTxt("First input is not a row vector.");
        }
        if (mxGetM(prhs[1]) != mxGetN(prhs[0]) || mxGetN(prhs[1]) != mxGetN(prhs[0])) {
                mexErrMsgTxt("Second input is not an LxL matrix.");
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
        if (nrhs == 3) {
                m = (size_t)*mxGetPr(prhs[2]) > 0 ? (size_t)*mxGetPr(prhs[2]) - 1 : 0;
                if (m > L-1) {
                        m = L-1;
                }
        }
        else {
                m = L-1;
        }

        plhs[0] = callPrombs(prhs, L, m);

        return;
}
