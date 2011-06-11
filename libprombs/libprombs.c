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

#include <bayes_datatypes.h>
#include <bayes_logarithmetic.h>

static prob_t prombsExt_epsilon = 0.0001;

void prombs_init(prob_t epsilon) {
        prombsExt_epsilon = epsilon;
}

static
void logproduct(prob_t *result, Matrix *ak, size_t L, size_t i)
{
        prob_t tmp[L], elem;
        size_t j, k;

        for (j = 0; j < L; j++) {
                tmp[j] = result[j];
        }

        for (j = i; j < L+i; j++) {
                result[j-i] = -HUGE_VAL;
                for (k = i; k <= j; k++) {
                        if (j < L && k < L) {
                                elem = ak->mat[k][j];
                        }
                        else if (k == j) {
                                elem = logl(1);
                        }
                        else  {
                                elem = -HUGE_VAL;
                        }
                        result[j-i] = logadd(result[j-i], tmp[k-i] + elem);
                }
        }
}

static inline
void init_f(Matrix *ak, prob_t (*f)(int, int, void*), size_t L, void *data)
{
        size_t i, j;

        // initialise A^1 = (a^1_ij)_LxL <- (f(i,j))_LxL
        for (j = 0; j < L; j++) {
                ak->mat[0][j] = (*f)(0, j, data);
        }
        // pr is now initialized to
        // pr = [f(0,1), f(0,2), ..., f(0,L)]
        for (i = 1; i < L; i++) {
                for (j = i; j < L; j++) {
                        ak->mat[i][j] = (*f)(i, j, data);
                }
        }
}

// result: array where the result is saved
// g: contains the prior P(m_B) for m_B = 1,...,L
// L: the number of inputs (maximal number of bins)
// m: the maximal number of bins in a multibin
void prombs(
        prob_t *result,
        Matrix *ak,
        prob_t *g,
        prob_t (*f)(int, int, void*),
        size_t L,
        size_t m,
        void *data)
{
        prob_t pr[L];
        size_t i, j;

        // init
        init_f(ak, f, L, data);
        for (j = 0; j < L; j++) {
                pr[j] = ak->mat[0][j];
        }

        // compute the products
        for (i = 0; i < m; i++) {
                logproduct(pr, ak, L, i+1);
        }
        // save result
        for (i = 0; i < L-m-1; i++) {
                // models with i>m were not computed, store a zero
                result[L-1-i] = -HUGE_VAL;
        }
        for (i = L-m-1; i < L; i++) {
                // the actual results are saved here
                if (g[L-1-i] == -HUGE_VAL) {
                        result[L-1-i] = -HUGE_VAL;
                }
                else {
                        result[L-1-i] = pr[i] + g[L-1-i];
                }
        }
}

static prob_t (*prombsExt_f)(int, int, void*);
static prob_t (*prombsExt_h)(int, int, void*);
static prob_t prombsExt_fprime(int i, int j, void *data) {
        return (*prombsExt_f)(i, j, data) + prombsExt_epsilon*(*prombsExt_h)(i, j, data);
}

void prombsExt(
        prob_t *result,
        Matrix *ak,
        prob_t *g,
        prob_t (*f)(int, int, void*), // on log scale
        prob_t (*h)(int, int, void*), // on normal scale
        size_t L,
        size_t m,
        void *data)
{
        size_t i;
        prob_t tmp[L];
        prombsExt_f = f;
        prombsExt_h = h;
        prombs(result, ak, g, &prombsExt_fprime, L, m, data);
        prombs(tmp, ak, g, f, L, m, data);

        for (i = 0; i < L; i++) {
                if (result[i] != tmp[i]) {
                        result[i] = logsub(result[i], tmp[i]) - logl(prombsExt_epsilon);
                }
                else {
                        // this can happen if all counts are zero
                        result[i] = -HUGE_VAL;
                }
        }
}
