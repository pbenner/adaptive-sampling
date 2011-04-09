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

void prombs(prob_t *result, prob_t *g, prob_t (*f)(int, int), size_t L, size_t m)
{
        Matrix *ak = allocMatrix(L, L);
        prob_t pr[L];
        size_t i, j;

        // initialise A^1 = (a^1_ij)_LxL
        for (j = 0; j < L; j++) {
                ak->mat[0][j] = (*f)(0, j);
                pr[j] = ak->mat[0][j];
        }
        for (i = 1; i < L; i++) {
                for (j = i; j < L; j++) {
                        ak->mat[i][j] = (*f)(i, j);
                }
        }
        // compute the products
        for (i = 0; i < m; i++) {
                logproduct(pr, ak, L, i+1);
        }
        // save result
        for (i = 0; i < L; i++) {
                result[L-1-i] = pr[i] + g[L-1-i];
        }

        freeMatrix(ak);
}

static prob_t prombsExt_epsilon;
static prob_t (*prombsExt_f)(int, int);
static prob_t (*prombsExt_h)(int, int);
static prob_t prombsExt_fprime(int i, int j) {
        return (*prombsExt_f)(i, j) + prombsExt_epsilon*(*prombsExt_h)(i, j);
}

void prombsExt(
        prob_t *result,
        prob_t *g,
        prob_t (*f)(int, int), // on log scale
        prob_t (*h)(int, int), // on normal scale
        prob_t epsilon,
        size_t L, size_t m)
{
        size_t i;
        prob_t tmp[L];
        prombsExt_f = f;
        prombsExt_h = h;
        prombsExt_epsilon = epsilon;
        prombs(result, g, &prombsExt_fprime, L, m);
        prombs(tmp, g, f, L, m);

        for (i = 0; i < L; i++) {
                result[i] = logsub(result[i], tmp[i]) - logl(epsilon);
        }
}
