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
void logproduct(prob_t *result, prob_t **ak, size_t L, size_t i)
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
                                elem = ak[k][j];
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

static
prob_t ** allocMatrix(size_t L) {
        prob_t **a = (prob_t **)malloc(sizeof(prob_t *) * L);
        int i;

        for (i = 0; i < L; i++) {
                a[i] = (prob_t *)malloc(sizeof(prob_t) * L);
        }
        return a;
}

static
void freeMatrix(prob_t **a, size_t L) {
        int i;

        for (i = 0; i < L; i++) {
                free(a[i]);
        }
        free(a);
}

void prombs(prob_t *result, prob_t *g, prob_t (*f)(int, int), size_t L, size_t m)
{
        prob_t **ak = allocMatrix(L);
        prob_t pr[L];
        size_t i, j;

        // initialise A^1 = (a^1_ij)_LxL
        for (j = 0; j < L; j++) {
                ak[0][j] = (*f)(0, j);
                pr[j] = ak[0][j];
        }
        for (i = 1; i < L; i++) {
                for (j = i; j < L; j++) {
                        ak[i][j] = (*f)(i, j);
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

        freeMatrix(ak, L);
}

void prombsExt(
        prob_t *result,
        prob_t *g,
        prob_t (*f)(int, int),
        prob_t (*h)(int, int),
        prob_t epsilon,
        size_t L, size_t m)
{
        prob_t fprime(int i, int j) {
                return (*f)(i, j) + epsilon*(*h)(i, j);
        }
        size_t i;
        prob_t tmp[L];
        prombs(result, g, &fprime, L, m);
        prombs(tmp, g, f, L, m);

        for (i = 0; i < L; i++) {
                result[i] = logsub(result[i], tmp[i]) - logl(epsilon);
        }
}
