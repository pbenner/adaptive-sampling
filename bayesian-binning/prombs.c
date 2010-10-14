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

#include <exception.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>

#include "datatypes.h"
#include "logadd.h"

static
void successor(long double **ak, size_t L)
{
        size_t i, j;

        for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                        if (i < L-1 && j < L-1) {
                                ak[i][j] = ak[i+1][j+1];
                        }
                        else if (i == L-1 && j == L-1) {
                                ak[i][j] = logl(1);
                        }
                        else {
                                ak[i][j] = -HUGE_VAL;
                        }
                }
        }
}

static
void logproduct(long double **result, long double **ak, size_t L)
{
        long double tmp[L][L];
        size_t i, j, k;

        for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                        tmp[i][j] = result[i][j];
                }
        }

        for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                        result[i][j] = -HUGE_VAL;
                        for (k = 0; k < L; k++) {
                                result[i][j] = logadd(result[i][j], tmp[i][k] + ak[k][j]);
                        }
                }
        }        
}

static
long double ** allocMatrix(size_t L) {
        long double **a = (long double **)malloc(sizeof(long double *) * L);
        int i;

        for (i = 0; i < L; i++) {
                a[i] = (long double *)malloc(sizeof(long double) * L);
        }
        return a;
}

static
void freeMatrix(long double **a, size_t L) {
        int i;

        for (i = 0; i < L; i++) {
                free(a[i]);
        }
        free(a);
}

void prombs(long double *result, long double (*f)(int, int), size_t L, size_t m)
{
        long double **ak = allocMatrix(L);
        long double **pr = allocMatrix(L);
        size_t i, j;

        // initialise A^1 = (a^1_ij)_LxL
        for (i = 0; i < L; i++) {
                for (j = 0; j < L; j++) {
                        ak[i][j] = (i <= j ? (*f)(i, j) : -HUGE_VAL);
                        pr[i][j] = ak[i][j];
                }
        }

        for (i = 0; i < m; i++) {
                successor(ak, L);
                logproduct(pr, ak, L);
        }

        for (i = 0; i < L; i++) {
                result[L-1-i] = pr[0][i];
        }

        freeMatrix(ak, L);
        freeMatrix(pr, L);
}
