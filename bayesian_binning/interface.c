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

#include <bayes_exception.h>
#include <bayes_linalg.h>

#include "interface.h"
#include "bayesian-binning.h"

// export functions for interface.py
Vector * _allocVector(int size)              { return allocVector(size); }
void     _freeVector(Vector *v)              { freeVector(v); }
Matrix * _allocMatrix(int rows, int columns) { return allocMatrix(rows, columns); }
void     _freeMatrix(Matrix *m)              { freeMatrix(m); }
void     _free(void *ptr)                    { free(ptr); }

BinningResult * binning(size_t events, Matrix **counts, Matrix **alpha, Vector *beta, Matrix *gamma, Options *options)
{
        size_t i;

        // allocate gsl matrices
        gsl_matrix **counts_m = (gsl_matrix **)malloc(events*sizeof(gsl_matrix *));
        gsl_matrix **alpha_m  = (gsl_matrix **)malloc(events*sizeof(gsl_matrix *));
        for (i = 0; i < events; i++) {
                counts_m[i] = toGslMatrix(counts[i]);
                alpha_m[i]  = toGslMatrix(alpha[i]);
        }
        gsl_vector *beta_v  = toGslVector(beta);
        gsl_matrix *gamma_m = toGslMatrix(gamma);
        BinningResultGSL *resultGsl;
        BinningResult    *result = (BinningResult *)malloc(sizeof(BinningResult));

        resultGsl = bin_log(events, counts_m, alpha_m, beta_v, gamma_m, options);
        if (options->n_moments > 0) {
                result->moments = fromGslMatrix(resultGsl->moments);
        }
        else {
                result->moments = NULL;
        }
        if (options->marginal) {
                result->marginals = fromGslMatrix(resultGsl->marginals);
        }
        else {
                result->marginals = NULL;
        }
        result->bprob   = fromGslVector(resultGsl->bprob);
        result->mpost   = fromGslVector(resultGsl->mpost);
        result->differential_gain = fromGslVector(resultGsl->differential_gain);
        result->effective_counts  = fromGslVector(resultGsl->effective_counts);
        result->multibin_entropy  = fromGslVector(resultGsl->multibin_entropy);

        // free gsl matrices
        for (i = 0; i < events; i++) {
                gsl_matrix_free(counts_m[i]);
                gsl_matrix_free(alpha_m[i]);
        }
        gsl_vector_free(beta_v);
        gsl_matrix_free(gamma_m);
        if (options->n_moments > 0) {
                gsl_matrix_free(resultGsl->moments);
        }
        if (options->marginal) {
                gsl_matrix_free(resultGsl->marginals);
        }
        gsl_vector_free(resultGsl->bprob);
        gsl_vector_free(resultGsl->mpost);
        gsl_vector_free(resultGsl->differential_gain);
        gsl_vector_free(resultGsl->effective_counts);
        gsl_vector_free(resultGsl->multibin_entropy);
        free(resultGsl);

        return result;
}
