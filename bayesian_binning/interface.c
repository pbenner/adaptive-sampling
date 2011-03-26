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
                
BinningResult * binning(Matrix *counts, Matrix *alpha, Vector *prior, Options *options)
{
        gsl_matrix *counts_m = toGslMatrix(counts);
        gsl_matrix *alpha_m  = toGslMatrix(alpha);
        gsl_vector *prior_v  = toGslVector(prior);
        BinningResultGSL *resultGsl;
        BinningResult    *result = (BinningResult *)malloc(sizeof(BinningResult));

        resultGsl = bin_log(counts_m, alpha_m, prior_v, options);
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

        gsl_matrix_free(counts_m);
        gsl_matrix_free(alpha_m);
        gsl_vector_free(prior_v);
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
