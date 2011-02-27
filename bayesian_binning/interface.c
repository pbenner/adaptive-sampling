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

#include "interface.h"
#include "bayesian-binning.h"

Vector * allocVector(int size) {
        Vector *v  = (Vector *)malloc(sizeof(Vector));
        v->vec     = (double *)malloc(sizeof(double) * size);
        v->size    = size;
        return v;
}

Matrix * allocMatrix(int rows, int columns) {
        Matrix *m  = (Matrix * )malloc(sizeof(Matrix));
        m->mat     = (double **)malloc(sizeof(double *) * rows);
        m->rows    = rows;
        m->columns = columns;
        int i;
        for (i = 0; i < rows; i++) {
                m->mat[i] = (double *)malloc(sizeof(double) * columns);
        }
        return m;
}

void freeVector(Vector *v) {
        free(v->vec);
        free(v);
}

void freeMatrix(Matrix *m) {
        int i;
        for (i = 0; i < m->rows; i++) {
                free(m->mat[i]);
        }
        free(m->mat);
        free(m);
}

gsl_vector * toGslVector(const Vector *vector)
{
        gsl_vector *v = gsl_vector_alloc(vector->size);
        int i;

        for(i = 0; i < vector->size; i++) {
                gsl_vector_set(v, i, vector->vec[i]);
        }
        return v;
}

Vector * fromGslVector(const gsl_vector * vector)
{
        int i;
        int size    = vector->size;
        Vector *v   = allocVector(size);

        for(i = 0; i < size; i++) {
                v->vec[i] = gsl_vector_get(vector, i);
        }
        return v;
}

gsl_matrix * toGslMatrix(const Matrix *matrix)
{
        gsl_matrix *m = gsl_matrix_alloc(matrix->rows, matrix->columns);
        int i, j;

        for(i = 0; i < matrix->rows; i++) {
                for(j = 0; j < matrix->columns; j++) {
                        gsl_matrix_set(m, i, j, matrix->mat[i][j]);
                }
        }
        return m;
}

Matrix * fromGslMatrix(const gsl_matrix * matrix)
{
        int i, j;
        int rows    = matrix->size1;
        int columns = matrix->size2;
        Matrix *m   = allocMatrix(rows, columns);

        for(i = 0; i < rows; i++) {
                for(j = 0; j < columns; j++) {
                        m->mat[i][j] = gsl_matrix_get(matrix, i, j);
                }
        }
        return m;
}

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
