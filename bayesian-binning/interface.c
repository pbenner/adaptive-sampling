#include <stdio.h>
#include <stdlib.h>

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

Matrix * binning(Vector *counts, unsigned int trials, Vector *prior)
{
        gsl_vector *tmp1 = toGslVector(counts);
        gsl_vector *tmp2 = toGslVector(prior);
        gsl_matrix *tmp3 = bin(tmp1, trials, tmp2);
        Matrix *m = fromGslMatrix(tmp3);

        gsl_vector_free(tmp1);
        gsl_vector_free(tmp2);
        gsl_matrix_free(tmp3);

        return m;
}

Matrix * test(Matrix *v)
{
        return v;
}
