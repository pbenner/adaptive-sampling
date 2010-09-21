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
}

Matrix * binning(Vector *v)
{
        gsl_matrix *tmp = bin(v);
        Matrix *m = fromGslMatrix(tmp);

        gsl_matrix_free(tmp);

        return m;
}

Matrix * test(Matrix *v)
{
        return v;
}
