#include <stdio.h>
#include <stdlib.h>

#include "example.h"

int cube(int n)
{
        return n * n * n;
}

gsl_matrix * toGslMatrix(Matrix *matrix)
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

Matrix * test(Matrix *v)
{
        return v;
} 

double product(Vector *l) {
        int i;
        if (l->size <= 0) {
                return 0;
        }
        double ret_val = 1;
        for (i = 0; i < l->size; i++) {
                ret_val *= l->vec[i];
        }
        return ret_val;
}
