
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

typedef struct vector {
        int size;
        double *vec;
} Vector;

typedef struct matrix {
        int rows;
        int columns;
        double **mat;
} Matrix;

extern void freeVector(Vector *v);
extern void freeMatrix(Matrix *m);
extern Vector * allocVector(int size);
extern Matrix * allocMatrix(int rows, int columns);

extern Matrix * binning(Vector *counts, unsigned int trials);
extern Matrix * test(Matrix *v);
