
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

struct vector {
        int size;
        double *vec;
};

struct matrix {
        int rows;
        int columns;
        double **mat;
};

typedef struct vector Vector;
typedef struct matrix Matrix;

extern void freeVector(Vector *v);
extern void freeMatrix(Matrix *m);
extern Vector * allocVector(int size);
extern Matrix * allocMatrix(int rows, int columns);

extern Matrix * binning(Vector *v);
extern Matrix * test(Matrix *v);
