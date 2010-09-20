
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

extern int cube( int n );
extern Matrix * test( Matrix *v);
extern double product(Vector *l);
