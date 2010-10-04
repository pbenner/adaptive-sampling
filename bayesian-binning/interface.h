
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "datatypes.h"

extern void freeVector(Vector *v);
extern void freeMatrix(Matrix *m);
extern Vector * allocVector(int size);
extern Matrix * allocMatrix(int rows, int columns);

extern Matrix * binning(Vector *counts, unsigned int trials, Vector *prior, Options *options);
extern Matrix * test(Matrix *v);
