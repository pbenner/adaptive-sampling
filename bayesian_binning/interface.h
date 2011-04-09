
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include <bayes_datatypes.h>

extern void freeVector(Vector *v);
extern void freeMatrix(Matrix *m);
extern Vector * allocVector(int size);
extern Matrix * allocMatrix(int rows, int columns);

extern Vector * fromGslVector(const gsl_vector * vector);
extern Matrix * fromGslMatrix(const gsl_matrix * matrix);
extern gsl_vector * toGslVector(const Vector *vector);
extern gsl_matrix * toGslMatrix(const Matrix *matrix);

extern BinningResult * binning(Matrix *counts, Matrix *alpha, Vector *prior, Options *options);
