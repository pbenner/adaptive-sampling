
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

gsl_matrix * bin(gsl_vector *v)
{
        gsl_matrix *m = gsl_matrix_alloc(2,2);

        gsl_matrix_set(m, 0, 0, 1);
        gsl_matrix_set(m, 0, 1, 2);
        gsl_matrix_set(m, 1, 0, 3);
        gsl_matrix_set(m, 1, 1, 4);

        return m;
}
