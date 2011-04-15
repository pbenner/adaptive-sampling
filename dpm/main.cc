
using namespace std;

#include <iostream>
#include <cstdlib>
#include <ctime>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "init.hh"
#include "data.hh"
#include "cluster.hh"
#include "statistics.hh"
#include "gaussian-data.hh"
#include "gaussian-dpm.hh"

int main(void) {
        __dpm_init__();

        gsl_matrix* cov = gsl_matrix_alloc(2,2);
        gsl_matrix_set(cov, 0, 0, 0.5);
        gsl_matrix_set(cov, 1, 1, 0.5);
        gsl_matrix_set(cov, 0, 1, 0.2);
        gsl_matrix_set(cov, 1, 0, 0.2);

        GaussianData data(cov, 10, 4);
        GaussianDPM gdpm(data);

        for (int i = 0; i < 100000; i++) {
                gdpm.gibbsSample(10);
                cout << i << endl;
        }

        return 0;
}
