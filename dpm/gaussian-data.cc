/* Copyright (C) 2011 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

#include "gaussian-data.hh"
#include "statistics.hh"

GaussianData::GaussianData(
        gsl_matrix* cov,
        gsl_matrix* cov_0,
        gsl_vector* mu_0,
        int n, double* pi, size_t k)
        : Data()
{
        gsl_vector* mu = gsl_vector_alloc(2);
        double sample_x, sample_y;
        int tag = 0;

        BivariateNormal bg_0(cov_0, mu_0);
        Data::x_t _mu;
        _mu.push_back(0);
        _mu.push_back(0);
        // for every cluster
        for (size_t i = 0; i < k; i++) {
                // generate a new mean
                bg_0.sample(&_mu[0], &_mu[1]);
                // save mean
                means.push_back(_mu);
        }

        gsl_ran_discrete_t* gdd = gsl_ran_discrete_preproc(k, pi);
        // generate n samples
        for (int j = 0; j < n; j++) {
                int i = gsl_ran_discrete(_r, gdd);
                gsl_vector_set(mu, 0, means[i][0]);
                gsl_vector_set(mu, 1, means[i][1]);
                BivariateNormal bg(cov, mu);
                bg.sample(&sample_x, &sample_y);
                vector<double> v;
                v.push_back(sample_x);
                v.push_back(sample_y);
                Data::element e = { v, tag++, i };
                elements.push_back(e);
        }
        gsl_ran_discrete_free(gdd);
}

GaussianData::~GaussianData() {

}
