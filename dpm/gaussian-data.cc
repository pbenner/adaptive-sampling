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

#include <gsl/gsl_matrix.h>

#include "gaussian-data.hh"
#include "statistics.hh"

GaussianData::GaussianData(gsl_matrix* cov, int _n, int _k) {
        gsl_vector* mu = gsl_vector_alloc(2);
        double sample_x, sample_y;
        int tag = 0;
        n = _k * _n;

        // for every cluster
        for (int i = 0; i < _k; i++) {
                // generate a new mean
                gsl_vector_set(mu, 0, 20.0*(double)rand()/RAND_MAX);
                gsl_vector_set(mu, 0, 20.0*(double)rand()/RAND_MAX);
                BivariateNormal bg(cov, mu);
                // generate n samples
                for (int j = 0; j < _n; j++) {
                        bg.sample(&sample_x, &sample_y);
                        vector<double> v;
                        v.push_back(sample_x);
                        v.push_back(sample_y);
                        Data::element e = { v, tag++, i };
                        elements.push_back(e);
                }
        }
}

GaussianData::~GaussianData() {

}
