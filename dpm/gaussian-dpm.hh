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

#ifndef GAUSSIAN_DPM_HH
#define GAUSSIAN_DPM_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>

#include "data.hh"
#include "dpm.hh"
#include "cluster.hh"
#include "statistics.hh"
#include "gaussian-data.hh"

using namespace std;

class GaussianDPM : public DPM {
public:
        GaussianDPM(GaussianData& data,
                    gsl_matrix* _cov, gsl_matrix* _cov_0, gsl_vector* _mu_0);
        ~GaussianDPM();

        Distribution& posteriorPredictive(const Cluster::cluster& cluster);
        Distribution& predictive();

        void inverse(gsl_matrix* src, gsl_matrix* dst);
        double likelihood();
        void compute_statistics();
        vector<Data::x_t>* cluster_means();

        vector<Data::x_t>* get_original_means() {
                return da.get_means();
        }

        vector<Data::x_t>& get_hist_means() {
                return hist_means;
        }

private:
        GaussianData da;

        // likelihood parameters
        gsl_matrix* cov;
        gsl_matrix* cov_inv;

        // prior parameters
        gsl_vector* mu_0;
        gsl_matrix* cov_0;
        gsl_matrix* cov_inv_0;

        // predictive distribution
        gsl_matrix* predictive_cov;

        // posterior predictive distribution
        gsl_matrix* cov_n;
        gsl_vector* mu_n;
        gsl_vector* _mean;

        // inverse matrices
        gsl_matrix*      __inv_tmp;
        gsl_permutation* __inv_perm;
        int              __inv_s;

        // means
        vector<Data::x_t> hist_means;

        // private methods
        void _computeMean(const Cluster::cluster& cluster);
};

#endif /* GAUSSIAN_DPM_HH */
