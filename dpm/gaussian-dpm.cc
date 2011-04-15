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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "gaussian-dpm.hh"

using namespace std;

GaussianDPM::GaussianDPM(
        GaussianData* data,
        gsl_matrix* _cov,
        gsl_matrix* _cov_0,
        gsl_vector* _mu_0)
        : DPM(data)
{
        cov            = gsl_matrix_alloc(2,2);
        cov_0          = gsl_matrix_alloc(2,2);
        cov_inv        = gsl_matrix_alloc(2,2);
        cov_inv_0      = gsl_matrix_alloc(2,2);
        predictive_cov = gsl_matrix_alloc(2,2);
        cov_n          = gsl_matrix_alloc(2,2);
        mu_0           = gsl_vector_alloc(2);
        mu_n           = gsl_vector_alloc(2);

        _mean          = gsl_vector_alloc(2);

        // likelihood
        gsl_matrix_memcpy(cov, _cov);

        // prior
        gsl_matrix_memcpy(cov_0, _cov_0);
        gsl_vector_memcpy(mu_0,  _mu_0);

        // predictive distribution
        gsl_matrix_memcpy(predictive_cov, cov);
        gsl_matrix_add(predictive_cov, cov_0);

        // inverse covariance matrices
        __inv_tmp  = gsl_matrix_alloc(2, 2);
        __inv_perm = gsl_permutation_alloc(2);

        inverse(cov,   cov_inv);
        inverse(cov_0, cov_inv_0);

        // initialize distributions
        predictiveDist          = new BivariateNormal(predictive_cov, mu_0);
        posteriorPredictiveDist = new BivariateNormal();
}

GaussianDPM::~GaussianDPM() {
        gsl_matrix_free(cov);
        gsl_matrix_free(cov_0);
        gsl_matrix_free(cov_inv);
        gsl_matrix_free(cov_inv_0);
        gsl_matrix_free(predictive_cov);
        gsl_matrix_free(cov_n);
        gsl_vector_free(mu_0);
        gsl_vector_free(mu_n);

        gsl_vector_free(_mean);

        gsl_matrix_free(__inv_tmp);
        gsl_permutation_free(__inv_perm);

        delete(predictiveDist);
        delete(posteriorPredictiveDist);
        delete(da);
}

void GaussianDPM::inverse(gsl_matrix* src, gsl_matrix* dst) {
        __inv_s = 0;
        gsl_matrix_memcpy(__inv_tmp, src);
        gsl_linalg_LU_decomp(__inv_tmp, __inv_perm, &__inv_s);
        gsl_linalg_LU_invert(__inv_tmp, __inv_perm, dst);
}

void GaussianDPM::_computeMean(const Cluster::cluster& cluster) {
        double num = cluster.elements.size();
        // set _mean to zero
        gsl_vector_set(_mean, 0, 0);
        gsl_vector_set(_mean, 1, 0);

        // compute cluster mean
        for (Cluster::elements_t::const_iterator it  = cluster.elements.begin();
             it != cluster.elements.end(); it++) {
                gsl_vector_set(_mean, 0, gsl_vector_get(_mean, 0) + (*it)->x[0]);
                gsl_vector_set(_mean, 1, gsl_vector_get(_mean, 1) + (*it)->x[1]);
        }
        gsl_vector_set(_mean, 0, gsl_vector_get(_mean, 0)/num);
        gsl_vector_set(_mean, 1, gsl_vector_get(_mean, 1)/num);
}

Distribution& GaussianDPM::posteriorPredictive(const Cluster::cluster& cluster) {
        double num = cluster.elements.size();

        _computeMean(cluster);

        // posterior covariance
        gsl_matrix_memcpy(cov_n, cov_inv);
        gsl_matrix_scale(cov_n, num);
        gsl_matrix_add(cov_n, cov_inv_0);
        inverse(cov_n, cov_n);

        // posterior mean
        gsl_vector* tmp1 = gsl_vector_alloc(2);
        gsl_vector* tmp2 = gsl_vector_alloc(2);

        gsl_blas_dgemv(CblasNoTrans, 1.0, cov_inv,  _mean, 0.0, tmp1);
        gsl_blas_dgemv(CblasNoTrans, 1.0, cov_inv_0, mu_0, 0.0, tmp2);
        gsl_vector_scale(tmp1, num);
        gsl_vector_add(tmp1, tmp2);
        gsl_blas_dgemv(CblasTrans, 1.0, cov_n, tmp1, 0.0, mu_n);

        gsl_vector_free(tmp1);
        gsl_vector_free(tmp2);

        gsl_matrix_add(cov_n, cov);
        posteriorPredictiveDist->update(cov_n, mu_n);

        return *posteriorPredictiveDist;
}

Distribution& GaussianDPM::predictive() {
        return *predictiveDist;
}

double GaussianDPM::likelihood() {
        double likelihood = 0;
        for (Cluster::iterator it = cl.begin(); it != cl.end(); it++) {
                _computeMean(**it);
                double cluster_size = (*it)->elements.size();
                BivariateNormal bg(cov, _mean);
                for (Cluster::elements_t::iterator is  = (*it)->elements.begin();
                     is != (*it)->elements.end(); is++) {
                        likelihood += bg.pdf((*is)->x)/cluster_size;
                }
        }
        return likelihood/cl.size();
}

vector<Data::x_t>* GaussianDPM::cluster_means() {
        vector<Data::x_t>* means = new vector<Data::x_t>();

        for (Cluster::iterator it = cl.begin(); it != cl.end(); it++) {
                _computeMean(**it);
                Data::x_t mean;
                mean.push_back(gsl_vector_get(_mean, 0));
                mean.push_back(gsl_vector_get(_mean, 1));
                means->push_back(mean);
        }
        return means;
}

void GaussianDPM::compute_statistics() {
        hist_likelihood.push_back(likelihood());
        hist_num_clusters.push_back(cl.size());
        for (Cluster::iterator it = cl.begin(); it != cl.end(); it++) {
                _computeMean(**it);
                Data::x_t mean;
                mean.push_back(gsl_vector_get(_mean, 0));
                mean.push_back(gsl_vector_get(_mean, 1));
                hist_means.push_back(mean);
        }
}
