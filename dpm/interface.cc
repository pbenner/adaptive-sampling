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

#include <gsl/gsl_matrix.h>

namespace Bayes {
        extern "C" {
#include <bayes_linalg.h>
        }
}

#include "init.hh"
#include "cluster.hh"
#include "gaussian-data.hh"
#include "gaussian-dpm.hh"
#include "interface.hh"

static GaussianDPM* _gdpm;

extern "C" {

Bayes::Vector * _allocVector(int size)              { return Bayes::allocVector(size); }
void            _freeVector(Bayes::Vector *v)       { Bayes::freeVector(v); }
Bayes::Matrix * _allocMatrix(int rows, int columns) { return Bayes::allocMatrix(rows, columns); }
void            _freeMatrix(Bayes::Matrix *m)       { Bayes::freeMatrix(m); }
void            _free(void *ptr)                    { free(ptr); }

void _dpm_init(unsigned int n, unsigned int k) {
        __dpm_init__();

        gsl_matrix* cov = gsl_matrix_alloc(2,2);
        gsl_matrix_set(cov, 0, 0, 0.5);
        gsl_matrix_set(cov, 1, 1, 0.5);
        gsl_matrix_set(cov, 0, 1, 0.2);
        gsl_matrix_set(cov, 1, 0, 0.2);

        GaussianData data(cov, n, k);
        _gdpm = new GaussianDPM(data);
}

unsigned int _dpm_num_clusters() {
        return _gdpm->num_clusters();
}

Bayes::Matrix* _dpm_cluster(unsigned int c) {
        Cluster::cluster& cl = (*_gdpm)[c];
        int n = cl.elements.size();
        int m = 2;

        Bayes::Matrix* result = Bayes::allocMatrix(n, m);
        int i = 0;
        for (Cluster::elements_t::iterator it = cl.elements.begin();
             it != cl.elements.end(); it++) {
                result->mat[i][0] = (*it)->x[0];
                result->mat[i][1] = (*it)->x[1];
                i++;
        }

        return result;
}

Bayes::Vector* _dpm_original_tags(unsigned int c) {
        Cluster::cluster& cl = (*_gdpm)[c];
        int n = cl.elements.size();

        Bayes::Vector* result = Bayes::allocVector(n);
        int i = 0;
        for (Cluster::elements_t::iterator it = cl.elements.begin();
             it != cl.elements.end(); it++) {
                result->vec[i] = (double)((*it)->original_cluster);
                i++;
        }

        return result;
}

Bayes::Matrix* _dpm_hist_means() {
        vector<Data::x_t>& mean = _gdpm->get_hist_means();
        unsigned int len = mean.size();
        Bayes::Matrix* result = Bayes::allocMatrix(len, 2);

        for (unsigned int i = 0; i < len; i++) {
                result->mat[i][0] = mean[i][0];
                result->mat[i][1] = mean[i][1];
        }

        return result;
}

Bayes::Matrix* _dpm_means() {
        vector<Data::x_t>* means = _gdpm->cluster_means();
        unsigned int len = means->size();
        Bayes::Matrix* result = Bayes::allocMatrix(len, 2);

        for (unsigned int i = 0; i < len; i++) {
                result->mat[i][0] = (*means)[i][0];
                result->mat[i][1] = (*means)[i][1];
        }
        delete means;

        return result;
}

void _dpm_print() {
        cout << *_gdpm << endl;
}

void _dpm_sample(unsigned int n) {
        _gdpm->gibbsSample(n);
}

void _dpm_free() {
        delete(_gdpm);
}

}
