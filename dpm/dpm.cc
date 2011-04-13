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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>

#include "cluster.hh"
#include "dpm.hh"
#include "statistics.hh"

using namespace std;

DPM::DPM(Data& data) : da(data), cl(da), alpha(1.0) {
}

DPM::~DPM() {
}

void DPM::sample(Data::element& element) {
        cl.release(element);
        Distribution& pred     = predictive();
        Distribution& postPred = posteriorPredictive();
        Cluster::size_type num_clusters = cl.num_clusters();
        double weights[num_clusters];
        Cluster::size_type labels[num_clusters];

        for (Cluster::iterator it = cl.begin(); it != cl.end(); it++) {
                
        }

        gsl_ran_discrete_t* gdd = gsl_ran_discrete_preproc(num_clusters, weights);
        gsl_ran_discrete(_r, gdd);
        gsl_ran_discrete_free(gdd);
}
