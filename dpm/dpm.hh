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

#ifndef DPM_HH
#define DPM_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

using namespace std;

#include "data.hh"
#include "cluster.hh"
#include "statistics.hh"

class DPM {
public:
        // constructors and destructors
        DPM(Data& data);
        ~DPM();

        // operators
        friend ostream& operator<< (std::ostream& o, DPM const& dpm);

        // methods
        bool sample(Data::element& element);
        void gibbsSample(unsigned int steps);

        virtual Distribution& posteriorPredictive(const Cluster::cluster& cluster) {
                Distribution* dist = new Distribution;
                return *dist;
        }
        virtual Distribution& predictive() {
                Distribution* dist = new Distribution;
                return *dist;
        }

protected:
        Data da;
        Cluster cl;
        float alpha;
};

#endif /* DPM_HH */
