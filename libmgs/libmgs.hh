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

#ifndef LIBMGS_HH
#define LIBMGS_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <iostream>
#include <cstdlib>
#include <list>
#include <vector>

namespace Bayes {
        extern "C" {
#include <bayes/datatypes.h>
#include <bayes/logarithmetic.h>
        }
}

using namespace std;

typedef struct {
        size_t from;
        size_t to;
} bin_t;

class Multibin {
public:
         Multibin(size_t L);
         Multibin(Multibin* mb, uint32_t* breaks);
        ~Multibin();

        Multibin* copy();
        void insert_break(size_t i);
        void remove_break(size_t i);
        void switch_break(size_t i);

        void print();
        void print_bins();

        list<bin_t>* get_bins();
        size_t get_n_breaks();
        size_t get_n_bins();

private:
        size_t n;        // length of breaks[]
        size_t n_breaks; // number of possible breaks
        size_t n_bins;   // current number of bins
        uint32_t* breaks;
};


#endif /* LIBMGS_HH */
