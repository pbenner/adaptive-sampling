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

#include <config.h>

#include <stdlib.h>
#include <math.h>

#include "libmgs.hh"

extern "C" {

static
void sample_bin(
        size_t pos,
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        list<bin_t>* bins1 = mb->get_bins(); mb->switch_break(pos);
        list<bin_t>* bins2 = mb->get_bins();
        Bayes::prob_t r = (Bayes::prob_t)rand()/(RAND_MAX+1.0);

        Bayes::prob_t sum1 = 0;
        for (list<bin_t>::iterator it = bins1->begin(); it != bins1->end(); it++) {
                sum1 += (*f)((*it).from, (*it).to, data);
        }
        sum1 += g[bins1->size()-1];

        Bayes::prob_t sum2 = 0;
        for (list<bin_t>::iterator it = bins2->begin(); it != bins2->end(); it++) {
                sum2 += (*f)((*it).from, (*it).to, data);
        }
        sum2 += g[bins2->size()-1];

        Bayes::prob_t post = expl(sum1 - Bayes::logadd(sum1, sum2));

        if (r < post) {
                // use multibin 1
                mb->switch_break(pos);
                if (result) {
                        result[bins1->size()-1] =
                                Bayes::logadd(result[bins1->size()-1], sum1);
                }
        }
        else {
                // use multibin 2
                if (result) {
                        result[bins2->size()-1] =
                                Bayes::logadd(result[bins2->size()-1], sum2);
                }
        }
        free(bins1);
        free(bins2);
}

static
void sample_multibin(
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        size_t pos;

        for (pos = 0; pos < mb->get_n_breaks(); pos++) {
                sample_bin(pos, result, g, f, data, mb);
        }
}

void mgs(
        size_t N,
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        size_t L,
        void *data)
{
        Multibin *mb = new Multibin(L);
        size_t i;

        for (i = 0; i < L; i++) {
                result[i] = -HUGE_VAL;
        }

        // burn in
        for (i = 0; i < 100; i++) {
                sample_multibin(NULL, g, f, data, mb);
        }

        // sample
        for (i = 0; i < N; i++) {
                sample_multibin(result, g, f, data, mb);
        }

        free(mb);
}

}
