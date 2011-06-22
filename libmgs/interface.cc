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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <gsl/gsl_permutation.h>

#include "libmgs.hh"

extern "C" {

static Multibin** __multibins__;

static
void evaluate(
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        list<bin_t>* bins = mb->get_bins();
        Bayes::prob_t sum = 0;
        for (list<bin_t>::iterator it = bins->begin(); it != bins->end(); it++) {
                sum += (*f)((*it).from, (*it).to, data);
        }
        sum += g[bins->size()-1];

        result[bins->size()-1] =
                Bayes::logadd(result[bins->size()-1], sum);

        delete(bins);
}

static
void sample_bin(
        size_t pos,
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
                printf("Hello1: %lu:%lu\n", (long unsigned int)(*it).from, (long unsigned int)(*it).to);
                sum1 += (*f)((*it).from, (*it).to, data);
                printf("Hello2\n");
        }
        sum1 += g[bins1->size()-1];

        Bayes::prob_t sum2 = 0;
        for (list<bin_t>::iterator it = bins2->begin(); it != bins2->end(); it++) {
                sum2 += (*f)((*it).from, (*it).to, data);
        }
        sum2 += g[bins2->size()-1];

        // sample
        Bayes::prob_t post = expl(sum1 - Bayes::logadd(sum1, sum2));
        if (r < post) {
                // use multibin 1
                mb->switch_break(pos);
        }
        // else:
        // use multibin 2

        delete(bins1);
        delete(bins2);
}

static
void sample_multibin(
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        size_t pos;

        for (pos = 0; pos < mb->get_n_breaks(); pos++) {
                sample_bin(pos, g, f, data, mb);
        }
}

void mgs(
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        size_t L,
        void *data)
{
        size_t N;
        size_t i;

        for (i = 0; i < L; i++) {
                result[i] = -HUGE_VAL;
        }

        // evaluate samples
        for (N = 0; __multibins__[N]; N++) {
                evaluate(result, g, f, data, __multibins__[i]);
        }
        for (i = 0; i < L; i++) {
                if (result[i] != -HUGE_VAL) {
                        result[i] -= logl(N);
                }
        }
}

void mgs_init(
        size_t N,
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        size_t L,
        void *data)
{
        __multibins__    = (Multibin**)malloc((N+1)*sizeof(Multibin*));
        __multibins__[N] = (Multibin*)NULL;

        L = 16;
        printf("L: %u\n", (unsigned int)L);
        printf("N: %u\n", (unsigned int)N);

        size_t i;

        // burn in
        __multibins__[0] = new Multibin(L);
        for (i = 0; i < 200; i++) {
                sample_multibin(g, f, data, __multibins__[0]);
        }

        // sample
        for (i = 1; i < N; i++) {
                __multibins__[i] = __multibins__[i-1]->copy();
                sample_multibin(g, f, data, __multibins__[i]);
        }
}

void mgs_free()
{
        size_t i;

        for (i = 0; __multibins__[i]; i++) {
                delete(__multibins__[i]);
        }
        free(__multibins__);
}

}
