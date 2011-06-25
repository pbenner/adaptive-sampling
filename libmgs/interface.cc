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

#include <sys/time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

#include "libmgs.hh"

extern "C" {

#include <bayes/exception.h>

static Multibin** __multibins__;
static gsl_rng * __r__;

static
void sample_bin(
        size_t pos,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        Bayes::prob_t post;
        Bayes::prob_t sum1 = 0;
        Bayes::prob_t sum2 = 0;
        list<bin_t>* bins1 = mb->get_bins(); mb->switch_break(pos);
        list<bin_t>* bins2 = mb->get_bins();
        Bayes::prob_t r = (Bayes::prob_t)rand()/RAND_MAX;

        if (g[bins1->size()-1] > -HUGE_VAL) {
                for (list<bin_t>::iterator it = bins1->begin(); it != bins1->end(); it++) {
                        sum1 += (*f)((*it).from, (*it).to, data);
                }
                sum1 += g[bins1->size()-1];
        }
        else {
                goto err;
        }

        if (g[bins2->size()-1] > -HUGE_VAL) {
                for (list<bin_t>::iterator it = bins2->begin(); it != bins2->end(); it++) {
                        sum2 += (*f)((*it).from, (*it).to, data);
                }
                sum2 += g[bins2->size()-1];
        }
        else {
                mb->switch_break(pos);
                goto err;
        }

        // sample
        post = expl(sum1 - Bayes::logadd(sum1, sum2));
        if (r < post) {
                // use multibin 1
                mb->switch_break(pos);
        }
        // else:
        // use multibin 2
err:
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
        size_t N = mb->get_n_breaks();
        size_t pos;
        gsl_permutation * p = gsl_permutation_alloc(N);
        gsl_permutation_init(p);
        gsl_ran_shuffle (__r__, p->data, N, sizeof(size_t));

        for (pos = 0; pos < N; pos++) {
                sample_bin(gsl_permutation_get(p, pos), g, f, data, mb);
        }
        gsl_permutation_free(p);
}

void mgs_init(
        size_t R,
        size_t N,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        size_t L,
        void *data)
{
        const gsl_rng_type * T = gsl_rng_default;
        gsl_rng_env_setup();
        __r__ = gsl_rng_alloc (T);

        struct timeval tv;
        gettimeofday(&tv, NULL);
        time_t seed = tv.tv_sec*tv.tv_usec;
        srand(seed);


        __multibins__    = (Multibin**)malloc((N+1)*sizeof(Multibin*));
        __multibins__[N] = (Multibin* )NULL;

        size_t i;

        // burn in
        __multibins__[0] = new Multibin(L);
        for (i = 0; i < R; i++) {
                sample_multibin(g, f, data, __multibins__[0]);
        }

        // sample
        for (i = 1; i < N; i++) {
                if (i%100 == 0) {
                        notice(NONE, "Generating samples... %.1f%%", (float)100*(i+1)/N);
                }
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
        gsl_rng_free (__r__);
}

static
void evaluate(
        Bayes::prob_t *result,
        Bayes::prob_t *g,
        Bayes::prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        list<bin_t>* bins = mb->get_bins();
        if (g[bins->size()-1] > -HUGE_VAL) {
                Bayes::prob_t sum = 0;
                for (list<bin_t>::iterator it = bins->begin(); it != bins->end(); it++) {
                        sum += (*f)((*it).from, (*it).to, data);
                }
//                sum += g[bins->size()-1];

                result[bins->size()-1] =
                        Bayes::logadd(result[bins->size()-1], sum);
        }
        delete(bins);
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
                evaluate(result, g, f, data, __multibins__[N]);
        }
        for (i = 0; i < L; i++) {
                if (result[i] > -HUGE_VAL) {
                        result[i] -= logl(N);
                }
        }
}

}
