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

#include <adaptive-sampling/exception.h>
#include <adaptive-sampling/linalg.h>

static Multibin** __multibins__;
static size_t __N__;
static gsl_rng* __r__;
static size_t* __counts__;

static
void sample_bin(
        size_t pos,
        prob_t *g,
        prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        prob_t post;
        prob_t sum1 = 0;
        prob_t sum2 = 0;
        list<bin_t>* bins1 = mb->get_bins(); mb->switch_break(pos);
        list<bin_t>* bins2 = mb->get_bins();
        prob_t r = (prob_t)rand()/RAND_MAX;

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
        post = expl(sum1 - logadd(sum1, sum2));
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
        prob_t *g,
        prob_t (*f)(int, int, void*),
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
        prob_t *g,
        prob_t (*f)(int, int, void*),
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

        __N__ = N;
        __multibins__    = (Multibin**)malloc((N+1)*sizeof(Multibin*));
        __multibins__[N] = (Multibin* )NULL;
        __counts__ = (size_t*)malloc(L*sizeof(size_t));

        size_t i;
        // initialize counts
        for (i = 0; i < L; i++) {
                __counts__[i] = 0;
        }

        // burn in
        __multibins__[0] = new Multibin(L);
        for (i = 0; i < R; i++) {
                sample_multibin(g, f, data, __multibins__[0]);
        }

        __counts__[__multibins__[0]->get_n_bins()-1]++;
        // sample
        for (i = 1; i < N; i++) {
                if (i%100 == 0) {
                        notice(NONE, "Generating samples... %.1f%%", (float)100*(i+1)/N);
                }
                __multibins__[i] = __multibins__[i-1]->copy();
                sample_multibin(g, f, data, __multibins__[i]);

                __counts__[__multibins__[i]->get_n_bins()-1]++;
        }
}

void mgs_free()
{
        size_t i;

        for (i = 0; __multibins__[i]; i++) {
                delete(__multibins__[i]);
        }
        free(__multibins__);
        free(__counts__);
        gsl_rng_free (__r__);
}

static
void evaluate(
        prob_t *result,
        prob_t *g,
        prob_t (*f)(int, int, void*),
        void *data,
        Multibin* mb)
{
        list<bin_t>* bins = mb->get_bins();
        if (g[bins->size()-1] > -HUGE_VAL) {
                prob_t sum = 0;
                for (list<bin_t>::iterator it = bins->begin(); it != bins->end(); it++) {
                        sum += (*f)((*it).from, (*it).to, data);
                }

                result[bins->size()-1] =
                        logadd(sum, result[bins->size()-1]);
        }
        delete(bins);
}

size_t *
mgs_get_counts()
{
        return __counts__;
}

void
mgs_get_bprob(vector_t *bprob, size_t L)
{
        size_t breaks[L];
        size_t i;

        for (i = 0; i < L; i++) {
                breaks[i] = 0;
        }
        for (i = 0; i < __N__; i++) {
                __multibins__[i]->get_breaks(breaks);
        }
        for (i = 0; i < L; i++) {
                bprob->content[i] = (prob_t)breaks[i]/__N__;
        }
}

void mgs(
        prob_t *result,
        prob_t *g,
        prob_t (*f)(int, int, void*),
        size_t L,
        void *data)
{
        size_t i;

        for (i = 0; i < L; i++) {
                result[i] = -HUGE_VAL;
        }

        // evaluate samples
        for (i = 0; __multibins__[i]; i++) {
                evaluate(result, g, f, data, __multibins__[i]);
        }
        for (i = 0; i < L; i++) {
                if (result[i] > -HUGE_VAL) {
                        result[i] -= LOG(__N__);
                }
        }
}

}
