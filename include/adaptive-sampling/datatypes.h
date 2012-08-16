/* Copyright (C) 2010 Philipp Benner
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

#ifndef ADAPTIVE_SAMPLING_DATATYPES_H
#define ADAPTIVE_SAMPLING_DATATYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

/******************************************************************************
 * Data structures
 ******************************************************************************/

#include <adaptive-sampling/linalg.h>

typedef struct _options_ {
        float epsilon;
        int verbose;
        int prombsTest;
        /* break probabilities */
        int bprob;
        int threads;
        int stacksize;
        /* specify components of the utility function */
        int kl_psi;
        int kl_multibin;
        int effective_counts;
        int effective_posterior_counts;
        /* which event */
        int which;
        /* binning algorithm */
        int algorithm;
        /* burnin length and number of samples */
        int samples[2];
        /* compute marginal */
        int density;
        float density_step;
        struct {
                float from;
                float to;
        } density_range;
        int n_moments;
        int n_density;
        int model_posterior;
        /* use hidden markov model instead of prombs */
        int hmm;
        /* hmm parameters */
        float rho;
} options_t;

typedef struct _marginal_ {
        matrix_t *moments;
        matrix_t *density;
        vector_t *bprob;
        vector_t *mpost;
} marginal_t;

typedef struct _utility_ {
        matrix_t *expectation;
        vector_t *utility;
} utility_t;

#endif /* ADAPTIVE_SAMPLING_DATATYPES_H */
