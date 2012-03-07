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

typedef struct options {
        float epsilon;
        int verbose;
        int prombsTest;
        /* break probabilities */
        int bprob;
        int threads;
        int stacksize;
        /* compute sampling utility */
        int utility;
        /* specify components of the utility function */
        int kl_component;
        int kl_multibin;
        int effective_counts;
        /* which event */
        int which;
        /* binning algorithm */
        int algorithm;
        /* burnin length and number of samples */
        int samples[2];
        /* compute marginal */
        int marginal;
        float marginal_step;
        struct {
                float from;
                float to;
        } marginal_range;
        int n_moments;
        int n_marginals;
        int model_posterior;
} Options;

typedef struct binningResult{
        matrix_t *moments;
        matrix_t *marginals;
        vector_t *bprob;
        vector_t *mpost;
        vector_t *utility;
} BinningResult;

#endif /* ADAPTIVE_SAMPLING_DATATYPES_H */
