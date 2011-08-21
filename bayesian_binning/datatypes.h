/* Copyright (C) 2010, 2011 Philipp Benner
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

#ifndef DATATYPES_H
#define DATATYPES_H

#include <config.h>

#include <bayes/linalg.h>
#include <bayes/datatypes.h>

#include <gsl/gsl_matrix.h>

////////////////////////////////////////////////////////////////////////////////
// Data structures
////////////////////////////////////////////////////////////////////////////////

typedef struct options {
        float epsilon;
        int verbose;
        int prombsTest;
        // break probabilities
        int bprob;
        int threads;
        int stacksize;
        // compute sampling utility
        int utility;
        // specify components of the utility function
        int differential_entropy;
        int multibin_entropy;
        int effective_counts;
        // which event
        int which;
        // binning algorithm
        int algorithm;
        // burnin length and number of samples
        int samples[2];
        // compute marginal
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

typedef struct binningResultGSL {
        gsl_matrix *moments;
        gsl_matrix *marginals;
        gsl_vector *bprob;
        gsl_vector *mpost;
        gsl_vector *utility;
} BinningResultGSL;

typedef struct binningResult{
        Matrix *moments;
        Matrix *marginals;
        Vector *bprob;
        Vector *mpost;
        Vector *utility;
} BinningResult;

// data that has to be immutable
typedef struct {
        Options *options;
        // number of timesteps
        unsigned int L;
        unsigned int events;
        prob_t *prior_log;     // P(p,B|m_B)
        // counts and parameters
        gsl_matrix **counts;
        gsl_matrix **alpha;
        gsl_vector  *beta;     // P(m_B)
        gsl_matrix  *gamma;
} binData;

// mutable data, local to each thread
typedef struct {
        binData* bd;
        // temporary memory for prombs
        Matrix *ak;
        // break probability
        int bprob_pos;
        // effective counts
        int counts_pos;
        // moments
        struct {
                int pos;
                int n;
                int which;
        } add_event;
        // marginals
        struct {
                int pos;
                prob_t val;
                int which;
        } fix_prob;
} binProblem;

#endif /* DATATYPES_H */
