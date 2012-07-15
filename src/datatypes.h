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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/probtype.h>

/******************************************************************************
 * Data structures
 ******************************************************************************/

/* data that has to be immutable */
typedef struct {
        Options *options;
        /* number of timesteps */
        size_t L;
        size_t events;
        prob_t *prior_log;     /* P(p,B|m_B) */
        /* counts and parameters */
        matrix_t **counts;
        matrix_t **alpha;
        vector_t  *beta;       /* P(m_B) */
        matrix_t  *gamma;
        /* hmm parameters */
        prob_t rho;
} binData;

/* mutable data, local to each thread */
typedef struct {
        binData* bd;
        /* temporary memory for prombs */
        matrix_t *ak;
        /* break probability */
        int bprob_pos;
        /* effective counts */
        int counts_pos;
        /* moments */
        struct {
                int pos;
                int n;
                int which;
        } add_event;
        /* marginals */
        struct {
                int pos;
                prob_t val;
                int which;
        } fix_prob;
} binProblem;

#endif /* DATATYPES_H */
