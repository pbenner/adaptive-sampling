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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

#include <adaptive-sampling/exception.h>
#include <adaptive-sampling/logarithmetic.h>
#include <adaptive-sampling/mgs.h>
#include <adaptive-sampling/prombs.h>
#include <adaptive-sampling/datatypes.h>

#include <datatypes.h>
#include <model.h>
#include <threading.h>
#include <tools.h>

/******************************************************************************
 * Prombs moment functions
 ******************************************************************************/

static
prob_t moment(
        int nth,
        int pos,
        int which,
        prob_t evidence_ref,
        binProblem *bp)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bp->bd->L];

        bp->add_event.pos   = pos;
        bp->add_event.n     = nth;
        bp->add_event.which = which;
        evidence_log        = evidence(evidence_log_tmp, bp);
        bp->add_event.pos   = -1;
        bp->add_event.n     = 0;

        return EXP(evidence_log - evidence_ref);
}

/******************************************************************************
 * HMM moment function
 ******************************************************************************/

static
prob_t hmm_he(int from, int to, binProblem* bp)
{
        size_t i;
        prob_t ca[bp->bd->events];
        prob_t cb[bp->bd->events];

        for (i = 0; i < bp->bd->events; i++) {
                ca[i] = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
                cb[i] = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
        }
        ca[bp->add_event.which] += bp->add_event.n;

        return mbeta_log(ca, bp) - mbeta_log(cb, bp);
}

void hmm_computeMoments(
        matrix_t *moments,
        vector_t *forward,
        vector_t *backward,
        binProblem *bp)
{
        size_t i, j;
        vector_t* tmp = alloc_vector(bp->bd->L);

        for (i = 0; i < bp->bd->options->n_moments; i++) {
                bp->add_event.n     = i+1;
                bp->add_event.which = bp->bd->options->which;
                hmm_fb(tmp, forward, backward, &hmm_he, bp);
                bp->add_event.pos   = -1;
                bp->add_event.n     = 0;

                for (j = 0; j < bp->bd->L; j++) {
                        moments->content[i][j] = exp(tmp->content[j]);
                }
        }
        free_vector(tmp);
}

/******************************************************************************
 * Main
 ******************************************************************************/

static
void * computeMoments_thread(void* data_)
{
        pthread_data_t *data  = (pthread_data_t *)data_;
        binProblem *bp = data->bp;
        int i = data->i, j;
        matrix_t *moments = (matrix_t *)data->result;
        prob_t evidence_ref = data->evidence_ref;

        /* Moments */
        for (j = 0; j < bp->bd->options->n_moments; j++) {
                moments->content[j][i] = moment(j+1, i, bp->bd->options->which, evidence_ref, bp);
        }
        return NULL;
}

void computeMoments(
        matrix_t *moments,
        prob_t evidence_ref,
        binData *bd)
{
        threaded_computation((void *)moments, evidence_ref, bd, computeMoments_thread,
                             "Computing moments: %.1f%%");
}
