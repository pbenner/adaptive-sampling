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
#include <adaptive-sampling/uthash.h>

#include <datatypes.h>
#include <model.h>
#include <threading.h>
#include <utility.h>

/******************************************************************************
 * Prombs marginal functions
 ******************************************************************************/

static
prob_t marginal(
        int pos,
        prob_t val,
        int which,
        prob_t evidence_ref,
        binProblem *bp)
{
        prob_t evidence_log;
        prob_t evidence_log_tmp[bp->bd->L];

        bp->fix_prob.pos   = pos;
        bp->fix_prob.val   = val;
        bp->fix_prob.which = which;
        evidence_log       = evidence(evidence_log_tmp, bp);
        bp->fix_prob.pos   = -1;
        bp->fix_prob.val   =  0;
        bp->fix_prob.which =  0;

        return EXP(evidence_log - evidence_ref);
}

/******************************************************************************
 * Loop through all X
 ******************************************************************************/

static
void * computeMarginal_thread(void* data_)
{
        pthread_data_t *data  = (pthread_data_t *)data_;
        binProblem *bp = data->bp;
        int i = data->i, j;
        matrix_t *marginals = (matrix_t *)data->result;
        prob_t evidence_ref = data->evidence_ref;

        /* Moments */
        for (j = 0; j < bp->bd->options->n_marginals; j++) {
                prob_t p = j*bp->bd->options->marginal_step;
                if (bp->bd->options->marginal_range.from <= p &&
                    bp->bd->options->marginal_range.to   >= p &&
                    p != 0.0 && p != 1.0) {
                        marginals->content[i][j] = marginal(i, p, bp->bd->options->which, evidence_ref, bp);
                }
                else {
                        marginals->content[i][j] = 0;
                }
        }
        return NULL;
}

void computeMarginal(
        matrix_t *marginals,
        prob_t evidence_ref,
        binData *bd)
{
        threaded_computation((void *)marginals, evidence_ref, bd, computeMarginal_thread,
                             "Computing marginal: %.1f%%");
}
