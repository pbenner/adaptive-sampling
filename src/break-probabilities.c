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
#include <pthread.h>
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

////////////////////////////////////////////////////////////////////////////////
// Prombs entropy functions
////////////////////////////////////////////////////////////////////////////////

static
prob_t breakProb_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        // include only those bins that don't cover
        // position pos, which means, that only multi-bins
        // are included, that have a break at position
        // pos
        if (i < bp->bprob_pos && bp->bprob_pos <= j) {
                return -HUGE_VAL;
        }
        return iec_log(i, j, bp);
}
static
prob_t breakProb(size_t pos, prob_t evidence_ref, binProblem *bp)
{
        prob_t ev_log[bp->bd->L];

        bp->bprob_pos = pos;
        callBinningAlgorithm(&breakProb_f, ev_log, bp);

        return EXP(sumModels(ev_log, bp) - evidence_ref);
}

static
void * computeBreakProbabilities_thread(void* data_)
{
        pthread_data_t *data  = (pthread_data_t *)data_;
        binProblem *bp = data->bp;
        int i = data->i;
        vector_t *bprob = (vector_t *)data->result;
        prob_t evidence_ref = data->evidence_ref;

        bprob->content[i] = breakProb(i, evidence_ref, bp);
        return NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////

void computeBreakProbabilities(
        vector_t *bprob,
        prob_t evidence_ref,
        binData *bd)
{
        if (bd->options->algorithm == 2) {
                mgs_get_bprob(bprob, bd->L);
                return;
        }
        threaded_computation((void *)bprob, evidence_ref, bd, computeBreakProbabilities_thread);
}
