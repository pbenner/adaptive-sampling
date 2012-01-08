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

typedef struct {
        binProblem *bp;
        int i;
        vector_t *bprob;
        prob_t evidence_ref;
} pthread_data_bprob;

static
void * computeBreakProbabilities_thread(void* data_)
{
        pthread_data_bprob *data  = (pthread_data_bprob *)data_;
        binProblem *bp = data->bp;
        int i = data->i;
        vector_t *bprob = data->bprob;
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

        size_t i, j, rc;
        binProblem bp[bd->options->threads];
        pthread_t threads[bd->options->threads];
        pthread_data_bprob data[bd->options->threads];
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        if (bd->options->stacksize < PTHREAD_STACK_MIN) {
                if (pthread_attr_setstacksize (&attr, PTHREAD_STACK_MIN) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }
        else {
                if (pthread_attr_setstacksize (&attr, (size_t)bd->options->stacksize) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }

        for (j = 0; j < bd->options->threads && j < bd->L; j++) {
                binProblemInit(&bp[j], bd);
                data[j].bp = &bp[j];
                data[j].bprob = bprob;
                data[j].evidence_ref = evidence_ref;
        }
        for (i = 0; i < bd->L; i += bd->options->threads) {
                for (j = 0; j < bd->options->threads && i+j < bd->L; j++) {
                        notice(NONE, "Computing break probabilities: %.1f%%", (float)100*(i+j+1)/bd->L);
                        data[j].i = i+j;
                        rc = pthread_create(&threads[j], &attr, computeBreakProbabilities_thread, (void *)&data[j]);
                        if (rc) {
                                std_err(NONE, "Couldn't create thread.");
                        }
                }
                for (j = 0; j < bd->options->threads && i+j < bd->L; j++) {
                        rc = pthread_join(threads[j], NULL);
                        if (rc) {
                                std_err(NONE, "Couldn't join thread.");
                        }
                }
        }
        for (j = 0; j < bd->options->threads && j < bd->L; j++) {
                binProblemFree(&bp[j]);
        }
}
