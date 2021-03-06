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
#include <tools.h>

/******************************************************************************
 * HMM density functions
 ******************************************************************************/

static
prob_t hmm_hd(int from, int to, binProblem* bp)
{
        size_t i;
        prob_t counts;
        prob_t alpha [bp->bd->events];
        prob_t result = 0;

        for (i = 0; i < bp->bd->events; i++) {
                alpha[i] = countAlpha(i, from, from, bp); 
        }
        result -= mbeta_log(alpha, bp);
        for (i = 0; i < bp->bd->events; i++) {
                counts = countAlpha(i, from, from, bp) + countStatistic(i, from, to, bp);
                if (i == bp->fix_prob.which) {
                        result += (counts-1)*LOG(bp->fix_prob.val);
                }
                else {
                        result += (counts-1)*LOG(1-bp->fix_prob.val);
                }
        }

        return result;
}

void hmm_computeDensity(
       matrix_t *density,
       prob_t *forward,
       prob_t *backward,
       binProblem *bp)
{
        size_t i, j;
        prob_t tmp[bp->bd->L];

        for (j = 0; j < bp->bd->options->n_density; j++) {
                prob_t p = j*bp->bd->options->density_step;
                if (bp->bd->options->density_range.from <= p &&
                    bp->bd->options->density_range.to   >= p &&
                    p != 0.0 && p != 1.0) {

                        bp->fix_prob.pos   = -1;
                        bp->fix_prob.val   = p;
                        bp->fix_prob.which = bp->bd->options->which;

                        hmm_fb(tmp, forward, backward, &hmm_hd, bp);
                        for (i = 0; i < bp->bd->L; i++) {
                                density->content[i][j] = EXP(tmp[i]);
                        }

                        bp->fix_prob.pos   = -1;
                        bp->fix_prob.val   =  0;
                        bp->fix_prob.which =  0;
                }
                else {
                        for (i = 0; i < bp->bd->L; i++) {
                                density->content[i][j] = 0;
                        }
                }
        }
}

/******************************************************************************
 * Prombs density functions
 ******************************************************************************/

static
prob_t density(
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
void * computeDensity_thread(void* data_)
{
        pthread_data_t *data  = (pthread_data_t *)data_;
        binProblem *bp = data->bp;
        int i = data->i, j;
        matrix_t *result = (matrix_t *)data->result;
        prob_t evidence_ref = data->evidence_ref;

        /* Moments */
        for (j = 0; j < bp->bd->options->n_density; j++) {
                prob_t p = j*bp->bd->options->density_step;
                if (bp->bd->options->density_range.from <= p &&
                    bp->bd->options->density_range.to   >= p &&
                    p != 0.0 && p != 1.0) {
                        result->content[i][j] = density(i, p, bp->bd->options->which, evidence_ref, bp);
                }
                else {
                        result->content[i][j] = 0;
                }
        }
        return NULL;
}

void computeDensity(
        matrix_t *result,
        prob_t evidence_ref,
        binData *bd)
{
        threaded_computation((void *)result, evidence_ref, bd, computeDensity_thread,
                             "Computing density: %.1f%%");
}
