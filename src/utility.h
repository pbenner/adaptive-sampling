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

#ifndef UTILITY_H
#define UTILITY_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <datatypes.h>

#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/logarithmetic.h>
#include <adaptive-sampling/prombs.h>

////////////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////////////

static inline
prob_t sumModels(prob_t *ev_log, binProblem* bp)
{
        // sum up the vector returned by prombs
        prob_t sum = -HUGE_VAL;
        int i;

        for (i = 0; i < bp->bd->L; i++) {
                if (bp->bd->beta->content[i] > 0) {
                        sum = logadd(sum, ev_log[i]);
                }
        }

        return sum;
}

/* Find the smallest m_B for which the prior P(m_B) is nonzero. */
static inline
int minM(binProblem *bp)
{
        int i;
        for (i = bp->bd->L-1; i > 0; i--) {
                if (bp->bd->beta->content[i] > 0) {
                        return i;
                }
        }
        return i;
}

static inline
void binProblemInit(binProblem *bp, binData* bd)
{
        bp->bd              = bd;
        if (bd->options->algorithm == 0) {
                bp->ak      = alloc_matrix(bd->L, bd->L);
        }
        else {
                bp->ak      = NULL;
        }
        bp->bprob_pos       = -1;
        bp->counts_pos      = -1;
        bp->add_event.pos   = -1;
        bp->add_event.n     = 0;
        bp->add_event.which = bd->options->which;
        bp->fix_prob.pos    = -1;
        bp->fix_prob.val    = 0;
        bp->fix_prob.which  = bd->options->which;
}

static inline
void binProblemFree(binProblem *bp)
{
        if (bp->ak) {
                free_matrix(bp->ak);
        }
}

////////////////////////////////////////////////////////////////////////////////
// Count statistics
////////////////////////////////////////////////////////////////////////////////

static inline
size_t countStatistic(size_t event, int ks, int ke, binProblem *bp)
{
        if (bp != NULL && bp->add_event.which == event &&
            ks <= bp->add_event.pos && bp->add_event.pos <= ke) {
                return bp->bd->counts[event]->content[ks][ke] +
                        bp->add_event.n;
        }
        else if (ks <= ke) {
                return bp->bd->counts[event]->content[ks][ke];
        }
        else {
                return 0;
        }
}

static inline
prob_t countAlpha(size_t event, int ks, int ke, binProblem *bp)
{
        if (ks <= ke) {
                return bp->bd->alpha[event]->content[ks][ke];
        }
        else {
                return 0;
        }
}

////////////////////////////////////////////////////////////////////////////////
// Utility functions for Prombs
////////////////////////////////////////////////////////////////////////////////

static inline
void callBinningAlgorithm(
        prob_t (*f)(int, int, void*),
        prob_t *ev_log,
        binProblem *bp)
{
        switch (bp->bd->options->algorithm) {
        default:
        case 0:
                prombs(ev_log, bp->ak, bp->bd->prior_log, f, bp->bd->L, minM(bp), (void *)bp);
                break;
        case 1:
                prombs_tree(ev_log, bp->bd->prior_log, f, bp->bd->L, minM(bp), (void *)bp);
                break;
        case 2:
                mgs(ev_log, bp->bd->prior_log, f, bp->bd->L, (void *)bp);
                break;
        }
}

static inline
prob_t execPrombs_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        return iec_log(i, j, bp);
}

static inline
prob_t evidence(prob_t *ev_log, binProblem *bp)
{
        callBinningAlgorithm(execPrombs_f, ev_log, bp);

        return sumModels(ev_log, bp);
//        return prombs_rec(bp->bd->L, execPrombs_f, (void *)bp);
}

#endif /* UTILITY_H */
