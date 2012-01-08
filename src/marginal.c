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
// Prombs marginal functions
////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////
// Loop through all X
////////////////////////////////////////////////////////////////////////////////

void computeMarginal(
        prob_t **marginals,
        prob_t evidence_ref,
        binData *bd)
{
        binProblem bp; binProblemInit(&bp, bd);
        int i, j;

        for (i = 0; i < bd->L; i++) {
                notice(NONE, "Computing marginals... %.1f%%", (float)100*(i+1)/bd->L);
                marginals[i][0] = 0;
                for (j = 1; j < bd->options->n_marginals; j++) {
                        prob_t p = j*bd->options->marginal_step;
                        if (bd->options->marginal_range.from <= p &&
                            bd->options->marginal_range.to   >= p) {
                                marginals[i][j] = marginal(i, p, bd->options->which, evidence_ref, &bp);
                        }
                        else {
                                marginals[i][j] = 0;
                        }
                }
        }

        binProblemFree(&bp);
}
