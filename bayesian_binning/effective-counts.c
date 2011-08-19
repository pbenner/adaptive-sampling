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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <pthread.h>
#include <limits.h>
#include <sys/time.h>

#include <bayes/exception.h>
#include <bayes/logarithmetic.h>
#include <bayes/mgs.h>
#include <bayes/prombs.h>
#include <bayes/datatypes.h>
#include <bayes/uthash.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_psi.h>

#include <datatypes.h>
#include <model.h>
#include <utility.h>

////////////////////////////////////////////////////////////////////////////////
// Prombs effective counts functions
////////////////////////////////////////////////////////////////////////////////

static
prob_t effectiveCounts_f(int i, int j, void *data)
{
        binProblem *bp = (binProblem *)data;

        if (i <= bp->counts_pos && bp->counts_pos <= j) {
                int k;
                prob_t n = 0;
                for (k = 0; k < bd.events; k++) {
                        n += countStatistic(NULL, k, i, j) + countAlpha(k, i, j);
                }
                return log(n) + iec_log(NULL, i, j);
        }
        else {
                return iec_log(NULL, i, j);
        }
}
static
prob_t effectiveCounts(binProblem *bp, unsigned int pos, prob_t evidence_ref)
{
        prob_t ev_log[bd.L];

        bp->counts_pos = pos;
        prombs(ev_log, bp->ak, bd.prior_log, &effectiveCounts_f, bd.L, minM(), (void *)bp);

        return expl(sumModels(ev_log) - evidence_ref);
}

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////

void computeEffectiveCountsUtility(
        prob_t *result,
        prob_t evidence_ref)
{
        binProblem bp; binProblemInit(&bp);
        unsigned int i;

        for (i = 0; i < bd.L; i++) {
                notice(NONE, "Computing effective counts... %.1f%%", (float)100*(i+1)/bd.L);
                result[i] = -effectiveCounts(&bp, i, evidence_ref);
        }

        binProblemFree(&bp);
}
