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

#include <datatypes.h>

////////////////////////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////////////////////////

void computeModelPosteriors(
        prob_t *ev_log,
        prob_t *mpost,
        prob_t evidence_ref,
        binData *bd)
{
        unsigned int j;

        if (bd->options->algorithm == 2) {
                size_t *counts = mgs_get_counts();
                for (j = 0; j < bd->L; j++) {
                        if (bd->beta->content[j] > 0) {
                                mpost[j] = (prob_t)counts[j]/bd->options->samples[1];
                        }
                        else {
                                mpost[j] = 0;
                        }
                }
        }
        else {
                for (j = 0; j < bd->L; j++) {
                        if (bd->beta->content[j] > 0) {
                                mpost[j] = expl(ev_log[j] - evidence_ref);
                        }
                        else {
                                mpost[j] = 0;
                        }
                }
        }
}
