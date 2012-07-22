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

#ifndef LOCAL_UTILITY_H
#define LOCAL_UTILITY_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <adaptive-sampling/datatypes.h>

void computeKLUtility(vector_t *result, prob_t evidence_ref, binData *bd);
void hmm_computeUtility(
        vector_t *utility,
        vector_t *forward,
        vector_t *backward,
        binProblem *bp);
void hmm_computeNStepUtility(
        vector_t *utility,
        binProblem *bp);

#endif /* LOCAL_UTILITY_H */
