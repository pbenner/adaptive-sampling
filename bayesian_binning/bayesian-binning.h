/* Copyright (C) 2010 Philipp Benner
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

#ifndef BAYESIAN_BINNING_H
#define BAYESIAN_BINNING_H

#include <config.h>

#include <gsl/gsl_matrix.h>

extern void __init__(prob_t epsilon);
extern void __free__();

extern prob_t bin_entropy(size_t events, gsl_matrix **counts, gsl_matrix **alpha, gsl_vector *beta, gsl_matrix *gamma, Options *options);
extern BinningResultGSL * bin_log(size_t events, gsl_matrix **counts, gsl_matrix **alpha, gsl_vector *beta, gsl_matrix *gamma, Options *options);

#endif /* BAYESIAN_BINNING */
