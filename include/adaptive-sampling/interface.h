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

#ifndef ADAPTIVE_SAMPLING_INTERFACE_H
#define ADAPTIVE_SAMPLING_INTERFACE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <adaptive-sampling/datatypes.h>

void __init__(double epsilon);
void __free__();

marginal_t* posterior(
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t *beta,
        matrix_t *gamma,
        options_t *options);
utility_t* utility(
        int events,
        matrix_t **counts,
        matrix_t **alpha,
        vector_t  *beta,
        matrix_t  *gamma,
        options_t *options);

#endif /* ADAPTIVE_SAMPLING_INTERFACE */
