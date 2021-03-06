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

#ifndef _PROMBS_H_
#define _PROMBS_H_

#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/linalg.h>
#include <adaptive-sampling/probtype.h>

void __init_prombs__(prob_t epsilon);
void prombs(prob_t *result, prob_t **ak, prob_t *g, prob_t (*f)(int, int, void*), size_t L, size_t m, void *data);
prob_t prombs_rec(
        size_t j,
        prob_t (*f)(int, int, void*),
        void *data);
void prombsExt(
        prob_t *result,
        prob_t **ak,
        prob_t *g,
        prob_t (*f)(int, int, void*),
        prob_t (*h)(int, int, void*),
        size_t L, size_t m, void *data);

#endif /* _PROMBS_H_ */
