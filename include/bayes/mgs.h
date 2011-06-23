/* Copyright (C) 2011 Philipp Benner
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

#ifndef _MGS_H_
#define _MGS_H_

#include <bayes/datatypes.h>

extern void mgs(prob_t *result, prob_t *g, prob_t (*f)(int, int, void*), size_t L, void *data);
extern void mgs_init(size_t R, size_t N, prob_t *g, prob_t (*f)(int, int, void*), size_t L, void *data);
extern void mgs_free();

#endif /* _MGS_H_ */
