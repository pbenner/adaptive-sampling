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

#ifndef _LOGARITHMETIC_H_
#define _LOGARITHMETIC_H_

#include <math.h>

#include <bayes/datatypes.h>

/* Log Sum of Exponentials Algorithm */

static inline
prob_t logadd(prob_t a, prob_t b)
{
        if (a < b) return a == -HUGE_VAL ? b : b + log1pl(expl(a-b));
        else       return b == -HUGE_VAL ? a : a + log1pl(expl(b-a));
}

static inline
prob_t logsub(prob_t a, prob_t b)
{
        return b == -HUGE_VAL ? a : a + logl(1-expl(b-a));
}

#endif /* _LOGARITHMETIC_H_ */
