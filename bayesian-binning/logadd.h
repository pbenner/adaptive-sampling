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

#ifndef _LOGADD_H_
#define _LOGADD_H_

#include <math.h>

/* Log Sum of Exponentials Algorithm */

static inline
long double logadd(long double a, long double b)
{
        if (a < b) return a == -HUGE_VAL ? b : b + log1pl(expl(a-b));
        else       return b == -HUGE_VAL ? a : a + log1pl(expl(b-a));
}

#endif /* _LOGADD_H_ */
