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
// Model
////////////////////////////////////////////////////////////////////////////////

typedef struct {
        double key;
        double value;
        UT_hash_handle hh;
} lngamma_hash_t;

static lngamma_hash_t*  lngamma_map = NULL;
static pthread_rwlock_t lngamma_map_lock;

static
void lngamma_map_free() {
        lngamma_hash_t *current, *tmp;

        HASH_ITER(hh, lngamma_map, current, tmp) {
                HASH_DEL(lngamma_map, current);
                free(current);
        }
}

static
double hashed_lngamma(double p)
{
        lngamma_hash_t* s;
        if (pthread_rwlock_rdlock(&lngamma_map_lock) != 0) {
                fprintf(stderr, "Can't get lngamma_map_lock (r)\n");
                exit(EXIT_FAILURE);
        }
        HASH_FIND(hh, lngamma_map, &p, sizeof(double), s);
        pthread_rwlock_unlock(&lngamma_map_lock);

        if (s == NULL) {
                lngamma_hash_t* new = (lngamma_hash_t*)malloc(sizeof(lngamma_hash_t));
                new->key   = p;
                new->value = gsl_sf_lngamma(p);
                if (pthread_rwlock_wrlock(&lngamma_map_lock) != 0) {
                        fprintf(stderr, "Can't get lngamma_map_lock (w)\n");
                        exit(EXIT_FAILURE);
                }
                HASH_ADD(hh, lngamma_map, key, sizeof(double), new);
                pthread_rwlock_unlock(&lngamma_map_lock);
                return new->value;
        }
        else {
                return s->value;
        }
}

prob_t mbeta_log(prob_t *p)
{
        unsigned int i;
        prob_t sum1, sum2;

        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < bd.events; i++) {
                sum1 += p[i];
                sum2 += gsl_sf_lngamma(p[i]);
//                sum2 += hashed_lngamma(p[i]);
        }

        return sum2 - gsl_sf_lngamma(sum1);
//        return sum2 - hashed_lngamma(sum1);
}

/* P(E|B) */
prob_t iec_log(binProblem *bp, int kk, int k)
{
        unsigned int i;
        prob_t c[bd.events];
        prob_t alpha[bd.events];
        prob_t gamma = gsl_matrix_get(bd.gamma, kk, k);
        if (gamma == 0) {
                return -HUGE_VAL;
        }
        for (i = 0; i < bd.events; i++) {
                c[i]     = countStatistic(bp, i, kk, k) + countAlpha(i, kk, k);
                alpha[i] = countAlpha(i, kk, k);
        }
        if (bp != NULL && kk <= bp->fix_prob.pos && bp->fix_prob.pos <= k) {
                // compute marginals
                // TODO: extend to multinomial case
                if (bp->fix_prob.which == 0) {
                        return logl(gamma) + (c[0]-1)*logl(bp->fix_prob.val)
                                + (c[1]-1)*logl(1-bp->fix_prob.val)
                                - mbeta_log(alpha);
                }
                else {
                        return logl(gamma) + (c[0]-1)*log(1-bp->fix_prob.val)
                                + (c[1]-1)*log(bp->fix_prob.val)
                                - mbeta_log(alpha);
                }
        }
        else {
                return logl(gamma) + (mbeta_log(c) - mbeta_log(alpha));
        }
}

////////////////////////////////////////////////////////////////////////////////
// Model init
////////////////////////////////////////////////////////////////////////////////

void init_model() {
        if (pthread_rwlock_init(&lngamma_map_lock, NULL) != 0) {
                fprintf(stderr, "Can't create lngamma_map_lock\n");
                exit(EXIT_FAILURE);
        }
}

void free_model() {
        lngamma_map_free();
}
