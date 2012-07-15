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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <limits.h>
#include <sys/time.h>

#include <adaptive-sampling/exception.h>
#include <adaptive-sampling/logarithmetic.h>
#include <adaptive-sampling/mgs.h>
#include <adaptive-sampling/prombs.h>
#include <adaptive-sampling/datatypes.h>
#include <adaptive-sampling/uthash.h>

#include <gsl/gsl_sf_gamma.h>

#include <datatypes.h>
#include <model.h>
#include <tools.h>

#ifdef HAVE_LIB_PTHREAD
#include <pthread.h>
#endif /* HAVE_LIB_PTHREAD */

/******************************************************************************
 * Model
 ******************************************************************************/

typedef struct {
        double key;
        double value;
        UT_hash_handle hh;
} lngamma_hash_t;

static lngamma_hash_t*  lngamma_map = NULL;
#ifdef HAVE_LIB_PTHREAD
static pthread_rwlock_t lngamma_map_lock;
#endif /* HAVE_LIB_PTHREAD */

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
#ifdef HAVE_LIB_PTHREAD
        if (pthread_rwlock_rdlock(&lngamma_map_lock) != 0) {
                fprintf(stderr, "Can't get lngamma_map_lock (r)\n");
                exit(EXIT_FAILURE);
        }
#endif /* HAVE_LIB_PTHREAD */
        HASH_FIND(hh, lngamma_map, &p, sizeof(double), s);
#ifdef HAVE_LIB_PTHREAD
        pthread_rwlock_unlock(&lngamma_map_lock);
#endif /* HAVE_LIB_PTHREAD */

        if (s == NULL) {
                lngamma_hash_t* new = (lngamma_hash_t*)malloc(sizeof(lngamma_hash_t));
                new->key   = p;
                new->value = gsl_sf_lngamma(p);
#ifdef HAVE_LIB_PTHREAD
                if (pthread_rwlock_wrlock(&lngamma_map_lock) != 0) {
                        fprintf(stderr, "Can't get lngamma_map_lock (w)\n");
                        exit(EXIT_FAILURE);
                }
#endif /* HAVE_LIB_PTHREAD */
                HASH_ADD(hh, lngamma_map, key, sizeof(double), new);
#ifdef HAVE_LIB_PTHREAD
                pthread_rwlock_unlock(&lngamma_map_lock);
#endif /* HAVE_LIB_PTHREAD */
                return new->value;
        }
        else {
                return s->value;
        }
}

prob_t mbeta_log(prob_t *p, binProblem *bp)
{
        size_t i;
        prob_t sum1, sum2;

        sum1 = 0;
        sum2 = 0;
        for (i = 0; i < bp->bd->events; i++) {
                sum1 += p[i];
                sum2 += gsl_sf_lngamma(p[i]);
/*                sum2 += hashed_lngamma(p[i]); */
        }

        return sum2 - gsl_sf_lngamma(sum1);
/*        return sum2 - hashed_lngamma(sum1); */
}

/* P(E|B) */
prob_t iec_log(int kk, int k, binProblem *bp)
{
        size_t i;
        prob_t c[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t gamma = bp->bd->gamma->content[kk][k];
        if (gamma == 0) {
                return -HUGE_VAL;
        }
        for (i = 0; i < bp->bd->events; i++) {
                c[i]     = countStatistic(i, kk, k, bp) + countAlpha(i, kk, k, bp);
                alpha[i] = countAlpha(i, kk, k, bp);
        }
        if (bp != NULL && kk <= bp->fix_prob.pos && bp->fix_prob.pos <= k) {
                /* compute marginals
                 * TODO: extend to multinomial case */
                if (bp->fix_prob.which == 0) {
                        return LOG(gamma) + (c[0]-1)*LOG(bp->fix_prob.val)
                                + (c[1]-1)*LOG(1-bp->fix_prob.val)
                                - mbeta_log(alpha, bp);
                }
                else {
                        return LOG(gamma) + (c[0]-1)*LOG(1-bp->fix_prob.val)
                                + (c[1]-1)*LOG(bp->fix_prob.val)
                                - mbeta_log(alpha, bp);
                }
        }
        else {
                return LOG(gamma) + (mbeta_log(c, bp) - mbeta_log(alpha, bp));
        }
}

/******************************************************************************
 * Hidden Markov model
 ******************************************************************************/

static
prob_t hmm_hp(int from, int to, binProblem* bp)
{
        return 0.0;
}

static
prob_t hmm_forward_rec(vector_t *result, size_t j, size_t to, prob_t (*f)(int, int, binProblem*), binProblem* bp)
{
        size_t i, k;

        prob_t c[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t tmp;

        for (i = 0; i < bp->bd->events; i++) {
                c[i]     = countAlpha(i, 0, 0, bp) + countStatistic(i, 0, to, bp);
                alpha[i] = countAlpha(i, 0, 0, bp);
        }
        tmp = to*log(bp->bd->rho) + mbeta_log(c, bp) - mbeta_log(alpha, bp) + f(0, to, bp);
        for (k = 0; k < j; k++) {
                for (i = 0; i < bp->bd->events; i++) {
                        c[i]     = countAlpha(i, k+1, k+1, bp) + countStatistic(i, k+1, to, bp);
                        alpha[i] = countAlpha(i, k+1, k+1, bp);
                }
                tmp = logadd(tmp, (to-k-1)*log(bp->bd->rho) + log(1.0-bp->bd->rho) + result->content[k] + mbeta_log(c, bp) - mbeta_log(alpha, bp) + f(k+1, to, bp));
        }
        return tmp;
}

void hmm_forward(vector_t *result, binProblem* bp)
{
        size_t j;

        for (j = 0; j < bp->bd->L; j++) {
                result->content[j] = hmm_forward_rec(result, j, j, &hmm_hp, bp);
        }
}

static
prob_t hmm_backward_rec(vector_t *result, size_t j, prob_t (*f)(int, int, binProblem*), binProblem* bp)
{
        size_t i, k;

        prob_t c[bp->bd->events];
        prob_t alpha[bp->bd->events];
        prob_t tmp;

        for (i = 0; i < bp->bd->events; i++) {
                c[i]     = countAlpha(i, j, j, bp) + countStatistic(i, j, bp->bd->L-1, bp);
                alpha[i] = countAlpha(i, j, j, bp);
        }
        tmp = (bp->bd->L-j-1)*log(bp->bd->rho) + mbeta_log(c, bp) - mbeta_log(alpha, bp) + f(j, bp->bd->L-1, bp);
        for (k = j+1; k < bp->bd->L; k++) {
                for (i = 0; i < bp->bd->events; i++) {
                        c[i]     = countAlpha(i, j, j, bp) + countStatistic(i, j, k-1, bp);
                        alpha[i] = countAlpha(i, j, j, bp);
                }
                tmp = logadd(tmp, (k-1-j)*log(bp->bd->rho) + log(1.0-bp->bd->rho) + result->content[k] + mbeta_log(c, bp) - mbeta_log(alpha, bp) + f(j, k-1, bp));
        }
        return tmp;
}

void hmm_backward(vector_t *result, binProblem* bp)
{
        size_t j;

        for (j = 0; j < bp->bd->L; j++) {
                result->content[bp->bd->L-j-1] = hmm_backward_rec(result, bp->bd->L-j-1, &hmm_hp, bp);
        }
}

static
prob_t hmm_fb_rec(vector_t *forward, vector_t *backward, size_t j, prob_t (*f)(int, int, binProblem*), binProblem* bp)
{
        size_t k;
        prob_t tmp;

        tmp = hmm_forward_rec(forward, j, bp->bd->L-1, f, bp);
        for (k = j+1; k < bp->bd->L; k++) {
                tmp = logadd(tmp, log(1-bp->bd->rho) + hmm_forward_rec(forward, j, k-1, f, bp) + backward->content[k]);
        }
        return tmp;
}

void hmm_fb(vector_t *result, vector_t *forward, vector_t *backward, prob_t (*f)(int, int, binProblem*), binProblem* bp)
{
        size_t j;

        for (j = 0; j < bp->bd->L; j++) {
                result->content[j] = hmm_fb_rec(forward, backward, j, f, bp) - forward->content[bp->bd->L-1];
        }
}

/******************************************************************************
 * Model init
 ******************************************************************************/

void __init_model__() {
#ifdef HAVE_LIB_PTHREAD
        if (pthread_rwlock_init(&lngamma_map_lock, NULL) != 0) {
                fprintf(stderr, "Can't create lngamma_map_lock\n");
                exit(EXIT_FAILURE);
        }
#endif /* HAVE_LIB_PTHREAD */
}

void __free_model__() {
        lngamma_map_free();
}
