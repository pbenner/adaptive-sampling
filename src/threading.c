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

#include <adaptive-sampling/exception.h>

#include <datatypes.h>
#include <threading.h>
#include <utility.h>

#include <limits.h>
#include <pthread.h>

void threaded_computation(
        void *result,
        prob_t evidence_ref,
        binData *bd,
        void *(*f_thread)(void*),
        const char *msg)
{
        size_t i, j, rc;
        binProblem bp[bd->options->threads];
        pthread_t threads[bd->options->threads];
        pthread_data_t data[bd->options->threads];
        pthread_attr_t attr;
        pthread_attr_init(&attr);

        for (i = 0; i < bd->options->threads && i < bd->L; i++) {
                binProblemInit(&bp[i], bd);
                data[i].bp = &bp[i];
                data[i].result = result;
                data[i].evidence_ref = evidence_ref;
        }
        if (bd->options->stacksize < PTHREAD_STACK_MIN) {
                if (pthread_attr_setstacksize (&attr, PTHREAD_STACK_MIN) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }
        else {
                if (pthread_attr_setstacksize (&attr, (size_t)bd->options->stacksize) != 0) {
                        std_warn(NONE, "Couldn't set stack size.");
                }
        }
        for (i = 0; i < bd->L; i += bd->options->threads) {
                for (j = 0; j < bd->options->threads && i+j < bd->L; j++) {
                        notice(NONE, msg, (float)100*(i+j+1)/bd->L);
                        data[j].i = i+j;
                        rc = pthread_create(&threads[j], &attr, f_thread, (void *)&data[j]);
                        if (rc) {
                                std_err(NONE, "Couldn't create thread.");
                        }
                }
                for (j = 0; j < bd->options->threads && i+j < bd->L; j++) {
                        rc = pthread_join(threads[j], NULL);
                        if (rc) {
                                std_err(NONE, "Couldn't join thread.");
                        }
                }
        }
        for (j = 0; j < bd->options->threads && j < bd->L; j++) {
                binProblemFree(&bp[j]);
        }
}
