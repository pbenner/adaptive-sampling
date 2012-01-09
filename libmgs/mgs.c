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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <mgs.h>

multibin_t*
new_multibin(size_t L)
{
        size_t i;
        multibin_t* mb = (multibin_t*)malloc(sizeof(multibin_t));

        mb->n          = (L-2)%32 ? (L-2)/32+1 : (L-2)/32;
        mb->n_breaks   = L-1;
        mb->n_bins     = 1;
        mb->breaks     = (uint32_t*)malloc(mb->n*sizeof(uint32_t));

        for (i = 0; i < mb->n; i++) {
                mb->breaks[i] = 0;
        }
        return mb;
}

multibin_t*
clone_multibin(multibin_t* multibin)
{
        size_t i;
        multibin_t* mb = (multibin_t*)malloc(sizeof(multibin_t));

        mb->n          = multibin->n;
        mb->n_breaks   = multibin->n_breaks;
        mb->n_bins     = multibin->n_bins;
        mb->breaks     = (uint32_t*)malloc(mb->n*sizeof(uint32_t));

        for (i = 0; i < mb->n; i++) {
                mb->breaks[i] = multibin->breaks[i];
        }
        return mb;
}

void
free_multibin(multibin_t* multibin)
{
        free(multibin->breaks);
}

void
insert_break(multibin_t* multibin, size_t i)
{
        if (i < multibin->n_breaks) {
                size_t n = i/32;
                uint32_t m = 0x1 << (i%32);

                if ((multibin->breaks[n] & m) == 0) {
                        multibin->breaks[n] += m;
                        multibin->n_bins++;
                }
        }
}

void
remove_break(multibin_t* multibin, size_t i)
{
        if (i < multibin->n_breaks) {
                size_t n = i/32;
                uint32_t m = 0x1 << (i%32);

                if (multibin->breaks[n] & m) {
                        multibin->breaks[n] -= m;
                        multibin->n_bins--;
                }
        }
}

void
switch_break(multibin_t* multibin, size_t i)
{
        if (i < multibin->n_breaks) {
                size_t n = i/32;
                uint32_t m = 0x1 << (i%32);

                if (multibin->breaks[n] & m) {
                        multibin->breaks[n] -= m;
                        multibin->n_bins--;
                }
                else {
                        multibin->breaks[n] += m;
                        multibin->n_bins++;
                }
        }
}

void
print_mutlibin(multibin_t* multibin)
{
        size_t i, j;

        for (i = 0; i < multibin->n; i++) {
                uint32_t m = 0x1;
                for (j = 0; j < 32 && j*(i+1) < multibin->n_breaks; j++) {
                        if (multibin->breaks[i] & m) {
                                printf("1 ");
                        }
                        else {
                                printf("0 ");
                        }
                        m = m << 1;
                }
        }
        printf("\n\n");
}

void
get_bins(multibin_t* multibin, bin_t* bins)
{
        long int from  =  0;
        long int to    = -1;
        size_t i, j, k = 0;

        for (i = 0, j = 0; i < multibin->n; i++) {
                uint32_t p = multibin->breaks[i];
                uint32_t m;
                for (; j < multibin->n_bins-1 && p; j++) {
                        m    = p & (~p + 1);
                        from = to+1;
                        to   = 32*i+(size_t)log2(m);
                        bins[k].from = from;
                        bins[k].to   = to;
                        p -= m; k++;
                }
        }
        bins[k].from = to+1;
        bins[k].to   = multibin->n_breaks;
}

void
get_breaks(multibin_t* multibin, size_t *breaks)
{
        size_t i;
        bin_t bins[multibin->n_bins];

        get_bins(multibin, bins);
        for (i = 0; i < multibin->n_bins; i++) {
                breaks[bins[i].from]++;
        }
}
