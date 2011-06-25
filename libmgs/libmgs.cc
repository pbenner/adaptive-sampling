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

#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "libmgs.hh"

Multibin::Multibin(size_t L)
{
        size_t n = (L-2)%32 ? (L-2)/32+1 : (L-2)/32;
        size_t i;

        this->n = n;
        this->n_breaks = L-1;
        this->n_bins   = 1;
        this->breaks   = (uint32_t*)malloc(n*sizeof(uint32_t));

        for (i = 0; i < n; i++) {
                this->breaks[i] = 0;
        }
}

Multibin::Multibin(Multibin* mb, uint32_t* breaks)
{
        size_t L = mb->get_n_breaks()+1;
        size_t n = (L-2)%32 ? (L-2)/32+1 : (L-2)/32;
        size_t i;

        this->n = n;
        this->n_breaks = L-1;
        this->breaks   = (uint32_t*)malloc(n*sizeof(uint32_t));
        this->n_bins   = mb->get_n_bins();

        for (i = 0; i < n; i++) {
                this->breaks[i] = breaks[i];
        }
}

Multibin::~Multibin()
{
        free(this->breaks);
}

Multibin*
Multibin::copy()
{
        return new Multibin(this, this->breaks);
}

void
Multibin::insert_break(size_t i)
{
        if (i < this->n_breaks) {
                size_t n = i/32;
                uint32_t m = 0x1 << (i%32);

                if ((this->breaks[n] & m) == 0) {
                        this->breaks[n] += m;
                        this->n_bins++;
                }
        }
}

void
Multibin::remove_break(size_t i)
{
        if (i < this->n_breaks) {
                size_t n = i/32;
                uint32_t m = 0x1 << (i%32);

                if (this->breaks[n] & m) {
                        this->breaks[n] -= m;
                        this->n_bins--;
                }
        }
}

void
Multibin::switch_break(size_t i)
{
        if (i < this->n_breaks) {
                size_t n = i/32;
                uint32_t m = 0x1 << (i%32);

                if (this->breaks[n] & m) {
                        this->breaks[n] -= m;
                        this->n_bins--;
                }
                else {
                        this->breaks[n] += m;
                        this->n_bins++;
                }
        }
}

void
Multibin::print()
{
        size_t i, j;

        for (i = 0; i < this->n; i++) {
                uint32_t m = 0x1;
                for (j = 0; j < 32 && j*(i+1) < this->n_breaks; j++) {
                        if (this->breaks[i] & m) {
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

list<bin_t>*
Multibin::get_bins()
{
        long int from =  0;
        long int to   = -1;
        size_t i, j;
        list<bin_t>* bins = new list<bin_t>();

        for (i = 0, j = 0; i < this->n; i++) {
                uint32_t p = this->breaks[i];
                uint32_t m;
                for (; j < this->n_bins-1 && p; j++) {
                        m    = p & (~p + 1);
                        from = to+1;
                        to   = 32*i+(size_t)log2(m);
                        bin_t b = {from, to};
                        bins->push_back(b);
                        p -= m;
                }
        }
        from = to+1;
        to   = n_breaks;
        bin_t b = {from, to};
        bins->push_back(b);

        return bins;
}

size_t
Multibin::get_n_breaks()
{
        return this->n_breaks;
}

size_t
Multibin::get_n_bins()
{
        return this->n_bins;
}

void
Multibin::print_bins()
{
        list<bin_t>* bins = this->get_bins();

        for (list<bin_t>::iterator it = bins->begin(); it != bins->end(); it++) {
                printf("(%lu, %lu)\n", (unsigned long)(*it).from, (unsigned long)(*it).to);
        }

        delete(bins);
}
